//! `QualByDepth`, `GenotypeSummaries` and `LikelihoodRankSumTest`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.6.2.0).
//!
//! `QD`, `NCC`/`GQ_MEAN`/`GQ_STDDEV`, and `LikelihoodRankSum`: three annotations that read the
//! genotypes rather than the reads.
//!
//! # `QD` above 35 is **randomised**, and therefore not reproducible
//!
//! ```java
//! public static double fixTooHighQD(final double QD) {
//!     if ( QD < MAX_QD_BEFORE_FIXING ) { return QD; }
//!     return IDEAL_HIGH_QD + Utils.getRandomGenerator().nextGaussian() * JITTER_SIGMA;
//! }
//! ```
//!
//! *"The haplotype caller generates very high quality scores when multiple events are on the same
//! haplotype... VQSR will filter these out"*, so a `QD` at or above 35 is replaced by 30 plus a
//! Gaussian jitter. The value written to the VCF is then a draw from a random generator, and two
//! runs of the same tool on the same data agree only because the generator is seeded.
//!
//! This port **refuses** that branch rather than approximating it: `nextGaussian` goes through
//! `StrictMath.log`, which is fdlibm and not the correctly-rounded logarithm this crate has.
//! Measured on the jmath corpus, `Math.log` and `StrictMath.log` differ on 186 of 44,996 points,
//! so the ported logarithm is not a substitute. [`QualByDepthError::RandomisedAboveThreshold`] is
//! that refusal, and the suite measures where the boundary is rather than what is past it.
//!
//! # The depth `QD` divides by is not `DP`
//!
//! Only het and hom-var genotypes count. For each, the whole `AD` total is added, but only if the
//! **alternate** part of it is greater than one does the same total also go into an
//! "AD-restricted" tally, and if that tally is non-zero at the end it replaces the depth
//! entirely. So one sample with two alternate reads can discard the depth of every other sample.
//!
//! # `GQ_MEAN` and `GQ_STDDEV` are strings, and the standard deviation needs two genotypes
//!
//! Both are `String.format("%.2f", ...)`, and `GQ_STDDEV` is written only when more than one
//! genotype has a `GQ`, because a `DescriptiveStatistics` of one value has a standard deviation
//! of zero and the reference declines to report it.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::{Genotype, VariantContext};

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `GATKVCFConstants.QUAL_BY_DEPTH_KEY`.
pub const QUAL_BY_DEPTH_KEY: &str = "QD";
/// `GATKVCFConstants.NOCALL_CHROM_KEY`.
pub const NOCALL_CHROM_KEY: &str = "NCC";
/// `GATKVCFConstants.GQ_MEAN_KEY`.
pub const GQ_MEAN_KEY: &str = "GQ_MEAN";
/// `GATKVCFConstants.GQ_STDEV_KEY`.
pub const GQ_STDEV_KEY: &str = "GQ_STDDEV";
/// `GATKVCFConstants.LIKELIHOOD_RANK_SUM_KEY`.
pub const LIKELIHOOD_RANK_SUM_KEY: &str = "LikelihoodRankSum";

/// `QualByDepth.MAX_QD_BEFORE_FIXING`.
const MAX_QD_BEFORE_FIXING: f64 = 35.0;

/// What this port refuses rather than inventing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum QualByDepthError {
    /// `QD >= 35`, where the reference replaces the value with a random draw. See the module note.
    RandomisedAboveThreshold { raw: f64 },
}

/// Whether a genotype is het or hom-var, which is the only kind `QD` counts.
fn is_het_or_hom_var(genotype: &Genotype) -> bool {
    let called: Vec<&htsjdk_vcf::allele::Allele> = genotype
        .alleles
        .iter()
        .filter(|a| !a.is_no_call())
        .collect();
    if called.len() < 2 {
        return false;
    }
    let all_ref = called.iter().all(|a| a.is_reference());
    !all_ref
}

/// `QualByDepth.getDepth`.
pub fn qual_by_depth_depth(
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> i32 {
    let mut depth = 0i32;
    let mut ad_restricted_depth = 0i32;
    for genotype in &vc.genotypes {
        if !is_het_or_hom_var(genotype) {
            continue;
        }
        if let Some(ad) = &genotype.ad {
            let total: i32 = ad.iter().sum();
            if total != 0 {
                // The whole total goes into both tallies, but only the second is conditional.
                if total - ad[0] > 1 {
                    ad_restricted_depth += total;
                }
                depth += total;
                continue;
            }
        }
        if let Some(likelihoods) = likelihoods {
            if let Some(index) = likelihoods.index_of_sample(&genotype.sample_name) {
                depth += likelihoods.sample_evidence_count(index) as i32;
            }
        } else if let Some(dp) = genotype.dp {
            depth += dp;
        }
    }
    // One sample with two alternate reads can discard every other sample's depth.
    if ad_restricted_depth > 0 {
        depth = ad_restricted_depth;
    }
    depth
}

/// `QualByDepth.annotate`, with the randomised branch refused.
pub fn qual_by_depth(
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    raw_qual_approx: Option<i32>,
) -> Result<Option<String>, QualByDepthError> {
    // `vc.hasLog10PError()` is false for a QUAL of `.`, which htsjdk stores as
    // `NO_LOG10_PERROR`.
    let has_log10_perror = vc.log10_p_error != 1.0;
    if !has_log10_perror && raw_qual_approx.is_none() {
        return Ok(None);
    }
    if vc.genotypes.is_empty() {
        return Ok(None);
    }
    let depth = qual_by_depth_depth(vc, likelihoods);
    if depth == 0 {
        return Ok(None);
    }
    let qual = if has_log10_perror {
        -10.0 * vc.log10_p_error
    } else {
        raw_qual_approx.unwrap_or(0) as f64
    };
    let qd = qual / depth as f64;
    if qd >= MAX_QD_BEFORE_FIXING {
        return Err(QualByDepthError::RandomisedAboveThreshold { raw: qd });
    }
    Ok(Some(format_two_decimals(qd)))
}

/// `String.format("%.2f", value)`, half-up on the decimal expansion as Java rounds it.
pub fn format_two_decimals(value: f64) -> String {
    crate::rank_sum::format_decimals(value, 2)
}

/// `GenotypeSummaries`: `NCC`, `GQ_MEAN` and `GQ_STDDEV`.
pub struct GenotypeSummaries;

impl InfoFieldAnnotation for GenotypeSummaries {
    fn key_names(&self) -> Vec<&'static str> {
        vec![NOCALL_CHROM_KEY, GQ_MEAN_KEY, GQ_STDEV_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        _likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        if vc.genotypes.is_empty() {
            return Vec::new();
        }
        let mut out: Vec<(String, AnnotationValue)> = Vec::new();
        // `vc.getNoCallCount()`: no-call **alleles**, across every genotype, not no-call samples.
        let no_calls: usize = vc
            .genotypes
            .iter()
            .map(|g| g.alleles.iter().filter(|a| a.is_no_call()).count())
            .sum();
        out.push((
            NOCALL_CHROM_KEY.to_string(),
            AnnotationValue::Int(no_calls as i32),
        ));

        let values: Vec<f64> = vc
            .genotypes
            .iter()
            .filter_map(|g| g.gq.map(|gq| gq as f64))
            .collect();
        if values.is_empty() {
            return out;
        }
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        out.push((
            GQ_MEAN_KEY.to_string(),
            AnnotationValue::Str(format_two_decimals(mean)),
        ));
        if values.len() > 1 {
            // `Variance`, bias-corrected: the sum of squared deviations over n - 1, computed in
            // the two-pass form commons-math3 uses.
            let variance = values
                .iter()
                .map(|value| (value - mean) * (value - mean))
                .sum::<f64>()
                / (values.len() - 1) as f64;
            out.push((
                GQ_STDEV_KEY.to_string(),
                AnnotationValue::Str(format_two_decimals(variance.sqrt())),
            ));
        }
        out
    }
}

/// `LikelihoodRankSumTest`: the rank-sum test over the **likelihoods** themselves.
///
/// Its `getElementForRead` takes the best allele's likelihood, and the two-argument form the other
/// members implement answers empty, so this is the one member of the family whose value comes from
/// the matrix rather than from the read.
pub struct LikelihoodRankSumTest;

impl crate::rank_sum::RankSumTest for LikelihoodRankSumTest {
    fn vcf_key(&self) -> &'static str {
        LIKELIHOOD_RANK_SUM_KEY
    }

    fn element_for_read(&self, _read: &BamRecord, _vc: &VariantContext) -> Option<f64> {
        // The two-argument form, which the reference leaves empty with a comment saying it should
        // perhaps throw. Reached only if a caller bypasses the best-allele form.
        None
    }

    fn element_for_best_allele(&self, likelihood: f64) -> Option<f64> {
        Some(likelihood)
    }

    fn uses_best_allele(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_randomised_branch_is_refused_rather_than_drawn() {
        let mut vc = VariantContext::new(
            "chr1",
            100,
            vec![
                htsjdk_vcf::allele::Allele::from_str("A", true).unwrap(),
                htsjdk_vcf::allele::Allele::from_str("C", false).unwrap(),
            ],
        );
        vc.stop = 100;
        vc.log10_p_error = -100.0;
        let mut genotype = Genotype::new(
            "s1",
            vec![
                htsjdk_vcf::allele::Allele::from_str("A", true).unwrap(),
                htsjdk_vcf::allele::Allele::from_str("C", false).unwrap(),
            ],
        );
        genotype.ad = Some(vec![5, 5]);
        vc.genotypes.push(genotype);
        // 1000 / 10 = 100, well past the threshold.
        assert!(matches!(
            qual_by_depth(&vc, None, None),
            Err(QualByDepthError::RandomisedAboveThreshold { .. })
        ));
    }
}
