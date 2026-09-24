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
//! The draw is reproduced, not refused. `nextGaussian` goes through `StrictMath.log`, which is
//! FDLIBM and differs from the correctly rounded `Math.log` on 186 of the jmath corpus's 44,996
//! points; with `jmath::strict_math::log` (htsjdk-rs decision 0044) it is exact, and
//! [`gatk_engine::java_random::JavaRandom::next_gaussian`] is the polar method on top of it. What
//! the caller owns is the generator: `Utils.getRandomGenerator()` is one static stream per JVM, so
//! the value a site gets depends on how many Gaussians every earlier site drew. [`qual_by_depth`]
//! therefore takes the generator rather than making one, and [`fix_too_high_qd`] advances it only
//! when the ratio reaches 35.
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
use gatk_engine::java_random::JavaRandom;
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

/// `QualByDepth.IDEAL_HIGH_QD`.
const IDEAL_HIGH_QD: f64 = 30.0;

/// `QualByDepth.JITTER_SIGMA`.
const JITTER_SIGMA: f64 = 3.0;

/// `QualByDepth.fixTooHighQD`: a ratio at or above 35 becomes 30 plus a Gaussian with a standard
/// deviation of 3, drawn from the run's one generator. Below 35 the generator is not touched.
pub fn fix_too_high_qd(qd: f64, random: &mut JavaRandom) -> f64 {
    if qd < MAX_QD_BEFORE_FIXING {
        qd
    } else {
        IDEAL_HIGH_QD + random.next_gaussian() * JITTER_SIGMA
    }
}

/// Whether a genotype is het or hom-var, which is the only kind `QD` counts.
///
/// htsjdk's own `getType()`, which a haploid call answers too: a lone alternate is `HOM_VAR`, so
/// a haploid variant sample's depth counts, and a mixed call such as `./1` is neither.
fn is_het_or_hom_var(genotype: &Genotype) -> bool {
    genotype.is_het() || genotype.is_hom_var()
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

/// `QualByDepth.annotate`, drawing from `random` when the ratio reaches 35.
pub fn qual_by_depth(
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    raw_qual_approx: Option<i32>,
    random: &mut JavaRandom,
) -> Option<String> {
    // `vc.hasLog10PError()` is false for a QUAL of `.`, which htsjdk stores as
    // `NO_LOG10_PERROR`.
    let has_log10_perror = vc.log10_p_error != 1.0;
    if !has_log10_perror && raw_qual_approx.is_none() {
        return None;
    }
    if vc.genotypes.is_empty() {
        return None;
    }
    let depth = qual_by_depth_depth(vc, likelihoods);
    if depth == 0 {
        return None;
    }
    let qual = if has_log10_perror {
        -10.0 * vc.log10_p_error
    } else {
        raw_qual_approx.unwrap_or(0) as f64
    };
    let qd = qual / depth as f64;
    Some(format_two_decimals(fix_too_high_qd(qd, random)))
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
        // `vc.getNoCallCount()` counts **genotypes** that are no-call, not the no-call alleles
        // inside them, so a diploid no-call sample contributes one rather than two.
        let no_calls = vc
            .genotypes
            .iter()
            .filter(|g| g.alleles.iter().all(|a| a.is_no_call()))
            .count();
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

    fn site(qual: f64, ad: [i32; 2]) -> VariantContext {
        let alleles = vec![
            htsjdk_vcf::allele::Allele::from_str("A", true).unwrap(),
            htsjdk_vcf::allele::Allele::from_str("C", false).unwrap(),
        ];
        let mut vc = VariantContext::new("chr1", 100, alleles.clone());
        vc.stop = 100;
        vc.log10_p_error = qual / -10.0;
        let mut genotype = Genotype::new("s1", alleles);
        genotype.ad = Some(ad.to_vec());
        vc.genotypes.push(genotype);
        vc
    }

    #[test]
    fn a_ratio_past_35_is_replaced_by_a_draw_and_one_below_draws_nothing() {
        let mut random = JavaRandom::gatk();
        // 300 / 10 = 30: written as is, and the generator is not touched.
        assert_eq!(
            qual_by_depth(&site(300.0, [5, 5]), None, None, &mut random).as_deref(),
            Some("30.00")
        );
        assert_eq!(random, JavaRandom::gatk());
        // 1000 / 10 = 100: replaced by 30 plus three times the stream's first Gaussian.
        let expected = 30.0 + JavaRandom::gatk().next_gaussian() * 3.0;
        assert_eq!(
            qual_by_depth(&site(1000.0, [5, 5]), None, None, &mut random),
            Some(format_two_decimals(expected))
        );
        assert_ne!(random, JavaRandom::gatk());
    }
}
