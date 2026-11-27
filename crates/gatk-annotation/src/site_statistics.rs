//! `QualByDepth`, `GenotypeSummaries` and `LikelihoodRankSumTest`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.7.0.0; `QualByDepth.java`,
//! `GenotypeSummaries.java` and `LikelihoodRankSumTest.java` are unchanged since 4.6.2.0).
//!
//! `QD`, `NCC`/`GQ_MEAN`/`GQ_STDDEV`, and `LikelihoodRankSum`: three annotations that read the
//! genotypes rather than the reads.
//!
//! # `QD` above 35 is **randomised**, and therefore not reproducible
//!
//! ```java
//! public static double fixTooHighQD(final double QD) {
//!     if ( QD < MAX_QD_BEFORE_FIXING ) {
//!         return QD;
//!     } else {
//!         return IDEAL_HIGH_QD + Utils.getRandomGenerator().nextGaussian() * JITTER_SIGMA;
//!     }
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
//! The golden's `depth two-samples` row (`site_statistics.txt.gz`) is that: `AD`s of `9,1` and
//! `8,2` total twenty, and the depth is 10, the second sample's alone.
//!
//! # `GQ_MEAN` and `GQ_STDDEV` are strings, and the standard deviation needs two genotypes
//!
//! Both are `String.format("%.2f", ...)`, and `GQ_STDDEV` is written only when more than one
//! genotype has a `GQ`, because a `DescriptiveStatistics` of one value has a standard deviation
//! of zero and the reference declines to report it.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a `T` or nothing", checked by the compiler: Java's nullable reference, or an
//!   absent VCF field. `if let Some(x) = y { ... }` runs the block only when `y` holds a value,
//!   and binds it to `x`.
//! - `&mut JavaRandom` is a **writable borrow** of the caller's generator: a draw advances the
//!   caller's own stream, not a copy, which is what one static Java generator per JVM requires.
//! - `.iter().filter(..).count()` and `.iter().map(..).sum()` are Java streams. A closure is
//!   written `|g| expr`, Java's `g -> expr`.
//! - `x as f64` converts an integer to a `double`, as Java's widening does implicitly in
//!   `qual / depth`.

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
///
/// **Why the generator is a parameter**: the reference calls `Utils.getRandomGenerator()`, a
/// static. Taking it here makes the draw count visible to the caller, and the module's unit test
/// checks that a ratio of 30 leaves the generator unchanged and one of 100 advances it.
pub fn fix_too_high_qd(qd: f64, random: &mut JavaRandom) -> f64 {
    if qd < MAX_QD_BEFORE_FIXING {
        qd
    } else {
        IDEAL_HIGH_QD + random.next_gaussian() * JITTER_SIGMA
    }
}

/// Whether a genotype is het or hom-var, which is the only kind `QD` counts:
/// `if ( !genotype.isHet() && !genotype.isHomVar() ) { continue; }` in `getDepth`.
///
/// htsjdk's own `getType()`, which a haploid call answers too: a lone alternate is `HOM_VAR`, so
/// a haploid variant sample's depth counts, and a mixed call such as `./1` is neither.
fn is_het_or_hom_var(genotype: &Genotype) -> bool {
    genotype.is_het() || genotype.is_hom_var()
}

/// `QualByDepth.getDepth`.
///
/// **What**: the depth `QD` divides by, summed over het and hom-var genotypes only.
///
/// **How**, per genotype, in the reference's order: a non-zero `AD` total wins and the sample is
/// done (`continue`); otherwise the likelihood matrix's evidence count for the sample, or, when no
/// matrix was given, `DP`. A genotype with neither adds nothing.
///
/// **Not reproduced**: a sample the matrix does not hold. The reference passes
/// `indexOfSample`'s `-1` to `sampleEvidenceCount`, whose `Utils.validIndex` throws
/// `IllegalArgumentException`; here the sample adds nothing.
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
        // `genotype.hasAD()`. `&genotype.ad` borrows the field rather than moving it out.
        if let Some(ad) = &genotype.ad {
            // `(int) MathUtils.sum(AD)`, which sums in a `long` and narrows.
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
        // "no AD value or it is a dummy value": fall back to the matrix, then to `DP`.
        if let Some(likelihoods) = likelihoods {
            // `index_of_sample` answers `None` for an unknown sample where Java answers -1.
            // `as i32` narrows the unsigned count to Java's `int`.
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
///
/// **What**: `QD`, the formatted `QUAL / depth`, or `None` where the reference returns an empty
/// map: no `QUAL` and no `RAW_QUAL_APPROX`, no genotypes, or a depth of zero.
///
/// `raw_qual_approx` is the caller's reading of `vc.getAttributeAsInt(RAW_QUAL_APPROX_KEY, 0)`,
/// `None` when the attribute is absent. The reference's `NumberFormatException` branch, which
/// rethrows as a `GATKException`, is the caller's to model.
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
    // `vc.hasLog10PError()` decides again, so a site with both a `QUAL` and `RAW_QUAL_APPROX`
    // uses the `QUAL`.
    let qual = if has_log10_perror {
        -10.0 * vc.log10_p_error
    } else {
        raw_qual_approx.unwrap_or(0) as f64
    };
    // `qual / depth`: a `double` over an `int`, which Java widens to `double` before dividing.
    let qd = qual / depth as f64;
    Some(format_two_decimals(fix_too_high_qd(qd, random)))
}

/// `String.format("%.2f", value)`, half-up on the decimal expansion as Java rounds it.
///
/// Delegates to `crate::rank_sum::format_decimals`, which the rank-sum annotations also use.
pub fn format_two_decimals(value: f64) -> String {
    crate::rank_sum::format_decimals(value, 2)
}

/// `GenotypeSummaries`: `NCC`, `GQ_MEAN` and `GQ_STDDEV`.
///
/// `NCC` is written whenever there are genotypes, `GQ_MEAN` when at least one has a `GQ`, and
/// `GQ_STDDEV` when at least two do, which is the reference's `stats.getN() > 0L` and
/// `stats.getN() > 1L`. The golden's `summaries no-gq` and `summaries one-gq` rows
/// (`site_statistics.txt.gz`) show the first two cut-offs.
pub struct GenotypeSummaries;

impl InfoFieldAnnotation for GenotypeSummaries {
    fn key_names(&self) -> Vec<&'static str> {
        vec![NOCALL_CHROM_KEY, GQ_MEAN_KEY, GQ_STDEV_KEY]
    }

    /// `GenotypeSummaries.annotate`. The returned pairs are in the order of the reference's
    /// `LinkedHashMap` puts.
    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        _likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        // `if ( ! vc.hasGenotypes() )`.
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

        // The `GQ`s of the genotypes that have one, in genotype order: the `stats.addValue` loop.
        // `filter_map` keeps only the entries the closure answers `Some` for.
        let values: Vec<f64> = vc
            .genotypes
            .iter()
            .filter_map(|g| g.gq.map(|gq| gq as f64))
            .collect();
        if values.is_empty() {
            return out;
        }
        // `stats.getMean()`. `sum::<f64>()` names the type the sum is accumulated in.
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
/// the matrix rather than from the read:
///
/// ```java
/// return OptionalDouble.of(bestAllele.likelihood);
/// ```
///
/// The reference also throws `IllegalStateException` for a non-informative best allele. That
/// check is not here: `crate::rank_sum` skips non-informative best alleles before asking.
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
