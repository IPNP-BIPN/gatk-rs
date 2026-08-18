//! `GermlineFilter.germlineProbability`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.GermlineFilter` (GATK 4.6.2.0).
//!
//! The probability that an allele Mutect called somatic is really germline, from five numbers: the
//! normal's log odds, the two log odds of germline against somatic, the population allele frequency
//! and the log prior that the site is somatic.
//!
//! # The population frequency is three priors at once
//!
//! ```java
//! final double logPriorGermlineHet = Math.log(2*populationAF*(1-populationAF));
//! final double logPriorGermlineHomAlt = Math.log( MathUtils.square(populationAF));
//! final double logPriorNotGermline = Math.log(MathUtils.square(1 - populationAF));
//! ```
//!
//! One number decides all three, so a frequency of zero makes both germline hypotheses impossible
//! and the answer `0.0`, while a frequency of one makes the somatic hypothesis impossible and the
//! answer `1.0`. Neither is a special case in the code: both fall out of `log(0)`.
//!
//! # The answer is the first entry, not the second
//!
//! ```java
//! return NaturalLogUtils.normalizeLog(new double[] {logProbGermline, logProbSomatic}, false, true)[0];
//! ```
//!
//! [`crate::mutect_engine::posterior_probability_of_error`] has the same shape and returns index
//! **1**. These two functions answer opposite questions, and taking the wrong index turns a germline
//! filter into a somatic one that never fires.
//!
//! # The somatic prior enters twice
//!
//! Once as itself on the somatic side, and once through `log1mexp` on **both** germline sides. A
//! prior of one therefore answers `0.0` whatever the odds say, and a prior of negative infinity
//! answers `1.0`.

use crate::allele_filter::{sum_ads_over_samples, AlleleDepthTooShort, GenotypeData};
use crate::math_utils::{max_element_index, pow10};
use crate::natural_log_utils::{
    log1mexp, log_sum_exp, normalize_from_log_to_linear_space, NonFiniteSum,
};
use crate::somatic_clustering_model::{indel_length, AlternateAllele, SomaticClusteringModel};

/// `germlineProbability`.
///
/// `log_odds_of_germline_hom_alt_vs_somatic` is negative infinity when the caller judged the allele
/// fraction too low for a germline hom alt: the hypothesis is switched off by a value rather than by
/// a flag.
pub fn germline_probability(
    normal_log_odds: f64,
    log_odds_of_germline_het_vs_somatic: f64,
    log_odds_of_germline_hom_alt_vs_somatic: f64,
    population_af: f64,
    log_prior_somatic: f64,
) -> Result<f64, NonFiniteSum> {
    let log_prior_not_somatic = log1mexp(log_prior_somatic);
    let log_prior_germline_het = (2.0 * population_af * (1.0 - population_af)).ln();
    let log_prior_germline_hom_alt = (population_af * population_af).ln();
    let log_prior_not_germline = ((1.0 - population_af) * (1.0 - population_af)).ln();

    // Unnormalized, and the normal's odds reach both germline hypotheses and neither somatic one.
    let log_prob_germline_het = log_prior_germline_het
        + log_odds_of_germline_het_vs_somatic
        + normal_log_odds
        + log_prior_not_somatic;
    let log_prob_germline_hom_alt = log_prior_germline_hom_alt
        + log_odds_of_germline_hom_alt_vs_somatic
        + normal_log_odds
        + log_prior_not_somatic;
    let log_prob_germline = log_sum_exp(&[log_prob_germline_het, log_prob_germline_hom_alt])?;
    let log_prob_somatic = log_prior_not_germline + log_prior_somatic;

    let normalized = normalize_from_log_to_linear_space(&[log_prob_germline, log_prob_somatic])?;
    // The FIRST entry: germline.
    Ok(normalized[0])
}

/// `GermlineFilter`'s identity.
pub const FILTER_NAME: &str = "germline";

/// `phredScaledPosteriorAnnotationName`, `GERMLINE_QUAL_KEY`.
pub const ANNOTATION: &str = "GERMQ";

/// `GermlineFilter.EPSILON`, which brackets the population frequency.
pub const EPSILON: f64 = 1.0e-10;

/// `MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT`, below which the hom-alt hypothesis is switched off by
/// a value rather than by a flag.
pub const MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT: f64 = 0.9;

/// `NaturalLogUtils.LOG_ONE_HALF`.
fn log_one_half() -> f64 {
    jmath::math::log(0.5)
}

/// What the wrapper refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum GermlineError {
    /// A genotype's `AD` shorter than the record's allele count.
    AlleleDepth(AlleleDepthTooShort),
    /// `MathUtils.addToArrayInPlace`: the arrays must have the same length, and a genotype with no
    /// `AF` supplies a one-element default whatever the record's allele count is.
    ArrayLengthsDiffer,
    /// `logSumExp` or the normalisation over a sum that is not finite.
    NonFiniteSum(NonFiniteSum),
}

impl GermlineError {
    pub fn class(&self) -> &'static str {
        match self {
            GermlineError::AlleleDepth(_) => "java.lang.ArrayIndexOutOfBoundsException",
            GermlineError::ArrayLengthsDiffer => "java.lang.IllegalArgumentException",
            GermlineError::NonFiniteSum(_) => "java.lang.IllegalArgumentException",
        }
    }
}

impl From<AlleleDepthTooShort> for GermlineError {
    fn from(error: AlleleDepthTooShort) -> Self {
        GermlineError::AlleleDepth(error)
    }
}

impl From<NonFiniteSum> for GermlineError {
    fn from(error: NonFiniteSum) -> Self {
        GermlineError::NonFiniteSum(error)
    }
}

/// `Mutect2FilteringEngine.weightedAverageOfTumorAFs`.
///
/// The weights are the tumour genotypes' total depths, and the division by their sum comes last, so
/// a record with no tumour depth divides by zero rather than refusing.
pub fn weighted_average_of_tumor_afs<T>(
    genotypes: &[GenotypeData<T>],
    allele_fractions: &[Vec<f64>],
    alternate_count: usize,
) -> Result<Vec<f64>, GermlineError> {
    let mut total_weight = 0.0;
    let mut averages = vec![0.0; alternate_count];
    for (index, genotype) in genotypes.iter().enumerate() {
        if !genotype.tumor {
            continue;
        }
        let weight: f64 = genotype.allele_depths.iter().map(|d| f64::from(*d)).sum();
        total_weight += weight;
        // `getAttributeAsDoubleArray(g, AF, () -> new double[] {0.0}, 0.0)`: a genotype with no `AF`
        // supplies ONE zero, whatever the record's allele count is.
        let sample: Vec<f64> = if allele_fractions[index].is_empty() {
            vec![0.0]
        } else {
            allele_fractions[index].clone()
        };
        // `MathArrays.scaleInPlace` then `MathUtils.addToArrayInPlace`, which validates the lengths.
        if sample.len() != averages.len() {
            return Err(GermlineError::ArrayLengthsDiffer);
        }
        for (average, fraction) in averages.iter_mut().zip(&sample) {
            *average += weight * fraction;
        }
    }
    for average in averages.iter_mut() {
        *average *= 1.0 / total_weight;
    }
    Ok(averages)
}

/// `GermlineFilter.computeMinorAlleleFraction`.
///
/// With no tumour segmentation table every sample's minor allele fraction is `0.5`, so this is a
/// depth-weighted average of one repeated constant. `minor_allele_fractions` is one value per
/// genotype, `0.5` where no segment overlaps the record; the denominator is the tumour-only depth
/// sum the caller already computed.
pub fn compute_minor_allele_fraction<T>(
    genotypes: &[GenotypeData<T>],
    minor_allele_fractions: &[f64],
    allele_counts: &[i32],
) -> f64 {
    let mut weighted_sum = 0.0;
    for (index, genotype) in genotypes.iter().enumerate() {
        if !genotype.tumor {
            continue;
        }
        let depth: f64 = genotype.allele_depths.iter().map(|d| f64::from(*d)).sum();
        weighted_sum += minor_allele_fractions[index] * depth;
    }
    let total: f64 = allele_counts.iter().map(|d| f64::from(*d)).sum();
    weighted_sum / total
}

/// `BinomialDistribution.logProbability(x)`.
fn binomial_log_probability(trials: i32, x: i32, p: f64) -> f64 {
    if trials == 0 {
        return if x == 0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if x < 0 || x > trials {
        return f64::NEG_INFINITY;
    }
    jmath::saddle_point::log_binomial_probability(x, trials, p, 1.0 - p)
}

/// `GermlineFilter.calculateErrorProbability`, with `Mutect2VariantFilter`'s copy around it.
///
/// `tumor_log_10_odds` is `TLOD` and `population_negative_log10_af` is `POPAF`, both `None` when the
/// annotation is absent, which the required-annotation check answers with `0.0` per allele.
/// `normal_log_10_odds` is `NLOD`, absent meaning zero rather than a skip.
#[allow(clippy::too_many_arguments)]
pub fn germline_error_probabilities<T>(
    model: &mut SomaticClusteringModel,
    tumor_log_10_odds: Option<&[f64]>,
    population_negative_log10_af: Option<&[f64]>,
    normal_log_10_odds: Option<&[f64]>,
    genotypes: &[GenotypeData<T>],
    allele_fractions: &[Vec<f64>],
    minor_allele_fractions: &[f64],
    alternates: &[AlternateAllele],
    reference_length: i32,
) -> Result<Vec<f64>, GermlineError> {
    let alternate_count = alternates.len();
    let (Some(tumor_log_10_odds), Some(population_negative_log10_af)) =
        (tumor_log_10_odds, population_negative_log10_af)
    else {
        return Ok(vec![0.0; alternate_count]);
    };

    // `getTumorLogOdds` converts to natural log before the maximum is taken.
    let somatic_log_odds: Vec<f64> = tumor_log_10_odds
        .iter()
        .map(|value| crate::allele_likelihoods::log10_to_log(*value))
        .collect();
    let max_lod_index = max_element_index(&somatic_log_odds, 0, somatic_log_odds.len());

    // `Math.pow(10, -POPAF[maxLodIndex])`, and the two brackets around it.
    let population_af = pow10(-population_negative_log10_af[max_lod_index]);
    if population_af < EPSILON {
        return Ok(vec![0.0; alternate_count]);
    } else if population_af > 1.0 - EPSILON {
        return Ok(vec![1.0; alternate_count]);
    }

    let allele_counts = sum_ads_over_samples(alternate_count + 1, genotypes, true, false)?;
    let total_count: i32 = allele_counts.iter().sum();
    if total_count == 0 {
        return Ok(vec![0.0; alternate_count]);
    }
    // The depth array carries the reference, so the alternate needs `+ 1`; the weighted fractions do
    // not, so they take the index as it is.
    let alt_count = allele_counts[max_lod_index + 1];
    let alt_allele_fraction =
        weighted_average_of_tumor_afs(genotypes, allele_fractions, alternate_count)?[max_lod_index];

    let maf = compute_minor_allele_fraction(genotypes, minor_allele_fractions, &allele_counts);

    // Alt minor and alt major, added and halved in log space.
    let log_germline_likelihood = log_one_half()
        + log_sum_exp(&[
            binomial_log_probability(total_count, alt_count, maf),
            binomial_log_probability(total_count, alt_count, 1.0 - maf),
        ])?;
    let log_somatic_likelihood = model.log_likelihood_given_somatic(total_count, alt_count)?;
    let log_odds_of_germline_het_vs_somatic = log_germline_likelihood - log_somatic_likelihood;
    // Switched off by a value, not by a flag.
    let log_odds_of_germline_hom_alt_vs_somatic =
        if alt_allele_fraction < MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT {
            f64::NEG_INFINITY
        } else {
            0.0
        };

    // `NLOD` is log10 and is converted; absent, it is zero rather than a skip.
    let normal_lod = match normal_log_10_odds {
        Some(odds) => crate::allele_likelihoods::log10_to_log(odds[max_lod_index]),
        None => 0.0,
    };

    let prior = model
        .log_prior_of_somatic_variant(indel_length(reference_length, alternates[max_lod_index]));
    // The minus sign is the reference's: `NLOD` is the odds of the allele NOT being in the normal.
    let probability = germline_probability(
        -normal_lod,
        log_odds_of_germline_het_vs_somatic,
        log_odds_of_germline_hom_alt_vs_somatic,
        population_af,
        prior,
    )?;
    Ok(vec![
        crate::mutect_engine::round_finite_precision_errors(
            probability
        );
        alternate_count
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    const NORMAL_LOG_ODDS: f64 = 5.0;
    const HET_VS_SOMATIC: f64 = 1.0;
    const HOM_ALT_VS_SOMATIC: f64 = 0.0;
    const POPULATION_AF: f64 = 0.001;
    const LOG_PRIOR_SOMATIC: f64 = -13.0;

    fn baseline(af: f64, prior: f64) -> f64 {
        germline_probability(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            HOM_ALT_VS_SOMATIC,
            af,
            prior,
        )
        .expect("finite")
    }

    #[test]
    fn the_population_frequency_decides_both_ends() {
        // Nothing can be germline.
        assert_eq!(baseline(0.0, LOG_PRIOR_SOMATIC), 0.0);
        // Nothing can be somatic.
        assert_eq!(baseline(1.0, LOG_PRIOR_SOMATIC), 1.0);
        // And in between it rises with the frequency.
        assert!(baseline(1.0e-8, LOG_PRIOR_SOMATIC) < baseline(1.0e-4, LOG_PRIOR_SOMATIC));
        assert!(baseline(1.0e-4, LOG_PRIOR_SOMATIC) < baseline(0.5, LOG_PRIOR_SOMATIC));
    }

    #[test]
    fn the_somatic_prior_enters_twice() {
        // A prior of one: certainly somatic, whatever the odds.
        assert_eq!(baseline(POPULATION_AF, 0.0), 0.0);
        // A prior of zero: certainly not.
        assert_eq!(baseline(POPULATION_AF, f64::NEG_INFINITY), 1.0);
        // And it moves the answer everywhere in between.
        assert!(baseline(POPULATION_AF, -1.0) < baseline(POPULATION_AF, -13.0));
    }

    #[test]
    fn the_normal_odds_move_it_monotonically_and_an_infinity_is_refused() {
        let mut previous = 0.0;
        for odds in [-50.0, -10.0, -1.0, 0.0, 1.0, 10.0, 50.0] {
            let value = germline_probability(
                odds,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                POPULATION_AF,
                LOG_PRIOR_SOMATIC,
            )
            .expect("finite");
            assert!(value >= previous, "{odds}");
            previous = value;
        }
        // Negative infinity is simply zero.
        assert_eq!(
            germline_probability(
                f64::NEG_INFINITY,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                POPULATION_AF,
                LOG_PRIOR_SOMATIC
            )
            .expect("finite"),
            0.0
        );
        // Positive infinity is a refusal from the normaliser.
        let error = germline_probability(
            f64::INFINITY,
            HET_VS_SOMATIC,
            HOM_ALT_VS_SOMATIC,
            POPULATION_AF,
            LOG_PRIOR_SOMATIC,
        )
        .expect_err("not finite");
        assert_eq!(
            error.message(),
            "logValues must be non-infinite and non-NAN"
        );
    }

    #[test]
    fn the_hom_alt_hypothesis_is_switched_off_by_a_value() {
        let on = baseline(POPULATION_AF, LOG_PRIOR_SOMATIC);
        let off = germline_probability(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            f64::NEG_INFINITY,
            POPULATION_AF,
            LOG_PRIOR_SOMATIC,
        )
        .expect("finite");
        // At a rare allele the hom alt hypothesis is worth almost nothing.
        assert!(off < on);
        assert!(on - off < 1.0e-6);
    }

    #[test]
    fn both_corners_normalise_rather_than_refusing() {
        // No germline hypothesis and a certain somatic one.
        assert_eq!(
            germline_probability(
                NORMAL_LOG_ODDS,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                0.0,
                0.0
            )
            .expect("one side survives"),
            0.0
        );
        // No somatic hypothesis and an impossible somatic prior.
        assert_eq!(
            germline_probability(
                NORMAL_LOG_ODDS,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                1.0,
                f64::NEG_INFINITY
            )
            .expect("one side survives"),
            1.0
        );
    }
}
