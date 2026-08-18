//! `PolymeraseSlippageFilter`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.PolymeraseSlippageFilter`
//! (GATK 4.6.2.0).
//!
//! "Is this indel a PCR artifact of a short tandem repeat?" A somatic likelihood is weighed against
//! the likelihood that a polymerase slipped, and the second of those is a regularized beta.
//!
//! # The beta is about one allele and the likelihood is about all of them
//!
//! ```java
//! final int depth = (int) MathUtils.sum(ADs);
//! final int altCount = (int) MathUtils.sum(ADs) - ADs[0];
//! final double logSomaticLikelihood = ...logLikelihoodGivenSomatic(depth, altCount);
//! likelihoodGivenSlippageArtifact = Beta.regularizedBeta(slippageRate, ADs[1] + 1, ADs[0] + 1);
//! ```
//!
//! `altCount` sums **every** alternate allele's depth; `ADs[1]` is the **first** alternate's alone.
//! On a biallelic record they are the same number. On a triallelic one the two halves of the same
//! log-odds are about different alleles, and nothing reconciles them.
//!
//! # The prior's allele index is hard-coded
//!
//! `filteringEngine.posteriorProbabilityOfError(vc, logOdds, 0)` reads
//! `getLogPriorOfSomaticVariant(vc, 0)`, so the prior comes from alternate **zero**'s indel length
//! whichever allele slipped. The one probability that comes out is then copied to every alternate
//! by `Mutect2VariantFilter.errorProbabilities`.
//!
//! # `RPA` is parsed, not read
//!
//! ```java
//! vc.getAttributeAsList(REPEATS_PER_ALLELE_KEY).stream()
//!         .mapToInt(o -> Integer.parseInt(String.valueOf(o))).toArray();
//! ```
//!
//! `Integer.parseInt` is handed the *string form* of whatever the attribute holds, so an `RPA` of
//! `10.0` is a `NumberFormatException` on the string `"10.0"` rather than a truncation to ten.
//! `rpa.length < 2` is the only length check there is.
//!
//! # The fallback cannot be reached
//!
//! The `catch (MaxCountExceededException)` arm computes a binomial probability instead. Three-argument
//! `regularizedBeta` passes `Integer.MAX_VALUE` iterations, so nothing short of a diverging continued
//! fraction reaches it, and no golden can exercise it. It is [`SlippageError::FallbackUnmeasured`]
//! here rather than a guess.

use crate::allele_filter::{sum_ads_over_samples, AlleleDepthTooShort, GenotypeData};
use crate::mutect_engine::{posterior_probability_of_error, round_finite_precision_errors};
use crate::natural_log_utils::NonFiniteSum;
use crate::somatic_clustering_model::{indel_length, AlternateAllele, SomaticClusteringModel};
use jmath::beta::{regularized_beta, BetaError};
use jmath::continued_fraction::ContinuedFractionError;

/// `PolymeraseSlippageFilter`'s identity.
pub const FILTER_NAME: &str = "slippage";

/// `phredScaledPosteriorAnnotationName`, `POLYMERASE_SLIPPAGE_QUAL_KEY`.
pub const ANNOTATION: &str = "STRQ";

/// `M2FiltersArgumentCollection.DEFAULT_MIN_SLIPPAGE_LENGTH`.
pub const DEFAULT_MIN_SLIPPAGE_LENGTH: i32 = 8;

/// `M2FiltersArgumentCollection.DEFAULT_SLIPPAGE_RATE`.
pub const DEFAULT_SLIPPAGE_RATE: f64 = 0.1;

/// What this filter refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum SlippageError {
    /// `Integer.parseInt` on an `RPA` entry that is not an integer.
    NumberFormat { input: String },
    /// A genotype's `AD` is shorter than the record's allele count.
    AlleleDepth(AlleleDepthTooShort),
    /// `logSumExp` over a sum that is not finite.
    NonFiniteSum(NonFiniteSum),
    /// The beta refused for a reason that is not the unreachable iteration cap.
    Beta(BetaError),
    /// The continued fraction ran out of iterations, which is the arm whose `BinomialDistribution`
    /// fallback no golden can reach.
    FallbackUnmeasured,
}

impl SlippageError {
    /// The exception class the reference throws, for the refusals that have one.
    pub fn class(&self) -> Option<&'static str> {
        match self {
            SlippageError::NumberFormat { .. } => Some("java.lang.NumberFormatException"),
            SlippageError::AlleleDepth(_) => Some("java.lang.ArrayIndexOutOfBoundsException"),
            _ => None,
        }
    }

    /// `NumberFormatException`'s message, which quotes the input.
    pub fn message(&self) -> Option<String> {
        match self {
            SlippageError::NumberFormat { input } => Some(format!("For input string: \"{input}\"")),
            SlippageError::AlleleDepth(error) => Some(error.message()),
            _ => None,
        }
    }
}

impl From<AlleleDepthTooShort> for SlippageError {
    fn from(error: AlleleDepthTooShort) -> Self {
        SlippageError::AlleleDepth(error)
    }
}

impl From<NonFiniteSum> for SlippageError {
    fn from(error: NonFiniteSum) -> Self {
        SlippageError::NonFiniteSum(error)
    }
}

impl From<BetaError> for SlippageError {
    fn from(error: BetaError) -> Self {
        match error {
            BetaError::ContinuedFraction(ContinuedFractionError::NotConvergent { .. }) => {
                SlippageError::FallbackUnmeasured
            }
            other => SlippageError::Beta(other),
        }
    }
}

/// `PolymeraseSlippageFilter.errorProbabilities(vc, engine, referenceContext)`.
///
/// `repeats_per_allele` is `RPA` as written, one string per entry, because the reference parses the
/// string form; `None` is the annotation being absent. `repeat_unit` is `RU`, likewise. The answer
/// is one probability per alternate allele, all of them the same number.
#[allow(clippy::too_many_arguments)]
pub fn slippage_error_probabilities<T>(
    model: &mut SomaticClusteringModel,
    repeats_per_allele: Option<&[String]>,
    repeat_unit: Option<&str>,
    genotypes: &[GenotypeData<T>],
    alternates: &[AlternateAllele],
    reference_length: i32,
    min_slippage_length: i32,
    slippage_rate: f64,
) -> Result<Vec<f64>, SlippageError> {
    let alternate_count = alternates.len();
    let (Some(repeats_per_allele), Some(repeat_unit)) = (repeats_per_allele, repeat_unit) else {
        return Ok(vec![round_finite_precision_errors(0.0); alternate_count]);
    };

    let probability = calculate_error_probability(
        model,
        repeats_per_allele,
        repeat_unit,
        genotypes,
        alternates,
        reference_length,
        min_slippage_length,
        slippage_rate,
    )?;
    Ok(vec![
        round_finite_precision_errors(probability);
        alternate_count
    ])
}

/// `calculateErrorProbability`, without the copy.
#[allow(clippy::too_many_arguments)]
fn calculate_error_probability<T>(
    model: &mut SomaticClusteringModel,
    repeats_per_allele: &[String],
    repeat_unit: &str,
    genotypes: &[GenotypeData<T>],
    alternates: &[AlternateAllele],
    reference_length: i32,
    min_slippage_length: i32,
    slippage_rate: f64,
) -> Result<f64, SlippageError> {
    let mut repeats = Vec::with_capacity(repeats_per_allele.len());
    for entry in repeats_per_allele {
        repeats.push(parse_int(entry)?);
    }
    if repeats.len() < 2 {
        return Ok(0.0);
    }

    // `ru.length()` counts UTF-16 code units, which for a repeat unit of DNA is its bases.
    let reference_str_base_count = repeat_unit.encode_utf16().count() as i32 * repeats[0];
    let number_of_pcr_slips = repeats[0] - repeats[1];
    if !(reference_str_base_count >= min_slippage_length && number_of_pcr_slips.abs() == 1) {
        return Ok(0.0);
    }

    let allele_depths = sum_ads_over_samples(alternates.len() + 1, genotypes, true, false)?;
    if allele_depths.len() < 2 {
        return Ok(0.0);
    }
    // `(int) MathUtils.sum(int[])`, which sums into a `long` and is then truncated.
    let depth = (allele_depths.iter().map(|d| i64::from(*d)).sum::<i64>()) as i32;
    // The reference sums the array a second time rather than reusing `depth`, which is the same
    // number by construction.
    let alt_count = depth - allele_depths[0];
    let log_somatic_likelihood = model.log_likelihood_given_somatic(depth, alt_count)?;

    // `ADs[1]`, the FIRST alternate's depth, against `alt_count` above, which summed them all.
    let likelihood_given_slippage_artifact = regularized_beta(
        slippage_rate,
        f64::from(allele_depths[1]) + 1.0,
        f64::from(allele_depths[0]) + 1.0,
    )?;

    let log_odds = log_somatic_likelihood - jmath::math::log(likelihood_given_slippage_artifact);
    // `posteriorProbabilityOfError(vc, logOdds, 0)`: alternate ZERO's indel length, always.
    let prior = model.log_prior_of_somatic_variant(indel_length(reference_length, alternates[0]));
    Ok(posterior_probability_of_error(log_odds, prior)?)
}

/// `Integer.parseInt(String.valueOf(o))`, refusing where the reference throws.
fn parse_int(text: &str) -> Result<i32, SlippageError> {
    text.parse::<i32>()
        .map_err(|_| SlippageError::NumberFormat {
            input: text.to_string(),
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::somatic_clustering_model::PriorArguments;

    fn genotypes() -> Vec<GenotypeData<i32>> {
        vec![
            GenotypeData {
                tumor: true,
                allele_depths: vec![80, 20],
                values: Vec::new(),
            },
            GenotypeData {
                tumor: false,
                allele_depths: vec![90, 1],
                values: Vec::new(),
            },
        ]
    }

    fn answer(repeats: &[&str], repeat_unit: &str) -> Result<Vec<f64>, SlippageError> {
        let mut model = SomaticClusteringModel::new(PriorArguments::new(), None);
        let repeats: Vec<String> = repeats.iter().map(|r| r.to_string()).collect();
        slippage_error_probabilities(
            &mut model,
            Some(&repeats),
            Some(repeat_unit),
            &genotypes(),
            &[AlternateAllele {
                length: 1,
                symbolic: false,
            }],
            2,
            DEFAULT_MIN_SLIPPAGE_LENGTH,
            DEFAULT_SLIPPAGE_RATE,
        )
    }

    /// The gate needs both halves, and the base count is the repeat unit's length times the
    /// reference's repeat count.
    #[test]
    fn both_halves_of_the_gate_are_needed() {
        // Ten bases and one slip: filtered.
        assert!(answer(&["10", "9"], "A").expect("answered")[0] > 0.0);
        // Ten bases and two slips.
        assert_eq!(answer(&["10", "8"], "A").expect("answered"), vec![0.0]);
        // One slip and seven bases.
        assert_eq!(answer(&["7", "6"], "A").expect("answered"), vec![0.0]);
        // A two-base repeat unit reaches the minimum with half the repeats.
        assert!(answer(&["5", "4"], "AT").expect("answered")[0] > 0.0);
        // And an empty repeat unit is zero bases however long the repeat.
        assert_eq!(answer(&["100", "99"], "").expect("answered"), vec![0.0]);
    }

    /// `Integer.parseInt` is handed the string form, so a decimal is a refusal rather than a
    /// truncation, and the message quotes what it was given.
    #[test]
    fn a_decimal_repeat_count_is_refused_with_the_string_it_was_given() {
        let error = answer(&["10.0", "9.0"], "A").expect_err("refused");
        assert_eq!(error.class(), Some("java.lang.NumberFormatException"));
        assert_eq!(
            error.message(),
            Some("For input string: \"10.0\"".to_string())
        );
    }
}
