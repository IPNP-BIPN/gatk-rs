//! `ContaminationFilter`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ContaminationFilter`
//! (GATK 4.6.2.0).
//!
//! "Could this allele have come from another sample rather than from this one?" One of only two
//! filters in the `NON_SOMATIC` error type, which is what decides which other filters can mask it.
//!
//! # It can answer NaN
//!
//! ```java
//! new IndexRange(0, vc.getNAlleles()-1).forEach(i -> depthsAndPosteriorsPerAllele.add(new ArrayList<>()));
//! ...
//! new IndexRange(0, alleleFrequencies.length).forEach(i -> { ... });
//! ...
//! return depthsAndPosteriorsPerAllele.stream()
//!         .map(alleleData -> alleleData.isEmpty() ? Double.NaN : weightedMedianPosteriorProbability(alleleData))
//! ```
//!
//! The list is sized from the **record's** alleles and filled from **`POPAF`'s** length. An allele
//! the annotation does not cover keeps an empty list, and an empty list is `Double.NaN`, which
//! `roundFinitePrecisionErrors` passes through unchanged: `Math.min` and `Math.max` propagate NaN.
//! A record with no tumour sample at all answers NaN for every allele, and `ErrorProbabilities`
//! counts that as a filter that ran rather than dropping it the way it drops an empty list.
//!
//! # The other direction is an exception
//!
//! The same loop indexes `altADs`, which is sized from the genotype's `AD`, so a `POPAF` **longer**
//! than the record is an `ArrayIndexOutOfBoundsException` rather than a truncation.
//!
//! # Two hypotheses, compared by maximum
//!
//! ```java
//! logContaminantLikelihoodPerAllele[i] = Math.log(Math.max(single, many));
//! ```
//!
//! One contaminating sample carrying the allele, against many contaminating samples at the
//! population frequency. The larger likelihood wins outright rather than the two being summed, so
//! which hypothesis explains a site can change with the depth alone.
//!
//! # Everything collapses at a contamination of zero
//!
//! The default estimate is `0.0`, and a binomial at `p = 0` is zero for any positive count, so the
//! contaminant likelihood is zero, its logarithm negative infinity, the log-odds positive infinity
//! and the posterior of error exactly `0.0`. A run given neither a table nor an estimate builds a
//! filter that answers zero to everything.
//!
//! # The prior's allele index is the loop's
//!
//! `posteriorProbabilityOfError(vc, logOdds, i)`, so each alternate is priced with its own indel
//! length. [`crate::slippage_filter`] hard-codes zero in the same call.

use crate::allele_filter::{weighted_median_posterior_probability, GenotypeData};
use crate::math_utils::pow10;
use crate::mutect_engine::{posterior_probability_of_error, round_finite_precision_errors};
use crate::natural_log_utils::NonFiniteSum;
use crate::somatic_clustering_model::{
    binomial_probability, indel_length, AlternateAllele, SomaticClusteringModel,
};

/// `ContaminationFilter`'s identity.
pub const FILTER_NAME: &str = "contamination";

/// `phredScaledPosteriorAnnotationName`, `CONTAMINATION_QUAL_KEY`.
pub const ANNOTATION: &str = "CONTQ";

/// `M2FiltersArgumentCollection.DEFAULT_CONTAMINATION`.
pub const DEFAULT_CONTAMINATION: f64 = 0.0;

/// `EPSILON`, which is an instance field of the filter rather than a constant.
pub const EPSILON: f64 = 1.0e-10;

/// What this filter refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum ContaminationError {
    /// `altADs[i]` past the end, which a `POPAF` longer than the record reaches.
    IndexOutOfRange { index: usize, length: usize },
    /// `logSumExp` over a sum that is not finite.
    NonFiniteSum(NonFiniteSum),
}

impl ContaminationError {
    pub fn class(&self) -> Option<&'static str> {
        match self {
            ContaminationError::IndexOutOfRange { .. } => {
                Some("java.lang.ArrayIndexOutOfBoundsException")
            }
            ContaminationError::NonFiniteSum(_) => None,
        }
    }

    pub fn message(&self) -> Option<String> {
        match self {
            ContaminationError::IndexOutOfRange { index, length } => {
                Some(format!("Index {index} out of bounds for length {length}"))
            }
            ContaminationError::NonFiniteSum(_) => None,
        }
    }
}

impl From<NonFiniteSum> for ContaminationError {
    fn from(error: NonFiniteSum) -> Self {
        ContaminationError::NonFiniteSum(error)
    }
}

/// `Math.max(0, Math.min(contaminationFromFile, 1 - EPSILON))`.
///
/// The comment beside it says the upper clamp is there "to handle file with contamination == 1", so
/// a contamination of one and a contamination of two are the same number here.
pub fn clamp_contamination(contamination: f64) -> f64 {
    // `Math.min`/`Math.max` propagate NaN where Rust's do not, and this is the reference's order:
    // written as `clamp`, the guard above reads like a redundant special case rather than the thing
    // that makes the two forms equivalent.
    if contamination.is_nan() {
        return contamination;
    }
    #[allow(clippy::manual_clamp)]
    contamination.min(1.0 - EPSILON).max(0.0)
}

/// `ContaminationFilter.errorProbabilities(vc, engine, referenceContext)`.
///
/// `population_negative_log10_allele_frequencies` is `POPAF` as written, one entry per alternate
/// allele; `None` is the annotation being absent, which the required-annotation check answers with
/// an empty list.
///
/// `contaminations` is one value per genotype, already resolved:
/// `contaminationBySample.getOrDefault(sampleName, defaultContamination)`. The lookup itself, and
/// the `IllegalStateException` two tables naming one sample provoke out of `Collectors.toMap`,
/// belong to the table reader rather than here; nothing measured reaches them.
pub fn contamination_error_probabilities<T>(
    model: &mut SomaticClusteringModel,
    population_negative_log10_allele_frequencies: Option<&[f64]>,
    genotypes: &[GenotypeData<T>],
    contaminations: &[f64],
    alternates: &[AlternateAllele],
    reference_length: i32,
) -> Result<Vec<f64>, ContaminationError> {
    let Some(negative_log10_frequencies) = population_negative_log10_allele_frequencies else {
        return Ok(Vec::new());
    };

    // One list per ALTERNATE allele of the record, whatever length the annotation turns out to be.
    let mut depths_and_posteriors: Vec<Vec<(i32, f64)>> = vec![Vec::new(); alternates.len()];

    for (sample, genotype) in genotypes.iter().enumerate() {
        // `if (filteringEngine.isNormal(tumorGenotype)) continue;`
        if !genotype.tumor {
            continue;
        }
        let contamination = clamp_contamination(contaminations[sample]);
        let allele_depths = &genotype.allele_depths;
        let total_ad = (allele_depths.iter().map(|d| i64::from(*d)).sum::<i64>()) as i32;
        let alt_allele_depths = &allele_depths[1..];

        // `MathUtils.applyToArray(negativeLog10AlleleFrequencies, x -> Math.pow(10, -x))`.
        let allele_frequencies: Vec<f64> = negative_log10_frequencies
            .iter()
            .map(|x| pow10(-x))
            .collect();

        // Over EVERY alternate depth, not over the frequencies.
        let mut log_somatic_likelihood = Vec::with_capacity(alt_allele_depths.len());
        for alt_count in alt_allele_depths {
            log_somatic_likelihood.push(model.log_likelihood_given_somatic(total_ad, *alt_count)?);
        }

        let mut log_odds = Vec::with_capacity(allele_frequencies.len());
        for (index, frequency) in allele_frequencies.iter().enumerate() {
            let alt_count = *at(alt_allele_depths, index)?;
            // One contaminating sample, heterozygous or homozygous for the allele.
            let single_contaminant = 2.0
                * frequency
                * (1.0 - frequency)
                * binomial_probability(total_ad, alt_count, contamination / 2.0)
                + (frequency * frequency)
                    * binomial_probability(total_ad, alt_count, contamination);
            // Many contaminating samples, each at the population frequency.
            let many_contaminant =
                binomial_probability(total_ad, alt_count, contamination * frequency);
            // `Math.max`, which propagates NaN where Rust's `f64::max` does not.
            let larger = if single_contaminant.is_nan() || many_contaminant.is_nan() {
                f64::NAN
            } else {
                single_contaminant.max(many_contaminant)
            };
            let log_contaminant_likelihood = jmath::math::log(larger);
            log_odds.push(log_somatic_likelihood[index] - log_contaminant_likelihood);
        }

        // A second pass, as the reference writes it: the priors are read after every log-odds.
        for (index, odds) in log_odds.iter().enumerate() {
            let prior = model
                .log_prior_of_somatic_variant(indel_length(reference_length, alternates[index]));
            let posterior = posterior_probability_of_error(*odds, prior)?;
            depths_and_posteriors[index].push((*at(alt_allele_depths, index)?, posterior));
        }
    }

    Ok(depths_and_posteriors
        .iter_mut()
        .map(|allele_data| {
            let value = if allele_data.is_empty() {
                f64::NAN
            } else {
                weighted_median_posterior_probability(allele_data)
            };
            round_finite_precision_errors(value)
        })
        .collect())
}

/// One alternate depth, refusing where the reference throws `ArrayIndexOutOfBoundsException`.
fn at(values: &[i32], index: usize) -> Result<&i32, ContaminationError> {
    values
        .get(index)
        .ok_or(ContaminationError::IndexOutOfRange {
            index,
            length: values.len(),
        })
}
