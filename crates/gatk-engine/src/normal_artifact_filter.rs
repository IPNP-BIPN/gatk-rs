//! `NormalArtifactFilter`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NormalArtifactFilter`
//! (GATK 4.6.2.0).
//!
//! "Does the matched normal carry this allele too?" Two answers are combined: the normal-artifact
//! log odds Mutect2 already wrote into the record, and a p-value computed here from the normal's
//! depth against the site's median base quality.
//!
//! # One index, read two different ways
//!
//! ```java
//! final int indexOfMaxTumorLod = MathUtils.maxElementIndex(tumorLods);
//! final int tumorAltDepth = tumorAlleleDepths[indexOfMaxTumorLod + 1];
//! final double[] normalArtifactNegativeLogOdds = ...;
//! ... normalArtifactNegativeLogOdds[indexOfMaxTumorLod] ...
//! ```
//!
//! `indexOfMaxTumorLod` is the **alternate** allele's index. The depth arrays are indexed by the
//! record's alleles and so need the `+ 1`; `TLOD` and `NALOD` are per-alternate and do not. An
//! off-by-one either way is a wrong number rather than a crash.
//!
//! # Only the normal side of the ratio gate is guarded
//!
//! ```java
//! final double tumorAlleleFraction = (double) tumorAltDepth / tumorDepth;
//! final double normalAlleleFraction = normalDepth == 0 ? 0 : (double) normalAltDepth / normalDepth;
//! if (normalAlleleFraction < MIN_NORMAL_ARTIFACT_RATIO * tumorAlleleFraction) { return 0.0; }
//! ```
//!
//! The tumour side is a bare division, so a record with no tumour depth makes the comparison
//! `0 < NaN`, which is **false**: instead of returning zero, such a record falls through to the
//! arithmetic below and can be filtered.
//!
//! # The imputed base quality is unreachable from a record that lacks the annotation
//!
//! ```java
//! final int medianRefBaseQuality = vc.getAttributeAsIntList(MEDIAN_BASE_QUALITY_KEY, IMPUTED_NORMAL_BASE_QUALITY).get(0);
//! ```
//!
//! `getAttributeAsIntList(key, default)` maps the default over the elements of a **present** list,
//! replacing a null or `"."`. An absent key is `Collections.emptyList()`, and `.get(0)` on that is
//! an `IndexOutOfBoundsException`. The constant named `IMPUTED_NORMAL_BASE_QUALITY` never stands in
//! for a missing `MBQ`; it stands in for a missing element of one that is there.
//!
//! And the quality it reads is the **site's** `MBQ`, whatever the comment beside it says about the
//! average base quality of reference reads in the normal.
//!
//! # A missing required annotation is zero here, and an empty list next door
//!
//! `Mutect2VariantFilter.errorProbabilities` copies one probability to every alternate allele, and
//! when a required annotation is absent that probability is `0.0`. The per-allele base class
//! ([`crate::allele_filter`]) answers an empty list in the same situation, and `ErrorProbabilities`
//! drops an empty list rather than counting it. The two base classes disagree.

use crate::allele_filter::{sum_ads_over_samples, AlleleDepthTooShort, GenotypeData};
use crate::allele_likelihoods::log10_to_log;
use crate::math_utils::{max_element_index, qual_to_error_prob};
use crate::mutect_engine::{posterior_probability_of_error, round_finite_precision_errors};
use crate::natural_log_utils::NonFiniteSum;
use crate::somatic_clustering_model::SomaticClusteringModel;
use jmath::binomial::{cumulative_probability, BinomialError};

/// `NormalArtifactFilter`'s identity.
pub const FILTER_NAME: &str = "normal_artifact";

/// `MIN_NORMAL_ARTIFACT_RATIO`: "don't call normal artifact if allele fraction in normal is much
/// smaller than allele fraction in tumor".
pub const MIN_NORMAL_ARTIFACT_RATIO: f64 = 0.1;

/// `IMPUTED_NORMAL_BASE_QUALITY`, "only used if normal base quality annotation fails somehow" --
/// which, per the module note, is not the same thing as the annotation being absent.
pub const IMPUTED_NORMAL_BASE_QUALITY: i32 = 30;

/// `M2FiltersArgumentCollection.DEFAULT_NORMAL_P_VALUE_THRESHOLD`.
pub const DEFAULT_NORMAL_P_VALUE_THRESHOLD: f64 = 0.001;

/// What this filter refuses, each of them something the reference throws or has not been measured.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalArtifactError {
    /// `.get(0)` on the empty list an absent `MBQ` answers.
    MedianBaseQualityMissing,
    /// A genotype's `AD` is shorter than the record's allele count.
    AlleleDepth(AlleleDepthTooShort),
    /// The winning alternate index is past the end of one of the per-allele arrays.
    IndexOutOfRange { index: usize, length: usize },
    /// `normalizeFromLogToLinearSpace` over a sum that is not finite.
    NonFiniteSum(NonFiniteSum),
    /// The binomial p-value, or the beta underneath it.
    Binomial(BinomialError),
    /// `maxElementIndex` over an empty array, which the reference refuses and nothing has measured.
    NoTumorLogOdds,
}

impl From<AlleleDepthTooShort> for NormalArtifactError {
    fn from(error: AlleleDepthTooShort) -> Self {
        NormalArtifactError::AlleleDepth(error)
    }
}

impl From<NonFiniteSum> for NormalArtifactError {
    fn from(error: NonFiniteSum) -> Self {
        NormalArtifactError::NonFiniteSum(error)
    }
}

impl From<BinomialError> for NormalArtifactError {
    fn from(error: BinomialError) -> Self {
        NormalArtifactError::Binomial(error)
    }
}

/// `NormalArtifactFilter.errorProbabilities(vc, engine, referenceContext)`.
///
/// `tumor_log_10_odds` and `normal_artifact_log_10_odds` are the record's `TLOD` and `NALOD` as
/// written, in log10; `None` is the annotation being absent, which is what the required-annotation
/// check looks at. `median_base_qualities` is the record's `MBQ`, empty when it is absent.
///
/// The answer is one probability per alternate allele, all of them the same number.
pub fn normal_artifact_error_probabilities<T>(
    model: &SomaticClusteringModel,
    tumor_log_10_odds: Option<&[f64]>,
    normal_artifact_log_10_odds: Option<&[f64]>,
    median_base_qualities: &[i32],
    genotypes: &[GenotypeData<T>],
    allele_count: usize,
    normal_pileup_p_value_threshold: f64,
) -> Result<Vec<f64>, NormalArtifactError> {
    let alternate_count = allele_count - 1;
    // `requiredInfoAnnotations().stream().allMatch(vc::hasAttribute) ? calculate... : 0.0`, and the
    // 0.0 goes through the same rounding as a computed probability before being copied.
    let (Some(tumor_log_10_odds), Some(normal_artifact_log_10_odds)) =
        (tumor_log_10_odds, normal_artifact_log_10_odds)
    else {
        return Ok(vec![round_finite_precision_errors(0.0); alternate_count]);
    };

    let probability = calculate_error_probability(
        model,
        tumor_log_10_odds,
        normal_artifact_log_10_odds,
        median_base_qualities,
        genotypes,
        allele_count,
        normal_pileup_p_value_threshold,
    )?;
    Ok(vec![
        round_finite_precision_errors(probability);
        alternate_count
    ])
}

/// `calculateErrorProbability`, without the rounding and without the copy.
#[allow(clippy::too_many_arguments)]
fn calculate_error_probability<T>(
    model: &SomaticClusteringModel,
    tumor_log_10_odds: &[f64],
    normal_artifact_log_10_odds: &[f64],
    median_base_qualities: &[i32],
    genotypes: &[GenotypeData<T>],
    allele_count: usize,
    normal_pileup_p_value_threshold: f64,
) -> Result<f64, NormalArtifactError> {
    if tumor_log_10_odds.is_empty() {
        return Err(NormalArtifactError::NoTumorLogOdds);
    }
    // `getTumorLogOdds` converts to natural log first. The conversion is monotone, so the index it
    // is taken over is the same one the log10 values would give; it is written the reference's way
    // rather than shortened, because the values are read again below.
    let tumor_lods: Vec<f64> = tumor_log_10_odds.iter().map(|x| log10_to_log(*x)).collect();
    let index_of_max_tumor_lod = max_element_index(&tumor_lods, 0, tumor_lods.len());

    let tumor_allele_depths = sum_ads_over_samples(allele_count, genotypes, true, false)?;
    // `(int) MathUtils.sum(int[])`, which sums into a `long` and is then truncated to an `int`.
    let tumor_depth = (tumor_allele_depths
        .iter()
        .map(|d| i64::from(*d))
        .sum::<i64>()) as i32;
    let tumor_alt_depth = at(&tumor_allele_depths, index_of_max_tumor_lod + 1)?;

    let normal_allele_depths = sum_ads_over_samples(allele_count, genotypes, false, true)?;
    let normal_depth = (normal_allele_depths
        .iter()
        .map(|d| i64::from(*d))
        .sum::<i64>()) as i32;
    let normal_alt_depth = at(&normal_allele_depths, index_of_max_tumor_lod + 1)?;

    // The tumour side is unguarded: an empty tumour is `0 / 0`, which is NaN.
    let tumor_allele_fraction = f64::from(tumor_alt_depth) / f64::from(tumor_depth);
    let normal_allele_fraction = if normal_depth == 0 {
        0.0
    } else {
        f64::from(normal_alt_depth) / f64::from(normal_depth)
    };

    // `x < NaN` is false, so a record with no tumour depth does NOT return here.
    if normal_allele_fraction < MIN_NORMAL_ARTIFACT_RATIO * tumor_allele_fraction {
        return Ok(0.0);
    }

    // `applyToArrayInPlace(getAttributeAsDoubleArray(vc, NALOD), x -> -log10ToLog(x))`, which
    // negates in the natural-log space rather than negating the log10 value first.
    let negative_log_odds =
        -log10_to_log(at_f64(normal_artifact_log_10_odds, index_of_max_tumor_lod)?);
    let normal_artifact_probability = posterior_probability_of_error(
        negative_log_odds,
        model.log_prior_of_variant_versus_artifact(),
    )?;

    // The empty list an absent `MBQ` answers, not the imputed quality.
    let median_ref_base_quality = *median_base_qualities
        .first()
        .ok_or(NormalArtifactError::MedianBaseQualityMissing)?;
    let normal_p_value = 1.0
        - cumulative_probability(
            normal_depth,
            qual_to_error_prob(f64::from(median_ref_base_quality)),
            normal_alt_depth - 1,
        )?;

    Ok(if normal_p_value < normal_pileup_p_value_threshold {
        1.0
    } else {
        normal_artifact_probability
    })
}

/// One element of a depth array, refusing where the reference throws `ArrayIndexOutOfBounds`.
fn at(values: &[i32], index: usize) -> Result<i32, NormalArtifactError> {
    values
        .get(index)
        .copied()
        .ok_or(NormalArtifactError::IndexOutOfRange {
            index,
            length: values.len(),
        })
}

/// The same, for the log-odds arrays.
fn at_f64(values: &[f64], index: usize) -> Result<f64, NormalArtifactError> {
    values
        .get(index)
        .copied()
        .ok_or(NormalArtifactError::IndexOutOfRange {
            index,
            length: values.len(),
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::somatic_clustering_model::PriorArguments;

    fn genotypes(tumor: &[i32], normal: &[i32]) -> Vec<GenotypeData<i32>> {
        vec![
            GenotypeData {
                tumor: true,
                allele_depths: tumor.to_vec(),
                values: Vec::new(),
            },
            GenotypeData {
                tumor: false,
                allele_depths: normal.to_vec(),
                values: Vec::new(),
            },
        ]
    }

    fn probabilities(
        tumor: &[i32],
        normal: &[i32],
        median_base_qualities: &[i32],
    ) -> Result<Vec<f64>, NormalArtifactError> {
        let model = SomaticClusteringModel::new(PriorArguments::new(), None);
        normal_artifact_error_probabilities(
            &model,
            Some(&[20.0, 6.0]),
            Some(&[2.0, 0.5]),
            median_base_qualities,
            &genotypes(tumor, normal),
            3,
            DEFAULT_NORMAL_P_VALUE_THRESHOLD,
        )
    }

    /// A tumour with no depth makes the gate `0 < NaN`, which is false: the record is filtered
    /// rather than passed. A gate written the obvious way would return zero here.
    #[test]
    fn a_tumour_with_no_depth_falls_through_the_ratio_gate() {
        let clean_normal =
            probabilities(&[80, 20, 5], &[90, 0, 0], &[30, 30, 30]).expect("answered");
        assert_eq!(clean_normal, vec![0.0, 0.0], "the gate returns zero");
        let no_tumour = probabilities(&[0, 0, 0], &[90, 0, 0], &[30, 30, 30]).expect("answered");
        assert!(no_tumour[0] > 0.0, "and NaN does not");
    }

    /// The imputed quality never stands in for an absent annotation.
    #[test]
    fn a_missing_median_base_quality_is_a_refusal_and_not_the_imputed_thirty() {
        assert_eq!(
            probabilities(&[80, 20, 5], &[90, 10, 2], &[]),
            Err(NormalArtifactError::MedianBaseQualityMissing)
        );
        // And the constant that would have been used is still what the reference declares.
        assert_eq!(IMPUTED_NORMAL_BASE_QUALITY, 30);
    }

    /// A missing required annotation is a zero per alternate allele, not an empty list.
    #[test]
    fn a_missing_annotation_answers_one_zero_per_alternate_allele() {
        let model = SomaticClusteringModel::new(PriorArguments::new(), None);
        let answer = normal_artifact_error_probabilities(
            &model,
            None,
            Some(&[2.0, 0.5]),
            &[30, 30, 30],
            &genotypes(&[80, 20, 5], &[90, 10, 2]),
            3,
            DEFAULT_NORMAL_P_VALUE_THRESHOLD,
        )
        .expect("answered");
        assert_eq!(answer, vec![0.0, 0.0]);
    }
}
