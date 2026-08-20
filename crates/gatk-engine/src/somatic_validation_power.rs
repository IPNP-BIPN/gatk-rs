//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.PowerCalculationUtils`
//! (GATK 4.6.2.0).
//!
//! How likely a validation pileup of a given depth is to show a variant the discovery pileup
//! called, and how many reads count as showing it.
//!
//! # The power is one minus a beta-binomial tail
//!
//! The distribution's shapes are the discovery counts PLUS ONE -- `alt + 1` and
//! `total - alt + 1` -- so they are never zero and the constructor's own refusals are unreachable
//! from here. The tail is `1 - cumulativeProbability(minCount - 1)`, and that cumulative
//! probability is a plain uncompensated loop over `exp(logProbability(i))`.
//!
//! # The minimum count is a quantile floored at two
//!
//! `max(inverseCumulativeProbability(0.99), 2)` over a binomial of the validation depth at the
//! noise ratio. The floor is what makes a pileup of zero reads and a pileup of 317 clean ones both
//! answer two.

use crate::beta_binomial::{BetaBinomialDistribution, BetaBinomialError};
use crate::pileup::PileupElement;
use crate::read_pileup::ReadPileup;
use crate::variant_context_utils::{
    choose_allele_for_read, does_read_contain_allele, Allele, PileupAlleleError, Trilean,
};

/// `PowerCalculationUtils.P_VALUE_FOR_NOISE`.
pub const P_VALUE_FOR_NOISE: f64 = 0.99;

/// `PowerCalculationUtils.MINIMUM_NUM_READS_FOR_SIGNAL_COUNT`.
pub const MINIMUM_NUM_READS_FOR_SIGNAL_COUNT: i32 = 2;

/// What the power calculation refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum PowerError {
    /// `ParamUtils.isPositiveOrZero` on the validation depth.
    NegativeTotalCount,
    /// `ParamUtils.inRange` on the noise ratio.
    RatioOutOfRange,
    /// The distribution refused.
    Distribution(BetaBinomialError),
    /// The quantile refused.
    Quantile(jmath::binomial::BinomialError),
    /// `ParamUtils.isPositiveOrZero` on the base-quality cutoff, or the allele choice underneath.
    Allele(PileupAlleleError),
}

impl PowerError {
    /// The exception class the reference throws.
    pub fn java_class(&self) -> &'static str {
        match self {
            PowerError::Allele(error) => error.java_class(),
            _ => "java.lang.IllegalArgumentException",
        }
    }

    /// The message, including the reference's own doubled word in the second.
    pub fn message(&self) -> String {
        match self {
            PowerError::NegativeTotalCount => "Cannot have a negative total count.".to_string(),
            PowerError::RatioOutOfRange => {
                "Cannot have have a ratio that is outside of 0.0 - 1.0.".to_string()
            }
            PowerError::Distribution(error) => format!("{error:?}"),
            PowerError::Quantile(error) => format!("{error:?}"),
            PowerError::Allele(error) => error.message(),
        }
    }
}

/// `calculatePower`.
pub fn calculate_power(
    validation_total_count: i32,
    discovery_alt_count: i32,
    discovery_total_count: i32,
    minimum_count_for_signal: i32,
) -> Result<f64, PowerError> {
    let distribution = BetaBinomialDistribution::new(
        f64::from(discovery_alt_count) + 1.0,
        f64::from(discovery_total_count - discovery_alt_count) + 1.0,
        validation_total_count,
    )
    .map_err(PowerError::Distribution)?;
    let cumulative = distribution
        .cumulative_probability(minimum_count_for_signal - 1)
        .map_err(PowerError::Distribution)?;
    Ok(1.0 - cumulative)
}

/// `calculateMinCountForSignal`: the binomial quantile at 0.99, floored at two.
pub fn calculate_min_count_for_signal(
    validation_total_count: i32,
    max_signal_ratio_in_normal: f64,
) -> Result<i32, PowerError> {
    if validation_total_count < 0 {
        return Err(PowerError::NegativeTotalCount);
    }
    // `ParamUtils.inRange`, inclusive at both ends.
    if !(0.0..=1.0).contains(&max_signal_ratio_in_normal) {
        return Err(PowerError::RatioOutOfRange);
    }
    let quantile = jmath::binomial::inverse_cumulative_probability(
        validation_total_count,
        max_signal_ratio_in_normal,
        P_VALUE_FOR_NOISE,
    )
    .map_err(PowerError::Quantile)?;
    Ok(quantile.max(MINIMUM_NUM_READS_FOR_SIGNAL_COUNT))
}

/// `retrievePileupElements`: not a deletion, and at or above the cutoff.
///
/// The quality here is the pileup element's own base quality, which for a deletion would be the
/// constant 16, so the deletion filter runs first and the constant never decides anything.
fn elements_passing_quality<'a>(
    pileup: &ReadPileup<'a>,
    minimum_base_quality: i32,
) -> Vec<PileupElement<'a>> {
    pileup
        .elements
        .iter()
        .filter(|element| !element.is_deletion())
        .filter(|element| i32::from(element.qual()) >= minimum_base_quality)
        .cloned()
        .collect()
}

/// `calculateMaxAltRatio`: the fraction of the pileup that is not the reference allele.
///
/// THE TWO FILTERS ARE NOT COMPLEMENTS. An element is alternate when it does not contain the
/// reference allele OR it precedes an indel; it is reference when it does contain the reference
/// allele AND precedes neither. An element that answers UNKNOWN, which is a read ending inside the
/// allele, is in neither count, so a pileup of nothing but such reads has a denominator of zero and
/// the ratio is the literal 0.0 rather than a NaN. The caller depends on that: a NaN ratio is what
/// `calculateBasicValidationResult` refuses on.
pub fn calculate_max_alt_ratio(
    pileup: &ReadPileup<'_>,
    reference: &Allele,
    minimum_base_quality: i32,
) -> Result<f64, PowerError> {
    if minimum_base_quality < 0 {
        return Err(PowerError::Allele(
            PileupAlleleError::NegativeMinimumBaseQualityRatio,
        ));
    }
    let passing = elements_passing_quality(pileup, minimum_base_quality);
    let alternate = passing
        .iter()
        .filter(|element| {
            does_read_contain_allele(element, reference) == Trilean::False
                || element.is_before_deletion_start()
                || element.is_before_insertion()
        })
        .count();
    let reference_count = passing
        .iter()
        .filter(|element| {
            does_read_contain_allele(element, reference) == Trilean::True
                && !element.is_before_deletion_start()
                && !element.is_before_insertion()
        })
        .count();
    if reference_count + alternate == 0 {
        return Ok(0.0);
    }
    Ok(alternate as f64 / (reference_count as f64 + alternate as f64))
}

/// `calculateNumReadsSupportingAllele`: the elements whose chosen allele is this alternate.
///
/// The choice is made against a list of exactly one alternate, so an element supporting a different
/// alternate is not counted here even when the pileup was built for several.
pub fn calculate_num_reads_supporting_allele(
    pileup: &ReadPileup<'_>,
    reference: &Allele,
    alternate: &Allele,
    minimum_base_quality: i32,
) -> Result<i64, PowerError> {
    if minimum_base_quality < 0 {
        return Err(PowerError::Allele(
            PileupAlleleError::NegativeMinimumBaseQualityRatio,
        ));
    }
    let alternates = [alternate.clone()];
    let mut count = 0;
    for element in elements_passing_quality(pileup, minimum_base_quality) {
        let chosen = choose_allele_for_read(&element, reference, &alternates, minimum_base_quality)
            .map_err(PowerError::Allele)?;
        if chosen.as_ref() == Some(alternate) {
            count += 1;
        }
    }
    Ok(count)
}
