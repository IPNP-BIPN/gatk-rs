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
}

impl PowerError {
    /// The exception class the reference throws for the two argument checks.
    pub fn java_class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
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
