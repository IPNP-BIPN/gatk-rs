//! `BetaBinomialDistribution`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup` (GATK 4.6.2.0).
//!
//! A binomial whose success probability is itself drawn from a beta. `SomaticClusteringModel` uses
//! it for every cluster it carries: the flat background at `alpha = beta = 1`, and the high-AF
//! cluster at `alpha = 10, beta = 1`.
//!
//! # Three commons-math calls
//!
//! ```java
//! return k > n ? Double.NEGATIVE_INFINITY :
//!         CombinatoricsUtils.binomialCoefficientLog(n, k) + Beta.logBeta(k + alpha, n - k + beta)
//!                 - Beta.logBeta(alpha, beta);
//! ```
//!
//! All three live in [`jmath`], where the branch structure that makes them what they are is
//! documented. What this module adds is the argument checking and the order of the three terms.
//!
//! # The flat beta is not exactly flat
//!
//! At `alpha = beta = 1` the distribution is uniform over `0..=n` in the mathematics, and the
//! reference does not quite say so in doubles: at `n = 10` the answer is `-2.3978952727983707` at
//! `k = 0`, `-2.3978952727983716` at `k = 1` and `-2.39789527279837` at `k = 5`. The cancellation
//! between the coefficient and the two beta terms is approximate, and it lands on three different
//! doubles. A port that special-cased the uniform case would be smoother than the reference.
//!
//! # `k > n` is a branch, not arithmetic
//!
//! Negative infinity comes from the ternary rather than from a coefficient of zero, and it is
//! checked *after* the negative-`k` refusal.

use jmath::beta::{log_beta, BetaError};
use jmath::combinatorics::{binomial_coefficient_log, CombinatoricsError};

/// What the distribution refuses, each of them a `ParamUtils` check in the reference.
#[derive(Debug, Clone, PartialEq)]
pub enum BetaBinomialError {
    /// `alpha must be greater than zero.`
    AlphaNotPositive { alpha: f64 },
    /// `beta must be greater than zero.`
    BetaNotPositive { beta: f64 },
    /// `number of trials must be greater than (or equal to) zero.`
    TrialsNegative { n: i32 },
    /// `Number of successes must be greater than or equal to zero.`
    SuccessesNegative { k: i32 },
    /// A beta the port has not measured.
    Beta(BetaError),
    /// A coefficient the port has not measured.
    Combinatorics(CombinatoricsError),
}

impl From<BetaError> for BetaBinomialError {
    fn from(error: BetaError) -> Self {
        BetaBinomialError::Beta(error)
    }
}

impl From<CombinatoricsError> for BetaBinomialError {
    fn from(error: CombinatoricsError) -> Self {
        BetaBinomialError::Combinatorics(error)
    }
}

/// The two shape parameters and a number of trials, validated once at construction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BetaBinomialDistribution {
    alpha: f64,
    beta: f64,
    n: i32,
}

impl BetaBinomialDistribution {
    /// The reference's constructor, with its three `ParamUtils` checks in order.
    pub fn new(alpha: f64, beta: f64, n: i32) -> Result<Self, BetaBinomialError> {
        // `isPositive` is `> 0`, so it refuses NaN as well as zero.
        if alpha.is_nan() || alpha <= 0.0 {
            return Err(BetaBinomialError::AlphaNotPositive { alpha });
        }
        if beta.is_nan() || beta <= 0.0 {
            return Err(BetaBinomialError::BetaNotPositive { beta });
        }
        if n < 0 {
            return Err(BetaBinomialError::TrialsNegative { n });
        }
        Ok(Self { alpha, beta, n })
    }

    /// `logProbability(k)`.
    pub fn log_probability(&self, k: i32) -> Result<f64, BetaBinomialError> {
        if k < 0 {
            return Err(BetaBinomialError::SuccessesNegative { k });
        }
        if k > self.n {
            return Ok(f64::NEG_INFINITY);
        }
        // `k + alpha` and `n - k + beta`, the integer part summed as an integer first.
        Ok(binomial_coefficient_log(self.n as i64, k as i64)?
            + log_beta(k as f64 + self.alpha, (self.n - k) as f64 + self.beta)?
            - log_beta(self.alpha, self.beta)?)
    }

    /// `probability(k)`, which is `Math.exp` of the log.
    ///
    /// Under decision 0014 `Math.exp` cannot be transcribed, so this is the platform's `exp` and is
    /// not asserted bit-exact anywhere. `log_probability` is the one the clustering model uses.
    pub fn probability(&self, k: i32) -> Result<f64, BetaBinomialError> {
        Ok(self.log_probability(k)?.exp())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat() -> BetaBinomialDistribution {
        BetaBinomialDistribution::new(1.0, 1.0, 10).expect("valid")
    }

    #[test]
    fn the_flat_cluster_lands_on_three_different_doubles() {
        assert_eq!(
            flat().log_probability(0).expect("valid"),
            -2.3978952727983707
        );
        assert_eq!(
            flat().log_probability(1).expect("valid"),
            -2.3978952727983716
        );
        assert_eq!(flat().log_probability(5).expect("valid"), -2.39789527279837);
        // And `k = n` comes back to where `k = 0` was.
        assert_eq!(
            flat().log_probability(10).expect("valid"),
            -2.3978952727983707
        );
    }

    // The last value here is one digit short of `LN_2` and is a different double from it, which is
    // the whole point of asserting it, so the lint that wants the constant is switched off.
    #[allow(clippy::approx_constant)]
    #[test]
    fn the_high_af_cluster_rises_with_the_alternate_count() {
        let high = BetaBinomialDistribution::new(10.0, 1.0, 10).expect("valid");
        assert_eq!(high.log_probability(0).expect("valid"), -12.126791314602453);
        assert_eq!(high.log_probability(1).expect("valid"), -9.824206221608407);
        assert_eq!(high.log_probability(5).expect("valid"), -4.5248893547272875);
        // Not `-0.6931471805599453`: one digit shorter, and a different double.
        assert_eq!(high.log_probability(10).expect("valid"), -0.693147180559945);
    }

    #[test]
    fn a_count_past_the_total_is_negative_infinity_and_a_negative_one_is_a_refusal() {
        assert_eq!(
            flat().log_probability(11).expect("valid"),
            f64::NEG_INFINITY
        );
        assert_eq!(
            flat().log_probability(-1),
            Err(BetaBinomialError::SuccessesNegative { k: -1 })
        );
    }

    #[test]
    fn the_constructor_checks_in_the_reference_s_order() {
        assert_eq!(
            BetaBinomialDistribution::new(0.0, 1.0, 10),
            Err(BetaBinomialError::AlphaNotPositive { alpha: 0.0 })
        );
        assert_eq!(
            BetaBinomialDistribution::new(1.0, -1.0, 10),
            Err(BetaBinomialError::BetaNotPositive { beta: -1.0 })
        );
        assert_eq!(
            BetaBinomialDistribution::new(1.0, 1.0, -1),
            Err(BetaBinomialError::TrialsNegative { n: -1 })
        );
        // A bad alpha is reported even when the trial count is bad too, and NaN is a bad alpha.
        assert!(matches!(
            BetaBinomialDistribution::new(f64::NAN, 1.0, -1),
            Err(BetaBinomialError::AlphaNotPositive { alpha }) if alpha.is_nan()
        ));
    }
}
