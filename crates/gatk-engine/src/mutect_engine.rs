//! `Mutect2FilteringEngine`'s static arithmetic, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine`
//! (GATK 4.6.2.0).
//!
//! Three functions and two constants that every filter in the engine goes through.
//!
//! # The posterior is a two-entry normalisation, and the second entry is the answer
//!
//! ```java
//! final double[] unweightedPosteriorOfRealAndError = new double[] {logOddsOfRealVersusError + logPriorOfReal,
//!         NaturalLogUtils.log1mexp(logPriorOfReal)};
//! final double[] posteriorOfRealAndError = NaturalLogUtils.normalizeFromLogToLinearSpace(unweightedPosteriorOfRealAndError);
//! return posteriorOfRealAndError[1];
//! ```
//!
//! The value returned is the probability of **error**, which is why it falls as the odds of being
//! real rise. A prior of one makes it impossible, `log1mexp(0)` being negative infinity — except at
//! log odds of `-Infinity`, where both entries are `-Infinity` and the answer is `NaN`. A prior of
//! zero against log odds of `+Infinity` is the other end of that: the sum is `NaN` and the
//! normaliser refuses it.
//!
//! # Rounding the error is a clamp, and it keeps NaN
//!
//! ```java
//! public static double roundFinitePrecisionErrors(final double probability) {
//!     return Math.max(Math.min(probability, 1.0), 0.0);
//! }
//! ```
//!
//! Not a rounding: a probability of `1.0000000001` becomes `1`, `-Infinity` becomes `0`, and **NaN
//! stays NaN**, because `Math.min` and `Math.max` in Java propagate it. Rust's `f64::min` and
//! `f64::max` do the opposite — they return the non-NaN operand — so the obvious clamp answers
//! `1.0` for NaN and would call every such probability a certain artifact, where the reference
//! leaves it NaN for the caller to deal with.

use crate::allele_likelihoods::log10_to_log;
use crate::natural_log_utils::{log1mexp, normalize_from_log_to_linear_space, NonFiniteSum};

/// `EPSILON`.
pub const EPSILON: f64 = 1.0e-10;

/// `MIN_REPORTABLE_ERROR_PROBABILITY`.
pub const MIN_REPORTABLE_ERROR_PROBABILITY: f64 = 0.1;

/// `posteriorProbabilityOfError`: the probability that a call of these odds is an error.
pub fn posterior_probability_of_error(
    log_odds_of_real_versus_error: f64,
    log_prior_of_real: f64,
) -> Result<f64, NonFiniteSum> {
    let unweighted = [
        log_odds_of_real_versus_error + log_prior_of_real,
        log1mexp(log_prior_of_real),
    ];
    let posterior = normalize_from_log_to_linear_space(&unweighted)?;
    // The second entry: the error, not the truth.
    Ok(posterior[1])
}

/// `roundFinitePrecisionErrors`: `max(min(p, 1), 0)` with **Java's** NaN semantics.
pub fn round_finite_precision_errors(probability: f64) -> f64 {
    // `Math.min` and `Math.max` return NaN if either argument is NaN; `f64::min`/`f64::max` do not.
    if probability.is_nan() {
        return probability;
    }
    // `clamp` would do here, but the reference is two calls and the guard above is what makes them
    // equivalent; written as `clamp` the guard reads like a redundant special case.
    #[allow(clippy::manual_clamp)]
    probability.min(1.0).max(0.0)
}

/// `getTumorLogOdds`: `TLOD` is written as log10 odds and read as natural-log odds.
///
/// `None` in, `None` out: a record with no annotation answers null rather than an empty array, and
/// every caller has to tell the two apart.
pub fn tumor_log_odds(tumor_log_10_odds: Option<&[f64]>) -> Option<Vec<f64>> {
    tumor_log_10_odds.map(|odds| odds.iter().map(|value| log10_to_log(*value)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_posterior_falls_as_the_odds_rise() {
        let prior = -1.0;
        let low = posterior_probability_of_error(0.0, prior).expect("finite");
        let high = posterior_probability_of_error(10.0, prior).expect("finite");
        assert!(high < low);
        assert_eq!(low, 0.6321205588285577);
        assert_eq!(
            posterior_probability_of_error(10.0, -10.0).expect("finite"),
            0.4999886497599093
        );
        assert_eq!(
            posterior_probability_of_error(100.0, -100.0).expect("finite"),
            0.5
        );
    }

    #[test]
    fn a_prior_of_one_makes_the_error_impossible_except_at_minus_infinity() {
        for odds in [-10.0, -1.0, 0.0, 1.0, 10.0, 100.0, f64::INFINITY] {
            assert_eq!(
                posterior_probability_of_error(odds, 0.0).expect("finite"),
                0.0,
                "{odds}"
            );
        }
        // Both entries are -Infinity, and the normalisation of that is not a number.
        assert!(posterior_probability_of_error(f64::NEG_INFINITY, 0.0)
            .expect("no refusal")
            .is_nan());
    }

    #[test]
    fn infinite_odds_against_a_prior_of_zero_are_refused() {
        let error = posterior_probability_of_error(f64::INFINITY, f64::NEG_INFINITY)
            .expect_err("the sum is not a number");
        assert_eq!(error.class(), "java.lang.IllegalArgumentException");
        assert_eq!(
            error.message(),
            "logValues must be non-infinite and non-NAN"
        );
        // The other direction is fine: an impossible call against an impossible prior is certain.
        assert_eq!(
            posterior_probability_of_error(f64::NEG_INFINITY, f64::NEG_INFINITY).expect("finite"),
            1.0
        );
    }

    #[test]
    fn rounding_the_error_is_a_clamp_that_keeps_nan() {
        assert_eq!(round_finite_precision_errors(-0.5), 0.0);
        assert_eq!(round_finite_precision_errors(1.0 + 1.0e-10), 1.0);
        assert_eq!(round_finite_precision_errors(0.5), 0.5);
        assert_eq!(round_finite_precision_errors(f64::NEG_INFINITY), 0.0);
        assert_eq!(round_finite_precision_errors(f64::INFINITY), 1.0);
        // The one that a Rust clamp would get wrong.
        assert!(round_finite_precision_errors(f64::NAN).is_nan());
        #[allow(clippy::manual_clamp)]
        let naive = f64::NAN.min(1.0).max(0.0);
        assert_eq!(
            naive, 1.0,
            "the naive clamp calls a NaN probability a certain artifact"
        );
    }

    #[test]
    fn the_tumour_odds_are_converted_and_a_missing_annotation_is_not_empty() {
        assert_eq!(tumor_log_odds(Some(&[6.0])), Some(vec![13.815510557964275]));
        assert_eq!(
            tumor_log_odds(Some(&[0.0, -3.0])),
            Some(vec![0.0, -6.907755278982138])
        );
        assert_eq!(tumor_log_odds(None), None);
        // An empty annotation is not the same thing as none at all.
        assert_eq!(tumor_log_odds(Some(&[])), Some(Vec::new()));
    }

    #[test]
    fn the_two_constants_are_the_references() {
        assert_eq!(EPSILON, 1.0e-10);
        assert_eq!(MIN_REPORTABLE_ERROR_PROBABILITY, 0.1);
    }
}
