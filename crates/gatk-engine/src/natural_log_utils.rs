//! Ported from `org.broadinstitute.hellbender.utils.NaturalLogUtils` (GATK 4.6.2.0).
//!
//! The natural-log arithmetic every somatic likelihood goes through, and the layer that decides
//! whether `AllelePseudoDepth` (G1.9, #96) can be reproduced at all.
//!
//! # Which `exp` this is
//!
//! The reference calls `Math.exp`, whose only exact port is a transcription of GPL2-only HotSpot
//! source (htsjdk-rs decision 0014, and htsjdk-rs #71 for why no route round it is open). This
//! port calls [`jmath::strict_math::exp`], which is FDLIBM: permissively licensed, exact for
//! `StrictMath.exp`, and measured at **1 ulp** worst-case against `Math.exp` (htsjdk-rs decision
//! 0025).
//!
//! So this module is **not** bit-identical to the reference on every input, and saying otherwise
//! would be the mistake this programme keeps catching. What it is, is bounded: every value it
//! produces differs from the reference's by at most the accumulated effect of a 1-ulp `exp`. The
//! suite measures that effect rather than assuming it, and G1.9's whole argument is that for
//! `AllelePseudoDepth` the effect cannot reach the output, because the output is rounded to two
//! decimals first.
//!
//! # `logSumExp` has a path that is exact whatever `exp` does
//!
//! ```java
//! double sum = 1.0;
//! for (int i = 0; i < logValues.length; i++) {
//!     if (i == maxElementIndex || curVal == Double.NEGATIVE_INFINITY) { continue; }
//!     sum += Math.exp(curVal - maxValue);
//! }
//! return maxValue + (sum != 1.0 ? Math.log(sum) : 0.0);
//! ```
//!
//! The accumulator starts at **1.0**, not at zero, because the maximum's own term — `exp(0)` — is
//! skipped in the loop and folded in as that 1. Two consequences, and the second is the useful one:
//!
//!  * a `-Infinity` entry is skipped rather than contributing `exp(-inf) = 0`, which is the same
//!    value by a different route and is the reason the guard exists at all;
//!  * `sum != 1.0` then skips the `log`. So an array with **one** non-infinite entry, or with a
//!    single maximum and every other entry at `-Infinity`, returns `maxValue` **untouched** — no
//!    `exp`, no `log`, bit-identical to the reference by construction rather than by measurement.
//!
//! # The refusal is on the accumulator, not on the inputs
//!
//! `IllegalArgumentException("logValues must be non-infinite and non-NAN")` is thrown after the
//! loop, on `sum`, so it fires for a `NaN` **input** (which poisons the sum) and for an input large
//! enough to overflow the accumulator — but **not** for `-Infinity`, which the loop skipped, and
//! not for an all-`-Infinity` array, which returned early. The message names the inputs; the test
//! is on their sum.

use jmath::math::log;
use jmath::strict_math::exp;

use crate::math_utils::max_element_index;

/// `LOG_ONE_HALF`, which is `Math.log(0.5)` and not a literal.
pub fn log_one_half() -> f64 {
    log(0.5)
}

/// `LOG1MEXP_THRESHOLD`, the same value, used for a different purpose.
fn log1mexp_threshold() -> f64 {
    log(0.5)
}

/// `IllegalArgumentException` from [`log_sum_exp`], kept as a value so a caller can report the
/// class the reference would have thrown.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NonFiniteSum;

impl NonFiniteSum {
    pub fn class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }

    pub fn message(&self) -> &'static str {
        "logValues must be non-infinite and non-NAN"
    }
}

/// `logSumExp(double...)`.
///
/// The operation order is the reference's: the maximum is found first, every other term is scaled
/// by it before being exponentiated, and the accumulator is summed in index order. Floating-point
/// addition is not associative, so reordering the loop would change the result even with an
/// identical `exp`.
pub fn log_sum_exp(log_values: &[f64]) -> Result<f64, NonFiniteSum> {
    // `maxElementIndex` refuses an empty array; the port's own version does the same, so an empty
    // slice would panic there rather than reach this.
    let max_index = max_element_index(log_values, 0, log_values.len());
    let max_value = log_values[max_index];
    if max_value == f64::NEG_INFINITY {
        return Ok(max_value);
    }

    let mut sum = 1.0;
    for (index, &value) in log_values.iter().enumerate() {
        if index == max_index || value == f64::NEG_INFINITY {
            continue;
        }
        sum += exp(value - max_value);
    }

    if sum.is_nan() || sum == f64::INFINITY {
        return Err(NonFiniteSum);
    }
    // The guard is what makes a single-term array exact: no `log` is called at all.
    Ok(max_value + if sum != 1.0 { log(sum) } else { 0.0 })
}

/// `normalizeLog(array, takeLogOfOutput, inPlace)`.
///
/// `in_place` is not modelled: it decides whether the reference writes into the caller's array, and
/// this port returns a new vector either way. The distinction is observable only through aliasing,
/// which Rust does not permit here, and the *values* are identical on both paths.
pub fn normalize_log(array: &[f64], take_log_of_output: bool) -> Result<Vec<f64>, NonFiniteSum> {
    let log_sum = log_sum_exp(array)?;
    let result: Vec<f64> = array.iter().map(|value| value - log_sum).collect();
    if take_log_of_output {
        Ok(result)
    } else {
        Ok(result.into_iter().map(exp).collect())
    }
}

/// `normalizeFromLogToLinearSpace(array)`: `normalizeLog(array, false, true)`.
///
/// Every element goes through `exp` here, so unlike [`log_sum_exp`] there is no exact path: this
/// function is 1-ulp-bounded rather than bit-identical, on every input.
pub fn normalize_from_log_to_linear_space(array: &[f64]) -> Result<Vec<f64>, NonFiniteSum> {
    normalize_log(array, false)
}

/// `posteriors(logPriors, logLikelihoods)`, which is what the Dirichlet fixed point calls per read.
///
/// `MathArrays.ebeAdd` requires equal lengths; the reference reaches a
/// `DimensionMismatchException` otherwise, and this returns `None` rather than modelling an
/// exception nothing in G1.9 can produce.
pub fn posteriors(log_priors: &[f64], log_likelihoods: &[f64]) -> Option<Vec<f64>> {
    if log_priors.len() != log_likelihoods.len() {
        return None;
    }
    let summed: Vec<f64> = log_priors
        .iter()
        .zip(log_likelihoods)
        .map(|(prior, likelihood)| prior + likelihood)
        .collect();
    normalize_from_log_to_linear_space(&summed).ok()
}

/// `log1mexp(a)`: `log(1 - exp(a))` without losing precision.
///
/// ```java
/// if (a > 0) return Double.NaN;
/// if (a == 0) return Double.NEGATIVE_INFINITY;
/// return (a < LOG1MEXP_THRESHOLD) ? Math.log1p(-Math.exp(a)) : Math.log(-Math.expm1(a));
/// ```
///
/// The branch is the whole function, and it is ported with the threshold intact even though nothing
/// in G1.9 reaches it: the two arms are accurate in different halves of the range, and collapsing
/// them to either one alone loses precision exactly where the other was chosen to keep it.
///
/// `log1p` and `expm1` are the host libm's, which are **not** exact against the JVM (their rates
/// are in the jmath table). A call site that reaches this function inherits that, not the 1-ulp
/// bound `exp` carries.
pub fn log1mexp(a: f64) -> f64 {
    if a > 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return f64::NEG_INFINITY;
    }
    if a < log1mexp_threshold() {
        (-exp(a)).ln_1p()
    } else {
        log(-(a.exp_m1()))
    }
}
