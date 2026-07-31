//! The handful of `org.broadinstitute.hellbender.utils.MathUtils` and `QualityUtils` entries the
//! annotations reach through (GATK 4.6.2.0).
//!
//! These are not a general numerics library. They are the exact functions the ported call sites
//! name, written so the **order of the floating-point operations** matches the reference's, because
//! that order is what the golden measures.
//!
//! # `normalizeFromLog10ToLinearSpace` subtracts a log-sum, it does not divide by a sum
//!
//! ```java
//! final double log10Sum = log10SumLog10(array);
//! final double[] result = applyToArrayInPlace(array, x -> x - log10Sum);
//! return takeLog10OfOutput ? result : applyToArrayInPlace(result, x -> Math.pow(10.0, x));
//! ```
//!
//! So every element goes through `Math.pow(10.0, x)` individually. A "normalised" vector here is
//! therefore not guaranteed to sum to exactly one, and the callers that then compare two of its
//! entries for equality (`GenotypeUtils.computeDiploidGenotypeCounts` does, twice) are comparing
//! the results of two separate `pow` calls.
//!
//! # `Math.pow` is still the deferred function
//!
//! [`pow10`] is the platform's `powf`, as `fisher_exact` already does, and for the same reason:
//! `Math.pow` is a HotSpot intrinsic that decision 0007 (in `htsjdk-rs`) records as unported. The
//! conformance golden is what decides whether the two agree on the inputs these annotations
//! produce, rather than an assumption made here.

/// `MathUtils.maxElementIndex(array, start, endIndex)`.
///
/// Strictly greater, so the **first** maximum wins when two entries tie. That tie is reachable:
/// a genotype with PLs `[0, 0, X]` normalises to two equal likelihoods.
pub fn max_element_index(array: &[f64], start: usize, end: usize) -> usize {
    let mut max_i = start;
    for i in (start + 1)..end {
        if array[i] > array[max_i] {
            max_i = i;
        }
    }
    max_i
}

/// `Math.pow(10.0, x)`. See the module note: the deferred intrinsic, measured rather than assumed.
pub fn pow10(x: f64) -> f64 {
    10.0f64.powf(x)
}

/// `MathUtils.log10sumLog10(array)`.
///
/// The sum is `1.0 + sum over i != max of pow(10, a[i] - max)`, left to right with the maximum's
/// own term added as a literal zero, so the accumulation order is the array's.
pub fn log10_sum_log10(log10p: &[f64]) -> f64 {
    log10_sum_log10_range(log10p, 0, log10p.len())
}

/// `MathUtils.log10sumLog10(array, start, finish)`.
pub fn log10_sum_log10_range(log10p: &[f64], start: usize, finish: usize) -> f64 {
    if finish.saturating_sub(start) < 2 {
        return if finish == start {
            f64::NEG_INFINITY
        } else {
            log10p[start]
        };
    }
    let max_index = max_element_index(log10p, start, finish);
    let max_value = log10p[max_index];
    if max_value == f64::NEG_INFINITY {
        return max_value;
    }
    let mut sum = 0.0;
    for (i, value) in log10p.iter().enumerate().take(finish).skip(start) {
        // The maximum contributes a literal 0.0 rather than being skipped, so the number of
        // additions is the range's length whatever the maximum is.
        sum += if i == max_index {
            0.0
        } else {
            pow10(value - max_value)
        };
    }
    let sum = 1.0 + sum;
    // `Utils.validateArg(!isNaN(sum) && sum != POSITIVE_INFINITY)`: an IllegalArgumentException
    // there, which no ported caller can reach because the inputs are PLs.
    max_value + jmath::math::log10(sum)
}

/// `MathUtils.normalizeFromLog10ToLinearSpace`.
pub fn normalize_from_log10_to_linear_space(array: &[f64]) -> Vec<f64> {
    let log10_sum = log10_sum_log10(array);
    array.iter().map(|x| pow10(x - log10_sum)).collect()
}

/// `MathUtils.normalizeSumToOne`, which divides by the sum whatever the sum is.
///
/// There is no guard for a zero sum, so an all-zero input yields `NaN` rather than zeros. The
/// negative-sum check is an `IllegalArgumentException`, modelled here as `None`.
pub fn normalize_sum_to_one(array: &[f64]) -> Option<Vec<f64>> {
    if array.is_empty() {
        return Some(Vec::new());
    }
    let sum: f64 = array.iter().sum();
    if sum < 0.0 {
        return None;
    }
    Some(array.iter().map(|x| x / sum).collect())
}

/// `MathUtils.fastRound`, which is a **truncating cast**, not a rounding function.
///
/// ```java
/// return (d > 0.0) ? (int) (d + 0.5d) : (int) (d - 0.5d);
/// ```
///
/// Half-away-from-zero for ordinary input, but the cast truncates towards zero, so it rounds twice
/// on a value whose sum with a half is not representable, exactly the way `Math.round` stopped
/// doing in Java 7. Java's `(int)` narrowing clamps rather than wrapping and answers zero for NaN.
pub fn fast_round(d: f64) -> i32 {
    let shifted = if d > 0.0 { d + 0.5 } else { d - 0.5 };
    if shifted.is_nan() {
        return 0;
    }
    if shifted >= i32::MAX as f64 {
        return i32::MAX;
    }
    if shifted <= i32::MIN as f64 {
        return i32::MIN;
    }
    shifted as i32
}

/// `Math.min(double, double)`, which is **not** `f64::min`.
///
/// Rust's `min` returns the non-NaN argument; Java's propagates the NaN. It also distinguishes the
/// two zeros, which `f64::min` explicitly leaves unspecified. Both differences are reachable from
/// `Math.max(0., Math.min(1., pval))` in `ExcessHet`, where a NaN p-value would silently become 1.
pub fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    if a == 0.0 && b == 0.0 && b.is_sign_negative() {
        return b;
    }
    if a <= b {
        a
    } else {
        b
    }
}

/// `Math.max(double, double)`. See [`java_min`].
pub fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    if a == 0.0 && b == 0.0 && a.is_sign_negative() {
        return b;
    }
    if a >= b {
        a
    } else {
        b
    }
}

/// `QualityUtils.qualToProb(double)`: `1 - pow(10, qual / -10)`.
///
/// The `double` overload, not the cached `byte` one: an `int` GQ widens to `double` at the call
/// site in `GenotypeUtils`, so the cache is bypassed and the `Math.pow` is real.
pub fn qual_to_prob(qual: f64) -> f64 {
    1.0 - pow10(qual / -10.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_first_maximum_wins_a_tie() {
        assert_eq!(max_element_index(&[1.0, 1.0, 0.0], 0, 3), 0);
    }

    #[test]
    fn fast_round_is_half_away_from_zero() {
        assert_eq!(fast_round(2.5), 3);
        assert_eq!(fast_round(-2.5), -3);
        assert_eq!(fast_round(0.4), 0);
        // The truncating cast is what makes this not `Math.round`: the double below a half plus a
        // half is exactly one, so it rounds twice and answers one where half-up answers zero.
        assert_eq!(fast_round(0.499_999_999_999_999_94), 1);
    }

    #[test]
    fn an_all_zero_vector_normalises_to_nan() {
        let out = normalize_sum_to_one(&[0.0, 0.0]).expect("a non-negative sum");
        assert!(out.iter().all(|value| value.is_nan()));
    }
}
