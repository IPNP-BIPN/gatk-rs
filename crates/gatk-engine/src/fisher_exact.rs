//! `FisherExactTest` and the two `MathUtils` sums under it, ported from
//! `org.broadinstitute.hellbender.utils` (GATK 4.7.0.0; `FisherExactTest.java` is unchanged since
//! 4.6.2.0).
//!
//! The two-sided p-value of a 2x2 contingency table, which is what `FS`, the Fisher-strand
//! annotation, reports as a phred score.
//!
//! # It is R's algorithm, not a closed form
//!
//! ```java
//! //Note: this implementation follows the one in R base package
//! final double[] logds = support.mapToDouble(dist::logProbability);
//! final double threshold = logds[x[0][0] - lo] * REL_ERR;
//! final double[] log10ds = DoubleStream.of(logds).filter(d -> d <= threshold)...
//! ```
//!
//! Every point of the hypergeometric support is evaluated, the observed point's log density is
//! multiplied by `REL_ERR = 1 - 10e-7` to make a threshold, and every point at or below it is
//! summed. The relative slack is what keeps a point that ties the observed density inside the sum,
//! and it is applied by **multiplication on a negative number**, so the threshold sits slightly
//! *above* the observed density where the density is negative and slightly below it where it is
//! not.
//!
//! # The sum goes through base 10, and through `Math.pow`
//!
//! ```java
//! public static double sumLog10(final double[] log10values) {
//!     return Math.pow(10.0, log10SumLog10(Utils.nonNull(log10values)));
//! }
//! ```
//!
//! The densities arrive as natural logs, are converted to base 10 one by one, summed in log space
//! by a max-shifted `log10SumLog10`, and exponentiated back. `Math.pow` is the one function
//! decision 0007 deferred, so this port reaches the platform's `powf` and the golden is what says
//! whether that agrees. Reordering the conversion and the sum would change the result even if
//! `pow` agreed.
//!
//! # A support of one point is exactly 1.0, without computing anything
//!
//! ```java
//! if (support.size() <= 1){     //special case, support has only one value
//!     return 1.0;
//! }
//! ```
//!
//! A table whose margins allow a single arrangement has a p-value of exactly 1.0, returned before
//! any density is evaluated, so no `log`, `exp` or `pow` rounding can reach it.
//!
//! # Reading this file without reading Rust
//!
//! `docs/COMMENTING.md` has the general table. The pieces this file leans on:
//!
//! - `&[f64]` is a borrowed, read-only view of a `double[]`; the caller keeps the array;
//! - `Vec<f64>` is a growable array that owns its elements, here standing in for `double[]`;
//! - `|x| expr` is a closure, Java's lambda `x -> expr`;
//! - a chain such as `.iter().filter(..).map(..).collect()` is Java's `DoubleStream` pipeline
//!   ending in `toArray()`;
//! - a function's last expression is its return value; `return` only marks an early exit.

use jmath::saddle_point::hypergeometric_log_probability;

/// `FisherExactTest.REL_ERR`.
///
/// ```java
/// private static final double REL_ERR = 1 - 10e-7;
/// ```
///
/// `10e-7` is ten times `1e-7`, so this is `1 - 1e-6` (0.999999), not `1 - 1e-7`. The literal is
/// kept exactly as the Java writes it. Java promotes the integer `1` to `1.0` before subtracting;
/// Rust does not convert implicitly, so it writes `1.0`, and the double is the same.
const REL_ERR: f64 = 1.0 - 10e-7;

/// `MathUtils.LOG10_E`, which is `Math.log10(Math.E)`.
///
/// The constant is computed, not taken from Rust's own `std::f64::consts::LOG10_E`: the reference
/// calls `Math.log10` on the double nearest e, and `jmath::math::log10` is jmath's correctly
/// rounded port of `Math.log10`, so this is the same call on the same argument. A literal rounded
/// from the true log10(e) need not be the same double. It is a function rather than a Rust
/// `const` because a `const` must be computable at compile time and the jmath function is not
/// declared that way.
fn log10_e() -> f64 {
    jmath::math::log10(std::f64::consts::E)
}

/// `MathUtils.logToLog10`: a natural logarithm converted to base ten by one multiplication.
///
/// ```java
/// return ln * LOG10_E;
/// ```
///
/// A multiplication by the rounded constant, kept as one: dividing by the natural log of ten
/// instead need not give the same double. Also called by `allele_frequency_calculator`.
pub fn log_to_log10(ln: f64) -> f64 {
    ln * log10_e()
}

/// `MathUtils.log10SumLog10(values)`: the max-shifted sum.
///
/// The shift means the largest term contributes exactly 1 and never goes through `pow`, so the
/// answer depends on which element is the maximum and not only on the multiset of values.
///
/// The Java is `MathUtils.log10SumLog10(log10Values, start, finish)`, reached with `start = 0`
/// and `finish = log10Values.length`. This is a second copy of
/// `math_utils::log10_sum_log10_range`, which carries the full commentary; the steps here are the
/// same in the same order, with the range fixed to the whole slice. The loop it matches:
///
/// ```java
/// double sum = 1.0;
/// for (int i = start; i < finish; i++) {
///     final double curVal = log10Values[i];
///     if (i == maxElementIndex || curVal == Double.NEGATIVE_INFINITY) {
///         continue;
///     } else {
///         final double scaled_val = curVal - maxValue;
///         sum += Math.pow(10.0, scaled_val);
///     }
/// }
/// ```
pub fn log10_sum_log10(values: &[f64]) -> f64 {
    // `if (start >= finish) return Double.NEGATIVE_INFINITY;` with `start = 0`.
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }
    // `MathUtils.maxElementIndex`, inlined. `.enumerate()` pairs each value with its position,
    // and `.fold(0, |best, (index, value)| ...)` walks the pairs carrying `best` from one to the
    // next, starting at 0: a loop with an accumulator, Java's `reduce` with an identity. The test
    // is strictly greater, as in the reference, so the first of two equal maxima is kept.
    let max_index =
        values.iter().enumerate().fold(
            0usize,
            |best, (index, value)| {
                if *value > values[best] {
                    index
                } else {
                    best
                }
            },
        );
    let max_value = values[max_index];
    if max_value == f64::NEG_INFINITY {
        return max_value;
    }
    // `let mut` marks a variable that may be reassigned; a plain `let` cannot be.
    let mut sum = 1.0;
    for (index, value) in values.iter().enumerate() {
        // `value` is a reference into the slice and `*value` reads the `f64` through it. The
        // subtraction below works on the reference directly, because Rust defines arithmetic on
        // references to numbers.
        if index == max_index || *value == f64::NEG_INFINITY {
            continue;
        }
        sum += pow10(value - max_value);
    }
    if sum.is_nan() || sum == f64::INFINITY {
        // `IllegalArgumentException("log10 p: Values must be non-infinite and non-NAN")`. Nothing
        // in the Fisher path can reach it, since every term is at most the maximum.
        return f64::NAN;
    }
    // The `if` is an expression whose value is added to `max_value`, where Java writes a ternary
    // on `sum != 1.0`. A sum still exactly 1.0 skips the logarithm and adds `0.0`.
    max_value
        + if sum != 1.0 {
            jmath::math::log10(sum)
        } else {
            0.0
        }
}

/// `MathUtils.sumLog10`: the sum of the numbers whose base-ten logs are given.
///
/// ```java
/// return Math.pow(10.0, log10SumLog10(Utils.nonNull(log10values)));
/// ```
///
/// `Utils.nonNull` has nothing to check: a Rust slice cannot be null.
pub fn sum_log10(values: &[f64]) -> f64 {
    pow10(log10_sum_log10(values))
}

/// `Math.pow(10.0, x)`.
///
/// Decision 0007 (htsjdk-rs) deferred `Math.pow` because HotSpot's intrinsic is seeded by the
/// approximate `rcpps` instruction and is not correctly rounded, so no independent implementation
/// is guaranteed to agree; and the intrinsic cannot be transcribed, the JDK being GPL2 (decision
/// 0023). Decision 0027 later bounded it at 1 ulp, which is still not exact.
/// This calls the platform's, and the conformance suite is what says whether that is the same
/// double on the inputs the strand-bias annotations produce.
///
/// `docs/numeric-functions-a-ported-call-site-reaches.md` lists this function as one of the
/// host-libm call sites, with the same status as `math_utils::pow10`: 99.9378% agreement with
/// `Math.pow`, unbounded per call, and the bounded `strict_math::pow` (decision 0027) measured
/// worse at that site.
fn pow10(x: f64) -> f64 {
    10.0f64.powf(x)
}

/// `FisherExactTest.twoSidedPValue(normalizedTable)`.
///
/// The table is `[[a, b], [c, d]]`, and the reference refuses anything that is not 2x2.
/// The Java does that with `Utils.nonNull` and three `Utils.validateArg` checks, each failing
/// with `"input must be 2x2 "`. They have no counterpart here because the type does the work:
/// `[[i32; 2]; 2]` is an array of exactly two arrays of exactly two `int`s, fixed at compile
/// time and never null, so no other shape can be passed. The one caller, `strand_bias`, passes
/// the table through `normalize_contingency_table` first, as the Java's javadoc expects ("a
/// normalized table that ensures that the sum of all four entries is less than 2 * 200").
pub fn two_sided_p_value(table: [[i32; 2]; 2]) -> f64 {
    // `m` and `n` are the two row sums and `k` the first column's sum. `x[0][0]` follows a
    // hypergeometric law over `m + n` items with `m` successes and `k` draws.
    let m = table[0][0] + table[0][1];
    let n = table[1][0] + table[1][1];
    let k = table[0][0] + table[1][0];
    // `Math.max(0, k - n)` and `Math.min(k, m)`, written as methods on the integer.
    let lo = (k - n).max(0);
    let hi = k.min(m);

    // `IndexRange(lo, hi + 1)`: a support of one point is exactly 1.0, with nothing computed.
    if (hi + 1 - lo) <= 1 {
        return 1.0;
    }

    // `new HypergeometricDistribution(null, m+n, m, k)` then
    // `support.mapToDouble(dist::logProbability)`. The jmath function takes the same three
    // parameters in the same order (population, successes, sample) plus the point. `lo..=hi`
    // includes `hi`: the same points as the half-open `IndexRange(lo, hi + 1)`, in ascending
    // order, so `log_densities[i]` is the density at `lo + i`, as in the Java's `logds`.
    let log_densities: Vec<f64> = (lo..=hi)
        .map(|x| hypergeometric_log_probability(m + n, m, k, x))
        .collect();

    // The observed point's density, slackened by a relative error. Multiplying a negative log
    // density by a factor below 1 *raises* it, which is what keeps ties inside the sum.
    //
    // `(table[0][0] - lo) as usize` converts the `int` offset to Rust's index type. `as` does not
    // check: a negative value would wrap to a huge index and the bounds check would then panic,
    // where Java throws `ArrayIndexOutOfBoundsException`. Neither can happen for non-negative
    // counts, since `x[0][0] >= 0` and `x[0][0] = k - x[1][0] >= k - n`.
    let threshold = log_densities[(table[0][0] - lo) as usize] * REL_ERR;
    // `filter(d -> d <= threshold).map(MathUtils::logToLog10)`: compared in natural-log space,
    // converted to base ten only after. `**density` dereferences twice because `.iter()` yields
    // references to the elements and `.filter` hands its closure a reference to each of those.
    let log10_densities: Vec<f64> = log_densities
        .iter()
        .filter(|density| **density <= threshold)
        .map(|density| log_to_log10(*density))
        .collect();

    let p_value = sum_log10(&log10_densities);
    // "min is necessary as numerical precision can result in pValue being slightly greater than 1"
    //
    // `f64::min` is not `Math.min` (see `math_utils::java_min`): given a NaN it returns the other
    // argument where Java returns the NaN. No NaN reaches this line: a NaN density fails the
    // `<=` filter above, and `log10_sum_log10` returns NaN only for a NaN or infinite sum, which
    // terms that are each at most one, one per point of an `int` range, cannot make.
    p_value.min(1.0)
}

// `#[cfg(test)]` compiles this module only for `cargo test`; none of it ships.
#[cfg(test)]
mod tests {
    // Brings the parent module's items into scope, so the tests call them by their bare names.
    use super::*;

    #[test]
    fn a_support_of_one_point_is_exactly_one() {
        // Every read on one strand and one allele: the support cannot move.
        assert_eq!(two_sided_p_value([[0, 0], [0, 5]]), 1.0);
    }

    #[test]
    fn a_balanced_table_is_near_one_and_a_skewed_one_is_small() {
        let balanced = two_sided_p_value([[10, 10], [10, 10]]);
        let skewed = two_sided_p_value([[20, 0], [0, 20]]);
        assert!(balanced > 0.9, "balanced = {balanced}");
        assert!(skewed < 1e-8, "skewed = {skewed}");
    }
}
