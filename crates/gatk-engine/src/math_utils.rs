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
//!
//! # Reading this file without reading Rust
//!
//! Every function below is commented to the standard in `docs/COMMENTING.md`. Three pieces of
//! syntax recur and are worth learning once:
//!
//! - `&[f64]` is a **borrowed view** of a list of 64-bit floating-point numbers. "Borrowed" means
//!   the function may read it but does not own it and will not free it; the caller keeps the list.
//!   Java would write `double[]`, with the difference that Java's array could be modified through
//!   the reference and this one cannot, because it is not marked `mut` (mutable);
//! - `f64` is Java's `double`, `i32` is Java's `int`, `i64` is Java's `long`, `usize` is an
//!   unsigned integer wide enough to index memory (Java has no equivalent and uses `int`);
//! - a function's last expression is its return value; `return` is only written for early exits.

/// `MathUtils.maxElementIndex(array, start, endIndex)`.
///
/// **What**: the index of the largest value in `array[start..end]`.
///
/// **How**: walk once from `start + 1`, keeping the index of the best seen so far.
///
/// **Why it is written out rather than using a library maximum**: the comparison is strictly
/// greater-than, so the **first** of two equal maxima wins. A library routine is free to return
/// either. That tie is reachable: a genotype with PLs `[0, 0, X]` normalises to two equal
/// likelihoods, and which index is reported decides which pair of columns
/// `computeDiploidGenotypeCounts` reads.
pub fn max_element_index(array: &[f64], start: usize, end: usize) -> usize {
    // `let mut` declares a variable that may be reassigned. Without `mut` a Rust binding is
    // constant, which is the opposite default from Java's.
    let mut max_i = start;
    // `(a..b)` is a half-open range: it yields a, a+1, ... b-1, never b. Java's
    // `for (int i = a; i < b; i++)`.
    for i in (start + 1)..end {
        // Strictly greater, never `>=`. See the note above: this is the tie-breaking rule.
        if array[i] > array[max_i] {
            max_i = i;
        }
    }
    // No `return` keyword: the last expression of the function body is what it evaluates to.
    max_i
}

/// `Math.pow(10.0, x)`. See the module note: the deferred intrinsic, measured rather than assumed.
///
/// **What**: ten raised to the power `x`.
///
/// **How**: the platform's own `powf`. `10.0f64` writes the literal ten as a 64-bit float so that
/// the method call resolves to the 64-bit version rather than the 32-bit one.
///
/// **Why it is a named function rather than inlined at each call site**: so that there is exactly
/// one place to change if decision 0007 is ever closed by porting `Math.pow`, and so that a reader
/// grepping for the deferred function finds every use of it.
pub fn pow10(x: f64) -> f64 {
    10.0f64.powf(x)
}

/// `MathUtils.log10SumLog10(array)`, with a **capital S**.
///
/// **What**: given an array of base-ten logarithms, the logarithm of the sum of the numbers they
/// represent, computed without ever forming those numbers (which would overflow or underflow).
///
/// **How**: factor out the largest element. If `m` is the largest, then
/// `log10(sum of 10^a_i) = m + log10(sum of 10^(a_i - m))`, and every term of that inner sum is at
/// most one, so nothing overflows. The term for `m` itself is exactly one, which is why the
/// accumulator starts at `1.0` rather than at zero.
///
/// **Why the capital S matters**: `MathUtils` has two log-sums whose names differ by that one
/// letter, and they are not the same function:
///
/// | | accumulation | `-Infinity` entries | one element |
/// |---|---|---|---|
/// | `log10sumLog10` | `1.0 + (sum of terms)` | contribute `pow(10, -inf) = 0` | returned as is |
/// | `log10SumLog10` | `sum = 1.0` then `sum += term` | **skipped**, no addition at all | still summed |
///
/// The accumulation order is observable. For PLs of `[60, 0, 60]` the terms are `1e-6` twice, and
/// `(1.0 + 1e-6) + 1e-6` is two ulp away from `1.0 + (1e-6 + 1e-6)`. ("ulp" is one unit in the last
/// place: the distance between a floating-point number and the next one representable.)
/// `normalizeLog10`, and so every genotype count `ExcessHet` and `InbreedingCoeff` rest on, calls
/// the capital-S one. The golden caught the port calling the other, on the `equilibrium` cohort's
/// het count and nowhere else.
///
/// The last line is the third difference: a sum still exactly `1.0` skips `Math.log10` rather than
/// taking the logarithm of one, so a single-element array never reaches the logarithm at all.
pub fn log10_sum_log10(log10_values: &[f64]) -> f64 {
    // `.len()` is the number of elements, so this asks for the whole array.
    log10_sum_log10_range(log10_values, 0, log10_values.len())
}

/// `MathUtils.log10SumLog10(array, start, finish)`.
///
/// **What**: as [`log10_sum_log10`], but over the half-open slice `[start, finish)` only.
///
/// **How and why**: see [`log10_sum_log10`]. The three early exits below are each a branch of the
/// reference, reproduced in its order because the order decides which one fires first.
pub fn log10_sum_log10_range(log10_values: &[f64], start: usize, finish: usize) -> f64 {
    // An empty or inverted range. The reference answers negative infinity, which is the logarithm
    // of a sum of nothing, and it does this test **before** looking at the array, so an
    // out-of-bounds `start` never causes an index error on this path.
    if start >= finish {
        return f64::NEG_INFINITY;
    }
    let max_index = max_element_index(log10_values, start, finish);
    let max_value = log10_values[max_index];
    // If the largest entry is negative infinity then every entry is, so every number is zero and
    // the sum is zero, whose logarithm is negative infinity. Returning `max_value` rather than the
    // constant is the reference's own wording and gives the same bits.
    if max_value == f64::NEG_INFINITY {
        return max_value;
    }
    // Starts at one, not zero: the largest element's own term is `10^(m - m)`, which is one, and
    // the loop below skips it rather than adding it.
    let mut sum = 1.0f64;
    // `.iter()` walks the array by reference; `.enumerate()` pairs each item with its position;
    // `.take(finish)` stops after position `finish - 1`; `.skip(start)` discards the first `start`
    // pairs. Together they are Java's `for (int i = start; i < finish; i++)`, and the positions
    // reported by `enumerate` are absolute because `skip` comes after it.
    for (i, value) in log10_values.iter().enumerate().take(finish).skip(start) {
        // Two skips, not one. The maximum's term is omitted because it was folded into the
        // starting `1.0`; a negative-infinity term is omitted because the reference omits it, and
        // that is **not** the same as adding the zero it would have produced. Adding zero is an
        // operation and can change nothing here, but the count of additions changes which
        // intermediate sums exist, and this is the file where that matters.
        if i == max_index || *value == f64::NEG_INFINITY {
            continue;
        }
        // `*value` reads through the borrow: `value` is a reference to an `f64`, `*value` is the
        // `f64` itself. Java has no equivalent because Java's `double` is never a reference.
        sum += pow10(value - max_value);
    }
    // `throw new IllegalArgumentException("log10 p: Values must be non-infinite and non-NAN")`,
    // which no ported caller can reach because the inputs are PLs.
    //
    // The `if` is an expression here, not a statement: it evaluates to one of the two branches and
    // that value is added to `max_value`. Java would need a ternary.
    max_value
        + if sum != 1.0 {
            jmath::math::log10(sum)
        } else {
            0.0
        }
}

/// `MathUtils.normalizeFromLog10ToLinearSpace`.
///
/// **What**: turn an array of base-ten logarithms into the numbers themselves, scaled so they sum
/// to (approximately) one.
///
/// **How**: subtract the log-sum from every entry, then raise ten to each. Subtracting a logarithm
/// is dividing, so this divides every number by their total.
///
/// **Why "approximately"**: each element goes through its own `pow`, and a sum of separately
/// rounded quotients need not be exactly one. Callers that test two entries for equality are
/// comparing two `pow` results, not two exact fractions.
pub fn normalize_from_log10_to_linear_space(array: &[f64]) -> Vec<f64> {
    let log10_sum = log10_sum_log10(array);
    // `.map(...)` applies the closure to each element and `.collect()` gathers the results into a
    // new `Vec<f64>`, which is Java's `ArrayList<Double>` without the boxing. The closure
    // `|x| ...` takes one argument; `x` here is a reference, and `x - log10_sum` reads through it
    // automatically because subtraction is defined on references to numbers.
    array.iter().map(|x| pow10(x - log10_sum)).collect()
}

/// `MathUtils.normalizeSumToOne`, which divides by the sum whatever the sum is.
///
/// **What**: scale an array so its entries sum to one.
///
/// **How**: add everything up, then divide each entry by that total.
///
/// **Why there is no guard for a zero sum**: because the reference has none. An all-zero input
/// divides zero by zero and yields `NaN` ("not a number", the floating-point result of an
/// undefined operation) rather than zeros. A site with no informative reads gets an allele fraction
/// of `NaN`, and a consumer reading it as a number has to survive that.
///
/// **Why the return type is `Option`**: the reference's negative-sum check throws
/// `IllegalArgumentException`. `Option<T>` is Rust's "a value or nothing", checked by the compiler,
/// which is how this codebase models a reference exception that a caller might legitimately hit.
/// `None` here means "the reference would have thrown".
pub fn normalize_sum_to_one(array: &[f64]) -> Option<Vec<f64>> {
    if array.is_empty() {
        // `Some(x)` wraps a real value; `Vec::new()` is an empty list. The reference returns the
        // input array unchanged for an empty input, which is the same thing observed from outside.
        return Some(Vec::new());
    }
    // `.sum()` needs to know what type it is accumulating into, which the `: f64` annotation
    // supplies. It adds left to right, which is the order the reference's own loop uses.
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
/// **What**: round a number to the nearest integer, halves going away from zero.
///
/// **How**: add or subtract a half depending on the sign, then throw away the fractional part.
///
/// **Why that is not the same as rounding**: the addition happens first and is itself rounded. On a
/// value whose sum with a half is not exactly representable, the sum rounds **up** to the next
/// integer and the truncation then keeps it, so the function rounds twice. `0.49999999999999994` is
/// the largest double below a half; adding a half gives exactly `1.0`, and this function answers
/// one where half-up answers zero. `Math.round` stopped doing this in Java 7; `fastRound` never
/// did.
///
/// **Why the three guards**: Java's `(int)` narrowing clamps to the extremes rather than wrapping
/// around, and answers zero for `NaN`. Rust's `as` conversion has had exactly those semantics since
/// 1.45, so the guards are belt and braces and document the intent.
pub fn fast_round(d: f64) -> i32 {
    let shifted = if d > 0.0 { d + 0.5 } else { d - 0.5 };
    if shifted.is_nan() {
        return 0;
    }
    // `i32::MAX as f64` converts the largest 32-bit integer to a float so the two can be compared.
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
/// **What**: the smaller of two numbers, with Java's exact treatment of the two odd cases.
///
/// **How**: test for `NaN` first, then for the two zeros, then compare.
///
/// **Why it cannot be the built-in**: Rust's `min` returns the **non**-`NaN` argument, so
/// `min(1.0, NaN)` is `1.0`; Java's propagates the `NaN`. Rust's also leaves the choice between
/// `0.0` and `-0.0` explicitly unspecified; Java's prefers the negative zero. Both differences are
/// reachable from `Math.max(0., Math.min(1., pval))` in `ExcessHet`, where a `NaN` p-value would
/// silently become one under the built-in and stay `NaN` under this.
pub fn java_min(a: f64, b: f64) -> f64 {
    // `a != a` is Java's own idiom for "a is NaN", since NaN is the only value unequal to itself.
    // Rust spells it `is_nan()`, which is the same test.
    if a.is_nan() {
        return a;
    }
    // Both are zero **and** the second is the negative one: Java returns the negative zero. The two
    // zeros compare equal, so `a <= b` below could not distinguish them and this case must come
    // first.
    if a == 0.0 && b == 0.0 && b.is_sign_negative() {
        return b;
    }
    if a <= b {
        a
    } else {
        b
    }
}

/// `Math.max(double, double)`. See [`java_min`]: same three cases, mirrored.
pub fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    // Mirrored from `java_min`: here it is the **first** argument being the negative zero that
    // makes the reference return the other one.
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
/// **What**: turn a Phred-scaled quality score into the probability that the call it describes is
/// correct. Phred 30 means one error in a thousand, so the answer is 0.999.
///
/// **How**: `10^(-q/10)` is the error probability by definition of the scale; one minus it is the
/// probability of being right.
///
/// **Why the argument is `f64` and not a byte**: `QualityUtils` has two overloads. The byte one
/// reads a precomputed cache; the double one calls `Math.pow` for real. An `int` genotype quality
/// widens to `double` at the call site in `GenotypeUtils`, so Java picks the double overload, the
/// cache is bypassed, and the `pow` is the deferred intrinsic. Binding to the cached one here would
/// be a different function.
pub fn qual_to_prob(qual: f64) -> f64 {
    1.0 - pow10(qual / -10.0)
}

// `#[cfg(test)]` compiles this module only when running tests, so none of it ships.
#[cfg(test)]
mod tests {
    // Brings everything from the parent module into scope, so the tests can call the functions
    // above by their bare names.
    use super::*;

    #[test]
    fn the_first_maximum_wins_a_tie() {
        // Two equal maxima at positions 0 and 1. The strictly-greater comparison keeps the first.
        assert_eq!(max_element_index(&[1.0, 1.0, 0.0], 0, 3), 0);
    }

    #[test]
    fn fast_round_is_half_away_from_zero() {
        assert_eq!(fast_round(2.5), 3);
        assert_eq!(fast_round(-2.5), -3);
        assert_eq!(fast_round(0.4), 0);
        // The truncating cast is what makes this not `Math.round`: the double below a half plus a
        // half is exactly one, so it rounds twice and answers one where half-up answers zero. The
        // underscores in the literal are digit separators and have no effect on the value.
        assert_eq!(fast_round(0.499_999_999_999_999_94), 1);
    }

    #[test]
    fn an_all_zero_vector_normalises_to_nan() {
        let out = normalize_sum_to_one(&[0.0, 0.0]).expect("a non-negative sum");
        // `.all(...)` is true when the closure holds for every element. Zero divided by zero is
        // NaN, so both entries are.
        assert!(out.iter().all(|value| value.is_nan()));
    }
}
