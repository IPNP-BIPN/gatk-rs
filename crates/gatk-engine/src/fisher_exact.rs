//! `FisherExactTest` and the two `MathUtils` sums under it, ported from
//! `org.broadinstitute.hellbender.utils` (GATK 4.6.2.0).
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
//!     return Math.pow(10.0, log10SumLog10(log10values));
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

use jmath::saddle_point::hypergeometric_log_probability;

/// `FisherExactTest.REL_ERR`.
const REL_ERR: f64 = 1.0 - 10e-7;

/// `MathUtils.LOG10_E`, which is `Math.log10(Math.E)`.
fn log10_e() -> f64 {
    jmath::math::log10(std::f64::consts::E)
}

/// `MathUtils.logToLog10`.
pub fn log_to_log10(ln: f64) -> f64 {
    ln * log10_e()
}

/// `MathUtils.log10SumLog10(values)`: the max-shifted sum.
///
/// The shift means the largest term contributes exactly 1 and never goes through `pow`, so the
/// answer depends on which element is the maximum and not only on the multiset of values.
pub fn log10_sum_log10(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }
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
    let mut sum = 1.0;
    for (index, value) in values.iter().enumerate() {
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
    max_value
        + if sum != 1.0 {
            jmath::math::log10(sum)
        } else {
            0.0
        }
}

/// `MathUtils.sumLog10`.
pub fn sum_log10(values: &[f64]) -> f64 {
    pow10(log10_sum_log10(values))
}

/// `Math.pow(10.0, x)`.
///
/// Decision 0007 deferred `Math.pow`: its HotSpot intrinsic is GPL2 and cannot be transcribed
/// here, and it is not correctly rounded, so no independent implementation is guaranteed to agree.
/// This calls the platform's, and the conformance suite is what says whether that is the same
/// double on the inputs the strand-bias annotations produce.
fn pow10(x: f64) -> f64 {
    10.0f64.powf(x)
}

/// `FisherExactTest.twoSidedPValue(normalizedTable)`.
///
/// The table is `[[a, b], [c, d]]`, and the reference refuses anything that is not 2x2.
pub fn two_sided_p_value(table: [[i32; 2]; 2]) -> f64 {
    let m = table[0][0] + table[0][1];
    let n = table[1][0] + table[1][1];
    let k = table[0][0] + table[1][0];
    let lo = (k - n).max(0);
    let hi = k.min(m);

    // `IndexRange(lo, hi + 1)`: a support of one point is exactly 1.0, with nothing computed.
    if (hi + 1 - lo) <= 1 {
        return 1.0;
    }

    let log_densities: Vec<f64> = (lo..=hi)
        .map(|x| hypergeometric_log_probability(m + n, m, k, x))
        .collect();

    // The observed point's density, slackened by a relative error. Multiplying a negative log
    // density by a factor below 1 *raises* it, which is what keeps ties inside the sum.
    let threshold = log_densities[(table[0][0] - lo) as usize] * REL_ERR;
    let log10_densities: Vec<f64> = log_densities
        .iter()
        .filter(|density| **density <= threshold)
        .map(|density| log_to_log10(*density))
        .collect();

    let p_value = sum_log10(&log10_densities);
    // "min is necessary as numerical precision can result in pValue being slightly greater than 1"
    p_value.min(1.0)
}

#[cfg(test)]
mod tests {
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
