//! `MannWhitneyU`, ported from `org.broadinstitute.hellbender.utils.MannWhitneyU` (GATK 4.6.2.0).
//!
//! The rank-sum test behind `BaseQRankSum`, `MQRankSum`, `ReadPosRankSum` and
//! `ClippingRankSum`. Four of its decisions are not what "a Mann-Whitney U test" suggests.
//!
//! # The ranks are **floats**
//!
//! ```java
//! private static final class Rank { final double value; float rank; final int series; }
//! ...
//! float r1 = 0, r2 = 0;
//! ```
//!
//! The values are doubles and the ranks are single precision, and so are the sums of ranks that
//! the statistic is computed from. With more than about 2^24 rank units the sum stops being exact,
//! and long before that the *averaged* rank of a tie band is a float division. A port that ranked
//! in double precision would agree on small inputs and drift on realistic ones.
//!
//! # Two different tests, chosen by size
//!
//! ```java
//! if (n1 >= MINIMUM_NORMAL_N || n2 >= MINIMUM_NORMAL_N) { ... normal approximation ... }
//! else { p = permutationTest(...); z = NORMAL.inverseCumulativeProbability(p); }
//! ```
//!
//! With **either** series at 10 or more it is the normal approximation; otherwise it enumerates
//! every permutation of the group labels. So the reported Z comes from a cumulative probability in
//! one case and from an inverse cumulative probability in the other, and those two are not
//! inverses of each other in commons-math3 (see `jmath::normal`). The two regimes therefore do not
//! meet smoothly at the boundary, and that discontinuity is in the reference.
//!
//! # The continuity correction disappears when everything is tied
//!
//! ```java
//! if (nties == 0) { correction = 0; }
//! ```
//!
//! `transformTies` deliberately reports **zero** ties when every single value is tied, because the
//! sigma formula breaks down there; the correction is then dropped so the answer comes out at
//! exactly p = 0.5. Two different mechanisms conspire to produce one intended number.
//!
//! # The permutation p-value counts half of its own bin
//!
//! ```java
//! double sumOfAllSmallerBins = histo.get(round(2 * testStatU)).getValue() / 2.0;
//! ```
//!
//! Half the observed bin plus everything more extreme, rather than the cumulative distribution,
//! *"which gives a p-value of 1 in the most extreme case and doesn't result in a usable
//! z-score"*. The histogram is of **twice** U, because U is integer or half-integer.

use jmath::normal::NormalDistribution;

/// `MannWhitneyU.TestType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TestType {
    FirstDominates,
    SecondDominates,
    TwoSided,
}

/// `MannWhitneyU.MINIMUM_NORMAL_N`: the length at which the normal approximation takes over.
pub const MINIMUM_NORMAL_N: usize = 10;

/// `MannWhitneyU.Result`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MannWhitneyResult {
    pub u: f64,
    pub z: f64,
    pub p: f64,
    pub median_shift: f64,
}

/// One ranked observation. The rank is a `float`, as in the reference.
#[derive(Debug, Clone, Copy)]
struct Rank {
    value: f64,
    rank: f32,
    series: u8,
}

/// `calculateRank`: merge two sorted series, then average the rank of each tie band.
///
/// The merge is stable towards series 1 (`series1[i] <= series2[j]` keeps the first), which
/// decides nothing about the statistic but does decide the `series` labels inside a tie band, and
/// therefore which side each averaged rank is added to.
#[allow(clippy::if_same_then_else)]
fn calculate_rank(series1: &[f64], series2: &[f64]) -> (Vec<Rank>, Vec<usize>) {
    let mut a = series1.to_vec();
    let mut b = series2.to_vec();
    // `Arrays.sort(double[])`: -0.0 before 0.0, NaN last.
    a.sort_by(|x, y| x.total_cmp(y));
    b.sort_by(|x, y| x.total_cmp(y));

    let mut ranks: Vec<Rank> = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0usize, 0usize);
    while ranks.len() < a.len() + b.len() {
        // The rank handed to the constructor is the *incremented* counter, so ranks start at 1.
        let r = ranks.len() as f32 + 1.0;
        if i >= a.len() {
            ranks.push(Rank {
                value: b[j],
                rank: r,
                series: 2,
            });
            j += 1;
        // The two middle branches are the same push with the same label, and they are kept
        // apart because the reference keeps them apart: one is "series 2 is exhausted" and the
        // other is "series 1 compares less or equal", and only the second decides a tie.
        } else if j >= b.len() {
            ranks.push(Rank {
                value: a[i],
                rank: r,
                series: 1,
            });
            i += 1;
        } else if a[i] <= b[j] {
            ranks.push(Rank {
                value: a[i],
                rank: r,
                series: 1,
            });
            i += 1;
        } else {
            ranks.push(Rank {
                value: b[j],
                rank: r,
                series: 2,
            });
            j += 1;
        }
    }

    let mut num_of_ties: Vec<usize> = Vec::new();
    let mut index = 0usize;
    while index < ranks.len() {
        let mut rank = ranks[index].rank;
        let mut count = 1usize;
        let mut j = index + 1;
        while j < ranks.len() && ranks[j].value == ranks[index].value {
            rank += ranks[j].rank;
            count += 1;
            j += 1;
        }
        if count > 1 {
            // A float division, and the sum above was a float sum.
            rank /= count as f32;
            for slot in ranks.iter_mut().skip(index).take(count) {
                slot.rank = rank;
            }
            num_of_ties.push(count);
        }
        index += count;
    }

    (ranks, num_of_ties)
}

/// `transformTies`: the tie term of sigma, and zero when *everything* is tied.
fn transform_ties(num_of_ranks: usize, num_of_ties: &[usize]) -> f64 {
    let mut total = 0.0;
    for &count in num_of_ties {
        if count != num_of_ranks {
            // `Math.pow(count, 3) - count`. The exponent is an integer and the result is exactly
            // representable for every count a read pileup can produce, so the deferred `Math.pow`
            // of decision 0007 is not reached in a way that could differ.
            total += (count as f64).powi(3) - count as f64;
        }
    }
    total
}

/// `calculateU1andU2`.
fn u1_and_u2(series1: &[f64], series2: &[f64]) -> (f64, f64, f64) {
    let (ranks, num_of_ties) = calculate_rank(series1, series2);
    let ties_for_sigma = transform_ties(ranks.len(), &num_of_ties);

    // Single precision, as in the reference.
    let (mut r1, mut r2) = (0.0f32, 0.0f32);
    for rank in &ranks {
        if rank.series == 1 {
            r1 += rank.rank;
        } else {
            r2 += rank.rank;
        }
    }

    let n1 = series1.len() as f64;
    let n2 = series2.len() as f64;
    let u1 = r1 as f64 - ((n1 * (n1 + 1.0)) / 2.0);
    let u2 = r2 as f64 - ((n2 * (n2 + 1.0)) / 2.0);
    (u1, u2, ties_for_sigma)
}

/// `calculateZ`, including the continuity correction and the rule that drops it.
pub fn calculate_z(u: f64, n1: usize, n2: usize, nties: f64, which_side: TestType) -> f64 {
    let n1f = n1 as f64;
    let n2f = n2 as f64;
    let m = (n1 * n2) as f64 / 2.0;

    let mut correction = match which_side {
        TestType::TwoSided => {
            if (u - m) >= 0.0 {
                0.5
            } else {
                -0.5
            }
        }
        TestType::FirstDominates => -0.5,
        TestType::SecondDominates => 0.5,
    };
    if nties == 0.0 {
        correction = 0.0;
    }

    let sigma = ((n1 * n2) as f64 / 12.0
        * ((n1f + n2f + 1.0) - nties / ((n1f + n2f) * (n1f + n2f - 1.0))))
        .sqrt();
    (u - m - correction) / sigma
}

/// `median(double[])`, which assumes the array is already sorted and does **not** sort it.
///
/// This is not `MathUtils.median`: no commons-math3, no interpolation beyond the two-element
/// average, and an odd count takes the upper of the two middles.
pub fn median(sorted: &[f64]) -> f64 {
    let len = sorted.len();
    let mid = len / 2;
    #[allow(clippy::manual_is_multiple_of)]
    if len % 2 == 0 {
        (sorted[mid] + sorted[mid - 1]) / 2.0
    } else {
        sorted[mid]
    }
}

/// Every distinct permutation of a multiset of labels, in the order the reference's
/// next-permutation walk produces them.
///
/// The reference collects them into a `HashSet<List<Integer>>` first, which changes nothing here:
/// the walk already yields each permutation once, and the histogram they feed does not care about
/// order.
fn permutations(labels: &[u8]) -> Vec<Vec<u8>> {
    let mut current = labels.to_vec();
    let mut out = vec![current.clone()];
    loop {
        let mut k = None;
        for i in (0..current.len().saturating_sub(1)).rev() {
            if current[i] < current[i + 1] {
                k = Some(i);
                break;
            }
        }
        let Some(k) = k else { break };
        let mut l = k + 1;
        for i in (k + 1..current.len()).rev() {
            if current[k] < current[i] {
                l = i;
                break;
            }
        }
        current.swap(k, l);
        current[k + 1..].reverse();
        out.push(current.clone());
    }
    out
}

/// `permutationTest`: the exact test, used when both series are shorter than ten.
fn permutation_test(series1: &[f64], series2: &[f64], test_stat_u: f64) -> f64 {
    let n1 = series1.len();
    let n2 = series2.len();
    let (ranks, _) = calculate_rank(series1, series2);

    // The first permutation is `n1` zeroes then `n2` ones, which is the smallest arrangement, so
    // the walk enumerates all of them.
    let mut labels = vec![0u8; n1 + n2];
    for label in labels.iter_mut().skip(n1) {
        *label = 1;
    }

    // A histogram of **twice** U, because U is integer or half-integer.
    let mut histogram: Vec<(i64, f64)> = Vec::new();
    for permutation in permutations(&labels) {
        let mut sum = 0.0f64;
        for (index, group) in permutation.iter().enumerate() {
            if *group == 0 {
                // `MathUtils.sum` over the ranks assigned to group one, in permutation order.
                sum += ranks[index].rank as f64;
            }
        }
        let new_u = sum - ((n1 * (n1 + 1)) as f64 / 2.0);
        let key = jmath::fast_math::round(2.0 * new_u);
        match histogram.iter_mut().find(|(id, _)| *id == key) {
            Some((_, value)) => *value += 1.0,
            None => histogram.push((key, 1.0)),
        }
    }

    let observed = jmath::fast_math::round(2.0 * test_stat_u);
    // Half of the observed bin, plus everything more extreme. `histo.get(key)` would raise a
    // NullPointerException for a key the histogram does not hold, and the observed statistic is
    // always one of the permutations, so the bin is always there.
    let mut sum_of_all_smaller_bins = histogram
        .iter()
        .find(|(id, _)| *id == observed)
        .map(|(_, value)| value / 2.0)
        .unwrap_or(f64::NAN);
    for (id, value) in &histogram {
        if *id < observed {
            sum_of_all_smaller_bins += value;
        }
    }
    let total: f64 = histogram.iter().map(|(_, value)| value).sum();
    sum_of_all_smaller_bins / total
}

/// `MannWhitneyU.test(series1, series2, whichSide)`.
///
/// The arrays are sorted in place by the reference; here they are copied, which is invisible to
/// the result and visible to a caller that reused the array afterwards.
pub fn test(series1: &[f64], series2: &[f64], which_side: TestType) -> MannWhitneyResult {
    let n1 = series1.len();
    let n2 = series2.len();
    if n1 == 0 || n2 == 0 {
        // `new Result(Float.NaN, ...)`: a float NaN widened to double, which is the same bits as
        // a double NaN here.
        return MannWhitneyResult {
            u: f64::NAN,
            z: f64::NAN,
            p: f64::NAN,
            median_shift: f64::NAN,
        };
    }

    let (u1, u2, nties) = u1_and_u2(series1, series2);
    let u = match which_side {
        TestType::TwoSided => u1.min(u2),
        TestType::FirstDominates => u1,
        TestType::SecondDominates => u2,
    };

    let normal = NormalDistribution::default();
    let (z, p);
    if n1 >= MINIMUM_NORMAL_N || n2 >= MINIMUM_NORMAL_N {
        z = calculate_z(u, n1, n2, nties, which_side);
        let two_sided = 2.0 * normal.cumulative_probability(0.0 + z * 1.0);
        p = if which_side == TestType::TwoSided {
            two_sided
        } else {
            two_sided / 2.0
        };
    } else {
        // Only the one-sided exact test exists; the reference logs a warning for the others and
        // runs it anyway.
        p = permutation_test(series1, series2, u);
        z = normal.inverse_cumulative_probability(p).unwrap_or(f64::NAN);
    }

    let mut sorted1 = series1.to_vec();
    let mut sorted2 = series2.to_vec();
    sorted1.sort_by(|a, b| a.total_cmp(b));
    sorted2.sort_by(|a, b| a.total_cmp(b));
    MannWhitneyResult {
        u,
        z,
        p,
        median_shift: (median(&sorted1) - median(&sorted2)).abs(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_empty_series_is_all_nan() {
        let result = test(&[], &[1.0], TestType::FirstDominates);
        assert!(result.u.is_nan() && result.z.is_nan() && result.p.is_nan());
    }

    #[test]
    fn everything_tied_drops_the_continuity_correction() {
        let tied = vec![5.0; 12];
        let result = test(&tied, &tied, TestType::FirstDominates);
        // Half the mass on each side, exactly, which is what the two mechanisms are for.
        assert!((result.p - 0.5).abs() < 1e-12, "p = {}", result.p);
    }

    #[test]
    fn the_permutation_walk_produces_every_distinct_arrangement() {
        let labels = [0u8, 0, 1, 1];
        assert_eq!(permutations(&labels).len(), 6);
    }
}
