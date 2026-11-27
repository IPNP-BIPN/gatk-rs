//! `MannWhitneyU`, ported from `org.broadinstitute.hellbender.utils.MannWhitneyU` (GATK 4.7.0.0;
//! `MannWhitneyU.java` is unchanged since 4.6.2.0).
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
//! and the *averaged* rank of a tie band is a float division. Ranks are whole or half numbers, so a
//! float sum of them stays exact up to about 4,095 observations; a port that ranked in double
//! precision would agree below that and could drift above it.
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
//! # The continuity correction disappears when everything is tied, and when nothing is
//!
//! ```java
//! if (nties == 0) { correction = 0; }
//! ```
//!
//! `transformTies` deliberately reports **zero** ties when every single value is tied, because the
//! sigma formula breaks down there; the correction is then dropped so the answer comes out at
//! exactly p = 0.5. Two different mechanisms conspire to produce one intended number. The same
//! zero is also what `transformTies` returns for data with no ties at all, so tie-free input (the
//! golden's `large-ramp`) gets no continuity correction either.
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
//!
//! # Reading this file without reading Rust
//!
//! - `&[f64]` is a borrowed, read-only view of a `double[]`; the caller keeps the array. The
//!   reference sorts the caller's arrays in place (inside `calculateRank`); a read-only view
//!   cannot be sorted, so each function that needs sorted data copies it first (`to_vec()`).
//! - `x as f32` and `x as f64` are Java's `(float) x` and `(double) x`.
//! - A tuple such as `(Vec<Rank>, Vec<usize>)` is an anonymous record. Tuples stand in for the
//!   reference's `RankedData` and `TestStatistic` holder classes, which carry nothing else.
//! - `Option<T>` is "a value or nothing" (`Some(x)` or `None`), a nullable reference that the
//!   compiler makes the code check. `let Some(k) = k else { break };` unwraps it or leaves the
//!   loop.
//! - A `match` over `TestType` must name every variant (the compiler rejects a missing one), and
//!   like most Rust constructs it is an expression: its value is the value of the chosen arm.
//! - `usize` is an unsigned index type Java lacks. Java's `int` arithmetic on lengths wraps past
//!   2^31 - 1; `usize` arithmetic does not (see `calculate_z` for where that matters).

use jmath::normal::NormalDistribution;

/// `MannWhitneyU.TestType`.
///
/// Which series the one-sided null hypothesis says dominates. Every rank-sum annotation passes
/// `FirstDominates` with the alt values first; see `test`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TestType {
    FirstDominates,
    SecondDominates,
    TwoSided,
}

/// `MannWhitneyU.MINIMUM_NORMAL_N`: the length at which the normal approximation takes over.
pub const MINIMUM_NORMAL_N: usize = 10;

/// `MannWhitneyU.Result`.
///
/// The Java getters are public fields here. `median_shift` is
/// `Math.abs(median(series1) - median(series2))` (MannWhitneyU.java, `test`), so it does not
/// depend on which series comes first; for the one-sided tests `u`, `z` and `p` do.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MannWhitneyResult {
    pub u: f64,
    pub z: f64,
    pub p: f64,
    pub median_shift: f64,
}

/// One ranked observation. The rank is a `float`, as in the reference.
///
/// ```java
/// final double value;
/// float rank;
/// final int series;
/// ```
///
/// `series` is 1 or 2. The reference's `compareTo` is never called (the `Rank[]` is built in
/// order and never sorted), so it is not ported.
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
///
/// Each rank is the merge position plus one. In `ranks[r++] = new Rank(series2[j++], r, 2);` Java
/// evaluates the index `r++` before the argument `r`, so the argument is the incremented counter.
/// A tie band's ranks are consecutive integers, summed in `float` (exact below 2^24) and divided
/// by the band size; the true average is an integer or half-integer, so the division is exact.
/// `#[allow(clippy::if_same_then_else)]` silences the lint that flags the two equal branches.
#[allow(clippy::if_same_then_else)]
fn calculate_rank(series1: &[f64], series2: &[f64]) -> (Vec<Rank>, Vec<usize>) {
    let mut a = series1.to_vec();
    let mut b = series2.to_vec();
    // `Arrays.sort(double[])`: -0.0 before 0.0, NaN last.
    // `|x, y| ...` is a closure, Java's lambda. `total_cmp` agrees with `Arrays.sort` on every
    // input except a NaN with its sign bit set, which it orders first where Java puts all NaN last.
    a.sort_by(|x, y| x.total_cmp(y));
    b.sort_by(|x, y| x.total_cmp(y));

    let mut ranks: Vec<Rank> = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0usize, 0usize);
    // Java's write index `r` is `ranks.len()`: each `push` appends one element.
    while ranks.len() < a.len() + b.len() {
        // The rank handed to the constructor is the *incremented* counter, so ranks start at 1.
        // Java widens the `int` once; this widens, then adds 1.0. Same value below 2^24 positions.
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
        // `<=` and `==` on f64 are Java's: -0.0 equals 0.0, and every comparison with NaN is
        // false, so a NaN in series 1 yields to series 2 here and never joins a tie band below.
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
            // `iter_mut().skip(index).take(count)` yields writable references to
            // `ranks[index..index + count]`: Java's `for (int j = i; j < i + count; ++j)`.
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
///
/// ```java
/// if (count != numOfRanks) {
///     transformedTies.add((Math.pow(count, 3)) - count);
/// }
/// ```
///
/// `count^3 - count` is zero only for a count of 0 or 1, and every recorded band has at least
/// two, so the result is zero in exactly two cases: no tie band at all, or one band that is
/// everything. `calculate_z` cannot tell them apart. The reference collects the terms into a
/// list and sums it in a second loop, in the same order; one loop makes the same additions.
fn transform_ties(num_of_ranks: usize, num_of_ties: &[usize]) -> f64 {
    let mut total = 0.0;
    // `for &count in` copies each element out of the borrowed slice.
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
///
/// The tuple is `(U1, U2, transformTies(..))`, the reference's `TestStatistic`.
///
/// ```java
/// float r1 = 0, r2 = 0;
/// ```
///
/// R1 and R2 are accumulated in `float` and only then widened to `double`. Every rank is an
/// integer or half-integer, which a `float` holds exactly below 2^23, so the sums are exact while
/// N(N+1)/2 < 2^23, that is for up to 4,095 observations in both series together. Past that a sum
/// can round, and a port that summed in double could then report a different U and Z. No golden row
/// reaches that size: the largest fixture in `tests/data/mann_whitney.txt.gz`,
/// `large-interleaved`, has 1,200.
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
///
/// ```java
/// if (nties == 0) {
///     correction = 0;
/// }
/// ```
///
/// `nties` is `transform_ties`, which is zero in two cases, not one: every value tied (the case
/// the reference's own comment names, where dropping the correction gives p = 0.5) and **no ties
/// at all**. So the continuity correction is applied only when the data has at least one tie band
/// and is not entirely tied; tie-free data (`large-ramp` in `tests/data/mann_whitney.txt.gz`)
/// gets none.
///
/// The reference does the length arithmetic in `int`:
///
/// ```java
/// double m = (n1 * n2) / 2d;
/// ```
///
/// and likewise `n1 * n2 / 12d` and `(n1 + n2) * (n1 + n2 - 1)` inside sigma, all of which wrap
/// past 2^31 - 1. Here they are `usize` and `f64` and do not wrap, so this port departs from the
/// reference once `n1 * n2 >= 2^31`, or, when `nties` is not zero, once `n1 + n2 >= 46,342`. The
/// golden's `z` rows stop at `n1 + n2 = 103`.
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

    // `sqrt` is correctly rounded under IEEE 754 in Rust and in Java, so it is not a libm gap.
    let sigma = ((n1 * n2) as f64 / 12.0
        * ((n1f + n2f + 1.0) - nties / ((n1f + n2f) * (n1f + n2f - 1.0))))
        .sqrt();
    (u - m - correction) / sigma
}

/// `median(double[])`, which assumes the array is already sorted and does **not** sort it.
///
/// This is not `MathUtils.median`: no commons-math3, no interpolation beyond the two-element
/// average: an even count averages the two middle values, an odd count takes the middle one.
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
///
/// The walk is `calculatePermutations`, the lexicographic next-permutation algorithm: find the
/// last `k` with `current[k] < current[k + 1]`, the last `l > k` with `current[k] < current[l]`,
/// swap them, reverse the tail. Started from the sorted labels it reaches every distinct
/// arrangement once.
/// The reference's static `PERMUTATIONS` cache is dropped; it only avoids recomputing the set.
fn permutations(labels: &[u8]) -> Vec<Vec<u8>> {
    let mut current = labels.to_vec();
    let mut out = vec![current.clone()];
    loop {
        let mut k = None;
        // `i` from `len - 2` down to 0, Java's `for (int i = temp.length - 2; i >= 0; i--)`.
        // `saturating_sub` stops at 0 instead of wrapping the unsigned `usize`, so a list of
        // fewer than two labels runs no iterations.
        for i in (0..current.len().saturating_sub(1)).rev() {
            if current[i] < current[i + 1] {
                k = Some(i);
                break;
            }
        }
        let Some(k) = k else { break };
        // Java starts `l` at -1. Since `current[k] < current[k + 1]`, the search below always
        // finds an `l`, so the starting value is never used.
        let mut l = k + 1;
        for i in (k + 1..current.len()).rev() {
            if current[k] < current[i] {
                l = i;
                break;
            }
        }
        current.swap(k, l);
        // `current[k + 1..]` is the tail from `k + 1`, reversed in place: Java's swap loop.
        current[k + 1..].reverse();
        out.push(current.clone());
    }
    out
}

/// `permutationTest`: the exact test, used when both series are shorter than ten.
///
/// The reference's `Histogram<Long>` (htsjdk) is a `TreeMap` from twice U to a count, read in key
/// order; here it is a list of `(key, count)` pairs in first-seen order. Every count and every
/// partial sum is a multiple of one half far below 2^53 (at most C(18, 9) = 48,620
/// permutations), so the order of the additions cannot change the result.
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
        // `find` returns an `Option` holding a writable reference to the matching bin, if any:
        // increment it, or add the bin. Java's `histo.increment(FastMath.round(2 * newU))`.
        match histogram.iter_mut().find(|(id, _)| *id == key) {
            Some((_, value)) => *value += 1.0,
            None => histogram.push((key, 1.0)),
        }
    }

    let observed = jmath::fast_math::round(2.0 * test_stat_u);
    // Half of the observed bin, plus everything more extreme. `histo.get(key)` would raise a
    // NullPointerException for a key the histogram does not hold. For `FirstDominates`, the only
    // side any caller asks for, the observed statistic is one of the permutations, so the bin is
    // there; for another side that argument does not hold, and where Java would raise, this
    // returns NaN.
    let mut sum_of_all_smaller_bins = histogram
        .iter()
        .find(|(id, _)| *id == observed)
        .map(|(_, value)| value / 2.0)
        .unwrap_or(f64::NAN);
    // Strictly smaller keys only, the FIRST_DOMINATES direction ("in the FIRST_DOMINATES case
    // the smaller bins", says the reference), whatever side the caller asked for.
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
///
/// Every caller passes the **alt** values first and `FirstDominates`. `RankSumTest.java` (one line
/// in the source, wrapped here):
///
/// ```java
/// mannWhitneyU.test(Doubles.toArray(altQuals), Doubles.toArray(refQuals),
/// MannWhitneyU.TestType.FIRST_DOMINATES);
/// ```
///
/// `AS_RankSumTest.java` does the same, as do `rank_sum.rs` and `allele_specific_rank_sum.rs` in
/// `gatk-annotation`. A negative Z therefore means the alt values rank low. Swapping the two
/// arguments is, in the normal regime, not an exact negation of Z: the correction stays -0.5
/// for `FirstDominates` whichever series is first, so in exact arithmetic the swapped Z is
/// `-Z + 1/sigma` whenever the correction applies, and exactly `-Z` only when `nties` is zero.
/// In the permutation regime the `tiny-separated` and `tiny-reversed` rows of
/// `tests/data/mann_whitney.txt.gz` (the same data, swapped) give Z = -1.95996... and
/// +1.95996... Asking `SecondDominates` of the same order is an exact negation in the normal
/// regime (its `qualities` rows: U 104 and 40, Z +2.236 and -2.236).
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

    // The reference's `calculateTwoSidedU` and `calculateOneSidedU` are the `match` below.
    let (u1, u2, nties) = u1_and_u2(series1, series2);
    let u = match which_side {
        TestType::TwoSided => u1.min(u2),
        TestType::FirstDominates => u1,
        TestType::SecondDominates => u2,
    };

    let normal = NormalDistribution::default();
    // Declared without a value and assigned once in each branch, as Java's `double z;` and
    // `double p;`.
    let (z, p);
    if n1 >= MINIMUM_NORMAL_N || n2 >= MINIMUM_NORMAL_N {
        z = calculate_z(u, n1, n2, nties, which_side);
        // `0.0 + z * 1.0` is the reference's `NORMAL_MEAN + z * NORMAL_SD`, kept as written.
        let two_sided = 2.0 * normal.cumulative_probability(0.0 + z * 1.0);
        p = if which_side == TestType::TwoSided {
            two_sided
        } else {
            two_sided / 2.0
        };
    } else {
        // Only the one-sided exact test exists; the reference logs a warning for the others and
        // runs it anyway. The warning's guard is `if (whichSide != TestType.FIRST_DOMINATES) {`,
        // so it covers `SecondDominates` too; the log line is not ported.
        p = permutation_test(series1, series2, u);
        // `None` would be the reference's `OutOfRangeException` for p outside [0, 1]. With
        // `FirstDominates` the permutation p is (half the observed bin + the smaller bins) / all
        // bins: above 0, and below 1 because half the observed bin is left out, so the fallback
        // is unreachable. For the other sides `u` need not be a U1 key: `([1], [1, 2])` with
        // `SecondDominates` looks up key 3 in keys {1, 4}, where Java throws a
        // NullPointerException and `permutation_test` returns NaN.
        z = normal.inverse_cumulative_probability(p).unwrap_or(f64::NAN);
    }

    // The reference's `median(series1)` gets arrays that `calculateRank` already sorted in place,
    // and `median` does not sort, so the copies here are sorted to match.
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
