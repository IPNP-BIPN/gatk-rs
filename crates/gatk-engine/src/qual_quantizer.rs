//! `QualQuantizer` and `QuantizationInfo`, ported from
//! `org.broadinstitute.hellbender.utils.recalibration` (GATK 4.6.2.0).
//!
//! The map from every original quality score to the quantized one that replaces it. It is the last
//! table a recalibration report carries, `BaseRecalibrator` computes it and `ApplyBQSR` reads it
//! back, so it is the interface between the two tools.
//!
//! The algorithm is a greedy agglomeration. Every quality score starts as its own interval; the
//! adjacent pair whose merge costs least is merged; repeat until only `n_levels` intervals remain.
//! The cost is a penalty in log space, weighted by observations.
//!
//! # Every leaf carries a fixed quality, and that changes what an error rate means
//!
//! `QualInterval.getErrorRate` reads:
//!
//! ```java
//! if ( hasFixedQual() ) return QualityUtils.qualToErrorProb((byte) fixedQual);
//! else if ( nObservations == 0 ) return 0.0;
//! else return (nErrors+1) / (1.0 * (nObservations+1));
//! ```
//!
//! and `quantize()` builds every leaf through the constructor that **sets** `fixedQual` to the
//! quality score itself. So the third branch, the one that looks like the definition, is only ever
//! reached by a merged interval. A leaf's error rate is the theoretical one its Phred score
//! declares, not the one its counts imply, and its quantized quality is that score unchanged.
//!
//! That is why a quantized quality of **zero** is reachable even though
//! [`error_prob_to_qual`] clamps to one: a leaf that was never merged never goes through it.
//!
//! # And its error count saturates
//!
//! ```java
//! final double nErrors = nObs * errorRate;
//! new QualInterval(qStart, qStart, nObs, (int) Math.floor(nErrors), 0, (byte)qStart);
//! ```
//!
//! A `long` count of observations, a `double` product, and an `int` in the middle. Java's narrowing
//! cast **clamps** rather than wrapping, so three billion observations at quality zero give
//! `Integer.MAX_VALUE` errors and not three billion. Rust's `as` clamps the same way, which is the
//! one place where the obvious port is also the faithful one.

use crate::math_utils::{pow10, qual_to_error_prob};

/// `QualityUtils.MIN_USABLE_Q_SCORE`, the default below which merges are free.
pub const MIN_USABLE_Q_SCORE: i32 = 6;

/// `QualityUtils.MAX_SAM_QUAL_SCORE`.
pub const MAX_SAM_QUAL_SCORE: i32 = 93;

/// What the quantizer refuses, and the two ends it does not refuse but cannot survive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QuantizerError {
    /// `GATKException`, whose message lists the whole histogram.
    NegativeCounts(String),
    /// `GATKException("nLevels must be >= 0")`.
    NegativeLevels,
    /// `GATKException("minInterestingQual must be >= 0")`.
    NegativeMinInterestingQual,
    /// The reference's `NullPointerException`: with `nLevels` of zero the merge loop runs once more
    /// than there are pairs and dereferences the minimum it never found.
    NoPairToMerge,
    /// The reference's `NoSuchElementException`, with no message: the minimum of an empty histogram.
    EmptyHistogram,
}

impl QuantizerError {
    pub fn message(&self) -> String {
        match self {
            QuantizerError::NegativeCounts(histogram) => {
                format!("Quality score histogram has negative values at: {histogram}")
            }
            QuantizerError::NegativeLevels => "nLevels must be >= 0".to_string(),
            QuantizerError::NegativeMinInterestingQual => {
                "minInterestingQual must be >= 0".to_string()
            }
            QuantizerError::NoPairToMerge => {
                "Cannot read field \"subIntervals\" because \"minMerge\" is null".to_string()
            }
            // The reference's exception carries no message at all.
            QuantizerError::EmptyHistogram => "null".to_string(),
        }
    }
}

/// `QualityUtils.errorProbToQual(errorRate)`.
///
/// Round the Phred score, then clamp to `[1, maxQual]`. Two ends worth naming: an error rate of one
/// is Phred zero and comes back as **1**, and an error rate of zero is Phred infinity, which
/// `Math.round` turns into `Long.MAX_VALUE`, the narrowing cast clamps to `Integer.MAX_VALUE` and
/// the bound brings back to 93.
///
/// `None` where the reference throws: the argument must be a probability, which NaN is not.
pub fn error_prob_to_qual(error_rate: f64) -> Option<u8> {
    // `MathUtils.isValidProbability`: `0.0 <= p <= 1.0`, which is false for NaN.
    if !(0.0..=1.0).contains(&error_rate) {
        return None;
    }
    // `final double d = Math.round(...)` then `(int) d`. The round returns a `long` and is
    // **widened to a double** before the narrowing, so the cast that reaches the bound is
    // double-to-int, which SATURATES. Going straight from the long would be a long-to-int
    // narrowing, which wraps, and an error rate of zero would come back as 1 instead of 93.
    let rounded = jmath::math::round(-10.0 * jmath::math::log10(error_rate)) as f64;
    Some(bound_qual(rounded as i32, MAX_SAM_QUAL_SCORE))
}

/// `QualityUtils.boundQual(qual, maxQual)`: clamped to `[1, maxQual]`, never to zero.
pub fn bound_qual(qual: i32, max_qual: i32) -> u8 {
    (qual.min(max_qual).max(1) & 0xFF) as u8
}

/// `QualQuantizer.QualInterval`: a contiguous run of quality scores, with both ends inclusive.
#[derive(Debug, Clone, PartialEq)]
pub struct QualInterval {
    pub q_start: i32,
    pub q_end: i32,
    pub n_observations: i64,
    pub n_errors: i64,
    pub level: i32,
    /// -1 for a merged interval. See the module note: a leaf always has one, and that decides both
    /// its error rate and its quality.
    pub fixed_qual: i32,
    /// The two intervals this one was merged from, left first. Empty for a leaf.
    pub sub_intervals: Vec<QualInterval>,
}

impl QualInterval {
    /// `getName()`: `qStart-qEnd`, which is what the golden identifies an interval by.
    pub fn name(&self) -> String {
        format!("{}-{}", self.q_start, self.q_end)
    }

    /// `hasFixedQual()`.
    pub fn has_fixed_qual(&self) -> bool {
        self.fixed_qual != -1
    }

    /// `getErrorRate()`. See the module note: the third branch is only for merged intervals.
    pub fn error_rate(&self) -> f64 {
        if self.has_fixed_qual() {
            // The byte overload, which reads a cache the double one bypasses. Both are
            // `10^(-q/10)`, and the cache was filled from the same function.
            qual_to_error_prob(self.fixed_qual as f64)
        } else if self.n_observations == 0 {
            0.0
        } else {
            (self.n_errors + 1) as f64 / (1.0 * (self.n_observations + 1) as f64)
        }
    }

    /// `getQual()`: the fixed quality if there is one, otherwise the Phred of the error rate.
    pub fn qual(&self) -> u8 {
        if self.has_fixed_qual() {
            // The reference casts the int field to a byte with no clamp, which is how a quantized
            // quality of zero exists.
            self.fixed_qual as u8
        } else {
            error_prob_to_qual(self.error_rate()).unwrap_or(1)
        }
    }

    /// `merge(toMerge)`: contiguous only, and the level is one above the higher of the two.
    ///
    /// The reference throws when the two are not adjacent. That is unreachable from `quantize`,
    /// which only ever pairs neighbours in the sorted set, so it is not ported as an error.
    pub fn merge(&self, other: &QualInterval) -> QualInterval {
        let (left, right) = if self.q_start < other.q_start {
            (self, other)
        } else {
            (other, self)
        };
        QualInterval {
            q_start: left.q_start,
            q_end: right.q_end,
            n_observations: left.n_observations + right.n_observations,
            n_errors: left.n_errors + right.n_errors,
            level: left.level.max(right.level) + 1,
            fixed_qual: -1,
            sub_intervals: vec![left.clone(), right.clone()],
        }
    }

    /// `getPenalty()`: the cost of approximating everything under this interval with its own error
    /// rate.
    pub fn penalty(&self, min_interesting_qual: i32) -> f64 {
        self.calc_penalty(self.error_rate(), min_interesting_qual)
    }

    /// `calcPenalty(globalErrorRate)`.
    ///
    /// Three things decide the answer and only one of them is arithmetic:
    ///
    ///  * **a global error rate of zero is a penalty of zero**, tested first, so an empty histogram
    ///    never reaches the leaves and every merge costs the same;
    ///  * **a leaf at or below `minInterestingQual` contributes nothing**, which is what "free to
    ///    merge" means;
    ///  * and the sum is over **leaves**, recursively, each against the same global rate rather than
    ///    against its parent's.
    fn calc_penalty(&self, global_error_rate: f64, min_interesting_qual: i32) -> f64 {
        if global_error_rate == 0.0 {
            return 0.0;
        }
        if self.sub_intervals.is_empty() {
            if self.q_end <= min_interesting_qual {
                0.0
            } else {
                (log10(self.error_rate()) - log10(global_error_rate)).abs()
                    * self.n_observations as f64
            }
        } else {
            self.sub_intervals
                .iter()
                .map(|interval| interval.calc_penalty(global_error_rate, min_interesting_qual))
                .sum()
        }
    }
}

/// `Math.log10`, kept named so the deferred intrinsic has one call site here as elsewhere.
fn log10(x: f64) -> f64 {
    jmath::math::log10(x)
}

/// `QualQuantizer`: the histogram, the levels, and the map it produced.
#[derive(Debug, Clone, PartialEq)]
pub struct QualQuantizer {
    pub n_levels: i32,
    pub min_interesting_qual: i32,
    /// The final forest, in `q_start` order.
    pub quantized_intervals: Vec<QualInterval>,
    /// `getOriginalToQuantizedMap()`: one entry per bin of the histogram.
    pub original_to_quantized_map: Vec<u8>,
}

impl QualQuantizer {
    /// The constructor, whose three checks run in the reference's order.
    pub fn new(
        n_observations_per_qual: &[i64],
        n_levels: i32,
        min_interesting_qual: i32,
    ) -> Result<QualQuantizer, QuantizerError> {
        // `Collections.min` of an empty list is a NoSuchElementException, and it happens before
        // either of the two checks below.
        let minimum = n_observations_per_qual
            .iter()
            .min()
            .ok_or(QuantizerError::EmptyHistogram)?;
        if *minimum < 0 {
            let histogram = n_observations_per_qual
                .iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join(", ");
            return Err(QuantizerError::NegativeCounts(histogram));
        }
        if n_levels < 0 {
            return Err(QuantizerError::NegativeLevels);
        }
        if min_interesting_qual < 0 {
            return Err(QuantizerError::NegativeMinInterestingQual);
        }

        let quantized_intervals =
            quantize(n_observations_per_qual, n_levels, min_interesting_qual)?;
        let original_to_quantized_map =
            intervals_to_map(&quantized_intervals, n_observations_per_qual.len());

        Ok(QualQuantizer {
            n_levels,
            min_interesting_qual,
            quantized_intervals,
            original_to_quantized_map,
        })
    }
}

/// `quantize()`: one interval per quality score, then merge until `n_levels` remain.
fn quantize(
    histogram: &[i64],
    n_levels: i32,
    min_interesting_qual: i32,
) -> Result<Vec<QualInterval>, QuantizerError> {
    let mut intervals: Vec<QualInterval> = histogram
        .iter()
        .enumerate()
        .map(|(q_start, n_obs)| {
            let q_start = q_start as i32;
            let error_rate = qual_to_error_prob(q_start as f64);
            let n_errors = *n_obs as f64 * error_rate;
            QualInterval {
                q_start,
                q_end: q_start,
                n_observations: *n_obs,
                // `(int) Math.floor(nErrors)` widened back to a long. See the module note: the
                // narrowing clamps, so a huge count saturates at Integer.MAX_VALUE.
                n_errors: n_errors.floor() as i32 as i64,
                level: 0,
                fixed_qual: q_start,
                sub_intervals: Vec::new(),
            }
        })
        .collect();

    while intervals.len() as i32 > n_levels {
        merge_lowest_penalty_intervals(&mut intervals, min_interesting_qual)?;
    }
    Ok(intervals)
}

/// `mergeLowestPenaltyIntervals`: merge the adjacent pair whose merge costs least.
///
/// **The first minimum wins.** The comparison is strictly less-than over a walk in `q_start` order,
/// so a tie goes to the leftmost pair, and ties are the rule rather than the exception: an empty
/// histogram gives every pair a penalty of zero.
///
/// **With one interval left there is no pair.** The reference builds two iterators, skips one on the
/// second, and then dereferences a minimum it never found. That is the `nLevels = 0` end.
fn merge_lowest_penalty_intervals(
    intervals: &mut Vec<QualInterval>,
    min_interesting_qual: i32,
) -> Result<(), QuantizerError> {
    let mut min_merge: Option<(usize, QualInterval, f64)> = None;
    for index in 0..intervals.len().saturating_sub(1) {
        let merged = intervals[index].merge(&intervals[index + 1]);
        let penalty = merged.penalty(min_interesting_qual);
        match &min_merge {
            Some((_, _, best)) if penalty >= *best => {}
            _ => min_merge = Some((index, merged, penalty)),
        }
    }
    let (index, merged, _) = min_merge.ok_or(QuantizerError::NoPairToMerge)?;
    // `intervals.removeAll(minMerge.subIntervals)` over a TreeSet ordered by qStart, then `add`,
    // which puts the merged interval where the left one was.
    intervals.remove(index + 1);
    intervals[index] = merged;
    Ok(())
}

/// `intervalsToMap`: every quality score in an interval takes that interval's quality.
///
/// The reference fills the map with `Byte.MIN_VALUE` first and throws if any survives. That cannot
/// happen, because the intervals partition the histogram, so it is not ported as an error.
fn intervals_to_map(intervals: &[QualInterval], size: usize) -> Vec<u8> {
    let mut map = vec![0u8; size];
    for interval in intervals {
        let qual = interval.qual();
        for q in interval.q_start..=interval.q_end {
            if let Some(slot) = map.get_mut(q as usize) {
                *slot = qual;
            }
        }
    }
    map
}

/// `QuantizationInfo`: the map, the counts it came from, and how many levels it turned out to have.
#[derive(Debug, Clone, PartialEq)]
pub struct QuantizationInfo {
    pub quantized_quals: Vec<u8>,
    pub empirical_qual_counts: Vec<i64>,
    pub quantization_levels: i32,
}

impl QuantizationInfo {
    /// The two-argument constructor, which counts the levels out of the map it is given.
    pub fn new(quantized_quals: Vec<u8>, empirical_qual_counts: Vec<i64>) -> QuantizationInfo {
        let quantization_levels = calculate_quantization_levels(&quantized_quals);
        QuantizationInfo {
            quantized_quals,
            empirical_qual_counts,
            quantization_levels,
        }
    }

    /// `quantizeQualityScores(nLevels)`: replace the map, leaving the level count alone.
    ///
    /// The reference does **not** update `quantizationLevels` here. Only the constructor and
    /// [`QuantizationInfo::no_quantization`] set it, so after this call the count can describe a map
    /// that no longer exists.
    pub fn quantize_quality_scores(&mut self, n_levels: i32) -> Result<(), QuantizerError> {
        let quantizer =
            QualQuantizer::new(&self.empirical_qual_counts, n_levels, MIN_USABLE_Q_SCORE)?;
        self.quantized_quals = quantizer.original_to_quantized_map;
        Ok(())
    }

    /// `noQuantization()`: the identity map, for the **first 93 entries only**.
    ///
    /// The loop is `for (i = 0; i < quantizationLevels; i++)` after `quantizationLevels` was set to
    /// 93, so entry 93 keeps whatever it had. The golden carries the whole map, so the last entry is
    /// visible.
    pub fn no_quantization(&mut self) {
        self.quantization_levels = MAX_SAM_QUAL_SCORE;
        for i in 0..self.quantization_levels as usize {
            if let Some(slot) = self.quantized_quals.get_mut(i) {
                *slot = i as u8;
            }
        }
    }
}

/// `calculateQuantizationLevels`: a count of **changes**, not of distinct values.
///
/// It starts from -1, which no unsigned quality can equal, so the first entry always counts. A map
/// that returns to a value it already used counts that as another level: `2,2,10,10,2,2,30` is four.
pub fn calculate_quantization_levels(quantized_quals: &[u8]) -> i32 {
    let mut last: i32 = -1;
    let mut levels = 0;
    for qual in quantized_quals {
        if *qual as i32 != last {
            levels += 1;
            last = *qual as i32;
        }
    }
    levels
}

/// `QualityUtils.qualToProbLog10`, kept beside the rest of the quantizer's arithmetic.
pub fn qual_to_prob_log10(qual: u8) -> f64 {
    log10(1.0 - pow10(qual as f64 / -10.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_leaf_reads_its_error_rate_off_its_fixed_quality() {
        let leaf = QualInterval {
            q_start: 30,
            q_end: 30,
            // Counts that imply a completely different rate, and are ignored.
            n_observations: 1000,
            n_errors: 500,
            level: 0,
            fixed_qual: 30,
            sub_intervals: Vec::new(),
        };
        assert_eq!(leaf.error_rate(), qual_to_error_prob(30.0));
        assert_eq!(leaf.qual(), 30);
        // The same interval merged loses the fixed quality and the counts start to matter.
        let merged = QualInterval {
            fixed_qual: -1,
            sub_intervals: vec![leaf.clone(), leaf.clone()],
            ..leaf
        };
        assert_eq!(merged.error_rate(), 501.0 / 1001.0);
    }

    #[test]
    fn the_error_count_of_a_leaf_saturates_at_an_int() {
        // Three billion observations at quality zero: every base an error, and the count clamps.
        let quantizer = QualQuantizer::new(&[3_000_000_000, 3_000_000_000], 1, 6).unwrap();
        let interval = &quantizer.quantized_intervals[0];
        assert_eq!(interval.n_observations, 6_000_000_000);
        assert_eq!(interval.n_errors, 2 * i32::MAX as i64);
    }

    #[test]
    fn a_quantized_quality_of_zero_is_reachable() {
        // Five bins and sixteen levels: nothing merges, so every leaf keeps its own quality,
        // including zero, which `error_prob_to_qual` would have clamped to one.
        let quantizer = QualQuantizer::new(&[100; 5], 16, 6).unwrap();
        assert_eq!(quantizer.original_to_quantized_map, vec![0, 1, 2, 3, 4]);
        assert_eq!(error_prob_to_qual(1.0), Some(1));
    }

    #[test]
    fn the_two_ends_of_error_prob_to_qual() {
        assert_eq!(error_prob_to_qual(0.0), Some(93));
        assert_eq!(error_prob_to_qual(1.0), Some(1));
        assert_eq!(error_prob_to_qual(0.001), Some(30));
        assert_eq!(error_prob_to_qual(-0.1), None);
        assert_eq!(error_prob_to_qual(f64::NAN), None);
    }

    #[test]
    fn no_levels_is_the_missing_pair_and_no_histogram_is_the_missing_minimum() {
        assert_eq!(
            QualQuantizer::new(&[1; 10], 0, 6).unwrap_err(),
            QuantizerError::NoPairToMerge
        );
        assert_eq!(
            QualQuantizer::new(&[], 4, 6).unwrap_err(),
            QuantizerError::EmptyHistogram
        );
        assert_eq!(
            QualQuantizer::new(&[1, 1, -1], 4, 6).unwrap_err().message(),
            "Quality score histogram has negative values at: 1, 1, -1"
        );
    }

    #[test]
    fn the_level_count_counts_changes_and_not_values() {
        assert_eq!(calculate_quantization_levels(&[2, 2, 10, 10, 2, 2, 30]), 4);
        assert_eq!(calculate_quantization_levels(&[0, 0, 0]), 1);
        assert_eq!(calculate_quantization_levels(&[]), 0);
    }
}
