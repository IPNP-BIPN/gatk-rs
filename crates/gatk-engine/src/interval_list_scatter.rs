//! Ported from `picard.util.IntervalList.IntervalListScatterMode` and the six classes under it
//! (Picard 3.5.0, as GATK 4.7.0.0 bundles it; the nine files of `picard.util.IntervalList` are
//! unchanged since Picard 3.4.0, which GATK 4.6.2.0 bundled).
//!
//! The five ways an interval list is divided into shards, which is what `SplitIntervals` and every
//! scatter-gather pipeline built on it inherits.
//!
//! # Only one mode uniques its input
//!
//! `INTERVAL_SUBDIVISION`'s `preprocessIntervalList` is `uniqued()`; every other mode's is
//! `sorted()`. The same overlapping list therefore comes out as different INTERVALS under the
//! other four, not merely in different shards, and the names of a merged run are joined with a
//! pipe by the uniquing.
//!
//! # The ideal weight is taken from the unique base count either way
//!
//! ```java
//! return Math.max(1, (int) Math.floorDiv(intervalList.getUniqueBaseCount(), nCount));
//! ```
//!
//! `getUniqueBaseCount` uniques the list before summing, so a mode that never uniques its own
//! input still divides by a number that assumes the overlaps are gone: its shards then carry more
//! bases than the ideal says. The no-subdivision modes raise the ideal to the widest interval on
//! top of that, which is what makes a large scatter count come back with fewer shards than asked.
//!
//! # The last shard takes everything left
//!
//! The loop offers intervals only while `intervalsReturned < scatterCount`, and then flushes the
//! whole queue into the list it is building. The last shard is therefore unbounded, and a scatter
//! count of one is the whole list in one shard whatever the weights say.
//!
//! # The projection is a double division by the shards still to come
//!
//! The reference's `next()` computes, before each interval it offers,
//! `(weightRemaining - scatterer.listWeight(runningIntervalList))` divided by
//! `((double) scatterCount - intervalsReturned)`: the weight not yet placed, spread over the shards
//! still to come.
//!
//! Only the overflow mode and the distributed-remainder mode read it, and reading it is what
//! makes their decisions depend on how far into the scatter they are rather than on the list
//! alone.
//!
//! # Reading this file without reading Rust
//!
//! - The reference has one `IntervalListScatterer` class per mode; here the five are one `enum`,
//!   and each method `match`es on the mode where the reference dispatches on the class.
//!   `matches!(x, A | B)` is "x is one of A, B".
//! - `i64` is Java's `long`, `i32` its `int`. `i64::from(x)` widens an `int` exactly; `x as i32`
//!   truncates a `long` to its low 32 bits, as Java's `(int)` does.
//! - `(Option<Interval>, Option<Interval>)` is `takeSome`'s two-element list, with `None` for
//!   each `null`.

use htsjdk_bam::interval::{Interval, IntervalList};

/// `IntervalListScatterMode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScatterMode {
    /// `IntervalListScattererWithSubdivision`.
    IntervalSubdivision,
    /// `IntervalListScattererWithoutSubdivision`.
    BalancingWithoutIntervalSubdivision,
    /// `IntervalListScattererWithoutSubdivisionWithOverflow`.
    BalancingWithoutIntervalSubdivisionWithOverflow,
    /// `IntervalListScattererByIntervalCount`.
    IntervalCount,
    /// `IntervalListScattererByIntervalCountWithDistributedRemainder`.
    IntervalCountWithDistributedRemainder,
}

impl ScatterMode {
    /// The enum constant's own name.
    pub fn name(&self) -> &'static str {
        match self {
            ScatterMode::IntervalSubdivision => "INTERVAL_SUBDIVISION",
            ScatterMode::BalancingWithoutIntervalSubdivision => {
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION"
            }
            ScatterMode::BalancingWithoutIntervalSubdivisionWithOverflow => {
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
            }
            ScatterMode::IntervalCount => "INTERVAL_COUNT",
            ScatterMode::IntervalCountWithDistributedRemainder => {
                "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER"
            }
        }
    }

    /// The five modes in the order the enum declares them.
    pub const ALL: [ScatterMode; 5] = [
        ScatterMode::IntervalSubdivision,
        ScatterMode::BalancingWithoutIntervalSubdivision,
        ScatterMode::BalancingWithoutIntervalSubdivisionWithOverflow,
        ScatterMode::IntervalCount,
        ScatterMode::IntervalCountWithDistributedRemainder,
    ];

    /// Whether the weight of an interval is its bases or the interval itself: the three modes
    /// that extend `IntervalListScattererByBaseCount`, against the two by interval count.
    fn counts_bases(&self) -> bool {
        !matches!(
            self,
            ScatterMode::IntervalCount | ScatterMode::IntervalCountWithDistributedRemainder
        )
    }

    /// `preprocessIntervalList`: uniqued for the subdividing mode, sorted for every other.
    pub fn preprocess(&self, list: &IntervalList) -> IntervalList {
        match self {
            // `uniqued()` is `uniqued(true)`, which concatenates the names it merged.
            ScatterMode::IntervalSubdivision => list.uniqued(true),
            _ => list.sorted(),
        }
    }

    /// Whether the shards this mode produces carry `SO:coordinate` on their `@HD` line.
    ///
    /// A shard's header is the PREPROCESSED list's header, and the two preprocessing calls differ:
    /// `sorted()` stamps `SortOrder.coordinate` on the copy it returns, while `uniqued()` builds a
    /// fresh list around a clone of the ORIGINAL header and stamps nothing. So the four modes that
    /// only sort write a sort order into every file, and the one that uniques does not, from lists
    /// that are sorted either way.
    pub fn stamps_sort_order(&self) -> bool {
        !matches!(self, ScatterMode::IntervalSubdivision)
    }

    /// `intervalWeight`: `interval.length()` by bases, 1 by count.
    fn interval_weight(&self, interval: &Interval) -> i64 {
        if self.counts_bases() {
            i64::from(interval.length())
        } else {
            1
        }
    }

    /// `listWeight`: the base count, which is NOT the unique base count, or the interval count.
    /// (`intervalList.getBaseCount()` and `intervalList.size()`.)
    pub fn list_weight(&self, list: &IntervalList) -> i64 {
        if self.counts_bases() {
            list.intervals
                .iter()
                .map(|interval| i64::from(interval.length()))
                .sum()
        } else {
            list.intervals.len() as i64
        }
    }

    /// `deduceIdealSplitWeight`, on the list the mode has already preprocessed.
    ///
    /// The by-base-count form is the one the module note quotes. Its `(int)` cast sits *inside*
    /// the `max`: a quotient above `Integer.MAX_VALUE` keeps only its low 32 bits there, which for
    /// a quotient below 2^32 is negative, and the `max` then answers 1. This keeps the quotient in
    /// `i64` without the cast, so the two part for such a quotient. The interval-count modes cast
    /// after the `max` (`(int) Math.max(1, ...)`), on a count that cannot reach that size.
    pub fn ideal_split_weight(&self, list: &IntervalList, scatter_count: i32) -> i64 {
        let base = if self.counts_bases() {
            // `getUniqueBaseCount()`, which uniques first whatever the mode does.
            let unique: i64 = list
                .uniqued(true)
                .intervals
                .iter()
                .map(|interval| i64::from(interval.length()))
                .sum();
            // `1.max(x)` is `Math.max(1, x)`.
            1.max(floor_div(unique, i64::from(scatter_count)))
        } else {
            1.max(floor_div(self.list_weight(list), i64::from(scatter_count)))
        };
        match self {
            ScatterMode::BalancingWithoutIntervalSubdivision
            | ScatterMode::BalancingWithoutIntervalSubdivisionWithOverflow => {
                // `IntervalListScattererWithoutSubdivision.deduceIdealSplitWeight`:
                // `Math.max(widestIntervalLength, splitWidth)`.
                // `orElse(1)`: an empty list's widest interval is one, not zero.
                let widest = list
                    .intervals
                    .iter()
                    .map(|interval| i64::from(interval.length()))
                    .max()
                    .unwrap_or(1);
                widest.max(base)
            }
            _ => base,
        }
    }

    /// `takeSome`: what goes into the shard being built, and what goes back on the queue.
    ///
    /// Per mode, as the five classes have it:
    ///
    ///  * subdivision: all of it if it fits in `idealSplitWeight - currentSize`, nothing if that is
    ///    zero, otherwise the first `amount` bases and the rest;
    ///  * without subdivision: all or nothing by `shouldIncludeInterval`, which the overflow mode
    ///    widens with `|| idealSplitWeight < projectedSizeOfRemaining`;
    ///  * by interval count: one more while `amount > 0`;
    ///  * with distributed remainder: one more while `projectSizeOfRemaining > currentSize`, which
    ///    ignores the ideal weight entirely.
    fn take_some(
        &self,
        interval: &Interval,
        ideal: i64,
        current_size: i64,
        projected_remaining: f64,
    ) -> (Option<Interval>, Option<Interval>) {
        match self {
            ScatterMode::IntervalSubdivision => {
                let amount = ideal - current_size;
                if amount >= i64::from(interval.length()) {
                    return (Some(interval.clone()), None);
                }
                if amount == 0 {
                    return (None, Some(interval.clone()));
                }
                // The cut keeps the strand and the name on both halves.
                // `interval.getStart() + (int) amount - 1`, and `+ (int) amount` for the right.
                let mut left = interval.clone();
                left.end = interval.start + amount as i32 - 1;
                let mut right = interval.clone();
                right.start = interval.start + amount as i32;
                (Some(left), Some(right))
            }
            ScatterMode::BalancingWithoutIntervalSubdivision
            | ScatterMode::BalancingWithoutIntervalSubdivisionWithOverflow => {
                let projected = current_size + self.interval_weight(interval);
                let overflow = matches!(
                    self,
                    ScatterMode::BalancingWithoutIntervalSubdivisionWithOverflow
                );
                let include =
                    projected <= ideal || (overflow && (ideal as f64) < projected_remaining);
                if include {
                    (Some(interval.clone()), None)
                } else {
                    (None, Some(interval.clone()))
                }
            }
            ScatterMode::IntervalCount => {
                if ideal - current_size > 0 {
                    (Some(interval.clone()), None)
                } else {
                    (None, Some(interval.clone()))
                }
            }
            ScatterMode::IntervalCountWithDistributedRemainder => {
                if projected_remaining > current_size as f64 {
                    (Some(interval.clone()), None)
                } else {
                    (None, Some(interval.clone()))
                }
            }
        }
    }
}

/// `Math.floorDiv`, which rounds toward negative infinity rather than toward zero.
///
/// Rust's `/` on integers truncates toward zero, as Java's does; the correction subtracts one when
/// there is a remainder and the signs differ. The counts here are not negative, so it never fires
/// on the inputs this file passes, but it keeps the function the one it is named for.
fn floor_div(numerator: i64, denominator: i64) -> i64 {
    let quotient = numerator / denominator;
    if (numerator % denominator != 0) && ((numerator < 0) != (denominator < 0)) {
        quotient - 1
    } else {
        quotient
    }
}

/// What scattering refuses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScatterError {
    /// `IntervalListScatter`'s constructor.
    ScatterCountBelowOne,
}

impl ScatterError {
    /// The class of the exception the reference throws.
    pub fn java_class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }

    /// `throw new IllegalArgumentException("scatterCount < 1");`
    pub fn message(&self) -> String {
        "scatterCount < 1".to_string()
    }
}

/// `IntervalListScatterer.scatter`, which is the iterator run to exhaustion.
///
/// The reference preprocesses twice, in the `IntervalListScatter` constructor and again in
/// `ScatterState` on the already preprocessed list; this preprocesses once. Each pass of the outer
/// loop is one `ScatterState.next()`.
pub fn scatter(
    list: &IntervalList,
    mode: ScatterMode,
    scatter_count: i32,
) -> Result<Vec<IntervalList>, ScatterError> {
    if scatter_count < 1 {
        return Err(ScatterError::ScatterCountBelowOne);
    }
    let processed = mode.preprocess(list);
    let ideal = mode.ideal_split_weight(&processed, scatter_count);
    // `new ArrayDeque<>(processedIntervals.getIntervals())`: `pop_front` is `pollFirst` and
    // `push_front` is `addFirst`.
    let mut queue: std::collections::VecDeque<Interval> =
        processed.intervals.iter().cloned().collect();
    let mut weight_remaining = mode.list_weight(&processed);
    let mut intervals_returned: i64 = 0;
    let mut shards = Vec::new();

    while !queue.is_empty() {
        intervals_returned += 1;
        let mut running = IntervalList::new(processed.dictionary.clone());
        let mut closed = false;
        while !queue.is_empty() && intervals_returned < i64::from(scatter_count) {
            let interval = queue.pop_front().expect("the queue is not empty");
            let current_size = mode.list_weight(&running);
            // The `-1` in the reference's denominator is this `intervals_returned`, which was
            // incremented before the shard was started.
            let projected_remaining = (weight_remaining - current_size) as f64
                / (f64::from(scatter_count) - intervals_returned as f64);
            let (take, give_back) =
                mode.take_some(&interval, ideal, current_size, projected_remaining);
            if let Some(back) = give_back {
                queue.push_front(back);
            }
            // A `None` first half closes the shard: the reference's `return runningIntervalList`,
            // which `closed` and the `break` stand for.
            match take {
                None => {
                    weight_remaining -= mode.list_weight(&running);
                    closed = true;
                    break;
                }
                Some(interval) => running.intervals.push(interval),
            }
        }
        if closed {
            shards.push(running);
            continue;
        }
        // The flush: everything left goes into this shard, whatever its weight.
        while let Some(interval) = queue.pop_front() {
            running.intervals.push(interval);
        }
        // The reference throws `NoSuchElementException` for an empty shard here. The queue was
        // not empty when this pass began, so the shard cannot be empty and the test always holds.
        if !running.intervals.is_empty() {
            shards.push(running);
        }
    }
    Ok(shards)
}
