//! Ported from `org.broadinstitute.hellbender.engine.FeatureIntervalIterator` and the traversal
//! path of `org.broadinstitute.hellbender.engine.FeatureDataSource` (GATK 4.6.2.0).
//!
//! Which variants a `VariantWalker` sees, and in what order. Two behaviours decide it and neither
//! is what "restrict the traversal to `-L`" suggests.
//!
//! # De-duplication remembers **one** interval, not all of them
//!
//! A feature is emitted when it is *novel*, and novel means "does not overlap the **previous**
//! interval":
//!
//! ```java
//! private boolean featureIsNovel( final T feature ) {
//!     return previousInterval == null || ! previousInterval.overlaps(new SimpleInterval(feature));
//! }
//! ```
//!
//! One interval back, not a set of everything seen. So a variant covered by intervals 1 and 3 but
//! not by 2 is handed to `apply` **twice**, and the same variant covered by two adjacent intervals
//! is handed over once. The class documents that intervals must be non-overlapping and sorted, and
//! nothing enforces it: a caller that breaks the precondition gets duplicates rather than an error.
//!
//! # `previousInterval` lags by one query, including at the end
//!
//! It is assigned from `currentInterval` at the *start* of each query, so during the first
//! interval it is null and during interval *n* it is interval *n-1*. When the intervals run out
//! both are cleared, which matters only if the iterator is restarted.
//!
//! # An empty interval list is not an empty traversal
//!
//! `setIntervalsForTraversal` maps both null and an empty list to null, which means *no*
//! restriction: passing an empty list traverses the whole file rather than nothing.

use crate::interval::SimpleInterval;

/// The location of one feature, which is all the traversal looks at.
pub trait Located {
    fn contig(&self) -> &str;
    fn start(&self) -> i32;
    fn stop(&self) -> i32;
}

/// `FeatureDataSource.setIntervalsForTraversal`: null and empty are the same thing, and that thing
/// is "no restriction".
pub fn intervals_for_traversal(intervals: Option<&[SimpleInterval]>) -> Option<&[SimpleInterval]> {
    match intervals {
        Some(list) if !list.is_empty() => Some(list),
        _ => None,
    }
}

/// `FeatureDataSource.iterator()`: the whole file, or `FeatureIntervalIterator` over the intervals.
///
/// `features` stands in for the indexed reader: it is the file's records in file order, and the
/// per-interval query is the subset overlapping that interval, in the same order. That is what a
/// Tribble query returns for a sorted, indexed file, and the index itself is htsjdk's problem
/// rather than the traversal's.
pub fn traverse<'a, T: Located>(
    features: &'a [T],
    intervals: Option<&[SimpleInterval]>,
) -> Vec<&'a T> {
    let Some(intervals) = intervals_for_traversal(intervals) else {
        return features.iter().collect();
    };

    let mut emitted: Vec<&'a T> = Vec::new();
    let mut previous: Option<&SimpleInterval> = None;

    for (index, current) in intervals.iter().enumerate() {
        // `previousInterval = currentInterval` happens *before* the query, so the first interval
        // has no previous and every later one has exactly the interval before it.
        previous = if index == 0 {
            None
        } else {
            Some(&intervals[index - 1])
        };

        for feature in features {
            if !current.overlaps(feature.contig(), feature.start(), feature.stop()) {
                continue;
            }
            // The novelty test: one interval of memory, so a feature can be emitted more than once
            // when the intervals are not the sorted, non-overlapping list the class asks for.
            let novel = match previous {
                None => true,
                Some(before) => !before.overlaps(feature.contig(), feature.start(), feature.stop()),
            };
            if novel {
                emitted.push(feature);
            }
        }
    }
    let _ = previous;

    emitted
}
