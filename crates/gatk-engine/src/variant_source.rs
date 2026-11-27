//! Ported from `org.broadinstitute.hellbender.engine.FeatureIntervalIterator` and the traversal
//! path of `org.broadinstitute.hellbender.engine.FeatureDataSource` (GATK 4.7.0.0; both files are
//! unchanged since 4.6.2.0).
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
//! One interval back, not a set of everything seen, so a variant covered by intervals 1 and 3 but
//! not by 2 would be handed to `apply` twice.
//!
//! **The golden says it is not.** `-L chr1:100-160 -L chr1:500-600 -L chr1:100-160` produces two
//! `apply` calls, not four, and `-L chr1:200-210 -L chr1:100-160` traverses in sorted order. The
//! argument layer ([`crate::interval_args`]) sorts and merges before the traversal ever sees the
//! list, so the precondition the iterator documents is always satisfied through a command line and
//! the one-interval memory is unreachable from there.
//!
//! That is worth stating exactly, because the two claims are different: the behaviour is real in
//! the class and unreachable in the tool. A caller that builds the list itself and skips the
//! argument layer still gets duplicates, and this port still reproduces them.
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
//!
//! # Reading this file without reading Rust
//!
//! - `Option<&[SimpleInterval]>` is "a borrowed, read-only list of intervals, or nothing": Java's
//!   nullable `List<SimpleInterval>`;
//! - `Vec<&'a T>` is a list of references into the caller's `features`, not copies. `'a` is a
//!   lifetime: it ties the result to the slice it points into and changes nothing at run time;
//! - `let Some(x) = y else { ... }` unwraps `y` or takes the early exit in the braces;
//! - a `trait` is a Java interface.

use crate::interval::SimpleInterval;

/// The location of one feature, which is all the traversal looks at: Java's `Locatable`, with
/// `stop` for `getEnd`. Coordinates are 1-based and inclusive, as Tribble returns them.
pub trait Located {
    fn contig(&self) -> &str;
    fn start(&self) -> i32;
    fn stop(&self) -> i32;
}

// A reference to something located is located too, so a caller can traverse a list of
// references without copying the records. `**self` follows the reference to the record; `?Sized`
// admits types whose size is not known at compile time (a `dyn Located`, for instance).
impl<T: Located + ?Sized> Located for &T {
    fn contig(&self) -> &str {
        (**self).contig()
    }
    fn start(&self) -> i32 {
        (**self).start()
    }
    fn stop(&self) -> i32 {
        (**self).stop()
    }
}

/// `FeatureDataSource.setIntervalsForTraversal`: null and empty are the same thing, and that thing
/// is "no restriction".
///
/// ```java
/// // Treat null and empty interval lists the same
/// intervalsForTraversal = (intervals != null && !intervals.isEmpty()) ? intervals : null;
/// ```
///
/// The reference's next check, a `UserException` when the source does not support random
/// access, has no counterpart: `traverse` can always filter the records it is handed.
pub fn intervals_for_traversal(intervals: Option<&[SimpleInterval]>) -> Option<&[SimpleInterval]> {
    // `Some(list) if ...` is a guarded match arm: it matches only when the guard holds, and `_`
    // catches everything else (null, or an empty list).
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
///
/// The reference is a lazy iterator; this collects the same sequence eagerly: interval by
/// interval in the order given, and within one interval in the file's order.
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
        // has no previous and every later one has exactly the interval before it:
        // `previousInterval = currentInterval;` then `currentInterval = intervalIterator.next();`
        // in `queryNextInterval`. `.enumerate()` pairs each interval with its index.
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
    // Discards `previous` and has no effect on the result. The reference's end-of-list reset,
    // `currentInterval = previousInterval = null;`, has no counterpart because this function
    // returns rather than being restarted.
    let _ = previous;

    // The final expression, with no `return`, is the function's value.
    emitted
}
