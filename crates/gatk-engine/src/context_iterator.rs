//! Ported from `org.broadinstitute.hellbender.utils.locusiterator.IntervalAlignmentContextIterator`,
//! `org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator` and
//! `org.broadinstitute.hellbender.utils.iterators.IntervalOverlappingIterator` (GATK 4.6.2.0).
//!
//! What [`crate::locus_iterator`] yields is one context per *covered* locus. A tool that asked for
//! empty loci wants one per *requested* locus, covered or not, and a tool that asked for intervals
//! without empty loci wants the covered ones filtered to those intervals. Those are two different
//! iterators wrapped around the same source, and which one a walker gets is decided by
//! `AlignmentContextIteratorBuilder`.
//!
//! Three behaviours here are decisions:
//!
//!  * **an empty locus carries an empty pileup, not an absent one.** `createEmptyAlignmentContext`
//!    builds a real context at the requested position whose pileup has no elements, so a tool sees
//!    a position with zero depth rather than not seeing the position;
//!  * **the source is advanced past, not rewound.** `advanceAlignmentContextToCurrentInterval`
//!    walks the covered contexts forward while they sit before the requested locus. A covered
//!    locus that no requested interval contains is therefore consumed and dropped, and it cannot
//!    come back if a later interval would have wanted it, which is why the intervals must be
//!    sorted and merged before they get here;
//!  * **the filter and the emitter disagree about what to do at the end.** When the source runs
//!    out, `IntervalOverlappingIterator` stops, while the empty-loci iterator keeps producing
//!    empty contexts until the requested loci run out.

use crate::interval::SimpleInterval;
use crate::locus_iterator::AlignmentContext;
use crate::read_pileup::ReadPileup;
use htsjdk_bam::header::SamHeader;

/// `IntervalUtils.compareLocatables`: contig index in the dictionary, then start, then end.
fn compare_locatables(
    left: (&str, i32, i32),
    right: (&str, i32, i32),
    header: &SamHeader,
) -> std::cmp::Ordering {
    let index = |contig: &str| {
        header
            .sequences
            .iter()
            .position(|s| s.name == contig)
            .unwrap_or(usize::MAX)
    };
    index(left.0)
        .cmp(&index(right.0))
        .then(left.1.cmp(&right.1))
        .then(left.2.cmp(&right.2))
}

/// `IntervalLocusIterator`: each interval split into its single-base loci, in order.
pub fn interval_loci(intervals: &[SimpleInterval]) -> Vec<SimpleInterval> {
    let mut loci = Vec::new();
    for interval in intervals {
        for position in interval.start..=interval.end {
            loci.push(SimpleInterval {
                contig: interval.contig.clone(),
                start: position,
                end: position,
            });
        }
    }
    loci
}

/// `IntervalOverlappingIterator`: the covered contexts that overlap the requested intervals.
///
/// The interval cursor only advances, as the source cursor does, so this is a merge of two sorted
/// streams rather than a containment test per context.
pub fn overlapping<'a>(
    contexts: Vec<AlignmentContext<'a>>,
    intervals: &[SimpleInterval],
    header: &SamHeader,
) -> Vec<AlignmentContext<'a>> {
    let mut kept = Vec::new();
    let mut interval_index = 0usize;
    for context in contexts {
        loop {
            let Some(interval) = intervals.get(interval_index) else {
                return kept;
            };
            let here = (context.contig.as_str(), context.position, context.position);
            let there = (interval.contig.as_str(), interval.start, interval.end);
            if interval.contig == context.contig
                && interval.start <= context.position
                && context.position <= interval.end
            {
                kept.push(context);
                break;
            }
            match compare_locatables(there, here, header) {
                // The interval is behind the context: take the next interval.
                std::cmp::Ordering::Less => interval_index += 1,
                // The interval is ahead: drop this context and take the next one.
                _ => break,
            }
        }
    }
    kept
}

/// `IntervalAlignmentContextIterator`: one context per requested locus, empty where uncovered.
pub fn with_empty_loci<'a>(
    contexts: Vec<AlignmentContext<'a>>,
    intervals: &[SimpleInterval],
    header: &SamHeader,
) -> Vec<AlignmentContext<'a>> {
    let loci = interval_loci(intervals);
    let mut out = Vec::with_capacity(loci.len());
    let mut source = contexts.into_iter().peekable();

    for locus in loci {
        // `advanceAlignmentContextToCurrentInterval`: consume every covered context that sits
        // before this locus. They are dropped, not buffered.
        while let Some(context) = source.peek() {
            let here = (context.contig.as_str(), context.position, context.position);
            let there = (locus.contig.as_str(), locus.start, locus.end);
            if compare_locatables(there, here, header) == std::cmp::Ordering::Greater {
                source.next();
            } else {
                break;
            }
        }

        let covered = source.peek().is_some_and(|context| {
            context.contig == locus.contig
                && locus.start <= context.position
                && context.position <= locus.end
        });
        if covered {
            out.push(source.next().expect("peeked"));
        } else {
            out.push(AlignmentContext {
                contig: locus.contig.clone(),
                position: locus.start,
                pileup: ReadPileup::new(&locus.contig, locus.start, Vec::new()),
            });
        }
    }
    out
}

/// `AlignmentContextIteratorBuilder.createAlignmentContextIterator`, as a routing decision.
///
/// `areIntervalsSpecified` is `intervals != null`, so `Some(vec![])` is *specified* while `None` is
/// not. That distinction is the whole reason this is an enum rather than a slice.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Route {
    /// No intervals and no empty loci: the covered contexts, unfiltered.
    Unbounded,
    /// Intervals without empty loci: filtered to them.
    Overlapping,
    /// Empty loci: one context per requested locus.
    EmptyLoci,
    /// `Utils.nonEmpty` in `IntervalOverlappingIterator`'s constructor: intervals were specified
    /// as an empty list, and empty loci were not asked for.
    RejectedEmptyIntervalList,
}

/// Which iterator the builder returns.
pub fn route(emit_empty_loci: bool, intervals: Option<&[SimpleInterval]>) -> Route {
    if emit_empty_loci {
        // With no intervals at all the builder substitutes the whole reference, so the route is
        // the same either way.
        return Route::EmptyLoci;
    }
    match intervals {
        None => Route::Unbounded,
        Some([]) => Route::RejectedEmptyIntervalList,
        Some(_) => Route::Overlapping,
    }
}
