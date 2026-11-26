//! Ported from `org.broadinstitute.hellbender.engine.IntervalWalker` and the parts of
//! `org.broadinstitute.hellbender.engine.GATKTool` it runs on (GATK 4.6.2.0).
//!
//! The traversal itself is a loop over `userIntervals`, so the port is thin. What it stands on is
//! [`gatk_engine::interval_args`], which turns `-L`, `-XL` and the padding, set and merging rules
//! into that list.
//!
//! Two facts belong here rather than there, because they are the walker's and not the parser's:
//!
//!  * **`requiresIntervals()` is true**, which makes `-L` a *required* Barclay argument. So the
//!    `parseIntervals` branch that treats an absent `-L` as "the whole reference, minus `-XL`" is
//!    unreachable from this walker: `-XL chr1` on its own is rejected as a missing argument before
//!    any interval is parsed. The golden records exactly that, and the branch is still ported
//!    because tools with an optional interval collection do reach it;
//!  * **each interval gets a windowless `ReferenceContext`**, built from the interval rather than
//!    from a read, so a tool that wants flanking bases has to widen the window itself.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, E>` is "a value or an error"; [`TraversalError`] is how this file models the
//!   exceptions the Java would throw, and `?` after a call hands its error straight back to the
//!   caller (Java's `throws`, at that one call site);
//! - `Option<T>` is "a value or nothing", Java's nullable reference checked by the compiler;
//! - `&T` is a read-only borrow (the caller keeps ownership) and `&mut T` a writable one, of
//!   which only one may exist at a time.

use gatk_engine::context::ReferenceContext;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::interval_args::{self, IntervalArgumentError, IntervalArguments};
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_engine::reference::ReferenceFileSource;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// One `apply` call: the interval, the reads overlapping it, and the reference under it.
///
/// The Java hands these to the tool's `apply` as it goes; here they are returned as a list so the
/// runner decides what each tool does with them. No `FeatureContext` is built.
pub struct Applied {
    pub interval: SimpleInterval,
    pub reads: Vec<BamRecord>,
    pub context: ReferenceContext,
}

/// What stopped a traversal before it began.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TraversalError {
    /// Barclay's `CommandLineException.MissingArgument`: `requiresIntervals()` and no `-L`.
    MissingIntervalArgument,
    /// The interval arguments did not resolve.
    Intervals(IntervalArgumentError),
    /// The reads could not be read.
    Reads(ReadsError),
}

/// `IntervalArgumentCollection` for a walker whose `requiresIntervals()` is true.
///
/// ```java
/// public boolean requiresIntervals() {
///     return true;
/// }
/// ```
///
/// The check is Barclay's, not the parser's, and it fires on `-L` being absent whatever `-XL`
/// says: an empty required argument is a missing argument.
pub fn traversal_intervals(
    arguments: &IntervalArguments,
    header: &SamHeader,
) -> Result<Vec<SimpleInterval>, TraversalError> {
    if arguments.include.is_empty() {
        return Err(TraversalError::MissingIntervalArgument);
    }
    let parameters =
        interval_args::parse_intervals(arguments, header).map_err(TraversalError::Intervals)?;
    // `traverseUnmappedReads` has nowhere to go here: IntervalWalker walks intervals, and the
    // unmapped request was already separated out of the list. `-L unmapped` alone therefore runs
    // a traversal of zero intervals rather than failing.
    Ok(parameters.intervals)
}

/// `IntervalWalker.traverse`: one `apply` per interval, in the parsed order.
///
/// ```java
/// for ( final SimpleInterval interval : userIntervals ) {
///     apply(interval,
///           new ReadsContext(reads, interval, readFilter),
///           new ReferenceContext(reference, interval),
///           new FeatureContext(features, interval));
/// ```
///
/// `filter` is `makeReadFilter()`, passed in as a closure (`&dyn Fn`, a borrowed function value,
/// Java's `Predicate<GATKRead>`). The reference is optional because the Java's `reference` is
/// null when no `-R` was given.
pub fn traverse(
    source: &ReadsDataSource,
    reference: Option<&mut ReferenceFileSource>,
    arguments: &IntervalArguments,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<Vec<Applied>, TraversalError> {
    let header = source.header().clone();
    let intervals = traversal_intervals(arguments, &header)?;

    let mut applied = Vec::with_capacity(intervals.len());
    // Rebound as `mut` so each iteration can reborrow it: `as_deref_mut` below lends the
    // reference out for one interval and gets it back, where moving it would use it up.
    let mut reference = reference;
    for interval in intervals {
        // `new ReadsContext(reads, interval, readFilter)`: the reads overlapping this interval
        // alone, filtered, and re-queried per interval rather than streamed once. The walker's own
        // documentation calls the lack of caching out: "The current implementation uses no
        // caching". `map_err` converts the reader's error into this file's; `.filter(..)
        // .collect()` is Java's `stream().filter(..).collect(toList())`.
        let reads: Vec<BamRecord> = source
            .query(std::slice::from_ref(&interval))
            .map_err(TraversalError::Reads)?
            .into_iter()
            .filter(|read| filter(read))
            .collect();

        // `new ReferenceContext(reference, interval)`: leading and trailing window of zero.
        // `unwrap_or_else` turns a failed context into an empty one rather than an error.
        let context = match reference.as_deref_mut() {
            Some(reference) => ReferenceContext::new(reference, Some(interval.clone()), 0, 0),
            None => ReferenceContext::without_source(Some(interval.clone()), 0, 0),
        }
        .unwrap_or_else(|_| ReferenceContext::empty());

        applied.push(Applied {
            interval,
            reads,
            context,
        });
    }
    Ok(applied)
}
