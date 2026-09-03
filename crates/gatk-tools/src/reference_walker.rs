//! Ported from `org.broadinstitute.hellbender.engine.ReferenceWalker` and the parts of
//! `org.broadinstitute.hellbender.engine.GATKTool` it runs on (GATK 4.6.2.0).
//!
//! Thirty lines of traversal standing on the interval machinery. What makes it worth its own module
//! is the branch it reaches that no other walker in this port does.
//!
//! # An absent `-L` is the whole reference
//!
//! `getTraversalIntervals` is `hasUserSuppliedIntervals() ? userIntervals : hasReference() ?
//! getAllIntervalsForReference(dictionary) : null`. [`crate::interval_walker`] cannot get there:
//! its `requiresIntervals()` is true, so Barclay rejects a run without `-L` before any interval is
//! parsed. `ReferenceWalker`'s is false, so all three arms are live here, and the golden measures
//! the first two.
//!
//! `-XL` with no `-L` is therefore legal as well, and goes through
//! [`gatk_engine::interval_args::parse_intervals`], whose include set is the whole reference for
//! exactly that case.
//!
//! # Every locus is one base, and the window is a method
//!
//! The interval list is cut by `IntervalLocusIterator`, which is `ShardedIntervalIterator` at shard
//! size one, so an eleven-base interval is eleven `apply` calls. Each call gets a
//! `ReferenceContext` built from the locus and from `getReferenceWindow(locus)` -- a METHOD, not an
//! argument, so a tool that wants flanking bases overrides it and still sees one locus at a time.
//! The window it returns is what decides which bases, reads and features `apply` sees.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::interval_args::{self, IntervalArgumentError, IntervalArguments};
use gatk_engine::locus_shards::interval_loci;
use gatk_engine::reference::{ReferenceError, ReferenceFileSource};
use htsjdk_bam::header::{SamHeader, SequenceRecord};

/// One `apply` call: the window it was given, and the bases under it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Applied {
    /// `referenceContext.getWindow()`, which is the locus unless the walker widened it.
    pub window: SimpleInterval,
    /// `referenceContext.getBases()`, already upper-cased with IUPAC codes flattened to `N`.
    pub bases: Vec<u8>,
}

/// What stopped a traversal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TraversalError {
    /// The interval arguments did not resolve.
    Intervals(IntervalArgumentError),
    /// The reference could not answer.
    Reference(ReferenceError),
}

/// `getReferenceDictionary()` for a tool whose only data source is the reference.
///
/// A `SamHeader` because that is what the interval parser reads a dictionary out of; nothing else
/// on it is set, and nothing else is looked at.
pub fn dictionary(reference: &ReferenceFileSource) -> SamHeader {
    let mut header = SamHeader::new();
    for (name, length) in reference.sequences() {
        header
            .sequences
            .push(SequenceRecord::new(name, *length as i32));
    }
    header
}

/// `GATKTool.getTraversalIntervals`.
///
/// With no interval argument at all this is `getAllIntervalsForReference`, one interval per contig
/// covering all of it, in dictionary order. With any argument -- including `-XL` alone -- it is the
/// parser's answer.
pub fn traversal_intervals(
    arguments: &IntervalArguments,
    header: &SamHeader,
) -> Result<Vec<SimpleInterval>, TraversalError> {
    if !arguments.specified() {
        return Ok(header
            .sequences
            .iter()
            .map(|sequence| {
                SimpleInterval::new(&sequence.name, 1, sequence.length)
                    .expect("a contig length is at least one")
            })
            .collect());
    }
    let parameters =
        interval_args::parse_intervals(arguments, header).map_err(TraversalError::Intervals)?;
    Ok(parameters.intervals)
}

/// `ReferenceWalker.traverse`: one `apply` per base, in interval order.
///
/// `window` is `getReferenceWindow`, whose default is the locus itself. It is passed in rather than
/// overridden because the reference makes it a method on the walker and the tools that widen it are
/// the exception: `CountBasesInReference` takes the default.
pub fn traverse(
    reference: &mut ReferenceFileSource,
    arguments: &IntervalArguments,
    window: impl Fn(&SimpleInterval) -> SimpleInterval,
) -> Result<Vec<Applied>, TraversalError> {
    let header = dictionary(reference);
    let intervals = traversal_intervals(arguments, &header)?;
    traverse_intervals(reference, &intervals, window)
}

/// The same traversal over intervals somebody else resolved.
///
/// It exists because the dictionary the intervals resolve against is NOT always the reference's:
/// `getBestAvailableSequenceDictionary` prefers a `--sequence-dictionary` over it, so a run with
/// both resolves `-L` against the master and then queries the FASTA, which is how the reference
/// answers `Given reference file does not have data at the requested contig` rather than refusing
/// the interval. A traversal that resolved its own intervals could not reach that.
pub fn traverse_intervals(
    reference: &mut ReferenceFileSource,
    intervals: &[SimpleInterval],
    window: impl Fn(&SimpleInterval) -> SimpleInterval,
) -> Result<Vec<Applied>, TraversalError> {
    let intervals = intervals.to_vec();
    let mut applied = Vec::new();
    for locus in interval_loci(&intervals) {
        let reference_window = window(&locus);
        let bases = reference
            .query(
                &reference_window.contig,
                reference_window.start,
                reference_window.end,
            )
            .map_err(TraversalError::Reference)?;
        applied.push(Applied {
            window: reference_window,
            bases,
        });
    }
    Ok(applied)
}
