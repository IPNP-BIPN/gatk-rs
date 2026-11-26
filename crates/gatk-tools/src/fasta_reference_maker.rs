//! Ported from `org.broadinstitute.hellbender.tools.walkers.fasta.FastaReferenceMaker`
//! (GATK 4.7.0.0; `FastaReferenceMaker.java` is unchanged since 4.6.2.0).
//!
//! A `ReferenceWalker` whose `apply` appends one base, and whose output is a FASTA written through
//! [`htsjdk_bam::fasta_writer`] with its index and dictionary beside it.
//!
//! # The output sequences are numbered, and the contig is a description
//!
//! `appendSequence(String.valueOf(contigCount), description, basesPerLine, bases)`: the NAME is a
//! counter starting at one and the contig and coordinates are the description, so a FASTA made from
//! `chr1:1-12` holds a sequence called `1` whose header reads `>1 chr1:1-12`. Anything downstream
//! that expected the contig's name gets an integer.
//!
//! # A gap starts a new sequence, and so does a contig boundary
//!
//! `advancePosition` compares each locus with the last through
//! `lastPosition.withinDistanceOf(interval, 1)`, which is false when the contig differs as well as
//! when the positions are more than one apart. So two abutting intervals are one output sequence,
//! two a base apart are two, and a run with no `-L` over a two-contig reference writes two.
//!
//! # The description is the span written, not the interval asked for
//!
//! `currentSequenceStartPosition` is the first locus of the sequence and the end is the last locus
//! seen, so the header describes what was written. With intervals that is the same thing; with a
//! gap in the middle it is not.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, MakerError>` is "a value or an error"; it is how a Java exception is modelled.
//!   `?` after a call returns its error to the caller at once, a `throws` at one call site, and
//!   `.map_err(MakerError::Writer)` first wraps the callee's error in this file's error type.
//! - `Option<SimpleInterval>` is a value or nothing (`Some(x)` or `None`); it stands for Java's
//!   nullable `lastPosition`, with the compiler forcing every use to say what happens on `None`.
//! - `&` before a type or value is a borrowed reference: the callee reads it and the caller keeps
//!   it. `&mut` is a borrowed reference the callee may write through.
//! - A body's last expression, with no `;`, is its return value.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::interval_args::IntervalArguments;
use gatk_engine::reference::ReferenceFileSource;
use htsjdk_bam::fasta_writer::{
    FastaOutputs, FastaReferenceWriter, FastaWriterError, DEFAULT_BASES_PER_LINE,
};

use crate::reference_walker::{self, TraversalError};

/// What stopped the tool.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MakerError {
    /// The traversal did not start, or the reference refused a window.
    Traversal(TraversalError),
    /// The writer refused, which is where a line width of zero lands.
    Writer(FastaWriterError),
}

/// `doWork`: the FASTA, the `.fai` and the `.dict`.
///
/// **How**: resolve `-L`/`-XL` against the reference's own dictionary, then [`run_over`].
///
/// `bases_per_line` is `--line-width`, whose default is the writer's own.
pub fn run(
    reference: &mut ReferenceFileSource,
    arguments: &IntervalArguments,
    bases_per_line: usize,
) -> Result<FastaOutputs, MakerError> {
    let intervals =
        reference_walker::traversal_intervals(arguments, &reference_walker::dictionary(reference))
            .map_err(MakerError::Traversal)?;
    run_over(reference, &intervals, bases_per_line)
}

/// What the writer leaves behind when the traversal refuses.
///
/// The writer is built in `onTraversalStart` and closed in `closeTool`, so a run that dies in the
/// middle of the traversal still writes its three files: the FASTA and the `.fai` are empty, and
/// the dictionary holds nothing but its `@HD` line. This is that state, produced by the same
/// writer rather than spelled out, so it cannot drift from what a real close would emit.
///
/// `closeTool` swallows the writer's `IllegalStateException`, which is why the files survive:
///
/// ```java
/// } catch (IllegalStateException e){
///     //sink this
/// ```
pub fn empty_outputs(bases_per_line: usize) -> Result<FastaOutputs, MakerError> {
    let writer = FastaReferenceWriter::new(bases_per_line, true).map_err(MakerError::Writer)?;
    // `close_streams` and not `close`: an empty reference is refused by `close()`, and what is
    // wanted here is the three files the reference's `finally` clause leaves behind, not that
    // refusal. The tool has already refused for its own reason by the time this is called.
    writer.close_streams().map_err(MakerError::Writer)
}

/// The same run over intervals somebody else resolved.
///
/// It exists for the reason [`reference_walker::traverse_intervals`] does: a command line resolves
/// its intervals against the BEST available dictionary, which a `--sequence-dictionary` outranks
/// the reference in, and the traversal then queries the FASTA with them.
///
/// **How**: the reference's `apply`/`advancePosition` state machine, run over the list of loci the
/// traversal returns rather than as a callback: the closure passed to `traverse_intervals`
/// (`|locus| locus.clone()`, Java's lambda) only records each locus, and the loop below then does
/// what `addToReference` does to each, in the same order.
pub fn run_over(
    reference: &mut ReferenceFileSource,
    intervals: &[SimpleInterval],
    bases_per_line: usize,
) -> Result<FastaOutputs, MakerError> {
    // The writer is built in `onTraversalStart`, before a single locus is read, so a width of zero
    // is refused before the reference is touched.
    let mut writer = FastaReferenceWriter::new(bases_per_line, true).map_err(MakerError::Writer)?;

    let applied =
        reference_walker::traverse_intervals(reference, intervals, |locus: &SimpleInterval| {
            locus.clone()
        })
        .map_err(MakerError::Traversal)?;

    let mut count = 0usize;
    let mut last: Option<SimpleInterval> = None;
    let mut start_position = 0;
    let mut sequence: Vec<u8> = Vec::new();

    for call in &applied {
        let interval = &call.window;
        // `advancePosition`: a first locus opens a sequence, and a locus that is not within one of
        // the last closes the sequence before opening the next.
        let new_sequence = match &last {
            None => true,
            Some(previous) => !within_distance_of(previous, interval, 1),
        };
        if new_sequence {
            if last.is_some() {
                finalize(
                    &mut writer,
                    count,
                    &last,
                    start_position,
                    &sequence,
                    bases_per_line,
                )?;
            }
            count += 1;
            start_position = interval.start;
            sequence.clear();
        }
        last = Some(interval.clone());
        // `getBase()`, which at the default window is the window's only byte.
        sequence.push(call.bases[0]);
    }

    // `onTraversalSuccess` finalizes whatever is open. A traversal that applied nothing still calls
    // it, and the reference's `finalizeSequence` then dereferences a null `lastPosition`
    // (`lastPosition.getContig()`), which cannot happen here, because a reference with no bases
    // has no contigs either.
    if last.is_some() {
        finalize(
            &mut writer,
            count,
            &last,
            start_position,
            &sequence,
            bases_per_line,
        )?;
    }

    writer.close().map_err(MakerError::Writer)
}

/// `finalizeSequence`: the numbered name, the description, and the bases collected so far.
///
/// The description, broken here to fit:
///
/// ```java
/// final String description = lastPosition.getContig() + ":" + currentSequenceStartPosition
///         + "-" + lastPosition.getEnd();
/// ```
///
/// `.expect(..)` unwraps the `Option` and panics (aborts, like an unchecked Java exception) on
/// `None`; both callers check `last.is_some()` first.
fn finalize(
    writer: &mut FastaReferenceWriter,
    count: usize,
    last: &Option<SimpleInterval>,
    start_position: i32,
    sequence: &[u8],
    bases_per_line: usize,
) -> Result<(), MakerError> {
    let last = last.as_ref().expect("a sequence is open");
    let description = format!("{}:{}-{}", last.contig, start_position, last.end);
    writer
        .start_sequence_with(&count.to_string(), &description, bases_per_line)
        .map_err(MakerError::Writer)?;
    writer.append_bases(sequence).map_err(MakerError::Writer)
}

/// `Locatable.withinDistanceOf`, which is false across contigs whatever the distance.
///
/// htsjdk 5.0.0 `Locatable.java` (the second line is broken here to fit):
///
/// ```java
/// return contigsMatch(other)
///         && CoordMath.overlaps(getStart(), getEnd(),
///                 other.getStart() - distance, other.getEnd() + distance);
/// ```
///
/// `CoordMath.overlaps` tests whether either end of the second span falls inside the first or
/// the second encloses the first; for the single-base loci compared here (start equals end) that
/// is the two inequalities below, with `distance` moved to the other side of the second.
fn within_distance_of(left: &SimpleInterval, right: &SimpleInterval, distance: i32) -> bool {
    left.contig == right.contig
        && left.start <= right.end + distance
        && right.start <= left.end + distance
}

/// `FastaReferenceMaker.basesPerLine`'s default, which is the writer's.
pub const DEFAULT_LINE_WIDTH: usize = DEFAULT_BASES_PER_LINE;
