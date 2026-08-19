//! Ported from `org.broadinstitute.hellbender.tools.walkers.fasta.FastaReferenceMaker`
//! (GATK 4.6.2.0).
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
/// `bases_per_line` is `--line-width`, whose default is the writer's own.
pub fn run(
    reference: &mut ReferenceFileSource,
    arguments: &IntervalArguments,
    bases_per_line: usize,
) -> Result<FastaOutputs, MakerError> {
    // The writer is built in `onTraversalStart`, before a single locus is read, so a width of zero
    // is refused before the reference is touched.
    let mut writer = FastaReferenceWriter::new(bases_per_line, true).map_err(MakerError::Writer)?;

    let applied =
        reference_walker::traverse(reference, arguments, |locus: &SimpleInterval| locus.clone())
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
    // it, and the reference then writes a sequence with a null name -- which cannot happen here,
    // because a reference with no bases has no contigs either.
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
fn within_distance_of(left: &SimpleInterval, right: &SimpleInterval, distance: i32) -> bool {
    left.contig == right.contig
        && left.start <= right.end + distance
        && right.start <= left.end + distance
}

/// `FastaReferenceMaker.basesPerLine`'s default, which is the writer's.
pub const DEFAULT_LINE_WIDTH: usize = DEFAULT_BASES_PER_LINE;
