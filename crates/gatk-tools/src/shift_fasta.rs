//! Ported from `org.broadinstitute.hellbender.tools.walkers.fasta.ShiftFasta` (GATK 4.6.2.0).
//!
//! Rotates every contig of a reference and writes four files: the shifted FASTA (with its index and
//! dictionary, from [`htsjdk_bam::fasta_writer`]), a chain file that shifts coordinates back, and a
//! pair of interval lists naming the region each end of the shift covers.
//!
//! It is not a walker. `traverse` is overridden outright and loops over the sequence dictionary, so
//! `-L` never reaches it and every contig is read whole.
//!
//! # A contig that is not shifted is dropped
//!
//! `shiftContig` does its work only when `0 < offset < contigLength`. Otherwise it logs and returns,
//! so that contig appears in NO output file: not unshifted in the FASTA, not in the index, not in
//! the dictionary, not in the chain. A shifted reference can therefore have fewer contigs than the
//! reference it came from, which is the sort of thing a caller notices later rather than sooner.
//!
//! # The two halves are appended separately
//!
//! `appendBases(basesAtEnd).appendBases(basesAtStart)` is two calls, and the writer breaks lines on
//! a running count rather than per call, so the join between the halves falls mid-line. Writing the
//! concatenation in one call would produce the same bytes; writing a newline between them would
//! not, and this records which of those the reference does.

use htsjdk_bam::fasta_writer::{FastaOutputs, FastaReferenceWriter, FastaWriterError};

use gatk_engine::reference::{ReferenceError, ReferenceFileSource};

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShiftError {
    /// `UserException.BadInput`: an offset list that is not one per contig.
    BadOffsetList { given: usize, contigs: usize },
    /// The reference could not answer.
    Reference(ReferenceError),
    /// The writer refused.
    Writer(FastaWriterError),
}

impl ShiftError {
    /// The exception class the reference throws.
    pub fn java_class(&self) -> &'static str {
        match self {
            ShiftError::BadOffsetList { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
            ShiftError::Reference(_) => "org.broadinstitute.hellbender.exceptions.UserException",
            ShiftError::Writer(error) => error.java_class(),
        }
    }

    /// The message, which for the bad list carries the `Bad input: ` prefix the exception adds.
    pub fn message(&self) -> String {
        match self {
            ShiftError::BadOffsetList { given, contigs } => format!(
                "Bad input: Shift offset list size {given} must equal number of contigs in the \
                 reference {contigs}"
            ),
            ShiftError::Reference(error) => format!("{error:?}"),
            ShiftError::Writer(error) => error.message(),
        }
    }
}

/// The four outputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShiftOutputs {
    /// The shifted FASTA and the two files the writer produces beside it.
    pub reference: FastaOutputs,
    /// The shift-back chain file.
    pub chain: String,
    /// The unshifted interval list.
    pub intervals: String,
    /// The shifted interval list.
    pub shifted_intervals: String,
}

/// `doWork`, on a reference already open.
///
/// `offsets` is `--shift-offset-list`: empty means half of each contig's own length.
pub fn run(
    reference: &mut ReferenceFileSource,
    offsets: &[i32],
    bases_per_line: usize,
) -> Result<ShiftOutputs, ShiftError> {
    let contigs: Vec<(String, usize)> = reference.sequences().to_vec();
    if !offsets.is_empty() && offsets.len() != contigs.len() {
        return Err(ShiftError::BadOffsetList {
            given: offsets.len(),
            contigs: contigs.len(),
        });
    }

    let mut writer = FastaReferenceWriter::new(bases_per_line, true).map_err(ShiftError::Writer)?;
    let mut outputs = ShiftOutputs {
        reference: FastaOutputs::default(),
        chain: String::new(),
        intervals: String::new(),
        shifted_intervals: String::new(),
    };
    // `chainId` starts at one and counts across contigs, so the second contig's records are 3 and 4.
    let mut chain_id = 1;

    for (index, (name, length)) in contigs.iter().enumerate() {
        let contig_length = *length as i32;
        let offset = if offsets.is_empty() {
            contig_length / 2
        } else {
            offsets[index]
        };
        // The whole of `shiftContig` is inside this test; a contig outside it is skipped, not
        // copied.
        if offset <= 0 || offset >= contig_length {
            continue;
        }

        let bases = reference
            .query(name, 1, contig_length)
            .map_err(ShiftError::Reference)?;
        let (start, end) = bases.split_at(offset as usize);
        let shift_back_offset = bases.len() as i32 - offset;

        // `addToShiftedReference`: the tail, then the head, as two appends.
        writer
            .start_sequence_with(name, "", bases_per_line)
            .map_err(ShiftError::Writer)?;
        writer.append_bases(end).map_err(ShiftError::Writer)?;
        writer.append_bases(start).map_err(ShiftError::Writer)?;

        // `addToChainFile`: two records, each followed by its own length and a blank line.
        outputs.chain.push_str(&chain_string(
            name,
            shift_back_offset,
            contig_length,
            offset,
            bases.len() as i32,
            0,
            shift_back_offset,
            chain_id,
        ));
        chain_id += 1;
        outputs
            .chain
            .push_str(&format!("\n{shift_back_offset}\n\n"));
        outputs.chain.push_str(&chain_string(
            name,
            offset - 1,
            contig_length,
            0,
            offset,
            shift_back_offset,
            bases.len() as i32,
            chain_id,
        ));
        chain_id += 1;
        outputs.chain.push_str(&format!("\n{offset}\n\n"));

        // `addToIntervalFiles`: the same start for both, and an end that differs by the parity.
        let interval_start = offset / 2;
        let interval_end = interval_start + contig_length / 2 - 1;
        outputs
            .intervals
            .push_str(&format!("{name}:{interval_start}-{interval_end}\n"));
        outputs.shifted_intervals.push_str(&format!(
            "{name}:{interval_start}-{}\n",
            interval_end + contig_length % 2
        ));
    }

    outputs.reference = writer.close().map_err(ShiftError::Writer)?;
    Ok(outputs)
}

/// `createChainString`: thirteen tab-separated fields, the contig named twice.
#[allow(clippy::too_many_arguments)]
fn chain_string(
    name: &str,
    score: i32,
    length: i32,
    start: i32,
    end: i32,
    shift_back_start: i32,
    shift_back_end: i32,
    id: i32,
) -> String {
    [
        "chain".to_string(),
        score.to_string(),
        name.to_string(),
        length.to_string(),
        "+".to_string(),
        shift_back_start.to_string(),
        shift_back_end.to_string(),
        name.to_string(),
        length.to_string(),
        "+".to_string(),
        start.to_string(),
        end.to_string(),
        id.to_string(),
    ]
    .join("\t")
}
