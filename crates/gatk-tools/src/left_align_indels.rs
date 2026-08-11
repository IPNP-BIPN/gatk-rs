//! Ported from `org.broadinstitute.hellbender.tools.LeftAlignIndels` (GATK 4.6.2.0).
//!
//! The eighth whole tool of the record-transform archetype, and the first whose
//! `requiresReference()` is true.
//!
//! ```java
//! public void apply(GATKRead read, ReferenceContext ref, FeatureContext featureContext) {
//!     if (read.isUnmapped() || read.numCigarElements() <= 1) {
//!         outputWriter.addRead(read);
//!         return;
//!     }
//!     final CigarBuilder.Result result =
//!         AlignmentUtils.leftAlignIndels(read.getCigar(), ref.getBases(), read.getBases(), 0);
//!     read.setCigar(result.getCigar());
//!     if (result.getLeadingDeletionBasesRemoved() > 0) {
//!         read.setPosition(read.getContig(), read.getStart() + result.getLeadingDeletionBasesRemoved());
//!     }
//!     outputWriter.addRead(read);
//! }
//! ```
//!
//! # The window is the read, not the contig
//!
//! A `ReadWalker` builds `new ReferenceContext(reference, new SimpleInterval(read))`, so the bases
//! `apply` receives are the read's own span with no padding, and the read start passed into
//! [`gatk_engine::alignment_utils::left_align_indels`] is 0 because the window already starts at
//! the read. An indel can therefore only be moved as far as the read itself reaches. A port that
//! queried the contig instead would left-align further than the reference does and produce a file
//! that looks healthier and is wrong.
//!
//! # A deletion moves the read; an insertion does not
//!
//! Measured: `4M1D5M` at `chr1:6` over a homopolymer comes back as `9M` at `chr1:7`. The deletion
//! walked to the front of the window, [`gatk_engine::cigar_builder::CigarBuilder`] dropped it, and
//! the tool moved the read right by the one reference base it removed. `4M2I4M` over an `AC`
//! repeat becomes `2I8M` at the same start: an insertion that walks to the front is kept, and only
//! a deletion is dropped, which is why the tool has a line for one and not the other.
//!
//! This is the second tool of this archetype that changes a read's position, after
//! [`crate::print_distant_mates`], and it does it for an unrelated reason.
//!
//! # Two classes of read never reach the call
//!
//! An unmapped read, and a read whose cigar has one element or none, are written straight out. The
//! guard is `read.isUnmapped()`, which is the adapter's three criteria rather than the flag alone,
//! and it runs before the reference is touched.

use htsjdk_bam::record::BamRecord;

use gatk_engine::alignment_utils::{self, AlignmentError};
use gatk_engine::read;
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_engine::reference::ReferenceFileSource;

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK LeftAlignIndels";

/// What a run produces: the output BAM and its index, or the refusal the util raised.
pub type RunResult = Result<Result<(Vec<u8>, Option<Vec<u8>>), AlignmentError>, ReadsError>;

/// `apply` for one read: the left-aligned cigar, and the position it moved to.
///
/// `reference_bases` is the window the walker handed this read, which is the read's own span.
/// Returns `false` when the read was passed through untouched, which is the first branch of the
/// reference's `apply` rather than a no-op outcome of the call.
pub fn left_align(read: &mut BamRecord, reference_bases: &[u8]) -> Result<bool, AlignmentError> {
    // We cannot deal with screwy records, and a read with one cigar element is a trivial case.
    if read::is_unmapped(read) || read.cigar.num_elements() <= 1 {
        return Ok(false);
    }

    let result =
        alignment_utils::left_align_indels(&read.cigar, reference_bases, &read.read_bases, 0)?;
    read.cigar = result.cigar;
    if result.leading_deletion_bases_removed > 0 {
        // `setPosition(read.getContig(), read.getStart() + removed)`: the contig does not change,
        // so only the start moves.
        read.alignment_start += result.leading_deletion_bases_removed as i32;
    }
    Ok(true)
}

/// `LeftAlignIndels`: every read the traversal reaches, with its indels moved as far left as the
/// read's own window allows.
pub fn left_align_indels(
    source: &ReadsDataSource,
    reference: &mut ReferenceFileSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> RunResult {
    let applied = crate::read_walker::traverse_with_reference(
        source,
        Some(reference),
        &options.intervals,
        false,
        filter,
    )?;

    let mut records = Vec::with_capacity(applied.len());
    for mut entry in applied {
        let bases = entry
            .context
            .bases(reference)
            .map_err(|error| ReadsError::Malformed(format!("{error:?}")))?;
        if let Err(error) = left_align(&mut entry.read, &bases) {
            return Ok(Err(error));
        }
        records.push(entry.read);
    }

    let header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    Ok(Ok(write_records(
        &header,
        &records,
        options.create_output_bam_index,
    )?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::text_parse::parse_cigar;

    fn read(start: i32, cigar: &str, bases: &str) -> BamRecord {
        BamRecord {
            read_name: "r".to_string(),
            flags: 0,
            reference_index: 0,
            alignment_start: start,
            mapping_quality: 60,
            cigar: parse_cigar(cigar).unwrap(),
            read_bases: bases.as_bytes().to_vec(),
            ..BamRecord::default()
        }
    }

    #[test]
    fn a_deletion_that_walks_off_the_front_moves_the_read() {
        let mut record = read(6, "4M1D5M", "AAAATTTTT");
        assert!(left_align(&mut record, b"AAAAATTTTT").unwrap());
        assert_eq!(record.cigar.to_text(), "9M");
        assert_eq!(record.alignment_start, 7);
    }

    #[test]
    fn an_insertion_that_walks_off_the_front_does_not() {
        let mut record = read(21, "4M2I4M", "ACACACACAC");
        assert!(left_align(&mut record, b"ACACACAC").unwrap());
        assert_eq!(record.cigar.to_text(), "2I8M");
        assert_eq!(record.alignment_start, 21, "only a deletion moves the read");
    }

    #[test]
    fn an_indel_with_no_repeat_to_move_into_stays() {
        let mut record = read(17, "4M1D5M", "TTTTCACAC");
        assert!(left_align(&mut record, b"TTTTACACAC").unwrap());
        assert_eq!(record.cigar.to_text(), "4M1D5M");
        assert_eq!(record.alignment_start, 17);
    }

    #[test]
    fn a_single_element_cigar_never_reaches_the_call() {
        let mut record = read(11, "10M", "TTTTTTTTTT");
        assert!(
            !left_align(&mut record, b"TTTTTTTTTT").unwrap(),
            "passed through"
        );
        assert_eq!(record.cigar.to_text(), "10M");
    }

    #[test]
    fn an_unmapped_read_never_reaches_the_call() {
        let mut record = read(41, "4M1D5M", "AAAATTTTT");
        record.flags |= 0x4;
        assert!(!left_align(&mut record, b"").unwrap(), "passed through");
        assert_eq!(record.cigar.to_text(), "4M1D5M", "and keeps its cigar");
    }
}
