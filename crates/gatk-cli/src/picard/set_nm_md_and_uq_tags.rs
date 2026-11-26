//! `SetNmMdAndUqTags.doWork`, over [`gatk_tools::set_nm_md_and_uq_tags`].
//!
//! The header's sort order is checked before the writer exists, every record is fixed against the
//! contig the walker hands over BY INDEX, and the records are written presorted under the input's
//! header. `SetNmAndUqTags` is a deprecated subclass that overrides nothing, so it runs this too.
//!
//! # Reading this file without reading Rust
//!
//! - `?` returns the error to the caller at once, Java's `throws`;
//! - `for record in &mut records` walks the list handing out a writable borrow of each record,
//!   so `fix_record` edits it in place, as Java's `peek(rec -> fixRecord(rec, refSeqWalker))`;
//! - the chain on `reference_index` is an `Option`, a value or nothing, passed along: `.ok()`
//!   turns a failed conversion (a negative index) into nothing, `.and_then` looks the contig up
//!   only when there is an index, and `.ok_or_else` turns nothing into the exception.

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::set_nm_md_and_uq_tags::{check_sort_order, fix_record, Arguments};

use super::gather_bam_files::{read_reads, write_reads};
use super::{assert_readable, assert_writable, thrown};
use crate::runners::{argument, flag};

/// `SetNmMdAndUqTags.doWork`, in the reference's order: INPUT and OUTPUT asserted, the reader
/// opened, the sort order refused before any writer exists, then every record fixed and written.
pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    let reference = argument(parser, "REFERENCE_SEQUENCE").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let (header, mut records) = read_reads(&input)?;
    // `if (reader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate)`, with
    // a missing `SO` read as `unsorted`, the header's default.
    check_sort_order(header.attributes.get("SO"))
        .map_err(|error| thrown(error.java_class(), error.message()))?;

    let contigs = htsjdk_bam::fasta::read_fasta_file(&reference).map_err(|error| {
        thrown(
            "htsjdk.samtools.SAMException",
            format!("{reference}: {error:?}"),
        )
    })?;
    let arguments = Arguments {
        is_bisulfite_sequence: flag(parser, "IS_BISULFITE_SEQUENCE"),
        set_only_uq: flag(parser, "SET_ONLY_UQ"),
    };
    for record in &mut records {
        // `if (!record.getReadUnmappedFlag())`: 0x4 is the unmapped flag, and an unmapped read
        // never asks the walker for a contig.
        if record.flags & 0x4 != 0 {
            continue;
        }
        // `ReferenceSequenceFileWalker.get(index)`, which walks the FASTA by position. The port
        // looks the contig up for every mapped read; the reference's `fixUq` only does so
        // `if (record.getBaseQualities() != SAMRecord.NULL_QUALS)`, which matters under
        // `SET_ONLY_UQ` alone.
        let bases = usize::try_from(record.reference_index)
            .ok()
            .and_then(|index| contigs.get(index))
            .map(|contig| contig.bases.as_slice())
            .ok_or_else(|| {
                thrown(
                    "htsjdk.samtools.SAMException",
                    format!("No reference sequence for index {}", record.reference_index),
                )
            })?;
        fix_record(record, bases, &arguments);
    }
    write_reads(
        parser,
        &output,
        &header,
        &records,
        flag(parser, "CREATE_INDEX"),
        flag(parser, "CREATE_MD5_FILE"),
    )?;
    Ok(0)
}
