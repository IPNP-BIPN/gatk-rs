//! `SetNmMdAndUqTags.doWork`, over [`gatk_tools::set_nm_md_and_uq_tags`].
//!
//! The header's sort order is checked before the writer exists, every record is fixed against the
//! contig the walker hands over BY INDEX, and the records are written presorted under the input's
//! header. `SetNmAndUqTags` is a deprecated subclass that overrides nothing, so it runs this too.

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::set_nm_md_and_uq_tags::{check_sort_order, fix_record, Arguments};

use super::gather_bam_files::{read_reads, write_reads};
use super::{assert_readable, assert_writable, thrown};
use crate::runners::{argument, flag};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    let reference = argument(parser, "REFERENCE_SEQUENCE").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let (header, mut records) = read_reads(&input)?;
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
        if record.flags & 0x4 != 0 {
            continue;
        }
        // `ReferenceSequenceFileWalker.get(index)`, which walks the FASTA by position.
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
