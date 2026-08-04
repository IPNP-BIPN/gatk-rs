//! Ported from `org.broadinstitute.hellbender.tools.PrintReads` and the header handling in
//! `GATKTool.getHeaderForSAMWriter` (GATK 4.6.2.0).
//!
//! The first whole tool here, and the first whose claim is bytes. `PrintReads` is a `ReadWalker`
//! whose `apply` writes the read out, so the interesting part is not the traversal, which is
//! [`crate::read_walker`], but what happens to the header on the way.
//!
//! Two things there decide the output's bytes:
//!
//!  * a `@PG` record is appended, with `ID` = the **tool name**, `VN` = the GATK version, `CL` =
//!    the whole expanded command line and `PN` = the tool name again. The tool name is
//!    `GATK PrintReads`, with a space, while the command line begins `PrintReads`: they are not
//!    the same string, which is measurable and easy to get wrong;
//!  * the ID is made unique by appending `.1`, `.2` and so on until it is free, so running the
//!    tool over its own output produces `GATK PrintReads.1` rather than replacing the first
//!    record or emitting a duplicate.
//!
//! `getHeaderForSAMWriter` mutates the reads header in place rather than copying it, which is not
//! observable in one run's output and is reproduced anyway.

use gatk_engine::reads::{ReadsDataSource, ReadsError};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

pub use crate::sam_output::Options;

/// `GATKTool.getToolName()` for this tool, which is not the string the command line starts with.
pub const TOOL_NAME: &str = "GATK PrintReads";

/// `GATKTool.createProgramGroupID` for this tool, kept as a named re-export because the suite and
/// the sibling tools reach it by this path.
pub fn create_program_group_id(header: &SamHeader, tool_name: &str) -> String {
    crate::sam_output::create_program_group_id(header, tool_name)
}

/// `GATKTool.getHeaderForSAMWriter` with this tool's name.
pub fn header_for_sam_writer(header: &SamHeader, options: &Options) -> SamHeader {
    crate::sam_output::header_for_sam_writer(header, TOOL_NAME, options)
}

/// `PrintReads`: the reads that survive the traversal, written back out.
///
/// Returns the BAM's bytes and, when asked for, its `.bai`.
pub fn print_reads(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ReadsError> {
    let records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    let header = header_for_sam_writer(source.header(), options);
    crate::sam_output::write_records(&header, &records, options.create_output_bam_index)
}
