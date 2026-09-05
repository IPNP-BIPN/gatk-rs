//! Ported from `org.broadinstitute.hellbender.tools.walkers.UnmarkDuplicates` (GATK 4.6.2.0).
//!
//! The second whole tool here, and the point of it is the measurement rather than the tool: G2's
//! calibration gate asks what the *second* member of the largest archetype costs once the first has
//! paid for the engine. The answer is this file, and it is short because
//! [`crate::sam_output`] holds everything `PrintReads` already paid for.
//!
//! ```java
//! public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
//!     read.setIsDuplicate(false);
//!     outputWriter.addRead(read);
//! }
//! ```
//!
//! # It is not `PrintReads` with one line changed
//!
//! One override makes the difference, and it is not in `apply`:
//!
//! ```java
//! public List<ReadFilter> getDefaultReadFilters() {
//!     return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
//! }
//! ```
//!
//! `PrintReads` takes `GATKTool`'s default, which is `WellformedReadFilter` — six predicates that
//! drop reads with an inconsistent length, a missing sequence, an invalid alignment start and so
//! on. This tool replaces the whole list with one filter that keeps everything. So on a file
//! containing a malformed read the two tools emit **different sets of reads**, and it is the
//! filter rather than the transform that decides it.
//!
//! That is the kind of difference an archetype hides. Two tools in the same archetype, both
//! `ReadWalker`s whose `apply` is two lines, and their default traversal is not the same.

use gatk_engine::reads::{ReadsDataSource, ReadsError};
use htsjdk_bam::record::BamRecord;

use crate::sam_output::{header_for_sam_writer, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK UnmarkDuplicates";

/// `SAMFlag.DUPLICATE_READ`, the 0x400 bit this tool exists to clear.
pub const DUPLICATE_READ_FLAG: u16 = 0x400;

/// `read.setIsDuplicate(false)`.
///
/// Unconditional: a read that was not flagged is written back unchanged rather than skipped, which
/// matters because the traversal still counts it and the writer still re-encodes it.
pub fn unmark(read: &mut BamRecord) {
    read.flags &= !DUPLICATE_READ_FLAG;
}

/// `UnmarkDuplicates`: every read, with 0x400 cleared, written back out.
///
/// The filter is a parameter rather than a default because the tool's own default is
/// `ALLOW_ALL_READS`; passing anything else models `--read-filter` on the command line.
pub fn unmark_duplicates(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ReadsError> {
    unmark_duplicates_with(
        source,
        options,
        filter,
        htsjdk_bgzf::DEFAULT_COMPRESSION_LEVEL,
        htsjdk_bgzf::Deflater::Jdk,
    )
}

/// The same, with the BGZF compression named: see [`crate::sam_output::write_records_with`].
///
/// A command line needs this rather than the default: `GATKConfig` sets `samjdk.compression_level`
/// to TWO and routes it through GKL, so a BAM written at the library default is a different file.
pub fn unmark_duplicates_with(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ReadsError> {
    let mut records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    for record in &mut records {
        unmark(record);
    }
    let header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    crate::sam_output::write_records_with(
        &header,
        &records,
        options.create_output_bam_index,
        level,
        deflater,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_the_duplicate_bit_is_touched() {
        let mut read = BamRecord {
            // Paired, first of pair, reverse strand, duplicate, and supplementary: everything
            // except 0x400 must survive.
            flags: 0x1 | 0x40 | 0x10 | DUPLICATE_READ_FLAG | 0x800,
            ..BamRecord::default()
        };
        unmark(&mut read);
        assert_eq!(read.flags, 0x1 | 0x40 | 0x10 | 0x800);
    }

    #[test]
    fn a_read_that_was_never_flagged_is_unchanged() {
        let mut read = BamRecord {
            flags: 0x1 | 0x40,
            ..BamRecord::default()
        };
        unmark(&mut read);
        assert_eq!(read.flags, 0x1 | 0x40);
    }
}
