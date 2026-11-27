//! Ported from `org.broadinstitute.hellbender.tools.FixMisencodedBaseQualityReads` and
//! `MisencodedBaseQualityReadTransformer` (GATK 4.7.0.0; `FixMisencodedBaseQualityReads.java` and
//! `MisencodedBaseQualityReadTransformer.java` are unchanged since 4.6.2.0).
//!
//! The fourth whole tool of the record-transform archetype, and the first whose transform can
//! refuse the file it was given.
//!
//! ```java
//!         final byte[] quals = read.getBaseQualities();
//!         for ( int i = 0; i < quals.length; i++ ) {
//!             quals[i] -= ILLUMINA_ENCODING_FIX_VALUE;
//!             if ( quals[i] < 0 )
//! ```
//!
//! and the `if` throws `UserException.BadInput` with the message [`MisencodedError::message`]
//! returns; after the loop, `read.setBaseQualities(quals);` and `return read;`.
//!
//! # It keeps the default read filters, where the other two replace them
//!
//! `UnmarkDuplicates` and `RevertBaseQualityScores` both override `getDefaultReadFilters` with
//! `ALLOW_ALL_READS`. This one does not, so it takes `GATKTool`'s default of
//! `WellformedReadFilter`. Three tools in one archetype and two different default traversals; this
//! is the one that keeps the default, which is why naming `AllowAllReadsReadFilter` changes what it
//! emits and would not change what the other two emit.
//!
//! # The refusal is a property of the traversal, not of the file
//!
//! A quality below 31 aborts the whole run. An interval that excludes the read carrying it does
//! not: the same file succeeds or fails depending on what the traversal reaches, which the corpus
//! carries as two cases over one fixture.
//!
//! # A read with no qualities passes through
//!
//! The loop does not run, the empty array is set back, and nothing refuses. `*` in a SAM is not a
//! quality of zero, and this is the difference showing.
//!
//! # Reading this file without reading Rust
//!
//! - `u8` is an unsigned byte, 0 to 255, where Java's `byte` is signed, -128 to 127. That is the
//!   one place this file's arithmetic differs in kind from the Java's: see [`fix`].
//! - `Result<T, E>` is "a value or an error". [`RunResult`] nests two: the outer error is a read
//!   failure, the inner one the tool's own refusal, so a caller can tell them apart.
//! - `?` after a call propagates its error to the caller, Java's `throws` at one call site.
//! - `&mut BamRecord` is a borrowed, writable reference: [`fix`] changes the caller's record in
//!   place, as the Java changes the read's array.

use gatk_engine::reads::{ReadsDataSource, ReadsError};
use htsjdk_bam::record::BamRecord;

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK FixMisencodedBaseQualityReads";

/// `ILLUMINA_ENCODING_FIX_VALUE`: Illumina's 64 less Phred's 33.
pub const ILLUMINA_ENCODING_FIX_VALUE: u8 = 31;

/// Why the transform refused the file: the `UserException.BadInput` that
/// `MisencodedBaseQualityReadTransformer.apply` throws. It carries nothing, because the Java
/// message names neither the read nor the base.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MisencodedError;

impl MisencodedError {
    /// The reference's message, which says the read "was correctly encoded": the tool has been
    /// asked to fix an encoding the read does not have.
    pub fn message(&self) -> String {
        "while fixing mis-encoded base qualities we encountered a read that was correctly encoded; \
         we cannot handle such a mixture of reads so unfortunately the BAM must be fixed with some \
         other tool"
            .to_string()
    }

    /// The Java class, which is the nested one rather than `UserException` itself.
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    }
}

/// The transform: subtract 31 from every quality, refusing the first one that would go below zero.
///
/// The check is **per base and after the subtraction**, so the bases before the offending one are
/// already rewritten when it throws. That is only observable through the reads written before the
/// one that failed, which is why the corpus puts the offending base in the middle of a read in the
/// middle of a file.
///
/// A quality of exactly 31 becomes 0, which is not below zero: the boundary is inclusive on the
/// side that succeeds.
///
/// Java's `byte` is signed and this `u8` is not, so the two agree only for stored qualities up to
/// 158. A stored byte of 159 to 255 is -97 to -1 in Java, goes below zero by the subtraction and
/// is refused; here it is at least 31 and is rewritten. The SAM text format cannot carry a
/// quality above 93, so only a BAM can reach that range.
pub fn fix(read: &mut BamRecord) -> Result<(), MisencodedError> {
    // `iter_mut()` yields a writable reference to each quality in turn; `*quality` is the value
    // it points at.
    for quality in read.base_qualities.iter_mut() {
        // Java subtracts on a signed byte and compares against zero. A quality below 31 wraps
        // there and is negative; here it would wrap the other way, so the comparison is made
        // before the subtraction rather than after it.
        if *quality < ILLUMINA_ENCODING_FIX_VALUE {
            return Err(MisencodedError);
        }
        *quality -= ILLUMINA_ENCODING_FIX_VALUE;
    }
    Ok(())
}

/// What a run produces: the output BAM and its index, or the refusal, or a failure to read.
pub type RunResult = Result<Result<(Vec<u8>, Option<Vec<u8>>), MisencodedError>, ReadsError>;

/// `FixMisencodedBaseQualityReads`: every read the traversal reaches, with 31 taken off every
/// quality.
///
/// The filter is a parameter because this tool's default is `GATKTool`'s, not `ALLOW_ALL_READS`:
/// the caller passes whichever one the command line resolved to. `&dyn Fn(&BamRecord) -> bool` is
/// any function or closure from a record to a boolean, Java's `Predicate<GATKRead>`.
///
/// The Java transforms and writes one read at a time,
/// `outputWriter.addRead(transform.apply(read));`; this traverses first and transforms the whole
/// list, and on a refusal writes nothing at all.
pub fn fix_misencoded_base_quality_reads(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> RunResult {
    let mut records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    // `if let Err(error) = ...` runs the block only when `fix` refused, binding the error.
    for record in &mut records {
        if let Err(error) = fix(record) {
            return Ok(Err(error));
        }
    }
    fixed_records_to_bam(source, options, records)
}

/// The same, with the BGZF compression named: see [`crate::sam_output::write_records_with`].
///
/// A real `gatk` run writes at level TWO through GKL rather than at htsjdk's default of five.
pub fn fix_misencoded_base_quality_reads_with(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
) -> RunResult {
    let mut records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    for record in &mut records {
        if let Err(error) = fix(record) {
            return Ok(Err(error));
        }
    }
    let header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    Ok(Ok(crate::sam_output::write_records_with(
        &header,
        &records,
        options.create_output_bam_index,
        level,
        deflater,
    )?))
}

/// Writes the fixed records under the header `GATKTool` would write, at htsjdk's default
/// compression.
fn fixed_records_to_bam(
    source: &ReadsDataSource,
    options: &Options,
    records: Vec<BamRecord>,
) -> RunResult {
    let header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    Ok(Ok(write_records(
        &header,
        &records,
        options.create_output_bam_index,
    )?))
}
