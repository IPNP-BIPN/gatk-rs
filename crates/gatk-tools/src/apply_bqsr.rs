//! Ported from `org.broadinstitute.hellbender.tools.walkers.bqsr.ApplyBQSR` (GATK 4.6.2.0).
//!
//! The eighth whole tool of the record-transform archetype, and the one every brick of
//! [`gatk_engine::recal_datum`], [`gatk_engine::covariates`],
//! [`gatk_engine::recalibration_tables`], [`gatk_engine::qual_quantizer`],
//! [`gatk_engine::recalibration_report`] and [`gatk_engine::bqsr_transformer`] was built for.
//!
//! The tool's own body is four lines:
//!
//! ```java
//! public ReadTransformer makePostReadFilterTransformer(){
//!     return new BQSRReadTransformer(getHeaderForReads(), bqsrRecalFile, bqsrArgs);
//! }
//! public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
//!     outputWriter.addRead(read);
//! }
//! ```
//!
//! # The transformer runs after the read filters
//!
//! `makePostReadFilterTransformer`, not `makePreReadFilterTransformer`. A read the filters drop is
//! therefore **never recalibrated and never written**: it does not pass through unrecalibrated, it
//! disappears. This tool takes `GATKTool`'s default filter list, `WellformedReadFilter`, so a read
//! whose cigar disagrees with its length is gone from the output, and the golden carries a run with
//! that filter disabled to show the same file producing two reads instead of one.
//!
//! # And the recalibration file decides the read group keys, not the BAM
//!
//! The covariates come from the report's own `RecalTable0`. A BAM whose read group the report does
//! not name is refused, and `--allow-missing-read-group` makes those reads **quantized but not
//! recalibrated** rather than passed through untouched.

use gatk_engine::bqsr_transformer::{ApplyBqsrArguments, BqsrError, BqsrReadTransformer};
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_engine::recalibration_report::{RecalibrationReport, RecalibrationReportError};
use htsjdk_bam::record::BamRecord;

use crate::sam_output::{header_for_sam_writer, Options};

/// `GATKTool.getToolName()` for this tool, which is not the string the command line starts with.
pub const TOOL_NAME: &str = "GATK ApplyBQSR";

/// What the tool stops with.
#[derive(Debug, Clone, PartialEq)]
pub enum ApplyBqsrError {
    /// The recalibration file could not be read.
    Report(RecalibrationReportError),
    /// A read the transformer refused, which is a read group the table has nothing for.
    Transform(BqsrError),
    /// Anything reading or writing the BAM refused.
    Reads(ReadsError),
}

impl ApplyBqsrError {
    pub fn message(&self) -> String {
        match self {
            ApplyBqsrError::Report(error) => error.message(),
            ApplyBqsrError::Transform(error) => error.message(),
            ApplyBqsrError::Reads(error) => format!("{error:?}"),
        }
    }

    /// The Java class the reference throws, which the non-user handler prints.
    ///
    /// One class for all three: the report's own refusals and the transformer's are
    /// `UserException`s, and a read the reader cannot make sense of is one too.
    pub fn java_class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }
}

/// `ApplyBQSR`: every read that survives the filters, recalibrated and written back out.
///
/// Returns the BAM's bytes and, when asked for, its `.bai`. The writer is created presorted, as the
/// reference's `createSAMWriter(output, true)` is, so the traversal's order is the output's.
///
/// `recal_text` is the recalibration table's own text, because the report is read whole before the
/// traversal starts.
pub fn apply_bqsr(
    source: &ReadsDataSource,
    recal_text: &str,
    arguments: &ApplyBqsrArguments,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ApplyBqsrError> {
    apply_bqsr_with(
        source,
        recal_text,
        arguments,
        options,
        filter,
        htsjdk_bgzf::DEFAULT_COMPRESSION_LEVEL,
        htsjdk_bgzf::Deflater::Jdk,
    )
}

/// The same, with the BGZF compression named: see [`crate::sam_output::write_records_with`].
#[allow(clippy::too_many_arguments)]
pub fn apply_bqsr_with(
    source: &ReadsDataSource,
    recal_text: &str,
    arguments: &ApplyBqsrArguments,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ApplyBqsrError> {
    let mut report = RecalibrationReport::parse(recal_text).map_err(ApplyBqsrError::Report)?;

    // The filters run first and the transformer second, which is the whole point of
    // `makePostReadFilterTransformer`. A read dropped here never reaches the recalibration.
    let records = crate::read_walker::traverse(source, &options.intervals, filter)
        .map_err(ApplyBqsrError::Reads)?;

    let header = source.header().clone();
    let mut transformer = BqsrReadTransformer::new(
        &header,
        &mut report.tables,
        &mut report.quantization_info,
        &report.covariates,
        arguments,
    )
    .map_err(ApplyBqsrError::Transform)?;

    let mut recalibrated = Vec::with_capacity(records.len());
    for record in &records {
        recalibrated.push(
            transformer
                .apply(record)
                .map_err(ApplyBqsrError::Transform)?,
        );
    }

    let out_header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    crate::sam_output::write_records_with(
        &out_header,
        &recalibrated,
        options.create_output_bam_index,
        level,
        deflater,
    )
    .map_err(ApplyBqsrError::Reads)
}
