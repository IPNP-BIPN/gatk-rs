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
    /// NOT one class for all three. The transformer's refusals are `GATKException`s: a
    /// recalibration table written from a different BAM answers `Read group rg1 not found in the
    /// recalibration table.` and leaves **exit 3**, where a `UserException` would leave 2 and be
    /// printed as `A USER ERROR has occurred`. Measured on row 10 of this tool's covering array,
    /// which is a row the corpus only reached once `--bqsr-recal-file` had a second value.
    ///
    /// The report's own refusals and the reader's are user errors, and stay one.
    pub fn java_class(&self) -> &'static str {
        match self {
            ApplyBqsrError::Transform(error) => error.java_class(),
            ApplyBqsrError::Report(_) | ApplyBqsrError::Reads(_) => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
        }
    }

    /// Whether `Main` prints it as a user error, which is what decides the exit code.
    pub fn is_user(&self) -> bool {
        self.java_class() == "org.broadinstitute.hellbender.exceptions.UserException"
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
    let report = RecalibrationReport::parse(recal_text).map_err(ApplyBqsrError::Report)?;

    // The filters run first and the transformer second, which is the whole point of
    // `makePostReadFilterTransformer`. A read dropped here never reaches the recalibration.
    let records = crate::read_walker::traverse(source, &options.intervals, filter)
        .map_err(ApplyBqsrError::Reads)?;

    let recalibrated = recalibrate_with(report, source.header(), arguments, &records)?;

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

/// [`apply_bqsr_with`] straight to `output`, a read at a time: the reads are recalibrated as the
/// traversal hands them over and written as they are recalibrated, so neither the reads nor the
/// output are held. A refusal, the traversal's or the transformer's, comes at the read that causes
/// it, as in the reference's own loop, and leaves `output` as it was: see
/// [`crate::sam_output::write_bam_file`]. The inner error is the file's.
#[allow(clippy::too_many_arguments)]
pub fn apply_bqsr_to(
    output: &std::path::Path,
    source: &ReadsDataSource,
    recal_text: &str,
    arguments: &ApplyBqsrArguments,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
) -> Result<Result<crate::sam_output::Written, std::io::Error>, ApplyBqsrError> {
    // The report and the transformer are refused here, before a read is read, as the reference's
    // constructor refuses them; the workers below build their own, which then cannot fail.
    let mut report = RecalibrationReport::parse(recal_text).map_err(ApplyBqsrError::Report)?;
    let transformer = BqsrReadTransformer::new(
        source.header(),
        &mut report.tables,
        &mut report.quantization_info,
        &report.covariates,
        arguments,
    )
    .map_err(ApplyBqsrError::Transform)?;
    let records = crate::read_walker::traverse_iter(source, &options.intervals, false, filter)
        .map_err(ApplyBqsrError::Reads)?;
    let out_header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    let write = |recalibrated: &mut dyn Iterator<Item = Result<BamRecord, ReadsError>>| {
        crate::sam_output::write_bam_file(
            output,
            &out_header,
            recalibrated,
            options.create_output_bam_index,
            level,
            deflater,
        )
    };
    // The writer takes reads errors only; a transformer refusal is kept here and stops it.
    let mut refused = None;
    let workers = std::thread::available_parallelism()
        .map(|threads| threads.get().saturating_sub(1).min(MAX_WORKERS))
        .unwrap_or(0);
    let written = if workers < 2 {
        let mut transformer = transformer;
        let mut recalibrated = records.map(|record| {
            let record = record?;
            transformer.apply(&record).map_err(|error| {
                refused = Some(error);
                ReadsError::Malformed(String::new())
            })
        });
        write(&mut recalibrated)
    } else {
        drop(transformer);
        recalibrate_in_parallel(
            records,
            workers,
            source.header(),
            recal_text,
            arguments,
            &mut refused,
            write,
        )
    };
    match (written, refused) {
        (_, Some(error)) => Err(ApplyBqsrError::Transform(error)),
        (written, None) => written.map_err(ApplyBqsrError::Reads),
    }
}

/// At most this many recalibrating threads: past it the traversal and the writer are the limit.
const MAX_WORKERS: usize = 6;
/// Reads per batch handed to a worker.
const BATCH: usize = 2048;

/// The recalibration of `records` on `workers` threads, each with its own report and transformer,
/// handed to `write` in the traversal's order. Each datum computes its empirical quality once,
/// from a prior that depends only on its own keys, so a thread's copy answers what the one
/// transformer would. Refusals come in the traversal's order too: a batch's transformer refusal
/// before a later batch's, and a read the traversal refuses after the batches before it.
fn recalibrate_in_parallel<R>(
    mut records: impl Iterator<Item = Result<BamRecord, ReadsError>>,
    workers: usize,
    header: &htsjdk_bam::header::SamHeader,
    recal_text: &str,
    arguments: &ApplyBqsrArguments,
    refused: &mut Option<BqsrError>,
    write: impl FnOnce(&mut dyn Iterator<Item = Result<BamRecord, ReadsError>>) -> R,
) -> R {
    use std::sync::mpsc;
    std::thread::scope(|scope| {
        let lanes: Vec<_> = (0..workers)
            .map(|_| {
                let (to_worker, batches) = mpsc::sync_channel::<Vec<BamRecord>>(1);
                let (to_writer, results) = mpsc::sync_channel(1);
                scope.spawn(move || {
                    let Ok(mut report) = RecalibrationReport::parse(recal_text) else {
                        return;
                    };
                    let Ok(mut transformer) = BqsrReadTransformer::new(
                        header,
                        &mut report.tables,
                        &mut report.quantization_info,
                        &report.covariates,
                        arguments,
                    ) else {
                        return;
                    };
                    for batch in batches {
                        let recalibrated: Result<Vec<BamRecord>, BqsrError> =
                            batch.iter().map(|read| transformer.apply(read)).collect();
                        let failed = recalibrated.is_err();
                        if to_writer.send(recalibrated).is_err() || failed {
                            return;
                        }
                    }
                });
                (to_worker, results)
            })
            .collect();
        let mut in_flight = std::collections::VecDeque::new();
        let mut next_lane = 0;
        let mut traversal_refusal = None;
        let mut exhausted = false;
        let mut current = Vec::new().into_iter();
        let mut stopped = false;
        let mut recalibrated = std::iter::from_fn(|| loop {
            if stopped {
                return None;
            }
            if let Some(read) = current.next() {
                return Some(Ok(read));
            }
            while !exhausted && in_flight.len() < workers {
                let batch: Vec<BamRecord> = records
                    .by_ref()
                    .map_while(|record| match record {
                        Ok(record) => Some(record),
                        Err(error) => {
                            traversal_refusal = Some(error);
                            None
                        }
                    })
                    .take(BATCH)
                    .collect();
                exhausted = traversal_refusal.is_some() || batch.len() < BATCH;
                if batch.is_empty() {
                    break;
                }
                if lanes[next_lane].0.send(batch).is_err() {
                    stopped = true;
                    return Some(Err(ReadsError::Malformed(
                        "a recalibrating thread stopped".to_string(),
                    )));
                }
                in_flight.push_back(next_lane);
                next_lane = (next_lane + 1) % workers;
            }
            let Some(lane) = in_flight.pop_front() else {
                stopped = true;
                return traversal_refusal.take().map(Err);
            };
            match lanes[lane].1.recv() {
                Ok(Ok(batch)) => current = batch.into_iter(),
                Ok(Err(error)) => {
                    *refused = Some(error);
                    stopped = true;
                    return Some(Err(ReadsError::Malformed(String::new())));
                }
                Err(_) => {
                    stopped = true;
                    return Some(Err(ReadsError::Malformed(
                        "a recalibrating thread stopped".to_string(),
                    )));
                }
            }
        });
        let written = write(&mut recalibrated);
        // Closing the lanes ends the workers' loops, so the scope can join them.
        drop(lanes);
        written
    })
}

/// `BQSRReadTransformer` over reads already filtered, which `ApplyBQSR` shares with the two Spark
/// tools that apply a table. `ApplyBQSRSparkFn` builds one transformer per partition, and one
/// partition is one transformer.
pub fn recalibrate(
    header: &htsjdk_bam::header::SamHeader,
    recal_text: &str,
    arguments: &ApplyBqsrArguments,
    records: &[BamRecord],
) -> Result<Vec<BamRecord>, ApplyBqsrError> {
    let report = RecalibrationReport::parse(recal_text).map_err(ApplyBqsrError::Report)?;
    recalibrate_with(report, header, arguments, records)
}

fn recalibrate_with(
    mut report: RecalibrationReport,
    header: &htsjdk_bam::header::SamHeader,
    arguments: &ApplyBqsrArguments,
    records: &[BamRecord],
) -> Result<Vec<BamRecord>, ApplyBqsrError> {
    let mut transformer = BqsrReadTransformer::new(
        header,
        &mut report.tables,
        &mut report.quantization_info,
        &report.covariates,
        arguments,
    )
    .map_err(ApplyBqsrError::Transform)?;

    let mut recalibrated = Vec::with_capacity(records.len());
    for record in records {
        recalibrated.push(
            transformer
                .apply(record)
                .map_err(ApplyBqsrError::Transform)?,
        );
    }
    Ok(recalibrated)
}
