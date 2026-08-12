//! Ported from `org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator` (GATK 4.6.2.0).
//!
//! The tool that writes the table [`crate::apply_bqsr`] reads, and the one that closes the BQSR
//! cycle. Its body is three calls: count every read that survives the filters, finalise the tables,
//! and write the report.
//!
//! # Its default read filters are seven, not one
//!
//! ```java
//! public static List<ReadFilter> getBQSRSpecificReadFilterList() {
//!     filters.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
//!     filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
//!     filters.add(ReadFilterLibrary.MAPPED);
//!     filters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
//!     filters.add(ReadFilterLibrary.NOT_DUPLICATE);
//!     filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
//! }
//! ```
//!
//! plus `WellformedReadFilter`. That is a fifth pattern of `getDefaultReadFilters` across the ported
//! tools, after taking the engine's default, replacing it with `ALLOW_ALL_READS`, extending it with
//! four more, and replacing it with a single non-Wellformed filter. It is also what decides which
//! reads are counted at all, so it is part of the table's contents and not a detail of the
//! traversal.
//!
//! # And the known sites must be indexed
//!
//! `--known-sites` is required and takes any feature file, VCF or BED. Without an index the
//! reference refuses before reading a single read, because the traversal queries it by interval.

use gatk_engine::base_recalibration_engine::{
    BaseRecalibrationEngine, EngineArguments, EngineError,
};
use gatk_engine::interval::SimpleInterval;
use gatk_engine::qual_quantizer::{QualQuantizer, QuantizationInfo, MIN_USABLE_Q_SCORE};
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_engine::recal_utils::output_recalibration_report;
use htsjdk_bam::record::BamRecord;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK BaseRecalibrator";

/// `RecalibrationArgumentCollection.QUANTIZING_LEVELS`.
pub const QUANTIZING_LEVELS: i32 = 16;

/// What the tool stops with.
#[derive(Debug, Clone, PartialEq)]
pub enum BaseRecalibratorError {
    Engine(EngineError),
    Reads(ReadsError),
}

impl BaseRecalibratorError {
    pub fn message(&self) -> String {
        match self {
            BaseRecalibratorError::Engine(error) => error.message(),
            BaseRecalibratorError::Reads(error) => format!("{error:?}"),
        }
    }
}

/// `BaseRecalibrator`: the recalibration table this run produced, as text.
///
/// `contig_bases` is the whole reference contig, because the counting pass needs the read's span for
/// its comparison and a wider window for BAQ.
///
/// The quantization is computed **after** `finalizeData`, so it depends on every read the traversal
/// kept, and the report is written from the finalised tables.
pub fn base_recalibrator(
    source: &ReadsDataSource,
    contig_bases: &[u8],
    known_sites: &[SimpleInterval],
    arguments: &EngineArguments,
    quantizing_levels: i32,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<String, BaseRecalibratorError> {
    let header = source.header().clone();
    let records =
        crate::read_walker::traverse(source, &[], filter).map_err(BaseRecalibratorError::Reads)?;

    let mut engine = BaseRecalibrationEngine::new(arguments.clone(), &header)
        .map_err(BaseRecalibratorError::Engine)?;
    for record in &records {
        engine
            .process_read(record, &header, contig_bases, known_sites)
            .map_err(BaseRecalibratorError::Engine)?;
    }
    engine
        .finalize_data()
        .map_err(BaseRecalibratorError::Engine)?;

    let quantization = quantization_info(&engine, quantizing_levels);
    Ok(output_recalibration_report(
        &arguments.covariates,
        quantizing_levels,
        &quantization,
        &engine.tables,
        &engine.covariates,
    ))
}

/// `new QuantizationInfo(finalTables, QUANTIZING_LEVELS)`: the empirical quality histogram of the
/// quality score table, quantized.
///
/// It is built from the **finalised** tables, so the empirical qualities it bins are the ones the
/// rounding left behind.
fn quantization_info(engine: &BaseRecalibrationEngine, levels: i32) -> QuantizationInfo {
    let mut histogram = vec![0i64; 94];
    for (_, datum) in engine.tables.quality_score_table().all_leaves() {
        let empirical =
            gatk_engine::math_utils::fast_round(datum.borrow_mut().empirical_quality()) as usize;
        if empirical < histogram.len() {
            histogram[empirical] += datum.borrow().num_observations();
        }
    }
    match QualQuantizer::new(&histogram, levels, MIN_USABLE_Q_SCORE) {
        Ok(quantizer) => QuantizationInfo::new(quantizer.original_to_quantized_map, histogram),
        // `nLevels = 0` is the reference's own null dereference; the tool never passes it, because
        // its argument's default is sixteen.
        Err(_) => QuantizationInfo::new(vec![0; 94], histogram),
    }
}
