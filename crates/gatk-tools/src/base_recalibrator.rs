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
/// `contig_bases` answers a whole reference contig by name, because the counting pass needs the
/// read's span for its comparison and a wider window for BAQ. Each read is compared against ITS
/// contig: taking the reads' first one compared a chr2 read against chr1's bases.
///
/// The quantization is computed **after** `finalizeData`, so it depends on every read the traversal
/// kept, and the report is written from the finalised tables.
/// `contig_bases` is a CLOSURE and not a slice, because the order matters: the reference is opened
/// at startup and queried only once the traversal has produced reads, so a file whose reads cannot
/// be read fails on the reads and not on the reference. A row that hands a BAM the wrong index and
/// a reference the reads' contig is missing from answers the read failure, which is what the
/// reference answers.
#[allow(clippy::too_many_arguments)]
pub fn base_recalibrator(
    source: &ReadsDataSource,
    contig_bases: &mut dyn FnMut(&str) -> Result<Vec<u8>, BaseRecalibratorError>,
    known_sites: &[SimpleInterval],
    arguments: &EngineArguments,
    quantizing_levels: i32,
    filter: &dyn Fn(&BamRecord) -> bool,
    intervals: &[SimpleInterval],
) -> Result<String, BaseRecalibratorError> {
    let header = source.header().clone();
    // The traversal is bounded like any other walker's. Hardcoding no intervals here counted the
    // reads a command line had excluded: `--intervals chr1:1-6000 --exclude-intervals chr1:1-500`
    // leaves seven of the corpus's eight reads, and the eighth's nine countable bases were nine
    // observations the recalibration table carried too many.
    let records = crate::read_walker::traverse(source, intervals, filter)
        .map_err(BaseRecalibratorError::Reads)?;

    // The bases come after the traversal, for the reason in the doc comment above, and once per
    // contig.
    let mut fetched: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();

    let mut engine = BaseRecalibrationEngine::new(arguments.clone(), &header)
        .map_err(BaseRecalibratorError::Engine)?;
    for record in &records {
        let contig = usize::try_from(record.reference_index)
            .ok()
            .and_then(|index| header.sequences.get(index))
            .map(|sequence| sequence.name.clone());
        let bases: &[u8] = match contig {
            Some(contig) => {
                if !fetched.contains_key(&contig) {
                    let bases = contig_bases(&contig)?;
                    fetched.insert(contig.clone(), bases);
                }
                &fetched[&contig]
            }
            None => &[],
        };
        let overlapping = sites_overlapping(record, &header, known_sites);
        engine
            .process_read(record, &header, bases, &overlapping)
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

/// The known sites the engine is handed for one read: those on its contig overlapping
/// `new SimpleInterval(read)`, its aligned span, which is the walker's `FeatureContext` and
/// `JoinReadsWithVariants`' query alike. A site in a soft clip alone is not one of them, and a read
/// whose span is no valid interval is handed none.
fn sites_overlapping(
    record: &BamRecord,
    header: &htsjdk_bam::header::SamHeader,
    known_sites: &[SimpleInterval],
) -> Vec<SimpleInterval> {
    let Some(contig) = usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
    else {
        return Vec::new();
    };
    let start = record.alignment_start;
    let end = record.alignment_end();
    if start < 1 || end < start {
        return Vec::new();
    }
    known_sites
        .iter()
        .filter(|site| site.contig == contig.name && site.start <= end && site.end >= start)
        .cloned()
        .collect()
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
