//! `GatherBQSRReports`, ported from `RecalibrationReport.gatherReports` (GATK 4.6.2.0).
//!
//! The gather that ends a scattered BQSR run: several recalibration reports summed into one.
//!
//! # It is not a concatenation
//!
//! ```java
//! final RecalibrationReport result = inputs.stream()
//!         .map(i -> new RecalibrationReport(new GATKReport(i), allReadGroups))
//!         .reduce(RecalibrationReport::combine)
//!         .filter(r -> !r.isEmpty())
//!         .orElseThrow(() -> new GATKException("there is no usable data in any input file"));
//! result.quantizationInfo = new QuantizationInfo(result.recalibrationTables, result.RAC.QUANTIZING_LEVELS);
//! ```
//!
//! Four things happen there that a copy would not do. The READ GROUPS ARE THE UNION of every
//! input's, computed before any report is parsed, so every shard numbers its keys the same way. The
//! COMBINE SUMS the observations and the errors and leaves the empirical qualities to be
//! recomputed, which the writer then does per row. The QUANTIZATION IS RECOMPUTED from the summed
//! quality-score table, at the FIRST report's level count. And an EMPTY REPORT IS SKIPPED by the
//! combine, so a gather of nothing but empty shards falls out of the `filter` and refuses.
//!
//! The arguments table written at the end is the first input's. What it holds is the recalibration
//! ARGUMENT COLLECTION and not the command line, so every shard of one scattered run carries the
//! same one and gathering two shards in either order gives the same bytes.

use gatk_engine::qual_quantizer::{QualQuantizer, QuantizationInfo, MIN_USABLE_Q_SCORE};
use gatk_engine::recal_utils::output_recalibration_report;
use gatk_engine::recalibration_report::{RecalibrationReport, RecalibrationReportError};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK GatherBQSRReports";

/// What the gather refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum GatherError {
    /// `Utils.nonEmpty(inputs, ...)`.
    NoInputs,
    /// Every input was empty, so the reduce's `filter` left nothing.
    NoUsableData,
    /// One of the reports could not be read.
    Report(RecalibrationReportError),
    /// Two shards' tables could not be summed, which needs tables of different shapes.
    Combine(gatk_engine::recalibration_tables::CombineError),
}

impl GatherError {
    pub fn java_class(&self) -> &'static str {
        match self {
            GatherError::NoInputs => "java.lang.IllegalArgumentException",
            GatherError::NoUsableData => "org.broadinstitute.hellbender.exceptions.GATKException",
            // Every reader refusal this port has met is one of the report reader's own, whose
            // class its suite pins; the gather adds nothing to it.
            GatherError::Report(_) | GatherError::Combine(_) => {
                "org.broadinstitute.hellbender.exceptions.GATKException"
            }
        }
    }

    pub fn message(&self) -> String {
        match self {
            GatherError::NoInputs => "Cannot gather an empty list of inputs".to_string(),
            GatherError::NoUsableData => "there is no usable data in any input file".to_string(),
            GatherError::Report(error) => error.message(),
            GatherError::Combine(error) => error.message(),
        }
    }
}

/// `new QuantizationInfo(tables, levels)`, on the summed tables.
fn quantization_info(
    tables: &gatk_engine::recalibration_tables::RecalibrationTables,
    levels: i32,
) -> QuantizationInfo {
    let mut histogram = vec![0i64; 94];
    for (_, datum) in tables.quality_score_table().all_leaves() {
        let empirical =
            gatk_engine::math_utils::fast_round(datum.borrow_mut().empirical_quality()) as usize;
        if empirical < histogram.len() {
            histogram[empirical] += datum.borrow().num_observations();
        }
    }
    match QualQuantizer::new(&histogram, levels, MIN_USABLE_Q_SCORE) {
        Ok(quantizer) => QuantizationInfo::new(quantizer.original_to_quantized_map, histogram),
        Err(_) => QuantizationInfo::new(vec![0; 94], histogram),
    }
}

/// `gatherReports`, as the text of the gathered report.
pub fn gather(reports: &[&str]) -> Result<String, GatherError> {
    if reports.is_empty() {
        return Err(GatherError::NoInputs);
    }
    // The union, sorted, taken from every input before any of them is parsed properly.
    let mut all_read_groups: Vec<String> = Vec::new();
    for text in reports {
        // `new GATKReport(input).getReadGroups()`, which is the `RecalTable0` column and not a
        // parsed report: the reference reads the read groups of every input before it builds one.
        let report = gatk_engine::gatk_report::Report::parse(text)
            .map_err(|error| GatherError::Report(RecalibrationReportError::Report(error)))?;
        let groups = report
            .read_groups()
            .map_err(|error| GatherError::Report(RecalibrationReportError::Report(error)))?;
        for group in groups {
            if !all_read_groups.contains(&group) {
                all_read_groups.push(group);
            }
        }
    }
    all_read_groups.sort();

    let mut combined: Option<RecalibrationReport> = None;
    for text in reports {
        let report = RecalibrationReport::parse_with_read_groups(text, Some(&all_read_groups))
            .map_err(GatherError::Report)?;
        match &mut combined {
            // `combine` returns `this` untouched when the OTHER report is empty.
            Some(total) => {
                if !report.is_empty() {
                    total
                        .tables
                        .combine(&report.tables)
                        .map_err(GatherError::Combine)?;
                }
            }
            None => combined = Some(report),
        }
    }
    let result = combined.expect("at least one input");
    if result.is_empty() {
        return Err(GatherError::NoUsableData);
    }
    let quantization = quantization_info(&result.tables, result.quantizing_levels);
    Ok(output_recalibration_report(
        &result.arguments,
        result.quantizing_levels,
        &quantization,
        &result.tables,
        &result.covariates,
    ))
}
