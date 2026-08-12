//! The writing side of `RecalUtils`, ported from
//! `org.broadinstitute.hellbender.utils.recalibration` (GATK 4.6.2.0).
//!
//! [`crate::recalibration_report`] reads a recalibration table; this writes one. The two are the
//! same five tables in the same order, and `BaseRecalibrator` writing a file that `ApplyBQSR` reads
//! back is the cycle both halves have to agree on.
//!
//! # The two additional covariate tables share one report table
//!
//! ```java
//! } else if (allCovsReportTable == null && recalibrationTables.isAdditionalCovariateTable(table)) {
//!     reportTable = makeNewTableWithColumns(columnNames, reportTableName, sort);
//!     rowIndex = 0;
//!     allCovsReportTable = reportTable;
//!     addToList = true;
//! } else {
//!     reportTable = allCovsReportTable;
//!     addToList = false;
//! }
//! ```
//!
//! The context and cycle tables both write into one `RecalTable2`, and `rowIndex` is **not** reset
//! between them. The reference's own comment calls that "knowledge about the ordering of tables".
//! The rows are then sorted by their values, so the two covariates come out interleaved.
//!
//! # And the row keys never reach the file
//!
//! Every table here sorts `SORT_BY_COLUMN`, which orders by the row's **data** and not by its key,
//! so the integer row index the writer counts with is invisible in the output. It matters only
//! because two rows written under the same index would overwrite each other.

use crate::covariates::{CovariateKind, RecalibrationArguments, StandardCovariateList};
use crate::gatk_report::{Report, Sorting, Table, Value};
use crate::qual_quantizer::QuantizationInfo;
use crate::recal_datum::EventType;
use crate::recalibration_report::{
    ALL_COVARIATES_REPORT_TABLE_TITLE, ARGUMENT_REPORT_TABLE_TITLE,
    QUALITY_SCORE_REPORT_TABLE_TITLE, QUANTIZED_REPORT_TABLE_TITLE, READGROUP_REPORT_TABLE_TITLE,
};
use crate::recalibration_tables::RecalibrationTables;

/// `RecalibrationArgumentCollection.generateReportTable(covariateNames)`.
///
/// Seventeen rows, in the order the reference adds them, sorted by the first column when written.
/// **The value column has an empty format**, so it is `Unknown` and every value is written through
/// `String.valueOf`, which is why an integer argument has no padding rule of its own.
pub fn generate_argument_table(
    arguments: &RecalibrationArguments,
    quantizing_levels: i32,
    covariate_names: &str,
) -> Table {
    let mut table = Table::new(
        ARGUMENT_REPORT_TABLE_TITLE,
        "Recalibration argument collection values used in this run",
        Sorting::SortByColumn,
    );
    table.add_column("Argument", "%s");
    // The empty format, which becomes `%s` with the `Unknown` data type.
    table.add_column("Value", "");

    let mut row = |name: &str, value: Value| {
        // `addRowID(name, true)` populates the first column with the row's own id.
        table.set(name, "Argument", Value::Str(name.to_string()));
        table.set(name, "Value", value);
    };

    row("covariate", Value::Str(covariate_names.to_string()));
    row("no_standard_covs", Value::Bool(false));
    row("run_without_dbsnp", Value::Bool(false));
    row("solid_recal_mode", Value::Str("SET_Q_ZERO".to_string()));
    row(
        "solid_nocall_strategy",
        Value::Str("THROW_EXCEPTION".to_string()),
    );
    row(
        "mismatches_context_size",
        Value::Int(arguments.mismatches_context_size as i64),
    );
    row(
        "indels_context_size",
        Value::Int(arguments.indels_context_size as i64),
    );
    row("mismatches_default_quality", Value::Int(-1));
    row("deletions_default_quality", Value::Int(45));
    row("insertions_default_quality", Value::Int(45));
    row(
        "maximum_cycle_value",
        Value::Int(arguments.maximum_cycle_value as i64),
    );
    row(
        "low_quality_tail",
        Value::Int(arguments.low_qual_tail as i64),
    );
    row("default_platform", Value::Str("null".to_string()));
    row("force_platform", Value::Str("null".to_string()));
    row("quantizing_levels", Value::Int(quantizing_levels as i64));
    row("recalibration_report", Value::Str("null".to_string()));
    row("binary_tag_name", Value::Str("null".to_string()));
    table
}

/// `QuantizationInfo.generateReportTable()`: one row per quality score, all 94 of them.
pub fn generate_quantization_table(info: &QuantizationInfo) -> Table {
    let mut table = Table::new(
        QUANTIZED_REPORT_TABLE_TITLE,
        "Quality quantization map",
        Sorting::SortByColumn,
    );
    table.add_column("QualityScore", "%d");
    table.add_column("Count", "%d");
    table.add_column("QuantizedScore", "%d");
    for qual in 0..info.quantized_quals.len() {
        let key = qual.to_string();
        table.set(&key, "QualityScore", Value::Int(qual as i64));
        table.set(
            &key,
            "Count",
            Value::Int(info.empirical_qual_counts.get(qual).copied().unwrap_or(0)),
        );
        table.set(
            &key,
            "QuantizedScore",
            Value::Int(info.quantized_quals[qual] as i64),
        );
    }
    table
}

/// `RecalUtils.generateReportTables(recalibrationTables, covariates)`.
///
/// Three tables out of four: the read group table, the quality score table, and **one** table for
/// both additional covariates. See the module note.
pub fn generate_report_tables(
    tables: &RecalibrationTables,
    covariates: &StandardCovariateList,
) -> Vec<Table> {
    let mut result = Vec::new();
    let mut all_covariates: Option<Table> = None;
    let mut row_index = 0usize;

    for (index, table) in tables.all_tables.iter().enumerate() {
        let is_read_group = index == 0;
        let is_additional = index >= 2;

        // The columns, in the order the reference adds them, which is the order they are written.
        let mut columns: Vec<(&str, &str)> = vec![("ReadGroup", "%s")];
        if !is_read_group {
            columns.push(("QualityScore", "%d"));
            if is_additional {
                columns.push(("CovariateValue", "%s"));
                columns.push(("CovariateName", "%s"));
            }
        }
        columns.push(("EventType", "%s"));
        columns.push(("EmpiricalQuality", "%.4f"));
        if is_read_group {
            // The read group table alone carries the estimated reported quality.
            columns.push(("EstimatedQReported", "%.4f"));
        }
        columns.push(("Observations", "%d"));
        columns.push(("Errors", "%.2f"));

        let name = if is_read_group {
            READGROUP_REPORT_TABLE_TITLE
        } else if index == 1 {
            QUALITY_SCORE_REPORT_TABLE_TITLE
        } else {
            ALL_COVARIATES_REPORT_TABLE_TITLE
        };

        // The first additional covariate makes the shared table; the second writes into it.
        let mut report = if is_additional && all_covariates.is_some() {
            all_covariates.take().unwrap()
        } else {
            row_index = 0;
            let mut fresh = Table::new(name, "", Sorting::SortByColumn);
            for (column, format) in &columns {
                fresh.add_column(column, format);
            }
            fresh
        };

        let kind = covariates.kinds()[index];
        for (keys, datum) in table.all_leaves() {
            let key = row_index.to_string();
            let mut key_index = 0;
            report.set(
                &key,
                "ReadGroup",
                Value::Str(
                    covariates
                        .read_group
                        .format_key(keys[key_index])
                        .unwrap_or("null")
                        .to_string(),
                ),
            );
            key_index += 1;
            if !is_read_group {
                report.set(&key, "QualityScore", Value::Int(keys[key_index] as i64));
                key_index += 1;
                if is_additional {
                    let value = match kind {
                        CovariateKind::Context => covariates
                            .context
                            .format_key(keys[key_index])
                            .ok()
                            .flatten()
                            .unwrap_or_else(|| "null".to_string()),
                        CovariateKind::Cycle => covariates.cycle.format_key(keys[key_index]),
                        _ => keys[key_index].to_string(),
                    };
                    report.set(&key, "CovariateValue", Value::Str(value));
                    report.set(
                        &key,
                        "CovariateName",
                        Value::Str(kind.parsed_name().to_string()),
                    );
                    key_index += 1;
                }
            }
            let event =
                EventType::from_index(keys[key_index]).unwrap_or(EventType::BaseSubstitution);
            report.set(
                &key,
                "EventType",
                Value::Str(event.representation().to_string()),
            );
            // This call computes and caches the empirical quality, which is why writing a table
            // changes the datums in it.
            report.set(
                &key,
                "EmpiricalQuality",
                Value::Double(datum.borrow_mut().empirical_quality()),
            );
            if is_read_group {
                report.set(
                    &key,
                    "EstimatedQReported",
                    Value::Double(datum.borrow().reported_quality()),
                );
            }
            report.set(
                &key,
                "Observations",
                Value::Int(datum.borrow().num_observations()),
            );
            report.set(
                &key,
                "Errors",
                Value::Double(datum.borrow().num_mismatches()),
            );
            row_index += 1;
        }

        if is_additional {
            all_covariates = Some(report);
        } else {
            result.push(report);
        }
    }
    if let Some(shared) = all_covariates {
        result.push(shared);
    }
    result
}

/// `RecalUtils.outputRecalibrationReport`: the whole five-table report as text.
pub fn output_recalibration_report(
    arguments: &RecalibrationArguments,
    quantizing_levels: i32,
    quantization_info: &QuantizationInfo,
    tables: &RecalibrationTables,
    covariates: &StandardCovariateList,
) -> String {
    let mut report = Report::new();
    report.add_table(generate_argument_table(
        arguments,
        quantizing_levels,
        &covariates.covariate_names(),
    ));
    report.add_table(generate_quantization_table(quantization_info));
    for table in generate_report_tables(tables, covariates) {
        report.add_table(table);
    }
    report.write()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_argument_table_has_seventeen_rows_and_an_untyped_value_column() {
        let table = generate_argument_table(
            &RecalibrationArguments::default(),
            16,
            "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate",
        );
        assert_eq!(table.rows.len(), 17);
        // The empty format becomes `%s`, and its data type stays `Unknown`.
        assert_eq!(table.columns[1].format, "%s");
        assert_eq!(
            table.columns[1].data_type,
            crate::gatk_report::DataType::Unknown
        );
    }

    #[test]
    fn the_additional_covariates_share_one_table() {
        let covariates = StandardCovariateList::new(
            &RecalibrationArguments::default(),
            &["unit-rg1".to_string()],
        )
        .unwrap();
        let tables = RecalibrationTables::new(&covariates).unwrap();
        let report = generate_report_tables(&tables, &covariates);
        // Four recalibration tables, three report tables.
        assert_eq!(tables.all_tables.len(), 4);
        assert_eq!(report.len(), 3);
        assert_eq!(report[2].name, ALL_COVARIATES_REPORT_TABLE_TITLE);
        assert_eq!(report[2].columns.len(), 8);
    }
}
