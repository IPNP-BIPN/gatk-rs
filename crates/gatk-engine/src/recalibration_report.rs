//! `RecalibrationReport`, ported from `org.broadinstitute.hellbender.utils.recalibration`
//! (GATK 4.6.2.0).
//!
//! A recalibration table read off disk: five named tables in a `GATKReport`, turned back into a
//! [`RecalibrationTables`], a [`StandardCovariateList`] and a [`QuantizationInfo`]. Every part is
//! measured on its own already; this is the assembly, and the assembly is where the keys are chosen.
//!
//! | table | what it holds |
//! |---|---|
//! | `Arguments` | the covariate sizes the tables were built with |
//! | `Quantized` | the map from every quality to its quantized value |
//! | `RecalTable0` | one datum per (read group, event) |
//! | `RecalTable1` | one per (read group, reported quality, event) |
//! | `RecalTable2` | one per (read group, reported quality, covariate value, event) |
//!
//! # The read group keys are the sorted order, not the file's
//!
//! `getReadGroups` returns a `TreeSet` and the covariate list is built from it, so a report written
//! with `zebra` first comes back with `alpha` as key 0. That is what makes a gathered report number
//! its groups the same way whichever input was read first, and it means **the numbering in the file
//! is not the numbering in memory**.
//!
//! It also reads **only `RecalTable0`**. A read group named in `RecalTable1` but not in
//! `RecalTable0` gets the missing key -1 and is written at a negative index, which the reference
//! reports as `ArrayIndexOutOfBoundsException: Index -1 out of bounds for length 1`.
//!
//! # And the empirical quality column is thrown away
//!
//! Every datum is built with a reported quality of one and then corrected through the setter, and
//! its empirical quality is left uncomputed on purpose, so `ApplyBQSR` can recompute it against
//! whatever prior it is run with. A report whose `EmpiricalQuality` column is nonsense parses
//! without complaint.

use std::cell::RefCell;
use std::rc::Rc;

use crate::covariates::{CovariateKind, RecalibrationArguments, StandardCovariateList};
use crate::gatk_report::{Report, ReportReadError, Table, Value};
use crate::qual_quantizer::QuantizationInfo;
use crate::recal_datum::{EventType, RecalDatum, RecalDatumError};
use crate::recalibration_tables::{NestedArrayError, RecalibrationTables, SharedDatum};

/// `RecalUtils`' table titles.
pub const ARGUMENT_REPORT_TABLE_TITLE: &str = "Arguments";
pub const QUANTIZED_REPORT_TABLE_TITLE: &str = "Quantized";
pub const READGROUP_REPORT_TABLE_TITLE: &str = "RecalTable0";
pub const QUALITY_SCORE_REPORT_TABLE_TITLE: &str = "RecalTable1";
pub const ALL_COVARIATES_REPORT_TABLE_TITLE: &str = "RecalTable2";

/// `RecalUtils`' column names.
const ARGUMENT_COLUMN_NAME: &str = "Argument";
const ARGUMENT_VALUE_COLUMN_NAME: &str = "Value";
const QUANTIZED_VALUE_COLUMN_NAME: &str = "QuantizedScore";
const QUANTIZED_COUNT_COLUMN_NAME: &str = "Count";
const READGROUP_COLUMN_NAME: &str = "ReadGroup";
const EVENT_TYPE_COLUMN_NAME: &str = "EventType";
const ESTIMATED_Q_REPORTED_COLUMN_NAME: &str = "EstimatedQReported";
const QUALITY_SCORE_COLUMN_NAME: &str = "QualityScore";
const COVARIATE_VALUE_COLUMN_NAME: &str = "CovariateValue";
const COVARIATE_NAME_COLUMN_NAME: &str = "CovariateName";
const NUMBER_OBSERVATIONS_COLUMN_NAME: &str = "Observations";
const NUMBER_ERRORS_COLUMN_NAME: &str = "Errors";

/// Everything reading a recalibration report can refuse.
#[derive(Debug, Clone, PartialEq)]
pub enum RecalibrationReportError {
    /// Anything the underlying `GATKReport` reader refuses.
    Report(ReportReadError),
    /// The reference's `ArrayIndexOutOfBoundsException`, from a read group the covariate does not
    /// know reaching a table as the key -1. See the module note.
    Nested(NestedArrayError),
    /// A datum the report describes that `RecalDatum` refuses.
    Datum(RecalDatumError),
    /// `UserException`: the `covariate` argument names something other than the four standard ones.
    NonStandardCovariates { supported: String, found: String },
    /// `IllegalArgumentException` from `EventType.eventFrom(String)`.
    UnknownEventType(String),
    /// The reference's `NullPointerException`: a covariate name `RecalTable2` carries that the
    /// standard list does not hold.
    UnknownCovariateName(String),
    /// `GATKException` from `asLong` or `asDouble` on a cell of a type they do not accept.
    UnexpectedCellType { column: String },
}

impl RecalibrationReportError {
    pub fn message(&self) -> String {
        match self {
            RecalibrationReportError::Report(error) => error.message(),
            RecalibrationReportError::Nested(error) => error.message(),
            RecalibrationReportError::Datum(error) => error.message(),
            RecalibrationReportError::NonStandardCovariates { supported, found } => format!(
                "Non-standard covariates are not supported. Only the following are supported \
                 {supported} but was {found}"
            ),
            RecalibrationReportError::UnknownEventType(name) => {
                format!("Event {name} does not exist.")
            }
            RecalibrationReportError::UnknownCovariateName(name) => {
                format!("no covariate named {name}")
            }
            RecalibrationReportError::UnexpectedCellType { column } => {
                format!("Object in column {column} is not of the expected type")
            }
        }
    }
}

/// `RecalibrationReport`.
#[derive(Debug)]
pub struct RecalibrationReport {
    pub arguments: RecalibrationArguments,
    /// `QUANTIZING_LEVELS`, kept because re-quantizing needs it.
    pub quantizing_levels: i32,
    pub quantization_info: QuantizationInfo,
    pub covariates: StandardCovariateList,
    pub tables: RecalibrationTables,
}

impl RecalibrationReport {
    /// `new RecalibrationReport(InputStream)`.
    ///
    /// The order matters and is the reference's. `getReadGroups` runs **first**, because the public
    /// constructor is `this(report, report.getReadGroups())` and Java evaluates the argument before
    /// the body, so a report with no `RecalTable0` is refused by that name rather than by the
    /// argument table. Then the arguments, which size the covariates; then the quantization map;
    /// then the three data tables, keyed by the read groups already fixed.
    pub fn parse(text: &str) -> Result<RecalibrationReport, RecalibrationReportError> {
        RecalibrationReport::parse_with_read_groups(text, None)
    }

    /// `new RecalibrationReport(report, allReadGroups)`, the constructor the gather uses.
    ///
    /// `read_groups` is the UNION of every input's read groups, so that a shard which never saw one
    /// still numbers its keys the same way as the shard that did. `None` is the public
    /// constructor, which takes the report's own.
    pub fn parse_with_read_groups(
        text: &str,
        read_groups: Option<&[String]>,
    ) -> Result<RecalibrationReport, RecalibrationReportError> {
        let report = Report::parse(text).map_err(RecalibrationReportError::Report)?;

        // `RecalTable0` is asked for FIRST, before the arguments, because the public constructor is
        // `this(report, report.getReadGroups())` and Java evaluates the argument before the body. A
        // report with no read group table is refused by name here and not by the argument table.
        let own_read_groups = report
            .read_groups()
            .map_err(RecalibrationReportError::Report)?;
        let read_groups: Vec<String> = read_groups
            .map(<[String]>::to_vec)
            .unwrap_or(own_read_groups);

        let argument_table = report
            .table_named(ARGUMENT_REPORT_TABLE_TITLE)
            .map_err(RecalibrationReportError::Report)?;
        let (arguments, quantizing_levels) = parse_arguments(argument_table)?;

        let quantized_table = report
            .table_named(QUANTIZED_REPORT_TABLE_TITLE)
            .map_err(RecalibrationReportError::Report)?;
        let quantization_info = parse_quantization(quantized_table)?;

        let covariates = StandardCovariateList::new(&arguments, &read_groups).map_err(|error| {
            RecalibrationReportError::Report(ReportReadError::NotANumber(error.message()))
        })?;
        let mut tables = RecalibrationTables::with_read_groups(&covariates, read_groups.len())
            .map_err(RecalibrationReportError::Nested)?;

        parse_read_group_table(
            report
                .table_named(READGROUP_REPORT_TABLE_TITLE)
                .map_err(RecalibrationReportError::Report)?,
            &covariates,
            &mut tables,
        )?;
        parse_quality_score_table(
            report
                .table_named(QUALITY_SCORE_REPORT_TABLE_TITLE)
                .map_err(RecalibrationReportError::Report)?,
            &covariates,
            &mut tables,
        )?;
        parse_all_covariates_table(
            report
                .table_named(ALL_COVARIATES_REPORT_TABLE_TITLE)
                .map_err(RecalibrationReportError::Report)?,
            &covariates,
            &mut tables,
        )?;

        Ok(RecalibrationReport {
            arguments,
            quantizing_levels,
            quantization_info,
            covariates,
            tables,
        })
    }

    /// `isEmpty()`: no table holds a datum.
    pub fn is_empty(&self) -> bool {
        self.tables.is_empty()
    }
}

/// The cell at one row of a named column, or the reference's missing-table error.
fn cell<'a>(table: &'a Table, row: usize, column: &str) -> Option<&'a Value> {
    let index = table.columns.iter().position(|c| c.name == column)?;
    table.rows.get(row).and_then(|values| values.get(index))
}

/// `asLong(Object)`: a `Long`, an `Integer` widened, or a `Double` **truncated**.
fn as_long(value: Option<&Value>, column: &str) -> Result<i64, RecalibrationReportError> {
    match value {
        Some(Value::Int(number)) => Ok(*number),
        // `((Double) o).longValue()`, which truncates toward zero.
        Some(Value::Double(number)) => Ok(*number as i64),
        _ => Err(RecalibrationReportError::UnexpectedCellType {
            column: column.to_string(),
        }),
    }
}

/// `decodeByte(Object)`: a `Byte`, a `String` parsed, or a `Long` narrowed.
///
/// The `QualityScore` column is written `%s`, so a parsed report always takes the `String` branch.
fn decode_byte(table: &Table, row: usize, column: &str) -> Result<i8, RecalibrationReportError> {
    match cell(table, row, column) {
        Some(Value::Str(text)) => {
            text.parse::<i8>()
                .map_err(|_| RecalibrationReportError::UnexpectedCellType {
                    column: column.to_string(),
                })
        }
        // `((Long) value).byteValue()`, which truncates to the low eight bits.
        Some(Value::Int(number)) => Ok(*number as i8),
        _ => Err(RecalibrationReportError::UnexpectedCellType {
            column: column.to_string(),
        }),
    }
}

/// `asDouble(Object)`: a `Double`, or an `Integer` or `Long` widened.
fn as_double(value: Option<&Value>, column: &str) -> Result<f64, RecalibrationReportError> {
    match value {
        Some(Value::Double(number)) => Ok(*number),
        Some(Value::Int(number)) => Ok(*number as f64),
        _ => Err(RecalibrationReportError::UnexpectedCellType {
            column: column.to_string(),
        }),
    }
}

/// `getRecalDatum(table, row, hasEstimatedQReportedColumn)`.
///
/// The datum is built with a reported quality of **one** and then corrected, because the constructor
/// refuses a negative quality and the setter is what the report's own value goes through. The
/// empirical quality is deliberately left uncomputed.
fn recal_datum(
    table: &Table,
    row: usize,
    has_estimated_q_reported: bool,
) -> Result<SharedDatum, RecalibrationReportError> {
    let observations = as_long(
        cell(table, row, NUMBER_OBSERVATIONS_COLUMN_NAME),
        NUMBER_OBSERVATIONS_COLUMN_NAME,
    )?;
    let errors = as_double(
        cell(table, row, NUMBER_ERRORS_COLUMN_NAME),
        NUMBER_ERRORS_COLUMN_NAME,
    )?;
    let reported = if has_estimated_q_reported {
        // The read group table alone carries this column, and it is a Double.
        as_double(
            cell(table, row, ESTIMATED_Q_REPORTED_COLUMN_NAME),
            ESTIMATED_Q_REPORTED_COLUMN_NAME,
        )?
    } else {
        // `decodeByte(get(row, QUALITY_SCORE_COLUMN_NAME))`, and **not** `asDouble`. The writer
        // declares that column `%s`, so a parsed report hands back a String and the decoder parses
        // it as a byte; a `%d` column would have handed back a Long and taken the other branch.
        decode_byte(table, row, QUALITY_SCORE_COLUMN_NAME)? as f64
    };

    let mut datum =
        RecalDatum::new(observations, errors, 1).map_err(RecalibrationReportError::Datum)?;
    datum
        .set_reported_quality(reported)
        .map_err(RecalibrationReportError::Datum)?;
    Ok(Rc::new(RefCell::new(datum)))
}

/// `QualityScoreCovariate.keyFromValue` on a row's `QualityScore` cell.
fn quality_key(table: &Table, row: usize) -> Result<i32, RecalibrationReportError> {
    decode_byte(table, row, QUALITY_SCORE_COLUMN_NAME).map(|quality| quality as i32)
}

/// The `EventType` a row names, by its one-letter representation.
fn event_of(table: &Table, row: usize) -> Result<EventType, RecalibrationReportError> {
    let name = match cell(table, row, EVENT_TYPE_COLUMN_NAME) {
        Some(Value::Str(text)) => text.clone(),
        other => format!("{other:?}"),
    };
    EventType::from_representation(&name).ok_or(RecalibrationReportError::UnknownEventType(name))
}

/// The text of a cell, which is what every key lookup is given.
fn text_of(table: &Table, row: usize, column: &str) -> String {
    match cell(table, row, column) {
        Some(Value::Str(text)) => text.clone(),
        Some(Value::Int(number)) => number.to_string(),
        Some(Value::Double(number)) => number.to_string(),
        Some(Value::Bool(flag)) => flag.to_string(),
        Some(Value::Char(character)) => character.to_string(),
        _ => "null".to_string(),
    }
}

fn parse_read_group_table(
    table: &Table,
    covariates: &StandardCovariateList,
    tables: &mut RecalibrationTables,
) -> Result<(), RecalibrationReportError> {
    for row in 0..table.rows.len() {
        let group =
            covariates
                .read_group
                .key_from_value(&text_of(table, row, READGROUP_COLUMN_NAME));
        let event = event_of(table, row)?.ordinal() as i32;
        let datum = recal_datum(table, row, true)?;
        tables
            .read_group_table_mut()
            .put(datum, &[group, event])
            .map_err(RecalibrationReportError::Nested)?;
    }
    Ok(())
}

fn parse_quality_score_table(
    table: &Table,
    covariates: &StandardCovariateList,
    tables: &mut RecalibrationTables,
) -> Result<(), RecalibrationReportError> {
    for row in 0..table.rows.len() {
        let group =
            covariates
                .read_group
                .key_from_value(&text_of(table, row, READGROUP_COLUMN_NAME));
        // `QualityScoreCovariate.keyFromValue`, whose String branch is `Byte.parseByte`.
        let quality = quality_key(table, row)?;
        let event = event_of(table, row)?.ordinal() as i32;
        let datum = recal_datum(table, row, false)?;
        tables
            .quality_score_table_mut()
            .put(datum, &[group, quality, event])
            .map_err(RecalibrationReportError::Nested)?;
    }
    Ok(())
}

fn parse_all_covariates_table(
    table: &Table,
    covariates: &StandardCovariateList,
    tables: &mut RecalibrationTables,
) -> Result<(), RecalibrationReportError> {
    for row in 0..table.rows.len() {
        let group =
            covariates
                .read_group
                .key_from_value(&text_of(table, row, READGROUP_COLUMN_NAME));
        let quality = quality_key(table, row)?;
        let name = text_of(table, row, COVARIATE_NAME_COLUMN_NAME);
        let kind = covariates
            .covariate_by_parsed_name(&name)
            .ok_or_else(|| RecalibrationReportError::UnknownCovariateName(name.clone()))?;
        let value = text_of(table, row, COVARIATE_VALUE_COLUMN_NAME);
        let key = match kind {
            CovariateKind::Context => covariates.context.key_from_value(&value),
            CovariateKind::Cycle => covariates
                .cycle
                .key_from_value(value.parse::<i32>().unwrap_or(0))
                .map_err(|error| RecalibrationReportError::UnknownCovariateName(error.message()))?,
            // The special covariates never appear in this table.
            CovariateKind::ReadGroup => covariates.read_group.key_from_value(&value),
            CovariateKind::QualityScore => value.parse::<i32>().unwrap_or(0),
        };
        let event = event_of(table, row)?.ordinal() as i32;
        let datum = recal_datum(table, row, false)?;
        let index = covariates.index_by_class(kind) as usize;
        tables.all_tables[index]
            .put(datum, &[group, quality, key, event])
            .map_err(RecalibrationReportError::Nested)?;
    }
    Ok(())
}

/// `initializeQuantizationTable`: read **by row index**, not by the `QualityScore` column.
///
/// A table whose rows are out of order is therefore read as though they were in order. The values
/// go through `toString` and then a parse, so a `%d` column that came back as a `Long` is turned
/// into a byte the long way round.
fn parse_quantization(table: &Table) -> Result<QuantizationInfo, RecalibrationReportError> {
    let mut quals = vec![0u8; 94];
    let mut counts = vec![0i64; 94];
    for row in 0..table.rows.len() {
        if row >= quals.len() {
            break;
        }
        quals[row] = text_of(table, row, QUANTIZED_VALUE_COLUMN_NAME)
            .parse::<i16>()
            .map_err(|_| RecalibrationReportError::UnexpectedCellType {
                column: QUANTIZED_VALUE_COLUMN_NAME.to_string(),
            })? as u8;
        counts[row] = text_of(table, row, QUANTIZED_COUNT_COLUMN_NAME)
            .parse::<i64>()
            .map_err(|_| RecalibrationReportError::UnexpectedCellType {
                column: QUANTIZED_COUNT_COLUMN_NAME.to_string(),
            })?;
    }
    Ok(QuantizationInfo::new(quals, counts))
}

/// `initializeArgumentCollectionTable`: a chain of name comparisons over the `Arguments` table.
///
/// `null` is turned into a real absence first, so `binary_tag_name` of `null` is nothing rather than
/// the four characters. The `covariate` argument is checked against the four standard class names
/// and anything else is refused.
fn parse_arguments(
    table: &Table,
) -> Result<(RecalibrationArguments, i32), RecalibrationReportError> {
    let mut arguments = RecalibrationArguments::default();
    let mut quantizing_levels = 16;
    let standard = "[ReadGroupCovariate, QualityScoreCovariate, ContextCovariate, CycleCovariate]";

    for row in 0..table.rows.len() {
        let name = text_of(table, row, ARGUMENT_COLUMN_NAME);
        let raw = text_of(table, row, ARGUMENT_VALUE_COLUMN_NAME);
        // "if (value.equals("null")) value = null;", before any comparison.
        let value = if raw == "null" { None } else { Some(raw) };

        match (name.as_str(), value.as_deref()) {
            ("covariate", Some(value)) => {
                let found = format!("[{}]", value.split(',').collect::<Vec<_>>().join(", "));
                if found != standard {
                    return Err(RecalibrationReportError::NonStandardCovariates {
                        supported: standard.to_string(),
                        found,
                    });
                }
            }
            ("mismatches_context_size", Some(value)) => {
                arguments.mismatches_context_size = value.parse().unwrap_or(2);
            }
            ("indels_context_size", Some(value)) => {
                arguments.indels_context_size = value.parse().unwrap_or(3);
            }
            ("maximum_cycle_value", Some(value)) => {
                arguments.maximum_cycle_value = value.parse().unwrap_or(500);
            }
            ("low_quality_tail", Some(value)) => {
                arguments.low_qual_tail = value.parse().unwrap_or(2);
            }
            ("quantizing_levels", Some(value)) => {
                quantizing_levels = value.parse().unwrap_or(16);
            }
            _ => {}
        }
    }
    Ok((arguments, quantizing_levels))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_null_argument_is_an_absence_and_not_four_characters() {
        // The `covariate` row is skipped when its value is `null`, so a report that names no
        // covariates at all is not refused.
        let mut table = Table::new("Arguments", "", crate::gatk_report::Sorting::DoNotSort);
        table.add_column("Argument", "%s");
        table.add_column("Value", "%s");
        table.set("0", "Argument", Value::Str("covariate".into()));
        table.set("0", "Value", Value::Str("null".into()));
        assert!(parse_arguments(&table).is_ok());
    }

    #[test]
    fn a_non_standard_covariate_list_is_refused_with_both_lists() {
        let mut table = Table::new("Arguments", "", crate::gatk_report::Sorting::DoNotSort);
        table.add_column("Argument", "%s");
        table.add_column("Value", "%s");
        table.set("0", "Argument", Value::Str("covariate".into()));
        table.set(
            "0",
            "Value",
            Value::Str("ReadGroupCovariate,NonesuchCovariate".into()),
        );
        let error = parse_arguments(&table).unwrap_err();
        assert!(
            error.message().contains("NonesuchCovariate"),
            "{}",
            error.message()
        );
        assert!(
            error.message().contains("ContextCovariate"),
            "{}",
            error.message()
        );
    }
}
