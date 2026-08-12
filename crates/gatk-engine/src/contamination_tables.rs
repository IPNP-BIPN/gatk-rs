//! `ContaminationRecord` and `MinorAlleleFractionRecord`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.contamination` (GATK 4.6.2.0).
//!
//! The two tables `CalculateContamination` writes: the answer, and the segments of minor allele
//! fraction it was computed over. Both sit on the same tsv layer as [`crate::pileup_summary`], and
//! they carry the sample differently.
//!
//! # The same sample name is quoted twice over, in two different places
//!
//! The contamination table has **no metadata at all**: its first line is the header and its sample
//! is a column. The segments table writes `#<METADATA>SAMPLE=` before its own. So a sample name
//! holding a tab is quoted as a value in one file and as a whole comment line in the other, for
//! the same name and the same writer.
//!
//! # The interval validates, and it does so late
//!
//! ```java
//! return new MinorAlleleFractionRecord(new SimpleInterval(contig, start, end), maf);
//! ```
//!
//! A row whose end is before its start parses as a table without complaint and is refused when the
//! **record** is built, so the failure is `IllegalArgumentException: Invalid interval. Contig:chr1
//! start:100 end:50` rather than the table's `BadInput`. A start of zero goes the same way, since
//! the interval is one-based.

use crate::interval::SimpleInterval;
use crate::pileup_summary::SAMPLE_METADATA_TAG;
use crate::tsv_table::{java_double_to_string, write_table, Table, TableError};

/// The contamination table's columns, whose sample is one of them.
pub const CONTAMINATION_COLUMNS: [&str; 3] = ["sample", "contamination", "error"];
/// The segments table's columns.
pub const SEGMENT_COLUMNS: [&str; 4] = ["contig", "start", "end", "minor_allele_fraction"];

/// One estimate: a sample, its contamination and the error on it.
#[derive(Debug, Clone, PartialEq)]
pub struct ContaminationRecord {
    pub sample: String,
    pub contamination: f64,
    pub error: f64,
}

/// One segment and the minor allele fraction over it.
#[derive(Debug, Clone, PartialEq)]
pub struct MinorAlleleFractionRecord {
    pub segment: SimpleInterval,
    pub minor_allele_fraction: f64,
}

/// What the two readers refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContaminationTableError {
    /// The table itself was malformed, or a column was missing, which is the tsv layer's refusal.
    Table(TableError),
    /// A segment the interval refused, which happens when the record is built and not when the
    /// file is parsed.
    InvalidInterval {
        contig: String,
        start: i32,
        end: i32,
    },
}

impl ContaminationTableError {
    /// The message the reference carries.
    pub fn message(&self) -> String {
        match self {
            ContaminationTableError::Table(error) => error.message(),
            ContaminationTableError::InvalidInterval { contig, start, end } => {
                format!("Invalid interval. Contig:{contig} start:{start} end:{end}")
            }
        }
    }

    /// The Java class, which is the interval's own for a bad segment.
    pub fn java_class(&self) -> &'static str {
        match self {
            ContaminationTableError::Table(error) => error.java_class(),
            ContaminationTableError::InvalidInterval { .. } => "java.lang.IllegalArgumentException",
        }
    }
}

impl From<TableError> for ContaminationTableError {
    fn from(error: TableError) -> ContaminationTableError {
        ContaminationTableError::Table(error)
    }
}

/// `ContaminationRecord.writeToFile`: no metadata line, so the header is the first line.
pub fn write_contamination(records: &[ContaminationRecord]) -> String {
    let rows: Vec<Vec<String>> = records
        .iter()
        .map(|record| {
            vec![
                record.sample.clone(),
                java_double_to_string(record.contamination),
                java_double_to_string(record.error),
            ]
        })
        .collect();
    write_table(&CONTAMINATION_COLUMNS, &rows, &[])
}

/// `ContaminationRecord.readFromFile`, which takes its columns by name and ignores any metadata.
pub fn read_contamination(
    text: &str,
    source: &str,
) -> Result<Vec<ContaminationRecord>, ContaminationTableError> {
    let table = Table::parse(text, source)?;
    let mut records = Vec::with_capacity(table.rows.len());
    for (index, row) in table.rows.iter().enumerate() {
        let line = table.row_lines[index];
        records.push(ContaminationRecord {
            sample: table.get(row, "sample")?.to_string(),
            contamination: get_double(&table, row, "contamination", source, line)?,
            error: get_double(&table, row, "error", source, line)?,
        });
    }
    Ok(records)
}

/// `MinorAlleleFractionRecord.writeToFile`: the sample as metadata, then the segments.
pub fn write_segments(sample: &str, records: &[MinorAlleleFractionRecord]) -> String {
    let rows: Vec<Vec<String>> = records
        .iter()
        .map(|record| {
            vec![
                record.segment.contig.clone(),
                record.segment.start.to_string(),
                record.segment.end.to_string(),
                java_double_to_string(record.minor_allele_fraction),
            ]
        })
        .collect();
    write_table(&SEGMENT_COLUMNS, &rows, &[(SAMPLE_METADATA_TAG, sample)])
}

/// `MinorAlleleFractionRecord.readFromFile`: the sample, which may be missing, and the segments.
pub fn read_segments(
    text: &str,
    source: &str,
) -> Result<(Option<String>, Vec<MinorAlleleFractionRecord>), ContaminationTableError> {
    let table = Table::parse(text, source)?;
    let mut records = Vec::with_capacity(table.rows.len());
    for (index, row) in table.rows.iter().enumerate() {
        let line = table.row_lines[index];
        let contig = table.get(row, "contig")?.to_string();
        let start = table.get_int(row, "start", source, line)?;
        let end = table.get_int(row, "end", source, line)?;
        let fraction = get_double(&table, row, "minor_allele_fraction", source, line)?;
        // The interval is what refuses a backwards or zero-based segment, once every column has
        // already been read.
        let segment = SimpleInterval::new(&contig, start, end).ok_or(
            ContaminationTableError::InvalidInterval {
                contig: contig.clone(),
                start,
                end,
            },
        )?;
        records.push(MinorAlleleFractionRecord {
            segment,
            minor_allele_fraction: fraction,
        });
    }
    Ok((table.metadata.get(SAMPLE_METADATA_TAG).cloned(), records))
}

/// `SimpleInterval.toString()`, which is `contig:start-end` and not what the columns hold.
pub fn interval_to_string(segment: &SimpleInterval) -> String {
    format!("{}:{}-{}", segment.contig, segment.start, segment.end)
}

/// `DataLine.getDouble`, whose refusal is the int getter's message. See
/// [`crate::pileup_summary`], where the same copy is described.
fn get_double(
    table: &Table,
    row: &[String],
    column: &str,
    source: &str,
    line: u64,
) -> Result<f64, TableError> {
    let value = table.get(row, column)?;
    parse_java_double(value).ok_or_else(|| TableError::NotAnInteger {
        source: source.to_string(),
        line,
        column: column.to_string(),
        value: value.to_string(),
    })
}

fn parse_java_double(value: &str) -> Option<f64> {
    match value {
        "NaN" => Some(f64::NAN),
        "Infinity" | "+Infinity" => Some(f64::INFINITY),
        "-Infinity" => Some(f64::NEG_INFINITY),
        _ => value.parse::<f64>().ok(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_contamination_table_has_no_metadata_line() {
        let text = write_contamination(&[ContaminationRecord {
            sample: "sample-a".to_string(),
            contamination: 0.0,
            error: 0.0,
        }]);
        assert_eq!(text, "sample\tcontamination\terror\nsample-a\t0.0\t0.0\n");
    }

    /// The same name, quoted as a value in one table and as a comment line in the other.
    #[test]
    fn a_tab_in_the_sample_is_quoted_in_two_different_places() {
        let contamination = write_contamination(&[ContaminationRecord {
            sample: "has\ta tab".to_string(),
            contamination: 0.5,
            error: f64::NAN,
        }]);
        assert!(
            contamination.contains("\"has\ta tab\"\t0.5\tNaN"),
            "{contamination}"
        );

        let segments = write_segments("has\ta tab", &[]);
        assert!(
            segments.starts_with("\"#<METADATA>SAMPLE=has\ta tab\"\n"),
            "{segments}"
        );
    }

    #[test]
    fn a_backwards_segment_is_the_intervals_refusal_and_not_the_tables() {
        let text = "#<METADATA>SAMPLE=sample-a\ncontig\tstart\tend\tminor_allele_fraction\nchr1\t100\t50\t0.5\n";
        let error = read_segments(text, "x").unwrap_err();
        assert_eq!(
            error.message(),
            "Invalid interval. Contig:chr1 start:100 end:50"
        );
        assert_eq!(error.java_class(), "java.lang.IllegalArgumentException");
    }

    #[test]
    fn the_columns_are_taken_by_name_and_a_missing_one_names_itself() {
        let reordered =
            read_contamination("error\tsample\tcontamination\n0.01\tsample-a\t0.1\n", "x")
                .expect("the order does not matter");
        assert_eq!(reordered[0].contamination, 0.1);
        assert_eq!(reordered[0].error, 0.01);

        let short = read_contamination("sample\tcontamination\nsample-a\t0.1\n", "x").unwrap_err();
        assert_eq!(short.message(), "there is no such column: error");
    }
}
