//! `PileupSummary`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.contamination.PileupSummary` (GATK 4.6.2.0).
//!
//! The record every contamination tool passes around: six columns, a sample carried as metadata,
//! and a handful of quantities derived from the counts. Nothing here is a walker; it is the table
//! `GetPileupSummaries` writes, `GatherPileupSummaries` concatenates and `CalculateContamination`
//! reads.
//!
//! # Three refusals for one mistake
//!
//! ```java
//! sample = reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG);
//! writer.writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
//! ...
//! final String thisSample = reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG);
//! if (! thisSample.equals(sample)) {
//! ```
//!
//! A file with no `SAMPLE` metadata fails differently depending on where it sits in the list. First
//! in the list, the missing value reaches `writeMetadata`, whose `Utils.nonNull` gives a plain
//! `IllegalArgumentException` reading "Null object is not allowed here."; anywhere after, it
//! reaches `thisSample.equals`, which is a `NullPointerException`. Only two files of two different
//! samples get the message the code meant to give. A reader on its own takes the same file in its
//! stride and returns no sample beside the records.
//!
//! # The comparator orders by an index that can be -1
//!
//! ```java
//! final int contigIndex1 = contigsInOrder.indexOf(ps1.getContig());
//! ```
//!
//! `indexOf` is -1 for a contig the dictionary does not have, and -1 sorts before 0, so an unplaced
//! contig comes out **first** rather than last or refused.
//!
//! # Nothing bounds the allele frequency
//!
//! The reference frequency is `1 - alleleFrequency` with no check, so a table holding 2.0 gives
//! -1.0, and a `NaN` propagates. Only the alt fraction is guarded, and only against a total of
//! zero, where it gives 0 rather than `NaN`.

use crate::base_utils::simple_base_to_base_index;
use crate::tsv_table::{java_double_to_string, write_table, Table, TableError};
use std::cmp::Ordering;

/// `TableUtils.SAMPLE_METADATA_TAG`, which is upper case where the column names are not.
pub const SAMPLE_METADATA_TAG: &str = "SAMPLE";

/// The columns, in the order the writer emits them.
pub const COLUMNS: [&str; 6] = [
    "contig",
    "position",
    "ref_count",
    "alt_count",
    "other_alt_count",
    "allele_frequency",
];

/// One site's counts and the frequency of its alternate allele.
#[derive(Debug, Clone, PartialEq)]
pub struct PileupSummary {
    pub contig: String,
    pub position: i32,
    pub ref_count: i32,
    pub alt_count: i32,
    pub other_alt_count: i32,
    /// `refCount + altCount + otherAltsCount`, derived here and never written to the table.
    pub total_count: i32,
    pub allele_frequency: f64,
}

impl PileupSummary {
    /// The constructor the table reader uses, whose total is the sum of the three counts.
    ///
    /// The other constructor, from a `VariantContext` and a pileup, sums the four base counts
    /// instead and derives `otherAlts` from that total, so a base the pileup could not index is
    /// counted there and not here.
    pub fn new(
        contig: &str,
        position: i32,
        ref_count: i32,
        alt_count: i32,
        other_alt_count: i32,
        allele_frequency: f64,
    ) -> PileupSummary {
        PileupSummary {
            contig: contig.to_string(),
            position,
            ref_count,
            alt_count,
            other_alt_count,
            total_count: ref_count + alt_count + other_alt_count,
            allele_frequency,
        }
    }

    /// The other constructor: `PileupSummary(VariantContext, ReadPileup)`.
    ///
    /// The counts come from `getBaseCounts`, which counts `A`, `C`, `G` and `T` alone: a deletion
    /// at the site is skipped and an `N` is not counted, so `totalCount` is over those four bases
    /// and is NOT the pileup's depth. `otherAltsCount` is what is left after the reference and the
    /// first alternate, so a site whose pileup holds nothing but `N`s summarises as all zeroes.
    ///
    /// A base that is not `ACGT` in the record's own alleles indexes at -1 upstream, which is an
    /// array access with a negative index and therefore an exception; nothing that reaches this
    /// constructor through `GetPileupSummaries` can, because the record has to be a SNP first.
    pub fn from_base_counts(
        contig: &str,
        position: i32,
        allele_frequency: f64,
        reference_base: u8,
        alternate_base: u8,
        base_counts: [i32; 4],
    ) -> Option<PileupSummary> {
        let alt_index = simple_base_to_base_index(alternate_base);
        let ref_index = simple_base_to_base_index(reference_base);
        if alt_index < 0 || ref_index < 0 {
            return None;
        }
        let alt_count = base_counts[alt_index as usize];
        let ref_count = base_counts[ref_index as usize];
        let total_count: i32 = base_counts.iter().sum();
        Some(PileupSummary {
            contig: contig.to_string(),
            position,
            ref_count,
            alt_count,
            other_alt_count: total_count - alt_count - ref_count,
            total_count,
            allele_frequency,
        })
    }

    /// `getAltFraction`, whose only guard is against an empty site.
    pub fn alt_fraction(&self) -> f64 {
        if self.total_count == 0 {
            0.0
        } else {
            f64::from(self.alt_count) / f64::from(self.total_count)
        }
    }

    /// `getMinorAlleleFraction`: the alt fraction against everything else, not the rarer base.
    pub fn minor_allele_fraction(&self) -> f64 {
        let fraction = self.alt_fraction();
        fraction.min(1.0 - fraction)
    }

    /// `getRefFrequency`, which nothing bounds.
    pub fn ref_frequency(&self) -> f64 {
        1.0 - self.allele_frequency
    }
}

/// What the table reader and the gatherer refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PileupSummaryError {
    /// The table itself was malformed, which is the tsv layer's refusal.
    Table(TableError),
    /// Two files of two different samples, which is the message the gatherer meant to give.
    DifferentSamples { first: String, other: String },
    /// The **first** file had no sample, so the missing value reached `writeMetadata`.
    NoSampleToWrite,
    /// A **later** file had no sample, so the missing value reached `equals`.
    NoSampleToCompare,
}

impl PileupSummaryError {
    /// The message the reference carries.
    pub fn message(&self) -> String {
        match self {
            PileupSummaryError::Table(error) => error.message(),
            PileupSummaryError::DifferentSamples { first, other } => format!(
                "Combining PileupSummaryTables from different samples is not supported. Got samples {first} and {other}"
            ),
            PileupSummaryError::NoSampleToWrite => "Null object is not allowed here.".to_string(),
            PileupSummaryError::NoSampleToCompare => {
                "Cannot invoke \"String.equals(Object)\" because \"thisSample\" is null".to_string()
            }
        }
    }

    /// The Java class, which is a different one for each of the three ways a sample can be missing.
    pub fn java_class(&self) -> &'static str {
        match self {
            PileupSummaryError::Table(error) => error.java_class(),
            PileupSummaryError::DifferentSamples { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
            PileupSummaryError::NoSampleToWrite => "java.lang.IllegalArgumentException",
            PileupSummaryError::NoSampleToCompare => "java.lang.NullPointerException",
        }
    }
}

impl From<TableError> for PileupSummaryError {
    fn from(error: TableError) -> PileupSummaryError {
        PileupSummaryError::Table(error)
    }
}

/// `writeToFile(sample, records, outputTable)`: the sample as metadata, then the header, then the
/// rows.
///
/// The header is written whether or not a record follows, so an empty run still leaves a readable
/// table.
pub fn write_to_file(sample: &str, records: &[PileupSummary]) -> String {
    let rows: Vec<Vec<String>> = records
        .iter()
        .map(|record| {
            vec![
                record.contig.clone(),
                record.position.to_string(),
                record.ref_count.to_string(),
                record.alt_count.to_string(),
                record.other_alt_count.to_string(),
                // The one column that is not an integer, and the one that keeps its point.
                java_double_to_string(record.allele_frequency),
            ]
        })
        .collect();
    write_table(&COLUMNS, &rows, &[(SAMPLE_METADATA_TAG, sample)])
}

/// `readFromFile(tableFile)`: the sample the metadata carried, which may be missing, and the
/// records.
pub fn read_from_file(
    text: &str,
    source: &str,
) -> Result<(Option<String>, Vec<PileupSummary>), PileupSummaryError> {
    let table = Table::parse(text, source)?;
    let mut records = Vec::with_capacity(table.rows.len());
    for (index, row) in table.rows.iter().enumerate() {
        let line = table.row_lines[index];
        let contig = table.get(row, "contig")?.to_string();
        let position = table.get_int(row, "position", source, line)?;
        let ref_count = table.get_int(row, "ref_count", source, line)?;
        let alt_count = table.get_int(row, "alt_count", source, line)?;
        let other_alt_count = table.get_int(row, "other_alt_count", source, line)?;
        let allele_frequency = get_double(&table, row, "allele_frequency", source, line)?;
        records.push(PileupSummary::new(
            &contig,
            position,
            ref_count,
            alt_count,
            other_alt_count,
            allele_frequency,
        ));
    }
    Ok((table.metadata.get(SAMPLE_METADATA_TAG).cloned(), records))
}

/// `DataLine.getDouble`, whose refusal is the **int** getter's message.
///
/// ```java
/// throw formatErrorFactory.apply(String.format("expected int value for column %s but found %s", ...));
/// ```
///
/// The line is copied twice into the double getter, once for each branch, so a malformed frequency
/// is reported as a bad integer. It reaches the user, so the port says the same thing.
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

/// `Double.parseDouble`, for the spellings the writer produces and the ones a hand-written table
/// may hold.
fn parse_java_double(value: &str) -> Option<f64> {
    match value {
        "NaN" => Some(f64::NAN),
        "Infinity" | "+Infinity" => Some(f64::INFINITY),
        "-Infinity" => Some(f64::NEG_INFINITY),
        _ => value.parse::<f64>().ok(),
    }
}

/// `writeToFile(inputFiles, output)`: several tables concatenated, in the order they are given.
///
/// Each input is its text and the name the messages use. The sample is the **first** file's, and
/// every other file is compared against it, which is where the two missing-sample failures live.
pub fn gather(inputs: &[(&str, &str)]) -> Result<String, PileupSummaryError> {
    let mut sample: Option<String> = None;
    let mut records = Vec::new();
    for (index, (text, source)) in inputs.iter().enumerate() {
        let (this_sample, mut these_records) = read_from_file(text, source)?;
        if index == 0 {
            // `writeMetadata` is called before anything is compared, so a missing sample here is
            // the writer's refusal and not the comparison's.
            sample = Some(
                this_sample
                    .clone()
                    .ok_or(PileupSummaryError::NoSampleToWrite)?,
            );
        }
        let this_sample = this_sample.ok_or(PileupSummaryError::NoSampleToCompare)?;
        let first = sample.clone().unwrap_or_default();
        if this_sample != first {
            return Err(PileupSummaryError::DifferentSamples {
                first,
                other: this_sample,
            });
        }
        records.append(&mut these_records);
    }
    Ok(write_to_file(&sample.unwrap_or_default(), &records))
}

/// `PileupSummaryComparator`: the dictionary's order, and -1 for a contig it does not have.
pub fn compare(dictionary: &[String], first: &PileupSummary, second: &PileupSummary) -> Ordering {
    let index = |contig: &str| {
        dictionary
            .iter()
            .position(|name| name == contig)
            .map_or(-1i64, |position| position as i64)
    };
    let (left, right) = (index(&first.contig), index(&second.contig));
    if left != right {
        left.cmp(&right)
    } else {
        first.position.cmp(&second.position)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_empty_site_has_a_fraction_of_zero_and_not_a_nan() {
        let empty = PileupSummary::new("chr1", 300, 0, 0, 0, 0.0);
        assert_eq!(empty.total_count, 0);
        assert_eq!(empty.alt_fraction(), 0.0);
        assert_eq!(empty.minor_allele_fraction(), 0.0);
    }

    #[test]
    fn nothing_bounds_the_reference_frequency() {
        assert_eq!(
            PileupSummary::new("chr1", 600, 2, 2, 0, 2.0).ref_frequency(),
            -1.0
        );
        assert!(PileupSummary::new("chr1", 700, 3, 1, 0, f64::NAN)
            .ref_frequency()
            .is_nan());
    }

    #[test]
    fn an_unknown_contig_sorts_before_every_known_one() {
        let dictionary: Vec<String> = ["chr1", "chr2", "chr10"]
            .iter()
            .map(|name| name.to_string())
            .collect();
        let mut records = [
            PileupSummary::new("chr10", 5, 1, 1, 0, 0.5),
            PileupSummary::new("chr2", 50, 1, 1, 0, 0.5),
            PileupSummary::new("chrUn", 1, 1, 1, 0, 0.5),
            PileupSummary::new("chr1", 900, 1, 1, 0, 0.5),
            PileupSummary::new("chr1", 100, 1, 1, 0, 0.5),
        ];
        records.sort_by(|first, second| compare(&dictionary, first, second));
        let places: Vec<String> = records
            .iter()
            .map(|record| format!("{}:{}", record.contig, record.position))
            .collect();
        assert_eq!(
            places.join(","),
            "chrUn:1,chr1:100,chr1:900,chr2:50,chr10:5"
        );
    }

    #[test]
    fn a_missing_sample_fails_differently_at_each_position() {
        let nameless = "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\nchr1\t100\t10\t5\t0\t0.5\n".to_string();
        let named = write_to_file("sample-a", &[PileupSummary::new("chr1", 100, 1, 1, 0, 0.5)]);

        let first = gather(&[(&nameless, "nameless.table"), (&named, "part.table")]).unwrap_err();
        assert_eq!(first, PileupSummaryError::NoSampleToWrite);
        assert_eq!(first.java_class(), "java.lang.IllegalArgumentException");

        let second = gather(&[(&named, "part.table"), (&nameless, "nameless.table")]).unwrap_err();
        assert_eq!(second, PileupSummaryError::NoSampleToCompare);
        assert_eq!(second.java_class(), "java.lang.NullPointerException");

        // On its own, the same file is read without complaint.
        let (sample, records) = read_from_file(&nameless, "nameless.table").expect("reads");
        assert_eq!(sample, None);
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn a_bad_frequency_is_reported_as_a_bad_integer() {
        let text = "#<METADATA>SAMPLE=sample-a\ncontig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\nchr1\t100\t10\t5\t0\tnot-a-number\n";
        let error = read_from_file(text, "broken.table").unwrap_err();
        assert_eq!(
            error.message(),
            "format error in 'broken.table' at line 3: expected int value for column allele_frequency but found not-a-number"
        );
    }
}
