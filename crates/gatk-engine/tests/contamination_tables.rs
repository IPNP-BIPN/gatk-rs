//! Conformance for the two contamination output tables against GATK 4.6.2.0, compared as written
//! text, as the records read back, and as the refusals.
//!
//! Golden from `tools/readfilter-conformance/ContaminationTablesDump.java`.
//!
//! # What this suite is for
//!
//!  * **one table carries no metadata and the other does**, so the same sample name is quoted as a
//!    value in one and as a whole comment line in the other;
//!  * **both double columns are `Double.toString`'s spelling**;
//!  * **the interval validates late**, so a backwards segment is its refusal and not the table's;
//!  * **and the columns are taken by name**, so another order reads and a missing one names itself.

use gatk_corpus as corpus;
use gatk_engine::contamination_tables::{
    interval_to_string, read_contamination, read_segments, write_contamination, write_segments,
    ContaminationRecord, MinorAlleleFractionRecord,
};
use gatk_engine::interval::SimpleInterval;
use gatk_engine::tsv_table::java_double_to_string as java;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/contamination_tables.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// The reverse of the dump's `escape`, scanning once so a real backslash is never read as a tab.
fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn labelled(text: &str, kind: &str, label: &str) -> String {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no {kind} row for {label}"))[1]
        .to_string()
}

fn written(text: &str, label: &str) -> String {
    unescape(&labelled(text, "written", label))
}

fn contamination_records() -> Vec<ContaminationRecord> {
    let record = |sample: &str, contamination: f64, error: f64| ContaminationRecord {
        sample: sample.to_string(),
        contamination,
        error,
    };
    vec![
        record("sample-a", 0.0, 0.0),
        record("sample-a", 0.05, 0.001),
        record("sample-a", 1.0 / 3.0, 1e-7),
        record("has\ta tab", 0.5, f64::NAN),
    ]
}

fn segment_records() -> Vec<MinorAlleleFractionRecord> {
    let record = |contig: &str, start: i32, end: i32, fraction: f64| MinorAlleleFractionRecord {
        segment: SimpleInterval::new(contig, start, end).expect("a valid segment"),
        minor_allele_fraction: fraction,
    };
    vec![
        record("chr1", 1, 1000, 0.5),
        record("chr1", 1001, 2000, 0.0),
        record("chr2", 5, 5, 1.0 / 3.0),
    ]
}

#[test]
fn every_written_table_is_the_reference() {
    let text = golden();
    assert_eq!(
        write_contamination(&contamination_records()),
        written(&text, "contamination")
    );
    assert_eq!(
        write_contamination(&[]),
        written(&text, "contamination-empty")
    );
    assert_eq!(
        write_segments("sample-a", &segment_records()),
        written(&text, "segments")
    );
    assert_eq!(
        write_segments("sample-a", &[]),
        written(&text, "segments-empty")
    );
}

#[test]
fn every_contamination_record_read_back_is_the_reference() {
    let text = golden();
    for label in [
        "contamination",
        "contamination-with-metadata",
        "contamination-reordered",
    ] {
        let records = read_contamination(&written(&text, label), "x").expect("the table is read");
        let ours: Vec<String> = records
            .iter()
            .map(|record| {
                format!(
                    "{},{},{}",
                    record.sample,
                    java(record.contamination),
                    java(record.error)
                )
            })
            .collect();
        let theirs: Vec<String> = rows(&text, "read")
            .into_iter()
            .filter(|row| row[0] == label && row.len() > 3 && !row[3].is_empty())
            .map(|row| unescape(row[3]))
            .collect();
        assert_eq!(ours, theirs, "read/{label}");
    }

    // The empty table has its header and no record.
    assert!(
        read_contamination(&written(&text, "contamination-empty"), "x")
            .expect("reads")
            .is_empty()
    );
}

#[test]
fn every_segment_read_back_is_the_reference() {
    let text = golden();
    for label in ["segments", "segments-nameless"] {
        let (sample, records) =
            read_segments(&written(&text, label), "x").expect("the table reads");
        let expected: Vec<Vec<&str>> = rows(&text, "read")
            .into_iter()
            .filter(|row| row[0] == label)
            .collect();
        assert_eq!(
            sample.clone().unwrap_or_else(|| "null".to_string()),
            expected[0][1],
            "sample/{label}"
        );
        let ours: Vec<String> = records
            .iter()
            .map(|record| {
                format!(
                    "{},{},{},{},{}",
                    interval_to_string(&record.segment),
                    record.segment.contig,
                    record.segment.start,
                    record.segment.end,
                    java(record.minor_allele_fraction)
                )
            })
            .collect();
        let theirs: Vec<String> = expected.iter().map(|row| row[3].to_string()).collect();
        assert_eq!(ours, theirs, "read/{label}");
    }

    // The empty segments table still carries its sample.
    let (sample, records) = read_segments(&written(&text, "segments-empty"), "x").expect("reads");
    assert_eq!(sample.as_deref(), Some("sample-a"));
    assert!(records.is_empty());
}

/// The interval's two refusals, and the missing column, each with its own class.
#[test]
fn every_refusal_is_the_reference() {
    let text = golden();
    let expected = |label: &str| labelled(&text, "error", label);

    for label in ["segments-backwards", "segments-zero-start"] {
        let error = read_segments(&written(&text, label), "x").expect_err("the interval refuses");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            expected(label),
            "error/{label}"
        );
    }

    let error =
        read_contamination(&written(&text, "contamination-short"), "x").expect_err("no column");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        expected("contamination-short")
    );
}
