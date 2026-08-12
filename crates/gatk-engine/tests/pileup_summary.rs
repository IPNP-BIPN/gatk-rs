//! Conformance for `PileupSummary` against GATK 4.6.2.0, compared as written tables, as the
//! values read back, and as the refusals.
//!
//! Golden from `tools/readfilter-conformance/PileupSummaryDump.java`.
//!
//! # What this suite is for
//!
//!  * **the allele frequency keeps `Double.toString`'s spelling**, because the rounding branch of
//!    `DataLine.set(int, double)` is overwritten by the `return` under it;
//!  * **a sample name with a tab quotes the whole comment line**, and the reader parses it back;
//!  * **three refusals for one mistake**, decided by where the file with no sample sits;
//!  * **an unknown contig sorts first**, because the dictionary's index is -1;
//!  * **and a malformed frequency is reported as a bad integer**, in the int getter's words.

use gatk_corpus as corpus;
use gatk_engine::pileup_summary::{
    compare, gather, read_from_file, write_to_file, PileupSummary, PileupSummaryError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/pileup_summary.txt.gz"),
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

/// The file each label was written as, from the golden itself.
fn written(text: &str, label: &str) -> String {
    unescape(&labelled(text, "written", label))
}

/// The records each written label was built from, mirroring the dump.
fn source_records(label: &str) -> (&'static str, Vec<PileupSummary>) {
    match label {
        "frequencies" => (
            "sample-a",
            vec![
                PileupSummary::new("chr1", 100, 10, 5, 0, 1.0),
                PileupSummary::new("chr1", 200, 10, 5, 2, 0.5),
                PileupSummary::new("chr1", 300, 0, 0, 0, 0.0),
                PileupSummary::new("chr1", 400, 1, 3, 0, 1.0 / 3.0),
                PileupSummary::new("chr1", 500, 7, 0, 0, 1e-7),
                PileupSummary::new("chr1", 600, 2, 2, 0, 2.0),
                PileupSummary::new("chr1", 700, 3, 1, 0, f64::NAN),
                PileupSummary::new("chr1", 800, 3, 1, 0, f64::INFINITY),
            ],
        ),
        "empty" => ("sample-a", vec![]),
        "odd-sample" => (
            "has\ta tab",
            vec![PileupSummary::new("chr1", 1, 1, 1, 1, 0.25)],
        ),
        "part-one" => (
            "sample-a",
            vec![
                PileupSummary::new("chr1", 100, 10, 5, 0, 0.5),
                PileupSummary::new("chr1", 200, 10, 5, 0, 0.5),
            ],
        ),
        "part-two" => (
            "sample-a",
            vec![PileupSummary::new("chr2", 100, 1, 1, 0, 0.25)],
        ),
        "other-sample" => (
            "sample-b",
            vec![PileupSummary::new("chr1", 300, 1, 1, 0, 0.25)],
        ),
        other => panic!("no written case {other}"),
    }
}

#[test]
fn every_written_table_is_the_reference() {
    let text = golden();
    for label in [
        "frequencies",
        "empty",
        "odd-sample",
        "part-one",
        "part-two",
        "other-sample",
    ] {
        let (sample, records) = source_records(label);
        assert_eq!(
            write_to_file(sample, &records),
            written(&text, label),
            "written/{label}"
        );
    }
}

/// Every quantity the record derives, against the reference's own `String.valueOf`.
#[test]
fn every_derived_quantity_is_the_reference() {
    let text = golden();
    let (_, records) = source_records("frequencies");
    let expected: Vec<String> = rows(&text, "derived")
        .into_iter()
        .filter(|row| row[0] == "frequencies")
        .map(|row| row[2].to_string())
        .collect();
    let ours: Vec<String> = records
        .iter()
        .map(|record| {
            format!(
                "{},{},{},{}",
                record.total_count,
                java(record.alt_fraction()),
                java(record.minor_allele_fraction()),
                java(record.ref_frequency())
            )
        })
        .collect();
    assert_eq!(ours, expected);
}

/// `String.valueOf(double)`, which is `Double.toString`.
fn java(value: f64) -> String {
    gatk_engine::tsv_table::java_double_to_string(value)
}

#[test]
fn every_read_record_is_the_reference() {
    let text = golden();
    for label in ["frequencies", "odd-sample", "nameless"] {
        let file = written(&text, label);
        let (sample, records) = read_from_file(&file, "x").expect("the table is read");
        let expected: Vec<Vec<&str>> = rows(&text, "read")
            .into_iter()
            .filter(|row| row[0] == label)
            .collect();

        // The sample the golden holds, which is the string "null" where there was none.
        let theirs = unescape(expected[0][1]);
        assert_eq!(
            sample.clone().unwrap_or_else(|| "null".to_string()),
            theirs,
            "sample/{label}"
        );

        // An empty table has one row in the golden and no record at all.
        if expected[0].len() < 4 || expected[0][3].is_empty() {
            assert!(records.is_empty(), "records/{label}");
            continue;
        }
        let ours: Vec<String> = records
            .iter()
            .map(|record| {
                format!(
                    "{},{},{},{},{},{}",
                    record.contig,
                    record.position,
                    record.ref_count,
                    record.alt_count,
                    record.other_alt_count,
                    java(record.allele_frequency)
                )
            })
            .collect();
        let theirs: Vec<String> = expected.iter().map(|row| row[3].to_string()).collect();
        assert_eq!(ours, theirs, "read/{label}");
    }

    // The empty table carries its sample and no record.
    let (sample, records) = read_from_file(&written(&text, "empty"), "x").expect("read");
    assert_eq!(sample.as_deref(), Some("sample-a"));
    assert!(records.is_empty());
}

#[test]
fn gathering_keeps_the_order_the_files_were_given_in() {
    let text = golden();
    let one = written(&text, "part-one");
    let two = written(&text, "part-two");
    let gathered = gather(&[(&two, "part-two.table"), (&one, "part-one.table")]).expect("gathers");
    assert_eq!(
        gathered,
        unescape(&labelled(&text, "gathered", "same-sample"))
    );
}

#[test]
fn a_missing_sample_fails_by_position_and_a_second_sample_by_name() {
    let text = golden();
    let expected = |label: &str| labelled(&text, "error", label);
    let one = written(&text, "part-one");
    let other = written(&text, "other-sample");
    let nameless = written(&text, "nameless");

    let two_samples = gather(&[(&one, "part-one.table"), (&other, "other-sample.table")])
        .expect_err("two samples are refused");
    assert_eq!(
        format!(
            "{}:Bad input: {}",
            two_samples.java_class(),
            two_samples.message()
        ),
        expected("two-samples")
    );

    let first = gather(&[(&nameless, "nameless.table"), (&one, "part-one.table")])
        .expect_err("the writer refuses a null sample");
    assert_eq!(first, PileupSummaryError::NoSampleToWrite);
    assert_eq!(
        format!("{}:{}", first.java_class(), first.message()),
        expected("nameless-first")
    );

    let second = gather(&[(&one, "part-one.table"), (&nameless, "nameless.table")])
        .expect_err("the comparison refuses a null sample");
    assert_eq!(second, PileupSummaryError::NoSampleToCompare);
    assert_eq!(
        format!("{}:{}", second.java_class(), second.message()),
        expected("nameless-second")
    );
}

#[test]
fn a_contig_the_dictionary_lacks_sorts_first() {
    let text = golden();
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
        labelled(&text, "sorted", "dictionary-order")
    );
}

/// Both malformed fields, of which the frequency is reported in the int getter's words.
#[test]
fn every_malformed_field_is_the_reference() {
    let text = golden();
    let expected = |label: &str| labelled(&text, "error", label);
    let header = "#<METADATA>SAMPLE=sample-a\n\
                  contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n";

    for (label, row) in [
        ("broken", "chr1\tx\t10\t5\t0\t0.5\n"),
        ("broken-frequency", "chr1\t100\t10\t5\t0\tnot-a-number\n"),
    ] {
        let file = format!("{header}{row}");
        let error = read_from_file(&file, &format!("pileupsummary-dump/{label}.table"))
            .expect_err("this row is refused");
        assert_eq!(
            format!("{}:Bad input: {}", error.java_class(), error.message()),
            expected(label),
            "error/{label}"
        );
    }
}
