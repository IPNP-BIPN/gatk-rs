//! Conformance for `CalculateAverageCombinedAnnotations` against GATK 4.6.2.0, compared as the
//! whole output file of every run.
//!
//! Golden from `tools/readfilter-conformance/CalculateAverageCombinedAnnotationsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the averages**, down to the encoder's two decimals;
//!  * **a divisor of zero leaving the record untouched**, with no `AVERAGE_` field at all;
//!  * **an annotation a record does not carry being skipped** for that record only;
//!  * **the header gaining a line per requested annotation**, including one the input never
//!    declared;
//!  * **and the refusal for a record with no `RAW_GT_COUNT`**, which comes after the earlier
//!    records were already written.

use gatk_corpus as corpus;
use gatk_tools::calculate_average_combined_annotations::{
    apply, header_with_averages, AverageError,
};
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::vcf_file::write_vcf;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/calculate_average_combined_annotations.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

/// The reference's own `##GATKCommandLine` line, which the golden masks and this port does not
/// write at all: it is dropped from both sides before they are compared.
fn without_command_line(file: &str) -> String {
    file.lines()
        .filter(|line| !line.starts_with("##GATKCommandLine"))
        .collect::<Vec<&str>>()
        .join("\n")
        + "\n"
}

fn run(text: &str, label: &str, annotations: &[&str]) -> Result<String, AverageError> {
    let annotations: Vec<String> = annotations.iter().map(|a| a.to_string()).collect();
    let file = read_vcf(&value(text, "input", label)).expect("a vcf the reader accepts");
    let header = header_with_averages(&file.header, &annotations)?;
    let mut records = Vec::new();
    for record in &file.records {
        records.push(apply(record, &annotations)?);
    }
    Ok(write_vcf(&header, &records).expect("a file the writer accepts"))
}

#[test]
fn every_averaged_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, annotations) in [
        ("two-annotations", vec!["SUMMED", "AS_QD"]),
        ("one-annotation", vec!["SUMMED"]),
        ("undeclared", vec!["NOT_THERE"]),
        ("no-records", vec!["SUMMED"]),
    ] {
        let ours = run(&text, label, &annotations).expect("a run the tool allows");
        assert_eq!(
            without_command_line(&ours),
            without_command_line(&value(&text, "averaged", label)),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 4, "the golden's outputs");
}

/// The three shapes one file produces: an average beside the source, a record left alone because
/// its divisor is zero, and a record whose annotation is absent.
#[test]
fn a_divisor_of_zero_leaves_the_record_untouched() {
    let text = golden();
    let written = value(&text, "averaged", "two-annotations");
    let info: Vec<&str> = written
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').nth(7).expect("an INFO column"))
        .collect();
    assert_eq!(
        info,
        vec![
            "AS_QD=9.0;AVERAGE_AS_QD=3.00;AVERAGE_SUMMED=10.00;RAW_GT_COUNT=0,2,1;SUMMED=30.0",
            // Divisor zero: nothing added, though SUMMED is right there.
            "RAW_GT_COUNT=5,0,0;SUMMED=30.0",
            // AS_QD absent from this record, so only one average.
            "AVERAGE_SUMMED=7.00;RAW_GT_COUNT=0,1,0;SUMMED=7.0",
            "RAW_GT_COUNT=0,3,0",
            "AS_QD=4.0;AVERAGE_AS_QD=2.00;AVERAGE_SUMMED=2.50;RAW_GT_COUNT=0,1,1;SUMMED=5.0",
        ]
    );
}

/// The line is added for an annotation the input never declared.
#[test]
fn the_header_declares_every_requested_average() {
    let text = golden();
    let written = value(&text, "averaged", "undeclared");
    assert!(written.lines().any(|line| line.starts_with(
        "##INFO=<ID=AVERAGE_NOT_THERE,Number=1,Type=Float,Description=\"Average of NOT_THERE"
    )));
    assert!(!written
        .lines()
        .any(|line| line.starts_with("##INFO=<ID=NOT_THERE,")));
}

#[test]
fn a_record_with_no_counts_is_refused_by_its_site() {
    let text = golden();
    let error = run(&text, "missing-counts", &["SUMMED"]).expect_err("the refusal");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "missing-counts")
    );
    // The record before it was written, which the golden keeps as the partial output.
    let partial = value(&text, "averaged", "missing-counts-partial");
    assert_eq!(
        partial
            .lines()
            .filter(|line| !line.starts_with('#'))
            .count(),
        1
    );
    assert!(partial.contains("AVERAGE_SUMMED=10.00"));
}
