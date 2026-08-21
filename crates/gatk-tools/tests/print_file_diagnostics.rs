//! Conformance for `PrintFileDiagnostics` against GATK 4.6.2.0, compared as the whole report.
//!
//! Golden from `tools/readfilter-conformance/PrintFileDiagnosticsDump.java`, whose index files
//! travel as base64.
//!
//! # What this suite is for
//!
//!  * **the BAI report, whole**: bin summaries, hexadecimal chunk offsets, the metadata bin, the
//!    linear index and the no-coordinate count;
//!  * **the empty reference's spacing**, `n_bin=0` where every other line writes `n_bin= 4`;
//!  * **and the refusal for an extension no analyzer claims**, which quotes the raw argument.
//!
//! # What this port does not do
//!
//! The CRAM analyzer, and the CRAI one. The golden records the CRAI report, which is three lines
//! of `CRAIEntry.toString` under a header; reading a `.crai` needs the CRAM crate and the report
//! itself carries no arithmetic, so it waits for the CRAM bricks.

use gatk_corpus as corpus;
use gatk_tools::print_file_diagnostics::{analyzer_for, bai_report, Analyzer, DiagnosticsError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/print_file_diagnostics.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn fixture(text: &str, label: &str) -> Vec<u8> {
    let prefix = format!("fixture\t{label}\t");
    corpus::decode_base64(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries fixture/{label}")),
    )
}

fn report(text: &str, label: &str) -> String {
    let prefix = format!("report\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries report/{label}")),
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

#[test]
fn the_bai_report_matches_the_golden() {
    let text = golden();
    let ours = bai_report(&fixture(&text, "bai")).expect("a bai the tool reads");
    assert_eq!(ours, report(&text, "bai"));
}

/// The reference with no reads is written with its own spacing, which is the difference a port is
/// most likely to smooth over.
#[test]
fn an_empty_reference_is_printed_without_the_space() {
    let text = golden();
    let ours = bai_report(&fixture(&text, "bai")).expect("a bai the tool reads");
    assert!(ours.contains("Reference 1 has n_bin=0\n"));
    assert!(ours.contains("Reference 1 has n_intv=0\n"));
    // While the references that do have bins carry one.
    assert!(ours.contains("Reference 0 has n_bin= 4\n"));
    assert!(ours.contains("Reference 0 has n_intv= 4\n"));
}

/// The pseudo-bin is counted in `n_bin` and printed after the real bins, out of numeric order,
/// with the counts in the same hexadecimal as the offsets.
#[test]
fn the_metadata_bin_is_counted_and_printed_apart() {
    let text = golden();
    let ours = bai_report(&fixture(&text, "bai")).expect("a bai the tool reads");
    let lines: Vec<&str> = ours.lines().collect();
    let bins: Vec<&&str> = lines
        .iter()
        .filter(|line| line.starts_with("  Ref 0 bin "))
        .collect();
    // Three real bins and the pseudo-bin last, though 37450 is the largest number of the four.
    assert_eq!(bins.len(), 4);
    assert!(bins[3].starts_with("  Ref 0 bin 37450 has n_chunk= 2"));
    // Two unmapped reads, so the second pseudo-chunk holds 3 aligned and 0 unaligned for
    // reference 0.
    assert!(ours.contains("     Chunk:  start: 3 end: 0\n"));
    assert_eq!(ours.lines().last(), Some("No Coordinate Count=2"));
}

#[test]
fn the_analyzer_is_chosen_by_the_name() {
    let text = golden();
    assert_eq!(analyzer_for("reads.bai"), Ok(Analyzer::Bai));
    assert_eq!(analyzer_for("reads.cram.crai"), Ok(Analyzer::Crai));
    assert_eq!(analyzer_for("reads.cram"), Ok(Analyzer::Cram));

    let error = analyzer_for("<dir>/notes.txt").expect_err("the refusal");
    assert_eq!(
        error,
        DiagnosticsError::Unsupported {
            raw: "<dir>/notes.txt".to_string()
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "unsupported")
    );
}

/// The CRAI report is in the golden for a later brick: its shape is recorded here so the port that
/// writes it has something to answer to.
#[test]
fn the_crai_report_is_measured_and_not_yet_written() {
    let text = golden();
    let crai = report(&text, "crai");
    assert!(crai.starts_with(
        "\nSeqId AlignmentStart AlignmentSpan ContainerOffset SliceOffset SliceSize\n"
    ));
    assert_eq!(
        crai.lines()
            .filter(|line| line.starts_with('0') || line.starts_with('1'))
            .count(),
        3
    );
}
