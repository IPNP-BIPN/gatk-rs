//! Conformance for the tranches file against GATK 4.6.2.0, compared as the FILTER header lines
//! `ApplyVQSR` writes and as every refusal reading one can produce.
//!
//! Golden from `tools/readfilter-conformance/ApplyVqsrTranchesDump.java`.
//!
//! # What this suite is for
//!
//!  * **the two optional columns are not optional**, `numKnown` defaulting to `-1` and then being
//!    refused for being negative;
//!  * **the same reader words its refusals differently**, two of them naming no file at all;
//!  * **a `numNovel` past an int cannot be read back**, though the field is a `long`;
//!  * **the last tranche never becomes a filter and the first becomes two**;
//!  * **and the FILTER IDs are the `filterName` column**, not anything synthesized.

use gatk_corpus as corpus;
use gatk_engine::tranches::read_tranches;
use gatk_tools::apply_vqsr::{keep, low_vqslod_filter_line, tranche_filter_lines};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/apply_vqsr_tranches.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The whole text of one tranches file, and the path the reference read it from.
fn tranches_file(text: &str, label: &str) -> (String, String) {
    let body = unescape(
        rows(text, "tranches")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no tranches file {label}"))[1],
    );
    (format!("apply-vqsr-tranches-dump/{label}.tranches"), body)
}

/// The `##FILTER` lines of one run, in the order the golden holds them.
fn filters(text: &str, run: &str) -> Vec<String> {
    rows(text, "filter")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .collect()
}

/// The class and the message of one run's refusal.
fn refusal(text: &str, run: &str) -> (String, String) {
    let row = rows(text, "error")
        .into_iter()
        .find(|row| row[0] == run)
        .unwrap_or_else(|| panic!("no refusal {run}"));
    let (class, message) = row[1].split_once(':').expect("class and message");
    (class.to_string(), unescape(message))
}

/// Every run that reached the end: its tranches file, its level, and the lines it wrote.
///
/// htsjdk adds a `##FILTER=<ID=PASS,...>` line of its own, which is not the tool's and is dropped.
fn lines_of(text: &str, file: &str, level: Option<f64>) -> Vec<String> {
    let (path, body) = tranches_file(text, file);
    let tranches = read_tranches(&path, &body).expect("this file reads");
    match level {
        Some(level) => {
            let kept = keep(&tranches, level);
            tranche_filter_lines(&kept, level)
                .expect("this level keeps something")
                .iter()
                .map(|line| line.to_line())
                .collect()
        }
        None => vec![low_vqslod_filter_line(None).to_line()],
    }
}

fn written(text: &str, run: &str) -> Vec<String> {
    let mut lines: Vec<String> = filters(text, run)
        .into_iter()
        .filter(|line| !line.starts_with("##FILTER=<ID=PASS,"))
        .collect();
    lines.sort();
    lines
}

fn ours(text: &str, file: &str, level: Option<f64>) -> Vec<String> {
    let mut lines = lines_of(text, file, level);
    // The reference collects them into a HashSet and htsjdk writes them sorted, so only the set is
    // the measurement here.
    lines.sort();
    lines
}

#[test]
fn every_filter_line_matches_the_golden_byte_for_byte() {
    let text = golden();
    for (run, file, level) in [
        ("level-99", "tranches", Some(99.0)),
        ("level-0", "tranches", Some(0.0)),
        ("custom-names", "custom-names", Some(0.0)),
        ("no-level", "tranches", None),
    ] {
        assert_eq!(ours(&text, file, level), written(&text, run), "{run}");
    }
}

#[test]
fn the_last_tranche_never_becomes_a_filter_and_the_first_becomes_two() {
    let text = golden();
    let lines = written(&text, "level-0");
    // Three tranches in the file, three lines, and two names between them.
    assert_eq!(lines.len(), 3);
    assert!(lines
        .iter()
        .any(|line| line.contains("ID=VQSRTrancheSNP99.00to100.00+,")));
    assert!(lines
        .iter()
        .any(|line| line.contains("ID=VQSRTrancheSNP99.00to100.00,")));
    // The tranche kept whole is named by nothing.
    assert!(!lines
        .iter()
        .any(|line| line.contains("VQSRTrancheSNP0.00to90.00")));
}

#[test]
fn the_filter_ids_are_the_files_own_names() {
    let text = golden();
    let lines = written(&text, "custom-names");
    assert!(lines.iter().any(|line| line.contains("ID=tight+,")));
    assert!(lines.iter().any(|line| line.contains("ID=middling,")));
    // `loose` is the tranche kept whole.
    assert!(!lines.iter().any(|line| line.contains("loose")));
}

#[test]
fn every_refusal_carries_the_references_class_and_words() {
    let text = golden();
    for (run, file) in [
        ("missing-required", "missing-required"),
        ("missing-optional", "missing-optional"),
        ("bad-model", "bad-model"),
        ("missing-model", "missing-model"),
        ("invalid-value", "invalid-value"),
        ("novel-past-int", "novel-past-int"),
        ("unreasonable", "unreasonable"),
        ("short-row", "short-row"),
        ("short-header", "short-header"),
    ] {
        let (path, body) = tranches_file(&text, file);
        let error = read_tranches(&path, &body).expect_err("this file does not read");
        let (class, message) = refusal(&text, run);
        assert_eq!(error.class(), class, "{run}");
        assert_eq!(error.message(), message, "{run}");
    }
}

#[test]
fn two_of_the_refusals_name_no_file_at_all() {
    let text = golden();
    let (_, missing) = refusal(&text, "missing-required");
    assert!(
        missing.starts_with("Unknown file is malformed:"),
        "{missing}"
    );
    let (_, invalid) = refusal(&text, "invalid-value");
    assert!(
        invalid.starts_with("Unknown file is malformed:"),
        "{invalid}"
    );
    // While the two length checks do name it.
    let (_, short) = refusal(&text, "short-row");
    assert!(
        short.starts_with("File apply-vqsr-tranches-dump/short-row.tranches is malformed:"),
        "{short}"
    );
}

#[test]
fn a_num_novel_past_an_int_is_the_same_refusal_as_a_word() {
    let text = golden();
    let (_, big) = refusal(&text, "novel-past-int");
    let (_, word) = refusal(&text, "invalid-value");
    assert_eq!(big, word);
    assert!(big.ends_with("Invalid value for key numNovel"));
}

#[test]
fn a_level_above_every_tranche_is_a_refusal_of_the_tools_own() {
    let text = golden();
    let (path, body) = tranches_file(&text, "tranches");
    let tranches = read_tranches(&path, &body).expect("a good file");
    let kept = keep(&tranches, 100.1);
    let error = tranche_filter_lines(&kept, 100.1).expect_err("nothing survives");
    let (class, message) = refusal(&text, "level-above-everything");
    assert_eq!(error.class(), class);
    assert_eq!(error.message(), message);
}
