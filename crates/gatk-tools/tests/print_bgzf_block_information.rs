//! Conformance for `PrintBGZFBlockInformation` against GATK 4.6.2.0, compared as the whole report
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/PrintBgzfBlockInformationDump.java`, whose eight input
//! files travel as base64: a test that built its own bgzf would be inventing the framing this suite
//! is about.
//!
//! # What this suite is for
//!
//!  * **every block's offset and both its sizes**, read from the framing and never decompressed;
//!  * **a premature terminator reported twice**, once above the block that follows it and once in
//!    a summary joined with a comma and no space;
//!  * **the number reported being the terminator's own**, one less than the block being printed;
//!  * **the missing-terminator banner** for a file whose last block carries data;
//!  * **a file that is only a terminator**, which is accepted;
//!  * **the two startup refusals**, a regular gzip and an uncompressed file earning the same
//!    message and a missing file its own;
//!  * **and a truncated block**, whose partial report is kept and whose run still fails.

use gatk_corpus as corpus;
use gatk_tools::print_bgzf_block_information::{is_block_compressed, report, Refusal};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/print_bgzf_block_information.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The `report\t<label>=` rows.
fn expected_report(text: &str, label: &str) -> String {
    let prefix = format!("report\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries report/{label}")),
    )
}

/// The `error\t<label>\t` rows.
fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

fn file(text: &str, label: &str) -> Vec<u8> {
    let prefix = format!("file\t{label}\t");
    corpus::decode_base64(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries file/{label}")),
    )
}

/// The file each label was written under, which the report's first line quotes.
fn name(label: &str) -> &'static str {
    match label {
        "whole" => "whole.gz",
        "no-terminator" => "no-terminator.gz",
        "premature" => "premature.gz",
        "premature-twice" => "premature-twice.gz",
        "terminator-only" => "terminator-only.gz",
        "regular-gzip" => "regular.gz",
        "plain-text" => "plain.txt",
        "truncated" => "truncated.gz",
        other => panic!("{other} is in the golden but not named here"),
    }
}

fn path(label: &str) -> String {
    format!("<dir>/{}", name(label))
}

#[test]
fn every_report_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "whole",
        "no-terminator",
        "premature",
        "premature-twice",
        "terminator-only",
    ] {
        let (ours, refused) = report(&file(&text, label), name(label), &path(label));
        assert_eq!(refused, None, "{label}");
        assert_eq!(ours, expected_report(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 5, "the golden's reports");
}

#[test]
fn the_two_startup_refusals_match_the_golden() {
    let text = golden();
    for label in ["regular-gzip", "plain-text"] {
        assert!(!is_block_compressed(&file(&text, label)), "{label}");
        let error = Refusal::NotBlockCompressed { path: path(label) };
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
    }
    // The well formed files pass the same check.
    for label in ["whole", "premature", "terminator-only"] {
        assert!(is_block_compressed(&file(&text, label)), "{label}");
    }
}

#[test]
fn a_missing_file_is_refused_by_its_path() {
    let text = golden();
    let error = Refusal::DoesNotExist {
        path: "<dir>/absent.gz".to_string(),
    };
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "absent")
    );
}

/// The truncated run keeps the blocks it managed to print and fails afterwards, so both halves are
/// compared: the partial report on disk and the exception.
#[test]
fn a_truncated_block_keeps_the_report_it_wrote() {
    let text = golden();
    let (ours, refused) = report(
        &file(&text, "truncated"),
        name("truncated"),
        &path("truncated"),
    );
    assert_eq!(ours, expected_report(&text, "truncated"));
    let error = refused.expect("the parse failure");
    assert_eq!(
        error,
        Refusal::PrematureEndOfFile {
            path: path("truncated")
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "truncated")
    );
}
