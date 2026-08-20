//! Conformance for `GatherBQSRReports` against GATK 4.6.2.0, compared as the whole gathered report
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/GatherBQSRReportsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the gather is a sum, not a concatenation**: `same-twice` doubles every observation and
//!    every error, and its empirical qualities are the qualities of the doubled counts;
//!  * **the quantization is recomputed** from the summed quality-score table;
//!  * **the arguments table is the first input's**, and holds the recalibration argument
//!    collection rather than the command line, so `two-shards` and `reversed` are the same bytes;
//!  * **an empty shard is skipped** by the combine;
//!  * **and a gather of nothing but empty shards refuses**.
//!
//! The shard tables are the golden's own: they come from `BaseRecalibrator`, whose suite already
//! pins them, and are fed to the gather here as text.

use gatk_corpus as corpus;
use gatk_tools::gather_bqsr_reports::{gather, GatherError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/gather_bqsr_reports.txt.gz"),
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

/// The shards of each gather, by the golden's own shard rows.
fn shards(text: &str, label: &str) -> Vec<String> {
    let first = value(text, "shard", "first");
    let second = value(text, "shard", "second");
    let empty = value(text, "shard", "empty");
    match label {
        "two-shards" => vec![first, second],
        "reversed" => vec![second, first],
        "one-shard" => vec![first],
        "same-twice" => vec![first.clone(), first],
        "with-empty" => vec![first, empty],
        "all-empty" => vec![empty.clone(), empty],
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_gather_writes_what_the_reference_writes() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "two-shards",
        "reversed",
        "one-shard",
        "same-twice",
        "with-empty",
    ] {
        let owned = shards(&text, label);
        let inputs: Vec<&str> = owned.iter().map(String::as_str).collect();
        let ours = gather(&inputs).expect("a gather with usable data");
        assert_eq!(ours, value(&text, "gathered", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 5, "the golden's gathers");
}

#[test]
fn a_gather_of_empty_shards_refuses() {
    let text = golden();
    let owned = shards(&text, "all-empty");
    let inputs: Vec<&str> = owned.iter().map(String::as_str).collect();
    let error = gather(&inputs).expect_err("no usable data");
    assert_eq!(error, GatherError::NoUsableData);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "all-empty")
    );
}

#[test]
fn the_order_of_the_shards_does_not_reach_the_output() {
    let text = golden();
    // The arguments table is the first input's, but what it holds is the recalibration ARGUMENT
    // COLLECTION and not the command line, and every shard of one scattered run has the same one.
    // So gathering the same two shards in either order gives the same bytes, and a port that
    // carried the inputs into the output would not.
    assert_eq!(
        value(&text, "gathered", "two-shards"),
        value(&text, "gathered", "reversed")
    );
}
