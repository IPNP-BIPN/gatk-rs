//! Conformance for the three Mutect gathers against GATK 4.6.2.0, compared as the whole output
//! file of every run.
//!
//! Golden from `tools/readfilter-conformance/MutectGathersDump.java`.
//!
//! # What this suite is for
//!
//!  * **`MergeMutectStats` writes its sum as a double**, so shards of `1` and `2` give `3.0`;
//!  * **and refuses any statistic outside its aggregation map**;
//!  * **`GatherPileupSummaries` sorts its files by their first record**, not by the order given;
//!  * **and drops the files with no records before sorting**;
//!  * **`GatherNormalArtifactData` concatenates in the order given**, which the `reversed` run
//!    shows by producing the same two rows the other way round.

use gatk_corpus as corpus;
use gatk_tools::mutect_gathers::{
    gather_normal_artifact_data, gather_pileup_summaries, merge_stats, GatherError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mutect_gathers.txt.gz"),
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

fn dictionary() -> Vec<String> {
    vec!["chr1".to_string(), "chr2".to_string()]
}

const SUMMARY_COLUMNS: &str =
    "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n";

/// The shards of each stats run.
fn stats_run(label: &str) -> Vec<&'static str> {
    match label {
        "three-shards" => vec![
            "statistic\tvalue\ncallable\t100.0\n",
            "statistic\tvalue\ncallable\t250.0\n",
            "statistic\tvalue\ncallable\t0.0\n",
        ],
        "integers" => vec![
            "statistic\tvalue\ncallable\t1\n",
            "statistic\tvalue\ncallable\t2\n",
        ],
        "one-shard" => vec!["statistic\tvalue\ncallable\t7.5\n"],
        "empty-shard" => vec!["statistic\tvalue\ncallable\t7.5\n", "statistic\tvalue\n"],
        "unknown-statistic" => vec!["statistic\tvalue\ncallable\t1.0\nother\t2.0\n"],
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The shards of each pileup run, each with the name the messages use.
fn pileup_run(label: &str) -> Vec<(String, String)> {
    let with =
        |sample: &str, rows: &str| format!("#<METADATA>SAMPLE={sample}\n{SUMMARY_COLUMNS}{rows}");
    match label {
        "out-of-order" => vec![
            (with("sample", "chr2\t10\t10\t5\t0\t0.5\n"), "0".to_string()),
            (
                with(
                    "sample",
                    "chr1\t10\t20\t2\t0\t0.1\nchr1\t50\t18\t4\t1\t0.2\n",
                ),
                "1".to_string(),
            ),
        ],
        "empty-shard" => vec![
            (with("sample", ""), "0".to_string()),
            (with("sample", "chr1\t10\t20\t2\t0\t0.1\n"), "1".to_string()),
        ],
        "two-samples" => vec![
            (with("first", "chr1\t10\t20\t2\t0\t0.1\n"), "0".to_string()),
            (with("second", "chr1\t50\t18\t4\t1\t0.2\n"), "1".to_string()),
        ],
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The shards of each normal-artifact run.
fn artifact_run(label: &str) -> Vec<&'static str> {
    let columns = "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n";
    let snv =
        "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n1\t20\t5\t30\t1.0\tSNV\n";
    let indel = "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n0\t25\t7\t35\t1.0\tINDEL\n";
    match label {
        "two-shards" => vec![snv, indel],
        "reversed" => vec![indel, snv],
        "empty-shard" => vec![columns, snv],
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_gather_writes_what_the_reference_writes() {
    let text = golden();
    let mut compared = 0;

    for label in ["three-shards", "integers", "one-shard", "empty-shard"] {
        let ours = merge_stats(&stats_run(label)).expect("a statistic the map knows");
        assert_eq!(ours, value(&text, "stats", label), "stats/{label}");
        compared += 1;
    }

    for label in ["out-of-order", "empty-shard"] {
        let shards = pileup_run(label);
        let inputs: Vec<(&str, &str)> = shards
            .iter()
            .map(|(text, source)| (text.as_str(), source.as_str()))
            .collect();
        let ours = gather_pileup_summaries(&inputs, &dictionary()).expect("one sample");
        assert_eq!(ours, value(&text, "pileup", label), "pileup/{label}");
        compared += 1;
    }

    for label in ["two-shards", "reversed", "empty-shard"] {
        let ours = gather_normal_artifact_data(&artifact_run(label));
        assert_eq!(ours, value(&text, "artifact", label), "artifact/{label}");
        compared += 1;
    }

    assert_eq!(compared, 9, "the golden's outputs");
}

#[test]
fn the_refusals_carry_the_references_messages() {
    let text = golden();

    let unknown = merge_stats(&stats_run("unknown-statistic")).expect_err("an unknown statistic");
    assert_eq!(unknown, GatherError::UnknownStatistic("other".to_string()));
    assert_eq!(
        format!("{}:{}", unknown.java_class(), unknown.message()),
        refusal(&text, "stats-unknown-statistic")
    );

    let shards = pileup_run("two-samples");
    let inputs: Vec<(&str, &str)> = shards
        .iter()
        .map(|(text, source)| (text.as_str(), source.as_str()))
        .collect();
    let samples = gather_pileup_summaries(&inputs, &dictionary()).expect_err("two samples");
    assert_eq!(
        format!("{}:{}", samples.java_class(), samples.message()),
        refusal(&text, "pileup-two-samples")
    );
}
