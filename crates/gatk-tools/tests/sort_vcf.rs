//! Conformance for `SortVcf` against Picard 3.4.0, compared as the whole output file of every run.
//!
//! Golden from `tools/readfilter-conformance/SortVcfDump.java`, which carries every input of every
//! run as well as its output. The dictionary declares `chr10` before `chr2`, so the order is
//! visibly the dictionary's rather than the alphabet's.
//!
//! # What this suite is for
//!
//!  * **ties keep the order the inputs were given**, which `two-files` and `two-files-reversed`
//!    separate;
//!  * **the header is the smart merge of every input's**, so a line only one file declared is in
//!    the output;
//!  * **the sample check is on the sorted names**, so two files whose columns are ordered
//!    differently agree, and the output's columns are the sorted names;
//!  * **and a record on an undeclared contig is a null pointer**, not a message.
//!
//! Two refusals name a path in the reference and are masked in the golden, so this port's messages
//! carry `<masked>` in the same place.

use gatk_corpus as corpus;
use gatk_tools::sort_vcf::sort;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sort_vcf.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// Every input of one run, in the order the harness wrote them.
fn inputs(text: &str, label: &str) -> Vec<String> {
    let mut found = Vec::new();
    for index in 0..8 {
        let prefix = format!("input\t{label}/{index}=");
        match text
            .lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
        {
            Some(value) => found.push(unescape(value)),
            None => break,
        }
    }
    assert!(!found.is_empty(), "the golden carries inputs for {label}");
    found
}

fn output(text: &str, label: &str) -> String {
    let prefix = format!("sorted\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries sorted/{label}")),
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
fn every_sorted_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "one-file",
        "two-files",
        "two-files-reversed",
        "merged-header",
        "conflicting-header",
        "two-samples",
        "sample-order-differs",
        "no-records",
    ] {
        let ours = sort(&inputs(&text, label)).expect("a run the tool allows");
        assert_eq!(ours, output(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 8, "the golden's outputs");
}

#[test]
fn the_three_refusals_are_the_dictionary_twice_and_the_undeclared_contig() {
    let text = golden();
    for label in [
        "different-dictionaries",
        "no-dictionary",
        "undeclared-contig",
    ] {
        let error = sort(&inputs(&text, label)).expect_err("a refusal");
        let theirs = refusal(&text, label);
        // The dictionary mismatch names both dictionaries in full; the port carries the sentence
        // that opens it and nothing more, so that row is compared by its beginning.
        if label == "different-dictionaries" {
            assert!(
                theirs.starts_with(&format!("{}:{}", error.java_class(), error.message())),
                "{label}: {theirs}"
            );
        } else {
            assert_eq!(
                format!("{}:{}", error.java_class(), error.message()),
                theirs,
                "{label}"
            );
        }
    }
}
