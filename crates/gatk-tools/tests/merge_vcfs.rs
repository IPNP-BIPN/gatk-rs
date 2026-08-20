//! Conformance for `MergeVcfs` against Picard 3.4.0, compared as the whole output file of every
//! run.
//!
//! Golden from `tools/readfilter-conformance/MergeVcfsDump.java`, which carries every input of
//! every run as well as its output.
//!
//! # What this suite is for
//!
//!  * **a tie is not decided by the input order**, so `two-files` and `reversed` write the same
//!    file, which is the opposite of `sort-vcf`;
//!  * **an input whose own records are out of order is a refusal**, not a repair;
//!  * **the contig check is about indices**, so a subset that shifts one is refused as surely as a
//!    reordering;
//!  * **two comments collapse into one**, sharing the key `MergeVcfs.comment`;
//!  * **two identical input headers collapse** before the merge sees them;
//!  * **and the sample check is on the sorted names**, so two files whose columns are ordered
//!    differently agree.

use gatk_corpus as corpus;
use gatk_tools::merge_vcfs::merge;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/merge_vcfs.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

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
    let prefix = format!("merged\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries merged/{label}")),
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

/// The comments each run was given.
fn comments(label: &str) -> Vec<String> {
    if label == "comments" {
        vec!["first note".to_string(), "second note".to_string()]
    } else {
        Vec::new()
    }
}

#[test]
fn every_merged_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "two-files",
        "reversed",
        "one-file",
        "three-files",
        "identical-files",
        "sample-order-differs",
        "comments",
        "empty-file",
    ] {
        let ours = merge(&inputs(&text, label), &comments(label)).expect("a run the tool allows");
        assert_eq!(ours, output(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 8, "the golden's outputs");
}

#[test]
fn a_tie_is_written_the_same_way_whichever_order_the_files_came_in() {
    let text = golden();
    assert_eq!(
        output(&text, "two-files"),
        output(&text, "reversed"),
        "the reference writes the same file both ways"
    );
    assert_eq!(
        merge(&inputs(&text, "two-files"), &[]).expect("a run"),
        merge(&inputs(&text, "reversed"), &[]).expect("a run"),
        "and so does the port"
    );
}

#[test]
fn the_five_refusals_carry_the_references_messages() {
    let text = golden();
    for label in [
        "unsorted-input",
        "subset-contigs",
        "reordered-contigs",
        "different-samples",
        "no-contigs",
    ] {
        let error = merge(&inputs(&text, label), &[]).expect_err("a refusal");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
    }
}
