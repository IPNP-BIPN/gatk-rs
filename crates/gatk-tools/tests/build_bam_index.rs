//! Conformance for `BuildBamIndex` against Picard 3.4.0, compared as the whole `.bai` of every run
//! and as the path each index landed on.
//!
//! Golden from `tools/readfilter-conformance/BuildBamIndexDump.java`, which carries every input BAM
//! and every index it wrote as base64.
//!
//! # What this suite is for
//!
//!  * **the index is byte for byte the reference's**, unmapped reads at the end included;
//!  * **an empty BAM still produces one**, a bin-less reference per sequence;
//!  * **the default output lands beside the process and not beside the input**, and replaces the
//!    extension only for a name ending `.bam`;
//!  * **a queryname header and an unsorted one are refused by the same message**, read from `SO`
//!    alone;
//!  * **and a sam file is refused by its type**, whatever its name.

use gatk_corpus as corpus;
use gatk_tools::build_bam_index::{build, default_output, is_bam};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/build_bam_index.txt.gz"),
    )
}

fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .to_string()
}

fn bytes(text: &str, kind: &str, label: &str) -> Vec<u8> {
    corpus::decode_base64(&row(text, kind, label))
}

#[test]
fn every_index_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, fixture) in [
        ("sorted", "sorted"),
        ("empty", "empty"),
        ("default-output", "sorted"),
        ("default-output-odd-name", "sorted"),
    ] {
        let ours = build(&bytes(&text, "fixture", fixture)).expect("a run the tool allows");
        assert_eq!(ours, bytes(&text, "index", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 4, "the golden's indexes");
}

#[test]
fn the_default_output_lands_beside_the_process() {
    let text = golden();
    // The input is `build-bam-index-dump/sorted.bam` and the index landed on `sorted.bai`, in the
    // working directory: the directory is dropped, not kept.
    assert_eq!(
        row(&text, "wrote", "default-output"),
        default_output("sorted.bam")
    );
    assert_eq!(default_output("sorted.bam"), "sorted.bai");
    // A name that does not end `.bam` keeps all of itself.
    assert_eq!(
        row(&text, "wrote", "default-output-odd-name"),
        default_output("sorted.bam.copy")
    );
    assert_eq!(default_output("sorted.bam.copy"), "sorted.bam.copy.bai");
}

#[test]
fn the_three_refusals_match_the_golden() {
    let text = golden();
    for (label, fixture) in [
        ("queryname", "queryname"),
        ("unsorted", "unsorted"),
        ("plain-sam", "plain-sam"),
    ] {
        let error = build(&bytes(&text, "fixture", fixture)).expect_err("a refusal");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            row(&text, "error", label),
            "{label}"
        );
    }
}

/// The sam file is refused for what it is and not for what it is called: the dump named it
/// `plain.sam`, and the port is handed the bytes with no name at all.
#[test]
fn a_sam_file_is_not_a_bam_by_its_bytes() {
    let text = golden();
    assert!(!is_bam(&bytes(&text, "fixture", "plain-sam")));
    assert!(is_bam(&bytes(&text, "fixture", "sorted")));
    assert!(is_bam(&bytes(&text, "fixture", "queryname")));
}
