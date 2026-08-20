//! Conformance for `SplitVcfs` against Picard 3.4.0, compared as both output files of every run.
//!
//! Golden from `tools/readfilter-conformance/SplitVcfsDump.java`, which carries each run's input as
//! well as its two outputs.
//!
//! # What this suite is for
//!
//!  * **a MIXED record, an MNP, a symbolic alternate and a monomorphic record are each in neither
//!    file**;
//!  * **a spanning-deletion alternate beside a SNP leaves the record a SNP**, the star being one
//!    base long;
//!  * **`STRICT` is on by default**, so the same input refuses unless it is turned off;
//!  * **both files carry the whole input header** however few records they hold;
//!  * **and `CREATE_INDEX` is on by default**, so a header declaring no contigs is a refusal while
//!    the same file with the index off writes.

use gatk_corpus as corpus;
use gatk_tools::split_vcfs::split;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/split_vcfs.txt.gz"),
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

/// `STRICT` and `CREATE_INDEX` per run, both of which default to true.
fn arguments(label: &str) -> (bool, bool) {
    match label {
        "every-type" => (false, true),
        "no-contigs-no-index" => (true, false),
        _ => (true, true),
    }
}

#[test]
fn both_files_match_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "every-type",
        "snps-only",
        "indels-only",
        "no-records",
        "strict-all-snps",
        "no-contigs-no-index",
    ] {
        let (strict, create_index) = arguments(label);
        let input = value(&text, "input", label);
        let ours = split(&input, strict, create_index).expect("a run the tool allows");
        assert_eq!(ours.snps, value(&text, "snps", label), "{label}: the SNPs");
        assert_eq!(
            ours.indels,
            value(&text, "indels", label),
            "{label}: the indels"
        );
        compared += 1;
    }
    assert_eq!(compared, 6, "the golden's runs");
}

#[test]
fn strict_and_the_index_are_the_two_refusals() {
    let text = golden();
    for label in ["strict", "no-contigs"] {
        let (strict, create_index) = arguments(label);
        let input = value(&text, "input", label);
        let error = split(&input, strict, create_index).expect_err("a refusal");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
    }
}
