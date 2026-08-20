//! Conformance for `MakeSitesOnlyVcf` against Picard 3.4.0, compared as the whole output file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/MakeSitesOnlyVcfDump.java`, which carries each run's
//! input as well as its output.
//!
//! # What this suite is for
//!
//!  * **the requested samples are sorted and de-duplicated**, so `zeta S=alpha` comes back
//!    `alpha zeta` and the same name twice is one column;
//!  * **a name the input never carried becomes a column of `./.`**;
//!  * **the INFO fields are not recomputed**, so a one-sample output still carries the whole
//!    file's AC and AN;
//!  * **the default drops the FORMAT column entirely** rather than leaving it empty;
//!  * **and `CREATE_INDEX` is on by default**, so a header declaring no contigs is a refusal while
//!    the same file with the index off writes.

use gatk_corpus as corpus;
use gatk_tools::make_sites_only_vcf::make_sites_only;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/make_sites_only_vcf.txt.gz"),
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

/// The samples asked for and whether the index was on, per run.
fn arguments(label: &str) -> (Vec<String>, bool) {
    let names = |list: &[&str]| list.iter().map(|s| (*s).to_string()).collect();
    match label {
        "sites-only" | "already-sites-only" | "no-contigs" => (Vec::new(), true),
        "one-sample" => (names(&["alpha"]), true),
        "one-sample-twice" => (names(&["alpha", "alpha"]), true),
        "two-samples-unsorted" => (names(&["zeta", "alpha"]), true),
        "all-samples" => (names(&["zeta", "alpha", "middle"]), true),
        "absent-sample" => (names(&["absent"]), true),
        "absent-and-present" => (names(&["absent", "alpha"]), true),
        "no-contigs-no-index" => (Vec::new(), false),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_output_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "sites-only",
        "one-sample",
        "one-sample-twice",
        "two-samples-unsorted",
        "all-samples",
        "absent-sample",
        "absent-and-present",
        "already-sites-only",
        "no-contigs-no-index",
    ] {
        let (samples, create_index) = arguments(label);
        let input = value(&text, "input", label);
        let ours = make_sites_only(&input, &samples, create_index).expect("a run the tool allows");
        assert_eq!(ours, value(&text, "sites", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 9, "the golden's outputs");
}

#[test]
fn a_header_with_no_contigs_cannot_be_indexed() {
    let text = golden();
    let (samples, create_index) = arguments("no-contigs");
    let input = value(&text, "input", "no-contigs");
    let error = make_sites_only(&input, &samples, create_index).expect_err("the index refusal");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-contigs")
    );
}
