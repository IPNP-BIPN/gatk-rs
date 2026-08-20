//! Conformance for `RenameSampleInVcf` against Picard 3.4.0, compared as the whole output file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/RenameSampleInVcfDump.java`, which carries each run's
//! input as well as its output, so this reads the reference's own bytes.
//!
//! # What this suite is for
//!
//!  * **a sites-only VCF is renamed rather than refused**, and comes back with a sample column and
//!    a missing genotype on every record;
//!  * **the new name is not validated**, so a space and a bare number both land in the header;
//!  * **`OLD_SAMPLE_NAME` is checked against the first sample only**, and names what was there;
//!  * **the records go back out through the writer**, so `50.00` comes back `50` and `10.5` comes
//!    back `10.50`;
//!  * **and a record using an undeclared INFO key is the writer's refusal**, not the tool's.

use gatk_corpus as corpus;
use gatk_tools::rename_sample_in_vcf::rename;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/rename_sample_in_vcf.txt.gz"),
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

/// The new name and the asserted old name of each run.
fn arguments(label: &str) -> (&'static str, Option<&'static str>) {
    match label {
        "plain" | "sites-only" | "two-samples" | "undeclared-info" | "no-records" => {
            ("renamed", None)
        }
        "old-name-right" => ("renamed", Some("NA12878")),
        "old-name-wrong" => ("renamed", Some("OTHER")),
        "name-with-space" => ("two words", None),
        "name-that-is-a-number" => ("12345", None),
        "same-name" => ("NA12878", None),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_renamed_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "plain",
        "old-name-right",
        "name-with-space",
        "name-that-is-a-number",
        "same-name",
        "sites-only",
        "no-records",
    ] {
        let (new_name, old_name) = arguments(label);
        let input = value(&text, "input", label);
        let ours = rename(&input, new_name, old_name).expect("a run the tool allows");
        assert_eq!(ours, value(&text, "renamed", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 7, "the golden's outputs");
}

#[test]
fn the_three_refusals_are_two_of_the_tools_and_one_of_the_writers() {
    let text = golden();
    for label in ["old-name-wrong", "two-samples", "undeclared-info"] {
        let (new_name, old_name) = arguments(label);
        let input = value(&text, "input", label);
        let error = rename(&input, new_name, old_name).expect_err("a refusal");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
    }
}
