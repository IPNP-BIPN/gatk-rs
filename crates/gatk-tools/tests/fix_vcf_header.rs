//! Conformance for `FixVcfHeader` against Picard 3.4.0, compared as the whole output file of every
//! run.
//!
//! Golden from `tools/readfilter-conformance/FixVcfHeaderDump.java`, which carries each run's input
//! and, where there is one, its replacement header.
//!
//! # What this suite is for
//!
//!  * **an invented line is always `Number=.` and `Type=String`**, whatever the value looked like;
//!  * **the six standard FORMAT lines arrive whatever the file uses**, which `no-records` shows by
//!    coming back with more declarations than it started with;
//!  * **a record limit can leave the header wrong**, and the write then refuses at the record that
//!    uses the undeclared key;
//!  * **`ENFORCE_SAME_SAMPLES` names the first index that differs**;
//!  * **and with it off the input's samples are kept**, so a sites-only replacement header still
//!    writes the input's columns.

use gatk_corpus as corpus;
use gatk_tools::fix_vcf_header::{fix, fix_with_header};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fix_vcf_header.txt.gz"),
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

/// One run, whichever of the two paths it takes.
fn run(text: &str, label: &str) -> Result<String, gatk_tools::fix_vcf_header::FixError> {
    let input = value(text, "input", label);
    match label {
        "undeclared" | "nothing-missing" | "no-records" => fix(&input, -1),
        "first-one-record" => fix(&input, 1),
        "replacement-header" | "different-sample" | "incomplete-header" => {
            fix_with_header(&input, &value(text, "header", label), true)
        }
        "different-sample-unenforced" | "sites-only-header" => {
            fix_with_header(&input, &value(text, "header", label), false)
        }
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_fixed_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "undeclared",
        "nothing-missing",
        "no-records",
        "replacement-header",
        "different-sample-unenforced",
        "sites-only-header",
    ] {
        let ours = run(&text, label).expect("a run the tool allows");
        assert_eq!(ours, value(&text, "fixed", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 6, "the golden's outputs");
}

#[test]
fn the_three_refusals_carry_the_references_messages() {
    let text = golden();
    for label in ["first-one-record", "different-sample", "incomplete-header"] {
        let error = run(&text, label).expect_err("a refusal");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
    }
}
