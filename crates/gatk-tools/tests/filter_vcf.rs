//! Conformance for `FilterVcf` against Picard 3.4.0, compared as the whole output file of every
//! run.
//!
//! Golden from `tools/readfilter-conformance/FilterVcfDump.java`, which carries each run's input as
//! well as its output.
//!
//! # What this suite is for
//!
//!  * **a passing genotype's `FT` never reaches the file**, so FT appears only on the records
//!    where something was filtered;
//!  * **a record whose genotypes were all filtered is `AllGtsFiltered`**;
//!  * **`LowQD` skips a record with no QD** however high the threshold, and `StrandBias` reaches
//!    one with no FS only when the threshold is negative;
//!  * **the allele-balance filter groups by the genotype's alleles**, so two samples with the same
//!    het call share one tally;
//!  * **and the four FILTER lines and the FT line are added to the header** whatever the file did.

use gatk_corpus as corpus;
use gatk_tools::filter_vcf::{filter, Thresholds};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/filter_vcf.txt.gz"),
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

/// The thresholds of each run.
fn thresholds(label: &str) -> Thresholds {
    let base = Thresholds::default();
    match label {
        "defaults" | "no-records" | "no-contigs" => base,
        "min-ab" => Thresholds {
            min_ab: 0.4,
            ..base
        },
        "max-fs" => Thresholds {
            max_fs: 0.5,
            ..base
        },
        "min-qd" => Thresholds {
            min_qd: 25.0,
            ..base
        },
        "min-gq" => Thresholds { min_gq: 50, ..base },
        "min-dp" => Thresholds { min_dp: 10, ..base },
        "negative-fs" => Thresholds {
            max_fs: -1.0,
            ..base
        },
        "everything" => Thresholds {
            min_ab: 0.45,
            max_fs: 0.5,
            min_qd: 25.0,
            min_gq: 50,
            min_dp: 10,
        },
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_filtered_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "defaults",
        "min-ab",
        "max-fs",
        "min-qd",
        "min-gq",
        "min-dp",
        "negative-fs",
        "everything",
        "no-records",
    ] {
        let ours = filter(&value(&text, "input", label), &thresholds(label))
            .expect("a run the tool allows");
        assert_eq!(ours, value(&text, "filtered", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 9, "the golden's outputs");
}

#[test]
fn a_header_with_no_contigs_is_the_one_refusal() {
    let text = golden();
    let error = filter(
        &value(&text, "input", "no-contigs"),
        &thresholds("no-contigs"),
    )
    .expect_err("the dictionary refusal");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-contigs")
    );
}
