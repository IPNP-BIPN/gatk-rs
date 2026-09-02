//! Conformance for `SequenceDictionaryUtils` against GATK 4.6.2.0, over pairs of dictionaries.
//!
//! Golden from `tools/readfilter-conformance/SequenceDictionaryValidationDump.java`. Four rows of
//! `CountVariants`' covering array disagreed on this and nothing else: they pass a BAM to
//! `--input` alongside a VCF to `--variant`, and the two dictionaries are compared before the
//! traversal (#1038).
//!
//! # What this suite is for
//!
//!  * **three of the eight outcomes needing to be asked for**, so an argument decides whether the
//!    same pair is accepted;
//!  * **a length of zero being equivalent to any length**;
//!  * **two EMPTY dictionaries having no common contigs**, which is a refusal;
//!  * **reversing the common contigs being a SUPERSET without the ordering check**;
//!  * **the same relative order at different absolute positions being its own case**;
//!  * **a common subset being accepted unless a superset was required**, the refusal then naming
//!    the missing contigs;
//!  * **and the human-order check needing chr1, chr2 and chr10 by name AND by length.**
//!
//! The golden is committed and re-derived by the `sequence-dictionary-validation` suite on every
//! run; the dump can still be overridden with an environment variable while a harness change is
//! being checked.

use gatk_tools::sequence_dictionary::{compare, validate};
use htsjdk_bam::header::SequenceRecord;

/// The pairs the harness runs, as `name:length` lists in its own order.
const CASES: &[(&str, &[&str], &[&str])] = &[
    (
        "identical",
        &["chr1:100", "chr2:200"],
        &["chr1:100", "chr2:200"],
    ),
    (
        "superset",
        &["chr1:100", "chr2:200", "chr3:300"],
        &["chr1:100", "chr2:200"],
    ),
    ("common-subset", &["chr1:100"], &["chr1:100", "chr2:200"]),
    ("no-common-contigs", &["chr1:100"], &["chrA:100"]),
    ("unequal-lengths", &["chr1:100"], &["chr1:200"]),
    ("zero-length", &["chr1:0"], &["chr1:200"]),
    ("zero-length-both-sides", &["chr1:0"], &["chr1:0"]),
    (
        "reversed",
        &["chr1:100", "chr2:200"],
        &["chr2:200", "chr1:100"],
    ),
    (
        "different-indices",
        &["chrX:10", "chr1:100", "chr2:200"],
        &["chr1:100", "chr2:200"],
    ),
    ("empty-first", &[], &["chr1:100"]),
    ("empty-both", &[], &[]),
    (
        "lexicographic-human-first",
        &["chr1:249250621", "chr10:135534747", "chr2:243199373"],
        &["chr1:249250621", "chr2:243199373", "chr10:135534747"],
    ),
    (
        "lexicographic-human-second",
        &["chr1:249250621", "chr2:243199373", "chr10:135534747"],
        &["chr1:249250621", "chr10:135534747", "chr2:243199373"],
    ),
    (
        "lexicographic-not-human",
        &["chr1:100", "chr10:300", "chr2:200"],
        &["chr1:100", "chr2:200", "chr10:300"],
    ),
];

fn dictionary(records: &[&str]) -> Vec<SequenceRecord> {
    records
        .iter()
        .map(|record| {
            let (name, length) = record.rsplit_once(':').expect("a name and a length");
            SequenceRecord::new(name, length.parse().expect("a length"))
        })
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t")
        .replace("\\n", "\n")
        .replace("\\\\", "\\")
}

fn field(dump: &str, kind: &str, key: &str) -> String {
    let prefix = format!("{kind}\t{key}\t");
    dump.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
        .unwrap_or_else(|| panic!("{kind}/{key} is not in the dump"))
}

#[test]
fn every_pair_compares_as_the_reference_compares_it() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `SEQUENCE_DICTIONARY_DUMP` still overrides it, which is how a harness change is checked
    // before CI sees it.
    let dump = match std::env::var("SEQUENCE_DICTIONARY_DUMP") {
        Ok(path) => {
            std::fs::read_to_string(path).expect("the dump named by SEQUENCE_DICTIONARY_DUMP")
        }
        Err(_) => gatk_corpus::read_golden(
            &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/data/sequence_dictionary_validation.txt.gz"),
        ),
    };

    let mut compared = 0;
    for (case, a, b) in CASES {
        let (a, b) = (dictionary(a), dictionary(b));
        for ordering in [false, true] {
            assert_eq!(
                compare(&a, &b, ordering).name(),
                field(&dump, "compare", &format!("{case}\t{ordering}")),
                "{case}, ordering checked: {ordering}"
            );
            compared += 1;
        }
        for superset in [false, true] {
            for ordering in [false, true] {
                let ours = match validate("reads", &a, "features", &b, superset, ordering) {
                    Ok(()) => "ok".to_string(),
                    Err(refusal) => format!("{}: {}", refusal.java_class(), refusal.message()),
                };
                assert_eq!(
                    ours,
                    field(
                        &dump,
                        "validate",
                        &format!("{case}\t{superset}\t{ordering}")
                    ),
                    "{case}, superset required: {superset}, ordering checked: {ordering}"
                );
                compared += 1;
            }
        }
    }
    assert_eq!(compared, CASES.len() * 6);

    // Every row of the dump is answered: a pair added to the harness and not here would pass.
    let rows = dump
        .lines()
        .filter(|line| line.starts_with("compare\t") || line.starts_with("validate\t"))
        .count();
    assert_eq!(rows, compared, "the dump carries a case this test does not");
    println!("rows={compared}");
}
