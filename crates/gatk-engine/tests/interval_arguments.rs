//! Conformance for `IntervalArgumentCollection` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/IntervalArgumentsDump.java`. The runners read
//! `--intervals` and ignored the other four arguments of the collection every walker carries, so
//! `--interval-padding`, `--interval-exclusion-padding`, `--interval-set-rule` and
//! `--exclude-intervals` changed nothing in the port and changed the answer in the reference.
//!
//! # What this suite is for
//!
//!  * **padding being applied per `-L` argument**, before the set operator, and each padded batch
//!    being sorted and merged with `ALL` whatever the merging rule is;
//!  * **the fold short-circuiting on an empty side**, so the first `-L` is never intersected with
//!    anything and `INTERSECTION` over three arguments is `((a ∩ b) ∩ c)`;
//!  * **an empty intersection being a refusal**, and an exclusion that removes everything being
//!    another one, each quoting the raw argument strings;
//!  * **`-XL` with no `-L` meaning the whole reference**, contig by contig, in the dictionary's
//!    own order;
//!  * **padding being clamped to the contig at both ends**;
//!  * **and `unmapped` being a traversal flag rather than an interval**, accepted on `-L` and
//!    refused on `-XL`.
//!
//! The golden is committed and re-derived by the `interval-arguments` suite on every run; the
//! dump can still be overridden with an environment variable while a harness change is being
//! checked.

use gatk_engine::interval::MergingRule;
use gatk_engine::interval_arguments::{traversal_parameters, SetRule};
use htsjdk_bam::header::{SamHeader, SequenceRecord};

/// The dictionary every case resolves against, whose order is not alphabetical.
fn dictionary() -> SamHeader {
    let mut header = SamHeader::default();
    for (name, length) in [("chr1", 1000), ("chr2", 500), ("chr10", 200)] {
        header.sequences.push(SequenceRecord::new(name, length));
    }
    header
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_string()).collect()
}

/// One command line: include, exclude, set rule, merging rule, padding, exclusion padding.
type Case = (
    &'static str,
    &'static [&'static str],
    &'static [&'static str],
    SetRule,
    MergingRule,
    i32,
    i32,
);

const UNION: SetRule = SetRule::Union;
const INTERSECTION: SetRule = SetRule::Intersection;
const ALL: MergingRule = MergingRule::All;
const OVERLAPPING: MergingRule = MergingRule::OverlappingOnly;

const CASES: &[Case] = &[
    ("one-interval", &["chr1:100-200"], &[], UNION, ALL, 0, 0),
    ("padded", &["chr1:100-200"], &[], UNION, ALL, 50, 0),
    (
        "padded-past-the-start",
        &["chr1:10-20"],
        &[],
        UNION,
        ALL,
        50,
        0,
    ),
    (
        "padded-past-the-end",
        &["chr2:480-490"],
        &[],
        UNION,
        ALL,
        50,
        0,
    ),
    (
        "padding-merges-two",
        &["chr1:100-200", "chr1:260-300"],
        &[],
        UNION,
        ALL,
        30,
        0,
    ),
    (
        "adjacent-all",
        &["chr1:100-200", "chr1:201-300"],
        &[],
        UNION,
        ALL,
        0,
        0,
    ),
    (
        "adjacent-overlapping-only",
        &["chr1:100-200", "chr1:201-300"],
        &[],
        UNION,
        OVERLAPPING,
        0,
        0,
    ),
    (
        "union-two",
        &["chr1:100-200", "chr1:150-300"],
        &[],
        UNION,
        ALL,
        0,
        0,
    ),
    (
        "intersection-two",
        &["chr1:100-200", "chr1:150-300"],
        &[],
        INTERSECTION,
        ALL,
        0,
        0,
    ),
    (
        "intersection-three",
        &["chr1:100-400", "chr1:200-500", "chr1:300-600"],
        &[],
        INTERSECTION,
        ALL,
        0,
        0,
    ),
    (
        "intersection-empty",
        &["chr1:100-200", "chr1:300-400"],
        &[],
        INTERSECTION,
        ALL,
        0,
        0,
    ),
    (
        "intersection-across-contigs",
        &["chr1:100-200", "chr2:100-200"],
        &[],
        INTERSECTION,
        ALL,
        0,
        0,
    ),
    (
        "exclude-middle",
        &["chr1:100-300"],
        &["chr1:150-200"],
        UNION,
        ALL,
        0,
        0,
    ),
    (
        "exclude-prefix",
        &["chr1:100-300"],
        &["chr1:50-150"],
        UNION,
        ALL,
        0,
        0,
    ),
    (
        "exclude-everything",
        &["chr1:100-200"],
        &["chr1:1-1000"],
        UNION,
        ALL,
        0,
        0,
    ),
    (
        "exclusion-padded",
        &["chr1:100-300"],
        &["chr1:180-190"],
        UNION,
        ALL,
        0,
        20,
    ),
    (
        "both-paddings",
        &["chr1:100-300"],
        &["chr1:180-190"],
        UNION,
        ALL,
        10,
        20,
    ),
    ("exclude-only", &[], &["chr1:1-900"], UNION, ALL, 0, 0),
    (
        "exclude-only-whole-contig",
        &[],
        &["chr2"],
        UNION,
        ALL,
        0,
        0,
    ),
    (
        "unmapped-included",
        &["chr1:100-200", "unmapped"],
        &[],
        UNION,
        ALL,
        0,
        0,
    ),
    ("unmapped-only", &["unmapped"], &[], UNION, ALL, 0, 0),
    (
        "unmapped-excluded",
        &["chr1:100-200"],
        &["unmapped"],
        UNION,
        ALL,
        0,
        0,
    ),
    ("whole-contig", &["chr2"], &[], UNION, ALL, 0, 0),
    (
        "three-contigs-out-of-order",
        &["chr10:1-10", "chr2:1-10", "chr1:1-10"],
        &[],
        UNION,
        ALL,
        0,
        0,
    ),
];

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t")
        .replace("\\n", "\n")
        .replace("\\\\", "\\")
}

fn field(dump: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    dump.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

#[test]
fn every_command_line_resolves_as_the_reference_resolves_it() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `INTERVAL_ARGUMENTS_DUMP` still overrides it, which is how a harness change is checked
    // before CI sees it.
    let dump = match std::env::var("INTERVAL_ARGUMENTS_DUMP") {
        Ok(path) => {
            std::fs::read_to_string(path).expect("the dump named by INTERVAL_ARGUMENTS_DUMP")
        }
        Err(_) => gatk_corpus::read_golden(
            &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/data/interval_arguments.txt.gz"),
        ),
    };
    let header = dictionary();

    for (case, include, exclude, set_rule, merging_rule, padding, exclusion_padding) in CASES {
        let answer = traversal_parameters(
            &strings(include),
            &strings(exclude),
            &header,
            *set_rule,
            *merging_rule,
            *padding,
            *exclusion_padding,
        );
        match answer {
            Ok(parameters) => {
                let rendered = if parameters.intervals.is_empty() {
                    "(empty)".to_string()
                } else {
                    parameters
                        .intervals
                        .iter()
                        .map(|interval| {
                            format!("{}:{}-{}", interval.contig, interval.start, interval.end)
                        })
                        .collect::<Vec<_>>()
                        .join(" ")
                };
                assert_eq!(
                    rendered,
                    field(&dump, "intervals", case)
                        .unwrap_or_else(|| panic!("{case}: the golden refused, the port did not")),
                    "{case}"
                );
                assert_eq!(
                    parameters.traverse_unmapped.to_string(),
                    field(&dump, "unmapped", case).expect("the unmapped row"),
                    "{case}: unmapped"
                );
            }
            Err(error) => {
                assert_eq!(
                    format!("{}: {}", error.java_class(), error.message()),
                    field(&dump, "error", case)
                        .unwrap_or_else(|| panic!("{case}: the port refused, the golden did not")),
                    "{case}"
                );
            }
        }
    }

    // Every case in the dump is answered here.
    let cases: std::collections::BTreeSet<&str> = dump
        .lines()
        .filter_map(|line| line.split('\t').nth(1))
        .collect();
    assert_eq!(
        cases.len(),
        CASES.len(),
        "the dump carries a case this test does not"
    );
}
