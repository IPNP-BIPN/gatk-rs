//! Conformance for Mutect's engine-free hard filters against GATK 4.6.2.0, compared as the answer
//! each filter gives each record, and as all four refusals.
//!
//! Golden from `tools/readfilter-conformance/MutectHardFiltersDump.java`.
//!
//! # What this suite is for
//!
//!  * **only a long insertion is judged by the reference's mapping quality**, a deletion's indel
//!    length being negative;
//!  * **a negative median read position is never an artifact**;
//!  * **strict strand bias answers an empty list** when it is switched off;
//!  * **the fragment-length filter looks at one allele only**;
//!  * **and four filters break four ways on a record with no annotations**.

use gatk_corpus as corpus;
use gatk_engine::mutect_hard_filters::{
    base_quality_artifacts, clustered_events_is_artifact, fragment_length_is_artifact,
    mapping_quality_artifacts, multiallelic_is_artifact, read_position_artifacts,
    strict_strand_artifacts,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mutect_hard_filters.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

/// The answer the reference gave for one filter and one record.
fn answer(text: &str, label: &str) -> String {
    rows(text, "filter")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no answer {label}"))[1]
        .to_string()
}

/// The refusal the reference gave, as its class and message.
fn refusal(text: &str, label: &str) -> (String, String) {
    let row = rows(text, "error")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no refusal {label}"));
    let (class, message) = row[1].split_once(':').expect("class and message");
    (class.to_string(), message.replace("\\\"", "\""))
}

/// A list of booleans, printed the way Java prints one.
fn list(values: &[bool]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    )
}

/// The annotations of one of the dump's records.
struct Record {
    median_base_qualities: Vec<i32>,
    median_mapping_qualities: Vec<i32>,
    median_fragment_lengths: Vec<i32>,
    median_read_positions: Vec<i32>,
    haplotype_event_counts: Vec<i32>,
    region_event_count: i32,
    tumour_log_odds: Vec<f64>,
    strand_counts: Vec<Vec<i32>>,
    /// `getIndelLengths()`: the alt length minus the ref length, or nothing for a SNP.
    indel_lengths: Option<Vec<i32>>,
}

fn record(name: &str) -> Record {
    let common = |mbq: Vec<i32>,
                  mmq: Vec<i32>,
                  mfrl: Vec<i32>,
                  mpos: Vec<i32>,
                  ecnth: Vec<i32>,
                  ecnt: i32,
                  tlod: Vec<f64>,
                  sb: Vec<Vec<i32>>,
                  indels: Option<Vec<i32>>| Record {
        median_base_qualities: mbq,
        median_mapping_qualities: mmq,
        median_fragment_lengths: mfrl,
        median_read_positions: mpos,
        haplotype_event_counts: ecnth,
        region_event_count: ecnt,
        tumour_log_odds: tlod,
        strand_counts: sb,
        indel_lengths: indels,
    };
    match name {
        "poor-snp" => common(
            vec![30, 10],
            vec![60, 20],
            vec![300, 380],
            vec![2],
            vec![1],
            1,
            vec![6.0],
            vec![vec![5, 5], vec![0, 7]],
            None,
        ),
        "good-snp" => common(
            vec![30, 35],
            vec![60, 60],
            vec![300, 305],
            vec![20],
            vec![1],
            1,
            vec![6.0],
            vec![vec![5, 5], vec![4, 3]],
            None,
        ),
        "triallelic" => common(
            vec![30, 10, 35],
            vec![60, 20, 60],
            vec![300, 380, 302],
            vec![2, 20],
            vec![1, 4],
            3,
            vec![6.0, 4.0],
            vec![vec![5, 5], vec![0, 7], vec![3, 4]],
            None,
        ),
        // A seven-base deletion: ALT length minus REF length is -7.
        "long-deletion" => common(
            vec![30, 35],
            vec![60, 20],
            vec![300, 305],
            vec![20],
            vec![1],
            1,
            vec![6.0],
            vec![vec![5, 5], vec![4, 3]],
            Some(vec![-7]),
        ),
        // The same annotations on a seven-base insertion.
        "long-insertion" => common(
            vec![30, 35],
            vec![60, 20],
            vec![300, 305],
            vec![20],
            vec![1],
            1,
            vec![6.0],
            vec![vec![5, 5], vec![4, 3]],
            Some(vec![7]),
        ),
        "negative-position" => common(
            vec![30, 35],
            vec![60, 60],
            vec![300, 305],
            vec![-1],
            vec![1],
            1,
            vec![6.0],
            vec![vec![5, 5], vec![4, 3]],
            None,
        ),
        other => panic!("no record {other}"),
    }
}

const RECORDS: [&str; 6] = [
    "poor-snp",
    "good-snp",
    "triallelic",
    "long-deletion",
    "long-insertion",
    "negative-position",
];

#[test]
fn every_answer_matches_the_golden() {
    let text = golden();
    for name in RECORDS {
        let record = record(name);
        assert_eq!(
            list(&base_quality_artifacts(&record.median_base_qualities, 20.0)),
            answer(&text, &format!("base-quality-{name}")),
            "base-quality-{name}"
        );
        assert_eq!(
            list(
                &mapping_quality_artifacts(
                    &record.median_mapping_qualities,
                    record.indel_lengths.as_deref(),
                    30.0,
                    5
                )
                .expect("a ref entry")
            ),
            answer(&text, &format!("mapping-quality-{name}")),
            "mapping-quality-{name}"
        );
        assert_eq!(
            list(&read_position_artifacts(&record.median_read_positions, 5.0)),
            answer(&text, &format!("read-position-{name}")),
            "read-position-{name}"
        );
        assert_eq!(
            list(&strict_strand_artifacts(&record.strand_counts, 1)),
            answer(&text, &format!("strict-strand-{name}")),
            "strict-strand-{name}"
        );
        assert_eq!(
            fragment_length_is_artifact(&record.median_fragment_lengths, 50.0)
                .expect("two entries")
                .to_string(),
            answer(&text, &format!("fragment-length-{name}")),
            "fragment-length-{name}"
        );
        assert_eq!(
            clustered_events_is_artifact(
                &record.haplotype_event_counts,
                record.region_event_count,
                2,
                2
            )
            .expect("a maximum")
            .to_string(),
            answer(&text, &format!("clustered-events-{name}")),
            "clustered-events-{name}"
        );
        assert_eq!(
            multiallelic_is_artifact(Some(&record.tumour_log_odds), 1)
                .expect("an array")
                .to_string(),
            answer(&text, &format!("multiallelic-{name}")),
            "multiallelic-{name}"
        );
    }
}

#[test]
fn only_a_long_insertion_is_judged_by_the_reference() {
    let text = golden();
    // The same MMQ of 60,20 on both, and two different answers.
    assert_eq!(answer(&text, "mapping-quality-long-deletion"), "[true]");
    assert_eq!(answer(&text, "mapping-quality-long-insertion"), "[false]");
    assert_eq!(
        record("long-deletion").median_mapping_qualities,
        record("long-insertion").median_mapping_qualities
    );
}

#[test]
fn a_negative_median_read_position_is_never_an_artifact() {
    let text = golden();
    // A position of 2 is under the minimum and is an artifact; -1 is under it and is not.
    assert_eq!(answer(&text, "read-position-poor-snp"), "[true]");
    assert_eq!(answer(&text, "read-position-negative-position"), "[false]");
}

#[test]
fn strict_strand_bias_answers_an_empty_list_when_it_is_switched_off() {
    let text = golden();
    assert_eq!(answer(&text, "strict-strand-poor-snp"), "[true]");
    assert_eq!(answer(&text, "strict-strand-off-poor-snp"), "[]");
    assert_eq!(
        list(&strict_strand_artifacts(
            &record("poor-snp").strand_counts,
            0
        )),
        "[]"
    );
}

#[test]
fn the_annotation_list_is_the_same_before_and_after() {
    let text = golden();
    assert_eq!(answer(&text, "mutation-before"), "[60, 20]");
    assert_eq!(answer(&text, "mutation-after"), "[60, 20]");
    // The port takes a slice and copies, so the caller's list cannot be touched at all.
    let qualities = [60, 20];
    mapping_quality_artifacts(&qualities, None, 30.0, 5).expect("a ref entry");
    assert_eq!(qualities, [60, 20]);
}

#[test]
fn four_filters_break_four_ways_on_a_record_with_no_annotations() {
    let text = golden();
    for (label, error) in [
        (
            "mapping-quality-no-annotations",
            mapping_quality_artifacts(&[], None, 30.0, 5).expect_err("nothing to remove"),
        ),
        (
            "fragment-length-no-annotations",
            fragment_length_is_artifact(&[], 50.0).expect_err("nothing to get"),
        ),
        (
            "clustered-events-no-annotations",
            clustered_events_is_artifact(&[], 0, 2, 2).expect_err("no maximum"),
        ),
        (
            "multiallelic-no-annotations",
            multiallelic_is_artifact(None, 1).expect_err("no array"),
        ),
    ] {
        let (class, message) = refusal(&text, label);
        assert_eq!(error.class(), class, "{label}");
        assert_eq!(error.message(), message, "{label}");
    }
    // And the three that answer an empty list instead of refusing.
    assert_eq!(answer(&text, "base-quality-no-annotations"), "[]");
    assert_eq!(answer(&text, "read-position-no-annotations"), "[]");
    assert_eq!(answer(&text, "strict-strand-no-annotations"), "[]");
}
