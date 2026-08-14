//! Conformance for `ErrorProbabilities` against GATK 4.6.2.0, compared as the combined probability
//! of every record and the number of filters that survived.
//!
//! Golden from `tools/readfilter-conformance/MutectErrorProbabilitiesDump.java`.
//!
//! # What this suite is for
//!
//!  * **a filter that answers an empty list is dropped**, so four filters leave three;
//!  * **one error type is a maximum**, so failing two filters is no worse than failing one;
//!  * **ragged lists are a refusal**;
//!  * **and a site-level filter's one answer is copied to every allele**.
//!
//! Every probability the dump can reach is 0 or 1, because each filter that answers a fraction needs
//! the somatic clustering model or a contamination table. What is compared here is therefore the
//! shape of the combination; its arithmetic is asserted in the module's own unit tests, on values
//! the golden cannot produce.

use gatk_corpus as corpus;
use gatk_engine::error_probabilities::{combined, kept, ErrorType, FilterAnswer};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mutect_error_probabilities.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn expected(text: &str, kind: &str, label: &str) -> String {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no {kind} {label}"))[1]
        .to_string()
}

/// A list of probabilities, printed the way Java prints one.
fn list(values: &[f64]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| gatk_engine::tsv_table::java_double_to_string(*value))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn artifact(probabilities: &[f64]) -> FilterAnswer {
    FilterAnswer {
        error_type: ErrorType::Artifact,
        probabilities: probabilities.to_vec(),
    }
}

/// The four filters' answers for one of the dump's records: base quality, read position, clustered
/// events, and the strand filter that is switched off and therefore empty.
fn answers(record: &str) -> Vec<FilterAnswer> {
    match record {
        // MBQ 30,10 fails base quality; MPOS 20 passes; ECNT 1 passes.
        "one-allele-one-failure" => vec![
            artifact(&[1.0]),
            artifact(&[0.0]),
            artifact(&[0.0]),
            artifact(&[]),
        ],
        // The same allele failing two of the three.
        "one-allele-two-failures" => vec![
            artifact(&[1.0]),
            artifact(&[1.0]),
            artifact(&[0.0]),
            artifact(&[]),
        ],
        "one-allele-no-failure" => vec![
            artifact(&[0.0]),
            artifact(&[0.0]),
            artifact(&[0.0]),
            artifact(&[]),
        ],
        // Two alternates: the first fails both allele filters, the second neither.
        "two-alleles" => vec![
            artifact(&[1.0, 0.0]),
            artifact(&[1.0, 0.0]),
            artifact(&[0.0, 0.0]),
            artifact(&[]),
        ],
        // A site-level filter firing: its one answer copied to both alleles.
        "two-alleles-site-filter" => vec![
            artifact(&[0.0, 0.0]),
            artifact(&[0.0, 0.0]),
            artifact(&[1.0, 1.0]),
            artifact(&[]),
        ],
        // The symbolic alternate has been removed, leaving one allele.
        "symbolic-allele" => vec![
            artifact(&[1.0]),
            artifact(&[1.0]),
            artifact(&[0.0]),
            artifact(&[]),
        ],
        other => panic!("no record {other}"),
    }
}

const RECORDS: [&str; 6] = [
    "one-allele-one-failure",
    "one-allele-two-failures",
    "one-allele-no-failure",
    "two-alleles",
    "two-alleles-site-filter",
    "symbolic-allele",
];

#[test]
fn every_combined_probability_matches_the_golden() {
    let text = golden();
    for record in RECORDS {
        let answers = answers(record);
        assert_eq!(
            list(&combined(&answers).expect("equal lengths")),
            expected(&text, "combined", record),
            "{record}"
        );
    }
}

#[test]
fn a_filter_that_answers_an_empty_list_is_dropped() {
    let text = golden();
    // Four filters offered, three surviving, in every record of the dump.
    for record in RECORDS {
        assert_eq!(expected(&text, "filters", record), "3", "{record}");
        assert_eq!(kept(&answers(record)).len(), 3, "{record}");
    }
}

#[test]
fn one_error_type_is_a_maximum() {
    let text = golden();
    // Failing two filters is the same answer as failing one.
    assert_eq!(
        expected(&text, "combined", "one-allele-one-failure"),
        expected(&text, "combined", "one-allele-two-failures")
    );
    assert_eq!(
        expected(&text, "combined", "one-allele-no-failure"),
        "[0.0]"
    );
}

#[test]
fn a_site_filter_copies_its_answer_to_every_allele() {
    let text = golden();
    // One allele fails the allele filters, the other does not.
    assert_eq!(expected(&text, "combined", "two-alleles"), "[1.0, 0.0]");
    // The site filter reaches both.
    assert_eq!(
        expected(&text, "combined", "two-alleles-site-filter"),
        "[1.0, 1.0]"
    );
}

#[test]
fn a_symbolic_alternate_is_removed_before_the_combination() {
    let text = golden();
    // The record has two alternates and one of them is symbolic: one probability comes back.
    assert_eq!(expected(&text, "combined", "symbolic-allele"), "[1.0]");
    assert_eq!(
        combined(&answers("symbolic-allele"))
            .expect("equal lengths")
            .len(),
        1
    );
}

#[test]
fn ragged_lists_are_a_refusal() {
    let text = golden();
    let row = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "ragged-lists")
        .expect("a refusal");
    let (class, message) = row[1].split_once(':').expect("class and message");
    // Two filters answering two alleles and one answering one.
    let error = combined(&[artifact(&[1.0, 0.0]), artifact(&[1.0])]).expect_err("two and one");
    assert_eq!(error.class(), class);
    assert_eq!(error.message(), message);
}
