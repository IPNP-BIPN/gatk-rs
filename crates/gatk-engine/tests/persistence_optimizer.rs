//! Conformance for `PersistenceOptimizer` against GATK 4.6.2.0, compared as the indices it returns
//! and the persistences beside them.
//!
//! Golden from `tools/readfilter-conformance/PersistenceOptimizerDump.java`.
//!
//! # What this suite is for
//!
//!  * **the ordering is `Double.compare`**, so `-0.0` sorts below `0.0` and a `NaN` above
//!    everything, which moves the global minimum and can make a persistence negative zero;
//!  * **the sort is stable**, so a plateau's minimum is its leftmost point;
//!  * **the global minimum is prepended**, and its persistence is the whole range;
//!  * **and every data set is compared as text**, so a spelling that differs is a failure even
//!    where the value would compare equal.

use gatk_corpus as corpus;
use gatk_engine::persistence_optimizer::persistence_optimizer;
use gatk_engine::tsv_table::java_double_to_string as java;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/persistence_optimizer.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn labelled(text: &str, kind: &str, label: &str) -> String {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no {kind} row for {label}"))
        .get(1)
        .copied()
        .unwrap_or("")
        .to_string()
}

/// `Double.parseDouble` for the spellings the dump prints.
fn parse(value: &str) -> f64 {
    match value {
        "NaN" => f64::NAN,
        "Infinity" => f64::INFINITY,
        "-Infinity" => f64::NEG_INFINITY,
        other => other.parse().expect("a double"),
    }
}

/// The data of one label, taken from the golden itself, so the port is fed exactly what the
/// reference was fed, including the negative zeroes.
fn data(text: &str, label: &str) -> Vec<f64> {
    let row = labelled(text, "data", label);
    if row.is_empty() {
        return Vec::new();
    }
    row.split(',').map(parse).collect()
}

fn labels(text: &str) -> Vec<String> {
    rows(text, "data")
        .into_iter()
        .map(|row| row[0].to_string())
        .collect()
}

#[test]
fn every_set_of_minima_is_the_reference() {
    let text = golden();
    for label in labels(&text) {
        let answer = persistence_optimizer(&data(&text, &label)).expect("the data is accepted");
        let ours: Vec<String> = answer
            .minima_indices
            .iter()
            .map(|index| index.to_string())
            .collect();
        assert_eq!(
            ours.join(","),
            labelled(&text, "minima", &label),
            "minima/{label}"
        );
    }
}

/// The persistences as text, which is how a negative zero and a NaN are told apart from a zero and
/// a number that merely compares equal.
#[test]
fn every_persistence_is_the_reference_to_the_digit() {
    let text = golden();
    for label in labels(&text) {
        let answer = persistence_optimizer(&data(&text, &label)).expect("the data is accepted");
        let ours: Vec<String> = answer
            .persistences
            .iter()
            .map(|value| java(*value))
            .collect();
        assert_eq!(
            ours.join(","),
            labelled(&text, "persistence", &label),
            "persistence/{label}"
        );
    }
}

/// The two orderings that are not `<`, pinned on their own.
#[test]
fn the_signed_zero_and_the_nan_are_where_the_reference_puts_them() {
    let text = golden();

    assert_eq!(labelled(&text, "minima", "signed-zeroes"), "1,2");
    let signed = persistence_optimizer(&data(&text, "signed-zeroes")).expect("data");
    assert!(
        signed.persistences[1].is_sign_negative() && signed.persistences[1] == 0.0,
        "the second persistence is a negative zero"
    );

    let with_nan = persistence_optimizer(&data(&text, "with-nan")).expect("data");
    assert!(with_nan.persistences[0].is_nan(), "the NaN is the maximum");
    assert_eq!(labelled(&text, "minima", "with-nan"), "2,0,4");
}

#[test]
fn empty_data_is_the_references_refusal() {
    let text = golden();
    let error = persistence_optimizer(&[]).unwrap_err();
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        labelled(&text, "error", "empty")
    );
}
