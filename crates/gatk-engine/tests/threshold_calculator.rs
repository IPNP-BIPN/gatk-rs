//! Conformance for `ThresholdCalculator` against GATK 4.6.2.0, compared as the threshold every
//! strategy answers for every list, the list the sort left behind, and both refusals.
//!
//! Golden from `tools/readfilter-conformance/ThresholdCalculatorDump.java`.
//!
//! # What this suite is for
//!
//!  * **the optimal F score keeps the last tie**, so four identical posteriors answer `1.0`;
//!  * **its answer is three-way**, and only the middle case reports a posterior;
//!  * **the false-discovery walk steps back one**, and answers `0.0` or `1.0` at the two ends;
//!  * **an empty list is answered at opposite ends by the two strategies**, which is what a second
//!    relearn computes from;
//!  * **and the list is sorted in place**.

use gatk_corpus as corpus;
use gatk_engine::threshold_calculator::{
    threshold_from_false_discovery_rate, threshold_from_optimal_f_score, Strategy,
    ThresholdCalculator, ThresholdError,
};
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/threshold_calculator.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

/// The threshold the reference wrote for one label.
fn expected(text: &str, label: &str) -> String {
    rows(text, "threshold")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no threshold {label}"))[1]
        .to_string()
}

/// The six lists the dump drives every strategy over.
fn list(name: &str) -> Vec<f64> {
    match name {
        "spread" => vec![0.9, 0.1, 0.5, 0.02, 0.3],
        "tied" => vec![0.2, 0.2, 0.2, 0.2],
        "hopeless" => vec![0.99, 0.98, 0.97],
        "clean" => vec![0.001, 0.002, 0.003],
        "single" => vec![0.4],
        "empty" => Vec::new(),
        other => panic!("no list {other}"),
    }
}

/// One run of the dump: a fresh calculator, one list, one relearn.
fn ours(strategy: Strategy, name: &str, max_false_discovery_rate: f64, beta: f64) -> String {
    let mut calculator = ThresholdCalculator::new(strategy, 0.123, max_false_discovery_rate, beta);
    calculator.add(&list(name));
    calculator.relearn().expect("valid parameters");
    java_double_to_string(calculator.threshold())
}

#[test]
fn every_threshold_matches_the_golden() {
    let text = golden();
    for strategy in [
        Strategy::Constant,
        Strategy::FalseDiscoveryRate,
        Strategy::OptimalFScore,
    ] {
        for name in ["spread", "tied", "hopeless", "clean", "single", "empty"] {
            let label = format!("{}-{name}", strategy.name());
            assert_eq!(
                ours(strategy, name, 0.05, 1.0),
                expected(&text, &label),
                "{label}"
            );
        }
    }
    // The two rates and the two betas.
    assert_eq!(
        ours(Strategy::FalseDiscoveryRate, "spread", 1.0, 1.0),
        expected(&text, "FALSE_DISCOVERY_RATE-loose")
    );
    assert_eq!(
        ours(Strategy::FalseDiscoveryRate, "spread", 0.001, 1.0),
        expected(&text, "FALSE_DISCOVERY_RATE-tight")
    );
    assert_eq!(
        ours(Strategy::OptimalFScore, "spread", 0.05, 0.0),
        expected(&text, "OPTIMAL_F_SCORE-beta-zero")
    );
    assert_eq!(
        ours(Strategy::OptimalFScore, "spread", 0.05, 10.0),
        expected(&text, "OPTIMAL_F_SCORE-beta-ten")
    );
}

#[test]
fn the_optimal_f_score_keeps_the_last_tie() {
    let text = golden();
    // Four identical posteriors of 0.2, and the answer is neither 0.2 nor 0.0.
    assert_eq!(expected(&text, "OPTIMAL_F_SCORE-tied"), "1.0");
    let mut tied = list("tied");
    assert_eq!(
        threshold_from_optimal_f_score(&mut tied, 1.0).expect("beta"),
        1.0
    );
    // The same list under the other strategy is at the far end.
    assert_eq!(expected(&text, "FALSE_DISCOVERY_RATE-tied"), "0.0");
}

#[test]
fn the_optimal_f_score_answers_three_ways() {
    let text = golden();
    // The last index, by two different routes.
    assert_eq!(expected(&text, "OPTIMAL_F_SCORE-clean"), "1.0");
    assert_eq!(expected(&text, "OPTIMAL_F_SCORE-single"), "1.0");
    // A posterior of its own.
    assert_eq!(expected(&text, "OPTIMAL_F_SCORE-hopeless"), "0.97");
    let mut hopeless = list("hopeless");
    assert_eq!(
        threshold_from_optimal_f_score(&mut hopeless, 1.0).expect("beta"),
        0.97
    );
    // And no index at all.
    assert_eq!(expected(&text, "OPTIMAL_F_SCORE-empty"), "0.0");
}

#[test]
fn a_second_relearn_computes_from_an_empty_list() {
    let text = golden();
    for (strategy, name) in [
        (Strategy::FalseDiscoveryRate, "FALSE_DISCOVERY_RATE"),
        (Strategy::OptimalFScore, "OPTIMAL_F_SCORE"),
        (Strategy::Constant, "CONSTANT"),
    ] {
        let mut calculator = ThresholdCalculator::new(strategy, 0.123, 0.05, 1.0);
        calculator.add(&list("spread"));
        calculator.relearn().expect("valid");
        assert_eq!(
            java_double_to_string(calculator.threshold()),
            expected(&text, &format!("{name}-first")),
            "{name} first"
        );
        assert!(calculator.accumulated().is_empty());
        calculator.relearn().expect("nothing left");
        assert_eq!(
            java_double_to_string(calculator.threshold()),
            expected(&text, &format!("{name}-second")),
            "{name} second"
        );
    }
    // The two strategies answer the emptied list at opposite ends.
    assert_eq!(expected(&text, "FALSE_DISCOVERY_RATE-second"), "1.0");
    assert_eq!(expected(&text, "OPTIMAL_F_SCORE-second"), "0.0");
}

#[test]
fn the_list_is_sorted_in_place() {
    let text = golden();
    let before = rows(&text, "sorted")
        .into_iter()
        .find(|row| row[0] == "before")
        .expect("before")[1]
        .to_string();
    let after = rows(&text, "sorted")
        .into_iter()
        .find(|row| row[0] == "after")
        .expect("after")[1]
        .to_string();
    assert_eq!(before, "[0.9, 0.1, 0.5]");
    assert_eq!(after, "[0.1, 0.5, 0.9]");

    let mut ours = vec![0.9, 0.1, 0.5];
    threshold_from_false_discovery_rate(&mut ours, 0.05).expect("a rate");
    assert_eq!(
        format!(
            "[{}]",
            ours.iter()
                .map(|value| java_double_to_string(*value))
                .collect::<Vec<_>>()
                .join(", ")
        ),
        after
    );
}

#[test]
fn both_refusals_carry_the_references_class_and_words() {
    let text = golden();
    let refusal = |label: &str| -> (String, String) {
        let row = rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no refusal {label}"));
        let (class, message) = row[1].split_once(':').expect("class and message");
        (class.to_string(), message.to_string())
    };

    let (class, message) = refusal("negative-beta");
    let error = threshold_from_optimal_f_score(&mut list("spread"), -1.0).expect_err("negative");
    assert_eq!(error.class(), class);
    assert_eq!(error.message(), message);
    assert_eq!(error, ThresholdError::NegativeBeta);

    let (class, message) = refusal("negative-rate");
    let error =
        threshold_from_false_discovery_rate(&mut list("spread"), -0.5).expect_err("negative");
    assert_eq!(error.class(), class);
    assert_eq!(error.message(), message);
}
