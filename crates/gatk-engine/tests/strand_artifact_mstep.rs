//! Conformance for `StrandArtifactFilter`'s M step and the Brent optimiser under it, against
//! GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/StrandArtifactMStepDump.java`.
//!
//! # What this suite is for
//!
//!  * **the optimiser's fifteen cases are analytic**, so they pin Brent's trajectory rather than one
//!    caller's use of it: the bounds are never reached, a flat objective answers the guess, and a
//!    NaN objective stops at its second point;
//!  * **the tolerances and the budget are refusals**, with messages whose numbers are formatted;
//!  * **a second pass with nothing accumulated does not keep the first pass's shape**, the guess
//!    being `INITIAL_ALPHA_STRAND` rather than the current alpha;
//!  * **the two mass sums are over different sets**, the artifact one over the sites above `0.1` and
//!    the non-artifact one over all of them.
//!
//! Every row is compared and every row is bit-identical, `DoubleStream.sum`'s compensation and the
//! beta binomial included.

use gatk_corpus as corpus;
use gatk_engine::strand_artifact_filter::{
    calculate_artifact_probabilities, learn_parameters, parse_strand_bias_table, EStep,
    LearnedParameters, INITIAL_ALPHA_STRAND, INITIAL_BETA_STRAND, INITIAL_STRAND_ARTIFACT_PRIOR,
};
use gatk_engine::tsv_table::java_double_to_string;
use jmath::brent::{maximize, BrentError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/strand_artifact_mstep.txt.gz"),
    )
}

fn rows() -> Vec<(String, String, String)> {
    golden()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let mut fields = line.splitn(3, '\t');
            (
                fields.next().expect("a kind").to_string(),
                fields.next().expect("a label").to_string(),
                fields.next().expect("a payload").to_string(),
            )
        })
        .collect()
}

/// Every arm of the match below hands back a different closure type, so each is boxed.
fn boxed(objective: impl Fn(f64) -> f64 + 'static) -> Box<dyn Fn(f64) -> f64> {
    Box::new(objective)
}

/// One optimiser case: its objective and the seven arguments the dump passed.
fn optimise(label: &str) -> Result<(f64, f64), BrentError> {
    let quadratic_at_three = |x: f64| -(x - 3.0) * (x - 3.0);
    let (objective, min, max, guess, relative, absolute, evaluations) = match label {
        "quadratic-at-three" => (boxed(quadratic_at_three), 0.01, 100.0, 1.0, 0.01, 0.01, 100),
        "quadratic-at-fifty" => (
            boxed(|x: f64| -(x - 50.0) * (x - 50.0)),
            0.01,
            100.0,
            1.0,
            0.01,
            0.01,
            100,
        ),
        "maximum-at-the-lower-bound" => (boxed(|x: f64| -x), 0.01, 100.0, 1.0, 0.01, 0.01, 100),
        "maximum-at-the-upper-bound" => (boxed(|x: f64| x), 0.01, 100.0, 1.0, 0.01, 0.01, 100),
        "flat" => (boxed(|_| 1.0), 0.01, 100.0, 1.0, 0.01, 0.01, 100),
        "two-maxima" => (boxed(f64::sin), 0.01, 20.0, 1.0, 0.01, 0.01, 100),
        "guess-at-the-minimum" => (
            boxed(quadratic_at_three),
            0.01,
            100.0,
            0.01,
            0.01,
            0.01,
            100,
        ),
        "guess-at-the-maximum" => (
            boxed(quadratic_at_three),
            0.01,
            100.0,
            100.0,
            0.01,
            0.01,
            100,
        ),
        "tight-tolerance" => (
            boxed(quadratic_at_three),
            0.01,
            100.0,
            1.0,
            1e-10,
            1e-10,
            1000,
        ),
        "reversed-interval" => (boxed(quadratic_at_three), 100.0, 0.01, 1.0, 0.01, 0.01, 100),
        "too-few-evaluations" => (boxed(quadratic_at_three), 0.01, 100.0, 1.0, 0.01, 0.01, 3),
        "relative-tolerance-too-small" => (
            boxed(quadratic_at_three),
            0.01,
            100.0,
            1.0,
            1e-17,
            0.01,
            100,
        ),
        "absolute-tolerance-zero" => (boxed(quadratic_at_three), 0.01, 100.0, 1.0, 0.01, 0.0, 100),
        "guess-outside-the-interval" => (
            boxed(quadratic_at_three),
            0.01,
            100.0,
            200.0,
            0.01,
            0.01,
            100,
        ),
        "nan-objective" => (boxed(|_| f64::NAN), 0.01, 100.0, 1.0, 0.01, 0.01, 100),
        other => panic!("no optimiser case named {other}"),
    };
    maximize(objective, min, max, guess, relative, absolute, evaluations)
        .map(|pair| (pair.point, pair.value))
}

/// The strand tables one learning case accumulated.
fn tables(label: &str) -> Vec<&'static str> {
    match label {
        "no-data" => vec![],
        "one-strong-artifact" | "learned-twice" => vec!["50,50|20,0"],
        "one-weak-site" => vec!["50,50|10,10"],
        "strong-and-weak" => vec!["50,50|20,0", "50,50|10,10"],
        "two-strong-artifacts" => vec!["50,50|20,0", "50,50|0,30"],
        "every-site-weak" => vec!["50,50|10,10", "50,50|9,11"],
        "deep-artifact" => vec!["2000,2000|400,0"],
        other => panic!("no learning case named {other}"),
    }
}

/// `accumulateDataForLearning` over one biallelic SNV per table, at the initial parameters.
fn accumulate(label: &str) -> Vec<EStep> {
    let mut steps = Vec::new();
    for table in tables(label) {
        let parsed = parse_strand_bias_table(table).expect("parsed");
        steps.extend(
            calculate_artifact_probabilities(
                &parsed,
                &[0],
                INITIAL_STRAND_ARTIFACT_PRIOR,
                INITIAL_ALPHA_STRAND,
                INITIAL_BETA_STRAND,
            )
            .expect("answered"),
        );
    }
    steps
}

fn render(learned: &LearnedParameters) -> String {
    format!(
        "{},{},{}",
        java_double_to_string(learned.strand_artifact_prior),
        java_double_to_string(learned.alpha_strand),
        java_double_to_string(learned.beta_strand)
    )
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 32, "the golden's row count");
    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "opt" => {
                let (_, expected) = payload.rsplit_once('=').expect("an answer");
                let (point, value) = optimise(label).expect("optimised");
                assert_eq!(
                    format!(
                        "{},{}",
                        java_double_to_string(point),
                        java_double_to_string(value)
                    ),
                    expected,
                    "opt {label}"
                );
            }
            "accumulated" => {
                // `learned-twice-after` is the list after a pass, which `learnParameters` clears.
                let expected: usize = payload.parse().expect("a count");
                let ours = if label.ends_with("-after") {
                    0
                } else {
                    accumulate(label).len()
                };
                assert_eq!(ours, expected, "accumulated {label}");
            }
            "learned" => {
                let steps = if label.ends_with("-second") {
                    // The second pass starts from an empty list.
                    Vec::new()
                } else {
                    accumulate(label.trim_end_matches("-first"))
                };
                assert_eq!(
                    render(&learn_parameters(&steps)),
                    *payload,
                    "learned {label}"
                );
            }
            "error" => {
                let error = optimise(label).expect_err("refused");
                assert_eq!(
                    *payload,
                    format!("{}:{}", error.class(), error.message()),
                    "error {label}"
                );
            }
            other => panic!("no row kind {other}"),
        }
    }
}
