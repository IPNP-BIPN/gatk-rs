//! Conformance for `Mutect2FilteringEngine`'s static arithmetic against GATK 4.6.2.0, compared as
//! every posterior of the grid, every clamp, and every conversion.
//!
//! Golden from `tools/readfilter-conformance/MutectEngineArithmeticDump.java`.
//!
//! # What this suite is for
//!
//!  * **the posterior is the probability of error**, and falls as the odds of being real rise;
//!  * **the two infinities do not agree**: one answers NaN, the other is refused;
//!  * **rounding the error is a clamp that keeps NaN**;
//!  * **and the tumour odds are converted from log10**.
//!
//! # One of the forty-eight is not bit-identical
//!
//! The posterior puts every entry through `exp`, whose bit-exact transcription is **withdrawn**
//! under htsjdk-rs decision 0014: the HotSpot source it would have to be transcribed from is GPL2
//! with no Classpath Exception. Forty-seven of the grid's forty-eight inputs agree to the bit with
//! the reference; `100.0,-1.0` differs by one ulp, `6.392138950083686E-44` against
//! `...687E-44`. The test asserts the bound and names the input, so the divergence cannot spread
//! unnoticed to a second one.

use gatk_corpus as corpus;
use gatk_engine::mutect_engine::{
    posterior_probability_of_error, round_finite_precision_errors, tumor_log_odds, EPSILON,
    MIN_REPORTABLE_ERROR_PROBABILITY,
};
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mutect_engine_arithmetic.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

/// The label the dump prints a double under, which is `Double.toString`.
fn label(value: f64) -> String {
    java_double_to_string(value)
}

const LOG_ODDS: [f64; 8] = [
    -10.0,
    -1.0,
    0.0,
    1.0,
    10.0,
    100.0,
    f64::NEG_INFINITY,
    f64::INFINITY,
];
const LOG_PRIORS: [f64; 6] = [0.0, -0.001, -1.0, -10.0, -100.0, f64::NEG_INFINITY];

#[test]
fn every_posterior_matches_the_golden() {
    let text = golden();
    let answers = rows(&text, "posterior");
    let refusals = rows(&text, "error");
    let mut compared = 0;
    let mut within_one_ulp: Vec<String> = Vec::new();
    for odds in LOG_ODDS {
        for prior in LOG_PRIORS {
            let key = format!("{},{}", label(odds), label(prior));
            match posterior_probability_of_error(odds, prior) {
                Ok(value) => {
                    let expected = answers
                        .iter()
                        .find(|row| row[0] == key)
                        .unwrap_or_else(|| panic!("no posterior {key}"))[1];
                    if java_double_to_string(value) != expected {
                        // Not bit-identical, and cannot be: every entry of the normalisation goes
                        // through `exp`, whose bit-exact transcription is withdrawn under
                        // htsjdk-rs decision 0014. What is asserted is the ulp bound the port
                        // claims, and the divergence is named so it cannot grow unnoticed.
                        let theirs: f64 = expected.parse().expect("a double");
                        let ulps = ((value.to_bits() as i64) - (theirs.to_bits() as i64)).abs();
                        assert!(ulps <= 1, "{key}: {value} against {expected}, {ulps} ulps");
                        within_one_ulp.push(key.clone());
                    }
                }
                Err(error) => {
                    let row = refusals
                        .iter()
                        .find(|row| row[0] == key)
                        .unwrap_or_else(|| panic!("no refusal {key}"));
                    let (class, message) = row[1].split_once(':').expect("class and message");
                    assert_eq!(error.class(), class, "{key}");
                    assert_eq!(error.message(), message, "{key}");
                }
            }
            compared += 1;
        }
    }
    assert_eq!(compared, 48);
    assert_eq!(answers.len() + refusals.len(), 48);
    // Exactly one input of the grid is not bit-identical, and it is this one.
    assert_eq!(within_one_ulp, vec!["100.0,-1.0".to_string()]);
}

#[test]
fn the_posterior_is_the_probability_of_error() {
    let text = golden();
    let of = |odds: f64, prior: f64| -> String {
        rows(&text, "posterior")
            .into_iter()
            .find(|row| row[0] == format!("{},{}", label(odds), label(prior)))
            .expect("a posterior")[1]
            .to_string()
    };
    // The same prior, rising odds, falling error.
    assert_eq!(of(0.0, -1.0), "0.6321205588285577");
    assert_eq!(of(1.0, -1.0), "0.38730016321971794");
    assert_eq!(of(10.0, -1.0), "7.800378925839786E-5");
    // A prior of one makes it impossible.
    assert_eq!(of(-10.0, 0.0), "0.0");
}

#[test]
fn the_two_infinities_do_not_agree() {
    let text = golden();
    // Log odds of -Infinity against a prior of one: both entries are -Infinity.
    let nan = rows(&text, "posterior")
        .into_iter()
        .find(|row| row[0] == "-Infinity,0.0")
        .expect("a posterior")[1]
        .to_string();
    assert_eq!(nan, "NaN");
    assert!(posterior_probability_of_error(f64::NEG_INFINITY, 0.0)
        .expect("no refusal")
        .is_nan());
    // Log odds of +Infinity against a prior of zero: the sum is NaN and the normaliser refuses.
    let row = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "Infinity,-Infinity")
        .expect("a refusal");
    let (class, message) = row[1].split_once(':').expect("class and message");
    let error =
        posterior_probability_of_error(f64::INFINITY, f64::NEG_INFINITY).expect_err("not a number");
    assert_eq!(error.class(), class);
    assert_eq!(error.message(), message);
}

#[test]
fn rounding_the_error_is_a_clamp_that_keeps_nan() {
    let text = golden();
    let clamps = rows(&text, "round");
    for value in [
        -0.5,
        -1.0e-12,
        0.0,
        0.5,
        1.0,
        1.0 + 1.0e-10,
        2.0,
        f64::NAN,
        f64::NEG_INFINITY,
        f64::INFINITY,
    ] {
        let key = label(value);
        let expected = clamps
            .iter()
            .find(|row| row[0] == key)
            .unwrap_or_else(|| panic!("no clamp {key}"))[1];
        assert_eq!(
            java_double_to_string(round_finite_precision_errors(value)),
            expected,
            "{key}"
        );
    }
    // The two constants, printed by the same dump.
    assert_eq!(
        clamps
            .iter()
            .find(|row| row[0] == "EPSILON")
            .expect("epsilon")[1],
        java_double_to_string(EPSILON)
    );
    assert_eq!(
        clamps
            .iter()
            .find(|row| row[0] == "MIN_REPORTABLE_ERROR_PROBABILITY")
            .expect("the minimum")[1],
        java_double_to_string(MIN_REPORTABLE_ERROR_PROBABILITY)
    );
}

#[test]
fn the_tumour_odds_are_converted_from_log10() {
    let text = golden();
    let of = |name: &str| -> String {
        rows(&text, "tumorlogodds")
            .into_iter()
            .find(|row| row[0] == name)
            .unwrap_or_else(|| panic!("no conversion {name}"))[1]
            .to_string()
    };
    let rendered = |values: Option<Vec<f64>>| -> String {
        match values {
            None => "null".to_string(),
            Some(values) => format!(
                "[{}]",
                values
                    .iter()
                    .map(|value| java_double_to_string(*value))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        }
    };
    assert_eq!(rendered(tumor_log_odds(Some(&[6.0]))), of("one"));
    assert_eq!(rendered(tumor_log_odds(Some(&[6.0, 4.0]))), of("two"));
    assert_eq!(
        rendered(tumor_log_odds(Some(&[0.0, -3.0]))),
        of("zero-and-negative")
    );
    // A record with no annotation is null, which is not an empty array.
    assert_eq!(rendered(tumor_log_odds(None)), of("absent"));
    assert_eq!(of("absent"), "null");
}
