//! Conformance for `LearnReadOrientationModel` against GATK 4.6.2.0, compared as every flat prior,
//! every responsibility vector and every refusal.
//!
//! Golden from `tools/readfilter-conformance/LearnReadOrientationModelDump.java`.
//!
//! The EM loop is not compared here: it needs a counts file, which is `CollectF1R2Counts`'s output
//! and has its own suite. What is compared is the two functions the fit is built from.
//!
//! # What this suite is for
//!
//!  * **the flat prior not being flat over twelve states**;
//!  * **at most two artefact states being reachable for one site**;
//!  * **the F1R2 count alone deciding which of the two takes the mass**;
//!  * **`givenNotHomRef` renormalising rather than rescaling**;
//!  * **a state with a zero prior staying at zero**;
//!  * **depth and alternate depth moving the answer separately**;
//!  * **and two of the constructor's three refusals not being its own.**

use gatk_corpus as corpus;
use gatk_tools::learn_read_orientation_model::{
    compute_responsibilities, flat_prior, is_canonical_kmer, ref_to_ref_artifacts,
    validate_context, Base, ModelError, State, NUM_STATES,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/learn_read_orientation_model.txt.gz"),
    )
}

fn line(text: &str, kind: &str, name: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{name}"))
        .to_string()
}

fn refusal(text: &str, label: &str) -> (String, String) {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = row.split_once(':').expect("a class and a message");
    (class.to_string(), message.to_string())
}

/// The dump prints every double to ten decimal places, so the comparison is on that text.
fn format(values: &[f64; NUM_STATES]) -> String {
    values
        .iter()
        .map(|value| format!("{value:.10}"))
        .collect::<Vec<String>>()
        .join(",")
}

fn measured(text: &str, kind: &str, name: &str) -> Vec<f64> {
    line(text, kind, name)
        .split(',')
        .map(|value| value.parse().expect("a probability"))
        .collect()
}

/// label, reference, alternate, alt depth, F1R2 alt count, depth, prior, given-not-hom-ref.
type Case = (
    &'static str,
    Base,
    Base,
    i32,
    i32,
    i32,
    [f64; NUM_STATES],
    bool,
);

fn cases() -> Vec<Case> {
    let flat_a = flat_prior(Base::A);
    let flat_c = flat_prior(Base::C);
    let mut without_artifacts = flat_a;
    without_artifacts[State::F1R2C.index()] = 0.0;
    without_artifacts[State::F2R1C.index()] = 0.0;
    let mut only_f1r2 = [0.0; NUM_STATES];
    only_f1r2[State::F1R2C.index()] = 1.0;
    vec![
        ("no-alt", Base::A, Base::C, 0, 0, 50, flat_a, false),
        (
            "no-alt-not-hom-ref",
            Base::A,
            Base::C,
            0,
            0,
            50,
            flat_a,
            true,
        ),
        ("artifact-f1r2", Base::A, Base::C, 5, 5, 50, flat_a, false),
        ("artifact-f2r1", Base::A, Base::C, 5, 0, 50, flat_a, false),
        ("balanced", Base::A, Base::C, 6, 3, 50, flat_a, false),
        ("het", Base::A, Base::C, 25, 12, 50, flat_a, false),
        ("hom-var", Base::A, Base::C, 50, 25, 50, flat_a, false),
        (
            "artifact-deep",
            Base::A,
            Base::C,
            50,
            50,
            500,
            flat_a,
            false,
        ),
        ("artifact-g", Base::A, Base::G, 5, 5, 50, flat_a, false),
        ("alt-is-ref", Base::A, Base::A, 5, 5, 50, flat_a, false),
        ("ref-c", Base::C, Base::A, 5, 5, 50, flat_c, false),
        (
            "prior-without-artifacts",
            Base::A,
            Base::C,
            5,
            5,
            50,
            without_artifacts,
            false,
        ),
        (
            "prior-only-f1r2",
            Base::A,
            Base::C,
            5,
            5,
            50,
            only_f1r2,
            false,
        ),
        (
            "prior-only-f1r2-wrong-way",
            Base::A,
            Base::C,
            5,
            0,
            50,
            only_f1r2,
            false,
        ),
    ]
}

#[test]
fn every_value_matches_the_golden() {
    let text = golden();
    // The state order every prior array is indexed in.
    assert_eq!(
        line(&text, "order", "states"),
        State::all()
            .iter()
            .map(|state| state.name())
            .collect::<Vec<&str>>()
            .join(",")
    );
    let mut compared = 1;
    for base in [Base::A, Base::C, Base::G, Base::T] {
        assert_eq!(
            format(&flat_prior(base)),
            line(&text, "flat", base.name()),
            "{}",
            base.name()
        );
        compared += 1;
    }
    for (label, reference, alternate, alt_depth, f1r2, depth, prior, not_hom_ref) in cases() {
        assert_eq!(
            format(&compute_responsibilities(
                reference,
                alternate,
                alt_depth,
                f1r2,
                depth,
                &prior,
                not_hom_ref
            )),
            line(&text, "resp", label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 1 + 4 + 14, "the values the golden carries");
}

/// Ten states share the mass, not twelve, and which two are zeroed depends on the reference base.
#[test]
fn the_flat_prior_is_not_flat() {
    for base in [Base::A, Base::C, Base::G, Base::T] {
        let prior = flat_prior(base);
        assert_eq!(prior.iter().filter(|value| **value == 0.0).count(), 2);
        assert_eq!(
            prior.iter().filter(|value| **value == 0.1).count(),
            10,
            "a tenth, not a twelfth"
        );
        assert!((prior.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        for state in ref_to_ref_artifacts(base) {
            assert_eq!(prior[state.index()], 0.0, "{}", state.name());
            assert_eq!(state.artifact_base(), Some(base));
        }
    }
    // A twelfth is what a genuinely flat prior would give, and no entry is one.
    assert!(!flat_prior(Base::A).contains(&(1.0 / NUM_STATES as f64)));
}

/// At most two of the eight can be non-zero, and the F1R2 count alone decides which.
#[test]
fn the_f1r2_count_decides_which_artifact_state() {
    let text = golden();
    let of = |label: &str| measured(&text, "resp", label);
    let f1r2 = of("artifact-f1r2");
    let f2r1 = of("artifact-f2r1");
    assert!(f1r2[State::F1R2C.index()] > 0.97);
    assert!(f2r1[State::F2R1C.index()] > 0.97);
    assert_eq!(
        f1r2[State::F1R2C.index()],
        f2r1[State::F2R1C.index()],
        "the same value, the other way round"
    );
    // The six other artefact states are zero on both, and so are the ref-to-ref two.
    for state in State::all().iter().filter(|state| state.is_artifact()) {
        if matches!(state, State::F1R2C | State::F2R1C) {
            continue;
        }
        assert_eq!(f1r2[state.index()], 0.0, "{}", state.name());
    }
    // Split evenly, neither of the two takes anything worth having.
    let balanced = of("balanced");
    assert!(balanced[State::F1R2C.index()] < 0.001);
    assert!(balanced[State::SomaticHet.index()] > 0.9);
    // A different observed alternate moves the mass to a different pair.
    let g = of("artifact-g");
    assert!(g[State::F1R2G.index()] > 0.97);
    assert_eq!(g[State::F1R2C.index()], 0.0);
    // A different reference base zeroes a different pair, and the answer shifts with it.
    let c = of("ref-c");
    assert!(c[State::F1R2A.index()] > 0.97);
    // And an alternate that IS the reference leaves every artefact state at zero.
    let same = of("alt-is-ref");
    for state in State::all().iter().filter(|state| state.is_artifact()) {
        assert_eq!(same[state.index()], 0.0, "{}", state.name());
    }
}

/// Hom ref is zeroed AFTER the posteriors, so what is left is renormalised.
#[test]
fn given_not_hom_ref_renormalises() {
    let text = golden();
    let with = measured(&text, "resp", "no-alt");
    let without = measured(&text, "resp", "no-alt-not-hom-ref");
    assert!(with[State::HomRef.index()] > 0.75);
    assert_eq!(without[State::HomRef.index()], 0.0);
    // The two artefact states keep their RATIO to each other and to the rest.
    let scale = with[State::HomRef.index()];
    assert!(
        (without[State::F1R2C.index()] * (1.0 - scale) - with[State::F1R2C.index()]).abs() < 1e-9,
        "renormalised over what is left, not rescaled"
    );
    assert!((without.iter().sum::<f64>() - 1.0).abs() < 1e-12);
}

/// It stays at zero, so the flat prior's zeros stick, and a degenerate prior forces the answer.
#[test]
fn a_state_with_a_zero_prior_stays_at_zero() {
    let text = golden();
    let ruled_out = measured(&text, "resp", "prior-without-artifacts");
    assert_eq!(ruled_out[State::F1R2C.index()], 0.0);
    assert_eq!(ruled_out[State::F2R1C.index()], 0.0);
    // The same counts that gave F1R2_C 0.97 under the flat prior give it nothing at all here.
    assert!(measured(&text, "resp", "artifact-f1r2")[State::F1R2C.index()] > 0.97);

    // A prior on ONE state returns that state at 1.0 whatever the counts say, in either direction.
    for label in ["prior-only-f1r2", "prior-only-f1r2-wrong-way"] {
        let forced = measured(&text, "resp", label);
        assert_eq!(forced[State::F1R2C.index()], 1.0, "{label}");
        assert_eq!(forced.iter().filter(|value| **value != 0.0).count(), 1);
    }
}

/// The same alternate fraction at ten times the depth is more certain.
#[test]
fn depth_and_alt_depth_move_the_answer_separately() {
    let text = golden();
    let shallow = measured(&text, "resp", "artifact-f1r2")[State::F1R2C.index()];
    let deep = measured(&text, "resp", "artifact-deep")[State::F1R2C.index()];
    assert!(deep > shallow);
    assert!(shallow > 0.97 && shallow < 0.98);
    assert!(deep > 0.9999997);

    // And a half fraction is a germline het rather than an artefact, at either orientation.
    let het = measured(&text, "resp", "het");
    assert!(het[State::GermlineHet.index()] > 0.7);
    assert_eq!(het[State::F1R2C.index()], 0.0);
    let hom = measured(&text, "resp", "hom-var");
    assert!(hom[State::HomVar.index()] > 0.99);
}

/// Two are validated with their own messages; the third is Apache Commons refusing a matrix with
/// zero rows, which the validations never reach.
#[test]
fn two_of_the_three_refusals_are_not_the_constructors_own() {
    let text = golden();

    let (class, message) = refusal(&text, "bad-context");
    assert_eq!(class, "java.lang.IllegalStateException");
    let produced = validate_context("ACGT", 1).expect_err("the wrong length");
    assert_eq!(
        produced,
        ModelError::ContextLength {
            context: "ACGT".to_string()
        }
    );
    assert_eq!(produced.message(), message);

    let (class, message) = refusal(&text, "non-canonical");
    assert_eq!(class, "java.lang.IllegalStateException");
    let produced = validate_context("AGT", 1).expect_err("a non-canonical kmer");
    assert_eq!(
        produced,
        ModelError::NonCanonicalKmer {
            context: "AGT".to_string()
        }
    );
    assert_eq!(produced.message(), message);

    // The one that is not a validation at all.
    let (class, message) = refusal(&text, "empty-design-matrix");
    assert_eq!(
        class,
        "org.apache.commons.math3.exception.NotStrictlyPositiveException"
    );
    let produced = validate_context("TCA", 0).expect_err("an empty design matrix");
    assert_eq!(produced, ModelError::EmptyDesignMatrix);
    assert_eq!(produced.message(), message);
    // TCA passes both validations, which is what makes the crash reachable.
    assert!(is_canonical_kmer("TCA"));
    assert!(!is_canonical_kmer("AGT"), "the middle base decides");
    assert!(validate_context("TCA", 1).is_ok());
}
