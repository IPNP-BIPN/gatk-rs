//! Conformance for `OptimizationUtils.max`, which is Apache Commons Math 3's `BrentOptimizer`,
//! against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/BrentOptimizerDump.java`.
//!
//! This is an engine primitive rather than a tool. It is here because `CreateSomaticPanelOfNormals`
//! fits its BETA field through it, and that field cannot be reproduced without it.
//!
//! # What this suite is for
//!
//!  * **an exact parabola being solved in one step**, so the tolerances never reach it;
//!  * **the tolerances deciding where the search stops on anything else**, badly when loose;
//!  * **the search being local to its INTERVAL rather than to its guess**;
//!  * **a symmetric function stopping a hair off centre**, which prints as a negative zero;
//!  * **and three things being refused rather than clamped.**

use gatk_corpus as corpus;
use gatk_engine::brent_optimizer::{maximize, OptimizerError, GOLDEN_SECTION};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/brent_optimizer.txt.gz"),
    )
}

fn measured(text: &str, label: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("max\t{label}=")))
        .unwrap_or_else(|| panic!("the golden carries max/{label}"))
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

/// The dump prints the point and the value to fourteen decimal places.
fn format(point: f64, value: f64) -> String {
    format!("{point:.14},{value:.14}")
}

fn parabola(x: f64) -> f64 {
    -(x - 3.0) * (x - 3.0) + 10.0
}

/// A peak the interpolation cannot fit, and piecewise linear so no transcendental enters it.
fn kinked(s: f64) -> f64 {
    -(s - 7.0).abs() - s / 50.0
}

/// Two maxima, at 4 and at 16, both at zero.
fn two_peaks(x: f64) -> f64 {
    -(x - 4.0) * (x - 4.0) * (x - 16.0) * (x - 16.0)
}

/// label, function, min, max, guess, relative tolerance, absolute tolerance, budget.
type Case = (&'static str, fn(f64) -> f64, f64, f64, f64, f64, f64, usize);

fn cases() -> Vec<Case> {
    vec![
        (
            "parabola-default",
            parabola as fn(f64) -> f64,
            0.0,
            10.0,
            1.0,
            0.001,
            0.001,
            1000,
        ),
        (
            "parabola-tight",
            parabola,
            0.0,
            10.0,
            1.0,
            1e-12,
            1e-12,
            1000,
        ),
        ("parabola-loose", parabola, 0.0, 10.0, 1.0, 0.1, 0.1, 1000),
        (
            "parabola-guess-high",
            parabola,
            0.0,
            10.0,
            9.0,
            0.001,
            0.001,
            1000,
        ),
        (
            "parabola-guess-exact",
            parabola,
            0.0,
            10.0,
            3.0,
            0.001,
            0.001,
            1000,
        ),
        (
            "parabola-guess-at-min",
            parabola,
            0.0,
            10.0,
            0.0,
            0.001,
            0.001,
            1000,
        ),
        (
            "parabola-guess-at-max",
            parabola,
            0.0,
            10.0,
            10.0,
            0.001,
            0.001,
            1000,
        ),
        ("kinked-settings", kinked, 0.01, 100.0, 1.0, 0.01, 0.1, 100),
        ("kinked-tight", kinked, 0.01, 100.0, 1.0, 1e-12, 1e-12, 1000),
        ("kinked-loose", kinked, 0.01, 100.0, 1.0, 0.5, 0.5, 1000),
        (
            "two-peaks-low-guess",
            two_peaks,
            0.0,
            20.0,
            1.0,
            0.001,
            0.001,
            1000,
        ),
        (
            "two-peaks-high-guess",
            two_peaks,
            0.0,
            20.0,
            14.0,
            0.001,
            0.001,
            1000,
        ),
        (
            "two-peaks-low-interval",
            two_peaks,
            0.0,
            10.0,
            1.0,
            0.001,
            0.001,
            1000,
        ),
        (
            "two-peaks-high-interval",
            two_peaks,
            10.0,
            20.0,
            14.0,
            0.001,
            0.001,
            1000,
        ),
        (
            "symmetric",
            (|x: f64| -x * x) as fn(f64) -> f64,
            -5.0,
            5.0,
            2.0,
            0.001,
            0.001,
            1000,
        ),
    ]
}

#[test]
fn every_maximum_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, function, min, max, guess, relative, absolute, budget) in cases() {
        let result = maximize(function, min, max, guess, relative, absolute, budget)
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));
        assert_eq!(
            format(result.point, result.value),
            measured(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 15, "the searches that returned a point");
}

/// The interpolation fits it exactly, so every setting returns the vertex and the tolerances never
/// reach it.
#[test]
fn an_exact_parabola_is_solved_in_one_step() {
    let text = golden();
    let vertex = measured(&text, "parabola-default");
    for label in [
        "parabola-tight",
        "parabola-loose",
        "parabola-guess-high",
        "parabola-guess-exact",
        "parabola-guess-at-min",
        "parabola-guess-at-max",
    ] {
        assert_eq!(measured(&text, label), vertex, "{label}");
    }
    assert_eq!(vertex, "3.00000000000000,10.00000000000000");
}

/// On a function the interpolation cannot fit they decide where it stops, and a loose one stops
/// well short of the maximum.
#[test]
fn the_tolerances_decide_where_it_stops() {
    let text = golden();
    let point = |label: &str| -> f64 {
        measured(&text, label)
            .split(',')
            .next()
            .expect("a point")
            .parse()
            .expect("a point")
    };
    let settings = point("kinked-settings");
    let tight = point("kinked-tight");
    let loose = point("kinked-loose");
    assert!((tight - 7.0).abs() < 1e-10, "the real peak");
    assert!(
        (settings - 7.0).abs() > 0.03,
        "the panel's own settings stop short"
    );
    assert!(
        (loose - 7.0).abs() > 0.28,
        "and a loose tolerance stops further short still"
    );
    // Each looser setting is a worse answer, which is what makes the tolerance matter.
    assert!(kinked(loose) < kinked(settings));
    assert!(kinked(settings) < kinked(tight));
}

/// The search is local to its INTERVAL rather than to its guess: over the whole range both guesses
/// land on the same peak, and bracketing is what separates them.
#[test]
fn the_search_is_local_to_its_interval() {
    let text = golden();
    let point = |label: &str| -> f64 {
        measured(&text, label)
            .split(',')
            .next()
            .expect("a point")
            .parse()
            .expect("a point")
    };
    let low_guess = point("two-peaks-low-guess");
    let high_guess = point("two-peaks-high-guess");
    assert!(
        low_guess > 15.9 && low_guess < 16.1,
        "a guess of 1 still finds 16"
    );
    assert!(high_guess > 15.9 && high_guess < 16.1);

    // Bracketing the interval around each peak is what actually separates them.
    let low_interval = point("two-peaks-low-interval");
    let high_interval = point("two-peaks-high-interval");
    assert!(low_interval > 3.9 && low_interval < 4.1);
    assert!(high_interval > 15.9 && high_interval < 16.1);
    // Both really are maxima, so neither is a failure.
    assert!(two_peaks(low_interval) > -0.001 && two_peaks(high_interval) > -0.01);
}

/// It returns the best point it evaluated, which for a symmetric function is a hair off centre: the
/// printed `-0.00000000000000` is a small NEGATIVE number rounded, not a negative zero.
#[test]
fn a_symmetric_function_stops_a_hair_off_centre() {
    let text = golden();
    assert_eq!(
        measured(&text, "symmetric"),
        "-0.00000000000000,-0.00000000000000"
    );
    let result = maximize(|x: f64| -x * x, -5.0, 5.0, 2.0, 0.001, 0.001, 1000).expect("a maximum");
    assert_eq!(result.point, -1.1102230246251565e-16);
    assert_ne!(result.point, 0.0, "not a zero of either sign");
    assert_eq!(format!("{:.14}", result.point), "-0.00000000000000");
}

/// A budget too small, a guess outside the interval and an interval the wrong way round.
#[test]
fn three_things_are_refused_rather_than_clamped() {
    let text = golden();

    let (class, message) = refusal(&text, "too-few-evaluations");
    assert_eq!(
        class,
        "org.apache.commons.math3.exception.TooManyEvaluationsException"
    );
    let produced = maximize(parabola, 0.0, 10.0, 1.0, 1e-14, 1e-14, 5).expect_err("no budget");
    assert_eq!(
        produced,
        OptimizerError::TooManyEvaluations { max_evaluations: 5 }
    );
    assert_eq!(produced.message(), message);
    // The same search with room converges, so it is the budget and not the function.
    assert!(maximize(parabola, 0.0, 10.0, 1.0, 1e-14, 1e-14, 1000).is_ok());

    for (label, guess) in [("guess-below-min", -1.0), ("guess-above-max", 11.0)] {
        let (class, message) = refusal(&text, label);
        assert_eq!(
            class, "org.apache.commons.math3.exception.OutOfRangeException",
            "{label}"
        );
        let produced = maximize(parabola, 0.0, 10.0, guess, 0.001, 0.001, 1000).expect_err(label);
        assert_eq!(
            produced,
            OptimizerError::OutOfRange {
                value: guess,
                min: 0.0,
                max: 10.0
            },
            "{label}"
        );
        assert_eq!(produced.message(), message, "{label}");
    }

    let (class, message) = refusal(&text, "inverted-interval");
    assert_eq!(
        class,
        "org.apache.commons.math3.exception.NumberIsTooLargeException"
    );
    let produced = maximize(parabola, 10.0, 0.0, 5.0, 0.001, 0.001, 1000).expect_err("inverted");
    assert_eq!(
        produced,
        OptimizerError::IntervalInverted {
            min: 10.0,
            max: 0.0
        }
    );
    assert_eq!(produced.message(), message);
    // The interval is checked before the guess, so an inverted interval with a bad guess reports
    // the interval.
    assert_eq!(
        maximize(parabola, 10.0, 0.0, -1.0, 0.001, 0.001, 1000).expect_err("inverted"),
        OptimizerError::IntervalInverted {
            min: 10.0,
            max: 0.0
        }
    );
}

/// The constant the golden-section step is scaled by.
#[test]
fn the_golden_section_constant() {
    assert_eq!(GOLDEN_SECTION, 0.5 * (3.0 - 5.0_f64.sqrt()));
}
