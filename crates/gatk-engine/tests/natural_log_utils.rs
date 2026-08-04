//! Conformance for `NaturalLogUtils` against the oracle.
//!
//! Golden from `tools/annotation-conformance/NaturalLogUtilsDump.java`, compared as **raw bit
//! patterns**: the question this suite exists to answer is whether a port built on a
//! permissively-licensed `exp` lands on the same `double`, and decimal rendering discards exactly
//! the last-place difference that question is about.
//!
//! # This suite is allowed to diverge, and says where
//!
//! The reference calls `Math.exp`; the port calls FDLIBM, which htsjdk-rs decision 0025 measured at
//! **1 ulp** worst case against it. So a row-by-row match is not guaranteed by construction, and
//! asserting one would be asserting something this programme has already measured to be false in
//! general.
//!
//! What the test asserts instead is the shape of the answer:
//!
//!  * every row that reaches **no** `exp` and **no** `log` must match **exactly**, because nothing
//!    in it can differ. Those are the rows where `logSumExp`'s accumulator is still `1.0` at the
//!    end — a single non-infinite entry, or a maximum with every other entry at `-Infinity`, or a
//!    difference so large that `1 + exp(diff)` rounds back to `1`;
//!  * every other row must be within **1 ulp**, and the test reports the worst distance it saw.
//!
//! A row that diverges by more than 1 ulp means the port's arithmetic is wrong, not that `exp` is
//! approximate, and that is the distinction the suite is built to make.

use std::io::Read;

use gatk_engine::natural_log_utils::{
    log1mexp, log_sum_exp, normalize_from_log_to_linear_space, posteriors,
};

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/natural_log_utils.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

/// `Double.doubleToRawLongBits` in reverse.
fn from_bits(hex: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(hex, 16).expect("hex bits"))
}

fn to_bits(value: f64) -> String {
    format!("{:x}", value.to_bits())
}

fn parse_list(text: &str) -> Vec<f64> {
    text.split(',').map(from_bits).collect()
}

/// The distance in representable doubles, which is the honest unit for "how far apart".
///
/// Only meaningful for two finite values of the same sign, which is every divergence this suite can
/// produce: the exact rows are exact, and the rest are ordinary finite results.
fn ulps_apart(ours: f64, theirs: f64) -> i64 {
    if ours.to_bits() == theirs.to_bits() {
        return 0;
    }
    if ours.is_nan() && theirs.is_nan() {
        return 0;
    }
    if !ours.is_finite() || !theirs.is_finite() {
        return i64::MAX;
    }
    (ours.to_bits() as i64 - theirs.to_bits() as i64).abs()
}

/// A row whose answer cannot depend on `exp` or `log`, so it must match to the bit.
///
/// Not a list of labels: the property is computed from the inputs, the same way the reference's
/// own control flow computes it. If the accumulator is still exactly `1.0` after the loop, no
/// `log` was called either, and the result is `maxValue` unchanged.
fn is_exact_by_construction(inputs: &[f64]) -> bool {
    if inputs.iter().all(|value| *value == f64::NEG_INFINITY) {
        return true;
    }
    let max =
        inputs.iter().copied().fold(
            f64::NEG_INFINITY,
            |best, value| if value > best { value } else { best },
        );
    if !max.is_finite() {
        return false;
    }
    // The same accumulation the port performs, asking only whether it moved off 1.0.
    let mut seen_max = false;
    let mut sum = 1.0f64;
    for &value in inputs {
        if !seen_max && value == max {
            seen_max = true;
            continue;
        }
        if value == f64::NEG_INFINITY {
            continue;
        }
        sum += (value - max).exp();
    }
    sum == 1.0
}

#[test]
fn the_arithmetic_is_the_references_within_one_ulp() {
    let text = golden();
    let mut compared = 0usize;
    let mut exact_rows = 0usize;
    let mut worst = 0i64;
    let mut worst_row = String::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[0] {
            "lse" => {
                let (label, expected, inputs) = (fields[1], fields[2], parse_list(fields[3]));
                let produced = log_sum_exp(&inputs);
                if let Some(class) = expected.strip_prefix("E:") {
                    // The refusal is on the accumulator, so the port must refuse the same rows.
                    let error = produced.expect_err(&format!("{label} must be refused"));
                    assert!(
                        class.starts_with(error.class()),
                        "{label}: expected {class}, got {}",
                        error.class()
                    );
                    compared += 1;
                    continue;
                }
                let ours = produced.unwrap_or_else(|_| panic!("{label} must not be refused"));
                let theirs = from_bits(expected);
                let distance = ulps_apart(ours, theirs);
                if is_exact_by_construction(&inputs) {
                    exact_rows += 1;
                    assert_eq!(
                        to_bits(ours),
                        expected,
                        "{label} reaches neither exp nor log, so it must be bit-identical"
                    );
                }
                assert!(
                    distance <= 1,
                    "{label}: {distance} ulp apart, ours {} theirs {expected}",
                    to_bits(ours)
                );
                if distance > worst {
                    worst = distance;
                    worst_row = label.to_string();
                }
                compared += 1;
            }
            "norm" => {
                let (label, expected, inputs) = (fields[1], fields[2], parse_list(fields[3]));
                let ours = normalize_from_log_to_linear_space(&inputs)
                    .unwrap_or_else(|_| panic!("{label} must not be refused"));
                let theirs: Vec<f64> = parse_list(expected);
                assert_eq!(ours.len(), theirs.len(), "{label}: length");
                for (index, (ours, theirs)) in ours.iter().zip(&theirs).enumerate() {
                    let distance = ulps_apart(*ours, *theirs);
                    assert!(distance <= 1, "{label}[{index}]: {distance} ulp apart");
                    if distance > worst {
                        worst = distance;
                        worst_row = format!("{label}[{index}]");
                    }
                    compared += 1;
                }
            }
            "post" => {
                let (label, expected) = (fields[1], fields[2]);
                let (priors, likelihoods) = fields[3].split_once('|').expect("two input lists");
                let ours = posteriors(&parse_list(priors), &parse_list(likelihoods))
                    .unwrap_or_else(|| panic!("{label} must produce a result"));
                let theirs: Vec<f64> = parse_list(expected);
                assert_eq!(ours.len(), theirs.len(), "{label}: length");
                for (index, (ours, theirs)) in ours.iter().zip(&theirs).enumerate() {
                    let distance = ulps_apart(*ours, *theirs);
                    assert!(distance <= 1, "{label}[{index}]: {distance} ulp apart");
                    if distance > worst {
                        worst = distance;
                        worst_row = format!("{label}[{index}]");
                    }
                    compared += 1;
                }
            }
            "l1me" => {
                let (input, expected) = (from_bits(fields[1]), fields[2]);
                let ours = log1mexp(input);
                let distance = ulps_apart(ours, from_bits(expected));
                assert!(
                    distance <= 1,
                    "log1mexp({input}): {distance} ulp apart, ours {}",
                    to_bits(ours)
                );
                if distance > worst {
                    worst = distance;
                    worst_row = format!("log1mexp({input})");
                }
                compared += 1;
            }
            _ => {}
        }
    }

    // 55 is what the golden carries today, not a round number: 22 logSumExp rows, 12 normalize
    // values, 10 posterior values and 11 log1mexp calls. The assertion is against the corpus
    // silently shrinking, so it tracks the real count rather than a guess at it.
    assert_eq!(compared, 55, "the golden changed size");
    assert!(
        exact_rows >= 6,
        "only {exact_rows} rows were exact by construction; the suite has lost the cases that \
         pin the algorithm rather than the exponential"
    );
    println!(
        "NaturalLogUtils: {compared} values compared, {exact_rows} exact by construction, \
         worst divergence {worst} ulp{}",
        if worst == 0 {
            String::new()
        } else {
            format!(" at {worst_row}")
        }
    );
}

/// The accumulator starting at `1.0` is what makes a single-term array exact, and it is worth an
/// assertion of its own because it is the property #96 rests on.
#[test]
fn a_single_term_returns_its_input_untouched() {
    for value in [0.0, -3.5, 700.0, -745.0, 1e-300] {
        assert_eq!(
            log_sum_exp(&[value]).expect("finite"),
            value,
            "a one-element logSumExp must return its input, with no exp and no log"
        );
    }
    // A maximum with everything else at negative infinity is the same case by a different route.
    assert_eq!(
        log_sum_exp(&[-1.5, f64::NEG_INFINITY]).expect("finite"),
        -1.5
    );
    assert_eq!(
        log_sum_exp(&[f64::NEG_INFINITY, -1.5]).expect("finite"),
        -1.5
    );
    // All negative infinity returns early, before the loop.
    assert_eq!(
        log_sum_exp(&[f64::NEG_INFINITY, f64::NEG_INFINITY]).expect("finite"),
        f64::NEG_INFINITY
    );
}
