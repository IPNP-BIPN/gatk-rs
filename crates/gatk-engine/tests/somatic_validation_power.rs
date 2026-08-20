//! Conformance for the power a validation pileup has, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/SomaticValidationPowerDump.java`.
//!
//! # What this suite is for
//!
//!  * **the cumulative probability is a plain uncompensated loop**, `+=` over `probability(i)`,
//!    which is not the `DoubleStream.sum` every other accumulation in this port has had to be;
//!  * **it stops growing past the trial count**, because `k > n` short-circuits to negative
//!    infinity before any log term is evaluated, so the row at `n + 1` repeats the row at `n`;
//!  * **a negative `k` is a refusal**, not a zero, which is what makes `minCount - 1` a
//!    precondition rather than a convenience;
//!  * **the minimum count is a binomial quantile floored at two**, so a pileup of no reads and a
//!    clean pileup of 317 both answer two;
//!  * **and the power's shapes are the discovery counts plus one**, never zero, which keeps the
//!    distribution constructor's own refusals out of reach.
//!
//! Every row is reproduced, and every row but ten is bit-identical. Each term goes through
//! `Math.exp(logProbability(k))`, which under decision 0014 is the platform's `exp` rather than a
//! transcription, and the ten rows that need an allowance are the tail of one uniform distribution
//! where an ulp in one term is carried by the uncompensated sum into every partial sum after it.

use gatk_corpus as corpus;
use gatk_engine::beta_binomial::{BetaBinomialDistribution, BetaBinomialError};
use gatk_engine::somatic_validation_power::{
    calculate_min_count_for_signal, calculate_power, PowerError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/somatic_validation_power.txt.gz"),
    )
}

/// A shape field, written as the raw bits of the double.
fn bits(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap_or_else(|_| panic!("bits: {text}")))
}

fn integer(text: &str) -> i32 {
    text.parse()
        .unwrap_or_else(|_| panic!("an integer: {text}"))
}

/// The dump's `%016x`, which is the raw bits and not the canonical NaN.
fn hex(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

/// The rows that need the decision 0014 allowance, and how many ulps each needs.
///
/// One family: the flat beta-binomial (`alpha = beta = 1`) over twenty trials, from `k = 12` on.
/// Its twenty-one terms are all near `1/21` but not bit-identical to each other, and one of them
/// comes back from the platform's `exp` an ulp below the reference's. The sum is a plain `+=` with
/// no compensation, so that ulp is not absorbed: it stands in every partial sum after it, and a
/// second term compounds it to two from `k = 15` before the last two rows round back to one.
///
/// Nothing else in the golden needs an allowance. The three deeper shapes agree bit for bit at
/// every `k`, and every `power` row does too, because a power's cumulative sum stops at
/// `minCount - 1`, which is four at the most.
const ALLOWANCE: &[(&str, i64)] = &[
    ("3ff0000000000000,3ff0000000000000,20,12", 1),
    ("3ff0000000000000,3ff0000000000000,20,13", 1),
    ("3ff0000000000000,3ff0000000000000,20,14", 1),
    ("3ff0000000000000,3ff0000000000000,20,15", 2),
    ("3ff0000000000000,3ff0000000000000,20,16", 2),
    ("3ff0000000000000,3ff0000000000000,20,17", 2),
    ("3ff0000000000000,3ff0000000000000,20,18", 2),
    ("3ff0000000000000,3ff0000000000000,20,19", 2),
    ("3ff0000000000000,3ff0000000000000,20,20", 1),
    ("3ff0000000000000,3ff0000000000000,20,21", 1),
];

fn distribution(alpha: &str, beta: &str, n: &str) -> BetaBinomialDistribution {
    BetaBinomialDistribution::new(bits(alpha), bits(beta), integer(n)).expect("shapes in range")
}

/// The three refusals, reported as the reference's exception and message.
fn refusal(label: &str) -> String {
    let text = match label {
        "negative-k-cumulative" => {
            match BetaBinomialDistribution::new(1.0, 1.0, 5)
                .expect("shapes in range")
                .cumulative_probability(-1)
            {
                Err(BetaBinomialError::SuccessesNegative { .. }) => {
                    "java.lang.IllegalArgumentException:Number of successes must be greater than \
                     or equal to zero."
                        .to_string()
                }
                other => panic!("{label}: {other:?}"),
            }
        }
        "ratio-above-one" | "negative-total" => {
            let answer = if label == "ratio-above-one" {
                calculate_min_count_for_signal(10, 1.5)
            } else {
                calculate_min_count_for_signal(-1, 0.1)
            };
            match answer {
                Err(error @ (PowerError::RatioOutOfRange | PowerError::NegativeTotalCount)) => {
                    format!("{}:{}", error.java_class(), error.message())
                }
                other => panic!("{label}: {other:?}"),
            }
        }
        other => panic!("an unexpected refusal: {other}"),
    };
    format!("error\t{label}\t{text}")
}

/// One row's value against the golden's, exactly unless the allowance names it.
fn agrees(kind: &str, label: &str, ours: f64, theirs: &str) -> bool {
    if hex(ours) == theirs {
        return true;
    }
    let Some((_, allowed)) = ALLOWANCE.iter().find(|(row, _)| *row == label) else {
        return false;
    };
    assert_eq!(
        kind, "cumulative",
        "{label}: an allowance on the wrong kind"
    );
    let reference = bits(theirs);
    let apart = ((ours.to_bits() as i64) - (reference.to_bits() as i64)).abs();
    apart <= *allowed
}

#[test]
fn every_row_matches_the_golden() {
    let text = golden();
    let mut rows = 0;
    let mut allowed_rows = 0;
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.splitn(3, '\t');
        let kind = fields.next().expect("a kind");
        let rest = fields.next().expect("a label");
        if kind == "error" {
            assert_eq!(refusal(rest), line);
            rows += 1;
            continue;
        }
        let (label, theirs) = rest.split_once('=').expect("a value");
        let parts: Vec<&str> = label.split(',').collect();
        match kind {
            "moments" => {
                let it = distribution(parts[0], parts[1], parts[2]);
                let ours = format!(
                    "{},{}",
                    hex(it.numerical_mean()),
                    hex(it.numerical_variance())
                );
                assert_eq!(ours, theirs, "moments {label}");
            }
            "cumulative" => {
                let ours = distribution(parts[0], parts[1], parts[2])
                    .cumulative_probability(integer(parts[3]))
                    .expect("a count in range");
                assert!(
                    agrees(kind, label, ours, theirs),
                    "cumulative {label}: {} against {theirs}",
                    hex(ours)
                );
                if hex(ours) != theirs {
                    allowed_rows += 1;
                }
            }
            "power" => {
                let ours = calculate_power(
                    integer(parts[0]),
                    integer(parts[1]),
                    integer(parts[2]),
                    integer(parts[3]),
                )
                .expect("a pileup in range");
                assert_eq!(hex(ours), theirs, "power {label}");
            }
            "mincount" => {
                let count = calculate_min_count_for_signal(integer(parts[0]), bits(parts[1]))
                    .expect("a ratio in range");
                assert_eq!(count.to_string(), theirs, "mincount {label}");
            }
            other => panic!("an unexpected row: {other}"),
        }
        rows += 1;
    }
    assert_eq!(rows, 304, "the golden's row count");
    assert_eq!(allowed_rows, ALLOWANCE.len(), "the rows that needed an ulp");
}
