//! Conformance for the beta-binomial against GATK 4.6.2.0, compared as every row of the dump.
//!
//! Golden from `tools/readfilter-conformance/BetaBinomialDump.java`.
//!
//! # What this suite is for
//!
//!  * **`logBeta` is the four-branch NSWC expansion**, not `logGamma(p) + logGamma(q) -
//!    logGamma(p + q)`, and the fifty-two shapes cross every branch boundary;
//!  * **`logBeta(1, 1)` is negative zero**, which the comparison catches only because it is made on
//!    the formatted string rather than on the double;
//!  * **`binomialCoefficientLog` picks its route by `n`**, an exact integer coefficient at ten and a
//!    rounded floating product at a hundred and a thousand;
//!  * **the flat cluster is not bit-identically uniform**, landing on three different doubles at
//!    `n = 10`;
//!  * **and a count past the total is negative infinity while a negative one is a refusal**.
//!
//! Every row is reproduced, and every row is bit-identical: unlike the sibling Mutect suites nothing
//! here goes through `exp`, so decision 0014 does not apply and there is no ulp allowance.

use gatk_corpus as corpus;
use gatk_engine::beta_binomial::{BetaBinomialDistribution, BetaBinomialError};
use gatk_engine::tsv_table::java_double_to_string;
use jmath::beta::log_beta;
use jmath::combinatorics::binomial_coefficient_log;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/beta_binomial.txt.gz"),
    )
}

/// A label field, parsed the way `Double.toString` wrote it.
fn double(text: &str) -> f64 {
    text.parse().unwrap_or_else(|_| panic!("a double: {text}"))
}

fn integer(text: &str) -> i32 {
    text.parse()
        .unwrap_or_else(|_| panic!("an integer: {text}"))
}

/// The dump's line for one beta-binomial, refusal included.
fn beta_binomial_line(label: &str) -> String {
    let fields: Vec<&str> = label.split(',').collect();
    let answer =
        BetaBinomialDistribution::new(double(fields[0]), double(fields[1]), integer(fields[2]))
            .and_then(|distribution| distribution.log_probability(integer(fields[3])));
    match answer {
        Ok(value) => format!("betabinom\t{label}\t{}", java_double_to_string(value)),
        // The one refusal reachable here, reported as the reference's exception and message.
        Err(BetaBinomialError::SuccessesNegative { .. }) => format!(
            "error\t{label}\tjava.lang.IllegalArgumentException:Number of successes must be \
             greater than or equal to zero."
        ),
        Err(other) => panic!("{label}: {other:?}"),
    }
}

#[test]
fn every_row_matches_the_golden() {
    let text = golden();
    let mut rows = 0;
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.splitn(3, '\t');
        let kind = fields.next().expect("a kind");
        let label = fields.next().expect("a label");
        let ours = match kind {
            "logbeta" => {
                let (p, q) = label.split_once(',').expect("two shapes");
                let value = log_beta(double(p), double(q)).expect("in range");
                format!("logbeta\t{label}\t{}", java_double_to_string(value))
            }
            "binomlog" => {
                let (n, k) = label.split_once(',').expect("a total and a count");
                let value = binomial_coefficient_log(integer(n) as i64, integer(k) as i64)
                    .expect("in range");
                format!("binomlog\t{label}\t{}", java_double_to_string(value))
            }
            "betabinom" | "error" => beta_binomial_line(label),
            other => panic!("an unexpected row: {other}"),
        };
        assert_eq!(ours, line);
        rows += 1;
    }
    assert_eq!(rows, 89, "the golden's row count");
}
