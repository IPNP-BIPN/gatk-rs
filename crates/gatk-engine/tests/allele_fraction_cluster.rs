//! Conformance for the two allele-fraction clusters against GATK 4.6.2.0, compared as every row of
//! the dump.
//!
//! Golden from `tools/readfilter-conformance/AlleleFractionClusterDump.java`.
//!
//! # What this suite is for
//!
//!  * **the corrected log likelihood is a TLOD plus four Dirichlet normalisations**, two of them at
//!    the flat shape, which do not cancel in doubles;
//!  * **`logDirichletNormalization` on a single `1.0` is negative zero** and on a single `0.5` is
//!    positive zero, so the comparison is made on the formatted string;
//!  * **a parameter of zero is `NaN`**, commons-math's `logGamma` answering `NaN` at zero rather
//!    than diverging;
//!  * **the fuzzy binomial clamps its mean first**, so `1.0` and `2.0` produce the same shape as
//!    `0.99`;
//!  * **and the two shape refusals are two differently-worded messages.**
//!
//! Every row is reproduced, and every row is bit-identical: nothing here goes through `exp`, so
//! decision 0014 does not apply and there is no ulp allowance.

use gatk_corpus as corpus;
use gatk_engine::allele_fraction_cluster::{
    corrected_log_likelihood, fuzzy_binomial, log_likelihood, AlleleFractionCluster,
    BetaDistributionShape, Datum, ShapeError,
};
use gatk_engine::somatic_likelihoods::log_dirichlet_normalization;
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/allele_fraction_cluster.txt.gz"),
    )
}

fn double(text: &str) -> f64 {
    text.parse().unwrap_or_else(|_| panic!("a double: {text}"))
}

fn integer(text: &str) -> i32 {
    text.parse()
        .unwrap_or_else(|_| panic!("an integer: {text}"))
}

/// The dump's refusal line for a shape, whose two messages come from two different call paths.
fn shape_refusal(label: &str, error: &ShapeError) -> String {
    let (name, value) = match error {
        ShapeError::Alpha { alpha } => ("alpha", *alpha),
        ShapeError::Beta { beta } => ("beta", *beta),
    };
    format!(
        "error\t{label}\tjava.lang.IllegalArgumentException:{name} must be greater than 0 but got {}",
        java_double_to_string(value)
    )
}

/// The cluster a `loglike` or `corrected` label names, and the counts after it.
///
/// Labels are `binomial-<mean>-<total>,<alt>` and
/// `betabinomial-<alpha>,<beta>-<total>,<alt>`, so the split is on the last dash.
fn cluster_and_counts(label: &str) -> (AlleleFractionCluster, i32, i32) {
    let (head, counts) = label.rsplit_once('-').expect("a counts suffix");
    let (total, alt) = counts.split_once(',').expect("a total and an alt count");
    let cluster = if let Some(mean) = head.strip_prefix("binomial-") {
        AlleleFractionCluster::binomial(double(mean)).expect("a shape")
    } else if let Some(shape) = head.strip_prefix("betabinomial-") {
        let (alpha, beta) = shape.split_once(',').expect("two shape parameters");
        AlleleFractionCluster::beta_binomial(
            BetaDistributionShape::new(double(alpha), double(beta)).expect("a shape"),
        )
    } else {
        panic!("an unexpected cluster: {head}")
    };
    (cluster, integer(total), integer(alt))
}

/// The `tlod-<odds>-<total>,<alt>` rows, which are the flat cluster with a moving TLOD.
fn tlod_row(label: &str) -> String {
    let (head, counts) = label.rsplit_once('-').expect("a counts suffix");
    let (total, alt) = counts.split_once(',').expect("a total and an alt count");
    let odds = double(head.strip_prefix("tlod-").expect("a tlod label"));
    let datum = Datum::new(odds, 0.0, 0.0, integer(alt), integer(total), 0);
    let value = corrected_log_likelihood(&datum, BetaDistributionShape::FLAT);
    format!("corrected\t{label}\t{}", java_double_to_string(value))
}

/// The `datum` rows, whose label names the pair of probabilities rather than carrying them.
fn datum_row(label: &str) -> String {
    let (artifact, non_somatic) = match label {
        "both-zero" => (0.0, 0.0),
        "artifact-only" => (0.3, 0.0),
        "non-somatic-only" => (0.0, 0.3),
        "both" => (0.3, 0.3),
        "artifact-certain" => (1.0, 0.5),
        "tiny" => (1.0e-10, 1.0e-10),
        other => panic!("an unexpected datum: {other}"),
    };
    let datum = Datum::new(0.0, artifact, non_somatic, 5, 10, 0);
    format!(
        "datum\t{label}\t{}",
        java_double_to_string(datum.non_sequencing_error_prob())
    )
}

/// The `fuzzy` and `error` rows, whose labels are either a mean or one of the four named refusals.
fn shape_row(kind: &str, label: &str) -> String {
    let named: Option<(f64, f64)> = match label {
        "alpha-zero" => Some((0.0, 1.0)),
        "beta-zero" => Some((1.0, 0.0)),
        "alpha-negative" => Some((-1.0, 1.0)),
        "beta-nan" => Some((1.0, f64::NAN)),
        _ => None,
    };
    if let Some((alpha, beta)) = named {
        return match BetaDistributionShape::new(alpha, beta) {
            Ok(shape) => format!(
                "fuzzy\t{label}\t{},{}",
                java_double_to_string(shape.alpha()),
                java_double_to_string(shape.beta())
            ),
            Err(error) => shape_refusal(label, &error),
        };
    }
    assert_eq!(kind, "fuzzy", "only a mean reaches here");
    let shape = fuzzy_binomial(double(label)).expect("a shape");
    format!(
        "fuzzy\t{label}\t{},{}",
        java_double_to_string(shape.alpha()),
        java_double_to_string(shape.beta())
    )
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
            "dirichlet" => {
                let parameters: Vec<f64> = label.split(',').map(double).collect();
                format!(
                    "dirichlet\t{label}\t{}",
                    java_double_to_string(log_dirichlet_normalization(&parameters))
                )
            }
            "fuzzy" | "error" => shape_row(kind, label),
            "loglike" => {
                let (cluster, total, alt) = cluster_and_counts(label);
                let value = log_likelihood(cluster.shape(), total, alt).expect("in range");
                format!("loglike\t{label}\t{}", java_double_to_string(value))
            }
            "corrected" if label.starts_with("tlod-") => tlod_row(label),
            "corrected" => {
                let (cluster, total, alt) = cluster_and_counts(label);
                let datum = Datum::new(0.0, 0.0, 0.0, alt, total, 0);
                let value = cluster.corrected_log_likelihood(&datum);
                format!("corrected\t{label}\t{}", java_double_to_string(value))
            }
            "datum" => datum_row(label),
            other => panic!("an unexpected row: {other}"),
        };
        assert_eq!(ours, line);
        rows += 1;
    }
    assert_eq!(rows, 215, "the golden's row count");
}
