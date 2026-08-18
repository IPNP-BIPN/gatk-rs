//! Conformance for `StrandArtifactFilter`'s E step against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/StrandArtifactEStepDump.java`.
//!
//! # What this suite is for
//!
//!  * **the strand counts come out of a string**, so an empty field and a word are both refusals and
//!    a one-entry table drops the whole filter;
//!  * **two branches answer a hard zero** rather than abstaining;
//!  * **the sequencing-error prior takes three steps in the indel size**, and an insertion and a
//!    deletion of the same length take the same one;
//!  * **the symbolic allele's data is removed before the totals are summed**;
//!  * **a prior of one answers a forward responsibility of exactly `1.0` beside a reverse one of
//!    `1.3E-25`**, which sum to more than one in the reals and back to exactly `1.0` in doubles, and
//!    a prior of zero answers two hard zeros.
//!
//! Every row reaches the beta binomial and `normalizeLog10`'s `Math.pow(10, x)`, so a divergence
//! would be named in [`ULPS_APART`].

use gatk_corpus as corpus;
use gatk_engine::somatic_clustering_model::AlternateAllele;
use gatk_engine::strand_artifact_filter::{
    calculate_artifact_probabilities, error_probabilities, indel_size, parse_strand_bias_table,
    strand_artifact_probability, EStep, StrandArtifactError, ANNOTATION, FILTER_NAME,
    INITIAL_ALPHA_STRAND, INITIAL_BETA_STRAND, INITIAL_STRAND_ARTIFACT_PRIOR,
};
use gatk_engine::tsv_table::java_double_to_string;

/// The rows that need an allowance, and how many ulps each needs.
const ULPS_APART: [(&str, i64); 0] = [];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/strand_artifact_estep.txt.gz"),
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

/// One record: its strand table, and one alternate per non-symbolic alternate allele.
struct Case {
    table: &'static str,
    reference_length: i32,
    /// The alternates whose data survives the symbolic removal, paired with the table's entries.
    alternates: Vec<i32>,
}

fn snv(table: &'static str) -> Case {
    Case {
        table,
        reference_length: 1,
        alternates: vec![1],
    }
}

fn case(label: &str) -> Case {
    match label {
        "balanced" => snv("50,50|10,10"),
        "forward-only" | "snv" => snv("50,50|20,0"),
        "reverse-only" => snv("50,50|0,20"),
        "one-forward-read" => snv("50,50|1,0"),
        "no-alt-reads" => snv("50,50|0,0"),
        "shallow" => snv("5,5|3,0"),
        "deep" => snv("2000,2000|400,0"),
        "one-base-deletion" => Case {
            table: "50,50|20,0",
            reference_length: 2,
            alternates: vec![1],
        },
        "two-base-insertion" => Case {
            table: "50,50|20,0",
            reference_length: 1,
            alternates: vec![3],
        },
        "three-base-insertion" => Case {
            table: "50,50|20,0",
            reference_length: 1,
            alternates: vec![4],
        },
        "four-base-deletion" => Case {
            table: "50,50|20,0",
            reference_length: 5,
            alternates: vec![1],
        },
        "five-base-deletion" => Case {
            table: "50,50|20,0",
            reference_length: 6,
            alternates: vec![1],
        },
        "two-alternates" => Case {
            table: "50,50|20,0|5,5",
            reference_length: 1,
            alternates: vec![1, 1],
        },
        // The symbolic alternate is last, so the removal leaves the first alternate's size in
        // place: one entry survives.
        "symbolic-alternate" => Case {
            table: "50,50|20,0",
            reference_length: 1,
            alternates: vec![1],
        },
        "one-entry-table" => snv("50,50"),
        "empty-table" | "no-table" => snv(""),
        "bracketed-table" => snv("[50,50|20,0]"),
        "spaced-table" => snv("50, 50 | 20, 0"),
        "non-integer-table" => snv("50,50|twenty,0"),
        "empty-field" => snv("50,50|,0"),
        other => panic!("no case named {other}"),
    }
}

/// A record's E steps, at the initial parameters.
fn steps(label: &str) -> Result<Vec<EStep>, StrandArtifactError> {
    let case = case(label);
    let table = parse_strand_bias_table(case.table)?;
    let sizes: Vec<i32> = case
        .alternates
        .iter()
        .map(|length| {
            indel_size(
                case.reference_length,
                AlternateAllele {
                    length: *length,
                    symbolic: false,
                },
            )
        })
        .collect();
    calculate_artifact_probabilities(
        &table,
        &sizes,
        INITIAL_STRAND_ARTIFACT_PRIOR,
        INITIAL_ALPHA_STRAND,
        INITIAL_BETA_STRAND,
    )
}

/// The dump's six direct calls to the package-private `strandArtifactProbability`.
fn direct(label: &str) -> Result<EStep, StrandArtifactError> {
    let (prior, forward, reverse, forward_alt, reverse_alt, size) = match label {
        "direct-default-prior" => (0.001, 50, 50, 20, 0, 0),
        "direct-high-prior" => (0.5, 50, 50, 20, 0, 0),
        "direct-zero-prior" => (0.0, 50, 50, 20, 0, 0),
        "direct-prior-of-one" => (1.0, 50, 50, 20, 0, 0),
        "direct-long-indel" => (0.001, 50, 50, 20, 0, 3),
        "direct-no-reads" => (0.001, 0, 0, 0, 0, 0),
        other => panic!("no direct case named {other}"),
    };
    strand_artifact_probability(
        prior,
        forward,
        reverse,
        forward_alt,
        reverse_alt,
        size,
        INITIAL_ALPHA_STRAND,
        INITIAL_BETA_STRAND,
    )
}

fn render(step: &EStep) -> String {
    format!(
        "{},{},{},{},{},{}",
        java_double_to_string(step.forward_artifact_responsibility),
        java_double_to_string(step.reverse_artifact_responsibility),
        step.forward_count,
        step.reverse_count,
        step.forward_alt_count,
        step.reverse_alt_count
    )
}

fn printed(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|v| java_double_to_string(*v)).collect();
    format!("[{}]", parts.join(", "))
}

/// Compare as text, with the allowance named row by row.
fn same(ours: &str, payload: &str, label: &str) {
    if ours == payload {
        return;
    }
    let Some(ulps) = ULPS_APART
        .iter()
        .find(|(name, _)| *name == label)
        .map(|(_, ulps)| *ulps)
    else {
        panic!("{label}: ours {ours}, reference {payload}");
    };
    for (ours, expected) in ours.split(',').zip(payload.split(',')) {
        if ours == expected {
            continue;
        }
        let (ours, theirs): (f64, f64) = (
            ours.parse().expect("a double"),
            expected.parse().expect("a double"),
        );
        let apart = (ours.to_bits() as i64 - theirs.to_bits() as i64).abs();
        assert!(apart <= ulps, "{label}: {apart} ulps apart, allowed {ulps}");
    }
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 50, "the golden's row count");
    // Each label's `estep` rows are consumed in order.
    let mut consumed: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "default" => {
                let ours = match label.as_str() {
                    "initialAlphaStrand" => java_double_to_string(INITIAL_ALPHA_STRAND),
                    "initialBetaStrand" => java_double_to_string(INITIAL_BETA_STRAND),
                    "initialStrandArtifactPrior" => {
                        java_double_to_string(INITIAL_STRAND_ARTIFACT_PRIOR)
                    }
                    other => panic!("no default named {other}"),
                };
                assert_eq!(*payload, ours, "default {label}");
            }
            "name" => assert_eq!(
                *payload,
                format!("{FILTER_NAME},ARTIFACT,{ANNOTATION},AS_SB_TABLE")
            ),
            "estep" if label.starts_with("direct-") => {
                same(&render(&direct(label).expect("answered")), payload, label)
            }
            "estep" => {
                let ours = steps(label).expect("answered");
                let index = consumed.entry(label.clone()).or_insert(0);
                same(&render(&ours[*index]), payload, label);
                *index += 1;
            }
            "prob" => {
                let ours = error_probabilities(&steps(label).expect("answered"));
                assert_eq!(printed(&ours), *payload, "prob {label}");
            }
            "error" => {
                let error = steps(label).expect_err("refused");
                assert_eq!(
                    *payload,
                    format!(
                        "{}:{}",
                        error.class().expect("a class"),
                        error.message().expect("a message")
                    ),
                    "error {label}"
                );
            }
            other => panic!("no row kind {other}"),
        }
    }
}
