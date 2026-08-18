//! Conformance for `NormalArtifactFilter` and the commons-math under it, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/NormalArtifactFilterDump.java`.
//!
//! # What this suite is for
//!
//!  * **`Beta.regularizedBeta` is two arrangements of the same fraction**, and which one a call
//!    takes decides its last digits: one ulp either side of the branch boundary answers the same
//!    double by two different routes, and the boundary itself answers another;
//!  * **the binomial CDF's guards are answers**, not refusals: the filter asks about `-1` whenever
//!    the normal is clean, and about `x >= trials` whenever every read supports the allele;
//!  * **a missing `MBQ` is an out-of-bounds, not the imputed 30**;
//!  * **only the normal side of the ratio gate is guarded**, so a record with no tumour depth
//!    compares against NaN and falls through instead of returning zero;
//!  * **and a missing required annotation is `0.0` per allele here**, where the per-allele base
//!    class answers an empty list.
//!
//! Every row is compared. The `posterior` and `prob` rows reach the clustering model's `exp`, whose
//! bit-exact transcription is withdrawn under htsjdk-rs decision 0014, so they carry an allowance
//! named row by row in [`ULPS_APART`]; everything else, the twenty-two regularized betas included,
//! is bit-identical.

use gatk_corpus as corpus;
use gatk_engine::allele_filter::GenotypeData;
use gatk_engine::math_utils::qual_to_error_prob;
use gatk_engine::mutect_engine::posterior_probability_of_error;
use gatk_engine::normal_artifact_filter::{
    normal_artifact_error_probabilities, NormalArtifactError, DEFAULT_NORMAL_P_VALUE_THRESHOLD,
    FILTER_NAME,
};
use gatk_engine::somatic_clustering_model::{PriorArguments, SomaticClusteringModel};
use gatk_engine::tsv_table::java_double_to_string;
use jmath::beta::regularized_beta;
use jmath::binomial::cumulative_probability;

/// The rows that need the decision 0014 allowance, and how many ulps each needs.
///
/// One row, and it needs **three**. `posterior nalod-2.0` asks `logSumExp` for
/// `exp(-2.4079456086518722)`, where `StrictMath.exp` answers `0x3fb70a3d70a3d708` and `Math.exp`
/// answers `0x3fb70a3d70a3d709`: one ulp, the bound htsjdk-rs decision 0025 measured. That single
/// ulp moves the accumulator from `1.09` to `1.0899999999999999`, which moves the log-sum by an
/// ulp, and the final `exp` of a value shifted by that lands three ulps from the reference's
/// `0.0825688073394495`. The other six posteriors and every one of the sixteen `prob` rows are
/// bit-identical, so this is one input finding the seam rather than the seam being everywhere.
const ULPS_APART: [(&str, i64); 1] = [("nalod-2.0", 3)];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/normal_artifact_filter.txt.gz"),
    )
}

/// The golden's rows as `(kind, label, payload)`.
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

/// `<inputs>=<answer>`, which is how every arithmetic row is written.
fn split(payload: &str) -> (Vec<f64>, &str) {
    let (inputs, answer) = payload.rsplit_once('=').expect("an answer");
    (
        inputs
            .split(',')
            .map(|value| value.parse::<f64>().expect("a double"))
            .collect(),
        answer,
    )
}

/// Compare as text, so a negative zero and a NaN are compared as the reference wrote them.
///
/// A row that diverges is reported rather than panicked on, so one run names every divergence.
fn same(ours: f64, expected: &str, what: &str, complaints: &mut Vec<String>) {
    let printed = java_double_to_string(ours);
    if printed == expected {
        return;
    }
    let theirs: f64 = expected.parse().expect("a double");
    let apart = (ours.to_bits() as i64 - theirs.to_bits() as i64).abs();
    match ULPS_APART
        .iter()
        .find(|(label, _)| *label == what)
        .map(|(_, ulps)| *ulps)
    {
        Some(ulps) if apart <= ulps => {}
        Some(ulps) => complaints.push(format!(
            "{what}: {apart} ulps apart, allowed {ulps} (ours {printed}, reference {expected})"
        )),
        None => complaints.push(format!(
            "{what}: ours {printed}, reference {expected} ({apart} ulps)"
        )),
    }
}

/// One tumour sample and, usually, one normal.
fn genotypes(tumor: &[i32], normal: Option<&[i32]>) -> Vec<GenotypeData<i32>> {
    let mut samples = vec![GenotypeData {
        tumor: true,
        allele_depths: tumor.to_vec(),
        values: Vec::new(),
    }];
    if let Some(normal) = normal {
        samples.push(GenotypeData {
            tumor: false,
            allele_depths: normal.to_vec(),
            values: Vec::new(),
        });
    }
    samples
}

/// The dump's record for one `prob` label: the arguments the filter is handed.
struct Case {
    tumor_log_10_odds: Option<Vec<f64>>,
    normal_artifact_log_10_odds: Option<Vec<f64>>,
    median_base_qualities: Vec<i32>,
    tumor: Vec<i32>,
    normal: Option<Vec<i32>>,
    allele_count: usize,
}

fn case(label: &str) -> Case {
    // The base record: a triallelic site whose normal carries the allele.
    let mut case = Case {
        tumor_log_10_odds: Some(vec![20.0, 6.0]),
        normal_artifact_log_10_odds: Some(vec![2.0, 0.5]),
        median_base_qualities: vec![30, 30, 30],
        tumor: vec![80, 20, 5],
        normal: Some(vec![90, 10, 2]),
        allele_count: 3,
    };
    match label {
        "normal-carries-the-allele" => {}
        "second-allele" => case.tumor_log_10_odds = Some(vec![6.0, 20.0]),
        "normal-is-clean" => case.normal = Some(vec![90, 0, 0]),
        "normal-at-the-ratio" => case.normal = Some(vec![80, 2, 0]),
        "no-normal-sample" => case.normal = None,
        "no-tumor-depth" => case.tumor = vec![0, 0, 0],
        "no-tumor-depth-clean-normal" => {
            case.tumor = vec![0, 0, 0];
            case.normal = Some(vec![90, 0, 0]);
        }
        "normal-is-all-alt" => case.normal = Some(vec![0, 20, 0]),
        "low-median-base-quality" => case.median_base_qualities = vec![2, 2, 2],
        "high-median-base-quality" => case.median_base_qualities = vec![60, 60, 60],
        "no-mbq" => case.median_base_qualities = Vec::new(),
        "negative-nalod" => {
            case.median_base_qualities = vec![2, 2, 2];
            case.normal_artifact_log_10_odds = Some(vec![-5.0, -1.0]);
        }
        "large-nalod" => {
            case.median_base_qualities = vec![2, 2, 2];
            case.normal_artifact_log_10_odds = Some(vec![20.0, 20.0]);
        }
        "no-nalod" => case.normal_artifact_log_10_odds = None,
        "no-tlod" => case.tumor_log_10_odds = None,
        "biallelic" => {
            case.tumor_log_10_odds = Some(vec![20.0]);
            case.normal_artifact_log_10_odds = Some(vec![2.0]);
            case.median_base_qualities = vec![30, 30];
            case.tumor = vec![80, 20];
            case.normal = Some(vec![90, 10]);
            case.allele_count = 2;
        }
        other => panic!("no case named {other}"),
    }
    case
}

fn answer(label: &str) -> Result<Vec<f64>, NormalArtifactError> {
    let case = case(label);
    let model = SomaticClusteringModel::new(PriorArguments::new(), None);
    normal_artifact_error_probabilities(
        &model,
        case.tumor_log_10_odds.as_deref(),
        case.normal_artifact_log_10_odds.as_deref(),
        &case.median_base_qualities,
        &genotypes(&case.tumor, case.normal.as_deref()),
        case.allele_count,
        DEFAULT_NORMAL_P_VALUE_THRESHOLD,
    )
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 65, "the golden's row count");
    let mut complaints: Vec<String> = Vec::new();
    let mut seen = 0;
    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "regbeta" => {
                let (inputs, expected) = split(payload);
                let ours = regularized_beta(inputs[0], inputs[1], inputs[2]).expect("answered");
                same(ours, expected, label, &mut complaints);
            }
            "cdf" => {
                let (inputs, expected) = split(payload);
                let ours = cumulative_probability(inputs[0] as i32, inputs[1], inputs[2] as i32)
                    .expect("answered");
                same(ours, expected, label, &mut complaints);
            }
            "qualerr" => {
                let (inputs, expected) = split(payload);
                same(
                    qual_to_error_prob(inputs[0]),
                    expected,
                    label,
                    &mut complaints,
                );
            }
            "prior" => {
                let model = SomaticClusteringModel::new(PriorArguments::new(), None);
                same(
                    model.log_prior_of_variant_versus_artifact(),
                    payload,
                    label,
                    &mut complaints,
                );
            }
            "posterior" => {
                let (inputs, expected) = split(payload);
                let model = SomaticClusteringModel::new(PriorArguments::new(), None);
                let ours = posterior_probability_of_error(
                    inputs[0],
                    model.log_prior_of_variant_versus_artifact(),
                )
                .expect("finite");
                same(ours, expected, label, &mut complaints);
            }
            "name" => {
                // `<filterName>,<errorType>,<annotation>,<required annotations>`.
                assert_eq!(*payload, format!("{FILTER_NAME},ARTIFACT,none,NALOD,TLOD"));
            }
            "prob" => {
                let probabilities = answer(label).expect("answered");
                let printed: Vec<String> = probabilities
                    .iter()
                    .map(|value| java_double_to_string(*value))
                    .collect();
                let ours = format!("[{}]", printed.join(", "));
                if ours != *payload {
                    // Fall back to the per-value comparison, which is where an allowance applies.
                    let expected: Vec<&str> = payload
                        .trim_start_matches('[')
                        .trim_end_matches(']')
                        .split(", ")
                        .collect();
                    assert_eq!(probabilities.len(), expected.len(), "{label}: allele count");
                    for (ours, expected) in probabilities.iter().zip(expected) {
                        same(*ours, expected, label, &mut complaints);
                    }
                }
            }
            "error" => {
                // `IndexOutOfBoundsException: Index 0 out of bounds for length 0`, which is
                // `.get(0)` on the empty list an absent `MBQ` answers.
                assert_eq!(label, "prob-no-mbq");
                assert_eq!(
                    payload,
                    "java.lang.IndexOutOfBoundsException:Index 0 out of bounds for length 0"
                );
                assert_eq!(
                    answer("no-mbq"),
                    Err(NormalArtifactError::MedianBaseQualityMissing)
                );
            }
            other => panic!("no row kind {other}"),
        }
        seen += 1;
    }
    assert_eq!(seen, rows.len(), "every row compared");
    assert!(complaints.is_empty(), "{}", complaints.join("\n"));
}
