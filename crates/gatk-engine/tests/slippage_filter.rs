//! Conformance for `PolymeraseSlippageFilter` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/PolymeraseSlippageFilterDump.java`.
//!
//! # What this suite is for
//!
//!  * **the beta is about `ADs[1]` and the likelihood beside it is about every alternate**, which
//!    only agree on a biallelic record;
//!  * **the prior's allele index is hard-coded to zero**, whichever allele slipped;
//!  * **`RPA` is parsed from its string form**, so `10.0` is a `NumberFormatException` rather than a
//!    truncation;
//!  * **the gate needs `ru.length() * rpa[0] >= minSlippageLength` and exactly one slip**, so an
//!    empty repeat unit and a two-repeat contraction both answer zero;
//!  * **and two calls through one model are two identical answers**, even though the prior inserts
//!    before it reads.
//!
//! Every row is compared. The probabilities reach the clustering model's `exp`, whose bit-exact
//! transcription is withdrawn under htsjdk-rs decision 0014, so any divergence would be named in
//! [`ULPS_APART`]; none is.

use gatk_corpus as corpus;
use gatk_engine::allele_filter::GenotypeData;
use gatk_engine::slippage_filter::{
    slippage_error_probabilities, SlippageError, ANNOTATION, DEFAULT_MIN_SLIPPAGE_LENGTH,
    DEFAULT_SLIPPAGE_RATE, FILTER_NAME,
};
use gatk_engine::somatic_clustering_model::{
    AlternateAllele, PriorArguments, SomaticClusteringModel,
};
use gatk_engine::tsv_table::java_double_to_string;

/// The rows that need the decision 0014 allowance. None of them does.
const ULPS_APART: [(&str, i64); 0] = [];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/slippage_filter.txt.gz"),
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

/// The record's reference is `AA`, its first alternate the deletion `A`, its second the insertion
/// `AAA`.
const REFERENCE_LENGTH: i32 = 2;

fn deletion() -> AlternateAllele {
    AlternateAllele {
        length: 1,
        symbolic: false,
    }
}

fn insertion() -> AlternateAllele {
    AlternateAllele {
        length: 3,
        symbolic: false,
    }
}

fn genotype(tumor: bool, allele_depths: &[i32]) -> GenotypeData<i32> {
    GenotypeData {
        tumor,
        allele_depths: allele_depths.to_vec(),
        values: Vec::new(),
    }
}

/// The dump's record for one label.
struct Case {
    repeats: Option<Vec<String>>,
    repeat_unit: Option<String>,
    genotypes: Vec<GenotypeData<i32>>,
    alternates: Vec<AlternateAllele>,
    slippage_rate: f64,
}

fn strings(values: &[&str]) -> Option<Vec<String>> {
    Some(values.iter().map(|v| v.to_string()).collect())
}

fn case(label: &str) -> Case {
    let mut case = Case {
        repeats: strings(&["10", "9"]),
        repeat_unit: Some("A".to_string()),
        genotypes: vec![genotype(true, &[80, 20]), genotype(false, &[90, 1])],
        alternates: vec![deletion()],
        slippage_rate: DEFAULT_SLIPPAGE_RATE,
    };
    match label {
        "contracted-by-one" | "called-twice-first" | "called-twice-second" => {}
        "expanded-by-one" => case.repeats = strings(&["9", "10"]),
        "contracted-by-two" => case.repeats = strings(&["10", "8"]),
        "two-base-repeat-unit" => {
            case.repeats = strings(&["5", "4"]);
            case.repeat_unit = Some("AT".to_string());
        }
        "below-the-minimum" => case.repeats = strings(&["7", "6"]),
        "at-the-minimum" => case.repeats = strings(&["8", "7"]),
        "empty-repeat-unit" => case.repeat_unit = Some(String::new()),
        "one-entry-in-rpa" => case.repeats = strings(&["10"]),
        "non-integer-rpa" => case.repeats = strings(&["ten", "nine"]),
        "decimal-rpa" => case.repeats = strings(&["10.0", "9.0"]),
        "no-rpa" => case.repeats = None,
        "no-ru" => case.repeat_unit = None,
        "no-alt-reads" => case.genotypes[0] = genotype(true, &[100, 0]),
        "all-alt-reads" => case.genotypes[0] = genotype(true, &[0, 100]),
        "shallow" => case.genotypes[0] = genotype(true, &[4, 2]),
        "deep" => case.genotypes[0] = genotype(true, &[800, 200]),
        "triallelic" => {
            case.repeats = strings(&["10", "9", "11"]);
            case.genotypes = vec![genotype(true, &[80, 20, 40]), genotype(false, &[90, 1, 0])];
            case.alternates = vec![deletion(), insertion()];
        }
        "slippage-rate-one" => case.slippage_rate = 1.0,
        "slippage-rate-zero" => case.slippage_rate = 0.0,
        "slippage-rate-tiny" => case.slippage_rate = 1.0e-12,
        other => panic!("no case named {other}"),
    }
    case
}

fn answer(label: &str, model: &mut SomaticClusteringModel) -> Result<Vec<f64>, SlippageError> {
    let case = case(label);
    slippage_error_probabilities(
        model,
        case.repeats.as_deref(),
        case.repeat_unit.as_deref(),
        &case.genotypes,
        &case.alternates,
        REFERENCE_LENGTH,
        DEFAULT_MIN_SLIPPAGE_LENGTH,
        case.slippage_rate,
    )
}

/// A fresh model per case, as the dump uses a fresh engine per case.
fn alone(label: &str) -> Result<Vec<f64>, SlippageError> {
    let mut model = SomaticClusteringModel::new(PriorArguments::new(), None);
    answer(label, &mut model)
}

fn printed(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|v| java_double_to_string(*v)).collect();
    format!("[{}]", parts.join(", "))
}

/// Compare as text, with the decision 0014 allowance named row by row.
fn same(ours: &[f64], payload: &str, label: &str) {
    let mine = printed(ours);
    if mine == payload {
        return;
    }
    let allowance = ULPS_APART
        .iter()
        .find(|(name, _)| *name == label)
        .map(|(_, ulps)| *ulps);
    let Some(ulps) = allowance else {
        panic!("{label}: ours {mine}, reference {payload}");
    };
    let expected: Vec<&str> = payload
        .trim_start_matches('[')
        .trim_end_matches(']')
        .split(", ")
        .collect();
    assert_eq!(ours.len(), expected.len(), "{label}: allele count");
    for (ours, expected) in ours.iter().zip(expected) {
        let theirs: f64 = expected.parse().expect("a double");
        let apart = (ours.to_bits() as i64 - theirs.to_bits() as i64).abs();
        assert!(apart <= ulps, "{label}: {apart} ulps apart, allowed {ulps}");
    }
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 25, "the golden's row count");
    // The two `called-twice` rows share one model, which the first call inserted into.
    let mut shared = SomaticClusteringModel::new(PriorArguments::new(), None);
    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "default" => {
                let ours = match label.as_str() {
                    "minSlippageLength" => DEFAULT_MIN_SLIPPAGE_LENGTH.to_string(),
                    "slippageRate" => java_double_to_string(DEFAULT_SLIPPAGE_RATE),
                    other => panic!("no default named {other}"),
                };
                assert_eq!(*payload, ours, "default {label}");
            }
            "name" => assert_eq!(
                *payload,
                format!("{FILTER_NAME},ARTIFACT,{ANNOTATION},RPA,RU"),
                "name {label}"
            ),
            "prob" => {
                let ours = if label.starts_with("called-twice") {
                    answer(label, &mut shared).expect("answered")
                } else {
                    alone(label).expect("answered")
                };
                same(&ours, payload, label);
            }
            "error" => {
                let error = alone(label).expect_err("refused");
                let rendered = format!(
                    "{}:{}",
                    error.class().expect("a class"),
                    error.message().expect("a message")
                );
                assert_eq!(*payload, rendered, "error {label}");
            }
            other => panic!("no row kind {other}"),
        }
    }
}
