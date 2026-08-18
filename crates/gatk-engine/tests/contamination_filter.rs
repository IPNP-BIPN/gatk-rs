//! Conformance for `ContaminationFilter` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/ContaminationFilterDump.java`.
//!
//! # What this suite is for
//!
//!  * **this filter can answer NaN**, for an allele `POPAF` does not cover and for every allele of a
//!    record with no tumour sample;
//!  * **a `POPAF` longer than the record is an exception**, not a truncation;
//!  * **the clamp makes a contamination of one and a contamination of two the same number**;
//!  * **the two contamination hypotheses are compared by maximum**, not summed;
//!  * **and a contamination of zero, a negative one and an infinite `POPAF` all answer exactly
//!    `0.0`**, which is what the default estimate does to every record.
//!
//! Every row reaches the clustering model's `exp`, whose bit-exact transcription is withdrawn under
//! htsjdk-rs decision 0014, so any divergence would be named in [`ULPS_APART`].

use gatk_corpus as corpus;
use gatk_engine::allele_filter::GenotypeData;
use gatk_engine::contamination_filter::{
    contamination_error_probabilities, ContaminationError, ANNOTATION, DEFAULT_CONTAMINATION,
    FILTER_NAME,
};
use gatk_engine::somatic_clustering_model::{
    AlternateAllele, PriorArguments, SomaticClusteringModel,
};
use gatk_engine::tsv_table::java_double_to_string;

/// The rows that need the decision 0014 allowance, and how many ulps each needs.
///
/// Two rows, one ulp each: `deep`'s first allele is `7.695503150241148E-54` upstream and `...149`
/// here, and `rare-allele`'s is `1.0505352307751092E-12` against `...094`. Every intermediate on
/// the way to `deep` is bit-identical, checked against the reference in the pinned container:
/// `binomialCoefficientLog(1050, 200)` is `507.79546692773414` on both sides,
/// `binomialProbability(1050, 200, p)` agrees for both contamination hypotheses, and
/// `logLikelihoodGivenSomatic(1050, 200)` is `-6.947547001280557` on both. What is left is the
/// final `posteriorProbabilityOfError`, whose `normalizeFromLogToLinearSpace` reaches
/// `StrictMath.exp` here and `Math.exp` upstream, bounded at 1 ulp by htsjdk-rs decision 0025.
/// The other eighteen probability rows, four of which are below `1e-50`, are exact.
const ULPS_APART: [(&str, i64); 2] = [("deep", 1), ("rare-allele", 1)];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/contamination_filter.txt.gz"),
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

/// The record's reference is one base and both alternates are one base: no indel either way.
const REFERENCE_LENGTH: i32 = 1;

fn substitution() -> AlternateAllele {
    AlternateAllele {
        length: 1,
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

struct Case {
    population_frequencies: Option<Vec<f64>>,
    genotypes: Vec<GenotypeData<i32>>,
    contamination: f64,
}

fn case(label: &str) -> Case {
    let mut case = Case {
        population_frequencies: Some(vec![2.0, 3.0]),
        genotypes: vec![genotype(true, &[80, 20, 5]), genotype(false, &[90, 1, 0])],
        contamination: 0.05,
    };
    match label {
        "contamination-default" => case.contamination = DEFAULT_CONTAMINATION,
        "contamination-five-percent" => {}
        "contamination-half" => case.contamination = 0.5,
        "contamination-one" => case.contamination = 1.0,
        "contamination-above-one" => case.contamination = 2.0,
        "contamination-negative" => case.contamination = -0.5,
        "common-allele" => case.population_frequencies = Some(vec![0.5, 0.5]),
        "rare-allele" => case.population_frequencies = Some(vec![9.0, 9.0]),
        "frequency-of-one" => case.population_frequencies = Some(vec![0.0, 0.0]),
        "infinite-popaf" => case.population_frequencies = Some(vec![f64::INFINITY, f64::INFINITY]),
        "short-popaf" => case.population_frequencies = Some(vec![2.0]),
        "long-popaf" => case.population_frequencies = Some(vec![2.0, 3.0, 4.0]),
        "no-popaf" => case.population_frequencies = None,
        "normal-only" => case.genotypes = vec![genotype(false, &[90, 1, 0])],
        "two-tumours" => {
            case.genotypes = vec![
                genotype(true, &[80, 20, 5]),
                genotype(true, &[40, 1, 30]),
                genotype(false, &[90, 1, 0]),
            ]
        }
        "no-alt-reads" => case.genotypes[0] = genotype(true, &[100, 0, 0]),
        "all-alt-reads" => case.genotypes[0] = genotype(true, &[0, 100, 0]),
        "shallow" => case.genotypes[0] = genotype(true, &[4, 2, 1]),
        "deep" => case.genotypes[0] = genotype(true, &[800, 200, 50]),
        "one-alt-read-deep" => case.genotypes[0] = genotype(true, &[999, 1, 0]),
        other => panic!("no case named {other}"),
    }
    case
}

fn answer(label: &str) -> Result<Vec<f64>, ContaminationError> {
    let case = case(label);
    let mut model = SomaticClusteringModel::new(PriorArguments::new(), None);
    // No contamination tables are read, so every sample takes the same estimate.
    let contaminations = vec![case.contamination; case.genotypes.len()];
    contamination_error_probabilities(
        &mut model,
        case.population_frequencies.as_deref(),
        &case.genotypes,
        &contaminations,
        &[substitution(), substitution()],
        REFERENCE_LENGTH,
    )
}

fn printed(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|v| java_double_to_string(*v)).collect();
    format!("[{}]", parts.join(", "))
}

fn same(ours: &[f64], payload: &str, label: &str, complaints: &mut Vec<String>) {
    let mine = printed(ours);
    if mine == payload {
        return;
    }
    let allowance = ULPS_APART
        .iter()
        .find(|(name, _)| *name == label)
        .map(|(_, ulps)| *ulps);
    let Some(ulps) = allowance else {
        complaints.push(format!("{label}: ours {mine}, reference {payload}"));
        return;
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
        if apart > ulps {
            complaints.push(format!("{label}: {apart} ulps apart, allowed {ulps}"));
        }
    }
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 22, "the golden's row count");
    let mut complaints: Vec<String> = Vec::new();
    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "default" => {
                assert_eq!(*label, "contaminationEstimate".to_string());
                assert_eq!(*payload, java_double_to_string(DEFAULT_CONTAMINATION));
            }
            "name" => assert_eq!(
                *payload,
                format!("{FILTER_NAME},NON_SOMATIC,{ANNOTATION},POPAF")
            ),
            "prob" => same(
                &answer(label).expect("answered"),
                payload,
                label,
                &mut complaints,
            ),
            "error" => {
                let error = answer(label).expect_err("refused");
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
    assert!(complaints.is_empty(), "{}", complaints.join("\n"));
}
