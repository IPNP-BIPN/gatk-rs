//! Conformance for `GermlineFilter.calculateErrorProbability` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/GermlineWrapperDump.java`.
//!
//! # What this suite is for
//!
//!  * **the hom-alt switch is worth six orders of magnitude on one hundredth of allele fraction**;
//!  * **two early returns bracket the population frequency**, and a `POPAF` of zero is a hard one;
//!  * **one index, three conventions**: `POPAF` per alternate, the depths with `+ 1`, the weighted
//!    fractions per alternate;
//!  * **a record with no `NLOD` uses zero** rather than skipping the normal;
//!  * **and the two helpers weight the same depths for different purposes.**
//!
//! Every row reaches the clustering model's `exp`, so a divergence would be named in
//! [`ULPS_APART`].

use gatk_corpus as corpus;
use gatk_engine::allele_filter::GenotypeData;
use gatk_engine::germline_filter::{
    compute_minor_allele_fraction, germline_error_probabilities, weighted_average_of_tumor_afs,
    ANNOTATION, FILTER_NAME,
};
use gatk_engine::somatic_clustering_model::{
    AlternateAllele, PriorArguments, SomaticClusteringModel,
};
use gatk_engine::tsv_table::java_double_to_string;

/// The rows that need the decision 0014 allowance, and how many ulps each needs.
///
/// Two rows, one ulp each. Both go through `germlineProbability`, whose
/// `normalizeFromLogToLinearSpace` reaches `StrictMath.exp` here and `Math.exp` upstream, bounded at
/// 1 ulp by htsjdk-rs decision 0025. The other eleven probabilities are exact, including the two
/// hard brackets and the pair either side of the hom-alt switch, which is what the suite is for.
const ULPS_APART: [(&str, i64); 2] = [("rare-allele", 1), ("normal-says-somatic", 1)];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/germline_wrapper.txt.gz"),
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

fn substitution() -> AlternateAllele {
    AlternateAllele {
        length: 1,
        symbolic: false,
    }
}

fn tumour(depths: &[i32]) -> GenotypeData<i32> {
    GenotypeData {
        tumor: true,
        allele_depths: depths.to_vec(),
        values: Vec::new(),
    }
}

fn normal(depths: &[i32]) -> GenotypeData<i32> {
    GenotypeData {
        tumor: false,
        allele_depths: depths.to_vec(),
        values: Vec::new(),
    }
}

/// One record's arguments, as the dump built them.
struct Case {
    tumor_log_10_odds: Option<Vec<f64>>,
    population_af: Option<Vec<f64>>,
    normal_log_10_odds: Option<Vec<f64>>,
    genotypes: Vec<GenotypeData<i32>>,
    allele_fractions: Vec<Vec<f64>>,
    alternates: Vec<AlternateAllele>,
}

/// The dump's `record(...)`: one tumour and one normal, one alternate.
fn biallelic(
    population_af: f64,
    normal_log_10_odds: Option<f64>,
    tumor_depths: &[i32],
    allele_fraction: f64,
) -> Case {
    Case {
        tumor_log_10_odds: Some(vec![20.0]),
        population_af: Some(vec![population_af]),
        normal_log_10_odds: normal_log_10_odds.map(|value| vec![value]),
        genotypes: vec![tumour(tumor_depths), normal(&[90, 1])],
        allele_fractions: vec![vec![allele_fraction], Vec::new()],
        alternates: vec![substitution()],
    }
}

fn case(label: &str) -> Case {
    match label {
        "rare-allele" | "one-tumour" => biallelic(6.0, Some(-4.0), &[80, 20], 0.2),
        "common-allele" => biallelic(1.0, Some(-4.0), &[80, 20], 0.2),
        "fixed-in-the-population" => biallelic(0.0, Some(-4.0), &[80, 20], 0.2),
        "below-the-epsilon" => biallelic(400.0, Some(-4.0), &[80, 20], 0.2),
        "hom-alt-off" => biallelic(6.0, Some(-4.0), &[10, 80], 0.89),
        "hom-alt-on" => biallelic(6.0, Some(-4.0), &[10, 80], 0.9),
        "normal-says-germline" => biallelic(6.0, Some(-8.0), &[80, 20], 0.2),
        "normal-says-somatic" => biallelic(6.0, Some(8.0), &[80, 20], 0.2),
        "no-nlod" => biallelic(6.0, None, &[80, 20], 0.2),
        "no-depth" => biallelic(6.0, Some(-4.0), &[0, 0], 0.0),
        "second-allele-wins" => Case {
            tumor_log_10_odds: Some(vec![6.0, 20.0]),
            population_af: Some(vec![6.0, 2.0]),
            normal_log_10_odds: Some(vec![-4.0, -4.0]),
            genotypes: vec![tumour(&[60, 20, 20]), normal(&[90, 1, 0])],
            allele_fractions: vec![vec![0.2, 0.2], Vec::new()],
            alternates: vec![substitution(), substitution()],
        },
        "no-tlod" => Case {
            tumor_log_10_odds: None,
            ..biallelic(6.0, Some(-4.0), &[80, 20], 0.2)
        },
        "no-popaf" => Case {
            population_af: None,
            ..biallelic(6.0, Some(-4.0), &[80, 20], 0.2)
        },
        "two-tumours" => Case {
            genotypes: vec![tumour(&[80, 20]), tumour(&[20, 30]), normal(&[90, 1])],
            allele_fractions: vec![vec![0.2], vec![0.6], Vec::new()],
            ..biallelic(6.0, Some(-4.0), &[80, 20], 0.2)
        },
        other => panic!("no case named {other}"),
    }
}

/// The filter's answer, as `Mutect2VariantFilter` copies it to every alternate.
fn probabilities(label: &str) -> Vec<f64> {
    let case = case(label);
    let mut model = SomaticClusteringModel::new(PriorArguments::new(), None);
    // No tumour segmentation table, so every sample's minor allele fraction is 0.5.
    let minor = vec![0.5; case.genotypes.len()];
    germline_error_probabilities(
        &mut model,
        case.tumor_log_10_odds.as_deref(),
        case.population_af.as_deref(),
        case.normal_log_10_odds.as_deref(),
        &case.genotypes,
        &case.allele_fractions,
        &minor,
        &case.alternates,
        1,
    )
    .expect("answered")
}

fn printed(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|v| java_double_to_string(*v)).collect();
    format!("[{}]", parts.join(", "))
}

/// Compare as text, with the decision 0014 allowance named row by row; every divergence is
/// reported rather than panicked on, so one run names them all.
fn same(ours: &str, payload: &str, label: &str, complaints: &mut Vec<String>) {
    if ours == payload {
        return;
    }
    let allowance = ULPS_APART
        .iter()
        .find(|(name, _)| *name == label)
        .map(|(_, ulps)| *ulps);
    let ours_values: Vec<&str> = ours
        .trim_matches(|c| c == '[' || c == ']')
        .split(", ")
        .collect();
    let theirs: Vec<&str> = payload
        .trim_matches(|c| c == '[' || c == ']')
        .split(", ")
        .collect();
    if ours_values.len() != theirs.len() {
        complaints.push(format!("{label}: ours {ours}, reference {payload}"));
        return;
    }
    let mut worst = 0i64;
    for (ours, expected) in ours_values.iter().zip(&theirs) {
        if ours == expected {
            continue;
        }
        let (a, b): (f64, f64) = (
            ours.parse().expect("a double"),
            expected.parse().expect("a double"),
        );
        worst = worst.max((a.to_bits() as i64 - b.to_bits() as i64).abs());
    }
    match allowance {
        Some(ulps) if worst <= ulps => {}
        _ => complaints.push(format!(
            "{label}: {worst} ulps apart (ours {ours}, reference {payload})"
        )),
    }
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 18, "the golden's row count");
    let mut complaints: Vec<String> = Vec::new();
    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "name" => assert_eq!(
                *payload,
                format!("{FILTER_NAME},NON_SOMATIC,{ANNOTATION},TLOD,POPAF")
            ),
            "prob" => same(
                &printed(&probabilities(label)),
                payload,
                label,
                &mut complaints,
            ),
            "weightedaf" => {
                let case = case(label);
                let ours = weighted_average_of_tumor_afs(
                    &case.genotypes,
                    &case.allele_fractions,
                    case.alternates.len(),
                )
                .expect("averaged");
                // `Arrays.toString(double[])`, which is the same shape as a list's.
                same(&printed(&ours), payload, label, &mut complaints);
            }
            "maf" => {
                let case = case(label);
                let minor = vec![0.5; case.genotypes.len()];
                let allele_counts = gatk_engine::allele_filter::sum_ads_over_samples(
                    case.alternates.len() + 1,
                    &case.genotypes,
                    true,
                    false,
                )
                .expect("summed");
                let ours = compute_minor_allele_fraction(&case.genotypes, &minor, &allele_counts);
                same(
                    &java_double_to_string(ours),
                    payload,
                    label,
                    &mut complaints,
                );
            }
            other => panic!("no row kind {other}"),
        }
    }
    assert!(complaints.is_empty(), "{}", complaints.join("\n"));
}
