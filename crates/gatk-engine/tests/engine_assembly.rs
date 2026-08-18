//! Conformance for the assembled filtering engine against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/EngineAssemblyDump.java`.
//!
//! # What this suite is for
//!
//!  * **seventeen of the eighteen filters answer a fully annotated record**, `StrictStrandBiasFilter`
//!    being switched off by its default and dropped;
//!  * **an empty list, a zero and a NaN are three different answers**;
//!  * **the mode reaches the arithmetic**: the same record's tumour-evidence probability differs
//!    between default and mitochondrial mode;
//!  * **the per-type maximum comes first and the independence product second**;
//!  * **and a site-level filter can fire while every allele passes.**
//!
//! Every row is compared. The allowances are named in [`ULPS_APART`], all of them the
//! `StrictMath.exp` seam of htsjdk-rs decision 0025.

use gatk_corpus as corpus;
use gatk_engine::accumulate_data::AccumulationAllele;
use gatk_engine::allele_filter::GenotypeData;
use gatk_engine::apply_filters::{apply_filters, FilterKind};
use gatk_engine::error_probabilities::{by_type, combined, kept, ErrorType};
use gatk_engine::filtering_engine::{error_probabilities_by_filter, EngineArguments, Record};
use gatk_engine::haplotype_filter::FilteredHaplotypeFilter;
use gatk_engine::mutect_filter_list::FilterArguments;
use gatk_engine::somatic_clustering_model::{
    AlternateAllele, PriorArguments, SomaticClusteringModel,
};
use gatk_engine::strand_artifact_filter::LearnedParameters;
use gatk_engine::tsv_table::java_double_to_string;

/// The rows that need the decision 0014 allowance, and how many ulps each needs.
///
/// Two rows, one ulp each, and the same number in both: the first alternate's tumour-evidence
/// probability on the records that carry two alternates, `2.970295757934691E-14` upstream and
/// `...6914` here. `probabilityOfSequencingError` reaches `StrictMath.exp` here and `Math.exp`
/// upstream, bounded at 1 ulp by htsjdk-rs decision 0025. The biallelic records' tumour-evidence
/// probabilities, both error-type maxima, every combined probability and all eight applied results
/// are exact.
const ULPS_APART: [(&str, i64); 2] = [
    ("two-alternates/TumorEvidenceFilter", 1),
    ("with-non-ref/TumorEvidenceFilter", 1),
];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/engine_assembly.txt.gz"),
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

fn real() -> AccumulationAllele {
    AccumulationAllele {
        allele: AlternateAllele {
            length: 1,
            symbolic: false,
        },
        non_ref: false,
    }
}

fn non_ref() -> AccumulationAllele {
    AccumulationAllele {
        allele: AlternateAllele {
            length: 0,
            symbolic: true,
        },
        non_ref: true,
    }
}

/// The three parallel lists one call to the dump's `genotypes(alleleCount)` produces.
struct Samples {
    genotypes: Vec<GenotypeData<i32>>,
    allele_fractions: Vec<Vec<f64>>,
    phasing: Vec<(Option<String>, Option<String>)>,
}

/// The dump's `genotypes(alleleCount)`: one tumour and one normal, both hom-ref calls.
fn genotypes(allele_count: usize) -> Samples {
    let mut tumour = vec![0; allele_count];
    let mut normal = vec![0; allele_count];
    tumour[0] = 80;
    normal[0] = 90;
    for index in 1..allele_count {
        tumour[index] = 30 - 10 * index as i32;
        normal[index] = 1;
    }
    let fractions: Vec<f64> = (0..allele_count - 1)
        .map(|index| 0.2 / (index as f64 + 1.0))
        .collect();
    Samples {
        genotypes: vec![
            GenotypeData {
                tumor: true,
                allele_depths: tumour,
                values: Vec::new(),
            },
            GenotypeData {
                tumor: false,
                allele_depths: normal,
                values: Vec::new(),
            },
        ],
        allele_fractions: vec![fractions, Vec::new()],
        phasing: vec![
            (Some("0|1".to_string()), Some("100_A_C".to_string())),
            (None, None),
        ],
    }
}

fn repeat_double(value: f64, count: usize) -> Vec<f64> {
    vec![value; count]
}

fn repeat_int(value: i32, count: usize) -> Vec<i32> {
    vec![value; count]
}

/// The dump's `annotated(alleles, alternateCount)`.
fn annotated(alternates: Vec<AccumulationAllele>) -> Record {
    let alternate_count = alternates.len();
    let allele_count = alternate_count + 1;
    let samples = genotypes(allele_count);
    let mut table = String::from("40,40");
    for _ in 0..alternate_count {
        table.push_str("|10,10");
    }
    Record {
        start: 100,
        reference_length: 1,
        alternates,
        genotypes: samples.genotypes,
        allele_fractions: samples.allele_fractions,
        phasing: samples.phasing,
        tumor_log_10_odds: Some(repeat_double(20.0, alternate_count)),
        normal_artifact_log_10_odds: Some(repeat_double(2.0, alternate_count)),
        normal_log_10_odds: None,
        population_af: Some(repeat_double(6.0, alternate_count)),
        median_base_quality: Some(repeat_int(30, allele_count)),
        median_mapping_quality: Some(repeat_int(60, allele_count)),
        median_fragment_length: Some(repeat_int(300, allele_count)),
        median_read_position: Some(repeat_int(25, alternate_count)),
        unique_alt_read_count: Some(repeat_int(8, alternate_count)),
        strand_bias_table: Some(table),
        n_count: Some(0),
        event_count_in_region: Some(1),
        event_count_in_haplotype: Some(1),
        repeats_per_allele: Some(
            repeat_int(10, allele_count)
                .into_iter()
                .map(|value| value.to_string())
                .collect(),
        ),
        repeat_unit: Some("A".to_string()),
        in_panel_of_normals: false,
        indel_lengths: None,
    }
}

fn case(label: &str) -> (Record, bool) {
    match label {
        "fully-annotated" => (annotated(vec![real()]), false),
        "fully-annotated-mitochondria" => (annotated(vec![real()]), true),
        "two-alternates" => (annotated(vec![real(), real()]), false),
        "with-non-ref" => (annotated(vec![real(), non_ref()]), false),
        "bare" => {
            let samples = genotypes(2);
            (
                Record {
                    start: 100,
                    reference_length: 1,
                    alternates: vec![real()],
                    genotypes: samples.genotypes,
                    allele_fractions: samples.allele_fractions,
                    phasing: samples.phasing,
                    ..Record::default()
                },
                false,
            )
        }
        "in-panel-of-normals" => {
            let mut record = annotated(vec![real()]);
            record.in_panel_of_normals = true;
            (record, false)
        }
        "weak-evidence" => {
            let mut record = annotated(vec![real()]);
            record.tumor_log_10_odds = Some(vec![0.5]);
            (record, false)
        }
        "poor-base-quality" => {
            let mut record = annotated(vec![real()]);
            record.median_base_quality = Some(vec![30, 2]);
            (record, false)
        }
        other => panic!("no case named {other}"),
    }
}

/// The engine's answers for one label, in construction order.
fn answers(label: &str) -> Vec<gatk_engine::filtering_engine::EngineAnswer> {
    let (record, mitochondria) = case(label);
    let arguments = EngineArguments {
        list: FilterArguments {
            mitochondria,
            ..FilterArguments::default()
        },
        ..EngineArguments::default()
    };
    // The mode reaches the priors too: `getLogSnvPrior()` answers the mitochondrial value when the
    // flag is set and the argument is still the default.
    let priors = PriorArguments {
        mitochondria,
        ..PriorArguments::new()
    };
    let mut model = SomaticClusteringModel::new(priors, None);
    let haplotype = FilteredHaplotypeFilter::new(100);
    let strand = LearnedParameters::default();
    error_probabilities_by_filter(&mut model, &haplotype, &strand, &arguments, &record)
        .expect("answered")
}

fn printed(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|v| java_double_to_string(*v)).collect();
    format!("[{}]", parts.join(", "))
}

fn same(ours: &str, payload: &str, what: &str, complaints: &mut Vec<String>) {
    if ours == payload {
        return;
    }
    let allowance = ULPS_APART
        .iter()
        .find(|(name, _)| *name == what)
        .map(|(_, ulps)| *ulps);
    let mine: Vec<&str> = ours
        .trim_matches(|c| c == '[' || c == ']')
        .split(", ")
        .collect();
    let theirs: Vec<&str> = payload
        .trim_matches(|c| c == '[' || c == ']')
        .split(", ")
        .collect();
    if mine.len() != theirs.len() {
        complaints.push(format!("{what}: ours {ours}, reference {payload}"));
        return;
    }
    let mut worst = 0i64;
    for (a, b) in mine.iter().zip(&theirs) {
        if a == b {
            continue;
        }
        let (x, y): (f64, f64) = (a.parse().unwrap_or(f64::NAN), b.parse().unwrap_or(f64::NAN));
        worst = worst.max((x.to_bits() as i64 - y.to_bits() as i64).abs());
    }
    match allowance {
        Some(ulps) if worst <= ulps => {}
        _ => complaints.push(format!(
            "{what}: {worst} ulps apart (ours {ours}, reference {payload})"
        )),
    }
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 155, "the golden's row count");
    let mut complaints: Vec<String> = Vec::new();
    let mut consumed: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    for (kind, label, payload) in &rows {
        let answers = answers(label);
        match kind.as_str() {
            "filter" => {
                let (class, values) = payload.split_once('=').expect("a value");
                let index = consumed.entry(label.clone()).or_insert(0);
                let ours = answers
                    .get(*index)
                    .unwrap_or_else(|| panic!("{label}: no answer {index}"));
                assert_eq!(ours.class, class, "{label} answer {index}");
                same(
                    &printed(&ours.probabilities),
                    values,
                    &format!("{label}/{class}"),
                    &mut complaints,
                );
                *index += 1;
            }
            "type" => {
                let (name, values) = payload.split_once('=').expect("a value");
                let error_type = match name {
                    "ARTIFACT" => ErrorType::Artifact,
                    "NON_SOMATIC" => ErrorType::NonSomatic,
                    other => panic!("no error type {other}"),
                };
                let all: Vec<_> = answers.iter().map(|a| a.as_error_probability()).collect();
                let ours = by_type(&kept(&all), error_type).expect("combined");
                same(
                    &printed(&ours),
                    values,
                    &format!("{label}/{name}"),
                    &mut complaints,
                );
            }
            "combined" => {
                let all: Vec<_> = answers.iter().map(|a| a.as_error_probability()).collect();
                let ours = combined(&all).expect("combined");
                same(
                    &printed(&ours),
                    payload,
                    &format!("{label}/combined"),
                    &mut complaints,
                );
            }
            "applied" => {
                let (record, _) = case(label);
                let applied: Vec<_> = answers
                    .iter()
                    .map(|answer| {
                        let annotation = match answer.kind {
                            FilterKind::PerSite => Some(answer.name.to_string()),
                            FilterKind::PerAllele => None,
                        };
                        answer.as_applied(annotation, true)
                    })
                    .collect();
                let alleles: Vec<AlternateAllele> =
                    record.alternates.iter().map(|a| a.allele).collect();
                // The threshold is the initial posterior threshold, which no pass has moved.
                let result = apply_filters(&applied, &alleles, 0.1).expect("applied");
                let mut names = result.filters.clone();
                names.sort();
                let column = if names.is_empty() {
                    "PASS".to_string()
                } else {
                    names.join(";")
                };
                let ours = format!("{column}|{}", result.as_filter_status);
                if ours != *payload {
                    complaints.push(format!("applied {label}: ours {ours}, reference {payload}"));
                }
            }
            other => panic!("no row kind {other}"),
        }
    }
    assert!(complaints.is_empty(), "{}", complaints.join("\n"));
}
