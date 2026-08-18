//! Conformance for the four remaining hard filters against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/RemainingHardFiltersDump.java`.
//!
//! # What this suite is for
//!
//!  * **none of the four can fire with the default arguments**, and the defaults are compared as
//!    values rather than assumed;
//!  * **`DuplicatedAltReadFilter`'s answer is as long as its annotation**, not as long as the
//!    record: four counts on a two-alternate record answer four probabilities;
//!  * **the two base classes disagree on a missing required annotation**, in the same golden: the
//!    per-allele one answers `[]`, which `ErrorProbabilities` drops, and the per-site one answers
//!    `[0.0, 0.0]`, which it counts;
//!  * **`MinAlleleFractionFilter` answers "not an artifact" from an absence**, an allele with no
//!    data being `orElse(1.0)`;
//!  * **and `PanelOfNormalsFilter` reads presence rather than value**, so `PON=false` is filtered.
//!
//! Every row is compared and every row is bit-identical. Nothing here is arithmetic beyond one
//! division.

use gatk_corpus as corpus;
use gatk_engine::allele_filter::alt_data_by_allele;
use gatk_engine::allele_filter::GenotypeData;
use gatk_engine::mutect_engine::round_finite_precision_errors;
use gatk_engine::mutect_filter_list::FILTERS;
use gatk_engine::mutect_hard_filters::{
    duplicated_alt_read_artifacts, error_probability, min_allele_fraction_artifacts,
    n_ratio_is_artifact, panel_of_normals_is_artifact, DEFAULT_MAX_N_RATIO, DEFAULT_MIN_AF,
    DEFAULT_MIN_UNIQUE_ALT_READS,
};
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/remaining_hard_filters.txt.gz"),
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

/// The record's three alleles, and its two alternates.
const ALLELES: [&str; 3] = ["A", "C", "G"];

/// One genotype: its allele depths, and the allele fractions it carries (empty when it has none).
fn genotype(tumor: bool, allele_depths: &[i32], fractions: &[f64]) -> GenotypeData<f64> {
    GenotypeData {
        tumor,
        allele_depths: allele_depths.to_vec(),
        values: fractions.to_vec(),
    }
}

/// The base record's two samples.
fn base_genotypes() -> Vec<GenotypeData<f64>> {
    vec![
        genotype(true, &[80, 20, 5], &[0.2, 0.05]),
        genotype(false, &[90, 1, 0], &[0.01, 0.0]),
    ]
}

/// `Mutect2AlleleFilter.errorProbabilities`: an empty list when a required annotation is missing,
/// and one rounded probability per element of whatever the filter answered otherwise.
fn allele_probabilities(present: bool, artifacts: Vec<bool>) -> Vec<f64> {
    if !present {
        return Vec::new();
    }
    artifacts
        .into_iter()
        .map(|artifact| round_finite_precision_errors(error_probability(artifact)))
        .collect()
}

/// `Mutect2VariantFilter.errorProbabilities`: one probability copied to every alternate allele, and
/// that probability is `0.0` when a required annotation is missing.
fn site_probabilities(present: bool, artifact: bool, alternate_count: usize) -> Vec<f64> {
    let probability = round_finite_precision_errors(if present {
        error_probability(artifact)
    } else {
        0.0
    });
    vec![probability; alternate_count]
}

/// `sumADsOverSamples(vc, true, true)`: both samples, indexed by the record's allele count.
fn all_allele_depths(genotypes: &[GenotypeData<f64>]) -> Vec<i32> {
    let mut totals = vec![0; ALLELES.len()];
    for genotype in genotypes {
        for (total, depth) in totals.iter_mut().zip(&genotype.allele_depths) {
            *total += depth;
        }
    }
    totals
}

/// `getAltDataByAllele(vc, g -> g.hasExtendedAttribute(AF) && isTumor(g), ...)`.
fn fractions_by_alt_allele(genotypes: &[GenotypeData<f64>]) -> Vec<Vec<f64>> {
    let alleles: Vec<String> = ALLELES.iter().map(|a| a.to_string()).collect();
    alt_data_by_allele(&alleles, genotypes, |g| g.tumor && !g.values.is_empty())
        .into_iter()
        .map(|(_, values)| values)
        .collect()
}

fn probabilities(label: &str) -> Vec<f64> {
    match label {
        // `DuplicatedAltReadFilter`, whose answer is as long as its annotation.
        "duplicate-default" => allele_probabilities(
            true,
            duplicated_alt_read_artifacts(&[3, 1], DEFAULT_MIN_UNIQUE_ALT_READS),
        ),
        "duplicate-threshold-two" => {
            allele_probabilities(true, duplicated_alt_read_artifacts(&[3, 1], 2))
        }
        "duplicate-short-list" => {
            allele_probabilities(true, duplicated_alt_read_artifacts(&[1], 2))
        }
        "duplicate-long-list" => {
            allele_probabilities(true, duplicated_alt_read_artifacts(&[1, 1, 1, 1], 2))
        }
        "duplicate-no-annotation" => allele_probabilities(false, Vec::new()),

        // `NRatioFilter`, over both samples' depths.
        "nratio-default" => site_probabilities(
            true,
            n_ratio_is_artifact(
                &all_allele_depths(&base_genotypes()),
                4,
                DEFAULT_MAX_N_RATIO,
            ),
            2,
        ),
        "nratio-threshold-half" => site_probabilities(
            true,
            n_ratio_is_artifact(&all_allele_depths(&base_genotypes()), 4, 0.5),
            2,
        ),
        "nratio-no-alt-reads" => {
            let genotypes = vec![
                genotype(true, &[80, 0, 0], &[0.0, 0.0]),
                genotype(false, &[90, 0, 0], &[0.0, 0.0]),
            ];
            site_probabilities(
                true,
                n_ratio_is_artifact(&all_allele_depths(&genotypes), 4, 0.5),
                2,
            )
        }
        "nratio-no-annotation" => site_probabilities(false, false, 2),
        "nratio-at-the-threshold" => site_probabilities(
            true,
            n_ratio_is_artifact(&all_allele_depths(&base_genotypes()), 13, 0.5),
            2,
        ),

        // `MinAlleleFractionFilter`, which requires no annotation at all.
        "minaf-default" => allele_probabilities(
            true,
            min_allele_fraction_artifacts(
                &fractions_by_alt_allele(&base_genotypes()),
                DEFAULT_MIN_AF,
            ),
        ),
        "minaf-threshold-tenth" => allele_probabilities(
            true,
            min_allele_fraction_artifacts(&fractions_by_alt_allele(&base_genotypes()), 0.1),
        ),
        "minaf-no-allele-fraction" => {
            let genotypes = vec![genotype(true, &[80, 20, 5], &[])];
            allele_probabilities(
                true,
                min_allele_fraction_artifacts(&fractions_by_alt_allele(&genotypes), 0.1),
            )
        }
        "minaf-normal-is-low" => {
            let genotypes = vec![
                genotype(true, &[80, 20, 5], &[0.5, 0.5]),
                genotype(false, &[90, 1, 0], &[0.001, 0.001]),
            ];
            allele_probabilities(
                true,
                min_allele_fraction_artifacts(&fractions_by_alt_allele(&genotypes), 0.1),
            )
        }
        "minaf-full-length-list" => {
            let genotypes = vec![genotype(true, &[80, 20, 5], &[0.9, 0.05, 0.9])];
            allele_probabilities(
                true,
                min_allele_fraction_artifacts(&fractions_by_alt_allele(&genotypes), 0.1),
            )
        }
        "minaf-at-the-threshold" => {
            let genotypes = vec![genotype(true, &[80, 20, 5], &[0.1, 0.1])];
            allele_probabilities(
                true,
                min_allele_fraction_artifacts(&fractions_by_alt_allele(&genotypes), 0.1),
            )
        }

        // `PanelOfNormalsFilter`, which reads presence rather than value.
        "pon-absent" => site_probabilities(true, panel_of_normals_is_artifact(false), 2),
        "pon-present" | "pon-false" | "pon-empty-string" => {
            site_probabilities(true, panel_of_normals_is_artifact(true), 2)
        }
        other => panic!("no case named {other}"),
    }
}

fn printed(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|v| java_double_to_string(*v)).collect();
    format!("[{}]", parts.join(", "))
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 27, "the golden's row count");
    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "default" => {
                let ours = match label.as_str() {
                    "uniqueAltReadCount" => DEFAULT_MIN_UNIQUE_ALT_READS.to_string(),
                    "nRatio" => java_double_to_string(DEFAULT_MAX_N_RATIO),
                    "minAf" => java_double_to_string(DEFAULT_MIN_AF),
                    other => panic!("no default named {other}"),
                };
                assert_eq!(*payload, ours, "default {label}");
            }
            "name" => {
                let filter = FILTERS
                    .iter()
                    .find(|f| f.class == label)
                    .unwrap_or_else(|| panic!("no filter class {label}"));
                // `phredScaledPosteriorAnnotationName` is empty for every hard filter.
                assert_eq!(
                    *payload,
                    format!("{},{},none", filter.filter_name, filter.error_type.name()),
                    "name {label}"
                );
            }
            "filter" => assert_eq!(printed(&probabilities(label)), *payload, "filter {label}"),
            other => panic!("no row kind {other}"),
        }
    }
}
