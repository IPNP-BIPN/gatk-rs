//! Conformance for the EM iteration and the quantile initialisation against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/SomaticClusteringLearnDump.java`.
//!
//! # What this suite is for
//!
//!  * **`learn` is sequential within an epoch**, so the same data in two orders learn two shapes
//!    that `%.2f` prints identically and the probes tell apart;
//!  * **the floors are applied every step**, so a cluster with no responsibility cannot move and one
//!    with none of its own lands at `mean = 0.990`;
//!  * **the initialisation splits peaks until the BIC stops improving**, at most five times;
//!  * **the prior map is rewritten only when there is a callable-site count**;
//!  * **and every number in `clusteringMetadata` is formatted before anyone sees it**, HALF_UP.
//!
//! The whole learning path runs through `exp`, whose bit-exact transcription is withdrawn under
//! htsjdk-rs decision 0014: the binomial density is `exp` of a saddle point, and every
//! responsibility is a normalisation of logs. Divergences are therefore bounded rather than absent,
//! and every row that needs the allowance is named in `ULPS_APART` with the size it needs.

use gatk_corpus as corpus;
use gatk_engine::allele_fraction_cluster::{
    learn_beta_binomial, learn_binomial, AlleleFractionCluster, BetaDistributionShape, Datum,
};
use gatk_engine::java_format::format_decimals;
use gatk_engine::somatic_clustering_model::{
    binomial_probability, AlternateAllele, PriorArguments, SomaticClusteringModel,
};
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/somatic_clustering_learn.txt.gz"),
    )
}

/// The probe counts the dump uses.
const PROBES: [(i32, i32); 3] = [(10, 5), (100, 3), (100, 50)];

/// The rows that need the decision 0014 allowance, and how many ulps each needs.
const ULPS_APART: [(&str, i64); 0] = [];

fn datum(total: i32, alt: i32) -> Datum {
    Datum::new(20.0, 0.0, 0.0, alt, total, 0)
}

fn data(counts: &[(i32, i32)]) -> Vec<Datum> {
    counts
        .iter()
        .map(|(total, alt)| datum(*total, *alt))
        .collect()
}

fn describe(cluster: &AlleleFractionCluster) -> String {
    match cluster {
        AlleleFractionCluster::BetaBinomial(shape) => format!(
            "alpha = {}, beta = {}",
            format_decimals(shape.alpha(), 2),
            format_decimals(shape.beta(), 2)
        ),
        AlleleFractionCluster::Binomial(shape) => format!(
            "mean = {}",
            format_decimals(shape.alpha() / (shape.alpha() + shape.beta()), 3)
        ),
    }
}

/// One `learn` call, printed as the dump prints it.
fn learned(
    label: &str,
    cluster: AlleleFractionCluster,
    data: &[Datum],
    responsibilities: &[f64],
) -> Vec<String> {
    let learned = match cluster {
        AlleleFractionCluster::BetaBinomial(shape) => AlleleFractionCluster::BetaBinomial(
            learn_beta_binomial(shape, data, responsibilities).expect("a shape"),
        ),
        AlleleFractionCluster::Binomial(_) => AlleleFractionCluster::Binomial(
            learn_binomial(data, responsibilities).expect("a shape"),
        ),
    };
    let mut rows = vec![format!("learn\t{label}\t{}", describe(&learned))];
    for (total, alt) in PROBES {
        let value = learned.log_likelihood(total, alt).expect("in range");
        rows.push(format!(
            "probe\t{label}-{total},{alt}\t{}",
            java_double_to_string(value)
        ));
    }
    rows
}

fn flat() -> AlleleFractionCluster {
    AlleleFractionCluster::beta_binomial(BetaDistributionShape::FLAT)
}

fn high_af() -> AlleleFractionCluster {
    AlleleFractionCluster::beta_binomial(BetaDistributionShape::new(10.0, 1.0).expect("a shape"))
}

fn binomial(mean: f64) -> AlleleFractionCluster {
    AlleleFractionCluster::binomial(mean).expect("a shape")
}

/// A model given the data one datum at a time, then told to learn.
fn model(counts: &[Datum], callable_sites: Option<f64>) -> SomaticClusteringModel {
    let mut model = SomaticClusteringModel::new(PriorArguments::new(), callable_sites);
    feed(&mut model, counts);
    model.learn_and_clear_accumulated_data().expect("learned");
    model
}

fn feed(model: &mut SomaticClusteringModel, counts: &[Datum]) {
    for datum in counts {
        let mut ads = [datum.total_count() - datum.alt_count(), datum.alt_count()];
        model
            .record(
                &mut ads,
                &[datum.tumor_log_odds()],
                &[datum.artifact_prob()],
                &[0.0],
                &[AlternateAllele {
                    length: 1,
                    symbolic: false,
                }],
                1,
            )
            .expect("recorded");
    }
}

/// The `metadata`, `prior` and `artifactprior` rows one learned model produces.
fn model_rows(label: &str, model: &mut SomaticClusteringModel) -> Vec<String> {
    let mut rows: Vec<String> = model
        .clustering_metadata()
        .into_iter()
        .map(|(key, value)| format!("metadata\t{label}\t{key}={value}"))
        .collect();
    for length in [0, 1, -1, 40] {
        rows.push(format!(
            "prior\t{label}-{length}\t{}",
            java_double_to_string(model.log_prior_of_somatic_variant(length))
        ));
    }
    rows.push(format!(
        "artifactprior\t{label}\t{}",
        java_double_to_string(model.log_prior_of_variant_versus_artifact())
    ));
    rows
}

fn ours() -> Vec<String> {
    let mut rows = Vec::new();
    for n in [10, 100] {
        for k in [0, 1, 5, 10] {
            if k <= n {
                for f in [0.0, 0.01, 0.1, 0.5, 0.99, 1.0] {
                    rows.push(format!(
                        "binomprob\t{n},{k},{}\t{}",
                        java_double_to_string(f),
                        java_double_to_string(binomial_probability(n, k, f))
                    ));
                }
            }
        }
    }

    let clonal = data(&[(100, 50), (100, 48), (100, 52), (100, 51)]);
    let ones = vec![1.0; clonal.len()];
    rows.extend(learned("betabinomial-flat-clonal", flat(), &clonal, &ones));
    rows.extend(learned(
        "betabinomial-highaf-clonal",
        high_af(),
        &clonal,
        &ones,
    ));
    rows.extend(learned("binomial-clonal", binomial(0.5), &clonal, &ones));

    let mut reversed = clonal.clone();
    reversed.reverse();
    rows.extend(learned(
        "betabinomial-flat-reversed",
        flat(),
        &reversed,
        &ones,
    ));
    rows.extend(learned(
        "binomial-reversed",
        binomial(0.5),
        &reversed,
        &ones,
    ));

    let halved = vec![0.5; clonal.len()];
    rows.extend(learned(
        "betabinomial-flat-halved",
        flat(),
        &clonal,
        &halved,
    ));
    rows.extend(learned("binomial-halved", binomial(0.5), &clonal, &halved));

    let zeroes = vec![0.0; clonal.len()];
    rows.extend(learned("betabinomial-flat-zero", flat(), &clonal, &zeroes));
    rows.extend(learned("binomial-zero", binomial(0.5), &clonal, &zeroes));

    rows.extend(learned("betabinomial-flat-empty", flat(), &[], &[]));
    rows.extend(learned("binomial-empty", binomial(0.5), &[], &[]));

    let subclonal = data(&[(100, 5), (100, 7), (100, 4), (100, 6)]);
    let subclonal_ones = vec![1.0; subclonal.len()];
    rows.extend(learned(
        "betabinomial-flat-subclonal",
        flat(),
        &subclonal,
        &subclonal_ones,
    ));
    rows.extend(learned(
        "binomial-subclonal",
        binomial(0.5),
        &subclonal,
        &subclonal_ones,
    ));

    let bimodal = data(&[
        (100, 45),
        (100, 48),
        (100, 50),
        (100, 14),
        (100, 16),
        (100, 15),
        (100, 47),
        (100, 13),
    ]);
    for (label, data, callable) in [
        ("clonal", clonal.clone(), Some(10000.0)),
        ("subclonal", subclonal.clone(), Some(10000.0)),
        ("bimodal", bimodal, Some(10000.0)),
        ("one-datum", data(&[(100, 50)]), Some(10000.0)),
        ("no-data", Vec::new(), Some(10000.0)),
        ("clonal-no-callable-sites", clonal.clone(), None),
        ("clonal-zero-callable-sites", clonal.clone(), Some(0.0)),
    ] {
        let mut learned_model = model(&data, callable);
        rows.extend(model_rows(label, &mut learned_model));
    }

    // Learned twice, with new data between the rounds.
    let mut twice = model(&clonal, Some(10000.0));
    rows.extend(model_rows("twice-first", &mut twice));
    feed(&mut twice, &subclonal);
    twice.learn_and_clear_accumulated_data().expect("learned");
    rows.extend(model_rows("twice-second", &mut twice));
    rows
}

#[test]
fn every_row_matches_the_golden() {
    let text = golden();
    let expected: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(expected.len(), 357, "the golden's row count");

    let mine = ours();
    assert_eq!(mine.len(), expected.len(), "every row is accounted for");

    for (ours, theirs) in mine.iter().zip(&expected) {
        if ours == theirs {
            continue;
        }
        let label = theirs.split('\t').nth(1).unwrap_or_default();
        let allowance = ULPS_APART
            .iter()
            .find(|(name, _)| *name == label)
            .map(|(_, ulps)| *ulps);
        let allowed = allowance.unwrap_or_else(|| panic!("{label}: {ours} against {theirs}"));
        let value: f64 = ours
            .rsplit('\t')
            .next()
            .and_then(|text| text.parse().ok())
            .unwrap_or_else(|| panic!("{label}: {ours} against {theirs}"));
        let reference: f64 = theirs
            .rsplit('\t')
            .next()
            .and_then(|text| text.parse().ok())
            .unwrap_or_else(|| panic!("{label}: {ours} against {theirs}"));
        let ulps = ((value.to_bits() as i64) - (reference.to_bits() as i64)).abs();
        assert!(
            ulps <= allowed,
            "{label}: {ours} against {theirs}, {ulps} ulps"
        );
    }
}
