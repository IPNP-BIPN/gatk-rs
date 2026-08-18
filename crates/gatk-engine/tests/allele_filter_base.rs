//! Conformance for the allele filter's shared machinery against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/AlleleFilterBaseDump.java`.
//!
//! # What this suite is for
//!
//!  * **the gather zips two iterators and stops at the shorter one**, so a short list shifts the
//!    whole gather rather than shortening one allele's list;
//!  * **the alternate gather still starts at the caller's first element**, so full-length data gives
//!    the first alternate the reference's value;
//!  * **the weighted median sorts the caller's list in place** and leans to the lower half;
//!  * **a short AD array is an out-of-bounds, not a skip**;
//!  * **and a missing annotation is an empty list, not a zero.**
//!
//! The `evidence` rows go through the clustering model's `exp`, whose bit-exact transcription is
//! withdrawn under htsjdk-rs decision 0014, so they are compared with an allowance. One row needs
//! it and needs **one ulp**: `strong`, whose second allele is `0.74812961299286` upstream and
//! `0.7481296129928601` here. `ULPS_APART` names it, so a second row cannot start diverging
//! quietly. Everything else in this suite is exact.

use gatk_corpus as corpus;
use gatk_engine::allele_filter::{
    alt_data_by_allele, data_by_allele, sum_ads_over_samples, tumor_evidence_error_probabilities,
    weighted_median_posterior_probability, AlleleDepthTooShort, GenotypeData,
    TUMOR_EVIDENCE_ANNOTATION, TUMOR_EVIDENCE_ERROR_TYPE, TUMOR_EVIDENCE_FILTER_NAME,
};
use gatk_engine::somatic_clustering_model::{
    AlternateAllele, PriorArguments, SomaticClusteringModel,
};
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/allele_filter_base.txt.gz"),
    )
}

/// The rows that need the decision 0014 allowance, and how many ulps each needs.
///
/// `strong` alone, and one ulp: its second allele's probability is `0.74812961299286` upstream and
/// `0.7481296129928601` here. Every other `evidence` row is bit-identical, including the two whose
/// probabilities sit against 1.0, where an ulp would not show.
const ULPS_APART: [(&str, i64); 1] = [("strong", 1)];

/// The record every gather row is built from: two tumour samples and a normal.
fn genotypes(values: [Vec<i32>; 3]) -> Vec<GenotypeData<i32>> {
    let [first, second, normal] = values;
    vec![
        GenotypeData {
            tumor: true,
            allele_depths: vec![80, 20, 5],
            values: first,
        },
        GenotypeData {
            tumor: true,
            allele_depths: vec![70, 30, 10],
            values: second,
        },
        GenotypeData {
            tumor: false,
            allele_depths: vec![99, 1, 0],
            values: normal,
        },
    ]
}

fn alleles() -> Vec<String> {
    vec!["A".to_string(), "C".to_string(), "G".to_string()]
}

/// One gather, printed the way the dump prints a `LinkedHashMap`'s entries.
fn gather_rows(kind: &str, label: &str, gathered: &[(String, Vec<i32>)]) -> Vec<String> {
    gathered
        .iter()
        .map(|(allele, values)| {
            let list = values
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<String>>()
                .join(", ");
            format!("{kind}\t{label}\t{allele}=[{list}]")
        })
        .collect()
}

/// `ImmutablePair.toString`, which is `(left,right)` with no space.
fn pair_list(pairs: &[(i32, f64)]) -> String {
    let inner = pairs
        .iter()
        .map(|(depth, posterior)| format!("({depth},{})", java_double_to_string(*posterior)))
        .collect::<Vec<String>>()
        .join(", ");
    format!("[{inner}]")
}

fn median_rows(label: &str, input: &[(i32, f64)]) -> Vec<String> {
    let mut mutable = input.to_vec();
    let answer = weighted_median_posterior_probability(&mut mutable);
    vec![
        format!("median\t{label}\t{}", java_double_to_string(answer)),
        format!("sorted\t{label}\t{}", pair_list(&mutable)),
    ]
}

fn ads_row(label: &str, genotypes: &[GenotypeData<i32>], tumor: bool, normal: bool) -> String {
    match sum_ads_over_samples(3, genotypes, tumor, normal) {
        Ok(totals) => {
            let list = totals
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<String>>()
                .join(", ");
            format!("ads\t{label}\t[{list}]")
        }
        Err(AlleleDepthTooShort { index, length }) => format!(
            "error\tads-{label}\tjava.lang.ArrayIndexOutOfBoundsException:Index {index} out of \
             bounds for length {length}"
        ),
    }
}

/// The dump's two alternate alleles, both single-base substitutions.
fn alternates() -> Vec<AlternateAllele> {
    vec![
        AlternateAllele {
            length: 1,
            symbolic: false,
        },
        AlternateAllele {
            length: 1,
            symbolic: false,
        },
    ]
}

fn evidence_row(label: &str, tumor_log_odds: Option<&[f64]>) -> String {
    // A model built the way the dump builds its engine: default arguments, no stats table.
    let mut model = SomaticClusteringModel::new(PriorArguments::new(), None);
    let probabilities = tumor_evidence_error_probabilities(
        &mut model,
        tumor_log_odds,
        &genotypes([vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]]),
        &alternates(),
        1,
    );
    let list = probabilities
        .iter()
        .map(|value| java_double_to_string(*value))
        .collect::<Vec<String>>()
        .join(", ");
    format!("evidence\t{label}\t[{list}]")
}

/// `getTumorLogOdds`, which is `log10ToLog` of the annotation.
fn log_odds(log10: &[f64]) -> Vec<f64> {
    log10
        .iter()
        .map(|value| value * std::f64::consts::LN_10)
        .collect()
}

fn ours() -> Vec<String> {
    let mut rows = Vec::new();
    let full = genotypes([vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]]);
    rows.extend(gather_rows(
        "byallele",
        "three-alleles-three-values",
        &data_by_allele(&alleles(), &full, |g| g.tumor),
    ));
    rows.extend(gather_rows(
        "altbyallele",
        "three-alleles-three-values",
        &alt_data_by_allele(&alleles(), &full, |g| g.tumor),
    ));

    let short = genotypes([vec![1], vec![4, 5, 6], vec![7, 8, 9]]);
    rows.extend(gather_rows(
        "byallele",
        "short-list",
        &data_by_allele(&alleles(), &short, |g| g.tumor),
    ));
    rows.extend(gather_rows(
        "altbyallele",
        "short-list",
        &alt_data_by_allele(&alleles(), &short, |g| g.tumor),
    ));

    // The overlong run keeps one tumour sample and the normal.
    let overlong = vec![
        GenotypeData {
            tumor: true,
            allele_depths: vec![80, 20, 5],
            values: vec![1, 2, 3, 4, 5],
        },
        GenotypeData {
            tumor: false,
            allele_depths: vec![99, 1, 0],
            values: vec![7, 8, 9],
        },
    ];
    rows.extend(gather_rows(
        "byallele",
        "long-list",
        &data_by_allele(&alleles(), &overlong, |g| g.tumor),
    ));

    let normal_only = vec![GenotypeData {
        tumor: false,
        allele_depths: vec![99, 1, 0],
        values: vec![7, 8, 9],
    }];
    rows.extend(gather_rows(
        "byallele",
        "no-tumor",
        &data_by_allele(&alleles(), &normal_only, |g| g.tumor),
    ));

    for (label, pairs) in [
        ("even-split", vec![(10, 0.1), (10, 0.9)]),
        ("one-dominant", vec![(1, 0.1), (99, 0.9)]),
        ("out-of-order", vec![(5, 0.9), (5, 0.2), (5, 0.5)]),
        ("all-equal", vec![(4, 0.5), (4, 0.5), (4, 0.5)]),
        ("zero-depth", vec![(0, 0.3), (0, 0.7)]),
        ("single", vec![(7, 0.42)]),
        ("empty", Vec::new()),
    ] {
        rows.extend(median_rows(label, &pairs));
    }

    rows.push(ads_row("tumor-only", &full, true, false));
    rows.push(ads_row("normal-only", &full, false, true));
    rows.push(ads_row("both", &full, true, true));
    rows.push(ads_row("neither", &full, false, false));
    rows.push(ads_row(
        "short-ad",
        &[
            GenotypeData {
                tumor: true,
                allele_depths: vec![80, 20],
                values: vec![1, 2, 3],
            },
            GenotypeData {
                tumor: false,
                allele_depths: vec![99, 1, 0],
                values: vec![7, 8, 9],
            },
        ],
        true,
        false,
    ));

    rows.push(format!(
        "name\ttumor-evidence\t{TUMOR_EVIDENCE_FILTER_NAME},{},{TUMOR_EVIDENCE_ANNOTATION}",
        TUMOR_EVIDENCE_ERROR_TYPE.name()
    ));
    rows.push(evidence_row("strong", Some(&log_odds(&[20.0, 6.0]))));
    rows.push(evidence_row("weak", Some(&log_odds(&[1.0, 0.5]))));
    rows.push(evidence_row("negative", Some(&log_odds(&[-3.0, -0.1]))));
    rows.push(evidence_row("no-tlod", None));
    rows.push(evidence_row("short-tlod", Some(&log_odds(&[20.0]))));
    rows
}

#[test]
fn every_row_matches_the_golden() {
    let text = golden();
    let expected: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(expected.len(), 41, "the golden's row count");

    let mine = ours();
    assert_eq!(mine.len(), expected.len(), "every row is accounted for");

    for (ours, theirs) in mine.iter().zip(&expected) {
        if ours == theirs {
            continue;
        }
        let label = theirs.split('\t').nth(1).unwrap_or_default();
        let allowed = ULPS_APART
            .iter()
            .find(|(name, _)| *name == label)
            .map(|(_, ulps)| *ulps)
            .unwrap_or_else(|| panic!("{label}: {ours} against {theirs}"));
        // A list of probabilities: compare entry by entry.
        let mine_values = ours.rsplit('\t').next().expect("a value");
        let their_values = theirs.rsplit('\t').next().expect("a value");
        let parse = |text: &str| -> Vec<f64> {
            text.trim_matches(['[', ']'])
                .split(", ")
                .filter(|piece| !piece.is_empty())
                .map(|piece| piece.parse().expect("a double"))
                .collect()
        };
        let (mine_values, their_values) = (parse(mine_values), parse(their_values));
        assert_eq!(mine_values.len(), their_values.len(), "{label}");
        for (value, reference) in mine_values.iter().zip(&their_values) {
            let ulps = ((value.to_bits() as i64) - (reference.to_bits() as i64)).abs();
            assert!(
                ulps <= allowed,
                "{label}: {ours} against {theirs}, {ulps} ulps"
            );
        }
    }
}
