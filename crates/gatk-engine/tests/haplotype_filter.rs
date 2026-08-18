//! Conformance for `FilteredHaplotypeFilter` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/FilteredHaplotypeFilterDump.java`.
//!
//! # What this suite is for
//!
//!  * **the interning changes the answer.** Three accumulation rows differ only in whether a
//!    filter's name is the constant or an equal copy of it, and they accumulate `0.1`, `0.9` and
//!    `0.1`;
//!  * **the first pass answers zero to everything**, and learning empties the accumulating map
//!    while the learned one keeps every locus;
//!  * **the record filters itself**, the distance test including zero;
//!  * **a tie in `AF` keeps the first tumour genotype**, which answers `0.2` where the greater
//!    fraction answers `0.8`;
//!  * **and a record with no tumour sample is a refusal**, not a zero.
//!
//! Every row is compared and every row is bit-identical. Nothing here is arithmetic.

use gatk_corpus as corpus;
use gatk_engine::error_probabilities::ErrorType;
use gatk_engine::haplotype_filter::{
    FilterAnswer, FilterIdentity, FilteredHaplotypeFilter, PhasedGenotype,
    ARTIFACT_IN_NORMAL_FILTER_NAME, DEFAULT_MAX_INTRA_HAPLOTYPE_DISTANCE, FILTER_NAME,
    GERMLINE_RISK_FILTER_NAME,
};
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/haplotype_filter.txt.gz"),
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

fn tumour(pgt: Option<&str>, pid: Option<&str>, allele_fraction: f64) -> PhasedGenotype {
    PhasedGenotype {
        tumor: true,
        allele_fractions: vec![allele_fraction],
        phasing_gt: pgt.map(str::to_string),
        phasing_id: pid.map(str::to_string),
    }
}

fn normal() -> PhasedGenotype {
    PhasedGenotype {
        tumor: false,
        allele_fractions: Vec::new(),
        phasing_gt: None,
        phasing_id: None,
    }
}

/// The dump's `record(start, pgt, pid, af)`: one tumour and one normal.
fn record(pgt: Option<&str>, pid: Option<&str>, allele_fraction: f64) -> Vec<PhasedGenotype> {
    vec![tumour(pgt, pid, allele_fraction), normal()]
}

/// The dump's `twoTumours`, both phased against the same `PID`.
fn two_tumours(
    first_af: f64,
    first_pgt: &str,
    second_af: f64,
    second_pgt: &str,
) -> Vec<PhasedGenotype> {
    vec![
        tumour(Some(first_pgt), Some("100_A_C"), first_af),
        tumour(Some(second_pgt), Some("100_A_C"), second_af),
    ]
}

fn answer(name: &str, error_type: ErrorType, probability: f64) -> FilterAnswer {
    // The identity a name has when it IS the constant. `interned` below is the other case.
    let identity = if name == GERMLINE_RISK_FILTER_NAME {
        FilterIdentity::Germline
    } else if name == ARTIFACT_IN_NORMAL_FILTER_NAME {
        FilterIdentity::NormalArtifact
    } else {
        FilterIdentity::Other
    };
    FilterAnswer {
        name: name.to_string(),
        identity,
        error_type,
        probabilities: vec![probability],
    }
}

/// The same name, not the same object: `new String(...)` upstream.
fn copied(name: &str, error_type: ErrorType, probability: f64) -> FilterAnswer {
    FilterAnswer {
        name: name.to_string(),
        identity: FilterIdentity::Other,
        error_type,
        probabilities: vec![probability],
    }
}

fn artifact(name: &str, probability: f64) -> FilterAnswer {
    answer(name, ErrorType::Artifact, probability)
}

/// The filter with the golden's four accumulated loci learned.
fn learned_filter() -> FilteredHaplotypeFilter {
    let mut filter = FilteredHaplotypeFilter::new(100);
    let base = artifact("base_qual", 0.8);
    filter.accumulate_data_for_learning(&[base], &record(Some("0|1"), Some("100_A_C"), 0.3), 100);
    filter.accumulate_data_for_learning(
        &[artifact("base_qual", 0.4)],
        &record(Some("0|1"), Some("100_A_C"), 0.3),
        150,
    );
    filter.accumulate_data_for_learning(
        &[artifact("base_qual", 0.9)],
        &record(Some("0|1"), Some("100_A_C"), 0.3),
        500,
    );
    filter.accumulate_data_for_learning(
        &[artifact("base_qual", 0.7)],
        &record(Some("1|0"), Some("100_A_C"), 0.3),
        120,
    );
    filter
}

/// The two-tumour filter, whose haplotypes carry different probabilities.
fn two_tumour_filter() -> FilteredHaplotypeFilter {
    let mut filter = FilteredHaplotypeFilter::new(100);
    filter.accumulate_data_for_learning(
        &[artifact("base_qual", 0.2)],
        &record(Some("0|1"), Some("100_A_C"), 0.1),
        100,
    );
    filter.accumulate_data_for_learning(
        &[artifact("base_qual", 0.8)],
        &record(Some("1|0"), Some("100_A_C"), 0.9),
        100,
    );
    filter
}

/// One accumulation over a chosen set of filter answers, as the dump's `accumulated` does.
fn accumulation(label: &str) -> FilteredHaplotypeFilter {
    let answers: Vec<FilterAnswer> = match label {
        "artifact-only" => vec![artifact("base_qual", 0.3), artifact("map_qual", 0.7)],
        "non-somatic-excluded" => vec![
            artifact("base_qual", 0.3),
            answer("contamination", ErrorType::NonSomatic, 0.99),
        ],
        "self-excluded" => vec![artifact("base_qual", 0.3), artifact(FILTER_NAME, 0.99)],
        "germline-drops-normal-artifact" => vec![
            answer(GERMLINE_RISK_FILTER_NAME, ErrorType::NonSomatic, 0.5),
            artifact(ARTIFACT_IN_NORMAL_FILTER_NAME, 0.9),
            artifact("base_qual", 0.1),
        ],
        "germline-below-the-threshold" => vec![
            answer(GERMLINE_RISK_FILTER_NAME, ErrorType::NonSomatic, 0.1),
            artifact(ARTIFACT_IN_NORMAL_FILTER_NAME, 0.9),
            artifact("base_qual", 0.1),
        ],
        "germline-at-the-threshold" => vec![
            answer(GERMLINE_RISK_FILTER_NAME, ErrorType::NonSomatic, 0.25),
            artifact(ARTIFACT_IN_NORMAL_FILTER_NAME, 0.9),
            artifact("base_qual", 0.1),
        ],
        // The name says `germline` and is not the constant, so `==` does not match it.
        "germline-name-not-interned" => vec![
            copied(GERMLINE_RISK_FILTER_NAME, ErrorType::NonSomatic, 0.5),
            artifact(ARTIFACT_IN_NORMAL_FILTER_NAME, 0.9),
            artifact("base_qual", 0.1),
        ],
        "normal-artifact-name-not-interned" => vec![
            answer(GERMLINE_RISK_FILTER_NAME, ErrorType::NonSomatic, 0.5),
            copied(ARTIFACT_IN_NORMAL_FILTER_NAME, ErrorType::Artifact, 0.9),
            artifact("base_qual", 0.1),
        ],
        // And this one IS matched, because that comparison is `.equals`.
        "self-name-not-interned" => vec![
            copied(FILTER_NAME, ErrorType::Artifact, 0.99),
            artifact("base_qual", 0.1),
        ],
        other => panic!("no accumulation named {other}"),
    };
    let mut filter = FilteredHaplotypeFilter::new(100);
    filter.accumulate_data_for_learning(&answers, &record(Some("0|1"), Some("100_A_C"), 0.3), 100);
    filter
}

/// `[(locus,probability), ...]` as `ImmutablePair`'s `toString` renders a list of them.
fn render(values: &[(i32, f64)]) -> String {
    let parts: Vec<String> = values
        .iter()
        .map(|(locus, probability)| format!("({locus},{})", java_double_to_string(*probability)))
        .collect();
    format!("[{}]", parts.join(", "))
}

/// The `accumulated`/`learned` rows a filter's map produces, sorted by key as the dump sorts them.
fn map_rows(map: &[(String, Vec<(i32, f64)>)]) -> Vec<String> {
    if map.is_empty() {
        return vec!["(empty)".to_string()];
    }
    let mut keys: Vec<&(String, Vec<(i32, f64)>)> = map.iter().collect();
    keys.sort_by(|a, b| a.0.cmp(&b.0));
    keys.iter()
        .map(|(key, values)| format!("{key}={}", render(values)))
        .collect()
}

fn printed(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|v| java_double_to_string(*v)).collect();
    format!("[{}]", parts.join(", "))
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 32, "the golden's row count");

    let learned = learned_filter();
    // The dump shows the accumulating map, then learns, then shows both.
    let mut after = learned.clone();
    after.learn_parameters_and_clear_accumulated_data();
    let two_tumour = two_tumour_filter();
    let mut two_tumour_learned = two_tumour.clone();
    two_tumour_learned.learn_parameters_and_clear_accumulated_data();
    let mut shared = FilteredHaplotypeFilter::new(100);
    shared.accumulate_data_for_learning(
        &[artifact("base_qual", 0.6)],
        &two_tumours(0.1, "0|1", 0.9, "0|1"),
        100,
    );

    // Every `accumulated`/`learned` row of one label, in the order the dump prints them.
    let mut expected_map_rows: std::collections::HashMap<(&str, &str), Vec<String>> =
        std::collections::HashMap::new();
    expected_map_rows.insert(
        ("accumulated", "three-loci"),
        map_rows(learned.accumulated()),
    );
    expected_map_rows.insert(("learned", "three-loci"), map_rows(after.learned()));
    expected_map_rows.insert(
        ("accumulated", "after-learning"),
        map_rows(after.accumulated()),
    );
    expected_map_rows.insert(
        ("accumulated", "two-tumours"),
        map_rows(two_tumour.accumulated()),
    );
    expected_map_rows.insert(
        ("accumulated", "shared-haplotype"),
        map_rows(shared.accumulated()),
    );
    for label in [
        "artifact-only",
        "non-somatic-excluded",
        "self-excluded",
        "germline-drops-normal-artifact",
        "germline-below-the-threshold",
        "germline-at-the-threshold",
        "germline-name-not-interned",
        "normal-artifact-name-not-interned",
        "self-name-not-interned",
    ] {
        expected_map_rows.insert(
            ("accumulated", label),
            map_rows(accumulation(label).accumulated()),
        );
    }
    // Each label's rows are consumed in order as the golden presents them.
    let mut consumed: std::collections::HashMap<(&str, &str), usize> =
        std::collections::HashMap::new();

    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "default" => {
                assert_eq!(
                    *label,
                    "maxDistanceToFilteredCallOnSameHaplotype".to_string()
                );
                assert_eq!(*payload, DEFAULT_MAX_INTRA_HAPLOTYPE_DISTANCE.to_string());
            }
            "name" => assert_eq!(*payload, format!("{FILTER_NAME},ARTIFACT,none,none")),
            "accumulated" | "learned" => {
                let key = (kind.as_str(), label.as_str());
                let (key, ours) = expected_map_rows
                    .get_key_value(&key)
                    .unwrap_or_else(|| panic!("no map for {kind} {label}"));
                let index = consumed.entry(*key).or_insert(0);
                assert_eq!(ours[*index], *payload, "{kind} {label} row {index}");
                *index += 1;
            }
            "prob" => {
                let (genotypes, start, filter) = match label.as_str() {
                    "first-pass" => (
                        record(Some("0|1"), Some("100_A_C"), 0.3),
                        100,
                        FilteredHaplotypeFilter::new(100),
                    ),
                    "at-its-own-locus" => (
                        record(Some("0|1"), Some("100_A_C"), 0.3),
                        100,
                        after.clone(),
                    ),
                    "within-of-two" => (
                        record(Some("0|1"), Some("100_A_C"), 0.3),
                        150,
                        after.clone(),
                    ),
                    "out-of-range" => (
                        record(Some("0|1"), Some("100_A_C"), 0.3),
                        300,
                        after.clone(),
                    ),
                    "exactly-at-the-distance" => (
                        record(Some("0|1"), Some("100_A_C"), 0.3),
                        200,
                        after.clone(),
                    ),
                    "one-past-the-distance" => (
                        record(Some("0|1"), Some("100_A_C"), 0.3),
                        201,
                        after.clone(),
                    ),
                    "other-haplotype" => (
                        record(Some("1|0"), Some("100_A_C"), 0.3),
                        120,
                        after.clone(),
                    ),
                    "unknown-haplotype" => (
                        record(Some("1|1"), Some("999_A_C"), 0.3),
                        100,
                        after.clone(),
                    ),
                    "no-pgt" => (record(None, Some("100_A_C"), 0.3), 100, after.clone()),
                    "no-pid" => (record(Some("0|1"), None, 0.3), 100, after.clone()),
                    "greatest-af-wins" => (
                        two_tumours(0.1, "0|1", 0.9, "1|0"),
                        100,
                        two_tumour_learned.clone(),
                    ),
                    "tie-keeps-the-first" => (
                        two_tumours(0.5, "0|1", 0.5, "1|0"),
                        100,
                        two_tumour_learned.clone(),
                    ),
                    other => panic!("no probability case named {other}"),
                };
                let ours = filter
                    .error_probabilities(&genotypes, start, 1)
                    .expect("answered");
                assert_eq!(printed(&ours), *payload, "prob {label}");
            }
            "error" => {
                assert_eq!(*label, "no-tumour-sample".to_string());
                let error = after
                    .error_probabilities(&[normal()], 100, 1)
                    .expect_err("refused");
                assert_eq!(*payload, format!("{}:{}", error.class(), error.message()));
            }
            other => panic!("no row kind {other}"),
        }
    }
}
