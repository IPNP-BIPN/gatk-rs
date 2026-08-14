//! Conformance for `Mutect2FilteringEngine` as constructed against GATK 4.6.2.0, compared as the
//! sample partition and the constants a fresh engine holds.
//!
//! Golden from `tools/readfilter-conformance/MutectEngineConstructionDump.java`.
//!
//! # What this suite is for
//!
//!  * **every sample not named normal is a tumour sample**, including one under a differently-cased
//!    key and one not in the VCF at all;
//!  * **the threshold starts at the argument collection's default**;
//!  * **a clustering model with no data still has priors**, and only the SNV one carries `log(1/3)`.
//!
//! The two posterior rows of the golden are not recomputed here: they go through the somatic
//! clustering model, which is not ported. What is compared is everything the port has.

use gatk_corpus as corpus;
use gatk_engine::mutect_engine::{
    default_log_indel_prior, default_log_prior_of_variant_versus_artifact, default_log_snv_prior,
    is_normal, is_tumor, log_prior_of_somatic_variant, normal_samples,
    DEFAULT_INITIAL_POSTERIOR_THRESHOLD,
};
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/mutect_engine_construction.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn expected(text: &str, kind: &str, label: &str) -> String {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no {kind} {label}"))[1]
        .to_string()
}

/// The header the dump built: two normal samples and a third under a differently-cased key.
fn header_lines() -> Vec<(String, String)> {
    vec![
        ("normal_sample".to_string(), "N1".to_string()),
        ("normal_sample".to_string(), "N2".to_string()),
        ("Normal_Sample".to_string(), "N3".to_string()),
    ]
}

#[test]
fn every_sample_matches_the_golden() {
    let text = golden();
    let normals = normal_samples(&header_lines());
    for sample in ["T1", "N1", "N2", "N3", "never-mentioned"] {
        let ours = if is_normal(&normals, sample) {
            "normal"
        } else {
            "tumour"
        };
        assert_eq!(ours, expected(&text, "sample", sample), "{sample}");
    }
}

#[test]
fn every_sample_not_named_normal_is_a_tumour_sample() {
    let text = golden();
    // Declared under a key differing only in case.
    assert_eq!(expected(&text, "sample", "N3"), "tumour");
    // Not in the VCF at all.
    assert_eq!(expected(&text, "sample", "never-mentioned"), "tumour");
    let normals = normal_samples(&header_lines());
    assert!(is_tumor(&normals, "N3"));
    assert!(is_tumor(&normals, "never-mentioned"));
    assert_eq!(normals, vec!["N1".to_string(), "N2".to_string()]);
}

#[test]
fn the_threshold_starts_at_the_default() {
    let text = golden();
    assert_eq!(
        java_double_to_string(DEFAULT_INITIAL_POSTERIOR_THRESHOLD),
        expected(&text, "value", "initial-threshold")
    );
}

#[test]
fn a_model_with_no_data_still_has_priors() {
    let text = golden();
    assert_eq!(
        java_double_to_string(default_log_prior_of_variant_versus_artifact()),
        expected(&text, "value", "log-prior-variant-versus-artifact")
    );
    // A SNV: the SNV prior plus log(1/3).
    assert_eq!(
        java_double_to_string(log_prior_of_somatic_variant(0)),
        expected(&text, "value", "log-somatic-prior-snp")
    );
    // A two-base deletion: the indel prior, and nothing added.
    assert_eq!(
        java_double_to_string(log_prior_of_somatic_variant(-2)),
        expected(&text, "value", "log-somatic-prior-indel")
    );
    // The log(1/3) is the whole of the difference between the priors and their defaults.
    let snp: f64 = expected(&text, "value", "log-somatic-prior-snp")
        .parse()
        .expect("a double");
    let indel: f64 = expected(&text, "value", "log-somatic-prior-indel")
        .parse()
        .expect("a double");
    assert_eq!(indel, default_log_indel_prior());
    assert!((snp - default_log_snv_prior() - (1.0f64 / 3.0).ln()).abs() < 1.0e-15);
    // And the two priors are nearer each other than the two defaults they come from.
    assert!((snp - indel).abs() < (default_log_snv_prior() - default_log_indel_prior()).abs());
}
