//! Conformance for `buildFiltersList` and the mode-dependent argument getters against
//! GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/FilterListByModeDump.java`.
//!
//! # What this suite is for
//!
//!  * **eighteen filters by default, twelve in mitochondrial mode, thirteen in microbial**, and the
//!    slippage filter moves to the end of the microbial list;
//!  * **`NormalArtifactFilter` is in the mitochondrial list**, under a comment saying it is not;
//!  * **the mapping-quality getter remembers its first answer**, because it writes to the field it
//!    reads;
//!  * **and passing the default prior explicitly is indistinguishable from not passing it.**
//!
//! Every row is compared and every row is bit-identical.

use gatk_corpus as corpus;
use gatk_engine::mutect_filter_list::{build_filters_list, FilterArguments};
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/filter_list_by_mode.txt.gz"),
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

fn arguments(mitochondria: bool, microbial: bool) -> FilterArguments {
    FilterArguments {
        mitochondria,
        microbial,
        ..FilterArguments::default()
    }
}

fn mode(label: &str) -> FilterArguments {
    match label {
        "default" => arguments(false, false),
        "mitochondria" => arguments(true, false),
        "microbial" => arguments(false, true),
        "both-modes" => arguments(true, true),
        other => panic!("no mode named {other}"),
    }
}

/// The `argument` rows, each label naming what the dump did to the collection first.
fn argument_value(label: &str, name: &str) -> String {
    match name {
        "minMedianMappingQuality" => {
            let mut collection = match label {
                "mapping-quality-default" => arguments(false, false),
                "mapping-quality-mitochondria" => arguments(true, false),
                "mapping-quality-microbial" => arguments(false, true),
                "mapping-quality-explicit" => FilterArguments {
                    min_median_mapping_quality: 42,
                    ..arguments(false, true)
                },
                // The pair that shows the memoisation: asked before the flag is set, and again
                // after.
                "mapping-quality-asked-first" | "mapping-quality-asked-again" => {
                    let mut remembered = arguments(false, false);
                    remembered.min_median_mapping_quality();
                    if label.ends_with("-again") {
                        remembered.microbial = true;
                    }
                    remembered
                }
                "mapping-quality-flag-set-first" => arguments(false, true),
                other => panic!("no mapping-quality case named {other}"),
            };
            collection.min_median_mapping_quality().to_string()
        }
        "logSnvPrior" | "logIndelPrior" => {
            let collection = match label {
                "priors-default" => arguments(false, false),
                "priors-mitochondria" => arguments(true, false),
                "priors-microbial" => arguments(false, true),
                // Exactly the default, passed explicitly, under mitochondrial mode.
                "priors-explicitly-the-default" => arguments(true, false),
                "priors-explicitly-other" => FilterArguments {
                    log_snv_prior: -12.0,
                    log_indel_prior: -13.0,
                    ..arguments(true, false)
                },
                other => panic!("no prior case named {other}"),
            };
            java_double_to_string(if name == "logSnvPrior" {
                collection.log_snv_prior()
            } else {
                collection.log_indel_prior()
            })
        }
        other => panic!("no argument named {other}"),
    }
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 21, "the golden's row count");
    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "list" => {
                let classes = build_filters_list(&mode(label));
                assert_eq!(
                    format!("{}={}", classes.len(), classes.join(",")),
                    *payload,
                    "list {label}"
                );
            }
            "argument" => {
                let (name, expected) = payload.split_once('=').expect("a value");
                assert_eq!(
                    argument_value(label, name),
                    expected,
                    "argument {label} {name}"
                );
            }
            other => panic!("no row kind {other}"),
        }
    }
}
