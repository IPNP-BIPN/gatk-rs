//! Conformance for `applyFiltersAndAccumulateOutputStats` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/ApplyFiltersDump.java`.
//!
//! # What this suite is for
//!
//!  * **the symbolic-allele removal is applied twice**, and a record whose symbolic allele comes
//!    first is `FAIL`ed although its one real allele passed;
//!  * **a record with no per-allele filter at all is a `NoSuchElementException`**;
//!  * **a filter can fire without being named**, the reporting floor being
//!    `min(maxErrorProb, 0.1)`;
//!  * **`SITE` is a placeholder**, so a record can be `PASS` with a filtered allele beside it;
//!  * **and the phred-scaled annotation is written when the filter did not fire**, and not written
//!    when its required annotations are missing.
//!
//! Every row is compared and every row is bit-identical.

use gatk_corpus as corpus;
use gatk_engine::apply_filters::{apply_filters, ApplyError, FilterAnswer, FilterKind};
use gatk_engine::somatic_clustering_model::AlternateAllele;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/apply_filters.txt.gz"),
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

fn real() -> AlternateAllele {
    AlternateAllele {
        length: 1,
        symbolic: false,
    }
}

fn symbolic() -> AlternateAllele {
    AlternateAllele {
        length: 0,
        symbolic: true,
    }
}

fn allele(name: &str, probabilities: &[f64]) -> FilterAnswer {
    FilterAnswer {
        name: name.to_string(),
        kind: FilterKind::PerAllele,
        probabilities: probabilities.to_vec(),
        annotation: None,
        required_annotations_present: true,
    }
}

fn site(name: &str, probability: f64) -> FilterAnswer {
    FilterAnswer {
        name: name.to_string(),
        kind: FilterKind::PerSite,
        probabilities: vec![probability],
        annotation: None,
        required_annotations_present: true,
    }
}

fn annotated(name: &str, probability: f64, annotation: &str, required: bool) -> FilterAnswer {
    FilterAnswer {
        annotation: Some(annotation.to_string()),
        required_annotations_present: required,
        ..site(name, probability)
    }
}

/// The dump's records. A per-site filter's list is one entry per alternate upstream, because
/// `Mutect2VariantFilter.errorProbabilities` copies it; the port takes the one probability.
struct Case {
    answers: Vec<FilterAnswer>,
    alternates: Vec<AlternateAllele>,
    threshold: f64,
}

fn case(label: &str) -> Case {
    let biallelic = vec![real()];
    let triallelic = vec![real(), real()];
    let passing = allele("base_qual", &[0.0]);
    match label {
        // One alternate, and a filter that answered for two: the extra string is never consumed.
        "everything-passes" => Case {
            answers: vec![allele("base_qual", &[0.1, 0.1])],
            alternates: biallelic,
            threshold: 0.5,
        },
        "one-allele-filtered" => Case {
            answers: vec![allele("base_qual", &[0.9, 0.1])],
            alternates: triallelic,
            threshold: 0.5,
        },
        "every-allele-same-filter" => Case {
            answers: vec![allele("base_qual", &[0.9, 0.9])],
            alternates: triallelic,
            threshold: 0.5,
        },
        "every-allele-different-filters" => Case {
            answers: vec![
                allele("base_qual", &[0.9, 0.1]),
                allele("map_qual", &[0.1, 0.9]),
            ],
            alternates: triallelic,
            threshold: 0.5,
        },
        "one-allele-two-filters" => Case {
            answers: vec![
                allele("base_qual", &[0.9, 0.1]),
                allele("map_qual", &[0.9, 0.1]),
            ],
            alternates: triallelic,
            threshold: 0.5,
        },
        "site-and-allele-filters" => Case {
            answers: vec![
                allele("base_qual", &[0.9, 0.1]),
                allele("map_qual", &[0.1, 0.9]),
                site("germline", 0.95),
            ],
            alternates: triallelic,
            threshold: 0.5,
        },
        "only-site-filters" => Case {
            answers: vec![site("germline", 0.9)],
            alternates: biallelic,
            threshold: 0.5,
        },
        "below-the-reporting-floor" => Case {
            answers: vec![passing, site("germline", 0.99), site("contamination", 0.05)],
            alternates: biallelic,
            threshold: 0.01,
        },
        "below-the-floor-alone" => Case {
            answers: vec![passing, site("contamination", 0.05)],
            alternates: biallelic,
            threshold: 0.01,
        },
        "annotation-written" => Case {
            answers: vec![passing, annotated("germline", 0.9, "GERMQ", true)],
            alternates: biallelic,
            threshold: 0.5,
        },
        // `Mutect2VariantFilter.errorProbabilities` has ALREADY turned the filter's 0.9 into 0.0,
        // its required annotation being absent; the second check here only decides the annotation.
        "annotation-required-annotation-missing" => Case {
            answers: vec![passing, annotated("germline", 0.0, "GERMQ", false)],
            alternates: biallelic,
            threshold: 0.5,
        },
        "annotation-without-the-filter" => Case {
            answers: vec![passing, annotated("germline", 0.1, "GERMQ", true)],
            alternates: biallelic,
            threshold: 0.5,
        },
        // `ErrorProbabilities` already removed the symbolic allele's entry, so the list is short.
        "symbolic-last" => Case {
            answers: vec![allele("base_qual", &[0.9])],
            alternates: vec![real(), symbolic()],
            threshold: 0.5,
        },
        "symbolic-first" => Case {
            answers: vec![allele("base_qual", &[0.1])],
            alternates: vec![symbolic(), real()],
            threshold: 0.5,
        },
        "empty-list" => Case {
            answers: vec![allele("base_qual", &[]), allele("map_qual", &[0.9, 0.9])],
            alternates: triallelic,
            threshold: 0.5,
        },
        "no-filters" => Case {
            answers: vec![],
            alternates: triallelic,
            threshold: 0.5,
        },
        "threshold-zero" => Case {
            answers: vec![allele("base_qual", &[0.0])],
            alternates: biallelic,
            threshold: 0.0,
        },
        "threshold-one" => Case {
            answers: vec![allele("base_qual", &[1.0])],
            alternates: biallelic,
            threshold: 1.0,
        },
        other => panic!("no case named {other}"),
    }
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 20, "the golden's row count");
    for (kind, label, payload) in &rows {
        let case = case(label);
        let answer = apply_filters(&case.answers, &case.alternates, case.threshold);
        match kind.as_str() {
            "applied" => {
                let applied = answer.expect("applied");
                // The dump prints the FILTER set sorted, and `PASS` when it is empty.
                let mut names = applied.filters.clone();
                names.sort();
                let column = if names.is_empty() {
                    "PASS".to_string()
                } else {
                    names.join(";")
                };
                assert_eq!(
                    format!("{column}|{}", applied.as_filter_status),
                    *payload,
                    "applied {label}"
                );
            }
            "attribute" => {
                let applied = answer.expect("applied");
                let (key, value) = payload.split_once('=').expect("a value");
                let ours = applied
                    .annotations
                    .iter()
                    .find(|(name, _)| name == key)
                    .unwrap_or_else(|| panic!("no annotation {key} for {label}"));
                assert_eq!(ours.1.to_string(), value, "attribute {label} {key}");
            }
            "error" => {
                let error: ApplyError = answer.expect_err("refused");
                assert_eq!(
                    format!("{}:{}", error.class(), error.message()),
                    *payload,
                    "error {label}"
                );
            }
            other => panic!("no row kind {other}"),
        }
    }
}
