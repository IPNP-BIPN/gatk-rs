//! Conformance for `CalculateGenotypePosteriors` against GATK 4.6.2.0, compared as the prior and
//! every recalled genotype of every site of every run.
//!
//! Golden from `tools/readfilter-conformance/CalculateGenotypePosteriorsDump.java`.
//!
//! Only the population-prior half is ported. No pedigree is passed in the measured runs, so
//! `FamilyLikelihoods` never runs; it is a second algorithm and gets its own brick.
//!
//! # What this suite is for
//!
//!  * **MLEAC being preferred over AC**, and a panel with no MLEAC falling through to it;
//!  * **the reference count being AN minus the alternates**;
//!  * **the flat-prior fallback**, and the ten-sample rule that leads to it;
//!  * **the reference samples entering as chromosomes, and only where the panel is silent**;
//!  * **the SNP pseudocount being chosen by allele length**;
//!  * **and the genotypes being recalled from the posteriors**.

use gatk_corpus as corpus;
use gatk_tools::calculate_genotype_posteriors::{
    apply, dirichlet_multinomial, get_dirichlet_prior, gls_to_pls, Allele, Genotype, Options,
    Record,
};
use std::collections::BTreeMap;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/calculate_genotype_posteriors.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{name}")),
    )
}

fn info(column: &str) -> BTreeMap<String, String> {
    column
        .split(';')
        .filter_map(|part| part.split_once('='))
        .map(|(key, value)| (key.to_string(), value.to_string()))
        .collect()
}

fn alleles(reference: &str, alternates: &str) -> Vec<Allele> {
    let mut out = vec![Allele {
        bases: reference.to_string(),
        is_ref: true,
    }];
    for alternate in alternates.split(',') {
        out.push(Allele {
            bases: alternate.to_string(),
            is_ref: false,
        });
    }
    out
}

/// The records of a VCF the golden carries, genotypes included when there are any.
fn records(vcf: &str) -> Vec<Record> {
    let mut samples: Vec<String> = Vec::new();
    let mut out = Vec::new();
    for line in vcf.lines() {
        if line.starts_with("#CHROM") {
            samples = line.split('\t').skip(9).map(str::to_string).collect();
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        let alleles = alleles(columns[3], columns[4]);
        let mut genotypes = Vec::new();
        if columns.len() > 8 {
            let keys: Vec<&str> = columns[8].split(':').collect();
            for (index, sample) in samples.iter().enumerate() {
                let values: Vec<&str> = columns[9 + index].split(':').collect();
                let get = |key: &str| {
                    keys.iter()
                        .position(|name| *name == key)
                        .map(|position| values[position])
                };
                genotypes.push(Genotype {
                    sample: sample.clone(),
                    alleles: get("GT")
                        .expect("a GT")
                        .split(['/', '|'])
                        .map(|value| value.parse().expect("an allele index"))
                        .collect(),
                    depth: get("DP").map(|value| value.parse().expect("a depth")),
                    likelihoods: get("PL").map(|value| {
                        value
                            .split(',')
                            .map(|part| part.parse().expect("a likelihood"))
                            .collect()
                    }),
                    posteriors: None,
                });
            }
        }
        out.push(Record {
            id: columns[2].to_string(),
            start: columns[1].parse().expect("a position"),
            alleles,
            attributes: info(columns[7]),
            genotypes,
        });
    }
    out
}

/// One sample's column: the recalled `GT`, the `GQ` and the `PP`.
type Called = (Vec<usize>, Option<i32>, Option<Vec<i32>>);

/// What one run wrote at one site: the `PG`, and every sample's column.
#[derive(Debug, PartialEq)]
struct Row {
    prior: Option<Vec<i32>>,
    genotypes: Vec<Called>,
}

fn measured(text: &str, label: &str) -> Vec<Row> {
    let body = section(text, "out", label);
    let mut samples = 0;
    let mut out = Vec::new();
    for line in body.lines() {
        if line.starts_with("#CHROM") {
            samples = line.split('\t').skip(9).count();
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        let attributes = info(columns[7]);
        let keys: Vec<&str> = columns[8].split(':').collect();
        let mut genotypes = Vec::new();
        for index in 0..samples {
            let values: Vec<&str> = columns[9 + index].split(':').collect();
            let get = |key: &str| {
                keys.iter()
                    .position(|name| *name == key)
                    .map(|position| values[position])
            };
            genotypes.push((
                get("GT")
                    .expect("a GT")
                    .split(['/', '|'])
                    .map(|value| value.parse().expect("an allele index"))
                    .collect::<Vec<usize>>(),
                get("GQ").map(|value| value.parse().expect("a GQ")),
                get("PP").map(|value| {
                    value
                        .split(',')
                        .map(|part| part.parse().expect("a posterior"))
                        .collect::<Vec<i32>>()
                }),
            ));
        }
        out.push(Row {
            prior: attributes.get("PG").map(|value| {
                value
                    .split(',')
                    .map(|part| part.parse().expect("a prior"))
                    .collect()
            }),
            genotypes,
        });
    }
    out
}

/// The panel with every MLEAC stripped, which is what the dump wrote as `ac-only.vcf`.
fn panel_without_mleac(panel: &str) -> String {
    panel
        .lines()
        .filter(|line| !line.starts_with("##INFO=<ID=MLEAC"))
        .map(|line| {
            let mut kept = line.to_string();
            for part in line.split('\t').nth(7).unwrap_or("").split(';') {
                if part.starts_with("MLEAC=") {
                    kept = kept.replace(&format!(";{part}"), "");
                }
            }
            kept
        })
        .collect::<Vec<String>>()
        .join("\n")
}

/// label, resources, reference samples, options, skip.
fn runs(
    panel: &[Record],
    ac_only: &[Record],
) -> Vec<(&'static str, Vec<Record>, i32, Options, bool)> {
    let default = Options::default();
    vec![
        ("no-panel", vec![], 0, default, false),
        ("panel", panel.to_vec(), 0, default, false),
        ("ref-samples", vec![], 1000, default, false),
        ("panel-and-ref", panel.to_vec(), 1000, default, false),
        (
            "default-to-ac",
            panel.to_vec(),
            0,
            Options {
                use_mleac: false,
                ..default
            },
            false,
        ),
        (
            "ignore-input",
            panel.to_vec(),
            0,
            Options {
                use_input_samples_allele_counts: false,
                ..default
            },
            false,
        ),
        (
            "discovered-off",
            panel.to_vec(),
            1000,
            Options {
                ignore_input_samples_for_missing_resources: true,
                ..default
            },
            false,
        ),
        (
            "flat-indels",
            panel.to_vec(),
            0,
            Options {
                use_flat_priors_for_indels: true,
                ..default
            },
            false,
        ),
        (
            "priors",
            panel.to_vec(),
            0,
            Options {
                snp_prior_dirichlet: 0.01,
                indel_prior_dirichlet: 0.0001,
                ..default
            },
            false,
        ),
        ("panel-ac-only", ac_only.to_vec(), 0, default, false),
    ]
}

fn produced(input: &[Record], resources: &[Record], num_ref: i32, options: &Options) -> Vec<Row> {
    input
        .iter()
        .map(|record| {
            let posteriors =
                apply(record, resources, num_ref, options, false).expect("population priors");
            Row {
                prior: posteriors.prior,
                genotypes: posteriors
                    .genotypes
                    .into_iter()
                    .map(|called| (called.alleles, called.gq, called.posteriors))
                    .collect(),
            }
        })
        .collect()
}

#[test]
fn every_prior_and_recalled_genotype_matches_the_golden() {
    let text = golden();
    let input = records(&section(&text, "vcf", "input"));
    let panel = records(&section(&text, "vcf", "panel"));
    let ac_only = records(&panel_without_mleac(&section(&text, "vcf", "panel")));

    let mut compared = 0;
    for (label, resources, num_ref, options, _) in runs(&panel, &ac_only) {
        assert_eq!(
            produced(&input, &resources, num_ref, &options),
            measured(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 10, "the runs that applied population priors");
}

/// The same panel, read two ways: MLEAC=20 gives one prior and AC=200 another, and a panel with no
/// MLEAC at all falls through to AC and matches the second exactly.
#[test]
fn mleac_is_preferred_and_ac_is_the_fallback() {
    let text = golden();
    let panel = records(&section(&text, "vcf", "panel"));
    assert_eq!(
        panel[0].attributes.get("MLEAC").map(String::as_str),
        Some("20")
    );
    assert_eq!(
        panel[0].attributes.get("AC").map(String::as_str),
        Some("200")
    );

    let with_mleac = measured(&text, "panel")[0].prior.clone().expect("a prior");
    let with_ac = measured(&text, "default-to-ac")[0]
        .prior
        .clone()
        .expect("a prior");
    let no_mleac = measured(&text, "panel-ac-only")[0]
        .prior
        .clone()
        .expect("a prior");
    assert_ne!(with_mleac, with_ac);
    assert_eq!(with_ac, no_mleac, "no MLEAC is the same as preferring AC");

    // 20 alternate chromosomes out of 2000, plus the input's own three of each.
    assert_eq!(
        gls_to_pls(&get_dirichlet_prior(
            &[1980.0 + 3.0 + 1e-3, 20.0 + 3.0 + 1e-3],
            2,
            false
        )),
        with_mleac
    );
    // And 200 out of 2000.
    assert_eq!(
        gls_to_pls(&get_dirichlet_prior(
            &[1800.0 + 3.0 + 1e-3, 200.0 + 3.0 + 1e-3],
            2,
            false
        )),
        with_ac
    );
}

/// Nothing supplies counts, so the prior is 1.0 everywhere and PP comes out equal to PL.
#[test]
fn a_site_nothing_covers_gets_a_flat_prior() {
    let text = golden();
    let input = records(&section(&text, "vcf", "input"));
    for row in measured(&text, "no-panel") {
        if let Some(prior) = &row.prior {
            assert!(prior.iter().all(|value| *value == 0), "{prior:?}");
        }
    }
    // The last site is absent from the panel, so even the panel run is flat there.
    let panel_rows = measured(&text, "panel");
    let absent = panel_rows.last().expect("the absent site");
    assert!(absent
        .prior
        .as_ref()
        .expect("a prior")
        .iter()
        .all(|v| *v == 0));

    // Three samples is below the ten the rule asks for.
    assert_eq!(input[0].genotypes.len(), 3);
    // And the PP equals the PL wherever the prior is flat.
    let no_panel = measured(&text, "no-panel");
    for (row, record) in no_panel.iter().zip(&input) {
        for (called, genotype) in row.genotypes.iter().zip(&record.genotypes) {
            if let (Some(posteriors), Some(likelihoods)) = (&called.2, &genotype.likelihoods) {
                assert_eq!(posteriors, likelihoods);
            }
        }
    }
}

/// The count is of samples and what enters the prior is chromosomes, and it only applies where the
/// panel is silent.
#[test]
fn the_reference_samples_are_doubled_and_only_used_where_the_panel_is_silent() {
    let text = golden();
    let with_ref = measured(&text, "ref-samples")[0]
        .prior
        .clone()
        .expect("a prior");
    // 2 * 1000 reference chromosomes, plus the input's own three of each.
    assert_eq!(
        gls_to_pls(&get_dirichlet_prior(
            &[2000.0 + 3.0 + 1e-3, 3.0 + 1e-3],
            2,
            false
        )),
        with_ref
    );

    // With the panel, the first five sites take the panel's prior and only the sixth takes the
    // padding.
    let both = measured(&text, "panel-and-ref");
    let panel_only = measured(&text, "panel");
    for index in 0..4 {
        assert_eq!(both[index].prior, panel_only[index].prior, "site {index}");
    }
    assert_eq!(both.last().expect("the absent site").prior, Some(with_ref));
}

/// An allele the same length as the reference takes the SNP pseudocount and every other allele the
/// indel one, so one site can mix them.
#[test]
fn the_pseudocount_is_chosen_by_allele_length() {
    let text = golden();
    let input = records(&section(&text, "vcf", "input"));
    let mixed = input
        .iter()
        .find(|record| record.id == "mixed")
        .expect("the mixed site");
    assert_eq!(mixed.reference().bases, "AT");
    assert_eq!(
        mixed.alleles[1].bases, "GT",
        "the same length as the reference"
    );
    assert_eq!(mixed.alleles[2].bases, "A", "shorter");

    // At the mixed site the panel's counts are in the hundreds, so the pseudocounts cannot be
    // seen there. The rule is asserted where nothing else is: with no counts at all, an alternate
    // of the reference's length takes the SNP pseudocount and a shorter one the indel one, and
    // making the two differ makes the two priors differ.
    let apart = Options {
        snp_prior_dirichlet: 1.0,
        indel_prior_dirichlet: 0.5,
        ..Options::default()
    };
    // No genotypes, so the prior is the whole of the answer, and one assumed reference sample so
    // that it is not the flat one.
    let same_length = Record {
        alleles: alleles("AT", "GT"),
        genotypes: vec![],
        ..mixed.clone()
    };
    let shorter = Record {
        alleles: alleles("AT", "A"),
        genotypes: vec![],
        ..mixed.clone()
    };
    let of = |record: &Record| {
        apply(record, &[], 1, &apart, false)
            .expect("population priors")
            .prior
    };
    assert_ne!(of(&same_length), of(&shorter));
    // And with the two pseudocounts equal they agree again.
    let together = Options {
        snp_prior_dirichlet: 1.0,
        indel_prior_dirichlet: 1.0,
        ..Options::default()
    };
    let of = |record: &Record| {
        apply(record, &[], 1, &together, false)
            .expect("population priors")
            .prior
    };
    assert_eq!(of(&same_length), of(&shorter));
}

/// A genotype can change, and its GQ is recomputed from the posteriors rather than kept.
#[test]
fn the_genotypes_are_recalled_from_the_posteriors() {
    let text = golden();
    let input = records(&section(&text, "vcf", "input"));
    let strong = measured(&text, "discovered-off");
    let mut changed = 0;
    for (row, record) in strong.iter().zip(&input) {
        for (called, genotype) in row.genotypes.iter().zip(&record.genotypes) {
            if called.0 != genotype.alleles {
                changed += 1;
            }
        }
    }
    assert!(changed > 0, "a strong reference prior recalls a genotype");

    // And a GQ moved with it.
    let flat = measured(&text, "no-panel");
    assert_ne!(flat[0].genotypes[0].1, strong[0].genotypes[0].1);
}

/// The Dirichlet-multinomial itself, on the two-allele shape the priors are built from.
#[test]
fn the_dirichlet_multinomial_is_the_reference_formula() {
    // Equal counts make the heterozygote the most likely genotype, which is what a small symmetric
    // prior looks like.
    let prior = get_dirichlet_prior(&[3.001, 3.001], 2, false);
    assert!(prior[1] > prior[0] && prior[1] > prior[2]);
    assert_eq!(prior[0], prior[2]);

    // A flat prior is literally 1.0, not a normalised probability.
    assert_eq!(
        get_dirichlet_prior(&[3.0, 3.0], 2, true),
        vec![1.0, 1.0, 1.0]
    );

    // And the three genotypes of a biallelic site sum to one in probability space.
    let total: f64 = prior.iter().map(|value| 10f64.powf(*value)).sum();
    assert!((total - 1.0).abs() < 1e-9, "{total}");

    // The formula, spelled out for one genotype.
    assert!((dirichlet_multinomial(&[3.001, 3.001], &[1, 1]) - prior[1]).abs() < 1e-12);
}
