//! Conformance for `CalculateGenotypePosteriors`' family priors against GATK 4.6.2.0, compared as
//! every trio member's posteriors and joint tags at every site of every run.
//!
//! Golden from `tools/readfilter-conformance/FamilyPriorsDump.java`.
//!
//! # What this suite is for
//!
//!  * **a Mendelian violation being overturned at the default de novo prior**, and standing at a
//!    higher one;
//!  * **the non-violation coefficient not being one**;
//!  * **an uncalled parent being a uniform third and an uncalled child stopping the trio**;
//!  * **the joint likelihood being read at the posterior's argmax**;
//!  * **and family priors never reaching a triallelic site**.

use gatk_corpus as corpus;
use gatk_tools::calculate_genotype_posteriors::{genotype_alleles, gls_to_pls};
use gatk_tools::family_priors::{
    combination_mv_count, configuration_likelihoods, likelihoods_as_map, updated_genotypes,
    GenotypeType, Member, DEFAULT_DE_NOVO_PRIOR, LOG10_OF_ONE_THIRD,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/family_priors.txt.gz"),
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

/// One sample's column of one row, by key.
fn field<'a>(keys: &[&str], values: &'a [&str], key: &str) -> Option<&'a str> {
    keys.iter()
        .position(|name| *name == key)
        .map(|index| values[index])
        .filter(|value| *value != ".")
}

/// The genotype type a `GT` implies, over the diploid calls this fixture uses.
fn genotype_type(gt: &str) -> GenotypeType {
    let alleles: Vec<&str> = gt.split(['/', '|']).collect();
    if alleles.contains(&".") {
        return GenotypeType::Uncalled;
    }
    match (alleles[0], alleles[1]) {
        ("0", "0") => GenotypeType::HomRef,
        (a, b) if a == b => GenotypeType::HomVar,
        _ => GenotypeType::Het,
    }
}

/// Every site of one run, as a map from sample name to its column.
struct Site {
    id: String,
    alternates: usize,
    columns: Vec<(String, Vec<String>)>,
    keys: Vec<String>,
}

fn sites(text: &str, label: &str) -> Vec<Site> {
    let body = section(text, "out", label);
    let mut samples: Vec<String> = Vec::new();
    let mut out = Vec::new();
    for line in body.lines() {
        if line.starts_with("#CHROM") {
            samples = line.split('\t').skip(9).map(str::to_string).collect();
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        out.push(Site {
            id: columns[2].to_string(),
            alternates: columns[4].split(',').count(),
            keys: columns[8].split(':').map(str::to_string).collect(),
            columns: samples
                .iter()
                .enumerate()
                .map(|(index, sample)| {
                    (
                        sample.clone(),
                        columns[9 + index].split(':').map(str::to_string).collect(),
                    )
                })
                .collect(),
        });
    }
    out
}

impl Site {
    fn get(&self, sample: &str, key: &str) -> Option<String> {
        let keys: Vec<&str> = self.keys.iter().map(String::as_str).collect();
        let column = self
            .columns
            .iter()
            .find(|(name, _)| name == sample)
            .map(|(_, values)| values)?;
        let values: Vec<&str> = column.iter().map(String::as_str).collect();
        // A column can be shorter than the format when the genotype is a bare no-call.
        if values.len() < keys.len() {
            return None;
        }
        field(&keys, &values, key).map(str::to_string)
    }

    fn member(&self, sample: &str) -> Member {
        let gt = self.get(sample, "GT").unwrap_or_else(|| {
            self.columns
                .iter()
                .find(|(n, _)| n == sample)
                .expect("a sample")
                .1[0]
                .clone()
        });
        Member {
            sample: sample.to_string(),
            genotype_type: genotype_type(&gt),
            likelihoods: self.get(sample, "PL").map(|value| {
                value
                    .split(',')
                    .map(|part| part.parse().expect("a likelihood"))
                    .collect()
            }),
            posteriors: None,
        }
    }
}

const TRIO: [&str; 3] = ["mom", "dad", "kid"];

/// label, de novo prior.
const RUNS: &[(&str, f64)] = &[
    ("family-only", DEFAULT_DE_NOVO_PRIOR),
    ("denovo-high", 0.001),
    ("denovo-low", 1e-9),
];

/// Every trio member's posteriors and joint tags, at every site of every family-prior run.
#[test]
fn every_trio_matches_the_golden() {
    let text = golden();
    // The likelihoods go in from the run that applied no family priors at all, so the port is
    // never handed its own output.
    let inputs = sites(&text, "skip-family");
    let mut compared = 0;
    for (label, de_novo) in RUNS {
        let measured = sites(&text, label);
        for (input, output) in inputs.iter().zip(&measured) {
            // Family priors never reach a site with more than one alternate.
            if input.alternates > 1 {
                for sample in TRIO {
                    assert_eq!(output.get(sample, "PP"), None, "{label}/{}", input.id);
                }
                continue;
            }
            let members: Vec<Member> = TRIO.iter().map(|sample| input.member(sample)).collect();
            let configuration = configuration_likelihoods(
                Some(&members[0]),
                Some(&members[1]),
                Some(&members[2]),
                *de_novo,
            );
            let Some(configuration) = configuration else {
                // The child was uncalled, so nothing was written.
                for sample in TRIO {
                    assert_eq!(output.get(sample, "PP"), None, "{label}/{}", input.id);
                }
                continue;
            };
            let updated = updated_genotypes(
                Some(&members[0]),
                Some(&members[1]),
                Some(&members[2]),
                &configuration,
            );
            for produced in &updated {
                let Some(posteriors) = &produced.posteriors else {
                    assert_eq!(
                        output.get(&produced.sample, "PP"),
                        None,
                        "{label}/{}/{}",
                        input.id,
                        produced.sample
                    );
                    continue;
                };
                let expected: Vec<i32> = output
                    .get(&produced.sample, "PP")
                    .unwrap_or_else(|| panic!("{label}/{}/{} has a PP", input.id, produced.sample))
                    .split(',')
                    .map(|part| part.parse().expect("a posterior"))
                    .collect();
                assert_eq!(
                    *posteriors, expected,
                    "{label}/{}/{}",
                    input.id, produced.sample
                );

                for (tag, value) in [
                    ("JL", produced.joint_likelihood),
                    ("JP", produced.joint_posterior),
                ] {
                    let measured_tag: Option<i32> = output
                        .get(&produced.sample, tag)
                        .map(|value| value.parse().expect("a joint tag"));
                    assert_eq!(
                        value, measured_tag,
                        "{tag} {label}/{}/{}",
                        input.id, produced.sample
                    );
                }

                // And the genotype the posteriors recall.
                let best = posteriors
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, value)| **value)
                    .map(|(index, _)| index)
                    .expect("a best genotype");
                let recalled = genotype_alleles(2, best);
                let measured_gt = output
                    .get(&produced.sample, "GT")
                    .expect("a GT")
                    .split(['/', '|'])
                    .map(|value| value.parse::<usize>().expect("an allele"))
                    .collect::<Vec<usize>>();
                assert_eq!(
                    recalled, measured_gt,
                    "GT {label}/{}/{}",
                    input.id, produced.sample
                );
                compared += 1;
            }
        }
    }
    assert_eq!(
        compared,
        3 * 3 * 3 + 2 * 3,
        "the trio members that were updated"
    );
}

/// At the default prior the violation is overturned; at a thousand times that it stands.
#[test]
fn a_violation_is_overturned_at_the_default_prior() {
    let text = golden();
    let default = sites(&text, "family-only");
    let high = sites(&text, "denovo-high");
    let violation = default
        .iter()
        .position(|site| site.id == "violation")
        .expect("the site");

    // The input had two hom-ref parents and a hom-var child.
    let input = &sites(&text, "skip-family")[violation];
    assert_eq!(input.get("mom", "GT").as_deref(), Some("0/0"));
    assert_eq!(input.get("kid", "GT").as_deref(), Some("1/1"));

    // At 1e-6 all three come out het.
    for sample in TRIO {
        assert_eq!(
            default[violation].get(sample, "GT").as_deref(),
            Some("0/1"),
            "{sample}"
        );
    }
    // At 1e-3 the violation stands.
    assert_eq!(high[violation].get("mom", "GT").as_deref(), Some("0/0"));
    assert_eq!(high[violation].get("kid", "GT").as_deref(), Some("1/1"));
}

/// The consistent combinations are scaled down by ten times the de novo prior, which is a term
/// nothing documents.
#[test]
fn the_non_violation_coefficient_is_not_one() {
    let de_novo = 0.25;
    let consistent = 1.0 - 10.0 * de_novo - de_novo * de_novo;
    assert!(consistent < 1.0);
    // At a quarter it is actually NEGATIVE, which log10 turns into NaN.
    assert!(consistent < 0.0, "{consistent}");

    // At the default it is just under one.
    let default =
        1.0 - 10.0 * DEFAULT_DE_NOVO_PRIOR - DEFAULT_DE_NOVO_PRIOR * DEFAULT_DE_NOVO_PRIOR;
    // 1 - 10e-6 - 1e-12, which is a hair BELOW 0.99999 rather than above it.
    assert!(default < 0.99999 && default > 0.9999);
    assert!(default.log10() < 0.0);
}

/// A missing parent is a uniform third and the pair is still judged; a missing child stops the
/// trio before the matrix is filled.
#[test]
fn an_uncalled_parent_is_a_third_and_an_uncalled_child_stops_the_trio() {
    let text = golden();
    let inputs = sites(&text, "skip-family");
    let measured = sites(&text, "family-only");

    let no_father = inputs
        .iter()
        .position(|site| site.id == "no-father")
        .expect("the site");
    let members: Vec<Member> = TRIO.iter().map(|s| inputs[no_father].member(s)).collect();
    assert_eq!(
        members[1].genotype_type,
        GenotypeType::Uncalled,
        "the father"
    );
    assert_eq!(
        likelihoods_as_map(Some(&members[1])),
        [LOG10_OF_ONE_THIRD; 3],
        "an uncalled member is three equal thirds"
    );
    // The mother and the child were still recomputed.
    assert!(measured[no_father].get("mom", "PP").is_some());
    assert!(measured[no_father].get("kid", "PP").is_some());
    // And their joint tags are -1, because not all three were called.
    assert_eq!(measured[no_father].get("mom", "JL").as_deref(), Some("-1"));

    let no_child = inputs
        .iter()
        .position(|site| site.id == "no-child")
        .expect("the site");
    let members: Vec<Member> = TRIO.iter().map(|s| inputs[no_child].member(s)).collect();
    assert!(
        configuration_likelihoods(
            Some(&members[0]),
            Some(&members[1]),
            Some(&members[2]),
            1e-6
        )
        .is_none(),
        "an uncalled child stops the trio"
    );
    for sample in TRIO {
        assert_eq!(measured[no_child].get(sample, "PP"), None, "{sample}");
    }
}

/// The Mendelian violation table, which counts one parent when only one is called.
#[test]
fn the_violation_count_is_the_reference_table() {
    use GenotypeType::{Het, HomRef, HomVar, Uncalled};
    // Two hom-ref parents and a hom-var child: both parents violated.
    assert_eq!(combination_mv_count(HomRef, HomRef, HomVar), 2);
    // One of each: consistent.
    assert_eq!(combination_mv_count(Het, Het, HomVar), 0);
    assert_eq!(combination_mv_count(HomRef, HomVar, Het), 0);
    // A het child with two hom-ref parents violates once.
    assert_eq!(combination_mv_count(HomRef, HomRef, Het), 1);
    // With ONE parent, a het child is always consistent.
    assert_eq!(combination_mv_count(HomRef, Uncalled, Het), 0);
    assert_eq!(combination_mv_count(Uncalled, HomVar, Het), 0);
    // And a single hom-ref parent still violates a hom-var child, once.
    assert_eq!(combination_mv_count(HomRef, Uncalled, HomVar), 1);
    // No child, or no parents at all: nothing to violate.
    assert_eq!(combination_mv_count(HomRef, HomRef, Uncalled), 0);
    assert_eq!(combination_mv_count(Uncalled, Uncalled, HomVar), 0);
}

/// `log10(1/3)` is written with seven digits, and that rounded value reaches the output.
#[test]
fn the_uniform_third_is_a_rounded_constant() {
    let exact = (1.0f64 / 3.0).log10();
    assert_ne!(LOG10_OF_ONE_THIRD, exact);
    assert!((LOG10_OF_ONE_THIRD - exact).abs() < 1e-7);
    // Three of them normalise back to a flat prior, whichever value is used.
    assert_eq!(
        gls_to_pls(&[LOG10_OF_ONE_THIRD; 3]),
        vec![0, 0, 0],
        "a uniform member contributes nothing on its own"
    );
}
