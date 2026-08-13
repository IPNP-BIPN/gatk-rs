//! Conformance for `AlleleSubsettingUtils.subsetAlleles` against GATK 4.6.2.0, compared as the
//! whole genotype that comes back from every call.
//!
//! Golden from `tools/readfilter-conformance/SubsetAllelesDump.java`.
//!
//! # What this suite is for
//!
//!  * **the PLs are permuted and rescaled**, and **the GQ is recomputed even when nothing moved**;
//!  * **a PL array of the wrong length is dropped**, GQ with it;
//!  * **the no-data cases come before the assignment method**;
//!  * **and the two no-call methods differ**, one keeping everything and one keeping the depth.

use gatk_corpus as corpus;
use gatk_engine::subset_alleles::{subset_alleles, AssignmentMethod, Genotype};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/subset_alleles.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn row(text: &str, kind: &str, label: &str) -> Vec<String> {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no {kind} row for {label}"))
        .into_iter()
        .map(|field| field.to_string())
        .collect()
}

fn method_of(name: &str) -> AssignmentMethod {
    match name {
        "BEST_MATCH_TO_ORIGINAL" => AssignmentMethod::BestMatchToOriginal,
        "USE_PLS_TO_ASSIGN" => AssignmentMethod::UsePlsToAssign,
        "SET_TO_NO_CALL" => AssignmentMethod::SetToNoCall,
        "SET_TO_NO_CALL_NO_ANNOTATIONS" => AssignmentMethod::SetToNoCallNoAnnotations,
        other => panic!("no method {other}"),
    }
}

/// The three alleles the dump used, so a base can be turned back into an index.
const ALLELES: [&str; 3] = ["A", "C", "G"];

fn numbers(text: &str) -> Option<Vec<i32>> {
    if text.is_empty() {
        return None;
    }
    Some(
        text.split(',')
            .map(|value| value.parse().expect("a number"))
            .collect(),
    )
}

/// `alleles|PL|GQ|AD|DP|attributes` back into a genotype.
fn parse_genotype(text: &str) -> Genotype {
    let fields: Vec<&str> = text.split('|').collect();
    let alleles = fields[0]
        .split('/')
        .map(|base| {
            if base.is_empty() || base == "." {
                None
            } else {
                Some(
                    ALLELES
                        .iter()
                        .position(|allele| *allele == base)
                        .unwrap_or_else(|| panic!("allele {base}")),
                )
            }
        })
        .collect();
    Genotype {
        alleles,
        pl: numbers(fields[1]),
        gq: fields[2].parse().ok(),
        ad: numbers(fields[3]),
        dp: fields[4].parse().ok(),
        attributes: if fields[5].is_empty() {
            Vec::new()
        } else {
            fields[5]
                .split(';')
                .map(|entry| {
                    let (key, value) = entry.split_once('=').expect("an attribute");
                    (key.to_string(), value.to_string())
                })
                .collect()
        },
    }
}

/// How the dump prints a genotype, so the two can be compared as text.
///
/// The calls are indices into the KEPT list, not the original one, which is what the subsetting
/// translated them to: with `kept = [0, 2, 1]` an index of 2 is the allele that was 1.
fn rendered(genotype: &Genotype, kept: &[usize]) -> String {
    let alleles: Vec<String> = genotype
        .alleles
        .iter()
        .map(|allele| match allele {
            None => String::new(),
            Some(index) => ALLELES[kept[*index]].to_string(),
        })
        .collect();
    let list = |values: &Option<Vec<i32>>| match values {
        None => String::new(),
        Some(values) => values
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(","),
    };
    let attributes: Vec<String> = genotype
        .attributes
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect();
    format!(
        "{}|{}|{}|{}|{}|{}",
        alleles.join("/"),
        list(&genotype.pl),
        genotype.gq.map(|gq| gq.to_string()).unwrap_or_default(),
        list(&genotype.ad),
        genotype.dp.map(|dp| dp.to_string()).unwrap_or_default(),
        attributes.join(";")
    )
}

fn labels(text: &str) -> Vec<String> {
    rows(text, "in")
        .into_iter()
        .map(|row| row[0].to_string())
        .collect()
}

#[test]
fn every_subsetting_gives_the_reference_s_genotype() {
    let text = golden();
    let all = labels(&text);
    assert!(
        all.len() >= 15,
        "every call is in the golden: {}",
        all.len()
    );

    for label in &all {
        let input = parse_genotype(&row(&text, "in", label)[1]);
        let out = row(&text, "out", label);
        let method = method_of(&out[1]);
        let kept: Vec<usize> = out[2]
            .split(',')
            .map(|value| value.parse().expect("an allele index"))
            .collect();

        let ours = subset_alleles(&[input], 2, 3, &kept, method).expect("a genotype");
        assert_eq!(rendered(&ours[0], &kept), out[3], "out/{label}");
    }
}

/// Keeping every allele leaves the PLs alone and still changes the GQ.
#[test]
fn the_gq_is_recomputed_even_when_nothing_moved() {
    let text = golden();
    let input = parse_genotype(&row(&text, "in", "het-keep-all")[1]);
    assert_eq!(input.gq, Some(50));
    let ours = subset_alleles(
        std::slice::from_ref(&input),
        2,
        3,
        &[0, 1, 2],
        AssignmentMethod::BestMatchToOriginal,
    )
    .expect("a genotype");
    assert_eq!(ours[0].pl, input.pl, "the PLs did not move");
    assert_eq!(ours[0].gq, Some(30), "and the GQ did");
}

/// The genotype with no data at all, which is cleared rather than subset.
#[test]
fn a_genotype_with_no_data_keeps_nothing() {
    let text = golden();
    let input = parse_genotype(&row(&text, "in", "gq-and-dp-zero")[1]);
    let ours = subset_alleles(
        &[input],
        2,
        3,
        &[0, 1],
        AssignmentMethod::BestMatchToOriginal,
    )
    .expect("a genotype");
    assert_eq!(rendered(&ours[0], &[0, 1]), "/|||||");

    // The same genotype with a depth keeps its PLs and its AD, and is still no-called.
    let with_depth = parse_genotype(&row(&text, "in", "gq-zero")[1]);
    let ours = subset_alleles(
        &[with_depth],
        2,
        3,
        &[0, 1],
        AssignmentMethod::BestMatchToOriginal,
    )
    .expect("a genotype");
    assert_eq!(rendered(&ours[0], &[0, 1]), row(&text, "out", "gq-zero")[3]);
}

/// Two names one word apart, two different genotypes out.
#[test]
fn the_two_no_call_methods_differ() {
    let text = golden();
    let plain = row(&text, "out", "het-no-call")[3].clone();
    let no_annotations = row(&text, "out", "het-no-annotations")[3].clone();
    assert_ne!(plain, no_annotations);
    // One keeps the likelihoods, the depth and the AD; the other keeps only the depth.
    assert!(plain.starts_with("/|50,0,60|50|10,12|30|"));
    assert_eq!(no_annotations, "/||||30|");
}
