//! Conformance for splitting a record whose genotypes carry likelihoods, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/SplitWithLikelihoodsDump.java`.
//!
//! # What this suite is for
//!
//!  * **each output record carries its own subset of the PLs**, rescaled on its own;
//!  * **the call does not follow the likelihoods** under `BEST_MATCH_TO_ORIGINAL`;
//!  * **AC and AF are recomputed from the subset calls**;
//!  * **one het-non-ref empties everything**, keeping only the depth;
//!  * **and the trimming still happens afterwards**, so a record can move.

use gatk_corpus as corpus;
use gatk_engine::subset_alleles::Genotype;
use gatk_engine::variant_context_utils::{split_variant_context_to_biallelics, Allele, Variant};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/split_with_likelihoods.txt.gz"),
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

fn place(text: &str) -> (i32, i32) {
    let (start, stop) = text.split_once('-').expect("a span");
    (
        start.parse().expect("a start"),
        stop.parse().expect("a stop"),
    )
}

fn alleles(text: &str) -> Vec<Allele> {
    text.split(',')
        .map(|field| match field.strip_suffix("(ref)") {
            Some(bases) => Allele::new(bases.as_bytes(), true),
            None => Allele::new(field.as_bytes(), false),
        })
        .collect()
}

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

/// `s0=A/C|50,0,60|50|10,12|30;s1=...` back into genotypes, calls as indices into `alleles`.
fn genotypes(text: &str, alleles: &[Allele]) -> Vec<Genotype> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(';')
        .map(|entry| {
            let (_, rest) = entry.split_once('=').expect("a sample");
            let fields: Vec<&str> = rest.split('|').collect();
            Genotype {
                alleles: fields[0]
                    .split('/')
                    .map(|bases| {
                        if bases.is_empty() {
                            return None;
                        }
                        Some(
                            alleles
                                .iter()
                                .position(|allele| allele.bases == bases.as_bytes())
                                .unwrap_or_else(|| panic!("allele {bases} is not in the record")),
                        )
                    })
                    .collect(),
                pl: numbers(fields[1]),
                gq: fields[2].parse().ok(),
                ad: numbers(fields[3]),
                dp: fields[4].parse().ok(),
                attributes: Vec::new(),
            }
        })
        .collect()
}

/// The sample names, which the golden carries in every genotype field.
fn samples(text: &str) -> Vec<String> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(';')
        .map(|entry| entry.split_once('=').expect("a sample").0.to_string())
        .collect()
}

fn variant(fields: &[String]) -> Variant {
    let (start, stop) = place(&fields[1]);
    let alleles = alleles(&fields[2]);
    let genotypes = genotypes(fields.get(3).map(String::as_str).unwrap_or(""), &alleles);
    Variant {
        contig: "chr1".to_string(),
        start,
        stop,
        alleles,
        genotypes,
        attributes: Vec::new(),
    }
}

/// How the dump prints an output record.
fn rendered(variant: &Variant, names: &[String]) -> String {
    let alleles: Vec<String> = variant
        .alleles
        .iter()
        .map(|allele| {
            format!(
                "{}{}",
                String::from_utf8_lossy(&allele.bases),
                if allele.is_reference { "(ref)" } else { "" }
            )
        })
        .collect();
    let attributes: Vec<String> = variant
        .attributes
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect();
    let list = |values: &Option<Vec<i32>>| match values {
        None => String::new(),
        Some(values) => values
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(","),
    };
    let genotypes: Vec<String> = variant
        .genotypes
        .iter()
        .zip(names)
        .map(|(genotype, name)| {
            let called: Vec<String> = genotype
                .alleles
                .iter()
                .map(|allele| match allele {
                    None => String::new(),
                    Some(index) => {
                        String::from_utf8_lossy(&variant.alleles[*index].bases).to_string()
                    }
                })
                .collect();
            format!(
                "{name}={}|{}|{}|{}|{}",
                called.join("/"),
                list(&genotype.pl),
                genotype.gq.map(|gq| gq.to_string()).unwrap_or_default(),
                list(&genotype.ad),
                genotype.dp.map(|dp| dp.to_string()).unwrap_or_default()
            )
        })
        .collect();
    format!(
        "{}-{}\t{}\t{}\t{}",
        variant.start,
        variant.stop,
        alleles.join(","),
        attributes.join(";"),
        genotypes.join(";")
    )
}

fn labels(text: &str) -> Vec<String> {
    rows(text, "in")
        .into_iter()
        .map(|row| row[0].to_string())
        .collect()
}

#[test]
fn every_split_is_the_reference_s() {
    let text = golden();
    let all = labels(&text);
    assert!(
        all.len() >= 6,
        "every record is in the golden: {}",
        all.len()
    );

    for label in &all {
        let fields = row(&text, "in", label);
        let input = variant(&fields);
        let names = samples(fields.get(3).map(String::as_str).unwrap_or(""));
        let ours = split_variant_context_to_biallelics(&input, false)
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));

        let expected: Vec<Vec<String>> = rows(&text, "out")
            .into_iter()
            .filter(|row| row[0] == label)
            .map(|row| row.into_iter().map(|field| field.to_string()).collect())
            .collect();
        assert_eq!(ours.len(), expected.len(), "count/{label}");
        for (index, record) in ours.iter().enumerate() {
            let want = &expected[index];
            assert_eq!(
                rendered(record, &names),
                format!("{}\t{}\t{}\t{}", want[2], want[3], want[4], want[5]),
                "out/{label}/{index}"
            );
        }
    }
}

/// The call and the likelihoods disagree, and the port writes the disagreement.
#[test]
fn the_call_does_not_follow_the_likelihoods() {
    let text = golden();
    let fields = row(&text, "in", "two-samples");
    let input = variant(&fields);
    let ours = split_variant_context_to_biallelics(&input, false).expect("two records");

    // The record that kept the first alternate: s1 was A/G, and G is gone.
    let s1 = &ours[0].genotypes[1];
    assert_eq!(s1.alleles, vec![Some(0), Some(0)], "called hom-ref");
    let pl = s1.pl.as_ref().expect("likelihoods");
    // Yet its own subset PLs make the heterozygote the most likely genotype by 10 phred.
    assert_eq!(pl[1], 0, "the het is the best likelihood");
    assert!(pl[0] > 0 && pl[2] > 0);
}

/// One 1/2 call in one sample, and every sample of every record keeps only its depth.
#[test]
fn one_het_non_ref_empties_every_record() {
    let text = golden();
    let input = variant(&row(&text, "in", "het-non-ref"));
    let ours = split_variant_context_to_biallelics(&input, false).expect("two records");
    for record in &ours {
        assert!(record.attributes.is_empty());
        for genotype in &record.genotypes {
            assert!(genotype.alleles.iter().all(Option::is_none));
            assert!(genotype.pl.is_none() && genotype.gq.is_none() && genotype.ad.is_none());
            assert!(genotype.dp.is_some(), "the depth survives");
        }
    }
}

/// The counts follow the subset calls, so a hom-alt is everything in one record and nothing in the
/// other.
#[test]
fn the_counts_are_recomputed_per_record() {
    let text = golden();
    let input = variant(&row(&text, "in", "hom-alt"));
    let ours = split_variant_context_to_biallelics(&input, false).expect("two records");
    let of = |record: &Variant, key: &str| {
        record
            .attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.clone())
            .unwrap_or_default()
    };
    assert_eq!(
        (of(&ours[0], "AC"), of(&ours[0], "AF")),
        ("0".into(), "0.0".into())
    );
    assert_eq!(
        (of(&ours[1], "AC"), of(&ours[1], "AF")),
        ("2".into(), "1.0".into())
    );
}
