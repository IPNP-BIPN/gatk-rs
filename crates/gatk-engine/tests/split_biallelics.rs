//! Conformance for `GATKVariantContextUtils.splitVariantContextToBiallelics` against GATK 4.6.2.0,
//! compared as the list of records that comes back from every call.
//!
//! Golden from `tools/readfilter-conformance/SplitBiallelicsDump.java`.
//!
//! # What this suite is for
//!
//!  * **a non-variant record becomes an empty list** and a biallelic one comes back as it stands;
//!  * **one het-non-ref call empties every record**, calls and attributes both;
//!  * **AC, AF and AN survive the filter and are then recomputed**, so a record with no genotypes
//!    loses them anyway;
//!  * **and every output is right trimmed on its own**, so two alternates land in two places.

use gatk_corpus as corpus;
use gatk_engine::subset_alleles::Genotype;
use gatk_engine::variant_context_utils::{split_variant_context_to_biallelics, Allele, Variant};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/split_biallelics.txt.gz"),
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

/// Whether each label asked for left trimming, which the golden does not carry.
fn trim_left(label: &str) -> bool {
    label == "left-trim"
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

/// `AC=[1, 2];AN=4;DP=30` back into pairs, in the order the dump sorted them.
fn attributes(text: &str) -> Vec<(String, String)> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(';')
        .map(|entry| {
            let (key, value) = entry.split_once('=').expect("an attribute");
            (key.to_string(), value.to_string())
        })
        .collect()
}

/// `s0=A/C;s1=C/G` back into allele indices, an empty call being a no-call.
fn genotypes(text: &str, alleles: &[Allele]) -> Vec<Genotype> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(';')
        .map(|entry| {
            let (_, called) = entry.split_once('=').expect("a sample");
            Genotype {
                alleles: called
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
                pl: None,
                gq: None,
                ad: None,
                dp: None,
                attributes: Vec::new(),
            }
        })
        .collect()
}

/// An `in` row as a record.
fn variant(fields: &[String]) -> Variant {
    let (start, stop) = place(&fields[1]);
    let alleles = alleles(&fields[2]);
    let genotypes = genotypes(fields.get(4).map(String::as_str).unwrap_or(""), &alleles);
    Variant {
        contig: "chr1".to_string(),
        start,
        stop,
        attributes: attributes(fields.get(3).map(String::as_str).unwrap_or("")),
        alleles,
        genotypes,
    }
}

/// How the dump prints one output record: place, alleles, attributes, genotypes.
fn rendered(variant: &Variant, samples: &[String]) -> String {
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
    let genotypes: Vec<String> = variant
        .genotypes
        .iter()
        .zip(samples)
        .map(|(genotype, sample)| {
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
            format!("{sample}={}", called.join("/"))
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

/// The sample names, which the golden carries in every genotype field.
fn samples(text: &str, label: &str) -> Vec<String> {
    let fields = row(text, "in", label);
    let genotypes = fields.get(4).map(String::as_str).unwrap_or("");
    if genotypes.is_empty() {
        return Vec::new();
    }
    genotypes
        .split(';')
        .map(|entry| entry.split_once('=').expect("a sample").0.to_string())
        .collect()
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
        all.len() >= 12,
        "every call is in the golden: {}",
        all.len()
    );

    for label in &all {
        let input = variant(&row(&text, "in", label));
        let ours = split_variant_context_to_biallelics(&input, trim_left(label))
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));

        let count: usize = row(&text, "count", label)[1].parse().expect("a count");
        assert_eq!(ours.len(), count, "count/{label}");

        let names = samples(&text, label);
        let expected: Vec<Vec<String>> = rows(&text, "out")
            .into_iter()
            .filter(|row| row[0] == label)
            .map(|row| row.into_iter().map(|field| field.to_string()).collect())
            .collect();
        for (index, record) in ours.iter().enumerate() {
            let fields = &expected[index];
            assert_eq!(
                rendered(record, &names),
                format!("{}\t{}\t{}\t{}", fields[2], fields[3], fields[4], fields[5]),
                "out/{label}/{index}"
            );
        }
    }
}

/// The two shapes that never split, one of which disappears.
#[test]
fn a_record_with_no_alternate_disappears() {
    let text = golden();
    let none = variant(&row(&text, "in", "no-alternate"));
    assert!(split_variant_context_to_biallelics(&none, false)
        .expect("an empty list")
        .is_empty());

    // And a biallelic record comes back as it stands, untrimmed.
    let biallelic = variant(&row(&text, "in", "biallelic"));
    let ours = split_variant_context_to_biallelics(&biallelic, false).expect("one record");
    assert_eq!(ours, vec![biallelic.clone()]);
    assert_eq!(row(&text, "same", "biallelic")[1], "true");
}

/// One 1/2 call, and every sample in every record loses its call and its attributes.
#[test]
fn one_het_non_ref_call_empties_everything() {
    let text = golden();
    let input = variant(&row(&text, "in", "het-non-ref-with-attributes"));
    assert!(!input.attributes.is_empty(), "the input carries attributes");
    let ours = split_variant_context_to_biallelics(&input, false).expect("two records");
    for record in &ours {
        assert!(record.attributes.is_empty());
        assert!(record
            .genotypes
            .iter()
            .all(|genotype| genotype.alleles.iter().all(Option::is_none)));
    }
}

/// The three counts survive the filter and are recomputed, or vanish for want of genotypes.
#[test]
fn the_counts_are_recomputed_and_a_record_without_genotypes_loses_them() {
    let text = golden();
    let counted = variant(&row(&text, "in", "no-het-non-ref"));
    let ours = split_variant_context_to_biallelics(&counted, false).expect("two records");
    assert_eq!(
        ours[0].attributes,
        vec![
            ("AC".to_string(), "1".to_string()),
            ("AF".to_string(), "0.16666666666666666".to_string()),
            ("AN".to_string(), "6".to_string()),
        ]
    );

    let no_genotypes = variant(&row(&text, "in", "attributes"));
    assert!(!no_genotypes.attributes.is_empty());
    let ours = split_variant_context_to_biallelics(&no_genotypes, false).expect("two records");
    assert!(ours.iter().all(|record| record.attributes.is_empty()));
}
