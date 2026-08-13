//! Conformance for `GATKVariantContextUtils.trimAlleles` against GATK 4.6.2.0, compared as the
//! record that comes back from every call.
//!
//! Golden from `tools/readfilter-conformance/TrimAllelesDump.java`.
//!
//! # What this suite is for
//!
//!  * **one allele of length one turns the whole thing off**, and the spanning deletion is the one
//!    exception;
//!  * **symbolic alleles and `*` are kept but not compared**;
//!  * **an allele trimmed to nothing gets one base back**, at whichever end is left;
//!  * **and each direction can be asked for on its own**, which the same record shows three ways.

use gatk_corpus as corpus;
use gatk_engine::subset_alleles::Genotype;
use gatk_engine::variant_context_utils::{trim_alleles, Allele, Variant};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/trim_alleles.txt.gz"),
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

/// Which directions each label was called with, which the golden does not carry.
fn arguments(label: &str) -> (bool, bool) {
    match label {
        "shared-suffix-forward-only" | "shared-prefix-forward-only" => (true, false),
        "shared-suffix-reverse-only" | "shared-prefix-reverse-only" => (false, true),
        "no-trim-requested" => (false, false),
        _ => (true, true),
    }
}

/// `100-103` back into its two numbers.
fn place(text: &str) -> (i32, i32) {
    let (start, stop) = text.split_once('-').expect("a span");
    (
        start.parse().expect("a start"),
        stop.parse().expect("a stop"),
    )
}

/// `ACGT(ref),AGT,*` back into alleles.
fn alleles(text: &str) -> Vec<Allele> {
    text.split(',')
        .map(|field| match field.strip_suffix("(ref)") {
            Some(bases) => Allele::new(bases.as_bytes(), true),
            None => Allele::new(field.as_bytes(), false),
        })
        .collect()
}

/// `s1=ACGT/AGT` back into allele indices.
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

/// How the dump prints a record, so the two can be compared as text.
fn rendered(variant: &Variant) -> String {
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
    format!("{}-{}\t{}", variant.start, variant.stop, alleles.join(","))
}

fn labels(text: &str) -> Vec<String> {
    rows(text, "in")
        .into_iter()
        .map(|row| row[0].to_string())
        .collect()
}

#[test]
fn every_call_gives_the_reference_s_record() {
    let text = golden();
    let all = labels(&text);
    assert!(
        all.len() >= 17,
        "every call is in the golden: {}",
        all.len()
    );

    for label in &all {
        let input = variant(&row(&text, "in", label));
        let (forward, reverse) = arguments(label);
        let ours = trim_alleles(&input, forward, reverse)
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));

        let expected = row(&text, "out", label);
        assert_eq!(
            rendered(&ours),
            format!("{}\t{}", expected[1], expected[2]),
            "out/{label}"
        );
    }
}

/// The record the reference handed straight back, which the golden marks.
#[test]
fn what_comes_back_untouched_is_untouched_here_too() {
    let text = golden();
    for label in labels(&text) {
        let same = row(&text, "same", &label)[1] == "true";
        let input = variant(&row(&text, "in", &label));
        let (forward, reverse) = arguments(&label);
        let ours = trim_alleles(&input, forward, reverse).expect("a record");
        assert_eq!(ours == input, same, "same/{label}");
    }
}

/// The test that turns the whole thing off, and the allele that is exempt from it.
#[test]
fn one_short_allele_stops_the_trim_and_the_star_does_not() {
    let text = golden();
    for label in ["one-base-allele", "snp-beside-indel"] {
        let input = variant(&row(&text, "in", label));
        assert_eq!(
            trim_alleles(&input, true, true).expect("untouched"),
            input,
            "{label}"
        );
    }

    // The same shape with `*` in place of the short allele is trimmed, and the `*` survives.
    let star = variant(&row(&text, "in", "spanning-deletion"));
    let trimmed = trim_alleles(&star, true, true).expect("trimmed");
    assert_ne!(trimmed, star);
    assert!(trimmed.alleles.iter().any(|allele| allele.is_span_del()));
    assert_eq!(rendered(&trimmed).split('\t').next().unwrap(), "100-101");
}

/// Each direction on its own, over one record, which is three different answers.
#[test]
fn each_direction_can_be_asked_for_alone() {
    let text = golden();
    let both = row(&text, "out", "shared-suffix")[1].clone();
    let forward = row(&text, "out", "shared-suffix-forward-only")[1].clone();
    let reverse = row(&text, "out", "shared-suffix-reverse-only")[1].clone();
    // A shared suffix is only the reverse trim's business.
    assert_eq!(both, "100-101");
    assert_eq!(forward, "100-103");
    assert_eq!(reverse, "100-101");

    let both = row(&text, "out", "shared-prefix")[1].clone();
    let forward = row(&text, "out", "shared-prefix-forward-only")[1].clone();
    let reverse = row(&text, "out", "shared-prefix-reverse-only")[1].clone();
    assert_eq!(both, "103-103");
    assert_eq!(forward, "103-103");
    assert_eq!(reverse, "100-103");
}
