//! Conformance for `GATKVariantContextUtils.leftAlignAndTrim` against GATK 4.6.2.0, compared as
//! the record that comes back from every call.
//!
//! Golden from `tools/readfilter-conformance/LeftAlignAndTrimDump.java`.
//!
//! # What this suite is for
//!
//!  * **the window widens and can stop short**, so the same deletion lands in four places;
//!  * **three things come back untouched**, a non-indel, a window of zero, and a shift of zero;
//!  * **trimming decides whether the record shrinks**;
//!  * **and the genotypes are remapped**, which is visible only because the alleles moved.

use gatk_corpus as corpus;
use gatk_engine::variant_context_utils::{left_align_and_trim, Allele, Variant};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/left_align_and_trim.txt.gz"),
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

/// The reference the dump aligned against, which is the whole contig.
fn reference(text: &str) -> Vec<u8> {
    text.lines()
        .find_map(|line| line.strip_prefix("reference\t"))
        .expect("the reference row")
        .as_bytes()
        .to_vec()
}

/// The arguments each label was called with, which the golden does not carry.
fn arguments(label: &str) -> (i32, bool) {
    match label {
        "deletion-narrow-window" => (2, true),
        "deletion-exact-window" => (7, true),
        "deletion-in-long-run-narrow" => (10, true),
        "deletion-in-long-run-twenty" => (20, true),
        "zero-window" => (0, true),
        "negative-window" => (-5, true),
        "no-trim" | "no-trim-untrimmed-alleles" => (1000, false),
        _ => (1000, true),
    }
}

/// `chr1:17-18` back into its three pieces.
fn place(text: &str) -> (String, i32, i32) {
    let (contig, span) = text.split_once(':').expect("a contig");
    let (start, stop) = span.split_once('-').expect("a span");
    (
        contig.to_string(),
        start.parse().expect("a start"),
        stop.parse().expect("a stop"),
    )
}

/// `AA*,A` back into alleles, the `*` marking the reference.
fn alleles(text: &str) -> Vec<Allele> {
    text.split(',')
        .map(|field| match field.strip_suffix('*') {
            Some(bases) => Allele::new(bases.as_bytes(), true),
            None => Allele::new(field.as_bytes(), false),
        })
        .collect()
}

/// `s1=AA/A;s2=A/A` back into allele indices, which is how the port carries a genotype.
fn genotypes(text: &str, alleles: &[Allele]) -> Vec<Vec<usize>> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(';')
        .map(|entry| {
            let (_, called) = entry.split_once('=').expect("a sample");
            called
                .split('/')
                .map(|bases| {
                    alleles
                        .iter()
                        .position(|allele| allele.bases == bases.as_bytes())
                        .unwrap_or_else(|| panic!("allele {bases} is not in the record"))
                })
                .collect()
        })
        .collect()
}

fn variant(fields: &[String]) -> Variant {
    let (contig, start, stop) = place(&fields[1]);
    let alleles = alleles(&fields[2]);
    let genotypes = genotypes(fields.get(3).map(String::as_str).unwrap_or(""), &alleles);
    Variant {
        contig,
        start,
        stop,
        alleles,
        genotypes,
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
                if allele.is_reference { "*" } else { "" }
            )
        })
        .collect();
    format!(
        "{}:{}-{}\t{}",
        variant.contig,
        variant.start,
        variant.stop,
        alleles.join(",")
    )
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
    let bases = reference(&text);
    let all = labels(&text);
    assert!(
        all.len() >= 19,
        "every call is in the golden: {}",
        all.len()
    );

    for label in &all {
        let input = variant(&row(&text, "in", label));
        let (max_leading_bases, trim) = arguments(label);
        let ours = left_align_and_trim(&input, &bases, max_leading_bases, trim)
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
    let bases = reference(&text);
    for label in labels(&text) {
        let same = row(&text, "same", &label)[1] == "true";
        let input = variant(&row(&text, "in", &label));
        let (max_leading_bases, trim) = arguments(&label);
        let ours = left_align_and_trim(&input, &bases, max_leading_bases, trim).expect("a record");
        assert_eq!(ours == input, same, "same/{label}");
    }
}

/// One deletion, four windows, four answers, and only the widest is left aligned.
#[test]
fn the_window_decides_where_the_deletion_lands() {
    let text = golden();
    let places: Vec<String> = [
        "deletion-narrow-window",
        "deletion-in-long-run-narrow",
        "deletion-in-long-run-twenty",
        "deletion-in-long-run",
    ]
    .iter()
    .map(|label| row(&text, "out", label)[1].clone())
    .collect();
    assert_eq!(
        places,
        vec!["chr1:15-16", "chr1:49-50", "chr1:39-40", "chr1:30-31"]
    );
}

/// The genotypes follow the alleles, which only shows because the alleles moved.
#[test]
fn the_genotypes_are_remapped_with_the_alleles() {
    let text = golden();
    let bases = reference(&text);
    let input = variant(&row(&text, "in", "with-genotypes"));
    assert_eq!(input.genotypes, vec![vec![0, 1], vec![1, 1]]);
    let ours = left_align_and_trim(&input, &bases, 1000, true).expect("aligned");
    assert_eq!(ours.genotypes, input.genotypes);

    // Which the golden spells out with the new bases: s1=GA/G;s2=G/G.
    let expected = row(&text, "out", "with-genotypes");
    let rendered: Vec<String> = ours
        .genotypes
        .iter()
        .zip(["s1", "s2"])
        .map(|(genotype, sample)| {
            let called: Vec<String> = genotype
                .iter()
                .map(|index| String::from_utf8_lossy(&ours.alleles[*index].bases).to_string())
                .collect();
            format!("{sample}={}", called.join("/"))
        })
        .collect();
    assert_eq!(rendered.join(";"), expected[3]);
}
