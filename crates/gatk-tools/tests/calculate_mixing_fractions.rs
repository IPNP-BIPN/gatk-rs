//! Conformance for `CalculateMixingFractions` against GATK 4.6.2.0, compared as the whole output
//! table of every run and as the class and message of the refusal.
//!
//! Golden from `tools/readfilter-conformance/CalculateMixingFractionsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the row order is a `HashMap` bucket order**, so a header of `zebra, alpha, mike` writes
//!    `zebra, mike, alpha`;
//!  * **one uncounted sample makes every row NaN**, the normalizer being a sum of fractions;
//!  * **singleton is either of two tests**, and `AC=2` with one het is counted;
//!  * **the sample is the first het**, and a singleton with no het is dropped;
//!  * **a read passing the site inside a deletion counts towards the total alone**;
//!  * **and a vendor-failed read is skipped by the tool**.

use gatk_corpus as corpus;
use gatk_tools::calculate_mixing_fractions::{
    is_biallelic_singleton_het_snp, mixing_fractions, site_counts, table, variant_sample,
    AltAndTotalReadCounts, CalculateMixingFractionsError,
};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};
use std::collections::HashMap;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/calculate_mixing_fractions.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// The header's sample list, which is the map's insertion order.
const SAMPLES: [&str; 3] = ["zebra", "alpha", "mike"];

/// The BAM the dump built, decoded from its base64 fixture row.
fn reads(text: &str) -> Vec<BamRecord> {
    let encoded = rows(text, "fixture")
        .into_iter()
        .find(|row| row[0] == "reads")
        .expect("the pooled bam")[1]
        .to_string();
    let bytes = corpus::decode_base64(&encoded);
    let decompressed = htsjdk_bgzf::read::decompress_all(&bytes).expect("the fixture is BGZF");
    let reader = htsjdk_bam::reader::BamReader::new(&decompressed).expect("the fixture opens");
    reader.map(|record| record.expect("a record")).collect()
}

/// The records of one input, decoded as far as the two tests look at them.
fn variants(text: &str, label: &str) -> Vec<VariantContext> {
    let whole = unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    );
    whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let mut alleles = vec![Allele::create(field[3].as_bytes(), true).expect("a reference")];
            for alternate in field[4].split(',') {
                alleles.push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
            }
            let mut variant =
                VariantContext::new(field[0], field[1].parse().expect("a position"), alleles);
            if field[7] != "." {
                for entry in field[7].split(';') {
                    if let Some((key, value)) = entry.split_once('=') {
                        variant
                            .attributes
                            .push((key.to_string(), Value::Str(value.to_string())));
                    }
                }
            }
            let called = variant.alleles.clone();
            variant.genotypes = SAMPLES
                .iter()
                .enumerate()
                .map(|(index, sample)| {
                    let call = field[9 + index];
                    Genotype::new(
                        sample,
                        call.split(['/', '|'])
                            .map(|allele| match allele.parse::<usize>() {
                                Ok(at) => called[at].clone(),
                                Err(_) => Allele::no_call(),
                            })
                            .collect(),
                    )
                })
                .collect();
            variant
        })
        .collect()
}

/// One run: which file it reads and which sites the intervals leave in the traversal.
fn setup(label: &str) -> (&'static str, Option<(i64, i64)>, bool) {
    match label {
        "every-shape" => ("every-shape", None, true),
        "one-site" => ("one-site", None, true),
        "every-sample-counted" => ("every-sample-counted", None, true),
        // No -I at all, so there is no pileup to count.
        "no-reads" => ("every-sample-counted", None, false),
        // `-L chr1:75-85` leaves the site at 80 alone.
        "one-interval" => ("every-sample-counted", Some((75, 85)), true),
        other => panic!("no setup for {other}"),
    }
}

/// The traversal and both branches of `apply`, as the tool runs them.
fn counts(text: &str, label: &str) -> HashMap<String, AltAndTotalReadCounts> {
    let (file, interval, has_reads) = setup(label);
    let pooled = reads(text);
    let mut buckets: HashMap<String, AltAndTotalReadCounts> = HashMap::new();

    for variant in variants(text, file) {
        if let Some((start, end)) = interval {
            if variant.stop < start || variant.start > end {
                continue;
            }
        }
        if !is_biallelic_singleton_het_snp(&variant) {
            continue;
        }
        let Some(sample) = variant_sample(&variant) else {
            continue;
        };
        let alt = variant.alleles[1].base_string().as_bytes()[0];
        let overlapping: Vec<BamRecord> = if has_reads {
            pooled
                .iter()
                .filter(|read| {
                    let start = read.alignment_start;
                    let end = start + read.cigar.reference_length() as i32 - 1;
                    start <= variant.start as i32 && variant.start as i32 <= end
                })
                .cloned()
                .collect()
        } else {
            Vec::new()
        };
        let site = site_counts(&overlapping, variant.start as i32, alt);
        let bucket = buckets.entry(sample).or_default();
        bucket.alt += site.alt;
        bucket.total += site.total;
    }
    buckets
}

#[test]
fn every_table_matches_the_golden_byte_for_byte() {
    let text = golden();
    let tables = rows(&text, "table");
    assert_eq!(tables.len(), 5, "five of the six runs write a table");

    let samples: Vec<String> = SAMPLES.iter().map(|name| name.to_string()).collect();
    for row in tables {
        let (label, expected) = (row[0], unescape(row[1]));
        let fractions =
            mixing_fractions(&samples, &counts(&text, label)).expect("no treeified bucket");
        assert_eq!(table(&fractions), expected, "the table of {label}");
    }
}

#[test]
fn the_rows_are_in_hash_order_and_not_the_headers() {
    let text = golden();
    let written = unescape(
        rows(&text, "table")
            .into_iter()
            .find(|row| row[0] == "every-sample-counted")
            .expect("the run where every sample is counted")[1],
    );
    let order: Vec<&str> = written
        .lines()
        .skip(1)
        .map(|line| line.split('\t').next().expect("a sample"))
        .collect();
    assert_eq!(order, vec!["zebra", "mike", "alpha"]);
    assert_ne!(order, SAMPLES.to_vec(), "not the header's order");
}

#[test]
fn one_uncounted_sample_makes_every_row_nan() {
    let text = golden();
    for label in ["one-site", "no-reads", "one-interval"] {
        let written = unescape(
            rows(&text, "table")
                .into_iter()
                .find(|row| row[0] == label)
                .unwrap_or_else(|| panic!("no run {label}"))[1],
        );
        assert!(
            written.lines().skip(1).all(|line| line.ends_with("NaN")),
            "{label} is NaN throughout"
        );
    }
}

#[test]
fn a_read_passing_the_site_inside_a_deletion_counts_towards_the_total_alone() {
    let text = golden();
    let bucket = counts(&text, "every-shape");
    // `mike` is counted only at the site the deletion spans.
    assert_eq!(
        bucket.get("mike").copied().expect("mike is counted"),
        AltAndTotalReadCounts { alt: 0, total: 1 }
    );
}

#[test]
fn the_vendor_failed_read_is_skipped_by_the_tool() {
    let text = golden();
    let bucket = counts(&text, "every-shape");
    // At the second site three reads overlap and one of them fails the vendor check.
    assert_eq!(
        bucket.get("alpha").copied().expect("alpha is counted"),
        AltAndTotalReadCounts { alt: 1, total: 2 }
    );
}

#[test]
fn the_refusal_carries_the_references_class_and_words() {
    let text = golden();
    let row = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "output-is-a-directory")
        .expect("the refusal");
    let (class, message) = row[1].split_once(':').expect("class and message");

    let refused = CalculateMixingFractionsError::CouldNotCreateOutputFile {
        path: "calculatemixingfractions-dump/.".to_string(),
    };
    assert_eq!(refused.class(), class);
    assert_eq!(refused.message(), message);
}
