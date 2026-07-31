//! Conformance for `F1R2`, `F2R1`, `FAD` and the fragment grouping under them, against the oracle.
//!
//! Golden from `tools/annotation-conformance/FragmentCountsDump.java`.
//!
//! ```text
//! frag    three-reads-one-name  s1  100-119:p/65,p/129;100-119:b/65
//! grouped three-reads-one-name  A=-4609434218613702656,...
//! f1r2    unpaired-forward      F1R2=[0, 0][[I];F2R1=[1, 1][[I]
//! ```
//!
//! The `three-reads-one-name` pair of rows is the sharpest: the `Fragment` holds **two** reads,
//! because the supplementary alignment was dropped, while its likelihood is the sum over **three**,
//! because the grouping summed the whole group before the fragment was built. The trimming applies
//! to the evidence object and not to the arithmetic.

use std::io::Read;

use gatk_annotation::fragment_counts::{fragment_allele_depths, orientation_bias_counts};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

const START: i64 = 105;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fragment_counts.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn reference() -> Allele {
    Allele::from_str("A", true).expect("an allele")
}

fn alternate() -> Allele {
    Allele::from_str("C", false).expect("an allele")
}

/// name, flags, mapping quality, base quality, best allele index.
fn composition(label: &str) -> Vec<(&'static str, u16, u8, u8, usize)> {
    match label {
        "paired-f1r2" => vec![
            ("p", 0x41, 60, 30, 0),
            ("p", 0x81, 60, 30, 0),
            ("q", 0x51, 60, 30, 1),
            ("q", 0x91, 60, 30, 1),
        ],
        "paired-f2r1" => vec![
            ("p", 0x51, 60, 30, 0),
            ("p", 0x91, 60, 30, 0),
            ("q", 0x41, 60, 30, 1),
            ("q", 0x81, 60, 30, 1),
        ],
        "unpaired-forward" => vec![("a", 0, 60, 30, 0), ("b", 0, 60, 30, 1)],
        "unpaired-reverse" => vec![("a", 0x10, 60, 30, 0), ("b", 0x10, 60, 30, 1)],
        "mixed" => vec![
            ("p", 0x41, 60, 30, 0),
            ("p", 0x81, 60, 30, 0),
            ("q", 0x51, 60, 30, 1),
            ("r", 0x41, 60, 30, 1),
            ("s", 0x91, 60, 30, 0),
        ],
        "low-mapping-quality" => vec![("a", 0x41, 0, 30, 0), ("b", 0x41, 60, 30, 1)],
        "unavailable-mapping-quality" => vec![("a", 0x41, 255, 30, 0), ("b", 0x41, 60, 30, 1)],
        "low-base-quality" => vec![("a", 0x41, 60, 5, 0), ("b", 0x41, 60, 30, 1)],
        "singleton-and-pair" => vec![
            ("p", 0x41, 60, 30, 0),
            ("p", 0x81, 60, 30, 0),
            ("a", 0x41, 60, 30, 1),
        ],
        "three-reads-one-name" => vec![
            ("p", 0x41, 60, 30, 0),
            ("p", 0x81, 60, 30, 0),
            ("p", 0x841, 60, 30, 0),
            ("b", 0x41, 60, 30, 1),
        ],
        "second-read-low-quality" => vec![
            ("p", 0x41, 60, 30, 0),
            ("p", 0x81, 60, 2, 0),
            ("b", 0x41, 60, 30, 1),
        ],
        "empty" => vec![],
        other => panic!("{other} has no composition"),
    }
}

fn read(name: &str, flags: u16, mapping_quality: u8, base_quality: u8) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        flags,
        reference_index: 0,
        alignment_start: 100,
        mapping_quality,
        cigar: htsjdk_bam::text_parse::parse_cigar("20M").expect("a cigar"),
        read_bases: vec![b'A'; 20],
        base_qualities: vec![base_quality; 20],
        ..Default::default()
    }
}

fn read_likelihoods(label: &str) -> AlleleLikelihoods<BamRecord> {
    let composition = composition(label);
    let reads: Vec<BamRecord> = composition
        .iter()
        .map(|(name, flags, mq, bq, _)| read(name, *flags, *mq, *bq))
        .collect();
    let values: Vec<Vec<f64>> = (0..2)
        .map(|a| {
            composition
                .iter()
                .map(|(_, _, _, _, best)| if *best == a { -1.0 } else { -10.0 })
                .collect()
        })
        .collect();
    AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        AlleleList::new(&[reference(), alternate()]),
        vec![reads],
        vec![values],
    )
    .expect("a matrix")
}

fn variant_context() -> VariantContext {
    let mut vc = VariantContext::new("chr1", START, vec![reference(), alternate()]);
    vc.stop = START;
    vc.genotypes
        .push(Genotype::new("s1", vec![reference(), alternate()]));
    vc
}

/// `Arrays.toString(int[])` as the dump renders it, with the `[I` class name Java gives an
/// `int[]`.
fn rendered_ints(key: &str, values: &[i32]) -> String {
    format!(
        "{key}=[{}][[I]",
        values
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    )
}

#[test]
fn every_fragment_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("frag\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.splitn(3, '\t').collect();
        let label = fields[0];
        let expected = fields.get(2).copied().unwrap_or("");
        let grouped = read_likelihoods(label)
            .group_by_fragment()
            .expect("a grouping");
        let ours = grouped
            .sample_evidence(0)
            .unwrap_or(&[])
            .iter()
            .map(|fragment| {
                let names = fragment
                    .reads
                    .iter()
                    .map(|read| format!("{}/{}", read.read_name, read.flags))
                    .collect::<Vec<_>>()
                    .join(",");
                format!("{}-{}:{}", fragment.start, fragment.end, names)
            })
            .collect::<Vec<_>>()
            .join(";");
        assert_eq!(ours, expected, "fragments of {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no fragment rows");
    println!("{count} fragment lists identical");
}

#[test]
fn every_grouped_likelihood_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("grouped\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("a label and the values");
        let grouped = read_likelihoods(label)
            .group_by_fragment()
            .expect("a grouping");
        let evidence_count = grouped.sample_evidence(0).map(|e| e.len()).unwrap_or(0);
        let ours = (0..grouped.number_of_alleles())
            .map(|allele| {
                let values = (0..evidence_count)
                    .map(|e| (grouped.value(0, allele, e).to_bits() as i64).to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                format!(
                    "{}={values}",
                    grouped
                        .get_allele(allele)
                        .map(|a| a.display_string())
                        .unwrap_or_default()
                )
            })
            .collect::<Vec<_>>()
            .join(";");
        assert_eq!(ours, expected, "grouped likelihoods of {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no grouped rows");
    println!("{count} grouped likelihood tables bit-identical");
}

#[test]
fn every_orientation_count_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("f1r2\t") else {
            continue;
        };
        let (label, expected) = match rest.split_once('\t') {
            Some(pair) => pair,
            None => (rest, ""),
        };
        let vc = variant_context();
        let grouped = read_likelihoods(label)
            .group_by_fragment()
            .expect("a grouping");
        let ours = match orientation_bias_counts(&vc, "s1", Some(&grouped)).expect("counts") {
            None => String::new(),
            Some((f1r2, f2r1)) => format!(
                "{};{}",
                rendered_ints("F1R2", &f1r2),
                rendered_ints("F2R1", &f2r1)
            ),
        };
        assert_eq!(ours, expected, "orientation counts of {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no orientation rows");
    println!("{count} orientation count pairs identical");
}

#[test]
fn every_fragment_depth_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("fad\t") else {
            continue;
        };
        let (label, expected) = match rest.split_once('\t') {
            Some(pair) => pair,
            None => (rest, ""),
        };
        let vc = variant_context();
        let grouped = read_likelihoods(label)
            .group_by_fragment()
            .expect("a grouping");
        let ours = match fragment_allele_depths(&vc, "s1", true, Some(&grouped)).expect("depths") {
            None => String::new(),
            Some(depths) => rendered_ints("FAD", &depths),
        };
        assert_eq!(ours, expected, "fragment depths of {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no depth rows");
    println!("{count} fragment depth vectors identical");
}

#[test]
fn the_trimming_applies_to_the_evidence_and_not_to_the_arithmetic() {
    // Three reads share a name; the supplementary one is dropped from the fragment but its
    // likelihood was already summed in.
    let grouped = read_likelihoods("three-reads-one-name")
        .group_by_fragment()
        .expect("a grouping");
    let fragment = &grouped.sample_evidence(0).expect("evidence")[0];
    assert_eq!(fragment.reads.len(), 2);
    assert_eq!(grouped.value(0, 0, 0), -3.0);
}
