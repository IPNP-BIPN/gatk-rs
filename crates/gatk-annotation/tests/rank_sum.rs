//! Conformance for the four rank-sum annotations, against the oracle.
//!
//! Golden from `tools/annotation-conformance/RankSumAnnotationDump.java`.
//!
//! ```text
//! anno  BaseQualityRankSumTest  three-and-three  BaseQRankSum=-1.960[java.lang.String]
//! anno  ReadPosRankSumTest      three-and-three  ReadPosRankSum=0.000[java.lang.String]
//! anno  BaseQualityRankSumTest  ref-only         (absent)
//! ```
//!
//! The value is a `java.lang.String` and not a Double, formatted to three decimals. A site whose
//! alternate series is empty has **no** key at all rather than one reading zero.

use std::io::Read;

use gatk_annotation::info_annotation::AnnotationValue;
use gatk_annotation::rank_sum::{
    annotate, BaseQualityRankSumTest, ClippingRankSumTest, MappingQualityRankSumTest,
    ReadPosRankSumTest,
};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

const START: i64 = 105;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/rank_sum.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn allele(bases: &str, is_ref: bool) -> Allele {
    Allele::from_str(bases, is_ref).expect("an allele")
}

fn site(with_genotypes: bool) -> VariantContext {
    let mut vc = VariantContext::new("chr1", START, vec![allele("A", true), allele("C", false)]);
    vc.stop = START;
    if with_genotypes {
        vc.genotypes.push(htsjdk_vcf::variant::Genotype::new(
            "s1",
            vec![allele("A", true), allele("C", false)],
        ));
    }
    vc
}

fn read(name: &str, mapping_quality: u8, base_quality: u8, start: i32, cigar: &str) -> BamRecord {
    let parsed = htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar");
    let length: usize = parsed
        .elements
        .iter()
        .filter(|e| {
            matches!(
                e.op,
                htsjdk_bam::cigar::Op::M | htsjdk_bam::cigar::Op::I | htsjdk_bam::cigar::Op::S
            )
        })
        .map(|e| e.length as usize)
        .sum();
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality,
        cigar: parsed,
        read_bases: vec![b'A'; length],
        base_qualities: vec![base_quality; length],
        ..Default::default()
    }
}

fn likelihoods_for(reads: Vec<BamRecord>, values: Vec<Vec<f64>>) -> AlleleLikelihoods<BamRecord> {
    AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        AlleleList::new(&[allele("A", true), allele("C", false)]),
        vec![reads],
        vec![values],
    )
    .expect("a matrix")
}

/// The dump's `twoSided`: `ref_count` reads on the reference, `alt_count` on the alternate, each
/// with a descending quality ramp.
fn two_sided(
    ref_count: usize,
    alt_count: usize,
    ref_base: i32,
    alt_base: i32,
    ref_mapping: i32,
    alt_mapping: i32,
) -> AlleleLikelihoods<BamRecord> {
    let mut reads = Vec::new();
    for i in 0..ref_count {
        reads.push(read(
            &format!("r{i}"),
            (ref_mapping - i as i32) as u8,
            (ref_base - i as i32) as u8,
            100,
            "10M",
        ));
    }
    for i in 0..alt_count {
        reads.push(read(
            &format!("a{i}"),
            (alt_mapping - i as i32) as u8,
            (alt_base - i as i32) as u8,
            100,
            "10M",
        ));
    }
    let values = (0..2)
        .map(|a| {
            (0..reads.len())
                .map(|e| {
                    let is_ref = e < ref_count;
                    if (a == 0) == is_ref {
                        -1.0
                    } else {
                        -10.0
                    }
                })
                .collect()
        })
        .collect();
    likelihoods_for(reads, values)
}

/// A matrix where the first `ref_count` reads support the reference and the rest the alternate.
fn split(reads: Vec<BamRecord>, ref_count: usize) -> AlleleLikelihoods<BamRecord> {
    let values = (0..2)
        .map(|a| {
            (0..reads.len())
                .map(|e| {
                    let is_ref = e < ref_count;
                    if (a == 0) == is_ref {
                        -1.0
                    } else {
                        -10.0
                    }
                })
                .collect()
        })
        .collect();
    likelihoods_for(reads, values)
}

fn case(label: &str) -> (VariantContext, Option<AlleleLikelihoods<BamRecord>>) {
    match label {
        "null-likelihoods" => (site(true), None),
        "no-genotypes" => (site(false), Some(two_sided(6, 6, 30, 20, 60, 40))),
        "empty-matrix" => (
            site(true),
            Some(likelihoods_for(Vec::new(), vec![Vec::new(), Vec::new()])),
        ),
        "three-and-three" => (site(true), Some(two_sided(3, 3, 30, 20, 60, 40))),
        "nine-and-nine" => (site(true), Some(two_sided(9, 9, 30, 20, 60, 40))),
        "ten-and-nine" => (site(true), Some(two_sided(10, 9, 30, 20, 60, 40))),
        "twelve-and-twelve" => (site(true), Some(two_sided(12, 12, 30, 20, 60, 40))),
        "identical-groups" => (site(true), Some(two_sided(12, 12, 30, 30, 60, 60))),
        "ref-only" => (site(true), Some(two_sided(12, 0, 30, 20, 60, 40))),
        "alt-only" => (site(true), Some(two_sided(0, 12, 30, 20, 60, 40))),
        "mapq-zero-and-unavailable" => (
            site(true),
            Some(split(
                vec![
                    read("r0", 0, 30, 100, "10M"),
                    read("r1", 255, 30, 100, "10M"),
                    read("r2", 60, 30, 100, "10M"),
                    read("a0", 0, 20, 100, "10M"),
                    read("a1", 255, 20, 100, "10M"),
                    read("a2", 40, 20, 100, "10M"),
                ],
                3,
            )),
        ),
        "reads-past-the-variant" => (
            site(true),
            Some(split(
                vec![
                    read("r0", 60, 30, 100, "10M"),
                    read("r1", 60, 30, 200, "10M"),
                    read("a0", 40, 20, 100, "10M"),
                    read("a1", 40, 20, 200, "10M"),
                ],
                2,
            )),
        ),
        "hard-clipped" => (
            site(true),
            Some(split(
                vec![
                    read("r0", 60, 30, 100, "10M"),
                    read("r1", 60, 30, 100, "3H10M"),
                    read("a0", 40, 20, 100, "5H10M5H"),
                    read("a1", 40, 20, 100, "10M2H"),
                ],
                2,
            )),
        ),
        "leading-insertion" => (
            site(true),
            Some(split(
                vec![
                    read("r0", 60, 30, 100, "10M"),
                    read("a0", 40, 20, START as i32 + 1, "3I7M"),
                ],
                1,
            )),
        ),
        other => panic!("{other} has no fixture"),
    }
}

fn rendered(entries: Vec<(String, AnnotationValue)>) -> String {
    entries
        .iter()
        .map(|(key, value)| {
            let AnnotationValue::Str(text) = value else {
                panic!("{key} is not a String");
            };
            format!("{key}={text}[java.lang.String]")
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn one(
    name: &str,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> String {
    match name {
        "BaseQualityRankSumTest" => {
            rendered(annotate(&BaseQualityRankSumTest, None, vc, likelihoods))
        }
        "MappingQualityRankSumTest" => {
            rendered(annotate(&MappingQualityRankSumTest, None, vc, likelihoods))
        }
        "ReadPosRankSumTest" => rendered(annotate(&ReadPosRankSumTest, None, vc, likelihoods)),
        "ClippingRankSumTest" => rendered(annotate(&ClippingRankSumTest, None, vc, likelihoods)),
        other => panic!("unknown annotation {other}"),
    }
}

#[test]
fn every_annotation_answers_as_the_reference_answers() {
    let text = golden();
    let mut count = 0;
    let mut absent = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("anno\t") else {
            continue;
        };
        let mut fields = rest.splitn(3, '\t');
        let name = fields.next().expect("an annotation");
        let label = fields.next().expect("a label");
        let expected = fields.next().unwrap_or("");
        if expected.is_empty() {
            absent += 1;
        }
        let (vc, likelihoods) = case(label);
        assert_eq!(
            one(name, &vc, likelihoods.as_ref()),
            expected,
            "{name} on {label}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no annotation rows");
    println!("{count} annotation answers identical, {absent} of them absent keys");
}
