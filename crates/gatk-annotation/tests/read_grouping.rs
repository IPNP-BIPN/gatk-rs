//! Conformance for `AS_UNIQ_ALT_READ_COUNT`, `BQHIST` and `REF_BASES`, against the oracle.
//!
//! Golden from `tools/annotation-conformance/ReadGroupingAnnotationDump.java`.
//!
//! ```text
//! anno  UniqueAltReadCount    three-duplicates  AS_UNIQ_ALT_READ_COUNT=1   three reads, one fragment
//! anno  UniqueAltReadCount    two-alternates    AS_UNIQ_ALT_READ_COUNT=1|1 joined with a pipe
//! anno  BaseQualityHistogram  ref-and-alt       BQHIST=[30, 1, 0, 31, 1, 0, 32, 0, 1, 33, 0, 1]
//! refbases  short-window                        ACGTACGTACGTNNNNNNNNN     padded on the right
//! ```

use std::io::Read;

use gatk_annotation::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_annotation::read_grouping::{BaseQualityHistogram, ReferenceBases, UniqueAltReadCount};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

const START: i64 = 105;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_grouping.txt.gz");
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

fn site(triallelic: bool) -> VariantContext {
    let alleles = if triallelic {
        vec![allele("A", true), allele("C", false), allele("G", false)]
    } else {
        vec![allele("A", true), allele("C", false)]
    };
    let mut vc = VariantContext::new("chr1", START, alleles);
    vc.stop = START;
    vc
}

fn read(name: &str, start: i32, fragment: i32, base_quality: u8, mapping_quality: u8) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality,
        cigar: htsjdk_bam::text_parse::parse_cigar("20M").expect("a cigar"),
        read_bases: vec![b'A'; 20],
        base_qualities: vec![base_quality; 20],
        inferred_insert_size: fragment,
        ..Default::default()
    }
}

fn matrix(
    reads: Vec<BamRecord>,
    alleles: Vec<Allele>,
    values: Vec<Vec<f64>>,
) -> AlleleLikelihoods<BamRecord> {
    AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        AlleleList::new(&alleles),
        vec![reads],
        vec![values],
    )
    .expect("a matrix")
}

/// Every read supports the alternate.
fn on_alt(pairs: &[(i32, i32)]) -> AlleleLikelihoods<BamRecord> {
    let reads: Vec<BamRecord> = pairs
        .iter()
        .enumerate()
        .map(|(i, (start, fragment))| read(&format!("r{i}"), *start, *fragment, 30, 60))
        .collect();
    let count = reads.len();
    matrix(
        reads,
        vec![allele("A", true), allele("C", false)],
        vec![vec![-10.0; count], vec![-1.0; count]],
    )
}

fn case(label: &str) -> (VariantContext, AlleleLikelihoods<BamRecord>) {
    match label {
        "empty-matrix" => (
            site(false),
            matrix(
                Vec::new(),
                vec![allele("A", true), allele("C", false)],
                vec![Vec::new(), Vec::new()],
            ),
        ),
        "three-distinct-fragments" => (site(false), on_alt(&[(100, 300), (101, 300), (102, 300)])),
        "three-duplicates" => (site(false), on_alt(&[(100, 300), (100, 300), (100, 300)])),
        "same-start-different-length" => (site(false), on_alt(&[(100, 300), (100, 301)])),
        "mixed" => (
            site(false),
            on_alt(&[(100, 300), (100, 300), (101, 300), (101, 300), (102, 400)]),
        ),
        "ref-and-alt" => {
            let reads = vec![
                read("r0", 100, 300, 30, 60),
                read("r1", 100, 300, 31, 60),
                read("a0", 101, 300, 32, 60),
                read("a1", 102, 300, 33, 60),
            ];
            (
                site(false),
                matrix(
                    reads,
                    vec![allele("A", true), allele("C", false)],
                    vec![
                        vec![-1.0, -1.0, -10.0, -10.0],
                        vec![-10.0, -10.0, -1.0, -1.0],
                    ],
                ),
            )
        }
        "two-alternates" => {
            let reads = vec![
                read("r0", 100, 300, 30, 60),
                read("a0", 101, 300, 31, 60),
                read("b0", 102, 300, 32, 60),
            ];
            let values = (0..3)
                .map(|a| (0..3).map(|e| if a == e { -1.0 } else { -10.0 }).collect())
                .collect();
            (
                site(true),
                matrix(
                    reads,
                    vec![allele("A", true), allele("C", false), allele("G", false)],
                    values,
                ),
            )
        }
        "varied-qualities" => {
            let reads = vec![
                read("r0", 100, 300, 20, 60),
                read("r1", 100, 301, 30, 60),
                read("a0", 101, 300, 20, 60),
                read("a1", 101, 301, 40, 60),
            ];
            (
                site(false),
                matrix(
                    reads,
                    vec![allele("A", true), allele("C", false)],
                    vec![
                        vec![-1.0, -1.0, -10.0, -10.0],
                        vec![-10.0, -10.0, -1.0, -1.0],
                    ],
                ),
            )
        }
        "mapq-zero-dropped" => {
            let reads = vec![
                read("r0", 100, 300, 30, 0),
                read("a0", 101, 300, 31, 0),
                read("a1", 102, 300, 32, 60),
            ];
            (
                site(false),
                matrix(
                    reads,
                    vec![allele("A", true), allele("C", false)],
                    vec![vec![-1.0, -10.0, -10.0], vec![-10.0, -1.0, -1.0]],
                ),
            )
        }
        other => panic!("{other} has no fixture"),
    }
}

fn rendered(entries: Vec<(String, AnnotationValue)>) -> String {
    entries
        .iter()
        .map(|(key, value)| match value {
            AnnotationValue::Str(text) => format!("{key}={text}[java.lang.String]"),
            AnnotationValue::List(values) => {
                let numbers: Vec<String> = values
                    .iter()
                    .map(|v| match v {
                        AnnotationValue::Int(n) => n.to_string(),
                        other => panic!("{other:?} is not an int"),
                    })
                    .collect();
                // `ArrayList.toString`.
                format!("{key}=[{}][java.util.ArrayList]", numbers.join(", "))
            }
            other => panic!("{other:?} has no rendering"),
        })
        .collect::<Vec<_>>()
        .join(";")
}

#[test]
fn every_annotation_answers_as_the_reference_answers() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("anno\t") else {
            continue;
        };
        let mut fields = rest.splitn(3, '\t');
        let name = fields.next().expect("an annotation");
        let label = fields.next().expect("a label");
        let expected = fields.next().unwrap_or("");
        let (vc, likelihoods) = case(label);
        let ours = match name {
            "UniqueAltReadCount" => {
                rendered(UniqueAltReadCount.annotate(None, &vc, Some(&likelihoods)))
            }
            "BaseQualityHistogram" => {
                rendered(BaseQualityHistogram.annotate(None, &vc, Some(&likelihoods)))
            }
            other => panic!("unknown annotation {other}"),
        };
        assert_eq!(ours, expected, "{name} on {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no annotation rows");
    println!("{count} annotation answers identical");
}

#[test]
fn every_reference_window_produces_the_same_string() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("refbases\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let _label = fields.next().expect("a label");
        let window_start: i64 = fields.next().expect("a start").parse().expect("a number");
        let bases = fields.next().expect("the bases");
        let variant_start: i64 = fields.next().expect("a start").parse().expect("a number");
        let expected = fields.next().expect("the result");
        let mut vc = site(false);
        vc.start = variant_start;
        vc.stop = variant_start;
        assert_eq!(
            ReferenceBases::local_bases(window_start, bases.as_bytes(), &vc),
            expected,
            "window at {window_start} for a variant at {variant_start}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no refbases rows");
    println!("{count} reference windows identical");
}
