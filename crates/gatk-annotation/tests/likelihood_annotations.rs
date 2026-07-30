//! Conformance for the three annotations whose whole input is the likelihood matrix.
//!
//! Golden from `tools/annotation-conformance/LikelihoodAnnotationDump.java`.
//!
//! The rows that separate them, all three on the same input:
//!
//! ```text
//! anno  Coverage            empty-matrix      (nothing)
//! anno  MappingQualityZero  empty-matrix      MQ0=0[java.lang.String]
//! anno  CountNs             empty-matrix      NCount=0[java.lang.Long]
//! ```
//!
//! Only `Coverage` tests the evidence count, only `MappingQualityZero` tests the site, and the
//! Java type of the value differs between two annotations that both write a count.

use std::io::Read;

use gatk_annotation::coverage::{CountNs, Coverage, MappingQualityZero};
use gatk_annotation::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

/// The dump's variant start: ten bases into a read that starts at 100.
const VARIANT_START: i64 = 105;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/likelihood_annotations.txt.gz");
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

fn variant_site() -> VariantContext {
    let mut vc = VariantContext::new(
        "chr1",
        VARIANT_START,
        vec![allele("A", true), allele("C", false)],
    );
    vc.stop = VARIANT_START;
    vc
}

fn reference_only_site() -> VariantContext {
    let mut vc = VariantContext::new("chr1", VARIANT_START, vec![allele("A", true)]);
    vc.stop = VARIANT_START;
    vc
}

fn read(name: &str, mapping_quality: u8, bases: &str, start: i32, cigar_text: &str) -> BamRecord {
    let cigar = htsjdk_bam::text_parse::parse_cigar(cigar_text).expect("a parsable cigar");
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality,
        read_bases: bases.as_bytes().to_vec(),
        base_qualities: vec![30; bases.len()],
        cigar,
        ..Default::default()
    }
}

/// A matrix over one sample, with the likelihoods left at zero: none of the three reads them.
fn likelihoods(reads: Vec<BamRecord>) -> AlleleLikelihoods<BamRecord> {
    matrix(vec![("s1", reads)])
}

fn matrix(samples: Vec<(&str, Vec<BamRecord>)>) -> AlleleLikelihoods<BamRecord> {
    let names: Vec<String> = samples
        .iter()
        .map(|(name, _)| (*name).to_string())
        .collect();
    let alleles = AlleleList::new(&[allele("A", true), allele("C", false)]);
    let values: Vec<Vec<Vec<f64>>> = samples
        .iter()
        .map(|(_, reads)| vec![vec![0.0; reads.len()], vec![0.0; reads.len()]])
        .collect();
    let evidence: Vec<Vec<BamRecord>> = samples.into_iter().map(|(_, reads)| reads).collect();
    AlleleLikelihoods::new(SampleList::new(&names), alleles, evidence, values)
        .expect("a well-formed matrix")
}

struct Case {
    label: &'static str,
    vc: VariantContext,
    likelihoods: Option<AlleleLikelihoods<BamRecord>>,
}

fn cases() -> Vec<Case> {
    let case = |label, vc, likelihoods| Case {
        label,
        vc,
        likelihoods,
    };
    vec![
        case("null-likelihoods", variant_site(), None),
        case("empty-matrix", variant_site(), Some(likelihoods(vec![]))),
        case(
            "sample-without-evidence",
            variant_site(),
            Some(matrix(vec![("s1", vec![]), ("s2", vec![])])),
        ),
        case(
            "three-plain-reads",
            variant_site(),
            Some(likelihoods(vec![
                read("r0", 60, "ACGTACGTAC", 100, "10M"),
                read("r1", 60, "ACGTACGTAC", 100, "10M"),
                read("r2", 60, "ACGTACGTAC", 100, "10M"),
            ])),
        ),
        case(
            "two-mapq-zero",
            variant_site(),
            Some(likelihoods(vec![
                read("r0", 0, "ACGTACGTAC", 100, "10M"),
                read("r1", 0, "ACGTACGTAC", 100, "10M"),
                read("r2", 60, "ACGTACGTAC", 100, "10M"),
            ])),
        ),
        case(
            "one-n-at-start",
            variant_site(),
            Some(likelihoods(vec![
                read("r0", 60, "ACGTANGTAC", 100, "10M"),
                read("r1", 60, "ACGTACGTAC", 100, "10M"),
                read("r2", 60, "ACGTACGTAC", 100, "10M"),
            ])),
        ),
        case(
            "lower-case-n",
            variant_site(),
            Some(likelihoods(vec![read("r0", 60, "ACGTAnGTAC", 100, "10M")])),
        ),
        case(
            "n-beside-the-start",
            variant_site(),
            Some(likelihoods(vec![read("r0", 60, "ACGTNCGTAC", 100, "10M")])),
        ),
        case(
            "n-in-soft-clip",
            variant_site(),
            Some(likelihoods(vec![read("r0", 60, "NNNNNNACGT", 106, "6S4M")])),
        ),
        case(
            "read-past-the-variant",
            variant_site(),
            Some(likelihoods(vec![read("r0", 60, "ACGTACGTAC", 200, "10M")])),
        ),
        case(
            "deletion-over-start",
            variant_site(),
            Some(likelihoods(vec![read(
                "r0",
                60,
                "ACGTACGTAC",
                100,
                "5M3D5M",
            )])),
        ),
        case(
            "two-samples",
            variant_site(),
            Some(matrix(vec![
                (
                    "s1",
                    vec![
                        read("a0", 0, "ACGTANGTAC", 100, "10M"),
                        read("a1", 60, "ACGTACGTAC", 100, "10M"),
                    ],
                ),
                ("s2", vec![read("b0", 0, "ACGTACGTAC", 100, "10M")]),
            ])),
        ),
        case(
            "monomorphic-site",
            reference_only_site(),
            Some(likelihoods(vec![read("r0", 0, "ACGTANGTAC", 100, "10M")])),
        ),
    ]
}

/// The dump's rendering of one map.
fn render(entries: &[(String, AnnotationValue)]) -> String {
    entries
        .iter()
        .map(|(key, value)| {
            format!(
                "{key}={}[{}]",
                value.to_java_string().expect("no Doubles here"),
                value.java_class()
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn value(text: &str, annotation: &str, label: &str) -> String {
    let needle = format!("anno\t{annotation}\t{label}\t");
    text.lines()
        .find(|line| line.starts_with(&needle))
        .unwrap_or_else(|| panic!("no row for {annotation}/{label}"))[needle.len()..]
        .to_string()
}

#[test]
fn coverage_matches_the_reference() {
    let text = golden();
    for case in cases() {
        let ours = render(&Coverage.annotate(None, &case.vc, case.likelihoods.as_ref()));
        assert_eq!(
            ours,
            value(&text, "Coverage", case.label),
            "Coverage on {}",
            case.label
        );
    }
}

#[test]
fn mapping_quality_zero_matches_the_reference() {
    let text = golden();
    for case in cases() {
        let ours = render(&MappingQualityZero.annotate(None, &case.vc, case.likelihoods.as_ref()));
        assert_eq!(
            ours,
            value(&text, "MappingQualityZero", case.label),
            "MappingQualityZero on {}",
            case.label
        );
    }
}

#[test]
fn count_ns_matches_the_reference() {
    let text = golden();
    for case in cases() {
        let ours = render(&CountNs.annotate(None, &case.vc, case.likelihoods.as_ref()));
        assert_eq!(
            ours,
            value(&text, "CountNs", case.label),
            "CountNs on {}",
            case.label
        );
    }
}

/// The rows where the three disagree on the same input.
#[test]
fn the_three_disagree_on_the_same_input() {
    let text = golden();

    // An empty matrix: only Coverage tests the evidence count.
    assert_eq!(value(&text, "Coverage", "empty-matrix"), "");
    assert_eq!(
        value(&text, "MappingQualityZero", "empty-matrix"),
        "MQ0=0[java.lang.String]"
    );
    assert_eq!(
        value(&text, "CountNs", "empty-matrix"),
        "NCount=0[java.lang.Long]"
    );

    // A non-variant site: only MappingQualityZero tests the site.
    assert_eq!(value(&text, "MappingQualityZero", "monomorphic-site"), "");
    assert!(value(&text, "Coverage", "monomorphic-site").starts_with("DP=1"));
    assert!(value(&text, "CountNs", "monomorphic-site").starts_with("NCount=1"));

    // Two counts, two Java types, both rendering the same way.
    assert!(value(&text, "Coverage", "three-plain-reads").ends_with("[java.lang.String]"));
    assert!(value(&text, "CountNs", "three-plain-reads").ends_with("[java.lang.Long]"));

    // The N tests: upper case only, at the start only, and not through a soft clip.
    assert_eq!(
        value(&text, "CountNs", "one-n-at-start"),
        "NCount=1[java.lang.Long]"
    );
    assert_eq!(
        value(&text, "CountNs", "lower-case-n"),
        "NCount=0[java.lang.Long]"
    );
    assert_eq!(
        value(&text, "CountNs", "n-beside-the-start"),
        "NCount=0[java.lang.Long]"
    );
    assert_eq!(
        value(&text, "CountNs", "n-in-soft-clip"),
        "NCount=0[java.lang.Long]"
    );
}
