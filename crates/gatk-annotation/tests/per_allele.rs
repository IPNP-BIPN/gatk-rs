//! Conformance for `MBQ`, `MMQ`, `MPOS` and `MFRL`, against the oracle.
//!
//! Golden from `tools/annotation-conformance/PerAlleleAnnotationDump.java`.
//!
//! The golden corrected the port on one point, and it is the shape of the data structure rather
//! than the arithmetic:
//!
//! ```text
//! anno  BaseQuality  allele-missing-from-matrix  E:java.lang.NullPointerException:...
//! ```
//!
//! The buckets are keyed by the **matrix's** alleles and the loop runs over the **variant's**, so
//! a variant allele the matrix never held looks up `null` and `aggregate` dereferences it. The
//! port had it as an empty list, which would have produced the annotation's value for no reads:
//! a plausible number where the reference produces no record at all.
//!
//! The rows that carry the rest:
//!
//! ```text
//! anno  MappingQuality  empty-matrix       MMQ=[60, 60]   the invented value, not zero
//! anno  MappingQuality  median-on-a-half   MMQ=[60, 21]   20 and 21 round up, not to even
//! anno  ReadPosition    hard-clipped       MPOS=[8]       hard clips count as bases
//! anno  ReadPosition    three-alleles      MPOS=[4, 50]   two numbers, not three
//! ```

use std::io::Read;

use gatk_annotation::info_annotation::AnnotationValue;
use gatk_annotation::per_allele::{
    annotate, BaseQuality, FragmentLength, MappingQuality, MissingBucket, PerAlleleAnnotation,
    ReadPosition,
};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

const START: i64 = 105;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/per_allele.txt.gz");
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

fn site(alleles: Vec<Allele>) -> VariantContext {
    let mut vc = VariantContext::new("chr1", START, alleles);
    vc.stop = START;
    vc
}

/// The dump's `read`, with the bases and qualities its cigar implies.
fn read(
    name: &str,
    mapping_quality: u8,
    base_quality: u8,
    start: i32,
    cigar: &str,
    fragment_length: i32,
) -> BamRecord {
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
        inferred_insert_size: fragment_length,
        ..Default::default()
    }
}

fn likelihoods_for(
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

/// Every read supports one named allele strongly, as the dump's `supporting` builds it.
fn supporting(target: usize, reads: Vec<BamRecord>) -> AlleleLikelihoods<BamRecord> {
    let alleles = vec![allele("A", true), allele("C", false)];
    let values = (0..alleles.len())
        .map(|a| {
            reads
                .iter()
                .map(|_| if a == target { -1.0 } else { -10.0 })
                .collect()
        })
        .collect();
    likelihoods_for(reads, alleles, values)
}

/// The matrix each label was measured against, and the variant it was measured at.
fn case(label: &str) -> (VariantContext, Option<AlleleLikelihoods<BamRecord>>) {
    let biallelic = vec![allele("A", true), allele("C", false)];
    let triallelic = vec![allele("A", true), allele("C", false), allele("G", false)];
    match label {
        "null-likelihoods" => (site(biallelic), None),
        "empty-matrix" => (
            site(biallelic.clone()),
            Some(likelihoods_for(
                Vec::new(),
                biallelic,
                vec![Vec::new(), Vec::new()],
            )),
        ),
        "one-read-each" => (
            site(biallelic.clone()),
            Some(likelihoods_for(
                vec![
                    read("ref0", 30, 20, 100, "10M", 200),
                    read("alt0", 50, 40, 100, "10M", 600),
                ],
                biallelic,
                vec![vec![-1.0, -10.0], vec![-10.0, -1.0]],
            )),
        ),
        "three-on-alt" => (
            site(biallelic),
            Some(supporting(
                1,
                vec![
                    read("r0", 20, 30, 100, "10M", 300),
                    read("r1", 40, 35, 100, "10M", 400),
                    read("r2", 60, 40, 100, "10M", 500),
                ],
            )),
        ),
        "four-on-alt" => (
            site(biallelic),
            Some(supporting(
                1,
                vec![
                    read("r0", 20, 30, 100, "10M", 300),
                    read("r1", 30, 31, 100, "10M", 301),
                    read("r2", 40, 32, 100, "10M", 302),
                    read("r3", 50, 33, 100, "10M", 303),
                ],
            )),
        ),
        "median-on-a-half" => (
            site(biallelic),
            Some(supporting(
                1,
                vec![
                    read("r0", 20, 30, 100, "10M", 300),
                    read("r1", 21, 31, 100, "10M", 301),
                ],
            )),
        ),
        "mapq-zero" => (
            site(biallelic),
            Some(supporting(
                1,
                vec![
                    read("r0", 0, 30, 100, "10M", 300),
                    read("r1", 40, 35, 100, "10M", 400),
                ],
            )),
        ),
        "mapq-unavailable" => (
            site(biallelic),
            Some(supporting(
                1,
                vec![
                    read("r0", 255, 30, 100, "10M", 300),
                    read("r1", 40, 35, 100, "10M", 400),
                ],
            )),
        ),
        "uninformative" => (
            site(biallelic.clone()),
            Some(likelihoods_for(
                vec![read("r0", 40, 35, 100, "10M", 400)],
                biallelic,
                vec![vec![-1.0], vec![-1.1]],
            )),
        ),
        "read-past-the-variant" => (
            site(biallelic),
            Some(supporting(1, vec![read("r0", 40, 35, 200, "10M", 400)])),
        ),
        "deletion-over-start" => (
            site(biallelic),
            Some(supporting(1, vec![read("r0", 40, 35, 100, "5M3D5M", 400)])),
        ),
        "hard-clipped" => (
            site(biallelic),
            Some(supporting(1, vec![read("r0", 40, 35, 100, "3H10M5H", 400)])),
        ),
        "soft-clipped-over-start" => (
            site(biallelic),
            Some(supporting(1, vec![read("r0", 40, 35, 106, "6S4M", 400)])),
        ),
        "negative-fragment-length" => (
            site(biallelic),
            Some(supporting(
                1,
                vec![
                    read("r0", 40, 35, 100, "10M", -400),
                    read("r1", 40, 35, 100, "10M", -300),
                ],
            )),
        ),
        "three-alleles" => (
            site(triallelic.clone()),
            Some(likelihoods_for(
                vec![read("r0", 40, 35, 100, "10M", 400)],
                triallelic,
                vec![vec![-10.0], vec![-1.0], vec![-10.0]],
            )),
        ),
        "allele-missing-from-matrix" => (
            site(triallelic),
            Some(likelihoods_for(
                vec![read("r0", 40, 35, 100, "10M", 400)],
                vec![allele("A", true), allele("C", false)],
                vec![vec![-10.0], vec![-1.0]],
            )),
        ),
        "two-samples" => (
            site(biallelic.clone()),
            Some(
                AlleleLikelihoods::new(
                    SampleList::new(&["s1".to_string(), "s2".to_string()]),
                    AlleleList::new(&biallelic),
                    vec![
                        vec![read("a0", 40, 30, 100, "10M", 300)],
                        vec![read("b0", 50, 40, 100, "10M", 500)],
                    ],
                    vec![vec![vec![-10.0], vec![-1.0]], vec![vec![-10.0], vec![-1.0]]],
                )
                .expect("a matrix"),
            ),
        ),
        other => panic!("{other} has no fixture"),
    }
}

/// The dump's rendering of one annotation's result.
fn rendered(result: Result<Vec<(String, AnnotationValue)>, MissingBucket>) -> String {
    match result {
        Err(error) => format!("E:{}:{}", error.class(), error.message()),
        Ok(entries) => entries
            .iter()
            .map(|(key, value)| {
                let AnnotationValue::List(values) = value else {
                    panic!("{key} is not an int[]");
                };
                let numbers: Vec<String> = values
                    .iter()
                    .map(|v| match v {
                        AnnotationValue::Int(n) => n.to_string(),
                        other => panic!("{other:?} is not an int"),
                    })
                    .collect();
                // `Arrays.toString(int[])` and the class of a Java int array.
                format!("{key}=[{}][[I]", numbers.join(", "))
            })
            .collect::<Vec<_>>()
            .join(";"),
    }
}

fn one(
    name: &str,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> String {
    match name {
        "BaseQuality" => rendered(annotate(&BaseQuality, None, vc, likelihoods)),
        "MappingQuality" => rendered(annotate(&MappingQuality, None, vc, likelihoods)),
        "ReadPosition" => rendered(annotate(&ReadPosition, None, vc, likelihoods)),
        "FragmentLength" => rendered(annotate(&FragmentLength, None, vc, likelihoods)),
        other => panic!("unknown annotation {other}"),
    }
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
        assert_eq!(
            one(name, &vc, likelihoods.as_ref()),
            expected,
            "{name} on {label}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no annotation rows");
    println!("{count} annotation answers identical");
}

/// The four values for no reads, which are four different numbers and the reason the parent
/// cannot be factored.
#[test]
fn each_member_invents_its_own_value_for_an_allele_with_no_reads() {
    assert_eq!(BaseQuality.value_for_no_reads(), 0);
    assert_eq!(MappingQuality.value_for_no_reads(), 60);
    assert_eq!(ReadPosition.value_for_no_reads(), 50);
    assert_eq!(FragmentLength.value_for_no_reads(), 0);
    // And only one of the four leaves the reference allele out.
    assert!(BaseQuality.include_ref_allele());
    assert!(MappingQuality.include_ref_allele());
    assert!(!ReadPosition.include_ref_allele());
    assert!(FragmentLength.include_ref_allele());
}
