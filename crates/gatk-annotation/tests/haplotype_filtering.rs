//! Conformance for `HaplotypeFilteringAnnotation`, against the oracle.
//!
//! Golden from `tools/annotation-conformance/HaplotypeFilteringAnnotationDump.java`.
//!
//! ```text
//! keys                                       ASSEMBLED_HAPS,FILTERED_HAPS
//! anno  duplicate-bases-same-uniqueness      ASSEMBLED_HAPS=1   one entry: equal haplotypes
//! anno  duplicate-bases-different-uniqueness ASSEMBLED_HAPS=2   two: the field is in equals
//! anno  negative-filtered-count              FILTERED_HAPS=-1   the getter does not clamp
//! ```

use std::io::Read;

use gatk_annotation::haplotype_filtering::{
    HaplotypeFilteringAnnotation, HaplotypeLikelihoods, JumboInfoAnnotation,
};
use gatk_annotation::info_annotation::AnnotationValue;
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use gatk_engine::fragment::Fragment;
use gatk_engine::haplotype::Haplotype;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

const VARIANT_START: i64 = 105;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/haplotype_filtering.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn hap(bases: &str, is_ref: bool, uniqueness: i32) -> Haplotype {
    let mut haplotype = Haplotype::new(bases.as_bytes(), is_ref).expect("a haplotype");
    haplotype.set_uniqueness_value(uniqueness);
    haplotype
}

fn site() -> VariantContext {
    let alleles = vec![
        Allele::from_str("A", true).expect("an allele"),
        Allele::from_str("C", false).expect("an allele"),
    ];
    let mut vc = VariantContext::new("chr1", VARIANT_START, alleles);
    vc.stop = VARIANT_START;
    vc
}

fn read(name: &str, start: i32) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality: 60,
        cigar: htsjdk_bam::text_parse::parse_cigar("20M").expect("a cigar"),
        read_bases: vec![b'A'; 20],
        base_qualities: vec![30; 20],
        inferred_insert_size: 300,
        ..Default::default()
    }
}

/// The likelihood values the dump sets: allele 0 at -1 and the rest at -10. Nothing this
/// annotation reads depends on them, which is one of the things the suite pins.
fn values(allele_count: usize, evidence_count: usize) -> Vec<Vec<Vec<f64>>> {
    vec![(0..allele_count)
        .map(|allele| vec![if allele == 0 { -1.0 } else { -10.0 }; evidence_count])
        .collect()]
}

fn by_read(haplotypes: Vec<Haplotype>, filtered: i32) -> AlleleLikelihoods<BamRecord, Haplotype> {
    let reads = vec![read("r0", 100), read("r1", 101)];
    let alleles = AlleleList::new(&haplotypes);
    let count = alleles.number_of_alleles();
    let mut likelihoods = AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        alleles,
        vec![reads],
        values(count, 2),
    )
    .expect("a matrix");
    likelihoods.set_filtered_haplotype_count(filtered);
    likelihoods
}

fn by_fragment(
    haplotypes: Vec<Haplotype>,
    filtered: i32,
) -> AlleleLikelihoods<Fragment, Haplotype> {
    let fragments = vec![
        Fragment::create_and_avoid_failure(&[read("r0", 100)]).expect("a fragment"),
        Fragment::create_and_avoid_failure(&[read("r1", 101)]).expect("a fragment"),
    ];
    let alleles = AlleleList::new(&haplotypes);
    let count = alleles.number_of_alleles();
    let mut likelihoods = AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        alleles,
        vec![fragments],
        values(count, 2),
    )
    .expect("a matrix");
    likelihoods.set_filtered_haplotype_count(filtered);
    likelihoods
}

fn no_evidence(
    haplotypes: Vec<Haplotype>,
    filtered: i32,
) -> AlleleLikelihoods<BamRecord, Haplotype> {
    let alleles = AlleleList::new(&haplotypes);
    let count = alleles.number_of_alleles();
    let mut likelihoods = AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        alleles,
        vec![Vec::new()],
        values(count, 0),
    )
    .expect("a matrix");
    likelihoods.set_filtered_haplotype_count(filtered);
    likelihoods
}

fn three() -> Vec<Haplotype> {
    vec![
        hap("ACGT", true, 0),
        hap("ACGA", false, 0),
        hap("ACGC", false, 0),
    ]
}

/// `anno\t<label>\t<key>=<value>[<class>];...`, with the two keys sorted as the dump sorts them.
fn row(label: &str, likelihoods: &HaplotypeLikelihoods<'_>) -> String {
    let mut result = HaplotypeFilteringAnnotation.annotate(None, &site(), None, None, likelihoods);
    result.sort_by(|left, right| left.0.cmp(&right.0));
    let rendered: Vec<String> = result
        .iter()
        .map(|(key, value)| {
            format!(
                "{key}={}[{}]",
                value.to_java_string().expect("an integer renders"),
                value.java_class()
            )
        })
        .collect();
    format!("anno\t{label}\t{}", rendered.join(";"))
}

#[test]
fn matches_the_oracle() {
    let golden = golden();
    let expected: Vec<&str> = golden
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();

    let mut actual: Vec<String> = Vec::new();
    actual.push(format!(
        "keys\t{}",
        HaplotypeFilteringAnnotation.key_names().join(",")
    ));

    let cases: Vec<(&str, Vec<Haplotype>, i32)> = vec![
        ("no-haplotypes", Vec::new(), 0),
        ("one-haplotype", vec![hap("ACGT", false, 0)], 0),
        ("three-haplotypes", three(), 0),
        (
            "duplicate-bases-same-uniqueness",
            vec![hap("ACGT", false, 0), hap("ACGT", false, 0)],
            0,
        ),
        (
            "duplicate-bases-different-uniqueness",
            vec![hap("ACGT", false, 0), hap("ACGT", false, 1)],
            0,
        ),
        (
            "duplicate-bases-different-ref-flag",
            vec![hap("ACGT", true, 0), hap("ACGT", false, 0)],
            0,
        ),
        (
            "different-lengths",
            vec![hap("ACGT", false, 0), hap("ACGTACGT", false, 0)],
            0,
        ),
        ("unfiltered", three(), 0),
        ("two-filtered", three(), 2),
        ("filtered-exceeds-remaining", vec![hap("ACGT", true, 0)], 5),
        (
            "negative-filtered-count",
            vec![hap("ACGT", true, 0), hap("ACGA", false, 0)],
            -1,
        ),
    ];
    for (label, haplotypes, filtered) in cases {
        let likelihoods = by_read(haplotypes, filtered);
        actual.push(row(label, &HaplotypeLikelihoods::ByRead(&likelihoods)));
    }

    // The two branches of the engine's ternary, over the same haplotypes.
    let read_typed = by_read(three(), 1);
    actual.push(row("by-read", &HaplotypeLikelihoods::ByRead(&read_typed)));
    let fragment_typed = by_fragment(three(), 1);
    actual.push(row(
        "by-fragment",
        &HaplotypeLikelihoods::ByFragment(&fragment_typed),
    ));

    let empty = no_evidence(three(), 1);
    actual.push(row("no-evidence", &HaplotypeLikelihoods::ByRead(&empty)));

    assert_eq!(actual.len(), expected.len(), "row count");
    for (produced, oracle) in actual.iter().zip(expected.iter()) {
        assert_eq!(produced, oracle);
    }
}

/// The two counts a caller reads, without going through the dump's rendering.
#[test]
fn duplicate_haplotypes_collapse_only_when_equal() {
    let same = by_read(vec![hap("ACGT", false, 0), hap("ACGT", false, 0)], 0);
    let different = by_read(vec![hap("ACGT", false, 0), hap("ACGT", false, 1)], 0);
    assert_eq!(same.number_of_alleles(), 1);
    assert_eq!(different.number_of_alleles(), 2);
    // The hash is the bases alone, so the pair that equality separates still collides.
    assert_eq!(
        hap("ACGT", false, 0).java_hash_code(),
        hap("ACGT", false, 1).java_hash_code()
    );
}

/// An unset count is a zero rather than a missing key: both keys are always written.
#[test]
fn an_unfiltered_matrix_reports_zero() {
    let likelihoods = by_read(three(), 0);
    let result = HaplotypeFilteringAnnotation.annotate(
        None,
        &site(),
        None,
        None,
        &HaplotypeLikelihoods::ByRead(&likelihoods),
    );
    assert_eq!(result.len(), 2);
    assert_eq!(result[1].1, AnnotationValue::Int(0));
}
