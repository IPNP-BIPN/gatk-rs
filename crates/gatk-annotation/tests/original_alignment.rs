//! Conformance for `OriginalAlignment` and the `OA` tag it reads, against the oracle.
//!
//! Golden from `tools/annotation-conformance/OriginalAlignmentDump.java`.
//!
//! Two rows carry the annotation's actual behaviour:
//!
//! ```text
//! anno  unmapped-oa   OCM=1[java.lang.Long]
//! anno  tlod-missing  OCM=1[java.lang.Long]
//! ```
//!
//! An unmapped read's tag is `*,0,*,*,0,0;`, so its original contig is the string `*`, which
//! differs from every real contig and is therefore counted as a mismatch. And a `TLOD` of `.` is
//! not treated as missing: it becomes `-1`, then `-ln(10)`, which wins the maximum of a
//! one-element array, so the annotation still picks an allele and still counts.

use std::io::Read;

use gatk_annotation::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_annotation::original_alignment::{oa_contig, OriginalAlignment, OA_TAG};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use gatk_engine::context::ReferenceContext;
use gatk_engine::interval::SimpleInterval;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Value, VariantContext};

const START: i64 = 105;

/// Likelihoods making the single read informative for the first alternate allele, and for the
/// second, and for neither.
const FOR_ALT1: [f64; 3] = [-5.0, 0.0, -5.0];
const FOR_ALT2: [f64; 3] = [-5.0, -5.0, 0.0];
const UNINFORMATIVE: [f64; 3] = [0.0, 0.0, 0.0];

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/original_alignment.txt.gz");
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

fn site(tlod: Option<&str>) -> VariantContext {
    let mut vc = VariantContext::new(
        "chr1",
        START,
        vec![allele("A", true), allele("C", false), allele("G", false)],
    );
    vc.stop = START;
    if let Some(tlod) = tlod {
        vc.attributes
            .push(("TLOD".to_string(), Value::Str(tlod.to_string())));
    }
    vc
}

fn read(name: &str, oa: Option<&str>) -> BamRecord {
    let mut record = BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: 100,
        mapping_quality: 60,
        read_bases: b"ACGTACGTAC".to_vec(),
        base_qualities: vec![30; 10],
        cigar: htsjdk_bam::text_parse::parse_cigar("10M").expect("a cigar"),
        ..Default::default()
    };
    if let Some(oa) = oa {
        record
            .tags
            .insert(Tag::new(&OA_TAG), TagValue::Str(oa.to_string()));
    }
    record
}

fn matrix(reads: Vec<BamRecord>, per_allele: [f64; 3]) -> AlleleLikelihoods<BamRecord> {
    let count = reads.len();
    let values = vec![vec![
        vec![per_allele[0]; count],
        vec![per_allele[1]; count],
        vec![per_allele[2]; count],
    ]];
    AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        AlleleList::new(&[allele("A", true), allele("C", false), allele("G", false)]),
        vec![reads],
        values,
    )
    .expect("a well-formed matrix")
}

/// The dump builds `new ReferenceContext((ReferenceDataSource) null, interval)`, which is the
/// source-less context: it knows where it is and has no bases, which is all this annotation asks.
fn reference() -> ReferenceContext {
    ReferenceContext::without_source(
        Some(SimpleInterval::new("chr1", START as i32, START as i32).expect("an interval")),
        0,
        0,
    )
    .expect("a source-less context")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let needle = format!("{kind}\t{label}\t");
    text.lines()
        .find(|line| line.starts_with(&needle))
        .unwrap_or_else(|| panic!("no {kind} row for {label}"))[needle.len()..]
        .to_string()
}

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

/// One fixture: a label, the site's TLOD, the reads, and the per-allele likelihood every read
/// gets.
struct Case {
    label: &'static str,
    tlod: Option<&'static str>,
    reads: Vec<BamRecord>,
    per_allele: [f64; 3],
}

/// The dump's fixtures, in its order, minus the one whose answer is a throw.
fn cases() -> Vec<Case> {
    let case = |label, tlod, reads, per_allele| Case {
        label,
        tlod,
        reads,
        per_allele,
    };
    vec![
        case(
            "no-tlod",
            None,
            vec![read("r0", Some("chr2,100,+,10M,60,0;"))],
            FOR_ALT1,
        ),
        case(
            "one-mismatch",
            Some("10.0"),
            vec![read("r0", Some("chr2,100,+,10M,60,0;"))],
            FOR_ALT1,
        ),
        case(
            "same-contig",
            Some("10.0"),
            vec![read("r0", Some("chr1,100,+,10M,60,0;"))],
            FOR_ALT1,
        ),
        case("no-oa-tag", Some("10.0"), vec![read("r0", None)], FOR_ALT1),
        case(
            "unmapped-oa",
            Some("10.0"),
            vec![read("r0", Some("*,0,*,*,0,0;"))],
            FOR_ALT1,
        ),
        case(
            "informative-for-alt2",
            Some("10.0,1.0"),
            vec![read("r0", Some("chr2,100,+,10M,60,0;"))],
            FOR_ALT2,
        ),
        case(
            "second-alt-wins",
            Some("1.0,10.0"),
            vec![read("r0", Some("chr2,100,+,10M,60,0;"))],
            FOR_ALT2,
        ),
        case(
            "tlod-tie",
            Some("10.0,10.0"),
            vec![read("r0", Some("chr2,100,+,10M,60,0;"))],
            FOR_ALT1,
        ),
        case(
            "tlod-missing",
            Some("."),
            vec![read("r0", Some("chr2,100,+,10M,60,0;"))],
            FOR_ALT1,
        ),
        case(
            "uninformative",
            Some("10.0"),
            vec![read("r0", Some("chr2,100,+,10M,60,0;"))],
            UNINFORMATIVE,
        ),
        case(
            "three-reads",
            Some("10.0"),
            vec![
                read("r0", Some("chr2,100,+,10M,60,0;")),
                read("r1", Some("chr1,100,+,10M,60,0;")),
                read("r2", Some("chr3,100,+,10M,60,0;")),
            ],
            FOR_ALT1,
        ),
    ]
}

#[test]
fn every_count_matches_the_reference() {
    let text = golden();
    for case in cases() {
        let label = case.label;
        let vc = site(case.tlod);
        let likelihoods = matrix(case.reads, case.per_allele);
        let ours = render(&OriginalAlignment.annotate(Some(&reference()), &vc, Some(&likelihoods)));
        assert_eq!(ours, value(&text, "anno", label), "OCM on {label}");
    }
}

#[test]
fn every_oa_contig_matches_the_reference() {
    let text = golden();
    let cases = [
        ("ordinary", "chr2,100,+,10M,60,0;"),
        ("unmapped", "*,0,*,*,0,0;"),
        ("contig-with-underscore", "chr_2,100,+,10M,60,0;"),
        ("no-comma", "chr2"),
        ("empty", ""),
        ("leading-comma", ",100,+,10M,60,0;"),
    ];
    for (label, tag) in cases {
        let record = read("probe", Some(tag));
        assert_eq!(
            oa_contig(&record).expect("a tag"),
            value(&text, "oacontig", label),
            "getOAContig on {label}"
        );
    }
}

/// The rows that decide what this annotation actually counts.
#[test]
fn the_rows_that_decide_what_is_counted() {
    let text = golden();

    // An unmapped read's original contig is the string "*", which differs from every real contig,
    // so a read that was unmapped before realignment counts as a mismatch.
    assert_eq!(value(&text, "oacontig", "unmapped"), "*");
    assert_eq!(value(&text, "anno", "unmapped-oa"), "OCM=1[java.lang.Long]");

    // A TLOD of "." is -1 and then -ln(10), an ordinary number that wins a one-element maximum,
    // so the annotation still picks an allele and still counts.
    assert_eq!(
        value(&text, "anno", "tlod-missing"),
        "OCM=1[java.lang.Long]"
    );
    // Whereas no TLOD at all is the one case that writes nothing.
    assert_eq!(value(&text, "anno", "no-tlod"), "");

    // A tie between two TLODs goes to the earliest alternate allele, so the read informative for
    // the first one is counted.
    assert_eq!(value(&text, "anno", "tlod-tie"), "OCM=1[java.lang.Long]");
    // And a read informative for the other alternate is not.
    assert_eq!(
        value(&text, "anno", "informative-for-alt2"),
        "OCM=0[java.lang.Long]"
    );

    // Uninformative reads are not counted however well their contigs mismatch.
    assert_eq!(
        value(&text, "anno", "uninformative"),
        "OCM=0[java.lang.Long]"
    );

    // The value is a Long, like CountNs and unlike Coverage.
    assert!(value(&text, "anno", "one-mismatch").ends_with("[java.lang.Long]"));

    // Null likelihoods throw here, where the other three likelihood-reading annotations answer an
    // empty map. The port's `count` takes the matrix by reference, so the state is unrepresentable
    // rather than guarded.
    assert!(value(&text, "anno", "null-likelihoods")
        .starts_with("E:java.lang.IllegalArgumentException"));
}
