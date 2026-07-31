//! Conformance for `FS`, `SOR` and `SB`, against the oracle.
//!
//! Golden from `tools/annotation-conformance/StrandBiasAnnotationDump.java`.
//!
//! ```text
//! anno    FisherStrand      skewed       FS=49.656[java.lang.String]
//! anno    FisherStrand      deep-skewed  FS=580.560     after the table was normalised
//! anno    StrandBiasBySample monomorphic SB=[5, 5, 0, 0][java.util.ArrayList]
//! ```
//!
//! `SB` is written at a **monomorphic** site where `FS` and `SOR` write nothing: the two info
//! annotations return early on `!vc.isVariant()` and the genotype one does not test the site.

use std::io::Read;

use gatk_annotation::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_annotation::strand_bias::{
    calculate_sor, contingency_table, p_value_for_contingency_table, FisherStrand,
    StrandBiasBySample, StrandOddsRatio,
};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

const START: i64 = 105;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/strand_bias.txt.gz");
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

fn site(sb_field: Option<Value>) -> VariantContext {
    let mut vc = VariantContext::new("chr1", START, vec![allele("A", true), allele("C", false)]);
    vc.stop = START;
    let mut genotype = Genotype::new("s1", vec![allele("A", true), allele("C", false)]);
    if let Some(value) = sb_field {
        genotype.extended.push(("SB".to_string(), value));
    }
    vc.genotypes.push(genotype);
    vc
}

fn monomorphic() -> VariantContext {
    let mut vc = VariantContext::new("chr1", START, vec![allele("A", true)]);
    vc.stop = START;
    vc
}

fn read(name: &str, reverse: bool) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        flags: if reverse { 0x10 } else { 0 },
        reference_index: 0,
        alignment_start: 100,
        mapping_quality: 60,
        cigar: htsjdk_bam::text_parse::parse_cigar("10M").expect("a cigar"),
        read_bases: vec![b'A'; 10],
        base_qualities: vec![30; 10],
        ..Default::default()
    }
}

/// The dump's `balanced`: the four counts, in one sample.
fn balanced(
    ref_fwd: usize,
    ref_rev: usize,
    alt_fwd: usize,
    alt_rev: usize,
) -> AlleleLikelihoods<BamRecord> {
    let mut reads = Vec::new();
    let mut is_ref = Vec::new();
    for i in 0..ref_fwd {
        reads.push(read(&format!("rf{i}"), false));
        is_ref.push(true);
    }
    for i in 0..ref_rev {
        reads.push(read(&format!("rr{i}"), true));
        is_ref.push(true);
    }
    for i in 0..alt_fwd {
        reads.push(read(&format!("af{i}"), false));
        is_ref.push(false);
    }
    for i in 0..alt_rev {
        reads.push(read(&format!("ar{i}"), true));
        is_ref.push(false);
    }
    let values = (0..2)
        .map(|a| {
            is_ref
                .iter()
                .map(|r| if (a == 0) == *r { -1.0 } else { -10.0 })
                .collect()
        })
        .collect();
    AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        AlleleList::new(&[allele("A", true), allele("C", false)]),
        vec![reads],
        vec![values],
    )
    .expect("a matrix")
}

/// Two samples with one read each: below the `FS` threshold apart, above it pooled.
fn two_samples() -> AlleleLikelihoods<BamRecord> {
    AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string(), "s2".to_string()]),
        AlleleList::new(&[allele("A", true), allele("C", false)]),
        vec![vec![read("a", false)], vec![read("b", true)]],
        vec![vec![vec![-1.0], vec![-10.0]], vec![vec![-1.0], vec![-10.0]]],
    )
    .expect("a matrix")
}

fn case(label: &str) -> (VariantContext, Option<AlleleLikelihoods<BamRecord>>) {
    match label {
        "null-likelihoods" => (site(None), None),
        "monomorphic" => (monomorphic(), Some(balanced(5, 5, 5, 5))),
        "balanced" => (site(None), Some(balanced(5, 5, 5, 5))),
        "skewed" => (site(None), Some(balanced(10, 0, 0, 10))),
        "one-read" => (site(None), Some(balanced(1, 0, 0, 0))),
        "two-reads" => (site(None), Some(balanced(1, 1, 0, 0))),
        "three-reads" => (site(None), Some(balanced(2, 1, 0, 0))),
        "deep-balanced" => (site(None), Some(balanced(200, 200, 200, 200))),
        "deep-skewed" => (site(None), Some(balanced(400, 1, 1, 400))),
        "ref-only" => (site(None), Some(balanced(5, 5, 0, 0))),
        "alt-only" => (site(None), Some(balanced(0, 0, 5, 5))),
        "empty-matrix" => (site(None), Some(balanced(0, 0, 0, 0))),
        "sb-field-string" => (
            site(Some(Value::Str("1,2,3,4".to_string()))),
            Some(balanced(50, 50, 50, 50)),
        ),
        "sb-field-list" => (
            site(Some(Value::List(vec![
                Value::Int(4),
                Value::Int(3),
                Value::Int(2),
                Value::Int(1),
            ]))),
            Some(balanced(50, 50, 50, 50)),
        ),
        "sb-field-below-threshold" => (
            site(Some(Value::Str("0,1,0,0".to_string()))),
            Some(balanced(50, 50, 50, 50)),
        ),
        "two-samples-each-small" => (site(None), Some(two_samples())),
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
            "FisherStrand" => rendered(FisherStrand.annotate(None, &vc, likelihoods.as_ref())),
            "StrandOddsRatio" => {
                rendered(StrandOddsRatio.annotate(None, &vc, likelihoods.as_ref()))
            }
            "StrandBiasBySample" => {
                match StrandBiasBySample.counts(&vc, "s1", likelihoods.as_ref()) {
                    // `ArrayList.toString`, which is what the dump printed.
                    Some(counts) => format!(
                        "SB=[{}][java.util.ArrayList]",
                        counts
                            .iter()
                            .map(|c| c.to_string())
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                    None => String::new(),
                }
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
fn every_table_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("table\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let label = fields.next().expect("a label");
        let min_count: i32 = fields.next().expect("a count").parse().expect("a number");
        let expected = fields.next().expect("the cells");
        let (_, likelihoods) = case(label);
        let likelihoods = likelihoods.expect("a matrix");
        let table = contingency_table(&likelihoods, &site(None), min_count);
        assert_eq!(
            format!(
                "{},{},{},{}",
                table[0][0], table[0][1], table[1][0], table[1][1]
            ),
            expected,
            "{label} at minCount {min_count}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no table rows");
    println!("{count} contingency tables identical");
}

#[test]
fn every_statistic_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut fisher = 0;
    let mut sor = 0;
    for line in text.lines() {
        let (kind, rest) = if let Some(rest) = line.strip_prefix("fisher\t") {
            ("fisher", rest)
        } else if let Some(rest) = line.strip_prefix("sor\t") {
            ("sor", rest)
        } else {
            continue;
        };
        let (cells, expected) = rest.split_once('\t').expect("cells and a result");
        let numbers: Vec<i32> = cells
            .split(',')
            .map(|c| c.parse().expect("a number"))
            .collect();
        let table = [[numbers[0], numbers[1]], [numbers[2], numbers[3]]];
        let ours = if kind == "fisher" {
            fisher += 1;
            p_value_for_contingency_table(table)
        } else {
            sor += 1;
            calculate_sor(table)
        };
        let want = f64::from_bits(expected.parse::<i64>().expect("bits") as u64);
        assert_eq!(
            ours.to_bits(),
            want.to_bits(),
            "{kind}({cells}) = {ours:e}, reference {want:e}"
        );
    }
    assert!(
        fisher > 0 && sor > 0,
        "the golden carries no statistic rows"
    );
    println!("{fisher} Fisher p-values and {sor} odds ratios bit-identical");
}
