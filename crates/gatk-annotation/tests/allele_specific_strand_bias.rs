//! Conformance for `AS_FS`, `AS_SOR` and the strand table they share, against the oracle.
//!
//! Golden from `tools/annotation-conformance/AlleleSpecificStrandBiasDump.java`.
//!
//! ```text
//! sbraw     two-reads    AS_SB_TABLE=0,0|0,0[java.lang.String]
//! sbcombine empty-entry  AS_SB_TABLE=10,8|[java.lang.String]
//! phred     0.0          4659327993935762513
//! ```
//!
//! The `sbraw` rows for one, two and three reads are the minimum-threshold rule: a sample with two
//! informative reads contributes nothing at all, so its counts vanish from the table even for the
//! allele that had both of them.

use std::io::Read;

use gatk_annotation::allele_specific_strand_bias::{
    annotate_direct, annotate_raw_data, combine_raw_data, finalize_raw_data,
    phred_scale_error_rate, AsStrandBias, AsStrandBiasError, AS_SB_TABLE_KEY,
};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

const START: i64 = 105;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/allele_specific_strand_bias.txt.gz");
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

fn reference() -> Allele {
    allele("A", true)
}

fn alternate() -> Allele {
    allele("C", false)
}

fn second_alternate() -> Allele {
    allele("G", false)
}

/// Per sample, the forward/reverse counts for the reference and for each alternate.
fn composition(label: &str) -> Option<Vec<Vec<[i32; 2]>>> {
    let table: Vec<Vec<[i32; 2]>> = match label {
        "balanced" => vec![vec![[5, 5], [5, 5]]],
        "skewed" => vec![vec![[10, 0], [0, 10]]],
        "ref-only" => vec![vec![[5, 5], [0, 0]]],
        "alt-only" => vec![vec![[0, 0], [5, 5]]],
        "one-read" => vec![vec![[1, 0], [0, 0]]],
        "two-reads" => vec![vec![[1, 1], [0, 0]]],
        "three-reads" => vec![vec![[2, 1], [0, 0]]],
        "two-samples-each-small" => vec![vec![[1, 1], [0, 0]], vec![[1, 1], [0, 0]]],
        "two-samples-one-large" => vec![vec![[1, 1], [0, 0]], vec![[5, 5], [5, 5]]],
        "multiallelic" => vec![vec![[6, 6], [3, 3], [2, 2]]],
        "all-forward" => vec![vec![[8, 0], [6, 0]]],
        "all-reverse" => vec![vec![[0, 8], [0, 6]]],
        "empty" => vec![vec![[0, 0], [0, 0]]],
        "null-likelihoods" => return None,
        other => panic!("{other} has no composition"),
    };
    Some(table)
}

fn alleles_for(label: &str) -> Vec<Allele> {
    if label == "multiallelic" {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    }
}

fn variant_context(label: &str) -> VariantContext {
    let mut vc = VariantContext::new("chr1", START, alleles_for(label));
    vc.stop = START;
    let samples = composition(label).map(|c| c.len()).unwrap_or(1);
    for s in 0..samples {
        vc.genotypes.push(Genotype::new(
            &format!("s{s}"),
            vec![reference(), alternate()],
        ));
    }
    vc
}

fn read(name: &str, reverse: bool) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        flags: if reverse { 0x10 } else { 0 },
        reference_index: 0,
        alignment_start: 100,
        mapping_quality: 60,
        cigar: htsjdk_bam::text_parse::parse_cigar("20M").expect("a cigar"),
        read_bases: vec![b'A'; 20],
        base_qualities: vec![30; 20],
        ..Default::default()
    }
}

fn likelihoods(label: &str) -> Option<AlleleLikelihoods<BamRecord>> {
    let composition = composition(label)?;
    let alleles = alleles_for(label);
    let mut samples = Vec::new();
    let mut evidence: Vec<Vec<BamRecord>> = Vec::new();
    let mut values: Vec<Vec<Vec<f64>>> = Vec::new();
    for (s, per_allele) in composition.iter().enumerate() {
        let mut reads = Vec::new();
        let mut best = Vec::new();
        for (a, counts) in per_allele.iter().enumerate() {
            for (strand, count) in counts.iter().enumerate() {
                for i in 0..*count {
                    reads.push(read(&format!("s{s}a{a}d{strand}i{i}"), strand == 1));
                    best.push(a);
                }
            }
        }
        let matrix: Vec<Vec<f64>> = (0..alleles.len())
            .map(|a| {
                best.iter()
                    .map(|b| if *b == a { -1.0 } else { -10.0 })
                    .collect()
            })
            .collect();
        samples.push(format!("s{s}"));
        evidence.push(reads);
        values.push(matrix);
    }
    Some(
        AlleleLikelihoods::new(
            SampleList::new(&samples),
            AlleleList::new(&alleles),
            evidence,
            values,
        )
        .expect("a matrix"),
    )
}

fn raw_strings(label: &str) -> Vec<String> {
    let raws: Vec<&str> = match label {
        "one-source" => vec!["10,8|3,4"],
        "two-sources" => vec!["10,8|3,4", "2,3|1,1"],
        "three-alleles" => vec!["10,8|3,4|2,2"],
        "empty-entry" => vec!["10,8|"],
        "zero-entry" => vec!["10,8|0,0"],
        "bracketed" => vec!["[10,8|3, 4]"],
        "spaced" => vec!["10, 8|3, 4"],
        "wrong-count" => vec!["10,8"],
        "extreme" => vec!["4000,1|1,4000"],
        other => panic!("{other} has no raw fixture"),
    };
    raws.into_iter().map(|s| s.to_string()).collect()
}

fn combine_alleles(label: &str) -> Vec<Allele> {
    if label == "three-alleles" {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    }
}

fn annotation_for(name: &str) -> AsStrandBias {
    match name {
        "AS_FisherStrand" => AsStrandBias::Fisher,
        "AS_StrandOddsRatio" => AsStrandBias::OddsRatio,
        other => panic!("unknown annotation {other}"),
    }
}

fn rendered(entries: &[(String, String)]) -> String {
    let mut sorted: Vec<&(String, String)> = entries.iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));
    sorted
        .iter()
        .map(|(key, value)| format!("{key}={value}[java.lang.String]"))
        .collect::<Vec<_>>()
        .join(";")
}

#[test]
fn every_raw_table_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("sbraw\t") else {
            continue;
        };
        let (label, expected) = match rest.split_once('\t') {
            Some(pair) => pair,
            None => (rest, ""),
        };
        let vc = variant_context(label);
        let ours = match annotate_raw_data(AsStrandBias::Fisher, &vc, likelihoods(label).as_ref()) {
            None => String::new(),
            Some((key, value)) => rendered(&[(key, value)]),
        };
        assert_eq!(ours, expected, "AS_SB_TABLE on {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no raw rows");
    println!("{count} allele-specific strand tables identical");
}

#[test]
fn every_direct_answer_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("sbdirect\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.splitn(3, '\t').collect();
        let name = fields[0];
        let label = fields[1];
        let expected = fields.get(2).copied().unwrap_or("");
        let vc = variant_context(label);
        let matrix = likelihoods(label);
        // The direct path is the pooled one, but with the allele-specific minimum count of two
        // rather than the plain StrandOddsRatio's zero.
        let ours = annotate_direct(annotation_for(name), &vc, matrix.as_ref())
            .iter()
            .map(|(key, value)| {
                let text = value.to_java_string().expect("a renderable value");
                format!("{key}={text}[java.lang.String]")
            })
            .collect::<Vec<_>>()
            .join(";");
        assert_eq!(ours, expected, "{name} on {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no direct rows");
    println!("{count} direct answers identical");
}

#[test]
fn every_combination_and_finalisation_matches_the_reference() {
    let text = golden();
    let mut combined_count = 0;
    let mut finalised = 0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("sbcombine\t") {
            let (label, expected) = match rest.split_once('\t') {
                Some(pair) => pair,
                None => (rest, ""),
            };
            let alleles = combine_alleles(label);
            let ours = match combine_raw_data(&alleles, &raw_strings(label)) {
                Ok(text) => rendered(&[(AS_SB_TABLE_KEY.to_string(), text)]),
                Err(AsStrandBiasError::AlleleCountMismatch { .. }) => {
                    "E:java.lang.IllegalStateException".to_string()
                }
                Err(other) => panic!("{label}: {other:?}"),
            };
            assert_eq!(ours, expected, "combination of {label}");
            combined_count += 1;
        } else if let Some(rest) = line.strip_prefix("sbfinal\t") {
            let fields: Vec<&str> = rest.splitn(3, '\t').collect();
            let name = fields[0];
            let label = fields[1];
            let expected = fields.get(2).copied().unwrap_or("");
            let alleles = combine_alleles(label);
            let combined = combine_raw_data(&alleles, &raw_strings(label)).expect("a combination");
            let ours =
                match finalize_raw_data(annotation_for(name), &alleles, &alleles, Some(&combined))
                    .expect("a finalisation")
                {
                    None => String::new(),
                    Some((key, reduced, raw_key, raw_value)) => {
                        rendered(&[(key, reduced), (raw_key, raw_value)])
                    }
                };
            assert_eq!(ours, expected, "{name} finalising {label}");
            finalised += 1;
        }
    }
    assert!(
        combined_count > 0 && finalised > 0,
        "the golden carries no combination rows"
    );
    println!("{combined_count} combinations and {finalised} finalisations identical");
}

#[test]
fn every_phred_scale_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("phred\t") else {
            continue;
        };
        let (rate, expected) = rest.split_once('\t').expect("a rate and a result");
        let rate: f64 = rate.parse().expect("a double");
        let want = expected.parse::<i64>().expect("bits") as u64;
        let ours = phred_scale_error_rate(rate);
        assert_eq!(
            ours.to_bits(),
            want,
            "phredScaleErrorRate({rate:e}) = {ours}, reference {}",
            f64::from_bits(want)
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no phred rows");
    println!("{count} phred scalings bit-identical");
}
