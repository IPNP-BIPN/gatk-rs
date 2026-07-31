//! Conformance for `AD`, `AF`, `DP` and the marginalisation under them, against the oracle.
//!
//! Golden from `tools/annotation-conformance/DepthPerAlleleDump.java`.
//!
//! The rows that carry the claim are the order rows, and they say the allele order is **not** the
//! variant's:
//!
//! ```text
//! order  A*,C        C,A*
//! order  A*,C,G,T    T,C,A*,G
//! ```
//!
//! `Collectors.toMap` builds a `HashMap` and `marginalize` takes its key set, so the new matrix's
//! allele order follows `Allele.hashCode`. `searchBestAllele` breaks a tie by keeping the first
//! index, so that order decides which allele a tied read is counted for.

use std::io::Read;

use gatk_annotation::depth_per_allele::{
    allele_depths, allele_fractions, informative_depth, marginalisation_order,
};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

const START: i64 = 105;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/depth_per_allele.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

/// The dump renders an allele as its bases with a `*` for the reference.
fn shown(allele: &Allele) -> String {
    format!(
        "{}{}",
        allele.display_string(),
        if allele.is_reference() { "*" } else { "" }
    )
}

fn parse_allele(text: &str) -> Allele {
    match text.strip_suffix('*') {
        Some(bases) => Allele::from_str(bases, true).expect("an allele"),
        None => Allele::from_str(text, false).expect("an allele"),
    }
}

fn read(name: &str, quality: u8) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: 100,
        mapping_quality: 60,
        cigar: htsjdk_bam::text_parse::parse_cigar("20M").expect("a cigar"),
        read_bases: vec![b'A'; 20],
        base_qualities: vec![quality; 20],
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

fn biallelic(
    ref_reads: usize,
    alt_reads: usize,
    uninformative: bool,
) -> AlleleLikelihoods<BamRecord> {
    let mut reads = Vec::new();
    for i in 0..ref_reads {
        reads.push(read(&format!("r{i}"), 30));
    }
    for i in 0..alt_reads {
        reads.push(read(&format!("a{i}"), 30));
    }
    let total = reads.len();
    let values = (0..2)
        .map(|a| {
            (0..total)
                .map(|e| {
                    let is_ref = e < ref_reads;
                    if uninformative {
                        if a == 0 {
                            -1.0
                        } else {
                            -1.1
                        }
                    } else if (a == 0) == is_ref {
                        -1.0
                    } else {
                        -10.0
                    }
                })
                .collect()
        })
        .collect();
    matrix(
        reads,
        vec![
            Allele::from_str("A", true).unwrap(),
            Allele::from_str("C", false).unwrap(),
        ],
        values,
    )
}

fn case(label: &str) -> (VariantContext, AlleleLikelihoods<BamRecord>) {
    let biallelic_site = || {
        let mut vc = VariantContext::new(
            "chr1",
            START,
            vec![
                Allele::from_str("A", true).unwrap(),
                Allele::from_str("C", false).unwrap(),
            ],
        );
        vc.stop = START;
        vc
    };
    match label {
        "two-and-two" => (biallelic_site(), biallelic(2, 2, false)),
        "all-ref" => (biallelic_site(), biallelic(4, 0, false)),
        "all-alt" => (biallelic_site(), biallelic(0, 4, false)),
        "empty" => (biallelic_site(), biallelic(0, 0, false)),
        "uninformative" => (biallelic_site(), biallelic(2, 2, true)),
        "triallelic" => {
            let reads = vec![
                read("r0", 30),
                read("a0", 30),
                read("b0", 30),
                read("t0", 30),
            ];
            let alleles = vec![
                Allele::from_str("A", true).unwrap(),
                Allele::from_str("C", false).unwrap(),
                Allele::from_str("G", false).unwrap(),
            ];
            let values = (0..3)
                .map(|a| {
                    (0..4)
                        .map(|e| {
                            if e == 3 {
                                if a == 1 || a == 2 {
                                    -1.0
                                } else {
                                    -10.0
                                }
                            } else if a == e {
                                -1.0
                            } else {
                                -10.0
                            }
                        })
                        .collect()
                })
                .collect();
            let mut vc = VariantContext::new("chr1", START, alleles.clone());
            vc.stop = START;
            (vc, matrix(reads, alleles, values))
        }
        "matrix-has-extra-allele" => {
            let reads = vec![read("r0", 30), read("a0", 30), read("x0", 30)];
            let alleles = vec![
                Allele::from_str("A", true).unwrap(),
                Allele::from_str("C", false).unwrap(),
                Allele::from_str("T", false).unwrap(),
            ];
            let values = (0..3)
                .map(|a| (0..3).map(|e| if a == e { -1.0 } else { -10.0 }).collect())
                .collect();
            (biallelic_site(), matrix(reads, alleles, values))
        }
        other => panic!("{other} has no fixture"),
    }
}

#[test]
fn the_marginalisation_order_is_the_reference_hash_order() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("order\t") else {
            continue;
        };
        let (input, expected) = rest.split_once('\t').expect("in and out");
        let alleles: Vec<Allele> = input.split(',').map(parse_allele).collect();
        let ours: Vec<String> = marginalisation_order(&alleles).iter().map(shown).collect();
        assert_eq!(ours.join(","), expected, "order of {input}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no order rows");
    println!("{count} marginalisation orders identical");
}

#[test]
fn every_depth_and_fraction_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    let mut nan_sign_exemptions = 0;
    for line in text.lines() {
        let (kind, rest) = if let Some(rest) = line.strip_prefix("ad\t") {
            ("ad", rest)
        } else if let Some(rest) = line.strip_prefix("af\t") {
            ("af", rest)
        } else if let Some(rest) = line.strip_prefix("dp\t") {
            ("dp", rest)
        } else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("a label and a value");
        let (vc, likelihoods) = case(label);
        let ours = match kind {
            "ad" => allele_depths(&vc, Some(&likelihoods), "s1", true)
                .map(|counts| {
                    counts
                        .iter()
                        .map(|c| c.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default(),
            "af" => allele_fractions(&vc, None, Some(&likelihoods), "s1", true)
                .map(|fractions| {
                    fractions
                        .iter()
                        .map(|f| (f.to_bits() as i64).to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default(),
            _ => informative_depth(Some(&likelihoods), "s1", true)
                .map(|d| d.to_string())
                .unwrap_or_default(),
        };
        if ours != expected && kind == "af" {
            // Decision 0012: an all-zero AD divides zero by zero, and the NaN's sign is the FPU's.
            let parse = |text: &str| {
                text.split(',')
                    .filter_map(|part| part.parse::<i64>().ok())
                    .map(|bits| f64::from_bits(bits as u64))
                    .collect::<Vec<f64>>()
            };
            let (a, b) = (parse(&ours), parse(expected));
            if a.len() == b.len()
                && a.iter().zip(&b).all(|(x, y)| {
                    x.to_bits() == y.to_bits()
                        || (x.is_nan() && y.is_nan() && !cfg!(target_arch = "x86_64"))
                })
            {
                nan_sign_exemptions += 1;
                count += 1;
                continue;
            }
        }
        assert_eq!(ours, expected, "{kind} on {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no rows");
    if cfg!(target_arch = "x86_64") {
        assert_eq!(nan_sign_exemptions, 0, "nothing to exempt on x86-64");
    } else {
        assert_eq!(
            nan_sign_exemptions, 2,
            "the NaN-sign exemption count changed"
        );
    }
    println!("{count} depths and fractions identical, {nan_sign_exemptions} NaN-sign exemptions");
}

#[test]
fn every_marginalised_likelihood_is_bit_identical() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("marginal\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let label = fields.next().expect("a label");
        let _sample = fields.next().expect("a sample");
        let allele = parse_allele(fields.next().expect("an allele"));
        let expected = fields.next().unwrap_or("");
        let (vc, likelihoods) = case(label);
        let order = marginalisation_order(&vc.alleles);
        let new_to_old: Vec<(Allele, Vec<Allele>)> =
            order.iter().map(|a| (a.clone(), vec![a.clone()])).collect();
        let marginal = likelihoods.marginalize(&new_to_old).expect("a matrix");
        let index = marginal.index_of_allele(&allele).expect("the allele");
        let values: Vec<String> = (0..marginal.sample_evidence_count(0))
            .map(|e| (marginal.value(0, index, e).to_bits() as i64).to_string())
            .collect();
        assert_eq!(values.join(","), expected, "{label} allele {allele:?}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no marginal rows");
    println!("{count} marginalised rows bit-identical");
}
