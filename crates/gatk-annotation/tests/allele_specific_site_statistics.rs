//! Conformance for `AS_QD`, `AS_MQ`, `AS_InbreedingCoeff` and the heterozygosity calculator,
//! against the oracle.
//!
//! Golden from `tools/annotation-conformance/AlleleSpecificSiteStatisticsDump.java`.
//!
//! ```text
//! asqd       one-alt-read  AS_QD=25.36[java.lang.String]     <- a random draw, refused here
//! asmq       no-ad         AS_MQ=Infinity[java.lang.String]
//! asmqfinal  3600.00|      AS_MQ=.[...];AS_RAW_MQ=3600.00|0.00[...]
//! ```
//!
//! Three `AS_QD` rows are values the reference randomised, for the same reason `QD` does: the
//! replacement at thirty-five goes through `StrictMath.log`, which is fdlibm. The suite asserts
//! that the port refuses exactly those rows, that the raw ratio really was at or past the
//! threshold, and that what the reference wrote is not that ratio.
//!
//! The heterozygosity counts are compared as raw bits, because they are sums of normalised
//! likelihoods and a decimal rendering would hide the last ulp.

use std::io::Read;

use gatk_annotation::allele_specific_site_statistics::{
    ad_counts, allele_depths, as_inbreeding_coefficient, as_qual_by_depth, as_rms_data,
    as_rms_finalized_string, as_rms_parse_raw, as_rms_raw_string, encode_value_list,
    AsSiteStatisticError, AS_QUAL_BY_DEPTH_KEY, AS_RAW_RMS_MAPPING_QUALITY_KEY,
    AS_RMS_MAPPING_QUALITY_KEY,
};
use gatk_annotation::heterozygosity::heterozygosity_counts;
use gatk_annotation::site_statistics::QualByDepthError;
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

const START: i64 = 105;
const MAX_QD_BEFORE_FIXING: f64 = 35.0;

const HOM_REF: [i32; 3] = [0, 60, 600];
const HET: [i32; 3] = [60, 0, 60];
const HOM_VAR: [i32; 3] = [600, 60, 0];

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/allele_specific_site_statistics.txt.gz");
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

// ---- AS_QualByDepth ---------------------------------------------------------------------------

fn qd_triallelic(label: &str) -> bool {
    matches!(label, "triallelic" | "mixed-qd" | "empty-slot")
}

/// The dump's `AS_QUAL` attribute, as `getAttributeAsList` would see it: one element for a String.
fn qd_as_qual(label: &str) -> Option<Vec<String>> {
    match label {
        "as-qual-approx" | "empty-slot" | "no-qual" => None,
        "as-vardp" => Some(vec!["300".to_string()]),
        "high-qd" => Some(vec!["5000".to_string()]),
        // Set programmatically as one String, so the list has **one** element and the triallelic
        // count check fails.
        "mixed-qd" => Some(vec!["70,5000".to_string()]),
        "triallelic" => Some(vec!["70,140".to_string()]),
        _ => Some(vec!["300".to_string()]),
    }
}

fn qd_as_qual_approx(label: &str) -> Option<&'static str> {
    match label {
        "as-qual-approx" => Some("0|300"),
        "empty-slot" => Some("0|300|"),
        _ => None,
    }
}

fn qd_as_vardp(label: &str) -> Option<&'static str> {
    if label == "as-vardp" {
        Some("10|8")
    } else {
        None
    }
}

fn qd_context(label: &str) -> VariantContext {
    let triallelic = qd_triallelic(label);
    let alleles = if triallelic {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    };
    let mut vc = VariantContext::new("chr1", START, alleles);
    vc.stop = START;
    if label == "no-genotypes" {
        return vc;
    }
    let mut genotype = match label {
        "hom-ref-only" => Genotype::new("s0", vec![reference(), reference()]),
        _ => Genotype::new("s0", vec![reference(), alternate()]),
    };
    match label {
        "hom-ref-only" => {
            genotype.ad = Some(if triallelic {
                vec![10, 0, 0]
            } else {
                vec![10, 0]
            })
        }
        "no-ad" => genotype.dp = Some(17),
        "one-alt-read" => genotype.ad = Some(vec![10, 1]),
        _ => {
            genotype.ad = Some(if triallelic {
                vec![10, 4, 6]
            } else {
                vec![10, 8]
            })
        }
    }
    vc.genotypes.push(genotype);
    vc
}

// ---- AS_RMSMappingQuality ---------------------------------------------------------------------

fn mq_composition(label: &str) -> Option<Vec<Vec<u8>>> {
    let table: Vec<Vec<u8>> = match label {
        "biallelic" => vec![vec![60, 60, 60], vec![30, 30]],
        "triallelic" => vec![vec![60, 60], vec![30, 30], vec![20, 20]],
        "ref-only" => vec![vec![60, 60, 60], vec![]],
        "alt-only" => vec![vec![], vec![30, 30]],
        "all-unavailable" => vec![vec![255, 255], vec![255, 255]],
        "mixed-unavailable" => vec![vec![60, 255], vec![30, 255]],
        "no-ad" => vec![vec![60, 60], vec![30, 30]],
        "null-likelihoods" => return None,
        other => panic!("{other} has no composition"),
    };
    Some(table)
}

fn mq_context(label: &str) -> VariantContext {
    let triallelic = label == "triallelic";
    let alleles = if triallelic {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    };
    let mut vc = VariantContext::new("chr1", START, alleles);
    vc.stop = START;
    let mut genotype = Genotype::new("s1", vec![reference(), alternate()]);
    if label != "no-ad" {
        genotype.ad = Some(if triallelic {
            vec![2, 2, 2]
        } else {
            vec![3, 2]
        });
    }
    vc.genotypes.push(genotype);
    vc
}

fn read(name: &str, mapping_quality: u8) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: 100,
        mapping_quality,
        cigar: htsjdk_bam::text_parse::parse_cigar("20M").expect("a cigar"),
        read_bases: vec![b'A'; 20],
        base_qualities: vec![30; 20],
        ..Default::default()
    }
}

fn mq_likelihoods(label: &str) -> Option<AlleleLikelihoods<BamRecord>> {
    let composition = mq_composition(label)?;
    let alleles = if label == "triallelic" {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    };
    let mut reads = Vec::new();
    let mut best = Vec::new();
    for (a, qualities) in composition.iter().enumerate() {
        for (i, mq) in qualities.iter().enumerate() {
            reads.push(read(&format!("a{a}i{i}"), *mq));
            best.push(a);
        }
    }
    let values: Vec<Vec<f64>> = (0..alleles.len())
        .map(|a| {
            best.iter()
                .map(|b| if *b == a { -1.0 } else { -10.0 })
                .collect()
        })
        .collect();
    Some(
        AlleleLikelihoods::new(
            SampleList::new(&["s1".to_string()]),
            AlleleList::new(&alleles),
            vec![reads],
            vec![values],
        )
        .expect("a matrix"),
    )
}

/// The dump's `asmqFinal` context, whose allele count follows the raw string's slot count.
fn mq_final_context(raw: &str) -> VariantContext {
    let slots = raw.split('|').count();
    let alleles = if slots > 2 {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    };
    let mut vc = VariantContext::new("chr1", START, alleles);
    vc.stop = START;
    let mut genotype = Genotype::new("s1", vec![reference(), alternate()]);
    genotype.ad = Some(if slots > 2 { vec![2, 2, 2] } else { vec![3, 2] });
    vc.genotypes.push(genotype);
    vc
}

// ---- AS_InbreedingCoeff -----------------------------------------------------------------------

fn ic_context(label: &str) -> VariantContext {
    let triallelic = label == "triallelic";
    let alleles = if triallelic {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    };
    let mut vc = VariantContext::new("chr1", START, alleles);
    vc.stop = START;
    let mut push = |name: &str, called: Vec<Allele>, pl: Option<Vec<i32>>, gq: Option<i32>| {
        let mut genotype = Genotype::new(name, called);
        genotype.pl = pl;
        genotype.gq = gq;
        vc.genotypes.push(genotype);
    };
    match label {
        "ten-het" => {
            for i in 0..10 {
                push(
                    &format!("s{i}"),
                    vec![reference(), alternate()],
                    Some(HET.to_vec()),
                    None,
                );
            }
        }
        "ten-hom-ref" => {
            for i in 0..10 {
                push(
                    &format!("s{i}"),
                    vec![reference(), reference()],
                    Some(HOM_REF.to_vec()),
                    None,
                );
            }
        }
        "nine-samples" => {
            for i in 0..9 {
                push(
                    &format!("s{i}"),
                    vec![reference(), alternate()],
                    Some(HET.to_vec()),
                    None,
                );
            }
        }
        "twenty-mixed" => {
            for i in 0..5 {
                push(
                    &format!("r{i}"),
                    vec![reference(), reference()],
                    Some(HOM_REF.to_vec()),
                    None,
                );
            }
            for i in 0..10 {
                push(
                    &format!("h{i}"),
                    vec![reference(), alternate()],
                    Some(HET.to_vec()),
                    None,
                );
            }
            for i in 0..5 {
                push(
                    &format!("v{i}"),
                    vec![alternate(), alternate()],
                    Some(HOM_VAR.to_vec()),
                    None,
                );
            }
        }
        "triallelic" => {
            for i in 0..12 {
                push(
                    &format!("s{i}"),
                    vec![reference(), alternate()],
                    Some(vec![60, 0, 60, 60, 60, 600]),
                    None,
                );
            }
        }
        "gq-only" => {
            for i in 0..12 {
                push(
                    &format!("s{i}"),
                    vec![reference(), alternate()],
                    None,
                    Some(30),
                );
            }
        }
        "no-genotypes" => {}
        other => panic!("{other} has no cohort fixture"),
    }
    vc
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
fn every_quality_by_depth_matches_or_is_refused_as_randomised() {
    let text = golden();
    let mut identical = 0;
    let mut refused = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("asqd\t") else {
            continue;
        };
        let (label, expected) = match rest.split_once('\t') {
            Some(pair) => pair,
            None => (rest, ""),
        };
        let vc = qd_context(label);
        let quals = qd_as_qual(label);
        let ours = as_qual_by_depth(
            &vc,
            quals.as_deref(),
            qd_as_qual_approx(label),
            qd_as_vardp(label),
        );
        match ours {
            Err(AsSiteStatisticError::QualCountMismatch { .. })
            | Err(AsSiteStatisticError::QualApproxCountMismatch { .. }) => {
                assert_eq!(expected, "E:java.lang.IllegalStateException", "{label}");
                identical += 1;
            }
            Err(other) => panic!("{label}: {other:?}"),
            Ok(None) => {
                assert_eq!(expected, "", "{label}");
                identical += 1;
            }
            Ok(Some(entries)) => {
                if entries.iter().any(|entry| entry.is_err()) {
                    // Every refused entry must be one the reference randomised, and what it wrote
                    // must not be the raw ratio.
                    let written = expected
                        .strip_prefix("AS_QD=")
                        .and_then(|rest| rest.strip_suffix("[java.lang.String]"))
                        .expect("a written value");
                    for entry in &entries {
                        let Err(QualByDepthError::RandomisedAboveThreshold { raw }) = entry else {
                            continue;
                        };
                        assert!(
                            *raw >= MAX_QD_BEFORE_FIXING,
                            "{label} refused a ratio of {raw}, which is below the threshold"
                        );
                        let as_written = format!("{raw:.2}");
                        assert_ne!(
                            written, as_written,
                            "{label}: the reference wrote the raw ratio, so nothing was randomised"
                        );
                    }
                    refused += 1;
                } else {
                    let values: Vec<f64> = entries
                        .iter()
                        .map(|entry| *entry.as_ref().unwrap())
                        .collect();
                    let ours = rendered(&[(
                        AS_QUAL_BY_DEPTH_KEY.to_string(),
                        encode_value_list(&values, 2),
                    )]);
                    assert_eq!(ours, expected, "AS_QD on {label}");
                    identical += 1;
                }
            }
        }
    }
    assert!(identical > 0, "the golden carries no AS_QD rows");
    assert!(refused > 0, "the golden carries no randomised AS_QD rows");
    println!("{identical} AS_QD answers identical, {refused} refused as randomised");
}

#[test]
fn every_allele_depth_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("asdepths\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("a label and a result");
        let vc = qd_context(label);
        let ours = match allele_depths(&vc) {
            None => "null".to_string(),
            Some(depths) => depths
                .iter()
                .map(|d| d.to_string())
                .collect::<Vec<_>>()
                .join(","),
        };
        assert_eq!(ours, expected, "depths of {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no depth rows");
    println!("{count} allele-depth vectors identical");
}

#[test]
fn every_mapping_quality_answer_matches_the_reference() {
    let text = golden();
    let mut direct = 0;
    let mut raw = 0;
    let mut finalised = 0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("asmqraw\t") {
            let (label, expected) = match rest.split_once('\t') {
                Some(pair) => pair,
                None => (rest, ""),
            };
            let vc = mq_context(label);
            let ours = match mq_likelihoods(label) {
                None => String::new(),
                Some(matrix) => {
                    let data = as_rms_data(&matrix);
                    rendered(&[(
                        AS_RAW_RMS_MAPPING_QUALITY_KEY.to_string(),
                        as_rms_raw_string(&vc.alleles, &data),
                    )])
                }
            };
            assert_eq!(ours, expected, "AS_RAW_MQ on {label}");
            raw += 1;
        } else if let Some(rest) = line.strip_prefix("asmqfinal\t") {
            let (raw_string, expected) = rest.split_once('\t').expect("a raw string and a result");
            let vc = mq_final_context(raw_string);
            let parsed = as_rms_parse_raw(&vc.alleles, raw_string).expect("a parse");
            let finalised_string =
                as_rms_finalized_string(&vc, &parsed).expect("genotypes to divide by");
            let ours = rendered(&[
                (AS_RMS_MAPPING_QUALITY_KEY.to_string(), finalised_string),
                (
                    AS_RAW_RMS_MAPPING_QUALITY_KEY.to_string(),
                    as_rms_raw_string(&vc.alleles, &parsed),
                ),
            ]);
            assert_eq!(ours, expected, "finalising {raw_string}");
            finalised += 1;
        } else if let Some(rest) = line.strip_prefix("asmq\t") {
            let (label, expected) = match rest.split_once('\t') {
                Some(pair) => pair,
                None => (rest, ""),
            };
            let vc = mq_context(label);
            let ours = match mq_likelihoods(label) {
                None => String::new(),
                Some(matrix) => {
                    let data = as_rms_data(&matrix);
                    let text = as_rms_finalized_string(&vc, &data).expect("genotypes to divide by");
                    rendered(&[(AS_RMS_MAPPING_QUALITY_KEY.to_string(), text)])
                }
            };
            assert_eq!(ours, expected, "AS_MQ on {label}");
            direct += 1;
        }
    }
    assert!(
        direct > 0 && raw > 0 && finalised > 0,
        "the golden carries no mapping-quality rows"
    );
    println!("{direct} AS_MQ, {raw} AS_RAW_MQ and {finalised} finalisations identical");
}

#[test]
fn every_heterozygosity_count_and_coefficient_matches_the_reference() {
    let text = golden();
    let mut coefficients = 0;
    let mut counts = 0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("asic\t") {
            let (label, expected) = match rest.split_once('\t') {
                Some(pair) => pair,
                None => (rest, ""),
            };
            let vc = ic_context(label);
            let ours = match as_inbreeding_coefficient(&vc) {
                None => String::new(),
                Some(value) => format!("AS_InbreedingCoeff={value}[java.lang.String]"),
            };
            assert_eq!(ours, expected, "AS_InbreedingCoeff on {label}");
            coefficients += 1;
        } else if let Some(rest) = line.strip_prefix("hetcounts\t") {
            let fields: Vec<&str> = rest.split('\t').collect();
            let label = fields[0];
            let vc = ic_context(label);
            let ours = heterozygosity_counts(&vc);
            assert_eq!(
                ours.sample_count.to_string(),
                fields[1],
                "sample count of {label}"
            );
            for entry in fields[2].split(';') {
                let (name, bits) = entry.split_once('=').expect("an allele and its bits");
                let want = bits.parse::<i64>().expect("bits") as u64;
                let got = ours.het_count(&allele(name, false));
                assert_eq!(got.to_bits(), want, "het count of {name} in {label}");
            }
            for entry in fields[3].split(';') {
                let (name, bits) = entry.split_once('=').expect("an allele and its bits");
                let want = bits.parse::<i64>().expect("bits") as u64;
                let is_reference = name == "A";
                let got = ours.allele_count(&allele(name, is_reference));
                assert_eq!(got.to_bits(), want, "allele count of {name} in {label}");
            }
            counts += 1;
        }
    }
    assert!(
        coefficients > 0 && counts > 0,
        "the golden carries no heterozygosity rows"
    );
    println!("{coefficients} coefficients and {counts} count tables identical");
}

#[test]
fn the_two_depth_definitions_disagree_where_the_reference_says_they_do() {
    // `AS_QualByDepth.getAlleleDepths` drops a genotype whose alternate depth is one; the
    // `AS_RMSMappingQuality.getADcounts` of the same site keeps it.
    let vc = qd_context("one-alt-read");
    assert_eq!(allele_depths(&vc), Some(vec![0, 0]));
    let counts = ad_counts(&vc).expect("genotypes");
    assert_eq!(counts[1].1, 1);
}
