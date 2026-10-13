//! Conformance for the DRAGEN genotyping models (BQD and FRD) against the oracle.
//!
//! Golden from `tools/readfilter-conformance/DragenGenotypesDump.java`. Every row carries its own
//! inputs or follows rows that do (`reference`, then `case`, `read`, `lk` for a case's model runs);
//! this test recomputes each row in turn and requires the golden's line exactly.

use std::collections::HashMap;

use gatk_corpus as corpus;
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use gatk_engine::dragen_genotypes::{
    indel_prior, map_mapping_quality_to_phred_likelihood_score, transform_mapping_quality,
    DragenGenotypesModel, DragenSample,
};
use gatk_engine::dragstr::{DragstrParams, DragstrReferenceAnalyzer};
use gatk_engine::genotype_likelihood_calculator::gls_to_pls;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_vcf::allele::Allele;

const REVERSE: u16 = 0x10;

struct Case {
    alleles: Vec<Allele>,
    ploidy: usize,
    variant: i32,
    padding: i32,
    offset: i32,
    str_on: bool,
    max_effective_depth: i32,
    samples: Vec<String>,
    reads: Vec<Vec<BamRecord>>,
    /// `[sample][allele][read]`.
    lk: Vec<Vec<Vec<f64>>>,
}

fn hex(value: f64) -> String {
    format!("{:x}", value.to_bits())
}

fn parse_bits(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).expect("hex bits"))
}

fn record(name: &str, start: i32, cigar: &str, reverse: bool, mq: u8, quals: &[u8]) -> BamRecord {
    let cigar = htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar");
    BamRecord {
        read_name: name.to_string(),
        flags: if reverse { REVERSE } else { 0 },
        reference_index: 0,
        alignment_start: start,
        mapping_quality: mq,
        read_bases: vec![b'A'; quals.len()],
        base_qualities: quals.to_vec(),
        cigar,
        ..BamRecord::default()
    }
}

/// The model over one case and configuration, rendered as the dump's `filtered` and `gl` rows.
fn run_case(
    label: &str,
    case: &Case,
    config: &str,
    reference: &[u8],
) -> (Vec<String>, Vec<String>) {
    let mut filtered_rows = Vec::new();
    let mut gl_rows = Vec::new();
    let samples = SampleList::new(&case.samples);
    let alleles = AlleleList::new(&case.alleles);
    let evidence: Vec<Vec<usize>> = case.reads.iter().map(|r| (0..r.len()).collect()).collect();
    let mut likelihoods =
        AlleleLikelihoods::new(samples, alleles, evidence, case.lk.clone()).expect("likelihoods");
    likelihoods
        .filter_poorly_modeled_evidence(|_, _, _| -20.0)
        .expect("filtered");
    let params = DragstrParams::default_params();
    let dragstrs = case.str_on.then(|| {
        DragstrReferenceAnalyzer::of(
            reference,
            0,
            reference.len() as i32,
            params.maximum_period(),
        )
        .expect("an analyzer")
    });
    let model = DragenGenotypesModel::new(
        config != "FRD",
        config != "BQD",
        case.padding,
        case.max_effective_depth,
    )
    .expect("a model");
    let api = indel_prior(dragstrs.as_ref(), &params, case.offset).expect("a prior");
    for s in 0..case.samples.len() {
        let kept: Vec<BamRecord> = likelihoods
            .sample_evidence(s)
            .unwrap()
            .iter()
            .map(|&i| case.reads[s][i].clone())
            .collect();
        let filtered: Vec<BamRecord> = likelihoods
            .filtered_evidence(s)
            .iter()
            .map(|&i| case.reads[s][i].clone())
            .collect();
        filtered_rows.push(format!(
            "filtered\t{label}\t{s}\t{}",
            filtered
                .iter()
                .map(|r| r.read_name.as_str())
                .collect::<Vec<_>>()
                .join(",")
        ));
        let matrix: Vec<Vec<f64>> = (0..case.alleles.len())
            .map(|a| {
                (0..kept.len())
                    .map(|e| likelihoods.value(s, a, e))
                    .collect()
            })
            .collect();
        let sample = DragenSample {
            likelihoods: &matrix,
            reads: &kept,
            filtered: &filtered,
        };
        match model.sample_log10_likelihoods(
            case.ploidy,
            &case.alleles,
            &sample,
            case.variant - case.padding,
            reference,
            case.offset,
            api,
        ) {
            Ok(values) => gl_rows.push(format!(
                "gl\t{label}\t{config}\t{s}\t{}\t{}",
                values.iter().map(|v| hex(*v)).collect::<Vec<_>>().join(","),
                gls_to_pls(&values)
                    .iter()
                    .map(i32::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )),
            Err(message) => gl_rows.push(format!(
                "error\t{label}\t{config}\tGATKException: {message}"
            )),
        }
    }
    (filtered_rows, gl_rows)
}

/// `bad-<start>-<end>`, either of which may be negative.
fn bad_bounds(label: &str) -> (i32, i32) {
    let rest = &label["bad-".len()..];
    let split = rest[1..].find('-').expect("two bounds") + 1;
    (
        rest[..split].parse().unwrap(),
        rest[split + 1..].parse().unwrap(),
    )
}

#[test]
fn every_row_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/dragen_genotypes.txt.gz"),
    );
    let mut reference: Vec<u8> = Vec::new();
    let mut case: Option<(String, Case)> = None;
    let mut computed: HashMap<String, (Vec<String>, Vec<String>)> = HashMap::new();
    let mut compared = 0;
    let mut sample_index: HashMap<String, usize> = HashMap::new();
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        let got = match f[0] {
            "reference" => {
                reference = f[1].as_bytes().to_vec();
                continue;
            }
            "mq" => {
                let mq: i32 = f[1].parse().unwrap();
                match map_mapping_quality_to_phred_likelihood_score(mq) {
                    Ok(v) => format!("mq\t{mq}\t{}", hex(v)),
                    Err(m) => format!("mq\t{mq}\terror\tGATKException: {m}"),
                }
            }
            "xq" => {
                let (mq, xq) = match f[1] {
                    "tagged" => (60, Some(187)),
                    "plain" => (42, None),
                    _ => (60, Some(0)),
                };
                let mut r = record(f[1], 1000, "10M", false, mq, &[30; 10]);
                if let Some(xq) = xq {
                    r.tags.insert(Tag::new(b"XQ"), TagValue::Int(xq));
                }
                transform_mapping_quality(&mut r).unwrap();
                format!("xq\t{}\t{}", f[1], r.mapping_quality)
            }
            "str" => {
                let (start, end, max_period, position): (i32, i32, usize, i32) = (
                    f[3].parse().unwrap(),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                );
                let a =
                    DragstrReferenceAnalyzer::of(f[2].as_bytes(), start, end, max_period).unwrap();
                format!(
                    "str\t{}\t{}\t{start}\t{end}\t{max_period}\t{position}\t{}\t{}",
                    f[1],
                    f[2],
                    a.period(position).unwrap(),
                    a.repeat_length(position).unwrap()
                )
            }
            "strerror" => {
                let (start, end) = bad_bounds(f[1]);
                match DragstrReferenceAnalyzer::of(&reference, start, end, 8) {
                    Ok(_) => format!("strerror\t{}\tnone", f[1]),
                    Err(m) => format!("strerror\t{}\tIndexOutOfBoundsException: {m}", f[1]),
                }
            }
            "case" => {
                let samples: Vec<String> = f[9].split(',').map(str::to_string).collect();
                sample_index = samples
                    .iter()
                    .enumerate()
                    .map(|(i, s)| (s.clone(), i))
                    .collect();
                let alleles: Vec<Allele> = f[2]
                    .split(',')
                    .enumerate()
                    .map(|(i, a)| Allele::create(a.as_bytes(), i == 0).unwrap())
                    .collect();
                let n = samples.len();
                let a = alleles.len();
                case = Some((
                    f[1].to_string(),
                    Case {
                        alleles,
                        ploidy: f[3].parse().unwrap(),
                        variant: f[4].parse().unwrap(),
                        padding: f[5].parse().unwrap(),
                        offset: f[6].parse().unwrap(),
                        str_on: f[7] == "true",
                        max_effective_depth: f[8].parse().unwrap(),
                        samples,
                        reads: vec![Vec::new(); n],
                        lk: vec![vec![Vec::new(); a]; n],
                    },
                ));
                computed.clear();
                continue;
            }
            "read" => {
                let (_, c) = case.as_mut().unwrap();
                let quals: Vec<u8> = f[8].split(',').map(|q| q.parse().unwrap()).collect();
                let r = record(
                    f[3],
                    f[4].parse().unwrap(),
                    f[5],
                    f[6] == "true",
                    f[7].parse().unwrap(),
                    &quals,
                );
                c.reads[sample_index[f[2]]].push(r);
                continue;
            }
            "lk" => {
                let (_, c) = case.as_mut().unwrap();
                let s: usize = f[2].parse().unwrap();
                let a: usize = f[3].parse().unwrap();
                c.lk[s][a] = if f[4].is_empty() {
                    Vec::new()
                } else {
                    f[4].split(',').map(parse_bits).collect()
                };
                continue;
            }
            "filtered" => {
                let (label, c) = case.as_ref().unwrap();
                let entry = computed
                    .entry("BQD".to_string())
                    .or_insert_with(|| run_case(label, c, "BQD", &reference));
                let s: usize = f[2].parse().unwrap();
                entry.0[s].clone()
            }
            "gl" | "error" => {
                let (label, c) = case.as_ref().unwrap();
                let config = f[2].to_string();
                let entry = computed
                    .entry(config.clone())
                    .or_insert_with(|| run_case(label, c, &config, &reference));
                if f[0] == "error" {
                    entry
                        .1
                        .iter()
                        .find(|r| r.starts_with("error"))
                        .cloned()
                        .unwrap_or_default()
                } else {
                    let s: usize = f[3].parse().unwrap();
                    entry.1[s].clone()
                }
            }
            other => panic!("unknown row kind {other}"),
        };
        assert_eq!(got, line, "row {compared}");
        compared += 1;
    }
    assert!(compared > 500, "only {compared} rows compared");
}
