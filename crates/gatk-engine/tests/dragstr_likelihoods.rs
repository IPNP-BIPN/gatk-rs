//! Conformance for DRAGstr in the read likelihoods against the oracle.
//!
//! Golden from `tools/readfilter-conformance/DragstrDump.java`. Every row carries its own inputs or
//! follows rows that do (`paramtext` for the parsed tables, `haplotype` and `read` for the engine);
//! this test recomputes each row in turn and requires the golden's line exactly.

use std::collections::HashMap;

use gatk_corpus as corpus;
use gatk_engine::dragstr::{
    dragstr_impute, non_symmetrical_impute, parse, DragstrError, DragstrParams,
    DragstrReadStrAnalyzer,
};
use gatk_engine::haplotype::Haplotype;
use gatk_engine::pair_hmm_likelihood_engine::{
    LikelihoodEngineArguments, PairHmmLikelihoodEngine, PcrErrorModel,
};
use htsjdk_bam::record::BamRecord;

fn csv(values: &[u8]) -> String {
    values
        .iter()
        .map(u8::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn render_error(label: &str, e: &DragstrError) -> String {
    let message = match e {
        DragstrError::BadInput(m) => format!("Bad input: {m}"),
        DragstrError::IllegalArgument(m) | DragstrError::IllegalState(m) => m.clone(),
        DragstrError::NullPointer => "null".to_string(),
        other => format!("{other:?}"),
    };
    format!("error\t{label}\t{}: {message}", e.class())
}

fn param_row(
    label: &str,
    params: &DragstrParams,
    table: &str,
    period: i32,
    repeats: i32,
) -> String {
    let value = match table {
        "gop" => params.gop(period, repeats),
        "gcp" => params.gcp(period, repeats),
        _ => params.api(period, repeats),
    }
    .unwrap();
    format!(
        "param\t{label}\t{table}\t{period}\t{repeats}\t{:x}",
        value.to_bits()
    )
}

#[test]
fn every_row_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/dragstr_likelihoods.txt.gz"),
    );
    let default = DragstrParams::default_params();
    let mut parsed: HashMap<String, Result<DragstrParams, DragstrError>> = HashMap::new();
    let mut haplotypes: Vec<Haplotype> = Vec::new();
    let mut reads: Vec<BamRecord> = Vec::new();
    let mut likelihoods: Option<HashMap<(usize, String), String>> = None;
    let mut compared = 0;
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        let got = match f[0] {
            "str" => {
                let max_period: usize = f[1].parse().unwrap();
                let pos: usize = f[3].parse().unwrap();
                let a = DragstrReadStrAnalyzer::of(f[2].as_bytes(), max_period).unwrap();
                let repeats: Vec<String> = (1..=max_period as i32)
                    .map(|p| a.number_of_repeats(pos, p).unwrap().to_string())
                    .collect();
                format!(
                    "str\t{}\t{}\t{pos}\t{}\t{}\t{}",
                    f[1],
                    f[2],
                    a.most_repeated_period(pos).unwrap(),
                    a.number_of_most_repeats(pos).unwrap(),
                    repeats.join(",")
                )
            }
            "param" => {
                let params = if f[1] == "default" {
                    &default
                } else {
                    parsed[f[1]].as_ref().unwrap()
                };
                param_row(
                    f[1],
                    params,
                    f[2],
                    f[3].parse().unwrap(),
                    f[4].parse().unwrap(),
                )
            }
            "paramtext" => {
                let text = f[2].replace('|', "\n");
                parsed.insert(f[1].to_string(), parse(&text, &format!("{}.txt", f[1])));
                line.to_string()
            }
            "parsed" => {
                let p = parsed[f[1]].as_ref().unwrap();
                format!(
                    "parsed\t{}\t{}\t{}",
                    f[1],
                    p.maximum_period(),
                    p.maximum_repeats()
                )
            }
            "error" => render_error(f[1], parsed[f[1]].as_ref().unwrap_err()),
            "impute" => {
                let params = if f[1] == "default" {
                    &default
                } else {
                    parsed[f[1]].as_ref().unwrap()
                };
                let (gop, gcp) = dragstr_impute(params, f[2].as_bytes()).unwrap();
                format!("impute\t{}\t{}\t{}\t{}", f[1], f[2], csv(&gop), csv(&gcp))
            }
            "nonsym" => {
                let n = |i: usize| -> u8 { f[i].parse().unwrap() };
                let (del, ins, gcp) = non_symmetrical_impute(
                    n(1),
                    n(2),
                    n(3),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                );
                format!(
                    "nonsym\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    f[1],
                    f[2],
                    f[3],
                    f[4],
                    f[5],
                    csv(&del),
                    csv(&ins),
                    csv(&gcp)
                )
            }
            "haplotype" => {
                haplotypes.push(Haplotype::new(f[2].as_bytes(), f[1] == "0").unwrap());
                line.to_string()
            }
            "read" => {
                reads.push(BamRecord {
                    read_name: f[1].to_string(),
                    reference_index: 0,
                    alignment_start: 1,
                    mapping_quality: f[4].parse().unwrap(),
                    cigar: htsjdk_bam::text_parse::parse_cigar(&format!("{}M", f[2].len()))
                        .unwrap(),
                    read_bases: f[2].as_bytes().to_vec(),
                    base_qualities: f[3].bytes().map(|b| b - 33).collect(),
                    ..Default::default()
                });
                line.to_string()
            }
            "lk" => {
                let table = likelihoods.get_or_insert_with(|| {
                    let engine = PairHmmLikelihoodEngine::new(LikelihoodEngineArguments {
                        gap_continuation_penalty: 10,
                        log10_global_read_mismapping_rate: -4.5,
                        pcr_error_model: PcrErrorModel::None,
                        modify_soft_clipped_bases: true,
                        dragstr_params: Some(DragstrParams::default_params()),
                        ..LikelihoodEngineArguments::default()
                    })
                    .unwrap();
                    let result = engine
                        .compute_read_likelihoods(
                            &haplotypes,
                            &["s1".to_string()],
                            &[reads.clone()],
                        )
                        .unwrap();
                    let mut table = HashMap::new();
                    for a in 0..result.number_of_alleles() {
                        for (r, read) in result.sample_evidence(0).unwrap().iter().enumerate() {
                            table.insert(
                                (a, read.read_name.clone()),
                                format!("{:x}", result.value(0, a, r).to_bits()),
                            );
                        }
                    }
                    table
                });
                let a: usize = f[1].parse().unwrap();
                format!("lk\t{a}\t{}\t{}", f[2], table[&(a, f[2].to_string())])
            }
            other => panic!("unexpected row {other}"),
        };
        assert_eq!(got, line);
        compared += 1;
    }
    assert!(compared > 2000, "{compared} rows");
}
