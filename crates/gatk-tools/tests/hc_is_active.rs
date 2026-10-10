//! Conformance for `HaplotypeCallerEngine.isActive` against the oracle.
//!
//! Golden from `tools/readfilter-conformance/HcIsActiveDump.java`. Each case gives its samples,
//! its arguments and its reads; this test runs them through the locus iterator as the assembly
//! region walker does, evaluates every locus, and renders the `active` rows, which must be the
//! golden's in order.

use gatk_corpus as corpus;
use gatk_engine::locus_iterator::{self, LocusIteratorOptions};
use gatk_engine::read_states::{DownsamplingInfo, ReadStateManager};
use gatk_tools::haplotype_caller_engine::{ActiveRegionEvaluator, IsActiveArguments};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue, Tags};

fn bits(value: f64) -> String {
    format!("{:x}", value.to_bits())
}

struct Case<'a> {
    fields: Vec<&'a str>,
    reads: Vec<Vec<&'a str>>,
}

fn run(reference: &[u8], case: &Case) -> Vec<String> {
    let label = case.fields[1];
    let samples: Vec<String> = case.fields[2].split(',').map(str::to_string).collect();
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", reference.len() as _));
    for s in &samples {
        let mut group = ReadGroup::new(&format!("rg{s}"));
        group.attributes.set("SM", s);
        header.read_groups.push(group);
    }
    let reads: Vec<BamRecord> = case
        .reads
        .iter()
        .map(|f| {
            let mut tags = Tags::new();
            tags.insert(Tag::new(b"RG"), TagValue::Str(format!("rg{}", f[3])));
            BamRecord {
                read_name: f[2].to_string(),
                reference_index: 0,
                alignment_start: f[4].parse().unwrap(),
                mapping_quality: 60,
                cigar: htsjdk_bam::text_parse::parse_cigar(f[5]).unwrap(),
                read_bases: f[6].as_bytes().to_vec(),
                base_qualities: f[7].bytes().map(|b| b - 33).collect(),
                tags,
                ..Default::default()
            }
        })
        .collect();
    let arguments = IsActiveArguments {
        sample_ploidy: case.fields[3].parse().unwrap(),
        snp_heterozygosity: f64::from_bits(u64::from_str_radix(case.fields[4], 16).unwrap()),
        min_base_quality_score: case.fields[5].parse().unwrap(),
        ref_model_deletion_quality: case.fields[6].parse().unwrap(),
        ..IsActiveArguments::default()
    };
    let mut evaluator = ActiveRegionEvaluator::new(arguments, samples.clone());
    let sample_options: Vec<Option<String>> = samples.iter().cloned().map(Some).collect();
    let states = ReadStateManager::new(sample_options.clone(), DownsamplingInfo::NONE).unwrap();
    let contexts = locus_iterator::contexts(
        &reads,
        sample_options,
        &header,
        LocusIteratorOptions {
            include_deletions: true,
            include_ns: false,
        },
        states,
    )
    .unwrap();
    let mut out = Vec::new();
    for context in &contexts {
        let ref_base = reference[context.position as usize - 1];
        let state = evaluator.is_active(context, ref_base, &header).unwrap();
        out.push(format!(
            "active\t{label}\t{}\t{}\t{}\t{}\t{}\t{}",
            context.position,
            context.pileup.size(),
            bits(state.prob),
            state.kind.name(),
            state.result_value.map_or("null".to_string(), bits),
            bits(state.original_active_prob)
        ));
    }
    out
}

#[test]
fn every_locus_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/hc_is_active.txt.gz"),
    );
    let mut reference = Vec::new();
    let mut cases: Vec<Case> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "reference" => reference = f[1].as_bytes().to_vec(),
            "case" => cases.push(Case {
                fields: f,
                reads: Vec::new(),
            }),
            "read" => cases.last_mut().unwrap().reads.push(f),
            _ => expected.push(line),
        }
    }
    let got: Vec<String> = cases.iter().flat_map(|c| run(&reference, c)).collect();
    assert_eq!(got.len(), expected.len());
    for (g, e) in got.iter().zip(&expected) {
        assert_eq!(g, e);
    }
}
