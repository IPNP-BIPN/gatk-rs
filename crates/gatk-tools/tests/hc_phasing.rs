//! Conformance for HaplotypeCaller's physical phasing (`AssemblyBasedCallerUtils.phaseCalls`)
//! against the oracle.
//!
//! Golden from `tools/readfilter-conformance/HcPhasingDump.java`: the calls of
//! `assignGenotypeLikelihoods` with phasing on and no annotations. Each case gives the haplotypes,
//! the reads and the likelihoods; this test genotypes and phases them and encodes every call as
//! the VCF line the golden holds.

use gatk_corpus as corpus;
use gatk_engine::allele_frequency_calculator::{AlleleFrequencyCalculator, Priors};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_tools::genotyping_engine::{Configuration, GenotypingEngine, SubsetMethod};
use gatk_tools::hc_genotyping::{assign_genotype_likelihoods, HcGenotypingArguments};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::encoder::{MissingFields, VcfEncoder};
use htsjdk_vcf::header::VcfHeader;

const REF: &str = "GATTACACGTTGCAAGTCCGATAGCTTAGGCATCGATCGGATCCATGCAAGTCAGTAAAAAGTCCATGAC";
const START: i32 = 1000;

fn from_bits(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}

#[derive(Default)]
struct Case<'a> {
    fields: Vec<&'a str>,
    haplotypes: Vec<Vec<&'a str>>,
    reads: Vec<Vec<&'a str>>,
    lk: Vec<Vec<&'a str>>,
}

fn run(case: &Case) -> Vec<String> {
    let f = &case.fields;
    let label = f[1];
    let samples: Vec<String> = f[2].split(',').map(str::to_string).collect();
    let ploidy: usize = f[3].parse().unwrap();
    let window = SimpleInterval::new("chr1", f[4].parse().unwrap(), f[5].parse().unwrap()).unwrap();
    let ref_loc = SimpleInterval::new("chr1", START, START + REF.len() as i32 - 1).unwrap();
    let mut haplotypes: Vec<Haplotype> = case
        .haplotypes
        .iter()
        .map(|h| {
            let mut hap = Haplotype::new(h[4].as_bytes(), h[3] == "true").unwrap();
            hap.set_cigar(&htsjdk_bam::text_parse::parse_cigar(h[5]).unwrap())
                .unwrap();
            hap.set_alignment_start_hap_wrt_ref(0);
            hap.set_genome_location(ref_loc.clone());
            hap.set_score(from_bits(h[6]));
            hap
        })
        .collect();
    let mut evidence: Vec<Vec<BamRecord>> = vec![Vec::new(); samples.len()];
    for r in &case.reads {
        let s = samples.iter().position(|x| x == r[2]).unwrap();
        evidence[s].push(BamRecord {
            read_name: r[3].to_string(),
            reference_index: 0,
            alignment_start: r[4].parse().unwrap(),
            mapping_quality: 60,
            cigar: htsjdk_bam::text_parse::parse_cigar(&format!("{}M", r[5])).unwrap(),
            read_bases: vec![b'A'; r[5].parse().unwrap()],
            base_qualities: vec![30; r[5].parse().unwrap()],
            ..Default::default()
        });
    }
    let mut values: Vec<Vec<Vec<f64>>> = vec![Vec::new(); samples.len()];
    for row in &case.lk {
        let s: usize = row[2].parse().unwrap();
        let parsed: Vec<f64> = if row.len() > 4 && !row[4].is_empty() {
            row[4].split(',').map(from_bits).collect()
        } else {
            Vec::new()
        };
        values[s].push(parsed);
    }
    let likelihoods = AlleleLikelihoods::new(
        SampleList::new(&samples),
        AlleleList::new(&haplotypes),
        evidence,
        values,
    )
    .unwrap();
    let priors = Priors {
        sample_ploidy: ploidy,
        ..Priors::default()
    };
    let mut engine = GenotypingEngine::new(
        Configuration {
            standard_confidence_for_calling: from_bits(f[6]),
            max_alternate_alleles: f[7].parse().unwrap(),
            sample_ploidy: ploidy,
            annotate_number_of_alleles_discovered: false,
            emit_all_active_sites: false,
            allele_specific: false,
            emit_all_confident_sites: false,
            annotate_all_sites_with_pls: false,
            force_keep_all_alleles: false,
            assignment_method: SubsetMethod::UsePlsToAssign,
        },
        AlleleFrequencyCalculator::make_calculator(&priors),
    );
    let arguments = HcGenotypingArguments {
        sample_ploidy: ploidy,
        disable_spanning_event_genotyping: f[8] != "true",
        do_physical_phasing: true,
        ..HcGenotypingArguments::default()
    };
    let calls = assign_genotype_likelihoods(
        &mut engine,
        &arguments,
        &mut haplotypes,
        &likelihoods,
        &samples,
        REF.as_bytes(),
        &ref_loc,
        &window,
        5000,
    )
    .unwrap();
    let header = VcfHeader {
        lines: Vec::new(),
        samples: samples.clone(),
    };
    let encoder = VcfEncoder::new(&header).with_missing_fields(MissingFields::Allow);
    let mut out: Vec<String> = calls
        .iter()
        .map(|vc| format!("call\t{label}\t{}", encoder.encode(vc).unwrap()))
        .collect();
    out.push(format!("calls\t{label}\t{}", calls.len()));
    out
}

#[test]
fn every_phased_call_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/hc_phasing.txt.gz"),
    );
    let mut cases: Vec<Case> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "case" => cases.push(Case {
                fields: f,
                ..Default::default()
            }),
            "haplotype" => cases.last_mut().unwrap().haplotypes.push(f),
            "read" => cases.last_mut().unwrap().reads.push(f),
            "lk" => cases.last_mut().unwrap().lk.push(f),
            _ => expected.push(line),
        }
    }
    let got: Vec<String> = cases.iter().flat_map(run).collect();
    assert_eq!(got.len(), expected.len(), "{}", got.join("\n"));
    for (g, e) in got.iter().zip(&expected) {
        assert_eq!(g, e);
    }
}
