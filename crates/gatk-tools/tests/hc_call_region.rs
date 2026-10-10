//! Conformance for `HaplotypeCallerEngine.callRegion` against the oracle.
//!
//! Golden from `tools/readfilter-conformance/HcCallRegionDump.java`: whole active regions called
//! with HaplotypeCaller's defaults and default annotations, the random generator reset per region.
//! This test rebuilds each region from the golden's reads and encodes every call as its VCF line.

use gatk_annotation::catalogue;
use gatk_corpus as corpus;
use gatk_engine::allele_frequency_calculator::{AlleleFrequencyCalculator, Priors};
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::assembly_region_trimmer::{AssemblyRegionTrimmer, TrimmerArguments};
use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_random::JavaRandom;
use gatk_engine::pair_hmm_likelihood_engine::{LikelihoodEngineArguments, PairHmmLikelihoodEngine};
use gatk_engine::smith_waterman::SmithWatermanJavaAligner;
use gatk_tools::genotyping_engine::{Configuration, GenotypingEngine, SubsetMethod};
use gatk_tools::hc_genotyping::{CallAnnotator, HcGenotypingArguments};
use gatk_tools::hc_region::{
    call_region, haplotype_caller_assembler, CallRegionContext, RegionAssemblyArguments,
};
use gatk_tools::variant_annotator_engine::Engine;
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue, Tags};
use htsjdk_vcf::encoder::{MissingFields, VcfEncoder};
use htsjdk_vcf::header::VcfHeader;

const ANNOTATIONS: [&str; 13] = [
    "BaseQualityRankSumTest",
    "ChromosomeCounts",
    "Coverage",
    "DepthPerAlleleBySample",
    "DepthPerSampleHC",
    "ExcessHet",
    "FisherStrand",
    "InbreedingCoeff",
    "MappingQualityRankSumTest",
    "QualByDepth",
    "RMSMappingQuality",
    "ReadPosRankSumTest",
    "StrandOddsRatio",
];

fn run(reference: &[u8], case: &[&str], reads: &[Vec<&str>]) -> Vec<String> {
    let label = case[1];
    let samples: Vec<String> = case[5].split(',').map(str::to_string).collect();
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", reference.len() as _));
    for s in &samples {
        let mut group = ReadGroup::new(&format!("rg{s}"));
        group.attributes.set("SM", s);
        header.read_groups.push(group);
    }
    let active =
        SimpleInterval::new("chr1", case[2].parse().unwrap(), case[3].parse().unwrap()).unwrap();
    let mut region =
        AssemblyRegion::with_padding(active, true, case[4].parse().unwrap(), &header).unwrap();
    for r in reads {
        let mut tags = Tags::new();
        tags.insert(Tag::new(b"RG"), TagValue::Str(format!("rg{}", r[3])));
        let record = BamRecord {
            read_name: r[2].to_string(),
            flags: r[8].parse().unwrap(),
            reference_index: 0,
            alignment_start: r[4].parse().unwrap(),
            mapping_quality: r[9].parse().unwrap(),
            cigar: htsjdk_bam::text_parse::parse_cigar(r[5]).unwrap(),
            read_bases: r[6].as_bytes().to_vec(),
            base_qualities: r[7].bytes().map(|b| b - 33).collect(),
            tags,
            ..Default::default()
        };
        region.add(record, &header).unwrap();
    }
    let assembler = haplotype_caller_assembler(10);
    let trimmer =
        AssemblyRegionTrimmer::new(TrimmerArguments::default(), reference.len() as i32).unwrap();
    let likelihood_engine = PairHmmLikelihoodEngine::new(LikelihoodEngineArguments {
        modify_soft_clipped_bases: true,
        ..LikelihoodEngineArguments::default()
    })
    .unwrap();
    let genotyping = HcGenotypingArguments {
        do_physical_phasing: true,
        ..HcGenotypingArguments::default()
    };
    let context = CallRegionContext {
        header: &header,
        samples: &samples,
        contig_bases: reference,
        assembly: &RegionAssemblyArguments::default(),
        assembler: &assembler,
        aligner: &SmithWatermanJavaAligner,
        trimmer: &trimmer,
        likelihood_engine: &likelihood_engine,
        genotyping: &genotyping,
        mapping_quality_threshold: 20,
    };
    let mut engine = GenotypingEngine::new(
        Configuration {
            standard_confidence_for_calling: 30.0,
            max_alternate_alleles: 6,
            sample_ploidy: 2,
            annotate_number_of_alleles_discovered: false,
            emit_all_active_sites: false,
            allele_specific: false,
            emit_all_confident_sites: false,
            annotate_all_sites_with_pls: false,
            force_keep_all_alleles: false,
            assignment_method: SubsetMethod::UsePlsToAssign,
        },
        AlleleFrequencyCalculator::make_calculator(&Priors::default()),
    );
    let annotator = Engine {
        resolved: ANNOTATIONS
            .iter()
            .map(|n| catalogue::entry(n).unwrap())
            .collect(),
        overlap_names: Vec::new(),
        expressions: Vec::new(),
        allele_concordance: false,
    };
    let mut random = JavaRandom::new(47382911);
    let calls = call_region(
        region,
        &context,
        &mut engine,
        Some(CallAnnotator {
            engine: &annotator,
            random: &mut random,
        }),
    )
    .unwrap();
    let vcf_header = VcfHeader {
        lines: Vec::new(),
        samples: samples.clone(),
    };
    let encoder = VcfEncoder::new(&vcf_header).with_missing_fields(MissingFields::Allow);
    let mut out: Vec<String> = calls
        .iter()
        .map(|vc| format!("call\t{label}\t{}", encoder.encode(vc).unwrap()))
        .collect();
    out.push(format!("calls\t{label}\t{}", calls.len()));
    out
}

#[test]
fn every_region_call_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/hc_call_region.txt.gz"),
    );
    let mut reference: Vec<u8> = Vec::new();
    let mut cases: Vec<(Vec<&str>, Vec<Vec<&str>>)> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "reference" => reference = f[1].as_bytes().to_vec(),
            "case" => cases.push((f, Vec::new())),
            "read" => cases.last_mut().unwrap().1.push(f),
            _ => expected.push(line),
        }
    }
    let got: Vec<String> = cases
        .iter()
        .flat_map(|(c, r)| run(&reference, c, r))
        .collect();
    for (g, e) in got.iter().zip(&expected) {
        assert_eq!(g, e);
    }
    assert_eq!(got.len(), expected.len());
}
