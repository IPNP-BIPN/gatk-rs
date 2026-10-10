//! Conformance for `ReadThreadingAssembler.runLocalAssembly` against the oracle.
//!
//! Golden from `tools/readfilter-conformance/ReadThreadingAssemblerDump.java`. Each case lays its
//! inputs out before its outputs: the contig, the region, the assembler's settings and every read.
//! This test rebuilds the region from those rows, runs the port, and renders the `result`,
//! `haplotype`, `set` and `error` rows the dump prints, which must be the golden's, in order.

use gatk_corpus as corpus;
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::read_threading_assembler::{
    AssemblerError, AssemblerSettings, ReadThreadingAssembler,
};
use gatk_engine::smith_waterman::{SmithWatermanJavaAligner, NEW_SW_PARAMETERS, STANDARD_NGS};
use gatk_engine::tsv_table::java_double_to_string;
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue, Tags};

const CONTIG: &str = "chr1";

/// The rows of one case, inputs and outputs.
#[derive(Default)]
struct Case {
    label: String,
    contig: String,
    active_start: i32,
    active_end: i32,
    padding: i32,
    reference_padding: i32,
    settings: Vec<(String, String)>,
    reads: Vec<Vec<String>>,
    outputs: Vec<String>,
}

fn cases(golden: &str) -> Vec<Case> {
    let mut cases: Vec<Case> = Vec::new();
    for line in golden.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[0] {
            "case" => cases.push(Case {
                label: fields[1].to_string(),
                contig: fields[2].to_string(),
                active_start: fields[3].parse().unwrap(),
                active_end: fields[4].parse().unwrap(),
                padding: fields[5].parse().unwrap(),
                reference_padding: fields[6].parse().unwrap(),
                ..Default::default()
            }),
            "settings" => {
                let case = cases.last_mut().unwrap();
                for pair in &fields[2..] {
                    let (key, value) = pair.split_once('=').unwrap();
                    case.settings.push((key.to_string(), value.to_string()));
                }
            }
            "read" => cases
                .last_mut()
                .unwrap()
                .reads
                .push(fields[2..].iter().map(|f| f.to_string()).collect()),
            _ => cases.last_mut().unwrap().outputs.push(line.to_string()),
        }
    }
    cases
}

impl Case {
    fn setting(&self, key: &str) -> &str {
        &self
            .settings
            .iter()
            .find(|(k, _)| k == key)
            .unwrap_or_else(|| panic!("{}: no setting {key}", self.label))
            .1
    }

    fn flag(&self, key: &str) -> bool {
        self.setting(key).parse().unwrap()
    }

    fn int(&self, key: &str) -> i32 {
        self.setting(key).parse().unwrap()
    }

    fn bits(&self, key: &str) -> f64 {
        f64::from_bits(u64::from_str_radix(self.setting(key), 16).unwrap())
    }
}

fn header(contig_length: usize) -> SamHeader {
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new(CONTIG, contig_length as _));
    for sample in ["s1", "s2"] {
        let mut group = ReadGroup::new(&format!("rg-{sample}"));
        group.attributes.set("SM", sample);
        header.read_groups.push(group);
    }
    header
}

fn record(fields: &[String]) -> BamRecord {
    let mut tags = Tags::new();
    tags.insert(Tag::new(b"RG"), TagValue::Str(format!("rg-{}", fields[1])));
    BamRecord {
        read_name: fields[0].clone(),
        reference_index: 0,
        alignment_start: fields[2].parse().unwrap(),
        mapping_quality: 60,
        cigar: htsjdk_bam::text_parse::parse_cigar(&fields[3]).expect("a cigar"),
        read_bases: fields[4].as_bytes().to_vec(),
        base_qualities: fields[5].bytes().map(|q| q - 33).collect(),
        tags,
        ..Default::default()
    }
}

fn render_error(label: &str, stage: &str, error: &AssemblerError) -> String {
    format!("error\t{label}\t{stage}\t{error:?}")
}

fn render_location(location: Option<&SimpleInterval>) -> String {
    location.map_or("null".to_string(), |l| {
        format!("{}:{}-{}", l.contig, l.start, l.end)
    })
}

fn render_cigar(cigar: Option<&Cigar>) -> String {
    cigar.map_or("null".to_string(), |c| c.to_text())
}

/// `ReferenceConfidenceModel.createReferenceHaplotype`.
fn reference_haplotype(
    region: &AssemblyRegion,
    ref_bases: &[u8],
    ref_loc: &SimpleInterval,
) -> Haplotype {
    let mut h = Haplotype::new(ref_bases, true).unwrap();
    h.set_genome_location(region.padded_span().clone());
    h.set_alignment_start_hap_wrt_ref(region.padded_span().start - ref_loc.start);
    h.set_cigar(&Cigar::new(vec![CigarElement {
        length: ref_bases.len() as u32,
        op: Op::M,
    }]))
    .unwrap();
    h
}

fn run(case: &Case) -> Vec<String> {
    let label = &case.label;
    let header = header(case.contig.len());
    let active = SimpleInterval::new(CONTIG, case.active_start, case.active_end).unwrap();
    let mut region = AssemblyRegion::with_padding(active, true, case.padding, &header).unwrap();
    for fields in &case.reads {
        region.add(record(fields), &header).unwrap();
    }
    let padded = region.padded_span().clone();
    let ref_loc = SimpleInterval::new(
        CONTIG,
        (padded.start - case.reference_padding).max(1),
        (padded.end + case.reference_padding).min(case.contig.len() as i32),
    )
    .unwrap();
    let contig = case.contig.as_bytes();
    let full_ref = &contig[ref_loc.start as usize - 1..ref_loc.end as usize];
    let ref_bases = &contig[padded.start as usize - 1..padded.end as usize];

    let kmer_sizes = case
        .setting("kmers")
        .split(',')
        .map(|k| k.parse().unwrap())
        .collect();
    let settings = AssemblerSettings {
        max_allowed_paths: case.int("maxPaths") as usize,
        kmer_sizes,
        dont_increase_kmer_sizes_for_cycles: case.flag("dontIncrease"),
        allow_non_unique_kmers_in_ref: case.flag("allowNonUnique"),
        num_pruning_samples: case.int("pruningSamples") as usize,
        prune_factor: case.int("pruneFactor"),
        use_adaptive_pruning: case.flag("adaptive"),
        initial_error_rate_for_pruning: case.bits("initialErrorRate"),
        pruning_log_odds_threshold: case.bits("pruningLogOdds"),
        pruning_seeding_log_odds_threshold: case.bits("seedingLogOdds"),
        max_unpruned_variants: case.int("maxUnprunedVariants"),
        enable_legacy_graph_cycle_detection: case.flag("legacyCycles"),
        min_matching_bases_to_dangling_end_recovery: case.int("minMatching"),
    };
    let mut assembler = ReadThreadingAssembler::new(&settings).unwrap();
    assembler.set_recover_dangling_branches(case.flag("recoverDangling"));
    assembler.set_recover_all_dangling_branches(case.flag("recoverAll"));
    assembler.set_min_dangling_branch_length(case.int("minDanglingLength"));
    assembler.set_min_base_quality_to_use_in_assembly(case.int("minBaseQuality") as u8);
    let aligner = SmithWatermanJavaAligner;

    let mut out = Vec::new();
    let reads = assembler.hard_clipped_reads(&region, &header).unwrap();
    match assembler.assemble(
        &reads,
        &reference_haplotype(&region, ref_bases, &ref_loc),
        &aligner,
        &STANDARD_NGS,
    ) {
        Ok(results) => {
            for (i, r) in results.iter().enumerate() {
                let graph = r.seq_graph.as_ref();
                out.push(format!(
                    "result\t{label}\t{i}\t{}\t{}\t{}\t{}",
                    r.kmer_size.map_or("-".to_string(), |k| k.to_string()),
                    r.status.name(),
                    graph.map_or("-".to_string(), |g| g.vertex_count().to_string()),
                    graph.map_or("-".to_string(), |g| g.edge_count().to_string()),
                ));
            }
        }
        Err(e) => out.push(render_error(label, "assemble", &e)),
    }

    match assembler.run_local_assembly(
        &region,
        reference_haplotype(&region, ref_bases, &ref_loc),
        full_ref,
        &ref_loc,
        &header,
        &aligner,
        &STANDARD_NGS,
        &NEW_SW_PARAMETERS,
    ) {
        Ok(set) => {
            for (i, h) in set.haplotype_list().iter().enumerate() {
                out.push(format!(
                    "haplotype\t{label}\t{i}\t{}\t{}\t{}\t{}\t{:x},{}\t{}\t{}",
                    h.is_reference(),
                    String::from_utf8(h.bases()).unwrap(),
                    render_cigar(h.cigar()),
                    h.alignment_start_hap_wrt_ref(),
                    h.score().to_bits(),
                    java_double_to_string(h.score()),
                    h.kmer_size(),
                    render_location(h.genome_location()),
                ));
            }
            out.push(format!(
                "set\t{label}\t{}\t{}\t{}",
                set.haplotype_count(),
                set.is_variation_present(),
                set.reference_index().map_or(-1, |i| i as i64),
            ));
        }
        Err(e) => out.push(render_error(label, "run", &e)),
    }
    out
}

#[test]
fn every_case_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/read_threading_assembler.txt.gz"),
    );
    let cases = cases(&golden);
    assert!(cases.len() >= 40, "{} cases", cases.len());
    let mut compared = 0;
    for case in &cases {
        let got = run(case);
        assert_eq!(
            got.len(),
            case.outputs.len(),
            "{}: rows\n{}",
            case.label,
            got.join("\n")
        );
        for (g, e) in got.iter().zip(&case.outputs) {
            assert_eq!(g, e, "{}", case.label);
            compared += 1;
        }
    }
    assert!(compared >= 364, "{compared} rows compared");
}
