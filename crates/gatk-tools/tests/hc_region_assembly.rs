//! Conformance for HaplotypeCaller's region assembly against the oracle.
//!
//! Golden from `tools/readfilter-conformance/HcRegionAssemblyDump.java`: `assembleReads`,
//! `getVariationEvents`, the trimmer and `trimTo` over a region whose reads the golden lists. This
//! test rebuilds each region, runs the same steps, and renders the `finalized`, `assembled`,
//! `event`, `trim`, `trimmed`, `set` and `error` rows, which must be the golden's in order.

use gatk_corpus as corpus;
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::assembly_region_trimmer::{AssemblyRegionTrimmer, TrimmerArguments};
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::smith_waterman::SmithWatermanJavaAligner;
use gatk_tools::hc_region::{assemble_reads, haplotype_caller_assembler, RegionAssemblyArguments};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue, Tags};

fn loc(i: &SimpleInterval) -> String {
    format!("{}:{}-{}", i.contig, i.start, i.end)
}

fn show(h: &Haplotype) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}",
        h.is_reference(),
        String::from_utf8(h.bases()).unwrap(),
        h.cigar().map_or("null".to_string(), |c| c.to_text()),
        h.alignment_start_hap_wrt_ref(),
        h.genome_location().map_or("null".to_string(), loc)
    )
}

fn run(reference: &[u8], case: &[&str], reads: &[Vec<&str>]) -> Vec<String> {
    let label = case[1];
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", reference.len() as _));
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "s1");
    header.read_groups.push(group);
    let active =
        SimpleInterval::new("chr1", case[2].parse().unwrap(), case[3].parse().unwrap()).unwrap();
    let mut region =
        AssemblyRegion::with_padding(active, true, case[4].parse().unwrap(), &header).unwrap();
    for r in reads {
        let mut tags = Tags::new();
        tags.insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
        let record = BamRecord {
            read_name: r[2].to_string(),
            flags: r[8].parse().unwrap(),
            reference_index: 0,
            alignment_start: r[4].parse().unwrap(),
            mapping_quality: 60,
            cigar: htsjdk_bam::text_parse::parse_cigar(r[5]).unwrap(),
            read_bases: r[6].as_bytes().to_vec(),
            base_qualities: r[7].bytes().map(|b| b - 33).collect(),
            tags,
            ..Default::default()
        };
        region.add(record, &header).unwrap();
    }
    let mut out = Vec::new();
    let assembler = haplotype_caller_assembler(10);
    let samples = vec!["s1".to_string()];
    let mut untrimmed = match assemble_reads(
        &mut region,
        &RegionAssemblyArguments::default(),
        &header,
        &samples,
        reference,
        &assembler,
        &SmithWatermanJavaAligner,
    ) {
        Ok(set) => set,
        Err(e) => {
            out.push(format!("error\t{label}\t{e:?}"));
            return out;
        }
    };
    for r in region.reads() {
        out.push(format!(
            "finalized\t{label}\t{}\t{}\t{}\t{}\t{}",
            r.read_name,
            r.alignment_start,
            r.cigar.to_text(),
            String::from_utf8_lossy(&r.read_bases),
            r.base_qualities
                .iter()
                .map(|q| (q + 33) as char)
                .collect::<String>()
        ));
    }
    for (i, h) in untrimmed.haplotype_list().iter().enumerate() {
        out.push(format!("assembled\t{label}\t{i}\t{}", show(h)));
    }
    let events = untrimmed.variation_events(0).unwrap();
    for e in &events {
        out.push(format!(
            "event\t{label}\t{}\t{}\t{}",
            e.start(),
            e.ref_allele().display_string(),
            e.alt_allele().display_string()
        ));
    }
    let trimmer =
        AssemblyRegionTrimmer::new(TrimmerArguments::default(), reference.len() as i32).unwrap();
    let padded = region.padded_span().clone();
    let padded_bases = &reference[padded.start as usize - 1..padded.end as usize];
    let result = trimmer
        .trim(&region, &events, padded_bases, padded.start)
        .unwrap();
    if result.variant_span.is_none() {
        out.push(format!("trim\t{label}\tnull\tnull"));
        return out;
    }
    let variant_region = trimmer.variant_region(&result, &region, &header).unwrap();
    out.push(format!(
        "trim\t{label}\t{}\t{}",
        loc(variant_region.span()),
        loc(variant_region.padded_span())
    ));
    let trimmed = untrimmed.trim_to(variant_region).unwrap();
    for (i, h) in trimmed.haplotype_list().iter().enumerate() {
        out.push(format!(
            "trimmed\t{label}\t{i}\t{}\t{}",
            show(h),
            h.kmer_size()
        ));
    }
    let genotyping = trimmed.region_for_genotyping().unwrap();
    out.push(format!(
        "set\t{label}\t{}\t{}\t{}\t{}\t{}",
        trimmed.haplotype_count(),
        trimmed.is_variation_present(),
        loc(genotyping.span()),
        loc(genotyping.padded_span()),
        genotyping
            .reads()
            .iter()
            .map(|r| r.read_name.as_str())
            .collect::<Vec<_>>()
            .join(",")
    ));
    out
}

#[test]
fn every_region_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/hc_region_assembly.txt.gz"),
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
