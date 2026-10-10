//! Conformance for `AssemblyBasedCallerUtils.finalizeRegion` against the oracle.
//!
//! Golden from `tools/readfilter-conformance/FinalizeRegionDump.java`. The golden carries the
//! contig, every read and each case's switches; this test rebuilds the region, finalizes it, and
//! renders the `kept`, `region` and `error` rows, which must be the golden's in order.

use gatk_corpus as corpus;
use gatk_engine::assembly_based_caller_utils::{finalize_region, FinalizeArguments};
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::interval::SimpleInterval;
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue, Tags};

fn phred(q: &[u8]) -> String {
    q.iter().map(|&b| (b + 33) as char).collect()
}

fn int_tag(read: &BamRecord, tag: &[u8; 2]) -> String {
    match read.tags.get(Tag::new(tag)) {
        Some(TagValue::Int(v)) => v.to_string(),
        Some(other) => format!("{other:?}"),
        None => "null".to_string(),
    }
}

fn record(f: &[&str]) -> BamRecord {
    let mut tags = Tags::new();
    tags.insert(Tag::new(b"RG"), TagValue::Str(format!("rg-{}", f[7])));
    let flags: u16 = f[2].parse().unwrap();
    let paired = flags & 1 != 0;
    BamRecord {
        read_name: f[1].to_string(),
        flags,
        reference_index: 0,
        alignment_start: f[3].parse().unwrap(),
        mapping_quality: 60,
        cigar: htsjdk_bam::text_parse::parse_cigar(f[4]).unwrap(),
        mate_reference_index: if paired { 0 } else { -1 },
        mate_alignment_start: f[5].parse().unwrap(),
        inferred_insert_size: f[6].parse().unwrap(),
        read_bases: f[8].as_bytes().to_vec(),
        base_qualities: f[9].bytes().map(|b| b - 33).collect(),
        tags,
    }
}

fn kept(label: &str, list: &str, reads: &[BamRecord], out: &mut Vec<String>) {
    for (i, r) in reads.iter().enumerate() {
        out.push(format!(
            "kept\t{label}\t{list}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.read_name,
            r.alignment_start,
            r.cigar.to_text(),
            String::from_utf8_lossy(&r.read_bases),
            phred(&r.base_qualities),
            int_tag(r, b"os"),
            int_tag(r, b"oe"),
        ));
    }
}

#[test]
fn every_finalization_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/finalize_region.txt.gz"),
    );
    let mut contig_length = 0usize;
    let mut reads: Vec<BamRecord> = Vec::new();
    let mut cases: Vec<Vec<&str>> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "contig" => contig_length = f[1].len(),
            "read" => reads.push(record(&f)),
            "case" => cases.push(f),
            _ => expected.push(line),
        }
    }
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", contig_length as _));
    for s in ["s1", "s2"] {
        let mut group = ReadGroup::new(&format!("rg-{s}"));
        group.attributes.set("SM", s);
        header.read_groups.push(group);
    }
    let samples = vec!["s1".to_string(), "s2".to_string()];

    let mut got = Vec::new();
    for case in &cases {
        let label = case[1];
        let flag = |i: usize| case[i] == "true";
        let mut region = AssemblyRegion::with_padding(
            SimpleInterval::new("chr1", 200, 300).unwrap(),
            true,
            50,
            &header,
        )
        .unwrap();
        let padded = region.padded_span().clone();
        for r in &reads {
            if padded.overlaps(
                "chr1",
                gatk_engine::read_utils::start(r),
                gatk_engine::read_utils::end(r),
            ) {
                region.add(r.clone(), &header).unwrap();
            }
        }
        if flag(9) {
            region.set_finalized(true);
        }
        let arguments = FinalizeArguments {
            error_correct_reads: flag(2),
            dont_use_soft_clipped_bases: flag(3),
            min_tail_quality: case[4].parse().unwrap(),
            correct_overlapping_base_qualities: flag(5),
            soft_clip_low_quality_ends: flag(6),
            override_softclip_fragment_check: flag(7),
            track_hardclipped_reads: flag(8),
        };
        match finalize_region(&mut region, &arguments, &header, &samples) {
            Ok(()) => {
                kept(label, "reads", region.reads(), &mut got);
                kept(
                    label,
                    "hardclipped",
                    region.hard_clipped_pileup_reads(),
                    &mut got,
                );
                got.push(format!(
                    "region\t{label}\t{}\t{}\t{}",
                    region.reads().len(),
                    region.hard_clipped_pileup_reads().len(),
                    region.is_finalized()
                ));
            }
            Err(e) => got.push(format!("error\t{label}\t{e:?}")),
        }
    }
    assert_eq!(got.len(), expected.len(), "{}", got.join("\n"));
    for (g, e) in got.iter().zip(&expected) {
        assert_eq!(g, e);
    }
}
