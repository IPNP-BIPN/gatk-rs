//! `traverse_with_pileups` walking a shard a window at a time against walking it whole.
//!
//! The windowed walk holds only the reads overlapping its window and replays the downsampler from a
//! copy of the generator to fill the regions, so everything a tool can see has to come out the
//! same: every evaluator call with its pileup, every region with its reads, and the generator's
//! state when the first region is handed out.

use gatk_engine::assembly_region_iterator::AssemblyRegionArgs;
use gatk_engine::assembly_region_walker::{traverse_with_pileups, Downsampler, WalkerReads};
use gatk_engine::downsampling::positional_downsample_with;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_random::JavaRandom;
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    header.sequences.push(SequenceRecord::new("chr1", 6000));
    header.sequences.push(SequenceRecord::new("chr2", 3000));
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "sample1");
    header.read_groups.push(group);
    header
}

/// Coordinate-sorted reads with piles at one start (so the downsampler draws), deletions,
/// insertions, soft clips, long spans and coverage gaps.
fn reads() -> Vec<BamRecord> {
    let mut random = JavaRandom::new(7);
    let mut out = Vec::new();
    for (reference_index, length) in [(0, 6000), (1, 3000)] {
        let mut start = 1;
        while start < length - 400 {
            let pile = if random.next_int_bound(9) == 0 {
                8
            } else {
                1 + random.next_int_bound(3)
            };
            for _ in 0..pile {
                let cigar = match random.next_int_bound(6) {
                    0 => "20M5D30M",
                    1 => "10S40M",
                    2 => "25M3I22M",
                    3 => "30M200N20M",
                    _ => "50M",
                };
                let cigar = htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar");
                let bases = cigar.read_length() as usize;
                let mut tags = htsjdk_bam::tag::Tags::new();
                tags.insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
                out.push(BamRecord {
                    read_name: format!("r{}", out.len()),
                    reference_index,
                    alignment_start: start,
                    mapping_quality: 60,
                    read_bases: (0..bases)
                        .map(|i| b"ACGT"[(i * 7 + start as usize) % 4])
                        .collect(),
                    base_qualities: vec![30; bases],
                    cigar,
                    tags,
                    ..Default::default()
                });
            }
            // Mostly dense, now and then a gap no read covers.
            start += if random.next_int_bound(40) == 0 {
                300
            } else {
                1 + random.next_int_bound(12)
            };
        }
    }
    out
}

/// Everything the tool sees: the evaluator's calls, then each region and its reads, then the
/// generator as the regions are drawn.
fn run(window: Option<i32>, intervals: &[SimpleInterval]) -> Vec<String> {
    let header = header();
    let reads = reads();
    let mut walker = WalkerReads::held(&reads, &header);
    walker.window = window;
    let args = AssemblyRegionArgs {
        max_reads_per_alignment_start: 3,
        ..Default::default()
    };
    let new_batch = |_: &BamRecord, previous: &BamRecord, next: &BamRecord| {
        (previous.reference_index, previous.alignment_start)
            != (next.reference_index, next.alignment_start)
    };
    let run = |reads: &[BamRecord], random: &mut JavaRandom| {
        positional_downsample_with(reads, 3, random, false)
    };
    let mut random = JavaRandom::new(42);
    let mut seen = Vec::new();
    let groups = gatk_engine::assembly_region_iterator::group_intervals_by_contig(intervals);
    for group in groups {
        let regions = traverse_with_pileups(
            &walker,
            &group,
            &[Some("sample1".to_string())],
            &args,
            &header,
            Downsampler {
                new_batch: &new_batch,
                run: &run,
            },
            &mut random,
            &mut |locus, context| {
                let names: Vec<&str> = context
                    .map(|c| {
                        c.pileup
                            .elements
                            .iter()
                            .map(|e| e.read.read_name.as_str())
                            .collect()
                    })
                    .unwrap_or_default();
                seen.push(format!(
                    "{}:{} {}",
                    locus.contig,
                    locus.start,
                    names.join(",")
                ));
                (
                    (names.len().is_multiple_of(5) && !names.is_empty()) as u8 as f64,
                    None,
                )
            },
        )
        .expect("a traversal");
        for region in regions.into_iter().flatten() {
            let region = region.expect("a region");
            let names: Vec<&str> = region
                .region
                .reads()
                .iter()
                .map(|r| r.read_name.as_str())
                .collect();
            seen.push(format!(
                "{:?} {:?} {} {} next={}",
                region.region.span(),
                region.region.padded_span(),
                region.is_active,
                names.join(","),
                // The tool's own draw while calling the region.
                random.next_int(),
            ));
        }
    }
    seen
}

#[test]
fn windows_walk_what_the_whole_shard_walks() {
    let intervals = [
        SimpleInterval::new("chr1", 1, 2200).unwrap(),
        SimpleInterval::new("chr1", 2500, 5000).unwrap(),
        SimpleInterval::new("chr2", 40, 2900).unwrap(),
    ];
    let whole = run(None, &intervals);
    assert!(
        whole.iter().any(|line| line.contains("true")),
        "some region is active"
    );
    for window in [1, 7, 50, 333, 5000] {
        let windowed = run(Some(window), &intervals);
        assert_eq!(windowed.len(), whole.len(), "window {window}");
        for (a, b) in windowed.iter().zip(&whole) {
            assert_eq!(a, b, "window {window}");
        }
    }
}
