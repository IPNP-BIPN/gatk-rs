//! Conformance for the assembly-region traversal against the oracle.
//!
//! Goldens from `tools/readfilter-conformance/AssemblyRegionIteratorDump.java`.
//!
//! The rows that carry the suite:
//!
//! ```text
//! shard  shard-adjacent           chr1:100-300              chr1:90-310
//! shard  shard-merged-by-padding  chr1:100-200,chr1:216-300 chr1:90-310
//! region trav-whole-contig  2  chr1:65-114|chr1:55-124|true|3
//! ```
//!
//! Two intervals that touch are **already one** in the unpadded list, because
//! `getIntervalsWithFlanks` sorts and merges whatever the padding is; two intervals fifteen bases
//! apart stay two there and become one once padded. And the first active region starts at 65 while
//! the first read starts at 101, because the band pass filter spreads probability backwards.

use gatk_corpus as corpus;
use gatk_engine::assembly_region_iterator::{
    assembly_regions, group_intervals_by_contig, AssemblyRegionArgs, ReadShard,
};
use gatk_engine::interval::SimpleInterval;
use gatk_engine::locus_iterator::{self, LocusIteratorOptions};
use gatk_engine::read_states::{DownsamplingInfo, ReadStateManager};
use gatk_engine::variant_source::Located;
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const CHR1: i32 = 1000;
const CHR2: i32 = 500;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/assembly_region_iterator.txt.gz"),
    )
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    header.sequences.push(SequenceRecord::new("chr1", CHR1));
    header.sequences.push(SequenceRecord::new("chr2", CHR2));
    for (id, sample) in [("rg1", "sampleA"), ("rg2", "sampleB")] {
        let mut group = ReadGroup::new(id);
        group.attributes.set("SM", sample);
        header.read_groups.push(group);
    }
    header
}

fn interval(contig: &str, start: i32, end: i32) -> SimpleInterval {
    SimpleInterval::new(contig, start, end).expect("a valid interval")
}

fn read(name: &str, group: &str, cigar: &str, start: i32) -> BamRecord {
    let cigar = htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar");
    let length = cigar.read_length() as usize;
    let mut tags = htsjdk_bam::tag::Tags::new();
    tags.insert(Tag::new(b"RG"), TagValue::Str(group.to_string()));
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality: 60,
        read_bases: (0..length).map(|i| b"ACGT"[i % 4]).collect(),
        base_qualities: vec![30; length],
        cigar,
        tags,
        ..Default::default()
    }
}

/// The dump's fixture: a covered block at 101-140, a hole, a second at 201-230, a lone read at 400.
fn reads() -> Vec<BamRecord> {
    vec![
        read("a1", "rg1", "20M", 101),
        read("a2", "rg1", "20M", 111),
        read("a3", "rg2", "20M", 121),
        read("b1", "rg1", "20M", 201),
        read("b2", "rg2", "20M", 211),
        read("c1", "rg1", "20M", 400),
    ]
}

/// Label, padding, intervals: the `shard` cases in the dump's order.
fn shard_case(label: &str) -> (i32, Vec<SimpleInterval>) {
    match label {
        "shard-one" => (10, vec![interval("chr1", 100, 200)]),
        "shard-adjacent" => (
            10,
            vec![interval("chr1", 100, 200), interval("chr1", 201, 300)],
        ),
        "shard-overlapping" => (
            10,
            vec![interval("chr1", 100, 220), interval("chr1", 200, 300)],
        ),
        "shard-unsorted" => (
            10,
            vec![interval("chr1", 300, 400), interval("chr1", 100, 200)],
        ),
        "shard-merged-by-padding" => (
            10,
            vec![interval("chr1", 100, 200), interval("chr1", 216, 300)],
        ),
        "shard-off-contig" => (10, vec![interval("chr1", 5, 20)]),
        "shard-zero-padding" => (
            0,
            vec![interval("chr1", 100, 200), interval("chr1", 201, 300)],
        ),
        "shard-negative-padding" => (-1, vec![interval("chr1", 100, 200)]),
        other => panic!("unknown shard case {other}"),
    }
}

/// Label, and the six validated fields in the dump's order.
const ARGS: &[(&str, i32, i32, i32, i32, i32, i32)] = &[
    ("args-default", 5, 50, 10, 50, 20, 75),
    ("args-zero-min", 0, 50, 10, 50, 20, 75),
    ("args-zero-max", 5, 0, 10, 50, 20, 75),
    ("args-min-above-max", 60, 50, 10, 50, 20, 75),
    ("args-negative-padding", 5, 50, -1, 50, 20, 75),
    ("args-negative-max-reads", 5, 50, 10, -1, 20, 75),
    ("args-negative-snp-padding", 5, 50, 10, 50, -1, 75),
    ("args-negative-indel-padding", 5, 50, 10, 50, 20, -1),
    ("args-two-wrong", 0, 50, -1, 50, 20, 75),
];

/// Label, padding, min size, max size, intervals.
fn traversal_case(label: &str) -> (i32, i32, i32, Vec<SimpleInterval>) {
    match label {
        "trav-whole-contig" => (10, 5, 50, vec![interval("chr1", 1, CHR1)]),
        "trav-covered-only" => (10, 5, 50, vec![interval("chr1", 101, 140)]),
        "trav-two-intervals" => (
            10,
            5,
            50,
            vec![interval("chr1", 101, 140), interval("chr1", 201, 230)],
        ),
        "trav-merged-intervals" => (
            40,
            5,
            50,
            vec![interval("chr1", 101, 140), interval("chr1", 201, 230)],
        ),
        "trav-small-max" => (10, 5, 20, vec![interval("chr1", 1, CHR1)]),
        "trav-large-min" => (10, 400, 500, vec![interval("chr1", 101, 140)]),
        "trav-no-padding" => (0, 5, 50, vec![interval("chr1", 101, 140)]),
        "trav-no-reads" => (10, 5, 50, vec![interval("chr2", 100, 140)]),
        "trav-lone-read" => (10, 5, 50, vec![interval("chr1", 395, 425)]),
        other => panic!("unknown traversal case {other}"),
    }
}

fn row<'a>(text: &'a str, prefix: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(prefix))
        .unwrap_or_else(|| panic!("the golden carries no row {prefix:?}"))
}

fn join(intervals: &[SimpleInterval]) -> String {
    intervals
        .iter()
        .map(|i| format!("{}:{}-{}", i.contig, i.start, i.end))
        .collect::<Vec<_>>()
        .join(",")
}

/// A locus of the wrapped iterator, as a `Located`.
struct LocatedContext {
    contig: String,
    position: i32,
    depth: usize,
}

impl Located for LocatedContext {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.position
    }
    fn stop(&self) -> i32 {
        self.position
    }
}

#[test]
fn every_shard_matches_the_reference() {
    let text = golden();
    let header = header();

    for label in [
        "shard-one",
        "shard-adjacent",
        "shard-overlapping",
        "shard-unsorted",
        "shard-merged-by-padding",
        "shard-off-contig",
        "shard-zero-padding",
        "shard-negative-padding",
    ] {
        let (padding, intervals) = shard_case(label);
        match ReadShard::new(&intervals, padding, &header) {
            None => {
                assert!(
                    text.lines()
                        .any(|line| line.starts_with(&format!("error\t{label}\t"))),
                    "{label}: we refused a padding the reference accepted"
                );
            }
            Some(shard) => {
                assert_eq!(
                    format!(
                        "{}\t{}",
                        join(&shard.intervals),
                        join(&shard.padded_intervals)
                    ),
                    row(&text, &format!("shard\t{label}\t")),
                    "{label}"
                );
            }
        }
    }
}

#[test]
fn every_validation_matches_the_reference() {
    let text = golden();

    for (label, min, max, padding, max_reads, snp, indel) in ARGS {
        let args = AssemblyRegionArgs {
            min_assembly_region_size: *min,
            max_assembly_region_size: *max,
            assembly_region_padding: *padding,
            max_reads_per_alignment_start: *max_reads,
            snp_padding_for_genotyping: *snp,
            indel_padding_for_genotyping: *indel,
            ..Default::default()
        };
        let ours = match args.validate() {
            Ok(()) => "ok".to_string(),
            Err(error) => format!("E:{}:{}", error.class(), error.message()),
        };
        assert_eq!(ours, row(&text, &format!("args\t{label}\t")), "{label}");
    }
}

#[test]
fn every_traversal_matches_the_reference() {
    let text = golden();
    let header = header();
    let records = reads();

    for label in [
        "trav-whole-contig",
        "trav-covered-only",
        "trav-two-intervals",
        "trav-merged-intervals",
        "trav-small-max",
        "trav-large-min",
        "trav-no-padding",
        "trav-no-reads",
        "trav-lone-read",
    ] {
        let (padding, min_size, max_size, intervals) = traversal_case(label);
        let args = AssemblyRegionArgs {
            assembly_region_padding: padding,
            min_assembly_region_size: min_size,
            max_assembly_region_size: max_size,
            ..Default::default()
        };
        let shard = ReadShard::new(&intervals, padding, &header).expect("a valid padding");

        // The dump's shard keeps only the reads overlapping its padded intervals, which is what
        // querying a BAM over those intervals returns.
        let shard_reads: Vec<BamRecord> = records
            .iter()
            .filter(|record| {
                shard.padded_intervals.iter().any(|interval| {
                    interval.overlaps(
                        "chr1",
                        gatk_engine::read_utils::start(record),
                        gatk_engine::read_utils::end(record),
                    )
                })
            })
            .cloned()
            .collect();

        let samples = vec![Some("sampleA".to_string()), Some("sampleB".to_string())];
        let states = ReadStateManager::new(samples.clone(), DownsamplingInfo::NONE)
            .expect("no downsampling");
        let contexts = locus_iterator::contexts(
            &shard_reads,
            samples,
            &header,
            LocusIteratorOptions {
                include_deletions: true,
                include_ns: true,
            },
            states,
        )
        .unwrap_or_else(|e| panic!("{label}: {e:?}"));
        let located: Vec<LocatedContext> = contexts
            .iter()
            .map(|context| LocatedContext {
                contig: context.contig.clone(),
                position: context.position,
                depth: context.pileup.size(),
            })
            .collect();

        // The dump's probe evaluator: a locus is active when at least two reads cover it.
        let traversed = assembly_regions(
            &located,
            &shard_reads,
            &shard,
            &args,
            &header,
            &|_, context| match context {
                Some(context) if context.depth >= 2 => 1.0,
                _ => 0.0,
            },
        )
        .unwrap_or_else(|e| panic!("{label}: {e:?}"));

        let expected: usize = row(&text, &format!("count\t{label}\t"))
            .parse()
            .expect("a number");
        assert_eq!(traversed.len(), expected, "{label}: region count");

        for (index, entry) in traversed.iter().enumerate() {
            let region = &entry.region;
            assert_eq!(
                format!(
                    "{}:{}-{}|{}:{}-{}|{}|{}",
                    region.span().contig,
                    region.span().start,
                    region.span().end,
                    region.padded_span().contig,
                    region.padded_span().start,
                    region.padded_span().end,
                    region.is_active(),
                    region.size()
                ),
                row(&text, &format!("region\t{label}\t{index}\t")),
                "{label}, region {index}"
            );
            let names: Vec<String> = region
                .reads()
                .iter()
                .map(|read| read.read_name.clone())
                .collect();
            assert_eq!(
                names.join(","),
                row(&text, &format!("rread\t{label}\t{index}\t")),
                "{label}, reads of region {index}"
            );
        }
    }
}

/// Two intervals that touch are already one before any padding is applied, and two fifteen bases
/// apart are not.
#[test]
fn a_shard_merges_at_zero_padding() {
    let text = golden();
    assert_eq!(
        row(&text, "shard\tshard-adjacent\t"),
        "chr1:100-300\tchr1:90-310"
    );
    assert_eq!(
        row(&text, "shard\tshard-merged-by-padding\t"),
        "chr1:100-200,chr1:216-300\tchr1:90-310"
    );
    // And zero padding still merges: the two lists are equal and both are the merged one.
    assert_eq!(
        row(&text, "shard\tshard-zero-padding\t"),
        "chr1:100-300\tchr1:100-300"
    );
}

/// `groupIntervalsByContig` groups on a *change* of contig, so an unsorted list gives two groups
/// for one contig. Nothing upstream notices because the walker only hands it sorted intervals.
#[test]
fn grouping_by_contig_trusts_its_input_to_be_sorted() {
    let sorted = vec![
        interval("chr1", 1, 10),
        interval("chr1", 20, 30),
        interval("chr2", 1, 10),
    ];
    assert_eq!(group_intervals_by_contig(&sorted).len(), 2);

    let unsorted = vec![
        interval("chr1", 1, 10),
        interval("chr2", 1, 10),
        interval("chr1", 20, 30),
    ];
    assert_eq!(group_intervals_by_contig(&unsorted).len(), 3);
}
