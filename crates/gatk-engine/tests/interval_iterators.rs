//! Conformance for the iterators between a shard of reads and an activity profile.
//!
//! Goldens from `tools/readfilter-conformance/IntervalIteratorsDump.java`.
//!
//! The suite's reason for existing is one comparison. `LocusIteratorByState` emits nothing for a
//! locus no read covers, and `IntervalAlignmentContextIterator` emits **something** there: an
//! `AlignmentContext` with an empty pileup, manufactured on the spot. The activity profile needs
//! that zero in order to close a region over the gap, so a port that skipped the wrapper would move
//! every region boundary downstream of a coverage hole.

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::locus_iterator::{self, AlignmentContext, LocusIteratorOptions};
use gatk_engine::locus_shards::{interval_alignment_contexts, interval_loci, sharded_intervals};
use gatk_engine::read_states::{DownsamplingInfo, ReadStateManager};
use gatk_engine::variant_source::Located;
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const CHR1: i32 = 300;
const CHR2: i32 = 200;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/interval_iterators.txt.gz"),
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

/// The dump's fixture: nothing before 101, a hole at 115-119, nothing after 140. The hole is
/// 115-119 and not 111-119 because r2 is 10M at 105, so it covers through 114: coverage is the
/// union of the reads, not of their starts, and the first golden is what said so.
fn reads() -> Vec<BamRecord> {
    vec![
        read("r1", "rg1", "10M", 101),
        read("r2", "rg1", "10M", 105),
        read("r3", "rg1", "10M", 120),
        read("r4", "rg2", "10M", 131),
    ]
}

/// Label, shard size, intervals: the `shard` cases in the dump's order.
fn shard_case(label: &str) -> (i32, Vec<SimpleInterval>) {
    match label {
        "shard-10-1" => (1, vec![interval("chr1", 10, 20)]),
        "shard-10-2" => (2, vec![interval("chr1", 10, 20)]),
        "shard-10-3" => (3, vec![interval("chr1", 10, 20)]),
        "shard-10-4" => (4, vec![interval("chr1", 10, 20)]),
        "shard-10-11" => (11, vec![interval("chr1", 10, 20)]),
        "shard-10-100" => (100, vec![interval("chr1", 10, 20)]),
        "shard-one-base" => (3, vec![interval("chr1", 10, 10)]),
        "shard-two-intervals" => (4, vec![interval("chr1", 10, 20), interval("chr2", 5, 9)]),
        "shard-zero" => (0, vec![interval("chr1", 10, 20)]),
        "shard-negative" => (-1, vec![interval("chr1", 10, 20)]),
        "shard-empty" => (5, vec![]),
        other => panic!("unknown shard case {other}"),
    }
}

fn locus_case(label: &str) -> Vec<SimpleInterval> {
    match label {
        "loci-one" => vec![interval("chr1", 10, 14)],
        "loci-two-contigs" => vec![interval("chr1", 10, 12), interval("chr2", 5, 6)],
        "loci-adjacent" => vec![interval("chr1", 10, 12), interval("chr1", 13, 14)],
        "loci-empty" => vec![],
        other => panic!("unknown locus case {other}"),
    }
}

fn context_case(label: &str) -> Vec<SimpleInterval> {
    match label {
        "ctx-covered" => vec![interval("chr1", 101, 110)],
        "ctx-leading-gap" => vec![interval("chr1", 95, 106)],
        "ctx-trailing-gap" => vec![interval("chr1", 135, 145)],
        "ctx-interior-gap" => vec![interval("chr1", 108, 122)],
        "ctx-all-gap" => vec![interval("chr1", 115, 119)],
        "ctx-two-intervals" => vec![interval("chr1", 101, 103), interval("chr1", 130, 133)],
        "ctx-other-contig" => vec![interval("chr2", 10, 14)],
        "ctx-both-contigs" => vec![interval("chr1", 138, 141), interval("chr2", 10, 12)],
        "ctx-whole-contig" => vec![interval("chr1", 1, CHR1)],
        "ctx-empty" => vec![],
        other => panic!("unknown context case {other}"),
    }
}

/// The `AlignmentContext` of a single locus, as a `Located`: one base wide.
struct LocatedContext {
    contig: String,
    position: i32,
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

fn rows<'a>(text: &'a str, kind: &str, label: &str) -> Vec<&'a str> {
    text.lines()
        .filter_map(|line| line.strip_prefix(&format!("{kind}\t{label}\t")))
        .map(|rest| rest.split_once('\t').expect("an index and a value").1)
        .collect()
}

fn count(text: &str, label: &str) -> Option<usize> {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("count\t{label}\t")))
        .map(|value| value.parse().expect("a number"))
}

fn error<'a>(text: &'a str, label: &str) -> Option<&'a str> {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
}

#[test]
fn every_shard_matches_the_reference() {
    let text = golden();
    let labels: Vec<&str> = text
        .lines()
        .filter_map(|line| {
            line.strip_prefix("shard\t")
                .or(line.strip_prefix("error\t"))
        })
        .filter_map(|rest| rest.split('\t').next())
        .filter(|label| label.starts_with("shard-"))
        .collect();
    let mut seen: Vec<&str> = Vec::new();
    for label in labels {
        if seen.contains(&label) {
            continue;
        }
        seen.push(label);
    }
    assert!(!seen.is_empty(), "the golden carries no shard rows");

    for label in &seen {
        let (size, intervals) = shard_case(label);
        match sharded_intervals(&intervals, size) {
            None => {
                // A refused shard size is an `error` row and no `count` row.
                assert!(
                    error(&text, label).is_some(),
                    "{label}: we refused a shard size the reference accepted"
                );
                assert!(count(&text, label).is_none(), "{label}");
            }
            Some(shards) => {
                assert_eq!(
                    shards.len(),
                    count(&text, label).unwrap_or_else(|| panic!("{label}: no count row")),
                    "{label}: shard count"
                );
                let expected = rows(&text, "shard", label);
                for (index, shard) in shards.iter().enumerate() {
                    assert_eq!(
                        format!("{}:{}-{}", shard.contig, shard.start, shard.end),
                        expected[index],
                        "{label}, shard {index}"
                    );
                }
            }
        }
    }
    println!("{} shard cases identical", seen.len());
}

#[test]
fn every_locus_matches_the_reference() {
    let text = golden();
    for label in [
        "loci-one",
        "loci-two-contigs",
        "loci-adjacent",
        "loci-empty",
    ] {
        let loci = interval_loci(&locus_case(label));
        assert_eq!(
            loci.len(),
            count(&text, label).unwrap_or_else(|| panic!("{label}: no count row")),
            "{label}: locus count"
        );
        let expected = rows(&text, "locus", label);
        for (index, locus) in loci.iter().enumerate() {
            assert_eq!(
                format!("{}:{}-{}", locus.contig, locus.start, locus.end),
                expected[index],
                "{label}, locus {index}"
            );
        }
    }
}

#[test]
fn every_manufactured_context_matches_the_reference() {
    let text = golden();
    let header = header();
    let records = reads();

    for label in [
        "ctx-covered",
        "ctx-leading-gap",
        "ctx-trailing-gap",
        "ctx-interior-gap",
        "ctx-all-gap",
        "ctx-two-intervals",
        "ctx-other-contig",
        "ctx-both-contigs",
        "ctx-whole-contig",
        "ctx-empty",
    ] {
        let samples = vec![Some("sampleA".to_string()), Some("sampleB".to_string())];
        let states = ReadStateManager::new(samples.clone(), DownsamplingInfo::NONE)
            .expect("no downsampling");
        // `AssemblyRegionIterator` builds its LocusIteratorByState with includeDeletions and
        // includeNs both **true**, which is the `keepUniqueReadListInLibs` constructor's default
        // and not the LocusWalker's.
        let contexts: Vec<AlignmentContext> = locus_iterator::contexts(
            &records,
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
            })
            .collect();

        let emitted = interval_alignment_contexts(&located, &context_case(label), &header);
        assert_eq!(
            emitted.len(),
            count(&text, label).unwrap_or_else(|| panic!("{label}: no count row")),
            "{label}: context count"
        );

        let expected = rows(&text, "ctx", label);
        for (index, locus) in emitted.iter().enumerate() {
            let (size, names) = match locus.context {
                None => (0, String::new()),
                Some(context) => (
                    contexts[context].pileup.size(),
                    contexts[context]
                        .pileup
                        .elements
                        .iter()
                        .map(|element| element.read.read_name.clone())
                        .collect::<Vec<_>>()
                        .join(","),
                ),
            };
            assert_eq!(
                format!(
                    "{}:{}\t{}\t{}",
                    locus.interval.contig, locus.interval.start, size, names
                ),
                expected[index],
                "{label}, context {index}"
            );
        }
    }
}

/// The row that justifies the wrapper: a locus inside the coverage hole is emitted, with a pileup
/// of zero, rather than skipped.
#[test]
fn an_uncovered_locus_is_emitted_with_an_empty_pileup() {
    let text = golden();
    // The fixture's hole is 115-119 and `ctx-all-gap` asks for exactly it: five loci, none covered.
    assert_eq!(count(&text, "ctx-all-gap"), Some(5));
    for row in rows(&text, "ctx", "ctx-all-gap") {
        let (_, rest) = row.split_once('\t').expect("a position and a size");
        assert!(rest.starts_with('0'), "{row}");
    }
    // And `LocusIteratorByState` on its own emits nothing there, which is the comparison: the
    // covered case has as many contexts as it has loci, the gap case has none of its own.
    assert_eq!(count(&text, "ctx-covered"), Some(10));
}
