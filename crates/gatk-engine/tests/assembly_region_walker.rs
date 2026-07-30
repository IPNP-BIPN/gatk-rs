//! Conformance for `AssemblyRegionWalker` against the oracle.
//!
//! Goldens from `tools/readfilter-conformance/AssemblyRegionWalkerDump.java`, a real
//! `AssemblyRegionWalker` run through the real command line, so its defaults are measured rather
//! than transcribed.
//!
//! The pair of runs that settles what `--force-active` does:
//!
//! ```text
//! apply  force-active-small-regions          0  chr1:1-20|chr1:1-120|true|3|120
//! apply  threshold-above-all-small-regions   0  chr1:1-20|chr1:1-120|false|3|120
//! count  force-active-small-regions          10
//! count  threshold-above-all-small-regions   10
//! ```
//!
//! Same threshold, same region sizes, same ten regions with the same ten boundaries. Only the flag
//! differs. `--force-active` rewrites `isActive` **after** the regions have been cut, so it changes
//! what every region claims to be without changing where any of them is. A port that folded it into
//! the evaluator would return one region covering everything.

use gatk_corpus as corpus;
use gatk_engine::assembly_region_iterator::AssemblyRegionArgs;
use gatk_engine::assembly_region_walker::{traverse, LocusDepth, WalkerError};
use gatk_engine::interval::SimpleInterval;
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const CONTIG_LENGTH: i32 = 200;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/assembly_region_walker.txt.gz"),
    )
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    for contig in ["chr1", "chr2"] {
        header
            .sequences
            .push(SequenceRecord::new(contig, CONTIG_LENGTH));
    }
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "sample1");
    group.attributes.set("PL", "ILLUMINA");
    header.read_groups.push(group);
    header
}

fn interval(contig: &str, start: i32, end: i32) -> SimpleInterval {
    SimpleInterval::new(contig, start, end).expect("a valid interval")
}

fn read(name: &str, contig: &str, start: i32, cigar: &str, with_group: bool) -> BamRecord {
    let cigar = htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar");
    let length = cigar.read_length().max(10) as usize;
    let mut tags = htsjdk_bam::tag::Tags::new();
    if with_group {
        tags.insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
    }
    BamRecord {
        read_name: name.to_string(),
        reference_index: if contig == "chr1" { 0 } else { 1 },
        alignment_start: start,
        mapping_quality: 60,
        read_bases: (0..length).map(|i| b"ACGT"[i % 4]).collect(),
        base_qualities: vec![30; length],
        cigar,
        tags,
        ..Default::default()
    }
}

/// The read-walker fixture **after** the two default filters have run, which is what the shard's
/// iterator hands the traversal. `r005` has no read group and `r006` carries an `N` operator, both
/// of which `WellformedReadFilter` removes; `m001` is mapped with an empty cigar and goes the same
/// way; `u001` and `x001` are unmapped and `MappedReadFilter` removes them.
///
/// The filtering itself is oracle-backed in its own suite, so reproducing it here would be testing
/// the filters twice and the traversal not at all.
fn filtered_reads() -> Vec<BamRecord> {
    vec![
        read("r001", "chr1", 10, "10M", true),
        read("r002", "chr1", 65, "10M", true),
        read("r003", "chr1", 120, "10M", true),
        read("r004", "chr1", 140, "5M10D5M", true),
        read("r007", "chr1", 170, "10M", true),
        read("r101", "chr2", 10, "10M", true),
    ]
}

/// The same fixture with the default filters disabled, which is the `no-filters` run: `r005` has
/// no read group, so it has no sample, and the sample partition refuses it.
fn unfiltered_reads() -> Vec<BamRecord> {
    let mut reads = filtered_reads();
    reads.insert(4, read("r005", "chr1", 150, "10M", false));
    reads.sort_by_key(|record| (record.reference_index, record.alignment_start));
    reads
}

/// Label, intervals (`None` for the whole reference), and the arguments that differ from default.
struct Case {
    label: &'static str,
    intervals: Option<Vec<SimpleInterval>>,
    args: AssemblyRegionArgs,
    filtered: bool,
}

fn cases() -> Vec<Case> {
    let whole = || {
        vec![
            interval("chr1", 1, CONTIG_LENGTH),
            interval("chr2", 1, CONTIG_LENGTH),
        ]
    };
    let small = || AssemblyRegionArgs {
        min_assembly_region_size: 5,
        max_assembly_region_size: 20,
        ..Default::default()
    };
    vec![
        Case {
            label: "all",
            intervals: Some(whole()),
            args: Default::default(),
            filtered: true,
        },
        Case {
            label: "chr1",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: Default::default(),
            filtered: true,
        },
        Case {
            label: "chr2",
            intervals: Some(vec![interval("chr2", 1, CONTIG_LENGTH)]),
            args: Default::default(),
            filtered: true,
        },
        Case {
            label: "two-on-one-contig",
            intervals: Some(vec![interval("chr1", 10, 40), interval("chr1", 150, 190)]),
            args: Default::default(),
            filtered: true,
        },
        Case {
            label: "two-contigs",
            intervals: Some(vec![interval("chr1", 10, 40), interval("chr2", 10, 40)]),
            args: Default::default(),
            filtered: true,
        },
        Case {
            label: "narrow",
            intervals: Some(vec![interval("chr1", 100, 110)]),
            args: Default::default(),
            filtered: true,
        },
        Case {
            label: "small-regions",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: small(),
            filtered: true,
        },
        Case {
            label: "zero-padding",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: AssemblyRegionArgs {
                assembly_region_padding: 0,
                ..Default::default()
            },
            filtered: true,
        },
        Case {
            label: "large-padding",
            intervals: Some(vec![interval("chr1", 100, 110)]),
            args: AssemblyRegionArgs {
                assembly_region_padding: 500,
                ..Default::default()
            },
            filtered: true,
        },
        Case {
            label: "force-active",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: AssemblyRegionArgs {
                force_active: true,
                ..Default::default()
            },
            filtered: true,
        },
        Case {
            label: "threshold-above-all",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: AssemblyRegionArgs {
                active_prob_threshold: 2.0,
                ..Default::default()
            },
            filtered: true,
        },
        Case {
            label: "force-active-above-threshold",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: AssemblyRegionArgs {
                active_prob_threshold: 2.0,
                force_active: true,
                ..Default::default()
            },
            filtered: true,
        },
        Case {
            label: "force-active-small-regions",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: AssemblyRegionArgs {
                active_prob_threshold: 2.0,
                force_active: true,
                ..small()
            },
            filtered: true,
        },
        Case {
            label: "threshold-above-all-small-regions",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: AssemblyRegionArgs {
                active_prob_threshold: 2.0,
                ..small()
            },
            filtered: true,
        },
        Case {
            label: "zero-propagation",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: AssemblyRegionArgs {
                max_prob_propagation_distance: 0,
                ..Default::default()
            },
            filtered: true,
        },
        Case {
            label: "max-starts-1",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: AssemblyRegionArgs {
                max_reads_per_alignment_start: 1,
                ..Default::default()
            },
            filtered: true,
        },
        Case {
            label: "max-starts-0",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: AssemblyRegionArgs {
                max_reads_per_alignment_start: 0,
                ..Default::default()
            },
            filtered: true,
        },
        Case {
            label: "no-filters",
            intervals: Some(vec![interval("chr1", 1, CONTIG_LENGTH)]),
            args: Default::default(),
            filtered: false,
        },
    ]
}

fn row<'a>(text: &'a str, prefix: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(prefix))
        .unwrap_or_else(|| panic!("the golden carries no row {prefix:?}"))
}

/// The probe evaluator: a locus is active when at least one read covers it.
fn is_active(locus: &LocusDepth) -> f64 {
    if locus.depth >= 1 {
        1.0
    } else {
        0.0
    }
}

#[test]
fn every_traversal_matches_the_reference() {
    let text = golden();
    let header = header();
    let samples = vec![Some("sample1".to_string())];

    for case in cases() {
        let reads = if case.filtered {
            filtered_reads()
        } else {
            unfiltered_reads()
        };
        let intervals = case.intervals.clone().expect("every case names its intervals");

        let result = traverse(&reads, &intervals, &samples, &case.args, &header, &is_active);
        let label = case.label;

        let expected_summary = row(&text, &format!("summary\t{label}\t"));
        match &result {
            Err(error) => {
                assert_eq!(
                    format!("E:{}:{}", error.class(), error.message()),
                    expected_summary,
                    "{label}: the refusal"
                );
                assert_eq!(row(&text, &format!("count\t{label}\t")), "0", "{label}");
                continue;
            }
            Ok(_) => assert_eq!(expected_summary, "ok", "{label}: the reference refused"),
        }

        let regions = result.expect("just matched Ok");
        let expected_count: usize = row(&text, &format!("count\t{label}\t"))
            .parse()
            .expect("a number");
        assert_eq!(regions.len(), expected_count, "{label}: apply calls");

        for (index, entry) in regions.iter().enumerate() {
            let region = &entry.region;
            // The reference bases the walker hands `apply` are the padded span's, clipped to the
            // contig, so their count is the padded span's size.
            let reference_bases = region.padded_span().size();
            assert_eq!(
                format!(
                    "{}:{}-{}|{}:{}-{}|{}|{}|{}",
                    region.span().contig,
                    region.span().start,
                    region.span().end,
                    region.padded_span().contig,
                    region.padded_span().start,
                    region.padded_span().end,
                    region.is_active(),
                    region.size(),
                    reference_bases
                ),
                row(&text, &format!("apply\t{label}\t{index}\t")),
                "{label}, apply {index}"
            );
            let names: Vec<String> = region
                .reads()
                .iter()
                .map(|read| read.read_name.clone())
                .collect();
            assert_eq!(
                names.join(","),
                row(&text, &format!("aread\t{label}\t{index}\t")),
                "{label}, reads of region {index}"
            );
        }
    }
}

/// `--force-active` changes every flag and not one boundary.
#[test]
fn forcing_active_moves_no_boundary() {
    let text = golden();
    let spans = |label: &str| -> Vec<String> {
        text.lines()
            .filter_map(|line| line.strip_prefix(&format!("apply\t{label}\t")))
            .map(|rest| {
                let value = rest.split_once('\t').expect("an index and a value").1;
                // span|paddedSpan, dropping the flag and the counts.
                value
                    .split('|')
                    .take(2)
                    .collect::<Vec<_>>()
                    .join("|")
            })
            .collect()
    };
    let flags = |label: &str| -> Vec<String> {
        text.lines()
            .filter_map(|line| line.strip_prefix(&format!("apply\t{label}\t")))
            .map(|rest| {
                rest.split('|').nth(2).expect("a flag").to_string()
            })
            .collect()
    };

    let forced = spans("force-active-small-regions");
    let unforced = spans("threshold-above-all-small-regions");
    assert_eq!(forced.len(), 10, "the golden carries ten regions");
    assert_eq!(forced, unforced, "the boundaries moved");
    assert!(flags("force-active-small-regions")
        .iter()
        .all(|flag| flag == "true"));
    assert!(flags("threshold-above-all-small-regions")
        .iter()
        .all(|flag| flag == "false"));
}

/// Disabling the default read filters is refused, and by the sample partition rather than by the
/// walker: the read that gets through has no read group and therefore no sample.
#[test]
fn the_default_filters_are_load_bearing() {
    let text = golden();
    let summary = row(&text, "summary\tno-filters\t");
    assert!(
        summary.starts_with("E:java.lang.IllegalStateException:"),
        "{summary}"
    );
    assert!(summary.contains("SamplePartitioner"), "{summary}");
    // And the port refuses in the same place.
    let error = traverse(
        &unfiltered_reads(),
        &[interval("chr1", 1, CONTIG_LENGTH)],
        &[Some("sample1".to_string())],
        &Default::default(),
        &header(),
        &is_active,
    )
    .expect_err("the sample partition refuses a read with no sample");
    assert!(matches!(error, WalkerError::ReadState(_)));
}
