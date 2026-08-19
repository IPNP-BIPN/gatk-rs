//! Conformance for `ASEReadCounter` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/ASEReadCounterDump.java`.
//!
//! # What this suite is for
//!
//!  * **the three overlap modes**, which disagree at the site where the two mates carry different
//!    bases and agree everywhere else;
//!  * **the ordered cascade**, so the improperly paired read is charged to `improperPairs` and to
//!    nothing else whatever its qualities;
//!  * **raw depth after the overlap filter**, which is why it is 4 where the pair merged and 5
//!    where it did not;
//!  * **the depth threshold removing lines rather than zeroing them**;
//!  * **and the header being written before the traversal**, so a run with no sites still has one
//!    line.
//!
//! The pileups are built here from the dump's own reads: the locus walker has its own suites, and
//! what this tool decides is which elements survive and which column each one lands in.

use gatk_corpus as corpus;
use gatk_engine::pileup::PileupElement;
use gatk_engine::read_pileup::ReadPileup;
use gatk_tools::ase_read_counter::{
    self, CountType, OutputFormat, SiteCounts, DEFAULT_READ_FILTERS,
};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/ase_read_counter.txt.gz"),
    )
}

/// One case's table, as rows with the escaping undone.
fn table(text: &str, label: &str) -> Vec<String> {
    let prefix = format!("table\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {label}"))
        .split("\\n")
        .filter(|row| !row.is_empty())
        .map(|row| row.replace("\\t", "\t"))
        .collect()
}

/// The dump's five reads, at chr1:10 with twenty matched bases.
fn read(
    name: &str,
    bases: &[u8],
    mapping_quality: u8,
    flags: u16,
    low_quality_at: Option<usize>,
) -> BamRecord {
    let mut qualities = vec![30u8; bases.len()];
    if let Some(index) = low_quality_at {
        qualities[index] = 20;
    }
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: 10,
        mapping_quality,
        cigar: Cigar {
            elements: vec![CigarElement {
                op: Op::M,
                length: bases.len() as u32,
            }],
        },
        read_bases: bases.to_vec(),
        base_qualities: qualities,
        mate_reference_index: 0,
        mate_alignment_start: 10,
        flags,
        ..Default::default()
    }
}

/// The fixture, in the order the BAM holds it.
fn reads() -> Vec<BamRecord> {
    vec![
        // The pair: proper, first and second of pair, agreeing at 12 and disagreeing at 16.
        read("p001", b"ACGTACGTACGTACGTACGT", 60, 0x1 | 0x2 | 0x40, None),
        read("p001", b"ACGTACATACGTACGTACGT", 60, 0x1 | 0x2 | 0x80, None),
        // The reference base at 12, and a mapping quality of 30.
        read("s001", b"ACTTACGTACGTACGTACGT", 30, 0, None),
        // A base quality of 20 at position 20, which is offset 10.
        read("s002", b"ACGTACGTACGTACGTACGT", 60, 0, Some(10)),
        // Paired and NOT properly paired.
        read("s003", b"ACGTACGTACGTACGTACGT", 60, 0x1 | 0x40, None),
    ]
}

/// The pileup at one position, in the order the reads were added.
fn pileup_at(records: &[BamRecord], position: i32) -> ReadPileup<'_> {
    let offset = position - 10;
    let elements: Vec<PileupElement> = records
        .iter()
        .map(|record| {
            PileupElement::for_read_and_offset(record, offset).expect("the read covers the locus")
        })
        .collect();
    ReadPileup::new("chr1", position, elements)
}

/// One site of the golden, counted the port's way.
fn counted(
    records: &[BamRecord],
    position: i32,
    count_type: CountType,
    minimum_mapping_quality: i32,
    minimum_base_quality: u8,
) -> SiteCounts {
    let pileup = pileup_at(records, position);
    let filtered = ase_read_counter::filter_pileup(&pileup, count_type);
    // The reference base at 12, 16 and 20 is `T`; the alternate the VCF names is `G`.
    ase_read_counter::count_site(
        &filtered,
        b'T',
        b'G',
        minimum_mapping_quality,
        minimum_base_quality,
    )
}

fn row(counts: SiteCounts, position: i32, id: &str) -> String {
    ase_read_counter::line(
        "chr1",
        position,
        id,
        b'T',
        b'G',
        counts,
        1,
        OutputFormat::RTable,
    )
    .expect("the site is deep enough")
}

#[test]
fn the_three_overlap_modes_match_the_golden() {
    let text = golden();
    let records = reads();

    for (label, count_type) in [
        ("default", CountType::CountFragmentsRequireSameBase),
        ("count-reads", CountType::CountReads),
        ("count-fragments", CountType::CountFragments),
    ] {
        let expected = table(&text, label);
        assert_eq!(expected[0], ase_read_counter::header(OutputFormat::RTable));
        let mine: Vec<String> = [(12, "rs1"), (16, "rs2"), (20, "rs3")]
            .iter()
            .map(|(position, id)| {
                row(
                    counted(&records, *position, count_type, 0, 0),
                    *position,
                    id,
                )
            })
            .collect();
        assert_eq!(mine, expected[1..], "{label}");
    }
}

/// The site where the mates disagree is the one that separates the three modes.
#[test]
fn a_disagreeing_pair_is_discarded_only_by_the_default() {
    let records = reads();
    let discarded = counted(&records, 16, CountType::CountFragmentsRequireSameBase, 0, 0);
    let kept_one = counted(&records, 16, CountType::CountFragments, 0, 0);
    let kept_both = counted(&records, 16, CountType::CountReads, 0, 0);

    assert_eq!(discarded.raw, 3, "the pair is gone");
    assert_eq!(kept_one.raw, 4, "one of the two survives");
    assert_eq!(kept_both.raw, 5, "both survive");
    // Only `COUNT_READS` sees the disagreeing base at all, and it is neither the reference nor the
    // alternate.
    assert_eq!(kept_both.other_bases, 1);
    assert_eq!(kept_one.other_bases, 0);
    assert_eq!(discarded.other_bases, 0);
}

#[test]
fn the_quality_thresholds_match_the_golden() {
    let text = golden();
    let records = reads();

    // `--min-mapping-quality 40` moves the read at 30 out of the reference count.
    let expected = table(&text, "min-mapq");
    let mine: Vec<String> = [(12, "rs1"), (16, "rs2"), (20, "rs3")]
        .iter()
        .map(|(position, id)| {
            row(
                counted(
                    &records,
                    *position,
                    CountType::CountFragmentsRequireSameBase,
                    40,
                    0,
                ),
                *position,
                id,
            )
        })
        .collect();
    assert_eq!(mine, expected[1..]);

    // `--min-base-quality 25` moves the one low base at position 20.
    let expected = table(&text, "min-baseq");
    let mine: Vec<String> = [(12, "rs1"), (16, "rs2"), (20, "rs3")]
        .iter()
        .map(|(position, id)| {
            row(
                counted(
                    &records,
                    *position,
                    CountType::CountFragmentsRequireSameBase,
                    0,
                    25,
                ),
                *position,
                id,
            )
        })
        .collect();
    assert_eq!(mine, expected[1..]);
}

/// The improperly paired read is charged to one column and no other, whatever else is true of it.
#[test]
fn the_cascade_charges_each_read_once() {
    let records = reads();
    let counts = counted(
        &records,
        12,
        CountType::CountFragmentsRequireSameBase,
        40,
        25,
    );
    assert_eq!(counts.improper_pairs, 1);
    // Four elements after the pair merged; one improper, one low mapping quality, two counted.
    assert_eq!(counts.raw, 4);
    assert_eq!(counts.low_mapping_quality, 1);
    assert_eq!(
        counts.improper_pairs
            + counts.low_mapping_quality
            + counts.low_base_quality
            + counts.other_bases
            + counts.total,
        counts.raw,
        "the five buckets partition the pileup"
    );
}

/// A site below the threshold produces no line at all.
#[test]
fn the_depth_threshold_removes_lines() {
    let text = golden();
    let records = reads();
    let expected = table(&text, "min-depth");
    // Three sites, and the one whose pair was discarded has a total of 2.
    let mine: Vec<String> = [(12, "rs1"), (16, "rs2"), (20, "rs3")]
        .iter()
        .filter_map(|(position, id)| {
            ase_read_counter::line(
                "chr1",
                *position,
                id,
                b'T',
                b'G',
                counted(
                    &records,
                    *position,
                    CountType::CountFragmentsRequireSameBase,
                    0,
                    0,
                ),
                3,
                OutputFormat::RTable,
            )
        })
        .collect();
    assert_eq!(mine.len(), 2, "the shallow site is gone");
    assert_eq!(mine, expected[1..]);
}

/// The format changes the separator and nothing else, and the header is written either way.
#[test]
fn the_format_is_a_separator() {
    let text = golden();
    assert_eq!(
        ase_read_counter::header(OutputFormat::Table),
        ase_read_counter::header(OutputFormat::RTable)
    );
    assert_eq!(table(&text, "table"), table(&text, "default"));
    let csv = table(&text, "csv");
    assert_eq!(csv[0], ase_read_counter::header(OutputFormat::Csv));
    // A window with no sites still writes the header line, and nothing else.
    assert_eq!(table(&text, "no-sites").len(), 1);
}

/// The eight filters are named, and the walker's own are not among them.
#[test]
fn the_default_filter_set_is_eight_and_excludes_wellformed() {
    assert_eq!(DEFAULT_READ_FILTERS.len(), 8);
    assert!(!DEFAULT_READ_FILTERS.contains(&"WellformedReadFilter"));
}
