//! Conformance for `CollectAllelicCounts` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/CollectAllelicCountsDump.java`.
//!
//! # What this suite is for
//!
//!  * **every locus is a row**, including the ones with no reads;
//!  * **except where the reference base is not `ACGT`**, which is a missing row rather than an
//!    empty one, and the only hole the table has;
//!  * **the alternate count is total minus reference**, so three different non-reference bases are
//!    all in it while the alternate base names only the commonest;
//!  * **the base quality threshold is the collector's**, so one base at quality 19 changes a count
//!    without changing which reads were read;
//!  * **and the header is a SAM header** with a read group the tool invents.

use gatk_corpus as corpus;
use gatk_engine::pileup::PileupElement;
use gatk_engine::read_pileup::ReadPileup;
use gatk_tools::collect_allelic_counts::{
    self, AllelicCount, ADDITIONAL_READ_FILTERS, DEFAULT_MINIMUM_BASE_QUALITY,
};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/collect_allelic_counts.txt.gz"),
    )
}

/// One case's file, unescaped.
fn file(text: &str, label: &str) -> String {
    let prefix = format!("table\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {label}"))
        .replace("\\t", "\t")
        .replace("\\n", "\n")
}

/// The dump's reads, at chr1:10 with six matched bases.
fn read(name: &str, bases: &[u8], low_quality_at: Option<usize>) -> BamRecord {
    let mut qualities = vec![30u8; bases.len()];
    if let Some(index) = low_quality_at {
        qualities[index] = 19;
    }
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: 10,
        mapping_quality: 60,
        cigar: Cigar {
            elements: vec![CigarElement {
                op: Op::M,
                length: bases.len() as u32,
            }],
        },
        read_bases: bases.to_vec(),
        base_qualities: qualities,
        mate_reference_index: -1,
        ..Default::default()
    }
}

fn reads() -> Vec<BamRecord> {
    vec![
        read("r001", b"ACGTAC", None),
        read("r002", b"ACGTAC", None),
        read("r003", b"AAAAAA", None),
        read("r004", b"CCCCCC", None),
        read("r005", b"NNNNNN", None),
        read("r006", b"ACGTAC", Some(0)),
    ]
}

/// The reference bases of chr1, which repeat `ACGT`.
fn reference_base(position: i32) -> u8 {
    b"ACGT"[((position - 1) % 4) as usize]
}

/// The whole table for one window, at one base quality threshold.
fn table(records: &[BamRecord], window: std::ops::RangeInclusive<i32>, minimum: u8) -> String {
    let counts: Vec<AllelicCount> = window
        .filter_map(|position| {
            let offset = position - 10;
            let elements: Vec<PileupElement> = records
                .iter()
                .filter_map(|record| PileupElement::for_read_and_offset(record, offset))
                .collect();
            let pileup = ReadPileup::new("chr1", position, elements);
            collect_allelic_counts::collect_at_locus(
                reference_base(position),
                &pileup,
                "chr1",
                position,
                minimum,
            )
        })
        .collect();
    collect_allelic_counts::write(
        &[("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        "NA1",
        &counts,
    )
}

#[test]
fn the_default_table_matches_the_golden() {
    let text = golden();
    let records = reads();
    assert_eq!(
        table(&records, 10..=19, DEFAULT_MINIMUM_BASE_QUALITY),
        file(&text, "default")
    );
}

/// One base at quality 19 is the whole difference between the two thresholds.
#[test]
fn the_base_quality_threshold_is_the_collectors() {
    let text = golden();
    let records = reads();
    assert_eq!(
        table(&records, 10..=19, 0),
        file(&text, "base-quality-zero")
    );
    assert_eq!(
        table(&records, 10..=19, 30),
        file(&text, "base-quality-thirty")
    );
    // The read is present either way; only its base is dropped.
    let default = file(&text, "default");
    let zero = file(&text, "base-quality-zero");
    assert_ne!(default, zero);
    assert!(default.contains("chr1\t10\t1\t3\tC\tA"));
    assert!(zero.contains("chr1\t10\t1\t4\tC\tA"));
}

/// A window with no reads is all zeroes and an `N` alternate, not a missing table.
#[test]
fn a_window_with_no_reads_is_rows_of_zeroes() {
    let text = golden();
    assert_eq!(
        table(&[], 60..=64, DEFAULT_MINIMUM_BASE_QUALITY),
        file(&text, "no-reads")
    );
}

/// The reference's `N` run is the one hole in the table.
#[test]
fn an_n_in_the_reference_skips_the_locus() {
    let text = golden();
    let rows = file(&text, "reference-n");
    // 123, 124, 129 and 130 are present; 125 to 128 are the `N` run and are absent.
    for present in [123, 124, 129, 130] {
        assert!(rows.contains(&format!("chr1\t{present}\t")), "{present}");
    }
    for absent in [125, 126, 127, 128] {
        assert!(!rows.contains(&format!("chr1\t{absent}\t")), "{absent}");
    }

    // And the port refuses the same locus for the same reason.
    let pileup = ReadPileup::new("chr1", 125, Vec::new());
    assert!(collect_allelic_counts::collect_at_locus(
        b'N',
        &pileup,
        "chr1",
        125,
        DEFAULT_MINIMUM_BASE_QUALITY
    )
    .is_none());
}

/// The alternate count is every non-reference base, and the alternate base is only the commonest.
#[test]
fn the_alternate_count_is_total_minus_reference() {
    let records = reads();
    // Position 11's reference is `G`; the pileup holds four `C` and two `A`, all non-reference.
    let offset = 1;
    let elements: Vec<PileupElement> = records
        .iter()
        .filter_map(|record| PileupElement::for_read_and_offset(record, offset))
        .collect();
    let pileup = ReadPileup::new("chr1", 11, elements);
    let count = collect_allelic_counts::collect_at_locus(
        b'G',
        &pileup,
        "chr1",
        11,
        DEFAULT_MINIMUM_BASE_QUALITY,
    )
    .expect("a row");
    assert_eq!(count.reference_count, 0);
    assert_eq!(count.alternate_count, 5, "every ACGT read that is not G");
    assert_eq!(count.alternate_nucleotide, b'C', "the commonest of them");
}

/// The four filters this tool adds are named, with the mapping quality threshold among them.
#[test]
fn the_tool_adds_four_filters() {
    assert_eq!(ADDITIONAL_READ_FILTERS.len(), 4);
    assert_eq!(ADDITIONAL_READ_FILTERS[3], "MappingQualityReadFilter");
    assert_eq!(collect_allelic_counts::DEFAULT_MINIMUM_MAPPING_QUALITY, 30);
    // The window covered only by the read at mapping quality 20 has rows, and they are empty.
    let rows = file(&golden(), "low-mapping-quality");
    assert!(rows.contains("chr1\t30\t0\t0\t"));
}
