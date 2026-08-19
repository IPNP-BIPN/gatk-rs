//! Conformance for the `Pileup` tool against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/PileupDump.java`.
//!
//! # What this suite is for
//!
//!  * **the line's shape**, which is five space-separated fields and then whatever the optional
//!    columns add, each preceded by a space of its own;
//!  * **the trailing space** a line with no metadata ends in, which is the features column being
//!    printed empty rather than skipped;
//!  * **the deletion filter running first**, so a locus covered only by a deletion prints an empty
//!    base string and an empty quality string, and the verbose deletion count is zero;
//!  * **and the reference base being `N`** when the run has none.
//!
//! The pileups are built here from the fixture's reads: the locus walker and the pileup itself have
//! their own suites, and what this tool decides is how one line is written.

use gatk_corpus as corpus;
use gatk_engine::pileup::PileupElement;
use gatk_engine::read_pileup::ReadPileup;
use gatk_tools::pileup;
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/pileup_tool.txt.gz"),
    )
}

fn row(text: &str, label: &str) -> String {
    let prefix = format!("pileup\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {label}"))
        .to_string()
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
}

/// The fixture's `r001`: ten bases at chr1:10, all matched.
fn read(name: &str, start: i32, cigar: Cigar, bases: &[u8]) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality: 60,
        cigar,
        read_bases: bases.to_vec(),
        base_qualities: vec![30; bases.len()],
        mate_reference_index: -1,
        ..Default::default()
    }
}

fn matched(length: u32) -> Cigar {
    Cigar {
        elements: vec![CigarElement { op: Op::M, length }],
    }
}

/// The reference bases of chr1 in `ReadWalkerDump.FASTA`, which repeat `ACGT`.
fn reference_base(position: i32) -> char {
    b"ACGT"[((position - 1) % 4) as usize] as char
}

#[test]
fn the_plain_line_matches_the_golden() {
    let text = golden();
    let record = read("r001", 10, matched(10), b"ACGTACGTAC");

    let mut written = String::new();
    for position in 10..=19 {
        let offset = position - 10;
        let element = PileupElement::for_read_and_offset(&record, offset).expect("an element");
        let pileup = ReadPileup::new("chr1", position, vec![element]);
        written.push_str(&pileup::line(
            &pileup,
            reference_base(position),
            &[],
            false,
            false,
        ));
    }
    assert_eq!(escape(&written), row(&text, "plain"));
}

#[test]
fn both_optional_columns_match_the_golden() {
    let text = golden();
    let record = read("r001", 10, matched(10), b"ACGTACGTAC");

    let mut written = String::new();
    for position in 10..=19 {
        let offset = position - 10;
        let element = PileupElement::for_read_and_offset(&record, offset).expect("an element");
        let pileup = ReadPileup::new("chr1", position, vec![element]);
        written.push_str(&pileup::line(
            &pileup,
            reference_base(position),
            &[],
            true,
            true,
        ));
    }
    assert_eq!(escape(&written), row(&text, "both"));
}

/// Without a reference every middle field is `N`, and nothing else moves.
#[test]
fn a_run_without_a_reference_writes_n() {
    let text = golden();
    let record = read("r001", 10, matched(10), b"ACGTACGTAC");

    let mut written = String::new();
    for position in 10..=19 {
        let offset = position - 10;
        let element = PileupElement::for_read_and_offset(&record, offset).expect("an element");
        let pileup = ReadPileup::new("chr1", position, vec![element]);
        written.push_str(&pileup::line(&pileup, 'N', &[], false, false));
    }
    assert_eq!(escape(&written), row(&text, "no-reference"));
}

/// A locus covered only by a deletion prints an empty base string, and the verbose count is zero
/// because the deletions were filtered out before it was taken.
#[test]
fn a_deleted_locus_prints_nothing_and_counts_nothing() {
    let deletion = read(
        "r004",
        140,
        Cigar {
            elements: vec![
                CigarElement {
                    op: Op::M,
                    length: 5,
                },
                CigarElement {
                    op: Op::D,
                    length: 10,
                },
                CigarElement {
                    op: Op::M,
                    length: 5,
                },
            ],
        },
        b"ACGTACGTAC",
    );
    // Position 145 is inside the deletion.
    let element = PileupElement::for_read_and_offset(&deletion, 4).expect("an element");
    assert!(!element.is_deletion(), "offset 4 is the last matched base");

    let deleted = ReadPileup::new("chr1", 145, vec![deleted_element(&deletion)]);
    let line = pileup::line(&deleted, 'A', &[], false, true);
    assert_eq!(line, "chr1 145 A    0 \n");
    assert!(
        golden().contains("chr1 145 A    0 "),
        "the golden holds the same line"
    );
}

/// A pileup element sitting on the deletion itself.
fn deleted_element(read: &BamRecord) -> PileupElement<'_> {
    PileupElement::new(
        read,
        4,
        CigarElement {
            op: Op::D,
            length: 10,
        },
        1,
        0,
    )
}

/// The features column is printed either way, which is what leaves the trailing space.
#[test]
fn the_features_column_is_printed_empty() {
    assert_eq!(pileup::features_string(&[]), "");
    assert_eq!(
        pileup::features_string(&["one".to_string(), "two".to_string()]),
        "[Feature(s): one, two]"
    );
    let record = read("r001", 10, matched(10), b"ACGTACGTAC");
    let element = PileupElement::for_read_and_offset(&record, 0).expect("an element");
    let pileup = ReadPileup::new("chr1", 10, vec![element]);
    assert!(
        pileup::line(&pileup, 'C', &[], false, false).ends_with(" \n"),
        "a line with no metadata ends in a space"
    );
}

/// The three filters this tool adds are named, so dropping one fails here.
#[test]
fn the_tool_adds_three_filters() {
    assert_eq!(pileup::ADDITIONAL_READ_FILTERS.len(), 3);
    assert_eq!(pileup::ADDITIONAL_READ_FILTERS[0], "NotDuplicateReadFilter");
    // The duplicate window produced an empty file, which is the filter working.
    assert_eq!(row(&golden(), "duplicate"), "");
}
