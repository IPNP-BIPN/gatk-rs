//! Conformance for `GetNormalArtifactData` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/GetNormalArtifactDataDump.java`.
//!
//! # What this suite is for
//!
//!  * **the seeded draw**, which is only observable when the keep probability is at its floor: of
//!    thirty eligible loci the reference keeps ONE, and the port keeps the same one;
//!  * **the two rejection rules**, one before the draw and one after, so a locus rejected by the
//!    second still consumed a number;
//!  * **the allele being chosen in the normal and counted in the tumour**;
//!  * **and the error probability moving the keep probability and nothing else.**
//!
//! The pileups are built here: the locus walker has its own suites, and what this tool decides is
//! which loci become rows.

use gatk_corpus as corpus;
use gatk_engine::java_random::JavaRandom;
use gatk_engine::pileup::PileupElement;
use gatk_engine::read_pileup::ReadPileup;
use gatk_tools::get_normal_artifact_data::{
    self, NormalArtifactRecord, Outcome, DEFAULT_ERROR_PROBABILITY, MIN_READ_LENGTH,
    STANDARD_MUTECT2_READ_FILTERS,
};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/normal_artifact_data.txt.gz"),
    )
}

fn file(text: &str, label: &str) -> String {
    let prefix = format!("table\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {label}"))
        .replace("\\t", "\t")
        .replace("\\n", "\n")
}

/// A forty-base read at chr1:10.
fn read(name: &str, bases: &[u8]) -> BamRecord {
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
        base_qualities: vec![30; bases.len()],
        mate_reference_index: -1,
        ..Default::default()
    }
}

/// The fixture's two sequences: one matching the reference from position 10, one all `A`.
fn matching() -> Vec<u8> {
    "CGTA".repeat(10).into_bytes()
}

fn alternate() -> Vec<u8> {
    vec![b'A'; 40]
}

/// The reference bases of chr1, which repeat `ACGT`.
fn reference_base(position: i32) -> u8 {
    b"ACGT"[((position - 1) % 4) as usize]
}

/// The normal: nine matching reads and one carrying the alternate.
fn normal_reads() -> Vec<BamRecord> {
    let mut reads: Vec<BamRecord> = (0..9)
        .map(|index| read(&format!("n{index}"), &matching()))
        .collect();
    reads.push(read("n9", &alternate()));
    reads
}

/// The tumour of the first fixture: seven matching, three carrying the alternate.
fn tumor_reads() -> Vec<BamRecord> {
    let mut reads: Vec<BamRecord> = (0..7)
        .map(|index| read(&format!("t{index}"), &matching()))
        .collect();
    reads.extend((7..10).map(|index| read(&format!("t{index}"), &alternate())));
    reads
}

/// The tumour of the second fixture: ten matching reads and no alternate at all.
fn matching_tumor() -> Vec<BamRecord> {
    (0..10)
        .map(|index| read(&format!("t{index}"), &matching()))
        .collect()
}

fn pileup_at(records: &[BamRecord], position: i32) -> ReadPileup<'_> {
    let offset = position - 10;
    ReadPileup::new(
        "chr1",
        position,
        records
            .iter()
            .filter_map(|record| PileupElement::for_read_and_offset(record, offset))
            .collect(),
    )
}

/// Run the traversal the tool's way: one generator for the whole window, drawn on per candidate.
fn traverse(
    normal: &[BamRecord],
    tumor: &[BamRecord],
    window: std::ops::RangeInclusive<i32>,
    error_probability: f64,
) -> Vec<NormalArtifactRecord> {
    let mut generator = JavaRandom::gatk();
    let mut records = Vec::new();
    for position in window {
        let outcome = get_normal_artifact_data::apply(
            &pileup_at(normal, position),
            &pileup_at(tumor, position),
            reference_base(position),
            error_probability,
            || generator.next_double(),
        );
        if let Outcome::Kept(record) = outcome {
            records.push(*record);
        }
    }
    records
}

#[test]
fn every_table_matches_the_golden() {
    let text = golden();
    let normal = normal_reads();
    let tumor = tumor_reads();

    for (label, window, error) in [
        ("default", 10..=19, DEFAULT_ERROR_PROBABILITY),
        ("high-error", 10..=19, 0.1),
        ("wide", 10..=49, DEFAULT_ERROR_PROBABILITY),
        ("one-locus", 12..=12, DEFAULT_ERROR_PROBABILITY),
        // Position 13's reference is `A`, which the alternate read carries, so nothing is a
        // candidate and the generator is never drawn on.
        ("no-alternate", 13..=13, DEFAULT_ERROR_PROBABILITY),
    ] {
        let records = traverse(&normal, &tumor, window, error);
        assert_eq!(
            get_normal_artifact_data::write(&records),
            file(&text, label),
            "{label}"
        );
    }
}

/// The floor case: thirty eligible loci, a keep probability of 0.05, and three survivors.
///
/// Three is what the seed gives: draws 4, 8 and 11 of `Random(47382911)` are at or below 0.05. The
/// dump resets the generator before each case for exactly this reason -- without the reset the
/// answer would depend on how many candidates the earlier cases had, which is a property of the
/// harness rather than of the tool.
#[test]
fn the_seeded_draw_keeps_three_loci_of_thirty() {
    let text = golden();
    let records = traverse(
        &normal_reads(),
        &matching_tumor(),
        10..=49,
        DEFAULT_ERROR_PROBABILITY,
    );
    assert_eq!(records.len(), 3, "the generator keeps three of the thirty");
    assert_eq!(records[0].downsampling, 0.05, "the floor");
    assert_eq!(records[0].tumor_alt_count, 0);
    assert_eq!(
        get_normal_artifact_data::write(&records),
        file(&text, "floor")
    );
}

/// A normal the BAM does not name leaves an empty normal pileup and therefore no rows.
#[test]
fn an_unknown_normal_sample_produces_a_header_and_nothing_else() {
    let text = golden();
    let records = traverse(&[], &tumor_reads(), 10..=19, DEFAULT_ERROR_PROBABILITY);
    assert!(records.is_empty());
    assert_eq!(
        get_normal_artifact_data::write(&records),
        file(&text, "unknown-normal")
    );
}

/// The first rule draws nothing and the second draws before deciding.
#[test]
fn the_draw_happens_between_the_two_rejection_rules() {
    let normal = normal_reads();
    let tumor = tumor_reads();
    // Position 13's reference is `A`: no normal alternate, so no draw.
    let mut drawn = false;
    let outcome = get_normal_artifact_data::apply(
        &pileup_at(&normal, 13),
        &pileup_at(&tumor, 13),
        reference_base(13),
        DEFAULT_ERROR_PROBABILITY,
        || {
            drawn = true;
            0.5
        },
    );
    assert_eq!(outcome, Outcome::SkippedBeforeDraw);
    assert!(
        !drawn,
        "the first rule rejects before the generator is touched"
    );

    // A locus that passes the first rule always draws, even when the second rejects it.
    let mut drawn = false;
    let outcome = get_normal_artifact_data::apply(
        &pileup_at(&normal, 10),
        &pileup_at(&tumor, 10),
        reference_base(10),
        DEFAULT_ERROR_PROBABILITY,
        || {
            drawn = true;
            // A draw above the keep probability, so this locus is dropped by the draw itself.
            1.0
        },
    );
    assert_eq!(outcome, Outcome::SkippedAfterDraw);
    assert!(drawn);
}

/// Mutect2's filter list, and the read length that keeps a short read out of `apply` entirely.
#[test]
fn the_filter_list_is_twelve_and_holds_a_read_length() {
    assert_eq!(STANDARD_MUTECT2_READ_FILTERS.len(), 12);
    assert!(STANDARD_MUTECT2_READ_FILTERS.contains(&"ReadLengthReadFilter"));
    assert_eq!(MIN_READ_LENGTH, 30);
}
