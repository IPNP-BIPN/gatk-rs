//! Conformance for `PileupElement` against GATK 4.6.2.0.
//!
//! One row per element over 25 cigars, and one row per `createPileupForReadAndOffset` call at
//! every read offset, including the offsets inside a soft clip or an insertion where the reference
//! throws.
//!
//! Five of the cigars exist only to make the two cigar navigations disagree: `isBeforeInsertion`
//! looks at the adjacent element whatever it is, `isBeforeDeletionStart` skips to the next
//! on-genome one. On `3M2S3I3M` the golden shows the last base of the leading `3M` reporting
//! `isBeforeSoftClip` and *not* `isBeforeInsertion`, with `2S3I` between it and the next position.

use gatk_corpus as corpus;
use gatk_engine::alignment_state::AlignmentStateMachine;
use gatk_engine::pileup::PileupElement;
use gatk_engine::read_utils;
use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// The alignment start `PileupElementDump.makeRead` uses.
const ALIGNMENT_START: i32 = 101;

/// The cigars whose reads the harness gave `BI` and `BD` tags.
const TAGGED: [&str; 3] = ["10M", "5M3D5M", "3S7M"];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/pileup_element.txt.gz"),
    )
}

/// `SAMUtils.phredToFastq`: the tag the reference writes is FASTQ text, so the port has to write
/// the same text rather than the raw bytes.
fn phred_to_fastq(quals: &[u8]) -> String {
    quals.iter().map(|q| (q + 33) as char).collect()
}

fn read_for(cigar_text: &str) -> BamRecord {
    let cigar = htsjdk_bam::text_parse::parse_cigar(cigar_text).expect("a parsable cigar");
    let length = cigar.read_length() as usize;
    let mut tags = htsjdk_bam::tag::Tags::new();
    if TAGGED.contains(&cigar_text) {
        let insertion: Vec<u8> = (0..length).map(|i| 10 + i as u8).collect();
        let deletion: Vec<u8> = (0..length).map(|i| 30 + i as u8).collect();
        tags.insert(
            Tag::new(&read_utils::BQSR_BASE_INSERTION_QUALITIES),
            TagValue::Str(phred_to_fastq(&insertion)),
        );
        tags.insert(
            Tag::new(&read_utils::BQSR_BASE_DELETION_QUALITIES),
            TagValue::Str(phred_to_fastq(&deletion)),
        );
    }
    BamRecord {
        read_name: format!("read-{cigar_text}"),
        reference_index: 0,
        alignment_start: ALIGNMENT_START,
        mapping_quality: 60,
        read_bases: (0..length).map(|i| b"ACGT"[i % 4]).collect(),
        // Varied, as the harness varies them: a flat quality would hide an off-by-one.
        base_qualities: (0..length).map(|i| 20 + (i % 11) as u8).collect(),
        cigar,
        tags,
        ..Default::default()
    }
}

/// `CigarElement.toString`: length then the SAM character.
fn describe(element: Option<CigarElement>) -> String {
    match element {
        None => "null".to_string(),
        Some(element) => format!("{}{}", element.length, operator(element.op)),
    }
}

fn describe_all(elements: &[CigarElement]) -> String {
    if elements.is_empty() {
        return "-".to_string();
    }
    elements
        .iter()
        .map(|e| format!("{}{}", e.length, operator(e.op)))
        .collect()
}

fn operator(op: Op) -> String {
    match op {
        Op::Eq => "=".to_string(),
        other => format!("{other:?}").to_uppercase(),
    }
}

#[test]
fn every_element_matches_the_reference() {
    let text = golden();

    let mut cigars: Vec<String> = Vec::new();
    let mut rows: std::collections::HashMap<String, Vec<String>> = Default::default();
    let mut offsets: std::collections::HashMap<String, Vec<(i32, String)>> = Default::default();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("el\t") {
            let mut parts = rest.splitn(3, '\t');
            let cigar = parts.next().expect("a cigar").to_string();
            let _index = parts.next();
            rows.entry(cigar)
                .or_default()
                .push(parts.next().unwrap_or("").to_string());
        } else if let Some(rest) = line.strip_prefix("count\t") {
            let cigar = rest.split('\t').next().expect("a cigar").to_string();
            cigars.push(cigar);
        } else if let Some(rest) = line.strip_prefix("forOffset\t") {
            let mut parts = rest.split('\t');
            let cigar = parts.next().expect("a cigar").to_string();
            let offset: i32 = parts.next().expect("an offset").parse().expect("a number");
            let outcome = parts.next().expect("an outcome").to_string();
            offsets.entry(cigar).or_default().push((offset, outcome));
        }
    }
    assert!(!cigars.is_empty(), "the golden carries no count rows");

    let mut compared = 0;
    let mut probed = 0;
    for cigar in &cigars {
        let read = read_for(cigar);
        let expected = rows.get(cigar).cloned().unwrap_or_default();

        let mut ours: Vec<String> = Vec::new();
        let mut machine = AlignmentStateMachine::new(&read);
        while let Ok(Some(_)) = machine.step_forward_on_genome() {
            let element = PileupElement::from_state(&read, &machine).expect("not an edge");
            ours.push(format!(
                "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                element.offset,
                element.base() as char,
                element.qual(),
                element.base_insertion_qual(),
                element.base_deletion_qual(),
                element.is_deletion(),
                element.is_before_deletion_start(),
                element.is_after_deletion_end(),
                element.is_before_insertion(),
                element.is_after_insertion(),
                element.is_before_soft_clip(),
                element.is_after_soft_clip(),
                element.is_next_to_soft_clip(),
                element.at_start_of_current_cigar(),
                element.at_end_of_current_cigar(),
                element.length_of_immediately_following_indel(),
                match element.bases_of_immediately_following_insertion() {
                    None => "null".to_string(),
                    Some(bases) => String::from_utf8(bases).expect("ASCII bases"),
                },
                describe(element.previous_on_genome_cigar_element()),
                describe(element.next_on_genome_cigar_element()),
                describe_all(&element.between_prev_position()),
                describe_all(&element.between_next_position()),
                element.is_usable_base_for_annotation(),
            ));
        }
        assert_eq!(&ours, &expected, "{cigar}");
        compared += ours.len();

        // createPileupForReadAndOffset, including the offsets it refuses.
        for (offset, outcome) in offsets.get(cigar).into_iter().flatten() {
            let ours = match PileupElement::for_read_and_offset(&read, *offset) {
                None => "E:java.lang.IllegalStateException".to_string(),
                Some(element) => {
                    format!("ok:{}:{}", element.offset, element.current_cigar_offset)
                }
            };
            assert_eq!(&ours, outcome, "{cigar} at read offset {offset}");
            probed += 1;
        }
    }

    println!(
        "{compared} pileup elements over {} cigars, plus {probed} createPileupForReadAndOffset calls",
        cigars.len()
    );
}
