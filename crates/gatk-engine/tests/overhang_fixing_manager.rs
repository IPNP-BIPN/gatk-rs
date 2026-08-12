//! Conformance for `OverhangFixingManager` against GATK 4.6.2.0, compared as text.
//!
//! Golden from `tools/readfilter-conformance/OverhangFixingManagerDump.java`, which drives the
//! reference's manager directly and reaches its two protected predicates and its package-private
//! key builder by reflection.
//!
//! # What this suite is for
//!
//!  * **the span gate is `> readLength / 2`**, integer division against a strict comparison;
//!  * **there are two ways to call an overhang mismatched**, the tolerance and the half rule;
//!  * **the two overhang predicates mix strict and non-strict comparisons**;
//!  * **a splice already held returns nothing**, and so does a manager with fixing turned off;
//!  * **the queue flushes halfway under pressure**;
//!  * **the family is repaired on the way out**: NM, MD and NH cleared, SA tags written;
//!  * **a secondary alignment is written but not clipped** unless asked for;
//!  * **and the mate key carries the other end of the pair, and the read's old start**.
//!
//! The written reads are compared field by field against the golden's SAM lines rather than as
//! bytes: this is an in-memory class with no file of its own, so the SAM text the reference prints
//! is the only rendering there is, and the fields are what the port owns.

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::overhang_fixing_manager::{self as ofm, OverhangArguments, OverhangFixingManager};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/overhang_fixing_manager.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn reference_bases(text: &str) -> String {
    rows(text, "reference")[0][0].to_string()
}

/// The second contig of the dump's reference, which only the splice lifecycle needs.
const CHR2: &str =
    "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG";

fn header(chr1_length: i32) -> SamHeader {
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", chr1_length));
    header
        .sequences
        .push(SequenceRecord::new("chr2", CHR2.len() as i32));
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "s1");
    header.read_groups.push(group);
    header
}

/// A reference callback over the dump's two contigs, 1-based and inclusive.
fn reference(chr1: &str) -> impl FnMut(&str, i32, i32) -> Result<Vec<u8>, String> + '_ {
    move |contig: &str, start: i32, end: i32| {
        let bases = match contig {
            "chr1" => chr1,
            "chr2" => CHR2,
            other => return Err(format!("unknown contig {other}")),
        };
        Ok(bases.as_bytes()[(start - 1) as usize..end as usize].to_vec())
    }
}

fn parse_cigar(text: &str) -> Cigar {
    if text == "*" {
        return Cigar::new(Vec::new());
    }
    let mut elements = Vec::new();
    let mut length = 0u32;
    for byte in text.bytes() {
        if byte.is_ascii_digit() {
            length = length * 10 + u32::from(byte - b'0');
            continue;
        }
        let op = match byte {
            b'M' => Op::M,
            b'I' => Op::I,
            b'D' => Op::D,
            b'N' => Op::N,
            b'S' => Op::S,
            b'H' => Op::H,
            b'P' => Op::P,
            b'=' => Op::Eq,
            b'X' => Op::X,
            other => panic!("no cigar operator {}", other as char),
        };
        elements.push(CigarElement { length, op });
        length = 0;
    }
    Cigar::new(elements)
}

fn read(name: &str, start: i32, cigar: &str, bases: &str, flags: u16) -> BamRecord {
    let mut record = BamRecord {
        read_name: name.to_string(),
        flags,
        reference_index: 0,
        alignment_start: start,
        mapping_quality: 60,
        cigar: parse_cigar(cigar),
        read_bases: bases.as_bytes().to_vec(),
        base_qualities: vec![35; bases.len()],
        ..Default::default()
    };
    record
        .tags
        .insert(Tag::new(b"RG"), TagValue::Str("rg1".into()));
    // The three tags the family repair clears.
    record.tags.insert(Tag::new(b"NM"), TagValue::Int(1));
    record
        .tags
        .insert(Tag::new(b"MD"), TagValue::Str("20".into()));
    record.tags.insert(Tag::new(b"NH"), TagValue::Int(2));
    record
}

/// The dump's `mutate`: the given zero-based offsets changed to a base that is not there.
fn mutate(bases: &str, offsets: &[usize]) -> String {
    let mut chars: Vec<char> = bases.chars().collect();
    for &offset in offsets {
        chars[offset] = if chars[offset] == 'A' { 'C' } else { 'A' };
    }
    chars.into_iter().collect()
}

/// The eleven mandatory SAM columns plus the tags, rendered the way htsjdk prints them.
///
/// The tags come out ordered by their **binary** tag value, `first + (second << 8)`, which is the
/// order `SAMRecord` keeps them in and not alphabetical: `NH` sorts after `RG`.
fn sam_fields(record: &BamRecord, header: &SamHeader) -> Vec<String> {
    let contig = |index: i32| -> String {
        usize::try_from(index)
            .ok()
            .and_then(|i| header.sequences.get(i))
            .map(|s| s.name.clone())
            .unwrap_or_else(|| "*".to_string())
    };
    let rnext = if record.mate_reference_index < 0 {
        "*".to_string()
    } else if record.mate_reference_index == record.reference_index {
        "=".to_string()
    } else {
        contig(record.mate_reference_index)
    };

    let mut fields = vec![
        record.read_name.clone(),
        record.flags.to_string(),
        contig(record.reference_index),
        record.alignment_start.to_string(),
        record.mapping_quality.to_string(),
        record.cigar.to_text(),
        rnext,
        record.mate_alignment_start.to_string(),
        record.inferred_insert_size.to_string(),
        String::from_utf8(record.read_bases.clone()).expect("the bases are ASCII"),
        record
            .base_qualities
            .iter()
            .map(|&q| (q + 33) as char)
            .collect(),
    ];

    let mut tags: Vec<(Tag, &TagValue)> = record
        .tags
        .iter()
        .map(|(tag, value)| (*tag, value))
        .collect();
    tags.sort_by_key(|(tag, _)| {
        let name = tag.name();
        i32::from(name[0]) + (i32::from(name[1]) << 8)
    });
    for (tag, value) in tags {
        let name = tag.name();
        let name = format!("{}{}", name[0] as char, name[1] as char);
        fields.push(match value {
            TagValue::Str(text) => format!("{name}:Z:{text}"),
            TagValue::Int(number) => format!("{name}:i:{number}"),
            other => panic!("no rendering for {other:?}"),
        });
    }
    fields
}

fn expected_fields(line: &str) -> Vec<String> {
    line.split("\\t").map(|s| s.to_string()).collect()
}

#[test]
fn every_span_decision_is_the_reference() {
    let text = golden();
    let read = b"ACGTACGTAC";
    let same = b"ACGTACGTAC";
    let one_off = b"TCGTACGTAC";
    let two_off = b"TTGTACGTAC";
    let default = OverhangArguments::default();
    let narrow = OverhangArguments {
        max_bases_in_overhang: 3,
        ..OverhangArguments::default()
    };
    let strict = OverhangArguments {
        max_mismatches_in_overhang: 0,
        ..OverhangArguments::default()
    };

    let ours = |label: &str| -> bool {
        match label {
            "span-0" => ofm::overhanging_bases_mismatch(read, 0, 10, two_off, 0, 0, &default),
            "span-negative" => {
                ofm::overhanging_bases_mismatch(read, 0, 10, two_off, 0, -1, &default)
            }
            "span-6-of-10" => ofm::overhanging_bases_mismatch(read, 0, 10, two_off, 0, 6, &default),
            "span-5-of-10" => ofm::overhanging_bases_mismatch(read, 0, 10, two_off, 0, 5, &default),
            "span-4-of-10" => ofm::overhanging_bases_mismatch(read, 0, 10, two_off, 0, 4, &default),
            "span-5-of-9" => ofm::overhanging_bases_mismatch(read, 0, 9, two_off, 0, 5, &default),
            "identical" => ofm::overhanging_bases_mismatch(read, 0, 10, same, 0, 4, &default),
            "one-of-four" => ofm::overhanging_bases_mismatch(read, 0, 10, one_off, 0, 4, &default),
            "one-of-two" => ofm::overhanging_bases_mismatch(read, 0, 10, one_off, 0, 2, &default),
            "two-of-four" => ofm::overhanging_bases_mismatch(read, 0, 10, two_off, 0, 4, &default),
            "span-above-max-bases" => {
                ofm::overhanging_bases_mismatch(read, 0, 10, two_off, 0, 4, &narrow)
            }
            "zero-tolerance" => {
                ofm::overhanging_bases_mismatch(read, 0, 10, one_off, 0, 4, &strict)
            }
            other => panic!("no span case {other}"),
        }
    };

    let cases = rows(&text, "mismatch");
    assert_eq!(cases.len(), 12, "twelve span decisions");
    for row in cases {
        assert_eq!(ours(row[0]).to_string(), row[1], "mismatch/{}", row[0]);
    }
}

#[test]
fn every_overhang_verdict_is_the_reference() {
    let text = golden();
    let splice = SimpleInterval::new("chr1", 50, 60).expect("a valid splice");

    for row in rows(&text, "overhang") {
        let (start, end) = row[0].split_once('-').expect("a start and an end");
        let loc = SimpleInterval::new("chr1", start.parse().unwrap(), end.parse().unwrap())
            .expect("a valid read location");
        let ours = if ofm::is_left_overhang(&loc, &splice) {
            "left"
        } else if ofm::is_right_overhang(&loc, &splice) {
            "right"
        } else {
            "neither"
        };
        assert_eq!(ours, row[1], "overhang/{}", row[0]);
    }
}

#[test]
fn every_mate_key_is_the_reference() {
    let text = golden();
    for row in rows(&text, "key") {
        let ours = match row[0] {
            "first" => ofm::make_key("read1", true, 100),
            "second" => ofm::make_key("read1", false, 100),
            "zero-start" => ofm::make_key("read1", true, 0),
            other => panic!("no key case {other}"),
        };
        assert_eq!(ours, row[1], "key/{}", row[0]);
    }
}

#[test]
fn the_splice_lifecycle_is_the_reference() {
    let text = golden();
    let chr1 = reference_bases(&text);
    let header = header(chr1.len() as i32);
    let mut source = reference(&chr1);

    let mut manager = OverhangFixingManager::new(&header, OverhangArguments::default());
    let mut observed: Vec<(String, String, String)> = Vec::new();
    let record = |manager: &OverhangFixingManager,
                  observed: &mut Vec<(String, String, String)>,
                  label: &str,
                  added: Option<ofm::Splice>| {
        let held: Vec<String> = manager
            .splices()
            .iter()
            .map(|splice| {
                format!(
                    "{}:{}-{}",
                    splice.loc.contig, splice.loc.start, splice.loc.end
                )
            })
            .collect();
        observed.push((
            label.to_string(),
            if added.is_none() { "null" } else { "new" }.to_string(),
            held.join(","),
        ));
    };

    for (label, contig, start, end) in [
        ("first", "chr1", 50, 60),
        ("again", "chr1", 50, 60),
        ("second", "chr1", 20, 30),
        ("new-contig", "chr2", 10, 20),
        ("back-to-chr1", "chr1", 50, 60),
    ] {
        let added = manager
            .add_splice_position(contig, start, end, &mut source)
            .expect("the reference is readable");
        record(&manager, &mut observed, label, added);
    }

    // With fixing off nothing is remembered and every call returns nothing.
    let mut off = OverhangFixingManager::new(
        &header,
        OverhangArguments {
            do_not_fix_overhangs: true,
            ..OverhangArguments::default()
        },
    );
    let added = off
        .add_splice_position("chr1", 50, 60, &mut source)
        .expect("nothing is read at all");
    record(&off, &mut observed, "fixing-off", added);

    let expected: Vec<(String, String, String)> = rows(&text, "splice")
        .into_iter()
        .map(|row| {
            (
                row[0].to_string(),
                row[1].to_string(),
                row.get(2).copied().unwrap_or("").to_string(),
            )
        })
        .collect();
    assert_eq!(observed, expected);

    // And writing can only be activated once.
    manager.activate_writing().expect("the first call");
    let error = manager
        .activate_writing()
        .expect_err("the second is refused");
    let expected = rows(&text, "error")[0][1];
    assert_eq!(
        format!(
            "org.broadinstitute.hellbender.exceptions.GATKException:{}",
            error.message()
        ),
        expected
    );
}

/// One labelled run of the dump's `clipping` section.
fn run(label: &str, chr1: &str, header: &SamHeader) -> (usize, Vec<Vec<String>>) {
    let arguments = OverhangArguments {
        max_records_in_memory: if label == "pressure" { 2 } else { 100 },
        process_secondary_reads: label == "secondary-processed",
        ..OverhangArguments::default()
    };
    let mut manager = OverhangFixingManager::new(header, arguments);
    let mut source = reference(chr1);
    manager
        .activate_writing()
        .expect("writing is on from the start");

    let family: Vec<BamRecord> = match label {
        "left-matching" => vec![read("left", 55, "20M", &chr1[54..74], 0)],
        "left-mismatching" | "splice-after-read" => {
            vec![read(
                "left",
                55,
                "20M",
                &mutate(&chr1[54..74], &[0, 1, 2]),
                0,
            )]
        }
        "right-mismatching" => vec![read(
            "right",
            40,
            "15M",
            &mutate(&chr1[39..54], &[10, 11, 12]),
            0,
        )],
        "family" => vec![
            read("piece", 45, "10M10S", &chr1[44..64], 0),
            read("piece", 70, "10S10M", &chr1[59..79], 0),
        ],
        "pressure" => vec![
            read("a", 10, "10M", &chr1[9..19], 0),
            read("b", 20, "10M", &chr1[19..29], 0),
            read("c", 30, "10M", &chr1[29..39], 0),
            read("d", 40, "10M", &chr1[39..49], 0),
        ],
        "secondary" | "secondary-processed" => vec![read(
            "secondary",
            55,
            "20M",
            &mutate(&chr1[54..74], &[0, 1, 2]),
            256,
        )],
        other => panic!("no clipping case {other}"),
    };
    let splices: Vec<(i32, i32)> = if label == "pressure" {
        vec![]
    } else {
        vec![(50, 60)]
    };

    let splice_first = label != "splice-after-read";
    if splice_first {
        for (start, end) in &splices {
            manager
                .add_splice_position("chr1", *start, *end, &mut source)
                .expect("the reference is readable");
        }
    }

    if family.len() > 1 && family[0].read_name == family[1].read_name {
        manager.add_read_group(&family).expect("a non-empty family");
    } else {
        for record in &family {
            manager
                .add_read_group(std::slice::from_ref(record))
                .expect("a non-empty family");
        }
    }

    if !splice_first {
        for (start, end) in &splices {
            manager
                .add_splice_position("chr1", *start, *end, &mut source)
                .expect("the reference is readable");
        }
    }

    let written_before_flush = manager.written.len();
    manager.flush().expect("the flush");
    let written = manager
        .written
        .iter()
        .map(|record| sam_fields(record, header))
        .collect();
    (written_before_flush, written)
}

#[test]
fn every_written_read_is_the_reference() {
    let text = golden();
    let chr1 = reference_bases(&text);
    let header = header(chr1.len() as i32);

    let queues = rows(&text, "queue");
    let written = rows(&text, "written");
    assert!(!queues.is_empty(), "the golden has the clipping runs");

    for queue in &queues {
        let label = queue[0];
        let (before_flush, ours) = run(label, &chr1, &header);
        assert_eq!(before_flush.to_string(), queue[1], "queue/{label}");

        let expected: Vec<Vec<String>> = written
            .iter()
            .filter(|row| row[0] == label)
            .map(|row| expected_fields(row[1]))
            .collect();
        assert_eq!(ours, expected, "written/{label}");
    }
}

#[test]
fn the_mate_repair_is_the_reference() {
    let text = golden();
    let chr1 = reference_bases(&text);
    let header = header(chr1.len() as i32);
    let mut source = reference(&chr1);

    let mut manager = OverhangFixingManager::new(&header, OverhangArguments::default());

    // First pass: the clipped read records the key its mate will look up.
    let clipped = read(
        "pair",
        55,
        "20M",
        &mutate(&chr1[54..74], &[0, 1, 2]),
        0x1 | 0x40,
    );
    manager
        .add_splice_position("chr1", 50, 60, &mut source)
        .expect("the reference is readable");
    manager
        .add_read_group(std::slice::from_ref(&clipped))
        .expect("a non-empty family");
    manager.flush().expect("the flush");

    let mate = |flags: u16, mate_start: i32, with_mc: bool| -> BamRecord {
        let mut record = BamRecord {
            read_name: "pair".to_string(),
            flags,
            reference_index: 0,
            alignment_start: 100,
            mapping_quality: 60,
            cigar: parse_cigar("20M"),
            mate_reference_index: 0,
            mate_alignment_start: mate_start,
            read_bases: chr1.as_bytes()[99..119].to_vec(),
            base_qualities: vec![35; 20],
            ..Default::default()
        };
        record
            .tags
            .insert(Tag::new(b"RG"), TagValue::Str("rg1".into()));
        if with_mc {
            record
                .tags
                .insert(Tag::new(b"MC"), TagValue::Str("20M".into()));
        }
        record
    };

    let mut observed: Vec<(String, String, Vec<String>)> = Vec::new();
    let apply = |manager: &OverhangFixingManager,
                 observed: &mut Vec<(String, String, Vec<String>)>,
                 label: &str,
                 mut record: BamRecord| {
        let edited = manager.set_predicted_mate_information(&mut record);
        observed.push((
            label.to_string(),
            if edited { "edited" } else { "untouched" }.to_string(),
            sam_fields(&record, &header),
        ));
    };

    apply(
        &manager,
        &mut observed,
        "before-activation",
        mate(0x1 | 0x80, 55, true),
    );
    manager
        .activate_writing()
        .expect("writing is activated once");
    apply(
        &manager,
        &mut observed,
        "after-activation",
        mate(0x1 | 0x80, 55, true),
    );
    apply(
        &manager,
        &mut observed,
        "no-mc-tag",
        mate(0x1 | 0x80, 55, false),
    );
    apply(&manager, &mut observed, "unpaired", mate(0, 55, true));
    apply(
        &manager,
        &mut observed,
        "wrong-mate-start",
        mate(0x1 | 0x80, 56, true),
    );
    apply(
        &manager,
        &mut observed,
        "first-of-pair",
        mate(0x1 | 0x40, 55, true),
    );

    let expected: Vec<(String, String, Vec<String>)> = rows(&text, "mate")
        .into_iter()
        .map(|row| {
            (
                row[0].to_string(),
                row[1].to_string(),
                expected_fields(row[2]),
            )
        })
        .collect();
    assert_eq!(observed, expected);
}
