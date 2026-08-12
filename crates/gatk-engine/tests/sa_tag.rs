//! Conformance for `SATagBuilder` against GATK 4.6.2.0, compared as text.
//!
//! Golden from `tools/readfilter-conformance/SATagBuilderDump.java`. The fixtures are reads built
//! here rather than a BAM, because the class works on reads in memory and never touches a file.
//!
//! # What this suite is for
//!
//!  * **a unit is six fields split with a limit of -1**, so an empty trailing NM parses and five
//!    fields do not;
//!  * **the three validations are position, cigar and mapping quality**, each accepting `*`, and
//!    **NM is not validated at all**;
//!  * **the strand is normalised on the way out**, the one field that does not round trip;
//!  * **a non-supplementary read goes to the front of the list**, so the primary is the first unit
//!    of every tag;
//!  * **existing SA tags are preserved and stay first**;
//!  * **and a family of one writes nothing at all**.

use gatk_corpus as corpus;
use gatk_engine::read::flags;
use gatk_engine::sa_tag::{self, SaRead, SaTagBuilder};
use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const CONTIG_LENGTH: i32 = 5000;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sa_tag.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    for contig in ["chr1", "chr2", "chr9"] {
        header
            .sequences
            .push(SequenceRecord::new(contig, CONTIG_LENGTH));
    }
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "s1");
    header.read_groups.push(group);
    header
}

fn read(
    name: &str,
    contig: i32,
    start: i32,
    cigar: &str,
    flags: u16,
    mapq: u8,
    nm: Option<i64>,
) -> BamRecord {
    let mut record = BamRecord {
        read_name: name.to_string(),
        flags,
        reference_index: contig,
        alignment_start: start,
        mapping_quality: mapq,
        cigar: parse_cigar(cigar),
        read_bases: b"ACGTACGTAC".to_vec(),
        base_qualities: vec![35; 10],
        ..Default::default()
    };
    record
        .tags
        .insert(Tag::new(b"RG"), TagValue::Str("rg1".into()));
    if let Some(nm) = nm {
        record.tags.insert(Tag::new(b"NM"), TagValue::Int(nm));
    }
    record
}

fn parse_cigar(text: &str) -> Cigar {
    use htsjdk_bam::cigar::{CigarElement, Op};
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

/// The tag each labelled parse case starts from, mirroring the dump's table.
fn parse_input(label: &str) -> &'static str {
    match label {
        "one" => "chr1,100,+,10M,60,2;",
        "two" => "chr1,100,+,10M,60,2;chr2,200,-,5M5S,30,0;",
        "empty-nm" => "chr1,100,+,10M,60,;",
        "star-pos" => "chr1,*,+,10M,60,2;",
        "star-mapq" => "chr1,100,+,10M,*,2;",
        "star-cigar" => "chr1,100,+,*,60,2;",
        "odd-strand" => "chr1,100,x,10M,60,2;",
        "text-nm" => "chr1,100,+,10M,60,not-a-number;",
        "odd-cigar" => "chr1,100,+,1M1M1M,60,2;",
        "five-fields" => "chr1,100,+,10M,60;",
        "seven-fields" => "chr1,100,+,10M,60,2,extra;",
        "negative-pos" => "chr1,-1,+,10M,60,2;",
        "negative-mapq" => "chr1,100,+,10M,-1,2;",
        "bad-cigar" => "chr1,100,+,10Z,60,2;",
        "empty-cigar" => "chr1,100,+,,60,2;",
        other => panic!("no parse case {other}"),
    }
}

/// The read each labelled unit case is built from.
fn unit_input(label: &str) -> BamRecord {
    match label {
        "plain" => read("plain", 0, 100, "10M", 0, 60, None),
        "with-nm" => read("with-nm", 0, 100, "10M", 0, 60, Some(3)),
        "reverse" => read(
            "reverse",
            0,
            100,
            "10M",
            flags::READ_REVERSE_STRAND,
            60,
            Some(3),
        ),
        "zero-mapq" => read("zero-mapq", 0, 100, "10M", 0, 0, None),
        "unmapped" => {
            let mut record = read("unmapped", -1, 0, "*", flags::READ_UNMAPPED, 0, None);
            record.tags.remove(Tag::new(b"NM"));
            record
        }
        other => panic!("no unit case {other}"),
    }
}

/// The family each labelled group case is built from.
fn group_input(label: &str) -> Vec<BamRecord> {
    match label {
        "pair" => vec![
            read("piece", 0, 100, "5M5S", 0, 60, Some(1)),
            read("piece", 0, 200, "5S5M", 0, 60, Some(2)),
        ],
        "triple" => vec![
            read("piece", 0, 100, "3M7S", 0, 60, Some(1)),
            read("piece", 0, 200, "3S3M4S", 0, 60, Some(2)),
            read("piece", 0, 300, "6S4M", 0, 60, Some(3)),
        ],
        "existing-tag" => {
            let mut primary = read("piece", 0, 100, "5M5S", 0, 60, Some(1));
            primary.tags.insert(
                sa_tag::SA_TAG,
                TagValue::Str("chr9,999,-,4M,20,0;".to_string()),
            );
            vec![primary, read("piece", 0, 200, "5S5M", 0, 60, Some(2))]
        }
        "primary-already-supplementary" => vec![
            read(
                "piece",
                0,
                100,
                "5M5S",
                flags::SUPPLEMENTARY_ALIGNMENT,
                60,
                Some(1),
            ),
            read("piece", 0, 200, "5S5M", 0, 60, Some(2)),
        ],
        "single" => vec![read("alone", 0, 100, "10M", 0, 60, Some(1))],
        other => panic!("no group case {other}"),
    }
}

#[test]
fn every_parse_is_the_reference() {
    let text = golden();
    let header = header();

    let mut compared = 0;
    for row in rows(&text, "parse") {
        let (label, expected) = (row[0], row[1]);
        let mut record = read("carrier", 0, 1, "10M", 0, 60, None);
        record.tags.insert(
            sa_tag::SA_TAG,
            TagValue::Str(parse_input(label).to_string()),
        );
        let builder = SaTagBuilder::new(&record, &header).expect("the tag parses");
        builder.set_sa_tag(&mut record);
        let ours = match record.tags.get(sa_tag::SA_TAG) {
            Some(TagValue::Str(text)) => text.clone(),
            other => panic!("no SA tag: {other:?}"),
        };
        assert_eq!(ours, expected, "parse/{label}");
        compared += 1;
    }
    assert!(compared >= 9, "the golden has the round-tripping cases");
    println!("sa-tag: {compared} tags round tripped");
}

#[test]
fn every_refusal_is_the_reference() {
    let text = golden();
    let header = header();

    let mut compared = 0;
    for row in rows(&text, "error") {
        let (label, expected) = (row[0], row[1]);
        let mut record = read("carrier", 0, 1, "10M", 0, 60, None);
        record.tags.insert(
            sa_tag::SA_TAG,
            TagValue::Str(parse_input(label).to_string()),
        );
        let error = SaTagBuilder::new(&record, &header).expect_err("this tag is refused");
        // The class is the reference's, and so is the message.
        let ours = format!(
            "org.broadinstitute.hellbender.exceptions.GATKException:{}",
            error.message()
        );
        assert_eq!(ours, expected, "error/{label}");
        compared += 1;
    }
    assert_eq!(compared, 6, "six refusals");
}

#[test]
fn every_unit_is_the_reference() {
    let text = golden();
    let header = header();

    for row in rows(&text, "unit") {
        let (label, expected) = (row[0], row[1]);
        let unit = SaRead::from_record(&unit_input(label), &header);
        assert_eq!(unit.to_text(), expected, "unit/{label}");
    }
}

#[test]
fn every_family_is_the_reference() {
    let text = golden();
    let header = header();

    // The golden's rows are grouped by label and ordered primary first.
    let all = rows(&text, "group");
    let labels: Vec<&str> = {
        let mut seen = Vec::new();
        for row in &all {
            if !seen.contains(&row[0]) {
                seen.push(row[0]);
            }
        }
        seen
    };

    for label in labels {
        let mut family = group_input(label);
        let mut primary = family.remove(0);
        sa_tag::set_reads_as_supplemental(&mut primary, &mut family, &header)
            .expect("the family's own tags parse");

        let mut ours = Vec::new();
        for record in std::iter::once(&primary).chain(family.iter()) {
            let tag = match record.tags.get(sa_tag::SA_TAG) {
                Some(TagValue::Str(text)) => text.clone(),
                _ => "absent".to_string(),
            };
            ours.push(vec![
                format!("{}:{}", record.read_name, record.alignment_start),
                record.flags.to_string(),
                tag,
            ]);
        }

        let expected: Vec<Vec<String>> = all
            .iter()
            .filter(|row| row[0] == label)
            .map(|row| row[1..].iter().map(|s| s.to_string()).collect())
            .collect();
        assert_eq!(ours, expected, "group/{label}");
    }
}
