//! The conformance corpus, as it travels in every golden of this programme.
//!
//! One parser, used by every crate's tests. Two parsers would be two corpora: a golden says which
//! records the reference judged, and a second reading of those rows is a second set of records
//! that can differ in exactly the ways the goldens exist to catch.

use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;

/// Read and decompress a gzipped golden.
pub fn read_golden(path: &std::path::Path) -> String {
    use std::io::Read;
    let file = std::fs::File::open(path).unwrap_or_else(|e| {
        panic!(
            "{}: {e}. Regenerate with tools/conformance/run_suite.py --suites readfilters",
            path.display()
        )
    });
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden is not valid gzip");
    text
}

/// Decode a base64 field of a golden, standard alphabet, `=` padded.
///
/// Some fixtures are binary: the `ReadsDataSource` golden carries the BAM it queried and that
/// BAM's `.bai`, so the port reads the bytes the reference read rather than a file rebuilt to
/// match a description of them. Written out rather than pulled in, because a dependency to undo
/// `java.util.Base64` in a test is a dependency to audit for the life of the programme.
pub fn decode_base64(text: &str) -> Vec<u8> {
    fn value(byte: u8) -> Option<u32> {
        match byte {
            b'A'..=b'Z' => Some((byte - b'A') as u32),
            b'a'..=b'z' => Some((byte - b'a') as u32 + 26),
            b'0'..=b'9' => Some((byte - b'0') as u32 + 52),
            b'+' => Some(62),
            b'/' => Some(63),
            _ => None,
        }
    }

    let mut out = Vec::with_capacity(text.len() / 4 * 3);
    let mut accumulator = 0u32;
    let mut bits = 0u32;
    for byte in text.bytes() {
        let Some(value) = value(byte) else {
            continue; // padding and whitespace
        };
        accumulator = (accumulator << 6) | value;
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((accumulator >> bits) as u8);
        }
    }
    out
}

/// The header the reference judged the corpus against.
///
/// It travels in the golden because the resolved filters read the library, sample, platform and
/// contig lengths out of it: a port given a different header would be answering a different
/// question and could agree by accident.
pub fn header(text: &str) -> SamHeader {
    let mut header = SamHeader::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        match parts[0] {
            "sq" => header
                .sequences
                .push(SequenceRecord::new(parts[2], parts[3].parse().unwrap())),
            "rg" => {
                let mut group = ReadGroup::new(parts[1]);
                for field in &parts[2..] {
                    let (key, value) = field.split_once('=').expect("an @RG field is KEY=value");
                    // "null" is how the dump prints an absent attribute; setting it would make the
                    // port match on a string the reference never had.
                    if value != "null" {
                        group.attributes.set(key, value);
                    }
                }
                header.read_groups.push(group);
            }
            _ => {}
        }
    }
    header
}

/// The corpus, in the order the reference judged it.
pub fn records(text: &str) -> Vec<BamRecord> {
    let mut records = Vec::new();
    for line in text.lines() {
        let mut parts = line.splitn(3, '\t');
        if parts.next() != Some("record") {
            continue;
        }
        let index: usize = parts
            .next()
            .expect("a record row has an index")
            .parse()
            .unwrap();
        let fields: Vec<&str> = parts
            .next()
            .expect("a record row has fields")
            .split('|')
            .collect();
        assert_eq!(
            fields.len(),
            11,
            "record {index} has {} fields",
            fields.len()
        );

        // Fields rather than a SAM line, because the corpus contains a record whose flags say
        // mapped while its reference is absent, which is one of the three criteria of
        // GATKRead.isUnmapped and exactly what htsjdk's reader refuses to parse
        // ("RNAME is not specified but flags indicate mapped"). Routing the corpus through SAM
        // text would drop the case the filter most needs.
        let mut record = BamRecord {
            read_name: fields[0].to_string(),
            flags: fields[1].parse().unwrap(),
            reference_index: fields[2].parse().unwrap(),
            alignment_start: fields[3].parse().unwrap(),
            mapping_quality: fields[4].parse().unwrap(),
            mate_reference_index: fields[6].parse().unwrap(),
            mate_alignment_start: fields[7].parse().unwrap(),
            inferred_insert_size: fields[8].parse().unwrap(),
            read_bases: fields[9].as_bytes().to_vec(),
            base_qualities: if fields[10].is_empty() {
                Vec::new()
            } else {
                fields[10].split(',').map(|q| q.parse().unwrap()).collect()
            },
            ..BamRecord::default()
        };
        if fields[5] != "*" {
            record.cigar = htsjdk_bam::text_parse::parse_cigar(fields[5])
                .unwrap_or_else(|e| panic!("record {index} cigar does not parse: {e:?}"));
        }

        assert_eq!(
            records.len(),
            index,
            "records are out of order in the golden"
        );
        records.push(record);
    }

    // Tags travel on their own rows: an OA value ends with a semicolon, so any in-line separator
    // would collide with the data it carries. The type travels with them because `tp` is a byte
    // array, and printing one as text would carry its identity hash into the corpus.
    for line in text.lines() {
        let mut parts = line.splitn(5, '\t');
        if parts.next() != Some("tag") {
            continue;
        }
        let index: usize = parts.next().unwrap().parse().unwrap();
        let name = parts.next().expect("a tag row has a name").as_bytes();
        let kind = parts.next().expect("a tag row has a type");
        let value = parts.next().expect("a tag row has a value");
        let parsed = match kind {
            "bytes" => htsjdk_bam::tag::TagValue::ByteArray {
                values: value
                    .split(',')
                    .filter(|v| !v.is_empty())
                    .map(|v| v.parse().expect("a signed byte"))
                    .collect(),
                unsigned: false,
            },
            _ => htsjdk_bam::tag::TagValue::Str(value.to_string()),
        };
        records[index]
            .tags
            .insert(htsjdk_bam::Tag::new(&[name[0], name[1]]), parsed);
    }
    assert!(!records.is_empty(), "the golden carries no records");
    records
}
