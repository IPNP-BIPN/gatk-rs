//! Conformance for ITF8 and LTF8, against `htsjdk.samtools.cram.io.ITF8` and `LTF8`.
//!
//! Goldens from `tools/cram-conformance/CramVarintDump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the ones where two different streams read the same and
//! where a broken stream reads a number:
//!
//! ```text
//! itf8read  five-byte-agreeing          f000000112  18
//! itf8read  five-byte-nibble-disagrees  f0000001f2  18
//! itf8read  truncated-two               80          -1
//! itf8read  empty                       -           E:htsjdk.samtools.util.RuntimeEOFException:null
//! ```
//!
//! One byte differs between the first two and the answer does not, because the fifth byte's high
//! nibble is written and then discarded. And a stream that stops early is not refused: only an
//! empty one is.

use std::io::Read;

use htsjdk_cram::varint::{
    read_unsigned_itf8, read_unsigned_ltf8, write_unsigned_itf8, write_unsigned_ltf8,
};

/// The `itf8` round-trip values, in the dump's order.
const INTS: &[i32] = &[
    0,
    1,
    63,
    127,
    128,
    129,
    255,
    256,
    16383,
    16384,
    16385,
    2097151,
    2097152,
    2097153,
    268435455,
    268435456,
    268435457,
    i32::MAX - 1,
    i32::MAX,
    -1,
    -2,
    -128,
    -129,
    i32::MIN,
    i32::MIN + 1,
];

/// The `ltf8` round-trip values, in the dump's order.
const LONGS: &[i64] = &[
    0,
    1,
    127,
    128,
    16383,
    16384,
    2097151,
    2097152,
    268435455,
    268435456,
    34359738367,
    34359738368,
    4398046511103,
    4398046511104,
    562949953421311,
    562949953421312,
    72057594037927935,
    72057594037927936,
    i64::MAX - 1,
    i64::MAX,
    -1,
    -2,
    i64::MIN,
];

/// The hand-made ITF8 streams, in the dump's order.
const ITF8_STREAMS: &[(&str, &[u8])] = &[
    ("five-byte-agreeing", &[0xF0, 0x00, 0x00, 0x01, 0x12]),
    ("five-byte-nibble-disagrees", &[0xF0, 0x00, 0x00, 0x01, 0xF2]),
    ("five-byte-high-nibble-only", &[0xF0, 0x00, 0x00, 0x00, 0xF0]),
    ("one-byte-zero", &[0x00]),
    ("two-byte-zero", &[0x80, 0x00]),
    ("three-byte-zero", &[0xC0, 0x00, 0x00]),
    ("four-byte-zero", &[0xE0, 0x00, 0x00, 0x00]),
    ("five-byte-zero", &[0xF0, 0x00, 0x00, 0x00, 0x00]),
    ("truncated-two", &[0x80]),
    ("truncated-five", &[0xF0, 0x00]),
    ("empty", &[]),
];

const LTF8_STREAMS: &[(&str, &[u8])] = &[
    ("ltf8-eight-byte-zero", &[0xFE, 0, 0, 0, 0, 0, 0, 0]),
    ("ltf8-nine-byte-zero", &[0xFF, 0, 0, 0, 0, 0, 0, 0, 0]),
    ("ltf8-truncated", &[0xFF, 0]),
    ("ltf8-empty", &[]),
];

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_varint.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn hex(bytes: &[u8]) -> String {
    if bytes.is_empty() {
        return "-".to_string();
    }
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

#[test]
fn every_varint_encodes_and_decodes_as_the_reference_does() {
    let corpus = corpus();
    let golden: Vec<&str> = corpus.lines().filter(|line| !line.starts_with('#')).collect();
    let mut produced: Vec<String> = Vec::new();

    for value in INTS {
        let (bytes, bits) = write_unsigned_itf8(*value);
        let read = read_unsigned_itf8(&bytes).expect("what the writer wrote is readable").0;
        produced.push(format!("itf8\t{value}\t{}\t{bits}\t{read}", hex(&bytes)));
    }
    for value in LONGS {
        let (bytes, bits) = write_unsigned_ltf8(*value);
        let read = read_unsigned_ltf8(&bytes).expect("what the writer wrote is readable").0;
        produced.push(format!("ltf8\t{value}\t{}\t{bits}\t{read}", hex(&bytes)));
    }
    for (label, bytes) in ITF8_STREAMS {
        let outcome = match read_unsigned_itf8(bytes) {
            Ok((value, _)) => value.to_string(),
            Err(error) => format!("E:{}:{}", error.class(), error.message()),
        };
        produced.push(format!("itf8read\t{label}\t{}\t{outcome}", hex(bytes)));
    }
    for (label, bytes) in LTF8_STREAMS {
        let outcome = match read_unsigned_ltf8(bytes) {
            Ok((value, _)) => value.to_string(),
            Err(error) => format!("E:{}:{}", error.class(), error.message()),
        };
        produced.push(format!("ltf8read\t{label}\t{}\t{outcome}", hex(bytes)));
    }

    assert_eq!(
        produced.len(),
        golden.len(),
        "the port produced {} rows and the golden has {}",
        produced.len(),
        golden.len()
    );
    for (index, (mine, theirs)) in produced.iter().zip(golden.iter()).enumerate() {
        assert_eq!(mine, theirs, "row {index}");
    }
}

/// Stated on its own because it is the property a port is most likely to "fix".
#[test]
fn the_fifth_byte_carries_four_bits_and_the_other_four_are_ignored() {
    for high in 0u8..16 {
        let byte = (high << 4) | 0x02;
        assert_eq!(
            read_unsigned_itf8(&[0xF0, 0x00, 0x00, 0x01, byte]).unwrap().0,
            18,
            "high nibble {high:x} changed the answer"
        );
    }
}
