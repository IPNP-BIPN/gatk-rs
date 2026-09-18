//! Conformance for Kryo 4's `Output`, against GATK 4.6.2.0's own Kryo.
//!
//! Golden from `tools/readfilter-conformance/KryoStreamDump.java`, which wrote each case through a
//! real `com.esotericsoftware.kryo.io.Output` in the pinned container and printed the bytes.
//!
//! # What this suite is for
//!
//!  * **a fixed `int` being four bytes big-endian, where `writeInt(value, true)` is a varint**;
//!  * **a negative varint costing all five bytes**, which is what "optimize positive" means;
//!  * **a string taking the ASCII shape only above one character**, with no length and the top bit
//!    set on its last byte, and the length shape otherwise;
//!  * **and a null string being one byte, which is how the length tells it from an empty one.**
//!
//! # What is NOT compared here yet
//!
//! The golden also carries the containers and the two whole objects: `LongHopscotchSet`,
//! `LargeLongHopscotchSet`, `PSTree`, `PSKmerSet` and `PSTaxonomyDatabase`. Those bytes are the
//! containers' own layouts -- a hopscotch table's iteration order under FNV-1a, and a `HashMap`'s
//! over the accessions -- and the port does not carry them yet. They are measured so that the port
//! has something to be right against; IPNP-BIPN/gatk-rs#1181 tracks the half that remains, and the
//! rows are asserted here to be PRESENT so that a dump which stops producing them fails this suite
//! rather than passing it quietly.

use gatk_corpus as corpus;
use gatk_engine::kryo::Output;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/kryo_stream.txt.gz"),
    )
}

/// The hex the golden carries for one label.
fn stream(text: &str, label: &str) -> String {
    let prefix = format!("stream\t{label}=");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
        .unwrap_or_else(|| panic!("{label}"))
}

/// What the port writes for one case, as the same hex.
fn written(write: impl FnOnce(&mut Output)) -> String {
    let mut output = Output::new();
    write(&mut output);
    output.hex()
}

#[test]
fn a_fixed_int_is_four_bytes_and_a_varint_is_not() {
    let text = golden();
    assert_eq!(written(|out| out.write_int(0)), stream(&text, "int-zero"));
    assert_eq!(written(|out| out.write_int(1)), stream(&text, "int-one"));
    assert_eq!(
        written(|out| out.write_int(-1)),
        stream(&text, "int-minus-one")
    );
    assert_eq!(
        written(|out| out.write_int(i32::MAX)),
        stream(&text, "int-max")
    );
    assert_eq!(
        written(|out| out.write_int(i32::MIN)),
        stream(&text, "int-min")
    );
    // The same numbers through the varint, which is the encoding a port reaches for by mistake.
    assert_eq!(
        written(|out| out.write_varint(300)),
        stream(&text, "varint-300")
    );
    assert_eq!(
        written(|out| out.write_varint(-300)),
        stream(&text, "varint-minus-300")
    );
    assert_eq!(
        stream(&text, "varint-minus-300").len() / 2,
        5,
        "a negative varint is the full five bytes"
    );
}

#[test]
fn a_fixed_long_is_eight_bytes() {
    let text = golden();
    assert_eq!(written(|out| out.write_long(0)), stream(&text, "long-zero"));
    assert_eq!(
        written(|out| out.write_long(-1)),
        stream(&text, "long-minus-one")
    );
    assert_eq!(
        written(|out| out.write_long(i64::MAX)),
        stream(&text, "long-max")
    );
    // The mask a `PSKmerSet` carries, which is the one long this port writes for real.
    assert_eq!(
        written(|out| out.write_long(1_048_575)),
        stream(&text, "long-kmer-mask")
    );
}

#[test]
fn a_string_is_ascii_without_a_length_or_utf8_with_one() {
    let text = golden();
    assert_eq!(
        written(|out| out.write_string(None)),
        stream(&text, "string-null")
    );
    assert_eq!(
        written(|out| out.write_string(Some(""))),
        stream(&text, "string-empty")
    );
    // One character is not enough for the ASCII shape, which is the rule that surprises.
    assert_eq!(
        written(|out| out.write_string(Some("A"))),
        stream(&text, "string-one-char")
    );
    assert_eq!(
        written(|out| out.write_string(Some("1"))),
        stream(&text, "string-ascii")
    );
    assert_eq!(
        written(|out| out.write_string(Some("131567"))),
        stream(&text, "string-number")
    );
    assert_eq!(
        written(|out| out.write_string(Some("Homo sapiens neanderthalensis and more xx"))),
        stream(&text, "string-long-ascii")
    );
    assert_eq!(
        written(|out| out.write_string(Some("Crème"))),
        stream(&text, "string-utf8")
    );
}

#[test]
fn an_object_is_preceded_by_a_reference_marker() {
    let text = golden();
    // Kryo 4 has references on by default, so every object in these streams starts with `01`.
    for label in [
        "hopscotch-eight",
        "large-hopscotch",
        "kmer-set",
        "tree-three-nodes",
    ] {
        let bytes = stream(&text, label);
        assert!(bytes.starts_with("01"), "{label}: {bytes}");
    }
    assert_eq!(
        written(|out| out.write_object(true, |inner| inner.write_int(0))),
        "0100000000"
    );
}

/// The container rows exist and say what the port still owes, rather than being quietly absent.
#[test]
fn the_containers_are_measured_even_where_the_port_does_not_write_them_yet() {
    let text = golden();
    // The capacity a hopscotch table writes is the legal size ABOVE the one asked for, and the two
    // cases here both land on 251, so their streams are the same bytes. That is the measurement,
    // and a port that reproduces the requested capacity would be wrong.
    assert_eq!(
        stream(&text, "hopscotch-eight"),
        stream(&text, "hopscotch-sixty-four")
    );
    assert!(stream(&text, "hopscotch-empty").starts_with("01000000fb"));
    // The same three accessions in the two map kinds do NOT produce the same bytes: the tool builds
    // a `HashMap`, so its iteration order reaches the file and a port owes that order too.
    assert_ne!(
        stream(&text, "taxonomy-linked-map"),
        stream(&text, "taxonomy-hash-map")
    );
    let order: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with("key\ttaxonomy-hash-map\t"))
        .map(|line| line.rsplit_once('=').expect("a key").1)
        .collect();
    assert_eq!(
        order,
        vec!["NZ_CP009273.1", "NC_000913.3", "NC_002695.2"],
        "the HashMap's order is the file's"
    );
}
