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
//!  * **and the k-mer file whole**: `LongHopscotchSet`'s table order under FNV-1a,
//!    `LargeLongHopscotchSet`'s partitions, and the `PSKmerSet` `PathSeqBuildKmers` writes, byte
//!    for byte against the reference's own.
//!
//! # What is NOT compared here yet
//!
//! The taxonomy half: `PSTree` and `PSTaxonomyDatabase`. Their bytes carry a `HashMap`'s iteration
//! order -- over integer node ids in the first and over accession strings in the second -- which is
//! the obligation `docs/an-unspecified-order-that-reaches-the-output.md` already describes and
//! which this port does not meet yet. The rows are measured, and asserted here to be PRESENT, so a
//! dump that stopped producing them would fail this suite rather than pass it quietly.
//! IPNP-BIPN/gatk-rs#1181 tracks that half.

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

/// The containers, which the port now writes: the same bytes the reference's own tables produced.
#[test]
fn the_hopscotch_tables_land_where_the_reference_put_them() {
    use gatk_engine::hopscotch::{LargeLongHopscotchSet, LongHopscotchSet};
    use gatk_tools::pathseq_kryo;

    const KMERS: [i64; 8] = [1, 2, 3, 17, 1024, 65535, 1_048_577, 123_456_789];
    let text = golden();

    let mut small = LongHopscotchSet::with_capacity(8);
    for value in KMERS {
        small.add(value);
    }
    assert_eq!(
        written(
            |out| out.write_object(true, |inner| pathseq_kryo::write_hopscotch_set(
                inner, &small
            ))
        ),
        stream(&text, "hopscotch-eight")
    );

    // Asked for sixty-four, the table is the same 251 buckets, so the bytes are the same too.
    let mut large = LongHopscotchSet::with_capacity(64);
    for value in KMERS {
        large.add(value);
    }
    assert_eq!(
        written(
            |out| out.write_object(true, |inner| pathseq_kryo::write_hopscotch_set(
                inner, &large
            ))
        ),
        stream(&text, "hopscotch-sixty-four")
    );

    let empty = LongHopscotchSet::with_capacity(8);
    assert_eq!(
        written(
            |out| out.write_object(true, |inner| pathseq_kryo::write_hopscotch_set(
                inner, &empty
            ))
        ),
        stream(&text, "hopscotch-empty")
    );

    let mut partitioned = LargeLongHopscotchSet::new(KMERS.len() as i64);
    for value in KMERS {
        partitioned.add(value);
    }
    assert_eq!(
        written(
            |out| out.write_object(true, |inner| pathseq_kryo::write_large_hopscotch_set(
                inner,
                &partitioned
            ))
        ),
        stream(&text, "large-hopscotch")
    );

    // And the file itself: k of 31 under the mask the corpus uses.
    let file = pathseq_kryo::kmer_set_file(31, 1_048_575, &partitioned);
    let hex: String = file.iter().map(|byte| format!("{byte:02x}")).collect();
    assert_eq!(hex, stream(&text, "kmer-set"));
}

/// The taxonomy rows say what the port still owes, rather than being quietly absent.
#[test]
fn the_taxonomy_is_measured_even_where_the_port_does_not_write_it_yet() {
    let text = golden();
    // The tree and the database are measured and not yet reproduced; a dump that stopped writing
    // them would fail here.
    assert!(!stream(&text, "tree-three-nodes").is_empty());
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
