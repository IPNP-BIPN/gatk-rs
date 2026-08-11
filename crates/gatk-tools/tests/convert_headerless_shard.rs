//! Conformance for `ConvertHeaderlessHadoopBamShardToBam` against GATK 4.6.2.0, compared as
//! **bytes**.
//!
//! Golden from `tools/readfilter-conformance/ConvertHeaderlessShardDump.java`.
//!
//! # What this suite is for
//!
//! The first tool here that is not a `GATKTool` at all, and the only one whose central claim cannot
//! be checked by reading the output:
//!
//!  * **the shard is copied byte for byte**, so the assertion is a byte search and a three-part
//!    layout, not a read round trip. A port that decoded the records and wrote them back would pass
//!    every read-level check and fail this one;
//!  * **the header is encoded keeping the file's version**, which is the branch
//!    [`gatk_tools::print_reads_header`] does not take and the one htsjdk-rs#164 says the ordinary
//!    writer does not take either;
//!  * **the terminator appears exactly once**, because the header block is flushed rather than
//!    closed;
//!  * **a donor that is not a BAM is not refused**: the header comes out as `@HD VN:1.6` and the
//!    shard is appended to it anyway.

use gatk_corpus as corpus;
use gatk_tools::convert_headerless_shard::{self as tool, EMPTY_GZIP_BLOCK};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/convert_headerless_shard.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    text.lines()
        .find_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .unwrap_or_else(|| panic!("the golden lost its {kind} row"))
}

fn of_run<'a>(text: &'a str, kind: &str, label: &str) -> Vec<Vec<&'a str>> {
    rows(text, kind)
        .into_iter()
        .filter(|row| row[0] == label)
        .collect()
}

/// The header the reference read out of each run's `--bam-with-header`.
///
/// `plain` and `emptyshard` are given the donor; `badheader` is given the shard, which is not a BAM
/// and parses to an empty header.
fn header_for(label: &str, donor: &[u8]) -> SamHeader {
    match label {
        "plain" | "emptyshard" => {
            let decompressed =
                htsjdk_bgzf::read::decompress_all(donor).expect("the donor decompresses");
            BamReader::new(&decompressed)
                .expect("the donor is a BAM")
                .header
                .text
        }
        // SILENT validation on a headerless shard yields a header with nothing but its version.
        "badheader" => SamHeader::default(),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The shard each run was given.
fn shard_for(label: &str, shard: &[u8]) -> Vec<u8> {
    match label {
        "plain" | "badheader" => shard.to_vec(),
        "emptyshard" => Vec::new(),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();
    let donor = corpus::decode_base64(field(&text, "donor"));
    let shard = corpus::decode_base64(field(&text, "shard"));

    let outputs = rows(&text, "output");
    assert_eq!(outputs.len(), 3, "three runs, none refused");

    let mut compared = 0usize;
    for row in &outputs {
        let (label, expected_base64) = (row[0], row[1]);
        let ours =
            tool::convert_headerless_shard(&shard_for(label, &shard), &header_for(label, &donor))
                .expect("the header encodes");
        let expected = corpus::decode_base64(expected_base64);
        assert_eq!(ours.len(), expected.len(), "{label}: output length differs");
        if ours != expected {
            let at = ours
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            panic!("{label}: first byte difference at offset {at}");
        }
        compared += 1;
    }

    assert_eq!(compared, 3);
    println!("convert-headerless-shard: {compared} output files byte-identical");
}

/// The claim a read round trip cannot check: the shard's bytes are copied, not re-encoded.
#[test]
fn the_shard_is_copied_rather_than_re_encoded() {
    let text = golden();
    let shard = corpus::decode_base64(field(&text, "shard"));
    assert_eq!(
        shard.len().to_string(),
        field(&text, "shardlength"),
        "the golden's own length row"
    );

    // The reference reports the three parts of each output it measured.
    for row in rows(&text, "layout") {
        let (label, header_block, middle, terminator) = (
            row[0],
            row[1].parse::<usize>().unwrap(),
            row[2].parse::<usize>().unwrap(),
            row[3].parse::<usize>().unwrap(),
        );
        assert_eq!(terminator, EMPTY_GZIP_BLOCK.len(), "{label}");

        let expected = corpus::decode_base64(
            of_run(&text, "output", label)
                .first()
                .map(|row| row[1])
                .unwrap_or_else(|| panic!("no output for {label}")),
        );
        assert_eq!(
            expected.len(),
            header_block + middle + terminator,
            "{label}"
        );
        assert_eq!(
            &expected[header_block..header_block + middle],
            &shard_for(label, &shard)[..],
            "{label}: the middle is the shard, byte for byte"
        );
        assert_eq!(&expected[expected.len() - terminator..], &EMPTY_GZIP_BLOCK);
    }

    // And the reference said so itself.
    for row in rows(&text, "copiedverbatim") {
        assert_eq!(row[1], "true", "{}", row[0]);
    }
}

/// The terminator appears once, at the end, because the header block is flushed rather than closed.
#[test]
fn the_terminator_appears_exactly_once() {
    let text = golden();
    for row in rows(&text, "output") {
        let bytes = corpus::decode_base64(row[1]);
        let occurrences = bytes
            .windows(EMPTY_GZIP_BLOCK.len())
            .filter(|w| *w == EMPTY_GZIP_BLOCK)
            .count();
        assert_eq!(occurrences, 1, "{}", row[0]);
        assert_eq!(
            &bytes[bytes.len() - EMPTY_GZIP_BLOCK.len()..],
            &EMPTY_GZIP_BLOCK
        );
    }
}

/// A donor that is not a BAM produces a valid BAM with an empty header, rather than a refusal.
#[test]
fn a_donor_that_is_not_a_bam_is_not_refused() {
    let text = golden();
    assert!(
        rows(&text, "error").is_empty(),
        "the reference refuses nothing here"
    );
    let header = of_run(&text, "outputheader", "badheader")
        .first()
        .map(|row| row[1].to_string())
        .expect("the golden lost the badheader run");
    assert_eq!(header, "@HD\\tVN:1.6\\n");
    assert_eq!(SamHeader::default().encode(), "@HD\tVN:1.6\n");
}

/// Nothing is appended to a valid donor's header: there is no GATKTool to append a @PG.
#[test]
fn the_donors_header_survives_unchanged() {
    let text = golden();
    let header = of_run(&text, "outputheader", "plain")
        .first()
        .map(|row| row[1].to_string())
        .expect("the golden lost the plain run");
    assert!(header.contains("@PG\\tID:upstream"), "{header}");
    assert!(!header.contains("ConvertHeaderless"), "{header}");
    assert!(header.contains("@CO\\ta donor comment"), "{header}");
}
