//! Conformance for `SAMPileupCodec` against GATK 4.6.2.0, compared as text.
//!
//! Golden from `tools/readfilter-conformance/SAMPileupCodecDump.java`. Every case is a line of the
//! samtools mpileup format, decoded or refused.
//!
//! # What this suite is for
//!
//!  * **the bases column is a little language**, whose markers consume qualities at their own rate;
//!  * **a deletion becomes the letter `D`**, which the reads column itself would not accept;
//!  * **a coverage of zero returns before the columns after it are read**;
//!  * **the field count is 4 to 6, not the format's six**, so seven columns are refused and five
//!    fall out of the array with Java's own exception;
//!  * **two messages carry typos** that a port has to carry too;
//!  * **and `canDecode` is the extension and nothing else**.

use gatk_corpus as corpus;
use gatk_engine::sam_pileup::{can_decode, decode};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sam_pileup.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The line each labelled case decodes, taken from the golden itself.
fn line_of(text: &str, label: &str) -> String {
    for row in rows(text, "decode") {
        if row[0] == label {
            return unescape(row[1]);
        }
    }
    // A refused line is not in a `decode` row, so the cases are mirrored here.
    match label {
        "seven-columns" => "chr1\t10\tA\t2\t.,\tII\t~~".to_string(),
        "five-columns" => "chr1\t10\tA\t2\t.,".to_string(),
        "eight-columns" => "chr1\t10\tA\t2\t.,\tII\t~~\textra".to_string(),
        "coverage-mismatch" => "chr1\t10\tA\t5\t.,\tII".to_string(),
        "bad-position" => "chr1\tten\tA\t2\t.,\tII".to_string(),
        "bad-reference" => "chr1\t10\tZ\t2\t.,\tII".to_string(),
        "bad-base" => "chr1\t10\tA\t2\t.Z\tII".to_string(),
        "indel-without-length" => "chr1\t10\tA\t2\t.+AC,\tII".to_string(),
        "too-few-qualities" => "chr1\t10\tA\t2\t.,\tI".to_string(),
        "too-many-qualities" => "chr1\t10\tA\t2\t.,\tIII".to_string(),
        other => panic!("no case {other}"),
    }
}

#[test]
fn every_decoded_line_is_the_reference() {
    let text = golden();
    let cases = rows(&text, "decode");
    assert_eq!(cases.len(), 10, "ten lines decode");

    for row in cases {
        let label = row[0];
        let feature = decode(&unescape(row[1])).unwrap_or_else(|error| {
            panic!("{label} was refused: {}", error.message());
        });
        // `chr:pos ref cov`, then the bases, then the qualities as Java prints an array.
        let summary = format!(
            "{}:{} {} {}",
            feature.contig,
            feature.position,
            feature.reference_base as char,
            feature.size()
        );
        assert_eq!(summary, row[2], "summary/{label}");
        assert_eq!(feature.bases_string(), row[3], "bases/{label}");

        let quals = format!(
            "[{}]",
            feature
                .base_quals()
                .iter()
                .map(|q| q.to_string())
                .collect::<Vec<String>>()
                .join(",")
        );
        assert_eq!(quals, row[4], "quals/{label}");
    }
}

#[test]
fn every_refusal_is_the_reference_including_its_class() {
    let text = golden();
    let cases = rows(&text, "error");
    assert_eq!(cases.len(), 10, "ten lines are refused");

    for row in cases {
        let label = row[0];
        let error = decode(&line_of(&text, label)).expect_err("this line is refused");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            unescape(row[1]),
            "error/{label}"
        );
    }
}

#[test]
fn can_decode_is_the_reference() {
    let text = golden();
    for row in rows(&text, "candecode") {
        assert_eq!(
            can_decode(row[0]).to_string(),
            row[1],
            "candecode/{}",
            row[0]
        );
    }
}
