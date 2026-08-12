//! Conformance for `BAQ` against GATK 4.6.2.0, compared **bit for bit** where it is a double and
//! byte for byte where it is a quality.
//!
//! Golden from `tools/readfilter-conformance/BaqDump.java`. This is the hidden Markov model
//! `BaseRecalibrator` caps qualities with, so every number in a recalibration table rests on it.
//!
//! # What this suite is for
//!
//!  * **confidence is about placement, not about matching**: an exact match against a unique
//!    reference scores up to 91, the same match against `ACACAC...` scores 4 everywhere;
//!  * **most of the emission table is one**, so an `N` costs almost nothing where a real mismatch
//!    costs a lot;
//!  * **the emission's quality is floored, not the read's**;
//!  * **the reference window is not the read's span**, and an insertion moves it;
//!  * **an insertion or a soft clip keeps its raw quality**, through a `case S:` that falls through;
//!  * **the BQ tag is the difference, not the value**.

use gatk_corpus as corpus;
use gatk_engine::baq::{
    calc_baq_from_tag, encode_bq_tag, qual_to_prob_table, reference_window_for_read,
    state_aligned_position, state_is_indel, Baq, DEFAULT_BANDWIDTH, DEFAULT_GOP,
};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/baq.txt.gz"),
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

fn constant(text: &str, name: &str) -> String {
    rows(text, "const")
        .into_iter()
        .find(|row| row[0] == name)
        .unwrap_or_else(|| panic!("the golden has no constant {name}"))[1]
        .to_string()
}

fn bits(value: f64) -> String {
    format!("{:x}", value.to_bits())
}

fn join(values: &[u8]) -> String {
    values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

/// A cigar written the way the dump writes it, `4M2D6M`.
fn parse_cigar(text: &str) -> Vec<CigarElement> {
    let mut elements = Vec::new();
    let mut length = 0u32;
    for c in text.chars() {
        if c.is_ascii_digit() {
            length = length * 10 + c.to_digit(10).unwrap();
            continue;
        }
        let op = match c {
            'M' => Op::M,
            'I' => Op::I,
            'D' => Op::D,
            'N' => Op::N,
            'S' => Op::S,
            'H' => Op::H,
            'P' => Op::P,
            '=' => Op::Eq,
            'X' => Op::X,
            other => panic!("no cigar operator {other}"),
        };
        elements.push(CigarElement { length, op });
        length = 0;
    }
    elements
}

/// A read of ten bases with one cigar, which is every fixture the dump built.
fn read(cigar: &str, start: i32) -> BamRecord {
    BamRecord {
        read_name: "read".to_string(),
        reference_index: 0,
        alignment_start: start,
        cigar: Cigar::new(parse_cigar(cigar)),
        read_bases: b"ACGTACGTAC".to_vec(),
        base_qualities: vec![40; 10],
        mapping_quality: 60,
        ..BamRecord::default()
    }
}

#[test]
fn the_constants_are_the_references() {
    let text = golden();
    let baq = Baq::default();
    assert_eq!(constant(&text, "DEFAULT_GOP"), format!("{DEFAULT_GOP:.1}"));
    assert_eq!(
        constant(&text, "DEFAULT_BANDWIDTH"),
        DEFAULT_BANDWIDTH.to_string()
    );
    assert_eq!(constant(&text, "BAQ_TAG"), "BQ");
    assert_eq!(
        constant(&text, "minBaseQual"),
        baq.min_base_qual().to_string()
    );
    assert_eq!(constant(&text, "bandWidth"), baq.band_width().to_string());
    // The gap open probability is the Phred default converted, and the reference prints it the way
    // Java writes a double.
    assert_eq!(constant(&text, "gapOpenProb"), "1.0E-4");
    assert_eq!(baq.gap_open_prob(), 1.0e-4);
    assert_eq!(constant(&text, "gapExtensionProb"), "0.1");
    assert_eq!(baq.gap_extension_prob(), 0.1);
}

/// The quality cache and the emission table, as raw bits.
#[test]
fn the_emission_table_is_the_reference_bit_for_bit() {
    let text = golden();
    for row in rows(&text, "qual2prob") {
        let q: usize = row[0].parse().unwrap();
        assert_eq!(bits(qual_to_prob_table()[q]), row[1], "qual2prob[{q}]");
    }

    let baq = Baq::default();
    let entries = rows(&text, "epsilon");
    assert!(
        entries.len() >= 100,
        "a slice of the table, not a sample of one"
    );
    for row in entries {
        let reference = row[0].as_bytes()[0];
        let read = row[1].as_bytes()[0];
        let q: u8 = row[2].parse().unwrap();
        assert_eq!(
            bits(baq.calc_epsilon(reference, read, q)),
            row[3],
            "epsilon[{}][{}][{q}]",
            row[0],
            row[1]
        );
    }
}

/// The model itself, over the sequences the dump ran.
#[test]
fn every_run_of_the_model_is_the_reference() {
    let text = golden();
    let baq = Baq::default();
    let unique = b"GATTACAGGCTCTAGCAT";
    /// One run: a label, a reference, a query and its encoded qualities.
    type Run = (&'static str, &'static [u8], &'static [u8], &'static [u8]);
    let cases: Vec<Run> = vec![
        ("exact", unique, b"TTACAGGC", b"IIIIIIII"),
        ("mismatch", unique, b"TTACTGGC", b"IIIIIIII"),
        ("repeat", b"ACACACACACACACACACAC", b"ACACAC", b"IIIIII"),
        ("low-quality", unique, b"TTACAGGC", b"!!!!!!!!"),
        ("mixed-quality", unique, b"TTACAGGC", b"!\"#$%&'("),
        ("n-in-read", unique, b"TTANAGGC", b"IIIIIIII"),
        ("n-in-ref", b"GATTACANGCTCTAGCAT", b"TTACAGGC", b"IIIIIIII"),
        ("one-base", unique, b"A", b"I"),
        ("full-length", b"ACGTACGTAC", b"ACGTACGTAC", b"IIIIIIIIII"),
        ("query-longer-than-ref", b"ACGT", b"ACGTACGT", b"IIIIIIII"),
    ];
    let expected = |label: &str| -> Vec<String> {
        rows(&text, "hmm")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no hmm row {label}"))[1..]
            .iter()
            .map(|field| field.to_string())
            .collect()
    };

    for (label, reference, query, encoded) in cases {
        let quals: Vec<u8> = encoded.iter().map(|c| c - 33).collect();
        let mut state = vec![0i32; query.len()];
        let mut q = vec![0u8; query.len()];
        let returned = baq
            .hmm_glocal(
                reference,
                query,
                0,
                query.len() as i32,
                &quals,
                &mut state,
                &mut q,
            )
            .unwrap();
        let theirs = expected(label);
        assert_eq!(returned.to_string(), theirs[0], "{label}: return");
        assert_eq!(join(&q), theirs[1], "{label}: bq");
        assert_eq!(
            state
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(","),
            theirs[2],
            "{label}: state"
        );
    }
    println!("baq: 10 runs of the model compared");
}

#[test]
fn the_state_decoders_are_the_references() {
    let text = golden();
    for row in rows(&text, "state") {
        let value: i32 = row[0].parse().unwrap();
        assert_eq!(state_is_indel(value).to_string(), row[1], "state {value}");
        assert_eq!(
            state_aligned_position(value).to_string(),
            row[2],
            "state {value}"
        );
    }
}

/// The reference window and the per-read calculation, over the cigars that move them.
#[test]
fn every_read_window_and_calculation_is_the_reference() {
    let text = golden();
    let baq = Baq::default();
    let reference: Vec<u8> = "ACGT".repeat(25).into_bytes();
    let cases: Vec<(&str, &str, i32)> = vec![
        ("aligned", "10M", 50),
        ("leading-insertion", "3I7M", 50),
        ("trailing-insertion", "7M3I", 50),
        ("both-insertions", "2I6M2I", 50),
        ("soft-clipped", "3S7M", 50),
        ("deletion", "4M2D6M", 50),
        ("at-contig-start", "10M", 1),
        ("at-contig-end", "10M", 95),
        ("skipped-region", "4M2N4M", 50),
        ("all-soft-clipped", "10S", 50),
    ];

    for (label, cigar, start) in cases {
        let record = read(cigar, start);
        let window = reference_window_for_read(&record, "chr1", DEFAULT_BANDWIDTH).unwrap();
        let expected_window = rows(&text, "window")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no window row {label}"))[1]
            .to_string();
        assert_eq!(
            format!("{}:{}-{}", window.contig, window.start, window.end),
            expected_window,
            "{label}: window"
        );

        let expected_baq = rows(&text, "baq")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no baq row {label}"))[1]
            .to_string();
        let ours = baq.calc_baq_from_hmm(&record, &reference, 0).unwrap();
        match ours {
            Some(result) => assert_eq!(join(&result.bq), expected_baq, "{label}: bq"),
            None => assert_eq!(expected_baq, "null", "{label}: bq"),
        }
    }
}

/// The BQ tag, which carries the difference between the quality and the BAQ.
#[test]
fn the_tag_is_the_reference_in_both_directions() {
    let text = golden();
    let cases: Vec<Vec<u8>> = vec![
        vec![40, 40, 40, 40],
        vec![40, 39, 30, 0],
        vec![41, 45, 40, 40],
    ];
    for (index, values) in cases.iter().enumerate() {
        let quals = vec![40u8; values.len()];
        let tag = encode_bq_tag(&quals, values);
        let expected = rows(&text, "tag")
            .into_iter()
            .find(|row| row[0] == index.to_string())
            .unwrap_or_else(|| panic!("no tag row {index}"))[1]
            .to_string();
        assert_eq!(tag, expected, "tag {index}");

        for suffix in ["", "-overwrite"] {
            let label = format!("{index}{suffix}");
            let expected = rows(&text, "fromtag")
                .into_iter()
                .find(|row| row[0] == label)
                .unwrap_or_else(|| panic!("no fromtag row {label}"))[1]
                .to_string();
            assert_eq!(
                join(&calc_baq_from_tag("t", "chr1:1-4", &quals, Some(&tag), false).unwrap()),
                expected,
                "fromtag {label}"
            );
        }
    }

    // A read with no tag at all, asked for one both ways.
    let quals = vec![40u8; 4];
    let raw = rows(&text, "fromtag")
        .into_iter()
        .find(|row| row[0] == "bare-raw")
        .unwrap()[1]
        .to_string();
    assert_eq!(
        join(&calc_baq_from_tag("bare", "chr1:10-13", &quals, None, true).unwrap()),
        raw
    );
    let message = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "bare-strict")
        .unwrap()[2]
        .to_string();
    assert_eq!(
        calc_baq_from_tag("bare", "chr1:10-13", &quals, None, false)
            .unwrap_err()
            .message(),
        message
    );

    // A tag no encoder produces, whose difference takes the quality below zero.
    let malformed = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "malformed-tag")
        .unwrap()[2]
        .to_string();
    assert_eq!(
        calc_baq_from_tag("malformed", "chr1:10-13", &quals, Some("zzzz"), false)
            .unwrap_err()
            .message(),
        malformed
    );
}
