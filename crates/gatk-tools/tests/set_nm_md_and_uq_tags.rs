//! Conformance for `SetNmMdAndUqTags` against Picard 3.4.0, compared as the tags of every record
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/SetNmMdAndUqTagsDump.java`, which carries the
//! reference, the input BAM as base64 and each output both as base64 and as text.
//!
//! # What this suite is for
//!
//!  * **the three tags of every record**, replaced rather than filled in;
//!  * **an unmapped read left alone**, nonsense tags and all;
//!  * **UQ skipped for a record with no qualities**, which keeps the one it arrived with;
//!  * **`SET_ONLY_UQ`**, which leaves a wrong NM and a wrong MD where they were;
//!  * **the IUPAC comparison**, an `N` in the reference being a mismatch and an `N` in the read
//!    matching it;
//!  * **and the bisulfite run's disagreement**: `MD:Z:0C0C0C0C0C0C0C0C0` beside `NM:i:0`, the two
//!    written by different functions.
//!
//! The comparison is the records' tags rather than the output's bytes: the reference's writer
//! stamps a `@PG` line this port does not write, so the files differ in the header while every
//! record agrees. The text rows the dump carries are what the tags are read from.

use gatk_corpus as corpus;
use gatk_tools::set_nm_md_and_uq_tags::{check_sort_order, fix_record, Arguments, SetTagsError};
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bgzf::read::decompress_all;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/set_nm_md_and_uq_tags.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .to_string()
}

/// The `sam\t<label>=` rows, which are the output's records as text.
fn sam(text: &str, label: &str) -> String {
    let prefix = format!("sam\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries sam/{label}")),
    )
}

/// The reference's contigs, as the dump wrote the fasta.
fn reference(text: &str) -> Vec<(String, Vec<u8>)> {
    let fasta = unescape(
        text.lines()
            .find_map(|line| line.strip_prefix("fasta\t"))
            .expect("the golden carries the reference"),
    );
    let mut contigs: Vec<(String, Vec<u8>)> = Vec::new();
    for line in fasta.lines() {
        if let Some(name) = line.strip_prefix('>') {
            contigs.push((name.to_string(), Vec::new()));
        } else if let Some(last) = contigs.last_mut() {
            last.1.extend_from_slice(line.as_bytes());
        }
    }
    contigs
}

fn fixture(text: &str, label: &str) -> Vec<BamRecord> {
    let plain = decompress_all(&corpus::decode_base64(&row(text, "fixture", label)))
        .expect("a bgzf fixture");
    let reader = BamReader::new(&plain).expect("a bam fixture");
    reader.map(|record| record.expect("a record")).collect()
}

/// The tags of one record, rendered the way the golden's text rows carry them.
fn tags_of(record: &BamRecord) -> Vec<String> {
    let mut rendered: Vec<String> = Vec::new();
    for name in [b"MD", b"NM", b"UQ"] {
        match record.tags.get(Tag::new(name)) {
            Some(TagValue::Str(value)) => rendered.push(format!(
                "{}:Z:{value}",
                std::str::from_utf8(name).expect("a tag name")
            )),
            Some(TagValue::Int(value)) => rendered.push(format!(
                "{}:i:{value}",
                std::str::from_utf8(name).expect("a tag name")
            )),
            _ => {}
        }
    }
    rendered
}

/// The same three tags taken off one line of the golden's text, in the same order.
fn expected_tags(line: &str) -> Vec<String> {
    let mut rendered: Vec<String> = Vec::new();
    for name in ["MD", "NM", "UQ"] {
        if let Some(field) = line
            .split('\t')
            .find(|field| field.starts_with(&format!("{name}:")))
        {
            rendered.push(field.to_string());
        }
    }
    rendered
}

fn run(text: &str, label: &str, arguments: &Arguments) {
    let contigs = reference(text);
    let mut records = fixture(text, "sorted");
    for record in &mut records {
        let bases = if record.reference_index >= 0 {
            contigs[record.reference_index as usize].1.clone()
        } else {
            Vec::new()
        };
        fix_record(record, &bases, arguments);
    }
    let output = sam(text, label);
    let expected: Vec<&str> = output.lines().collect();
    assert_eq!(records.len(), expected.len(), "{label}: record count");
    for (record, line) in records.iter().zip(expected.iter()) {
        assert_eq!(
            line.split('\t').next().expect("a read name"),
            record.read_name,
            "{label}: the records line up"
        );
        assert_eq!(
            tags_of(record),
            expected_tags(line),
            "{label}: {}",
            record.read_name
        );
    }
}

#[test]
fn the_default_run_matches_the_golden() {
    let text = golden();
    run(&text, "defaults", &Arguments::default());
}

#[test]
fn set_only_uq_leaves_nm_and_md_alone() {
    let text = golden();
    run(
        &text,
        "only-uq",
        &Arguments {
            set_only_uq: true,
            ..Arguments::default()
        },
    );
}

#[test]
fn the_bisulfite_run_matches_the_golden() {
    let text = golden();
    run(
        &text,
        "bisulfite",
        &Arguments {
            is_bisulfite_sequence: true,
            ..Arguments::default()
        },
    );
    // The record the two functions disagree about, spelled out: MD lists eight mismatches and NM
    // counts none.
    let line = sam(&text, "bisulfite")
        .lines()
        .find(|line| line.starts_with("bisulfite\t"))
        .expect("the bisulfite record")
        .to_string();
    assert!(line.contains("MD:Z:0C0C0C0C0C0C0C0C0"));
    assert!(line.contains("NM:i:0"));
}

#[test]
fn a_queryname_sorted_input_is_refused() {
    let text = golden();
    let error = check_sort_order(Some("queryname")).expect_err("the sort order refusal");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        unescape(&row(&text, "error", "queryname"))
    );
    assert_eq!(check_sort_order(Some("coordinate")), Ok(()));
    // An absent SO is `unsorted`, which is refused by the same message.
    assert_eq!(
        check_sort_order(None),
        Err(SetTagsError::NotCoordinateSorted {
            found: "unsorted".to_string()
        })
    );
}
