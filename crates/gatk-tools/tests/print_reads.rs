//! Conformance for `PrintReads` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Every suite before this one compared decisions: which reads, which bases, which answer. This
//! one compares the file. The output BAM and its index travel in the golden in full, and the port
//! writes its own and must produce the same bytes.
//!
//! What that claim is about is declared rather than implied: the reference was run with
//! `--use-jdk-deflater`, because GATK's default is the Intel GKL deflater and htsjdk-rs
//! reproduces the JDK one. The command line lands in the `@PG` record's `CL`, so it is read out
//! of the golden and handed to the port rather than reconstructed: it carries the temporary paths
//! of the run that produced it, and inventing it would be inventing part of the file.

use gatk_corpus as corpus;
use gatk_engine::interval::{self, SimpleInterval};
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::{not_duplicate, with_header};
use gatk_tools::print_reads::{self, Options};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/print_reads.txt.gz"),
    )
}

fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .unwrap_or_else(|| panic!("the golden carries no {kind} row"))
}

/// The rows of one kind, by label.
fn rows<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str)> {
    text.lines()
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .filter_map(|rest| rest.split_once('\t'))
        .collect()
}

/// The arguments each labelled run was given. A label is a configuration, and there is nothing in
/// the row to derive it from.
fn configuration(label: &str) -> (Vec<&'static str>, bool, bool, bool) {
    // (intervals, default filters, not-duplicate filter, create index)
    match label {
        "all" => (vec![], true, false, true),
        "chr1" => (vec!["chr1"], true, false, true),
        "chr1:100-160" => (vec!["chr1:100-160"], true, false, true),
        "nofilter" => (vec![], false, false, true),
        "nodup" => (vec![], true, true, true),
        "noindex" => (vec![], true, false, false),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();

    let dir = std::env::temp_dir().join(format!("gatk-rs-printreads-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let bam = dir.join("reads.bam");
    let bai = dir.join("reads.bai");
    std::fs::write(&bam, corpus::decode_base64(field(&text, "bam"))).unwrap();
    std::fs::write(&bai, corpus::decode_base64(field(&text, "bai"))).unwrap();

    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();

    let command_lines: Vec<(&str, &str)> = rows(&text, "commandline");
    let outputs: Vec<(&str, &str)> = rows(&text, "output");
    let indexes: Vec<(&str, &str)> = rows(&text, "index");
    assert_eq!(command_lines.len(), outputs.len());

    let mut compared = 0;
    for (label, expected_base64) in &outputs {
        let (interval_strings, default_filters, no_duplicates, create_index) = configuration(label);
        let intervals: Vec<SimpleInterval> = interval_strings
            .iter()
            .map(|text| interval::parse_interval(text, &header).expect("a parsable interval"))
            .collect();

        let header_for_filter = header.clone();
        let filter: Box<dyn Fn(&BamRecord) -> bool> = if !default_filters {
            Box::new(|_: &BamRecord| true)
        } else if no_duplicates {
            Box::new(move |read: &BamRecord| {
                with_header::wellformed(read, &header_for_filter) && not_duplicate(read)
            })
        } else {
            Box::new(move |read: &BamRecord| with_header::wellformed(read, &header_for_filter))
        };

        let command_line = command_lines
            .iter()
            .find(|(l, _)| l == label)
            .map(|(_, cl)| *cl)
            .expect("a command line for every run");

        let options = Options {
            intervals,
            create_output_bam_index: create_index,
            command_line,
            ..Options::default()
        };
        let (ours, our_index) =
            print_reads::print_reads(&source, &options, filter.as_ref()).expect("the tool runs");

        let expected = corpus::decode_base64(expected_base64);
        assert_eq!(
            ours.len(),
            expected.len(),
            "{label}: output length differs ({} vs {})",
            ours.len(),
            expected.len()
        );
        if ours != expected {
            let at = ours
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            panic!("{label}: first byte difference at offset {at}");
        }

        let expected_index = indexes
            .iter()
            .find(|(l, _)| l == label)
            .map(|(_, value)| *value)
            .expect("an index row for every run");
        match (our_index, expected_index) {
            (None, "absent") => {}
            (Some(_), "absent") => panic!("{label}: the reference wrote no index and the port did"),
            (None, _) => panic!("{label}: the reference wrote an index and the port did not"),
            (Some(ours), expected) => {
                assert_eq!(ours, corpus::decode_base64(expected), "{label}: the .bai");
            }
        }
        compared += 1;
    }

    std::fs::remove_dir_all(&dir).ok();
    println!("{compared} PrintReads outputs, byte-identical with their indexes");
}
