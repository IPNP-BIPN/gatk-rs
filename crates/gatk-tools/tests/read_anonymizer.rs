//! Conformance for `ReadAnonymizer` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/ReadAnonymizerDump.java`. The output BAMs and their
//! indexes travel in full, base64, as every record-transform tool's do, and so does the input
//! fixture, its index and the reference.
//!
//! # What this suite is for
//!
//! The ninth whole tool of the archetype, and the first whose transform makes the read a different
//! length:
//!
//!  * **a deletion adds the reference bases to the read** and an insertion drops its own;
//!  * **every M, X and D becomes one operator**, so elements of different kinds collapse;
//!  * **a matching base keeps its quality and a replaced one takes `--ref-base-quality`**;
//!  * **every attribute but the read group is cleared**;
//!  * **its default read filters are a sixth pattern**, and the first without `WellformedReadFilter`.

use gatk_corpus as corpus;
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::with_header;
use gatk_tools::read_anonymizer::{self, AnonymizerArguments};
use gatk_tools::sam_output::Options;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_anonymizer.txt.gz"),
    )
}

fn pairs<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str)> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .filter_map(|rest| rest.split_once('\t'))
        .collect()
}

fn value<'a>(text: &'a str, kind: &str, label: &str) -> Option<&'a str> {
    pairs(text, kind)
        .into_iter()
        .find(|(name, _)| *name == label)
        .map(|(_, value)| value)
}

/// The arguments and options each labelled run was given.
struct Configuration {
    arguments: AnonymizerArguments,
    create_index: bool,
    program_record: bool,
}

fn configuration(label: &str) -> Configuration {
    let base = Configuration {
        arguments: AnonymizerArguments::default(),
        create_index: true,
        program_record: true,
    };
    match label {
        "plain" => base,
        "simple-cigar" => Configuration {
            arguments: AnonymizerArguments {
                use_simple_cigar: true,
                ..AnonymizerArguments::default()
            },
            ..base
        },
        "ref-qual-0" => Configuration {
            arguments: AnonymizerArguments {
                ref_base_quality: 0,
                ..AnonymizerArguments::default()
            },
            ..base
        },
        "ref-qual-60" => base,
        "no-index" => Configuration {
            create_index: false,
            ..base
        },
        "no-program-record" => Configuration {
            program_record: false,
            ..base
        },
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The seven default filters, of which `WellformedReadFilter` is **not** one.
fn default_filter(header: &SamHeader) -> impl Fn(&BamRecord) -> bool + '_ {
    move |read: &BamRecord| {
        gatk_readfilter::valid_alignment_start(read)
            && gatk_readfilter::valid_alignment_end(read)
            && gatk_readfilter::read_length_equals_cigar_length(read)
            && gatk_readfilter::seq_is_stored(read)
            && gatk_readfilter::matching_bases_and_quals(read)
            && gatk_readfilter::mapped(read)
            && with_header::alignment_agrees_with_header(read, header)
    }
}

fn install_fixtures(text: &str, dir: &std::path::Path) {
    std::fs::create_dir_all(dir).expect("a scratch directory");
    for (label, encoded) in pairs(text, "fixture") {
        std::fs::write(
            dir.join(format!("{label}.bam")),
            corpus::decode_base64(encoded),
        )
        .expect("the fixture bam");
    }
    for (label, encoded) in pairs(text, "fixtureindex") {
        std::fs::write(
            dir.join(format!("{label}.bai")),
            corpus::decode_base64(encoded),
        )
        .expect("the fixture index");
    }
}

#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-anonymizer-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let reference = text
        .lines()
        .find_map(|line| line.strip_prefix("reference\t"))
        .expect("the golden carries the reference");

    let outputs = pairs(&text, "output");
    assert_eq!(outputs.len(), 6, "six runs finished");
    assert_eq!(pairs(&text, "error").len(), 1, "and one was refused");

    let mut compared = 0usize;
    for (label, expected_base64) in &outputs {
        let config = configuration(label);
        let source = ReadsDataSource::open(&dir.join("input.bam"), &dir.join("input.bai"))
            .expect("the fixture opens");
        let header = source.header().clone();
        let filter = default_filter(&header);

        let options = Options {
            create_output_bam_index: config.create_index,
            add_output_sam_program_record: config.program_record,
            command_line: value(&text, "commandline", label)
                .expect("the golden carries the command line"),
            ..Options::default()
        };

        let (out, index) = read_anonymizer::read_anonymizer(
            &source,
            reference.as_bytes(),
            &config.arguments,
            &options,
            &filter,
        )
        .expect("the tool runs");

        let expected = corpus::decode_base64(expected_base64);
        assert_eq!(out.len(), expected.len(), "{label}: output length differs");
        if out != expected {
            let at = out
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            panic!("{label}: first byte difference at offset {at}");
        }

        let expected_index = value(&text, "index", label).expect("an index row");
        match index {
            Some(index) => assert_eq!(
                index,
                corpus::decode_base64(expected_index),
                "{label}: the .bai"
            ),
            None => assert_eq!(expected_index, "absent", "{label}: the index"),
        }
        compared += 1;
    }
    println!("read-anonymizer: {compared} outputs compared byte for byte");
}

/// The one argument the parser refuses, before the tool ever runs.
#[test]
fn a_reference_quality_above_sixty_is_refused_by_the_parser() {
    let text = golden();
    let row = pairs(&text, "error")
        .into_iter()
        .next()
        .expect("the golden lost the refusal");
    assert_eq!(row.0, "ref-qual-61");
    let (exception, message) = row.1.split_once('\t').expect("an error row has a message");
    assert_eq!(exception, "OutOfRangeArgumentValue");
    assert!(
        message.contains("ref-base-quality"),
        "the message names the argument: {message}"
    );
    // The port takes a `u8` and leaves the range check to the caller, which is where the reference
    // puts it too: this is Barclay refusing the value, not the tool.
    assert_eq!(
        AnonymizerArguments::default().ref_base_quality,
        gatk_tools::read_anonymizer::DEFAULT_REF_BASE_QUALITY
    );
}
