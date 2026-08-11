//! Conformance for `LeftAlignIndels` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/LeftAlignIndelsToolDump.java`. The output BAMs and
//! their indexes travel in full, base64, as the rest of this archetype's do, and so do the fixture
//! with its index and the reference with its `.fai`: this is the first tool here that needs a
//! reference, and a test that built its own would be inventing part of the input.
//!
//! # What this suite is for
//!
//! The eighth whole tool of the archetype, and the first whose `requiresReference()` is true. The
//! call it is built around has its own suite; what is here is what the tool does with the answer:
//!
//!  * **the window is the read, not the contig.** A `ReadWalker` hands `apply` the read's own span
//!    with no padding, so an indel can only move as far as the read reaches. The golden carries
//!    the window of every read, so a port that widened it would fail on the window rather than
//!    only on the bytes;
//!  * **a deletion that walks off the front moves the read**, and an insertion that does the same
//!    does not: `4M1D5M` at `chr1:6` becomes `9M` at `chr1:7`, while `4M2I4M` becomes `2I8M` where
//!    it stood;
//!  * **two classes of read never reach the call**, an unmapped one and one whose cigar has a
//!    single element, and both come out unchanged;
//!  * **a run with no `-R` is refused by Barclay**, not by the engine's `requiresReference`, and
//!    the message is an argument-parsing one.
//!
//! The command line lands in the `@PG` record's `CL`, so it is read out of the golden and handed
//! to the port rather than reconstructed: it carries the paths of the run that produced it.

use gatk_corpus as corpus;
use gatk_engine::interval::{self, SimpleInterval};
use gatk_engine::reads::ReadsDataSource;
use gatk_engine::reference::ReferenceFileSource;
use gatk_readfilter::with_header;
use gatk_tools::left_align_indels as tool;
use gatk_tools::sam_output::Options;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/left_align_indels_tool.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
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

/// A row of one kind that carries a single field.
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

/// What each labelled run was given. A label is a configuration and the row carries nothing to
/// derive it from, so it is written here beside the dump that produced it.
struct Configuration {
    intervals: &'static [&'static str],
    create_index: bool,
    program_record: bool,
}

fn configuration(label: &str) -> Configuration {
    let base = Configuration {
        intervals: &[],
        create_index: true,
        program_record: true,
    };
    match label {
        "all" => base,
        // Drops the two reads that start past base 20, one of which is the unmapped one.
        "chr1head" => Configuration {
            intervals: &["chr1:1-20"],
            ..base
        },
        "noindex" => Configuration {
            create_index: false,
            ..base
        },
        "nopg" => Configuration {
            program_record: false,
            ..base
        },
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The fixture and the reference, written out so the port can open them.
fn install(text: &str, dir: &std::path::Path) {
    std::fs::create_dir_all(dir).expect("a scratch directory");
    std::fs::write(dir.join("ref.fasta"), unescape(field(text, "fasta"))).expect("the fasta");
    std::fs::write(dir.join("ref.fasta.fai"), unescape(field(text, "fai"))).expect("the fai");
    for row in rows(text, "fixture") {
        std::fs::write(
            dir.join(format!("{}.bam", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture bam");
    }
    for row in rows(text, "fixtureindex") {
        std::fs::write(
            dir.join(format!("{}.bai", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture index");
    }
}

#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-leftalign-{}", std::process::id()));
    install(&text, &dir);

    let outputs = rows(&text, "output");
    let indexes = rows(&text, "index");
    assert_eq!(outputs.len(), 4, "four runs finish and one is refused");

    let mut compared = 0usize;
    for row in &outputs {
        let (label, expected_base64) = (row[0], row[1]);
        let config = configuration(label);
        let source = ReadsDataSource::open(&dir.join("plain.bam"), &dir.join("plain.bai"))
            .expect("the fixture opens");
        let mut reference =
            ReferenceFileSource::open(&dir.join("ref.fasta")).expect("the reference opens");
        let header = source.header().clone();

        let intervals: Vec<SimpleInterval> = config
            .intervals
            .iter()
            .map(|text| interval::parse_interval(text, &header).expect("a parsable interval"))
            .collect();

        // This tool does not override its read filters, so its default is GATKTool's.
        let filter = move |read: &BamRecord| with_header::wellformed(read, &header);

        let command_line = of_run(&text, "commandline", label)
            .first()
            .map(|row| row.get(1).copied().unwrap_or(""))
            .unwrap_or("");
        let options = Options {
            intervals,
            create_output_bam_index: config.create_index,
            add_output_sam_program_record: config.program_record,
            command_line,
            ..Options::default()
        };

        let (ours, our_index) = tool::left_align_indels(&source, &mut reference, &options, &filter)
            .expect("the source reads")
            .expect("this label does not refuse");

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

        let expected_index = indexes
            .iter()
            .find(|index| index[0] == label)
            .map(|index| index[1])
            .expect("an index row for every output");
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

    assert_eq!(compared, 4);
    println!("left-align-indels-tool: {compared} output files byte-identical");
}

/// The window `apply` is handed is the read's own span, which is the bound on how far an indel can
/// move.
///
/// A port that queried the contig would left-align further than the reference does and produce a
/// file that looks healthier and is wrong. This asserts the window itself rather than only the
/// bytes that follow from it.
#[test]
fn the_reference_window_is_the_read_and_nothing_more() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-leftalign-win-{}", std::process::id()));
    install(&text, &dir);

    let source = ReadsDataSource::open(&dir.join("plain.bam"), &dir.join("plain.bai"))
        .expect("the fixture opens");
    let mut reference =
        ReferenceFileSource::open(&dir.join("ref.fasta")).expect("the reference opens");

    let applied = gatk_tools::read_walker::traverse_with_reference(
        &source,
        Some(&mut reference),
        &[],
        false,
        &|_: &BamRecord| true,
    )
    .expect("the reads");

    let expected = rows(&text, "refwindow");
    assert_eq!(applied.len(), expected.len(), "one window per read");

    for (mut entry, row) in applied.into_iter().zip(&expected) {
        let name = row[0];
        assert_eq!(entry.read.read_name, name);
        match entry.context.window() {
            // The unmapped read has no interval at all, which the dump prints as a star.
            None => assert_eq!(row[1], "*", "{name}: window"),
            Some(window) => {
                assert_eq!(window.contig, row[1], "{name}: contig");
                assert_eq!(window.start.to_string(), row[2], "{name}: start");
                assert_eq!(window.end.to_string(), row[3], "{name}: end");
            }
        }
        let bases = entry.context.bases(&mut reference).expect("the bases");
        if row[4] != "-" {
            assert_eq!(String::from_utf8_lossy(&bases), row[4], "{name}: bases");
        }
    }
}

/// What the tool did to each read, which is the finding rather than the byte count.
#[test]
fn a_deletion_moves_the_read_and_an_insertion_does_not() {
    let text = golden();
    let expected = of_run(&text, "reads", "all");
    assert_eq!(expected.len(), 6);

    let by_name = |name: &str| -> Vec<String> {
        expected
            .iter()
            .find(|row| row[1] == name)
            .map(|row| row.iter().map(|field| field.to_string()).collect())
            .unwrap_or_else(|| panic!("the golden lost {name}"))
    };

    // r0 went in as 4M1D5M at 6.
    let moved = by_name("r0");
    assert_eq!((moved[4].as_str(), moved[5].as_str()), ("7", "9M"));
    // r4 went in as 4M2I4M at 21.
    let stayed = by_name("r4");
    assert_eq!((stayed[4].as_str(), stayed[5].as_str()), ("21", "2I8M"));
    // r2 is the single-element cigar and r5 the unmapped read: both untouched.
    assert_eq!(
        (by_name("r2")[4].as_str(), by_name("r2")[5].as_str()),
        ("11", "10M")
    );
    assert_eq!(
        (by_name("r5")[4].as_str(), by_name("r5")[5].as_str()),
        ("41", "10M")
    );
    // r3's deletion has no repeat to move into.
    assert_eq!(
        (by_name("r3")[4].as_str(), by_name("r3")[5].as_str()),
        ("17", "4M1D5M")
    );
}

/// The run that never starts, and what refuses it.
#[test]
fn a_run_with_no_reference_is_refused_by_the_argument_parser() {
    let text = golden();
    let refusals = rows(&text, "error");
    assert_eq!(refusals.len(), 1);
    let (label, message) = (refusals[0][0], refusals[0][1]);
    assert_eq!(label, "noreference");
    // Not the engine's requiresReference, and not a UserException: Barclay refuses the command
    // line before the tool is built at all, so the port has nothing to reproduce here beyond
    // recording which layer said no.
    assert!(
        message.starts_with(
            "org.broadinstitute.barclay.argparser.CommandLineException$MissingArgument:"
        ),
        "the refusal changed layer: {message}"
    );
    assert!(message.contains("Argument 'reference' is required"));
}
