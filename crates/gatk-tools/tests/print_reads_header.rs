//! Conformance for `PrintReadsHeader` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/PrintReadsHeaderDump.java`. The tool's whole output is
//! a text file, so it travels twice: as base64, because the reference writes through an
//! `OutputStreamWriter` with no charset and only the bytes say which one was used, and as escaped
//! text, so a reader can see that nothing was appended without decoding base64 by hand.
//!
//! # What this suite is for
//!
//! The smallest tool of the archetype, taken as much for the calibration gate as for itself:
//!
//!  * **the header it prints is not the header it was given.** `encode(writer, header)` is the
//!    two-argument overload, `keepExistingVersionNumber = false`, so `VN` becomes the current
//!    version and leads the `@HD` line whatever the file said;
//!  * **neither of those is observable on a file htsjdk wrote**, because the ordinary BAM writer
//!    goes through the same overload. The golden's `builtvn` row is the evidence: a header holding
//!    `VN:1.5` when the writer was handed it produces a BAM that reads `VN:1.6`. That is
//!    htsjdk-rs#164, not this tool;
//!  * **nothing is appended**, because the tool reads `getHeaderForReads()` rather than
//!    `getHeaderForSAMWriter()`. Every other tool of this archetype does the opposite;
//!  * **a header with no sequence dictionary still prints**, and a run with no `-I` never starts.

use gatk_corpus as corpus;
use gatk_engine::reads::ReadsDataSource;
use gatk_tools::print_reads_header as tool;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/print_reads_header.txt.gz"),
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

/// Which fixture each labelled run was given.
fn fixture_for(label: &str) -> &'static str {
    match label {
        "oldversion" => "old_version",
        "rich" => "rich",
        "bare" => "bare",
        other => panic!("{other} is in the golden but not configured here"),
    }
}

fn install(text: &str, dir: &std::path::Path) {
    std::fs::create_dir_all(dir).expect("a scratch directory");
    for row in rows(text, "fixture") {
        std::fs::write(
            dir.join(format!("{}.bam", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture bam");
    }
}

#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-readsheader-{}", std::process::id()));
    install(&text, &dir);

    let outputs = rows(&text, "output");
    assert_eq!(outputs.len(), 3, "three runs finish and one is refused");

    let mut compared = 0usize;
    for row in &outputs {
        let (label, expected_base64) = (row[0], row[1]);
        let source =
            ReadsDataSource::open_unindexed(&dir.join(format!("{}.bam", fixture_for(label))))
                .expect("the fixture opens with no index");

        let ours = tool::print_reads_header(&source);
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

        // The same thing as text, which is what a reader of a failure would want to see.
        let expected_text = of_run(&text, "outputtext", label)
            .first()
            .map(|row| unescape(row[1]))
            .expect("an outputtext row for every output");
        assert_eq!(String::from_utf8_lossy(&ours), expected_text, "{label}");

        compared += 1;
    }

    assert_eq!(compared, 3);
    println!("print-reads-header: {compared} output files byte-identical");
}

/// The finding the measurement turned up, which belongs to htsjdk-rs rather than to this tool.
///
/// If htsjdk-rs#164 is ever fixed and this assertion starts failing, the golden is the thing to
/// re-derive: it would mean the writer stopped normalising, which nothing in this repository can
/// decide on its own.
#[test]
fn no_bam_htsjdk_wrote_carries_a_non_current_version() {
    let text = golden();
    let current = field(&text, "currentversion");
    assert_eq!(current, "1.6");

    // The header held 1.5 at the moment the writer was handed it.
    let built = rows(&text, "builtvn")
        .into_iter()
        .find(|row| row[0] == "old_version.bam")
        .map(|row| row[1].to_string())
        .expect("the golden lost its builtvn row");
    assert_eq!(built, "1.5");

    // And every header read back out of a fixture says 1.6.
    for row in rows(&text, "inputheader") {
        let first = unescape(row[1])
            .lines()
            .next()
            .unwrap_or_default()
            .to_string();
        assert!(
            first.starts_with(&format!("@HD\tVN:{current}")),
            "{}: {first}",
            row[0]
        );
    }
}

/// The tool adds nothing, which is what separates it from every other tool of this archetype.
#[test]
fn the_program_chain_is_the_files_own() {
    let text = golden();
    for label in ["oldversion", "rich", "bare"] {
        let output = of_run(&text, "outputtext", label)
            .first()
            .map(|row| unescape(row[1]))
            .unwrap_or_else(|| panic!("the golden lost {label}"));
        let programs: Vec<&str> = output
            .lines()
            .filter(|line| line.starts_with("@PG"))
            .collect();
        assert_eq!(
            programs,
            vec![
                "@PG\tID:upstream\tVN:1.0\tCL:upstream --in a.bam",
                "@PG\tID:downstream\tPP:upstream",
            ],
            "{label}"
        );
        assert!(!output.contains("PrintReadsHeader"), "{label}");
    }
}

/// A header with no sequence dictionary prints, and a run with no input never starts.
#[test]
fn an_empty_dictionary_prints_and_a_missing_input_does_not() {
    let text = golden();
    let bare = of_run(&text, "outputtext", "bare")
        .first()
        .map(|row| unescape(row[1]))
        .expect("the golden lost the bare run");
    assert_eq!(
        bare.lines().filter(|line| line.starts_with("@SQ")).count(),
        0,
        "the codec iterates an empty list rather than refusing"
    );
    assert!(bare.starts_with("@HD\t"));

    let refusals = rows(&text, "error");
    assert_eq!(refusals.len(), 1);
    assert_eq!(refusals[0][0], "noinput");
    // Barclay, before the tool is built at all, rather than the engine's requiresReads.
    assert!(refusals[0][1]
        .starts_with("org.broadinstitute.barclay.argparser.CommandLineException$MissingArgument:"));
}

/// The charset the reference wrote in, checked against the bytes rather than assumed.
#[test]
fn the_output_is_utf8_because_that_is_the_containers_default() {
    let text = golden();
    assert_eq!(field(&text, "charset"), "UTF-8");
    let rich = corpus::decode_base64(
        of_run(&text, "output", "rich")
            .first()
            .map(|row| row[1])
            .expect("the rich output"),
    );
    // The fixture's non-ASCII comment, as UTF-8 rather than as anything else.
    assert!(
        rich.windows(2).any(|w| w == "é".as_bytes()),
        "the accented comment did not survive as UTF-8"
    );
    assert!(String::from_utf8(rich).is_ok());
}
