//! The first tool this port runs from a command line, end to end.
//!
//! The bytes of the index are compared against the golden by
//! `gatk-tools/tests/index_feature_file.rs`, which hands the ported function the fixture's text.
//! What this suite is for is the layer between that function and a command line: reading the named
//! file, deciding what to call the output, writing it, and returning the path the reference
//! returns.
//!
//! # What this suite is for
//!
//!  * **the file plumbing: a path in, a file out, and the name the reference chose**;
//!  * **the source the header records, which is the input's own path and timestamp**;
//!  * **the refusals reaching the dispatcher as the reference's statuses**;
//!  * **and `handleResult` printing the path the tool returned.**

use gatk_corpus as corpus;
use gatk_tools::index_feature_file::{build, default_output, Source};
use gatk_tools::main_entry::{self, Failure};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/index_feature_file.txt.gz"),
    )
}

fn bytes(text: &str, kind: &str, label: &str) -> Vec<u8> {
    let prefix = format!("{kind}\t{label}\t");
    let encoded = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
        .unwrap_or_else(|| panic!("{kind}/{label}"));
    corpus::decode_base64(&encoded)
}

fn args(list: &[&str]) -> Vec<String> {
    list.iter().map(|arg| (*arg).to_string()).collect()
}

/// A directory of this test's own, named after the case so two cases cannot collide.
fn scratch(case: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("gatk-cli-{case}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    dir
}

/// The tool writes the index the port builds, under the name the reference chose.
#[test]
fn the_index_is_written_beside_its_input() {
    let text = golden();
    let dir = scratch("plain");
    let input = dir.join("reads.vcf");
    std::fs::write(&input, bytes(&text, "input", "plain")).expect("the fixture");
    let path = input.to_string_lossy().to_string();

    let written = gatk_cli::run(&args(&["IndexFeatureFile", "-I", &path]));
    assert_eq!(written.status, 0, "{}", written.stderr);
    // `handleResult` prints what the tool returned, which is the index's path.
    let expected_output = default_output(&path);
    assert_eq!(
        written.stdout,
        format!("Tool returned:\n{expected_output}\n")
    );
    // Which is the same name the golden recorded, once the dump's directory is put back.
    let recorded = golden()
        .lines()
        .find(|line| line.starts_with("returned\tplain\t"))
        .map(|line| line["returned\tplain\t".len()..].to_string())
        .expect("the returned path");
    assert_eq!(recorded, "<dir>/reads.vcf.idx");
    assert!(expected_output.ends_with("reads.vcf.idx"));

    // The file exists, and it is the bytes the port builds for that source: the path and the
    // timestamp the header records are the input file's own, which is what the plumbing is for.
    let produced = std::fs::read(&expected_output).expect("the index");
    let mut source = Source::new(&path);
    source.timestamp = std::fs::metadata(&input)
        .and_then(|data| data.modified())
        .ok()
        .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|since| since.as_millis() as i64)
        .expect("a timestamp");
    let text_in = String::from_utf8(bytes(&text, "input", "plain")).expect("a text fixture");
    let expected = build(&text_in, &source, &path).expect("the index");
    assert_eq!(produced, expected);
    // And the timestamp is NOT zero, which is what says the plumbing read it rather than leaving
    // the default the suite's own fixtures use.
    assert!(source.timestamp > 0);
}

/// An explicit output is used as it is given.
#[test]
fn an_explicit_output_is_used_as_given() {
    let text = golden();
    let dir = scratch("explicit");
    let input = dir.join("reads.vcf");
    std::fs::write(&input, bytes(&text, "input", "plain")).expect("the fixture");
    let output = dir.join("elsewhere.idx");
    let written = gatk_cli::run(&args(&[
        "IndexFeatureFile",
        "-I",
        &input.to_string_lossy(),
        "-O",
        &output.to_string_lossy(),
    ]));
    assert_eq!(written.status, 0, "{}", written.stderr);
    assert!(output.exists());
    assert!(!dir.join("reads.vcf.idx").exists());
    assert!(written.stdout.ends_with("elsewhere.idx\n"));
}

/// The second tool the port runs prints its report rather than writing a file.
#[test]
fn a_report_is_returned_where_a_file_is_not_named() {
    let dir = scratch("bgzf");
    // A tiny BGZF file: an empty block and the terminator, which is what an empty gzip member
    // and the reference's own terminator look like on disk.
    let terminator: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b,
        0x00, 0x03, 0x00, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    let file = dir.join("empty.bgzf");
    std::fs::write(&file, terminator).expect("the fixture");
    let written = gatk_cli::run(&args(&[
        "PrintBGZFBlockInformation",
        "--bgzf-file",
        &file.to_string_lossy(),
    ]));
    assert_eq!(written.status, 0, "{}", written.stderr);
    // `handleResult` prints what the tool returned, and what this tool returns with no `--output`
    // is the report itself.
    assert!(
        written
            .stdout
            .starts_with("Tool returned:\nBGZF block information for file: empty.bgzf"),
        "{}",
        written.stdout
    );
    // With an output it writes the file instead and returns nothing.
    let report = dir.join("report.txt");
    let written = gatk_cli::run(&args(&[
        "PrintBGZFBlockInformation",
        "--bgzf-file",
        &file.to_string_lossy(),
        "--output",
        &report.to_string_lossy(),
    ]));
    assert_eq!(written.status, 0, "{}", written.stderr);
    assert!(written.stdout.is_empty(), "{}", written.stdout);
    let text = std::fs::read_to_string(&report).expect("the report");
    assert!(
        text.starts_with("BGZF block information for file: empty.bgzf"),
        "{text}"
    );
    // And a file that is not block compressed is refused by the framing check rather than by a
    // codec search: sixteen bytes of text have a header the reader can read and cannot believe.
    let plain = dir.join("plain.txt");
    std::fs::write(&plain, b"not a bgzf file\n").expect("the fixture");
    let refused = gatk_cli::run(&args(&[
        "PrintBGZFBlockInformation",
        "--bgzf-file",
        &plain.to_string_lossy(),
    ]));
    assert_eq!(refused.status, main_entry::exit_status(Failure::User));
    // The reference asks whether the file is block compressed BEFORE it walks it, so the refusal
    // is about what the file is and not about how its first bytes frame. A covering array run
    // against the binary is what found the port answering the other one.
    assert!(
        refused.stderr.contains("File is not a valid BGZF file"),
        "{}",
        refused.stderr
    );
}

/// A refusal reaches the dispatcher as the reference's own status and message.
#[test]
fn a_refusal_reaches_the_dispatcher() {
    let text = golden();
    let dir = scratch("absent");
    let missing = dir.join("absent.vcf");
    let written = gatk_cli::run(&args(&[
        "IndexFeatureFile",
        "-I",
        &missing.to_string_lossy(),
    ]));
    // `CouldNotReadInputFile` is a UserException, which is status two and not one.
    assert_eq!(written.status, main_entry::exit_status(Failure::User));
    let recorded = text
        .lines()
        .find(|line| line.starts_with("error\tabsent\t"))
        .expect("the golden's refusal");
    let message = recorded
        .split_once("CouldNotReadInputFile:")
        .expect("the message")
        .1;
    assert!(
        written
            .stderr
            .contains(&message.replace("<dir>", &dir.to_string_lossy())),
        "{}",
        written.stderr
    );
    // A file that exists and has no codec is a different refusal with a different message.
    let notes = dir.join("notes.txt");
    std::fs::write(&notes, b"nothing a codec knows\n").expect("the fixture");
    let written = gatk_cli::run(&args(&["IndexFeatureFile", "-I", &notes.to_string_lossy()]));
    assert_eq!(written.status, main_entry::exit_status(Failure::User));
    assert!(
        written.stderr.contains("because no suitable codecs found"),
        "{}",
        written.stderr
    );
    // And a command line the parser refuses is status one, before any file is opened.
    let refused = gatk_cli::run(&args(&["IndexFeatureFile", "--no-such-argument", "1"]));
    assert_eq!(
        refused.status,
        main_entry::exit_status(Failure::CommandLine)
    );
}
