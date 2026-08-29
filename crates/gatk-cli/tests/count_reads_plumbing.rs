//! Conformance for the file plumbing under `CountReads` against GATK 4.6.2.0.
//!
//! `count-reads-and-bases` compares the number. This suite is the layer between that number and a
//! command line: what `gatk-rs CountReads` returns, what `-O` receives, and what a missing input
//! earns.
//!
//! The fixture is the shared corpus's BAM, rebuilt here rather than read from disk: the golden's
//! cases were run over eight reads on one contig, seven hundred bases apart, one of them flagged a
//! duplicate, and what the port needs is those reads rather than that file.
//!
//! # What this suite is for
//!
//!  * **the tool returning the count, which `handleResult` prints**;
//!  * **`-O` receiving the digits and nothing else, with no trailing newline**;
//!  * **`-O` not suppressing the return**;
//!  * **`-L` counting what the interval holds**;
//!  * **`--read-filter` adding to the defaults and `--disable-tool-default-read-filters`
//!    replacing them**;
//!  * **and a missing input refused in htsjdk's own words.**

use gatk_corpus as corpus;
use gatk_tools::main_entry::{self, Failure};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/count_reads_plumbing.txt.gz"),
    )
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
}

fn args(list: &[&str]) -> Vec<String> {
    list.iter().map(|arg| (*arg).to_string()).collect()
}

/// A directory of this test's own, named after the case.
fn scratch(case: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("gatk-cli-count-{case}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    dir
}

/// The shared fixture's BAM, as `tools/coverage/MakeFixtures.java` writes it: eight reads on one
/// contig seven hundred bases apart, the last of them flagged a duplicate, in one read group.
///
/// Rebuilt here rather than read from disk. The bytes of the file are the writer's business and
/// are compared in htsjdk-rs; what this suite needs is the reads the golden's cases were run over.
fn fixture(dir: &std::path::Path) -> std::path::PathBuf {
    use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
    use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
    use htsjdk_bam::record::BamRecord;
    use htsjdk_bam::tag::{Tag, TagValue, Tags};
    use htsjdk_bam::writer::BamWriter;

    let mut header = SamHeader::default();
    header.sequences.push(SequenceRecord::new("chr1", 100_000));
    let mut group = ReadGroup::new("rg1");
    for (key, value) in [
        ("SM", "sample1"),
        ("LB", "lib1"),
        ("PU", "unit1"),
        ("PL", "ILLUMINA"),
    ] {
        group.attributes.set(key, value);
    }
    header.read_groups.push(group);

    let mut writer = BamWriter::new(Vec::new(), &header)
        .expect("the header")
        .with_index();
    for index in 0..8 {
        // The read group is what `WellformedReadFilter` asks for: a read without it is filtered
        // out, and the tool's default filter is that one.
        let mut tags = Tags::new();
        tags.insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
        let record = BamRecord {
            tags,
            read_name: format!("HWI:1:FC:1:1:{}:{}", index + 1, index + 1),
            flags: if index == 7 { 0x400 } else { 0 },
            reference_index: 0,
            alignment_start: 100 + index * 700,
            mapping_quality: 60,
            cigar: Cigar::new(vec![CigarElement {
                op: Op::M,
                length: 10,
            }]),
            read_bases: b"ACGTACGTAC".to_vec(),
            base_qualities: vec![40; 10],
            ..BamRecord::default()
        };
        writer.write(&record).expect("the record");
    }
    let (bam, bai) = writer.finish_with_index().expect("the fixture");
    let path = dir.join("reads.bam");
    std::fs::write(&path, bam).expect("the fixture");
    std::fs::write(dir.join("reads.bam.bai"), bai).expect("the index");
    path
}

/// Every case the golden ran, through the dispatcher.
#[test]
fn the_plumbing_answers_what_the_reference_answered() {
    let text = golden();
    let dir = scratch("plumbing");
    let bam = fixture(&dir);
    let input = bam.to_string_lossy().to_string();

    let cases: Vec<(&str, Vec<String>)> = vec![
        ("the-whole-file", args(&["--input", &input, "--output", ""])),
        ("no-output-file", args(&["--input", &input])),
        (
            "an-interval",
            args(&["--input", &input, "--output", "", "-L", "chr1:1-1000"]),
        ),
        (
            "an-interval-with-nothing-in-it",
            args(&["--input", &input, "--output", "", "-L", "chr1:50000-60000"]),
        ),
        (
            "a-second-filter",
            args(&[
                "--input",
                &input,
                "--output",
                "",
                "--read-filter",
                "NotDuplicateReadFilter",
            ]),
        ),
        (
            "the-defaults-disabled",
            args(&[
                "--input",
                &input,
                "--output",
                "",
                "--disable-tool-default-read-filters",
            ]),
        ),
        (
            "a-mapping-quality-filter",
            args(&[
                "--input",
                &input,
                "--output",
                "",
                "--read-filter",
                "MappingQualityReadFilter",
                "--minimum-mapping-quality",
                "70",
            ]),
        ),
    ];

    for (case, argv) in cases {
        let output = dir.join(format!("{case}.txt"));
        let mut argv = argv;
        if let Some(position) = argv.iter().position(|value| value.is_empty()) {
            argv[position] = output.to_string_lossy().to_string();
        }
        let mut whole = vec!["CountReads".to_string()];
        whole.extend(argv.clone());
        let run = gatk_cli::run(&whole);
        assert_eq!(run.status, 0, "{case}: {}", run.stderr);

        // The tool returns the count, which `handleResult` prints under its own line.
        let returned = field(&text, "returned", case).unwrap_or_else(|| panic!("{case}"));
        assert_eq!(
            run.stdout,
            format!("Tool returned:\n{returned}\n"),
            "{case}"
        );

        // And the file holds the digits and nothing else, where one was named.
        match field(&text, "file", case) {
            None => assert!(!output.exists(), "{case}"),
            Some(encoded) => {
                let expected = corpus::decode_base64(&encoded);
                let produced = std::fs::read(&output).unwrap_or_else(|_| panic!("{case}"));
                assert_eq!(produced, expected, "{case}");
                // Which is the number with no trailing newline: a `println` would be one byte
                // longer, and the golden says it is not.
                assert_eq!(produced, returned.as_bytes());
            }
        }
    }
}

/// The two refusals, in the reference's own words and with its own statuses.
#[test]
fn the_refusals_are_the_reference_ones() {
    let text = golden();
    let dir = scratch("refusals");

    let missing = dir.join("nowhere.bam");
    let run = gatk_cli::run(&args(&[
        "CountReads",
        "--input",
        &missing.to_string_lossy(),
    ]));
    assert_eq!(run.status, main_entry::exit_status(Failure::User));
    let recorded = field(&text, "error", "an-absent-input").expect("the golden's refusal");
    let message = recorded.split_once(':').expect("the message").1;
    // The golden's path is the container's; what is compared is the wording around it.
    assert!(
        message.starts_with("Cannot read non-existent file: file://"),
        "{message}"
    );
    assert!(
        run.stderr
            .contains("Cannot read non-existent file: file://"),
        "{}",
        run.stderr
    );

    // No input at all is the parser's refusal, which is status one and not two.
    let run = gatk_cli::run(&args(&["CountReads"]));
    assert_eq!(run.status, main_entry::exit_status(Failure::CommandLine));
    let recorded = field(&text, "error", "no-input-at-all").expect("the golden's refusal");
    assert!(
        recorded.ends_with("Argument input was missing: Argument 'input' is required"),
        "{recorded}"
    );
    assert!(
        run.stderr
            .contains("Argument input was missing: Argument 'input' is required"),
        "{}",
        run.stderr
    );
}
