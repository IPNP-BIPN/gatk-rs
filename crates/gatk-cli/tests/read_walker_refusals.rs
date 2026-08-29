//! Conformance for a read walker's refusals against GATK 4.6.2.0.
//!
//! Golden from `tools/argument-conformance/ReadWalkerRefusalDump.java`: eleven inputs handed to
//! `CountReads`, with the exception class, whether it is a `UserException`, and the message.
//!
//! # What this suite is for
//!
//!  * **the three refusals, which are three classes with two statuses between them**;
//!  * **the same file refused differently once an interval is given**;
//!  * **the two inputs that are NOT refused: a BAM, and an empty file**;
//!  * **and the dispatcher answering with the reference's status, which is three for a refusal
//!    that is nobody's fault.**

use gatk_corpus as corpus;
use gatk_tools::main_entry::{self, Failure};
use gatk_tools::read_walker_refusal::{cannot_read, refusal, Refusal};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/read_walker_refusals.txt.gz"),
    )
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
}

fn scratch(case: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("gatk-cli-refusal-{case}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    dir
}

/// The four files the corpus carries, as the fixture program writes them.
fn vcf() -> String {
    let mut text = String::from("##fileformat=VCFv4.2\n");
    text.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    text.push_str("##contig=<ID=chr1,length=100000>\n");
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n");
    let mut position = 100;
    while position <= 5000 {
        text.push_str(&format!(
            "chr1\t{position}\trs{position}\tA\tC\t100\tPASS\t.\tGT\t0/1\n"
        ));
        position += 700;
    }
    text
}

fn bed() -> String {
    let mut text = String::new();
    let mut start = 100;
    while start <= 5000 {
        text.push_str(&format!("chr1\t{start}\t{}\tregion{start}\n", start + 50));
        start += 700;
    }
    text
}

/// Each case's refusal, class, kind and message, against the golden's.
#[test]
fn the_refusals_are_the_goldens() {
    let text = golden();
    let dir = scratch("files");
    let plain = dir.join("reads.vcf");
    std::fs::write(&plain, vcf()).expect("the fixture");
    let regions = dir.join("regions.bed");
    std::fs::write(&regions, bed()).expect("the fixture");
    let notes = dir.join("notes.txt");
    std::fs::write(&notes, "not a bam\n").expect("the fixture");
    let empty = dir.join("empty.bam");
    std::fs::write(&empty, []).expect("the fixture");
    let missing = dir.join("nowhere.bam");
    let directory = dir.join("adirectory");
    std::fs::create_dir_all(&directory).expect("a directory");

    // A file that is not a BAM, with the name the golden's message carries.
    let case =
        |case: &str, path: &std::path::Path, name: &str, compressed: bool, intervals: bool| {
            let bytes = std::fs::read(path).ok();
            let produced = refusal(
                name,
                path.exists(),
                path.is_dir(),
                bytes.as_deref(),
                compressed,
                intervals,
            );
            match field(&text, "class", case)
                .expect("the golden's class")
                .as_str()
            {
                "none" => assert!(produced.is_none(), "{case}"),
                class => {
                    let produced = produced.unwrap_or_else(|| panic!("{case}"));
                    assert_eq!(produced.exception(), class, "{case}");
                    assert_eq!(
                        field(&text, "kind", case).as_deref(),
                        Some(if produced.is_user() { "user" } else { "other" }),
                        "{case}"
                    );
                    assert_eq!(
                        field(&text, "message", case).as_deref(),
                        Some(produced.message().as_str()),
                        "{case}"
                    );
                }
            }
        };

    case("a-plain-vcf", &plain, "fixtures/reads.vcf", false, false);
    case("a-bed", &regions, "fixtures/regions.bed", false, false);
    case("a-text-file", &notes, "<dir>/notes.txt", false, false);
    case("an-empty-file", &empty, "<dir>/empty.bam", false, false);
    case("a-directory", &directory, "<dir>/adirectory", false, false);
    case(
        "a-path-that-does-not-exist",
        &missing,
        "<dir>/nowhere.bam",
        false,
        false,
    );
    // The same two files with an interval, which is a different refusal from the same bytes.
    case(
        "a-plain-vcf-with-an-interval",
        &plain,
        "fixtures/reads.vcf",
        false,
        true,
    );
    case(
        "an-empty-file-with-an-interval",
        &empty,
        "<dir>/empty.bam",
        false,
        true,
    );
    // And the block-compressed one, whose message has lost the file's name by then.
    let produced = refusal(
        "fixtures/reads.vcf.gz",
        true,
        false,
        Some(vcf().as_bytes()),
        true,
        false,
    )
    .expect("the refusal");
    assert_eq!(
        field(&text, "message", "a-block-compressed-vcf").as_deref(),
        Some(produced.message().as_str())
    );
    assert!(matches!(produced, Refusal::NotSamText { file: None, .. }));
}

/// The dispatcher turns each into the status the reference exits with.
#[test]
fn the_statuses_are_the_reference_ones() {
    let dir = scratch("statuses");
    let plain = dir.join("reads.vcf");
    std::fs::write(&plain, vcf()).expect("the fixture");

    // A `SAMFormatException` is nobody's fault, which is status THREE and not two.
    let run = gatk_cli::run(&[
        "CountReads".to_string(),
        "--input".to_string(),
        plain.to_string_lossy().to_string(),
    ]);
    assert_eq!(run.status, main_entry::exit_status(Failure::Other));
    assert!(
        run.stderr.contains("Error parsing text SAM file"),
        "{}",
        run.stderr
    );

    // The same file with an interval is the dictionary's refusal, and also status three.
    let run = gatk_cli::run(&[
        "CountReads".to_string(),
        "--input".to_string(),
        plain.to_string_lossy().to_string(),
        "-L".to_string(),
        "chr1:1-1000".to_string(),
    ]);
    assert_eq!(run.status, main_entry::exit_status(Failure::Other));
    assert!(
        run.stderr.contains("Dictionary cannot have size zero"),
        "{}",
        run.stderr
    );

    // A path that does not exist is the user's, which is status two.
    let run = gatk_cli::run(&[
        "CountReads".to_string(),
        "--input".to_string(),
        dir.join("nowhere.bam").to_string_lossy().to_string(),
    ]);
    assert_eq!(run.status, main_entry::exit_status(Failure::User));
    assert!(
        run.stderr
            .contains("Cannot read non-existent file: file://"),
        "{}",
        run.stderr
    );
    // And htsjdk's wording is the wording, trailing slash for a directory included.
    assert!(cannot_read("/x", true).ends_with("file:///x/"));
    assert!(cannot_read("/x", false).ends_with("file:///x"));
}
