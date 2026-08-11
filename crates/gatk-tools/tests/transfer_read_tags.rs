//! Conformance for `TransferReadTags` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/TransferReadTagsDump.java`. The output BAMs travel in
//! full, base64, as the rest of this archetype's do. No index travels with any of them, on either
//! side: a queryname-sorted BAM cannot have one, and the reference writes none.
//!
//! # What this suite is for
//!
//! The tenth whole tool of the archetype, and the first that is not a walker:
//!
//!  * **the traversal is the tool's own**, so no read filter runs. The fixture's reads carry no
//!    read group, which `WellformedReadFilter` rejects, and every one of them comes out. A port
//!    that reached for the read walker would produce a shorter file;
//!  * **every tag is transferred as a string.** `XI:i:42` arrives as `XI:Z:42` and `XN:f:1.5` as
//!    `XN:Z:1.5`, so the golden prints the Java class of every tag value beside it: the characters
//!    alone would not show the change;
//!  * **an aligned read past the end of the unmapped file is silently dropped**, with no exception
//!    and no warning, because the catch-up loop is bounded by `hasNext()`;
//!  * **the writer is not told the reads are sorted**, so a name whose two records sit in the file
//!    in the order the queryname comparator swaps comes out the other way round;
//!  * **the refusals come from four layers**, and the golden carries the message of each.
//!
//! The command line lands in the `@PG` record's `CL`, so it is read out of the golden and handed to
//! the port rather than reconstructed: it carries the paths of the run that produced it.

use gatk_corpus as corpus;
use gatk_engine::reads::ReadsDataSource;
use gatk_tools::sam_output::Options;
use gatk_tools::transfer_read_tags::{self as tool, TransferError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/transfer_read_tags.txt.gz"),
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

fn of_run<'a>(text: &'a str, kind: &str, label: &str) -> Vec<Vec<&'a str>> {
    rows(text, kind)
        .into_iter()
        .filter(|row| row[0] == label)
        .collect()
}

/// What each labelled run was given: the aligned fixture, the unmapped one, and the tags asked for.
fn configuration(label: &str) -> (&'static str, &'static str, &'static [&'static str]) {
    match label {
        "rx" => ("aligned", "unmapped", &["RX"]),
        "alltypes" => ("aligned", "unmapped", &["RX", "XI", "XN"]),
        "tail" => ("aligned_tail", "unmapped", &["RX"]),
        "unsorted" => ("aligned_unsorted", "unmapped", &["RX"]),
        "emptyaligned" => ("aligned_empty", "unmapped", &["RX"]),
        "bothempty" => ("aligned_empty", "unmapped_empty", &["RX"]),
        "gap" => ("aligned", "unmapped_gap", &["RX"]),
        "before" => ("aligned_before", "unmapped", &["RX"]),
        "missingtag" => ("aligned", "unmapped_missing", &["RX"]),
        "coordinate" => ("aligned_coordinate", "unmapped", &["RX"]),
        "emptyunmapped" => ("aligned", "unmapped_empty", &["RX"]),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The fixtures, written out so the port can open them. None has an index.
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

fn run(
    text: &str,
    dir: &std::path::Path,
    label: &str,
) -> Result<(Vec<u8>, Option<Vec<u8>>), TransferError> {
    let (aligned_name, unmapped_name, tags) = configuration(label);
    let aligned = ReadsDataSource::open_unindexed(&dir.join(format!("{aligned_name}.bam")))
        .expect("the aligned fixture opens with no index");
    let unmapped = ReadsDataSource::open_unindexed(&dir.join(format!("{unmapped_name}.bam")))
        .expect("the unmapped fixture opens with no index");

    let command_line = of_run(text, "commandline", label)
        .first()
        .map(|row| row.get(1).copied().unwrap_or(""))
        .unwrap_or("");
    let options = Options {
        command_line,
        ..Options::default()
    };
    let tags: Vec<String> = tags.iter().map(|tag| tag.to_string()).collect();
    tool::transfer_read_tags(&aligned, &unmapped, &tags, &options).expect("the sources read")
}

#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-transfertags-{}", std::process::id()));
    install(&text, &dir);

    let outputs = rows(&text, "output");
    let indexes = rows(&text, "index");
    assert_eq!(outputs.len(), 6, "six runs finish and six are refused");

    let mut compared = 0usize;
    for row in &outputs {
        let (label, expected_base64) = (row[0], row[1]);
        let (ours, our_index) = run(&text, &dir, label).expect("this label does not refuse");

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
        assert_eq!(
            expected_index, "absent",
            "{label}: a queryname-sorted output has nothing to index"
        );
        assert!(our_index.is_none(), "{label}: the port wrote an index");

        compared += 1;
    }

    assert_eq!(compared, 6);
    println!("transfer-read-tags: {compared} output files byte-identical");
}

/// The type is lost and the characters are kept, which only the class beside the value shows.
#[test]
fn every_tag_is_transferred_as_a_string() {
    let text = golden();
    let tags = |label: &str, name: &str| -> String {
        of_run(&text, "reads", label)
            .iter()
            .find(|row| row[1] == name)
            .map(|row| row[5].to_string())
            .unwrap_or_else(|| panic!("the golden lost {label}/{name}"))
    };

    // What went in: an Integer and a Float beside the two Strings.
    assert_eq!(
        tags("in:unmapped", "a1"),
        "RG=rg1:String;XI=42:Integer;XN=1.5:Float;RX=AAA-CCC:String"
    );
    // What came out: every one of them a String.
    assert_eq!(
        tags("alltypes", "a1"),
        "RG=rg1:String;XI=42:String;XN=1.5:String;RX=AAA-CCC:String"
    );
}

/// The finding: a read the tool cannot match is dropped, not refused.
#[test]
fn an_aligned_read_past_the_unmapped_file_is_dropped_in_silence() {
    let text = golden();
    let names = |label: &str| -> Vec<String> {
        of_run(&text, "reads", label)
            .iter()
            .map(|row| row[1].to_string())
            .collect()
    };

    // The aligned file held a1 and a9; the unmapped one ends at a6.
    assert_eq!(names("tail"), ["a1"]);
    // And nothing said so: no error row carries this label.
    assert!(
        of_run(&text, "error", "tail").is_empty(),
        "the reference raises nothing here"
    );
}

/// The traversal makes a name-only comparison; the writer makes the full one.
#[test]
fn the_writer_sorts_what_the_traversal_did_not() {
    let text = golden();
    let flags = |label: &str| -> Vec<String> {
        of_run(&text, "reads", label)
            .iter()
            .map(|row| format!("{}:{}", row[1], row[2]))
            .collect()
    };

    // In the file, second of pair before first of pair.
    assert_eq!(
        flags("in:aligned_unsorted"),
        ["a1:129", "a1:65", "a3:0", "a5:0"]
    );
    // In the output, the other way round: the queryname comparator's second tie-break.
    assert_eq!(flags("unsorted"), ["a1:65", "a1:129", "a3:0", "a5:0"]);
}

/// No read filter runs, which the fixture's read groups are chosen to show.
#[test]
fn no_read_filter_runs_between_the_two_files() {
    let text = golden();
    // The unsorted fixture's reads carry no RG at all, which WellformedReadFilter rejects.
    for row in of_run(&text, "reads", "in:aligned_unsorted") {
        assert_eq!(row[5], "-", "no tags at all on the way in");
    }
    // And all four come out, each carrying only the tag the tool transferred.
    let out = of_run(&text, "reads", "unsorted");
    assert_eq!(
        out.len(),
        4,
        "a walker would have dropped every one of them"
    );
    for row in out {
        assert_eq!(row[5], "RX=AAA-CCC:String");
    }
}

/// Six refusals, from four layers, each with the message the reference raises.
#[test]
fn the_refusals_are_the_references_and_come_from_four_layers() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-transfertags-err-{}", std::process::id()));
    install(&text, &dir);

    let message = |label: &str| -> String {
        of_run(&text, "error", label)
            .first()
            .map(|row| row[1].to_string())
            .unwrap_or_else(|| panic!("the golden lost the {label} refusal"))
    };
    let errors = rows(&text, "error");
    assert_eq!(errors.len(), 6);

    // Barclay, before the tool is built at all: the `Utils.nonEmpty` written for this is
    // unreachable, so the port keeps it and the golden shows what actually refuses.
    assert!(message("notags")
        .starts_with("org.broadinstitute.barclay.argparser.CommandLineException$MissingArgument:"));
    assert!(message("notags").contains("Argument 'read-tags' is required"));

    // Utils.validate, on the aligned side only.
    let (class, text_of) = message("coordinate")
        .split_once(':')
        .map(|(a, b)| (a.to_string(), b.to_string()))
        .expect("a class and a message");
    assert_eq!(class, "java.lang.IllegalStateException");
    assert_eq!(text_of, TransferError::AlignedNotQueryNameSorted.message());
    assert_eq!(
        run(&text, &dir, "coordinate").expect_err("refused"),
        TransferError::AlignedNotQueryNameSorted
    );

    // A UserException, and the port's own.
    assert!(
        message("emptyunmapped").ends_with(&TransferError::UnmappedEmptyAndAlignedIsNot.message())
    );
    assert_eq!(
        run(&text, &dir, "emptyunmapped").expect_err("refused"),
        TransferError::UnmappedEmptyAndAlignedIsNot
    );

    // Utils.nonNull, naming the unmapped read rather than the aligned one.
    assert!(message("missingtag").ends_with("The attribute is empty: read a3"));
    assert_eq!(
        run(&text, &dir, "missingtag").expect_err("refused"),
        TransferError::AttributeEmpty {
            unmapped: "a3".to_string()
        }
    );

    // The two IllegalStateException sites, which carry the same message from different places: one
    // inside the catch-up loop, one before it is entered.
    for (label, aligned, unmapped) in [("gap", "a3", "a4"), ("before", "a0", "a1")] {
        assert!(
            message(label).ends_with(&format!(
                "aligned read = {aligned}, unmapped read = {unmapped}"
            )),
            "{label}: {}",
            message(label)
        );
        assert_eq!(
            run(&text, &dir, label).expect_err("refused"),
            TransferError::NotInUnmapped {
                aligned: aligned.to_string(),
                unmapped: unmapped.to_string(),
            }
        );
    }
}
