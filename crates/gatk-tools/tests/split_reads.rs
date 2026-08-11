//! Conformance for `SplitReads` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/SplitReadsDump.java`. Every file each run left in its
//! output directory travels in full, base64, indexes included, and so do the two input fixtures
//! with their indexes: which files exist is the measurement, so a test that named the files it
//! expected would be asserting its own guess.
//!
//! # What this suite is for
//!
//! The seventh whole tool of the record-transform archetype, and the first that opens more than
//! one writer. Five things follow from that, and each has its own test here:
//!
//!  * **every output file of a run has a different header.** `getHeaderForSAMWriter` adds a `@PG`
//!    record to the reads header in place and hands back the same object, so the nth writer's file
//!    carries n records for this one tool;
//!  * **the files come from the header, not from the reads**, as the cross product of each
//!    splitter's values over every read group, duplicates kept;
//!  * **a null value is spelled `null` from the header and `unknown` from a read**, so the file
//!    the header promises is not the file the reads go to;
//!  * **on three splitters that difference aborts the run**, because the on-demand writer accepts
//!    exactly one key;
//!  * **a read with no read group is a null pointer**, and what saves the tool is
//!    `WellformedReadFilter` rather than the tool.
//!
//! One row is measured and not ported: `missingdir`. `IOUtil.assertDirectoryIsWritable` is a
//! precondition of the reference's command line rather than of the transform, and this port does
//! not touch a filesystem. The row stays in the golden as the measurement it is.
//!
//! The command line lands in the `@PG` record's `CL`, so it is read out of the golden and handed
//! to the port rather than reconstructed: it carries the paths of the run that produced it.

use gatk_corpus as corpus;
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::with_header;
use gatk_tools::sam_output::Options;
use gatk_tools::split_reads::{self as tool, Splitter};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/split_reads.txt.gz"),
    )
}

/// Rows of one kind, split on tabs, with the kind dropped.
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

/// What each labelled run was given. A label is a configuration and the row carries nothing to
/// derive it from, so it is written here beside the dump that produced it.
struct Configuration {
    fixture: &'static str,
    splitters: &'static [Splitter],
    /// `--disable-read-filter WellformedReadFilter`, which is what lets a read with no read group
    /// reach a splitter at all.
    wellformed: bool,
    create_index: bool,
    program_record: bool,
}

fn configuration(label: &str) -> Configuration {
    let base = Configuration {
        fixture: "plain",
        splitters: &[Splitter::Sample],
        wellformed: true,
        create_index: true,
        program_record: true,
    };
    match label {
        "sample" => base,
        "readgroup" => Configuration {
            splitters: &[Splitter::ReadGroupId],
            ..base
        },
        "library" => Configuration {
            splitters: &[Splitter::LibraryName],
            ..base
        },
        "all3" => Configuration {
            splitters: &[
                Splitter::Sample,
                Splitter::ReadGroupId,
                Splitter::LibraryName,
            ],
            ..base
        },
        "none" => Configuration {
            splitters: &[],
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
        "norg-default" => Configuration {
            fixture: "norg",
            ..base
        },
        "norg-nofilter" => Configuration {
            fixture: "norg",
            wellformed: false,
            ..base
        },
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The fixtures, written out so the port's reader can open them.
fn install_fixtures(text: &str, dir: &std::path::Path) {
    std::fs::create_dir_all(dir).expect("a scratch directory");
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

fn run(
    text: &str,
    dir: &std::path::Path,
    label: &str,
) -> Result<Vec<tool::OutputFile>, tool::SplitError> {
    let config = configuration(label);
    let source = ReadsDataSource::open(
        &dir.join(format!("{}.bam", config.fixture)),
        &dir.join(format!("{}.bai", config.fixture)),
    )
    .expect("the fixture opens");
    let header = source.header().clone();

    let filter: Box<dyn Fn(&BamRecord) -> bool> = if config.wellformed {
        Box::new(move |read: &BamRecord| with_header::wellformed(read, &header))
    } else {
        Box::new(|_: &BamRecord| true)
    };

    let command_line = of_run(text, "commandline", label)
        .first()
        .map(|row| row.get(1).copied().unwrap_or(""))
        .unwrap_or("");
    let options = Options {
        create_output_bam_index: config.create_index,
        add_output_sam_program_record: config.program_record,
        command_line,
        ..Options::default()
    };

    tool::split_reads(
        &source,
        &options,
        config.splitters,
        config.fixture,
        ".bam",
        filter.as_ref(),
    )
    .expect("the source reads")
}

/// Every file of every run that finishes, by name, byte for byte.
#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-splitreads-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let labels: Vec<&str> = rows(&text, "outcount")
        .into_iter()
        .map(|row| row[0])
        .collect();
    assert_eq!(labels.len(), 7, "seven runs finish and three do not");

    let mut compared = 0usize;
    for label in labels {
        let ours = run(&text, &dir, label).expect("this label finishes");
        let expected = of_run(&text, "outfile", label);
        let indexes = of_run(&text, "outindex", label);
        let count: usize = of_run(&text, "outcount", label)[0][1]
            .parse()
            .expect("a file count");
        assert_eq!(ours.len(), count, "{label}: number of files");
        assert_eq!(expected.len(), count);

        // The reference lists its directory, so the golden is by name; the port returns the files
        // in the order their writers were created.
        let mut by_name: Vec<&tool::OutputFile> = ours.iter().collect();
        by_name.sort_by(|a, b| a.name.cmp(&b.name));

        for (ours, row) in by_name.iter().zip(&expected) {
            let name = row[1];
            assert_eq!(ours.name, name, "{label}: file name");
            let bytes = corpus::decode_base64(row[2]);
            assert_eq!(ours.bam.len(), bytes.len(), "{label}/{name}: output length");
            if ours.bam != bytes {
                let at = ours
                    .bam
                    .iter()
                    .zip(&bytes)
                    .position(|(a, b)| a != b)
                    .unwrap_or(0);
                panic!("{label}/{name}: first byte difference at offset {at}");
            }
            let expected_index = indexes
                .iter()
                .find(|index| index[1] == name)
                .map(|index| index[2])
                .expect("an index row for every file");
            match (&ours.index, expected_index) {
                (None, "absent") => {}
                (Some(_), "absent") => panic!("{label}/{name}: the reference wrote no index"),
                (None, _) => panic!("{label}/{name}: the reference wrote an index"),
                (Some(ours), expected) => assert_eq!(
                    *ours,
                    corpus::decode_base64(expected),
                    "{label}/{name}: the .bai"
                ),
            }
            compared += 1;
        }
    }
    assert_eq!(compared, 16, "sixteen files across seven runs");
    println!("split-reads: {compared} output files byte-identical");
}

/// The nth writer of a run carries n `@PG` records for this one tool.
///
/// A port that built one header and handed it to every writer would produce files that are each
/// plausible and none of them the reference's.
#[test]
fn each_file_carries_one_more_program_record_than_the_last() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-splitreads-pg-{}", std::process::id()));
    install_fixtures(&text, &dir);

    // Splitting by library is the run that shows it: four writers, and the fourth is made during
    // the traversal rather than before it.
    let expected: Vec<Vec<&str>> = of_run(&text, "programs", "library");
    let ids: Vec<Vec<&str>> = expected
        .iter()
        .map(|row| row[2].split(';').collect())
        .collect();
    assert_eq!(
        ids,
        [
            vec!["upstream", "GATK SplitReads"],
            vec!["upstream", "GATK SplitReads", "GATK SplitReads.1"],
            vec![
                "upstream",
                "GATK SplitReads",
                "GATK SplitReads.1",
                "GATK SplitReads.2"
            ],
            vec![
                "upstream",
                "GATK SplitReads",
                "GATK SplitReads.1",
                "GATK SplitReads.2",
                "GATK SplitReads.3"
            ],
        ],
        "the golden lost the header growth this suite is for"
    );

    let ours = run(&text, &dir, "library").expect("it finishes");
    let mut by_name: Vec<&tool::OutputFile> = ours.iter().collect();
    by_name.sort_by(|a, b| a.name.cmp(&b.name));
    for (file, row) in by_name.iter().zip(&expected) {
        // Read the header back off the bytes rather than trusting the builder.
        let decompressed =
            htsjdk_bgzf::read::decompress_all(&file.bam).expect("the port's own output is BGZF");
        let header = htsjdk_bam::reader::BamReader::new(&decompressed)
            .expect("the port's own output opens")
            .header
            .text
            .clone();
        let ids: Vec<&str> = header.programs.iter().map(|p| p.id.as_str()).collect();
        assert_eq!(ids.join(";"), row[2], "the @PG list of {}", file.name);
    }
}

/// The file the header promises is not the file the reads go to.
#[test]
fn a_null_library_makes_two_files_one_of_them_empty() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-splitreads-null-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let names: Vec<&str> = of_run(&text, "outfile", "library")
        .iter()
        .map(|row| row[1])
        .collect();
    assert_eq!(
        names,
        [
            "plain.lib1.bam",
            "plain.lib2.bam",
            "plain.null.bam",
            "plain.unknown.bam"
        ]
    );

    let ours = run(&text, &dir, "library").expect("it finishes");
    let empty = ours
        .iter()
        .find(|file| file.name == "plain.null.bam")
        .expect("the file the header promised");
    let reads = ours
        .iter()
        .find(|file| file.name == "plain.unknown.bam")
        .expect("the file the reads went to");
    // The empty one is a valid BAM with an index, not a missing file.
    assert!(empty.index.is_some());
    assert!(empty.bam.len() < reads.bam.len());
}

/// The three runs that do not finish, which are findings rather than edge cases.
#[test]
fn three_runs_abort_and_two_of_them_are_the_tool_s_own() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-splitreads-err-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let errors = rows(&text, "error");
    assert_eq!(errors.len(), 3);

    for row in &errors {
        let (label, expected) = (row[0], row[1]);
        let (class, message) = expected.split_once(':').expect("a class and a message");
        // The output directory has to exist before the tool runs, which is a precondition of the
        // reference's command line rather than of the transform: this port has no filesystem to
        // check. The row is measured and left as measurement.
        if label == "missingdir" {
            assert_eq!(class, "htsjdk.samtools.SAMException");
            continue;
        }
        let error = run(&text, &dir, label).expect_err("this label aborts");
        assert_eq!(class, error.class(), "the class for {label}");
        assert_eq!(message, error.message(), "the message for {label}");
    }

    // And the same fixture with the default filters finishes, which is what makes the null
    // pointer a property of the filter chain rather than of the tool.
    assert!(run(&text, &dir, "norg-default").is_ok());
}
