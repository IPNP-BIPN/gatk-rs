//! Conformance for `UnmarkDuplicates` and `RevertBaseQualityScores` against GATK 4.6.2.0,
//! compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/RecordTransformDump.java`. The output BAMs and their
//! indexes travel in full, base64, as `PrintReads`'s do, and so do the three input fixtures and
//! their indexes: a test that built its own index would be inventing part of the input, the same
//! mistake as reconstructing a command line.
//!
//! # What this suite is for
//!
//! G2's calibration gate. `PrintReads` was the first whole tool; these are the second and third,
//! and the question is what a member of the largest archetype costs once the engine is paid for.
//! The answer is not in the transform — both `apply` bodies are two lines — it is in three things
//! the archetype hides, and each has its own case here:
//!
//!  * **the default read filters are not `PrintReads`'s.** Both tools override
//!    `getDefaultReadFilters` to `ALLOW_ALL_READS`, where `PrintReads` takes `GATKTool`'s default
//!    of `WellformedReadFilter`. The `wellformed` case asks for that filter explicitly, so the
//!    difference between "the tool's default" and "the engine's default" is visible rather than
//!    argued;
//!  * **`RevertBaseQualityScores` aborts the run** on a read with no `OQ`. The golden carries an
//!    `error` row rather than an output, and this file asserts the class and the message;
//!  * **an empty `OQ` is the same as an absent one**, which is the reference's own conflation
//!    inside `getOriginalBaseQualities` and reaches the same exception.
//!
//! The command line lands in the `@PG` record's `CL`, so it is read out of the golden and handed
//! to the port rather than reconstructed: it carries the paths of the run that produced it.

use gatk_corpus as corpus;
use gatk_engine::interval::{self, SimpleInterval};
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::with_header;
use gatk_tools::revert_base_quality_scores::{self, RevertError};
use gatk_tools::sam_output::Options;
use gatk_tools::unmark_duplicates;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/record_transform.txt.gz"),
    )
}

/// Rows of one kind that carry `<label>\t<value>`.
fn pairs<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str)> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .filter_map(|rest| rest.split_once('\t'))
        .collect()
}

/// Rows of one kind that carry `<tool>\t<label>\t<value>`.
fn triples<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str, &'a str)> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .filter_map(|rest| {
            let (tool, rest) = rest.split_once('\t')?;
            let (label, value) = rest.split_once('\t')?;
            Some((tool, label, value))
        })
        .collect()
}

fn lookup<'a>(rows: &[(&'a str, &'a str, &'a str)], tool: &str, label: &str) -> Option<&'a str> {
    rows.iter()
        .find(|(t, l, _)| *t == tool && *l == label)
        .map(|(_, _, value)| *value)
}

/// What each labelled run was given. A label is a configuration and the row carries nothing to
/// derive it from, so it is written here beside the dump that produced it.
struct Configuration {
    fixture: &'static str,
    intervals: &'static [&'static str],
    /// `--read-filter WellformedReadFilter`, which *adds* to the tool's own default rather than
    /// replacing it. The tool's default keeps everything, so adding this one is the whole filter.
    wellformed: bool,
    create_index: bool,
    program_record: bool,
}

fn configuration(tool: &str, label: &str) -> Configuration {
    let base = Configuration {
        fixture: "full",
        intervals: &[],
        wellformed: false,
        create_index: true,
        program_record: true,
    };
    match (tool, label) {
        ("UnmarkDuplicates", "all") | ("RevertBaseQualityScores", "all") => base,
        (_, "chr1") => Configuration {
            intervals: &["chr1"],
            ..base
        },
        (_, "chr1:100-160") => Configuration {
            intervals: &["chr1:100-160"],
            ..base
        },
        (_, "wellformed") => Configuration {
            wellformed: true,
            ..base
        },
        (_, "noindex") => Configuration {
            create_index: false,
            ..base
        },
        (_, "nopg") => Configuration {
            program_record: false,
            ..base
        },
        (_, "partial-input") | (_, "missing-oq") => Configuration {
            fixture: "partial",
            ..base
        },
        (_, "empty-oq") => Configuration {
            fixture: "empty-oq",
            ..base
        },
        (tool, other) => panic!("{tool}/{other} is in the golden but not configured here"),
    }
}

/// The three fixtures, written out so the port's reader can open them.
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
    let dir = std::env::temp_dir().join(format!("gatk-rs-recordtransform-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let outputs = triples(&text, "output");
    let indexes = triples(&text, "index");
    let command_lines = triples(&text, "commandline");
    let errors = triples(&text, "error");
    assert_eq!(outputs.len(), command_lines.len());
    assert!(
        errors.len() >= 2,
        "the golden lost the runs that abort, which are two of the three findings this suite is \
         for"
    );

    let mut compared = 0usize;
    for (tool, label, expected_base64) in &outputs {
        let config = configuration(tool, label);
        let bam = dir.join(format!("{}.bam", config.fixture));
        let bai = dir.join(format!("{}.bai", config.fixture));
        let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
        let header = source.header().clone();

        let intervals: Vec<SimpleInterval> = config
            .intervals
            .iter()
            .map(|text| interval::parse_interval(text, &header).expect("a parsable interval"))
            .collect();

        // Both tools default to ALLOW_ALL_READS, which is why the base case keeps everything.
        let header_for_filter = header.clone();
        let filter: Box<dyn Fn(&BamRecord) -> bool> = if config.wellformed {
            Box::new(move |read: &BamRecord| with_header::wellformed(read, &header_for_filter))
        } else {
            Box::new(|_: &BamRecord| true)
        };

        let command_line =
            lookup(&command_lines, tool, label).expect("a command line for every output");
        let options = Options {
            intervals,
            create_output_bam_index: config.create_index,
            add_output_sam_program_record: config.program_record,
            command_line,
            ..Options::default()
        };

        let (ours, our_index) = match *tool {
            "UnmarkDuplicates" => {
                unmark_duplicates::unmark_duplicates(&source, &options, filter.as_ref())
                    .expect("the tool runs")
            }
            "RevertBaseQualityScores" => revert_base_quality_scores::revert_base_quality_scores(
                &source,
                &options,
                filter.as_ref(),
            )
            .expect("the source reads")
            .expect("this label does not abort"),
            other => panic!("{other} is in the golden but not wired here"),
        };

        let expected = corpus::decode_base64(expected_base64);
        assert_eq!(
            ours.len(),
            expected.len(),
            "{tool}/{label}: output length differs"
        );
        if ours != expected {
            let at = ours
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            panic!("{tool}/{label}: first byte difference at offset {at}");
        }

        let expected_index = lookup(&indexes, tool, label).expect("an index row for every output");
        match (our_index, expected_index) {
            (None, "absent") => {}
            (Some(_), "absent") => {
                panic!("{tool}/{label}: the reference wrote no index and the port did")
            }
            (None, _) => {
                panic!("{tool}/{label}: the reference wrote an index and the port did not")
            }
            (Some(ours), expected) => {
                assert_eq!(
                    ours,
                    corpus::decode_base64(expected),
                    "{tool}/{label}: the .bai"
                );
            }
        }
        compared += 1;
    }

    assert!(compared >= 10, "only {compared} outputs compared");
    println!(
        "record-transform: {compared} output files byte-identical, {} refusals reproduced",
        errors.len()
    );
}

/// The runs that abort, which are the finding rather than an edge case.
///
/// A port that skipped the offending read, or passed it through, would produce a larger and
/// healthier-looking file than the reference. This asserts that it produces no file at all, and
/// that the class and message are the reference's.
#[test]
fn a_read_without_original_qualities_aborts_the_run() {
    let text = golden();
    let dir =
        std::env::temp_dir().join(format!("gatk-rs-recordtransform-e-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let errors = triples(&text, "error");
    assert!(!errors.is_empty(), "the golden carries no error rows");

    let mut checked = 0usize;
    for (tool, label, expected) in &errors {
        assert_eq!(
            *tool, "RevertBaseQualityScores",
            "only the revert tool refuses a read"
        );
        let config = configuration(tool, label);
        let source = ReadsDataSource::open(
            &dir.join(format!("{}.bam", config.fixture)),
            &dir.join(format!("{}.bai", config.fixture)),
        )
        .expect("the fixture opens");

        let options = Options {
            command_line: "",
            ..Options::default()
        };
        let outcome = revert_base_quality_scores::revert_base_quality_scores(
            &source,
            &options,
            &|_: &BamRecord| true,
        )
        .expect("the source reads");
        let error = outcome.expect_err("the run must abort");

        let (class, message) = expected
            .split_once(':')
            .expect("the golden's error row is class:message");
        assert_eq!(error.class(), class, "{label}: exception class");
        assert_eq!(error.message(), message, "{label}: exception message");
        assert!(
            matches!(error, RevertError::NoOriginalQualities { .. }),
            "{label}: an absent and an empty OQ reach the same refusal"
        );
        checked += 1;
    }
    assert_eq!(
        checked, 2,
        "both ways to have no original qualities must be covered: absent, and present but empty"
    );
}
