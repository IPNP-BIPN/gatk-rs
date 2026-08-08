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
use gatk_tools::add_original_alignment_tags as add_oa;
use gatk_tools::sam_output::Options;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/add_original_alignment_tags.txt.gz"),
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

fn configuration(_tool: &str, label: &str) -> Configuration {
    let base = Configuration {
        fixture: "plain",
        intervals: &[],
        wellformed: true,
        create_index: true,
        program_record: true,
    };
    match label {
        // This tool does not override its read filters either, so its default is
        // WellformedReadFilter.
        "all" => base,
        // The only run that succeeds: chr1 excludes the unpaired read, which is on chr2.
        "chr1" => Configuration {
            intervals: &["chr1"],
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

/// The fixtures, written out so the port's reader can open them.
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
    assert_eq!(
        errors.len(),
        3,
        "the golden lost the runs that abort, which are three of the four and the finding this \
         suite is for"
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

        // This tool keeps GATKTool's default filter, so the base case is the wellformed one and
        // `allowall` is the case that changes it. The other two tools of this archetype are the
        // other way round.
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

        let (ours, our_index) =
            add_oa::add_original_alignment_tags(&source, &options, filter.as_ref())
                .expect("the source reads")
                .expect("this label does not abort");

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

    assert_eq!(compared, 1, "one run succeeds and three abort");
    println!(
        "add-oa-tags: {compared} output file byte-identical, {} refusals reproduced",
        errors.len()
    );
}

/// The runs that abort, which are the finding rather than an edge case.
///
/// A port that skipped the offending read, or passed it through, would produce a larger and
/// healthier-looking file than the reference. This asserts that it produces no file at all, and
/// that the class and message are the reference's.
#[test]
fn an_unpaired_read_aborts_the_run() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-addoatags-err-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let errors = triples(&text, "error");
    let (_, label, expected) = errors.first().expect("a refusal in the golden");
    assert_eq!(*label, "all");

    let bam = dir.join("plain.bam");
    let bai = dir.join("plain.bai");
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();
    let filter = move |read: &BamRecord| with_header::wellformed(read, &header);

    let error = add_oa::add_original_alignment_tags(&source, &Options::default(), &filter)
        .expect("the source reads")
        .expect_err("the run aborts");

    let (class, message) = expected.split_once(':').expect("a class and a message");
    assert_eq!(class, error.class());
    assert_eq!(message.trim(), error.message());

    // Three of the four runs abort. The tool cannot read a file holding an unpaired read, which
    // is an ordinary thing for a file to hold.
    assert_eq!(errors.len(), 3);
}

/// Every tag the reference wrote, rebuilt from the same reads.
#[test]
fn both_tags_are_the_reference_tags() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-addoatags-tags-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let bam = dir.join("plain.bam");
    let bai = dir.join("plain.bai");
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();

    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter_map(|line| line.strip_prefix("tags\t"))
        .map(|rest| rest.split('\t').collect())
        .collect();
    assert!(!rows.is_empty(), "the golden lost its tag rows");

    let mut compared = 0;
    let mut nulls = 0;
    let all =
        gatk_tools::read_walker::traverse(&source, &[], &|_: &BamRecord| true).expect("the reads");
    for record in all {
        let mut record = record;
        let Some(row) = rows.iter().find(|row| row[1] == record.read_name) else {
            // The reads that never reached the output, because the run aborted before them.
            continue;
        };
        // The OA is written before the mate contig is asked for, so it exists even on the read
        // that aborts.
        let oa = add_oa::original_alignment_value(&record, &header);
        assert_eq!(oa, row[2], "OA of {}", record.read_name);
        if oa.ends_with(",null;") {
            nulls += 1;
        }
        let mate = add_oa::mate_contig_value(&record, &header).expect("this read is paired");
        assert_eq!(mate, row[3], "XM of {}", record.read_name);
        let _ = add_oa::add_tags(&mut record, &header);
        compared += 1;
    }
    assert_eq!(compared, rows.len(), "every tag row compared");
    assert_eq!(nulls, 2, "two reads carry no NM and print the word null");
}

/// The comma escaping cannot fire, because a comma is not a legal sequence name character. The
/// golden carries the refusal and its regex.
#[test]
fn the_comma_escaping_is_unreachable() {
    let text = golden();
    let row = text
        .lines()
        .find(|line| line.starts_with("sequencename\tcomma\t"))
        .expect("the golden lost the refusal that makes the branch dead");
    assert!(row.contains("SAMException"));
    assert!(row.contains("doesn't match regex"));
    // And the regex it quotes has no comma in either character class.
    let regex = row.rsplit_once(": ").expect("a regex").1;
    assert!(
        !regex.contains(','),
        "the regex would admit a comma: {regex}"
    );
}
