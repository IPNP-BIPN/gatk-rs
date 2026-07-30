//! Conformance for `IntervalUtils.parseIntervalArguments` against GATK 4.6.2.0.
//!
//! One row per `-L` argument: what it resolved to, or which exception refused it.
//!
//! The golden corrected an assumption the port was written on. `FeatureManager.isFeatureFile`
//! asks every codec whether it `canDecode` the path, and the codecs answer by **extension** rather
//! than by content, so a `.list` file holding a BED body is not a Feature file: it falls through
//! to the interval-file reader and dies parsing `chr1\t0\t10` as a genome location. The row named
//! `bed-contents-list-extension` is that measurement.
//!
//! Two cases are not compared yet and are named rather than skipped quietly: `.interval_list` and
//! `.bed` are Feature files, and reading them needs the codecs of G1.3. The test asserts that the
//! port refuses them *today* for the documented reason, so this assertion fails the moment
//! `FeatureDataSource` lands and has to be removed with it.

use gatk_corpus as corpus;
use gatk_engine::interval_args::{self, IntervalArgumentError};
use htsjdk_bam::header::{SamHeader, SequenceRecord};

const CONTIG_LENGTH: i32 = 200;

/// The cases whose answer needs Feature codecs. Both are files whose names look like intervals
/// and whose contents are decoded by a codec instead.
/// `.interval_list` is still pending: it is a Feature file to the reference, and htsjdk's
/// IntervalList reader is its own slice. `bed` left this list when the BED codec landed.
const PENDING_FEATURE_SOURCES: [&str; 1] = ["picard-interval-list"];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/interval_file.txt.gz"),
    )
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    for name in ["chr1", "chr2"] {
        header
            .sequences
            .push(SequenceRecord::new(name, CONTIG_LENGTH));
    }
    header
}

/// The files the harness wrote, by the label of the case that used them. The contents are the
/// input, so they are here rather than in the golden.
fn fixture(label: &str) -> Option<(&'static str, &'static str)> {
    // (file name, contents)
    match label {
        "list" => Some(("a.list", "chr1:1-10\nchr1:50-60\nchr2\n")),
        "intervals" => Some(("b.intervals", "chr1:1-10\nchr1:50-60\nchr2\n")),
        "whitespace" => Some(("ws.list", "\n  chr1:1-10  \n\n\tchr2:5-6\n\n")),
        "blank-only" => Some(("blank.list", "\n\n   \n")),
        "empty" => Some(("empty.list", "")),
        "uppercase-extension" => Some(("c.LIST", "chr1:1-5\n")),
        "unknown-extension" => Some(("d.txt", "chr1:1-10\n")),
        "bed" => Some(("f.bed", "chr1\t0\t10\nchr2\t4\t6\n")),
        "bed-contents-list-extension" => Some(("g.list", "chr1\t0\t10\nchr2\t4\t6\n")),
        // Written by htsjdk's IntervalList writer, header included.
        "picard-interval-list" => Some((
            "e.interval_list",
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:200\n@SQ\tSN:chr2\tLN:200\n\
             chr1\t1\t10\t+\t.\nchr2\t5\t6\t+\t.\n",
        )),
        _ => None,
    }
}

/// The argument each case passed, given the directory the fixtures live in.
fn argument(label: &str, dir: &std::path::Path) -> String {
    match label {
        "missing-list" => dir.join("absent.list").display().to_string(),
        "missing-unknown-extension" => dir.join("absent.txt").display().to_string(),
        "semicolon" => "chr1:1-10;chr2:1-10".to_string(),
        "literal" => "chr1:1-10".to_string(),
        "literal-whole-contig" => "chr2".to_string(),
        other => {
            let (name, _) = fixture(other).unwrap_or_else(|| panic!("{other} has no fixture"));
            dir.join(name).display().to_string()
        }
    }
}

/// The exception class the reference raised, for the refusal the port produced.
fn class_of(error: &IntervalArgumentError) -> &'static str {
    match error {
        IntervalArgumentError::IntervalFileEmpty => {
            "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile"
        }
        IntervalArgumentError::IntervalFileMissing(_)
        | IntervalArgumentError::FileIsNeitherFeaturesNorIntervals(_) => {
            "org.broadinstitute.hellbender.exceptions.UserException$CouldNotReadInputFile"
        }
        IntervalArgumentError::LegacySemicolonSyntax(_) => {
            "org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue"
        }
        // Every parse failure surfaces as one class: an unknown contig and malformed positions
        // are different messages of the same exception.
        IntervalArgumentError::Parse(_) => {
            "org.broadinstitute.hellbender.exceptions.UserException$MalformedGenomeLoc"
        }
        other => panic!("{other:?} has no reference class"),
    }
}

#[test]
fn every_argument_resolves_the_way_the_reference_resolves_it() {
    let text = golden();
    let header = header();

    let dir = std::env::temp_dir().join(format!("gatk-rs-intervalfile-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for label in [
        "list",
        "intervals",
        "whitespace",
        "blank-only",
        "empty",
        "uppercase-extension",
        "unknown-extension",
        "bed",
        "bed-contents-list-extension",
        "picard-interval-list",
    ] {
        let (name, contents) = fixture(label).expect("a fixture");
        std::fs::write(dir.join(name), contents).unwrap();
    }

    let mut compared = 0;
    let mut pending = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("case\t") else {
            continue;
        };
        let mut parts = rest.split('\t');
        let label = parts.next().expect("a label");
        let outcome = parts.next().expect("an outcome");
        let count: usize = parts.next().expect("a count").parse().expect("a number");
        let expected = parts.next().unwrap_or("");

        let result = interval_args::parse_interval_arguments(
            &argument(label, &dir),
            &header,
            &gatk_engine::feature_intervals::BedFeatureSource,
        );

        if PENDING_FEATURE_SOURCES.contains(&label) {
            // Named, not skipped: the reference reads these through a codec, and the port has no
            // codecs yet. When G1.3 lands, this assertion fails and goes away with the seam.
            assert!(
                matches!(
                    result,
                    Err(IntervalArgumentError::FileIsNeitherFeaturesNorIntervals(_))
                ),
                "{label}: expected the Feature seam to be empty, got {result:?}"
            );
            assert_eq!(outcome, "ok", "{label} is only pending because it succeeds");
            pending += 1;
            continue;
        }

        match (result, outcome) {
            (Ok(intervals), "ok") => {
                let ours: Vec<String> = intervals
                    .iter()
                    .map(|i| format!("{}:{}-{}", i.contig, i.start, i.end))
                    .collect();
                assert_eq!(ours.len(), count, "{label}: interval count");
                assert_eq!(ours.join("|"), expected, "{label}");
            }
            (Err(error), outcome) if outcome.starts_with("E:") => {
                assert_eq!(
                    format!("E:{}", class_of(&error)),
                    outcome,
                    "{label}: the wrong refusal"
                );
            }
            (Ok(_), outcome) => panic!("{label}: the reference raised {outcome}, the port did not"),
            (Err(error), _) => {
                panic!("{label}: the port raised {error:?}, the reference did not")
            }
        }
        compared += 1;
    }

    std::fs::remove_dir_all(&dir).ok();
    assert!(compared > 0, "the golden carries no case rows");
    println!(
        "{compared} -L arguments resolved identically, {pending} pending Feature sources (G1.3)"
    );
}
