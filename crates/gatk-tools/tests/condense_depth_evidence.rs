//! Conformance for `CondenseDepthEvidence` against GATK 4.6.2.0, compared as the whole output file
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/CondenseDepthEvidenceDump.java`.
//!
//! # What this suite is for
//!
//!  * **the merged intervals and their summed counts**;
//!  * **the maximum firing one bin late**, so a limit of 150 and a limit of 200 give the same file
//!    over hundred-base bins;
//!  * **the minimum dropping records outright**, the trailing accumulator included;
//!  * **adjacency being `end + 1 == start` on one contig**;
//!  * **and the two refusals**, a minimum above the maximum and an output named for another
//!    feature type.

use gatk_corpus as corpus;
use gatk_tools::condense_depth_evidence::{
    check_lengths, check_output, condense, read, write, Arguments, CondenseError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/condense_depth_evidence.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

fn run(text: &str, label: &str, arguments: &Arguments) -> String {
    let (samples, records) = read(&value(text, "input", label));
    write(&samples, &condense(&records, arguments))
}

#[test]
fn every_condensed_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, arguments) in [
        ("defaults", Arguments::default()),
        (
            "max-200",
            Arguments {
                max_interval_length: 200,
                ..Arguments::default()
            },
        ),
        (
            "max-150",
            Arguments {
                max_interval_length: 150,
                ..Arguments::default()
            },
        ),
        (
            "min-300",
            Arguments {
                min_interval_length: 300,
                ..Arguments::default()
            },
        ),
        (
            "min-200-max-300",
            Arguments {
                min_interval_length: 200,
                max_interval_length: 300,
            },
        ),
        ("single", Arguments::default()),
        (
            "single-dropped",
            Arguments {
                min_interval_length: 200,
                ..Arguments::default()
            },
        ),
    ] {
        assert_eq!(
            run(&text, label, &arguments),
            value(&text, "condensed", label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 7, "the golden's outputs");
}

/// The limit is tested against what is already held, so it fires one bin late and the interval
/// written is longer than the limit.
#[test]
fn the_maximum_is_exceeded_by_one_bin() {
    let text = golden();
    let hundred_fifty = run(
        &text,
        "max-150",
        &Arguments {
            max_interval_length: 150,
            ..Arguments::default()
        },
    );
    let two_hundred = run(
        &text,
        "max-200",
        &Arguments {
            max_interval_length: 200,
            ..Arguments::default()
        },
    );
    assert_eq!(hundred_fifty, two_hundred);
    // And what both write is intervals of two hundred bases, over a limit of one hundred and
    // fifty.
    let lengths: Vec<i32> = hundred_fifty
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            fields[2].parse::<i32>().expect("an end") - fields[1].parse::<i32>().expect("a start")
        })
        .collect();
    assert_eq!(lengths, vec![200, 200, 200, 200, 200, 200, 200]);
}

/// A run shorter than the minimum is not written at all, the last accumulator included.
#[test]
fn the_minimum_drops_records_rather_than_merging_them() {
    let text = golden();
    let kept = run(
        &text,
        "min-300",
        &Arguments {
            min_interval_length: 300,
            ..Arguments::default()
        },
    );
    assert_eq!(
        kept.lines().filter(|line| !line.starts_with('#')).count(),
        1
    );
    // The two hundred-base runs at the end of the file are gone, and nothing says so.
    assert!(kept.contains("chr1\t0\t1000\t55\t955\n"));

    // A single bin shorter than the minimum leaves a file with a header and nothing else.
    let dropped = run(
        &text,
        "single-dropped",
        &Arguments {
            min_interval_length: 200,
            ..Arguments::default()
        },
    );
    assert_eq!(dropped, "#Chr\tStart\tEnd\tsampleA\tsampleB\n");
}

#[test]
fn the_two_refusals_match_the_golden() {
    let text = golden();
    let error = check_lengths(&Arguments {
        min_interval_length: 500,
        max_interval_length: 100,
    })
    .expect_err("the length refusal");
    assert_eq!(error, CondenseError::MinimumAboveMaximum);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "min-above-max")
    );

    let error = check_output("<dir>/condensed.baf.txt").expect_err("the output refusal");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "wrong-output-type")
    );
    assert_eq!(check_output("<dir>/condensed.rd.txt"), Ok(()));
}
