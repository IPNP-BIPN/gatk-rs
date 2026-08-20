//! Conformance for `PreprocessIntervals` against GATK 4.6.2.0, compared as the whole interval list
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/PreprocessIntervalsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the overlap is resolved at the midpoint of the originals**, by an integer division that
//!    leans left: `overlap` cuts at 30 and `overlap-odd` at 31;
//!  * **the pass is sequential and in place**, which `overlap-three` shows;
//!  * **padding clamps at one and at the contig length**;
//!  * **the bins are laid from the start**, so the short bin is the last one;
//!  * **a bin length of zero means no binning at all**;
//!  * **the N filter is `allMatch` and case-insensitive**, so a bin of one non-N base survives, a
//!    bin of lower-case n does not, and a contig can come back with no bins at all;
//!  * **and no intervals at all means the whole reference**, one interval per contig.
//!
//! # What is compared
//!
//! Every byte of every list, the sequence lines included. The `M5` and `UR` fields are the
//! harness's mask: both come from the `.dict` the reference indexer wrote, the tool copies the
//! dictionary through untouched, and one of them is a path that depends on where the run happened.

use gatk_corpus as corpus;
use gatk_tools::filter_intervals::Interval;
use gatk_tools::preprocess_intervals::{preprocess, PreprocessError, Sequence};

const CONTIG_LENGTH: i32 = 240;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/preprocess_intervals.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The row of one kind and label.
fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

/// The dictionary as the interval list carries it, with the harness's mask in place.
fn sequences() -> Vec<Sequence> {
    ["chr1", "chr2"]
        .iter()
        .map(|name| Sequence {
            name: (*name).to_string(),
            length: CONTIG_LENGTH,
            md5: Some("<masked>".to_string()),
            uri: Some("<masked>".to_string()),
        })
        .collect()
}

/// The harness's reference, contig by contig.
///
/// chr1 is upper-case `ACGT` to 60, lower-case to 120, sixty Ns to 180 and upper-case again to
/// 240; chr2 is all N but for an `AC` at 150 and 151, its second line in lower case.
fn bases(contig: &str) -> Vec<u8> {
    let mut sequence = Vec::with_capacity(CONTIG_LENGTH as usize);
    for position in 1..=CONTIG_LENGTH {
        let base = match (contig, position) {
            ("chr1", 1..=60) => b"ACGT"[((position - 1) % 4) as usize],
            ("chr1", 61..=120) => b"acgt"[((position - 1) % 4) as usize],
            ("chr1", 121..=180) => b'N',
            ("chr1", _) => b"ACGT"[((position - 1) % 4) as usize],
            ("chr2", 150) => b'A',
            ("chr2", 151) => b'C',
            ("chr2", 61..=120) => b'n',
            ("chr2", _) => b'N',
            (other, _) => panic!("an unexpected contig: {other}"),
        };
        sequence.push(base);
    }
    sequence
}

fn interval(contig: &str, start: i32, end: i32) -> Interval {
    Interval {
        contig: contig.to_string(),
        start,
        end,
    }
}

/// Every run of the dump: its intervals, its bin length and its padding.
fn run(label: &str) -> (Option<Vec<Interval>>, i32, i32) {
    match label {
        "whole-genome" => (None, 50, 0),
        "whole-genome-padded" => (None, 50, 250),
        "uneven-bins" => (Some(vec![interval("chr1", 1, 100)]), 30, 0),
        "short-interval" => (Some(vec![interval("chr1", 10, 20)]), 50, 0),
        "no-bins" => (
            Some(vec![interval("chr1", 10, 20), interval("chr1", 150, 170)]),
            0,
            0,
        ),
        "clamped" => (
            Some(vec![interval("chr1", 5, 10), interval("chr1", 230, 235)]),
            0,
            20,
        ),
        "overlap" => (
            Some(vec![interval("chr1", 10, 20), interval("chr1", 41, 50)]),
            0,
            20,
        ),
        "overlap-odd" => (
            Some(vec![interval("chr1", 10, 20), interval("chr1", 42, 50)]),
            0,
            20,
        ),
        "overlap-three" => (
            Some(vec![
                interval("chr1", 10, 20),
                interval("chr1", 41, 50),
                interval("chr1", 71, 80),
            ]),
            0,
            20,
        ),
        "overlap-three-binned" => (
            Some(vec![
                interval("chr1", 10, 20),
                interval("chr1", 41, 50),
                interval("chr1", 71, 80),
            ]),
            10,
            20,
        ),
        "n-run" => (Some(vec![interval("chr1", 101, 200)]), 20, 0),
        "almost-all-n" => (Some(vec![interval("chr2", 1, CONTIG_LENGTH)]), 10, 0),
        "lower-case-n" => (Some(vec![interval("chr2", 61, 120)]), 20, 0),
        "straddling-bin" => (Some(vec![interval("chr1", 120, 140)]), 0, 0),
        "heavy-padding" => (
            Some(vec![interval("chr1", 100, 101), interval("chr1", 102, 103)]),
            0,
            60,
        ),
        "two-contigs" => (
            Some(vec![interval("chr1", 230, 235), interval("chr2", 1, 5)]),
            0,
            20,
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_list_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "whole-genome",
        "whole-genome-padded",
        "uneven-bins",
        "short-interval",
        "no-bins",
        "clamped",
        "overlap",
        "overlap-odd",
        "overlap-three",
        "overlap-three-binned",
        "n-run",
        "almost-all-n",
        "lower-case-n",
        "straddling-bin",
        "heavy-padding",
        "two-contigs",
    ] {
        let (intervals, bin_length, padding) = run(label);
        let ours = preprocess(
            intervals.as_deref(),
            &sequences(),
            bin_length,
            padding,
            bases,
        )
        .expect("a run the arguments allow");
        assert_eq!(ours, row(&text, "list", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 16, "the golden's lists");
}

#[test]
fn the_parser_bounds_are_the_two_refusals_a_port_can_reach() {
    let text = golden();
    for (label, error) in [
        (
            "negative-bin-length",
            PreprocessError::BinLengthOutOfRange(-1),
        ),
        ("negative-padding", PreprocessError::PaddingOutOfRange(-1)),
    ] {
        let bin_length = if label == "negative-bin-length" {
            -1
        } else {
            0
        };
        let padding = if label == "negative-padding" { -1 } else { 250 };
        let refused = preprocess(
            Some(&[interval("chr1", 10, 20)]),
            &sequences(),
            bin_length,
            padding,
            bases,
        )
        .expect_err("the parser's bound");
        assert_eq!(refused, error);
        assert_eq!(
            format!("{}:{}", refused.java_class(), refused.message()),
            row(&text, "error", label)
        );
    }

    // The three argument checks the tool makes before it starts. They are decisions about the
    // common interval arguments, which this port does not model as arguments at all, so they are
    // carried as messages and compared against the golden's.
    for (label, error) in [
        ("merging-rule", PreprocessError::MergingRule),
        ("interval-padding", PreprocessError::IntervalPadding),
        (
            "interval-exclusion-padding",
            PreprocessError::IntervalExclusionPadding,
        ),
    ] {
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            row(&text, "error", label)
        );
    }
}
