//! Conformance for `CallCopyRatioSegments` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/CallCopyRatioSegmentsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the bounds are inclusive and in copy-ratio space**, which the two segments placed at
//!    exactly 0.9 and exactly 1.1 are there to show;
//!  * **a single copy-neutral segment divides by zero** and calls everything neutral rather than
//!    failing;
//!  * **an empty copy-neutral set does the same** from the other side;
//!  * **the weighting is by interval length**, so a 900-base segment outweighs two 10-base ones;
//!  * **and the legacy file reorders the columns**, putting the call before the mean.

use gatk_corpus as corpus;
use gatk_tools::call_copy_ratio_segments::{
    self, Call, CallerError, CopyRatioSegment, DEFAULT_CALLING_Z_SCORE,
    DEFAULT_NEUTRAL_LOWER_BOUND, DEFAULT_NEUTRAL_UPPER_BOUND, DEFAULT_OUTLIER_Z_SCORE,
};

/// The SAM header every input and output carries.
const HEADER: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:GATKCopyNumber\tSM:SAMPLE\n";

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/call_copy_ratio_segments.txt.gz"),
    )
}

fn file(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .replace("\\t", "\t")
        .replace("\\n", "\n")
}

fn segment(start: i32, end: i32, mean_log2: f64) -> CopyRatioSegment {
    CopyRatioSegment {
        contig: "chr1".to_string(),
        start,
        end,
        num_points: 10,
        mean_log2_copy_ratio: mean_log2,
    }
}

/// The dump's six-segment spread, with two sitting exactly on the bounds.
fn spread() -> Vec<CopyRatioSegment> {
    vec![
        segment(1, 100, 0.0),
        segment(101, 200, -0.15200309344504997),
        segment(201, 300, 0.13750352374993502),
        segment(301, 400, 0.5),
        segment(401, 500, 1.5),
        segment(501, 600, -1.5),
    ]
}

fn called(segments: &[CopyRatioSegment], lower: f64, upper: f64) -> (String, String) {
    let calls = call_copy_ratio_segments::make_calls(
        segments,
        lower,
        upper,
        DEFAULT_OUTLIER_Z_SCORE,
        DEFAULT_CALLING_Z_SCORE,
    )
    .expect("the arguments are valid");
    (
        call_copy_ratio_segments::write_called(HEADER, segments, &calls),
        call_copy_ratio_segments::write_legacy("SAMPLE", segments, &calls),
    )
}

#[test]
fn every_called_table_matches_the_golden() {
    let text = golden();

    for (label, segments, lower, upper) in [
        (
            "spread",
            spread(),
            DEFAULT_NEUTRAL_LOWER_BOUND,
            DEFAULT_NEUTRAL_UPPER_BOUND,
        ),
        (
            "all-neutral",
            vec![
                segment(1, 100, 0.0),
                segment(101, 200, 0.05),
                segment(201, 300, -0.05),
            ],
            DEFAULT_NEUTRAL_LOWER_BOUND,
            DEFAULT_NEUTRAL_UPPER_BOUND,
        ),
        (
            "one-neutral",
            vec![
                segment(1, 100, 0.0),
                segment(101, 200, 1.5),
                segment(201, 300, -1.5),
            ],
            DEFAULT_NEUTRAL_LOWER_BOUND,
            DEFAULT_NEUTRAL_UPPER_BOUND,
        ),
        (
            "none-neutral",
            vec![segment(1, 100, 1.5), segment(101, 200, -1.5)],
            DEFAULT_NEUTRAL_LOWER_BOUND,
            DEFAULT_NEUTRAL_UPPER_BOUND,
        ),
        (
            "weighted",
            vec![
                segment(1, 10, 0.0),
                segment(11, 910, 0.09),
                segment(911, 920, -0.09),
                segment(921, 1000, 0.6),
            ],
            DEFAULT_NEUTRAL_LOWER_BOUND,
            DEFAULT_NEUTRAL_UPPER_BOUND,
        ),
        ("wide-bounds", spread(), 0.5, 2.0),
    ] {
        let (table, legacy) = called(&segments, lower, upper);
        assert_eq!(table, file(&text, "called", label), "{label}: the table");
        assert_eq!(
            legacy,
            file(&text, "legacy", label),
            "{label}: the legacy file"
        );
    }
}

/// A single copy-neutral segment divides by zero and calls everything neutral.
#[test]
fn one_copy_neutral_segment_calls_everything_neutral() {
    let segments = vec![
        segment(1, 100, 0.0),
        segment(101, 200, 1.5),
        segment(201, 300, -1.5),
    ];
    let statistics = call_copy_ratio_segments::calling_statistics(
        &segments,
        DEFAULT_NEUTRAL_LOWER_BOUND,
        DEFAULT_NEUTRAL_UPPER_BOUND,
        DEFAULT_OUTLIER_Z_SCORE,
    );
    assert!(
        !statistics.standard_deviation.is_finite(),
        "the denominator is zero, {statistics:?}"
    );
    let calls = call_copy_ratio_segments::make_calls(
        &segments,
        DEFAULT_NEUTRAL_LOWER_BOUND,
        DEFAULT_NEUTRAL_UPPER_BOUND,
        DEFAULT_OUTLIER_Z_SCORE,
        DEFAULT_CALLING_Z_SCORE,
    )
    .expect("valid arguments");
    assert!(calls.iter().all(|call| *call == Call::Neutral));
}

/// An empty copy-neutral set gives a NaN mean, and every comparison against it is false.
#[test]
fn an_empty_copy_neutral_set_calls_everything_neutral() {
    let segments = vec![segment(1, 100, 1.5), segment(101, 200, -1.5)];
    let statistics = call_copy_ratio_segments::calling_statistics(
        &segments,
        DEFAULT_NEUTRAL_LOWER_BOUND,
        DEFAULT_NEUTRAL_UPPER_BOUND,
        DEFAULT_OUTLIER_Z_SCORE,
    );
    assert!(statistics.mean.is_nan());
    let calls = call_copy_ratio_segments::make_calls(
        &segments,
        DEFAULT_NEUTRAL_LOWER_BOUND,
        DEFAULT_NEUTRAL_UPPER_BOUND,
        DEFAULT_OUTLIER_Z_SCORE,
        DEFAULT_CALLING_Z_SCORE,
    )
    .expect("valid arguments");
    assert!(calls.iter().all(|call| *call == Call::Neutral));
}

/// The bounds are inclusive, and they are compared in copy-ratio space.
#[test]
fn a_segment_exactly_on_a_bound_is_neutral() {
    let calls = call_copy_ratio_segments::make_calls(
        &spread(),
        DEFAULT_NEUTRAL_LOWER_BOUND,
        DEFAULT_NEUTRAL_UPPER_BOUND,
        DEFAULT_OUTLIER_Z_SCORE,
        DEFAULT_CALLING_Z_SCORE,
    )
    .expect("valid arguments");
    // The second and third segments are at exactly 0.9 and exactly 1.1.
    assert_eq!(calls[1], Call::Neutral);
    assert_eq!(calls[2], Call::Neutral);
    assert_eq!(calls[4], Call::Amplification);
    assert_eq!(calls[5], Call::Deletion);
}

/// The two argument refusals the golden carries.
#[test]
fn both_argument_refusals_match_the_golden() {
    let text = golden();

    let error = call_copy_ratio_segments::make_calls(
        &spread(),
        DEFAULT_NEUTRAL_LOWER_BOUND,
        DEFAULT_NEUTRAL_UPPER_BOUND,
        DEFAULT_OUTLIER_Z_SCORE,
        0.0,
    )
    .expect_err("a threshold of zero");
    assert_eq!(error, CallerError::CallingThresholdNotPositive);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        file(&text, "error", "zero-threshold").trim_end()
    );

    let error = call_copy_ratio_segments::make_calls(
        &spread(),
        2.0,
        1.0,
        DEFAULT_OUTLIER_Z_SCORE,
        DEFAULT_CALLING_Z_SCORE,
    )
    .expect_err("inverted bounds");
    assert_eq!(error, CallerError::BoundsNotOrdered);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        file(&text, "error", "inverted-bounds").trim_end()
    );
}

/// The length is the interval's, inclusive at both ends.
#[test]
fn the_length_is_inclusive() {
    assert_eq!(segment(1, 100, 0.0).length(), 100);
    assert_eq!(segment(5, 5, 0.0).length(), 1);
}
