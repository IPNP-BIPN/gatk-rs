//! Conformance for the decomposition `CalculateContamination` waits on, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/KernelSegmentationDump.java`.
//!
//! # What this suite is for
//!
//!  * **the subsample is seeded**, so the matrix is deterministic and a port that drew differently
//!    would be decomposing something else entirely;
//!  * **the singular values are compared as raw bits**, because the segmenter uses them as
//!    `1 / (sqrt(s) + 1e-10)` and a last-bit difference near zero is a factor of forty;
//!  * **`U` is compared entry by entry**, including its sign convention, which reaches the reduced
//!    observation matrix column by column;
//!  * **and the rank deficiency is part of the claim**: two of the three series decompose to
//!    matrices with exact zeros among the singular values, which is where a different
//!    implementation would diverge first.
//!
//! `docs/what-the-kernel-segmenter-needs-from-the-decomposition.md` works out which of these a
//! *different* implementation could be held to. This port is a transcription, so it is held to all
//! of them.
//!
//! The changepoint rows are read too, and they are the point of the exercise: they are what
//! `ContaminationSegmenter` consumes, and the only thing the decomposition sends downstream.

use gatk_corpus as corpus;
use gatk_engine::kernel_segmenter::{find_changepoints, sub_kernel_matrix, ChangepointSortOrder};
use gatk_engine::singular_value_decomposition::SingularValueDecomposition;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/kernel_segmentation.txt.gz"),
    )
}

/// The dump's kernel: a Gaussian of the difference, through `Math.exp` rather than `FastMath.exp`.
fn gaussian(variance: f64) -> impl Fn(f64, f64) -> f64 + Copy {
    move |x: f64, y: f64| (-(x - y) * (x - y) / (2.0 * variance)).exp()
}

/// The three series the dump decomposes, in its order.
fn series(label: &str) -> Vec<f64> {
    match label {
        // A step function with two changes, which is what a segmenter is for.
        "two-steps" => (0..90)
            .map(|i| {
                if i < 30 {
                    0.1
                } else if i < 60 {
                    0.5
                } else {
                    0.2
                }
            })
            .collect(),
        // A flat series, which has no changepoint to find.
        "flat" => vec![0.3; 90],
        // A ramp, where every point is a little different from the last.
        "ramp" => (0..90).map(|i| i as f64 / 90.0).collect(),
        other => panic!("unknown series {other}"),
    }
}

fn decompose(label: &str) -> SingularValueDecomposition {
    let matrix = sub_kernel_matrix(&series(label), 6, gaussian(0.01));
    SingularValueDecomposition::new(&matrix)
}

/// `%016x` of the raw bits, which is how the dump prints a double.
fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

#[test]
fn every_singular_value_and_u_entry_matches_the_golden() {
    let text = golden();
    let (mut singular_rows, mut u_rows) = (0, 0);

    for label in ["two-steps", "flat", "ramp"] {
        let svd = decompose(label);

        for line in text.lines() {
            let Some(rest) = line.strip_prefix("singular\t") else {
                continue;
            };
            let (row_label, entry) = rest.split_once('\t').expect("a label");
            if row_label != label {
                continue;
            }
            let (index, expected) = entry.split_once('=').expect("an index");
            let index: usize = index.parse().expect("a number");
            let expected = expected.split(',').next().expect("the bits");
            assert_eq!(
                bits(svd.singular_values[index]),
                expected,
                "singular value {index} of {label}"
            );
            singular_rows += 1;
        }

        for line in text.lines() {
            let Some(rest) = line.strip_prefix("u\t") else {
                continue;
            };
            let (row_label, entry) = rest.split_once('\t').expect("a label");
            if row_label != label {
                continue;
            }
            let (position, expected) = entry.split_once('=').expect("a position");
            let (row, column) = position.split_once(',').expect("a row and a column");
            let row: usize = row.parse().expect("a number");
            let column: usize = column.parse().expect("a number");
            assert_eq!(
                bits(svd.u[row][column]),
                expected,
                "U[{row}][{column}] of {label}"
            );
            u_rows += 1;
        }
    }

    assert_eq!(singular_rows, 18, "the golden's singular value rows");
    assert_eq!(u_rows, 108, "the golden's U rows");
}

/// The rank deficiency is not incidental: it is the reason the decision document exists.
#[test]
fn the_flat_series_decomposes_to_a_rank_one_matrix() {
    let svd = decompose("flat");
    assert_eq!(svd.singular_values[0], 6.0, "every kernel value is one");
    for (index, value) in svd.singular_values.iter().enumerate().skip(1) {
        assert_eq!(*value, 0.0, "singular value {index} is an exact zero");
    }
}

/// The dump's changepoint cases: the series, the maximum, and the two penalty factors.
fn changepoint_case(label: &str) -> (&'static str, usize, f64, f64) {
    match label {
        "two-steps" => ("two-steps", 10, 1.0, 1.0),
        "flat" => ("flat", 10, 1.0, 1.0),
        "ramp" => ("ramp", 10, 1.0, 1.0),
        // The same series with no penalty at all, and with a penalty ten times the default.
        "two-steps-lenient" => ("two-steps", 10, 0.0, 0.0),
        "two-steps-strict" => ("two-steps", 10, 10.0, 10.0),
        other => panic!("unknown case {other}"),
    }
}

#[test]
fn every_changepoint_row_matches_the_golden() {
    let text = golden();
    let mut rows = 0;

    for line in text.lines() {
        let Some(rest) = line.strip_prefix("changepoints\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("a label");
        let (series_label, maximum, linear, log_linear) = changepoint_case(label);
        let found = find_changepoints(
            &series(series_label),
            maximum,
            gaussian(0.01),
            6,
            &[8, 16],
            linear,
            log_linear,
            ChangepointSortOrder::Index,
        );
        // The dump prints `(none)` rather than an empty field, so an empty answer is visible.
        let printed = if found.is_empty() {
            "(none)".to_string()
        } else {
            found
                .iter()
                .map(|index| index.to_string())
                .collect::<Vec<String>>()
                .join(",")
        };
        assert_eq!(printed, expected, "changepoints for {label}");
        rows += 1;
    }

    assert_eq!(rows, 5, "the golden's changepoint rows");
}
