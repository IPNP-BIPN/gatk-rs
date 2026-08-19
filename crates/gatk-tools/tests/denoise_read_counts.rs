//! Conformance for `DenoiseReadCounts` with no panel of normals, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/DenoiseReadCountsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the scale-freedom** the first step buys: two inputs differing by a factor of two produce
//!    the same file;
//!  * **the floor**, which is why a zero count reads -29.897353 and not minus infinity;
//!  * **the interpolating median**, which an even number of intervals is there to show;
//!  * **and the two files being identical** when there is no panel.

use gatk_corpus as corpus;
use gatk_tools::denoise_read_counts::{self, EPSILON};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/denoise_read_counts.txt.gz"),
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

const SEQUENCES: [(&str, i32); 1] = [("chr1", 1000)];

fn sequences() -> Vec<(String, i32)> {
    SEQUENCES
        .iter()
        .map(|(name, length)| (name.to_string(), *length))
        .collect()
}

/// The dump's intervals: one hundred bases each, from position one.
fn intervals(count: usize) -> Vec<(String, i32, i32)> {
    (0..count)
        .map(|index| {
            let start = 1 + 100 * index as i32;
            ("chr1".to_string(), start, start + 99)
        })
        .collect()
}

fn table(counts: &[f64]) -> String {
    let values = denoise_read_counts::standardize(counts).expect("a positive median");
    denoise_read_counts::write(&sequences(), "SAMPLE", &intervals(counts.len()), &values)
}

#[test]
fn every_standardized_table_matches_the_golden() {
    let text = golden();

    for (label, counts) in [
        ("plain", vec![10.0, 20.0, 30.0, 40.0, 50.0]),
        ("doubled", vec![20.0, 40.0, 60.0, 80.0, 100.0]),
        ("even", vec![10.0, 20.0, 30.0, 40.0]),
        ("with-zero", vec![0.0, 20.0, 30.0]),
        ("flat", vec![7.0, 7.0, 7.0]),
        ("single", vec![42.0]),
    ] {
        let mine = table(&counts);
        assert_eq!(mine, file(&text, "standardized", label), "{label}");
        // With no panel the denoised file is the standardized one.
        assert_eq!(mine, file(&text, "denoised", label), "{label}: denoised");
    }
}

/// Doubling every count changes nothing, because the first step divides by the sum.
#[test]
fn the_output_is_scale_free() {
    let plain = table(&[10.0, 20.0, 30.0, 40.0, 50.0]);
    let doubled = table(&[20.0, 40.0, 60.0, 80.0, 100.0]);
    assert_eq!(plain, doubled);
    let tenfold = table(&[100.0, 200.0, 300.0, 400.0, 500.0]);
    assert_eq!(plain, tenfold);
}

/// A zero count is floored rather than infinite, and the floor joins the median.
#[test]
fn a_zero_count_is_floored() {
    assert_eq!(
        denoise_read_counts::safe_log2(0.0),
        denoise_read_counts::ln2_epsilon()
    );
    assert!(denoise_read_counts::ln2_epsilon().is_finite());
    let values = denoise_read_counts::standardize(&[0.0, 20.0, 30.0]).expect("a positive median");
    assert!(
        (values[0] + 29.897_352_853_986_263).abs() < 1e-9,
        "{:?}",
        values[0]
    );
    // Just under the epsilon is the floor; just over it is not.
    assert_eq!(
        denoise_read_counts::safe_log2(EPSILON / 2.0),
        denoise_read_counts::ln2_epsilon()
    );
    assert_ne!(
        denoise_read_counts::safe_log2(EPSILON * 2.0),
        denoise_read_counts::ln2_epsilon()
    );
}

/// An even number of intervals interpolates the median, so no value is exactly zero.
#[test]
fn an_even_count_interpolates_the_median() {
    let values = denoise_read_counts::standardize(&[10.0, 20.0, 30.0, 40.0]).expect("a median");
    assert!(
        values.iter().all(|value| *value != 0.0),
        "the median is between two values, {values:?}"
    );
    // And an odd count has exactly one zero, at the middle.
    let values = denoise_read_counts::standardize(&[10.0, 20.0, 30.0]).expect("a median");
    assert_eq!(values[1], 0.0);
}

/// Every count equal makes every standardized value zero.
#[test]
fn a_flat_row_standardizes_to_zeroes() {
    let values = denoise_read_counts::standardize(&[7.0, 7.0, 7.0]).expect("a median");
    assert_eq!(values, vec![0.0, 0.0, 0.0]);
}

/// An all-zero row divides by a sum of zero, so the median is a NaN -- and `isPositive` refuses a
/// NaN, which `median <= 0` would not.
#[test]
fn an_all_zero_row_is_refused() {
    let text = golden();
    let error = denoise_read_counts::standardize(&[0.0, 0.0, 0.0]).expect_err("a NaN median");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        text.lines()
            .find_map(|line| line.strip_prefix("error\tall-zero\t"))
            .expect("the golden carries the refusal")
    );
    // The median really is a NaN, so the test that catches it has to be the negated `>`.
    let fractional = denoise_read_counts::transform_to_fractional_coverage(&[0.0, 0.0, 0.0]);
    assert!(denoise_read_counts::median(&fractional).is_nan());
}
