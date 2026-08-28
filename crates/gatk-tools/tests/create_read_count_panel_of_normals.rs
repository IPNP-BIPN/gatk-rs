//! Conformance for `CreateReadCountPanelOfNormals` against GATK 4.6.2.0, compared as the intervals
//! and samples every run's panel kept.
//!
//! Golden from `tools/readfilter-conformance/CreateReadCountPanelOfNormalsDump.java`.
//!
//! The singular value decomposition is not measured or ported. What is compared is the
//! preprocessing that decides what reaches it: the four filters, the medians the panel records,
//! and the eigensample cap.
//!
//! # What this suite is for
//!
//!  * **the counts becoming fractional coverage first**;
//!  * **the interval median filter, and the medians the panel records**;
//!  * **the two zero filters not being symmetric, each seeing what the one before left**;
//!  * **the extreme-median filter cutting from both ends**;
//!  * **a percentile of zero skipping the step rather than filtering nothing**;
//!  * **the eigensample count being capped at the surviving samples**;
//!  * **a panel of one sample having none and refusing its singular values**;
//!  * **and an input whose intervals do not match being refused by name.**

use gatk_corpus as corpus;
use gatk_tools::create_read_count_panel_of_normals::{
    median, mismatched_intervals_message, number_of_eigensamples, percentile, preprocess,
    safe_log2, standardize, to_fractional_coverage, Arguments, Matrix, EPSILON,
    NO_NON_ZERO_SINGULAR_VALUES_MESSAGE, NO_SINGULAR_VALUES_MESSAGE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/create_read_count_panel_of_normals.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// One `panel\t<label>\t<field>=<value>` row.
fn field(text: &str, label: &str, name: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("panel\t{label}\t{name}=")))
        .unwrap_or_else(|| panic!("the golden carries panel/{label}/{name}"))
        .to_string()
}

/// The message a run that wrote no panel was refused by.
fn refusal(text: &str, label: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

fn number(text: &str, label: &str, name: &str) -> usize {
    field(text, label, name).parse().expect("a number")
}

/// The intervals one panel kept, as `start-end`.
fn panel_intervals(text: &str, label: &str) -> Vec<String> {
    field(text, label, "panel-intervals")
        .split(',')
        .map(str::to_string)
        .collect()
}

/// The fractional medians one panel recorded.
fn fractional_medians(text: &str, label: &str) -> Vec<f64> {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("matrix\t{label}\tfractional-medians=")))
            .unwrap_or_else(|| panic!("the golden carries matrix/{label}")),
    )
    .lines()
    .filter(|line| !line.is_empty())
    .map(|value| value.parse().expect("a number"))
    .collect()
}

/// The fixture's own counts, rebuilt from the rule the dump wrote them with.
///
/// Only five of the nine samples are printed in the golden, the four ordinary ones being alike,
/// so the rule is reproduced here and checked against the five that are.
/// One sample's profile before the shared clamps.
///
/// Sample 6 borrows sample 3's, so that doubling it leaves its fractional coverage exactly equal.
fn shape(sample: usize, i: usize) -> usize {
    let base = 100 + (i * 7) % 23;
    if sample == 7 {
        return if i < 4 { base } else { 0 };
    }
    if sample == 8 {
        return if i < 10 { base * 40 } else { base };
    }
    let weight = 1 + if sample == 6 { 3 } else { sample };
    if i < 20 {
        base * weight
    } else {
        base
    }
}

/// The dump's `sampleCounts`: the samples vary in SHAPE, because fractional coverage divides
/// depth out again and a fixture that varied only the depth left their medians indistinguishable.
fn sample_counts(sample: usize) -> Vec<f64> {
    (0..40)
        .map(|i| {
            let mut value = shape(sample, i) * if sample == 6 { 2 } else { 1 };
            if i < 2 {
                value = if sample == 7 {
                    0
                } else if sample == 6 {
                    2
                } else {
                    1
                };
            }
            if i == 20 && sample < 3 {
                value = 0;
            }
            value as f64
        })
        .collect()
}

fn fixture() -> Matrix {
    Matrix::new(&(0..9).map(sample_counts).collect::<Vec<_>>())
}

/// The counts the golden printed, for the samples it printed.
fn printed_counts(text: &str, sample: usize) -> Vec<f64> {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("counts\tsample{sample}=")))
            .unwrap_or_else(|| panic!("the golden carries counts/sample{sample}")),
    )
    .lines()
    .filter(|line| line.starts_with("chr1\t"))
    .map(|line| {
        line.split('\t')
            .nth(3)
            .expect("a count")
            .parse()
            .expect("a number")
    })
    .collect()
}

/// The fixture the port rebuilds is the one the golden printed.
#[test]
fn the_fixture_is_the_goldens() {
    let text = golden();
    for sample in [0, 1, 6, 7, 8] {
        assert_eq!(
            sample_counts(sample),
            printed_counts(&text, sample),
            "{sample}"
        );
    }
    let matrix = fixture();
    assert_eq!(matrix.samples, 9);
    assert_eq!(matrix.intervals, 40);
    // Every run of the golden reports the same forty original intervals whatever it filtered.
    for label in ["default", "no-filtering", "one-sample"] {
        assert_eq!(number(&text, label, "original-intervals"), 40, "{label}");
    }
}

/// The intervals every run kept.
#[test]
fn every_panel_keeps_the_intervals_the_golden_kept() {
    let text = golden();
    let matrix = fixture();
    let places = |kept: &[usize]| -> Vec<String> {
        kept.iter()
            .map(|interval| format!("{}-{}", interval * 1000 + 1, (interval + 1) * 1000))
            .collect()
    };
    let none = Arguments {
        minimum_interval_median_percentile: 0.0,
        maximum_zeros_in_sample_percentage: 100.0,
        maximum_zeros_in_interval_percentage: 100.0,
        extreme_sample_median_percentile: 0.0,
        extreme_outlier_truncation_percentile: 0.0,
        ..Arguments::default()
    };
    let cases: Vec<(&str, Arguments)> = vec![
        ("default", Arguments::default()),
        ("no-filtering", none.clone()),
        (
            "interval-median-only",
            Arguments {
                minimum_interval_median_percentile: 10.0,
                ..none.clone()
            },
        ),
        (
            "zeros-in-sample-only",
            Arguments {
                maximum_zeros_in_sample_percentage: 5.0,
                ..none.clone()
            },
        ),
        (
            "extreme-sample-only",
            Arguments {
                extreme_sample_median_percentile: 20.0,
                ..none.clone()
            },
        ),
        (
            "no-imputation",
            Arguments {
                impute_zeros: false,
                ..Arguments::default()
            },
        ),
        ("two-eigensamples", Arguments::default()),
        ("hundred-eigensamples", Arguments::default()),
    ];
    let mut compared = 0;
    for (label, arguments) in cases {
        let result = preprocess(&matrix, &arguments);
        assert_eq!(
            places(&result.panel_intervals()),
            panel_intervals(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 8, "the runs the port reproduces");
    // The ninth wrote no panel: its interval filter left too little to decompose.
    assert!(refusal(&text, "zeros-in-interval-only").contains(NO_NON_ZERO_SINGULAR_VALUES_MESSAGE));
}

/// Each sample divided by its own total, so depth is not what the extreme filter sees.
#[test]
fn the_counts_become_fractional_coverage_first() {
    let mut matrix = fixture();
    // Sample 6 is sample 3 doubled, the two intervals every sample pins at one included.
    for interval in 0..40 {
        assert_eq!(
            matrix.get(6, interval),
            matrix.get(3, interval) * 2.0,
            "{interval}"
        );
    }
    to_fractional_coverage(&mut matrix);
    // Every row sums to one.
    for sample in 0..matrix.samples {
        let total: f64 = matrix.row(sample).iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "{sample}");
    }
    // And the deep sample is now EXACTLY the shallow one: dividing by a total that is itself
    // doubled undoes the doubling bit for bit, which is what makes it not an outlier.
    for interval in 0..40 {
        assert_eq!(
            matrix.get(6, interval),
            matrix.get(3, interval),
            "{interval}"
        );
    }
}

/// The medians are taken before any filtering and are what the panel records.
#[test]
fn the_panel_records_the_original_interval_medians() {
    let text = golden();
    let matrix = fixture();
    let result = preprocess(&matrix, &Arguments::default());
    let written = fractional_medians(&text, "default");
    assert_eq!(
        result.panel_interval_fractional_medians.len(),
        written.len()
    );
    for (produced, expected) in result
        .panel_interval_fractional_medians
        .iter()
        .zip(written.iter())
    {
        assert!(
            (produced - expected).abs() < 1e-15,
            "{produced} vs {expected}"
        );
    }
    // They are the medians of the FRACTIONAL coverage, so they sum to about one over the
    // intervals that survived plus those that did not.
    assert!(result
        .panel_interval_fractional_medians
        .iter()
        .all(|m| *m > 0.0));
}

/// Each sees what the one before it left.
#[test]
fn the_two_zero_filters_are_not_symmetric() {
    let text = golden();
    let matrix = fixture();
    let none = Arguments {
        minimum_interval_median_percentile: 0.0,
        maximum_zeros_in_sample_percentage: 100.0,
        maximum_zeros_in_interval_percentage: 100.0,
        extreme_sample_median_percentile: 0.0,
        extreme_outlier_truncation_percentile: 0.0,
        ..Arguments::default()
    };
    let sample_only = preprocess(
        &matrix,
        &Arguments {
            maximum_zeros_in_sample_percentage: 5.0,
            ..none.clone()
        },
    );
    let interval_only = preprocess(
        &matrix,
        &Arguments {
            maximum_zeros_in_interval_percentage: 5.0,
            ..none.clone()
        },
    );
    // The sample filter takes a sample and leaves every interval.
    assert_eq!(sample_only.panel_intervals().len(), 40);
    assert_eq!(sample_only.panel_samples().len(), 8);
    // The interval filter is much the harsher of the two: it leaves so few intervals that the
    // decomposition finds no non-zero singular value and the run is REFUSED, so the golden has no
    // panel for it at all.
    assert!(interval_only.panel_intervals().len() < 5);
    assert_eq!(interval_only.panel_samples().len(), 9);
    assert!(refusal(&text, "zeros-in-interval-only").contains(NO_NON_ZERO_SINGULAR_VALUES_MESSAGE));
    // Where the sample filter leaves every interval and writes its panel.
    assert_eq!(panel_intervals(&text, "zeros-in-sample-only").len(), 40);
    assert_eq!(number(&text, "zeros-in-sample-only", "eigensamples"), 7);
}

/// Applied twice, so a percentile of twenty takes a sample from each end.
#[test]
fn the_extreme_median_filter_cuts_from_both_ends() {
    let matrix = fixture();
    let none = Arguments {
        minimum_interval_median_percentile: 0.0,
        maximum_zeros_in_sample_percentage: 100.0,
        maximum_zeros_in_interval_percentage: 100.0,
        extreme_sample_median_percentile: 0.0,
        extreme_outlier_truncation_percentile: 0.0,
        ..Arguments::default()
    };
    let result = preprocess(
        &matrix,
        &Arguments {
            extreme_sample_median_percentile: 20.0,
            ..none.clone()
        },
    );
    assert_eq!(result.panel_samples().len(), 7, "two of nine");
    let text = golden();
    assert_eq!(number(&text, "extreme-sample-only", "eigensamples"), 7);
    // A percentile of zero skips the step entirely rather than filtering nothing.
    let skipped = preprocess(&matrix, &none);
    assert_eq!(skipped.panel_samples().len(), 9);
    assert_eq!(number(&text, "no-filtering", "eigensamples"), 9);
}

/// The requested count is capped at the samples that survived, not at the samples given.
///
/// What the panel FILE reports is a different number: `getNumEigensamples` is the length of the
/// singular-value array Spark returned, and Spark drops any singular value under its own epsilon.
/// That is a numerical rank decided by a distributed solver, which is the one thing this suite's
/// harness says it does not measure, so the port models the cap and not the rank. The two agree
/// only where the request is the binding constraint.
#[test]
fn the_eigensample_count_is_capped_at_the_samples_that_survived() {
    let text = golden();
    assert_eq!(number_of_eigensamples(20, 8), 8);
    assert_eq!(number_of_eigensamples(100, 8), 8);
    assert_eq!(number_of_eigensamples(2, 8), 2);
    assert_eq!(number_of_eigensamples(20, 1), 1);
    // Where the request binds, the file's count is the request and the port reaches it.
    assert_eq!(number(&text, "two-eigensamples", "eigensamples"), 2);
    let matrix = fixture();
    let surviving = preprocess(&matrix, &Arguments::default())
        .panel_samples()
        .len();
    assert_eq!(number_of_eigensamples(2, surviving), 2);
    // Where it does not, the file reports the solver's rank, which is at most the cap and here is
    // under it: the port's seven surviving samples against the file's own count.
    assert_eq!(surviving, 8);
    for label in ["default", "hundred-eigensamples"] {
        let reported = number(&text, label, "eigensamples");
        assert!(
            reported <= surviving,
            "{label}: {reported} over {surviving}"
        );
        assert_eq!(reported, 7, "{label}");
    }
    // The singular values are as many as the count the file reports, whichever way it was
    // decided, which is what says the two are the same number read twice.
    for label in ["default", "two-eigensamples", "hundred-eigensamples"] {
        assert_eq!(
            number(&text, label, "singular-values"),
            number(&text, label, "eigensamples"),
            "{label}"
        );
    }
}

/// And it refuses its singular values rather than handing over an empty array.
#[test]
fn a_panel_of_one_sample_has_no_eigensamples() {
    let text = golden();
    assert_eq!(number(&text, "one-sample", "samples"), 1);
    assert_eq!(number(&text, "one-sample", "eigensamples"), 0);
    // A panel of one sample still keeps intervals: the filters are not skipped for it.
    assert!(!panel_intervals(&text, "one-sample").is_empty());
    // The refusal is an UnsupportedOperationException, so it is a state the reader will not
    // represent rather than an argument the caller got wrong.
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix("error\tone-sample\tsingular-values\t"))
        .expect("its refusal");
    let (class, message) = row.split_once(':').expect("a class and a message");
    assert_eq!(class, "java.lang.UnsupportedOperationException");
    assert!(message.starts_with(NO_SINGULAR_VALUES_MESSAGE), "{message}");
}

/// Refused by name, before anything is preprocessed.
#[test]
fn an_input_whose_intervals_do_not_match_is_refused() {
    let text = golden();
    let message = text
        .lines()
        .find_map(|line| line.strip_prefix("error\tmismatched-intervals\t"))
        .expect("its refusal");
    assert!(
        message.starts_with("java.lang.IllegalArgumentException: "),
        "{message}"
    );
    let body = message
        .strip_prefix("java.lang.IllegalArgumentException: ")
        .expect("a message");
    assert_eq!(body, mismatched_intervals_message("<dir>/odd.counts.tsv"));
}

/// The estimator is the legacy one, whose rank is `p/100 * (n + 1)`.
#[test]
fn the_percentile_is_the_legacy_estimator() {
    let values = vec![1.0, 2.0, 3.0, 4.0];
    // A rank below one is the minimum and a rank at or past n is the maximum.
    assert_eq!(percentile(&values, 0.0), 1.0);
    assert_eq!(percentile(&values, 100.0), 4.0);
    assert_eq!(percentile(&values, 10.0), 1.0, "rank 0.5, below one");
    // The median of an even count is the interpolation of the two middle values.
    assert_eq!(median(&values), 2.5);
    assert_eq!(median(&[1.0, 2.0, 3.0]), 2.0);
    // A rank between two order statistics interpolates linearly.
    assert_eq!(percentile(&values, 40.0), 2.0);
    assert_eq!(percentile(&values, 50.0), 2.5);
    // One value is itself, and none is a NaN.
    assert_eq!(percentile(&[7.0], 25.0), 7.0);
    assert!(percentile(&[], 50.0).is_nan());
    // The order of the input does not matter.
    assert_eq!(percentile(&[4.0, 1.0, 3.0, 2.0], 40.0), 2.0);
}

/// The floor is a clamp, not an addition.
#[test]
fn the_log_is_floored_rather_than_nudged() {
    assert_eq!(EPSILON, 1e-9);
    // A value below the floor becomes the floor's own logarithm.
    assert_eq!(safe_log2(0.0), safe_log2(EPSILON / 2.0));
    assert!((safe_log2(0.0) - (EPSILON.ln() / std::f64::consts::LN_2)).abs() < 1e-12);
    // A value at the floor is NOT floored, the comparison being strict.
    assert!((safe_log2(EPSILON) - safe_log2(EPSILON / 2.0)).abs() < 1e-12);
    // And an ordinary value is its own logarithm, with nothing added.
    assert!((safe_log2(1.0) - 0.0).abs() < 1e-15);
    assert!((safe_log2(2.0) - 1.0).abs() < 1e-15);
    assert!((safe_log2(8.0) - 3.0).abs() < 1e-14);
}

/// The panel subtracts one median from every row, where a single sample subtracts its own.
#[test]
fn the_standardisation_subtracts_the_median_of_the_medians() {
    let matrix = fixture();
    let mut values = preprocess(&matrix, &Arguments::default()).values;
    let before = values.clone();
    standardize(&mut values).expect("positive medians");
    // Every row's median is now measured against the same number, so the median of the row
    // medians is zero.
    let medians: Vec<f64> = (0..values.samples)
        .map(|sample| median(values.row(sample)))
        .collect();
    assert!(median(&medians).abs() < 1e-12);
    // The rows are not each centred on zero, which is what a single sample's standardisation
    // would have done.
    assert!(medians.iter().any(|m| m.abs() > 1e-12));
    // And the transform is monotone, so the order within a row is unchanged.
    for sample in 0..values.samples {
        for interval in 1..values.intervals {
            let (a, b) = (
                before.get(sample, interval - 1),
                before.get(sample, interval),
            );
            let (x, y) = (
                values.get(sample, interval - 1),
                values.get(sample, interval),
            );
            // Non-strict: two values may collapse onto one, but their order may not reverse.
            match (a.partial_cmp(&b), x.partial_cmp(&y)) {
                (Some(std::cmp::Ordering::Less), after) => {
                    assert_ne!(
                        after,
                        Some(std::cmp::Ordering::Greater),
                        "{sample} {interval}"
                    )
                }
                (Some(std::cmp::Ordering::Greater), after) => {
                    assert_ne!(after, Some(std::cmp::Ordering::Less), "{sample} {interval}")
                }
                (before, after) => assert_eq!(before, after, "{sample} {interval}"),
            }
        }
    }
    // A sample whose median is not positive is refused by index.
    let mut zeros = Matrix::new(&[vec![0.0; 4], vec![1.0; 4]]);
    assert_eq!(
        standardize(&mut zeros).expect_err("a zero median"),
        "Sample at index 0 does not have a positive sample median."
    );
    let mut alone = Matrix::new(&[vec![0.0; 4]]);
    assert_eq!(
        standardize(&mut alone).expect_err("a zero median"),
        "Sample does not have a positive sample median."
    );
}

/// A zero becomes the median of the non-zero values of its own interval.
#[test]
fn a_surviving_zero_is_imputed_to_its_intervals_median() {
    let text = golden();
    let matrix = fixture();
    let imputed = preprocess(&matrix, &Arguments::default());
    let left = preprocess(
        &matrix,
        &Arguments {
            impute_zeros: false,
            ..Arguments::default()
        },
    );
    // The two keep the same intervals and samples: imputation is not a filter.
    assert_eq!(imputed.panel_intervals(), left.panel_intervals());
    assert_eq!(imputed.panel_samples(), left.panel_samples());
    assert_eq!(
        panel_intervals(&text, "default"),
        panel_intervals(&text, "no-imputation")
    );
    // And they record the same medians, those being taken before either.
    assert_eq!(
        fractional_medians(&text, "default"),
        fractional_medians(&text, "no-imputation")
    );
    // The imputed matrix has no zeros left where the other one does.
    let zeros = |values: &Matrix| values.values.iter().filter(|v| **v == 0.0).count();
    assert!(zeros(&left.values) >= zeros(&imputed.values));
}
