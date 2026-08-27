//! Conformance for `SeriesStats` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/SeriesStatsDump.java`.
//!
//! A suite for a PRIMITIVE rather than for a tool: this is what both ground-truth tools keep their
//! numbers in, and neither can be ported without it.
//!
//! # What this suite is for
//!
//!  * **a percentile being an observed value**, at a truncated index;
//!  * **-0.0 and 0.0 being separate bins**;
//!  * **a NaN poisoning the bounds and leaving the median alone**;
//!  * **the deviation dividing by the count**;
//!  * **the CSV's format following the ADD PATH rather than the values**;
//!  * **and the digest of an empty series claiming a minimum of nought.**

use gatk_corpus as corpus;
use gatk_tools::series_stats::{java_double_to_int, SeriesStats};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/series_stats.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn line(text: &str, kind: &str, name: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{name}"))
        .to_string()
}

/// The dump prints every double to ten decimal places, and Java writes a NaN as `NaN`.
fn d(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else {
        format!("{value:.10}")
    }
}

fn stat_line(stats: &SeriesStats) -> String {
    format!(
        "{},{},{},{},{},{},{},{}",
        stats.count(),
        stats.uniq(),
        d(stats.min()),
        d(stats.max()),
        d(stats.mean()),
        d(stats.median()),
        d(stats.std()),
        d(stats.last())
    )
}

fn pct_line(stats: &SeriesStats) -> String {
    [0.0, 10.0, 25.0, 50.0, 75.0, 90.0, 99.0, 100.0]
        .iter()
        .map(|p| d(stats.percentile(*p)))
        .collect::<Vec<String>>()
        .join(",")
}

/// `Double.toString` as the dump's `entry.getKey() + ":"` renders it.
fn bins_line(stats: &SeriesStats) -> String {
    stats
        .bins()
        .iter()
        .map(|(key, count)| {
            format!(
                "{}:{count}",
                gatk_engine::tsv_table::java_double_to_string(key.0)
            )
        })
        .collect::<Vec<String>>()
        .join(",")
}

/// The nine series the dump built, in the order it built them.
fn series() -> Vec<(&'static str, SeriesStats)> {
    let ints = |values: &[i32]| {
        let mut stats = SeriesStats::new();
        for value in values {
            stats.add_int(*value);
        }
        stats
    };
    let doubles = |values: &[f64]| {
        let mut stats = SeriesStats::new();
        for value in values {
            stats.add(*value);
        }
        stats
    };
    let mut mixed = SeriesStats::new();
    mixed.add_int(1);
    mixed.add_int(2);
    mixed.add(3.5);
    vec![
        ("empty", SeriesStats::new()),
        ("one", ints(&[7])),
        ("even", ints(&[1, 2, 8, 9])),
        ("odd", ints(&[1, 2, 8, 9, 100])),
        ("repeated", ints(&[5, 5, 5, 5, 5, 5, 5, 5, 5, 100])),
        ("zeros", doubles(&[0.0, -0.0, 1.0])),
        ("nan", doubles(&[1.0, f64::NAN, 2.0])),
        ("doubles", doubles(&[1.5, 2.25, 2.25, 10.0])),
        ("mixed", mixed),
    ]
}

#[test]
fn every_series_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, stats) in series() {
        assert_eq!(
            stat_line(&stats),
            line(&text, "stat", label),
            "{label}/stat"
        );
        assert_eq!(pct_line(&stats), line(&text, "pct", label), "{label}/pct");
        assert_eq!(
            bins_line(&stats),
            line(&text, "bins", label),
            "{label}/bins"
        );
        assert_eq!(
            stats.to_digest(),
            line(&text, "digest", label),
            "{label}/digest"
        );
        compared += 4;
    }
    assert_eq!(compared, 36, "the values the golden carries");
}

#[test]
fn every_csv_matches_the_golden() {
    let text = golden();
    let mut integers = SeriesStats::new();
    for value in [1, 2, 8, 9] {
        integers.add_int(value);
    }
    assert_eq!(integers.csv(), unescape(&line(&text, "csv", "int")));

    let mut doubles = SeriesStats::new();
    for value in [1.5, 2.25, 2.25, 10.0] {
        doubles.add(value);
    }
    assert_eq!(doubles.csv(), unescape(&line(&text, "csv", "double")));

    let mut mixed = SeriesStats::new();
    mixed.add_int(1);
    mixed.add_int(2);
    mixed.add(3.5);
    assert_eq!(mixed.csv(), unescape(&line(&text, "csv", "mixed")));

    // Two whole numbers added as DOUBLES, which are written `1.000000` and not `1`.
    let mut whole = SeriesStats::new();
    whole.add(1.0);
    whole.add(2.0);
    assert_eq!(whole.csv(), unescape(&line(&text, "csv", "whole")));
    assert!(!whole.is_int_keys(), "the add path, not the values");
    assert!(integers.is_int_keys());
}

/// The walk returns a bin key at a truncated index, so the median of four is the third smallest.
#[test]
fn a_percentile_is_an_observed_value() {
    let mut even = SeriesStats::new();
    for value in [1, 2, 8, 9] {
        even.add_int(value);
    }
    assert_eq!(even.median(), 8.0, "not the 5 an interpolation would give");
    // The index: (int)(4 * 50 / 100) = 2, the third bin.
    assert_eq!(even.percentile(25.0), 2.0);
    assert_eq!(even.percentile(0.0), 1.0);
    assert_eq!(even.percentile(100.0), 9.0, "past the last bin");

    // A bin holding nine of ten values is what the walk has to step over.
    let mut repeated = SeriesStats::new();
    for value in [5, 5, 5, 5, 5, 5, 5, 5, 5, 100] {
        repeated.add_int(value);
    }
    assert_eq!(repeated.median(), 5.0);
    assert_eq!(repeated.percentile(90.0), 100.0, "the tenth value");
    assert_eq!(repeated.uniq(), 2);

    // A single value short-circuits to the LAST added.
    let mut one = SeriesStats::new();
    one.add_int(7);
    assert_eq!(one.percentile(0.0), 7.0);
    assert_eq!(one.percentile(100.0), 7.0);
}

/// A TreeMap over Double keeps them apart, so `getUniq` counts three for two zeros and a one.
#[test]
fn a_negative_zero_is_its_own_bin() {
    let mut zeros = SeriesStats::new();
    zeros.add(0.0);
    zeros.add(-0.0);
    zeros.add(1.0);
    assert_eq!(zeros.uniq(), 3);
    // And -0.0 sorts FIRST, so it is the minimum the bins report.
    let first = zeros.bins().keys().next().expect("a bin").0;
    assert!(first.is_sign_negative() && first == 0.0);
    assert_eq!(zeros.percentile(0.0), 0.0);
    assert!(zeros.percentile(0.0).is_sign_negative());
    // The mean sums them as one third, because -0.0 + 0.0 is 0.0.
    assert!((zeros.mean() - 1.0 / 3.0).abs() < 1e-12);
}

/// It poisons min, max, mean and the deviation, and leaves the median a real value.
#[test]
fn a_nan_poisons_the_bounds_and_not_the_median() {
    let mut stats = SeriesStats::new();
    stats.add(1.0);
    stats.add(f64::NAN);
    stats.add(2.0);
    assert!(stats.min().is_nan());
    assert!(stats.max().is_nan());
    assert!(stats.mean().is_nan());
    assert!(stats.std().is_nan());
    assert_eq!(stats.median(), 2.0, "the walk found a real bin");
    // The NaN sorts after everything, which is why the median is the second bin and not the third.
    assert!(stats.bins().keys().next_back().expect("a bin").0.is_nan());
    assert_eq!(stats.uniq(), 3);
}

/// The count, not the count less one.
#[test]
fn the_deviation_is_the_population_one() {
    let mut stats = SeriesStats::new();
    for value in [1, 2, 8, 9] {
        stats.add_int(value);
    }
    let mean = stats.mean();
    let population: f64 = [1.0, 2.0, 8.0, 9.0]
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / 4.0;
    assert!((stats.std() - population.sqrt()).abs() < 1e-12);
    let sample = population * 4.0 / 3.0;
    assert!(
        (stats.std() - sample.sqrt()).abs() > 0.1,
        "and not the sample one"
    );
}

/// An empty series counts as integer-keyed, so the digest casts its NaNs to nought.
#[test]
fn the_digest_of_an_empty_series_claims_a_minimum_of_nought() {
    let text = golden();
    let empty = SeriesStats::new();
    assert!(empty.min().is_nan(), "every other reader is told NaN");
    assert!(empty.is_int_keys(), "both counts are zero");
    assert_eq!(java_double_to_int(f64::NAN), 0);
    assert_eq!(
        empty.to_digest(),
        "count=0, min=0, max=0, median=0, bin.count=0"
    );
    assert_eq!(empty.to_digest(), line(&text, "digest", "empty"));
    // A series with one double takes the other branch and reports its NaNs as NaN.
    let mut with_double = SeriesStats::new();
    with_double.add(1.0);
    assert!(!with_double.is_int_keys());
    assert!(with_double.to_digest().contains("1.000000"));
}
