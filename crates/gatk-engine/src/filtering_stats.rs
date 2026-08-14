//! `FilterStats`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterStats` (GATK 4.6.2.0).
//!
//! The `.filteringStats.tsv` that `FilterMutectCalls` writes beside its VCF: some metadata lines,
//! five columns, and one row per filter.
//!
//! # The two roundings disagree
//!
//! ```java
//! .set(M2FilterStatsTableColumn.FALSE_POSITIVE_COUNT.toString(), stats.getFalsePositiveCount(), 2)
//! ...
//! writer.writeMetadata(THRESHOLD_METADATA_TAG, Double.toString(round(threshold)));
//! ```
//!
//! where `round` is `roundToNDecimalPlaces(x, 3)`. The columns keep two decimals and the metadata
//! three, so the same number appears twice in one file at two precisions.
//!
//! # Zero, infinity and the sign
//!
//! Both roundings are `Math.round((x + ulp(x)) * mult) / mult`, which is not a rounding function in
//! the usual sense at the edges:
//!
//!  * `Math.round` answers `0` for NaN, so a run with no passing call divides by zero and writes
//!    `0.0` rather than `NaN`;
//!  * it saturates at `Long.MAX_VALUE` for infinity, so a false-positive count over zero calls
//!    writes `9.223372036854776E15`;
//!  * and it is `floor(x + 0.5)`, so it leans **upwards** rather than away from zero: `1.005` comes
//!    out `1.01` while `-1.005` comes out `-1.0`, and `-0.0005` loses its sign entirely.
//!
//! # The rows are the caller's
//!
//! This writer filters nothing and sorts nothing. The rule that drops a filter which accounted for
//! neither a false positive nor a false negative lives in `FilteringOutputStats`, which needs the
//! filtering engine and is not ported here.

use crate::base_recalibration_engine::round_to_n_decimal_places;
use crate::tsv_table::{java_double_to_string, write_table};

/// `THRESHOLD_METADATA_TAG`, `FDR_METADATA_TAG` and `SENSITIVITY_METADATA_TAG`.
pub const THRESHOLD_METADATA_TAG: &str = "threshold";
pub const FDR_METADATA_TAG: &str = "fdr";
pub const SENSITIVITY_METADATA_TAG: &str = "sensitivity";

/// `M2FilterStatsTableColumn`, in the order the enum declares.
pub const COLUMNS: [&str; 5] = ["filter", "FP", "FDR", "FN", "FNR"];

/// One row: a filter and the four numbers attributed to it.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterStats {
    pub filter_name: String,
    pub false_positive_count: f64,
    pub false_discovery_rate: f64,
    pub false_negative_count: f64,
    pub false_negative_rate: f64,
}

/// A column value: `DataLine.set(name, value, 2)`.
pub fn column(value: f64) -> String {
    java_double_to_string(
        round_to_n_decimal_places(value, 2).expect("two places is more than zero"),
    )
}

/// A metadata value: `Double.toString(roundToNDecimalPlaces(x, 3))`.
pub fn metadata(value: f64) -> String {
    java_double_to_string(
        round_to_n_decimal_places(value, 3).expect("three places is more than zero"),
    )
}

/// `writeM2FilterSummary`: the whole file.
///
/// The three totals are divided here rather than by the caller, which is what puts `0/0` and
/// `x/0` in front of the rounding.
#[allow(clippy::too_many_arguments)]
pub fn write_summary(
    filter_stats: &[FilterStats],
    clustering_metadata: &[(String, String)],
    threshold: f64,
    total_calls: f64,
    expected_true_positives: f64,
    expected_false_positives: f64,
    expected_false_negatives: f64,
) -> String {
    let mut pairs: Vec<(String, String)> = clustering_metadata.to_vec();
    pairs.push((THRESHOLD_METADATA_TAG.to_string(), metadata(threshold)));
    pairs.push((
        FDR_METADATA_TAG.to_string(),
        metadata(expected_false_positives / total_calls),
    ));
    pairs.push((
        SENSITIVITY_METADATA_TAG.to_string(),
        metadata(expected_true_positives / (expected_true_positives + expected_false_negatives)),
    ));

    let rows: Vec<Vec<String>> = filter_stats
        .iter()
        .map(|stats| {
            vec![
                stats.filter_name.clone(),
                column(stats.false_positive_count),
                column(stats.false_discovery_rate),
                column(stats.false_negative_count),
                column(stats.false_negative_rate),
            ]
        })
        .collect();

    let borrowed: Vec<(&str, &str)> = pairs
        .iter()
        .map(|(key, value)| (key.as_str(), value.as_str()))
        .collect();
    write_table(&COLUMNS, &rows, &borrowed)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stats(name: &str, values: [f64; 4]) -> FilterStats {
        FilterStats {
            filter_name: name.to_string(),
            false_positive_count: values[0],
            false_discovery_rate: values[1],
            false_negative_count: values[2],
            false_negative_rate: values[3],
        }
    }

    #[test]
    fn the_two_roundings_disagree() {
        // The same number, two decimals in a column and three in the metadata.
        assert_eq!(column(0.123456), "0.12");
        assert_eq!(metadata(0.123456), "0.123");
        // And the half goes up.
        assert_eq!(column(1.005), "1.01");
        assert_eq!(column(2.675), "2.68");
        assert_eq!(column(0.005), "0.01");
    }

    #[test]
    fn zero_infinity_and_the_sign() {
        // 0/0, computed rather than written, which is how the tool reaches it.
        let zero = std::hint::black_box(0.0_f64);
        assert!((zero / zero).is_nan());
        assert_eq!(metadata(zero / zero), "0.0");
        assert_eq!(column(f64::NAN), "0.0");
        // x/0.
        let one = std::hint::black_box(1.0_f64);
        assert!((one / zero).is_infinite());
        assert_eq!(metadata(one / zero), "9.223372036854776E15");
        // `floor(x + 0.5)` leans upwards rather than away from zero.
        assert_eq!(column(-1.005), "-1.0");
        assert_eq!(column(-0.125), "-0.12");
        assert_eq!(column(-2.5), "-2.5");
        // And a small negative loses its sign.
        assert_eq!(metadata(-0.0005), "0.0");
        assert_eq!(column(-0.0), "0.0");
    }

    #[test]
    fn the_metadata_comes_before_the_columns_in_the_order_it_was_given() {
        let text = write_summary(
            &[stats("weak_evidence", [3.0, 0.15, 2.0, 0.08])],
            &[("clustering".to_string(), "1".to_string())],
            0.234,
            20.0,
            25.0,
            3.0,
            5.0,
        );
        assert_eq!(
            text,
            "#<METADATA>clustering=1\n\
             #<METADATA>threshold=0.234\n\
             #<METADATA>fdr=0.15\n\
             #<METADATA>sensitivity=0.833\n\
             filter\tFP\tFDR\tFN\tFNR\n\
             weak_evidence\t3.0\t0.15\t2.0\t0.08\n"
        );
    }

    #[test]
    fn the_rows_are_the_callers_and_a_file_may_have_none() {
        let text = write_summary(&[], &[], 0.1, 10.0, 9.0, 1.0, 2.0);
        assert!(text.ends_with("filter\tFP\tFDR\tFN\tFNR\n"));
        // Nothing is dropped and nothing is sorted: two rows come back in the order given.
        let text = write_summary(
            &[
                stats("second", [1.0, 0.1, 1.0, 0.1]),
                stats("first", [0.0, 0.0, 0.0, 0.0]),
            ],
            &[],
            0.2,
            10.0,
            9.0,
            1.0,
            2.0,
        );
        let names: Vec<&str> = text
            .lines()
            .skip_while(|line| line.starts_with('#'))
            .skip(1)
            .filter_map(|line| line.split('\t').next())
            .collect();
        assert_eq!(names, vec!["second", "first"]);
    }

    #[test]
    fn a_name_that_needs_quoting_gets_it() {
        let text = write_summary(
            &[
                stats("has\ttab", [1.0, 0.1, 1.0, 0.1]),
                stats("has\"quote", [1.0, 0.1, 1.0, 0.1]),
                stats("has,comma", [1.0, 0.1, 1.0, 0.1]),
            ],
            &[],
            0.2,
            10.0,
            9.0,
            1.0,
            2.0,
        );
        assert!(text.contains("\"has\ttab\"\t1.0"), "{text}");
        assert!(
            text.contains("\"has\"\"quote\"") || text.contains("\"has\\\"quote\""),
            "{text}"
        );
        // A comma is not a separator here, so it passes through bare.
        assert!(text.contains("\nhas,comma\t1.0"), "{text}");
    }
}
