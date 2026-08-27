//! `GroundTruthScorer`: every read scored against the reference it aligns to, and the report that
//! summarises them.
//!
//! The flow-based scoring is not ported. What is ported is the report the tool builds out of the
//! scores: the accumulators, the four table shapes, the phred each row carries, and the bins the
//! deviation and the base are folded into.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.groundtruth.GroundTruthScorer` in
//! GATK 4.6.2.0.

use crate::series_stats::SeriesStats;

/// The bounds the report is allocated with.
pub const QUAL_VALUE_MAX: usize = 60;
pub const HMER_VALUE_MAX: usize = 100;
/// `FlowBasedRead.DEFAULT_FLOW_ORDER.length() - 1`, the flow order being `TGCA`.
pub const BASE_VALUE_MAX: usize = 3;
pub const DEFAULT_FLOW_ORDER: &str = "TGCA";

/// `NORMALIZED_SCORE_THRESHOLD_DEFAULT`, which is NEGATIVE: the scores it bounds are.
pub const NORMALIZED_SCORE_THRESHOLD_DEFAULT: f64 = -0.1;
pub const DEFAULT_RATIO_THRESHOLD: f64 = 0.003;

/// The percentile columns the report carries when none are asked for.
pub const DEFAULT_QUALITY_PERCENTILES: &str = "10,25,50,75,90";

/// One cell of the report: how many observations were true and how many false.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Accumulator {
    pub false_count: u64,
    pub true_count: u64,
}

impl Accumulator {
    pub fn add(&mut self, value: bool) {
        if value {
            self.true_count += 1;
        } else {
            self.false_count += 1;
        }
    }

    pub fn count(&self) -> u64 {
        self.false_count + self.true_count
    }

    /// The FALSE rate, which is zero for an empty cell rather than a NaN.
    pub fn false_rate(&self) -> f64 {
        if self.count() == 0 {
            0.0
        } else {
            self.false_count as f64 / self.count() as f64
        }
    }
}

/// `deviationToBin`: `0,-1,1,-2,2...` become `0,1,2,3,4...`.
pub fn deviation_to_bin(deviation: i32) -> usize {
    if deviation >= 0 {
        (deviation * 2) as usize
    } else {
        ((-deviation * 2) - 1) as usize
    }
}

/// `binToDeviation`, which is NOT the inverse of `deviation_to_bin` as a number: it is a STRING,
/// and a positive deviation is written with a leading `+` while zero and the negatives are not.
pub fn bin_to_deviation(bin: usize) -> String {
    if bin == 0 {
        "0".to_string()
    } else if bin.is_multiple_of(2) {
        format!("+{}", bin / 2)
    } else {
        format!("{}", -((bin as i64 + 1) / 2))
    }
}

/// `binToBase`, which indexes the flow order.
pub fn bin_to_base(bin: usize) -> char {
    DEFAULT_FLOW_ORDER.chars().nth(bin).expect("a base")
}

/// The phred a rate becomes in the one-level table.
///
/// A rate of zero over a NON-EMPTY cell is replaced by the probability threshold before the
/// logarithm, so a cell that saw observations and no errors reports the threshold's phred rather
/// than zero. A rate of zero over an EMPTY cell is left alone and reports zero.
pub fn phred(rate: f64, count: u64, probability_threshold: f64) -> i64 {
    let effective = if rate == 0.0 && count != 0 && probability_threshold != 0.0 {
        probability_threshold
    } else {
        rate
    };
    if effective != 0.0 {
        (-10.0 * effective.log10()).ceil() as i64
    } else {
        0
    }
}

/// The four table names the report carries, built from the column names.
pub fn one_level_table_name(name: &str) -> String {
    format!("{name}Report")
}

pub fn two_level_table_name(first: &str, second: &str) -> String {
    format!("{first}_{second}Report")
}

pub fn four_level_table_name(first: &str, second: &str, third: &str, fourth: &str) -> String {
    format!("{first}_{second}_{third}_{fourth}_Report")
}

/// Whether a row is written, given where it sits and what it holds.
///
/// The ORIGIN is always written: `omit_zeros` skips an empty row only when at least one of its
/// indices is non-zero, so the first row of every table survives however empty it is.
pub fn keeps_row(omit_zeros: bool, indices: &[usize], count: u64) -> bool {
    if !omit_zeros {
        return true;
    }
    let at_origin = indices.iter().all(|index| *index == 0);
    at_origin || count != 0
}

/// How many rows the four-level table has when its zeros are kept.
///
/// It is the product of the four allocated dimensions and not of what was observed, which is why
/// the option that omits the zeros is not really optional.
pub fn four_level_row_count() -> usize {
    (QUAL_VALUE_MAX + 1)
        * (HMER_VALUE_MAX + 1)
        * deviation_to_bin(HMER_VALUE_MAX as i32 + 1)
        * (BASE_VALUE_MAX + 1)
}

/// The percentile table's columns, which are fixed except for the percentiles themselves.
pub const PERCENTILE_TABLE_NAME: &str = "PhredBinAccumulator";
pub const PERCENTILE_FIXED_COLUMNS: [&str; 7] =
    ["flow", "count", "min", "max", "mean", "median", "std"];

/// The percentile table's column names for a given `--quality-percentiles`.
pub fn percentile_columns(quality_percentiles: &str) -> Vec<String> {
    let mut columns: Vec<String> = PERCENTILE_FIXED_COLUMNS
        .iter()
        .map(|name| name.to_string())
        .collect();
    for percentile in quality_percentiles.split(',') {
        columns.push(format!("p{percentile}"));
    }
    columns
}

/// One flow's percentile row, which is a `SeriesStats` fed in PHRED space.
///
/// `addProb` takes a probability and stores `-10 * log10(p)`, so the series holds phreds and the
/// percentiles are over those rather than over the probabilities they came from.
#[derive(Debug, Clone, Default)]
pub struct PercentileReport {
    pub stats: SeriesStats,
}

impl PercentileReport {
    pub fn add_probability(&mut self, probability: f64) {
        self.stats.add(-10.0 * probability.log10());
    }

    /// One row of the percentile table, in the column order above.
    pub fn row(&self, index: usize, quality_percentiles: &str) -> Vec<f64> {
        let mut values = vec![
            index as f64,
            self.stats.count() as f64,
            self.stats.min(),
            self.stats.max(),
            self.stats.mean(),
            self.stats.median(),
            self.stats.std(),
        ];
        for percentile in quality_percentiles.split(',') {
            values.push(
                self.stats
                    .percentile(percentile.parse().expect("a percentile")),
            );
        }
        values
    }
}

/// The arguments that decide what is scored and what is written.
#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub use_softclipped_bases: bool,
    pub normalized_score_threshold: f64,
    pub add_mean_call: bool,
    pub no_output: bool,
    pub omit_zeros_from_report: bool,
    pub quality_percentiles: String,
    pub exclude_zero_flows: bool,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            use_softclipped_bases: false,
            normalized_score_threshold: NORMALIZED_SCORE_THRESHOLD_DEFAULT,
            add_mean_call: false,
            no_output: false,
            omit_zeros_from_report: false,
            quality_percentiles: DEFAULT_QUALITY_PERCENTILES.to_string(),
            exclude_zero_flows: false,
        }
    }
}

/// Whether a read's normalized score keeps it.
///
/// The comparison is a strict less-than against a NEGATIVE default, so a read is dropped when its
/// score falls below the threshold rather than above it.
pub fn keeps_read(normalized_score: f64, arguments: &Arguments) -> bool {
    normalized_score >= arguments.normalized_score_threshold
}

/// The two columns `--add-mean-call` appends to the CSV, in the order it appends them.
///
/// The probabilities come FIRST and the mean call second, which is the other way round from the
/// argument's own name.
pub const MEAN_CALL_COLUMNS: [&str; 2] = ["ReadProbs", "ReadMeanCall"];
