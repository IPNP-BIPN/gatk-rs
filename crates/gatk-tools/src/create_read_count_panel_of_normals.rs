//! `CreateReadCountPanelOfNormals`: which intervals and which samples reach the panel.
//!
//! The singular value decomposition is a distributed solver's and is not ported. What is ported is
//! everything that decides what reaches it: the transform to fractional coverage, the four
//! filters in the order they run, the imputation, the truncation, and the standardisation.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.copynumber.denoising.SVDDenoisingUtils` and
//! `org.broadinstitute.hellbender.tools.copynumber.CreateReadCountPanelOfNormals` in GATK 4.6.2.0.

/// The defaults the tool ships.
pub const DEFAULT_MINIMUM_INTERVAL_MEDIAN_PERCENTILE: f64 = 10.0;
pub const DEFAULT_MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE: f64 = 5.0;
pub const DEFAULT_MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE: f64 = 5.0;
pub const DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE: f64 = 2.5;
pub const DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE: f64 = 0.1;
pub const DEFAULT_NUMBER_OF_EIGENSAMPLES: usize = 20;

/// The arguments that decide what survives.
#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub minimum_interval_median_percentile: f64,
    pub maximum_zeros_in_sample_percentage: f64,
    pub maximum_zeros_in_interval_percentage: f64,
    pub extreme_sample_median_percentile: f64,
    pub impute_zeros: bool,
    pub extreme_outlier_truncation_percentile: f64,
    pub number_of_eigensamples: usize,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            minimum_interval_median_percentile: DEFAULT_MINIMUM_INTERVAL_MEDIAN_PERCENTILE,
            maximum_zeros_in_sample_percentage: DEFAULT_MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE,
            maximum_zeros_in_interval_percentage: DEFAULT_MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE,
            extreme_sample_median_percentile: DEFAULT_EXTREME_SAMPLE_MEDIAN_PERCENTILE,
            impute_zeros: true,
            extreme_outlier_truncation_percentile: DEFAULT_EXTREME_OUTLIER_TRUNCATION_PERCENTILE,
            number_of_eigensamples: DEFAULT_NUMBER_OF_EIGENSAMPLES,
        }
    }
}

/// `org.apache.commons.math3.stat.descriptive.rank.Percentile`, in its default estimation type.
///
/// The default is the LEGACY type: the rank is `p/100 * (n + 1)`, a rank below one is the
/// minimum, a rank at or above `n` is the maximum, and anything between is interpolated linearly
/// between the two neighbouring order statistics. It is NOT the type most other libraries use.
pub fn percentile(values: &[f64], p: f64) -> f64 {
    let mut sorted: Vec<f64> = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return sorted[0];
    }
    let position = p / 100.0 * (n as f64 + 1.0);
    if position < 1.0 {
        return sorted[0];
    }
    if position >= n as f64 {
        return sorted[n - 1];
    }
    let floor = position.floor();
    let difference = position - floor;
    let lower = sorted[floor as usize - 1];
    let upper = sorted[floor as usize];
    lower + difference * (upper - lower)
}

/// `org.apache.commons.math3.stat.descriptive.rank.Median`, which is the fiftieth percentile of
/// the same estimator.
pub fn median(values: &[f64]) -> f64 {
    percentile(values, 50.0)
}

/// The read counts, samples by intervals.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    pub samples: usize,
    pub intervals: usize,
    /// Row-major: `values[sample * intervals + interval]`.
    pub values: Vec<f64>,
}

impl Matrix {
    pub fn new(rows: &[Vec<f64>]) -> Matrix {
        let samples = rows.len();
        let intervals = rows.first().map_or(0, |row| row.len());
        Matrix {
            samples,
            intervals,
            values: rows.iter().flatten().copied().collect(),
        }
    }

    pub fn get(&self, sample: usize, interval: usize) -> f64 {
        self.values[sample * self.intervals + interval]
    }

    pub fn set(&mut self, sample: usize, interval: usize, value: f64) {
        self.values[sample * self.intervals + interval] = value;
    }

    pub fn row(&self, sample: usize) -> &[f64] {
        &self.values[sample * self.intervals..(sample + 1) * self.intervals]
    }

    pub fn column(&self, interval: usize) -> Vec<f64> {
        (0..self.samples).map(|s| self.get(s, interval)).collect()
    }
}

/// `transformToFractionalCoverage`: each sample divided by its own total.
///
/// This is what makes a sample sequenced twice as deeply as another contribute the same, so
/// sequencing depth is not what the extreme-median filter is looking at.
pub fn to_fractional_coverage(counts: &mut Matrix) {
    for sample in 0..counts.samples {
        let total: f64 = counts.row(sample).iter().sum();
        for interval in 0..counts.intervals {
            let value = counts.get(sample, interval);
            counts.set(sample, interval, value / total);
        }
    }
}

/// What one run of the preprocessing filtered out and kept.
#[derive(Debug, Clone, PartialEq)]
pub struct Preprocessed {
    /// One flag per sample, true when filtered OUT.
    pub filtered_samples: Vec<bool>,
    /// One flag per interval, true when filtered OUT.
    pub filtered_intervals: Vec<bool>,
    /// The original median of each interval that survived, before the division.
    pub panel_interval_fractional_medians: Vec<f64>,
    /// The surviving submatrix, imputed and truncated.
    pub values: Matrix,
}

impl Preprocessed {
    pub fn panel_intervals(&self) -> Vec<usize> {
        (0..self.filtered_intervals.len())
            .filter(|i| !self.filtered_intervals[*i])
            .collect()
    }

    pub fn panel_samples(&self) -> Vec<usize> {
        (0..self.filtered_samples.len())
            .filter(|i| !self.filtered_samples[*i])
            .collect()
    }
}

/// `preprocessPanel`: the four filters, in the order they run, then the imputation and the
/// truncation.
///
/// The ORDER is the behaviour. The interval medians are taken BEFORE anything is filtered, and
/// the division by them uses those original medians whatever is dropped afterwards. The sample
/// zero filter then counts zeros over the intervals that survived the median filter, and the
/// interval zero filter counts zeros over the samples that survived the sample filter, so the two
/// are not symmetric: each sees what the one before it left.
pub fn preprocess(counts: &Matrix, arguments: &Arguments) -> Preprocessed {
    let mut counts = counts.clone();
    to_fractional_coverage(&mut counts);
    let samples = counts.samples;
    let intervals = counts.intervals;
    let mut filtered_samples = vec![false; samples];
    let mut filtered_intervals = vec![false; intervals];

    // The medians every later step divides by, taken before any filtering.
    let original_medians: Vec<f64> = (0..intervals)
        .map(|interval| median(&counts.column(interval)))
        .collect();

    // A percentile of zero SKIPS the step rather than filtering nothing, which is not the same:
    // a threshold of zero would still drop an interval whose median is zero.
    if arguments.minimum_interval_median_percentile != 0.0 {
        let threshold = percentile(
            &original_medians,
            arguments.minimum_interval_median_percentile,
        );
        for interval in 0..intervals {
            if original_medians[interval] <= threshold {
                filtered_intervals[interval] = true;
            }
        }
    }

    // The division happens whatever the filters did, and it uses the ORIGINAL medians.
    let unfiltered_samples: Vec<usize> = (0..samples)
        .filter(|sample| !filtered_samples[*sample])
        .collect();
    for interval in 0..intervals {
        if filtered_intervals[interval] {
            continue;
        }
        for sample in &unfiltered_samples {
            let value = counts.get(*sample, interval);
            counts.set(*sample, interval, value / original_medians[interval]);
        }
    }

    // A percentage of a hundred skips the step.
    if arguments.maximum_zeros_in_sample_percentage != 100.0 {
        let passing_intervals = filtered_intervals.iter().filter(|f| !**f).count();
        let candidates: Vec<usize> = (0..samples)
            .filter(|sample| !filtered_samples[*sample])
            .collect();
        for sample in candidates {
            let zeros = (0..intervals)
                .filter(|interval| {
                    !filtered_intervals[*interval] && counts.get(sample, *interval) == 0.0
                })
                .count();
            if zeros as f64 / passing_intervals as f64
                >= arguments.maximum_zeros_in_sample_percentage / 100.0
            {
                filtered_samples[sample] = true;
            }
        }
    }

    if arguments.maximum_zeros_in_interval_percentage != 100.0 {
        let passing_samples = filtered_samples.iter().filter(|f| !**f).count();
        let candidates: Vec<usize> = (0..intervals)
            .filter(|interval| !filtered_intervals[*interval])
            .collect();
        for interval in candidates {
            let zeros = (0..samples)
                .filter(|sample| !filtered_samples[*sample] && counts.get(*sample, interval) == 0.0)
                .count();
            if zeros as f64 / passing_samples as f64
                >= arguments.maximum_zeros_in_interval_percentage / 100.0
            {
                filtered_intervals[interval] = true;
            }
        }
    }

    if arguments.extreme_sample_median_percentile != 0.0 {
        // The medians are taken for EVERY sample, filtered or not, which the reference calls
        // unnecessary bookkeeping; the comparison then applies to every sample too, so a sample
        // already filtered can be filtered again to no effect.
        let sample_medians: Vec<f64> = (0..samples)
            .map(|sample| {
                let kept: Vec<f64> = (0..intervals)
                    .filter(|interval| !filtered_intervals[*interval])
                    .map(|interval| counts.get(sample, interval))
                    .collect();
                median(&kept)
            })
            .collect();
        let minimum = percentile(&sample_medians, arguments.extreme_sample_median_percentile);
        let maximum = percentile(
            &sample_medians,
            100.0 - arguments.extreme_sample_median_percentile,
        );
        // Strictly outside, so a sample sitting exactly on either threshold is kept.
        let extreme: Vec<usize> = sample_medians
            .iter()
            .enumerate()
            .filter(|(_, value)| **value < minimum || **value > maximum)
            .map(|(sample, _)| sample)
            .collect();
        for sample in extreme {
            filtered_samples[sample] = true;
        }
    }

    let panel_intervals: Vec<usize> = (0..intervals)
        .filter(|interval| !filtered_intervals[*interval])
        .collect();
    let panel_samples: Vec<usize> = (0..samples)
        .filter(|sample| !filtered_samples[*sample])
        .collect();
    let mut values = Matrix {
        samples: panel_samples.len(),
        intervals: panel_intervals.len(),
        values: panel_samples
            .iter()
            .flat_map(|sample| {
                panel_intervals
                    .iter()
                    .map(|interval| counts.get(*sample, *interval))
                    .collect::<Vec<_>>()
            })
            .collect(),
    };
    let panel_interval_fractional_medians: Vec<f64> = panel_intervals
        .iter()
        .map(|interval| original_medians[*interval])
        .collect();

    if arguments.impute_zeros {
        // The median a zero becomes is over the NON-ZERO values of its interval alone, so an
        // interval that is all zeros in the panel imputes a NaN.
        let non_zero_medians: Vec<f64> = (0..values.intervals)
            .map(|interval| {
                let non_zero: Vec<f64> = values
                    .column(interval)
                    .into_iter()
                    .filter(|value| *value > 0.0)
                    .collect();
                median(&non_zero)
            })
            .collect();
        for sample in 0..values.samples {
            for (interval, replacement) in non_zero_medians.iter().enumerate() {
                if values.get(sample, interval) == 0.0 {
                    values.set(sample, interval, *replacement);
                }
            }
        }
    }

    if arguments.extreme_outlier_truncation_percentile != 0.0 {
        let minimum = percentile(
            &values.values,
            arguments.extreme_outlier_truncation_percentile,
        );
        let maximum = percentile(
            &values.values,
            100.0 - arguments.extreme_outlier_truncation_percentile,
        );
        for value in values.values.iter_mut() {
            if *value < minimum {
                *value = minimum;
            } else if *value > maximum {
                *value = maximum;
            }
        }
    }

    Preprocessed {
        filtered_samples,
        filtered_intervals,
        panel_interval_fractional_medians,
        values,
    }
}

/// The floor `safeLog2` clamps at, below which the logarithm is not taken at all.
pub const EPSILON: f64 = 1e-9;

/// `INV_LOG_2`, which is what the reference multiplies a natural logarithm by.
const INV_LOG_2: f64 = 1.0 / std::f64::consts::LN_2;

/// `safeLog2`: the base-two logarithm, floored rather than allowed to run to negative infinity.
///
/// A value BELOW the epsilon becomes `log2(epsilon)` outright; the epsilon is not added to it, so
/// a value just above the floor is not nudged and a value at the floor exactly is not floored.
pub fn safe_log2(x: f64) -> f64 {
    if x < EPSILON {
        EPSILON.ln() * INV_LOG_2
    } else {
        x.ln() * INV_LOG_2
    }
}

/// The refusal a sample whose median is not positive produces.
pub fn non_positive_median_message(sample: usize, samples: usize) -> String {
    if samples == 1 {
        "Sample does not have a positive sample median.".to_string()
    } else {
        format!("Sample at index {sample} does not have a positive sample median.")
    }
}

/// `divideBySampleMedianAndTransformToLog2`, then the median of the sample medians subtracted.
///
/// The second step is what separates the PANEL's standardisation from a single sample's: a panel
/// subtracts one median from every row, while a sample subtracts its own row's.
pub fn standardize(values: &mut Matrix) -> Result<(), String> {
    let sample_medians: Vec<f64> = (0..values.samples)
        .map(|sample| median(values.row(sample)))
        .collect();
    for (sample, sample_median) in sample_medians.iter().enumerate() {
        if *sample_median <= 0.0 {
            return Err(non_positive_median_message(sample, values.samples));
        }
    }
    for (sample, sample_median) in sample_medians.iter().enumerate() {
        for interval in 0..values.intervals {
            let value = values.get(sample, interval);
            values.set(sample, interval, safe_log2(value / sample_median));
        }
    }
    let log2_medians: Vec<f64> = (0..values.samples)
        .map(|sample| median(values.row(sample)))
        .collect();
    let median_of_medians = median(&log2_medians);
    for value in values.values.iter_mut() {
        *value -= median_of_medians;
    }
    Ok(())
}

/// The number of eigensamples the panel ends up with.
///
/// It is capped at the number of samples that SURVIVED, not at the number given, so a panel that
/// filtered one sample out of nine has eight however many were asked for. A panel of one sample
/// has none at all.
pub fn number_of_eigensamples(requested: usize, surviving_samples: usize) -> usize {
    requested.min(surviving_samples)
}

/// The refusal a panel with no eigensamples gives when asked for its singular values.
pub const NO_SINGULAR_VALUES_MESSAGE: &str = "No singular values were available.";

/// `HDF5SVDReadCountPanelOfNormals.create`, where the decomposition found nothing to keep.
///
/// A panel of more than one sample must yield at least one singular value over the solver's own
/// epsilon. Filtering hard enough to leave a handful of intervals is what reaches this, and the
/// message suggests the opposite of the filter that caused it.
pub const NO_NON_ZERO_SINGULAR_VALUES_MESSAGE: &str =
    "No non-zero singular values were found.  It may be necessary to use stricter parameters for \
filtering.  For example, use a larger value of minimum-interval-median-percentile.";

/// The refusal an input whose intervals do not match the others' produces.
pub fn mismatched_intervals_message(path: &str) -> String {
    format!("Intervals for read-counts file {path} do not match those in other read-counts files.")
}
