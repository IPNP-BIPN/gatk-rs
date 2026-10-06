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

/// `GCBiasCorrector`'s bins: 0%, 1%, ... 100%.
const NUMBER_OF_GC_BINS: usize = 101;

/// `GCBiasCorrector.correlationDecayRatePerBin`, from a correlation length of 0.02.
const CORRELATION_DECAY_RATE_PER_BIN: f64 = 1.0 / (0.02 * NUMBER_OF_GC_BINS as f64);

/// `gcContentToBinIndex`: `(int) Math.round(gc * 100)`. GC content is a fraction, so the halves
/// `Math.round` and `f64::round` part on (negative ones) never reach it, and a NaN lands in bin 0
/// on both sides.
fn gc_bin(gc_content: f64) -> usize {
    (gc_content * (NUMBER_OF_GC_BINS - 1) as f64).round() as usize
}

/// `GCBiasCorrector.correctGCBias`, in place, over every sample of the fractional coverage.
///
/// Each sample gets its own curve: the median coverage in each GC bin, smoothed over the bins
/// by weights that decay with distance and count each bin's intervals, so an empty bin weighs
/// nothing (its median, a placeholder 1.0, is never read). A value is multiplied by its bin's
/// inverse, then the sample is rescaled so its total coverage is what it was.
pub fn correct_gc_bias(counts: &mut Matrix, gc_content: &[f64]) {
    let bins: Vec<usize> = gc_content.iter().map(|gc| gc_bin(*gc)).collect();
    let totals: Vec<f64> = (0..counts.samples)
        .map(|sample| counts.row(sample).iter().map(|v| v.abs()).sum())
        .collect();
    let factors: Vec<Vec<f64>> = (0..counts.samples)
        .map(|sample| {
            let mut by_bin: Vec<Vec<f64>> = vec![Vec::new(); NUMBER_OF_GC_BINS];
            for (interval, value) in counts.row(sample).iter().enumerate() {
                by_bin[bins[interval]].push(*value);
            }
            let medians: Vec<f64> = by_bin
                .iter()
                .map(|bin| if bin.is_empty() { 1.0 } else { median(bin) })
                .collect();
            (0..NUMBER_OF_GC_BINS)
                .map(|bin| {
                    let weights: Vec<f64> = (0..NUMBER_OF_GC_BINS)
                        .map(|n| {
                            let distance = (bin as f64 - n as f64).abs();
                            by_bin[n].len() as f64
                                * (-distance * CORRELATION_DECAY_RATE_PER_BIN).exp()
                        })
                        .collect();
                    let mut dot = 0.0;
                    for (weight, median) in weights.iter().zip(&medians) {
                        dot += weight * median;
                    }
                    let norm: f64 = weights.iter().map(|w| w.abs()).sum();
                    1.0 / (dot / norm)
                })
                .collect()
        })
        .collect();
    for (sample, sample_factors) in factors.iter().enumerate() {
        for interval in 0..counts.intervals {
            let value = counts.get(sample, interval);
            counts.set(sample, interval, sample_factors[bins[interval]] * value);
        }
    }
    for (sample, total) in totals.iter().enumerate() {
        let corrected: f64 = counts.row(sample).iter().map(|v| v.abs()).sum();
        let normalization = total / corrected;
        for interval in 0..counts.intervals {
            let value = counts.get(sample, interval);
            counts.set(sample, interval, value * normalization);
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
    preprocess_with_gc(counts, None, arguments).unwrap_or_else(|message| panic!("{message}"))
}

/// `countNumberPassingFilter`'s refusal, a `UserException.BadInput`.
pub const FILTERED_EVERYTHING_MESSAGE: &str =
    "Bad input: Filtering removed all samples or intervals.  Select less strict filtering criteria.";

/// `countNumberPassingFilter`: how many survive, or the refusal when none do. The reference counts
/// at every step that runs (in the log line after it, and at the head of the two zero filters),
/// so the step that empties the panel is where the run stops.
fn passing(filter: &[bool]) -> Result<usize, String> {
    match filter.iter().filter(|f| !**f).count() {
        0 => Err(FILTERED_EVERYTHING_MESSAGE.to_string()),
        n => Ok(n),
    }
}

/// `preprocessPanel` with `--annotated-intervals`: the GC-bias correction runs on the fractional
/// coverage, before the first filter takes its medians.
pub fn preprocess_with_gc(
    counts: &Matrix,
    gc_content: Option<&[f64]>,
    arguments: &Arguments,
) -> Result<Preprocessed, String> {
    let mut counts = counts.clone();
    to_fractional_coverage(&mut counts);
    if let Some(gc_content) = gc_content {
        correct_gc_bias(&mut counts, gc_content);
    }
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
        passing(&filtered_intervals)?;
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
        let passing_intervals = passing(&filtered_intervals)?;
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
        passing(&filtered_samples)?;
    }

    if arguments.maximum_zeros_in_interval_percentage != 100.0 {
        let passing_samples = passing(&filtered_samples)?;
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
        passing(&filtered_intervals)?;
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
        passing(&filtered_samples)?;
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

    Ok(Preprocessed {
        filtered_samples,
        filtered_intervals,
        panel_interval_fractional_medians,
        values,
    })
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

/// `HDF5SVDReadCountPanelOfNormals.EPSILON`, also the reciprocal condition number handed to Spark.
pub const SVD_EPSILON: f64 = 1e-9;

/// The truncated decomposition of the standardised panel: singular values, decreasing, and the
/// eigensample vectors, one row per panel interval and one column per singular value.
#[derive(Debug, Clone, PartialEq)]
pub struct Decomposition {
    pub singular_values: Vec<f64>,
    pub eigensample_vectors: Vec<Vec<f64>>,
}

/// `RowMatrix.computeSVD(k, true, 1e-9)` over the panel transposed to intervals by samples.
///
/// A panel has fewer than a hundred samples in every case this port meets, which is where Spark
/// takes its local path: the samples-by-samples Gramian, its eigendecomposition, the square roots
/// as singular values, kept while at least `rCond` times the first, and U as A V S^-1.
///
/// Which of those values survive is not a property of the input. The Gramian of a standardised
/// panel is rank deficient, its trailing eigenvalues are rounding noise, and whether their roots
/// clear the threshold moves with the order a cluster sums in: the reference itself answered 8,
/// 6 and 7 on one fixture (docs/a-rank-is-not-a-byte.md). This solver is cyclic Jacobi, which is
/// deterministic; the coverage harness compares the bound and the branch, not these values.
pub fn truncated_svd(panel: &Matrix, k: usize) -> Decomposition {
    let n = panel.samples;
    let mut gramian = vec![vec![0.0; n]; n];
    for (i, row) in gramian.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = panel
                .row(i)
                .iter()
                .zip(panel.row(j))
                .map(|(a, b)| a * b)
                .sum();
        }
    }
    let (eigenvalues, eigenvectors) = symmetric_eigen(gramian);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|a, b| eigenvalues[*b].total_cmp(&eigenvalues[*a]));
    let sigmas: Vec<f64> = order
        .iter()
        .map(|i| eigenvalues[*i].max(0.0).sqrt())
        .collect();
    let threshold = SVD_EPSILON * sigmas.first().copied().unwrap_or(0.0);
    let kept = sigmas
        .iter()
        .take(k)
        .take_while(|sigma| **sigma >= threshold)
        .count();
    let singular_values = sigmas[..kept].to_vec();
    let eigensample_vectors = (0..panel.intervals)
        .map(|interval| {
            (0..kept)
                .map(|column| {
                    let vector = order[column];
                    let projection: f64 = (0..n)
                        .map(|sample| panel.get(sample, interval) * eigenvectors[sample][vector])
                        .sum();
                    projection / singular_values[column]
                })
                .collect()
        })
        .collect();
    Decomposition {
        singular_values,
        eigensample_vectors,
    }
}

/// Cyclic Jacobi rotations until the off-diagonal is negligible: eigenvalues, and the
/// eigenvectors as COLUMNS of the returned matrix.
fn symmetric_eigen(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for _sweep in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |j| *j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum();
        if off <= f64::EPSILON * f64::EPSILON {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for row in a.iter_mut() {
                    let (akp, akq) = (row[p], row[q]);
                    row[p] = c * akp - s * akq;
                    row[q] = s * akp + c * akq;
                }
                let (low, high) = a.split_at_mut(q);
                for (apk, aqk) in low[p].iter_mut().zip(high[0].iter_mut()) {
                    let (x, y) = (*apk, *aqk);
                    *apk = c * x - s * y;
                    *aqk = s * x + c * y;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ((0..n).map(|i| a[i][i]).collect(), v)
}

/// `HDF5SVDReadCountPanelOfNormals.CURRENT_PON_VERSION`.
pub const PON_VERSION: f64 = 7.0;

/// `HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX`, `Integer.MAX_VALUE / Byte.SIZE`.
pub const MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX: i64 = i32::MAX as i64 / 8;

/// The tool's `--maximum-chunk-size` default, the limit over 16.
pub const DEFAULT_MAXIMUM_CHUNK_SIZE: i64 = MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / 16;

/// One genomic interval, one-based and closed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

/// `HDF5Utils.writeChunkedDoubleMatrix`: the row and column counts, the number of chunks, and the
/// rows in chunks of as many whole rows as `maximum_chunk_size` values hold, the last one partial.
/// A matrix of fewer rows than one chunk holds is a single partial chunk.
pub fn write_chunked_matrix(
    file: &mut crate::hdf5_writer::Hdf5File,
    path: &str,
    rows: &[Vec<f64>],
    maximum_chunk_size: i64,
) {
    let row_count = rows.len();
    let columns = rows.first().map_or(0, Vec::len);
    let per_chunk = (maximum_chunk_size / columns as i64) as usize;
    let filled = row_count.checked_div(per_chunk).unwrap_or(0);
    let partial = filled == 0 || !row_count.is_multiple_of(per_chunk);
    file.double_array(&format!("{path}/num_rows"), &[row_count as f64]);
    file.double_array(&format!("{path}/num_columns"), &[columns as f64]);
    file.double_array(
        &format!("{path}/num_chunks"),
        &[(if partial { filled + 1 } else { filled }) as f64],
    );
    for chunk in 0..filled {
        file.double_matrix(
            &format!("{path}/chunk_{chunk}"),
            &rows[chunk * per_chunk..(chunk + 1) * per_chunk],
        );
    }
    if partial {
        file.double_matrix(
            &format!("{path}/chunk_{filled}"),
            &rows[filled * per_chunk..],
        );
    }
}

/// `HDF5Utils.writeIntervals`: the contigs in order of first appearance, and a three-row matrix
/// of contig index, start and end.
pub fn write_intervals(
    file: &mut crate::hdf5_writer::Hdf5File,
    path: &str,
    intervals: &[Interval],
) {
    let mut contigs: Vec<String> = Vec::new();
    let mut matrix: Vec<Vec<f64>> = (0..3)
        .map(|_| Vec::with_capacity(intervals.len()))
        .collect();
    for interval in intervals {
        let index = match contigs.iter().position(|c| *c == interval.contig) {
            Some(index) => index,
            None => {
                contigs.push(interval.contig.clone());
                contigs.len() - 1
            }
        };
        matrix[0].push(index as f64);
        matrix[1].push(interval.start as f64);
        matrix[2].push(interval.end as f64);
    }
    file.double_matrix(&format!("{path}/transposed_index_start_end"), &matrix);
    file.string_array(&format!("{path}/indexed_contig_names"), &contigs);
}

/// Everything `HDF5SVDReadCountPanelOfNormals.create` writes, in its order.
pub struct Panel<'a> {
    pub command_line: &'a str,
    /// `SAMTextHeaderCodec`'s text of a header holding only the dictionary.
    pub sequence_dictionary: &'a str,
    /// Samples by intervals, as read.
    pub original_counts: &'a Matrix,
    pub sample_filenames: &'a [String],
    pub intervals: &'a [Interval],
    pub gc_content: Option<&'a [f64]>,
    pub preprocessed: &'a Preprocessed,
    pub decomposition: Option<&'a Decomposition>,
    pub maximum_chunk_size: i64,
}

/// The panel as an HDF5 tree. The singular values and the eigensample vectors are written only
/// when a decomposition ran, which is more than one panel sample and more than zero eigensamples.
pub fn write_panel(panel: &Panel) -> crate::hdf5_writer::Hdf5File {
    let mut file = crate::hdf5_writer::Hdf5File::new();
    file.double_array("/version/value", &[PON_VERSION]);
    file.string_array("/command_line/value", &[panel.command_line.to_string()]);
    file.string_array(
        "/sequence_dictionary/value",
        &[panel.sequence_dictionary.to_string()],
    );
    let original: Vec<Vec<f64>> = (0..panel.original_counts.samples)
        .map(|sample| panel.original_counts.row(sample).to_vec())
        .collect();
    write_chunked_matrix(
        &mut file,
        "/original_data/read_counts_samples_by_intervals",
        &original,
        panel.maximum_chunk_size,
    );
    file.string_array("/original_data/sample_filenames", panel.sample_filenames);
    write_intervals(&mut file, "/original_data/intervals", panel.intervals);
    if let Some(gc_content) = panel.gc_content {
        file.double_array("/original_data/interval_gc_content", gc_content);
    }
    let panel_samples: Vec<String> = panel
        .preprocessed
        .panel_samples()
        .into_iter()
        .map(|sample| panel.sample_filenames[sample].clone())
        .collect();
    file.string_array("/panel/sample_filenames", &panel_samples);
    let panel_intervals: Vec<Interval> = panel
        .preprocessed
        .panel_intervals()
        .into_iter()
        .map(|interval| panel.intervals[interval].clone())
        .collect();
    write_intervals(&mut file, "/panel/intervals", &panel_intervals);
    file.double_array(
        "/panel/interval_fractional_medians",
        &panel.preprocessed.panel_interval_fractional_medians,
    );
    if let Some(decomposition) = panel.decomposition {
        file.double_array("/panel/singular_values", &decomposition.singular_values);
        let transposed: Vec<Vec<f64>> = (0..decomposition.singular_values.len())
            .map(|column| {
                decomposition
                    .eigensample_vectors
                    .iter()
                    .map(|row| row[column])
                    .collect()
            })
            .collect();
        write_chunked_matrix(
            &mut file,
            "/panel/transposed_eigensamples_samples_by_intervals",
            &transposed,
            panel.maximum_chunk_size,
        );
    }
    file
}
