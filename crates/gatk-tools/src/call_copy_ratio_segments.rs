//! Ported from `org.broadinstitute.hellbender.tools.copynumber.CallCopyRatioSegments` and
//! `caller.SimpleCopyRatioCaller` (GATK 4.6.2.0).
//!
//! Copy-ratio segments in, called segments out. No reads, no reference, no walker: this is a table
//! transformed by arithmetic, and every decision in it is a comparison against a statistic computed
//! from the table itself.
//!
//! # The statistics are computed twice, and the two passes are not interchangeable
//!
//! `calculateCallingStatistics` takes the copy-neutral segments, computes a length-weighted mean
//! and standard deviation over ALL of them, uses THOSE to drop the outliers, and recomputes over
//! what is left. The calling threshold is compared against the second pair; the outlier threshold
//! against the first. A port that filtered with the recomputed statistics would drop a different
//! set of segments.
//!
//! # The standard deviation's denominator is its own
//!
//! `sum / (((n - 1) / n) * totalLength)`, with `n` the segment COUNT and `totalLength` the sum of
//! their lengths. It is neither the population form nor the sample form, and at `n = 1` it is zero:
//! the answer is then an infinity or a NaN, and every comparison against it is false, so a run with
//! one copy-neutral segment calls everything NEUTRAL rather than failing.
//!
//! # Everything is in copy-ratio space
//!
//! The bounds, the mean, the deviation and the thresholds are all `2^log2`, never the log2 value.
//! Only the OUTPUT keeps the log2, formatted to six decimals.

/// `CalledCopyRatioSegment.Call`, whose printed forms are one character each.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Call {
    /// `-`.
    Deletion,
    /// `0`.
    Neutral,
    /// `+`.
    Amplification,
}

impl Call {
    /// The character the table carries.
    pub fn as_str(&self) -> &'static str {
        match self {
            Call::Deletion => "-",
            Call::Neutral => "0",
            Call::Amplification => "+",
        }
    }
}

/// One input row.
#[derive(Debug, Clone, PartialEq)]
pub struct CopyRatioSegment {
    /// The contig.
    pub contig: String,
    /// One-based inclusive start.
    pub start: i32,
    /// One-based inclusive end.
    pub end: i32,
    /// `NUM_POINTS_COPY_RATIO`, which is carried through and decides nothing.
    pub num_points: i32,
    /// `MEAN_LOG2_COPY_RATIO`.
    pub mean_log2_copy_ratio: f64,
}

impl CopyRatioSegment {
    /// `getLengthOnReference()`, which is inclusive at both ends.
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }

    /// `Math.pow(2, meanLog2CopyRatio)`, which is the space every comparison happens in.
    pub fn copy_ratio(&self) -> f64 {
        2.0f64.powf(self.mean_log2_copy_ratio)
    }
}

/// `SimpleCopyRatioCaller.DEFAULT` bounds and thresholds.
pub const DEFAULT_NEUTRAL_LOWER_BOUND: f64 = 0.9;
/// The upper bound of the copy-neutral band.
pub const DEFAULT_NEUTRAL_UPPER_BOUND: f64 = 1.1;
/// The z-score above which a copy-neutral segment is an outlier.
pub const DEFAULT_OUTLIER_Z_SCORE: f64 = 2.0;
/// The z-score above which a segment is called.
pub const DEFAULT_CALLING_Z_SCORE: f64 = 2.0;

/// What the constructor refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CallerError {
    /// A negative lower bound.
    NegativeLowerBound,
    /// A lower bound at or above the upper one.
    BoundsNotOrdered,
    /// An outlier threshold at or below zero.
    OutlierThresholdNotPositive,
    /// A calling threshold at or below zero.
    CallingThresholdNotPositive,
}

impl CallerError {
    /// The exception class, which is the same for all four.
    pub fn java_class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }

    /// The message, which for the bounds names both arguments.
    pub fn message(&self) -> &'static str {
        match self {
            CallerError::NegativeLowerBound => "Copy-neutral lower bound must be non-negative.",
            CallerError::BoundsNotOrdered => {
                "Copy-neutral lower bound (neutral-segment-copy-ratio-lower-bound) must be less \
                 than upper bound (neutral-segment-copy-ratio-upper-bound)."
            }
            CallerError::OutlierThresholdNotPositive => {
                "Outlier z-score threshold must be positive."
            }
            CallerError::CallingThresholdNotPositive => {
                "Calling z-score threshold must be positive."
            }
        }
    }
}

/// A length-weighted mean and standard deviation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Statistics {
    /// The length-weighted mean, in copy-ratio space.
    pub mean: f64,
    /// The length-weighted standard deviation, with the reference's own denominator.
    pub standard_deviation: f64,
}

/// `calculateLengthWeightedStatistics`.
///
/// Every sum here is a `DoubleStream.sum`, which is Kahan-compensated, so it is not a plain loop.
/// The denominator is `((n - 1) / n) * totalLength` with the division done in `double`, which for
/// `n = 1` is exactly zero.
pub fn length_weighted_statistics(segments: &[&CopyRatioSegment]) -> Statistics {
    let lengths: Vec<f64> = segments
        .iter()
        .map(|segment| f64::from(segment.length()))
        .collect();
    let total_length = gatk_engine::allele_fraction_cluster::double_stream_sum(&lengths);
    let count = segments.len();

    let weighted: Vec<f64> = segments
        .iter()
        .map(|segment| f64::from(segment.length()) * segment.copy_ratio())
        .collect();
    let mean = gatk_engine::allele_fraction_cluster::double_stream_sum(&weighted) / total_length;

    let squares: Vec<f64> = segments
        .iter()
        .map(|segment| f64::from(segment.length()) * (segment.copy_ratio() - mean).powi(2))
        .collect();
    let standard_deviation = (gatk_engine::allele_fraction_cluster::double_stream_sum(&squares)
        / (((count as f64 - 1.0) / count as f64) * total_length))
        .sqrt();

    Statistics {
        mean,
        standard_deviation,
    }
}

/// `calculateCallingStatistics`: the two passes, and the filter between them.
pub fn calling_statistics(
    segments: &[CopyRatioSegment],
    lower_bound: f64,
    upper_bound: f64,
    outlier_threshold: f64,
) -> Statistics {
    let neutral: Vec<&CopyRatioSegment> = segments
        .iter()
        .filter(|segment| is_neutral(segment, lower_bound, upper_bound))
        .collect();
    let unfiltered = length_weighted_statistics(&neutral);
    // The outlier test is INCLUSIVE, and it uses the unfiltered statistics rather than any
    // recomputed ones.
    let filtered: Vec<&CopyRatioSegment> = neutral
        .into_iter()
        .filter(|segment| {
            (segment.copy_ratio() - unfiltered.mean).abs()
                <= unfiltered.standard_deviation * outlier_threshold
        })
        .collect();
    length_weighted_statistics(&filtered)
}

/// The copy-neutral test, inclusive at both ends and in copy-ratio space.
fn is_neutral(segment: &CopyRatioSegment, lower_bound: f64, upper_bound: f64) -> bool {
    let copy_ratio = segment.copy_ratio();
    lower_bound <= copy_ratio && copy_ratio <= upper_bound
}

/// `makeCalls`: one call per segment.
///
/// A segment inside the band is neutral without any statistic being consulted. Outside it, the
/// deviation from the filtered mean decides, and a NaN mean makes both comparisons false -- which
/// is why an empty or single-segment copy-neutral set calls everything neutral.
pub fn make_calls(
    segments: &[CopyRatioSegment],
    lower_bound: f64,
    upper_bound: f64,
    outlier_threshold: f64,
    calling_threshold: f64,
) -> Result<Vec<Call>, CallerError> {
    if lower_bound < 0.0 {
        return Err(CallerError::NegativeLowerBound);
    }
    if lower_bound >= upper_bound {
        return Err(CallerError::BoundsNotOrdered);
    }
    if outlier_threshold <= 0.0 {
        return Err(CallerError::OutlierThresholdNotPositive);
    }
    if calling_threshold <= 0.0 {
        return Err(CallerError::CallingThresholdNotPositive);
    }

    let statistics = calling_statistics(segments, lower_bound, upper_bound, outlier_threshold);
    Ok(segments
        .iter()
        .map(|segment| {
            if is_neutral(segment, lower_bound, upper_bound) {
                return Call::Neutral;
            }
            let deviation = segment.copy_ratio() - statistics.mean;
            if deviation < -statistics.standard_deviation * calling_threshold {
                Call::Deletion
            } else if deviation > statistics.standard_deviation * calling_threshold {
                Call::Amplification
            } else {
                Call::Neutral
            }
        })
        .collect())
}

/// `CopyNumberFormatsUtils.DOUBLE_FORMAT`, which is `%.6f`.
pub fn format_double(value: f64) -> String {
    gatk_engine::java_format::format_decimals(value, 6)
}

/// The called-segments table: the input's header, its columns plus `CALL`, then the rows.
pub fn write_called(header: &str, segments: &[CopyRatioSegment], calls: &[Call]) -> String {
    let mut text = String::from(header);
    text.push_str("CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\tCALL\n");
    for (segment, call) in segments.iter().zip(calls) {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            segment.contig,
            segment.start,
            segment.end,
            segment.num_points,
            format_double(segment.mean_log2_copy_ratio),
            call.as_str()
        ));
    }
    text
}

/// The legacy `.igv.seg` file, whose columns are IGV's and whose order is not the table's.
///
/// The call comes BEFORE the mean here, and the sample name is a column rather than a header line.
pub fn write_legacy(sample: &str, segments: &[CopyRatioSegment], calls: &[Call]) -> String {
    let mut text = String::from("Sample\tChromosome\tStart\tEnd\tNum_Probes\tCall\tSegment_Mean\n");
    for (segment, call) in segments.iter().zip(calls) {
        text.push_str(&format!(
            "{sample}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            segment.contig,
            segment.start,
            segment.end,
            segment.num_points,
            call.as_str(),
            format_double(segment.mean_log2_copy_ratio)
        ));
    }
    text
}
