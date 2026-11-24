//! Ported from `org.broadinstitute.hellbender.tools.copynumber.DenoiseReadCounts`,
//! `denoising.SVDDenoisingUtils` and the reading half of `denoising.HDF5SVDReadCountPanelOfNormals`
//! (GATK 4.7.0.0).
//!
//! With no panel of normals and no GC annotations the tool does what its documentation calls
//! standardization, and nothing else: fractional coverage, divide by the sample median, log2,
//! subtract the median again. With a panel ([`denoise`]) the sample is first corrected by the
//! panel's GC content, cut to the panel's intervals and divided by their medians, then standardized
//! the same way, and the projection onto the panel's first eigensamples is subtracted.
//!
//! # Scale-free after the first step
//!
//! `transformToFractionalCoverage` divides every count by the row's sum, so doubling every count
//! changes nothing downstream. The golden has both, and they are the same file.
//!
//! # The floor, not an infinity
//!
//! `safeLog2(x)` is `x < 1e-9 ? log2(1e-9) : log(x) * INV_LOG_2`. A zero count therefore reads
//! **-29.897353** rather than minus infinity, and -- more importantly -- that floored value takes
//! part in the median like any other number. A port emitting an infinity would differ on the value
//! and on every median beside it.
//!
//! The logarithm is `Math.log(x) * INV_LOG_2` rather than a base-two logarithm, which is a
//! different last bit and is written that way here.
//!
//! # The median twice, and not the same median
//!
//! Once as a division BEFORE the log, once as a subtraction AFTER it. The second is taken over the
//! log values, so it is not the logarithm of the first.
//!
//! # `isPositive` is not `> 0` written the other way round
//!
//! An all-zero row divides by a sum of zero, so every value is a NaN and so is the median.
//! `ParamUtils.isPositive` refuses it, because `NaN > 0` is false. A port testing `median <= 0`
//! would NOT refuse it, because that comparison is false for a NaN as well, and would write NaNs
//! where the reference writes nothing. The test here is the reference's way round for that reason,
//! and the golden's `all-zero` row is what says which way that is.

/// `SVDDenoisingUtils.EPSILON`.
pub const EPSILON: f64 = 1e-9;

/// `MathUtils.INV_LOG_2`, which is how the reference divides by `ln 2`.
const INV_LOG_2: f64 = std::f64::consts::LOG2_E;

/// `SVDDenoisingUtils.LN2_EPSILON`: the value every count below [`EPSILON`] collapses to.
pub fn ln2_epsilon() -> f64 {
    EPSILON.ln() * INV_LOG_2
}

/// `safeLog2`.
pub fn safe_log2(value: f64) -> f64 {
    if value < EPSILON {
        ln2_epsilon()
    } else {
        value.ln() * INV_LOG_2
    }
}

/// What the standardization refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DenoiseError {
    /// `ParamUtils.isPositive` on the row median. The message names an index only when there is
    /// more than one sample, and on this path there never is.
    NonPositiveSampleMedian,
}

impl DenoiseError {
    /// The exception class the reference throws.
    pub fn java_class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }

    /// The message for a single sample, which is the only shape this path produces.
    pub fn message(&self) -> &'static str {
        "Sample does not have a positive sample median."
    }
}

/// `transformToFractionalCoverage`: divide by the row's sum.
///
/// The sum is `MathUtils.sum`, a plain loop, and not a stream's compensated one.
pub fn transform_to_fractional_coverage(counts: &[f64]) -> Vec<f64> {
    let mut total = 0.0;
    for count in counts {
        total += count;
    }
    counts.iter().map(|count| count / total).collect()
}

/// `new Median().evaluate(row)`, which is commons-math's and therefore interpolates.
pub fn median(values: &[f64]) -> f64 {
    jmath::percentile::median(values)
}

/// `preprocessAndStandardizeSample` with no GC annotations.
pub fn standardize(counts: &[f64]) -> Result<Vec<f64>, DenoiseError> {
    standardize_preprocessed(&transform_to_fractional_coverage(counts))
}

/// `divideBySampleMedianAndTransformToLog2`, then the log median subtracted: the standardization
/// both `preprocessAndStandardizeSample`s end with, over values already preprocessed.
fn standardize_preprocessed(fractional: &[f64]) -> Result<Vec<f64>, DenoiseError> {
    // `ParamUtils.isPositive`, whose test is `x > 0` -- so a NaN median is refused. Written as a
    // negated `>` rather than as `<= 0`, which would let the NaN through: the golden's `all-zero`
    // row is a refusal, and `NaN <= 0` is false.
    let sample_median = median(fractional);
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(sample_median > 0.0) {
        return Err(DenoiseError::NonPositiveSampleMedian);
    }
    let logs: Vec<f64> = fractional
        .iter()
        .map(|value| safe_log2(value / sample_median))
        .collect();

    // The second median, over the LOG values.
    let log_median = median(&logs);
    Ok(logs.iter().map(|value| value - log_median).collect())
}

/// `HDF5SVDReadCountPanelOfNormals.CURRENT_PON_VERSION`.
pub const CURRENT_PON_VERSION: f64 = 7.0;

/// What [`denoise`] reads of a panel of normals.
#[derive(Debug, Clone, PartialEq)]
pub struct Panel {
    /// `/original_data/intervals`: the intervals the panel's samples were counted over.
    pub original_intervals: Vec<Interval>,
    /// `/original_data/interval_gc_content`, when the panel was built with GC annotations.
    pub gc_content: Option<Vec<f64>>,
    /// `/panel/intervals`: the original intervals the panel's filters kept, in their order.
    pub panel_intervals: Vec<Interval>,
    /// `/panel/interval_fractional_medians`, one per panel interval.
    pub fractional_medians: Vec<f64>,
    /// `getEigensampleVectors`: panel intervals by eigensamples, the transpose of what is stored.
    /// Empty when the panel stores none (a single sample, or zero eigensamples requested).
    pub eigensamples: Vec<Vec<f64>>,
}

pub use crate::create_read_count_panel_of_normals::Interval;

/// What reading a panel refuses, as `HDF5SVDReadCountPanelOfNormals.read` refuses it.
#[derive(Debug, Clone, PartialEq)]
pub enum PanelError {
    /// `getVersion`, when `/version/value` is not there.
    OutOfDate,
    /// `read`, when the version is below [`CURRENT_PON_VERSION`].
    Older(f64),
    /// A dataset the panel must hold is missing; the reference fails in the HDF5 library there.
    Missing(&'static str),
}

impl PanelError {
    /// The message of the `UserException.BadInput`, without its `Bad input: ` prefix, or the
    /// path that was missing.
    pub fn message(&self) -> String {
        match self {
            PanelError::OutOfDate => format!(
                "The panel of normals is out of date and incompatible.  Please use a panel of \
                 normals that was created by CreateReadCountPanelOfNormals and is version {:.1}.",
                CURRENT_PON_VERSION
            ),
            // `%f`, which is six decimals.
            PanelError::Older(version) => format!(
                "The version of the specified panel of normals ({version:.6}) is older than the \
                 current version ({CURRENT_PON_VERSION:.6})."
            ),
            PanelError::Missing(path) => path.to_string(),
        }
    }
}

/// `HDF5SVDReadCountPanelOfNormals.read` and the getters `denoise` calls.
///
/// `getNumEigensamples` reads `num_rows` under the eigensample matrix and takes a failure to
/// read it as zero, which is how a panel with no decomposition answers; the vectors are then
/// never read.
pub fn read_panel(reader: &crate::hdf5_reader::Hdf5Reader) -> Result<Panel, PanelError> {
    use crate::create_read_count_panel_of_normals::{read_chunked_matrix, read_intervals};
    let (_, version) = reader
        .doubles("/version/value")
        .ok_or(PanelError::OutOfDate)?;
    let version = version.first().copied().ok_or(PanelError::OutOfDate)?;
    if version < CURRENT_PON_VERSION {
        return Err(PanelError::Older(version));
    }
    let vectors = "/panel/transposed_eigensamples_samples_by_intervals";
    let stored = reader
        .doubles(&format!("{vectors}/num_rows"))
        .and_then(|(_, values)| values.first().map(|value| *value as i32))
        .unwrap_or(0);
    let eigensamples = if stored == 0 {
        Vec::new()
    } else {
        let transposed =
            read_chunked_matrix(reader, vectors).ok_or(PanelError::Missing(vectors))?;
        let intervals = transposed.first().map_or(0, Vec::len);
        (0..intervals)
            .map(|interval| transposed.iter().map(|row| row[interval]).collect())
            .collect()
    };
    Ok(Panel {
        original_intervals: read_intervals(reader, "/original_data/intervals")
            .ok_or(PanelError::Missing("/original_data/intervals"))?,
        gc_content: reader
            .doubles("/original_data/interval_gc_content")
            .map(|(_, values)| values),
        panel_intervals: read_intervals(reader, "/panel/intervals")
            .ok_or(PanelError::Missing("/panel/intervals"))?,
        fractional_medians: reader
            .doubles("/panel/interval_fractional_medians")
            .map(|(_, values)| values)
            .ok_or(PanelError::Missing("/panel/interval_fractional_medians"))?,
        eigensamples,
    })
}

impl Panel {
    /// `getNumEigensamples`.
    pub fn number_of_eigensamples(&self) -> usize {
        self.eigensamples.first().map_or(0, Vec::len)
    }
}

/// `Utils.validateArg`'s message when the sample was counted over other intervals than the panel.
pub const MISMATCHED_INTERVALS: &str =
    "Sample intervals must be identical to the original intervals used to build the panel of normals.";

/// `SVDDenoisingUtils.denoise`: the standardized and the denoised copy ratios, over the panel's
/// intervals. `Err` is the message of the `IllegalArgumentException` the reference throws.
///
/// `eigensamples` is already `min(requested, available)`, as `doWork` computes it.
pub fn denoise(
    panel: &Panel,
    intervals: &[Interval],
    counts: &[f64],
    eigensamples: usize,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    if panel.original_intervals != intervals {
        return Err(MISMATCHED_INTERVALS.to_string());
    }
    // `preprocessAndStandardizeSample(panelOfNormals, ...)`: fractional coverage, the panel's GC
    // correction, the subset to the panel's intervals, and the division by their medians.
    let mut fractional = transform_to_fractional_coverage(counts);
    if let Some(gc_content) = &panel.gc_content {
        let mut matrix = crate::create_read_count_panel_of_normals::Matrix::new(&[fractional]);
        crate::create_read_count_panel_of_normals::correct_gc_bias(&mut matrix, gc_content);
        fractional = matrix.row(0).to_vec();
    }
    let kept: std::collections::HashSet<&Interval> = panel.panel_intervals.iter().collect();
    let subset: Vec<f64> = panel
        .original_intervals
        .iter()
        .zip(&fractional)
        .filter(|(interval, _)| kept.contains(interval))
        .zip(&panel.fractional_medians)
        .map(|((_, value), median)| value / median)
        .collect();
    let standardized =
        standardize_preprocessed(&subset).map_err(|error| error.message().to_string())?;
    if eigensamples == 0 || panel.number_of_eigensamples() == 0 {
        return Ok((standardized.clone(), standardized));
    }
    // `subtractProjection`: s - (s U_k) U_k^T, each product commons-math's
    // `Array2DRowRealMatrix.multiply`, a sum from zero over the shared index in order.
    let projected: Vec<f64> = (0..eigensamples)
        .map(|k| {
            let mut sum = 0.0;
            for (value, row) in standardized.iter().zip(&panel.eigensamples) {
                sum += value * row[k];
            }
            sum
        })
        .collect();
    let denoised = standardized
        .iter()
        .zip(&panel.eigensamples)
        .map(|(value, row)| {
            let mut sum = 0.0;
            for (weight, vector) in projected.iter().zip(row) {
                sum += weight * vector;
            }
            value - sum
        })
        .collect();
    Ok((standardized, denoised))
}

/// The output table's header, which is the input's with one column renamed.
pub fn header(sequences: &[htsjdk_bam::header::SequenceRecord], sample: &str) -> String {
    let mut text = crate::model_segments_tool::header(sample, Some(sequences));
    text.push_str("CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n");
    text
}

/// The whole file: the header, then one row per interval, formatted to six decimals.
pub fn write(
    sequences: &[htsjdk_bam::header::SequenceRecord],
    sample: &str,
    intervals: &[(String, i32, i32)],
    values: &[f64],
) -> String {
    let mut text = header(sequences, sample);
    for ((contig, start, end), value) in intervals.iter().zip(values) {
        text.push_str(&format!(
            "{contig}\t{start}\t{end}\t{}\n",
            gatk_engine::java_format::format_decimals(*value, 6)
        ));
    }
    text
}
