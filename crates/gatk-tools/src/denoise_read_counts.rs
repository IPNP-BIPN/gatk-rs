//! Ported from `org.broadinstitute.hellbender.tools.copynumber.DenoiseReadCounts` and the
//! no-panel half of `denoising.SVDDenoisingUtils` (GATK 4.6.2.0).
//!
//! With no panel of normals and no GC annotations the tool does what its documentation calls
//! standardization, and nothing else: fractional coverage, divide by the sample median, log2,
//! subtract the median again.
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
    let fractional = transform_to_fractional_coverage(counts);

    // `ParamUtils.isPositive`, whose test is `x > 0` -- so a NaN median is refused. Written as a
    // negated `>` rather than as `<= 0`, which would let the NaN through: the golden's `all-zero`
    // row is a refusal, and `NaN <= 0` is false.
    let sample_median = median(&fractional);
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

/// The output table's header, which is the input's with one column renamed.
pub fn header(sequences: &[(String, i32)], sample: &str) -> String {
    let mut text = String::from("@HD\tVN:1.6\n");
    for (name, length) in sequences {
        text.push_str(&format!("@SQ\tSN:{name}\tLN:{length}\n"));
    }
    text.push_str(&format!("@RG\tID:GATKCopyNumber\tSM:{sample}\n"));
    text.push_str("CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\n");
    text
}

/// The whole file: the header, then one row per interval, formatted to six decimals.
pub fn write(
    sequences: &[(String, i32)],
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
