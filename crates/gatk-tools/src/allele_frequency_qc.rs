//! `AlleleFrequencyQC`: the chi-squared statistic over the allele-frequency bins, and the p-value
//! it gives.
//!
//! The traversal is `VariantEval`'s with every knob preset, and is not ported here. What is ported
//! is what the tool does with the report that traversal writes: the grouping, the statistic, the
//! distribution it is read against, and the two metrics it writes.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.varianteval.AlleleFrequencyQC` and
//! `org.apache.commons.math3.distribution.ChiSquaredDistribution`.

/// `threshold`, under which the tool logs an error and changes nothing else.
pub const DEFAULT_THRESHOLD: f64 = 0.05;
/// `allowedVariance`, which stands in for the expected count Pearson's statistic would divide by.
pub const DEFAULT_ALLOWED_VARIANCE: f64 = 0.01;

/// The `METRIC_TYPE` column, which is a constant.
pub const METRIC_TYPE: &str = "Allele Frequency";

/// The `Filter` value the rows are cut down to before anything is grouped.
pub const CALLED: &str = "called";

/// `calculateChiSquaredStatistic`.
///
/// Not Pearson's: the expected count in the denominator is replaced by a constant variance, and
/// that variance is SQUARED, so a variance ten times larger divides the statistic by a hundred.
/// A bin holding fewer than two entries contributes nothing; on the reference's own path every bin
/// holds exactly two, one per eval track, so that guard never fires.
pub fn chi_squared_statistic(bins: &[Vec<f64>], variance: f64) -> f64 {
    let sum: f64 = bins
        .iter()
        .map(|entries| {
            if entries.len() >= 2 {
                (entries[0] - entries[1]).powi(2)
            } else {
                0.0
            }
        })
        .sum();
    sum / variance.powi(2)
}

/// The degrees of freedom: the bin count less one, counted over the bins the report holds and not
/// over the ones the data reached.
///
/// The allele-frequency stratifier emits a fixed ladder, so this is the same number whatever the
/// file holds, and a bin no variant reached contributes a term of nought to the statistic while
/// still being counted here.
pub fn degrees_of_freedom(bins: usize) -> f64 {
    bins as f64 - 1.0
}

/// `1 - ChiSquaredDistribution(df).cumulativeProbability(x)`: the upper tail.
///
/// commons-math evaluates the gamma distribution's CDF as the regularized lower incomplete gamma
/// function, which htsjdk-rs already carries, so the p-value here is that function's and not an
/// approximation of it.
pub fn p_value(statistic: f64, degrees_of_freedom: f64) -> f64 {
    if statistic <= 0.0 {
        return 1.0;
    }
    let cumulative = jmath::gamma::regularized_gamma_p(
        degrees_of_freedom / 2.0,
        statistic / 2.0,
        1e-14,
        i32::MAX,
    )
    .unwrap_or(f64::NAN);
    1.0 - cumulative
}

/// The two numbers the metrics file carries, from the bins the report gave.
pub fn metrics(bins: &[Vec<f64>], variance: f64) -> (f64, f64) {
    let statistic = chi_squared_statistic(bins, variance);
    (
        statistic,
        p_value(statistic, degrees_of_freedom(bins.len())),
    )
}

/// Whether the run complains, which is the whole of what the threshold decides.
pub fn complains(p_value: f64, threshold: f64) -> bool {
    p_value < threshold
}

/// The message that complaint carries.
pub fn complaint(p_value: f64) -> String {
    format!(
        "Allele frequencies between your array VCF and the expected VCF do not match with a \
         significant pvalue of {p_value}"
    )
}
