//! `ModelSegments`' two models and the smoothing between them.
//!
//! Ported from `org.broadinstitute.hellbender.tools.copynumber.models` (GATK 4.6.2.0):
//! `CopyRatioModeller`, `CopyRatioSamplers`, `CopyRatioSegmentedData`, `AlleleFractionModeller`,
//! `AlleleFractionSamplers`, `AlleleFractionLikelihoods`, `AlleleFractionInitializer`,
//! `AlleleFractionSegmentedData` and `MultidimensionalModeller`, with the Gibbs sampler of
//! `org.broadinstitute.hellbender.utils.mcmc` written out for the two states.
//!
//! # The draws are the comparison
//!
//! Every chain restarts `new Random(42)` (see [`gatk_engine::copy_number_mcmc`]), and the Gibbs
//! update visits the parameters in the order each state's `LinkedHashMap` was filled: variance,
//! outlier probability, segment means and outlier indicators for the copy ratios; mean bias, bias
//! variance, outlier probability and minor fractions for the allele fractions. Each sampler sees the
//! state as the samplers before it in the same sweep left it.
//!
//! # Sums are streams where the reference streams
//!
//! A segment's allele-fraction log likelihood and the total over segments are `DoubleStream.sum`,
//! which is compensated; the initial minor fractions and the minibatch sums are plain loops. The two
//! are different doubles, and each call site keeps the one the reference wrote.
//!
//! # `Math.exp` inside `logSumExp`
//!
//! A het's likelihood is `NaturalLogUtils.logSumExp` of three terms, which calls `Math.exp`. That
//! function is not portable (htsjdk-rs decision 0014); [`gatk_engine::natural_log_utils::log_sum_exp`]
//! is the repository's measured stand-in, within one ulp. The likelihoods only decide whether a
//! proposal lies on the slice and where Brent's search steps next, and no likelihood is printed.

use gatk_engine::copy_number_mcmc::{
    beta_log_density, beta_sample, commons_mean, commons_variance, deciles, double_stream_average,
    double_stream_sum, java_max, java_min, IllegalArgument, MinibatchSliceSampler,
    GIBBS_RANDOM_SEED,
};
use gatk_engine::java_random::JavaRandom;

use crate::main_entry::Thrown;

/// A locatable point: a copy-ratio interval, a het site, or a segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Span {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Span {
    /// `SimpleInterval.toString`.
    pub fn render(&self) -> String {
        format!("{}:{}-{}", self.contig, self.start, self.end)
    }

    /// Whether a one-base site lies inside, on the same contig.
    pub fn contains_position(&self, contig: &str, position: i32) -> bool {
        self.contig == contig && self.start <= position && position <= self.end
    }
}

/// One allelic count: the site and the two read counts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Het {
    pub contig: String,
    pub position: i32,
    pub ref_count: i32,
    pub alt_count: i32,
}

/// `SimplePosteriorSummary`: the tenth, fiftieth and ninetieth percentiles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PosteriorSummary {
    pub decile10: f64,
    pub decile50: f64,
    pub decile90: f64,
}

impl PosteriorSummary {
    fn of(samples: &[f64]) -> Self {
        let all = deciles(samples);
        PosteriorSummary {
            decile10: all[0],
            decile50: all[4],
            decile90: all[8],
        }
    }
}

/// `ModeledSegment`.
#[derive(Debug, Clone, PartialEq)]
pub struct ModeledSegment {
    pub span: Span,
    pub num_points_copy_ratio: i32,
    pub num_points_allele_fraction: i32,
    pub log2_copy_ratio: PosteriorSummary,
    pub minor_allele_fraction: PosteriorSummary,
}

impl ModeledSegment {
    fn new(
        span: Span,
        num_points_copy_ratio: i32,
        num_points_allele_fraction: i32,
        log2_copy_ratio: PosteriorSummary,
        minor_allele_fraction: PosteriorSummary,
    ) -> Result<Self, Thrown> {
        if !(num_points_copy_ratio > 0 || num_points_allele_fraction > 0) {
            return Err(illegal(format!(
                "Number of copy-ratio points or number of allele-fraction points must be positive: {}",
                span.render()
            )));
        }
        Ok(ModeledSegment {
            span,
            num_points_copy_ratio,
            num_points_allele_fraction,
            log2_copy_ratio,
            minor_allele_fraction,
        })
    }
}

/// An `IllegalArgumentException`.
pub fn illegal(message: impl Into<String>) -> Thrown {
    Thrown::non_user("java.lang.IllegalArgumentException", message)
}

fn from_illegal(error: IllegalArgument) -> Thrown {
    illegal(error.0)
}

/// The nine deciles of one global parameter, named as `ParameterDecileCollection` prints it.
pub type ParameterDeciles = Vec<(&'static str, [f64; 9])>;

/// `CopyRatioSegmentedData`: the copy ratios whose midpoints fall in each segment, in order.
struct CopyRatioData {
    /// Every copy ratio, in record order, which is what the range and the point count come from.
    all_values: Vec<f64>,
    /// `(value, segment index)` for each point inside a segment, segment by segment.
    points: Vec<(f64, usize)>,
    /// The `[start, end)` range of `points` for each segment.
    ranges: Vec<(usize, usize)>,
}

impl CopyRatioData {
    fn new(copy_ratios: &[(Span, f64)], segments: &[Span]) -> Self {
        let mut points = Vec::new();
        let mut ranges = Vec::with_capacity(segments.len());
        for (segment_index, segment) in segments.iter().enumerate() {
            let start = points.len();
            for (span, value) in copy_ratios {
                let midpoint = (span.start + span.end) / 2;
                if segment.contains_position(&span.contig, midpoint) {
                    points.push((*value, segment_index));
                }
            }
            ranges.push((start, points.len()));
        }
        CopyRatioData {
            all_values: copy_ratios.iter().map(|(_, value)| *value).collect(),
            points,
            ranges,
        }
    }

    fn num_points(&self) -> usize {
        self.all_values.len()
    }

    fn segment(&self, index: usize) -> &[(f64, usize)] {
        let (start, end) = self.ranges[index];
        &self.points[start..end]
    }
}

/// `Stream.min(Double::compareTo)` and `max`, which order NaN last and -0.0 before 0.0.
fn java_extreme(values: &[f64], greatest: bool) -> f64 {
    let mut iter = values.iter().copied();
    let Some(mut best) = iter.next() else {
        return f64::NAN;
    };
    for value in iter {
        let ordering = gatk_engine::persistence_optimizer::java_compare(value, best);
        let better = if greatest {
            ordering == std::cmp::Ordering::Greater
        } else {
            ordering == std::cmp::Ordering::Less
        };
        if better {
            best = value;
        }
    }
    best
}

/// `CopyRatioSamplers.normalTerm`.
fn normal_term(quantity: f64, mean: f64, variance: f64) -> f64 {
    (quantity - mean) * (quantity - mean) / (2.0 * variance)
}

/// `NaturalLogUtils.logSumLog`, with `Math.log` and `FastMath.exp`.
fn log_sum_log(a: f64, b: f64) -> f64 {
    if a > b {
        a + jmath::math::log(1.0 + jmath::fast_math::exp(b - a))
    } else {
        b + jmath::math::log(1.0 + jmath::fast_math::exp(a - b))
    }
}

const LOG2_COPY_RATIO_MIN: f64 = -50.0;
const LOG2_COPY_RATIO_MAX: f64 = 10.0;
const VARIANCE_MIN: f64 = 1e-6;
const OUTLIER_PROBABILITY_INITIAL: f64 = 0.05;
const OUTLIER_PROBABILITY_PRIOR_ALPHA: f64 = 5.0;
const OUTLIER_PROBABILITY_PRIOR_BETA: f64 = 95.0;
const APPROX_THRESHOLD: f64 = 0.1;

#[derive(Debug, Clone)]
struct CopyRatioState {
    variance: f64,
    outlier_probability: f64,
    segment_means: Vec<f64>,
    outlier_indicators: Vec<bool>,
}

/// `CopyRatioModeller`, fit once.
struct CopyRatioModeller {
    samples: Vec<CopyRatioState>,
    num_burn_in: usize,
}

impl CopyRatioModeller {
    fn fit(
        copy_ratios: &[(Span, f64)],
        segments: &[Span],
        num_samples: usize,
        num_burn_in: usize,
    ) -> Result<Self, Thrown> {
        let data = CopyRatioData::new(copy_ratios, segments);
        let num_segments = segments.len();

        let range_or_nan =
            java_extreme(&data.all_values, true) - java_extreme(&data.all_values, false);
        let data_range = if range_or_nan.is_nan() {
            LOG2_COPY_RATIO_MAX - LOG2_COPY_RATIO_MIN
        } else {
            range_or_nan
        };
        let variances: Vec<f64> = (0..num_segments)
            .map(|s| commons_variance(&data.segment(s).iter().map(|p| p.0).collect::<Vec<f64>>()))
            .filter(|v| !v.is_nan())
            .collect();
        let variance_or_nan = double_stream_average(&variances).unwrap_or(f64::NAN);
        let variance_estimate = if variance_or_nan.is_nan() {
            VARIANCE_MIN
        } else {
            java_max(variance_or_nan, VARIANCE_MIN)
        };
        let variance_width = 2.0 * variance_estimate;
        let variance_max = java_max(10.0 * variance_estimate, data_range * data_range);
        let mean_width =
            (variance_estimate * num_segments as f64 / data.num_points() as f64).sqrt();
        let segment_means: Vec<f64> = (0..num_segments)
            .map(|s| commons_mean(&data.segment(s).iter().map(|p| p.0).collect::<Vec<f64>>()))
            .map(|m| java_max(LOG2_COPY_RATIO_MIN, java_min(LOG2_COPY_RATIO_MAX, m)))
            .collect();
        let outlier_uniform_log_likelihood = -jmath::math::log(data_range);

        let mut state = CopyRatioState {
            variance: variance_estimate,
            outlier_probability: OUTLIER_PROBABILITY_INITIAL,
            segment_means,
            outlier_indicators: vec![false; data.num_points()],
        };

        let mut rng = JavaRandom::new(GIBBS_RANDOM_SEED);
        let mut samples = Vec::with_capacity(num_samples);
        samples.push(state.clone());
        for _ in 1..num_samples {
            // Variance.
            let non_outliers: Vec<(f64, usize)> = data
                .points
                .iter()
                .enumerate()
                .filter(|(index, _)| !state.outlier_indicators[*index])
                .map(|(_, point)| *point)
                .collect();
            let means = state.segment_means.clone();
            state.variance = MinibatchSliceSampler::new(
                &non_outliers,
                |_| Ok(0.0),
                |point: &(f64, usize), new_variance: f64| {
                    -0.5 * jmath::fast_math::log(new_variance)
                        - normal_term(point.0, means[point.1], new_variance)
                },
                VARIANCE_MIN,
                variance_max,
                variance_width,
                1000,
                APPROX_THRESHOLD,
            )
            .and_then(|mut sampler| sampler.sample(&mut rng, state.variance))
            .map_err(from_illegal)?;

            // Outlier probability.
            let mut num_outliers = 0usize;
            for index in 0..data.num_points() {
                if indicator(&state.outlier_indicators, index)? {
                    num_outliers += 1;
                }
            }
            state.outlier_probability = beta_sample(
                &mut rng,
                OUTLIER_PROBABILITY_PRIOR_ALPHA + num_outliers as f64,
                OUTLIER_PROBABILITY_PRIOR_BETA + data.num_points() as f64 - num_outliers as f64,
            );

            // Segment means.
            let mut new_means = Vec::with_capacity(num_segments);
            for segment_index in 0..num_segments {
                let (start, end) = data.ranges[segment_index];
                if start == end {
                    new_means.push(f64::NAN);
                    continue;
                }
                let indexed: Vec<(f64, usize)> = (start..end)
                    .map(|index| (data.points[index].0, index))
                    .collect();
                let indicators = &state.outlier_indicators;
                let variance = state.variance;
                let mean = MinibatchSliceSampler::new(
                    &indexed,
                    |_| Ok(0.0),
                    |point: &(f64, usize), new_mean: f64| {
                        if indicators.get(point.1).copied().unwrap_or(false) {
                            0.0
                        } else {
                            -normal_term(point.0, new_mean, variance)
                        }
                    },
                    LOG2_COPY_RATIO_MIN,
                    LOG2_COPY_RATIO_MAX,
                    mean_width,
                    100,
                    APPROX_THRESHOLD,
                )
                .and_then(|mut sampler| {
                    sampler.sample(&mut rng, state.segment_means[segment_index])
                })
                .map_err(from_illegal)?;
                new_means.push(mean);
            }
            state.segment_means = new_means;

            // Outlier indicators.
            let outlier_log_probability =
                jmath::math::log(state.outlier_probability) + outlier_uniform_log_likelihood;
            let prefactor = jmath::math::log(
                (1.0 - state.outlier_probability)
                    / (2.0 * std::f64::consts::PI * state.variance).sqrt(),
            );
            let mut indicators = Vec::with_capacity(data.points.len());
            for segment_index in 0..num_segments {
                for point in data.segment(segment_index) {
                    let not_outlier = prefactor
                        - normal_term(point.0, state.segment_means[segment_index], state.variance);
                    let conditional = jmath::fast_math::exp(
                        outlier_log_probability - log_sum_log(outlier_log_probability, not_outlier),
                    );
                    indicators.push(rng.next_double() < conditional);
                }
            }
            state.outlier_indicators = indicators;

            samples.push(state.clone());
        }
        Ok(CopyRatioModeller {
            samples,
            num_burn_in,
        })
    }

    fn kept(&self) -> &[CopyRatioState] {
        &self.samples[self.num_burn_in..]
    }

    fn segment_means_summaries(&self) -> Vec<PosteriorSummary> {
        let num_segments = self.samples[0].segment_means.len();
        (0..num_segments)
            .map(|s| {
                PosteriorSummary::of(
                    &self
                        .kept()
                        .iter()
                        .map(|state| state.segment_means[s])
                        .collect::<Vec<f64>>(),
                )
            })
            .collect()
    }

    fn global_deciles(&self) -> ParameterDeciles {
        let variance: Vec<f64> = self.kept().iter().map(|s| s.variance).collect();
        let outlier: Vec<f64> = self.kept().iter().map(|s| s.outlier_probability).collect();
        vec![
            ("VARIANCE", deciles(&variance)),
            ("OUTLIER_PROBABILITY", deciles(&outlier)),
        ]
    }
}

/// `List<Boolean>.get(index)`, which throws past the end.
fn indicator(indicators: &[bool], index: usize) -> Result<bool, Thrown> {
    indicators
        .get(index)
        .copied()
        .ok_or_else(|| out_of_bounds(index, indicators.len()))
}

fn out_of_bounds(index: usize, length: usize) -> Thrown {
    Thrown::non_user(
        "java.lang.IndexOutOfBoundsException",
        format!("Index {index} out of bounds for length {length}"),
    )
}

/// `AlleleFractionGlobalParameters`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct GlobalParameters {
    mean_bias: f64,
    bias_variance: f64,
    outlier_probability: f64,
}

impl GlobalParameters {
    fn alpha(&self) -> f64 {
        self.mean_bias * self.mean_bias / self.bias_variance
    }

    fn beta(&self) -> f64 {
        self.mean_bias / self.bias_variance
    }
}

const AF_EPSILON: f64 = 1e-10;

/// `AlleleFractionLikelihoods.log`: `Math.log` clamped below at `EPSILON`.
fn af_log(x: f64) -> f64 {
    jmath::math::log(java_max(AF_EPSILON, x))
}

fn bias_posterior_mode(alpha: f64, beta: f64, f: f64, a: i32, r: i32) -> f64 {
    let w = (1.0 - f) * (f64::from(a) - alpha + 1.0) + beta * f;
    java_max(
        ((w * w + 4.0 * beta * f * (1.0 - f) * (f64::from(r) + alpha - 1.0)).sqrt() - w)
            / (2.0 * beta * (1.0 - f)),
        AF_EPSILON,
    )
}

fn bias_posterior_curvature(alpha: f64, f: f64, r: i32, n: i32, lambda0: f64) -> f64 {
    let y = (1.0 - f) / (f + (1.0 - f) * lambda0);
    f64::from(n) * y * y - (f64::from(r) + alpha - 1.0) / (lambda0 * lambda0)
}

/// `AlleleFractionLikelihoods.hetLogLikelihood`.
fn het_log_likelihood(parameters: &GlobalParameters, minor_fraction: f64, het: &Het) -> f64 {
    let alpha = parameters.alpha();
    let beta = parameters.beta();
    let pi = parameters.outlier_probability;

    let log_pi = af_log(pi);
    let log_not_pi = af_log((1.0 - pi) / 2.0);
    let logc_common = alpha * af_log(beta) - jmath::gamma::log_gamma(alpha);
    let major_fraction = 1.0 - minor_fraction;
    let log_minor = af_log(minor_fraction);
    let log_major = af_log(major_fraction);

    let a = het.alt_count;
    let r = het.ref_count;
    let n = a + r;
    let (af, rf, nf) = (f64::from(a), f64::from(r), f64::from(n));

    let lambda0_alt = bias_posterior_mode(alpha, beta, minor_fraction, a, r);
    let kappa_alt = bias_posterior_curvature(alpha, minor_fraction, r, n, lambda0_alt);
    let rho_alt = java_max(1.0 - kappa_alt * lambda0_alt * lambda0_alt, AF_EPSILON);
    let tau_alt = java_max(-kappa_alt * lambda0_alt, AF_EPSILON);
    let logc_alt = logc_common
        + af * log_minor
        + rf * log_major
        + (rf + alpha - rho_alt) * af_log(lambda0_alt)
        + (tau_alt - beta) * lambda0_alt
        - nf * af_log(minor_fraction + major_fraction * lambda0_alt);
    let alt_minor =
        log_not_pi + logc_alt + jmath::gamma::log_gamma(rho_alt) - rho_alt * af_log(tau_alt);

    let lambda0_ref = bias_posterior_mode(alpha, beta, major_fraction, a, r);
    let kappa_ref = bias_posterior_curvature(alpha, major_fraction, r, n, lambda0_ref);
    let rho_ref = java_max(1.0 - kappa_ref * lambda0_ref * lambda0_ref, AF_EPSILON);
    let tau_ref = java_max(-kappa_ref * lambda0_ref, AF_EPSILON);
    let logc_ref = logc_common
        + af * log_major
        + rf * log_minor
        + (rf + alpha - rho_ref) * af_log(lambda0_ref)
        + (tau_ref - beta) * lambda0_ref
        - nf * af_log(major_fraction + minor_fraction * lambda0_ref);
    let ref_minor =
        log_not_pi + logc_ref + jmath::gamma::log_gamma(rho_ref) - rho_ref * af_log(tau_ref);

    let binomial = binomial_coefficient_log(n, a);
    let outlier = log_pi - jmath::math::log(f64::from(a + r + 1)) - binomial;
    log_sum_exp(&[alt_minor, ref_minor, outlier])
}

/// `CombinatoricsUtils.binomialCoefficientLog(n, k)`.
///
/// jmath refuses the exact route's arm for `n` from 62 to 66, which a site of that depth reaches.
/// That arm computes the coefficient itself, exactly, as a `long` (it only orders the
/// multiplications so that no intermediate overflows), so the coefficient is formed exactly here
/// in 128 bits and logged the same way: `FastMath.log` of the `long` widened to `double`.
fn binomial_coefficient_log(n: i32, k: i32) -> f64 {
    let (n, k) = (i64::from(n), i64::from(k));
    let trivial = n == k || k == 0 || k == 1 || k == n - 1;
    if (62..=66).contains(&n) && !trivial && k >= 0 && k <= n {
        let k = k.min(n - k);
        let mut result: u128 = 1;
        for j in 1..=k {
            result = result * (n - k + j) as u128 / j as u128;
        }
        return jmath::fast_math::log(result as i64 as f64);
    }
    jmath::combinatorics::binomial_coefficient_log(n, k).unwrap_or(f64::NAN)
}

/// `NaturalLogUtils.logSumExp`, with the host's `exp` standing in for `Math.exp`.
///
/// Not [`gatk_engine::natural_log_utils::log_sum_exp`], whose `exp` is FDLIBM: that stand-in is
/// within one ulp of `Math.exp`, and one ulp is enough here. On a row segmented at a penalty of
/// 0.1 (121 segments), FDLIBM moved one Brent step of the allele-fraction initializer and every
/// allele-fraction decile after it, where the host's `exp` answered the reference's bytes. Like
/// `AllelePseudoDepth`'s `pow`, the host libm is the closer stand-in on the points this reaches;
/// the covering array is what says so on the platform CI runs.
///
/// The refusal on a non-finite sum is the reference's `IllegalArgumentException`; a likelihood
/// has nowhere to report it from, so it answers NaN, which no comparison accepts.
fn log_sum_exp(values: &[f64]) -> f64 {
    // `MathUtils.maxElementIndex`: the first index holding the maximum.
    let mut max_index = 0;
    for (index, value) in values.iter().enumerate().skip(1) {
        if *value > values[max_index] {
            max_index = index;
        }
    }
    let max_value = values[max_index];
    if max_value == f64::NEG_INFINITY {
        return max_value;
    }
    let mut sum = 1.0;
    for (index, value) in values.iter().enumerate() {
        if index == max_index || *value == f64::NEG_INFINITY {
            continue;
        }
        sum += (value - max_value).exp();
    }
    if sum.is_nan() || sum == f64::INFINITY {
        return f64::NAN;
    }
    max_value
        + if sum != 1.0 {
            jmath::math::log(sum)
        } else {
            0.0
        }
}

/// `segmentLogLikelihood`, a stream sum over the segment's hets.
fn segment_log_likelihood(parameters: &GlobalParameters, minor_fraction: f64, hets: &[Het]) -> f64 {
    let terms: Vec<f64> = hets
        .iter()
        .map(|het| het_log_likelihood(parameters, minor_fraction, het))
        .collect();
    double_stream_sum(&terms)
}

/// `logLikelihood`, a stream sum over the segments.
fn total_log_likelihood(
    parameters: &GlobalParameters,
    minor_fractions: &[f64],
    data: &AlleleFractionData,
) -> f64 {
    let terms: Vec<f64> = (0..data.ranges.len())
        .map(|s| segment_log_likelihood(parameters, minor_fractions[s], data.segment(s)))
        .collect();
    double_stream_sum(&terms)
}

/// `AlleleFractionSegmentedData`.
struct AlleleFractionData {
    hets: Vec<Het>,
    segment_of: Vec<usize>,
    ranges: Vec<(usize, usize)>,
}

impl AlleleFractionData {
    fn new(hets: &[Het], segments: &[Span]) -> Self {
        let mut ordered = Vec::new();
        let mut segment_of = Vec::new();
        let mut ranges = Vec::with_capacity(segments.len());
        for (segment_index, segment) in segments.iter().enumerate() {
            let start = ordered.len();
            for het in hets {
                if segment.contains_position(&het.contig, het.position) {
                    ordered.push(het.clone());
                    segment_of.push(segment_index);
                }
            }
            ranges.push((start, ordered.len()));
        }
        AlleleFractionData {
            hets: ordered,
            segment_of,
            ranges,
        }
    }

    fn segment(&self, index: usize) -> &[Het] {
        let (start, end) = self.ranges[index];
        &self.hets[start..end]
    }
}

const MAX_REASONABLE_OUTLIER_PROBABILITY: f64 = 0.15;
const MAX_REASONABLE_MEAN_BIAS: f64 = 5.0;
const MAX_REASONABLE_BIAS_VARIANCE: f64 = 0.5;
const MAX_MINOR_ALLELE_FRACTION: f64 = 0.5;
const MIN_MINOR_FRACTION_SAMPLING_WIDTH: f64 = 1e-3;

/// `OptimizationUtils.argmax`: commons-math's Brent at `0.001` both ways and 1000 evaluations.
fn argmax(objective: impl Fn(f64) -> f64, min: f64, max: f64, guess: f64) -> Result<f64, Thrown> {
    jmath::brent::maximize(objective, min, max, guess, 0.001, 0.001, 1000)
        .map(|pair| pair.point)
        .map_err(|error| Thrown::non_user(error.class(), error.message()))
}

/// `AlleleFractionInitializer`: the mode of the likelihood, by coordinate-wise Brent searches.
fn initialize(data: &AlleleFractionData) -> Result<(GlobalParameters, Vec<f64>), Thrown> {
    let mut global = GlobalParameters {
        mean_bias: 1.0,
        bias_variance: 0.05,
        outlier_probability: 0.01,
    };
    let mut minor_fractions: Vec<f64> = (0..data.ranges.len())
        .map(|s| {
            let mut minor_count = 0.0;
            let mut total_count = 0.0;
            for het in data.segment(s) {
                let a = het.alt_count;
                let r = het.ref_count;
                let responsibility = match jmath::beta::regularized_beta(
                    0.5,
                    f64::from(a) + 1.0,
                    f64::from(r) + 1.0,
                ) {
                    Ok(value) => value,
                    Err(_) => {
                        if a < r {
                            1.0
                        } else {
                            0.0
                        }
                    }
                };
                minor_count +=
                    responsibility * f64::from(a) + (1.0 - responsibility) * f64::from(r);
                total_count += f64::from(a + r);
            }
            (minor_count + 1.0) / (total_count + 2.0)
        })
        .collect();

    let mut next = f64::NEG_INFINITY;
    let mut iteration = 1;
    loop {
        let previous = next;
        let old = global;
        let mean_bias = argmax(
            |value| {
                total_log_likelihood(
                    &GlobalParameters {
                        mean_bias: value,
                        ..old
                    },
                    &minor_fractions,
                    data,
                )
            },
            0.0,
            MAX_REASONABLE_MEAN_BIAS,
            old.mean_bias,
        )?;
        let bias_variance = argmax(
            |value| {
                total_log_likelihood(
                    &GlobalParameters {
                        bias_variance: value,
                        ..old
                    },
                    &minor_fractions,
                    data,
                )
            },
            0.0,
            MAX_REASONABLE_BIAS_VARIANCE,
            old.bias_variance,
        )?;
        let outlier_probability = argmax(
            |value| {
                total_log_likelihood(
                    &GlobalParameters {
                        outlier_probability: value,
                        ..old
                    },
                    &minor_fractions,
                    data,
                )
            },
            0.0,
            MAX_REASONABLE_OUTLIER_PROBABILITY,
            old.outlier_probability,
        )?;
        global = GlobalParameters {
            mean_bias,
            bias_variance,
            outlier_probability,
        };
        let mut estimated = Vec::with_capacity(minor_fractions.len());
        for (segment, guess) in minor_fractions.iter().enumerate() {
            let hets = data.segment(segment);
            estimated.push(argmax(
                |f| segment_log_likelihood(&global, f, hets),
                0.0,
                MAX_MINOR_ALLELE_FRACTION,
                *guess,
            )?);
        }
        minor_fractions = estimated;
        next = total_log_likelihood(&global, &minor_fractions, data);

        iteration += 1;
        if !(iteration < 50 && next - previous > 0.5) {
            break;
        }
    }
    Ok((global, minor_fractions))
}

/// `approximatePosteriorWidthAtMode`.
fn width_at_mode(log_pdf: impl Fn(f64) -> f64, mode: f64) -> f64 {
    let abs_mode = mode.abs();
    let epsilon = java_min(1e-6, abs_mode / 2.0);
    let default_width = abs_mode / 10.0;
    let second_derivative = (log_pdf(mode + epsilon) - 2.0 * log_pdf(mode)
        + log_pdf(mode - epsilon))
        / (epsilon * epsilon);
    if second_derivative < 0.0 {
        (-1.0 / second_derivative).sqrt()
    } else {
        default_width
    }
}

#[derive(Debug, Clone)]
struct AlleleFractionState {
    global: GlobalParameters,
    minor_fractions: Vec<f64>,
}

/// `AlleleFractionModeller`, fit once.
struct AlleleFractionModeller {
    samples: Vec<AlleleFractionState>,
    num_burn_in: usize,
}

impl AlleleFractionModeller {
    fn fit(
        hets: &[Het],
        segments: &[Span],
        prior_alpha: f64,
        num_samples: usize,
        num_burn_in: usize,
    ) -> Result<Self, Thrown> {
        let data = AlleleFractionData::new(hets, segments);
        let (initial, initial_minor) = initialize(&data)?;

        let mean_bias_width = width_at_mode(
            |value| {
                total_log_likelihood(
                    &GlobalParameters {
                        mean_bias: value,
                        ..initial
                    },
                    &initial_minor,
                    &data,
                )
            },
            initial.mean_bias,
        );
        let bias_variance_width = width_at_mode(
            |value| {
                total_log_likelihood(
                    &GlobalParameters {
                        bias_variance: value,
                        ..initial
                    },
                    &initial_minor,
                    &data,
                )
            },
            initial.bias_variance,
        );
        let outlier_width = width_at_mode(
            |value| {
                total_log_likelihood(
                    &GlobalParameters {
                        outlier_probability: value,
                        ..initial
                    },
                    &initial_minor,
                    &data,
                )
            },
            initial.outlier_probability,
        );
        let minor_widths: Vec<f64> = (0..segments.len())
            .map(|s| {
                width_at_mode(
                    |f| segment_log_likelihood(&initial, f, data.segment(s)),
                    initial_minor[s],
                )
            })
            .map(|w| java_max(w, MIN_MINOR_FRACTION_SAMPLING_WIDTH))
            .collect();

        let mut state = AlleleFractionState {
            global: initial,
            minor_fractions: initial_minor,
        };
        let mut rng = JavaRandom::new(GIBBS_RANDOM_SEED);
        let mut samples = Vec::with_capacity(num_samples);
        samples.push(state.clone());
        let all: Vec<usize> = (0..data.hets.len()).collect();
        for _ in 1..num_samples {
            // Mean bias.
            {
                let current = state.clone();
                state.global.mean_bias = MinibatchSliceSampler::new(
                    &all,
                    |_| Ok(0.0),
                    |index: &usize, value: f64| {
                        het_log_likelihood(
                            &GlobalParameters {
                                mean_bias: value,
                                ..current.global
                            },
                            current.minor_fractions[data.segment_of[*index]],
                            &data.hets[*index],
                        )
                    },
                    0.0,
                    MAX_REASONABLE_MEAN_BIAS,
                    mean_bias_width,
                    1000,
                    APPROX_THRESHOLD,
                )
                .and_then(|mut sampler| sampler.sample(&mut rng, current.global.mean_bias))
                .map_err(from_illegal)?;
            }
            // Bias variance.
            {
                let current = state.clone();
                state.global.bias_variance = MinibatchSliceSampler::new(
                    &all,
                    |_| Ok(0.0),
                    |index: &usize, value: f64| {
                        het_log_likelihood(
                            &GlobalParameters {
                                bias_variance: value,
                                ..current.global
                            },
                            current.minor_fractions[data.segment_of[*index]],
                            &data.hets[*index],
                        )
                    },
                    1e-10,
                    MAX_REASONABLE_BIAS_VARIANCE,
                    bias_variance_width,
                    1000,
                    APPROX_THRESHOLD,
                )
                .and_then(|mut sampler| sampler.sample(&mut rng, current.global.bias_variance))
                .map_err(from_illegal)?;
            }
            // Outlier probability.
            {
                let current = state.clone();
                state.global.outlier_probability = MinibatchSliceSampler::new(
                    &all,
                    |_| Ok(0.0),
                    |index: &usize, value: f64| {
                        het_log_likelihood(
                            &GlobalParameters {
                                outlier_probability: value,
                                ..current.global
                            },
                            current.minor_fractions[data.segment_of[*index]],
                            &data.hets[*index],
                        )
                    },
                    0.0,
                    MAX_REASONABLE_OUTLIER_PROBABILITY,
                    outlier_width,
                    1000,
                    APPROX_THRESHOLD,
                )
                .and_then(|mut sampler| {
                    sampler.sample(&mut rng, current.global.outlier_probability)
                })
                .map_err(from_illegal)?;
            }
            // Minor fractions.
            {
                let global = state.global;
                let mut minor_fractions = Vec::with_capacity(segments.len());
                for (segment_index, minor_width) in
                    minor_widths.iter().enumerate().take(segments.len())
                {
                    let hets_in_segment = data.segment(segment_index);
                    if hets_in_segment.is_empty() {
                        minor_fractions.push(f64::NAN);
                        continue;
                    }
                    let value = MinibatchSliceSampler::new(
                        hets_in_segment,
                        |f| beta_log_density(2.0 * f, prior_alpha, 1.0),
                        |het: &Het, f: f64| het_log_likelihood(&global, f, het),
                        0.0,
                        0.5,
                        *minor_width,
                        10,
                        APPROX_THRESHOLD,
                    )
                    .and_then(|mut sampler| {
                        sampler.sample(&mut rng, state.minor_fractions[segment_index])
                    })
                    .map_err(from_illegal)?;
                    minor_fractions.push(value);
                }
                state.minor_fractions = minor_fractions;
            }
            samples.push(state.clone());
        }
        Ok(AlleleFractionModeller {
            samples,
            num_burn_in,
        })
    }

    fn kept(&self) -> &[AlleleFractionState] {
        &self.samples[self.num_burn_in..]
    }

    fn minor_fraction_summaries(&self) -> Vec<PosteriorSummary> {
        let num_segments = self.samples[0].minor_fractions.len();
        (0..num_segments)
            .map(|s| {
                PosteriorSummary::of(
                    &self
                        .kept()
                        .iter()
                        .map(|state| state.minor_fractions[s])
                        .collect::<Vec<f64>>(),
                )
            })
            .collect()
    }

    fn global_deciles(&self) -> ParameterDeciles {
        let pick = |f: fn(&GlobalParameters) -> f64| -> [f64; 9] {
            deciles(
                &self
                    .kept()
                    .iter()
                    .map(|s| f(&s.global))
                    .collect::<Vec<f64>>(),
            )
        };
        vec![
            ("MEAN_BIAS", pick(|g| g.mean_bias)),
            ("BIAS_VARIANCE", pick(|g| g.bias_variance)),
            ("OUTLIER_PROBABILITY", pick(|g| g.outlier_probability)),
        ]
    }
}

/// The chain lengths `MultidimensionalModeller` was built with.
#[derive(Debug, Clone, Copy)]
pub struct ChainLengths {
    pub num_samples_copy_ratio: usize,
    pub num_burn_in_copy_ratio: usize,
    pub num_samples_allele_fraction: usize,
    pub num_burn_in_allele_fraction: usize,
}

/// `MultidimensionalModeller`.
pub struct MultidimensionalModeller<'a> {
    copy_ratios: &'a [(Span, f64)],
    hets: &'a [Het],
    prior_alpha: f64,
    lengths: ChainLengths,
    current_segments: Vec<Span>,
    pub modeled_segments: Vec<ModeledSegment>,
    is_model_fit: bool,
    copy_ratio_deciles: ParameterDeciles,
    allele_fraction_deciles: ParameterDeciles,
}

impl<'a> MultidimensionalModeller<'a> {
    /// The constructor, which fits the initial model.
    pub fn new(
        segments: Vec<Span>,
        copy_ratios: &'a [(Span, f64)],
        hets: &'a [Het],
        prior_alpha: f64,
        lengths: ChainLengths,
    ) -> Result<Self, Thrown> {
        if segments.is_empty() {
            return Err(illegal("Number of segments must be positive."));
        }
        let mut modeller = MultidimensionalModeller {
            copy_ratios,
            hets,
            prior_alpha,
            lengths,
            current_segments: segments,
            modeled_segments: Vec::new(),
            is_model_fit: false,
            copy_ratio_deciles: Vec::new(),
            allele_fraction_deciles: Vec::new(),
        };
        modeller.fit_model()?;
        Ok(modeller)
    }

    fn fit_model(&mut self) -> Result<(), Thrown> {
        let copy_ratio = CopyRatioModeller::fit(
            self.copy_ratios,
            &self.current_segments,
            self.lengths.num_samples_copy_ratio,
            self.lengths.num_burn_in_copy_ratio,
        )?;
        let allele_fraction = AlleleFractionModeller::fit(
            self.hets,
            &self.current_segments,
            self.prior_alpha,
            self.lengths.num_samples_allele_fraction,
            self.lengths.num_burn_in_allele_fraction,
        )?;
        let means = copy_ratio.segment_means_summaries();
        let fractions = allele_fraction.minor_fraction_summaries();
        self.modeled_segments.clear();
        for (index, segment) in self.current_segments.iter().enumerate() {
            let num_copy_ratio = self
                .copy_ratios
                .iter()
                .filter(|(span, _)| {
                    segment.contains_position(&span.contig, (span.start + span.end) / 2)
                })
                .count() as i32;
            let num_allele_fraction = self
                .hets
                .iter()
                .filter(|het| segment.contains_position(&het.contig, het.position))
                .count() as i32;
            self.modeled_segments.push(ModeledSegment::new(
                segment.clone(),
                num_copy_ratio,
                num_allele_fraction,
                means[index],
                fractions[index],
            )?);
        }
        self.copy_ratio_deciles = copy_ratio.global_deciles();
        self.allele_fraction_deciles = allele_fraction.global_deciles();
        self.is_model_fit = true;
        Ok(())
    }

    /// `smoothSegments`.
    pub fn smooth_segments(
        &mut self,
        max_iterations: i32,
        iterations_per_fit: i32,
        threshold_copy_ratio: f64,
        threshold_allele_fraction: f64,
    ) -> Result<(), Thrown> {
        for iteration in 1..=max_iterations {
            let previous = self.modeled_segments.len();
            let refit = iterations_per_fit > 0 && iteration % iterations_per_fit == 0;
            let merged = merge_similar_segments(
                &self.modeled_segments,
                threshold_copy_ratio,
                threshold_allele_fraction,
            );
            self.current_segments = merged.iter().map(|s| s.span.clone()).collect();
            if refit {
                self.fit_model()?;
            } else {
                self.modeled_segments = merged;
                self.is_model_fit = false;
            }
            if self.modeled_segments.len() == previous {
                break;
            }
        }
        if !self.is_model_fit {
            self.fit_model()?;
        }
        Ok(())
    }

    /// The global parameters' deciles, copy ratio first. `ensureModelIsFit` refits first when the
    /// last smoothing left the model unfit, which `smoothSegments` never does.
    pub fn parameter_deciles(&mut self) -> Result<(ParameterDeciles, ParameterDeciles), Thrown> {
        if !self.is_model_fit {
            self.fit_model()?;
        }
        Ok((
            self.copy_ratio_deciles.clone(),
            self.allele_fraction_deciles.clone(),
        ))
    }
}

fn are_similar(first: &PosteriorSummary, second: &PosteriorSummary, threshold: f64) -> bool {
    if first.decile50.is_nan() || second.decile50.is_nan() {
        return true;
    }
    let difference = (first.decile50 - second.decile50).abs();
    difference < threshold * (first.decile90 - first.decile10)
        || difference < threshold * (second.decile90 - second.decile10)
}

fn merge_summaries(first: &PosteriorSummary, second: &PosteriorSummary) -> PosteriorSummary {
    if first.decile50.is_nan() && !second.decile50.is_nan() {
        return *second;
    }
    if (!first.decile50.is_nan() && second.decile50.is_nan())
        || (first.decile50.is_nan() && second.decile50.is_nan())
    {
        return *first;
    }
    let sd1 = (first.decile90 - first.decile10) / 2.0;
    let sd2 = (second.decile90 - second.decile10) / 2.0;
    // `Math.pow(sd, 2.)`, the rounded square.
    let variance = 1.0 / (1.0 / (sd1 * sd1) + 1.0 / (sd2 * sd2));
    let mean = (first.decile50 / (sd1 * sd1) + second.decile50 / (sd2 * sd2)) * variance;
    let sd = variance.sqrt();
    // `new SimplePosteriorSummary(mean, mean - standardDeviation, mean + standardDeviation)`, whose
    // parameters are `(decile10, decile50, decile90)`: the merged mean lands in the TENTH
    // percentile and `mean - sd` in the median. The next comparison reads them in those places.
    PosteriorSummary {
        decile10: mean,
        decile50: mean - sd,
        decile90: mean + sd,
    }
}

/// `SimilarSegmentUtils.mergeSimilarSegments`: one pass, each segment merged with its right
/// neighbour for as long as the two are similar in both summaries.
fn merge_similar_segments(
    segments: &[ModeledSegment],
    threshold_copy_ratio: f64,
    threshold_allele_fraction: f64,
) -> Vec<ModeledSegment> {
    let mut merged = segments.to_vec();
    let mut index: i64 = 0;
    while index < merged.len() as i64 - 1 {
        let i = index as usize;
        let first = &merged[i];
        let second = &merged[i + 1];
        if first.span.contig == second.span.contig
            && are_similar(
                &first.log2_copy_ratio,
                &second.log2_copy_ratio,
                threshold_copy_ratio,
            )
            && are_similar(
                &first.minor_allele_fraction,
                &second.minor_allele_fraction,
                threshold_allele_fraction,
            )
        {
            let joined = ModeledSegment {
                span: Span {
                    contig: first.span.contig.clone(),
                    start: first.span.start.min(second.span.start),
                    end: first.span.end.max(second.span.end),
                },
                num_points_copy_ratio: first.num_points_copy_ratio + second.num_points_copy_ratio,
                num_points_allele_fraction: first.num_points_allele_fraction
                    + second.num_points_allele_fraction,
                log2_copy_ratio: merge_summaries(&first.log2_copy_ratio, &second.log2_copy_ratio),
                minor_allele_fraction: merge_summaries(
                    &first.minor_allele_fraction,
                    &second.minor_allele_fraction,
                ),
            };
            merged[i] = joined;
            merged.remove(i + 1);
            index -= 1;
        }
        index += 1;
    }
    merged
}
