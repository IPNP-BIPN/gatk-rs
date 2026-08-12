//! `RecalDatum` and `EventType`, ported from
//! `org.broadinstitute.hellbender.utils.recalibration` (GATK 4.6.2.0).
//!
//! A recalibration table is an array of these. Every number `BaseRecalibrator` writes into a report
//! and every number `ApplyBQSR` reads back out of one comes from a `RecalDatum`, so this is where
//! the arithmetic of base quality score recalibration is settled.
//!
//! The four-tuple it holds is `(reported quality, empirical quality, observations, errors)`, and
//! the interesting part is that none of the four is stored the way it is read.
//!
//! # The mismatch count is stored multiplied
//!
//! ```java
//! private static final double MULTIPLIER = 100000.0;  //See discussion in numMismatches about what the multiplier is.
//! ```
//!
//! and the discussion it points at says:
//!
//! ```text
//! The value of the MULTIPLIER that we found to give consistent results insensitive to sorting is 10000.0;
//! ```
//!
//! One zero apart. The constant is the behaviour and the comment is stale, and the difference is
//! observable: the scaling exists precisely to move the sum into a range where addition is
//! associative, so a port that believed the comment would drift on large tables.
//!
//! Every setter and the constructor multiply on the way in; every getter divides on the way out.
//! [`RecalDatum::combine`] goes through both, because it calls the public getter of the other datum
//! and hands the result to the incrementer, so the value it adds has been divided and multiplied
//! again rather than copied.
//!
//! # The empirical quality is an integer
//!
//! It was a `double` until February 2025 and is now an `int`, because
//! [`bayesian_estimate_of_empirical_quality`] had always returned a whole quality score. Every
//! getter widens it back to `double`, which is why a `%.2f` of it prints two zeros after the point
//! and never anything else.
//!
//! # And it is cached, with every setter invalidating it
//!
//! `getEmpiricalQuality(prior)` computes on the first call and returns the stored value on every
//! call after that, **whatever prior is asked for**. Measured:
//!
//! ```text
//! cache  prior-then-prior  first-with-10           13.0
//! cache  prior-then-prior  then-with-45            13.0
//! cache  prior-then-prior  after-setter-with-45    44.0
//! ```
//!
//! The setter in the middle stored the same observation count it already had. It invalidated the
//! cache all the same, and the second prior only then took effect. This is why the getters here
//! take `&mut self`: in Rust the mutation has to be in the signature, and it is the behaviour.
//!
//! # The smoothing is applied twice, with different arithmetic each time
//!
//! One error and two observations, added as doubles in [`RecalDatum::empirical_error_rate`] and
//! added to `(long)(mismatches + 0.5)` in `calc_empirical_quality`. That cast truncates after
//! adding a half, which is not the same as rounding: it is what makes half a mismatch count as one
//! and 0.49999 count as none, and both were measured.

use std::sync::LazyLock;

use crate::java_format::format_decimals;
use crate::math_utils::{max_element_index, qual_to_error_prob};

/// `SAMUtils.MAX_PHRED_SCORE`, the ceiling every empirical quality is clamped to.
pub const MAX_RECALIBRATED_Q_SCORE: i32 = 93;

/// The largest difference between an empirical and a reported quality the prior is tabulated for.
pub const MAX_GATK_USABLE_Q_SCORE: i32 = 40;

/// `QualityUtils.MAX_REASONABLE_Q_SCORE`, the top of the range the posterior is maximised over.
///
/// Twenty higher than [`MAX_GATK_USABLE_Q_SCORE`], which is the point worth noticing: the search
/// runs over 61 bins and the prior has only 41 entries, so every difference past 40 shares one
/// prior value.
pub const MAX_REASONABLE_Q_SCORE: i32 = 60;

/// See the module note. The constant, not the number the comment beside it names.
const MULTIPLIER: f64 = 100_000.0;

/// One error and one non-error observation, added so a datum with no errors is not quality
/// infinity.
const SMOOTHING_CONSTANT: i64 = 1;

/// The sentinel that means "not computed yet", which every setter writes back.
const UNINITIALIZED_EMPIRICAL_QUALITY: i32 = -1;

/// `Integer.MAX_VALUE - 1`, above which the binomial's counts are rescaled to fit an `int`.
const MAX_NUMBER_OF_OBSERVATIONS: i64 = i32::MAX as i64 - 1;

/// Every argument `RecalDatum` refuses, with the words the reference refuses it in.
#[derive(Debug, Clone, PartialEq)]
pub enum RecalDatumError {
    /// `numObservations < 0`, from the constructor or from `setNumObservations`.
    NegativeObservations,
    /// `numMismatches < 0`, from the constructor or from `setNumMismatches`.
    NegativeMismatches,
    /// `reportedQuality < 0`, from the constructor only.
    NegativeReportedQuality,
    /// `estimatedQReported < 0`, from `setReportedQuality`. A different wording for the same idea,
    /// because the field was called `estimatedQReported` when the check was written.
    NegativeEstimatedQReported,
    /// `estimatedQReported is infinite`.
    InfiniteEstimatedQReported,
    /// `estimatedQReported is NaN`.
    NaNEstimatedQReported,
    /// `empiricalQuality < 0`, from `setEmpiricalQuality`.
    NegativeEmpiricalQuality,
    /// `Utils.validateArg(qual >= 0.0)` inside `QualityUtils.qualToErrorProb(double)`, reached from
    /// `combine` when the reported quality is already NaN. See [`RecalDatum::combine`].
    QualityNotAtLeastZero(f64),
}

impl RecalDatumError {
    /// The exact `IllegalArgumentException` message, which the golden compares character for
    /// character.
    pub fn message(&self) -> String {
        match self {
            RecalDatumError::NegativeObservations => "numObservations < 0".to_string(),
            RecalDatumError::NegativeMismatches => "numMismatches < 0".to_string(),
            RecalDatumError::NegativeReportedQuality => "reportedQuality < 0".to_string(),
            RecalDatumError::NegativeEstimatedQReported => "estimatedQReported < 0".to_string(),
            RecalDatumError::InfiniteEstimatedQReported => {
                "estimatedQReported is infinite".to_string()
            }
            RecalDatumError::NaNEstimatedQReported => "estimatedQReported is NaN".to_string(),
            RecalDatumError::NegativeEmpiricalQuality => "empiricalQuality < 0".to_string(),
            // `Utils.validateArg` builds this with string concatenation, so the number is written
            // the way `String.valueOf(double)` writes it.
            RecalDatumError::QualityNotAtLeastZero(qual) => {
                format!("qual must be >= 0.0 but got {}", java_double(*qual))
            }
        }
    }
}

/// `String.valueOf(double)` for the handful of values the message above can carry.
///
/// Only NaN is reachable: the guard is `qual >= 0.0`, which is false for NaN and for every negative
/// number, and a reported quality can only become negative through the setter that refuses it.
fn java_double(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else if value == value.trunc() && value.abs() < 1e7 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

/// The 41 entries of `logPriorCache`: a Gaussian log density with mean zero and sigma one half,
/// evaluated at each integer difference between an empirical and a reported quality.
///
/// The reference builds it in a static initialiser from commons-math's `NormalDistribution`, whose
/// `logDensity` is
///
/// ```java
/// final double x0 = x - mean;
/// final double x1 = x0 / standardDeviation;
/// return -0.5 * x1 * x1 - logStandardDeviationPlusHalfLog2Pi;
/// ```
///
/// with the constant folded once in the constructor as
/// `FastMath.log(sd) + 0.5 * FastMath.log(2 * FastMath.PI)`.
///
/// **`FastMath.log` and not `Math.log`.** commons-math's logarithm is table-driven and is not
/// correctly rounded, so it disagrees with the platform's in the last bits, and the golden carries
/// the raw bit pattern of all 41 entries because that is the only way to see it.
///
/// **`FastMath.PI` and not `Math.PI`.** The two happen to be the same double, which the suite
/// asserts rather than assumes.
static LOG_PRIOR_CACHE: LazyLock<[f64; (MAX_GATK_USABLE_Q_SCORE + 1) as usize]> =
    LazyLock::new(|| {
        // The reference names both, so this does too, even though the mean is zero.
        let mean = 0.0;
        let sigma = 0.5;
        let log_sd_plus_half_log_2pi =
            jmath::fast_math::log(sigma) + 0.5 * jmath::fast_math::log(2.0 * std::f64::consts::PI);
        let mut cache = [0.0; (MAX_GATK_USABLE_Q_SCORE + 1) as usize];
        for (i, slot) in cache.iter_mut().enumerate() {
            let x0 = i as f64 - mean;
            let x1 = x0 / sigma;
            *slot = -0.5 * x1 * x1 - log_sd_plus_half_log_2pi;
        }
        cache
    });

/// The prior cache, for the suite that compares every entry against the reference's.
pub fn log_prior_cache() -> &'static [f64] {
    &*LOG_PRIOR_CACHE
}

/// `RecalDatum.getLogPrior(qualityScore, priorQualityScore)`.
///
/// **The cast runs before the absolute value.** The reference is
/// `Math.min(Math.abs((int) (qualityScore - priorQualityScore)), MAX_GATK_USABLE_Q_SCORE)`, so a
/// difference of `-0.5` becomes `(int) -0.5 = 0` and then `0`, where taking the absolute value
/// first would have given `1` and a different prior. Measured: a quality of 29.5 against a prior of
/// 30 gets the same value as 30 against 30.
///
/// **A NaN prior gives difference zero.** `(int)` of NaN is zero in Java, and a saturating `as`
/// cast is zero in Rust, so both index the cache at its first entry. That is reachable:
/// [`RecalDatum::combine`] can leave the reported quality NaN, and the empirical quality is then
/// computed against a flat prior.
pub fn get_log_prior(quality_score: f64, prior_quality_score: f64) -> f64 {
    let difference = ((quality_score - prior_quality_score) as i32)
        // `Math.abs(int)` overflows silently on `Integer.MIN_VALUE` and stays negative, so this is
        // the wrapping one rather than the panicking one. Reaching it needs an infinite prior,
        // which throws out of the cache lookup in the reference too.
        .wrapping_abs()
        .min(MAX_GATK_USABLE_Q_SCORE);
    LOG_PRIOR_CACHE[difference as usize]
}

/// `RecalDatum.getLogBinomialLikelihood(qualityScore, nObservations, nErrors)`.
///
/// Three escapes, and each one was measured:
///
///  * **no observations returns exactly `0.0`**, a log probability of one, before any distribution
///    is built;
///  * **counts above `Integer.MAX_VALUE - 1` are rescaled**, errors included, by a `Math.round` of
///    the scaled count, because the binomial's implementation caches on `int` arguments;
///  * **an infinite or NaN likelihood becomes `-Double.MAX_VALUE`** rather than `-Infinity`, which
///    keeps it comparable in the argmax that follows.
pub fn get_log_binomial_likelihood(
    quality_score: f64,
    mut n_observations: i64,
    mut n_errors: i64,
) -> f64 {
    if n_observations == 0 {
        return 0.0;
    }

    if n_observations > MAX_NUMBER_OF_OBSERVATIONS {
        let fraction = MAX_NUMBER_OF_OBSERVATIONS as f64 / n_observations as f64;
        n_errors = jmath::math::round(n_errors as f64 * fraction);
        n_observations = MAX_NUMBER_OF_OBSERVATIONS;
    }

    let log_likelihood = log_binomial_probability(
        n_observations as i32,
        n_errors as i32,
        qual_to_error_prob(quality_score),
    );
    if log_likelihood.is_infinite() || log_likelihood.is_nan() {
        -f64::MAX
    } else {
        log_likelihood
    }
}

/// `BinomialDistribution.logProbability(k)` from commons-math3, which is the one line
/// `MathUtils.logBinomialProbability` is.
///
/// The saddle point expansion under it is already ported, in `jmath`, for the Fisher strand
/// annotation. What is here is the wrapper: the two guards the distribution puts in front of it.
fn log_binomial_probability(n: i32, k: i32, p: f64) -> f64 {
    if n == 0 {
        return if k == 0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if k < 0 || k > n {
        f64::NEG_INFINITY
    } else {
        jmath::saddle_point::log_binomial_probability(k, n, p, 1.0 - p)
    }
}

/// `RecalDatum.bayesianEstimateOfEmpiricalQuality`.
///
/// The maximum a posteriori quality score under a Gaussian prior on the difference from the
/// reported quality and a binomial likelihood on the counts. The search is exhaustive over the 61
/// integer quality scores from zero to [`MAX_REASONABLE_Q_SCORE`], and the winner is the **first**
/// maximum, because `MathUtils.maxElementIndex` compares strictly greater-than.
pub fn bayesian_estimate_of_empirical_quality(
    n_observations: i64,
    n_errors: i64,
    prior_mean_quality_score: f64,
) -> i32 {
    let bins = (MAX_REASONABLE_Q_SCORE + 1) as usize;
    let log_posteriors: Vec<f64> = (0..bins)
        .map(|q| {
            get_log_prior(q as f64, prior_mean_quality_score)
                + get_log_binomial_likelihood(q as f64, n_observations, n_errors)
        })
        .collect();
    max_element_index(&log_posteriors, 0, log_posteriors.len()) as i32
}

/// The four-tuple, for one set of covariates.
///
/// `Clone` is the copy constructor: it copies the **raw** mismatch field and the cached empirical
/// quality, without dividing and multiplying, so a copy of a datum whose quality was already
/// computed carries that quality.
#[derive(Debug, Clone, PartialEq)]
pub struct RecalDatum {
    /// Scaled by nothing: a quality score as reported by the sequencer, or estimated from several
    /// of them by [`RecalDatum::combine`].
    reported_quality: f64,
    /// [`UNINITIALIZED_EMPIRICAL_QUALITY`] until computed. See the module note on caching.
    empirical_quality: i32,
    num_observations: i64,
    /// **Multiplied by [`MULTIPLIER`]**. Read it with [`RecalDatum::num_mismatches`].
    num_mismatches: f64,
}

impl RecalDatum {
    /// `new RecalDatum(numObservations, numMismatches, reportedQuality)`.
    ///
    /// The quality is a `byte` in the reference and an `i8` here, which matters: the check is
    /// `reportedQuality < 0`, and a quality of 200 written as a byte is -56 and refused, while a
    /// quality of 127 is accepted and is well above [`MAX_RECALIBRATED_Q_SCORE`].
    ///
    /// The three checks run in this order and the first one to fail decides the message.
    pub fn new(
        num_observations: i64,
        num_mismatches: f64,
        reported_quality: i8,
    ) -> Result<RecalDatum, RecalDatumError> {
        if num_observations < 0 {
            return Err(RecalDatumError::NegativeObservations);
        }
        // `< 0.0` and not `!(>= 0.0)`: a NaN mismatch count passes this and is stored.
        if num_mismatches < 0.0 {
            return Err(RecalDatumError::NegativeMismatches);
        }
        if reported_quality < 0 {
            return Err(RecalDatumError::NegativeReportedQuality);
        }
        Ok(RecalDatum {
            num_observations,
            num_mismatches: num_mismatches * MULTIPLIER,
            reported_quality: reported_quality as f64,
            empirical_quality: UNINITIALIZED_EMPIRICAL_QUALITY,
        })
    }

    /// `combine(other)`: add another datum's counts in and re-estimate the reported quality.
    ///
    /// **The reported quality is recomputed, not averaged.** Each side contributes the number of
    /// errors its own reported quality predicts, and the combined quality is whatever the total
    /// implies:
    ///
    /// ```java
    /// reportedQuality = -10 * Math.log10(expectedNumErrors / getNumObservations());
    /// ```
    ///
    /// Two consequences the golden carries. For two empty datums that is `-10*log10(0/0)`, so the
    /// field becomes **NaN**, and nothing stops it because the assignment is direct rather than
    /// through `setReportedQuality`, which refuses NaN. And combining onto that NaN then fails,
    /// because `QualityUtils.qualToErrorProb(double)` validates `qual >= 0.0` and NaN is not. That
    /// is the only reachable error here and it is why this returns a `Result`.
    ///
    /// The counts go through the public getter and the incrementer, so the mismatch count is
    /// divided by [`MULTIPLIER`] and multiplied by it again rather than copied across.
    ///
    /// **The sign of that NaN belongs to the processor.** `0.0 / 0.0` is the floating-point
    /// indefinite on x86-64, whose sign bit is set, and the default NaN on AArch64, whose sign bit
    /// is clear. Reference and port agree wherever they are run together and differ across
    /// architectures, so it is the one value in the golden compared as a NaN rather than as bits.
    pub fn combine(&mut self, other: &RecalDatum) -> Result<(), RecalDatumError> {
        let expected_num_errors = self.calc_expected_errors()? + other.calc_expected_errors()?;
        self.increment(other.num_observations(), other.num_mismatches());
        self.reported_quality =
            -10.0 * jmath::math::log10(expected_num_errors / self.num_observations() as f64);
        self.empirical_quality = UNINITIALIZED_EMPIRICAL_QUALITY;
        Ok(())
    }

    /// `setReportedQuality`, whose three checks name the field by the name it had when they were
    /// written.
    pub fn set_reported_quality(&mut self, reported_quality: f64) -> Result<(), RecalDatumError> {
        if reported_quality < 0.0 {
            return Err(RecalDatumError::NegativeEstimatedQReported);
        }
        if reported_quality.is_infinite() {
            return Err(RecalDatumError::InfiniteEstimatedQReported);
        }
        if reported_quality.is_nan() {
            return Err(RecalDatumError::NaNEstimatedQReported);
        }
        self.reported_quality = reported_quality;
        self.empirical_quality = UNINITIALIZED_EMPIRICAL_QUALITY;
        Ok(())
    }

    pub fn reported_quality(&self) -> f64 {
        self.reported_quality
    }

    /// `getReportedQualityAsByte()`, which is `(byte)(int)(Math.round(x))`.
    ///
    /// Two narrowings, and the outer one is where a quality comes back **negative**: a reported
    /// quality of 200 answers -56, and 127.5 rounds half-up to 128 and answers -128.
    pub fn reported_quality_as_byte(&self) -> i8 {
        jmath::math::round(self.reported_quality) as i32 as i8
    }

    /// `getEmpiricalErrorRate()`: errors over observations, smoothed, or exactly zero when there is
    /// nothing to divide.
    ///
    /// The smoothing is added to the counts as **integers** before the widening, which is the
    /// reference's `numObservations + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT` on a `long`.
    pub fn empirical_error_rate(&self) -> f64 {
        if self.num_observations == 0 {
            0.0
        } else {
            let double_mismatches = self.num_mismatches / MULTIPLIER + SMOOTHING_CONSTANT as f64;
            let double_observations =
                (self.num_observations + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT) as f64;
            double_mismatches / double_observations
        }
    }

    /// `setEmpiricalQuality(int)`: write the cache directly, so the counts no longer explain it.
    ///
    /// The reference also tests the argument for infinity and NaN. It is an `int`, so neither test
    /// can ever be true, and they are not ported.
    pub fn set_empirical_quality(&mut self, empirical_quality: i32) -> Result<(), RecalDatumError> {
        if empirical_quality < 0 {
            return Err(RecalDatumError::NegativeEmpiricalQuality);
        }
        self.empirical_quality = empirical_quality;
        Ok(())
    }

    /// `getEmpiricalQuality()`, with the reported quality as the prior.
    pub fn empirical_quality(&mut self) -> f64 {
        self.empirical_quality_with_prior(self.reported_quality)
    }

    /// `getEmpiricalQuality(priorQualityScore)`.
    ///
    /// **The prior is used only on the first call.** See the module note: this reads a cache that
    /// nothing but a setter clears, so a second call with a different prior gets the first answer
    /// back. Taking `&mut self` is what makes that visible in Rust.
    pub fn empirical_quality_with_prior(&mut self, prior_quality_score: f64) -> f64 {
        if self.empirical_quality == UNINITIALIZED_EMPIRICAL_QUALITY {
            self.calc_empirical_quality(prior_quality_score);
        }
        self.empirical_quality as f64
    }

    /// `getEmpiricalQualityAsByte()`, which is `(byte)(Math.round(...))` with **one** narrowing,
    /// unlike [`RecalDatum::reported_quality_as_byte`].
    pub fn empirical_quality_as_byte(&mut self) -> i8 {
        jmath::math::round(self.empirical_quality()) as i8
    }

    /// `toString()`: `"%d,%.2f,%.2f"` of observations, mismatches and empirical quality.
    ///
    /// It calls the empirical quality getter, so **printing a datum computes and caches it**. That
    /// is why this takes `&mut self` and why the golden prints each datum before asking for its
    /// quality separately.
    pub fn to_text(&mut self) -> String {
        let observations = self.num_observations;
        let mismatches = format_decimals(self.num_mismatches(), 2);
        let quality = format_decimals(self.empirical_quality(), 2);
        format!("{observations},{mismatches},{quality}")
    }

    /// `stringForCSV()`: the text above, then the reported quality and the difference from it.
    pub fn string_for_csv(&mut self) -> String {
        let text = self.to_text();
        let delta = self.empirical_quality() - self.reported_quality();
        format!(
            "{},{},{}",
            text,
            format_decimals(self.reported_quality(), 2),
            format_decimals(delta, 2)
        )
    }

    pub fn num_observations(&self) -> i64 {
        self.num_observations
    }

    pub fn set_num_observations(&mut self, num_observations: i64) -> Result<(), RecalDatumError> {
        if num_observations < 0 {
            return Err(RecalDatumError::NegativeObservations);
        }
        self.num_observations = num_observations;
        self.empirical_quality = UNINITIALIZED_EMPIRICAL_QUALITY;
        Ok(())
    }

    /// `getNumMismatches()`, which divides the stored value back down. See the module note.
    pub fn num_mismatches(&self) -> f64 {
        self.num_mismatches / MULTIPLIER
    }

    pub fn set_num_mismatches(&mut self, num_mismatches: f64) -> Result<(), RecalDatumError> {
        if num_mismatches < 0.0 {
            return Err(RecalDatumError::NegativeMismatches);
        }
        self.num_mismatches = num_mismatches * MULTIPLIER;
        self.empirical_quality = UNINITIALIZED_EMPIRICAL_QUALITY;
        Ok(())
    }

    /// `incrementNumObservations(by)`. Unchecked, so a negative `by` takes the count negative.
    pub fn increment_num_observations(&mut self, by: i64) {
        self.num_observations += by;
        self.empirical_quality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    /// `incrementNumMismatches(by)`. Unchecked, like the one above.
    pub fn increment_num_mismatches(&mut self, by: f64) {
        self.num_mismatches += by * MULTIPLIER;
        self.empirical_quality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    /// `increment(incObservations, incMismatches)`.
    ///
    /// **Neither argument is checked**, unlike every setter, so this is the way a datum's counts
    /// can go negative. Measured: incrementing a datum of one observation by -5 leaves -4.
    pub fn increment(&mut self, inc_observations: i64, inc_mismatches: f64) {
        self.num_observations += inc_observations;
        self.num_mismatches += inc_mismatches * MULTIPLIER;
        self.empirical_quality = UNINITIALIZED_EMPIRICAL_QUALITY;
    }

    /// `increment(isError)`: one observation, and one error or none.
    pub fn increment_by_observation(&mut self, is_error: bool) {
        self.increment(1, if is_error { 1.0 } else { 0.0 });
    }

    /// `calcExpectedErrors()`: how many errors the reported quality predicts over these
    /// observations.
    ///
    /// The guard is `QualityUtils.qualToErrorProb`'s, hoisted to the one call site that can fail
    /// it. `qualToErrorProb(double)` opens with `Utils.validateArg(qual >= 0.0)` on every call, and
    /// a reported quality is non-negative everywhere except after a [`RecalDatum::combine`] that
    /// left it NaN.
    fn calc_expected_errors(&self) -> Result<f64, RecalDatumError> {
        // The reference's condition is `qual >= 0.0` and this is its complement written out, which
        // is not `qual < 0.0`: NaN fails the reference's test and fails neither comparison.
        if self.reported_quality.is_nan() || self.reported_quality < 0.0 {
            return Err(RecalDatumError::QualityNotAtLeastZero(
                self.reported_quality,
            ));
        }
        Ok(self.num_observations() as f64 * qual_to_error_prob(self.reported_quality))
    }

    /// `calcEmpiricalQuality(priorQualityScore)`: the maximum a posteriori estimate, capped.
    ///
    /// The smoothing here is not the smoothing in [`RecalDatum::empirical_error_rate`]. The
    /// mismatch count goes through `(long)(x + 0.5)`, a **truncating cast after adding a half**,
    /// which the reference itself marks `TODO: why add 0.5?`. Half a mismatch becomes one error and
    /// 0.49999 becomes none, and the golden carries both.
    fn calc_empirical_quality(&mut self, prior_quality_score: f64) {
        let mismatches = (self.num_mismatches() + 0.5) as i64 + SMOOTHING_CONSTANT;
        let observations = self.num_observations + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT;
        let empirical_qual =
            bayesian_estimate_of_empirical_quality(observations, mismatches, prior_quality_score);
        self.empirical_quality = empirical_qual.min(MAX_RECALIBRATED_Q_SCORE);
    }
}

/// `EventType`: the three kinds of error a recalibration table counts, and therefore the three
/// tables a recalibration report holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EventType {
    BaseSubstitution,
    BaseInsertion,
    BaseDeletion,
}

impl EventType {
    /// `EventType.values()`, in declaration order, which is the ordinal order the tables are
    /// indexed by.
    pub const VALUES: [EventType; 3] = [
        EventType::BaseSubstitution,
        EventType::BaseInsertion,
        EventType::BaseDeletion,
    ];

    /// `eventFrom(int)`, which indexes the cached values array directly. Out of range is
    /// `ArrayIndexOutOfBoundsException` in the reference and `None` here.
    pub fn from_index(index: i32) -> Option<EventType> {
        usize::try_from(index).ok().and_then(|i| {
            if i < EventType::VALUES.len() {
                Some(EventType::VALUES[i])
            } else {
                None
            }
        })
    }

    /// `eventFrom(String)`, matched against the one-letter representation and not against the enum
    /// name. Unknown is `IllegalArgumentException("Event %s does not exist.")` in the reference.
    pub fn from_representation(representation: &str) -> Option<EventType> {
        EventType::VALUES
            .into_iter()
            .find(|event| event.representation() == representation)
    }

    /// `toString()`: the single letter a report's table name is built from.
    pub fn representation(&self) -> &'static str {
        match self {
            EventType::BaseSubstitution => "M",
            EventType::BaseInsertion => "I",
            EventType::BaseDeletion => "D",
        }
    }

    /// `prettyPrint()`: the words the report writes in its `EventType` column.
    pub fn pretty_print(&self) -> &'static str {
        match self {
            EventType::BaseSubstitution => "Base Substitution",
            EventType::BaseInsertion => "Base Insertion",
            EventType::BaseDeletion => "Base Deletion",
        }
    }

    /// `name()`, the enum constant's own name, which is not what `toString` returns.
    pub fn name(&self) -> &'static str {
        match self {
            EventType::BaseSubstitution => "BASE_SUBSTITUTION",
            EventType::BaseInsertion => "BASE_INSERTION",
            EventType::BaseDeletion => "BASE_DELETION",
        }
    }

    /// The ordinal, which is how a recalibration table addresses its third dimension.
    pub fn ordinal(&self) -> usize {
        match self {
            EventType::BaseSubstitution => 0,
            EventType::BaseInsertion => 1,
            EventType::BaseDeletion => 2,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_multiplier_is_the_constant_and_not_the_comment() {
        // A hundred thousand, and the comment beside it says ten thousand. Divide by the other one
        // and this value comes back ten times too large.
        let datum = RecalDatum::new(1000, 10.0, 30).unwrap();
        assert_eq!(datum.num_mismatches(), 10.0);
        assert_eq!(datum.num_mismatches * 1.0, 10.0 * 100_000.0);
    }

    #[test]
    fn the_cast_in_the_prior_runs_before_the_absolute_value() {
        // Half a point below the prior. Casting first gives zero; taking the absolute value first
        // would have given one, and a different entry of the cache.
        assert_eq!(get_log_prior(29.5, 30.0), get_log_prior(30.0, 30.0));
        assert_ne!(get_log_prior(29.0, 30.0), get_log_prior(30.0, 30.0));
    }

    #[test]
    fn a_nan_prior_is_a_flat_prior() {
        // `(int)` of NaN is zero in Java and a saturating `as` is zero in Rust, so every quality
        // score gets the same prior and the argmax is decided by the likelihood alone.
        assert_eq!(get_log_prior(0.0, f64::NAN), get_log_prior(60.0, f64::NAN));
    }

    #[test]
    fn the_prior_is_flat_past_forty_and_the_search_runs_to_sixty() {
        assert_eq!(get_log_prior(0.0, 40.0), get_log_prior(0.0, 60.0));
        assert_eq!(log_prior_cache().len(), 41);
    }

    #[test]
    fn no_observations_is_a_log_probability_of_one() {
        assert_eq!(get_log_binomial_likelihood(30.0, 0, 0), 0.0);
        assert_eq!(get_log_binomial_likelihood(0.0, 0, 0), 0.0);
    }

    #[test]
    fn an_impossible_outcome_is_the_most_negative_double_and_not_infinity() {
        // More errors than observations. The distribution answers -Infinity and the method
        // substitutes -Double.MAX_VALUE, so the argmax that follows can still compare it.
        assert_eq!(get_log_binomial_likelihood(30.0, 100, 200), -f64::MAX);
    }

    #[test]
    fn the_smoothing_in_the_quality_truncates_after_adding_a_half() {
        // Half a mismatch counts as one error; the double below a half counts as none. The two
        // datums differ by 0.00001 mismatches and by nothing else.
        let mut half = RecalDatum::new(1000, 0.5, 30).unwrap();
        let mut under = RecalDatum::new(1000, 0.49999, 30).unwrap();
        assert_eq!(half.empirical_quality(), under.empirical_quality());
        // The error rate, which smooths differently, does tell them apart.
        assert_ne!(half.empirical_error_rate(), under.empirical_error_rate());
    }

    #[test]
    fn the_empirical_quality_is_cached_and_a_setter_clears_it() {
        let mut datum = RecalDatum::new(1000, 10.0, 30).unwrap();
        let first = datum.empirical_quality_with_prior(10.0);
        // A different prior, and the cached answer comes back anyway.
        assert_eq!(datum.empirical_quality_with_prior(45.0), first);
        // The same observation count it already had, but a setter all the same.
        datum.set_num_observations(1000).unwrap();
        assert_ne!(datum.empirical_quality_with_prior(45.0), first);
    }

    #[test]
    fn combining_two_empty_datums_leaves_a_nan_that_the_setter_would_have_refused() {
        let mut empty = RecalDatum::new(0, 0.0, 30).unwrap();
        empty
            .combine(&RecalDatum::new(0, 0.0, 30).unwrap())
            .unwrap();
        assert!(empty.reported_quality().is_nan());
        assert_eq!(
            empty.set_reported_quality(f64::NAN),
            Err(RecalDatumError::NaNEstimatedQReported)
        );
        // And combining onto it fails, from the guard inside qualToErrorProb.
        let error = empty
            .combine(&RecalDatum::new(1000, 10.0, 30).unwrap())
            .unwrap_err();
        assert_eq!(error.message(), "qual must be >= 0.0 but got NaN");
    }

    #[test]
    fn increment_is_the_one_unchecked_way_in() {
        let mut datum = RecalDatum::new(1, 0.0, 30).unwrap();
        datum.increment(-5, -5.0);
        assert_eq!(datum.num_observations(), -4);
        assert_eq!(datum.num_mismatches(), -5.0);
        // Where every setter refuses the same thing.
        assert_eq!(
            datum.set_num_observations(-1),
            Err(RecalDatumError::NegativeObservations)
        );
    }

    #[test]
    fn the_two_byte_getters_narrow_differently() {
        let mut datum = RecalDatum::new(1, 0.0, 30).unwrap();
        datum.set_reported_quality(200.0).unwrap();
        assert_eq!(datum.reported_quality_as_byte(), -56);
        datum.set_reported_quality(127.5).unwrap();
        assert_eq!(datum.reported_quality_as_byte(), -128);
        datum.set_empirical_quality(200).unwrap();
        assert_eq!(datum.empirical_quality_as_byte(), -56);
    }

    #[test]
    fn the_event_letter_is_not_the_enum_name() {
        assert_eq!(EventType::BaseSubstitution.representation(), "M");
        assert_eq!(EventType::BaseSubstitution.name(), "BASE_SUBSTITUTION");
        assert_eq!(
            EventType::from_representation("M"),
            Some(EventType::BaseSubstitution)
        );
        assert_eq!(EventType::from_representation("BASE_SUBSTITUTION"), None);
        assert_eq!(EventType::from_index(3), None);
        assert_eq!(EventType::from_index(-1), None);
    }
}
