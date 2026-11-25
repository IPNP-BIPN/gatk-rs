//! The two `AlleleFractionCluster` implementations, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.clustering` (GATK 4.7.0.0;
//! `AlleleFractionCluster.java`, `BinomialCluster.java`, `BetaBinomialCluster.java`, `Datum.java`
//! and `readorientation/BetaDistributionShape.java` are unchanged since 4.6.2.0), with the
//! `BetaDistributionShape` and `Datum` they carry.
//!
//! `SomaticClusteringModel` keeps a `BinomialCluster` per allele fraction it tracks and a
//! `BetaBinomialCluster` beside them, and asks each for two numbers.
//!
//! # `correctedLogLikelihood` is not a likelihood
//!
//! It is the record's TLOD corrected for the cluster's shape:
//!
//! ```java
//! return datum.getTumorLogOdds() + logOddsCorrection(BetaDistributionShape.FLAT_BETA, betaDistributionShape, altCount, refCount);
//! ```
//!
//! and the correction is four Dirichlet normalisations, two of them at the flat shape. The flat
//! terms are constants and they do not cancel in doubles, so a port that dropped them would agree
//! with the mathematics and disagree with the reference.
//!
//! # `BinomialCluster` has no binomial in it
//!
//! Its constructor turns a mean into a beta shape at a fixed standard-deviation-over-mean of
//! `0.01`, and clamps the mean at `1 - 0.01` on the way. A cluster asked for a mean of `1.0` and
//! one asked for `2.0` come out at the same shape as one asked for `0.99`, and only a mean of `0.5`
//! is symmetric. `alphaPlusBeta` is `((1 - mean) / (mean * 0.0001)) - 1`, which is nearly ten
//! million at a mean of `0.001`: these are the shapes that drive [`jmath::beta::log_beta`] into its
//! `a >= 10` branch.
//!
//! Both `learn` methods are ported here, as `learn_binomial` and `learn_beta_binomial`. They
//! return the new shape rather than overwriting a field, and `SomaticClusteringModel` stores it.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<BetaDistributionShape, ShapeError>` is "a shape or the reference's
//!   `IllegalArgumentException`"; `Err(..)` is the throw, `Ok(..)` the return. A `?` after a call
//!   passes an `Err` straight up to the caller, as an uncaught exception would.
//! - `enum AlleleFractionCluster { Binomial(..), BetaBinomial(..) }` stands for the Java interface
//!   and its two implementations: a value is exactly one of the two variants and carries its
//!   shape, and a `match` on it must handle both.
//! - `#[derive(Clone, Copy, PartialEq)]` asks the compiler to generate copying and `==`; a `Copy`
//!   type is passed by value, as a Java `double` is.
//! - `x as f64` on an `i32` is Java's implicit `int` to `double` widening, which is exact.
//! - `data.iter().enumerate()` yields `(index, item)` pairs, the `n` and `data.get(n)` of the
//!   reference's counted loop.

use crate::beta_binomial::{BetaBinomialDistribution, BetaBinomialError};
use crate::somatic_likelihoods::log_dirichlet_normalization;

/// `readorientation.BetaDistributionShape`, whose two arguments are refused by two different call
/// paths. The fields are private and read through `alpha()` and `beta()`, the reference's
/// `getAlpha()` and `getBeta()`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BetaDistributionShape {
    alpha: f64,
    beta: f64,
}

/// The two refusals, kept apart because the reference words them differently.
#[derive(Debug, Clone, PartialEq)]
pub enum ShapeError {
    /// `ParamUtils.isPositive`: `alpha must be greater than 0 but got <alpha>`.
    Alpha { alpha: f64 },
    /// `Utils.validateArg`: `beta must be greater than 0 but got <beta>`.
    Beta { beta: f64 },
}

impl BetaDistributionShape {
    /// `FLAT_BETA`, the shape every correction is measured against.
    pub const FLAT: BetaDistributionShape = BetaDistributionShape {
        alpha: 1.0,
        beta: 1.0,
    };

    /// `new BetaDistributionShape(alpha, beta)`:
    ///
    /// ```java
    /// ParamUtils.isPositive(alpha, "alpha must be greater than 0 but got " + alpha);
    /// Utils.validateArg(beta > 0, "beta must be greater than 0 but got " + beta);
    /// ```
    ///
    /// `ParamUtils.isPositive` is `Utils.validateArg(val > 0, message)`, so both are the same test
    /// and the same exception class. Alpha is checked first, so a shape bad in both reports alpha.
    pub fn new(alpha: f64, beta: f64) -> Result<Self, ShapeError> {
        // Both checks are `> 0`, so both refuse NaN.
        if alpha.is_nan() || alpha <= 0.0 {
            return Err(ShapeError::Alpha { alpha });
        }
        if beta.is_nan() || beta <= 0.0 {
            return Err(ShapeError::Beta { beta });
        }
        Ok(Self { alpha, beta })
    }

    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    pub fn beta(&self) -> f64 {
        self.beta
    }
}

/// `Datum`, whose non-sequencing-error probability is computed once in the constructor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Datum {
    tumor_log_odds: f64,
    artifact_prob: f64,
    non_sequencing_error_prob: f64,
    alt_count: i32,
    total_count: i32,
    indel_length: i32,
}

impl Datum {
    /// `new Datum(tumorLogOdds, artifactProb, nonSomaticProb, altCount, totalCount, indelLength)`.
    /// `nonSomaticProb` is not kept: only the combined probability below is.
    pub fn new(
        tumor_log_odds: f64,
        artifact_prob: f64,
        non_somatic_prob: f64,
        alt_count: i32,
        total_count: i32,
        indel_length: i32,
    ) -> Self {
        Self {
            tumor_log_odds,
            artifact_prob,
            // `1 - (1 - artifactProb) * (1 - nonSomaticProb)`, which is not either input: an
            // artifact probability of 0.3 alone comes back as 0.30000000000000004.
            non_sequencing_error_prob: 1.0 - (1.0 - artifact_prob) * (1.0 - non_somatic_prob),
            alt_count,
            total_count,
            indel_length,
        }
    }

    pub fn tumor_log_odds(&self) -> f64 {
        self.tumor_log_odds
    }

    pub fn artifact_prob(&self) -> f64 {
        self.artifact_prob
    }

    pub fn non_sequencing_error_prob(&self) -> f64 {
        self.non_sequencing_error_prob
    }

    pub fn alt_count(&self) -> i32 {
        self.alt_count
    }

    pub fn total_count(&self) -> i32 {
        self.total_count
    }

    pub fn indel_length(&self) -> i32 {
        self.indel_length
    }
}

/// `BinomialCluster`'s fixed standard deviation over mean.
const STD_DEV_OVER_MEAN: f64 = 0.01;

/// `BinomialCluster.getFuzzyBinomial(unboundedMean)`.
///
/// The clamp is `Math.min(unboundedMean, 1 - 0.01)`, which propagates NaN in Java and does not in
/// Rust, so it is written out rather than taken from `f64::min`.
///
/// The clamp is pinned by the `allele_fraction_cluster.txt.gz` golden: its `fuzzy 0.99`,
/// `fuzzy 1.0` and `fuzzy 2.0` rows are all `99.01000000000009,1.0001010101010053`.
pub fn fuzzy_binomial(unbounded_mean: f64) -> Result<BetaDistributionShape, ShapeError> {
    let bound = 1.0 - STD_DEV_OVER_MEAN;
    // `Math.min` answers NaN when either argument is NaN; Rust's `f64::min` answers the other one.
    let mean = if unbounded_mean.is_nan() {
        f64::NAN
    } else {
        unbounded_mean.min(bound)
    };
    // `mean * MathUtils.square(STD_DEV_OVER_MEAN)`: the square is formed first, and
    // `(mean * 0.01) * 0.01` is a different double.
    let alpha_plus_beta = ((1.0 - mean) / (mean * (STD_DEV_OVER_MEAN * STD_DEV_OVER_MEAN))) - 1.0;
    let alpha = mean * alpha_plus_beta;
    BetaDistributionShape::new(alpha, alpha_plus_beta - alpha)
}

/// `BetaBinomialCluster.logOddsCorrection(originalBeta, newBeta, altCount, refCount)`.
///
/// Four normalisations, in the reference's order and with its signs. `g(1, 1)` is `0.0` and
/// `g(1 + alt, 1 + ref)` is not, so the flat pair does not vanish.
///
/// The reference's `return` (`BetaBinomialCluster.java:64-65`) is `g(new) - g(new + counts)`,
/// then `- g(original) + g(original + counts)`. The chain is evaluated left to right in both
/// languages, `((a - b) - c) + d`, so the rounding of each partial result is the reference's.
/// `let g = log_dirichlet_normalization;` names the function, as the reference's private
/// `g(double...)` does, and `&[x, y]` is the two-element array the varargs call builds.
pub fn log_odds_correction(
    original: BetaDistributionShape,
    new: BetaDistributionShape,
    alt_count: i32,
    ref_count: i32,
) -> f64 {
    let g = log_dirichlet_normalization;
    g(&[new.alpha(), new.beta()])
        - g(&[
            new.alpha() + alt_count as f64,
            new.beta() + ref_count as f64,
        ])
        - g(&[original.alpha(), original.beta()])
        + g(&[
            original.alpha() + alt_count as f64,
            original.beta() + ref_count as f64,
        ])
}

/// `BetaBinomialCluster.correctedLogLikelihood(datum, betaDistributionShape)`, which both clusters
/// use: `BinomialCluster` delegates to it with the shape its mean produced.
pub fn corrected_log_likelihood(datum: &Datum, shape: BetaDistributionShape) -> f64 {
    let alt_count = datum.alt_count();
    let ref_count = datum.total_count() - alt_count;
    datum.tumor_log_odds()
        + log_odds_correction(BetaDistributionShape::FLAT, shape, alt_count, ref_count)
}

/// `logLikelihood(totalCount, altCount)`, the same for both clusters once the shape is fixed:
/// `new BetaBinomialDistribution(null, alpha, beta, totalCount).logProbability(altCount)`, whose
/// `null` generator is never drawn from. The `?` passes the distribution's refusal through.
pub fn log_likelihood(
    shape: BetaDistributionShape,
    total_count: i32,
    alt_count: i32,
) -> Result<f64, BetaBinomialError> {
    BetaBinomialDistribution::new(shape.alpha(), shape.beta(), total_count)?
        .log_probability(alt_count)
}

/// `BetaBinomialCluster.RATE`, the learning rate. The reference's `MAX_RATE = 0.1` is declared
/// and never read, so it is not ported.
const RATE: f64 = 0.01;
/// `BetaBinomialCluster.NUM_EPOCHS`, the passes over the data per call to `learn`.
const NUM_EPOCHS: usize = 10;

/// `BetaBinomialCluster.learn(data, responsibilities)`.
///
/// Ten epochs of gradient ascent over `Gamma.digamma`, and three things about it that a smoother
/// implementation would get wrong:
///
///  * **the update is sequential within an epoch.** Alpha and beta are read and written per datum,
///    so the same data in a different order learn a different shape;
///  * **the floors are applied every step**, not at the end, so a shape that would have dipped
///    below `1.0` or `0.5` and come back cannot;
///  * **`Math.max` propagates NaN** where Rust's `f64::max` answers the other argument, so the
///    clamps are written out.
pub fn learn_beta_binomial(
    shape: BetaDistributionShape,
    data: &[Datum],
    responsibilities: &[f64],
) -> Result<BetaDistributionShape, ShapeError> {
    let mut alpha = shape.alpha();
    let mut beta = shape.beta();
    for _ in 0..NUM_EPOCHS {
        for (index, datum) in data.iter().enumerate() {
            let alt = datum.alt_count() as f64;
            let reference = (datum.total_count() - datum.alt_count()) as f64;
            // `datum.getTotalCount() + alpha + beta` in the reference: the `int` is widened and
            // the sum runs left to right, `(total + alpha) + beta`, as it does here.
            let digamma_of_total_plus_alpha_plus_beta =
                digamma(datum.total_count() as f64 + alpha + beta);
            let digamma_of_alpha_plus_beta = digamma(alpha + beta);
            let alpha_gradient =
                digamma(alpha + alt) - digamma_of_total_plus_alpha_plus_beta - digamma(alpha)
                    + digamma_of_alpha_plus_beta;
            let beta_gradient =
                digamma(beta + reference) - digamma_of_total_plus_alpha_plus_beta - digamma(beta)
                    + digamma_of_alpha_plus_beta;
            alpha = java_max(alpha + RATE * alpha_gradient * responsibilities[index], 1.0);
            beta = java_max(beta + RATE * beta_gradient * responsibilities[index], 0.5);
        }
    }
    BetaDistributionShape::new(alpha, beta)
}

/// `Gamma.digamma`, which jmath carries and which refuses nothing this reaches. jmath's version
/// returns a `Result`; `unwrap_or(f64::NAN)` takes the value, or `NaN` on a refusal.
fn digamma(x: f64) -> f64 {
    jmath::gamma::digamma(x).unwrap_or(f64::NAN)
}

/// `Math.max`, which answers NaN when either argument is NaN.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

/// `BinomialCluster.learn(data, responsibilities)`.
///
/// A responsibility-weighted mean with `0.0001` added to **both** sums. That is not a tie-breaker:
/// it is what keeps a cluster given no responsibility at all from dividing zero by zero. Such a
/// cluster lands at `0.0001 / 0.0001`, which the fuzzy binomial's own clamp turns into `0.99`.
///
/// Both sums are `DoubleStream.sum`, which is compensated rather than a plain loop
/// (`BinomialCluster.java:32-33`, `IntStream.range(..).mapToDouble(..).sum() + 0.0001`). The port
/// collects the products into a `Vec` first and sums that, in index order, as the stream does.
pub fn learn_binomial(
    data: &[Datum],
    responsibilities: &[f64],
) -> Result<BetaDistributionShape, ShapeError> {
    let alt_count = double_stream_sum(
        &data
            .iter()
            .enumerate()
            .map(|(index, datum)| datum.alt_count() as f64 * responsibilities[index])
            .collect::<Vec<f64>>(),
    ) + 0.0001;
    let total_count = double_stream_sum(
        &data
            .iter()
            .enumerate()
            .map(|(index, datum)| datum.total_count() as f64 * responsibilities[index])
            .collect::<Vec<f64>>(),
    ) + 0.0001;
    fuzzy_binomial(alt_count / total_count)
}

/// `DoubleStream.sum()`, which is Kahan summation with a simple sum kept beside it.
///
/// ```java
/// double tmp = summands[0] - summands[1];
/// double simpleSum = summands[summands.length - 1];
/// if (Double.isNaN(tmp) && Double.isInfinite(simpleSum)) { return simpleSum; } else { return tmp; }
/// ```
///
/// A plain loop is a different double on almost any list of more than two values, which is why this
/// is not `MathUtils.sum`.
pub fn double_stream_sum(values: &[f64]) -> f64 {
    let mut sum = 0.0;
    let mut compensation = 0.0;
    let mut simple = 0.0;
    for value in values {
        let tmp = value - compensation;
        let velvel = sum + tmp;
        compensation = (velvel - sum) - tmp;
        sum = velvel;
        simple += value;
    }
    // The compensation is SUBTRACTED. JDK 17 as first released added it ("better error bounds to
    // add both terms as the final sum"), which is JDK-8214761; the oracle's 17.0.19 carries the
    // fix, and a probe of twenty thousand random lists in the pinned image found 65 on which the
    // two finishes differ, every one of them agreeing with the subtraction.
    let total = sum - compensation;
    if total.is_nan() && simple.is_infinite() {
        simple
    } else {
        total
    }
}

/// The `AlleleFractionCluster` interface and its two implementations, which differ only in where
/// their shape comes from. `toString` is not ported; `learn` is the two free functions above.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlleleFractionCluster {
    /// `BinomialCluster`, built from a mean.
    Binomial(BetaDistributionShape),
    /// `BetaBinomialCluster`, built from a shape directly.
    BetaBinomial(BetaDistributionShape),
}

impl AlleleFractionCluster {
    /// `new BinomialCluster(mean)`.
    pub fn binomial(mean: f64) -> Result<Self, ShapeError> {
        Ok(AlleleFractionCluster::Binomial(fuzzy_binomial(mean)?))
    }

    /// `new BetaBinomialCluster(shape)`.
    pub fn beta_binomial(shape: BetaDistributionShape) -> Self {
        AlleleFractionCluster::BetaBinomial(shape)
    }

    /// The shape either variant carries. The `|` in the `match` arm binds `shape` from whichever
    /// variant it is, since both carry the same type.
    pub fn shape(&self) -> BetaDistributionShape {
        match self {
            AlleleFractionCluster::Binomial(shape) | AlleleFractionCluster::BetaBinomial(shape) => {
                *shape
            }
        }
    }

    /// `correctedLogLikelihood(datum)`. Both implementations call the static
    /// `BetaBinomialCluster.correctedLogLikelihood` with their own shape.
    pub fn corrected_log_likelihood(&self, datum: &Datum) -> f64 {
        corrected_log_likelihood(datum, self.shape())
    }

    /// `logLikelihood(totalCount, altCount)`, identical in both implementations.
    pub fn log_likelihood(
        &self,
        total_count: i32,
        alt_count: i32,
    ) -> Result<f64, BetaBinomialError> {
        log_likelihood(self.shape(), total_count, alt_count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_mean_is_clamped_before_the_shape_is_built() {
        let capped = fuzzy_binomial(0.99).expect("a shape");
        for mean in [1.0, 2.0, f64::INFINITY] {
            assert_eq!(fuzzy_binomial(mean).expect("a shape"), capped, "{mean}");
        }
        assert_eq!(capped.alpha(), 99.01000000000009);
        assert_eq!(capped.beta(), 1.0001010101010053);
        // Only the half-way mean is symmetric.
        let half = fuzzy_binomial(0.5).expect("a shape");
        assert_eq!(half.alpha(), 4999.5);
        assert_eq!(half.beta(), 4999.5);
    }

    #[test]
    fn the_flat_correction_does_not_vanish() {
        // Both shapes flat: every term cancels and the answer is the TLOD alone.
        let datum = Datum::new(0.0, 0.0, 0.0, 5, 10, 0);
        assert_eq!(
            corrected_log_likelihood(&datum, BetaDistributionShape::FLAT),
            0.0
        );
        // And the TLOD is carried through whatever it is.
        for odds in [5.0, -5.0, f64::NEG_INFINITY] {
            let moved = Datum::new(odds, 0.0, 0.0, 5, 10, 0);
            assert_eq!(
                corrected_log_likelihood(&moved, BetaDistributionShape::FLAT),
                odds
            );
        }
        let not_a_number = Datum::new(f64::NAN, 0.0, 0.0, 5, 10, 0);
        assert!(corrected_log_likelihood(&not_a_number, BetaDistributionShape::FLAT).is_nan());
    }

    #[test]
    fn the_datum_combines_its_two_probabilities() {
        assert_eq!(
            Datum::new(0.0, 0.3, 0.0, 5, 10, 0).non_sequencing_error_prob(),
            0.30000000000000004
        );
        assert_eq!(
            Datum::new(0.0, 0.3, 0.3, 5, 10, 0).non_sequencing_error_prob(),
            0.51
        );
        assert_eq!(
            Datum::new(0.0, 1.0, 0.5, 5, 10, 0).non_sequencing_error_prob(),
            1.0
        );
    }

    #[test]
    fn the_shape_refuses_its_two_arguments_separately() {
        assert_eq!(
            BetaDistributionShape::new(0.0, 1.0),
            Err(ShapeError::Alpha { alpha: 0.0 })
        );
        assert_eq!(
            BetaDistributionShape::new(1.0, 0.0),
            Err(ShapeError::Beta { beta: 0.0 })
        );
        assert!(matches!(
            BetaDistributionShape::new(1.0, f64::NAN),
            Err(ShapeError::Beta { beta }) if beta.is_nan()
        ));
    }
}
