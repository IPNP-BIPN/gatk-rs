//! Ported from `org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine` and
//! the `Dirichlet` under it (GATK 4.6.2.0).
//!
//! The variational fixed point `AllelePseudoDepth` ends in, and the place G1.9's open question
//! lives. `NaturalLogUtils` (G1.9.1) established that one ulp is what enters here; what the
//! iteration does with it is what this slice measures.
//!
//! # Two risks, and they are different
//!
//! ```java
//! double[] dirichletPosterior = new IndexRange(0, numberOfAlleles).mapToDouble(n -> 1.0);
//! boolean converged = false;
//! while (!converged) {
//!     final double[] alleleCounts = getEffectiveCounts(logLikelihoods, dirichletPosterior, weights);
//!     final double[] newDirichletPosterior = MathArrays.ebeAdd(alleleCounts, priorPseudocounts);
//!     converged = MathArrays.distance1(dirichletPosterior, newDirichletPosterior)
//!                     / MathUtils.sum(newDirichletPosterior) < CONVERGENCE_THRESHOLD;
//!     dirichletPosterior = newDirichletPosterior;
//! }
//! ```
//!
//! **Amplification.** A 1-ulp difference in one iteration's `exp` feeds the next iteration's
//! input. Whether it grows, shrinks or stays put is a property of this particular map, not
//! something that can be asserted from the outside.
//!
//! **The iteration count.** Convergence is a *threshold* test on
//! `distance1 / sum < 0.001`. A difference far too small to see in the values can still put one
//! side of that comparison on the other side of the threshold, and then the two runs do a
//! different number of iterations. That is a much larger divergence than one ulp, and it is why
//! [`allele_fractions_posterior`] returns the count alongside the result: a divergence has to be
//! attributable to one cause or the other, not merely observed.
//!
//! # `while (!converged)` tests **after** the step
//!
//! The loop is not `do { } while`, but it behaves as one: `converged` is false on entry, so the
//! body always runs at least once, and the test uses the value the body just computed. So the
//! minimum count is one, never zero, even for an input that was already at the fixed point.
//!
//! # The summation order is the answer
//!
//! `sumArrayFunction` accumulates over reads in **index order**, adding each read's contribution
//! element-wise into a running array. Floating-point addition is not associative: summing the same
//! reads in a different order gives a different `double`. The reference's order is transcribed
//! rather than rewritten, and the same goes for `MathUtils.sum`, which the convergence test
//! divides by.

use jmath::gamma;

use crate::natural_log_utils::{posteriors, NonFiniteSum};

/// `SomaticLikelihoodsEngine.CONVERGENCE_THRESHOLD`.
pub const CONVERGENCE_THRESHOLD: f64 = 0.001;

/// What a fixed-point run produced, and how long it took to get there.
///
/// The count is not diagnostics. It is the observable that separates "the arithmetic drifted by a
/// last bit" from "the two runs did different work", and G1.9 cannot be settled without it.
#[derive(Debug, Clone, PartialEq)]
pub struct Posterior {
    /// `dirichletPosterior` at the point the loop stopped.
    pub values: Vec<f64>,
    /// How many times the body ran. Never zero: the loop tests after the step.
    pub iterations: usize,
}

/// `Dirichlet.effectiveLogMultinomialWeights()`.
///
/// ```java
/// final double digammaOfSum = Gamma.digamma(MathUtils.sum(alpha));
/// return MathUtils.applyToArray(alpha, a -> (Gamma.digamma(a) - digammaOfSum));
/// ```
///
/// `digamma` is commons-math3's, ported and oracle-backed in htsjdk-rs (#62), so this is
/// composition rather than new numerics — and unlike everything downstream of it, it carries no
/// `exp` and therefore no 1-ulp bound. `sum` accumulates in index order.
pub fn effective_log_multinomial_weights(alpha: &[f64]) -> Option<Vec<f64>> {
    let digamma_of_sum = gamma::digamma(sum(alpha)).ok()?;
    alpha
        .iter()
        .map(|a| gamma::digamma(*a).ok().map(|d| d - digamma_of_sum))
        .collect()
}

/// `SomaticLikelihoodsEngine.logDirichletNormalization(double...)`.
///
/// ```java
/// final double logNumerator = Gamma.logGamma(MathUtils.sum(dirichletParams));
/// final double logDenominator = MathUtils.sum(MathUtils.applyToArray(dirichletParams, Gamma::logGamma));
/// return logNumerator - logDenominator;
/// ```
///
/// Four of these make the allele-fraction clusters' correction, so the two sides of the subtraction
/// matter more than the value does.
///
/// # A single parameter of one is negative zero
///
/// `logGamma(1.0)` is `-0.0`, and the denominator's sum starts from `0.0` and stays there, so the
/// answer is `-0.0 - 0.0`. At `0.5` both sides are the same positive number and the answer is `0.0`.
/// The two differ only in a sign bit no `==` can see.
///
/// # A parameter of zero is NaN, not an infinity
///
/// commons-math's `logGamma` answers `NaN` at and below zero rather than diverging, so a zero
/// parameter poisons the denominator and the subtraction carries the `NaN` out.
pub fn log_dirichlet_normalization(dirichlet_params: &[f64]) -> f64 {
    let log_numerator = gamma::log_gamma(sum(dirichlet_params));
    let logs: Vec<f64> = dirichlet_params
        .iter()
        .map(|p| gamma::log_gamma(*p))
        .collect();
    log_numerator - sum(&logs)
}

/// `MathUtils.sum(double[])`, accumulating in index order.
pub fn sum(values: &[f64]) -> f64 {
    let mut total = 0.0;
    for value in values {
        total += value;
    }
    total
}

/// `MathArrays.distance1`: the sum of absolute differences, again in index order.
pub fn distance1(first: &[f64], second: &[f64]) -> f64 {
    let mut total = 0.0;
    for (a, b) in first.iter().zip(second) {
        total += (a - b).abs();
    }
    total
}

/// `getEffectiveCounts(logLikelihoods, dirichletPrior, weights)`.
///
/// `log_likelihoods` is `[allele][read]`, matching the reference's `RealMatrix` whose rows are
/// alleles and whose columns are reads: `getColumn(read)` is one read's likelihood per allele.
///
/// The accumulation is `sumArrayFunction` over reads in index order, and the per-read vector is
/// `NaturalLogUtils.posteriors`, which is where the 1-ulp `exp` enters.
pub fn effective_counts(
    log_likelihoods: &[Vec<f64>],
    dirichlet_prior: &[f64],
    weights: Option<&[f64]>,
) -> Option<Vec<f64>> {
    let effective_log_weights = effective_log_multinomial_weights(dirichlet_prior)?;
    let read_count = log_likelihoods.first().map_or(0, Vec::len);
    if read_count == 0 {
        return None;
    }

    let mut total: Option<Vec<f64>> = None;
    for read in 0..read_count {
        let column: Vec<f64> = log_likelihoods.iter().map(|row| row[read]).collect();
        let mut unweighted = posteriors(&effective_log_weights, &column)?;
        if let Some(weights) = weights {
            // `applyToArrayInPlace(unweighted, d -> d * weights[read])`.
            for value in unweighted.iter_mut() {
                *value *= weights[read];
            }
        }
        total = Some(match total {
            None => unweighted,
            Some(mut running) => {
                for (slot, value) in running.iter_mut().zip(&unweighted) {
                    *slot += value;
                }
                running
            }
        });
    }
    total
}

/// `alleleFractionsPosterior(logLikelihoods, priorPseudocounts, weights)`.
///
/// The initial posterior is flat at `1.0` per allele, whatever the prior is, and the prior only
/// enters through the addition inside the loop.
///
/// There is no iteration cap in the reference, and none is added here: adding one would make the
/// port terminate where the reference does not, which is a divergence dressed as prudence. A
/// non-converging input would hang both.
pub fn allele_fractions_posterior(
    log_likelihoods: &[Vec<f64>],
    prior_pseudocounts: &[f64],
    weights: Option<&[f64]>,
) -> Result<Posterior, NonFiniteSum> {
    let allele_count = log_likelihoods.len();
    assert_eq!(
        allele_count,
        prior_pseudocounts.len(),
        "Must have one pseudocount per allele."
    );

    let mut dirichlet_posterior = vec![1.0; allele_count];
    let mut iterations = 0usize;
    loop {
        let Some(allele_counts) = effective_counts(log_likelihoods, &dirichlet_posterior, weights)
        else {
            return Err(NonFiniteSum);
        };
        let new_posterior: Vec<f64> = allele_counts
            .iter()
            .zip(prior_pseudocounts)
            .map(|(count, prior)| count + prior)
            .collect();
        iterations += 1;

        // The test is on the step just taken, so the body always runs at least once.
        let converged = distance1(&dirichlet_posterior, &new_posterior) / sum(&new_posterior)
            < CONVERGENCE_THRESHOLD;
        dirichlet_posterior = new_posterior;
        if converged {
            return Ok(Posterior {
                values: dirichlet_posterior,
                iterations,
            });
        }
    }
}
