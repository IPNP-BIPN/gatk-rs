//! `ThresholdCalculator`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ThresholdCalculator`
//! (GATK 4.7.0.0; `ThresholdCalculator.java` is unchanged since 4.6.2.0).
//!
//! The piece of `FilterMutectCalls` that decides where the error probability is cut: a list of
//! per-variant posteriors goes in, and one threshold comes out. A variant whose combined error
//! probability is at or below the threshold is called; above it, it is filtered.
//!
//! # The optimal F score keeps the last tie
//!
//! ```java
//! if (F >= optimalFScore) {
//!     optimalIndexInclusive = n;
//!     optimalFScore = F;
//! }
//! ```
//!
//! `>=` rather than `>`, so a run of equally good cut points resolves to the **largest** of them.
//! Four identical posteriors therefore end on the last index, and the answer is not that posterior
//! at all:
//!
//! ```java
//! return optimalIndexInclusive == -1 ? 0 : (optimalIndexInclusive == N - 1 ? 1 : posteriors.get(optimalIndexInclusive));
//! ```
//!
//! Three answers, two of them constants. The last index means `1`, no index at all means `0`, and
//! only the middle case reports a posterior.
//!
//! # An empty list is answered differently by the two strategies
//!
//! The false-discovery walk never exceeds anything and falls out at `1.0`; the F score finds no
//! index and returns `0`. Since `relearnThresholdAndClearAcumulatedProbabilities` **clears** what
//! it learned from, relearning twice runs the second pass over nothing, and a tool that does so
//! either passes everything or filters everything depending on which strategy it was given.
//!
//! The conformance golden is `crates/gatk-engine/tests/data/threshold_calculator.txt.gz`, written
//! by `tools/readfilter-conformance/ThresholdCalculatorDump.java` over six fixed lists.
//!
//! # Reading this file without reading Rust
//!
//! - `&mut [f64]` is a writable borrowed view of a list of doubles: the function may reorder it in
//!   place, as the reference's `Collections.sort(posteriors)` reorders the caller's `List`, but
//!   it does not own it.
//! - `Result<f64, ThresholdError>` is "a threshold or a refusal", how this codebase models a
//!   method that throws; `?` after a call rethrows its refusal to the caller.
//! - `Option<usize>` is "an index or nothing", used where the reference keeps `-1` as a sentinel.
//!   `usize` is the unsigned index type and cannot hold `-1`.

/// `ThresholdCalculator.Strategy`: `CONSTANT, FALSE_DISCOVERY_RATE, OPTIMAL_F_SCORE`.
///
/// The reference's `switch` has a `default` that throws `ShouldNeverReachHereException`; the
/// port's `match` names all three, and the compiler refuses a fourth left out, so it has none.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strategy {
    Constant,
    FalseDiscoveryRate,
    OptimalFScore,
}

impl Strategy {
    /// The enum constant's own name, which is what an unexpected-strategy message would quote.
    pub fn name(&self) -> &'static str {
        match self {
            Strategy::Constant => "CONSTANT",
            Strategy::FalseDiscoveryRate => "FALSE_DISCOVERY_RATE",
            Strategy::OptimalFScore => "OPTIMAL_F_SCORE",
        }
    }
}

/// What `ParamUtils.isPositiveOrZero` refuses, with the caller's own wording.
///
/// The check is `Utils.validateArg(val >= 0, message)`, so the reference refuses a `NaN` beta or
/// rate. The port tests `beta < 0.0` and `rate < 0.0`, which a `NaN` passes: the two differ there.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ThresholdError {
    /// `calculateThresholdBasedOnOptimalFScore` given a negative beta.
    NegativeBeta,
    /// `calculateThresholdBasedOnFalseDiscoveryRate` given a negative rate.
    NegativeFalsePositiveRate,
}

impl ThresholdError {
    /// The exception class: `Utils.validateArg` throws `IllegalArgumentException`.
    pub fn class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }

    /// The message the caller passed to `isPositiveOrZero`, verbatim.
    pub fn message(&self) -> &'static str {
        match self {
            ThresholdError::NegativeBeta => "requested F-score beta must be non-negative",
            ThresholdError::NegativeFalsePositiveRate => "requested FPR must be non-negative",
        }
    }
}

/// `ThresholdCalculator`: a strategy, a threshold, and whatever has been accumulated.
///
/// `error_probabilities` is the reference's package-private `List<Double> errorProbabilities`.
#[derive(Debug, Clone, PartialEq)]
pub struct ThresholdCalculator {
    strategy: Strategy,
    threshold: f64,
    max_false_discovery_rate: f64,
    f_score_beta: f64,
    error_probabilities: Vec<f64>,
}

impl ThresholdCalculator {
    /// `new ThresholdCalculator(strategy, initialThreshold, maxFalseDiscoveryRate, fScoreBeta)`.
    /// The initial threshold stands until the first relearn, and for `Constant` for good.
    pub fn new(
        strategy: Strategy,
        initial_threshold: f64,
        max_false_discovery_rate: f64,
        f_score_beta: f64,
    ) -> Self {
        ThresholdCalculator {
            strategy,
            threshold: initial_threshold,
            max_false_discovery_rate,
            f_score_beta,
            error_probabilities: Vec::new(),
        }
    }

    /// `addCombinedErrorProbabilites` (the reference's spelling): append, as `addAll` does.
    pub fn add(&mut self, error_probabilities: &[f64]) {
        self.error_probabilities
            .extend_from_slice(error_probabilities);
    }

    /// `relearnThresholdAndClearAcumulatedProbabilities`.
    ///
    /// The clearing is not a detail: a second call learns from nothing, which the two strategies
    /// answer at opposite ends.
    pub fn relearn(&mut self) -> Result<(), ThresholdError> {
        match self.strategy {
            // Don't adjust.
            Strategy::Constant => {}
            Strategy::FalseDiscoveryRate => {
                self.threshold = threshold_from_false_discovery_rate(
                    &mut self.error_probabilities,
                    self.max_false_discovery_rate,
                )?;
            }
            Strategy::OptimalFScore => {
                self.threshold = threshold_from_optimal_f_score(
                    &mut self.error_probabilities,
                    self.f_score_beta,
                )?;
            }
        }
        self.error_probabilities.clear();
        Ok(())
    }

    /// `getThreshold()`.
    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// What has been accumulated and not yet learned from.
    pub fn accumulated(&self) -> &[f64] {
        &self.error_probabilities
    }
}

/// `calculateThresholdBasedOnOptimalFScore`, whose argument is sorted in place.
///
/// **What**: the posterior at which to cut so that the expected F-beta score of what passes is
/// largest, where a posterior `p` counts `1 - p` of a true positive and `p` of a false one.
///
/// **How**: sort ascending, then walk the cut point up one variant at a time, keeping the
/// running expected true positives, false positives and false negatives, and the best cut.
///
/// **Why the arithmetic is spelled out**: `F` is the reference's expression term for term and in
/// its order of additions, `(1+beta*beta)*truePositives.getValue() /
/// ((1+beta*beta)*truePositives.getValue() + beta*beta*falseNegatives.getValue() +
/// falsePositives.getValue())`. With the `>=` tie rule, the last bit of `F` can decide which of
/// two cut points wins.
pub fn threshold_from_optimal_f_score(
    posteriors: &mut [f64],
    beta: f64,
) -> Result<f64, ThresholdError> {
    if beta < 0.0 {
        return Err(ThresholdError::NegativeBeta);
    }
    // `Collections.sort` on `Double`s is `Double.compare`. `total_cmp` is the same order except
    // that it puts a `NaN` with its sign bit set first rather than last.
    posteriors.sort_by(f64::total_cmp);

    // The reference's `posteriors.stream().mapToDouble(prob -> 1 - prob).sum()`. `.iter()`,
    // `.map(..)` and `.sum()` are the same stream; this `sum` adds left to right.
    let expected_true_positives: f64 = posteriors.iter().map(|prob| 1.0 - prob).sum();

    let mut true_positives = 0.0;
    let mut false_positives = 0.0;
    let mut false_negatives = expected_true_positives;
    // -1 means filter everything, which is what an empty list leaves behind.
    let mut optimal_index: Option<usize> = None;
    let mut optimal_f_score = 0.0;

    let count = posteriors.len();
    // `.enumerate()` pairs each element with its index: Java's `for (int n = 0; n < N; n++)`
    // with `posteriors.get(n)` read inside. The three `MutableDouble`s are plain `f64`s here:
    // `add` and `subtract` are `+=` and `-=`.
    for (index, posterior) in posteriors.iter().enumerate() {
        true_positives += 1.0 - posterior;
        false_positives += posterior;
        false_negatives -= 1.0 - posterior;
        let f = (1.0 + beta * beta) * true_positives
            / ((1.0 + beta * beta) * true_positives
                + beta * beta * false_negatives
                + false_positives);
        // `>=`: the last of a run of equally good cut points wins.
        if f >= optimal_f_score {
            optimal_index = Some(index);
            optimal_f_score = f;
        }
    }

    // The reference's nested ternary. The arm `Some(index) if index == count - 1` matches only
    // when there is an index, so `count - 1` is never taken of an empty list.
    Ok(match optimal_index {
        None => 0.0,
        Some(index) if index == count - 1 => 1.0,
        Some(index) => posteriors[index],
    })
}

/// `calculateThresholdBasedOnFalseDiscoveryRate`, whose argument is sorted in place.
///
/// **What**: the largest posterior such that the expected false-discovery rate among the variants
/// at or below it does not exceed the requested rate.
///
/// **How**: sort ascending and walk up, where the expected rate of the first `i + 1` is their mean
/// posterior. The first index that overshoots answers the posterior one below it, or `0.0`
/// (`thresholdForFilteringAll`) at the first index; a walk that never overshoots answers `1.0`
/// (`thresholdForFilteringNone`).
///
/// `(cumulativeExpectedFPs + posterior) / (i + 1)` is kept as one division of the running sum, not
/// a running mean, and the sum is updated only after the test, as the reference does.
pub fn threshold_from_false_discovery_rate(
    posteriors: &mut [f64],
    requested_false_positive_rate: f64,
) -> Result<f64, ThresholdError> {
    if requested_false_positive_rate < 0.0 {
        return Err(ThresholdError::NegativeFalsePositiveRate);
    }
    posteriors.sort_by(f64::total_cmp);

    let mut cumulative_expected_false_positives = 0.0;
    for (index, posterior) in posteriors.iter().enumerate() {
        // Java's `i + 1` is an `int` promoted to `double`; `index as f64 + 1.0` is the same value
        // for any index a list can have.
        let expected = (cumulative_expected_false_positives + posterior) / (index as f64 + 1.0);
        if expected > requested_false_positive_rate {
            // Step back one: the last posterior that kept the rate acceptable, or filter
            // everything.
            return Ok(if index > 0 {
                posteriors[index - 1]
            } else {
                0.0
            });
        }
        cumulative_expected_false_positives += posterior;
    }
    // Never exceeded: let everything pass.
    Ok(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spread() -> Vec<f64> {
        vec![0.9, 0.1, 0.5, 0.02, 0.3]
    }

    #[test]
    fn the_list_is_sorted_in_place() {
        let mut posteriors = vec![0.9, 0.1, 0.5];
        threshold_from_false_discovery_rate(&mut posteriors, 0.05).expect("a rate");
        assert_eq!(posteriors, vec![0.1, 0.5, 0.9]);
    }

    #[test]
    fn the_optimal_f_score_keeps_the_last_tie() {
        // Four equal posteriors: every cut point scores the same, and the last one wins, which the
        // three-way return then reports as 1.0 rather than as 0.2.
        let mut tied = vec![0.2, 0.2, 0.2, 0.2];
        assert_eq!(
            threshold_from_optimal_f_score(&mut tied, 1.0).expect("beta"),
            1.0
        );
    }

    #[test]
    fn the_optimal_f_score_answers_three_ways() {
        // No index at all: an empty list filters everything.
        assert_eq!(
            threshold_from_optimal_f_score(&mut [], 1.0).expect("beta"),
            0.0
        );
        // The last index: everything passes.
        let mut clean = vec![0.001, 0.002, 0.003];
        assert_eq!(
            threshold_from_optimal_f_score(&mut clean, 1.0).expect("beta"),
            1.0
        );
        // Anything else: the posterior at that index.
        let mut hopeless = vec![0.99, 0.98, 0.97];
        assert_eq!(
            threshold_from_optimal_f_score(&mut hopeless, 1.0).expect("beta"),
            0.97
        );
        let mut posteriors = spread();
        assert_eq!(
            threshold_from_optimal_f_score(&mut posteriors, 1.0).expect("beta"),
            0.5
        );
    }

    #[test]
    fn beta_weighs_recall_against_precision() {
        // Zero is precision alone, and cuts hard.
        let mut posteriors = spread();
        assert_eq!(
            threshold_from_optimal_f_score(&mut posteriors, 0.0).expect("beta"),
            0.02
        );
        // Ten is recall almost alone, and keeps everything.
        let mut posteriors = spread();
        assert_eq!(
            threshold_from_optimal_f_score(&mut posteriors, 10.0).expect("beta"),
            1.0
        );
    }

    #[test]
    fn the_false_discovery_rate_steps_back_one() {
        let mut posteriors = spread();
        // Sorted: 0.02, 0.1, 0.3, 0.5, 0.9. The rate is first exceeded at 0.1, so 0.02 is the
        // threshold.
        assert_eq!(
            threshold_from_false_discovery_rate(&mut posteriors, 0.05).expect("a rate"),
            0.02
        );
        // Exceeded at the very first posterior: filter everything.
        let mut tied = vec![0.2, 0.2, 0.2, 0.2];
        assert_eq!(
            threshold_from_false_discovery_rate(&mut tied, 0.05).expect("a rate"),
            0.0
        );
        // Never exceeded: pass everything.
        let mut clean = vec![0.001, 0.002, 0.003];
        assert_eq!(
            threshold_from_false_discovery_rate(&mut clean, 0.05).expect("a rate"),
            1.0
        );
    }

    #[test]
    fn an_empty_list_is_answered_at_opposite_ends() {
        assert_eq!(
            threshold_from_false_discovery_rate(&mut [], 0.05).expect("a rate"),
            1.0
        );
        assert_eq!(
            threshold_from_optimal_f_score(&mut [], 1.0).expect("beta"),
            0.0
        );
    }

    #[test]
    fn relearning_clears_what_it_learned_from() {
        for (strategy, second) in [
            (Strategy::FalseDiscoveryRate, 1.0),
            (Strategy::OptimalFScore, 0.0),
            // Constant never looks at any of it.
            (Strategy::Constant, 0.123),
        ] {
            let mut calculator = ThresholdCalculator::new(strategy, 0.123, 0.05, 1.0);
            calculator.add(&spread());
            calculator.relearn().expect("valid parameters");
            assert!(calculator.accumulated().is_empty());
            calculator.relearn().expect("nothing left");
            assert_eq!(calculator.threshold(), second, "{}", strategy.name());
        }
    }

    #[test]
    fn a_negative_beta_or_rate_is_refused() {
        assert_eq!(
            threshold_from_optimal_f_score(&mut spread(), -1.0)
                .expect_err("negative")
                .message(),
            "requested F-score beta must be non-negative"
        );
        assert_eq!(
            threshold_from_false_discovery_rate(&mut spread(), -0.5)
                .expect_err("negative")
                .message(),
            "requested FPR must be non-negative"
        );
        assert_eq!(
            ThresholdError::NegativeBeta.class(),
            "java.lang.IllegalArgumentException"
        );
    }
}
