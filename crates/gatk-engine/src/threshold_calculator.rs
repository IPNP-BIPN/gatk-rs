//! `ThresholdCalculator`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ThresholdCalculator`
//! (GATK 4.6.2.0).
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
//! index and returns `0`. Since `relearnThresholdAndClearAcumulatedProbabilities` **clears** what it
//! learned from, relearning twice runs the second pass over nothing, and a tool that does so either
//! passes everything or filters everything depending on which strategy it was given.

/// `ThresholdCalculator.Strategy`.
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

/// What `ParamUtils` refuses, with the caller's own wording.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ThresholdError {
    NegativeBeta,
    NegativeFalsePositiveRate,
}

impl ThresholdError {
    pub fn class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }

    pub fn message(&self) -> &'static str {
        match self {
            ThresholdError::NegativeBeta => "requested F-score beta must be non-negative",
            ThresholdError::NegativeFalsePositiveRate => "requested FPR must be non-negative",
        }
    }
}

/// The calculator itself: a strategy, a threshold, and whatever has been accumulated.
#[derive(Debug, Clone, PartialEq)]
pub struct ThresholdCalculator {
    strategy: Strategy,
    threshold: f64,
    max_false_discovery_rate: f64,
    f_score_beta: f64,
    error_probabilities: Vec<f64>,
}

impl ThresholdCalculator {
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

    /// `addCombinedErrorProbabilites`.
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

    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// What has been accumulated and not yet learned from.
    pub fn accumulated(&self) -> &[f64] {
        &self.error_probabilities
    }
}

/// `calculateThresholdBasedOnOptimalFScore`, whose argument is sorted in place.
pub fn threshold_from_optimal_f_score(
    posteriors: &mut [f64],
    beta: f64,
) -> Result<f64, ThresholdError> {
    if beta < 0.0 {
        return Err(ThresholdError::NegativeBeta);
    }
    posteriors.sort_by(f64::total_cmp);

    let expected_true_positives: f64 = posteriors.iter().map(|prob| 1.0 - prob).sum();

    let mut true_positives = 0.0;
    let mut false_positives = 0.0;
    let mut false_negatives = expected_true_positives;
    // -1 means filter everything, which is what an empty list leaves behind.
    let mut optimal_index: Option<usize> = None;
    let mut optimal_f_score = 0.0;

    let count = posteriors.len();
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

    Ok(match optimal_index {
        None => 0.0,
        Some(index) if index == count - 1 => 1.0,
        Some(index) => posteriors[index],
    })
}

/// `calculateThresholdBasedOnFalseDiscoveryRate`, whose argument is sorted in place.
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
        let expected = (cumulative_expected_false_positives + posterior) / (index as f64 + 1.0);
        if expected > requested_false_positive_rate {
            // Step back one: the last posterior that kept the rate acceptable, or filter everything.
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
