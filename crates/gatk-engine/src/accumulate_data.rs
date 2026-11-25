//! `Mutect2FilteringEngine.accumulateData` and the pass schedule `FilterMutectCalls` drives it
//! with, ported from `org.broadinstitute.hellbender.tools.walkers.mutect.filtering` (GATK 4.7.0.0;
//! `Mutect2FilteringEngine.java` and `FilterMutectCalls.java` are unchanged since 4.6.2.0).
//!
//! # Four passes, not two
//!
//! ```java
//! private static final int NUMBER_OF_LEARNING_PASSES = 2;
//! protected int numberOfPasses() { return NUMBER_OF_LEARNING_PASSES + 2; }
//! ```
//!
//! Passes 0, 1 and 2 all accumulate; only 0 and 1 learn parameters afterwards. Pass 2 accumulates a
//! whole traversal's worth of data and uses it **only** to relearn the threshold, the filters'
//! parameters being deliberately frozen: "it's important for filter parameters to stay the same and
//! only learn the threshold in the final pass so that the final threshold used corresponds exactly
//! to the filters". Pass 3 applies the filters and writes.
//!
//! # A record whose only alternate is `<NON_REF>` is skipped
//!
//! ```java
//! if (vc.getAlleles().stream().noneMatch(a -> a.isNonReference() && !a.isNonRefAllele())) {
//!     return;
//! }
//! ```
//!
//! `isNonReference()` is true of `<NON_REF>` as well, so the guard reads "no alternate other than
//! `<NON_REF>`". A GVCF-mode site contributes to neither the clustering model nor the threshold,
//! and a record with no alternate at all is skipped by the same test. A symbolic alternate that is
//! **not** `<NON_REF>` is not skipped, and still accumulates nothing: the symbolic removal inside
//! `ErrorProbabilities` has already emptied every list. Two records reach the same three zeroes by
//! different routes.
//!
//! # The three accumulators do not move together
//!
//! `record` drops an alternate whose artifact or non-somatic probability is an obvious artifact,
//! counting the first case, while the threshold calculator is fed the combined probabilities
//! regardless. A triallelic record whose two alternates are both obvious artifacts accumulates no
//! clustering data, two probabilities and two obvious artifacts.
//!
//! # A record with no `TLOD` is a refusal
//!
//! The reference dereferences `tumorLogOdds.length` inside `record`, which is a
//! `NullPointerException` rather than a skip: a VCF missing that annotation crashes a non-final
//! pass. One difference the golden cannot see: by the time it throws, `record` has already zeroed
//! the symbolic entries of the caller's depth array. This port refuses before calling `record`, so
//! the caller's array is untouched.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing" (`Some(x)` or `None`), Java's nullable reference with the
//!   null check enforced by the compiler. The pass schedules answer `None` where the reference
//!   throws `ShouldNeverReachHereException`.
//! - `Result<T, E>` is "a value or an error" (`Ok(x)` or `Err(e)`): how this codebase models a
//!   Java exception. `?` after a call returns its error to the caller at once, Java's `throws` at
//!   one call site.
//! - `&mut T` is a borrowed, writable reference: the function may change what the caller passed,
//!   as Java may write through an array or a list it is handed. `&[T]` is the read-only form.

use crate::somatic_clustering_model::{AlternateAllele, RecordError, SomaticClusteringModel};

/// `FilterMutectCalls.NUMBER_OF_LEARNING_PASSES`.
///
/// `i32` is Java's `int`, kept so the pass number `n` compares with it as it does in the reference.
pub const NUMBER_OF_LEARNING_PASSES: i32 = 2;

/// `numberOfPasses()`: the learning passes, one for the threshold, and one for calling.
pub const NUMBER_OF_PASSES: i32 = NUMBER_OF_LEARNING_PASSES + 2;

/// What `nthPassApply` does in pass `n`.
///
/// The reference does the work inside the walker callback; here the schedule is returned as a
/// value and the driver acts on it, so the schedule can be tested apart from the traversal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PassAction {
    /// `filteringEngine.accumulateData(...)`, for `n <= NUMBER_OF_LEARNING_PASSES`.
    Accumulate,
    /// `vcfWriter.add(filteringEngine.applyFiltersAndAccumulateOutputStats(...))`.
    ApplyAndWrite,
}

/// What `afterNthPass` does after pass `n`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AfterPassAction {
    /// `learnParameters()`, for `n < NUMBER_OF_LEARNING_PASSES`.
    LearnParameters,
    /// `learnThreshold()` alone, at `n == NUMBER_OF_LEARNING_PASSES`: the filters are frozen.
    LearnThresholdOnly,
    /// `writeFilteringStats(...)`.
    WriteFilteringStats,
}

/// `nthPassApply`'s schedule. `None` is the `ShouldNeverReachHereException`.
///
/// ```java
/// if (n <= NUMBER_OF_LEARNING_PASSES) {
///     filteringEngine.accumulateData(variant, referenceContext);
/// } else if (n == NUMBER_OF_LEARNING_PASSES + 1) {
/// ```
///
/// The `pass < 0` arm has no counterpart in the reference, whose `n <= NUMBER_OF_LEARNING_PASSES`
/// would accept a negative `n`. The walker counts its passes from zero, so the arm only matters to
/// a driver that miscounts, which it then refuses rather than accumulates.
pub fn action_for_pass(pass: i32) -> Option<PassAction> {
    // `if ... else if ...` is an expression in Rust: the branch taken is the function's value, so
    // no `return` is written.
    if pass < 0 {
        None
    } else if pass <= NUMBER_OF_LEARNING_PASSES {
        Some(PassAction::Accumulate)
    } else if pass == NUMBER_OF_LEARNING_PASSES + 1 {
        Some(PassAction::ApplyAndWrite)
    } else {
        None
    }
}

/// `afterNthPass`'s schedule. `None` is the `ShouldNeverReachHereException`.
///
/// Note the `<` here against the `<=` in [`action_for_pass`]: pass 2 accumulates but only relearns
/// the threshold afterwards, which is the module note's "four passes, not two".
pub fn action_after_pass(pass: i32) -> Option<AfterPassAction> {
    if pass < 0 {
        None
    } else if pass < NUMBER_OF_LEARNING_PASSES {
        Some(AfterPassAction::LearnParameters)
    } else if pass == NUMBER_OF_LEARNING_PASSES {
        Some(AfterPassAction::LearnThresholdOnly)
    } else if pass == NUMBER_OF_LEARNING_PASSES + 1 {
        Some(AfterPassAction::WriteFilteringStats)
    } else {
        None
    }
}

/// One alternate allele, as the accumulation guard and the clustering model each see it.
///
/// The reference reads both facts off one htsjdk `Allele`; here the caller extracts them, because
/// `allele.symbolic` (what `record` skips) and `non_ref` (what the guard skips) are different sets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AccumulationAllele {
    pub allele: AlternateAllele,
    /// `Allele.isNonRefAllele()`, which is true of `<NON_REF>` alone and not of every symbolic
    /// allele.
    pub non_ref: bool,
}

/// What the accumulation refuses: each variant stands for the Java exception named on it.
#[derive(Debug, Clone, PartialEq)]
pub enum AccumulateError {
    /// The `NullPointerException` on `tumorLogOdds.length` when the record has no `TLOD`.
    NoTumorLogOdds,
    /// What `record` itself refuses.
    Record(RecordError),
}

impl AccumulateError {
    /// The fully qualified Java exception class the reference throws.
    pub fn class(&self) -> &'static str {
        match self {
            AccumulateError::NoTumorLogOdds => "java.lang.NullPointerException",
            AccumulateError::Record(_) => "java.lang.IllegalArgumentException",
        }
    }

    /// The exception's message, as the JVM renders it. The first is the JDK's helpful
    /// `NullPointerException` text; the second is the `Utils.validateArg` message in
    /// `SomaticClusteringModel.record`.
    pub fn message(&self) -> String {
        match self {
            AccumulateError::NoTumorLogOdds => {
                "Cannot read the array length because \"tumorLogOdds\" is null".to_string()
            }
            AccumulateError::Record(_) => {
                "tumorADs must have one entry per allele including the ref allele".to_string()
            }
        }
    }
}

/// Lets `?` turn the clustering model's refusal into this type's `Record` variant, so a refusal
/// inside `record` propagates as it does in the reference.
impl From<RecordError> for AccumulateError {
    fn from(error: RecordError) -> Self {
        AccumulateError::Record(error)
    }
}

/// `accumulateData`, without the filters' own accumulation, which each filter owns.
///
/// `combined_probabilities` is `errorProbabilities.getCombinedErrorProbabilities()`; the artifact
/// and non-somatic lists are the other two the model reads. `tumor_allele_depths` is written
/// through, as the reference writes through the array it is handed.
///
/// Returns `false` when the guard skipped the record without touching anything.
///
/// The order is the reference's (`Mutect2FilteringEngine.accumulateData`): guard, then
/// `somaticClusteringModel.record(...)`, then `thresholdCalculator.addCombinedErrorProbabilites`.
/// So a record that `record` refuses adds nothing to the threshold calculator either.
#[allow(clippy::too_many_arguments)]
pub fn accumulate_data(
    model: &mut SomaticClusteringModel,
    accumulated_error_probabilities: &mut Vec<f64>,
    alternates: &[AccumulationAllele],
    tumor_allele_depths: &mut [i32],
    tumor_log_odds: Option<&[f64]>,
    artifact_probabilities: &[f64],
    non_somatic_probabilities: &[f64],
    combined_probabilities: &[f64],
    reference_length: i32,
) -> Result<bool, AccumulateError> {
    // `noneMatch(a -> a.isNonReference() && !a.isNonRefAllele())`: no alternate other than
    // `<NON_REF>`, which a record with no alternate at all also satisfies.
    if !alternates.iter().any(|alternate| !alternate.non_ref) {
        return Ok(false);
    }

    // `let Some(x) = y else { ... }` unwraps the `Option` or takes the early exit; there is no Java
    // equivalent. Here the exit is the `NullPointerException` the module note describes.
    let Some(tumor_log_odds) = tumor_log_odds else {
        return Err(AccumulateError::NoTumorLogOdds);
    };

    // `.iter().map(..).collect()` is Java's `stream().map(..).collect(toList())`; `|a| a.allele` is
    // a lambda, `a -> a.allele`.
    let alleles: Vec<AlternateAllele> = alternates.iter().map(|a| a.allele).collect();
    model.record(
        tumor_allele_depths,
        tumor_log_odds,
        artifact_probabilities,
        non_somatic_probabilities,
        &alleles,
        reference_length,
    )?;
    // `thresholdCalculator.addCombinedErrorProbabilites(...)`, spelled as the reference spells it.
    accumulated_error_probabilities.extend_from_slice(combined_probabilities);
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The schedule, which is the thing a driver is most likely to get wrong.
    #[test]
    fn three_passes_accumulate_and_only_two_of_them_learn() {
        assert_eq!(NUMBER_OF_PASSES, 4);
        for pass in 0..=2 {
            assert_eq!(
                action_for_pass(pass),
                Some(PassAction::Accumulate),
                "pass {pass}"
            );
        }
        assert_eq!(action_for_pass(3), Some(PassAction::ApplyAndWrite));
        assert_eq!(action_for_pass(4), None, "past the last pass");

        assert_eq!(action_after_pass(0), Some(AfterPassAction::LearnParameters));
        assert_eq!(action_after_pass(1), Some(AfterPassAction::LearnParameters));
        // The third accumulating pass learns the threshold ALONE: the filters stay as they were.
        assert_eq!(
            action_after_pass(2),
            Some(AfterPassAction::LearnThresholdOnly)
        );
        assert_eq!(
            action_after_pass(3),
            Some(AfterPassAction::WriteFilteringStats)
        );
        assert_eq!(action_after_pass(4), None);
    }
}
