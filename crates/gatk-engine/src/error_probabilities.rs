//! `ErrorProbabilities`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ErrorProbabilities`
//! (GATK 4.6.2.0).
//!
//! Every filter's answer, combined into one probability per alternate allele.
//!
//! # A filter that answers an empty list is dropped
//!
//! ```java
//! .entrySet().stream().filter(entry -> !entry.getValue().isEmpty())
//! ```
//!
//! A switched-off filter is not a filter that answered zero: it is removed from the map before
//! anything is combined, and before `transpose` checks that the surviving lists agree in length.
//!
//! # One type is a maximum, the types are independent
//!
//! Within an error type the **maximum** over filters wins for each allele, so failing two filters of
//! a type is no worse than failing one. Across types the answer is `1 - prod(1 - p)`, and the result
//! goes through [`crate::mutect_engine::round_finite_precision_errors`].
//!
//! # Ragged lists are a refusal
//!
//! `transpose` validates that every list is the same size. A site-level filter copies its one answer
//! to every alternate allele, which is what keeps them equal in the ordinary case; a record whose
//! annotation is shorter than its allele list breaks it.

use crate::mutect_engine::round_finite_precision_errors;

/// `ErrorType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ErrorType {
    Artifact,
    NonSomatic,
}

/// What `transpose` refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RaggedLists;

impl RaggedLists {
    pub fn class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }

    pub fn message(&self) -> &'static str {
        "lists are not the same size"
    }
}

/// One filter's answer: its error type and one probability per alternate allele.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterAnswer {
    pub error_type: ErrorType,
    pub probabilities: Vec<f64>,
}

/// `transpose(list)`, which refuses lists of differing length.
pub fn transpose(lists: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, RaggedLists> {
    let Some(first) = lists.first() else {
        return Ok(Vec::new());
    };
    if lists.iter().any(|list| list.len() != first.len()) {
        return Err(RaggedLists);
    }
    Ok((0..first.len())
        .map(|index| lists.iter().map(|list| list[index]).collect())
        .collect())
}

/// The probability of each error type, per allele: the maximum over that type's filters.
///
/// The empty answers have already been dropped by [`kept`], which is why nothing here treats an
/// absent filter as a zero.
pub fn by_type(answers: &[FilterAnswer], error_type: ErrorType) -> Result<Vec<f64>, RaggedLists> {
    let lists: Vec<Vec<f64>> = answers
        .iter()
        .filter(|answer| answer.error_type == error_type)
        .map(|answer| answer.probabilities.clone())
        .collect();
    Ok(transpose(&lists)?
        .into_iter()
        .map(|per_allele| {
            per_allele.into_iter().fold(
                0.0f64,
                |worst, value| if value > worst { value } else { worst },
            )
        })
        .collect())
}

/// The filters that were applied at all: the ones whose answer is not empty.
pub fn kept(answers: &[FilterAnswer]) -> Vec<FilterAnswer> {
    answers
        .iter()
        .filter(|answer| !answer.probabilities.is_empty())
        .cloned()
        .collect()
}

/// `getCombinedErrorProbabilities`: the whole pipeline, from every filter's answer to one
/// probability per allele.
pub fn combined(answers: &[FilterAnswer]) -> Result<Vec<f64>, RaggedLists> {
    let applied = kept(answers);
    let mut per_type: Vec<Vec<f64>> = Vec::new();
    for error_type in [ErrorType::Artifact, ErrorType::NonSomatic] {
        let probabilities = by_type(&applied, error_type)?;
        if !probabilities.is_empty() {
            per_type.push(probabilities);
        }
    }
    Ok(transpose(&per_type)?
        .into_iter()
        .map(|per_allele| {
            // Independent: the probability that none of the types is an error, complemented.
            let not_an_error = per_allele
                .into_iter()
                .fold(1.0, |product, probability| product * (1.0 - probability));
            round_finite_precision_errors(1.0 - not_an_error)
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn artifact(probabilities: &[f64]) -> FilterAnswer {
        FilterAnswer {
            error_type: ErrorType::Artifact,
            probabilities: probabilities.to_vec(),
        }
    }

    fn non_somatic(probabilities: &[f64]) -> FilterAnswer {
        FilterAnswer {
            error_type: ErrorType::NonSomatic,
            probabilities: probabilities.to_vec(),
        }
    }

    #[test]
    fn a_filter_that_answers_an_empty_list_is_dropped() {
        let answers = [artifact(&[0.0]), artifact(&[]), artifact(&[1.0])];
        assert_eq!(kept(&answers).len(), 2);
        // And the empty one does not make the lists ragged, having been removed first.
        assert_eq!(combined(&answers).expect("two of three"), vec![1.0]);
    }

    #[test]
    fn one_type_is_a_maximum() {
        // Two filters failing is no worse than one.
        assert_eq!(
            combined(&[artifact(&[1.0]), artifact(&[1.0])]).expect("equal lengths"),
            vec![1.0]
        );
        // And the worst filter decides, not the last.
        assert_eq!(
            by_type(
                &[artifact(&[0.3]), artifact(&[0.7]), artifact(&[0.1])],
                ErrorType::Artifact
            )
            .expect("equal lengths"),
            vec![0.7]
        );
    }

    #[test]
    fn the_types_are_independent() {
        // 1 - (1 - 0.5)(1 - 0.5).
        assert_eq!(
            combined(&[artifact(&[0.5]), non_somatic(&[0.5])]).expect("equal lengths"),
            vec![0.75]
        );
        // A type with no filter at all contributes nothing rather than a zero-length list.
        assert_eq!(combined(&[artifact(&[0.5])]).expect("one type"), vec![0.5]);
    }

    #[test]
    fn ragged_lists_are_a_refusal() {
        let error = combined(&[artifact(&[0.0, 0.0]), artifact(&[0.0])]).expect_err("two and one");
        assert_eq!(error.class(), "java.lang.IllegalArgumentException");
        assert_eq!(error.message(), "lists are not the same size");
        // The same across types.
        assert!(combined(&[artifact(&[0.0, 0.0]), non_somatic(&[0.0])]).is_err());
    }

    #[test]
    fn a_site_filter_copies_its_answer_to_every_allele() {
        // Which is what keeps the lists equal in the ordinary case.
        let site = artifact(&[1.0, 1.0]);
        let per_allele = artifact(&[1.0, 0.0]);
        assert_eq!(
            combined(&[site, per_allele]).expect("equal lengths"),
            vec![1.0, 1.0]
        );
    }
}
