//! `Mutect2FilteringEngine.applyFiltersAndAccumulateOutputStats`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering` (GATK 4.6.2.0).
//!
//! The step above the ten filters: probabilities in, the FILTER column and the `AS_FilterStatus`
//! annotation out. No arithmetic beyond one comparison and one rounding.
//!
//! # The symbolic-allele removal is applied twice, and the second one is a defect
//!
//! `ErrorProbabilities`' constructor already calls `removeDataForSymbolicAltAlleles` on every
//! filter's probability list. `applyFilters` then calls it **again**, on the per-allele filter
//! strings derived from those already-shortened lists:
//!
//! ```java
//! List<List<String>> filtersNonSymbolicAlleles = GATKVariantContextUtils.removeDataForSymbolicAltAlleles(vc, distinctFiltersByAllele);
//! if (!filtersNonSymbolicAlleles.stream().anyMatch(filterList -> filterList.contains(SITE_LEVEL_FILTERS))) {
//!     siteFiltersWithErrorProb.put(GATKVCFConstants.FAIL, 1.0);
//! }
//! ```
//!
//! On a record whose symbolic allele comes **first**, the second removal deletes the only surviving
//! entry, `anyMatch` over an empty list is false, and the site is `FAIL`ed although its one real
//! allele passed every filter. The golden runs the same record twice with the symbolic allele last
//! and first, and only the second is failed. This port reproduces it.
//!
//! # A record with no per-allele filter at all is a refusal
//!
//! `orderedASFilterStrings` walks a `ListIterator` over the transposed per-allele strings, and a
//! symbolic allele takes `SITE` **without advancing it**. If no filter answered per allele, that
//! iterator is empty and the first real alternate is a `NoSuchElementException` -- which is one
//! record's worth of missing annotations away, since every per-allele filter answers an empty list
//! when its annotations are absent.
//!
//! # A filter can fire without being named
//!
//! ```java
//! final double maxErrorProb = siteFiltersWithErrorProb.values().stream().mapToDouble(p->p).max().orElse(1);
//! if (entry.getValue() >= Math.min(maxErrorProb, MIN_REPORTABLE_ERROR_PROBABILITY)) { vcb.filter(entry.getKey()); }
//! ```
//!
//! The comment beside it says this "will not change the status of whether a variant is actually
//! filtered or not". It changes which names appear: `contamination` at `0.05` beside `germline` at
//! `0.99` filters the record and is absent from the column.
//!
//! # `SITE` is a placeholder, not a filter
//!
//! Every allele that passes a filter is recorded as `SITE`; `getDistinctFiltersForAllele` removes it
//! when the allele has any real filter and adds it back when it has none. So a triallelic record can
//! be `PASS` with `base_qual|SITE` beside it: one allele filtered, the site not.

use crate::qual_quantizer::error_prob_to_qual;
use crate::somatic_clustering_model::AlternateAllele;

/// `Mutect2FilteringEngine.EPSILON`.
pub const EPSILON: f64 = 1.0e-10;

/// `Mutect2FilteringEngine.MIN_REPORTABLE_ERROR_PROBABILITY`.
pub const MIN_REPORTABLE_ERROR_PROBABILITY: f64 = 0.1;

/// `GATKVCFConstants.SITE_LEVEL_FILTERS`, which is the string `SITE` and not a filter.
pub const SITE_LEVEL_FILTERS: &str = "SITE";

/// `GATKVCFConstants.FAIL`.
pub const FAIL: &str = "FAIL";

/// `GATKVCFConstants.AS_FILTER_STATUS_KEY`.
pub const AS_FILTER_STATUS_KEY: &str = "AS_FilterStatus";

/// Which base class a filter came from, which is what `ErrorProbabilities` partitions on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterKind {
    /// `Mutect2AlleleFilter`: one probability per alternate allele.
    PerAllele,
    /// `Mutect2VariantFilter`: one probability, copied to every alternate.
    PerSite,
}

/// One filter's answer, as `ErrorProbabilities` holds it: **after** its own symbolic removal.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterAnswer {
    pub name: String,
    pub kind: FilterKind,
    /// Empty when the filter was not evaluated, which `ErrorProbabilities` drops entirely.
    pub probabilities: Vec<f64>,
    /// `phredScaledPosteriorAnnotationName`.
    pub annotation: Option<String>,
    /// Whether `requiredInfoAnnotations().stream().allMatch(vc::hasAttribute)`, which is checked a
    /// second time here for the annotation alone.
    pub required_annotations_present: bool,
}

/// What one record's application produced.
#[derive(Debug, Clone, PartialEq)]
pub struct AppliedRecord {
    /// The FILTER column, in the order the names were recorded. Empty is `PASS`.
    pub filters: Vec<String>,
    pub as_filter_status: String,
    /// The phred-scaled posterior annotations, in the order they were written.
    pub annotations: Vec<(String, u8)>,
}

/// What this step refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum ApplyError {
    /// `Utils.validateArg` in `transpose`: the per-allele lists are not all the same length.
    RaggedLists,
    /// `.next()` on the exhausted iterator over the per-allele filter strings.
    NoFilterStringForAllele,
}

impl ApplyError {
    pub fn class(&self) -> &'static str {
        match self {
            ApplyError::RaggedLists => "java.lang.IllegalArgumentException",
            ApplyError::NoFilterStringForAllele => "java.util.NoSuchElementException",
        }
    }

    /// The message, which for the second is literally `null`: the reference throws
    /// `new NoSuchElementException()` with none.
    pub fn message(&self) -> &'static str {
        match self {
            ApplyError::RaggedLists => "lists are not the same size",
            ApplyError::NoFilterStringForAllele => "null",
        }
    }
}

/// `applyFiltersAndAccumulateOutputStats`, without the statistics.
///
/// `alternates` is the record's alternate alleles in order, symbolic ones included; `answers` are
/// the filters' probabilities as `ErrorProbabilities` holds them.
pub fn apply_filters(
    answers: &[FilterAnswer],
    alternates: &[AlternateAllele],
    threshold: f64,
) -> Result<AppliedRecord, ApplyError> {
    // `Math.min(1 - EPSILON, Math.max(EPSILON, getThreshold()))`.
    let error_threshold = (1.0 - EPSILON).min(EPSILON.max(threshold));

    // A `LinkedHashMap`: insertion order, and a repeated key keeps its first position.
    let mut site_filters: Vec<(String, f64)> = Vec::new();

    // `addFilterStrings` over the per-allele filters that answered.
    let allele_status_by_filter: Vec<Vec<String>> = answers
        .iter()
        .filter(|answer| answer.kind == FilterKind::PerAllele)
        .filter(|answer| !answer.probabilities.is_empty())
        .map(|answer| {
            answer
                .probabilities
                .iter()
                .map(|value| {
                    if *value > error_threshold {
                        answer.name.clone()
                    } else {
                        SITE_LEVEL_FILTERS.to_string()
                    }
                })
                .collect()
        })
        .collect();

    let filters_by_allele = transpose(&allele_status_by_filter)?;
    let distinct_filters_by_allele: Vec<Vec<String>> = filters_by_allele
        .iter()
        .map(|filters| distinct_filters_for_allele(filters))
        .collect();

    // `AnnotationUtils.encodeStringList`, which joins with a comma.
    let merged: Vec<String> = distinct_filters_by_allele
        .iter()
        .map(|filters| filters.join(","))
        .collect();

    // The walk: a symbolic alternate takes `SITE` and does NOT advance the iterator.
    let mut next = merged.iter();
    let mut ordered = Vec::with_capacity(alternates.len());
    for alternate in alternates {
        if alternate.symbolic {
            ordered.push(SITE_LEVEL_FILTERS.to_string());
        } else {
            ordered.push(
                next.next()
                    .ok_or(ApplyError::NoFilterStringForAllele)?
                    .clone(),
            );
        }
    }
    // `encodeAnyASListWithRawDelim`, which joins with a pipe.
    let as_filter_status = ordered.join("|");

    // A site-level filter from the allele-level ones, only when every allele agrees.
    for status in &allele_status_by_filter {
        if !status.is_empty()
            && status
                .iter()
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                == 1
            && !status.iter().any(|name| name == SITE_LEVEL_FILTERS)
        {
            put(&mut site_filters, &status[0], 1.0);
        }
    }

    // The per-site filters, whose annotation is written under its own check.
    let mut annotations = Vec::new();
    for answer in answers
        .iter()
        .filter(|answer| answer.kind == FilterKind::PerSite)
        .filter(|answer| !answer.probabilities.is_empty())
    {
        let probability = answer.probabilities[0];
        if let Some(annotation) = &answer.annotation {
            if answer.required_annotations_present {
                if let Some(qual) = error_prob_to_qual(probability) {
                    annotations.push((annotation.clone(), qual));
                }
            }
        }
        if probability > error_threshold {
            put(&mut site_filters, &answer.name, probability);
        }
    }

    // Every allele filtered, for different reasons. The removal below is the second one.
    if site_filters.is_empty() && !distinct_filters_by_allele.iter().all(Vec::is_empty) {
        let non_symbolic: Vec<&Vec<String>> = distinct_filters_by_allele
            .iter()
            .enumerate()
            .filter(|(index, _)| !alternates.get(*index).map(|a| a.symbolic).unwrap_or(false))
            .map(|(_, filters)| filters)
            .collect();
        if !non_symbolic
            .iter()
            .any(|filters| filters.iter().any(|name| name == SITE_LEVEL_FILTERS))
        {
            put(&mut site_filters, FAIL, 1.0);
        }
    }

    // Only the entries that reach the floor are named.
    let max_error_probability = site_filters
        .iter()
        .map(|(_, value)| *value)
        .fold(None::<f64>, |best, value| {
            Some(match best {
                None => value,
                Some(best) => java_max(best, value),
            })
        })
        .unwrap_or(1.0);
    let floor = max_error_probability.min(MIN_REPORTABLE_ERROR_PROBABILITY);
    let filters = site_filters
        .iter()
        .filter(|(_, value)| *value >= floor)
        .map(|(name, _)| name.clone())
        .collect();

    Ok(AppliedRecord {
        filters,
        as_filter_status,
        annotations,
    })
}

/// `Map.put`, which overwrites the value and keeps the key's first position.
fn put(map: &mut Vec<(String, f64)>, key: &str, value: f64) {
    match map.iter_mut().find(|(name, _)| name == key) {
        Some(entry) => entry.1 = value,
        None => map.push((key.to_string(), value)),
    }
}

/// `Math.max`, which propagates NaN where `f64::max` does not.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

/// `ErrorProbabilities.transpose`, which refuses ragged input and answers the input unchanged when
/// it is empty.
fn transpose(lists: &[Vec<String>]) -> Result<Vec<Vec<String>>, ApplyError> {
    if lists.is_empty() {
        return Ok(Vec::new());
    }
    let length = lists[0].len();
    if lists.iter().any(|list| list.len() != length) {
        return Err(ApplyError::RaggedLists);
    }
    Ok((0..length)
        .map(|index| lists.iter().map(|list| list[index].clone()).collect())
        .collect())
}

/// `getDistinctFiltersForAllele`.
///
/// `List.remove(Object)` removes the **first** occurrence, and the list has been deduplicated by
/// then, so there is only one to remove.
fn distinct_filters_for_allele(filters: &[String]) -> Vec<String> {
    let mut results: Vec<String> = Vec::new();
    for filter in filters {
        if !results.contains(filter) {
            results.push(filter.clone());
        }
    }
    if results.len() > 1 {
        if let Some(index) = results.iter().position(|name| name == SITE_LEVEL_FILTERS) {
            results.remove(index);
        }
    }
    if results.is_empty() {
        results.push(SITE_LEVEL_FILTERS.to_string());
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;

    fn real() -> AlternateAllele {
        AlternateAllele {
            length: 1,
            symbolic: false,
        }
    }

    fn symbolic() -> AlternateAllele {
        AlternateAllele {
            length: 0,
            symbolic: true,
        }
    }

    fn allele(name: &str, probabilities: &[f64]) -> FilterAnswer {
        FilterAnswer {
            name: name.to_string(),
            kind: FilterKind::PerAllele,
            probabilities: probabilities.to_vec(),
            annotation: None,
            required_annotations_present: true,
        }
    }

    fn site(name: &str, probability: f64) -> FilterAnswer {
        FilterAnswer {
            kind: FilterKind::PerSite,
            ..allele(name, &[probability])
        }
    }

    /// The same record twice, with the symbolic allele last and first. Only the second is failed,
    /// and its one real allele passed every filter.
    #[test]
    fn where_the_symbolic_allele_sits_decides_whether_the_site_fails() {
        // Symbolic last: the surviving entry is the real allele's, which was filtered.
        let last = apply_filters(&[allele("base_qual", &[0.9])], &[real(), symbolic()], 0.5)
            .expect("applied");
        assert_eq!(last.filters, vec!["base_qual".to_string()]);
        assert_eq!(last.as_filter_status, "base_qual|SITE");

        // Symbolic first: the surviving entry is the real allele's, which PASSED, and the second
        // removal deletes it anyway.
        let first = apply_filters(&[allele("base_qual", &[0.1])], &[symbolic(), real()], 0.5)
            .expect("applied");
        assert_eq!(first.as_filter_status, "SITE|SITE", "every allele passed");
        assert_eq!(
            first.filters,
            vec![FAIL.to_string()],
            "and the site is failed"
        );
    }

    /// A filter above the threshold and below the floor fires without being named.
    #[test]
    fn a_filter_can_fire_without_being_named() {
        let applied = apply_filters(
            &[
                allele("base_qual", &[0.0]),
                site("germline", 0.99),
                site("contamination", 0.05),
            ],
            &[real()],
            0.01,
        )
        .expect("applied");
        assert_eq!(applied.filters, vec!["germline".to_string()]);

        // Alone, the same probability is the maximum and is named.
        let alone = apply_filters(
            &[allele("base_qual", &[0.0]), site("contamination", 0.05)],
            &[real()],
            0.01,
        )
        .expect("applied");
        assert_eq!(alone.filters, vec!["contamination".to_string()]);
    }

    /// No per-allele filter answered, so the walk over the filter strings has nothing to take.
    #[test]
    fn a_record_with_no_per_allele_filter_is_a_refusal() {
        assert_eq!(
            apply_filters(&[site("germline", 0.9)], &[real()], 0.5),
            Err(ApplyError::NoFilterStringForAllele)
        );
        assert_eq!(
            apply_filters(&[], &[real()], 0.5),
            Err(ApplyError::NoFilterStringForAllele)
        );
    }

    /// `SITE` is added back to an allele that has no filter, and removed from one that has any.
    #[test]
    fn site_is_a_placeholder() {
        assert_eq!(distinct_filters_for_allele(&[]), vec!["SITE".to_string()]);
        assert_eq!(
            distinct_filters_for_allele(&["SITE".to_string(), "base_qual".to_string()]),
            vec!["base_qual".to_string()]
        );
        assert_eq!(
            distinct_filters_for_allele(&["SITE".to_string(), "SITE".to_string()]),
            vec!["SITE".to_string()],
            "deduplicated to one, and one is not more than one"
        );
    }
}
