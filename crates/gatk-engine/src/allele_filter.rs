//! `Mutect2AlleleFilter` and `Mutect2Filter`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering` (GATK 4.6.2.0), with
//! `TumorEvidenceFilter` on top of them.
//!
//! What every per-allele filter stands on: gather one value per allele across the samples, answer a
//! probability per alternate allele, and refuse to answer at all when the record lacks what the
//! filter reads.
//!
//! # The gather zips two iterators and stops at the shorter one
//!
//! ```java
//! Iterator<T> alleleDataIterator = getData.apply(g).iterator();
//! Iterator<List<T>> dataByAlleleIterator = dataByAllele.values().iterator();
//! while (alleleDataIterator.hasNext() && dataByAlleleIterator.hasNext())
//!     dataByAlleleIterator.next().add(alleleDataIterator.next());
//! ```
//!
//! No exception and no padding. Two consequences the measurement pinned:
//!
//!  * **a short list shifts the gather rather than shortening one allele's list.** A genotype giving
//!    `[1]` against three alleles contributes to the first allele alone, so a second sample's
//!    `[4, 5, 6]` lands as `A=[1, 4]`, `C=[5]`, `G=[6]`;
//!  * **`getAltDataByAllele` does not skip the reference in the data.** Its map is keyed by the
//!    alternate alleles alone, but the walk still starts at the caller's first element, so a caller
//!    passing a full-length per-allele list gives the first alternate the *reference's* value and
//!    drops the last.
//!
//! # An unanswerable filter answers nothing
//!
//! `errorProbabilities` checks that every required annotation is present and returns an **empty
//! list** if one is not. `ErrorProbabilities` drops an empty list entirely rather than counting it
//! as zero, so a filter that cannot be evaluated is not a filter that found nothing.

use crate::allele_fraction_cluster::Datum;
use crate::error_probabilities::ErrorType;
use crate::mutect_engine::round_finite_precision_errors;
use crate::somatic_clustering_model::{indel_length, AlternateAllele, SomaticClusteringModel};

/// One genotype, as much of it as the gather and the depth sum look at.
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypeData<T> {
    /// Whether the engine calls this sample a tumour, which is the usual precondition.
    pub tumor: bool,
    /// `Genotype.getAD()`.
    pub allele_depths: Vec<i32>,
    /// The per-allele values the filter asked for, in the genotype's own order.
    pub values: Vec<T>,
}

/// `getDataByAllele`: one list per allele of the record, the reference included.
///
/// The keys are the alleles in the record's order, so the result is returned as a `Vec` of
/// `(allele, values)` rather than a map: it is a `LinkedHashMap` upstream, and its order is the
/// record's.
pub fn data_by_allele<T: Clone>(
    alleles: &[String],
    genotypes: &[GenotypeData<T>],
    precondition: impl Fn(&GenotypeData<T>) -> bool,
) -> Vec<(String, Vec<T>)> {
    // `Collectors.toMap(identity, .., (a, b) -> a, LinkedHashMap::new)`: a repeated allele collapses
    // to one key, and the first occurrence wins.
    let mut gathered: Vec<(String, Vec<T>)> = Vec::new();
    for allele in alleles {
        if !gathered.iter().any(|(key, _)| key == allele) {
            gathered.push((allele.clone(), Vec::new()));
        }
    }
    combine(&mut gathered, genotypes, precondition);
    gathered
}

/// `getAltDataByAllele`: one list per **alternate** allele.
///
/// The caller is expected to supply alternate-only data. Nothing checks that: the walk starts at the
/// first element either way.
pub fn alt_data_by_allele<T: Clone>(
    alleles: &[String],
    genotypes: &[GenotypeData<T>],
    precondition: impl Fn(&GenotypeData<T>) -> bool,
) -> Vec<(String, Vec<T>)> {
    data_by_allele(&alleles[1..], genotypes, precondition)
}

/// `combineDataByAllele`, which is the zip.
fn combine<T: Clone>(
    gathered: &mut [(String, Vec<T>)],
    genotypes: &[GenotypeData<T>],
    precondition: impl Fn(&GenotypeData<T>) -> bool,
) {
    for genotype in genotypes.iter().filter(|g| precondition(g)) {
        for (slot, value) in gathered.iter_mut().zip(&genotype.values) {
            slot.1.push(value.clone());
        }
    }
}

/// What `sumADsOverSamples` refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleDepthTooShort {
    pub index: usize,
    pub length: usize,
}

impl AlleleDepthTooShort {
    /// `ArrayIndexOutOfBoundsException`'s message, which names both numbers.
    pub fn message(&self) -> String {
        format!(
            "Index {} out of bounds for length {}",
            self.index, self.length
        )
    }
}

/// `Mutect2FilteringEngine.sumADsOverSamples(vc, includeTumor, includeNormal)`.
///
/// Every selected genotype is indexed by the **record's** allele count, so a genotype whose AD array
/// is shorter is an out-of-bounds rather than a skip.
pub fn sum_ads_over_samples<T>(
    allele_count: usize,
    genotypes: &[GenotypeData<T>],
    include_tumor: bool,
    include_normal: bool,
) -> Result<Vec<i32>, AlleleDepthTooShort> {
    let mut totals = vec![0; allele_count];
    for genotype in genotypes {
        if !((include_tumor && genotype.tumor) || (include_normal && !genotype.tumor)) {
            continue;
        }
        // `new IndexRange(0, vc.getNAlleles()).forEach(n -> ADs[n] += ad[n])`: the range is the
        // record's allele count, so it is the genotype's array that runs out.
        for (index, total) in totals.iter_mut().enumerate() {
            let depth = genotype
                .allele_depths
                .get(index)
                .ok_or(AlleleDepthTooShort {
                    index,
                    length: genotype.allele_depths.len(),
                })?;
            *total += depth;
        }
    }
    Ok(totals)
}

/// `Mutect2Filter.weightedMedianPosteriorProbability(depthsAndPosteriors)`.
///
/// The lowest posterior that accounts for half the total alternate depth. Three things a smoother
/// implementation would get wrong: it **sorts the caller's list in place**, which is why this takes
/// `&mut`; an empty list is `0.0` rather than a refusal; and the test is `>=`, so an even split
/// returns the lower half's element.
pub fn weighted_median_posterior_probability(depths_and_posteriors: &mut [(i32, f64)]) -> f64 {
    let total_alt_depth: i32 = depths_and_posteriors.iter().map(|(depth, _)| depth).sum();
    // `Comparator.comparingDouble`, a stable sort on the posterior alone.
    depths_and_posteriors.sort_by(|a, b| a.1.total_cmp(&b.1));
    let mut cumulative = 0;
    for (depth, posterior) in depths_and_posteriors.iter() {
        cumulative += depth;
        if cumulative * 2 >= total_alt_depth {
            return *posterior;
        }
    }
    0.0
}

/// `TumorEvidenceFilter`'s identity.
pub const TUMOR_EVIDENCE_FILTER_NAME: &str = "weak_evidence";
pub const TUMOR_EVIDENCE_ERROR_TYPE: ErrorType = ErrorType::Sequencing;
/// `phredScaledPosteriorAnnotationName`, `SEQUENCING_QUAL_KEY`.
pub const TUMOR_EVIDENCE_ANNOTATION: &str = "SEQQ";

/// `TumorEvidenceFilter.errorProbabilities(vc, engine, referenceContext)`.
///
/// `tumor_log_odds` is `getTumorLogOdds(vc)`, already converted from log10, and `None` is the record
/// that has no `TLOD` at all: the required-annotation check answers an empty list rather than a
/// probability.
///
/// The `Datum` is built with **both** probabilities zero, so the filter reads the clustering model
/// without contributing to it; and the list is as long as the log-odds, not as long as the alternate
/// alleles, so a `TLOD` shorter than the record answers for the alleles it covers and no more.
pub fn tumor_evidence_error_probabilities<T>(
    model: &mut SomaticClusteringModel,
    tumor_log_odds: Option<&[f64]>,
    genotypes: &[GenotypeData<T>],
    alleles: &[AlternateAllele],
    reference_length: i32,
) -> Vec<f64> {
    let Some(log_odds) = tumor_log_odds else {
        return Vec::new();
    };
    let allele_depths = match sum_ads_over_samples(alleles.len() + 1, genotypes, true, false) {
        Ok(depths) => depths,
        Err(_) => return Vec::new(),
    };
    let total_count: i32 = allele_depths.iter().sum();
    log_odds
        .iter()
        .enumerate()
        .map(|(index, odds)| {
            let datum = Datum::new(
                *odds,
                0.0,
                0.0,
                allele_depths[index + 1],
                total_count,
                indel_length(reference_length, alleles[index]),
            );
            let probability = model
                .probability_of_sequencing_error(&datum)
                .unwrap_or(f64::NAN);
            round_finite_precision_errors(probability)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tumor(values: Vec<i32>) -> GenotypeData<i32> {
        GenotypeData {
            tumor: true,
            allele_depths: vec![80, 20, 5],
            values,
        }
    }

    #[test]
    fn a_short_list_shifts_the_gather() {
        let alleles = ["A".to_string(), "C".to_string(), "G".to_string()];
        let genotypes = vec![tumor(vec![1]), tumor(vec![4, 5, 6])];
        let gathered = data_by_allele(&alleles, &genotypes, |g| g.tumor);
        assert_eq!(gathered[0], ("A".to_string(), vec![1, 4]));
        assert_eq!(gathered[1], ("C".to_string(), vec![5]));
        assert_eq!(gathered[2], ("G".to_string(), vec![6]));
    }

    #[test]
    fn the_alternate_gather_still_starts_at_the_callers_first_element() {
        let alleles = ["A".to_string(), "C".to_string(), "G".to_string()];
        let genotypes = vec![tumor(vec![1, 2, 3]), tumor(vec![4, 5, 6])];
        let gathered = alt_data_by_allele(&alleles, &genotypes, |g| g.tumor);
        // The first ALT takes what the reference's slot took in the full gather.
        assert_eq!(gathered[0], ("C".to_string(), vec![1, 4]));
        assert_eq!(gathered[1], ("G".to_string(), vec![2, 5]));
    }

    #[test]
    fn a_repeated_allele_collapses_to_one_key() {
        let alleles = ["A".to_string(), "C".to_string(), "C".to_string()];
        let genotypes = vec![tumor(vec![1, 2, 3])];
        let gathered = data_by_allele(&alleles, &genotypes, |g| g.tumor);
        assert_eq!(gathered.len(), 2);
        assert_eq!(gathered[1], ("C".to_string(), vec![2]));
    }

    #[test]
    fn the_median_sorts_in_place_and_leans_low() {
        let mut even = [(10, 0.1), (10, 0.9)];
        assert_eq!(weighted_median_posterior_probability(&mut even), 0.1);
        let mut unsorted = [(5, 0.9), (5, 0.2), (5, 0.5)];
        assert_eq!(weighted_median_posterior_probability(&mut unsorted), 0.5);
        assert_eq!(unsorted, [(5, 0.2), (5, 0.5), (5, 0.9)]);
        // Every depth zero: the first element already satisfies `0 >= 0`.
        let mut zero = [(0, 0.3), (0, 0.7)];
        assert_eq!(weighted_median_posterior_probability(&mut zero), 0.3);
        assert_eq!(weighted_median_posterior_probability(&mut []), 0.0);
    }

    #[test]
    fn a_short_allele_depth_array_is_out_of_bounds() {
        let genotypes = vec![GenotypeData::<i32> {
            tumor: true,
            allele_depths: vec![80, 20],
            values: Vec::new(),
        }];
        assert_eq!(
            sum_ads_over_samples(3, &genotypes, true, false),
            Err(AlleleDepthTooShort {
                index: 2,
                length: 2
            })
        );
    }

    #[test]
    fn the_flags_select_the_samples() {
        let genotypes = vec![
            GenotypeData::<i32> {
                tumor: true,
                allele_depths: vec![80, 20, 5],
                values: Vec::new(),
            },
            GenotypeData::<i32> {
                tumor: false,
                allele_depths: vec![99, 1, 0],
                values: Vec::new(),
            },
        ];
        assert_eq!(
            sum_ads_over_samples(3, &genotypes, true, false).expect("in range"),
            vec![80, 20, 5]
        );
        assert_eq!(
            sum_ads_over_samples(3, &genotypes, false, true).expect("in range"),
            vec![99, 1, 0]
        );
        assert_eq!(
            sum_ads_over_samples(3, &genotypes, true, true).expect("in range"),
            vec![179, 21, 5]
        );
        assert_eq!(
            sum_ads_over_samples(3, &genotypes, false, false).expect("in range"),
            vec![0, 0, 0]
        );
    }
}
