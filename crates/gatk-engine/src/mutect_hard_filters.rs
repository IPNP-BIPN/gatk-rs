//! Mutect's engine-free hard filters, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering` (GATK 4.6.2.0).
//!
//! Seven filters that read one INFO annotation and answer without consulting the filtering engine:
//! `BaseQualityFilter`, `MappingQualityFilter`, `ReadPositionFilter`, `StrictStrandBiasFilter`,
//! `FragmentLengthFilter`, `ClusteredEventsFilter` and `MultiallelicFilter`.
//!
//! # Only a long insertion is judged by the reference's mapping quality
//!
//! ```java
//! final int refQual = mappingQualityByAllele.remove(0);
//! new IndexRange(0, mappingQualityByAllele.size()).forEach(i -> {
//!     if (indelLengths != null && indelLengths.get(i) >= longIndelSize) {
//!         mappingQualityByAllele.set(i, refQual);
//!     }
//! });
//! ```
//!
//! The comment above it says an indel that maps uniquely still maps badly, so the region's
//! mappability is the better proxy. But `getIndelLengths()` is the **alt** length minus the **ref**
//! length, which is negative for a deletion, and the test is `>=`. A seven-base deletion is
//! therefore judged on its own poor mapping quality while a seven-base insertion with the same
//! annotation is rescued. The golden runs exactly that pair.
//!
//! # A negative median read position is never an artifact
//!
//! ```java
//! // a negative value is possible due to a bug: https://github.com/broadinstitute/gatk/issues/5492
//! .map(readPos -> readPos > -1 && readPos < minMedianReadPosition)
//! ```
//!
//! The guard is the workaround, and it is why `MPOS = -1` passes a filter that a position of 2
//! fails. `MPOS` also has no reference entry, unlike `MBQ` and `MMQ`, so this filter reads its whole
//! list where the other two skip the first element.
//!
//! # The multiallelic threshold is in other units than the annotation
//!
//! `getTumorLogOdds` reads `TLOD`, which is written as **log10** odds, and converts it with
//! `MathUtils.log10ToLog` before the filter compares it to its hard-coded `5.0`. A `TLOD` of 4.0 is
//! therefore over the threshold rather than under it, and the golden's triallelic record, whose two
//! alternates score 6.0 and 4.0, is multiallelic on both of them.
//!
//! # Four filters break four ways on a record with no annotations
//!
//! Nothing here guards the annotations: the guard is `requiredInfoAnnotations` one level up, in
//! `Mutect2VariantFilter.errorProbabilities`. Called directly, the mapping-quality filter and the
//! fragment-length filter each run off the end of an empty list, the clustered-events filter takes
//! the maximum of nothing, and the multiallelic filter dereferences a null array. [`HardFilterError`]
//! is those four, with the exceptions the reference throws.

/// What a hard filter does when the annotation it needs is not there.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HardFilterError {
    /// `List.remove(0)` or `List.get(1)` past the end.
    IndexOutOfBounds { index: usize, length: usize },
    /// `OptionalInt.getAsInt()` on nothing.
    NoSuchElement,
    /// The tumour log odds array, which was never there.
    NullArray,
}

impl HardFilterError {
    pub fn class(&self) -> &'static str {
        match self {
            HardFilterError::IndexOutOfBounds { .. } => "java.lang.IndexOutOfBoundsException",
            HardFilterError::NoSuchElement => "java.util.NoSuchElementException",
            HardFilterError::NullArray => "java.lang.NullPointerException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            HardFilterError::IndexOutOfBounds { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
            HardFilterError::NoSuchElement => "No value present".to_string(),
            HardFilterError::NullArray => {
                "Cannot read the array length because \"array\" is null".to_string()
            }
        }
    }
}

/// `BaseQualityFilter`: `MBQ` carries the reference first, which is skipped.
pub fn base_quality_artifacts(median_base_qualities: &[i32], minimum: f64) -> Vec<bool> {
    median_base_qualities
        .iter()
        .skip(1)
        .map(|quality| (*quality as f64) < minimum)
        .collect()
}

/// `MappingQualityFilter`: the reference's own quality stands in for a long **insertion**.
///
/// `indel_lengths` is `getIndelLengths()`, which is `None` for a record that is not an indel and
/// negative for a deletion.
pub fn mapping_quality_artifacts(
    median_mapping_qualities: &[i32],
    indel_lengths: Option<&[i32]>,
    minimum: f64,
    long_indel_size: i32,
) -> Result<Vec<bool>, HardFilterError> {
    let Some((reference_quality, alternates)) = median_mapping_qualities.split_first() else {
        // `remove(0)` on an empty list.
        return Err(HardFilterError::IndexOutOfBounds {
            index: 0,
            length: 0,
        });
    };
    let mut qualities: Vec<i32> = alternates.to_vec();
    for (index, quality) in qualities.iter_mut().enumerate() {
        if let Some(lengths) = indel_lengths {
            // `>=`, against a length that is negative for every deletion.
            if lengths
                .get(index)
                .is_some_and(|length| *length >= long_indel_size)
            {
                *quality = *reference_quality;
            }
        }
    }
    Ok(qualities
        .into_iter()
        .map(|quality| (quality as f64) < minimum)
        .collect())
}

/// `ReadPositionFilter`: `MPOS` has no reference entry, and `-1` is never an artifact.
pub fn read_position_artifacts(median_read_positions: &[i32], minimum: f64) -> Vec<bool> {
    median_read_positions
        .iter()
        .map(|position| *position > -1 && (*position as f64) < minimum)
        .collect()
}

/// `StrictStrandBiasFilter`: an allele with no read on one strand at all.
///
/// The answer is **empty**, not a list of falses, when the filter is switched off or the strand
/// table is missing or holds one entry, so a caller cannot assume one boolean per allele.
pub fn strict_strand_artifacts(
    strand_counts_by_allele: &[Vec<i32>],
    minimum_reads_on_each_strand: i32,
) -> Vec<bool> {
    if minimum_reads_on_each_strand == 0 || strand_counts_by_allele.len() <= 1 {
        return Vec::new();
    }
    strand_counts_by_allele[1..]
        .iter()
        .map(|counts| counts.contains(&0))
        .collect()
}

/// `FragmentLengthFilter`: the first alternate alone, against the reference.
pub fn fragment_length_is_artifact(
    median_fragment_lengths: &[i32],
    maximum_difference: f64,
) -> Result<bool, HardFilterError> {
    let Some(alternate) = median_fragment_lengths.get(1) else {
        return Err(HardFilterError::IndexOutOfBounds {
            index: 1,
            length: median_fragment_lengths.len(),
        });
    };
    let reference = median_fragment_lengths[0];
    Ok((alternate - reference).abs() as f64 > maximum_difference)
}

/// `ClusteredEventsFilter`: the worst haplotype, or the whole region.
pub fn clustered_events_is_artifact(
    haplotype_event_counts: &[i32],
    region_event_count: i32,
    max_events_in_region: i32,
    max_events_in_haplotype: i32,
) -> Result<bool, HardFilterError> {
    // `.max().getAsInt()` with nothing to take the maximum of.
    let Some(worst) = haplotype_event_counts.iter().max() else {
        return Err(HardFilterError::NoSuchElement);
    };
    Ok(*worst > max_events_in_haplotype || region_event_count > max_events_in_region)
}

/// `MultiallelicFilter`, whose LOD threshold is hard-coded and not the one it was given.
pub const MULTIALLELIC_LOD_THRESHOLD: f64 = 5.0;

/// `MultiallelicFilter`: how many alternates the tumour evidence really supports.
///
/// **The threshold is in natural-log units and the annotation is not.** `getTumorLogOdds` reads
/// `TLOD`, which `Mutect2` writes as log10 odds, and converts it with `MathUtils.log10ToLog` before
/// anything compares it to `5.0`. A `TLOD` of 4.0 is therefore over the threshold, not under it:
/// 4 * ln(10) is 9.2. Comparing the raw annotation would count half as many alleles.
pub fn multiallelic_is_artifact(
    tumour_log_10_odds: Option<&[f64]>,
    number_of_alt_alleles_threshold: usize,
) -> Result<bool, HardFilterError> {
    let Some(log_10_odds) = tumour_log_10_odds else {
        return Err(HardFilterError::NullArray);
    };
    let passing = log_10_odds
        .iter()
        .map(|odds| crate::allele_likelihoods::log10_to_log(*odds))
        .filter(|odds| *odds > MULTIALLELIC_LOD_THRESHOLD)
        .count();
    Ok(passing > number_of_alt_alleles_threshold)
}

/// `M2FiltersArgumentCollection.DEFAULT_MIN_UNIQUE_ALT_READS`, against `count <= threshold`.
///
/// Zero, and a unique-alt-read count is at least one, so the filter cannot fire as configured.
pub const DEFAULT_MIN_UNIQUE_ALT_READS: i32 = 0;

/// `M2FiltersArgumentCollection.DEFAULT_MAX_N_RATIO`, against `ratio >= maximum`.
///
/// Positive infinity, and no finite ratio reaches it.
pub const DEFAULT_MAX_N_RATIO: f64 = f64::INFINITY;

/// `M2FiltersArgumentCollection.DEFAULT_MIN_AF`, against `max < minimum`.
///
/// Zero, and no allele fraction is below it.
pub const DEFAULT_MIN_AF: f64 = 0.0;

/// `DuplicatedAltReadFilter`: PCR duplicates with unique UMIs amplifying one allele.
///
/// The list is as long as `AS_UNIQ_ALT_READ_COUNT`, **not** as long as the record's alternate
/// alleles: the filter maps over the annotation and nothing compares the two. A four-element
/// annotation on a two-alternate record answers four probabilities, and a one-element one answers
/// one.
pub fn duplicated_alt_read_artifacts(unique_alt_read_counts: &[i32], threshold: i32) -> Vec<bool> {
    unique_alt_read_counts
        .iter()
        .map(|count| *count <= threshold)
        .collect()
}

/// `NRatioFilter`: too many Ns beside the alternate reads.
///
/// `allele_depths` is `sumADsOverSamples(vc, true, true)`, tumour **and** normal, and the alternate
/// count is the total minus the reference's entry.
///
/// The comment above the guard says "if there is no NCount annotation or the altCount is 0, don't
/// apply the filter", and only the second half is written. A missing `NCount` is
/// `getAttributeAsInt(key, 0)`, a zero rather than a skip; the skip lives one level up, in
/// `requiredInfoAnnotations`. The comparison is `>=`, so a ratio exactly at the maximum is an
/// artifact.
pub fn n_ratio_is_artifact(allele_depths: &[i32], n_count: i32, max_n_ratio: f64) -> bool {
    let total: i32 = allele_depths.iter().map(|d| i64::from(*d)).sum::<i64>() as i32;
    let alt_count = total - allele_depths[0];
    if alt_count == 0 {
        return false;
    }
    f64::from(n_count) / f64::from(alt_count) >= max_n_ratio
}

/// `MinAlleleFractionFilter`: the tumour's allele fraction is too low to believe.
///
/// `fractions_by_alt_allele` is what `getAltDataByAllele` gathered over the tumour genotypes that
/// carry `AF`, one list per alternate allele.
///
/// Two things this keeps. **An allele with no data is `orElse(1.0)`**, which is below no threshold,
/// so the filter answers "not an artifact" from an absence; `requiredInfoAnnotations` is empty, so
/// a record with no `AF` anywhere still reaches this and still passes. And the reference's
/// `.filter(entry -> !vc.getReference().equals(entry.getKey()))` is **dead code**, since the map it
/// filters is keyed on the alternate alleles alone; there is nothing here to model.
pub fn min_allele_fraction_artifacts(
    fractions_by_alt_allele: &[Vec<f64>],
    minimum: f64,
) -> Vec<bool> {
    fractions_by_alt_allele
        .iter()
        .map(|fractions| {
            // `Stream.max(Double::compare)`, which is a total order: `-0.0 < 0.0`, and NaN is the
            // largest of all. `total_cmp` is the same order; `f64::max` is not.
            fractions
                .iter()
                .copied()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap_or(1.0)
                < minimum
        })
        .collect()
}

/// `PanelOfNormalsFilter`: `vc.hasAttribute(PON)`.
///
/// Presence, not value. A record annotated `PON=false` or `PON=""` is filtered exactly as one
/// annotated `PON=true` is.
pub fn panel_of_normals_is_artifact(has_pon_attribute: bool) -> bool {
    has_pon_attribute
}

/// `HardFilter.calculateErrorProbability`: a hard filter answers one or zero.
pub fn error_probability(is_artifact: bool) -> f64 {
    if is_artifact {
        1.0
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_a_long_insertion_is_judged_by_the_reference() {
        let qualities = [60, 20];
        // A seven-base deletion: the indel length is negative, so no substitution happens.
        assert_eq!(
            mapping_quality_artifacts(&qualities, Some(&[-7]), 30.0, 5).expect("a ref entry"),
            vec![true]
        );
        // A seven-base insertion with the same annotation: the reference's 60 stands in.
        assert_eq!(
            mapping_quality_artifacts(&qualities, Some(&[7]), 30.0, 5).expect("a ref entry"),
            vec![false]
        );
        // A SNP has no indel lengths at all.
        assert_eq!(
            mapping_quality_artifacts(&qualities, None, 30.0, 5).expect("a ref entry"),
            vec![true]
        );
    }

    #[test]
    fn a_negative_median_read_position_is_never_an_artifact() {
        // MPOS has no reference entry, so every value is an alternate's.
        assert_eq!(read_position_artifacts(&[2], 5.0), vec![true]);
        assert_eq!(read_position_artifacts(&[-1], 5.0), vec![false]);
        assert_eq!(read_position_artifacts(&[2, 20], 5.0), vec![true, false]);
    }

    #[test]
    fn the_base_quality_filter_skips_the_reference() {
        // MBQ carries the reference first: two entries, one answer.
        assert_eq!(base_quality_artifacts(&[30, 10], 20.0), vec![true]);
        assert_eq!(
            base_quality_artifacts(&[30, 10, 35], 20.0),
            vec![true, false]
        );
        // And a missing annotation is an empty answer rather than a refusal.
        assert!(base_quality_artifacts(&[], 20.0).is_empty());
    }

    #[test]
    fn strict_strand_bias_answers_an_empty_list_when_it_is_switched_off() {
        let counts = vec![vec![5, 5], vec![0, 7]];
        assert_eq!(strict_strand_artifacts(&counts, 1), vec![true]);
        // Switched off: empty, not `[false]`.
        assert!(strict_strand_artifacts(&counts, 0).is_empty());
        // And a table with nothing but the reference in it, or none at all.
        assert!(strict_strand_artifacts(&[vec![5, 5]], 1).is_empty());
        assert!(strict_strand_artifacts(&[], 1).is_empty());
    }

    #[test]
    fn the_fragment_length_filter_looks_at_one_allele_only() {
        // 380 - 300 is over the difference; the third allele is never consulted.
        assert!(fragment_length_is_artifact(&[300, 380, 302], 50.0).expect("two entries"));
        assert!(!fragment_length_is_artifact(&[300, 305, 999], 50.0).expect("two entries"));
        // The difference is absolute.
        assert!(fragment_length_is_artifact(&[380, 300], 50.0).expect("two entries"));
    }

    #[test]
    fn the_multiallelic_filter_counts_in_natural_log_units() {
        // A TLOD of 4.0 is 9.2 once converted, so both of these count and the site is multiallelic.
        assert!(multiallelic_is_artifact(Some(&[6.0, 4.0]), 1).expect("an array"));
        // The threshold in the annotation's own units is 5 / ln(10), about 2.17.
        assert!(multiallelic_is_artifact(Some(&[2.18, 2.18]), 1).expect("an array"));
        assert!(!multiallelic_is_artifact(Some(&[2.17, 2.17]), 1).expect("an array"));
        // One allele over it is not more than one.
        assert!(!multiallelic_is_artifact(Some(&[6.0, 1.0]), 1).expect("an array"));
    }

    #[test]
    fn four_filters_break_four_ways_on_a_record_with_no_annotations() {
        assert_eq!(
            mapping_quality_artifacts(&[], None, 30.0, 5)
                .expect_err("nothing to remove")
                .message(),
            "Index 0 out of bounds for length 0"
        );
        assert_eq!(
            fragment_length_is_artifact(&[], 50.0)
                .expect_err("nothing to get")
                .message(),
            "Index 1 out of bounds for length 0"
        );
        assert_eq!(
            clustered_events_is_artifact(&[], 0, 2, 2)
                .expect_err("no maximum")
                .class(),
            "java.util.NoSuchElementException"
        );
        assert_eq!(
            multiallelic_is_artifact(None, 1)
                .expect_err("no array")
                .message(),
            "Cannot read the array length because \"array\" is null"
        );
    }

    #[test]
    fn a_hard_filter_answers_one_or_zero() {
        assert_eq!(error_probability(true), 1.0);
        assert_eq!(error_probability(false), 0.0);
    }

    /// The three defaults are each on the wrong side of their own comparison.
    #[test]
    fn none_of_the_three_thresholds_can_be_met() {
        // `count <= 0`, and a unique-alt-read count is at least one.
        assert_eq!(
            duplicated_alt_read_artifacts(&[1, 5, 100], DEFAULT_MIN_UNIQUE_ALT_READS),
            vec![false, false, false]
        );
        // `ratio >= Infinity`, whatever the counts.
        assert!(!n_ratio_is_artifact(
            &[10, 10, 10],
            1_000_000,
            DEFAULT_MAX_N_RATIO
        ));
        // `max < 0`, and an allele fraction is not negative.
        assert_eq!(
            min_allele_fraction_artifacts(&[vec![0.0], vec![1.0]], DEFAULT_MIN_AF),
            vec![false, false]
        );
    }

    /// An allele with no fraction at all is `orElse(1.0)`: an absence answers "not an artifact".
    #[test]
    fn an_allele_with_no_data_passes_every_threshold() {
        assert_eq!(
            min_allele_fraction_artifacts(&[Vec::new(), vec![0.05]], 0.1),
            vec![false, true]
        );
    }

    /// The one guard the N-ratio filter has, and the `>=` beside it.
    #[test]
    fn the_n_ratio_guard_is_the_alternate_count_and_the_comparison_is_inclusive() {
        // Every read is reference: the alternate count is zero and nothing is computed.
        assert!(!n_ratio_is_artifact(&[100, 0, 0], 50, 0.5));
        // Exactly at the maximum, which `>=` calls an artifact.
        assert!(n_ratio_is_artifact(&[100, 20, 20], 20, 0.5));
        assert!(!n_ratio_is_artifact(&[100, 20, 20], 19, 0.5));
    }
}
