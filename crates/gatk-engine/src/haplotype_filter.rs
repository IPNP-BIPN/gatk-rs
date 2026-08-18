//! `FilteredHaplotypeFilter`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilteredHaplotypeFilter`
//! (GATK 4.6.2.0).
//!
//! "Did something else on this haplotype already look like an artifact?" The first filter here with
//! state that outlives a record: it accumulates one artifact probability per phased haplotype
//! across the whole callset, and answers from what the **previous** pass learned.
//!
//! # The first pass answers zero to everything
//!
//! ```java
//! private Map<String, List<Pair<Integer, Double>>> accumulatingPhasedProbabilities = new HashMap<>();
//! private Map<String, List<Pair<Integer, Double>>> phasedProbabilities = new HashMap<>();
//! ```
//!
//! One is written and the other is read, and `learnParameters` moves the first to the second
//! between passes. Before any learning the read map is empty, so this filter's output depends on
//! how many passes have run and not only on the record in front of it.
//!
//! # Two of the three name comparisons are `==`
//!
//! ```java
//! .filter(e -> e.getKey().filterName() == GATKVCFConstants.GERMLINE_RISK_FILTER_NAME)
//! .filter(e -> !(ignoreNormalArtifact && e.getKey().filterName() == GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME))
//! .filter(e -> !e.getKey().filterName().equals(filterName()))
//! ```
//!
//! The first two are object identity, and they work only because both sides are compile-time
//! constants and therefore the same interned `String`. A filter whose name is an equal but
//! separately allocated string is silently **not** matched: the golden accumulates `0.1` when the
//! germline filter's name is the constant and `0.9` when it is a copy of it, from the same
//! probabilities. The third comparison is `.equals`, and it matches a copy.
//!
//! Rust has no interning to reproduce, so the distinction is carried explicitly by
//! [`FilterIdentity`]: it says which constant a name **is**, while [`FilterAnswer::name`] carries
//! what the name *says*. A caller that conflates the two would lose exactly the behaviour the
//! golden pins.
//!
//! # The record filters itself
//!
//! The distance test is `|accumulated locus - vc.getStart()| <= max`, which includes zero. The
//! filter excludes its own *filter* from the accumulated probability, as its comment says, and does
//! not exclude its own *locus*: a site some other filter called an artifact in the first pass
//! filters itself in the second.
//!
//! # `.get()` on an empty `Optional`
//!
//! `vc.getGenotypes().stream().filter(isTumor).max(...).get()` has no guard: a record with no
//! tumour sample is a `NoSuchElementException` out of a filter, not a skip.

use crate::error_probabilities::ErrorType;
use crate::mutect_engine::round_finite_precision_errors;

/// `FilteredHaplotypeFilter`'s identity.
pub const FILTER_NAME: &str = "haplotype";

/// `GERMLINE_PROBABILITY_TO_IGNORE_NORMAL_ARTIFACT`.
pub const GERMLINE_PROBABILITY_TO_IGNORE_NORMAL_ARTIFACT: f64 = 0.25;

/// `M2FiltersArgumentCollection.DEFAULT_MAX_INTRA_HAPLOTYPE_DISTANCE`.
pub const DEFAULT_MAX_INTRA_HAPLOTYPE_DISTANCE: i32 = 100;

/// `GATKVCFConstants.GERMLINE_RISK_FILTER_NAME`.
pub const GERMLINE_RISK_FILTER_NAME: &str = "germline";

/// `GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME`.
pub const ARTIFACT_IN_NORMAL_FILTER_NAME: &str = "normal_artifact";

/// Which constant a filter's name **is**, as against what it says.
///
/// See the module note: two of the three comparisons are `==` on `String`, so a name that merely
/// equals the constant does not match. `Other` is that case, and it is not hypothetical -- the
/// golden runs it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterIdentity {
    /// The name is `GATKVCFConstants.GERMLINE_RISK_FILTER_NAME` itself.
    Germline,
    /// The name is `GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME` itself.
    NormalArtifact,
    /// Any other name, including a string equal to one of the two above.
    Other,
}

/// One filter's answer, as `ErrorProbabilities.getProbabilitiesByFilter()` holds it.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterAnswer {
    /// What `filterName()` returns, compared with `.equals` against this filter's own name.
    pub name: String,
    /// Which constant that name is, compared with `==` against the two the reference names.
    pub identity: FilterIdentity,
    pub error_type: ErrorType,
    /// One probability per alternate allele.
    pub probabilities: Vec<f64>,
}

/// One genotype, as much of it as this filter looks at.
#[derive(Debug, Clone, PartialEq)]
pub struct PhasedGenotype {
    pub tumor: bool,
    /// `AF`, empty when the annotation is absent: the reference's default is `{0.0}`, whose maximum
    /// is zero, so the two are the same number and not the same thing.
    pub allele_fractions: Vec<f64>,
    /// `PGT`.
    pub phasing_gt: Option<String>,
    /// `PID`.
    pub phasing_id: Option<String>,
}

impl PhasedGenotype {
    /// `makePhasingString`: `PGT + PID`, and nothing at all if either is missing.
    pub fn phasing_string(&self) -> Option<String> {
        match (&self.phasing_gt, &self.phasing_id) {
            (Some(gt), Some(id)) => Some(format!("{gt}{id}")),
            _ => None,
        }
    }

    /// `MathUtils.arrayMax(getAttributeAsDoubleArray(g, AF, () -> new double[] {0.0}, 0.0))`.
    fn greatest_allele_fraction(&self) -> f64 {
        if self.allele_fractions.is_empty() {
            return 0.0;
        }
        // `arrayMax` is `Arrays.stream(...).max().getAsDouble()`, whose comparison is `Math.max`
        // and therefore propagates NaN. Nothing measured supplies one.
        self.allele_fractions
            .iter()
            .copied()
            .reduce(java_max)
            .expect("the list is not empty")
    }
}

/// What this filter refuses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HaplotypeFilterError {
    /// `.get()` on the empty `Optional` a record with no tumour sample leaves.
    NoTumourGenotype,
}

impl HaplotypeFilterError {
    pub fn class(&self) -> &'static str {
        "java.util.NoSuchElementException"
    }

    pub fn message(&self) -> &'static str {
        "No value present"
    }
}

/// The filter, and the two maps it keeps.
///
/// The reference's maps are `HashMap`s whose iteration order is never observed: every read is a
/// lookup by key. They are kept here as association lists in insertion order, which makes the
/// accumulation order visible where the reference's would not be, and changes no answer.
#[derive(Debug, Clone)]
pub struct FilteredHaplotypeFilter {
    max_intra_haplotype_distance: i32,
    accumulating: Vec<(String, Vec<(i32, f64)>)>,
    learned: Vec<(String, Vec<(i32, f64)>)>,
}

impl FilteredHaplotypeFilter {
    pub fn new(max_intra_haplotype_distance: i32) -> Self {
        Self {
            max_intra_haplotype_distance,
            accumulating: Vec::new(),
            learned: Vec::new(),
        }
    }

    /// What has been accumulated but not yet learned.
    pub fn accumulated(&self) -> &[(String, Vec<(i32, f64)>)] {
        &self.accumulating
    }

    /// What the previous pass learned, which is what answers.
    pub fn learned(&self) -> &[(String, Vec<(i32, f64)>)] {
        &self.learned
    }

    /// `errorProbabilities`: one probability per alternate allele, all of them the same.
    pub fn error_probabilities(
        &self,
        genotypes: &[PhasedGenotype],
        start: i32,
        alternate_count: usize,
    ) -> Result<Vec<f64>, HaplotypeFilterError> {
        let probability = self.calculate_error_probability(genotypes, start)?;
        Ok(vec![
            round_finite_precision_errors(probability);
            alternate_count
        ])
    }

    /// `calculateErrorProbability`.
    fn calculate_error_probability(
        &self,
        genotypes: &[PhasedGenotype],
        start: i32,
    ) -> Result<f64, HaplotypeFilterError> {
        // `Stream.max` reduces with `BinaryOperator.maxBy`, which returns the LEFT argument when the
        // comparison is not negative: a tie keeps the first tumour genotype, not the last.
        let tumour = genotypes
            .iter()
            .filter(|g| g.tumor)
            .reduce(|a, b| {
                if a.greatest_allele_fraction()
                    .total_cmp(&b.greatest_allele_fraction())
                    .is_ge()
                {
                    a
                } else {
                    b
                }
            })
            .ok_or(HaplotypeFilterError::NoTumourGenotype)?;

        let Some(phasing) = tumour.phasing_string() else {
            return Ok(0.0);
        };
        let Some(phased) = self
            .learned
            .iter()
            .find(|(key, _)| *key == phasing)
            .map(|(_, values)| values)
        else {
            return Ok(0.0);
        };

        // `<=`, and the accumulated locus may be this very record's.
        //
        // `DoubleStream.max()` is `reduce(Math::max)`, so an empty stream is the `orElse(0.0)` and
        // a NaN anywhere in it propagates; `f64::max` would swallow one.
        let mut best: Option<f64> = None;
        for (_, probability) in phased.iter().filter(|(locus, _)| {
            f64::from((locus - start).abs()) <= f64::from(self.max_intra_haplotype_distance)
        }) {
            best = Some(match best {
                None => *probability,
                Some(current) => java_max(current, *probability),
            });
        }
        Ok(best.unwrap_or(0.0))
    }

    /// `accumulateDataForLearning`.
    ///
    /// `answers` is `ErrorProbabilities.getProbabilitiesByFilter()`, whose iteration order is a
    /// hash order upstream and does not matter: every use of it is a maximum.
    pub fn accumulate_data_for_learning(
        &mut self,
        answers: &[FilterAnswer],
        genotypes: &[PhasedGenotype],
        start: i32,
    ) {
        // `==` on `String`: the identity, not the characters.
        let germline_probability = answers
            .iter()
            .filter(|answer| answer.identity == FilterIdentity::Germline)
            .flat_map(|answer| answer.probabilities.iter().copied())
            .fold(None::<f64>, |best, value| Some(max_by_compare(best, value)))
            .unwrap_or(0.0);

        let ignore_normal_artifact =
            germline_probability > GERMLINE_PROBABILITY_TO_IGNORE_NORMAL_ARTIFACT;

        let artifact_probability = answers
            .iter()
            .filter(|answer| answer.error_type != ErrorType::NonSomatic)
            // Again `==`, so a name that merely equals the constant is not dropped.
            .filter(|answer| {
                !(ignore_normal_artifact && answer.identity == FilterIdentity::NormalArtifact)
            })
            // And here `.equals`, so a copy of this filter's own name IS dropped.
            .filter(|answer| answer.name != FILTER_NAME)
            .flat_map(|answer| answer.probabilities.iter().copied())
            .fold(None::<f64>, |best, value| Some(max_by_compare(best, value)))
            .unwrap_or(0.0);

        // One entry per tumour genotype, so two samples sharing a phasing string add the same locus
        // twice.
        for genotype in genotypes.iter().filter(|g| g.tumor) {
            let Some(phasing) = genotype.phasing_string() else {
                continue;
            };
            match self
                .accumulating
                .iter_mut()
                .find(|(key, _)| *key == phasing)
            {
                Some((_, values)) => values.push((start, artifact_probability)),
                None => self
                    .accumulating
                    .push((phasing, vec![(start, artifact_probability)])),
            }
        }
    }

    /// `clearAccumulatedData`, which replaces the map rather than emptying it.
    pub fn clear_accumulated_data(&mut self) {
        self.accumulating = Vec::new();
    }

    /// `learnParameters`, which moves the accumulating map to the learned one.
    ///
    /// The reference assigns the reference, so until `clearAccumulatedData` replaces the
    /// accumulating map the two names point at one object. Nothing reads it in between.
    pub fn learn_parameters(&mut self) {
        self.learned = self.accumulating.clone();
    }

    /// `learnParametersAndClearAccumulatedData`, in that order.
    pub fn learn_parameters_and_clear_accumulated_data(&mut self) {
        self.learn_parameters();
        self.clear_accumulated_data();
    }
}

/// `Math.max`, which propagates NaN where `f64::max` returns the other operand.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

/// `Stream.max(Double::compareTo)`: a total order, and the left argument wins a tie.
fn max_by_compare(best: Option<f64>, value: f64) -> f64 {
    match best {
        None => value,
        Some(best) => {
            if best.total_cmp(&value).is_ge() {
                best
            } else {
                value
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tumour(pgt: &str, allele_fraction: f64) -> PhasedGenotype {
        PhasedGenotype {
            tumor: true,
            allele_fractions: vec![allele_fraction],
            phasing_gt: Some(pgt.to_string()),
            phasing_id: Some("100_A_C".to_string()),
        }
    }

    fn answer(name: &str, identity: FilterIdentity, error_type: ErrorType, p: f64) -> FilterAnswer {
        FilterAnswer {
            name: name.to_string(),
            identity,
            error_type,
            probabilities: vec![p],
        }
    }

    /// The one behaviour a port is most likely to get wrong: two of the three comparisons are
    /// identity, and the third is equality, so an equal name is dropped by one and kept by two.
    #[test]
    fn an_equal_name_is_not_the_same_name_for_two_of_the_three_comparisons() {
        let genotypes = vec![tumour("0|1", 0.3)];
        let interned = vec![
            answer(
                GERMLINE_RISK_FILTER_NAME,
                FilterIdentity::Germline,
                ErrorType::NonSomatic,
                0.5,
            ),
            answer(
                ARTIFACT_IN_NORMAL_FILTER_NAME,
                FilterIdentity::NormalArtifact,
                ErrorType::Artifact,
                0.9,
            ),
            answer("base_qual", FilterIdentity::Other, ErrorType::Artifact, 0.1),
        ];
        let mut filter = FilteredHaplotypeFilter::new(100);
        filter.accumulate_data_for_learning(&interned, &genotypes, 100);
        assert_eq!(
            filter.accumulated()[0].1,
            vec![(100, 0.1)],
            "the constant matches"
        );

        // The same characters, a different object: the germline maximum is zero, so the
        // normal-artifact filter is kept and its 0.9 wins.
        let copied = vec![
            answer(
                GERMLINE_RISK_FILTER_NAME,
                FilterIdentity::Other,
                ErrorType::NonSomatic,
                0.5,
            ),
            answer(
                ARTIFACT_IN_NORMAL_FILTER_NAME,
                FilterIdentity::NormalArtifact,
                ErrorType::Artifact,
                0.9,
            ),
            answer("base_qual", FilterIdentity::Other, ErrorType::Artifact, 0.1),
        ];
        let mut filter = FilteredHaplotypeFilter::new(100);
        filter.accumulate_data_for_learning(&copied, &genotypes, 100);
        assert_eq!(
            filter.accumulated()[0].1,
            vec![(100, 0.9)],
            "a copy does not"
        );

        // And the filter's own name is compared with `.equals`, so a copy IS dropped.
        let own = vec![
            answer(
                FILTER_NAME,
                FilterIdentity::Other,
                ErrorType::Artifact,
                0.99,
            ),
            answer("base_qual", FilterIdentity::Other, ErrorType::Artifact, 0.1),
        ];
        let mut filter = FilteredHaplotypeFilter::new(100);
        filter.accumulate_data_for_learning(&own, &genotypes, 100);
        assert_eq!(
            filter.accumulated()[0].1,
            vec![(100, 0.1)],
            "equality, not identity"
        );
    }

    /// Nothing is learned until a pass ends, and the distance test includes zero.
    #[test]
    fn the_first_pass_answers_zero_and_the_second_filters_the_record_itself() {
        let genotypes = vec![tumour("0|1", 0.3)];
        let mut filter = FilteredHaplotypeFilter::new(100);
        filter.accumulate_data_for_learning(
            &[answer(
                "base_qual",
                FilterIdentity::Other,
                ErrorType::Artifact,
                0.8,
            )],
            &genotypes,
            100,
        );
        assert_eq!(
            filter
                .error_probabilities(&genotypes, 100, 1)
                .expect("answered"),
            vec![0.0],
            "before learning"
        );
        filter.learn_parameters_and_clear_accumulated_data();
        assert!(
            filter.accumulated().is_empty(),
            "and the accumulating map is replaced"
        );
        assert_eq!(
            filter
                .error_probabilities(&genotypes, 100, 1)
                .expect("answered"),
            vec![0.8],
            "the record filters itself"
        );
        // The boundary is inclusive, and one past it is out of range.
        assert_eq!(
            filter
                .error_probabilities(&genotypes, 200, 1)
                .expect("answered"),
            vec![0.8]
        );
        assert_eq!(
            filter
                .error_probabilities(&genotypes, 201, 1)
                .expect("answered"),
            vec![0.0]
        );
    }

    /// A record with no tumour sample is a refusal, and `.get()` is where the reference throws.
    #[test]
    fn a_record_with_no_tumour_sample_is_a_refusal() {
        let filter = FilteredHaplotypeFilter::new(100);
        let normal = PhasedGenotype {
            tumor: false,
            allele_fractions: Vec::new(),
            phasing_gt: Some("0|1".to_string()),
            phasing_id: Some("100_A_C".to_string()),
        };
        assert_eq!(
            filter.error_probabilities(&[normal], 100, 1),
            Err(HaplotypeFilterError::NoTumourGenotype)
        );
    }
}
