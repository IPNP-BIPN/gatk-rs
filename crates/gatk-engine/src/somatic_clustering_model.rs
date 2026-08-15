//! `SomaticClusteringModel` as constructed, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.clustering` (GATK 4.6.2.0).
//!
//! The state the constructor puts the model in, and the answers it gives from that state.
//!
//! # The prior getter mutates the map it reads
//!
//! ```java
//! if (!logVariantPriors.containsKey(indelLength)) {
//!     logVariantPriors.put(indelLength, logVariantPriors.values().stream().mapToDouble(d -> d).min().getAsDouble());
//! }
//! return logVariantPriors.get(indelLength) + (indelLength == 0 ? MathUtils.LOG_ONE_THIRD : 0);
//! ```
//!
//! The map is built for indel lengths in `-10..=10`, so anything outside that window is inserted on
//! first ask, at the minimum of what is already there. It takes `&mut self` here for that reason.
//! Nothing about it is a cache: after learning, the inserted value is one the EM iteration can
//! rewrite, and which lengths were asked about decides which ones it rewrites.
//!
//! # A SNV's prior is not an indel's prior with a different number
//!
//! Only a SNV gets `LOG_ONE_THIRD` added, the prior being per mutation and a SNV having had three
//! bases to choose from.
//!
//! # The mitochondrial defaults are decided by `==`
//!
//! `getLogSnvPrior` returns the mitochondrial default only while the field still holds the ordinary
//! default, compared by `==`. A mitochondrial run told its SNV prior explicitly, even at the same
//! value, gets the ordinary path.
//!
//! # The initial weights do not sum to one
//!
//! They are `log1p(0.01)` and `log(0.01)`, so the weights are `1.01` and `0.01`. It shows:
//! [`SomaticClusteringModel::log_likelihood_given_somatic`] at `(0, 0)` is a **positive** log
//! likelihood.
//!
//! # What is not here
//!
//! `learnAndClearAccumulatedData`, `initializeClusters` and `performEMIteration`. The EM iteration
//! and its quantile initialisation are their own slice; [`SomaticClusteringModel::record`]
//! accumulates the data they would consume.

use crate::allele_fraction_cluster::{AlleleFractionCluster, BetaDistributionShape, Datum};
use crate::mutect_engine::{log_one_third, posterior_probability_of_error};
use crate::natural_log_utils::{log_sum_exp, NonFiniteSum};

/// `MAX_INDEL_SIZE_IN_PRIOR_MAP`.
pub const MAX_INDEL_SIZE_IN_PRIOR_MAP: i32 = 10;

/// `INITIAL_HIGH_AF_WEIGHT`.
pub const INITIAL_HIGH_AF_WEIGHT: f64 = 0.01;

/// `OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD`, compared with `>` so a probability of exactly `0.9`
/// keeps its datum.
pub const OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD: f64 = 0.9;

/// The `M2FiltersArgumentCollection` fields this model reads, with the two getters over them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PriorArguments {
    pub log_snv_prior: f64,
    pub log_indel_prior: f64,
    pub initial_log_prior_of_variant_versus_artifact: f64,
    pub mitochondria: bool,
}

impl PriorArguments {
    /// `DEFAULT_LOG_SNV_PRIOR`, `log10ToLog(-6)`.
    pub fn default_log_snv_prior() -> f64 {
        crate::mutect_engine::default_log_snv_prior()
    }

    /// `DEFAULT_LOG_INDEL_PRIOR`, `log10ToLog(-7)`.
    pub fn default_log_indel_prior() -> f64 {
        crate::mutect_engine::default_log_indel_prior()
    }

    /// `DEFAULT_LOG_SNV_PRIOR_FOR_MITO`, `log10ToLog(-2.5)`.
    pub fn default_log_snv_prior_for_mito() -> f64 {
        crate::allele_likelihoods::log10_to_log(-2.5)
    }

    /// `DEFAULT_LOG_INDEL_PRIOR_FOR_MITO`, `log10ToLog(-3.75)`.
    pub fn default_log_indel_prior_for_mito() -> f64 {
        crate::allele_likelihoods::log10_to_log(-3.75)
    }

    /// `new M2FiltersArgumentCollection()`.
    pub fn new() -> Self {
        Self {
            log_snv_prior: Self::default_log_snv_prior(),
            log_indel_prior: Self::default_log_indel_prior(),
            initial_log_prior_of_variant_versus_artifact:
                crate::mutect_engine::default_log_prior_of_variant_versus_artifact(),
            mitochondria: false,
        }
    }

    /// `getLogSnvPrior()`, whose mitochondrial arm is guarded by `==` against the default.
    pub fn log_snv_prior(&self) -> f64 {
        if self.mitochondria && self.log_snv_prior == Self::default_log_snv_prior() {
            Self::default_log_snv_prior_for_mito()
        } else {
            self.log_snv_prior
        }
    }

    /// `getLogIndelPrior()`.
    pub fn log_indel_prior(&self) -> f64 {
        if self.mitochondria && self.log_indel_prior == Self::default_log_indel_prior() {
            Self::default_log_indel_prior_for_mito()
        } else {
            self.log_indel_prior
        }
    }
}

impl Default for PriorArguments {
    fn default() -> Self {
        Self::new()
    }
}

/// What `record` refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum RecordError {
    /// `Utils.validateArg`: `tumorADs must have one entry per allele including the ref allele`.
    AlleleDepthsMismatched,
}

/// One alternate allele of a record, as much of it as this model looks at.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlternateAllele {
    /// `Allele.length()`, which is **zero** for a symbolic allele rather than its text's length.
    pub length: i32,
    pub symbolic: bool,
}

/// `SomaticClusteringModel.indelLength(vc, altIndex)`.
///
/// A symbolic alternate's length is zero, so a symbolic allele against a one-base reference is an
/// indel length of `-1` rather than zero or a refusal.
pub fn indel_length(reference_length: i32, alternate: AlternateAllele) -> i32 {
    alternate.length - reference_length
}

/// The model, before any learning.
#[derive(Debug, Clone)]
pub struct SomaticClusteringModel {
    /// `logVariantPriors`, keyed by indel length. Grows when it is asked about a length it lacks.
    log_variant_priors: Vec<(i32, f64)>,
    log_variant_vs_artifact_prior: f64,
    /// `callableSites`, empty when the stats say there were fewer than one.
    callable_sites: Option<f64>,
    clusters: Vec<AlleleFractionCluster>,
    log_cluster_weights: Vec<f64>,
    /// The data `record` accumulated, which the EM iteration would consume.
    data: Vec<Datum>,
    /// `obviousArtifactCount`, incremented only on the artifact threshold and not on the other.
    obvious_artifact_count: i32,
}

impl SomaticClusteringModel {
    /// `new SomaticClusteringModel(MTFAC, mutectStats)`.
    ///
    /// `callable_sites` is the `callable` statistic if the stats carry one; a value below one is
    /// treated as absent, which is the reference's "something is seriously wrong" warning path.
    pub fn new(arguments: PriorArguments, callable_sites_from_stats: Option<f64>) -> Self {
        let mut log_variant_priors = Vec::new();
        for length in -MAX_INDEL_SIZE_IN_PRIOR_MAP..=MAX_INDEL_SIZE_IN_PRIOR_MAP {
            log_variant_priors.push((length, arguments.log_indel_prior()));
        }
        // Written after the loop, so zero holds the SNV prior rather than the indel one.
        for entry in log_variant_priors.iter_mut() {
            if entry.0 == 0 {
                entry.1 = arguments.log_snv_prior();
            }
        }
        let no_callable_sites = callable_sites_from_stats.is_some_and(|sites| sites < 1.0);
        Self {
            log_variant_priors,
            log_variant_vs_artifact_prior: arguments.initial_log_prior_of_variant_versus_artifact,
            callable_sites: if no_callable_sites {
                None
            } else {
                callable_sites_from_stats
            },
            clusters: vec![
                AlleleFractionCluster::beta_binomial(BetaDistributionShape::FLAT),
                AlleleFractionCluster::beta_binomial(
                    BetaDistributionShape::new(10.0, 1.0).expect("a valid initial shape"),
                ),
            ],
            // `log1p(0.01)` and `log(0.01)`: weights of 1.01 and 0.01, which do not sum to one.
            log_cluster_weights: vec![INITIAL_HIGH_AF_WEIGHT.ln_1p(), INITIAL_HIGH_AF_WEIGHT.ln()],
            data: Vec::new(),
            obvious_artifact_count: 0,
        }
    }

    /// `getLogPriorOfVariantVersusArtifact()`.
    pub fn log_prior_of_variant_versus_artifact(&self) -> f64 {
        self.log_variant_vs_artifact_prior
    }

    /// `callableSites`, present only when the stats named at least one.
    pub fn callable_sites(&self) -> Option<f64> {
        self.callable_sites
    }

    /// How many data `record` has kept, which is not how many it was given.
    pub fn accumulated(&self) -> usize {
        self.data.len()
    }

    /// `obviousArtifactCount`.
    pub fn obvious_artifact_count(&self) -> i32 {
        self.obvious_artifact_count
    }

    /// `getLogPriorOfSomaticVariant(indelLength)`, which inserts before it reads.
    pub fn log_prior_of_somatic_variant(&mut self, indel_length: i32) -> f64 {
        if !self
            .log_variant_priors
            .iter()
            .any(|(length, _)| *length == indel_length)
        {
            // `Stream.min` over the values, which is a plain minimum and not a NaN-aware one.
            let minimum = self
                .log_variant_priors
                .iter()
                .map(|(_, prior)| *prior)
                .fold(f64::INFINITY, f64::min);
            self.log_variant_priors.push((indel_length, minimum));
        }
        let prior = self
            .log_variant_priors
            .iter()
            .find(|(length, _)| *length == indel_length)
            .map(|(_, prior)| *prior)
            .expect("just inserted if it was missing");
        if indel_length == 0 {
            prior + log_one_third()
        } else {
            prior
        }
    }

    /// `logLikelihoodGivenSomatic(totalCount, altCount)`, the weighted sum over the clusters.
    ///
    /// The weights are not normalised, so this can be positive: at `(0, 0)` every cluster's
    /// likelihood is zero and the answer is `log(1.01 + 0.01)`'s worth of weight.
    pub fn log_likelihood_given_somatic(
        &self,
        total_count: i32,
        alt_count: i32,
    ) -> Result<f64, NonFiniteSum> {
        log_sum_exp(&self.cluster_log_likelihoods(total_count, alt_count))
    }

    fn cluster_log_likelihoods(&self, total_count: i32, alt_count: i32) -> Vec<f64> {
        self.clusters
            .iter()
            .enumerate()
            .map(|(index, cluster)| {
                self.log_cluster_weights[index]
                    + cluster
                        .log_likelihood(total_count, alt_count)
                        .expect("the initial shapes and non-negative counts are in range")
            })
            .collect()
    }

    /// `probabilityOfSequencingError(datum)`.
    ///
    /// The weighted sum here is over `correctedLogLikelihood`, not over `logLikelihood`: the TLOD
    /// enters, and the prior that follows depends on the datum's indel length, which is why this
    /// takes `&mut self`.
    pub fn probability_of_sequencing_error(&mut self, datum: &Datum) -> Result<f64, NonFiniteSum> {
        let log_cluster_likelihoods: Vec<f64> = self
            .clusters
            .iter()
            .enumerate()
            .map(|(index, cluster)| {
                self.log_cluster_weights[index] + cluster.corrected_log_likelihood(datum)
            })
            .collect();
        let variant_log_likelihood = log_sum_exp(&log_cluster_likelihoods)?;
        let prior = self.log_prior_of_somatic_variant(datum.indel_length());
        posterior_probability_of_error(variant_log_likelihood, prior)
    }

    /// `record(tumorADs, tumorLogOdds, artifactProbabilities, nonSomaticProbabilities, vc)`.
    ///
    /// `tumor_ads` is `&mut` because the reference mutates the caller's array: a symbolic alternate's
    /// depth is zeroed **in place** before the total is summed, and the engine sees the change.
    pub fn record(
        &mut self,
        tumor_ads: &mut [i32],
        tumor_log_odds: &[f64],
        artifact_probabilities: &[f64],
        non_somatic_probabilities: &[f64],
        alternates: &[AlternateAllele],
        reference_length: i32,
    ) -> Result<(), RecordError> {
        if tumor_ads.len() != alternates.len() + 1 {
            return Err(RecordError::AlleleDepthsMismatched);
        }
        for (index, alternate) in alternates.iter().enumerate() {
            if alternate.symbolic {
                tumor_ads[index + 1] = 0;
            }
        }
        let total_ad: i32 = tumor_ads.iter().sum();
        for (index, alternate) in alternates.iter().enumerate().take(tumor_log_odds.len()) {
            if alternate.symbolic {
                continue;
            }
            // The two thresholds are the same number and not the same behaviour.
            if artifact_probabilities[index] > OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD {
                self.obvious_artifact_count += 1;
                continue;
            } else if non_somatic_probabilities[index] > OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD {
                continue;
            }
            self.data.push(Datum::new(
                tumor_log_odds[index],
                artifact_probabilities[index],
                non_somatic_probabilities[index],
                tumor_ads[index + 1],
                total_ad,
                indel_length(reference_length, *alternate),
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn snp() -> AlternateAllele {
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

    #[test]
    fn the_prior_getter_grows_the_map() {
        let mut model = SomaticClusteringModel::new(PriorArguments::new(), None);
        let inside = model.log_prior_of_somatic_variant(1);
        let outside = model.log_prior_of_somatic_variant(50);
        // The minimum of a map holding one SNV prior and twenty indel priors is the indel one.
        assert_eq!(outside, inside);
        // And the length is in the map now, so asking again is a plain read.
        assert_eq!(model.log_prior_of_somatic_variant(50), outside);
    }

    #[test]
    fn only_a_snv_gets_the_third() {
        let mut model = SomaticClusteringModel::new(PriorArguments::new(), None);
        assert_eq!(
            model.log_prior_of_somatic_variant(0) - log_one_third(),
            PriorArguments::default_log_snv_prior()
        );
        assert_eq!(
            model.log_prior_of_somatic_variant(1),
            PriorArguments::default_log_indel_prior()
        );
    }

    #[test]
    fn the_mitochondrial_default_is_decided_by_an_equality() {
        let mut arguments = PriorArguments::new();
        arguments.mitochondria = true;
        assert_eq!(
            arguments.log_snv_prior(),
            PriorArguments::default_log_snv_prior_for_mito()
        );
        // Told a different value, it takes that value.
        arguments.log_snv_prior = -5.0;
        assert_eq!(arguments.log_snv_prior(), -5.0);
        // And told exactly the default, it goes back to the mitochondrial arm.
        arguments.log_snv_prior = PriorArguments::default_log_snv_prior();
        assert_eq!(
            arguments.log_snv_prior(),
            PriorArguments::default_log_snv_prior_for_mito()
        );
        // The indel prior is untouched by any of that.
        assert_eq!(
            arguments.log_indel_prior(),
            PriorArguments::default_log_indel_prior_for_mito()
        );
    }

    #[test]
    fn the_unnormalised_weights_make_a_positive_log_likelihood() {
        let model = SomaticClusteringModel::new(PriorArguments::new(), None);
        assert!(model.log_likelihood_given_somatic(0, 0).expect("finite") > 0.0);
    }

    #[test]
    fn record_zeroes_the_callers_array_and_drops_on_two_thresholds() {
        let mut model = SomaticClusteringModel::new(PriorArguments::new(), Some(1000.0));
        let mut ads = [80, 20, 5];
        model
            .record(
                &mut ads,
                &[6.0, 6.0],
                &[0.0, 0.0],
                &[0.0, 0.0],
                &[snp(), symbolic()],
                1,
            )
            .expect("recorded");
        assert_eq!(ads, [80, 20, 0], "the caller's array came back changed");
        assert_eq!(model.accumulated(), 1, "the symbolic alternate was skipped");
        // The total the datum kept is the sum after the zeroing, not before.
        assert_eq!(model.data[0].total_count(), 100);

        let mut artifact = SomaticClusteringModel::new(PriorArguments::new(), Some(1000.0));
        artifact
            .record(&mut [80, 20], &[6.0], &[0.95], &[0.0], &[snp()], 1)
            .expect("recorded");
        assert_eq!(artifact.accumulated(), 0);
        assert_eq!(artifact.obvious_artifact_count(), 1);

        let mut non_somatic = SomaticClusteringModel::new(PriorArguments::new(), Some(1000.0));
        non_somatic
            .record(&mut [80, 20], &[6.0], &[0.0], &[0.95], &[snp()], 1)
            .expect("recorded");
        assert_eq!(non_somatic.accumulated(), 0);
        // The other threshold does not count anything: this is what tells the two apart.
        assert_eq!(non_somatic.obvious_artifact_count(), 0);

        // And the threshold itself is `>`, so 0.9 exactly keeps the datum.
        let mut at_threshold = SomaticClusteringModel::new(PriorArguments::new(), Some(1000.0));
        at_threshold
            .record(&mut [80, 20], &[6.0], &[0.9], &[0.9], &[snp()], 1)
            .expect("recorded");
        assert_eq!(at_threshold.accumulated(), 1);
    }

    #[test]
    fn a_short_array_is_a_refusal_and_a_symbolic_length_is_negative() {
        let mut model = SomaticClusteringModel::new(PriorArguments::new(), None);
        assert_eq!(
            model.record(&mut [80], &[6.0], &[0.0], &[0.0], &[snp()], 1),
            Err(RecordError::AlleleDepthsMismatched)
        );
        assert_eq!(indel_length(1, symbolic()), -1);
        assert_eq!(
            indel_length(
                1,
                AlternateAllele {
                    length: 4,
                    symbolic: false
                }
            ),
            3
        );
    }

    #[test]
    fn a_callable_site_count_below_one_is_no_count_at_all() {
        assert_eq!(
            SomaticClusteringModel::new(PriorArguments::new(), Some(0.0)).callable_sites(),
            None
        );
        assert_eq!(
            SomaticClusteringModel::new(PriorArguments::new(), Some(1.0)).callable_sites(),
            Some(1.0)
        );
    }
}
