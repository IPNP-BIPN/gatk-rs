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

use crate::allele_fraction_cluster::{
    double_stream_sum, learn_beta_binomial, learn_binomial, AlleleFractionCluster,
    BetaDistributionShape, Datum, ShapeError,
};
use crate::java_format::format_decimals;
use crate::java_hash::hash_map_order;
use crate::math_utils::normalize_sum_to_one;
use crate::mutect_engine::{log_one_third, posterior_probability_of_error};
use crate::natural_log_utils::normalize_from_log_to_linear_space;
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
    /// `clustersHaveBeenInitialized`: the quantile initialisation runs once and never again.
    clusters_have_been_initialized: bool,
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
            clusters_have_been_initialized: false,
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

/// `NUM_ITERATIONS`, the rounds of EM per call to learn.
pub const NUM_ITERATIONS: usize = 5;

/// `MAX_BINOMIAL_CLUSTERS`.
pub const MAX_BINOMIAL_CLUSTERS: usize = 5;

/// `NUM_INITIALIZATION_QUANTILES`.
pub const NUM_INITIALIZATION_QUANTILES: usize = 50;

/// `MIN_QUANTILE_INDEX_FOR_MAKING_CLUSTER`, `(int) (0.1 * 50)`.
pub const MIN_QUANTILE_INDEX_FOR_MAKING_CLUSTER: usize = 5;

/// `MAX_FRACTION_OF_BACKGROUND_TO_SPLIT_OFF`.
pub const MAX_FRACTION_OF_BACKGROUND_TO_SPLIT_OFF: f64 = 0.9;

/// `REGULARIZING_PSEUDOCOUNT`.
pub const REGULARIZING_PSEUDOCOUNT: f64 = 1.0;

/// `MathUtils.binomialProbability(n, k, p)`, which is commons-math's
/// `BinomialDistribution.probability` over the saddle-point expansion.
///
/// The exponential is `FastMath.exp`, which is commons-math's own pure-Java one and therefore
/// portable: this path is **not** the one decision 0014 bounds. Rust's `f64::exp` is the system
/// libm and disagrees with it in the last bit often enough to show, which is what gatk-rs's
/// `contamination-filter` golden caught.
pub fn binomial_probability(n: i32, k: i32, p: f64) -> f64 {
    if n == 0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    if k < 0 || k > n {
        return 0.0;
    }
    let log_probability = jmath::saddle_point::log_binomial_probability(k, n, p, 1.0 - p);
    if log_probability == f64::NEG_INFINITY {
        0.0
    } else {
        jmath::fast_math::exp(log_probability)
    }
}

impl SomaticClusteringModel {
    /// `probabilityOfSomaticVariant(datum)`.
    fn probability_of_somatic_variant(&mut self, datum: &Datum) -> f64 {
        let artifact_prob = datum.artifact_prob();
        let non_somatic_prob = datum.non_sequencing_error_prob();
        let sequencing_error_prob = self
            .probability_of_sequencing_error(datum)
            .unwrap_or(f64::NAN);
        (1.0 - artifact_prob) * (1.0 - non_somatic_prob) * (1.0 - sequencing_error_prob)
    }

    /// `backgroundProbGivenSomatic(totalCount, altCount)`, the first cluster's normalised share.
    fn background_prob_given_somatic(&self, total_count: i32, alt_count: i32) -> f64 {
        normalize_from_log_to_linear_space(&self.cluster_log_likelihoods(total_count, alt_count))
            .map(|probabilities| probabilities[0])
            .unwrap_or(f64::NAN)
    }

    /// `learnAndClearAccumulatedData()`.
    pub fn learn_and_clear_accumulated_data(&mut self) -> Result<(), ShapeError> {
        if !self.clusters_have_been_initialized {
            self.initialize_clusters()?;
        }
        for _ in 0..NUM_ITERATIONS {
            self.perform_em_iteration(true)?;
        }
        self.data.clear();
        self.obvious_artifact_count = 0;
        Ok(())
    }

    /// `calculateAlleleFractionQuantiles()`.
    ///
    /// The sort is `Comparator.comparingDouble` on the allele fraction, which is stable, so equal
    /// fractions keep the order the data arrived in. The final `distinct()` keeps the first of each.
    fn allele_fraction_quantiles(&mut self) -> Vec<f64> {
        let mut fractions_and_probs: Vec<(f64, f64)> = Vec::with_capacity(self.data.len());
        for index in 0..self.data.len() {
            let datum = self.data[index];
            let fraction = datum.alt_count() as f64 / datum.total_count() as f64;
            fractions_and_probs.push((fraction, self.probability_of_somatic_variant(&datum)));
        }
        fractions_and_probs.sort_by(|a, b| a.0.total_cmp(&b.0));
        let total_somatic_prob: f64 = double_stream_sum(
            &fractions_and_probs
                .iter()
                .map(|p| p.1)
                .collect::<Vec<f64>>(),
        );
        let mut cumulative_prob = 0.0;
        let quantile_step = total_somatic_prob / NUM_INITIALIZATION_QUANTILES as f64;
        let mut quantile_prob = quantile_step;
        let mut quantiles: Vec<f64> = Vec::new();
        for (fraction, prob) in &fractions_and_probs {
            cumulative_prob += prob;
            if cumulative_prob > quantile_prob {
                quantiles.push(*fraction);
                while cumulative_prob > quantile_prob {
                    quantile_prob += quantile_step;
                }
            }
        }
        // `distinct()`, which is by equality and keeps the first.
        let mut distinct: Vec<f64> = Vec::new();
        for quantile in quantiles {
            if !distinct.contains(&quantile) {
                distinct.push(quantile);
            }
        }
        distinct
    }

    /// `calculateQuantileBackgroundResponsibilities(alleleFractionQuantiles, backgroundProbs)`.
    ///
    /// The density is the binomial one times `n + 1`, which is the flat-prior posterior density at
    /// that allele fraction.
    fn quantile_background_responsibilities(
        &self,
        quantiles: &[f64],
        background_probs: &[f64],
    ) -> Vec<f64> {
        let mut totals = vec![0.0; quantiles.len()];
        for (index, datum) in self.data.iter().enumerate() {
            let background_prob = background_probs[index];
            for (q, fraction) in quantiles.iter().enumerate() {
                let density =
                    binomial_probability(datum.total_count(), datum.alt_count(), *fraction);
                totals[q] += density * background_prob * (datum.total_count() as f64 + 1.0);
            }
        }
        totals
    }

    /// `calculatePeaksAndMasses(alleleFractionQuantiles, totalQuantileResponsibilities)`.
    ///
    /// Trapezoid quadrature between local minima, where the minimum test is on
    /// `Double.compare` rather than on `<`, so a NaN responsibility sorts above everything.
    fn peaks_and_masses(quantiles: &[f64], responsibilities: &[f64]) -> Vec<(f64, f64)> {
        let mut peaks: Vec<(f64, f64)> = Vec::new();
        let mut current_peak_mass = 0.0;
        let mut current_peak = 0.0;
        let mut current_peak_responsibility = 0.0;
        for q in 0..quantiles.len() {
            let left_responsibility = if q == 0 { 0.0 } else { responsibilities[q - 1] };
            let responsibility = responsibilities[q];
            let right_responsibility = if q == quantiles.len() - 1 {
                0.0
            } else {
                responsibilities[q + 1]
            };
            let left_fraction = if q == 0 { 0.0 } else { quantiles[q - 1] };
            let fraction = quantiles[q];
            current_peak_mass +=
                (fraction - left_fraction) * (left_responsibility + responsibility) / 2.0;
            if responsibility > current_peak_responsibility {
                current_peak = fraction;
                current_peak_responsibility = responsibility;
            }
            let left_compare = responsibility.total_cmp(&left_responsibility);
            let right_compare = responsibility.total_cmp(&right_responsibility);
            let local_min = (left_compare.is_lt() && right_compare.is_le())
                || (left_compare.is_le() && right_compare.is_lt());
            if (local_min && q > 0) || q == quantiles.len() - 1 {
                peaks.push((current_peak, current_peak_mass));
                current_peak_mass = 0.0;
                current_peak = fraction;
                current_peak_responsibility = responsibility;
            }
        }
        peaks
    }

    /// `initializeClusters()`.
    ///
    /// Splits the biggest peak off the background, runs five silent EM iterations, and keeps the
    /// split only while the BIC improves -- at most five times. A peak below the 0.1 quantile stops
    /// it outright.
    fn initialize_clusters(&mut self) -> Result<(), ShapeError> {
        let data = self.data.clone();
        let somatic_probs: Vec<f64> = data
            .iter()
            .map(|datum| self.probability_of_somatic_variant(datum))
            .collect();
        let mut previous_bic = f64::NEG_INFINITY;
        for _ in 0..MAX_BINOMIAL_CLUSTERS {
            let old_log_cluster_weights = self.log_cluster_weights.clone();
            let background_probs: Vec<f64> = data
                .iter()
                .enumerate()
                .map(|(index, datum)| {
                    somatic_probs[index]
                        * self.background_prob_given_somatic(datum.total_count(), datum.alt_count())
                })
                .collect();
            let quantiles = self.allele_fraction_quantiles();
            let responsibilities =
                self.quantile_background_responsibilities(&quantiles, &background_probs);
            let peaks = Self::peaks_and_masses(&quantiles, &responsibilities);
            if peaks.is_empty() {
                break;
            }
            // `sorted(comparingDouble(getRight).reversed()).findFirst()`: a stable sort, so the
            // first of equal masses wins.
            let mut sorted = peaks.clone();
            sorted.sort_by(|a, b| b.1.total_cmp(&a.1));
            let (biggest_peak, biggest_mass) = sorted[0];
            let floor_index = MIN_QUANTILE_INDEX_FOR_MAKING_CLUSTER.min(quantiles.len() - 1);
            if biggest_peak < quantiles[floor_index] {
                break;
            }
            let total_mass = double_stream_sum(&peaks.iter().map(|p| p.1).collect::<Vec<f64>>());
            let fraction_of_background_to_split =
                MAX_FRACTION_OF_BACKGROUND_TO_SPLIT_OFF.min(biggest_mass / total_mass);
            let new_cluster_log_weight =
                fraction_of_background_to_split.ln() + self.log_cluster_weights[0];
            // `log1p`, not `log(1 - x)`: the background keeps MORE weight than it had.
            let new_background_weight =
                fraction_of_background_to_split.ln_1p() + self.log_cluster_weights[0];
            self.clusters
                .push(AlleleFractionCluster::binomial(biggest_peak)?);
            self.log_cluster_weights.push(new_cluster_log_weight);
            self.log_cluster_weights[0] = new_background_weight;

            for _ in 0..NUM_ITERATIONS {
                self.perform_em_iteration(false)?;
            }

            let log_likelihoods: Vec<f64> = data
                .iter()
                .map(|datum| {
                    self.log_likelihood_given_somatic(datum.total_count(), datum.alt_count())
                        .unwrap_or(f64::NAN)
                })
                .collect();
            let weighted: Vec<f64> = log_likelihoods
                .iter()
                .enumerate()
                .map(|(index, value)| somatic_probs[index] * value)
                .collect();
            // `MathUtils.sum`, a plain loop, over a product formed first: not the compensated sum.
            let weighted_log_likelihood = crate::somatic_likelihoods::sum(&weighted);
            let effective_somatic_count = crate::somatic_likelihoods::sum(&somatic_probs);
            let num_parameters = 2.0 * self.clusters.len() as f64;
            let current_bic =
                weighted_log_likelihood - num_parameters * effective_somatic_count.ln();
            if current_bic < previous_bic {
                self.clusters.pop();
                self.log_cluster_weights = old_log_cluster_weights;
                break;
            }
            previous_bic = current_bic;
        }
        self.clusters_have_been_initialized = true;
        Ok(())
    }

    /// `performEMIteration(updateSomaticPriors)`.
    fn perform_em_iteration(&mut self, update_somatic_priors: bool) -> Result<(), ShapeError> {
        // `Collectors.toMap` over `-10..=10`, which is a HashMap: its iteration order is the one the
        // variant count is summed in, and a length outside the window is appended by `putIfAbsent`.
        let mut counts: Vec<(i32, f64)> = (-MAX_INDEL_SIZE_IN_PRIOR_MAP
            ..=MAX_INDEL_SIZE_IN_PRIOR_MAP)
            .map(|length| (length, 0.0))
            .collect();
        let data = self.data.clone();
        let mut responsibilities: Vec<Vec<f64>> = Vec::with_capacity(data.len());
        let mut total_cluster_responsibilities = vec![0.0; self.clusters.len()];
        for datum in &data {
            let somatic_prob = self.probability_of_somatic_variant(datum);
            let indel_length = datum.indel_length();
            if !counts.iter().any(|(length, _)| *length == indel_length) {
                counts.push((indel_length, 0.0));
            }
            if let Some(entry) = counts
                .iter_mut()
                .find(|(length, _)| *length == indel_length)
            {
                entry.1 += somatic_prob;
            }
            let cluster_log_likelihoods =
                self.cluster_log_likelihoods(datum.total_count(), datum.alt_count());
            let if_somatic = normalize_from_log_to_linear_space(&cluster_log_likelihoods)
                .unwrap_or_else(|_| vec![f64::NAN; cluster_log_likelihoods.len()]);
            let scaled: Vec<f64> = if_somatic
                .iter()
                .map(|value| somatic_prob * value)
                .collect();
            for (index, value) in scaled.iter().enumerate() {
                total_cluster_responsibilities[index] += value;
            }
            responsibilities.push(scaled);
        }
        for value in total_cluster_responsibilities.iter_mut() {
            *value += REGULARIZING_PSEUDOCOUNT;
        }
        self.log_cluster_weights = normalize_sum_to_one(&total_cluster_responsibilities)
            .expect("the pseudocount keeps the sum positive")
            .iter()
            .map(|value| value.ln())
            .collect();
        let technical_artifact_count = self.obvious_artifact_count as f64
            + double_stream_sum(
                &data
                    .iter()
                    .map(|datum| datum.artifact_prob())
                    .collect::<Vec<f64>>(),
            );
        // The map's values in the HashMap's own iteration order, summed the compensated way.
        let variant_count = double_stream_sum(&Self::hash_map_values(&counts));

        if update_somatic_priors {
            self.log_variant_vs_artifact_prior = ((variant_count + REGULARIZING_PSEUDOCOUNT)
                / (variant_count + technical_artifact_count + REGULARIZING_PSEUDOCOUNT * 2.0))
                .ln();
            if let Some(callable_sites) = self.callable_sites {
                for length in -MAX_INDEL_SIZE_IN_PRIOR_MAP..=MAX_INDEL_SIZE_IN_PRIOR_MAP {
                    let empirical_ratio = counts
                        .iter()
                        .find(|(key, _)| *key == length)
                        .map(|(_, value)| *value)
                        .unwrap_or(0.0)
                        / callable_sites;
                    let floor = if length == 0 { 1.0e-8 } else { 1.0e-9 };
                    let prior = empirical_ratio.max(floor).ln();
                    if let Some(entry) = self
                        .log_variant_priors
                        .iter_mut()
                        .find(|(key, _)| *key == length)
                    {
                        entry.1 = prior;
                    } else {
                        self.log_variant_priors.push((length, prior));
                    }
                }
            }
        }

        for index in 0..self.clusters.len() {
            let for_this_cluster: Vec<f64> = responsibilities
                .iter()
                .map(|values| values[index])
                .collect();
            let shape = self.clusters[index].shape();
            self.clusters[index] = match self.clusters[index] {
                AlleleFractionCluster::Binomial(_) => {
                    AlleleFractionCluster::Binomial(learn_binomial(&data, &for_this_cluster)?)
                }
                AlleleFractionCluster::BetaBinomial(_) => AlleleFractionCluster::BetaBinomial(
                    learn_beta_binomial(shape, &data, &for_this_cluster)?,
                ),
            };
        }
        Ok(())
    }

    /// The values of the reference's `HashMap<Integer, MutableDouble>`, in its iteration order.
    ///
    /// `Integer.hashCode` is the value itself, so the order is the table's and not the insertion
    /// order: it is what decides how the variant count's compensated sum accumulates.
    fn hash_map_values(counts: &[(i32, f64)]) -> Vec<f64> {
        // `Integer.hashCode` is the value itself.
        let entries: Vec<(i32, i32)> = counts
            .iter()
            .map(|(length, _)| (*length, *length))
            .collect();
        match hash_map_order(&entries) {
            Ok(order) => order
                .iter()
                .map(|key| {
                    counts
                        .iter()
                        .find(|(length, _)| length == key)
                        .map(|(_, value)| *value)
                        .unwrap_or(0.0)
                })
                .collect(),
            // An order the hash port has not measured: fall back to insertion order rather than
            // guessing, which a test will catch as a wrong sum rather than as a silent one.
            Err(_) => counts.iter().map(|(_, value)| *value).collect(),
        }
    }

    /// `clusteringMetadata()`, whose numbers are formatted before anyone sees them.
    ///
    /// `%.4f` on the weights and `%.2f`/`%.3f` on the shapes, all of them Java's HALF_UP. A cluster
    /// that moved in its last digits reads here as one that did not.
    pub fn clustering_metadata(&self) -> Vec<(String, String)> {
        let mut result = Vec::new();
        for length in -MAX_INDEL_SIZE_IN_PRIOR_MAP..=MAX_INDEL_SIZE_IN_PRIOR_MAP {
            let log_prior = self
                .log_variant_priors
                .iter()
                .find(|(key, _)| *key == length)
                .map(|(_, value)| *value)
                .unwrap_or(f64::NAN);
            let kind = if length == 0 {
                "SNV".to_string()
            } else if length < 0 {
                format!("deletion of length {}", length.abs())
            } else {
                format!("insertion of length {length}")
            };
            result.push((
                format!("Ln prior of {kind}"),
                crate::tsv_table::java_double_to_string(log_prior),
            ));
        }
        result.push((
            "Background beta-binomial cluster".to_string(),
            format!(
                "weight = {}, {}",
                format_decimals(self.log_cluster_weights[0].exp(), 4),
                Self::describe(&self.clusters[0])
            ),
        ));
        result.push((
            "High-AF beta-binomial cluster".to_string(),
            format!(
                "weight = {}, {}",
                format_decimals(self.log_cluster_weights[1].exp(), 4),
                Self::describe(&self.clusters[1])
            ),
        ));
        let mut rest: Vec<usize> = (2..self.clusters.len()).collect();
        // `sorted(comparingDouble(c -> -logClusterWeights[c]))`, a stable sort on the negated weight.
        rest.sort_by(|a, b| {
            (-self.log_cluster_weights[*a]).total_cmp(&-self.log_cluster_weights[*b])
        });
        for index in rest {
            result.push((
                "Binomial cluster".to_string(),
                format!(
                    "weight = {}, {}",
                    format_decimals(self.log_cluster_weights[index].exp(), 4),
                    Self::describe(&self.clusters[index])
                ),
            ));
        }
        result
    }

    /// The two `toString`s: `alpha = %.2f, beta = %.2f` and `mean = %.3f`.
    fn describe(cluster: &AlleleFractionCluster) -> String {
        match cluster {
            AlleleFractionCluster::BetaBinomial(shape) => format!(
                "alpha = {}, beta = {}",
                format_decimals(shape.alpha(), 2),
                format_decimals(shape.beta(), 2)
            ),
            AlleleFractionCluster::Binomial(shape) => format!(
                "mean = {}",
                format_decimals(shape.alpha() / (shape.alpha() + shape.beta()), 3)
            ),
        }
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
