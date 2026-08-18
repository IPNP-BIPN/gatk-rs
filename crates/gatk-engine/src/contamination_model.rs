//! Ported from `org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationModel`.
//!
//! The probabilistic model `CalculateContamination` learns and then uses as a genotyper. The
//! comment at the top of the reference is worth keeping: the model is not what produces the
//! contamination estimate. It is learned, used to decide which sites are hom alt or hom ref, and
//! the estimate then comes from the reads at those sites.
//!
//! # Three loops decide the answer, and all three are ordered
//!
//!  * **the model iterates three times** between minor allele fractions and contamination, each
//!    round feeding the previous contamination back in, so a port that solved the joint problem
//!    once lands somewhere else;
//!  * **the loss-of-heterozygosity search walks a threshold down** from 0.40 in steps of 0.04 and
//!    stops at the FIRST threshold that keeps more than a quarter of the sites. The subtraction is
//!    repeated rather than computed, so the thresholds are 0.4, 0.36, 0.32000000000000006 and so
//!    on, and a port that generated them any other way compares against different numbers;
//!  * **the strategy cascade is three strategies over one loop**, hom alt above 0.25, hom ref down
//!    to 0.20, then the unscrupulous hom ref, and it returns the first estimate whose standard
//!    error is small enough RELATIVE TO THE ESTIMATE ITSELF.
//!
//! # The standard error is a binary search
//!
//! The closed formula for the standard error depends on the true contamination rather than on the
//! estimate, and at low contamination substituting the estimate gives nonsense: an estimate of zero
//! would have an error of zero. The reference instead searches for the contaminations consistent
//! with the estimate, and only falls back to the formula when the search fails to bracket a zero.
//!
//! # Every sum here is one of two different sums
//!
//! `MathUtils.sum` is a plain loop; a stream's `.sum()` is Kahan-compensated
//! ([`crate::allele_fraction_cluster::double_stream_sum`]). The reference uses both, and which one
//! it uses is not a detail: on a list of a few hundred depths they are different doubles.

use crate::allele_fraction_cluster::double_stream_sum;
use crate::contamination_segmenter::find_segments;
use crate::pileup_summary::PileupSummary;

/// `INITIAL_MAF_THRESHOLD`.
pub const INITIAL_MAF_THRESHOLD: f64 = 0.40;
/// `MAF_TO_SWITCH_TO_HOM_REF`.
pub const MAF_TO_SWITCH_TO_HOM_REF: f64 = 0.25;
/// `MAF_TO_SWITCH_TO_UNSCRUPULOUS_HOM_REF`.
pub const MAF_TO_SWITCH_TO_UNSCRUPULOUS_HOM_REF: f64 = 0.20;
/// `UNSCRUPULOUS_HOM_REF_ALLELE_FRACTION`.
pub const UNSCRUPULOUS_HOM_REF_ALLELE_FRACTION: f64 = 0.15;
/// `UNSCRUPULOUS_HOM_REF_FRACTION_TO_REMOVE_FOR_POSSIBLE_LOH`.
pub const UNSCRUPULOUS_HOM_REF_FRACTION_TO_REMOVE_FOR_POSSIBLE_LOH: f64 = 0.1;
/// `UNSCRUPULOUS_HOM_REF_PERCENTILE`, computed the reference's way rather than written as 90.
pub const UNSCRUPULOUS_HOM_REF_PERCENTILE: f64 =
    100.0 * (1.0 - UNSCRUPULOUS_HOM_REF_FRACTION_TO_REMOVE_FOR_POSSIBLE_LOH);
/// `MINIMUM_UNSCRUPULOUS_HOM_REF_ALT_FRACTION_THRESHOLD`.
pub const MINIMUM_UNSCRUPULOUS_HOM_REF_ALT_FRACTION_THRESHOLD: f64 = 0.1;
/// `MAF_STEP_SIZE`.
pub const MAF_STEP_SIZE: f64 = 0.04;
/// `PRECISION_FOR_STANDARD_ERROR`.
const PRECISION_FOR_STANDARD_ERROR: f64 = 0.000001;
/// `HOM_REF`, an index into the genotype likelihoods.
pub const HOM_REF: usize = 0;
/// `HOM_ALT`, the same.
pub const HOM_ALT: usize = 3;
/// `NUM_ITERATIONS`.
const NUM_ITERATIONS: usize = 3;
/// `MIN_FRACTION_OF_SITES_TO_USE`.
const MIN_FRACTION_OF_SITES_TO_USE: f64 = 0.25;
/// `MIN_RELATIVE_ERROR`.
const MIN_RELATIVE_ERROR: f64 = 0.2;
/// `MIN_ABSOLUTE_ERROR`.
const MIN_ABSOLUTE_ERROR: f64 = 0.001;
/// `CONTAMINATION_INITIAL_GUESSES`.
const CONTAMINATION_INITIAL_GUESSES: [f64; 4] = [0.02, 0.05, 0.1, 0.2];

/// `MathUtils.binomialProbability(n, k, p)`, which is `BinomialDistribution.probability(k)`.
///
/// The distribution's `logProbability` is the saddle point expansion, and `probability` exponentiates
/// it unless it is negative infinity, where it answers an exact zero rather than an underflow.
fn binomial_probability(n: i32, k: i32, p: f64) -> f64 {
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

/// `MathUtils.sum(double[])`, a plain loop and not the stream's compensated sum.
fn plain_sum(values: &[f64]) -> f64 {
    let mut total = 0.0;
    for value in values {
        total += value;
    }
    total
}

/// `MathUtils.binarySearchFindZero`.
///
/// Bisects while the bracket is wider than `precision`, and gives up the moment the two ends have
/// the same sign, which is the reference's way of saying the function is not monotone here.
fn binary_search_find_zero(
    function: impl Fn(f64) -> f64,
    lower: f64,
    upper: f64,
    precision: f64,
) -> Option<f64> {
    let mut bottom = lower;
    let mut top = upper;
    while top - bottom > precision {
        let mid = (bottom + top) / 2.0;
        let bottom_value = function(bottom);
        let top_value = function(top);
        let mid_value = function(mid);
        // `FastMath.signum`, which answers the argument itself for NaN and keeps signed zeros
        // distinct, so two zeros of opposite sign do not compare equal here.
        if signum(bottom_value) == signum(top_value) {
            return None;
        }
        if signum(bottom_value) == signum(mid_value) {
            bottom = mid;
        } else {
            top = mid;
        }
    }
    Some((bottom + top) / 2.0)
}

/// `FastMath.signum(double)`.
fn signum(value: f64) -> f64 {
    if value < 0.0 {
        -1.0
    } else if value > 0.0 {
        1.0
    } else {
        // Both zeros and NaN come back as themselves.
        value
    }
}

/// The three strategies `calculateContaminationFromHoms` walks through.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Strategy {
    /// Ref reads in hom alt sites, which is the principled one.
    HomAlt,
    /// Alt reads in hom ref sites, which is more exposed to the population frequencies.
    HomRef,
    /// Hom refs picked by a percentile rather than by the model.
    UnscrupulousHomRef,
}

/// The learned model.
#[derive(Debug, Clone)]
pub struct ContaminationModel {
    /// The contamination the three iterations converged on.
    pub contamination: f64,
    /// `calculateErrorRate` over the sites the model was built from.
    pub error_rate: f64,
    /// One minor allele fraction per segment.
    pub minor_allele_fractions: Vec<f64>,
    /// The segments themselves, each carrying its sites.
    pub segments: Vec<Vec<PileupSummary>>,
}

/// One segment's interval and fraction, which is `MinorAlleleFractionRecord`.
#[derive(Debug, Clone, PartialEq)]
pub struct MinorAlleleFractionRecord {
    /// The contig.
    pub contig: String,
    /// The first site's start.
    pub start: i32,
    /// The last site's end.
    pub end: i32,
    /// The fraction.
    pub minor_allele_fraction: f64,
}

impl ContaminationModel {
    /// The constructor: segment, then alternate between fractions and contamination three times.
    pub fn new(sites: &[PileupSummary]) -> ContaminationModel {
        let error_rate = calculate_error_rate(sites);
        let segments = find_segments(sites);

        let mut fractions = vec![0.5; segments.len()];
        let mut contamination = 0.0;
        for _ in 0..NUM_ITERATIONS {
            for (index, segment) in segments.iter().enumerate() {
                fractions[index] =
                    calculate_minor_allele_fraction(contamination, error_rate, segment);
            }
            let (non_loh_segments, non_loh_fractions) = non_loh_segments(&segments, &fractions);
            contamination =
                calculate_contamination(error_rate, &non_loh_segments, &non_loh_fractions);
        }

        ContaminationModel {
            contamination,
            error_rate,
            minor_allele_fractions: fractions,
            segments,
        }
    }

    /// `segmentationRecords`.
    pub fn segmentation_records(&self) -> Vec<MinorAlleleFractionRecord> {
        self.segments
            .iter()
            .enumerate()
            .map(|(index, segment)| MinorAlleleFractionRecord {
                contig: segment[0].contig.clone(),
                start: segment[0].position,
                end: segment[segment.len() - 1].position,
                minor_allele_fraction: self.minor_allele_fractions[index],
            })
            .collect()
    }

    /// `calculateContaminationFromHoms`: the estimate and its standard error.
    ///
    /// The loop runs the threshold down to zero inclusive, and the strategy is chosen from the
    /// threshold rather than from the previous answer. A NaN estimate never satisfies the exit
    /// condition, because the comparison against it is false whichever way it is written.
    pub fn calculate_contamination_from_homs(&self, tumor_sites: &[PileupSummary]) -> (f64, f64) {
        let mut min_maf = INITIAL_MAF_THRESHOLD;
        while min_maf >= 0.0 {
            let strategy = if min_maf > MAF_TO_SWITCH_TO_HOM_REF {
                Strategy::HomAlt
            } else if min_maf > MAF_TO_SWITCH_TO_UNSCRUPULOUS_HOM_REF {
                Strategy::HomRef
            } else {
                Strategy::UnscrupulousHomRef
            };
            let (estimate, error) = self.calculate_contamination(strategy, tumor_sites, min_maf);
            if !estimate.is_nan() && error < (estimate * MIN_RELATIVE_ERROR + MIN_ABSOLUTE_ERROR) {
                return (estimate, error);
            }
            min_maf -= MAF_STEP_SIZE;
        }

        // The last resort, and the only place a fixed answer is returned.
        let (estimate, error) =
            self.calculate_contamination(Strategy::UnscrupulousHomRef, tumor_sites, 0.0);
        if estimate.is_nan() {
            (0.0, 1.0)
        } else {
            (estimate, error)
        }
    }

    /// One strategy's estimate and standard error.
    fn calculate_contamination(
        &self,
        strategy: Strategy,
        tumor_sites: &[PileupSummary],
        min_maf: f64,
    ) -> (f64, f64) {
        let use_hom_alt = strategy == Strategy::HomAlt;
        let genotyping_homs: Vec<PileupSummary> = match strategy {
            Strategy::HomAlt => self.get_type(HOM_ALT, min_maf),
            Strategy::HomRef => self.get_type(HOM_REF, min_maf),
            Strategy::UnscrupulousHomRef => {
                let candidates: Vec<PileupSummary> = tumor_sites
                    .iter()
                    .filter(|site| site.alt_fraction() < UNSCRUPULOUS_HOM_REF_ALLELE_FRACTION)
                    .cloned()
                    .collect();
                let fractions: Vec<f64> =
                    candidates.iter().map(|site| site.alt_fraction()).collect();
                // `new Percentile(90).evaluate(...)`, which is the legacy estimation type.
                let percentile = jmath::percentile::evaluate(
                    &fractions,
                    UNSCRUPULOUS_HOM_REF_PERCENTILE,
                    jmath::percentile::EstimationType::Legacy,
                );
                let threshold = java_max(
                    MINIMUM_UNSCRUPULOUS_HOM_REF_ALT_FRACTION_THRESHOLD,
                    percentile,
                );
                candidates
                    .into_iter()
                    .filter(|site| site.alt_fraction() <= threshold)
                    .collect()
            }
        };

        let homs = subset_sites(tumor_sites, &genotyping_homs);
        let tumor_error_rate = calculate_error_rate(tumor_sites);

        // Ref depth in hom alts, or alt depth in hom refs.
        let opposite_count = |site: &PileupSummary| {
            if use_hom_alt {
                site.ref_count
            } else {
                site.alt_count
            }
        };
        let opposite_frequency = |site: &PileupSummary| {
            if use_hom_alt {
                site.ref_frequency()
            } else {
                site.allele_frequency
            }
        };

        let total_depth: i64 = homs.iter().map(|site| i64::from(site.total_count)).sum();
        let opposite_depth: i64 = homs
            .iter()
            .map(|site| i64::from(opposite_count(site)))
            .sum();
        // `Math.round(double)`, which is `floor(x + 0.5)` and not a round-half-even.
        let error_depth = java_round(total_depth as f64 * tumor_error_rate / 3.0);
        let contamination_opposite_depth = (opposite_depth - error_depth).max(0);

        let weighted: Vec<f64> = homs
            .iter()
            .map(|site| f64::from(site.total_count) * opposite_frequency(site))
            .collect();
        let total_depth_weighted = double_stream_sum(&weighted);

        let contamination_estimate = contamination_opposite_depth as f64 / total_depth_weighted;

        let coefficient_one = double_stream_sum(
            &homs
                .iter()
                .map(|site| opposite_frequency(site) * f64::from(site.total_count))
                .collect::<Vec<f64>>(),
        );
        let coefficient_two = double_stream_sum(
            &homs
                .iter()
                .map(|site| {
                    let frequency = opposite_frequency(site);
                    frequency * (1.0 - frequency) * square(f64::from(site.total_count))
                })
                .collect::<Vec<f64>>(),
        );

        let empty = homs.is_empty();
        let error_function = |c: f64| {
            if empty {
                1.0
            } else {
                (coefficient_one * c * (1.0 - c) + coefficient_two * c * c).sqrt()
                    / total_depth_weighted
            }
        };

        let upper = binary_search_find_zero(
            |c| c - error_function(c) - contamination_estimate,
            contamination_estimate,
            1.0,
            PRECISION_FOR_STANDARD_ERROR,
        );
        let lower = binary_search_find_zero(
            |c| c + error_function(c) - contamination_estimate,
            0.0,
            contamination_estimate,
            PRECISION_FOR_STANDARD_ERROR,
        );

        let standard_error = match (upper, lower) {
            (Some(upper), Some(lower)) => java_max(
                upper - contamination_estimate,
                contamination_estimate - lower,
            ),
            // The reference says this should never happen, and falls back to the closed formula.
            _ => error_function(contamination_estimate),
        };

        (java_min(contamination_estimate, 1.0), standard_error)
    }

    /// `getType`: the sites in non-LoH segments whose posterior for one genotype is above a half.
    fn get_type(&self, genotype: usize, min_maf: f64) -> Vec<PileupSummary> {
        let indices: Vec<usize> = (0..self.segments.len())
            .filter(|index| self.minor_allele_fractions[*index] > min_maf)
            .collect();
        let mut sites = Vec::new();
        for index in indices {
            for site in &self.segments[index] {
                if probability(
                    site,
                    self.contamination,
                    self.error_rate,
                    self.minor_allele_fractions[index],
                    genotype,
                ) > 0.5
                {
                    sites.push(site.clone());
                }
            }
        }
        sites
    }
}

/// `Math.max(double, double)`, which propagates NaN where `f64::max` returns the other argument.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        return f64::NAN;
    }
    if a >= b {
        a
    } else {
        b
    }
}

/// `Math.min(double, double)`, for the same reason.
fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        return f64::NAN;
    }
    if a <= b {
        a
    } else {
        b
    }
}

/// `Math.round(double)`: `floor(x + 0.5)`, so a half always goes up rather than to even.
fn java_round(value: f64) -> i64 {
    (value + 0.5).floor() as i64
}

/// `MathUtils.square`.
fn square(value: f64) -> f64 {
    value * value
}

/// `subsetSites`: the sites that overlap any of the subset's loci.
///
/// The reference builds an `OverlapDetector` over the subset and asks it, which on single-base
/// pileup summaries is the same as asking whether contig and position both match.
fn subset_sites(sites: &[PileupSummary], subset: &[PileupSummary]) -> Vec<PileupSummary> {
    sites
        .iter()
        .filter(|site| {
            subset
                .iter()
                .any(|locus| locus.contig == site.contig && locus.position == site.position)
        })
        .cloned()
        .collect()
}

/// `calculateErrorRate`: one and a half times the other-alt fraction of every base.
///
/// The counts are summed as longs and divided as doubles, so a site with no coverage at all gives
/// a NaN rather than a zero.
pub fn calculate_error_rate(sites: &[PileupSummary]) -> f64 {
    let total_bases: i64 = sites.iter().map(|site| i64::from(site.total_count)).sum();
    let other_alt_bases: i64 = sites
        .iter()
        .map(|site| i64::from(site.other_alt_count))
        .sum();
    1.5 * (other_alt_bases as f64 / total_bases as f64)
}

/// `calculateMinorAlleleFraction`: Brent, between 0.1 and 0.5, started at 0.4.
fn calculate_minor_allele_fraction(
    contamination: f64,
    error_rate: f64,
    segment: &[PileupSummary],
) -> f64 {
    jmath::brent::maximize(
        |maf| segment_log_likelihood(segment, contamination, error_rate, maf),
        0.1,
        0.5,
        0.4,
        0.01,
        0.01,
        20,
    )
    .expect("the reference propagates a Brent failure rather than handling it")
    .point
}

/// `calculateContamination`: four starting points, and the best of the four optima.
///
/// `Collections.max` keeps the FIRST of equal values, so two starting points that converge to the
/// same likelihood are decided by the order of the guesses.
fn calculate_contamination(
    error_rate: f64,
    segments: &[Vec<PileupSummary>],
    fractions: &[f64],
) -> f64 {
    let optima: Vec<jmath::brent::PointValuePair> = CONTAMINATION_INITIAL_GUESSES
        .iter()
        .map(|initial| {
            jmath::brent::maximize(
                |c| model_log_likelihood(segments, c, error_rate, fractions),
                0.0,
                0.5,
                *initial,
                1.0e-4,
                1.0e-4,
                30,
            )
            .expect("the reference propagates a Brent failure rather than handling it")
        })
        .collect();

    let mut best = &optima[0];
    for candidate in &optima[1..] {
        if candidate.value > best.value {
            best = candidate;
        }
    }
    best.point
}

/// `getNonLOHSegments`: walk the threshold down until enough sites survive.
fn non_loh_segments(
    segments: &[Vec<PileupSummary>],
    fractions: &[f64],
) -> (Vec<Vec<PileupSummary>>, Vec<f64>) {
    let num_sites: usize = segments.iter().map(|segment| segment.len()).sum();
    let mut min_maf = INITIAL_MAF_THRESHOLD;
    while min_maf > 0.0 {
        let indices: Vec<usize> = (0..segments.len())
            .filter(|index| fractions[*index] > min_maf)
            .collect();
        let kept: usize = indices.iter().map(|index| segments[*index].len()).sum();
        if kept as f64 / num_sites as f64 > MIN_FRACTION_OF_SITES_TO_USE {
            return (
                indices
                    .iter()
                    .map(|index| segments[*index].clone())
                    .collect(),
                indices.iter().map(|index| fractions[*index]).collect(),
            );
        }
        min_maf -= MAF_STEP_SIZE;
    }
    (segments.to_vec(), fractions.to_vec())
}

/// `genotypeLikelihoods`: hom ref, alt minor, alt major, hom alt.
fn genotype_likelihoods(
    site: &PileupSummary,
    contamination: f64,
    error_rate: f64,
    minor_allele_fraction: f64,
) -> [f64; 4] {
    let f = site.allele_frequency;
    let k = site.alt_count;
    let n = k + site.ref_count;

    let priors = [(1.0 - f) * (1.0 - f), f * (1.0 - f), f * (1.0 - f), f * f];
    let allele_fractions = [
        error_rate / 3.0,
        minor_allele_fraction,
        1.0 - minor_allele_fraction,
        1.0 - error_rate,
    ];

    let mut likelihoods = [0.0; 4];
    for genotype in 0..4 {
        likelihoods[genotype] = priors[genotype]
            * binomial_probability(
                n,
                k,
                (1.0 - contamination) * allele_fractions[genotype] + contamination * f,
            );
    }
    likelihoods
}

/// `probability`: one genotype's share of the four likelihoods.
fn probability(
    site: &PileupSummary,
    contamination: f64,
    error_rate: f64,
    minor_allele_fraction: f64,
    genotype: usize,
) -> f64 {
    let likelihoods = genotype_likelihoods(site, contamination, error_rate, minor_allele_fraction);
    likelihoods[genotype] / plain_sum(&likelihoods)
}

/// `segmentLogLikelihood`, whose outer sum is a stream's and whose inner sum is not.
fn segment_log_likelihood(
    segment: &[PileupSummary],
    contamination: f64,
    error_rate: f64,
    minor_allele_fraction: f64,
) -> f64 {
    let terms: Vec<f64> = segment
        .iter()
        .map(|site| {
            plain_sum(&genotype_likelihoods(
                site,
                contamination,
                error_rate,
                minor_allele_fraction,
            ))
            .ln()
        })
        .collect();
    double_stream_sum(&terms)
}

/// `modelLogLikelihood`, whose sum is `IndexRange.sum` and therefore a plain loop.
fn model_log_likelihood(
    segments: &[Vec<PileupSummary>],
    contamination: f64,
    error_rate: f64,
    fractions: &[f64],
) -> f64 {
    let mut total = 0.0;
    for (index, segment) in segments.iter().enumerate() {
        total += segment_log_likelihood(segment, contamination, error_rate, fractions[index]);
    }
    total
}
