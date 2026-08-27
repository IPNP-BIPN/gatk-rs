//! `ModelSegments`' genotyping and segmentation, the two steps the multi-sample mode does.
//!
//! Ported from `org.broadinstitute.hellbender.tools.copynumber.ModelSegments`,
//! `org.broadinstitute.hellbender.tools.copynumber.utils.genotyping.NaiveHeterozygousPileupGenotypingUtils`
//! and `org.broadinstitute.hellbender.tools.copynumber.segmentation.MultisampleMultidimensionalKernelSegmenter`.
//!
//! The kernel segmenter itself is not here: it draws a subsample through a seeded
//! `java.util.Random`, which may not be transcribed. What is here is everything around it, which
//! is what decides the answer once the changepoints are known: which sites reach the genotyper,
//! which of them are called heterozygous, how many changepoints the cap allows, and how a list of
//! changepoint indices becomes a list of intervals.

use jmath::beta::{regularized_beta, BetaError};

/// `SomaticGenotypingArgumentCollection.minTotalAlleleCountCase`, which lets every site through.
pub const DEFAULT_MINIMUM_TOTAL_ALLELE_COUNT_CASE: i32 = 0;
/// `SomaticGenotypingArgumentCollection.genotypingHomozygousLogRatioThreshold`.
pub const DEFAULT_GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD: f64 = -10.0;
/// `SomaticGenotypingArgumentCollection.genotypingBaseErrorRate`.
pub const DEFAULT_GENOTYPING_BASE_ERROR_RATE: f64 = 5E-2;
/// `MultisampleMultidimensionalKernelSegmenter.MIN_NUM_POINTS_REQUIRED_PER_CHROMOSOME`.
pub const MINIMUM_POINTS_REQUIRED_PER_CHROMOSOME: usize = 10;
/// `MultisampleMultidimensionalKernelSegmenter.findSegmentation`, on intervals that disagree.
pub const MISMATCHED_COPY_RATIO_INTERVALS_MESSAGE: &str =
    "Copy-ratio intervals must be identical across all case samples.";

/// One pileup's counts at a site, which is an `AllelicCount`'s two numbers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AllelicCount {
    pub position: i32,
    pub reference_count: i32,
    pub alternate_count: i32,
}

impl AllelicCount {
    /// `AllelicCount.getTotalReadCount`.
    pub fn total_read_count(&self) -> i32 {
        self.reference_count + self.alternate_count
    }
}

/// `calculateHomozygousLogRatio`: the log of the ratio of the likelihood that the success
/// probability lies in the two error-rate tails to the likelihood that it lies between them.
///
/// The reference asks `Beta.regularizedBeta(1, r + 1, n - r + 1)` for the whole mass, which is one
/// up to the tolerance, and subtracts the two partial integrals from it rather than from one.
pub fn homozygous_log_ratio(count: AllelicCount, base_error_rate: f64) -> Result<f64, BetaError> {
    let r = f64::from(count.reference_count);
    let n = f64::from(count.total_read_count());
    let beta_all = regularized_beta(1.0, r + 1.0, n - r + 1.0)?;
    let beta_error = regularized_beta(base_error_rate, r + 1.0, n - r + 1.0)?;
    let beta_one_minus_error = regularized_beta(1.0 - base_error_rate, r + 1.0, n - r + 1.0)?;
    let beta_homozygous = beta_error + beta_all - beta_one_minus_error;
    let beta_heterozygous = beta_one_minus_error - beta_error;
    Ok(beta_homozygous.ln() - beta_heterozygous.ln())
}

/// `filterByHeterozygosity`'s predicate, which is STRICTLY below the threshold.
pub fn is_heterozygous(
    count: AllelicCount,
    homozygous_log_ratio_threshold: f64,
    base_error_rate: f64,
) -> Result<bool, BetaError> {
    Ok(homozygous_log_ratio(count, base_error_rate)? < homozygous_log_ratio_threshold)
}

/// `filterByTotalCount`, which keeps a site whose total reaches the floor.
///
/// A floor of zero returns the collection untouched rather than filtering it, which is the same
/// answer because no total is negative.
pub fn filter_by_total_count(counts: &[AllelicCount], minimum: i32) -> Vec<AllelicCount> {
    if minimum == 0 {
        return counts.to_vec();
    }
    counts
        .iter()
        .copied()
        .filter(|count| count.total_read_count() >= minimum)
        .collect()
}

/// `filterByOverlap`, whose EMPTY collection keeps nothing rather than everything.
pub fn filter_by_overlap(counts: &[AllelicCount], intervals: &[(i32, i32)]) -> Vec<AllelicCount> {
    if intervals.is_empty() {
        return Vec::new();
    }
    counts
        .iter()
        .copied()
        .filter(|count| {
            intervals
                .iter()
                .any(|(start, end)| *start <= count.position && count.position <= *end)
        })
        .collect()
}

/// `genotypeHets` in case-only mode, for one sample: the total-count floor, then the overlap with
/// the copy-ratio intervals, then the heterozygosity test, in that order.
pub fn genotype_hets(
    counts: &[AllelicCount],
    copy_ratio_intervals: &[(i32, i32)],
    minimum_total_allele_count: i32,
    homozygous_log_ratio_threshold: f64,
    base_error_rate: f64,
) -> Result<Vec<AllelicCount>, BetaError> {
    let mut filtered = filter_by_total_count(counts, minimum_total_allele_count);
    if !copy_ratio_intervals.is_empty() {
        filtered = filter_by_overlap(&filtered, copy_ratio_intervals);
    }
    let mut hets = Vec::new();
    for count in filtered {
        if is_heterozygous(count, homozygous_log_ratio_threshold, base_error_rate)? {
            hets.push(count);
        }
    }
    Ok(hets)
}

/// The intersection `genotypeHets` takes over the samples, which it reaches only for MORE than
/// one: the sites called heterozygous in every sample, in the first sample's order.
pub fn intersect_het_sites(per_sample: &[Vec<AllelicCount>]) -> Vec<i32> {
    if per_sample.len() <= 1 {
        return per_sample
            .first()
            .map(|counts| counts.iter().map(|count| count.position).collect())
            .unwrap_or_default();
    }
    per_sample[0]
        .iter()
        .map(|count| count.position)
        .filter(|position| {
            per_sample
                .iter()
                .all(|counts| counts.iter().any(|count| count.position == *position))
        })
        .collect()
}

/// `ImmutableSet.copyOf(windowSizes).asList()`: the window sizes are a SET, so a size named twice
/// segments exactly as it does named once, and the first occurrence decides the order.
///
/// The segmenter itself refuses a repeated size, which is why the tool deduplicates before calling
/// it rather than passing the argument list through.
pub fn window_sizes(named: &[i32]) -> Vec<i32> {
    let mut unique: Vec<i32> = Vec::new();
    for size in named {
        if !unique.contains(size) {
            unique.push(*size);
        }
    }
    unique
}

/// `maxNumSegmentsPerChromosome - 1`: the cap counts SEGMENTS, and one changepoint makes two.
pub fn maximum_changepoints_per_chromosome(maximum_segments_per_chromosome: i32) -> i32 {
    maximum_segments_per_chromosome - 1
}

/// The segments `findSegmentation` builds from the changepoint indices of one chromosome.
///
/// Each changepoint is the LAST index of its segment, so a segment runs from the point after the
/// previous changepoint to the changepoint itself. A chromosome with too few points is not
/// segmented at all: it becomes the one interval spanning its points.
///
/// The closing index is appended when the list does not already hold the point COUNT, which is a
/// number no changepoint index can be: the indices run to the count less one. The guard is
/// therefore always true and the closing index is always appended, which is what makes the last
/// segment close on the final point. Handing the final index in as a changepoint as well would
/// have it appended a second time and then read one point past the end.
pub fn segments_from_changepoints(
    points: &[(i32, i32)],
    changepoints: &[usize],
) -> Vec<(i32, i32)> {
    if points.is_empty() {
        return Vec::new();
    }
    if points.len() < MINIMUM_POINTS_REQUIRED_PER_CHROMOSOME {
        return vec![(points[0].0, points[points.len() - 1].1)];
    }
    let mut closed: Vec<usize> = changepoints.to_vec();
    if !closed.contains(&points.len()) {
        closed.push(points.len() - 1);
    }
    let mut segments = Vec::with_capacity(closed.len());
    let mut previous: i64 = -1;
    for changepoint in closed {
        let start = points[(previous + 1) as usize].0;
        let end = points[changepoint].1;
        segments.push((start, end));
        previous = changepoint as i64;
    }
    segments
}
