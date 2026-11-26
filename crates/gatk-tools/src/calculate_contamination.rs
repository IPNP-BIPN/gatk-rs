//! Ported from `org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination`
//! (GATK 4.7.0.0; `CalculateContamination.java` is unchanged since 4.6.2.0).
//!
//! The tool around [`gatk_engine::contamination_model`]. It is small, and everything in it that is
//! not bookkeeping is the coverage filter.
//!
//! # The coverage filter uses a median and a mean, and both are computed twice over
//!
//! The low threshold is a fraction of the MEDIAN coverage and the high threshold a multiple of the
//! MEAN, so a table with a long tail of deep sites is cut differently at each end. Both statistics
//! are computed AFTER dropping sites at or below `MIN_COVERAGE`, so an uncovered site that would
//! have dragged the median down is already gone when the median is taken.
//!
//! Both thresholds are strict: a site exactly at the low threshold is dropped, and so is one
//! exactly at the high threshold.
//!
//! # The two ratio arguments do nothing
//!
//! `--low-coverage-ratio-threshold` and `--high-coverage-ratio-threshold` are accepted and ignored
//! by the reference: their fields are constant variables under JLS 4.12.4, so the one read of each
//! is a compile-time constant and the parsed value is never looked at. [`run_from_command_line`]
//! reproduces that; [`run`] is the function a caller with its own thresholds wants.
//!
//! # Which model genotypes and which model segments are two questions
//!
//! With a matched normal the NORMAL's model genotypes the tumour, while the segmentation table, if
//! one is asked for, comes from a second model built on the TUMOUR. Without one, a single model
//! does both. That is why the tool can build two models over the same run.
//!
//! Java quotes are verbatim apart from line breaks added to fit 100 columns.
//!
//! # Reading this file without reading Rust
//!
//! - `&[PileupSummary]` is a borrowed, read-only view of a list; the caller keeps ownership.
//! - `Option<T>` is "a value or nothing": `None` is Java's `null` for an absent matched normal or
//!   segmentation, checked by the compiler.
//! - `.iter().filter(..).cloned().collect()` is Java's `stream().filter(..).collect(toList())`;
//!   `.cloned()` copies each element, because a filtered view would otherwise borrow the input.
//! - `match x { None => .., Some(v) => .. }` is a null check that must handle both cases.
//! - A parameter whose name starts with `_` is accepted and deliberately unused.

use gatk_engine::contamination_model::{ContaminationModel, MinorAlleleFractionRecord};
use gatk_engine::pileup_summary::PileupSummary;

/// `MIN_COVERAGE`.
pub const MIN_COVERAGE: i32 = 10;
/// `DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD`, written as the reference writes it.
pub const DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD: f64 = 1.0 / 2.0;
/// `DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD`.
pub const DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD: f64 = 3.0;

/// What the tool answers: the contamination record and, optionally, the segmentation table.
#[derive(Debug, Clone, PartialEq)]
pub struct Output {
    /// The contamination estimate.
    pub contamination: f64,
    /// Its standard error.
    pub error: f64,
    /// The tumour segmentation, present when one was asked for.
    pub segmentation: Option<Vec<MinorAlleleFractionRecord>>,
}

/// `Mean.evaluate(double[])`, which is a corrected two-pass mean and not a sum over a count.
///
/// The first pass is a plain `Sum`, the second accumulates the residuals against that estimate, and
/// the correction is added back. On a few hundred depths it differs from `sum / n` in the last
/// bits, which is enough to move a threshold and therefore a site.
pub fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let sample_size = values.len() as f64;
    let mut sum = 0.0;
    for value in values {
        sum += value;
    }
    let xbar = sum / sample_size;
    let mut correction = 0.0;
    for value in values {
        correction += value - xbar;
    }
    xbar + (correction / sample_size)
}

/// `filterSitesByCoverage`.
///
/// ```java
/// final double lowCoverageThreshold = medianCoverage * lowCoverageRatioThreshold;
/// final double highCoverageThreshold = meanCoverage * highCoverageRatioThreshold;
/// return coveredSites.stream()
///         .filter(ps -> ps.getTotalCount() > lowCoverageThreshold
///                 && ps.getTotalCount() < highCoverageThreshold)
///         .collect(Collectors.toList());
/// ```
///
/// The ratios are parameters here so [`run`] can take them; see [`run_from_command_line`] for why
/// the tool itself always passes the defaults.
pub fn filter_sites_by_coverage(
    all_sites: &[PileupSummary],
    low_coverage_ratio_threshold: f64,
    high_coverage_ratio_threshold: f64,
) -> Vec<PileupSummary> {
    // "Just in case the intervals given to GetPileupSummaries contained un-covered sites": these go
    // before the median so that a run of zeroes cannot move it.
    let covered: Vec<PileupSummary> = all_sites
        .iter()
        .filter(|site| site.total_count > MIN_COVERAGE)
        .cloned()
        .collect();
    let coverage: Vec<f64> = covered
        .iter()
        .map(|site| f64::from(site.total_count))
        .collect();
    // `new Median().evaluate(coverage)` and `new Mean().evaluate(coverage)`, commons-math3's.
    let median_coverage = jmath::percentile::median(&coverage);
    let mean_coverage = mean(&coverage);
    let low = median_coverage * low_coverage_ratio_threshold;
    let high = mean_coverage * high_coverage_ratio_threshold;
    // `into_iter()` consumes `covered`, so the kept sites move into the result without a copy.
    covered
        .into_iter()
        .filter(|site| {
            let count = f64::from(site.total_count);
            count > low && count < high
        })
        .collect()
}

/// `doWork` as the COMMAND LINE reaches it: the two ratio arguments are accepted and IGNORED.
///
/// `--low-coverage-ratio-threshold` and `--high-coverage-ratio-threshold` are declared, parsed and
/// written into their fields, and then never read. Both are `private final double` initialised
/// from a constant expression, which makes them constant variables under JLS 4.12.4 and their one
/// read inside `filterSitesByCoverage` a compile-time constant. A probe reading the fields back
/// after the parse shows the value that was asked for, and the answer does not move: a low ratio
/// of ten, which would drop every site in the table, changes nothing.
///
/// The port reproduces the reference, so this entry point takes the two ratios and drops them.
/// [`run`] is what a caller with its own thresholds wants, and no command line reaches it.
pub fn run_from_command_line(
    sites: &[PileupSummary],
    matched: Option<&[PileupSummary]>,
    segmentation: bool,
    _low_coverage_ratio_threshold: f64,
    _high_coverage_ratio_threshold: f64,
) -> Output {
    run(
        sites,
        matched,
        segmentation,
        DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD,
        DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD,
    )
}

/// `doWork`, on tables already read.
///
/// ```java
/// final ContaminationModel tumorModel =
///         matchedPileupSummariesTable == null ? genotypingModel : new ContaminationModel(sites);
/// ```
///
///
/// `matched` is the matched normal's sites, if there are any. `segmentation` asks for the tumour
/// segmentation table, which is what decides whether the second model is built at all.
pub fn run(
    sites: &[PileupSummary],
    matched: Option<&[PileupSummary]>,
    segmentation: bool,
    low_coverage_ratio_threshold: f64,
    high_coverage_ratio_threshold: f64,
) -> Output {
    let filtered = filter_sites_by_coverage(
        sites,
        low_coverage_ratio_threshold,
        high_coverage_ratio_threshold,
    );
    // The matched normal is filtered the same way, and it genotypes when it is there.
    let genotyping_sites = match matched {
        None => filtered.clone(),
        Some(matched) => filter_sites_by_coverage(
            matched,
            low_coverage_ratio_threshold,
            high_coverage_ratio_threshold,
        ),
    };

    let genotyping_model = ContaminationModel::new(&genotyping_sites);

    let segmentation = if segmentation {
        // Without a matched normal the genotyping model IS the tumour model, so it is not rebuilt.
        let records = match matched {
            None => genotyping_model.segmentation_records(),
            Some(_) => ContaminationModel::new(&filtered).segmentation_records(),
        };
        Some(records)
    } else {
        None
    };

    // `calculateContaminationFromHoms(sites)` returns a `Pair`; here a tuple, split into two names.
    let (contamination, error) = genotyping_model.calculate_contamination_from_homs(&filtered);
    Output {
        contamination,
        error,
        segmentation,
    }
}
