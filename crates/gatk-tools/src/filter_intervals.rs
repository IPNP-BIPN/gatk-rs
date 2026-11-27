//! Ported from `org.broadinstitute.hellbender.tools.copynumber.FilterIntervals` (GATK 4.7.0.0;
//! `FilterIntervals.java` is unchanged since 4.6.2.0).
//!
//! The intervals a copy-number panel may use. Three filters over one shared mask, and the answer
//! depends on the order they run in as much as on the thresholds.
//!
//! # A contig left with a single interval loses it
//!
//! After every other filter the tool counts the survivors PER CONTIG and removes any contig's only
//! one, with a warning. So a run that filters down to exactly one interval ends with NONE, and the
//! count check then refuses the run outright. It is the last thing that happens and the easiest to
//! miss: `filter-intervals`' `solitary` case exists because three probes were needed to find it.
//!
//! # The mask is shared, so the order is part of the arithmetic
//!
//! An interval failed by the GC filter is skipped by the count filters -- and, more to the point,
//! is not in the population the per-sample percentiles are computed over. Reordering the filters
//! would change the percentiles and therefore the answer.
//!
//! # Every bound here is inclusive, and the two count rules are not
//!
//! The annotation ranges and the percentile range are `min <= x && x <= max`, tested by negation so
//! that a NaN fails. The low-count and extreme-count rules are STRICTLY GREATER on the count of
//! offending samples, so at fifty per cent one sample of two is not enough.
//!
//! # Reading this file without reading Rust
//!
//! - The mask is `&mut [bool]`: a writable borrow of the caller's `boolean[] mask`, so each
//!   filter updates the one array in place, as the reference's methods do. `&[f64]` and
//!   `&[Vec<f64>]` are read-only borrows of a `double[]` and of a `double[][]`.
//! - `Result<usize, FilterError>` is "a count or the refusal"; `Err(..)` is the reference's
//!   `throw`, and the caller passes it on with `?`.
//! - `(a..b)` is the half-open range `a, a + 1, ..., b - 1`, Java's `IntStream.range(a, b)`;
//!   `.filter(..).count()` is the stream's `filter(..).count()`, and `|x| ...` is a lambda.
//! - `as f64` converts a count to `double`, which Java does implicitly in the comparison.
//! - `#[allow(clippy::neg_cmp_op_on_partial_ord)]` silences a lint that would rewrite
//!   `!(a <= b)` as `a > b`, which is not the same test when a NaN is involved.

/// `--minimum-gc-content`'s default.
pub const DEFAULT_MINIMUM_GC_CONTENT: f64 = 0.1;
/// `--maximum-gc-content`'s default.
pub const DEFAULT_MAXIMUM_GC_CONTENT: f64 = 0.9;
/// `--minimum-mappability`'s default.
pub const DEFAULT_MINIMUM_MAPPABILITY: f64 = 0.9;
/// `--maximum-mappability`'s default.
pub const DEFAULT_MAXIMUM_MAPPABILITY: f64 = 1.0;
/// `--minimum-segmental-duplication-content`'s default.
pub const DEFAULT_MINIMUM_SEGMENTAL_DUPLICATION_CONTENT: f64 = 0.0;
/// `--maximum-segmental-duplication-content`'s default.
pub const DEFAULT_MAXIMUM_SEGMENTAL_DUPLICATION_CONTENT: f64 = 0.5;
/// `--low-count-filter-count-threshold`'s default.
pub const DEFAULT_LOW_COUNT_THRESHOLD: i32 = 10;
/// `--low-count-filter-percentage-of-samples`'s default.
pub const DEFAULT_LOW_COUNT_PERCENTAGE: f64 = 50.0;
/// `--extreme-count-filter-minimum-percentile`'s default.
pub const DEFAULT_EXTREME_MINIMUM_PERCENTILE: f64 = 1.0;
/// `--extreme-count-filter-maximum-percentile`'s default.
pub const DEFAULT_EXTREME_MAXIMUM_PERCENTILE: f64 = 99.0;
/// `--extreme-count-filter-percentage-of-samples`'s default.
pub const DEFAULT_EXTREME_PERCENTAGE: f64 = 90.0;

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FilterError {
    /// Neither annotations nor counts were given.
    NoInputs,
    /// The intersection of the requested intervals with the inputs' was empty.
    EmptyIntersection,
    /// Every interval was filtered out.
    EverythingFiltered,
}

impl FilterError {
    /// The exception class the reference throws.
    pub fn java_class(&self) -> &'static str {
        match self {
            FilterError::NoInputs => "org.broadinstitute.hellbender.exceptions.UserException",
            FilterError::EmptyIntersection => "java.lang.IllegalArgumentException",
            FilterError::EverythingFiltered => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
        }
    }

    /// The message, with the `Bad input: ` prefix where the exception adds one. Note the TWO
    /// spaces in the last: the reference's own string has them.
    pub fn message(&self) -> &'static str {
        match self {
            FilterError::NoInputs => "Must provide annotated intervals or counts.",
            FilterError::EmptyIntersection => "At least one interval must remain after intersection.",
            FilterError::EverythingFiltered => {
                "Bad input: Filtering removed all intervals.  Select less strict filtering criteria."
            }
        }
    }
}

/// `countNumberPassing`, which THROWS when nothing passes rather than answering zero.
///
/// The reference's line, broken here to fit:
///
/// ```java
/// throw new UserException.BadInput(
///         "Filtering removed all intervals.  Select less strict filtering criteria.");
/// ```
///
/// The reference calls it inside the log message after every filter, so a filter that removes
/// the last interval refuses the run there, before any later filter runs; every `update_mask_by_`
/// function below ends with it for that reason.
fn count_passing(mask: &[bool]) -> Result<usize, FilterError> {
    let passing = mask.iter().filter(|filtered| !**filtered).count();
    if passing == 0 {
        return Err(FilterError::EverythingFiltered);
    }
    Ok(passing)
}

/// `updateMaskByAnnotationFilter`: fail every interval outside `[minimum, maximum]`.
///
/// The test is written as a negation, which is what makes a NaN annotation fail rather than pass:
/// `if (!(minValue <= value && value <= maxValue)) {`. An interval already failed is skipped, as
/// the reference's `.filter(i -> !mask[i])` skips it.
pub fn update_mask_by_annotation(
    mask: &mut [bool],
    values: &[f64],
    minimum: f64,
    maximum: f64,
) -> Result<usize, FilterError> {
    for (index, value) in values.iter().enumerate() {
        if mask[index] {
            continue;
        }
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(minimum <= *value && *value <= maximum) {
            mask[index] = true;
        }
    }
    count_passing(mask)
}

/// The low-count filter: fail an interval whose count is below the threshold in too many samples.
///
/// `counts[sample][interval]`. The comparison is STRICTLY GREATER against
/// `percentage * numSamples / 100`, so one sample of two at fifty per cent does not fail it.
///
/// "Below" is strictly less: `.filter(c -> c < lowCountFilterCountThreshold)`, with the `int`
/// threshold widened to `double` (`f64::from` here).
pub fn update_mask_by_low_counts(
    mask: &mut [bool],
    counts: &[Vec<f64>],
    threshold: i32,
    percentage: f64,
) -> Result<usize, FilterError> {
    let samples = counts.len();
    for index in 0..mask.len() {
        if mask[index] {
            continue;
        }
        let below = counts
            .iter()
            .filter(|row| row[index] < f64::from(threshold))
            .count();
        if below as f64 > percentage * samples as f64 / 100.0 {
            mask[index] = true;
        }
    }
    count_passing(mask)
}

/// The extreme-count filter: per sample, fail the intervals outside the percentile band.
///
/// The percentiles are taken over the SURVIVORS at this point, per sample, and a percentile of zero
/// is short-circuited to a threshold of zero rather than evaluated -- which is what makes
/// `--extreme-count-filter-minimum-percentile 0` mean "no lower bound".
///
/// The percentile is commons-math3's `new Percentile(p).evaluate(filteredCounts)`, whose default
/// estimation is `LEGACY`; `jmath::percentile` with `EstimationType::Legacy` stands for it.
///
/// Every interval of every sample is then tested against that sample's band, already failed or
/// not, as the reference's inner loop over all `numIntersectedIntervals` does; only the vote
/// that follows skips the intervals already failed.
pub fn update_mask_by_extreme_counts(
    mask: &mut [bool],
    counts: &[Vec<f64>],
    minimum_percentile: f64,
    maximum_percentile: f64,
    percentage: f64,
) -> Result<usize, FilterError> {
    let samples = counts.len();
    let intervals = mask.len();
    // `new boolean[numSamples][numIntersectedIntervals]`: `vec![x; n]` is an array of `n`
    // copies of `x`.
    let mut percentile_mask = vec![vec![false; intervals]; samples];

    for (sample, row) in counts.iter().enumerate() {
        let surviving: Vec<f64> = (0..intervals)
            .filter(|index| !mask[*index])
            .map(|index| row[index])
            .collect();
        let lower = if minimum_percentile == 0.0 {
            0.0
        } else {
            jmath::percentile::evaluate(
                &surviving,
                minimum_percentile,
                jmath::percentile::EstimationType::Legacy,
            )
        };
        let upper = if maximum_percentile == 0.0 {
            0.0
        } else {
            jmath::percentile::evaluate(
                &surviving,
                maximum_percentile,
                jmath::percentile::EstimationType::Legacy,
            )
        };
        for index in 0..intervals {
            let count = row[index];
            #[allow(clippy::neg_cmp_op_on_partial_ord)]
            if !(lower <= count && count <= upper) {
                percentile_mask[sample][index] = true;
            }
        }
    }

    for index in 0..intervals {
        if mask[index] {
            continue;
        }
        let offending = (0..samples)
            .filter(|sample| percentile_mask[*sample][index])
            .count();
        if offending as f64 > percentage * samples as f64 / 100.0 {
            mask[index] = true;
        }
    }
    count_passing(mask)
}

/// The last filter: a contig left with exactly one surviving interval loses it.
///
/// The reference counts the survivors per contig with
/// `Collectors.groupingBy(SimpleInterval::getContig, Collectors.counting())` and then fails each
/// survivor whose contig's count is one, with a warning. All counts are taken before any interval
/// is failed, so failing one cannot change another's count. The list of pairs stands for the
/// map; only lookups are made in it, so its order does not matter.
pub fn update_mask_by_solitary_intervals(
    mask: &mut [bool],
    contigs: &[String],
) -> Result<usize, FilterError> {
    let mut counts: Vec<(String, usize)> = Vec::new();
    for (index, contig) in contigs.iter().enumerate() {
        if mask[index] {
            continue;
        }
        match counts.iter_mut().find(|(name, _)| name == contig) {
            Some((_, total)) => *total += 1,
            None => counts.push((contig.clone(), 1)),
        }
    }
    for (index, contig) in contigs.iter().enumerate() {
        if mask[index] {
            continue;
        }
        let total = counts
            .iter()
            .find(|(name, _)| name == contig)
            .map(|(_, total)| *total)
            .expect("the contig was counted");
        if total == 1 {
            mask[index] = true;
        }
    }
    count_passing(mask)
}

/// One interval of the output, a `SimpleInterval`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    /// The contig.
    pub contig: String,
    /// One-based inclusive start.
    pub start: i32,
    /// One-based inclusive end.
    pub end: i32,
}

/// The Picard interval list: a SAM header, then five tab-separated columns per interval.
///
/// `PreprocessIntervals` writes the same list from the same struct, and so does this tool: the
/// dictionary is the input collection's METADATA, copied field for field. A port that rebuilt the
/// `@SQ` line from a name and a length alone dropped the `M5` the annotated intervals carried,
/// which is a differing byte in an otherwise identical file.
pub fn write(
    sequences: &[crate::preprocess_intervals::Sequence],
    intervals: &[Interval],
) -> String {
    crate::preprocess_intervals::write_list(sequences, intervals)
}
