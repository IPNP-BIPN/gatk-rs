//! Ported from `org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationSegmenter`.
//!
//! Seventy-four lines of Java that turn a list of pileup summaries into segments of constant minor
//! allele fraction, by way of [`crate::kernel_segmenter`] and the decomposition under it.
//!
//! # The changepoint convention is off by one on purpose
//!
//! The reference says it in a comment: a changepoint at index `n` means index `n` belongs to the
//! LEFT segment. That is why the list is bracketed with `-1` rather than `0`, and why every start
//! is `changepoint + 1`. A port that used the usual end-exclusive convention would shift every
//! boundary by one site.
//!
//! # What the segmenter sees is not what it returns
//!
//! Segmentation runs over the HET sites alone, those whose alt fraction is inside `[0.1, 0.9]`.
//! The intervals it produces are then used to look up **all** sites that overlap them, hom refs and
//! hom alts included, which is how a hom alt ends up inside a segment whose boundaries no hom alt
//! helped choose.

use crate::java_hash;
use crate::kernel_segmenter::{
    find_changepoints, segmentation_kernel, ChangepointSortOrder, SEGMENTATION_KERNEL_VARIANCE,
};
use crate::pileup_summary::PileupSummary;

/// `ALT_FRACTIONS_FOR_SEGMENTATION`, a `Range.of(0.1, 0.9)` and therefore closed at both ends.
pub const ALT_FRACTION_LOWER: f64 = 0.1;
/// The upper end of the same range.
pub const ALT_FRACTION_UPPER: f64 = 0.9;
/// `KERNEL_SEGMENTER_LINEAR_COST`.
pub const KERNEL_SEGMENTER_LINEAR_COST: f64 = 1.0;
/// `KERNEL_SEGMENTER_LOG_LINEAR_COST`.
pub const KERNEL_SEGMENTER_LOG_LINEAR_COST: f64 = 1.0;
/// `KERNEL_SEGMENTER_DIMENSION`.
pub const KERNEL_SEGMENTER_DIMENSION: usize = 100;
/// `POINTS_PER_SEGMENTATION_WINDOW`.
pub const POINTS_PER_SEGMENTATION_WINDOW: usize = 50;
/// `MAX_CHANGEPOINTS_PER_CHROMOSOME`.
pub const MAX_CHANGEPOINTS_PER_CHROMOSOME: usize = 10;

/// One segment's span, which is `SimpleInterval` reduced to what the caller reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Span {
    /// The contig the segment is on.
    pub contig: String,
    /// The first site's start.
    pub start: i32,
    /// The last site's end, which for a pileup summary is its position.
    pub end: i32,
}

/// `findSegments`: every contig's segments, each carrying the sites that overlap it.
///
/// The contigs are visited in `HashMap` order, not in the order they appear in the table:
/// `Collectors.groupingBy` builds a `HashMap` and the reference streams its values. That order is
/// what [`crate::java_hash::hash_set_order`] reproduces, and getting it wrong reorders the whole
/// answer even though each segment is unchanged.
pub fn find_segments(sites: &[PileupSummary]) -> Vec<Vec<PileupSummary>> {
    let mut contigs: Vec<String> = Vec::new();
    for site in sites {
        if !contigs.contains(&site.contig) {
            contigs.push(site.contig.clone());
        }
    }
    let ordered = java_hash::hash_set_order(&contigs).expect("the contig names hash into buckets");

    let mut segments = Vec::new();
    for contig in &ordered {
        let contig_sites: Vec<PileupSummary> = sites
            .iter()
            .filter(|site| &site.contig == contig)
            .cloned()
            .collect();
        for span in find_contig_segments(&contig_sites) {
            // `od.getOverlaps(segment)`, sorted by start. The detector holds every site rather than
            // only this contig's, so the contig has to be compared as well as the position.
            let mut overlapping: Vec<PileupSummary> = sites
                .iter()
                .filter(|site| {
                    site.contig == span.contig
                        && site.position >= span.start
                        && site.position <= span.end
                })
                .cloned()
                .collect();
            overlapping.sort_by_key(|site| site.position);
            segments.push(overlapping);
        }
    }
    segments
}

/// `findContigSegments`, on one contig's sites.
///
/// An empty het list is an empty answer rather than one segment: a contig with nothing between the
/// two fractions contributes no segment at all, and its sites therefore appear in no segment.
fn find_contig_segments(sites: &[PileupSummary]) -> Vec<Span> {
    let het_sites: Vec<&PileupSummary> = sites
        .iter()
        .filter(|site| {
            let fraction = site.alt_fraction();
            // `Range.contains`, which is closed at both ends. Written as two comparisons rather
            // than as a `RangeInclusive`, whose `contains` would answer the same here but reads as
            // a Rust range rather than as the reference's `Range.of(0.1, 0.9)`.
            #[allow(clippy::manual_range_contains)]
            {
                fraction >= ALT_FRACTION_LOWER && fraction <= ALT_FRACTION_UPPER
            }
        })
        .collect();
    if het_sites.is_empty() {
        return Vec::new();
    }

    let fractions: Vec<f64> = het_sites.iter().map(|site| site.alt_fraction()).collect();
    let found = find_changepoints(
        &fractions,
        MAX_CHANGEPOINTS_PER_CHROMOSOME,
        segmentation_kernel,
        KERNEL_SEGMENTER_DIMENSION,
        &[POINTS_PER_SEGMENTATION_WINDOW],
        KERNEL_SEGMENTER_LINEAR_COST,
        KERNEL_SEGMENTER_LOG_LINEAR_COST,
        ChangepointSortOrder::Index,
    );

    // The bracketing the reference explains in its comment: `-1` at the front so the first segment
    // starts at index zero, and the last index at the back so the last segment ends there.
    let mut changepoints: Vec<i64> = vec![-1];
    changepoints.extend(found.iter().map(|index| *index as i64));
    changepoints.push(het_sites.len() as i64 - 1);

    (0..changepoints.len() - 1)
        .map(|n| {
            let first = &het_sites[(changepoints[n] + 1) as usize];
            let last = &het_sites[changepoints[n + 1] as usize];
            Span {
                contig: first.contig.clone(),
                start: first.position,
                end: last.position,
            }
        })
        .collect()
}

/// The kernel variance the segmenter uses, re-exported so a caller does not reach into the
/// segmenter module for it.
pub const KERNEL_VARIANCE: f64 = SEGMENTATION_KERNEL_VARIANCE;
