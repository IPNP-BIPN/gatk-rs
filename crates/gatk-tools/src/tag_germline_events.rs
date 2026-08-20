//! `TagGermlineEvents`, ported from
//! `org.broadinstitute.hellbender.tools.copynumber.utils.TagGermlineEvents` and the
//! `SimpleGermlineTagger` under it (GATK 4.6.2.0).
//!
//! Each tumour segment is asked whether the matched normal carries the same event, and the answer
//! is written as one more column.
//!
//! # Two tests over two different lists, joined by OR
//!
//! A merged, non-neutral normal region tags when BOTH its breakpoints are seen within the padding
//! among the tumour segments, OR when one MERGED tumour run reciprocally overlaps it past the
//! threshold. The breakpoint search runs over the unmerged segments and the overlap over the merged
//! ones, so a normal whose breakpoints are out of reach still tags through the overlap.
//!
//! # And a third test decides which segments are tagged
//!
//! ```java
//! .filter(s -> ((Math.abs(s.getStart() - normalSeg.getStart()) <= paddingInBp) || (Math.abs(normalSeg.getEnd() - s.getEnd()) <= paddingInBp)
//!         || ((normalSeg.getStart() < s.getStart()) && (normalSeg.getEnd() > s.getEnd())))
//!         && (normalSeg.getInterval().intersect(s).size() > (s.getInterval().size() * reciprocalThreshold)))
//! ```
//!
//! One of the segment's own breakpoints within the padding, or STRICT containment by the normal
//! region, and then an intersection STRICTLY larger than the segment's own length times the
//! threshold. The strictness differs from `isReciprocalOverlap`, which is `>=` on both sides, and
//! at a threshold of zero the two therefore disagree about a segment the normal merely touches.
//!
//! # The merge that precedes all of it keeps only the call
//!
//! `mergedRegionsByAnnotation` builds each merged region with a map holding the call annotation and
//! nothing else, so whatever else those rows carried is gone before any comparison.

use crate::annotated_interval::{sort_by_dictionary, AnnotatedInterval};
use std::collections::BTreeMap;

/// `TagGermlineEvents.GERMLINE_TAG_HEADER`.
pub const GERMLINE_TAG_HEADER: &str = "POSSIBLE_GERMLINE";
/// `DEFAULT_GERMLINE_TAG_PADDING_IN_BP`, which is a thousand bases: the breakpoints are allowed to
/// be a long way apart before the padding stops seeing them.
pub const DEFAULT_PADDING_IN_BP: i32 = 1000;
/// The tool's `reciprocal-threshold` default.
pub const DEFAULT_RECIPROCAL_THRESHOLD: f64 = 0.75;
/// `CalledCopyRatioSegment.Call.NEUTRAL.getOutputString()`.
pub const NEUTRAL: &str = "0";
/// The three calls, as they are written.
pub const CALLS: [&str; 3] = ["+", "-", "0"];

/// What the tagger refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TagError {
    /// A tumour segment with no call.
    EmptyTumourCall(String),
    /// A normal segment with no call.
    EmptyNormalCall(String),
    /// `ParamUtils.isPositiveOrZero` on the padding.
    NegativePadding,
    /// `ParamUtils.inRange` on the threshold.
    ThresholdOutOfRange,
    /// `validateNoOverlappingIntervals`, which names the region and then every overlap of it.
    OverlappingIntervals { first: String, second: String },
}

impl TagError {
    pub fn java_class(&self) -> &'static str {
        match self {
            TagError::OverlappingIntervals { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
            _ => "java.lang.IllegalArgumentException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            TagError::EmptyTumourCall(column) => format!(
                "All tumor segments must have a call.  Call annotation (column header) must be: \
                 {column}"
            ),
            TagError::EmptyNormalCall(column) => format!(
                "All normal segments must have a call.  Call annotation (column header) name must \
                 be: {column}"
            ),
            TagError::NegativePadding => {
                "padding must be greater than or equal to zero.".to_string()
            }
            TagError::ThresholdOutOfRange => {
                "Reciprocal threshold must be between 0.0 and 1.0".to_string()
            }
            TagError::OverlappingIntervals { first, second } => {
                format!("Bad input: Overlap detected in input:  {first} overlapped {second}")
            }
        }
    }
}

/// `AnnotatedInterval.toString()`, which the overlap refusal names its rows with.
pub fn java_string(interval: &AnnotatedInterval) -> String {
    let annotations = interval
        .annotations
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect::<Vec<String>>()
        .join(", ");
    format!(
        "AnnotatedInterval{{interval={}:{}-{}, annotations={{{}}}}}",
        interval.contig, interval.start, interval.end, annotations
    )
}

fn overlaps(left: &AnnotatedInterval, right: &AnnotatedInterval) -> bool {
    left.contig == right.contig && left.start <= right.end && right.start <= left.end
}

/// The number of bases two regions share, which is zero when they do not overlap.
fn intersection(left: &AnnotatedInterval, right: &AnnotatedInterval) -> i64 {
    if !overlaps(left, right) {
        return 0;
    }
    i64::from(left.end.min(right.end) - left.start.max(right.start) + 1)
}

fn size(interval: &AnnotatedInterval) -> i64 {
    i64::from(interval.end - interval.start + 1)
}

/// `IntervalUtils.isReciprocalOverlap`, whose comparisons are both `>=` and which answers true for
/// any pair at all when the threshold is zero.
pub fn is_reciprocal_overlap(
    first: &AnnotatedInterval,
    second: &AnnotatedInterval,
    threshold: f64,
) -> bool {
    if threshold == 0.0 {
        return true;
    }
    overlaps(first, second)
        && intersection(first, second) as f64 >= size(second) as f64 * threshold
        && intersection(second, first) as f64 >= size(first) as f64 * threshold
}

/// `mergedRegionsByAnnotation`: neighbouring regions of the same contig and the same call, merged
/// into a region carrying ONLY that call.
pub fn merged_regions_by_annotation(
    regions: &[AnnotatedInterval],
    call_annotation: &str,
) -> Vec<AnnotatedInterval> {
    let mut merged = Vec::new();
    let mut index = 0;
    while index < regions.len() {
        let first = &regions[index];
        let call = first.annotations.get(call_annotation).cloned();
        let mut end = first.end;
        index += 1;
        while index < regions.len()
            && regions[index].contig == first.contig
            && regions[index].annotations.get(call_annotation) == call.as_ref()
        {
            end = regions[index].end;
            index += 1;
        }
        let mut annotations = BTreeMap::new();
        if let Some(call) = call {
            annotations.insert(call_annotation.to_string(), call);
        }
        merged.push(AnnotatedInterval {
            contig: first.contig.clone(),
            start: first.start,
            end,
            annotations,
        });
    }
    merged
}

/// `validateNoOverlappingIntervals`, which asks the overlap detector for each region's overlaps
/// and refuses when there is more than one.
///
/// The one is the region itself, so the message names it twice: `locatable` and then the whole set
/// of overlaps, itself included. That set is a `HashSet`, so its order is Java's hash order; with
/// the pairs this suite carries it is the input's order, which is what this port writes.
fn validate_no_overlaps(
    regions: &[AnnotatedInterval],
    _dictionary: &[String],
) -> Result<(), TagError> {
    for region in regions {
        let overlapping: Vec<&AnnotatedInterval> = regions
            .iter()
            .filter(|other| overlaps(region, other))
            .collect();
        if overlapping.len() > 1 {
            return Err(TagError::OverlappingIntervals {
                first: java_string(region),
                second: overlapping
                    .iter()
                    .map(|other| java_string(other))
                    .collect::<Vec<String>>()
                    .join(", "),
            });
        }
    }
    Ok(())
}

/// `tagTumorSegmentsWithGermlineActivity`.
pub fn tag_tumour_segments(
    tumour: &[AnnotatedInterval],
    normal: &[AnnotatedInterval],
    call_annotation: &str,
    dictionary: &[String],
    output_annotation: &str,
    padding_in_bp: i32,
    reciprocal_threshold: f64,
) -> Result<Vec<AnnotatedInterval>, TagError> {
    validate_no_overlaps(tumour, dictionary)?;
    validate_no_overlaps(normal, dictionary)?;
    if padding_in_bp < 0 {
        return Err(TagError::NegativePadding);
    }
    if !(0.0..=1.0).contains(&reciprocal_threshold) {
        return Err(TagError::ThresholdOutOfRange);
    }
    let mut tumour_segments = tumour.to_vec();
    let mut normal_segments = normal.to_vec();
    sort_by_dictionary(&mut tumour_segments, dictionary);
    sort_by_dictionary(&mut normal_segments, dictionary);

    let empty = |segment: &AnnotatedInterval| {
        segment
            .annotations
            .get(call_annotation)
            .map(String::is_empty)
            .unwrap_or(true)
    };
    // The normal is checked first, which is why a file with both faults names the normal.
    if normal_segments.iter().any(empty) {
        return Err(TagError::EmptyNormalCall(call_annotation.to_string()));
    }
    if tumour_segments.iter().any(empty) {
        return Err(TagError::EmptyTumourCall(call_annotation.to_string()));
    }

    let merged_normal = merged_regions_by_annotation(&normal_segments, call_annotation);
    let interesting: Vec<&AnnotatedInterval> = merged_normal
        .iter()
        .filter(|region| {
            region
                .annotations
                .get(call_annotation)
                .is_some_and(|call| !call.is_empty() && call != NEUTRAL)
        })
        .collect();

    // The tag of each tumour segment, by its position in the sorted list.
    let mut tags: Vec<Option<String>> = vec![None; tumour_segments.len()];
    for normal_region in interesting {
        let overlapping: Vec<usize> = tumour_segments
            .iter()
            .enumerate()
            .filter(|(_, segment)| overlaps(normal_region, segment))
            .map(|(index, _)| index)
            .collect();
        if overlapping.is_empty() {
            continue;
        }
        let overlapping_segments: Vec<AnnotatedInterval> = overlapping
            .iter()
            .map(|index| tumour_segments[*index].clone())
            .collect();
        let merged_tumour = merged_regions_by_annotation(&overlapping_segments, call_annotation);
        let reciprocal_seen = merged_tumour
            .iter()
            .any(|run| is_reciprocal_overlap(run, normal_region, reciprocal_threshold));
        let start_seen = overlapping_segments
            .iter()
            .any(|segment| (segment.start - normal_region.start).abs() <= padding_in_bp);
        let end_seen = overlapping_segments
            .iter()
            .any(|segment| (segment.end - normal_region.end).abs() <= padding_in_bp);
        if !((start_seen && end_seen) || reciprocal_seen) {
            continue;
        }
        let call = normal_region
            .annotations
            .get(call_annotation)
            .cloned()
            .unwrap_or_default();
        for index in overlapping {
            let segment = &tumour_segments[index];
            let breakpoint = (segment.start - normal_region.start).abs() <= padding_in_bp
                || (normal_region.end - segment.end).abs() <= padding_in_bp;
            let contained = normal_region.start < segment.start && normal_region.end > segment.end;
            // Strictly larger, unlike the reciprocal test's `>=`.
            let enough = intersection(normal_region, segment) as f64
                > size(segment) as f64 * reciprocal_threshold;
            if (breakpoint || contained) && enough {
                tags[index] = Some(call.clone());
            }
        }
    }

    Ok(tumour_segments
        .into_iter()
        .enumerate()
        .map(|(index, mut segment)| {
            let tag = tags[index].clone().unwrap_or_else(|| NEUTRAL.to_string());
            segment
                .annotations
                .insert(output_annotation.to_string(), tag);
            segment
        })
        .collect())
}
