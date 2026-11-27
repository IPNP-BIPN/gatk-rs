//! Ported from `org.broadinstitute.hellbender.tools.copynumber.AnnotateIntervals`
//! (GATK 4.7.0.0; `AnnotateIntervals.java` is unchanged since 4.6.2.0).
//!
//! GC content per interval, and two optional annotations read from BED tracks. Three small
//! functions, each with a rule that is not the obvious one.
//!
//! # The GC denominator is not the interval's length
//!
//! `gcCount / (double) (gcCount + atCount)`, over the bases the CACHING READER returned -- which
//! has already upper-cased them and turned every IUPAC code into `N`. So an interval holding
//! ambiguity codes has a smaller denominator than its length, and one holding none of `ACGT` at all
//! answers `NaN` rather than zero.
//!
//! # The BED track is zero-based and the interval is one-based
//!
//! `CoordMath.getOverlap(feature.getStart(), feature.getEnd() - 1, interval.getStart(),
//! interval.getEnd())`: the feature's end is decremented and its start is not, and the interval's
//! coordinates go in untouched. A BED line covering `0-12` over the interval `1-12` therefore
//! overlaps by ELEVEN, and the annotation is eleven twelfths rather than one.
//!
//! ```java
//! lengthWeightedSum += score *
//!         CoordMath.getOverlap(
//!                 feature.getStart(), feature.getEnd() - 1,       // zero-based
//!                 interval.getStart(), interval.getEnd());        // one-based
//! ```
//!
//! (`AnnotateIntervals.java:326-329`; the trailing comments are the reference's.) The "zero-based"
//! is the file's convention, not `getStart()`'s: htsjdk 5.0.0's `BEDCodec`, built with its default
//! `StartOffset.ONE`, adds one to the start column as it parses it:
//!
//! ```java
//! int start = Integer.parseInt(tokens[1]) + startOffsetValue;
//! ```
//!
//! So the reference's feature is one-based at both ends, and the `- 1` drops its last base. In the
//! golden every feature starts where an interval starts, so the interval's start is the one
//! `getOverlap` keeps and the start convention never shows.
//!
//! # A missing score is one
//!
//! `Double.isNaN(score) ? 1. : score`, so a BED file with no score column annotates every overlap
//! at full weight rather than at none.
//!
//! # Reading this file without reading Rust
//!
//! - `&[u8]` is a borrowed, read-only view of a byte array, and `&[&BedFeature]` one of
//!   references to features the caller owns; `b'C'` is the byte value of the character `C`.
//! - `i32`, `i64` and `f64` are Java's `int`, `long` and `double`; `f64::from(x)` and `x as f64`
//!   convert an integer to a double, exactly for every value here.
//! - `a.max(b)` and `a.min(b)` are `Math.max(a, b)` and `Math.min(a, b)`.
//! - `match` must cover every case; `_ => {}` is the "anything else, do nothing" arm, and
//!   `Some(..) if cond =>` is an arm that also has to pass `cond`.
//! - A body's last expression, with no `return` and no `;`, is its value.

/// One BED feature, reduced to what the annotators read.
#[derive(Debug, Clone, PartialEq)]
pub struct BedFeature {
    /// The contig.
    pub contig: String,
    /// Zero-based, inclusive.
    pub start: i32,
    /// Zero-based, exclusive -- which is why the annotator decrements it.
    pub end: i32,
    /// `Float.NaN` when the column is absent, which the annotator reads as one.
    pub score: f64,
}

/// `CopyNumberAnnotations`' keys, in the order the annotators are added.
pub const GC_CONTENT: &str = "GC_CONTENT";
/// The mappability annotation's column name.
pub const MAPPABILITY: &str = "MAPPABILITY";
/// The segmental-duplication annotation's column name.
pub const SEGMENTAL_DUPLICATION_CONTENT: &str = "SEGMENTAL_DUPLICATION_CONTENT";

/// `GCContentAnnotator.apply`.
///
/// ```java
/// final long gcCount = counter.get(Nucleotide.C) + counter.get(Nucleotide.G);
/// final long atCount = counter.get(Nucleotide.A) + counter.get(Nucleotide.T);
/// final long totalCount = gcCount + atCount;
/// return totalCount == 0 ? Double.NaN : gcCount / (double) totalCount;
/// ```
///
/// (`AnnotateIntervals.java:295-298`.) `bases` is what `ReferenceContext.getBases()` returned, so
/// it is already upper-cased with the IUPAC codes flattened. Anything that is not `ACGT` -- an
/// `N`, most of all -- counts in neither the numerator nor the denominator.
pub fn gc_content(bases: &[u8]) -> f64 {
    let mut gc = 0i64;
    let mut at = 0i64;
    for base in bases {
        match base.to_ascii_uppercase() {
            b'C' | b'G' => gc += 1,
            b'A' | b'T' => at += 1,
            _ => {}
        }
    }
    let total = gc + at;
    if total == 0 {
        f64::NAN
    } else {
        gc as f64 / total as f64
    }
}

/// `CoordMath.getOverlap(start1, end1, start2, end2)`.
///
/// `Math.max(0, Math.min(end1, end2) - Math.max(start1, start2) + 1)`, all inclusive.
fn overlap(first_start: i32, first_end: i32, second_start: i32, second_end: i32) -> i32 {
    0.max(first_end.min(second_end) - first_start.max(second_start) + 1)
}

/// `BEDLengthWeightedAnnotator.apply`: the score-weighted overlap over the interval's length.
///
/// `features[i].start` is compared as given; the reference compares the codec's start, which is
/// the BED column plus one (see the module note), so a caller holding the raw column passes a
/// start one less than the reference's.
///
/// The features are the ones the query returned for this interval; a feature that does not overlap
/// contributes zero rather than being an error. The sum is a plain `+=` loop rather than a stream,
/// which is one of the few places in this port where it is.
pub fn length_weighted_annotation(
    features: &[&BedFeature],
    interval_start: i32,
    interval_end: i32,
) -> f64 {
    let mut sum = 0.0;
    for feature in features {
        // A missing score arrives as NaN and is read as one.
        let score = if feature.score.is_nan() {
            1.0
        } else {
            feature.score
        };
        // The feature's end is decremented and the interval's is not: the BED half-open end becomes
        // an inclusive one, and the interval was inclusive already.
        sum += score
            * f64::from(overlap(
                feature.start,
                feature.end - 1,
                interval_start,
                interval_end,
            ));
    }
    sum / f64::from(interval_end - interval_start + 1)
}

/// What the tool refuses before it reads anything.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnnotateError {
    /// The engine's default merging rule, which this tool will not run under.
    MergingRuleNotOverlappingOnly,
    /// A track whose own features overlap.
    OverlappingTrack(String),
    /// An interval of zero length, which the traversal refuses per interval.
    ZeroLengthInterval(String),
}

impl AnnotateError {
    /// The exception class the reference throws.
    pub fn java_class(&self) -> &'static str {
        match self {
            AnnotateError::MergingRuleNotOverlappingOnly => "java.lang.IllegalArgumentException",
            AnnotateError::OverlappingTrack(_) | AnnotateError::ZeroLengthInterval(_) => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
        }
    }

    /// The message, with the `Bad input: ` prefix where the exception adds one. The `BadInput`
    /// texts are `AnnotateIntervals.java:217-218` and `:226`.
    pub fn message(&self) -> String {
        match self {
            AnnotateError::MergingRuleNotOverlappingOnly => {
                "Interval merging rule must be set to OVERLAPPING_ONLY.".to_string()
            }
            AnnotateError::OverlappingTrack(path) => format!(
                "Bad input: Feature track {path} contains overlapping intervals; these should be \
                 merged."
            ),
            AnnotateError::ZeroLengthInterval(interval) => {
                format!("Bad input: Interval cannot have zero length: {interval}")
            }
        }
    }
}

/// `checkForOverlaps`: the track's own features merged with `OVERLAPPING_ONLY` and counted.
///
/// Two features that TOUCH survive the merge as two, and two that overlap become one. The check is
/// a count comparison rather than a scan, and it ITERATES the track, which is why it fires on a
/// file the query path could not have opened.
///
/// ```java
/// final List<GenomeLoc> mergedGenomeLocs = IntervalUtils.mergeIntervalLocations(
///         genomeLocs, IntervalMergingRule.OVERLAPPING_ONLY);
/// if (genomeLocs.size() != mergedGenomeLocs.size()) {
/// ```
///
/// (`AnnotateIntervals.java:214-216`.) `mergeIntervalLocations` is one pass in file order, merging
/// each location into the previous one when they overlap; it does not sort. Here the same pass
/// keeps the running merged span in `previous` and counts the spans it starts. The test
/// `feature.start < end` compares a zero-based start against an exclusive end, which is the
/// reference's one-based overlap test shifted by one on both sides, but only in one direction:
/// the reference's `overlapsP` also asks whether the feature ends before the previous span.
pub fn track_has_overlaps(features: &[BedFeature]) -> bool {
    let mut merged = 0;
    // `Option<(&str, i32, i32)>`: nothing yet, or the current merged span as (contig, start, end),
    // its contig borrowed from the feature rather than copied.
    let mut previous: Option<(&str, i32, i32)> = None;
    for feature in features {
        match previous {
            Some((contig, _, end)) if contig == feature.contig && feature.start < end => {
                // Overlapping: the merge absorbs it, so the count does not grow.
                previous = Some((&feature.contig, feature.start, end.max(feature.end)));
            }
            _ => {
                merged += 1;
                previous = Some((&feature.contig, feature.start, feature.end));
            }
        }
    }
    merged != features.len()
}

/// The header of the output table, which is the columns the annotators that ran produced.
pub fn columns(annotations: &[&str]) -> String {
    let mut names = vec!["CONTIG".to_string(), "START".to_string(), "END".to_string()];
    names.extend(annotations.iter().map(|name| name.to_string()));
    names.join("\t")
}

/// One row: the interval, then each annotation formatted to six decimals.
///
/// A `NaN` is written as the three characters `NaN`, which is what `String.format("%.6f", ...)`
/// does with one.
pub fn row(contig: &str, start: i32, end: i32, annotations: &[f64]) -> String {
    let mut fields = vec![contig.to_string(), start.to_string(), end.to_string()];
    for value in annotations {
        fields.push(if value.is_nan() {
            "NaN".to_string()
        } else {
            gatk_engine::java_format::format_decimals(*value, 6)
        });
    }
    fields.join("\t")
}
