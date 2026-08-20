//! `CombineSegmentBreakpoints`, ported from
//! `org.broadinstitute.hellbender.tools.copynumber.utils.CombineSegmentBreakpoints` and
//! `IntervalUtils.combineAndSortBreakpoints` (GATK 4.6.2.0).
//!
//! Two segment files cut at every breakpoint either of them carries, each piece annotated from
//! both.
//!
//! # The cutting walks a sorted list of breakpoints
//!
//! Every segment contributes a START breakpoint at its start and an END breakpoint at its end. The
//! breakpoints of one contig are collected into a SET, sorted by position with STARTS BEFORE ENDS
//! at equal positions, and each consecutive pair becomes a piece:
//!
//!  * start then start: `[current, next - 1]`;
//!  * end then end: `[current + 1, next]`;
//!  * start then end: `[current, next]`;
//!  * end then start: `[current + 1, next - 1]`, and the piece is kept only while more starts than
//!    ends have been seen, which is what drops the gap between two segments neither file covers.
//!
//! The set is what makes a shared breakpoint appear once, and the start-before-end order is what
//! keeps a single-base segment from collapsing.
//!
//! # A piece takes its annotations from the first segment it overlaps
//!
//! One map per input file, from the file's column name to the output's. A column both files carry
//! is suffixed with that file's label, which defaults to `1` and `2`; a column only one file
//! carries keeps its name. A piece that no segment of a file overlaps takes EMPTY STRINGS for that
//! file's columns rather than leaving them out, so every row has every column.

use crate::annotated_interval::AnnotatedInterval;
use std::collections::BTreeMap;

/// The tool's default labels.
pub const DEFAULT_LABELS: [&str; 2] = ["1", "2"];

/// A breakpoint: a position and which end of a segment it is.
///
/// The order of the enum is the order the sort uses at equal positions: starts first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Breakpoint {
    Start,
    End,
}

/// `combineAndSortBreakpoints`, for one contig's worth of segments at a time.
fn pieces_for_contig(contig: &str, breakpoints: &[(i32, Breakpoint)]) -> Vec<(i32, i32)> {
    // A `HashSet` in the reference, so a position that is a start in both files appears once.
    let mut unique: Vec<(i32, Breakpoint)> = Vec::new();
    for point in breakpoints {
        if !unique.contains(point) {
            unique.push(*point);
        }
    }
    unique.sort();
    let _ = contig;
    let mut pieces = Vec::new();
    let mut starts_seen = 0;
    let mut ends_seen = 0;
    for index in 0..unique.len().saturating_sub(1) {
        let (current, current_kind) = unique[index];
        let (next, next_kind) = unique[index + 1];
        let is_current_start = current_kind == Breakpoint::Start;
        let is_next_start = next_kind == Breakpoint::Start;
        let mut start = if !is_current_start && !is_next_start {
            current + 1
        } else {
            current
        };
        let mut end = if is_current_start && is_next_start {
            next - 1
        } else {
            next
        };
        let between = !is_current_start && is_next_start;
        if between {
            start += 1;
            end -= 1;
        }
        if is_current_start {
            starts_seen += 1;
        } else {
            ends_seen += 1;
        }
        if (!between || starts_seen > ends_seen) && start <= end {
            pieces.push((start, end));
        }
    }
    pieces
}

/// Every piece the two files cut each other into, in dictionary order.
pub fn combine_and_sort_breakpoints(
    first: &[AnnotatedInterval],
    second: &[AnnotatedInterval],
    dictionary: &[String],
) -> Vec<(String, i32, i32)> {
    let mut out = Vec::new();
    for contig in dictionary {
        let breakpoints: Vec<(i32, Breakpoint)> = first
            .iter()
            .chain(second.iter())
            .filter(|segment| &segment.contig == contig)
            .flat_map(|segment| {
                [
                    (segment.start, Breakpoint::Start),
                    (segment.end, Breakpoint::End),
                ]
            })
            .collect();
        if breakpoints.is_empty() {
            continue;
        }
        for (start, end) in pieces_for_contig(contig, &breakpoints) {
            out.push((contig.clone(), start, end));
        }
    }
    out
}

/// The map from one file's column names to the output's.
///
/// A column both files carry is suffixed with that file's label; a column only one file carries
/// keeps its name.
pub fn output_header_map(
    own: &[String],
    other: &[String],
    label: &str,
) -> BTreeMap<String, String> {
    own.iter()
        .map(|name| {
            let output = if other.contains(name) {
                format!("{name}_{label}")
            } else {
                name.clone()
            };
            (name.clone(), output)
        })
        .collect()
}

/// `annotateCombinedIntervals`: one row per piece, annotated from the FIRST segment of each file
/// that overlaps it, and with empty strings where a file has nothing.
pub fn combine(
    first: &[AnnotatedInterval],
    second: &[AnnotatedInterval],
    dictionary: &[String],
    labels: [&str; 2],
    columns_of_interest: &[String],
) -> Vec<AnnotatedInterval> {
    let keep = |segment: &AnnotatedInterval| -> Vec<String> {
        segment
            .annotations
            .keys()
            .filter(|name| columns_of_interest.is_empty() || columns_of_interest.contains(name))
            .cloned()
            .collect()
    };
    let first_columns = first.first().map(keep).unwrap_or_default();
    let second_columns = second.first().map(keep).unwrap_or_default();
    let first_map = output_header_map(&first_columns, &second_columns, labels[0]);
    let second_map = output_header_map(&second_columns, &first_columns, labels[1]);

    combine_and_sort_breakpoints(first, second, dictionary)
        .into_iter()
        .map(|(contig, start, end)| {
            let mut annotations: BTreeMap<String, String> = BTreeMap::new();
            for (segments, map) in [(first, &first_map), (second, &second_map)] {
                let matching = segments.iter().find(|segment| {
                    segment.contig == contig && segment.start <= end && start <= segment.end
                });
                for (input, output) in map {
                    let value = matching
                        .and_then(|segment| segment.annotations.get(input))
                        .cloned()
                        .unwrap_or_default();
                    annotations.insert(output.clone(), value);
                }
            }
            AnnotatedInterval {
                contig,
                start,
                end,
                annotations,
            }
        })
        .collect()
}
