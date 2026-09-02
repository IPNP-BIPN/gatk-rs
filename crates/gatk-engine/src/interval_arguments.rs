//! `IntervalArgumentCollection`: the five arguments that turn `-L` and `-XL` into the list a
//! walker traverses.
//!
//! [`crate::interval`] resolves ONE query string. This is the layer above it, and it is where four
//! arguments live that a port reading only `--intervals` silently ignores:
//! `--interval-padding`, `--interval-exclusion-padding`, `--interval-set-rule` and
//! `--exclude-intervals`.
//!
//! # Padding happens per argument, before the set operator
//!
//! Each `-L` string is resolved, padded, and then **sorted and merged with `ALL`** whatever the
//! merging rule is; only then is the result folded into the accumulator. So padding can make two
//! intervals of the same argument into one before an intersection ever sees them, and the merging
//! rule the user chose does not apply to that step.
//!
//! # The fold short-circuits on an empty side
//!
//! `mergeListsBySetOperator` returns the other list whenever either is empty, so the FIRST `-L` is
//! never intersected with anything: it becomes the accumulator. And `INTERSECTION` over three
//! arguments is `((a ∩ b) ∩ c)` rather than the intersection of all three at once, which is the
//! same answer for intervals and not the same for the empty case.
//!
//! # Two of the outcomes are refusals rather than empty traversals
//!
//! An empty intersection and an exclusion that removes every included base are both
//! `CommandLineException$BadArgumentValue`, and each quotes the raw argument strings back.
//!
//! # `-XL` alone means the whole reference
//!
//! With no `-L`, the include set is built from the sequence dictionary, contig by contig, in the
//! dictionary's own order. It is not an empty set, and it is not "everything" as a special case
//! the traversal understands later.
//!
//! # `unmapped` is a flag, not an interval
//!
//! It is accepted on `-L`, where it sets the traversal's unmapped bit and is removed from the
//! interval list, and refused on `-XL`.
//!
//! Ported from `org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection`,
//! `org.broadinstitute.hellbender.utils.IntervalUtils` (`loadIntervals`, `mergeListsBySetOperator`,
//! `getIntervalsWithFlanks`) and
//! `org.broadinstitute.hellbender.utils.GenomeLocSortedSet.subtractRegions`.

use htsjdk_bam::header::SamHeader;

use crate::interval::{
    load_intervals, merge_interval_locations, parse_interval, MergingRule, ParseError,
    SimpleInterval,
};

/// The string `-L` and `-XL` accept in place of an interval.
pub const UNMAPPED: &str = "unmapped";

/// `IntervalSetRule`, whose default is `UNION`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SetRule {
    #[default]
    Union,
    Intersection,
}

impl SetRule {
    /// The constant's own name, which a refusal quotes.
    pub fn name(self) -> &'static str {
        match self {
            SetRule::Union => "UNION",
            SetRule::Intersection => "INTERSECTION",
        }
    }
}

/// What the collection refuses, each of which quotes the arguments back.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntervalArgumentError {
    /// A query string that does not resolve, which the parser refuses first.
    Parse(ParseError),
    /// `INTERSECTION` that left nothing.
    EmptyIntersection { include: Vec<String>, rule: SetRule },
    /// `-XL` that removed every included base.
    ExcludedEverything {
        include: Vec<String>,
        exclude: Vec<String>,
    },
    /// `-XL unmapped`, which the reference has never supported.
    UnmappedExcluded,
}

/// `List.toString` on a list of strings, which is what the messages embed.
fn java_list(values: &[String]) -> String {
    format!("[{}]", values.join(", "))
}

impl IntervalArgumentError {
    pub fn java_class(&self) -> &'static str {
        match self {
            IntervalArgumentError::Parse(_) => {
                "org.broadinstitute.hellbender.exceptions.UserException$MalformedGenomeLoc"
            }
            IntervalArgumentError::EmptyIntersection { .. }
            | IntervalArgumentError::ExcludedEverything { .. } => {
                "org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue"
            }
            IntervalArgumentError::UnmappedExcluded => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
        }
    }

    pub fn message(&self) -> String {
        match self {
            IntervalArgumentError::Parse(error) => format!("{error:?}"),
            IntervalArgumentError::EmptyIntersection { include, rule } => format!(
                "Argument -L, --interval-set-rule has a bad value: {},{}. The specified intervals \
                 had an empty intersection",
                java_list(include),
                rule.name()
            ),
            IntervalArgumentError::ExcludedEverything { include, exclude } => format!(
                "Argument -L,-XL has a bad value: {}, {}. The intervals specified for exclusion \
                 with -XL removed all territory specified by -L.",
                java_list(include),
                java_list(exclude)
            ),
            IntervalArgumentError::UnmappedExcluded => {
                "-XL unmapped is not currently supported".to_string()
            }
        }
    }
}

/// What a walker traverses: the intervals, and whether unmapped records are included.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraversalParameters {
    pub intervals: Vec<SimpleInterval>,
    pub traverse_unmapped: bool,
}

fn contig_index(header: &SamHeader, contig: &str) -> usize {
    header
        .sequences
        .iter()
        .position(|sequence| sequence.name == contig)
        .unwrap_or(usize::MAX)
}

fn contig_length(header: &SamHeader, contig: &str) -> i32 {
    header
        .sequences
        .iter()
        .find(|sequence| sequence.name == contig)
        .map(|sequence| sequence.length)
        .unwrap_or(0)
}

/// `sortAndMergeIntervals`: the dictionary's order, then the merging rule.
fn sort_and_merge(
    mut intervals: Vec<SimpleInterval>,
    header: &SamHeader,
    rule: MergingRule,
) -> Vec<SimpleInterval> {
    intervals.sort_by_key(|interval| {
        (
            contig_index(header, &interval.contig),
            interval.start,
            interval.end,
        )
    });
    merge_interval_locations(intervals, rule)
}

/// `mergeListsBySetOperator`, whose empty short-circuit is what makes the first argument special.
///
/// `Err(())` is `UserException.EmptyIntersection`, which the caller turns into the refusal that
/// quotes the argument strings.
fn merge_by_set_operator(
    one: Vec<SimpleInterval>,
    two: Vec<SimpleInterval>,
    rule: SetRule,
    header: &SamHeader,
) -> Result<Vec<SimpleInterval>, ()> {
    if one.is_empty() || two.is_empty() {
        return Ok(if one.is_empty() { two } else { one });
    }
    if rule == SetRule::Union {
        let mut all = one;
        all.extend(two);
        return Ok(all);
    }

    // The intersection walks both lists once, taking the overlap and dropping whichever interval
    // ends first. `isBefore` orders by CONTIG INDEX first, which is the dictionary's order and not
    // the contig name's, so two lists on different contigs advance rather than intersect.
    let mut result = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < one.len() && j < two.len() {
        let (a, b) = (&one[i], &two[j]);
        let (ai, bi) = (
            contig_index(header, &a.contig),
            contig_index(header, &b.contig),
        );
        let b_before_a = bi < ai || (bi == ai && b.end < a.start);
        let a_before_b = ai < bi || (ai == bi && a.end < b.start);
        if b_before_a {
            j += 1;
        } else if a_before_b {
            i += 1;
        } else {
            result.push(SimpleInterval {
                contig: a.contig.clone(),
                start: a.start.max(b.start),
                end: a.end.min(b.end),
            });
            if a.end < b.end {
                i += 1;
            } else {
                j += 1;
            }
        }
    }
    if result.is_empty() {
        return Err(());
    }
    Ok(result)
}

/// `loadIntervals` with a set rule and padding: resolve, pad, fold.
fn load_with(
    queries: &[String],
    header: &SamHeader,
    set_rule: SetRule,
    merging_rule: MergingRule,
    padding: i32,
) -> Result<Vec<SimpleInterval>, IntervalArgumentError> {
    let mut accumulated: Vec<SimpleInterval> = Vec::new();
    for query in queries {
        let mut resolved =
            vec![parse_interval(query, header).map_err(IntervalArgumentError::Parse)?];
        if padding > 0 {
            // `getIntervalsWithFlanks` pads and then sorts and merges with ALL, whatever the
            // caller's merging rule is.
            resolved = resolved
                .iter()
                .filter_map(|interval| {
                    interval.expand_within_contig(padding, contig_length(header, &interval.contig))
                })
                .collect();
            resolved = sort_and_merge(resolved, header, MergingRule::All);
        }
        accumulated =
            merge_by_set_operator(resolved, accumulated, set_rule, header).map_err(|()| {
                IntervalArgumentError::EmptyIntersection {
                    include: queries.to_vec(),
                    rule: set_rule,
                }
            })?;
    }
    Ok(sort_and_merge(accumulated, header, merging_rule))
}

/// `GenomeLocSortedSet.subtractRegions`, over two sorted lists.
fn subtract(
    include: &[SimpleInterval],
    exclude: &[SimpleInterval],
    header: &SamHeader,
) -> Vec<SimpleInterval> {
    let mut kept = Vec::new();
    for interval in include {
        let mut pieces = vec![interval.clone()];
        for hole in exclude {
            if hole.contig != interval.contig {
                continue;
            }
            let mut next = Vec::new();
            for piece in pieces {
                if piece.end < hole.start || hole.end < piece.start {
                    next.push(piece);
                    continue;
                }
                if piece.start < hole.start {
                    if let Some(left) =
                        SimpleInterval::new(&piece.contig, piece.start, hole.start - 1)
                    {
                        next.push(left);
                    }
                }
                if hole.end < piece.end {
                    if let Some(right) = SimpleInterval::new(&piece.contig, hole.end + 1, piece.end)
                    {
                        next.push(right);
                    }
                }
            }
            pieces = next;
        }
        kept.extend(pieces);
    }
    sort_and_merge(kept, header, MergingRule::OverlappingOnly)
}

/// `createSetFromSequenceDictionary`: the whole reference, contig by contig, in its own order.
pub fn whole_reference(header: &SamHeader) -> Vec<SimpleInterval> {
    header
        .sequences
        .iter()
        .filter_map(|sequence| SimpleInterval::new(&sequence.name, 1, sequence.length))
        .collect()
}

/// `parseIntervals`: everything the five arguments decide, in the reference's own order.
pub fn traversal_parameters(
    include: &[String],
    exclude: &[String],
    header: &SamHeader,
    set_rule: SetRule,
    merging_rule: MergingRule,
    padding: i32,
    exclusion_padding: i32,
) -> Result<TraversalParameters, IntervalArgumentError> {
    // `unmapped` is a flag rather than an interval, and the parser never sees it.
    let mut traverse_unmapped = false;
    let included: Vec<String> = include
        .iter()
        .filter(|query| {
            if query.as_str() == UNMAPPED {
                traverse_unmapped = true;
                false
            } else {
                true
            }
        })
        .cloned()
        .collect();
    if exclude.iter().any(|query| query == UNMAPPED) {
        return Err(IntervalArgumentError::UnmappedExcluded);
    }

    let include_set = if included.is_empty() {
        if include.is_empty() {
            // No `-L` at all: the include set is the whole reference.
            whole_reference(header)
        } else {
            // `-L unmapped` alone: every mapped interval is gone, and the flag is the answer.
            Vec::new()
        }
    } else {
        load_with(&included, header, set_rule, merging_rule, padding)?
    };

    let exclude_set = load_with(
        exclude,
        header,
        SetRule::Union,
        merging_rule,
        exclusion_padding,
    )?;
    if exclude_set.is_empty() {
        return Ok(TraversalParameters {
            intervals: include_set,
            traverse_unmapped,
        });
    }

    let intervals = subtract(&include_set, &exclude_set, header);
    if intervals.is_empty() {
        return Err(IntervalArgumentError::ExcludedEverything {
            include: include.to_vec(),
            exclude: exclude.to_vec(),
        });
    }
    Ok(TraversalParameters {
        intervals,
        traverse_unmapped,
    })
}

/// The one-argument form the runners used before the other four existed.
pub fn intervals_only(
    include: &[String],
    header: &SamHeader,
    merging_rule: MergingRule,
) -> Result<Vec<SimpleInterval>, ParseError> {
    load_intervals(include, header, merging_rule)
}
