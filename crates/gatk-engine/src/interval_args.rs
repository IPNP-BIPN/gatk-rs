//! Ported from `org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection`,
//! `org.broadinstitute.hellbender.utils.IntervalUtils`, `.GenomeLocParser` and
//! `.GenomeLocSortedSet` (GATK 4.6.2.0).
//!
//! `IntervalWalker.traverse` is five lines: it iterates `userIntervals` and calls `apply` on each.
//! Every decision the traversal makes was already taken before it starts, when `-L`, `-XL`,
//! `--interval-padding`, `--interval-set-rule` and `--interval-merging-rule` were turned into that
//! list. So this module is the walker, and the loop is a formality.
//!
//! Five behaviours here change which intervals a tool processes, and none is derivable from the
//! argument names:
//!
//!  * **padding is clamped, not extended.** `createGenomeLocOnContig` bounds the padded start at 1
//!    and the padded stop at the contig length, and returns *null* when both ends fall off the
//!    contig, so a padded interval can silently disappear rather than throw;
//!  * **padding is applied per `-L` argument, before the set rule**, and each padded batch is
//!    itself sorted and merged with `ALL` inside `getIntervalsWithFlanks`. Padding two arguments
//!    that then abut therefore merges them even when the merging rule is `OVERLAPPING_ONLY`,
//!    because the `ALL` merge already happened;
//!  * **`INTERSECTION` is a running fold, not an n-way intersection.** Each new argument is
//!    intersected against the accumulated result, and the accumulator starts empty, which the
//!    reference short-circuits: intersecting anything with an empty set returns the other set. So
//!    a single `-L` under `INTERSECTION` is the identity, not an empty result;
//!  * **an empty intersection is an error**, not an empty traversal;
//!  * **`-XL` subtracts territory, and its own padding is a separate argument.** The subtraction is
//!    a stack walk that can split one interval into two, and it always unions its own inputs
//!    whatever the include rule says.

use crate::interval::{self, MergingRule, ParseError, SimpleInterval};
use htsjdk_bam::header::SamHeader;

/// `IntervalSetRule`. Governs `-L` only: `-XL` is always unioned.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SetRule {
    #[default]
    Union,
    Intersection,
}

/// `TraversalParameters`: the intervals to walk, and whether unmapped records are wanted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraversalParameters {
    pub intervals: Vec<SimpleInterval>,
    /// `-L unmapped` was given. The string is separated out of the interval list rather than
    /// parsed into one, so it never reaches the sort or the merge.
    pub traverse_unmapped: bool,
}

/// What the interval arguments were, as parsed off the command line.
#[derive(Debug, Clone, Default)]
pub struct IntervalArguments {
    /// `-L`. Empty means "the whole reference" only when `-XL` is given; otherwise the tool has no
    /// user intervals at all.
    pub include: Vec<String>,
    /// `-XL`.
    pub exclude: Vec<String>,
    /// `--interval-padding`, `-ip`.
    pub padding: i32,
    /// `--interval-exclusion-padding`, `-ixp`.
    pub exclusion_padding: i32,
    /// `--interval-set-rule`, `-isr`.
    pub set_rule: SetRule,
    /// `--interval-merging-rule`, `-imr`.
    pub merging_rule: MergingRule,
}

impl IntervalArguments {
    /// `intervalsSpecified()`: either list being non-empty is enough.
    pub fn specified(&self) -> bool {
        !self.include.is_empty() || !self.exclude.is_empty()
    }
}

/// What can go wrong turning the arguments into a traversal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntervalArgumentError {
    /// A `-L` or `-XL` string did not parse.
    Parse(ParseError),
    /// `UserException.EmptyIntersection`, surfaced by Barclay as a bad argument value for
    /// `-L, --interval-set-rule`.
    EmptyIntersection,
    /// `-XL` removed every base `-L` asked for.
    ExclusionRemovedEverything,
    /// `-XL unmapped`, which the reference refuses outright.
    UnmappedExcluded,
    /// `getTraversalParameters` without any interval argument at all.
    NoIntervalsSpecified,
    /// The argument ends in `.list` or `.intervals` and there is no such file.
    IntervalFileMissing(String),
    /// `UserException.MalformedFile`: the interval file held no intervals.
    IntervalFileEmpty,
    /// The file exists and is neither a Feature file nor an interval file by extension.
    FileIsNeitherFeaturesNorIntervals(String),
    /// The removed `-L "a;b"` syntax, which is refused rather than split.
    LegacySemicolonSyntax(String),
}

impl From<ParseError> for IntervalArgumentError {
    fn from(error: ParseError) -> Self {
        IntervalArgumentError::Parse(error)
    }
}

/// The string `GenomeLocParser.isUnmappedGenomeLocString` accepts, case-insensitively after a trim.
fn is_unmapped_string(text: &str) -> bool {
    text.trim().eq_ignore_ascii_case("unmapped")
}

fn contig_index(header: &SamHeader, contig: &str) -> Option<usize> {
    header.sequences.iter().position(|s| s.name == contig)
}

fn contig_length(header: &SamHeader, contig: &str) -> Option<i32> {
    header
        .sequences
        .iter()
        .find(|s| s.name == contig)
        .map(|s| s.length)
}

/// `GenomeLocParser.createPaddedGenomeLoc` composed with `createGenomeLocOnContig`.
///
/// Returns `None` where the reference returns null: when the whole padded interval lies off the
/// contig. Clamping rather than extending is the point, and so is the fact that a padding of zero
/// short-circuits before any clamping, so an interval that is already off-contig survives `-ip 0`
/// and dies under `-ip 1`.
pub fn pad(interval: &SimpleInterval, padding: i32, header: &SamHeader) -> Option<SimpleInterval> {
    if padding == 0 {
        return Some(interval.clone());
    }
    let length = contig_length(header, &interval.contig)?;
    let bounded_start = (interval.start - padding).max(1);
    let bounded_stop = (interval.end + padding).min(length);
    if bounded_start > length || bounded_stop < 1 {
        return None;
    }
    Some(SimpleInterval {
        contig: interval.contig.clone(),
        start: bounded_start,
        end: bounded_stop,
    })
}

/// `IntervalUtils.sortAndMergeIntervals`: sort by contig index then coordinates, then merge.
fn sort_and_merge(
    mut intervals: Vec<SimpleInterval>,
    header: &SamHeader,
    rule: MergingRule,
) -> Vec<SimpleInterval> {
    intervals.sort_by_key(|i| {
        (
            contig_index(header, &i.contig).unwrap_or(usize::MAX),
            i.start,
            i.end,
        )
    });
    interval::merge_interval_locations(intervals, rule)
}

/// `IntervalUtils.getIntervalsWithFlanks`: pad every interval, then sort and merge with `ALL`.
///
/// The `ALL` is hard-coded upstream and is not the user's merging rule, which is why padding can
/// join two intervals that `OVERLAPPING_ONLY` would have kept apart.
fn with_flanks(
    intervals: Vec<SimpleInterval>,
    padding: i32,
    header: &SamHeader,
) -> Vec<SimpleInterval> {
    if intervals.is_empty() {
        return intervals;
    }
    let padded: Vec<SimpleInterval> = intervals
        .iter()
        .filter_map(|i| pad(i, padding, header))
        .collect();
    sort_and_merge(padded, header, MergingRule::All)
}

/// `GenomeLoc.isBefore`: strictly before, by contig index then by coordinate.
fn is_before(a: &SimpleInterval, b: &SimpleInterval, header: &SamHeader) -> bool {
    let (ia, ib) = (
        contig_index(header, &a.contig).unwrap_or(usize::MAX),
        contig_index(header, &b.contig).unwrap_or(usize::MAX),
    );
    ia < ib || (ia == ib && a.end < b.start)
}

/// `IntervalUtils.mergeListsBySetOperator`.
///
/// The empty short-circuit at the top is what makes a fold over one argument the identity under
/// `INTERSECTION`, and it returns the *other* list unchanged rather than an intersection of it
/// with nothing.
pub fn merge_lists_by_set_operator(
    one: Vec<SimpleInterval>,
    two: Vec<SimpleInterval>,
    rule: SetRule,
    header: &SamHeader,
) -> Result<Vec<SimpleInterval>, IntervalArgumentError> {
    if one.is_empty() || two.is_empty() {
        return Ok(if one.is_empty() { two } else { one });
    }
    if rule == SetRule::Union {
        let mut all = one;
        all.extend(two);
        return Ok(all);
    }

    let mut result = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while j < two.len() && i < one.len() {
        if is_before(&two[j], &one[i], header) {
            j += 1;
        } else if is_before(&one[i], &two[j], header) {
            i += 1;
        } else {
            // They overlap, so intersect and drop whichever ends first.
            result.push(SimpleInterval {
                contig: one[i].contig.clone(),
                start: one[i].start.max(two[j].start),
                end: one[i].end.min(two[j].end),
            });
            if one[i].end < two[j].end {
                i += 1;
            } else {
                j += 1;
            }
        }
    }
    if result.is_empty() {
        return Err(IntervalArgumentError::EmptyIntersection);
    }
    Ok(result)
}

/// `IntervalUtils.GATK_INTERVAL_FILE_EXTENSIONS`.
pub const GATK_INTERVAL_FILE_EXTENSIONS: [&str; 2] = [".list", ".intervals"];

/// `IntervalUtils.isGatkIntervalFile`: does this argument *look* like an interval file?
///
/// Extension only, lower-cased, and deliberately not "does it contain intervals": a contig may
/// contain a period, so the reference refuses to treat the mere presence of an extension as
/// evidence, and an argument with one of these two extensions is a file even when it is missing.
pub fn has_gatk_interval_file_extension(query: &str) -> bool {
    let lowered = query.to_lowercase();
    GATK_INTERVAL_FILE_EXTENSIONS
        .iter()
        .any(|extension| lowered.ends_with(extension))
}

/// `IntervalUtils.gatkIntervalFileToList`: one interval per non-blank line, trimmed.
///
/// A file with no intervals in it is a `MalformedFile`, not an empty list, which is a different
/// outcome from `-L` with no argument at all.
pub fn gatk_interval_file_to_list(
    text: &str,
    header: &SamHeader,
) -> Result<Vec<SimpleInterval>, IntervalArgumentError> {
    let mut intervals = Vec::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        intervals.push(interval::parse_interval(trimmed, header)?);
    }
    if intervals.is_empty() {
        return Err(IntervalArgumentError::IntervalFileEmpty);
    }
    Ok(intervals)
}

/// What a `-L` argument turned out to be. Reading a Feature file is the caller's job, because it
/// needs codecs that live above this module.
pub trait FeatureIntervals {
    /// `FeatureManager.isFeatureFile` composed with `IntervalUtils.featureFileToIntervals`.
    ///
    /// `None` means "not a Feature file", which is what sends the argument down the interval-file
    /// and then the literal branch. It is not an error.
    fn intervals_from_feature_file(
        &self,
        path: &std::path::Path,
        header: &SamHeader,
    ) -> Option<Result<Vec<SimpleInterval>, IntervalArgumentError>>;
}

/// The seam with nothing plugged into it: no argument is ever a Feature file.
pub struct NoFeatureSources;

impl FeatureIntervals for NoFeatureSources {
    fn intervals_from_feature_file(
        &self,
        _path: &std::path::Path,
        _header: &SamHeader,
    ) -> Option<Result<Vec<SimpleInterval>, IntervalArgumentError>> {
        None
    }
}

/// `IntervalUtils.parseIntervalArguments(parser, arg)`: one argument to a list of intervals.
///
/// The order of the tests is the reference's and is observable. A Feature file is recognised
/// first, so a `.list` that also parses as a Feature file goes down the Feature path; the
/// interval-file test comes second and *throws* when the extension matches but the file is
/// missing; only then does an existing file that is neither become an error; and only a
/// non-existent argument with neither extension is parsed as a literal interval.
pub fn parse_interval_arguments(
    query: &str,
    header: &SamHeader,
    features: &dyn FeatureIntervals,
) -> Result<Vec<SimpleInterval>, IntervalArgumentError> {
    if query.contains(';') {
        return Err(IntervalArgumentError::LegacySemicolonSyntax(
            query.to_string(),
        ));
    }
    let path = std::path::Path::new(query);
    if path.exists() {
        if let Some(result) = features.intervals_from_feature_file(path, header) {
            return result;
        }
    }
    if has_gatk_interval_file_extension(query) {
        if !path.exists() {
            return Err(IntervalArgumentError::IntervalFileMissing(
                query.to_string(),
            ));
        }
        let text = std::fs::read_to_string(path)
            .map_err(|_| IntervalArgumentError::IntervalFileMissing(query.to_string()))?;
        return gatk_interval_file_to_list(&text, header);
    }
    if path.exists() {
        return Err(IntervalArgumentError::FileIsNeitherFeaturesNorIntervals(
            query.to_string(),
        ));
    }
    Ok(vec![interval::parse_interval(query, header)?])
}

/// `IntervalUtils.loadIntervals`: parse each argument, pad it, fold it in under the set rule, then
/// sort and merge the whole.
///
/// `-L unmapped` is reported separately rather than parsed, because `GenomeLoc.UNMAPPED` compares
/// and merges by identity upstream and would be meaningless in a sorted set here.
///
/// Padding is applied to the whole batch one argument produced, which matters for a file: an
/// interval file of twenty lines is padded and `ALL`-merged as one unit before the set rule sees
/// it, exactly as if its twenty lines had been twenty `-L` arguments would *not* be.
pub fn load_intervals_with_features(
    queries: &[String],
    header: &SamHeader,
    set_rule: SetRule,
    merging_rule: MergingRule,
    padding: i32,
    features: &dyn FeatureIntervals,
) -> Result<(Vec<SimpleInterval>, bool), IntervalArgumentError> {
    let mut all: Vec<SimpleInterval> = Vec::new();
    let mut unmapped = false;
    for query in queries {
        if is_unmapped_string(query) {
            unmapped = true;
            continue;
        }
        let mut parsed = parse_interval_arguments(query, header, features)?;
        if padding > 0 {
            parsed = with_flanks(parsed, padding, header);
        }
        all = merge_lists_by_set_operator(parsed, all, set_rule, header)?;
    }
    Ok((sort_and_merge(all, header, merging_rule), unmapped))
}

/// [`load_intervals_with_features`] with no Feature sources plugged in.
pub fn load_intervals(
    queries: &[String],
    header: &SamHeader,
    set_rule: SetRule,
    merging_rule: MergingRule,
    padding: i32,
) -> Result<(Vec<SimpleInterval>, bool), IntervalArgumentError> {
    load_intervals_with_features(
        queries,
        header,
        set_rule,
        merging_rule,
        padding,
        &NoFeatureSources,
    )
}

/// `GenomeLocSortedSet.createSetFromSequenceDictionary`: one interval per contig, whole.
pub fn whole_reference(header: &SamHeader) -> Vec<SimpleInterval> {
    header
        .sequences
        .iter()
        .map(|s| SimpleInterval {
            contig: s.name.clone(),
            start: 1,
            end: s.length,
        })
        .collect()
}

/// `GenomeLocSortedSet.subtractRegions`, stack walk included.
///
/// Written as the reference wrote it rather than as a sweep, because the shape is observable: when
/// an excluded interval sits strictly inside an included one, `GenomeLoc.subtract` returns the
/// *after* piece before the *before* piece, both are pushed back onto the processing stack, and
/// the before piece is therefore reconsidered first. A sweep that emitted them in coordinate order
/// would agree here by accident and diverge as soon as a second exclusion overlaps the tail.
pub fn subtract_regions(
    include: &[SimpleInterval],
    exclude: &[SimpleInterval],
    header: &SamHeader,
) -> Vec<SimpleInterval> {
    let key = |i: &SimpleInterval| contig_index(header, &i.contig).unwrap_or(usize::MAX);
    let mut good: Vec<SimpleInterval> = Vec::new();
    // Reversed, so that `pop` yields the first interval: these are the reference's two stacks.
    let mut to_process: Vec<SimpleInterval> = include.to_vec();
    to_process.reverse();
    let mut to_exclude: Vec<SimpleInterval> = exclude.to_vec();
    to_exclude.reverse();

    while let Some(p) = to_process.last().cloned() {
        let Some(e) = to_exclude.last().cloned() else {
            to_process.reverse();
            good.extend(to_process);
            return good;
        };
        let (pk, ek) = (key(&p), key(&e));
        let overlaps = pk == ek && p.start <= e.end && e.start <= p.end;
        if overlaps {
            to_process.pop();
            // `GenomeLoc.subtract`: the after piece first, then the before piece.
            if p.start >= e.start && p.end <= e.end {
                // `e` contains `p` entirely, including the equal case: nothing survives.
            } else if e.start >= p.start && e.end <= p.end {
                // `afterStop - afterStart >= 0` upstream, with afterStart = e.end + 1.
                if p.end > e.end {
                    to_process.push(SimpleInterval {
                        contig: p.contig.clone(),
                        start: e.end + 1,
                        end: p.end,
                    });
                }
                // `beforeStop - beforeStart >= 0` upstream, with beforeStop = e.start - 1.
                if e.start > p.start {
                    to_process.push(SimpleInterval {
                        contig: p.contig.clone(),
                        start: p.start,
                        end: e.start - 1,
                    });
                }
            } else if e.start < p.start {
                to_process.push(SimpleInterval {
                    contig: p.contig.clone(),
                    start: e.end + 1,
                    end: p.end,
                });
            } else {
                to_process.push(SimpleInterval {
                    contig: p.contig.clone(),
                    start: p.start,
                    end: e.start - 1,
                });
            }
        } else if pk < ek {
            good.push(to_process.pop().unwrap());
        } else if pk > ek {
            to_exclude.pop();
        } else if p.end < e.start {
            good.push(to_process.pop().unwrap());
        } else {
            // `e` ends before `p` starts, so it can never affect anything still to come.
            to_exclude.pop();
        }
    }
    good
}

/// `IntervalArgumentCollection.parseIntervals`: the whole pipeline, from argument strings to the
/// list `IntervalWalker` iterates.
pub fn parse_intervals(
    arguments: &IntervalArguments,
    header: &SamHeader,
) -> Result<TraversalParameters, IntervalArgumentError> {
    if !arguments.specified() {
        return Err(IntervalArgumentError::NoIntervalsSpecified);
    }

    // No -L but a -XL: the include set is the entire reference, so that -XL has something to
    // subtract from. That set is built from the dictionary and is therefore already sorted.
    let (include, unmapped) = if arguments.include.is_empty() {
        (whole_reference(header), false)
    } else {
        load_intervals(
            &arguments.include,
            header,
            arguments.set_rule,
            arguments.merging_rule,
            arguments.padding,
        )?
    };

    let (exclude, exclude_unmapped) = load_intervals(
        &arguments.exclude,
        header,
        SetRule::Union,
        arguments.merging_rule,
        arguments.exclusion_padding,
    )?;
    if exclude_unmapped {
        return Err(IntervalArgumentError::UnmappedExcluded);
    }

    let intervals = if exclude.is_empty() {
        include
    } else {
        let mut remaining = subtract_regions(&include, &exclude, header);
        if remaining.is_empty() {
            return Err(IntervalArgumentError::ExclusionRemovedEverything);
        }
        // `createSetFromList` sorts what the stack walk produced and does **not** merge it: its
        // `add` throws on an overlap rather than merging one, and abutting intervals are inserted
        // side by side. Merging here would silently join two `-L` arguments that the user asked to
        // keep apart with `OVERLAPPING_ONLY`, on the sole grounds that a `-XL` elsewhere existed.
        remaining.sort_by_key(|i| {
            (
                contig_index(header, &i.contig).unwrap_or(usize::MAX),
                i.start,
                i.end,
            )
        });
        remaining
    };

    Ok(TraversalParameters {
        intervals,
        traverse_unmapped: unmapped,
    })
}
