//! Ported from `org.broadinstitute.hellbender.utils.SimpleInterval`, `GenomeLoc`,
//! `GenomeLocParser`, `IntervalUtils` and `htsjdk.samtools.util.OverlapDetector` (GATK 4.6.2.0).
//!
//! Intervals are 1-based and **closed at both ends**, which is the first thing to get wrong: a
//! half-open port agrees on every interval of length one and disagrees by a base everywhere else,
//! and a one-base error at an interval edge changes which reads a walker sees.
//!
//! # Why the parser is not a `split(':')`
//!
//! A contig name may itself contain a colon, so `HLA-A*01:01:01:01:1-100` is a real query, and
//! `chr1:100` is ambiguous the moment a contig is literally named `chr1:100`. GATK resolves a
//! query against the sequence dictionary and produces *every* valid interpretation, then refuses
//! the query if there is more than one. `IntervalUtils.getResolvedIntervals` is that resolution,
//! and it splits on the **last** colon, not the first.

use std::collections::HashMap;

use htsjdk_bam::header::SamHeader;

/// `SimpleInterval`: a contig with a 1-based, closed `[start, end]`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SimpleInterval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

/// `SimpleInterval.END_OF_CONTIG`.
const END_OF_CONTIG: char = '+';
const CONTIG_SEPARATOR: char = ':';
const START_END_SEPARATOR: char = '-';

impl SimpleInterval {
    /// `SimpleInterval.isValid`.
    pub fn is_valid(start: i32, end: i32) -> bool {
        start > 0 && end >= start
    }

    /// `SimpleInterval.parsePositionThrowOnFailure`: commas are stripped, so `1,000,000` parses.
    pub fn parse_position(text: &str) -> Option<i32> {
        text.replace(',', "").parse().ok()
    }

    pub fn overlaps(&self, contig: &str, start: i32, end: i32) -> bool {
        self.contig == contig && self.start <= end && start <= self.end
    }
}

/// What a query string can resolve to, kept apart from a plain failure.
#[derive(Debug, PartialEq, Eq)]
pub enum ParseError {
    /// The query names no contig the dictionary declares.
    UnknownContig(String),
    /// Two valid readings, for instance a contig literally named like an interval query.
    Ambiguous(String),
    /// The positions are not a valid interval, or did not parse as numbers.
    MalformedPositions(String),
}

/// `IntervalUtils.getResolvedIntervals`: every valid reading of one query string.
///
/// The order matters and is the reference's: the whole-string-as-a-contig reading comes first, so
/// a query that is both a contig name and an interval query resolves to the contig when the second
/// reading turns out to be malformed.
pub fn resolved_intervals(query: &str, header: &SamHeader) -> Vec<SimpleInterval> {
    let mut resolved = Vec::new();
    if let Some(sequence) = header.sequences.iter().find(|s| s.name == query) {
        resolved.push(SimpleInterval {
            contig: query.to_string(),
            start: 1,
            end: sequence.length,
        });
    }

    // The last colon, not the first: a contig name may contain colons of its own.
    let Some(last_colon) = query.rfind(CONTIG_SEPARATOR) else {
        return resolved;
    };
    let prefix = &query[..last_colon];
    let Some(prefix_sequence) = header.sequences.iter().find(|s| s.name == prefix) else {
        return resolved;
    };

    let last_dash = query.rfind(START_END_SEPARATOR);
    let positions = if let Some(rest) = query.strip_suffix(END_OF_CONTIG) {
        // "prefix:nnn+" runs to the end of the contig.
        SimpleInterval::parse_position(&rest[last_colon + 1..])
            .map(|start| (start, prefix_sequence.length))
    } else if last_dash.is_some_and(|dash| dash > last_colon) {
        let dash = last_dash.unwrap();
        SimpleInterval::parse_position(&query[last_colon + 1..dash])
            .zip(SimpleInterval::parse_position(&query[dash + 1..]))
    } else {
        SimpleInterval::parse_position(&query[last_colon + 1..]).map(|start| (start, start))
    };

    // A number that does not parse is only an error when there is no other reading: a query that
    // is also a contig name survives its own malformed suffix, with a warning upstream.
    if let Some((start, end)) = positions {
        if SimpleInterval::is_valid(start, end) {
            resolved.push(SimpleInterval {
                contig: prefix.to_string(),
                start,
                end,
            });
        }
    }
    resolved
}

/// `GenomeLocParser.parseGenomeLoc`: one query string to one interval, or an error.
pub fn parse_interval(query: &str, header: &SamHeader) -> Result<SimpleInterval, ParseError> {
    let query = query.trim();
    let resolved = resolved_intervals(query, header);
    match resolved.len() {
        // The reference distinguishes these two: an unknown contig and a malformed position are
        // different messages, and the second only surfaces when nothing else parsed.
        0 => {
            let contig = query.split(CONTIG_SEPARATOR).next().unwrap_or(query);
            if header.sequences.iter().any(|s| s.name == contig) {
                Err(ParseError::MalformedPositions(query.to_string()))
            } else {
                Err(ParseError::UnknownContig(query.to_string()))
            }
        }
        1 => Ok(resolved.into_iter().next().unwrap()),
        _ => Err(ParseError::Ambiguous(query.to_string())),
    }
}

/// `IntervalMergingRule`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MergingRule {
    /// Merge overlapping *and* adjacent intervals. `IntervalOverlapReadFilter` uses this one.
    All,
    /// Merge only intervals that actually overlap.
    OverlappingOnly,
}

/// `IntervalUtils.loadIntervals` with `UNION` and no padding: parse, sort, merge.
///
/// The sort is `GenomeLoc.compareTo`, which orders by *contig index* and not by contig name, so
/// the dictionary's order is what decides. Sorting by name would put `chr10` before `chr2` and
/// merge two intervals that the reference keeps apart.
pub fn load_intervals(
    queries: &[String],
    header: &SamHeader,
    rule: MergingRule,
) -> Result<Vec<SimpleInterval>, ParseError> {
    let mut intervals = Vec::new();
    for query in queries {
        intervals.push(parse_interval(query, header)?);
    }
    let index = |contig: &str| {
        header
            .sequences
            .iter()
            .position(|s| s.name == contig)
            .unwrap_or(usize::MAX)
    };
    intervals.sort_by_key(|i| (index(&i.contig), i.start, i.end));
    Ok(merge_interval_locations(intervals, rule))
}

/// `IntervalUtils.mergeIntervalLocations` over a sorted list.
pub fn merge_interval_locations(
    intervals: Vec<SimpleInterval>,
    rule: MergingRule,
) -> Vec<SimpleInterval> {
    if intervals.len() <= 1 {
        return intervals;
    }
    let mut merged: Vec<SimpleInterval> = Vec::new();
    let mut iter = intervals.into_iter();
    let mut previous = iter.next().unwrap();
    for current in iter {
        // `overlapsP`, then `contiguousP`, which is the one that treats [1,10] and [11,20] as one
        // interval: discontinuousP compares `start - 1 > stop`, so adjacency counts as contiguous.
        let same_contig = previous.contig == current.contig;
        let overlaps =
            same_contig && previous.start <= current.end && current.start <= previous.end;
        let contiguous = same_contig
            && (previous.start - 1) <= current.end
            && (current.start - 1) <= previous.end;
        if overlaps || (contiguous && rule == MergingRule::All) {
            previous = SimpleInterval {
                contig: previous.contig,
                start: previous.start.min(current.start),
                end: previous.end.max(current.end),
            };
        } else {
            merged.push(previous);
            previous = current;
        }
    }
    merged.push(previous);
    merged
}

/// `htsjdk.samtools.util.OverlapDetector`, reduced to what the read filter asks of it.
///
/// A locatable whose contig is absent from the map is not an error: the tree lookup misses and
/// the answer is false. That is what makes an unmapped read, whose contig is null, fail the
/// filter rather than crash it.
pub struct OverlapDetector {
    by_contig: HashMap<String, Vec<SimpleInterval>>,
}

impl OverlapDetector {
    pub fn create(intervals: Vec<SimpleInterval>) -> OverlapDetector {
        let mut by_contig: HashMap<String, Vec<SimpleInterval>> = HashMap::new();
        for interval in intervals {
            by_contig
                .entry(interval.contig.clone())
                .or_default()
                .push(interval);
        }
        OverlapDetector { by_contig }
    }

    /// `OverlapDetector.overlapsAny`.
    pub fn overlaps_any(&self, contig: Option<&str>, start: i32, end: i32) -> bool {
        let Some(contig) = contig else {
            return false;
        };
        let Some(intervals) = self.by_contig.get(contig) else {
            return false;
        };
        // The reference bails out before querying the tree when the locatable is empty, which a
        // read whose cigar consumes no reference is.
        if start > end {
            return false;
        }
        intervals
            .iter()
            .any(|interval| interval.overlaps(contig, start, end))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::SequenceRecord;

    fn header() -> SamHeader {
        let mut header = SamHeader::new();
        header.sequences.push(SequenceRecord::new("chr1", 2000));
        header.sequences.push(SequenceRecord::new("chr2", 1000));
        header
    }

    #[test]
    fn a_bare_contig_is_the_whole_contig() {
        assert_eq!(
            parse_interval("chr2", &header()),
            Ok(SimpleInterval {
                contig: "chr2".to_string(),
                start: 1,
                end: 1000
            })
        );
    }

    #[test]
    fn a_single_position_is_a_one_base_interval() {
        let parsed = parse_interval("chr1:100", &header()).unwrap();
        assert_eq!((parsed.start, parsed.end), (100, 100));
    }

    /// `chr1:1900+` runs to the contig's end, which the dictionary supplies.
    #[test]
    fn a_trailing_plus_runs_to_the_end_of_the_contig() {
        let parsed = parse_interval("chr1:1900+", &header()).unwrap();
        assert_eq!((parsed.start, parsed.end), (1900, 2000));
    }

    #[test]
    fn commas_are_stripped_from_positions() {
        let parsed = parse_interval("chr1:1,000-1,200", &header()).unwrap();
        assert_eq!((parsed.start, parsed.end), (1000, 1200));
    }

    /// A contig whose name contains a colon: the split is on the *last* colon.
    #[test]
    fn a_contig_name_may_contain_colons() {
        let mut header = header();
        header
            .sequences
            .push(SequenceRecord::new("HLA-A*01:01:01:01", 3503));
        let parsed = parse_interval("HLA-A*01:01:01:01:100-200", &header).unwrap();
        assert_eq!(parsed.contig, "HLA-A*01:01:01:01");
        assert_eq!((parsed.start, parsed.end), (100, 200));

        // And the bare name still resolves to the whole contig.
        let whole = parse_interval("HLA-A*01:01:01:01", &header).unwrap();
        assert_eq!((whole.start, whole.end), (1, 3503));
    }

    /// Two readings of one string is a refusal, not a preference.
    #[test]
    fn an_ambiguous_query_is_refused() {
        let mut header = header();
        header.sequences.push(SequenceRecord::new("chr1:100", 500));
        assert_eq!(
            parse_interval("chr1:100", &header),
            Err(ParseError::Ambiguous("chr1:100".to_string()))
        );
    }

    /// Adjacent intervals merge under ALL and stay apart under OVERLAPPING_ONLY.
    #[test]
    fn adjacency_is_a_merge_only_under_the_all_rule() {
        let queries = ["chr1:1-10".to_string(), "chr1:11-20".to_string()];
        let all = load_intervals(&queries, &header(), MergingRule::All).unwrap();
        assert_eq!(all.len(), 1);
        assert_eq!((all[0].start, all[0].end), (1, 20));

        let strict = load_intervals(&queries, &header(), MergingRule::OverlappingOnly).unwrap();
        assert_eq!(strict.len(), 2);
    }

    #[test]
    fn an_unmapped_locatable_overlaps_nothing() {
        let detector = OverlapDetector::create(vec![SimpleInterval {
            contig: "chr1".to_string(),
            start: 1,
            end: 100,
        }]);
        assert!(detector.overlaps_any(Some("chr1"), 50, 60));
        assert!(!detector.overlaps_any(Some("chr2"), 50, 60));
        assert!(!detector.overlaps_any(None, 50, 60));
        // An empty locatable: the reference checks this before touching the tree.
        assert!(!detector.overlaps_any(Some("chr1"), 60, 50));
    }
}
