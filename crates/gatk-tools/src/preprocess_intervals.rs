//! `PreprocessIntervals`, ported from
//! `org.broadinstitute.hellbender.tools.copynumber.PreprocessIntervals` (GATK 4.6.2.0).
//!
//! The bins every copy-number tool counts over: the input intervals padded, de-overlapped, split
//! into bins, and stripped of the bins that are all N.
//!
//! # The overlap is resolved at the midpoint of the originals
//!
//! ```java
//! final int newThisEnd = (originalThisEnd + originalNextStart) / 2;
//! final int newNextStart = newThisEnd + 1;
//! ```
//!
//! Not the midpoint of the padded intervals, and not a merge: two intervals whose padding runs into
//! one another are cut apart between the intervals they came from. The division is integer and
//! rounds toward zero, so an odd sum leans left: `(20 + 41) / 2` is 30, and the first interval ends
//! at 30 rather than at 31.
//!
//! The midpoint cannot fall outside either interval however large the padding is. The inputs
//! arrive sorted and non-overlapping, so the midpoint is at least the first's end and at most the
//! second's start, and the `SimpleInterval` the resolution builds never throws.
//!
//! # The pass is sequential and in place
//!
//! Each step compares interval `i` against interval `i + 1` and rewrites both, and the next step
//! reads the interval it just rewrote. Three intervals padded into one another therefore resolve
//! left to right, and the ORIGINAL of the middle one is still what the second midpoint uses.
//!
//! # The bins are laid from the start
//!
//! `for (binStart = start; binStart <= end; binStart += binLength)`, so the short bin is the last
//! one and an interval shorter than a bin is one bin of its own length. A bin length of zero means
//! no binning at all.
//!
//! # The N filter is `allMatch`
//!
//! A bin survives on one non-N base, and `Nucleotide.decode` is case-insensitive, so a stretch of
//! lower-case `n` is dropped exactly like an upper-case one. A whole contig can come back with no
//! bins at all, which is a list with a header and nothing under it.

use crate::filter_intervals::Interval;

/// The tool's defaults.
pub const DEFAULT_BIN_LENGTH: i32 = 1000;
pub const DEFAULT_PADDING: i32 = 250;

/// A contig, as the sequence dictionary carries it.
///
/// `md5` and `uri` are the `M5` and `UR` fields, which come from the `.dict` the reference indexer
/// wrote rather than from this tool: it copies the dictionary into the interval list untouched.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub length: i32,
    pub md5: Option<String>,
    pub uri: Option<String>,
}

impl Sequence {
    /// The `@SQ` line, with the two optional fields in the order htsjdk writes them.
    pub fn line(&self) -> String {
        let mut text = format!("@SQ\tSN:{}\tLN:{}", self.name, self.length);
        if let Some(md5) = &self.md5 {
            text.push_str(&format!("\tM5:{md5}"));
        }
        if let Some(uri) = &self.uri {
            text.push_str(&format!("\tUR:{uri}"));
        }
        text
    }
}

/// The Picard interval list this tool writes: the dictionary, then five columns per bin.
pub fn write_list(sequences: &[Sequence], intervals: &[Interval]) -> String {
    let mut text = String::from("@HD\tVN:1.6\n");
    for sequence in sequences {
        text.push_str(&sequence.line());
        text.push('\n');
    }
    for interval in intervals {
        text.push_str(&format!(
            "{}\t{}\t{}\t+\t.\n",
            interval.contig, interval.start, interval.end
        ));
    }
    text
}

/// What the tool refuses before it does anything.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PreprocessError {
    /// `CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection`.
    MergingRule,
    IntervalPadding,
    IntervalExclusionPadding,
    /// The parser's own bounds, which are not the tool's code at all.
    BinLengthOutOfRange(i32),
    PaddingOutOfRange(i32),
}

impl PreprocessError {
    pub fn java_class(&self) -> &'static str {
        match self {
            PreprocessError::BinLengthOutOfRange(_) | PreprocessError::PaddingOutOfRange(_) => {
                "org.broadinstitute.barclay.argparser.CommandLineException$OutOfRangeArgumentValue"
            }
            _ => "java.lang.IllegalArgumentException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            PreprocessError::MergingRule => {
                "Interval merging rule must be set to OVERLAPPING_ONLY.".to_string()
            }
            PreprocessError::IntervalPadding => "Interval padding must be set to 0.".to_string(),
            PreprocessError::IntervalExclusionPadding => {
                "Interval exclusion padding must be set to 0.".to_string()
            }
            PreprocessError::BinLengthOutOfRange(value) => {
                format!("Argument bin-length has a bad value: {value}. minimum allowed value 0")
            }
            PreprocessError::PaddingOutOfRange(value) => {
                format!("Argument padding has a bad value: {value}. minimum allowed value 0")
            }
        }
    }
}

/// `getAllIntervalsForReference`: one interval per contig, whole.
pub fn whole_reference(sequences: &[Sequence]) -> Vec<Interval> {
    sequences
        .iter()
        .map(|sequence| Interval {
            contig: sequence.name.clone(),
            start: 1,
            end: sequence.length,
        })
        .collect()
}

/// `padIntervals`: pad, then cut the overlaps apart at the originals' midpoints.
pub fn pad_intervals(inputs: &[Interval], padding: i32, sequences: &[Sequence]) -> Vec<Interval> {
    let length_of = |contig: &str| {
        sequences
            .iter()
            .find(|sequence| sequence.name == contig)
            .map(|sequence| sequence.length)
            .expect("the dictionary holds every contig the intervals name")
    };
    let mut padded: Vec<Interval> = inputs
        .iter()
        .map(|interval| Interval {
            contig: interval.contig.clone(),
            start: 1.max(interval.start - padding),
            end: (interval.end + padding).min(length_of(&interval.contig)),
        })
        .collect();
    for index in 0..padded.len().saturating_sub(1) {
        // `SimpleInterval.overlaps`, which is false across contigs.
        let overlaps = padded[index].contig == padded[index + 1].contig
            && padded[index].start <= padded[index + 1].end
            && padded[index + 1].start <= padded[index].end;
        if !overlaps {
            continue;
        }
        // The ORIGINALS, not the padded intervals, and not the interval this pass may already
        // have rewritten.
        let midpoint = (inputs[index].end + inputs[index + 1].start) / 2;
        padded[index].end = midpoint;
        padded[index + 1].start = midpoint + 1;
    }
    padded
}

/// `generateBins`.
pub fn generate_bins(intervals: &[Interval], bin_length: i32) -> Vec<Interval> {
    if bin_length == 0 {
        return intervals.to_vec();
    }
    let mut bins = Vec::new();
    for interval in intervals {
        let mut start = interval.start;
        while start <= interval.end {
            bins.push(Interval {
                contig: interval.contig.clone(),
                start,
                end: (start + bin_length - 1).min(interval.end),
            });
            start += bin_length;
        }
    }
    bins
}

/// `Nucleotide.decode(b) == Nucleotide.N`, which is case-insensitive and holds for nothing else.
pub fn is_n(base: u8) -> bool {
    base.eq_ignore_ascii_case(&b'N')
}

/// `filterBinsContainingOnlyNs`: a bin survives on one base that is not an N.
///
/// `bases` is the whole contig, one-based positions indexing from 1, which is what a reference
/// query hands back a slice of.
pub fn filter_bins_containing_only_ns(
    bins: &[Interval],
    bases: impl Fn(&str) -> Vec<u8>,
) -> Vec<Interval> {
    bins.iter()
        .filter(|bin| {
            let contig = bases(&bin.contig);
            let from = (bin.start - 1) as usize;
            let to = (bin.end as usize).min(contig.len());
            !contig[from..to].iter().all(|base| is_n(*base))
        })
        .cloned()
        .collect()
}

/// The whole run: the intervals in, the interval list out.
///
/// `inputs` is `None` for a run with no `-L` at all, which is the whole reference.
pub fn preprocess(
    inputs: Option<&[Interval]>,
    sequences: &[Sequence],
    bin_length: i32,
    padding: i32,
    bases: impl Fn(&str) -> Vec<u8>,
) -> Result<String, PreprocessError> {
    if bin_length < 0 {
        return Err(PreprocessError::BinLengthOutOfRange(bin_length));
    }
    if padding < 0 {
        return Err(PreprocessError::PaddingOutOfRange(padding));
    }
    let whole;
    let inputs = match inputs {
        Some(intervals) => intervals,
        None => {
            whole = whole_reference(sequences);
            &whole
        }
    };
    let padded = pad_intervals(inputs, padding, sequences);
    let bins = generate_bins(&padded, bin_length);
    let kept = filter_bins_containing_only_ns(&bins, bases);
    Ok(write_list(sequences, &kept))
}
