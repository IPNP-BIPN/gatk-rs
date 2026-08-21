//! `VcfToIntervalList`, ported from `picard.vcf.VcfToIntervalList` (Picard 3.4.0) and the two
//! pieces of htsjdk it leans on: `VCFFileReader.toIntervals` and
//! `IntervalList.IntervalMergerIterator`.
//!
//! A VCF read as a stream of intervals and written back out as a Picard interval list, merged on
//! the way through.
//!
//! # The merging is a stream and never sorts
//!
//! The iterator holds one interval and compares the next with it, so a file whose records are out
//! of order comes out out of order, and two intervals that would have merged do not when a third
//! stands between them. Nothing here sorts, uniques or checks: the tool inherits whatever order
//! the file had.
//!
//! # Abutting intervals merge
//!
//! ```java
//! current.overlaps(next) || (combineAbuttingIntervals && current.withinDistanceOf(next, 1))
//! ```
//!
//! `combineAbuttingIntervals` is true here, so 50-50 and 51-51 are one interval while 60-60 and
//! 62-62 are two.
//!
//! # The name counter counts unnamed records, and it counts them after the filtering
//!
//! An unnamed record is `interval-<n>` where `n` increments only when a record has no ID, so the
//! number is not the record's position in the file. And the stream filters before it maps, so
//! `INCLUDE_FILTERED` moves the numbers of every unnamed interval after a filtered one.
//!
//! # `INCLUDE_FILTERED` changes which intervals merge
//!
//! A filtered record between two others is the bridge that joins them, so keeping it is not merely
//! one interval more: it is one interval fewer where the bridge closed a gap.
//!
//! # The two ID methods do not only differ in the name
//!
//! `CONCAT_ALL` rebuilds the interval from the group with `min(start)` and `max(end)`; `USE_FIRST`
//! keeps `current.start`, which is the first member's start and is never lowered. Under a sorted
//! file the two agree. Under a file where a later record starts earlier than the one it overlaps
//! they do not, and this port follows the Java on both branches. The golden does not pin that
//! difference down: its unsorted records do not overlap.

use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::{Value, VariantContext};

/// One sequence of the dictionary, as far as an `@SQ` line reads it. `VCFContigHeaderLine`
/// carries the ID, the length and `assembly` into a `SAMSequenceRecord` and nothing else.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub length: i64,
    pub assembly: Option<String>,
}

impl Sequence {
    /// `SAMTextHeaderCodec.getSQLine`: SN and LN, then the attributes, of which there is at most
    /// one here.
    pub fn line(&self) -> String {
        let mut text = format!("@SQ\tSN:{}\tLN:{}", self.name, self.length);
        if let Some(assembly) = &self.assembly {
            text.push_str(&format!("\tAS:{assembly}"));
        }
        text
    }
}

/// One line of the interval list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    /// `null` prints as `.`, which this tool never produces: every interval is named.
    pub name: Option<String>,
}

/// `VARIANT_ID_METHOD`, which is a static field in the Java and so leaks between runs in one JVM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdMethod {
    ConcatAll,
    UseFirst,
}

impl IdMethod {
    /// `concatenate_ids`, the only thing the tool reads it for.
    pub fn concatenates(&self) -> bool {
        matches!(self, IdMethod::ConcatAll)
    }
}

/// What the run can fail with, neither of which the tool raises itself.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConvertError {
    /// `getSequenceDictionary()` is null when the header declares no contigs, and the codec walks
    /// into it.
    NullDictionary,
    /// The reader refused.
    Vcf(String, String),
}

impl ConvertError {
    pub fn java_class(&self) -> &str {
        match self {
            ConvertError::NullDictionary => "java.lang.NullPointerException",
            ConvertError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            ConvertError::NullDictionary => {
                "Cannot invoke \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because \
                 the return value of \"htsjdk.samtools.SAMFileHeader.getSequenceDictionary()\" is \
                 null"
                    .to_string()
            }
            ConvertError::Vcf(_, message) => message.clone(),
        }
    }
}

/// `UNKNOWN_SEQUENCE_LENGTH`, which a contig line with no length becomes.
pub const UNKNOWN_SEQUENCE_LENGTH: i64 = 0;

/// `VCFHeader.getSequenceDictionary()`, which is `None` when there are no contig lines at all.
pub fn dictionary(header: &VcfHeader) -> Option<Vec<Sequence>> {
    let field = |fields: &[(String, String)], key: &str| {
        fields
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.clone())
    };
    let sequences: Vec<Sequence> = header
        .lines
        .iter()
        .filter_map(|line| match line {
            HeaderLine::Contig { fields, .. } => Some(Sequence {
                name: field(fields, "ID").unwrap_or_default(),
                length: field(fields, "length")
                    .and_then(|value| value.parse().ok())
                    .unwrap_or(UNKNOWN_SEQUENCE_LENGTH),
                assembly: field(fields, "assembly"),
            }),
            _ => None,
        })
        .collect();
    if sequences.is_empty() {
        None
    } else {
        Some(sequences)
    }
}

/// `VariantContext.isFiltered()`: filters applied and not empty. A `PASS` is applied and empty, so
/// it is not filtered.
fn is_filtered(record: &VariantContext) -> bool {
    record
        .filters
        .as_ref()
        .is_some_and(|filters| !filters.is_empty())
}

/// `getCommonInfo().getAttributeAsInt(END, vc.getEnd())`. `stop` already carries the END the
/// decoder read, so this only matters for a value the decoder could not use.
fn end_of(record: &VariantContext) -> i64 {
    record
        .attributes
        .iter()
        .find(|(key, _)| key == "END")
        .and_then(|(_, value)| match value {
            Value::Int(number) => Some(*number),
            Value::Str(text) => text.parse().ok(),
            _ => None,
        })
        .unwrap_or(record.stop)
}

/// `VCFFileReader.toIntervals`, whose counter increments only on records with no ID and only on
/// records the filter let through.
pub fn to_intervals(records: &[VariantContext], include_filtered: bool) -> Vec<Interval> {
    let mut count = 0;
    records
        .iter()
        .filter(|record| include_filtered || !is_filtered(record))
        .map(|record| {
            let name = if record.id == "." || record.id.is_empty() {
                count += 1;
                format!("interval-{count}")
            } else {
                record.id.clone()
            };
            Interval {
                contig: record.contig.clone(),
                start: record.start,
                end: end_of(record),
                name: Some(name),
            }
        })
        .collect()
}

/// `IntervalList.merge(intervals, true)`: the names are a `LinkedHashSet` joined with a pipe, and
/// the bounds are the group's minimum and maximum.
fn merge_group(group: &[Interval]) -> Interval {
    let first = group.first().expect("a group is never empty");
    let mut names: Vec<String> = Vec::new();
    let mut start = first.start;
    let mut end = first.end;
    for interval in group {
        if let Some(name) = &interval.name {
            if !names.contains(name) {
                names.push(name.clone());
            }
        }
        start = start.min(interval.start);
        end = end.max(interval.end);
    }
    Interval {
        contig: first.contig.clone(),
        start,
        end,
        name: if names.is_empty() {
            None
        } else {
            Some(names.join("|"))
        },
    }
}

/// `IntervalList.IntervalMergerIterator(intervals, true, false, concatenate_names)`, drained.
pub fn merge_intervals(intervals: &[Interval], concatenate_names: bool) -> Vec<Interval> {
    let mut out: Vec<Interval> = Vec::new();
    // `current`, which is a MutableFeature: its start is the first member's and only its end moves.
    let mut current: Option<Interval> = None;
    let mut group: Vec<Interval> = Vec::new();

    let emit = |current: &Interval, group: &[Interval]| {
        if concatenate_names {
            merge_group(group)
        } else {
            current.clone()
        }
    };

    for next in intervals {
        match &mut current {
            None => {
                if concatenate_names {
                    group.push(next.clone());
                }
                current = Some(next.clone());
            }
            Some(held) => {
                // `overlaps` or `withinDistanceOf(next, 1)`, both of which need the same contig.
                let touching = held.contig == next.contig
                    && held.start <= next.end + 1
                    && next.start <= held.end + 1;
                if touching {
                    if concatenate_names {
                        group.push(next.clone());
                    }
                    held.end = held.end.max(next.end);
                } else {
                    out.push(emit(held, &group));
                    group.clear();
                    // `current.setAll(next)`, which takes the contig, the start and the end.
                    *held = next.clone();
                    if concatenate_names {
                        group.push(next.clone());
                    }
                }
            }
        }
    }
    if let Some(held) = &current {
        out.push(emit(held, &group));
    }
    out
}

/// The interval list the writer produces: the header the codec writes, then five columns per
/// interval.
pub fn write_list(sequences: &[Sequence], intervals: &[Interval]) -> String {
    let mut text = String::from("@HD\tVN:1.6\n");
    for sequence in sequences {
        text.push_str(&sequence.line());
        text.push('\n');
    }
    for interval in intervals {
        text.push_str(&format!(
            "{}\t{}\t{}\t+\t{}\n",
            interval.contig,
            interval.start,
            interval.end,
            interval.name.as_deref().unwrap_or(".")
        ));
    }
    text
}

/// `doWork()`: the whole run, text in and text out.
pub fn convert(
    input: &str,
    include_filtered: bool,
    method: IdMethod,
) -> Result<String, ConvertError> {
    let file = read_vcf(input).map_err(|failure| {
        ConvertError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    // The writer is opened with the dictionary before a single record is read, so a header with no
    // contigs fails before any interval is built.
    let sequences = dictionary(&file.header).ok_or(ConvertError::NullDictionary)?;
    let intervals = to_intervals(&file.records, include_filtered);
    let merged = merge_intervals(&intervals, method.concatenates());
    Ok(write_list(&sequences, &merged))
}
