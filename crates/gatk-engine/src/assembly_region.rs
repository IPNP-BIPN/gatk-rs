//! `AssemblyRegion`, ported from `org.broadinstitute.hellbender.engine.AssemblyRegion`
//! (GATK 4.6.2.0).
//!
//! [`crate::activity_profile`] decides *where* a region starts and stops. This is what the region
//! **is** once it exists, and it is the object `HaplotypeCaller` and `Mutect2` assemble: two spans
//! and the reads overlapping the wider one. The primary span is where variants get called; the
//! padded span is what gets assembled, because assembling over more territory improves the calls
//! inside the primary span.
//!
//! # Four things that are not what they look like
//!
//! **The padding constructor can fail with a message about null.** `makePaddedSpan` goes through
//! `IntervalUtils.trimIntervalToContig`, which returns `null` rather than throwing when the padded
//! interval cannot be placed on the contig at all. The null then reaches `Utils.nonNull` inside the
//! other constructor, so the reported failure is a null padded span and never mentions padding.
//! Here that is [`RegionError::NullPaddedSpan`].
//!
//! **`trim(span, padding)` is not what its own javadoc describes.** The javadoc works an example:
//! active 100-200 with padding 50, so the true span is 50-250; trimmed to 150-225, it says "here we
//! represent the assembly region as a region from 150-200 with 25 bp of padding". The code does
//! something else. It expands the **requested** span by the requested padding and intersects that
//! with the old padded span, so the padding is never recomputed to fit and the answer is a
//! different interval. The golden is the arbiter, and the conformance suite carries that exact
//! example so the disagreement is a row rather than a footnote.
//!
//! **Trimming reorders the reads.** Every read is re-clipped to the new padded span, the ones left
//! empty or no longer overlapping are dropped, and what survives is **sorted** with
//! `ReadCoordinateComparator`. So the read order of a trimmed region is not the order the reads
//! were added in, and the comparator is part of the region's observable output rather than an
//! internal detail. It lives in [`crate::read_utils::compare_read_coordinate`].
//!
//! **Trimming drops the hard-clipped pileup reads.** The new region is constructed empty and only
//! `addAll` is called on it, so the second read list does not survive a trim.
//!
//! # Why the header is an argument here and a field upstream
//!
//! Upstream the region holds the `SAMFileHeader`, and uses it for two things: contig lengths when
//! padding, and the comparator when trimming. Holding a borrow of the header inside the struct
//! would put a lifetime on every region and on everything that stores one, so it is passed to the
//! operations that need it instead. That is a Rust shape, not a behavioural difference: the same
//! header reaches the same three call sites.

use crate::clipping::{self, ClipError};
use crate::interval::{trim_interval_to_contig, SimpleInterval};
use crate::read;
use crate::read_utils;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// What the reference throws, kept apart so a port cannot collapse two different refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RegionError {
    /// `UserException.MissingContigInSequenceDictionary`, from the padding constructor.
    MissingContig { contig: String, dictionary: String },
    /// `Utils.nonNull` on a padded span that `trimIntervalToContig` returned as `null`.
    NullPaddedSpan,
    /// `Utils.validate(paddedSpan.contains(activeSpan), ...)`.
    PaddedDoesNotContainActive,
    /// `Utils.validateArg(padding >= 0, ...)` in `trim`.
    NegativePadding,
    /// `Utils.validateArg(paddedSpan.contains(span), ...)` in the two-span `trim`.
    RequestedPaddedDoesNotContain,
    /// `SimpleInterval::intersect` refusing two intervals that do not overlap.
    NoOverlapToIntersect { left: String, right: String },
    /// `new SimpleInterval(read)` refusing the read's own coordinates. This is how an unmapped
    /// read is rejected: on the interval's validation, before the overlap is ever tested, and the
    /// message never says the read is unmapped.
    InvalidReadInterval {
        contig: String,
        start: i32,
        end: i32,
    },
    /// The read does not overlap the padded span.
    ReadDoesNotOverlap { read_loc: String, padded: String },
    /// A read on a different contig from the ones already in the region.
    ReadOnDifferentContig { last: String, read: String },
    /// A read starting before the last one added.
    ReadOutOfOrder {
        last: String,
        last_start: i32,
        read: String,
        read_start: i32,
    },
    /// Clipping failed while trimming.
    Clip(ClipError),
}

impl RegionError {
    /// The Java class the reference throws, as the dump prints it.
    pub fn class(&self) -> &'static str {
        match self {
            RegionError::MissingContig { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$MissingContigInSequenceDictionary"
            }
            RegionError::Clip(_) => "java.lang.IllegalArgumentException",
            _ => "java.lang.IllegalArgumentException",
        }
    }

    /// The message the reference produces, verbatim.
    pub fn message(&self) -> String {
        match self {
            RegionError::MissingContig { contig, dictionary } => {
                format!("Contig {contig} not present in the sequence dictionary {dictionary} ")
            }
            RegionError::NullPaddedSpan => "Null object is not allowed here.".to_string(),
            RegionError::PaddedDoesNotContainActive => {
                "Padded span must contain active span.".to_string()
            }
            RegionError::NegativePadding => "the padding size must be 0 or greater".to_string(),
            RegionError::RequestedPaddedDoesNotContain => {
                "The requested padded span must fully contain the requested span".to_string()
            }
            RegionError::NoOverlapToIntersect { left, right } => format!(
                "SimpleInterval::intersect(): The two intervals need to overlap {left} {right}"
            ),
            RegionError::InvalidReadInterval { contig, start, end } => {
                format!("Invalid interval. Contig:{contig} start:{start} end:{end}")
            }
            RegionError::ReadDoesNotOverlap { read_loc, padded } => format!(
                "Read location {read_loc} doesn't overlap with active region padded span {padded}"
            ),
            RegionError::ReadOnDifferentContig { last, read } => format!(
                "Attempting to add a read to ActiveRegion not on the same contig as other reads: \
                 lastRead {last} attempting to add {read}"
            ),
            RegionError::ReadOutOfOrder {
                last,
                last_start,
                read,
                read_start,
            } => format!(
                "Attempting to add a read to ActiveRegion out of order w.r.t. other reads: \
                 lastRead {last} at {last_start} attempting to add {read} at {read_start}"
            ),
            RegionError::Clip(error) => format!("{error:?}"),
        }
    }
}

/// `SimpleInterval.toString()`.
pub fn render_interval(interval: &SimpleInterval) -> String {
    format!("{}:{}-{}", interval.contig, interval.start, interval.end)
}

/// `GATKRead.commonToString()`, which the ordering errors embed.
///
/// An unmapped read, or one with an empty cigar, renders as `<name> UNMAPPED` rather than as a
/// position: the comment upstream says `SAMRecord` blows up on `getAlignmentEnd` with a null cigar.
pub fn render_read(record: &BamRecord, header: &SamHeader) -> String {
    if read::is_unmapped(record) || record.cigar.elements.is_empty() {
        format!("{} UNMAPPED", record.read_name)
    } else {
        format!(
            "{} {}:{}-{}",
            record.read_name,
            contig_name(record, header),
            read_utils::start(record),
            read_utils::end(record)
        )
    }
}

/// `GATKRead.getContig()`: **null** for an unmapped read, whatever reference the record names.
///
/// The null is what reaches `SimpleInterval`'s message, so it is rendered as the literal `null`
/// Java prints rather than as an absent value.
fn contig_name(record: &BamRecord, header: &SamHeader) -> String {
    if read::is_unmapped(record) {
        return "null".to_string();
    }
    usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.clone())
        .unwrap_or_else(|| "null".to_string())
}

fn contig_length(header: &SamHeader, contig: &str) -> Option<i32> {
    header
        .sequences
        .iter()
        .find(|sequence| sequence.name == contig)
        .map(|sequence| sequence.length)
}

/// The sequence dictionary as `Arrays.deepToString` of its names, which is what the missing-contig
/// message embeds.
fn pretty_print_sequence_records(header: &SamHeader) -> String {
    let names: Vec<&str> = header
        .sequences
        .iter()
        .map(|sequence| sequence.name.as_str())
        .collect();
    format!("[{}]", names.join(", "))
}

/// A region of the genome the assembly engine works over.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyRegion {
    active_span: SimpleInterval,
    padded_span: SimpleInterval,
    is_active: bool,
    reads: Vec<BamRecord>,
    hard_clipped_pileup_reads: Vec<BamRecord>,
    has_been_finalized: bool,
}

impl AssemblyRegion {
    /// `new AssemblyRegion(activeSpan, isActive, padding, header)`.
    ///
    /// The padding is applied through `trimIntervalToContig`, so it is clamped to the contig rather
    /// than refused when it runs off an end, and it becomes [`RegionError::NullPaddedSpan`] when
    /// there is no part of the contig left to clamp to.
    pub fn with_padding(
        active_span: SimpleInterval,
        is_active: bool,
        padding: i32,
        header: &SamHeader,
    ) -> Result<AssemblyRegion, RegionError> {
        let Some(length) = contig_length(header, &active_span.contig) else {
            return Err(RegionError::MissingContig {
                contig: active_span.contig.clone(),
                dictionary: pretty_print_sequence_records(header),
            });
        };
        let padded = trim_interval_to_contig(
            &active_span.contig,
            active_span.start - padding,
            active_span.end + padding,
            length,
        );
        match padded {
            None => Err(RegionError::NullPaddedSpan),
            Some(padded) => AssemblyRegion::new(active_span, padded, is_active),
        }
    }

    /// `new AssemblyRegion(activeSpan, paddedSpan, isActive, header)`.
    ///
    /// The zero-size check upstream is unreachable: a `SimpleInterval` cannot be built with
    /// `end < start`, so its size is always at least 1 by the time it gets here.
    pub fn new(
        active_span: SimpleInterval,
        padded_span: SimpleInterval,
        is_active: bool,
    ) -> Result<AssemblyRegion, RegionError> {
        if !padded_span.contains(&active_span) {
            return Err(RegionError::PaddedDoesNotContainActive);
        }
        Ok(AssemblyRegion {
            active_span,
            padded_span,
            is_active,
            reads: Vec::new(),
            hard_clipped_pileup_reads: Vec::new(),
            has_been_finalized: false,
        })
    }

    pub fn span(&self) -> &SimpleInterval {
        &self.active_span
    }

    pub fn padded_span(&self) -> &SimpleInterval {
        &self.padded_span
    }

    pub fn is_active(&self) -> bool {
        self.is_active
    }

    /// Package-private upstream, and deliberately so: only the walker changes a region's activity
    /// after construction.
    pub fn set_is_active(&mut self, value: bool) {
        self.is_active = value;
    }

    pub fn is_finalized(&self) -> bool {
        self.has_been_finalized
    }

    pub fn set_finalized(&mut self, value: bool) {
        self.has_been_finalized = value;
    }

    pub fn reads(&self) -> &[BamRecord] {
        &self.reads
    }

    pub fn hard_clipped_pileup_reads(&self) -> &[BamRecord] {
        &self.hard_clipped_pileup_reads
    }

    /// `AssemblyRegion.size()`: the reads, not the span.
    pub fn size(&self) -> usize {
        self.reads.len()
    }

    /// `AssemblyRegion.clearReads()`, which clears **both** lists.
    pub fn clear_reads(&mut self) {
        self.reads.clear();
        self.hard_clipped_pileup_reads.clear();
    }

    /// `AssemblyRegion.add`.
    pub fn add(&mut self, record: BamRecord, header: &SamHeader) -> Result<(), RegionError> {
        let padded = self.padded_span.clone();
        validate_addition(&self.reads, &record, &padded, header)?;
        self.reads.push(record);
        Ok(())
    }

    /// `AssemblyRegion.addHardClippedPileupReads`, which validates against the same padded span but
    /// keeps its own ordering invariant.
    pub fn add_hard_clipped_pileup_read(
        &mut self,
        record: BamRecord,
        header: &SamHeader,
    ) -> Result<(), RegionError> {
        let padded = self.padded_span.clone();
        validate_addition(&self.hard_clipped_pileup_reads, &record, &padded, header)?;
        self.hard_clipped_pileup_reads.push(record);
        Ok(())
    }

    /// `AssemblyRegion.addAll`, which is `add` in a loop and therefore stops at the first refusal
    /// with the reads before it already in the region.
    pub fn add_all(
        &mut self,
        records: impl IntoIterator<Item = BamRecord>,
        header: &SamHeader,
    ) -> Result<(), RegionError> {
        for record in records {
            self.add(record, header)?;
        }
        Ok(())
    }

    /// `AssemblyRegion.trim(span, padding)`.
    ///
    /// Note what is padded: the **requested** span, not the region. The result is then intersected
    /// with the region's existing padded span, which is the only thing keeping it inside the
    /// original region. This is the entry point whose javadoc describes a different answer from the
    /// one the code produces.
    pub fn trim_with_padding(
        &self,
        span: &SimpleInterval,
        padding: i32,
        header: &SamHeader,
    ) -> Result<AssemblyRegion, RegionError> {
        if padding < 0 {
            return Err(RegionError::NegativePadding);
        }
        let length = contig_length(header, &span.contig).unwrap_or(i32::MAX);
        let Some(padded) = span.expand_within_contig(padding, length) else {
            return Err(RegionError::NullPaddedSpan);
        };
        self.trim(span, &padded, header)
    }

    /// `AssemblyRegion.trim(span, paddedSpan)`.
    ///
    /// The reads are hard-clipped to the new padded span, the empty and non-overlapping ones are
    /// dropped, and the survivors are sorted by `ReadCoordinateComparator`. The sort is the part a
    /// port is most likely to skip, and it is observable: a region trimmed to a span its reads
    /// already fit inside still comes back in comparator order rather than in the order it held.
    pub fn trim(
        &self,
        span: &SimpleInterval,
        padded_span: &SimpleInterval,
        header: &SamHeader,
    ) -> Result<AssemblyRegion, RegionError> {
        if !padded_span.contains(span) {
            return Err(RegionError::RequestedPaddedDoesNotContain);
        }
        let new_active =
            self.active_span
                .intersect(span)
                .ok_or_else(|| RegionError::NoOverlapToIntersect {
                    left: render_interval(&self.active_span),
                    right: render_interval(span),
                })?;
        let new_padded = self.padded_span.intersect(padded_span).ok_or_else(|| {
            RegionError::NoOverlapToIntersect {
                left: render_interval(&self.padded_span),
                right: render_interval(padded_span),
            }
        })?;

        let mut result = AssemblyRegion::new(new_active, new_padded.clone(), self.is_active)?;

        let mut trimmed: Vec<BamRecord> = Vec::new();
        for record in &self.reads {
            let clipped = clipping::hard_clip_to_region(
                record,
                Some(header),
                new_padded.start,
                new_padded.end,
            )
            .map_err(RegionError::Clip)?;
            // `GATKRead.isEmpty()` is "no bases", which is what hard-clipping a read out of the
            // region leaves behind, and the overlap test is redone because clipping moves the
            // start.
            if clipped.read_bases.is_empty() {
                continue;
            }
            if !overlaps(&clipped, &new_padded, header) {
                continue;
            }
            trimmed.push(clipped);
        }
        trimmed.sort_by(read_utils::compare_read_coordinate);

        result.add_all(trimmed, header)?;
        Ok(result)
    }
}

/// `AssemblyRegion.addToReadCollectionAndValidate`.
///
/// The order of the checks is the behaviour: the read's own interval is built **first**, so an
/// unmapped read fails on `SimpleInterval`'s validation and the message talks about an invalid
/// interval rather than about an unmapped read.
fn validate_addition(
    collection: &[BamRecord],
    record: &BamRecord,
    padded_span: &SimpleInterval,
    header: &SamHeader,
) -> Result<(), RegionError> {
    let contig = contig_name(record, header);
    let start = read_utils::start(record);
    let end = read_utils::end(record);
    if read::is_unmapped(record) || start < 1 || end < start {
        return Err(RegionError::InvalidReadInterval { contig, start, end });
    }
    if !padded_span.overlaps(&contig, start, end) {
        return Err(RegionError::ReadDoesNotOverlap {
            read_loc: format!("{contig}:{start}-{end}"),
            padded: render_interval(padded_span),
        });
    }

    if let Some(last) = collection.last() {
        if contig_name(last, header) != contig {
            return Err(RegionError::ReadOnDifferentContig {
                last: render_read(last, header),
                read: render_read(record, header),
            });
        }
        let last_start = read_utils::start(last);
        if start < last_start {
            return Err(RegionError::ReadOutOfOrder {
                last: render_read(last, header),
                last_start,
                read: render_read(record, header),
                read_start: start,
            });
        }
    }
    Ok(())
}

fn overlaps(record: &BamRecord, interval: &SimpleInterval, header: &SamHeader) -> bool {
    interval.overlaps(
        &contig_name(record, header),
        read_utils::start(record),
        read_utils::end(record),
    )
}
