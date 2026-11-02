//! `ReadClassifier`, `BreakpointEvidence`, `KSWindowFinder` and `IntervalCoverageFinder`
//! (GATK 4.6.2.0): what `ExtractSVEvidenceSpark` finds in a partition of reads, as text.
//!
//! # Windows are counted from where the reads are
//!
//! The K-S windows are tiled every 250 bases from the first testable read, not from the contig, and
//! a window is tested when the next testable read leaves it. So the last window is only tested by
//! the sentinel `FlatMapGluer` feeds in at the end of the partition, and a library's anomaly comes
//! out in the order of a `HashMap` of library names.
//!
//! # A read can be evidence several times over
//!
//! The order is the reference's: a soft clip at each end, a big indel per element, the pair's
//! discordance, then whatever the K-S finder had been waiting to say about the previous window.
//!
//! # High coverage is everything's first filter
//!
//! `findGenomewideHighCoverageIntervalsToIgnore` runs before anything is classified: a read wholly
//! inside a region that reached the peak depth is not evidence of anything.
//!
//! # `Integer.MAX_VALUE` is a mapping quality
//!
//! A discordant pair with no `MQ` tag takes its mate's mapping quality to be the largest there is,
//! so its distal target always passes `--min-evidence-mapq`.

use crate::read_metadata::{LibraryName, MetadataArguments, ReadMetadata};
use gatk_engine::java_hash::JavaHashMap;
use gatk_engine::read_utils;
use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// `FindBreakpointEvidenceSpark.DEPTH_WINDOW_SIZE`.
const DEPTH_WINDOW_SIZE: i32 = 100000;
/// `ReadClassifier`'s constants.
const MIN_SOFT_CLIP_LEN: i32 = 30;
const MIN_INDEL_LEN: u32 = 40;
const MIN_QUALITY: u8 = 15;
const MAX_LOW_QUALITY_SCORES: i32 = 3;
const MAX_NON_OUTLIER_ZISH_SCORE: f32 = 6.0;
/// `SplitRead.UNCERTAINTY`, and `LargeIndel`'s.
const SPLIT_UNCERTAINTY: i32 = 3;
const INDEL_UNCERTAINTY: i32 = 4;
/// `KSWindowFinder`'s.
const BLOCK_SIZE: i32 = 250;
const KS_SIGNIFICANCE: f32 = 3.1E-4;
const EVIDENCE_WEIGHT: i32 = 5;

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST: u16 = 0x40;
const SECOND: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const QC_FAIL: u16 = 0x200;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// What the reference throws inside the task, as the root cause's message.
pub type EvidenceError = String;

const UNPAIRED_MATE: &str = "Cannot get mate information for an unpaired read";

/// `FindBreakpointEvidenceSparkArgumentCollection`, the arguments this tool reads.
#[derive(Debug, Clone, Copy)]
pub struct EvidenceArguments {
    pub metadata: MetadataArguments,
    pub high_depth_coverage_peak_factor: i32,
    pub high_depth_coverage_factor: i32,
}

/// `SVInterval`, half open.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SvInterval {
    pub contig: i32,
    pub start: i32,
    pub end: i32,
}

impl SvInterval {
    pub fn overlaps(&self, that: &SvInterval) -> bool {
        self.contig == that.contig && self.start < that.end && that.start < self.end
    }

    pub fn overlap_len(&self, that: &SvInterval) -> i32 {
        if self.contig != that.contig {
            return 0;
        }
        (self.end.min(that.end) - self.start.max(that.start)).max(0)
    }

    pub fn is_upstream_of(&self, that: &SvInterval) -> bool {
        self.contig < that.contig || (self.contig == that.contig && self.end <= that.start)
    }

    pub fn gap_len(&self, that: &SvInterval) -> i32 {
        if self.contig != that.contig {
            return i32::MAX;
        }
        that.start - self.end
    }

    pub fn join(&self, that: &SvInterval) -> SvInterval {
        SvInterval {
            contig: self.contig,
            start: self.start.min(that.start),
            end: self.end.max(that.end),
        }
    }

    pub fn text(&self) -> String {
        format!("{}[{}:{}]", self.contig, self.start, self.end)
    }
}

/// `SVReadFilter.notJunk` and `isMapped`.
pub(crate) fn is_mapped(record: &BamRecord) -> bool {
    record.flags & (DUPLICATE | QC_FAIL | UNMAPPED) == 0
        && record.reference_index >= 0
        && record.alignment_start != 0
}

pub(crate) fn is_primary_line(record: &BamRecord) -> bool {
    record.flags & (SECONDARY | SUPPLEMENTARY) == 0
}

/// `CigarUtils.countAlignedBases`.
fn aligned_bases(record: &BamRecord) -> i64 {
    record
        .cigar
        .elements
        .iter()
        .filter(|element| matches!(element.op, Op::M | Op::Eq | Op::X))
        .map(|element| i64::from(element.length))
        .sum()
}

/// `isEvidence`.
fn is_evidence(record: &BamRecord, arguments: &MetadataArguments) -> bool {
    is_mapped(record)
        && record.mapping_quality >= arguments.min_evidence_mapq
        && aligned_bases(record) >= arguments.min_evidence_match_length
        && record.flags & SECONDARY == 0
}

/// `mateIsUnmapped`, which refuses to answer about an unpaired read.
fn mate_is_unmapped(record: &BamRecord) -> Result<bool, EvidenceError> {
    if record.flags & PAIRED == 0 {
        return Err(UNPAIRED_MATE.to_string());
    }
    Ok(record.flags & MATE_UNMAPPED != 0
        || record.mate_reference_index < 0
        || record.mate_alignment_start == 0)
}

/// `isTemplateLenTestable`.
fn is_template_length_testable(
    record: &BamRecord,
    arguments: &MetadataArguments,
) -> Result<bool, EvidenceError> {
    if !is_evidence(record, arguments) || !is_primary_line(record) {
        return Ok(false);
    }
    Ok(!mate_is_unmapped(record)?
        && record.flags & REVERSE == 0
        && record.flags & MATE_REVERSE != 0
        && record.reference_index == record.mate_reference_index
        && record.alignment_start - arguments.allowed_short_fragment_overhang
            <= record.mate_alignment_start)
}

fn read_group(record: &BamRecord) -> Option<&str> {
    match record.tags.get(Tag::new(b"RG")) {
        Some(TagValue::Str(group)) => Some(group),
        _ => None,
    }
}

fn string_tag<'a>(record: &'a BamRecord, tag: &[u8; 2]) -> Option<&'a str> {
    match record.tags.get(Tag::new(tag)) {
        Some(TagValue::Str(value)) => Some(value),
        _ => None,
    }
}

/// `ReadMetadata.getLibraryName`.
fn library_name(
    metadata: &ReadMetadata,
    group: Option<&str>,
) -> Result<LibraryName, EvidenceError> {
    let Some(group) = group else {
        return Ok(LibraryName(None));
    };
    match metadata.read_group_to_library.get(&group.to_string()) {
        Some(library) => Ok(LibraryName(library.clone())),
        None => Err(format!("No such read group in header: {group}")),
    }
}

/// `getFragmentLengthStatistics(...).getMaxNonOutlierFragmentSize()`, where a library nobody
/// counted is a null to dereference.
fn max_non_outlier_fragment_size(
    metadata: &ReadMetadata,
    group: Option<&str>,
) -> Result<i32, EvidenceError> {
    let library = library_name(metadata, group)?;
    match metadata.libraries.get(&library) {
        Some(statistics) => Ok(statistics.cdf.pop_stat(0.90)),
        None => Err("null".to_string()),
    }
}

/// `Cigar.getPaddedReferenceLength`.
fn padded_reference_length(elements: &[CigarElement]) -> i32 {
    elements
        .iter()
        .filter(|element| matches!(element.op, Op::M | Op::D | Op::N | Op::Eq | Op::X | Op::P))
        .map(|element| element.length as i32)
        .sum()
}

/// `TextCigarCodec.decode`, over what the tags hold.
fn text_cigar(text: &str) -> Result<Vec<CigarElement>, EvidenceError> {
    crate::pathseq_score::text_cigar(text).map_err(|error| match error {
        crate::pathseq_score::ScoreError::BadInput(message)
        | crate::pathseq_score::ScoreError::Other { message, .. } => message,
    })
}

/// `getLeadingMismatches`: the `MD` elements up to the first run of matches, through `TextMDCodec`'s
/// `\G(?:([0-9]+)|([ACTGNactgn])|(\^[ACTGNactgn]+))`, which stops at what it cannot read.
fn leading_mismatches(record: &BamRecord, from_start: bool) -> i32 {
    let Some(md) = string_tag(record, b"MD") else {
        return 0;
    };
    // (length, is a match)
    let mut elements = Vec::new();
    let bytes = md.as_bytes();
    let is_base = |byte: u8| b"ACTGNactgn".contains(&byte);
    let mut at = 0;
    while at < bytes.len() {
        if bytes[at].is_ascii_digit() {
            let from = at;
            while at < bytes.len() && bytes[at].is_ascii_digit() {
                at += 1;
            }
            elements.push((md[from..at].parse::<i32>().unwrap_or(0), true));
        } else if is_base(bytes[at]) {
            elements.push((1, false));
            at += 1;
        } else if bytes[at] == b'^' && bytes.get(at + 1).copied().is_some_and(is_base) {
            let from = at + 1;
            at += 1;
            while at < bytes.len() && is_base(bytes[at]) {
                at += 1;
            }
            elements.push(((at - from) as i32, false));
        } else {
            break;
        }
    }
    let mut count = 0;
    let order: Box<dyn Iterator<Item = &(i32, bool)>> = if from_start {
        Box::new(elements.iter())
    } else {
        Box::new(elements.iter().rev())
    };
    for (length, is_match) in order {
        if *is_match && *length > 0 {
            break;
        }
        count += length;
    }
    count
}

/// `restOfFragmentInterval`.
fn rest_of_fragment(
    record: &BamRecord,
    metadata: &ReadMetadata,
) -> Result<SvInterval, EvidenceError> {
    let template = max_non_outlier_fragment_size(metadata, read_group(record))?;
    let (start, width) = if record.flags & REVERSE != 0 {
        let read_start = record.alignment_start + leading_mismatches(record, true);
        let mut width = read_start - (read_utils::unclipped_end(record) + 1 - template);
        let mut start = read_start - width;
        if start < 1 {
            width += start - 1;
            start = 1;
        }
        (start, width)
    } else {
        let read_end = read_utils::end(record) + 1 - leading_mismatches(record, false);
        (
            read_end,
            read_utils::unclipped_start(record) + template - read_end,
        )
    };
    Ok(SvInterval {
        contig: record.reference_index,
        start,
        end: start + width,
    })
}

/// `fixedWidthInterval`.
pub fn fixed_width(contig: i32, offset: i32, uncertainty: i32) -> SvInterval {
    let mut width = 2 * uncertainty;
    let mut start = offset - uncertainty;
    if start < 1 {
        width += start - 1;
        start = 1;
    }
    SvInterval {
        contig,
        start,
        end: start + width,
    }
}

/// `TemplateFragmentOrdinal.forRead`.
pub(crate) fn fragment_ordinal(record: &BamRecord) -> &'static str {
    if record.flags & PAIRED == 0 {
        return "";
    }
    match (record.flags & FIRST != 0, record.flags & SECOND != 0) {
        (true, true) => "/0",
        (true, false) => "/1",
        (false, true) => "/2",
        (false, false) => "/?",
    }
}

/// `SplitRead`'s `SA` tag, parsed by the constructor whether or not it is looked at.
struct SaMapping {
    contig: String,
    start: i32,
    forward: bool,
    cigar: String,
    mapq: i32,
}

/// `parseSATag`.
fn parse_sa_tag(tag: &str) -> Result<Vec<SaMapping>, EvidenceError> {
    // `String.split` drops the empty strings at the end.
    let mut pieces: Vec<&str> = tag.split(';').collect();
    while pieces.last().is_some_and(|piece| piece.is_empty()) {
        pieces.pop();
    }
    let number = |text: &str| -> Result<i32, EvidenceError> {
        text.parse::<i32>()
            .map_err(|_| format!("For input string: \"{text}\""))
    };
    let mut mappings = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let fields: Vec<&str> = piece.split(',').collect();
        if fields.len() != 6 {
            return Err(format!("Could not parse SATag: {piece}"));
        }
        let start = number(fields[1])?;
        let mapq = number(fields[4])?;
        number(fields[5])?;
        mappings.push(SaMapping {
            contig: fields[0].to_string(),
            start,
            forward: fields[2] == "+",
            cigar: fields[3].to_string(),
            mapq,
        });
    }
    Ok(mappings)
}

/// `getContigID`.
fn contig_id(metadata: &ReadMetadata, name: &str) -> Result<i32, EvidenceError> {
    metadata
        .contig_ids
        .get(&name.to_string())
        .copied()
        .ok_or_else(|| format!("No such contig name: {name}"))
}

/// `SplitRead`'s distal targets: one supplementary alignment or none, that is good and elsewhere.
fn split_read_targets(
    record: &BamRecord,
    location: &SvInterval,
    clipped_at_start: bool,
    metadata: &ReadMetadata,
    min_mapq: i32,
) -> Result<Option<(SvInterval, bool)>, EvidenceError> {
    let Some(tag) = string_tag(record, b"SA") else {
        return Ok(None);
    };
    let mappings = parse_sa_tag(tag)?;
    if mappings.len() != 1 {
        return Ok(None);
    }
    let mapping = &mappings[0];
    let primary_forward = record.flags & REVERSE == 0;
    // `calculateDistalTargetStrand(saMapping, !clippedAtStart, primaryForward)`.
    let strand = if primary_forward == !clipped_at_start {
        !mapping.forward
    } else {
        mapping.forward
    };
    let contig = contig_id(metadata, &mapping.contig)?;
    let cigar = text_cigar(&mapping.cigar)?;
    let position = if strand {
        mapping.start + padded_reference_length(&cigar)
    } else {
        mapping.start
    };
    let interval = SvInterval {
        contig,
        start: (position - SPLIT_UNCERTAINTY).max(1),
        end: position + SPLIT_UNCERTAINTY + 1,
    };
    let good = mapping.mapq >= min_mapq
        && (interval.contig == location.contig
            || (!metadata.ignore_cross_contig(interval.contig)
                && !metadata.ignore_cross_contig(location.contig)))
        && !interval.overlaps(location);
    if !good {
        return Ok(None);
    }
    Ok(Some((interval, strand)))
}

/// `DiscordantReadPairEvidence`'s target, and that target as a distal one if it is.
fn discordant_pair_targets(
    record: &BamRecord,
    location: &SvInterval,
    metadata: &ReadMetadata,
    min_mapq: i32,
) -> Result<(SvInterval, Option<(SvInterval, bool)>), EvidenceError> {
    let mate_reverse = record.flags & MATE_REVERSE != 0;
    let max_fragment = max_non_outlier_fragment_size(metadata, read_group(record))?;
    let mate_length = match string_tag(record, b"MC") {
        Some(cigar) => padded_reference_length(&text_cigar(cigar)?),
        None => record.read_bases.len() as i32,
    };
    let mate_start = record.mate_alignment_start;
    let target = SvInterval {
        contig: record.mate_reference_index,
        start: if mate_reverse {
            mate_start - max_fragment + mate_length
        } else {
            mate_start + mate_length - 2
        }
        .max(0),
        end: if mate_reverse {
            mate_start + 2
        } else {
            mate_start + max_fragment
        },
    };
    let quality = match record.tags.get(Tag::new(b"MQ")) {
        Some(TagValue::Int(quality)) => *quality as i32,
        Some(TagValue::Str(text)) => text
            .parse::<i32>()
            .map_err(|_| format!("Attribute MQ not of integer type: {text}"))?,
        _ => i32::MAX,
    };
    let distal = location.is_upstream_of(&target)
        && quality >= min_mapq
        && !target.overlaps(location)
        && !metadata.ignore_cross_contig(target.contig);
    Ok((target, distal.then_some((target, !mate_reverse))))
}

/// Which subclass of `ReadEvidence`, with what its `toString` adds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadKind {
    SplitRead { sa_tag: Option<String> },
    LargeIndel,
    MateUnmapped,
    InterContigPair { target: SvInterval },
    SameStrandPair { target: SvInterval },
    OutiesPair { target: SvInterval },
    WeirdTemplateSize { target: SvInterval },
}

/// What `ReadEvidence` keeps of its read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadEvidence {
    pub kind: ReadKind,
    pub template_name: String,
    pub ordinal: &'static str,
    /// `evidenceUpstreamOfBreakpoint`.
    pub forward: bool,
    pub template_size: i32,
    pub cigar: String,
    pub mapping_quality: u8,
}

/// `BreakpointEvidence`'s class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvidenceKind {
    /// `BreakpointEvidence` itself, what a cluster is replaced with.
    Cluster,
    External,
    TemplateSizeAnomaly {
        read_count: i32,
    },
    Read(Box<ReadEvidence>),
}

/// `BreakpointEvidence`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Evidence {
    pub location: SvInterval,
    pub weight: i32,
    pub validated: bool,
    pub kind: EvidenceKind,
    /// `getDistalTargets(...)` where `hasDistalTargets`: the reference never finds more than one.
    pub distal: Option<(SvInterval, bool)>,
}

impl Evidence {
    /// `stringRep`.
    pub fn string_rep(&self) -> String {
        let class = match &self.kind {
            EvidenceKind::Cluster => "BreakpointEvidence",
            EvidenceKind::External => "ExternalEvidence",
            EvidenceKind::TemplateSizeAnomaly { .. } => "TemplateSizeAnomaly",
            EvidenceKind::Read(read) => match read.kind {
                ReadKind::SplitRead { .. } => "SplitRead",
                ReadKind::LargeIndel => "LargeIndel",
                ReadKind::MateUnmapped => "MateUnmapped",
                ReadKind::InterContigPair { .. } => "InterContigPair",
                ReadKind::SameStrandPair { .. } => "SameStrandPair",
                ReadKind::OutiesPair { .. } => "OutiesPair",
                ReadKind::WeirdTemplateSize { .. } => "WeirdTemplateSize",
            },
        };
        let distal = match &self.distal {
            Some((interval, strand)) => {
                format!("{}{}", interval.text(), if *strand { "1" } else { "0" })
            }
            None => String::new(),
        };
        let base = format!(
            "{}\t{}\t{class}\t{distal}",
            self.location.text(),
            self.weight
        );
        match &self.kind {
            EvidenceKind::Cluster | EvidenceKind::External => base,
            EvidenceKind::TemplateSizeAnomaly { read_count } => format!("{base}\t{read_count}"),
            EvidenceKind::Read(read) => format!(
                "{base}\t{}{}\t{}\t{}\t{}\t{}",
                read.template_name,
                read.ordinal,
                if read.forward { "1" } else { "0" },
                read.template_size,
                read.cigar,
                read.mapping_quality
            ),
        }
    }

    /// `toString`, which is what `--breakpoint-evidence-dir` is written with.
    pub fn to_text(&self) -> String {
        let base = format!("{}^{}", self.location.text(), self.weight);
        match &self.kind {
            EvidenceKind::Cluster => base,
            EvidenceKind::External => format!("{base}\tExternalEvidence"),
            EvidenceKind::TemplateSizeAnomaly { read_count } => {
                format!("{base}\tTemplateSizeAnomaly\t{read_count}")
            }
            EvidenceKind::Read(read) => {
                let base = format!("{base}\t{}{}", read.template_name, read.ordinal);
                match &read.kind {
                    ReadKind::SplitRead { sa_tag } => format!(
                        "{base}\tSplit\t{}\t{}",
                        read.cigar,
                        match sa_tag {
                            Some(tag) => format!(" SA: {tag}"),
                            None => " SA: None".to_string(),
                        }
                    ),
                    ReadKind::LargeIndel => format!("{base}\tIndel\t{}", read.cigar),
                    ReadKind::MateUnmapped => format!("{base}\tUnmappedMate"),
                    ReadKind::InterContigPair { target } => {
                        format!("{base}\tIntercontigPair\t{}", target.text())
                    }
                    ReadKind::OutiesPair { target } => {
                        format!("{base}\tOutiesPair\t{}", target.text())
                    }
                    ReadKind::SameStrandPair { target } => {
                        format!("{base}\tSameStrandPair\t{}", target.text())
                    }
                    ReadKind::WeirdTemplateSize { target } => format!(
                        "{base}\tTemplateSize\t{}\t{}",
                        target.text(),
                        read.template_size
                    ),
                }
            }
        }
    }
}

/// `ReadEvidence`'s constructor.
fn read_evidence(
    record: &BamRecord,
    location: SvInterval,
    kind: ReadKind,
    distal: Option<(SvInterval, bool)>,
    forward: bool,
) -> Evidence {
    Evidence {
        location,
        weight: 1,
        validated: false,
        kind: EvidenceKind::Read(Box::new(ReadEvidence {
            kind,
            template_name: record.read_name.clone(),
            ordinal: fragment_ordinal(record),
            forward,
            template_size: record.inferred_insert_size,
            cigar: record.cigar.to_text(),
            mapping_quality: record.mapping_quality,
        })),
        distal,
    }
}

/// `hasInitialSoftClip`.
fn has_initial_soft_clip(record: &BamRecord) -> bool {
    let elements = &record.cigar.elements;
    let Some(mut first) = elements.first() else {
        return false;
    };
    if first.op == Op::H && elements.len() > 1 {
        first = &elements[1];
    }
    let clip_start = first.length as i32 - MIN_SOFT_CLIP_LEN;
    first.op == Op::S && clip_start >= 0 && is_high_quality(&record.base_qualities, clip_start)
}

/// `hasFinalSoftClip`.
fn has_final_soft_clip(record: &BamRecord) -> bool {
    let elements = &record.cigar.elements;
    let Some(mut last) = elements.last() else {
        return false;
    };
    if last.op == Op::H && elements.len() > 1 {
        last = &elements[elements.len() - 2];
    }
    last.op == Op::S
        && last.length as i32 >= MIN_SOFT_CLIP_LEN
        && is_high_quality(
            &record.base_qualities,
            record.read_bases.len() as i32 - last.length as i32,
        )
}

/// `isHighQualityRegion`.
fn is_high_quality(qualities: &[u8], start: i32) -> bool {
    let mut low = 0;
    for at in start..start + MIN_SOFT_CLIP_LEN {
        if qualities[at as usize] < MIN_QUALITY {
            low += 1;
            if low > MAX_LOW_QUALITY_SCORES {
                return false;
            }
        }
    }
    true
}

/// `IntHistogram`, as a pair of `long[]` the finder toggles between.
struct HistogramPair {
    counts: [Vec<i64>; 2],
    totals: [i64; 2],
}

/// `KSWindowFinder`.
struct KsWindowFinder<'a> {
    metadata: &'a ReadMetadata,
    arguments: MetadataArguments,
    libraries: JavaHashMap<LibraryName, HistogramPair>,
    fill: usize,
    contig: Option<i32>,
    end: i32,
}

impl<'a> KsWindowFinder<'a> {
    fn new(metadata: &'a ReadMetadata, arguments: MetadataArguments) -> KsWindowFinder<'a> {
        let capacity = (metadata.libraries.len() as u64 * 4 / 3) as usize + 1;
        let mut libraries = JavaHashMap::with_capacity(capacity);
        for (name, statistics) in metadata.libraries.iter() {
            let size = statistics.cdf.fractions.len();
            libraries.insert(
                name.clone(),
                HistogramPair {
                    counts: [vec![0; size], vec![0; size]],
                    totals: [0, 0],
                },
            );
        }
        KsWindowFinder {
            metadata,
            arguments,
            libraries,
            fill: 0,
            contig: None,
            end: 0,
        }
    }

    /// `testReadAndGatherEvidence`.
    fn test_read(
        &mut self,
        record: &BamRecord,
        lines: &mut Vec<Evidence>,
    ) -> Result<(), EvidenceError> {
        if !is_template_length_testable(record, &self.arguments)? {
            return Ok(());
        }
        let same_block =
            self.contig == Some(record.reference_index) && record.alignment_start < self.end;
        if !same_block {
            self.check_histograms(lines);
            self.advance_block(record);
        }
        // `addObservation`.
        let library = library_name(self.metadata, read_group(record))?;
        let fill = self.fill;
        let Some(pair) = self.libraries.get_mut(&library) else {
            return Err("null".to_string());
        };
        let counts = &mut pair.counts[fill];
        let size = record.inferred_insert_size.unsigned_abs() as usize;
        let last = counts.len() - 1;
        counts[size.min(last)] += 1;
        pair.totals[fill] += 1;
        Ok(())
    }

    /// `checkHistograms`.
    fn check_histograms(&mut self, lines: &mut Vec<Evidence>) {
        let Some(contig) = self.contig else {
            return;
        };
        let old = self.fill ^ 1;
        let start = (self.end - 2 * BLOCK_SIZE).max(1);
        let interval = SvInterval {
            contig,
            start,
            end: self.end,
        };
        let names: Vec<LibraryName> = self.libraries.keys().cloned().collect();
        for name in names {
            let fill = self.fill;
            let pair = self.libraries.get_mut(&name).expect("a key just listed");
            let [first, second] = &mut pair.counts;
            let (old_counts, fill_counts) = if old == 0 {
                (first, &*second)
            } else {
                (second, &*first)
            };
            for (old_count, fill_count) in old_counts.iter_mut().zip(fill_counts) {
                *old_count += fill_count;
            }
            pair.totals[old] += pair.totals[fill];
            let observed = pair.totals[old];
            if observed > 0 {
                let cdf = &self
                    .metadata
                    .libraries
                    .get(&name)
                    .expect("the same keys")
                    .cdf;
                if is_different_by_ks(&cdf.fractions, cdf.counts, &pair.counts[old], observed) {
                    lines.push(Evidence {
                        location: interval,
                        weight: EVIDENCE_WEIGHT,
                        validated: false,
                        kind: EvidenceKind::TemplateSizeAnomaly {
                            read_count: observed as i32,
                        },
                        distal: None,
                    });
                }
            }
            pair.counts[old].iter_mut().for_each(|count| *count = 0);
            pair.totals[old] = 0;
        }
    }

    /// `advanceBlock`.
    fn advance_block(&mut self, record: &BamRecord) {
        let old = self.fill;
        self.fill ^= 1;
        self.end += BLOCK_SIZE;
        if self.contig != Some(record.reference_index) || record.alignment_start >= self.end {
            self.contig = Some(record.reference_index);
            self.end = record.alignment_start + BLOCK_SIZE;
            let names: Vec<LibraryName> = self.libraries.keys().cloned().collect();
            for name in names {
                let pair = self.libraries.get_mut(&name).expect("a key just listed");
                pair.counts[old].iter_mut().for_each(|count| *count = 0);
                pair.totals[old] = 0;
            }
        }
    }
}

/// `IntHistogram.CDF.isDifferentByKSStatistic`, at `KS_SIGNIFICANCE`.
fn is_different_by_ks(fractions: &[f32], n_counts: i64, sample: &[i64], m_counts: i64) -> bool {
    let c_of_alpha = (-0.5 * jmath::math::log(f64::from(KS_SIGNIFICANCE) / 2.0)).sqrt();
    let significant = (c_of_alpha
        * ((n_counts + m_counts) as f64 / n_counts as f64 / m_counts as f64).sqrt())
        as f32;
    let mut sum = 0i64;
    let mut index = 0;
    let mut current = 0.0f32;
    let mut previous_zero = false;
    let tracked_total = m_counts - sample[sample.len() - 1];
    while sum < tracked_total {
        let value = sample[index];
        if value == 0 {
            previous_zero = true;
        } else {
            if previous_zero && (fractions[index - 1] - current).abs() >= significant {
                return true;
            }
            previous_zero = false;
            sum += value;
            current = sum as f32 / m_counts as f32;
            if (fractions[index] - current).abs() >= significant {
                return true;
            }
        }
        index += 1;
    }
    false
}

/// `findGenomewideHighCoverageIntervalsToIgnore`, over one partition.
pub fn high_coverage_intervals(
    records: &[BamRecord],
    metadata: &ReadMetadata,
    contig_lengths: &[i32],
    arguments: &EvidenceArguments,
) -> Vec<SvInterval> {
    let mut windows = Vec::new();
    for (contig, length) in contig_lengths.iter().enumerate() {
        let mut start = 1;
        while start < *length {
            windows.push(SvInterval {
                contig: contig as i32,
                start,
                end: (*length).min(start + DEPTH_WINDOW_SIZE),
            });
            start += DEPTH_WINDOW_SIZE;
        }
    }
    let min_flanking = (arguments.high_depth_coverage_factor as f32 * metadata.coverage) as i32;
    let min_peak = (arguments.high_depth_coverage_peak_factor as f32 * metadata.coverage) as i32;

    // `IntervalCoverageFinder`.
    let mut coverage: Vec<Option<Vec<i32>>> = vec![None; windows.len()];
    let mut index = 0;
    for record in records {
        if !is_mapped(record) || metadata.ignore_cross_contig(record.reference_index) {
            continue;
        }
        let contig = record.reference_index;
        let start = record.alignment_start;
        while index < windows.len() {
            let window = &windows[index];
            if window.contig > contig || (window.contig == contig && window.end > start) {
                break;
            }
            index += 1;
        }
        if index >= windows.len() {
            break;
        }
        let read = SvInterval {
            contig,
            start,
            end: read_utils::end(record) + 1,
        };
        let mut at = index;
        while at < windows.len() && windows[at].overlaps(&read) {
            let window = windows[at];
            let counts =
                coverage[at].get_or_insert_with(|| vec![0; (window.end - window.start) as usize]);
            for position in read.start.max(window.start)..read.end.min(window.end) {
                counts[(position - window.start) as usize] += 1;
            }
            at += 1;
        }
    }

    // `getHighCoverageIntervalsInWindow`.
    let mut candidates: Vec<(SvInterval, bool)> = Vec::new();
    for (window, counts) in windows.iter().zip(&coverage) {
        let Some(counts) = counts else {
            continue;
        };
        let mut in_region = false;
        let mut region_start = -1;
        let mut found_high = false;
        for (at, count) in counts.iter().enumerate() {
            if *count >= min_flanking {
                if !in_region {
                    in_region = true;
                    found_high = false;
                    region_start = at as i32;
                }
                if *count >= min_peak {
                    found_high = true;
                }
            } else {
                if in_region && (found_high || region_start == 0) {
                    candidates.push((
                        SvInterval {
                            contig: window.contig,
                            start: window.start + region_start,
                            end: window.start + at as i32,
                        },
                        found_high,
                    ));
                }
                in_region = false;
            }
        }
        if in_region {
            candidates.push((
                SvInterval {
                    contig: window.contig,
                    start: window.start + region_start,
                    end: window.end,
                },
                found_high,
            ));
        }
    }
    candidates.sort_by_key(|(interval, _)| *interval);

    // Merge what a window boundary cut, then what a read could span.
    let mut bounded = Vec::new();
    let mut candidates = candidates.into_iter();
    if let Some(mut previous) = candidates.next() {
        for current in candidates {
            if previous.0.gap_len(&current.0) == 0 {
                previous = (previous.0.join(&current.0), previous.1 || current.1);
            } else {
                if previous.1 {
                    bounded.push(previous.0);
                }
                previous = current;
            }
        }
        if previous.1 {
            bounded.push(previous.0);
        }
    }
    let mut merged = Vec::new();
    let mut bounded = bounded.into_iter();
    if let Some(mut previous) = bounded.next() {
        for current in bounded {
            if previous.gap_len(&current) <= metadata.avg_read_len {
                previous = previous.join(&current);
            } else {
                merged.push(previous);
                previous = current;
            }
        }
        merged.push(previous);
    }
    merged
}

/// `--high-coverage-intervals`: a BED line per interval.
pub fn high_coverage_bed(intervals: &[SvInterval], metadata: &ReadMetadata) -> String {
    let mut text = String::new();
    for interval in intervals {
        text.push_str(&format!(
            "{}\t{}\t{}\n",
            metadata.contigs[interval.contig as usize],
            interval.start - 1,
            interval.end
        ));
    }
    text
}

/// `ReadClassifier` over one partition and its sentinel, each piece of evidence through `stringRep`.
pub fn classify(
    records: &[BamRecord],
    metadata: &ReadMetadata,
    ignored: &[SvInterval],
    arguments: &MetadataArguments,
) -> Result<Vec<String>, EvidenceError> {
    Ok(evidence(records, metadata, ignored, arguments)?
        .iter()
        .map(Evidence::string_rep)
        .collect())
}

fn inter_contig(target: SvInterval) -> ReadKind {
    ReadKind::InterContigPair { target }
}

fn same_strand(target: SvInterval) -> ReadKind {
    ReadKind::SameStrandPair { target }
}

fn outies(target: SvInterval) -> ReadKind {
    ReadKind::OutiesPair { target }
}

fn weird_template(target: SvInterval) -> ReadKind {
    ReadKind::WeirdTemplateSize { target }
}

/// `SVReadFilter.containedInRegionToIgnore`.
pub(crate) fn contained_in_region_to_ignore(interval: &SvInterval, ignored: &[SvInterval]) -> bool {
    let length = interval.end - interval.start;
    ignored
        .iter()
        .any(|region| region.overlaps(interval) && region.overlap_len(interval) == length)
}

/// `ReadClassifier` over one partition and its sentinel.
pub fn evidence(
    records: &[BamRecord],
    metadata: &ReadMetadata,
    ignored: &[SvInterval],
    arguments: &MetadataArguments,
) -> Result<Vec<Evidence>, EvidenceError> {
    let mut finder = KsWindowFinder::new(metadata, *arguments);
    let mut lines = Vec::new();
    let min_mapq = i32::from(arguments.min_evidence_mapq);
    for record in records {
        if !is_mapped(record) || metadata.ignore_cross_contig(record.reference_index) {
            continue;
        }
        if !is_evidence(record, arguments) {
            continue;
        }
        let span = SvInterval {
            contig: record.reference_index,
            start: record.alignment_start,
            end: read_utils::end(record),
        };
        if contained_in_region_to_ignore(&span, ignored) {
            continue;
        }
        let contig = record.reference_index;

        // `checkForSplitRead`.
        for at_start in [true, false] {
            let clipped = if at_start {
                has_initial_soft_clip(record)
            } else {
                has_final_soft_clip(record)
            };
            if clipped {
                let offset = if at_start {
                    record.alignment_start
                } else {
                    read_utils::end(record)
                };
                let location = fixed_width(contig, offset, SPLIT_UNCERTAINTY);
                let distal = split_read_targets(record, &location, at_start, metadata, min_mapq)?;
                let sa_tag = string_tag(record, b"SA").map(str::to_string);
                lines.push(read_evidence(
                    record,
                    location,
                    ReadKind::SplitRead { sa_tag },
                    distal,
                    !at_start,
                ));
            }
        }
        // `checkBigIndel`.
        let mut locus = record.alignment_start;
        for element in &record.cigar.elements {
            if element.length >= MIN_INDEL_LEN {
                let offset = match element.op {
                    Op::I => Some(locus),
                    Op::D => Some(locus + element.length as i32 / 2),
                    _ => None,
                };
                if let Some(offset) = offset {
                    let location = fixed_width(contig, offset, INDEL_UNCERTAINTY);
                    lines.push(read_evidence(
                        record,
                        location,
                        ReadKind::LargeIndel,
                        None,
                        true,
                    ));
                }
            }
            if element.op.consumes_reference_bases() {
                locus += element.length as i32;
            }
        }
        // `checkDiscordantPair`.
        if is_primary_line(record) {
            let forward = record.flags & REVERSE == 0;
            if mate_is_unmapped(record)? {
                let location = rest_of_fragment(record, metadata)?;
                lines.push(read_evidence(
                    record,
                    location,
                    ReadKind::MateUnmapped,
                    None,
                    forward,
                ));
            } else {
                let mate_contig = record.mate_reference_index;
                let kind = if contig != mate_contig {
                    (!metadata.ignore_cross_contig(contig)
                        && !metadata.ignore_cross_contig(mate_contig))
                    .then_some(inter_contig as fn(SvInterval) -> ReadKind)
                } else if forward == (record.flags & MATE_REVERSE == 0) {
                    Some(same_strand as fn(SvInterval) -> ReadKind)
                } else if if forward {
                    record.alignment_start - arguments.allowed_short_fragment_overhang
                        > record.mate_alignment_start
                } else {
                    record.alignment_start + arguments.allowed_short_fragment_overhang
                        < record.mate_alignment_start
                } {
                    Some(outies as fn(SvInterval) -> ReadKind)
                } else {
                    let library = library_name(metadata, read_group(record))?;
                    let statistics = metadata
                        .libraries
                        .get(&library)
                        .ok_or_else(|| "null".to_string())?;
                    (statistics.zish_score(record.inferred_insert_size.abs())
                        > MAX_NON_OUTLIER_ZISH_SCORE)
                        .then_some(weird_template as fn(SvInterval) -> ReadKind)
                };
                if let Some(kind) = kind {
                    let location = rest_of_fragment(record, metadata)?;
                    let (target, distal) =
                        discordant_pair_targets(record, &location, metadata, min_mapq)?;
                    lines.push(read_evidence(
                        record,
                        location,
                        kind(target),
                        distal,
                        forward,
                    ));
                }
            }
        }
        finder.test_read(record, &mut lines)?;
    }
    finder.check_histograms(&mut lines);
    Ok(lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_interval_prints_as_the_reference_does() {
        let interval = SvInterval {
            contig: 2,
            start: 10,
            end: 20,
        };
        assert_eq!(interval.text(), "2[10:20]");
        assert_eq!(fixed_width(0, 2, 3).text(), "0[1:5]");
    }

    #[test]
    fn the_ks_test_finds_a_shifted_sample() {
        let mut fractions = vec![0.0f32; 12];
        for (at, fraction) in fractions.iter_mut().enumerate() {
            *fraction = (at as f32 / 10.0).min(1.0);
        }
        let mut sample = vec![0i64; 12];
        sample[9] = 200;
        assert!(is_different_by_ks(&fractions, 10000, &sample, 200));
        let even = vec![10i64, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0];
        assert!(!is_different_by_ks(&fractions, 10000, &even, 100));
    }

    #[test]
    fn a_soft_clip_needs_quality() {
        assert!(is_high_quality(&[30; 40], 0));
        let mut bad = [30u8; 40];
        bad[..4].iter_mut().for_each(|q| *q = 2);
        assert!(!is_high_quality(&bad, 0));
    }
}
