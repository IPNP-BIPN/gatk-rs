//! `DiscoverVariantsFromContigAlignmentsSAMSpark`: assembled contigs' alignments read as chimeras,
//! and the chimeras as simple structural variants.
//!
//! The classes this is ported from are `AlignedContig`, `AlignmentInterval`,
//! `ContigAlignmentsModifier.GapSplitter`, `ContigChimericAlignmentIterativeInterpreter`,
//! `SimpleChimera`, `NovelAdjacencyAndAltHaplotype`, `BreakpointsInference`,
//! `BreakpointComplications`, `SimpleNovelAdjacencyAndChimericAlignmentEvidence`, `SimpleSVType`
//! and `SVVCFWriter`'s sort.
//!
//! Three things in them are not what they look like, and each is visible in the VCF:
//!
//! * **a contig's chimeras share one list of insertion mappings.** `parseOneContig` hands the same
//!   `ArrayList` to each `SimpleChimera` and keeps adding to it, so a chimera built before an
//!   inserted alignment was seen still reports it;
//! * **only the chimeras that are complete and not translocations reach the rest**, so the inverted
//!   duplication, the reference-order swap and the inter-chromosome paths of the classes above
//!   cannot be taken from this tool, only from [`assembly`]'s;
//! * **the order among records that sort equal is `AppendOnlyMap`'s**, the slot being the murmur3 of
//!   `NovelAdjacencyAndAltHaplotype.hashCode`.
//!
//! Everything thrown between the reads and the collected variants is thrown inside a Spark task.
//! The runner reports those, so a failure here is its message.

use gatk_engine::java_hash::{byte_array_hash_code, string_hash_code};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Value, VariantContext};

pub mod assembly;

/// `DiscoverVariantsFromContigAlignmentsArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY`.
const GAP_SENSITIVITY: i32 = 50;
/// `CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD`, which the tool passes rather than `--min-mq`.
const HIGH_MQ: i32 = 60;
/// `DEFAULT_MIN_ALIGNMENT_LENGTH`, which the tool passes rather than `--min-align-length`.
const MIN_ALIGNMENT_LENGTH: i32 = 50;
/// `AlignmentInterval.NO_NM` and `NO_AS`.
const NO_VALUE: i32 = -1;

/// A failure inside a task, as the message of what was thrown.
pub type Failed = String;

/// `SimpleInterval`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Interval {
    /// The constructor, which refuses what `SimpleInterval.isValid` does.
    pub fn new(contig: &str, start: i32, end: i32) -> Result<Self, Failed> {
        if start > 0 && end >= start {
            Ok(Self {
                contig: contig.to_string(),
                start,
                end,
            })
        } else {
            Err(format!(
                "Invalid interval. Contig:{contig} start:{start} end:{end}"
            ))
        }
    }

    fn size(&self) -> i32 {
        self.end - self.start + 1
    }

    fn contains(&self, other: &Interval) -> bool {
        self.contig == other.contig && self.start <= other.start && self.end >= other.end
    }

    fn overlaps(&self, other: &Interval) -> bool {
        self.contig == other.contig && self.start <= other.end && other.start <= self.end
    }

    fn hash_code(&self) -> i32 {
        self.start
            .wrapping_mul(31)
            .wrapping_add(self.end)
            .wrapping_mul(31)
            .wrapping_add(string_hash_code(&self.contig))
    }

    fn text(&self) -> String {
        format!("{}:{}-{}", self.contig, self.start, self.end)
    }
}

fn contig_index(header: &SamHeader, contig: &str) -> i32 {
    header
        .sequences
        .iter()
        .position(|sequence| sequence.name == contig)
        .map_or(-1, |index| index as i32)
}

/// `IntervalUtils.compareLocatables`.
fn compare_locatables(
    header: &SamHeader,
    one: (&str, i32, i32),
    two: (&str, i32, i32),
) -> std::cmp::Ordering {
    contig_index(header, one.0)
        .cmp(&contig_index(header, two.0))
        .then(one.1.cmp(&two.1))
        .then(one.2.cmp(&two.2))
}

/// `AlignmentInterval`.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Alignment {
    span: Interval,
    start_on_contig: i32,
    end_on_contig: i32,
    cigar: Vec<CigarElement>,
    forward: bool,
    map_quality: i32,
    mismatches: i32,
    score: i32,
    /// `alnModType`, as `toString` gives it.
    modification: char,
}

fn is_clipping(op: Op) -> bool {
    matches!(op, Op::S | Op::H)
}

fn is_alignment(op: Op) -> bool {
    matches!(op, Op::M | Op::Eq | Op::X)
}

fn is_indel(op: Op) -> bool {
    matches!(op, Op::I | Op::D)
}

fn consumes_read(op: Op) -> bool {
    matches!(op, Op::M | Op::I | Op::S | Op::Eq | Op::X)
}

fn consumes_reference(op: Op) -> bool {
    matches!(op, Op::M | Op::D | Op::N | Op::Eq | Op::X)
}

/// `TextCigarCodec.encode`.
fn cigar_text(elements: &[CigarElement]) -> String {
    if elements.is_empty() {
        return "*".to_string();
    }
    Cigar {
        elements: elements.to_vec(),
    }
    .to_text()
}

fn read_length(elements: &[CigarElement]) -> i32 {
    elements
        .iter()
        .filter(|element| consumes_read(element.op))
        .map(|element| element.length as i32)
        .sum()
}

fn reference_length(elements: &[CigarElement]) -> i32 {
    elements
        .iter()
        .filter(|element| consumes_reference(element.op))
        .map(|element| element.length as i32)
        .sum()
}

/// `CigarUtils.countClippedBases(cigar, tail)`, both clip types of one end.
fn clipped_bases(elements: &[CigarElement], from_left: bool) -> Result<i32, Failed> {
    if elements.len() < 2 {
        return match elements.first() {
            Some(only) if !is_clipping(only.op) => Ok(0),
            _ => Err("cigar is empty or completely clipped.".to_string()),
        };
    }
    let walk: Box<dyn Iterator<Item = &CigarElement>> = if from_left {
        Box::new(elements.iter())
    } else {
        Box::new(elements.iter().rev())
    };
    let mut clipped = 0;
    for element in walk {
        if !is_clipping(element.op) {
            return Ok(clipped);
        }
        clipped += element.length as i32;
    }
    Err(format!(
        "Input cigar {} is completely clipped.",
        cigar_text(elements)
    ))
}

/// `CigarUtils.computeAssociatedDistOnRead`.
fn associated_distance_on_read(
    cigar: &[CigarElement],
    start: i32,
    reference_distance: i32,
    backward: bool,
) -> Result<i32, Failed> {
    if reference_distance <= 0 || start <= 0 {
        return Err(format!(
            "start {start} or distance {reference_distance} is non-positive."
        ));
    }
    let mut elements = cigar.to_vec();
    if backward {
        elements.reverse();
    }
    let length = read_length(&elements);
    let to_skip = if backward { length - start } else { start - 1 };
    let mut read_consumed = 0;
    let mut reference_consumed = 0;
    for element in &elements {
        let before = read_consumed;
        if consumes_read(element.op) {
            read_consumed += element.length as i32;
        }
        if read_consumed <= to_skip {
            continue;
        }
        if consumes_reference(element.op) {
            reference_consumed += element.length as i32 - (to_skip - before).max(0);
        }
        if reference_consumed >= reference_distance {
            let excess = (reference_consumed - reference_distance).max(0);
            return Ok(read_consumed
                - to_skip
                - if consumes_read(element.op) { excess } else { 0 });
        }
    }
    Err(format!(
        "Cigar {} does not contain at least {reference_distance} reference bases past red start \
         {start}.",
        cigar_text(cigar)
    ))
}

fn optional_int(record: &BamRecord, tag: &[u8; 2]) -> i32 {
    match record.tags.get(htsjdk_bam::tag::Tag::new(tag)) {
        Some(htsjdk_bam::tag::TagValue::Int(value)) => *value as i32,
        _ => NO_VALUE,
    }
}

impl Alignment {
    /// `AlignmentInterval(SAMRecord)`.
    fn from_record(record: &BamRecord, header: &SamHeader) -> Result<Self, Failed> {
        if record.flags & 0x4 != 0 {
            return Err(format!(
                "sam record being used to construct AlignmentInterval is unmapped: {}",
                record_text(record, header)
            ));
        }
        let reverse = record.flags & 0x10 != 0;
        let reference = header
            .sequences
            .get(record.reference_index as usize)
            .map_or("*", |sequence| sequence.name.as_str());
        let start = record.alignment_start;
        let span = Interval::new(
            reference,
            start,
            start + reference_length(&record.cigar.elements) - 1,
        )?;
        let mut cigar = record.cigar.elements.clone();
        if reverse {
            cigar.reverse();
        }
        let left = clipped_bases(&cigar, true)?;
        let right = clipped_bases(&cigar, false)?;
        let unclipped: i32 = cigar
            .iter()
            .filter(|element| is_clipping(element.op) || consumes_read(element.op))
            .map(|element| element.length as i32)
            .sum();
        Ok(Self {
            span,
            start_on_contig: 1 + left,
            end_on_contig: unclipped - right,
            cigar,
            forward: !reverse,
            map_quality: i32::from(record.mapping_quality),
            mismatches: optional_int(record, b"NM"),
            score: optional_int(record, b"AS"),
            modification: 'O',
        })
    }

    /// `AlignmentInterval.toPackedString`.
    fn packed(&self) -> String {
        [
            self.start_on_contig.to_string(),
            self.end_on_contig.to_string(),
            self.span.text(),
            (if self.forward { "+" } else { "-" }).to_string(),
            cigar_text(&self.cigar),
            self.map_quality.to_string(),
            self.mismatches.to_string(),
            self.score.to_string(),
            self.modification.to_string(),
        ]
        .join("_")
    }

    /// `AlignmentInterval.toString`, the SA form.
    fn sa_text(&self) -> String {
        let mut cigar = self.cigar.clone();
        if !self.forward {
            cigar.reverse();
        }
        let mut text = format!(
            "{},{},{},{},{}",
            self.span.contig,
            self.span.start,
            if self.forward { "+" } else { "-" },
            cigar_text(&cigar),
            self.map_quality
        );
        if self.mismatches != NO_VALUE || self.score != NO_VALUE {
            text.push(',');
            if self.mismatches == NO_VALUE {
                text.push_str(&format!(".,{}", self.score));
            } else if self.score == NO_VALUE {
                text.push_str(&self.mismatches.to_string());
            } else {
                text.push_str(&format!("{},{}", self.mismatches, self.score));
            }
        }
        text
    }
}

fn overlap_on_contig(one: &Alignment, two: &Alignment) -> i32 {
    ((one.end_on_contig + 1).min(two.end_on_contig + 1)
        - one.start_on_contig.max(two.start_on_contig))
    .max(0)
}

/// `AlignmentInterval.overlapOnRefSpan`.
fn overlap_on_reference(one: &Alignment, two: &Alignment) -> i32 {
    if one.span.contig != two.span.contig {
        return 0;
    }
    ((one.span.end + 1).min(two.span.end + 1) - one.span.start.max(two.span.start)).max(0)
}

/// `ContigAlignmentsModifier.GapSplitter.splitGaps`.
fn split_gaps(one: Alignment, unclipped_length: i32) -> Result<Vec<Alignment>, Failed> {
    const NOT_SET: i32 = -1;
    let elements = one.cigar.clone();
    let mut count = elements.len();
    if count <= 1 {
        return Ok(vec![one]);
    }
    let mut contig_offset = 0;
    let mut index = 0;
    if elements[0].op == Op::H {
        contig_offset += elements[0].length as i32;
        index = 1;
    }
    if elements[count - 1].op == Op::H {
        count -= 1;
    }
    let mut reference_offset = if one.forward {
        one.span.start
    } else {
        -one.span.end
    };
    let (mut start_contig, mut start_index, mut start_reference) = (NOT_SET, 0, 0);
    let (mut end_contig, mut end_index, mut end_reference) = (NOT_SET, 0, 0);
    let mut chunks = Vec::new();
    while index < count {
        let element = elements[index];
        let length = element.length as i32;
        if is_alignment(element.op) {
            if start_contig == NOT_SET {
                start_contig = contig_offset;
                start_index = index;
                start_reference = reference_offset;
            }
            end_contig = contig_offset + length;
            end_index = index + 1;
            end_reference = reference_offset + length;
        } else if is_indel(element.op) && length >= GAP_SENSITIVITY && start_contig != NOT_SET {
            chunks.push((
                start_contig,
                start_index,
                start_reference,
                end_contig,
                end_index,
                end_reference,
            ));
            start_contig = NOT_SET;
        }
        if consumes_read(element.op) {
            contig_offset += length;
        }
        if consumes_reference(element.op) {
            reference_offset += length;
        }
        index += 1;
    }
    if start_contig != NOT_SET {
        chunks.push((
            start_contig,
            start_index,
            start_reference,
            end_contig,
            end_index,
            end_reference,
        ));
    }
    if chunks.len() < 2 {
        return Ok(vec![one]);
    }
    let soft = |length: i32| CigarElement {
        length: length as u32,
        op: Op::S,
    };
    let mut result = Vec::with_capacity(chunks.len());
    for (start_contig, start_index, start_reference, end_contig, end_index, end_reference) in chunks
    {
        let mut cigar = Vec::new();
        let mut initial_soft = start_contig;
        if elements[0].op == Op::H {
            cigar.push(elements[0]);
            initial_soft -= elements[0].length as i32;
        }
        if initial_soft > 0 {
            cigar.push(soft(initial_soft));
        }
        cigar.extend_from_slice(&elements[start_index..end_index]);
        let mut final_soft = unclipped_length - end_contig;
        let last = elements[elements.len() - 1];
        if last.op == Op::H {
            final_soft -= last.length as i32;
            if final_soft > 0 {
                cigar.push(soft(final_soft));
            }
            cigar.push(last);
        } else if final_soft > 0 {
            cigar.push(soft(final_soft));
        }
        let span = if one.forward {
            Interval::new(&one.span.contig, start_reference, end_reference - 1)?
        } else {
            Interval::new(&one.span.contig, -end_reference + 1, -start_reference)?
        };
        result.push(Alignment {
            span,
            start_on_contig: start_contig + 1,
            end_on_contig: end_contig,
            cigar,
            forward: one.forward,
            map_quality: one.map_quality,
            mismatches: NO_VALUE,
            score: NO_VALUE,
            modification: 'S',
        });
    }
    Ok(result)
}

/// `AlignedContig`.
#[derive(Clone)]
struct Contig {
    name: String,
    sequence: Vec<u8>,
    alignments: Vec<Alignment>,
}

/// `SequenceUtil.reverseComplement`.
fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => *other,
        })
        .collect()
}

/// `SAMRecord.toString`.
fn record_text(record: &BamRecord, header: &SamHeader) -> String {
    let mut text = record.read_name.clone();
    if record.flags & 0x1 != 0 {
        text.push_str(if record.flags & 0x40 != 0 {
            " 1/2"
        } else {
            " 2/2"
        });
    }
    text.push_str(&format!(" {}b", record.read_bases.len()));
    if record.flags & 0x4 != 0 {
        text.push_str(" unmapped read.");
    } else {
        let contig = header
            .sequences
            .get(record.reference_index as usize)
            .map_or("*", |sequence| sequence.name.as_str());
        let end = record.alignment_start + reference_length(&record.cigar.elements) - 1;
        text.push_str(&format!(
            " aligned to {contig}:{}-{end}.",
            record.alignment_start
        ));
    }
    text
}

/// `AlignedContig.parseReadsAndOptionallySplitGappedAlignments`.
fn parse_contig(
    reads: &[&BamRecord],
    header: &SamHeader,
    split_gapped: bool,
) -> Result<Contig, Failed> {
    let primary = reads
        .iter()
        .find(|read| read.flags & 0x800 == 0)
        .ok_or_else(|| format!("no primary alignment for read {}", reads[0].read_name))?;
    if primary
        .cigar
        .elements
        .iter()
        .any(|element| element.op == Op::H)
    {
        return Err(format!(
            "assumption that primary alignment does not contain hard clipping is invalid for \
             read: {}",
            record_text(primary, header)
        ));
    }
    let mut sequence = primary.read_bases.clone();
    let mut alignments = Vec::new();
    if primary.flags & 0x4 == 0 {
        if primary.flags & 0x10 != 0 {
            sequence = reverse_complement(&sequence);
        }
        let unclipped = primary.read_bases.len() as i32;
        for read in reads {
            let alignment = Alignment::from_record(read, header)?;
            if split_gapped {
                alignments.extend(split_gaps(alignment, unclipped)?);
            } else {
                alignments.push(alignment);
            }
        }
    }
    // `getAlignmentIntervalComparator`, in a stable sort.
    alignments.sort_by(|one, two| {
        one.start_on_contig
            .cmp(&two.start_on_contig)
            .then_with(|| one.span.contig.cmp(&two.span.contig))
            .then(one.span.start.cmp(&two.span.start))
    });
    Ok(Contig {
        name: primary.read_name.clone(),
        sequence,
        alignments,
    })
}

/// `StrandSwitch`, whose ordinal `hashCode` reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StrandSwitch {
    None,
    ForwardToReverse,
    ReverseToForward,
}

/// `SimpleChimera`.
#[derive(Debug, Clone)]
struct Chimera {
    contig: String,
    lower: Alignment,
    higher: Alignment,
    switch: StrandSwitch,
    forward_representation: bool,
    insertion_mappings: Vec<String>,
    /// `goodNonCanonicalMappingSATag`.
    non_canonical: String,
}

/// `TypeInferredFromSimpleChimera`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChimeraType {
    InterChromosomeStrandSwitch,
    InterChromosomeLeftMateFirst,
    InterChromosomeLeftMateSecond,
    ReferenceOrderSwap,
    StrandSwitch,
    SimpleDeletion,
    Replacement,
    SimpleInsertion,
    DuplicationExpansion,
    DuplicationComplex,
    DuplicationContraction,
}

/// `SimpleChimera.DistancesBetweenAlignmentsOnRefAndOnRead`.
struct Distances {
    gap_on_reference: i32,
    gap_on_contig: i32,
    left_end: i32,
    right_start: i32,
    first_contig_end: i32,
    second_contig_start: i32,
}

impl Chimera {
    fn new(
        lower: &Alignment,
        higher: &Alignment,
        mappings: &[String],
        contig: &str,
        non_canonical: &str,
        header: &SamHeader,
    ) -> Result<Self, Failed> {
        let switch = if lower.forward == higher.forward {
            StrandSwitch::None
        } else if lower.forward {
            StrandSwitch::ForwardToReverse
        } else {
            StrandSwitch::ReverseToForward
        };
        let same_contig = lower.span.contig == higher.span.contig;
        // `isForwardStrandRepresentation`.
        let forward_representation = match switch {
            StrandSwitch::None => lower.forward,
            _ if !same_contig => {
                contig_index(header, &lower.span.contig) < contig_index(header, &higher.span.contig)
            }
            StrandSwitch::ForwardToReverse => lower.span.end < higher.span.end,
            StrandSwitch::ReverseToForward => lower.span.start < higher.span.start,
        };
        Ok(Self {
            contig: contig.to_string(),
            lower: lower.clone(),
            higher: higher.clone(),
            switch,
            forward_representation,
            insertion_mappings: mappings.to_vec(),
            non_canonical: non_canonical.to_string(),
        })
    }

    fn switch_name(&self) -> &'static str {
        match self.switch {
            StrandSwitch::None => "NO_SWITCH",
            StrandSwitch::ForwardToReverse => "FORWARD_TO_REVERSE",
            StrandSwitch::ReverseToForward => "REVERSE_TO_FORWARD",
        }
    }

    /// `SimpleChimera.toString`.
    fn text(&self) -> String {
        format!(
            "SimpleChimera{{sourceContigName='{}', regionWithLowerCoordOnContig={}, \
             regionWithHigherCoordOnContig={}, strandSwitch={}, \
             isForwardStrandRepresentation={}, insertionMappings=[{}], \
             goodNonCanonicalMappingSATag='{}'}}",
            self.contig,
            self.lower.sa_text(),
            self.higher.sa_text(),
            self.switch_name(),
            self.forward_representation,
            self.insertion_mappings.join(", "),
            self.non_canonical
        )
    }

    /// `firstContigRegionRefSpanAfterSecond`.
    fn first_after_second(&self, header: &SamHeader) -> bool {
        compare_locatables(
            header,
            (
                &self.lower.span.contig,
                self.lower.span.start,
                self.lower.span.end,
            ),
            (
                &self.higher.span.contig,
                self.higher.span.start,
                self.higher.span.end,
            ),
        ) == std::cmp::Ordering::Greater
    }

    /// `getCoordinateSortedRefSpans`.
    fn sorted_spans(&self, header: &SamHeader) -> (&Interval, &Interval) {
        if self.first_after_second(header) {
            (&self.higher.span, &self.lower.span)
        } else {
            (&self.lower.span, &self.higher.span)
        }
    }

    /// `AssemblyContigWithFineTunedAlignments.hasIncompletePictureFromTwoAlignments`.
    fn incomplete(&self) -> bool {
        incomplete_from_two(&self.lower, &self.higher)
    }

    /// `isCandidateSimpleTranslocation`.
    fn translocation(&self) -> bool {
        if self.incomplete() {
            return false;
        }
        let (one, two) = (&self.lower.span, &self.higher.span);
        if one.contig != two.contig {
            return true;
        }
        if self.switch != StrandSwitch::None {
            return false;
        }
        if self.lower.forward {
            one.start > two.end
        } else {
            two.start > one.end
        }
    }

    /// `getDistancesBetweenAlignmentsOnRefAndOnRead`.
    fn distances(&self) -> Distances {
        let (left, right) = if self.forward_representation {
            (&self.lower.span, &self.higher.span)
        } else {
            (&self.higher.span, &self.lower.span)
        };
        Distances {
            gap_on_reference: right.start - left.end - 1,
            gap_on_contig: self.higher.start_on_contig - self.lower.end_on_contig - 1,
            left_end: left.end,
            right_start: right.start,
            first_contig_end: self.lower.end_on_contig,
            second_contig_start: self.higher.start_on_contig,
        }
    }

    /// `inferType`.
    fn infer_type(&self, header: &SamHeader) -> Result<ChimeraType, Failed> {
        if self.translocation() {
            return Ok(if self.lower.span.contig == self.higher.span.contig {
                ChimeraType::ReferenceOrderSwap
            } else if self.switch != StrandSwitch::None {
                ChimeraType::InterChromosomeStrandSwitch
            } else if self.forward_representation != self.first_after_second(header) {
                ChimeraType::InterChromosomeLeftMateFirst
            } else {
                ChimeraType::InterChromosomeLeftMateSecond
            });
        }
        if self.switch != StrandSwitch::None {
            return Ok(ChimeraType::StrandSwitch);
        }
        let distances = self.distances();
        let (on_reference, on_contig) = (distances.gap_on_reference, distances.gap_on_contig);
        Ok(if on_reference > 0 {
            if on_contig <= 0 {
                ChimeraType::SimpleDeletion
            } else {
                ChimeraType::Replacement
            }
        } else if on_reference < 0 {
            if on_contig >= 0 {
                ChimeraType::DuplicationExpansion
            } else {
                ChimeraType::DuplicationComplex
            }
        } else if on_contig > 0 {
            ChimeraType::SimpleInsertion
        } else if on_contig < 0 {
            ChimeraType::DuplicationContraction
        } else {
            return Err(format!(
                "Detected badly parsed chimeric alignment for identifying SV breakpoints; no \
                 rearrangement found: {}",
                self.text()
            ));
        })
    }
}

/// `AssemblyContigWithFineTunedAlignments.hasIncompletePictureFromTwoAlignments`.
fn incomplete_from_two(head: &Alignment, tail: &Alignment) -> bool {
    let (one, two) = (&head.span, &tail.span);
    if one.contig != two.contig {
        return false;
    }
    if one.contains(two) || two.contains(one) {
        return true;
    }
    if head.forward != tail.forward {
        one.overlaps(two)
    } else if head.forward {
        one.start > two.start && one.start <= two.end
    } else {
        two.start > one.start && two.start <= one.end
    }
}

/// `AlignedContig.simpleChimeraWithStichableAlignments`.
fn stitchable(one: &Alignment, two: &Alignment) -> Result<bool, Failed> {
    if one.start_on_contig > two.start_on_contig {
        return Err(format!(
            "Assumption that input intervals are sorted by their starts on read is \
             violated.\tFirst: {}\tSecond: {}",
            one.packed(),
            two.packed()
        ));
    }
    if one.span.contig != two.span.contig || one.forward != two.forward {
        return Ok(false);
    }
    let contains_on_read = |a: &Alignment, b: &Alignment| {
        a.start_on_contig <= b.start_on_contig && a.end_on_contig >= b.end_on_contig
    };
    if contains_on_read(one, two) || contains_on_read(two, one) {
        return Ok(false);
    }
    if one.span.contains(&two.span) || two.span.contains(&one.span) {
        return Ok(false);
    }
    if one.forward != (one.span.start < two.span.start) {
        return Ok(false);
    }
    let on_contig = overlap_on_contig(one, two);
    let on_reference = overlap_on_reference(one, two);
    Ok(if on_contig == 0 && on_reference == 0 {
        two.span.start - one.span.end == 1 && two.start_on_contig - one.end_on_contig == 1
    } else {
        on_contig == on_reference
    })
}

fn first_too_short(first: &Alignment, second: &Alignment) -> bool {
    first.span.size() - overlap_on_contig(first, second) < MIN_ALIGNMENT_LENGTH
}

/// `nextAlignmentMayBeInsertion`, wholly contained alignments filtered as the tool asks.
fn next_may_be_insertion(current: &Alignment, next: &Alignment) -> bool {
    next.map_quality < HIGH_MQ
        || first_too_short(next, current)
        || current.span.contains(&next.span)
        || next.span.contains(&current.span)
}

/// `parseOneContig`.
fn parse_one_contig(contig: &Contig, header: &SamHeader) -> Result<Vec<Chimera>, Failed> {
    if contig.alignments.len() < 2 {
        return Ok(Vec::new());
    }
    let mut position = 0;
    while contig.alignments[position].map_quality < HIGH_MQ
        && position + 1 < contig.alignments.len()
    {
        position += 1;
    }
    let mut current = &contig.alignments[position];
    let mut results: Vec<Chimera> = Vec::new();
    let mut mappings: Vec<String> = Vec::new();
    let mut walk = position + 1;
    while walk < contig.alignments.len() {
        let next = &contig.alignments[walk];
        walk += 1;
        if first_too_short(current, next) || stitchable(current, next)? {
            continue;
        }
        if next_may_be_insertion(current, next) {
            if walk < contig.alignments.len() {
                mappings.push(next.packed());
                continue;
            }
            break;
        }
        let chimera = Chimera::new(current, next, &mappings, &contig.name, "NONE", header)?;
        if !chimera.incomplete() && !chimera.translocation() {
            results.push(chimera);
        }
        current = next;
    }
    // One list, shared, read once it is whole.
    for chimera in &mut results {
        chimera.insertion_mappings = mappings.clone();
    }
    Ok(results)
}

/// `BreakpointComplications`, of the classes this tool's chimeras can reach.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Duplication {
    None,
    Precise { cigars: Vec<String> },
    Imprecise { affected: Interval },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ComplicationClass {
    SimpleOrDuplication,
    StrandSwitch,
    ReferenceOrderSwap,
    InterChromosome,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Complication {
    /// The class, which `equals` compares: the classes without a duplication have the same fields.
    class: ComplicationClass,
    homology: String,
    inserted: String,
    unit: Option<Interval>,
    copies_on_reference: i32,
    copies_on_contig: i32,
    duplication: Duplication,
}

/// `SimpleInterval.hashCode` and the rest, folded the way `hashCode` folds them.
fn fold(result: i32, value: i32) -> i32 {
    result.wrapping_mul(31).wrapping_add(value)
}

impl Complication {
    fn has_duplication(&self) -> bool {
        self.unit.is_some()
    }

    fn contraction(&self) -> bool {
        self.copies_on_reference > self.copies_on_contig
    }

    fn hash_code(&self) -> i32 {
        let mut result = fold(
            string_hash_code(&self.homology),
            string_hash_code(&self.inserted),
        );
        let Some(unit) = &self.unit else {
            return result;
        };
        result = fold(result, unit.hash_code());
        result = fold(result, self.copies_on_reference);
        result = fold(result, self.copies_on_contig);
        // Each `Strand.POSITIVE`, whose ordinal is zero.
        for _ in 0..self.copies_on_reference + self.copies_on_contig {
            result = fold(result, 0);
        }
        match &self.duplication {
            Duplication::Precise { cigars } => {
                let mut list = 1;
                for cigar in cigars {
                    list = fold(list, string_hash_code(cigar));
                }
                fold(result, list)
            }
            Duplication::Imprecise { affected } => fold(result, affected.hash_code()),
            Duplication::None => result,
        }
    }

    /// `toVariantAttributes`.
    fn attributes(&self) -> Vec<(&'static str, Value)> {
        let mut attributes = Vec::new();
        if !self.inserted.is_empty() {
            attributes.push(("INSSEQ", Value::Str(self.inserted.clone())));
            attributes.push(("INSLEN", Value::Int(self.inserted.len() as i64)));
        }
        if !self.homology.is_empty() {
            attributes.push(("HOMSEQ", Value::Str(self.homology.clone())));
            attributes.push(("HOMLEN", Value::Int(self.homology.len() as i64)));
        }
        if let Some(unit) = &self.unit {
            attributes.push(("DUP_REPEAT_UNIT_REF_SPAN", Value::Str(unit.text())));
            attributes.push((
                "DUP_NUM",
                Value::Str(format!(
                    "{},{}",
                    self.copies_on_reference, self.copies_on_contig
                )),
            ));
            attributes.push((
                "DUP_ORIENTATIONS",
                Value::Str("+".repeat(self.copies_on_contig as usize)),
            ));
            attributes.push((
                if self.copies_on_reference < self.copies_on_contig {
                    "EXPANSION"
                } else {
                    "CONTRACTION"
                },
                Value::Str(String::new()),
            ));
            match &self.duplication {
                Duplication::Precise { cigars } if !cigars.is_empty() => {
                    attributes.push(("DUP_SEQ_CIGARS", Value::Str(cigars.join(","))));
                }
                Duplication::Imprecise { affected } => {
                    attributes.push(("DUP_ANNOTATIONS_IMPRECISE", Value::Str(String::new())));
                    attributes.push(("DUP_IMPRECISE_AFFECTED_RANGE", Value::Str(affected.text())));
                }
                _ => {}
            }
        }
        attributes
    }
}

/// `Arrays.copyOfRange`, which pads past the end with zeros.
fn copy_of_range(bases: &[u8], from: i32, to: i32) -> Result<Vec<u8>, Failed> {
    if from > to {
        return Err(format!("{from} > {to}"));
    }
    if from < 0 || from as usize > bases.len() {
        return Err(format!("Array index out of range: {from}"));
    }
    let mut copied = vec![0u8; (to - from) as usize];
    let available = (bases.len() - from as usize).min(copied.len());
    copied[..available].copy_from_slice(&bases[from as usize..from as usize + available]);
    Ok(copied)
}

/// `reverseComplementIfNecessary`.
fn oriented(bases: Vec<u8>, first: &Alignment, second: &Alignment, first_after: bool) -> Vec<u8> {
    let flip = if first.forward == second.forward {
        !first.forward
    } else {
        first_after == first.forward
    };
    if flip {
        reverse_complement(&bases)
    } else {
        bases
    }
}

fn text_of(bases: Vec<u8>) -> String {
    String::from_utf8_lossy(&bases).into_owned()
}

/// `inferHomology`.
fn infer_homology(
    first: &Alignment,
    second: &Alignment,
    sequence: &[u8],
    first_after: bool,
) -> Result<String, Failed> {
    if first.end_on_contig >= second.start_on_contig {
        let bases = copy_of_range(sequence, second.start_on_contig - 1, first.end_on_contig)?;
        Ok(text_of(oriented(bases, first, second, first_after)))
    } else {
        Ok(String::new())
    }
}

/// `inferInsertedSequence`.
fn infer_inserted(
    first: &Alignment,
    second: &Alignment,
    sequence: &[u8],
    first_after: bool,
) -> Result<String, Failed> {
    if first.end_on_contig < second.start_on_contig - 1 {
        let bases = copy_of_range(sequence, first.end_on_contig, second.start_on_contig - 1)?;
        Ok(text_of(oriented(bases, first, second, first_after)))
    } else {
        Ok(String::new())
    }
}

/// `extractCigarForTandupExpansion`.
fn tandem_expansion_cigar(
    region: &Alignment,
    one_end: i32,
    two_begin: i32,
) -> Result<Vec<CigarElement>, Failed> {
    let mut result = Vec::with_capacity(region.cigar.len());
    let mut collecting = false;
    let mut position = if region.forward {
        region.span.start
    } else {
        region.span.end
    };
    for element in &region.cigar {
        if is_clipping(element.op) {
            continue;
        }
        let length = element.length as i32;
        if consumes_reference(element.op) {
            position += if region.forward { length } else { -length };
        }
        let (offset, overshoot) = if region.forward {
            (position - two_begin, position - one_end - 1)
        } else {
            (one_end - position, two_begin - position - 1)
        };
        if offset > 0 {
            if overshoot <= 0 {
                result.push(if collecting {
                    *element
                } else {
                    CigarElement {
                        length: offset as u32,
                        op: element.op,
                    }
                });
                collecting = true;
            } else {
                if length < overshoot {
                    return Err(format!(
                        "Cigar element being constructed with negative length: {}",
                        length - overshoot
                    ));
                }
                result.push(CigarElement {
                    length: (length - overshoot) as u32,
                    op: element.op,
                });
                break;
            }
        }
    }
    Ok(result)
}

/// `SmallDuplicationWithImpreciseDupRangeBreakpointComplications.TandemRepeatStructure`.
struct TandemRepeat {
    lower: i32,
    higher: i32,
    unit_length: i32,
    pseudo_homology: i32,
}

fn tandem_repeat(on_reference: i32, on_contig: i32) -> TandemRepeat {
    let expansion = on_reference < on_contig;
    let (on_lower, on_higher) = if expansion {
        (on_reference.abs(), on_contig.abs())
    } else {
        (on_contig.abs(), on_reference.abs())
    };
    let mut best = TandemRepeat {
        lower: 0,
        higher: 0,
        unit_length: 0,
        pseudo_homology: 0,
    };
    let mut error = i64::MAX;
    for lower in 1..10 {
        for higher in lower + 1..=2 * lower {
            let bound = if higher == 2 * lower {
                on_lower
            } else {
                on_higher
            };
            for unit in 2..=bound {
                for lambda in 0..unit {
                    let d1 = (2 * lower - higher) * unit + lambda;
                    let d2 = lower * unit + lambda;
                    let new_error = i64::from((on_higher - d1).abs() + (on_lower - d2).abs());
                    if new_error < error {
                        error = new_error;
                        best = TandemRepeat {
                            lower,
                            higher,
                            unit_length: unit,
                            pseudo_homology: lambda,
                        };
                    }
                    if error < 1 {
                        return best;
                    }
                }
            }
        }
    }
    best
}

/// `NovelAdjacencyAndAltHaplotype`.
#[derive(Debug, Clone)]
struct NovelAdjacency {
    left: Interval,
    right: Interval,
    switch: StrandSwitch,
    complication: Complication,
    /// `type`, which `equals` and `hashCode` leave out.
    kind: ChimeraType,
    alt_haplotype: Vec<u8>,
}

impl PartialEq for NovelAdjacency {
    fn eq(&self, other: &Self) -> bool {
        self.left == other.left
            && self.right == other.right
            && self.switch == other.switch
            && self.complication == other.complication
            && self.alt_haplotype == other.alt_haplotype
    }
}

impl NovelAdjacency {
    fn hash_code(&self) -> i32 {
        let mut result = self.left.hash_code();
        result = fold(result, self.right.hash_code());
        result = fold(
            result,
            match self.switch {
                StrandSwitch::None => 0,
                StrandSwitch::ForwardToReverse => 1,
                StrandSwitch::ReverseToForward => 2,
            },
        );
        result = fold(result, self.complication.hash_code());
        fold(result, byte_array_hash_code(&self.alt_haplotype))
    }

    /// `makeLocationString` of this adjacency.
    fn location(&self) -> String {
        location_string(
            &self.left.contig,
            self.left.start,
            &self.right.contig,
            self.right.end,
        )
    }

    /// `toString`, for a message.
    fn text(&self) -> String {
        format!(
            "{}\t{}\t{}\thomology: {}\tinserted sequence: {}",
            self.left.text(),
            self.right.text(),
            match self.switch {
                StrandSwitch::None => "NO_SWITCH",
                StrandSwitch::ForwardToReverse => "FORWARD_TO_REVERSE",
                StrandSwitch::ReverseToForward => "REVERSE_TO_FORWARD",
            },
            self.complication.homology,
            self.complication.inserted
        )
    }

    /// `getLengthForDupTandemExpansion`.
    fn expansion_length(&self) -> Result<i32, Failed> {
        if self.complication.contraction() {
            return Err(format!(
                "Trying to extract length from a duplication contraction: {}",
                self.text()
            ));
        }
        let unit = self.complication.unit.as_ref().map_or(0, Interval::size);
        Ok(self.complication.inserted.len() as i32
            + (self.complication.copies_on_contig - self.complication.copies_on_reference) * unit)
    }
}

fn location_string(one: &str, position_one: i32, two: &str, position_two: i32) -> String {
    if one == two {
        format!("{one}_{position_one}_{position_two}")
    } else {
        format!("{one}_{position_one}_{two}_{position_two}")
    }
}

/// `BreakpointsInference.validateInferredLocations`.
fn validate_locations(left: &Interval, right: &Interval, header: &SamHeader) -> Result<(), Failed> {
    let (left_index, right_index) = (
        contig_index(header, &left.contig),
        contig_index(header, &right.contig),
    );
    if right_index < left_index || (right_index == left_index && right.end < left.start) {
        return Err(
            "Inferred novel adjacency reference locations have left location after right \
             location."
                .to_string(),
        );
    }
    let length = |index: i32| header.sequences[index as usize].length;
    if left.end > length(left_index) || right.end > length(right_index) {
        return Err("Inferred breakpoint beyond reference sequence length.".to_string());
    }
    Ok(())
}

/// `new NovelAdjacencyAndAltHaplotype(chimera, contigSequence, dictionary)`, through the inference
/// and complication classes the chimera's type selects.
fn novel_adjacency(
    chimera: &Chimera,
    sequence: &[u8],
    header: &SamHeader,
) -> Result<NovelAdjacency, Failed> {
    let kind = chimera.infer_type(header)?;
    let first_after = chimera.first_after_second(header);
    let (first, second) = (&chimera.lower, &chimera.higher);
    let (left_span, right_span) = chimera.sorted_spans(header);
    let reference = first.span.contig.clone();
    let mut right_contig = reference.clone();
    let mut left_contig = reference.clone();
    let mut complication = Complication {
        class: ComplicationClass::SimpleOrDuplication,
        homology: String::new(),
        inserted: String::new(),
        unit: None,
        copies_on_reference: 0,
        copies_on_contig: 0,
        duplication: Duplication::None,
    };
    let (upstream, downstream, alt_haplotype);
    let mut validate = true;
    match kind {
        ChimeraType::SimpleDeletion | ChimeraType::Replacement | ChimeraType::SimpleInsertion => {
            let distances = chimera.distances();
            if distances.gap_on_reference > 0 && distances.gap_on_contig < 0 {
                complication.homology = infer_homology(first, second, sequence, first_after)?;
            } else {
                complication.inserted = infer_inserted(first, second, sequence, first_after)?;
            }
            upstream = left_span.end - complication.homology.len() as i32;
            downstream = right_span.start - 1;
            alt_haplotype = complication.inserted.clone().into_bytes();
            validate = false;
        }
        ChimeraType::DuplicationExpansion | ChimeraType::DuplicationContraction => {
            let distances = chimera.distances();
            if kind == ChimeraType::DuplicationExpansion {
                if distances.gap_on_contig != 0 {
                    complication.inserted = infer_inserted(first, second, sequence, first_after)?;
                }
                let unit = Interval::new(&reference, distances.right_start, distances.left_end)?;
                let mut cigars = Vec::with_capacity(2);
                let mut lengths = Vec::with_capacity(2);
                for region in [first, second] {
                    let mut cigar =
                        tandem_expansion_cigar(region, distances.left_end, distances.right_start)?;
                    if !first.forward {
                        cigar.reverse();
                    }
                    cigars.push(cigar_text(&cigar));
                    lengths.push(read_length(&cigar));
                }
                if !first.forward {
                    cigars.reverse();
                    lengths.reverse();
                }
                complication.unit = Some(unit);
                complication.copies_on_reference = 1;
                complication.copies_on_contig = 2;
                let (first_copy, second_copy) = if chimera.forward_representation {
                    (lengths[0], lengths[1])
                } else {
                    (lengths[1], lengths[0])
                };
                let from = distances.first_contig_end - first_copy;
                let to = distances.second_contig_start + second_copy - 1;
                let mut bases = copy_of_range(sequence, from, to)?;
                if !chimera.forward_representation {
                    bases = reverse_complement(&bases);
                }
                alt_haplotype = bases;
                let unit_size = complication.unit.as_ref().map_or(0, Interval::size);
                upstream = left_span.end - unit_size;
                downstream = right_span.start - 1;
                complication.duplication = Duplication::Precise { cigars };
            } else {
                complication.homology = infer_homology(first, second, sequence, first_after)?;
                complication.unit = Some(Interval::new(
                    &reference,
                    distances.left_end
                        - (distances.first_contig_end - distances.second_contig_start),
                    distances.left_end,
                )?);
                complication.copies_on_reference = 2;
                complication.copies_on_contig = 1;
                complication.duplication = Duplication::Precise { cigars: Vec::new() };
                upstream = left_span.end - complication.homology.len() as i32;
                downstream = right_span.start - 1;
                alt_haplotype = Vec::new();
            }
        }
        ChimeraType::DuplicationComplex => {
            let distances = chimera.distances();
            let repeat = tandem_repeat(distances.gap_on_reference, distances.gap_on_contig);
            let expansion = distances.gap_on_reference < distances.gap_on_contig;
            let unit_start =
                distances.left_end - repeat.pseudo_homology - repeat.unit_length * repeat.lower + 1;
            complication.homology = infer_homology(first, second, sequence, first_after)?;
            let unit = Interval::new(
                &first.span.contig,
                unit_start,
                unit_start + repeat.unit_length - 1,
            )?;
            let unit_size = unit.size();
            complication.unit = Some(unit);
            (
                complication.copies_on_reference,
                complication.copies_on_contig,
            ) = if expansion {
                (repeat.lower, repeat.higher)
            } else {
                (repeat.higher, repeat.lower)
            };
            let (left, right) = if chimera.forward_representation {
                (&first.span, &second.span)
            } else {
                (&second.span, &first.span)
            };
            let affected = if expansion {
                Interval::new(&left.contig, distances.right_start, distances.left_end)?
            } else {
                Interval::new(
                    &left.contig,
                    right.start - repeat.unit_length,
                    left.end + repeat.unit_length,
                )?
            };
            complication.duplication = Duplication::Imprecise { affected };
            let homology = complication.homology.len() as i32;
            let mut bases;
            if complication.contraction() {
                upstream = left_span.end - homology;
                downstream = right_span.start - 1;
                bases = copy_of_range(
                    sequence,
                    distances.second_contig_start - 1,
                    distances.first_contig_end,
                )?;
            } else {
                let difference = complication.copies_on_contig - complication.copies_on_reference;
                upstream = left_span.end - homology - difference * unit_size;
                downstream = right_span.start - 1;
                let hard_offset = |cigar: &[CigarElement]| match cigar.first() {
                    Some(element) if element.op == Op::H => element.length as i32,
                    _ => 0,
                };
                let walked = associated_distance_on_read(
                    &first.cigar,
                    distances.first_contig_end - hard_offset(&first.cigar),
                    -distances.gap_on_reference,
                    true,
                )?;
                let from = distances.first_contig_end - walked;
                let walked = associated_distance_on_read(
                    &second.cigar,
                    distances.second_contig_start - hard_offset(&second.cigar),
                    -distances.gap_on_reference,
                    false,
                )?;
                let to = distances.second_contig_start + walked - 1;
                bases = copy_of_range(sequence, from, to)?;
            }
            if !chimera.forward_representation {
                bases = reverse_complement(&bases);
            }
            alt_haplotype = bases;
        }
        ChimeraType::ReferenceOrderSwap => {
            complication.class = ComplicationClass::ReferenceOrderSwap;
            complication.homology = infer_homology(first, second, sequence, first_after)?;
            complication.inserted = infer_inserted(first, second, sequence, first_after)?;
            upstream = left_span.start;
            downstream = right_span.end - complication.homology.len() as i32;
            alt_haplotype = Vec::new();
        }
        ChimeraType::InterChromosomeStrandSwitch
        | ChimeraType::InterChromosomeLeftMateFirst
        | ChimeraType::InterChromosomeLeftMateSecond => {
            complication.class = ComplicationClass::InterChromosome;
            complication.homology = infer_homology(first, second, sequence, first_after)?;
            complication.inserted = infer_inserted(first, second, sequence, first_after)?;
            let homology = complication.homology.len() as i32;
            // `isFirstInPartner`.
            let first_in_partner = match chimera.switch {
                StrandSwitch::None => {
                    contig_index(header, &first.span.contig)
                        < contig_index(header, &second.span.contig)
                }
                _ => chimera.forward_representation,
            };
            let (one, two) = if first_in_partner {
                (first, second)
            } else {
                (second, first)
            };
            left_contig = one.span.contig.clone();
            right_contig = two.span.contig.clone();
            (upstream, downstream) = match chimera.switch {
                StrandSwitch::None => {
                    // `first` always keeps the homology; which end depends on the representation.
                    let first_end = if chimera.forward_representation {
                        first.span.end - homology
                    } else {
                        first.span.start
                    };
                    let second_end = if chimera.forward_representation {
                        second.span.start
                    } else {
                        second.span.end - homology
                    };
                    if first_in_partner {
                        (first_end, second_end)
                    } else {
                        (second_end, first_end)
                    }
                }
                StrandSwitch::ForwardToReverse => (one.span.end - homology, two.span.end),
                StrandSwitch::ReverseToForward => (one.span.start, two.span.start + homology),
            };
            alt_haplotype = Vec::new();
        }
        ChimeraType::StrandSwitch => {
            complication.class = ComplicationClass::StrandSwitch;
            complication.homology = infer_homology(first, second, sequence, first_after)?;
            complication.inserted = infer_inserted(first, second, sequence, first_after)?;
            let homology = complication.homology.len() as i32;
            if chimera.switch == StrandSwitch::ForwardToReverse {
                upstream = left_span.end - homology;
                downstream = right_span.end;
            } else {
                upstream = left_span.start;
                downstream = right_span.start + homology;
            }
            alt_haplotype = Vec::new();
        }
    }
    let left = Interval::new(&left_contig, upstream, upstream)?;
    let right = Interval::new(&right_contig, downstream, downstream)?;
    if validate {
        validate_locations(&left, &right, header)?;
    }
    Ok(NovelAdjacency {
        left,
        right,
        switch: chimera.switch,
        complication,
        kind,
        alt_haplotype,
    })
}

/// Spark's `AppendOnlyMap` walked in slot order, as `ExternalAppendOnlyMap.iterator` walks it when
/// nothing has spilled. A key keeps its first arrival's object.
pub fn append_only_group<K: PartialEq, V>(
    pairs: Vec<(K, V)>,
    hash_of: impl Fn(&K) -> i32,
) -> Vec<(K, Vec<V>)> {
    let slot_hash = |hash: i32| crate::pathseq_filter::murmur3_32_int(hash) as u32 as usize;
    let mut capacity = 64usize;
    let mut slots: Vec<Option<usize>> = vec![None; capacity];
    let mut groups: Vec<(K, Vec<V>, i32)> = Vec::new();
    for (key, value) in pairs {
        let hash = hash_of(&key);
        let mask = capacity - 1;
        let mut position = slot_hash(hash) & mask;
        let mut delta = 1;
        loop {
            match slots[position] {
                None => {
                    slots[position] = Some(groups.len());
                    groups.push((key, vec![value], hash));
                    if groups.len() > (0.7 * capacity as f64) as usize {
                        let doubled = capacity * 2;
                        let mut grown = vec![None; doubled];
                        for slot in slots.iter().flatten() {
                            let mut at = slot_hash(groups[*slot].2) & (doubled - 1);
                            let mut step = 1;
                            while grown[at].is_some() {
                                at = (at + step) & (doubled - 1);
                                step += 1;
                            }
                            grown[at] = Some(*slot);
                        }
                        slots = grown;
                        capacity = doubled;
                    }
                    break;
                }
                Some(group) if groups[group].2 == hash && groups[group].0 == key => {
                    groups[group].1.push(value);
                    break;
                }
                Some(_) => {
                    position = (position + delta) & mask;
                    delta += 1;
                }
            }
        }
    }
    let order: Vec<usize> = slots.into_iter().flatten().collect();
    let mut taken: Vec<Option<(K, Vec<V>)>> = groups
        .into_iter()
        .map(|(key, values, _)| Some((key, values)))
        .collect();
    order
        .into_iter()
        .map(|group| taken[group].take().expect("each group once"))
        .collect()
}

/// What `produceAnnotatedVcFromAssemblyEvidence` needs from outside.
pub trait Reference {
    /// `BasicReference.getBases`.
    fn bases(&mut self, contig: &str, start: i32, end: i32) -> Result<Vec<u8>, Failed>;
}

/// One external call of `--cnv-calls`, as `SVIntervalTree` holds it.
#[derive(Debug, Clone)]
pub struct CnvCall {
    pub contig: i32,
    pub start: i32,
    pub end: i32,
    /// `formatExternalCNVCallAnnotation`.
    pub annotation: String,
}

/// `inferSimpleTypeFromNovelAdjacency` and `produceAnnotatedVcFromAssemblyEvidence`.
fn annotated_variant(
    adjacency: &NovelAdjacency,
    evidence: &[Chimera],
    reference: &mut dyn Reference,
    header: &SamHeader,
    cnv_calls: Option<&[CnvCall]>,
) -> Result<VariantContext, Failed> {
    let start = adjacency.left.end;
    let end = adjacency.right.start;
    let has_inserted = !adjacency.complication.inserted.is_empty();
    let has_duplication = adjacency.complication.has_duplication();
    let contig = &adjacency.left.contig;
    let position = adjacency.left.start;
    let mut extra: Vec<(&'static str, Value)> = Vec::new();
    let (stop, id, symbol, svtype, length);
    if adjacency.switch == StrandSwitch::None {
        if start == end {
            if has_duplication {
                let unit = adjacency
                    .complication
                    .unit
                    .as_ref()
                    .expect("a duplication has its unit");
                id = format!(
                    "INS-DUPLICATION-TANDEM-EXPANSION_{}",
                    location_string(&unit.contig, unit.start, &unit.contig, unit.end)
                );
                length = adjacency.expansion_length()?;
                stop = position;
                symbol = "<DUP>";
                svtype = "DUP";
                extra.push(("EXPANSION", Value::Bool(true)));
            } else {
                if !has_inserted {
                    return Err(format!(
                        "Something went wrong in type inference, there's suspected insertion \
                         happening but no inserted sequence could be inferred {}",
                        adjacency.text()
                    ));
                }
                id = format!(
                    "INS_{}",
                    location_string(contig, position, contig, position)
                );
                length = adjacency.complication.inserted.len() as i32;
                stop = position;
                symbol = "<INS>";
                svtype = "INS";
            }
        } else {
            if has_duplication && has_inserted {
                return Err(format!(
                    "Something went wrong in novel adjacency interpretation:  inferring simple SV \
                     type from a novel adjacency between two different reference locations, but \
                     annotated with both inserted sequence and duplication, which is NOT \
                     simple.\n{}",
                    adjacency.text()
                ));
            }
            id = format!(
                "{}_{}",
                if has_duplication {
                    "DEL-DUPLICATION-TANDEM-CONTRACTION"
                } else {
                    "DEL"
                },
                adjacency.location()
            );
            length = -(adjacency.right.end - adjacency.left.start);
            stop = adjacency.right.end;
            symbol = "<DEL>";
            svtype = "DEL";
            if has_duplication {
                extra.push(("CONTRACTION", Value::Bool(true)));
            }
        }
    } else {
        let name = if adjacency.switch == StrandSwitch::ForwardToReverse {
            "INV55"
        } else {
            "INV33"
        };
        id = format!("{name}_{}", adjacency.location());
        length = adjacency.right.start - adjacency.left.end;
        stop = adjacency.right.end;
        symbol = "<INV>";
        svtype = "INV";
        extra.push((name, Value::Bool(true)));
    }
    let bases = reference.bases(contig, position, position)?;
    let reference_allele = Allele::create(&bases, true).map_err(|error| format!("{error:?}"))?;
    let alternate = Allele::from_str(symbol, false).map_err(|error| format!("{error:?}"))?;
    let mut variant = VariantContext::new(
        contig,
        i64::from(position),
        vec![reference_allele, alternate],
    );
    variant.stop = i64::from(stop);
    variant.id = id;

    set(&mut variant, "END", Value::Int(i64::from(stop)));
    set(&mut variant, "SVTYPE", Value::Str(svtype.to_string()));
    set(&mut variant, "SVLEN", Value::Int(i64::from(length)));
    for (key, value) in extra {
        set(&mut variant, key, value);
    }
    annotate_from_evidence(&mut variant, adjacency, evidence, header, cnv_calls, stop)?;
    Ok(variant)
}

/// A `VariantContextBuilder.attribute` on what is already built.
fn set(variant: &mut VariantContext, key: &str, value: Value) {
    match variant.attributes.iter_mut().find(|(name, _)| name == key) {
        Some(slot) => slot.1 = value,
        None => variant.attributes.push((key.to_string(), value)),
    }
}

/// What `produceAnnotatedVcFromAssemblyEvidence` adds to a type's basic information.
fn annotate_from_evidence(
    variant: &mut VariantContext,
    adjacency: &NovelAdjacency,
    evidence: &[Chimera],
    header: &SamHeader,
    cnv_calls: Option<&[CnvCall]>,
    stop: i32,
) -> Result<(), Failed> {
    let contig = variant.contig.clone();
    let position = variant.start as i32;
    for (key, value) in adjacency.complication.attributes() {
        set(variant, key, value);
    }

    // `getAssemblyEvidenceRelatedAnnotations`.
    let mut sorted: Vec<&Chimera> = evidence.iter().collect();
    sorted.sort_by(|one, two| one.contig.cmp(&two.contig));
    let minimum_quality =
        |chimera: &Chimera| chimera.lower.map_quality.min(chimera.higher.map_quality);
    let minimum_length = |chimera: &Chimera| {
        chimera.lower.span.size().min(chimera.higher.span.size())
            - overlap_on_contig(&chimera.lower, &chimera.higher)
    };
    let joined = |values: Vec<String>| values.join(",");
    set(variant, "TOTAL_MAPPINGS", Value::Int(sorted.len() as i64));
    set(
        variant,
        "HQ_MAPPINGS",
        Value::Int(
            sorted
                .iter()
                .filter(|chimera| minimum_quality(chimera) == HIGH_MQ)
                .count() as i64,
        ),
    );
    set(
        variant,
        "MAPPING_QUALITIES",
        Value::Str(joined(
            sorted
                .iter()
                .map(|chimera| minimum_quality(chimera).to_string())
                .collect(),
        )),
    );
    set(
        variant,
        "ALIGN_LENGTHS",
        Value::Str(joined(
            sorted
                .iter()
                .map(|chimera| minimum_length(chimera).to_string())
                .collect(),
        )),
    );
    set(
        variant,
        "MAX_ALIGN_LENGTH",
        Value::Int(i64::from(
            sorted
                .iter()
                .map(|chimera| minimum_length(chimera))
                .max()
                .unwrap_or(0),
        )),
    );
    set(
        variant,
        "CTG_NAMES",
        Value::Str(joined(
            sorted
                .iter()
                .map(|chimera| chimera.contig.clone())
                .collect(),
        )),
    );
    let mut mappings: Vec<String> = sorted
        .iter()
        .flat_map(|chimera| chimera.insertion_mappings.iter().cloned())
        .collect();
    mappings.sort();
    if !mappings.is_empty() {
        set(variant, "INSSEQ_MAP", Value::Str(mappings.join(",")));
    }
    let good_mappings: Vec<String> = sorted
        .iter()
        .filter(|chimera| chimera.non_canonical != "NONE")
        .map(|chimera| chimera.non_canonical.clone())
        .collect();
    if !good_mappings.is_empty() {
        set(
            variant,
            "CTG_GOOD_NONCANONICAL_MAPPING",
            Value::Str(good_mappings.join(",")),
        );
    }

    if !adjacency.alt_haplotype.is_empty() {
        set(
            variant,
            "SEQ_ALT_HAPLOTYPE",
            Value::Str(String::from_utf8_lossy(&adjacency.alt_haplotype).into_owned()),
        );
    }

    // `annotateWithExternalCNVCalls`, whose tree walks in interval order.
    if let Some(calls) = cnv_calls {
        // A breakend has no end, and `SVInterval` is built from the -1 it says so with.
        if stop < 0 {
            return Err(format!("provided end is negative: {stop}"));
        }
        let index = contig_index(header, &contig);
        let annotation: Vec<&str> = calls
            .iter()
            .filter(|call| call.contig == index && call.start < stop && position < call.end)
            .map(|call| call.annotation.as_str())
            .collect();
        if !annotation.is_empty() {
            set(
                variant,
                "EXTERNAL_CNV_CALLS",
                Value::Str(annotation.join(",")),
            );
        }
    }
    Ok(())
}

/// `SVVCFWriter.sortVariantsByCoordinate`, stable as `Stream.sorted` is.
fn sort_variants(variants: &mut [VariantContext], header: &SamHeader) {
    let inserted = |variant: &VariantContext| -> String {
        variant
            .attributes
            .iter()
            .find(|(name, _)| name == "INSSEQ")
            .and_then(|(_, value)| match value {
                Value::Str(text) => Some(text.clone()),
                _ => None,
            })
            .unwrap_or_default()
    };
    variants.sort_by(|one, two| {
        compare_locatables(
            header,
            (&one.contig, one.start as i32, one.stop as i32),
            (&two.contig, two.start as i32, two.stop as i32),
        )
        .then_with(|| inserted(one).cmp(&inserted(two)))
        .then_with(|| one.id.cmp(&two.id))
    });
}

/// The novel adjacencies and the chimeras behind each, as `groupByKey` leaves them cached.
pub struct Adjacencies(Vec<(NovelAdjacency, Vec<Chimera>)>);

/// The tool's first job: the reads to the novel adjacencies.
pub fn adjacencies(records: &[BamRecord], header: &SamHeader) -> Result<Adjacencies, Failed> {
    // `filter(!isSecondaryAlignment).groupBy(getName)`.
    let by_name = append_only_group(
        records
            .iter()
            .filter(|record| record.flags & 0x100 == 0)
            .map(|record| (record.read_name.clone(), record))
            .collect(),
        |name: &String| string_hash_code(name),
    );
    let mut contigs = Vec::with_capacity(by_name.len());
    for (_, reads) in &by_name {
        contigs.push(parse_contig(reads, header, true)?);
    }
    let contigs: Vec<&Contig> = contigs
        .iter()
        .filter(|contig| contig.alignments.len() > 1)
        .collect();
    adjacencies_of(&contigs, header)
}

/// `ContigChimericAlignmentIterativeInterpreter.discoverVariantsFromChimeras`, up to `groupByKey`.
fn adjacencies_of(contigs: &[&Contig], header: &SamHeader) -> Result<Adjacencies, Failed> {
    let mut found = Vec::new();
    for contig in contigs {
        for chimera in parse_one_contig(contig, header)? {
            let adjacency = novel_adjacency(&chimera, &contig.sequence, header)?;
            found.push((adjacency, chimera));
        }
    }
    Ok(Adjacencies(append_only_group(
        found,
        NovelAdjacency::hash_code,
    )))
}

impl Adjacencies {
    /// `readBreakpointsFromNarls`, whose `SVInterval`s are the only thing it can throw from.
    pub fn check_padded_breakpoints(&self, padding: i32) -> Result<(), Failed> {
        for (adjacency, _) in &self.0 {
            let padded = padding + adjacency.complication.homology.len() as i32;
            for location in [&adjacency.left, &adjacency.right] {
                check_padded_start(location.start - padded)?;
            }
        }
        Ok(())
    }

    /// The tool's second job, then `SVVCFWriter`'s sort.
    pub fn variants(
        &self,
        header: &SamHeader,
        reference: &mut dyn Reference,
        cnv_calls: Option<&[CnvCall]>,
    ) -> Result<Vec<VariantContext>, Failed> {
        let mut variants = Vec::with_capacity(self.0.len());
        for (adjacency, evidence) in &self.0 {
            variants.push(annotated_variant(
                adjacency, evidence, reference, header, cnv_calls,
            )?);
        }
        sort_variants(&mut variants, header);
        Ok(variants)
    }
}

/// `SVInterval`'s constructor, on the start a padding moved.
pub fn check_padded_start(start: i32) -> Result<(), Failed> {
    if start < 0 {
        return Err(format!("provided start is negative: {start}"));
    }
    Ok(())
}

/// `GATKSVVCFHeaderLines` and the `END` line `SVVCFWriter.getVcfHeader` adds.
pub fn header_lines() -> Vec<htsjdk_vcf::header::HeaderLine> {
    use htsjdk_vcf::header::{Cardinality as C, HeaderLine as H, LineType as T};
    let alt = |id: &str, description: &str| H::Structured {
        key: "ALT".to_string(),
        fields: vec![
            ("ID".to_string(), id.to_string()),
            ("Description".to_string(), description.to_string()),
        ],
    };
    vec![
        alt("INV", "Inversion of reference sequence"),
        alt("DEL", "Deletion relative to the reference"),
        alt(
            "INS",
            "Insertion of novel sequence relative to the reference",
        ),
        alt(
            "DUP",
            "Region of elevated copy number relative to the reference",
        ),
        alt(
            "DUP:INV",
            "Region of elevated copy number relative to the reference, with some copies inverted",
        ),
        alt("CPX", "Complex rearrangement of reference sequence"),
        H::info(
            "END",
            C::Fixed(1),
            T::Integer,
            "Stop position of the interval",
        ),
        H::info(
            "SVTYPE",
            C::Fixed(1),
            T::String,
            "Type of structural variant",
        ),
        H::info(
            "SVLEN",
            C::Unbounded,
            T::Integer,
            "Difference in length between REF and ALT alleles",
        ),
        H::info(
            "SEQ_ALT_HAPLOTYPE",
            C::A,
            T::Character,
            "Alt haplotype sequence, one per alt allele",
        ),
        H::info(
            "INSSEQ",
            C::Unbounded,
            T::String,
            "Inserted sequence at the breakpoint",
        ),
        H::info(
            "INSLEN",
            C::A,
            T::Integer,
            "Length of inserted sequence (note for duplication records, this does not count the \
             extra copies of the duplicated sequence)",
        ),
        H::info(
            "INSSEQ_MAP",
            C::Unbounded,
            T::String,
            "Alignments of inserted sequence",
        ),
        H::info(
            "HOMSEQ",
            C::Unbounded,
            T::String,
            "Sequence of base pair identical micro-homology at event breakpoints",
        ),
        H::info(
            "HOMLEN",
            C::Fixed(1),
            T::Integer,
            "Length of base pair identical micro-homology at event breakpoints",
        ),
        H::info(
            "READ_PAIR_SUPPORT",
            C::Fixed(1),
            T::Integer,
            "Number of discordant read pairs supporting the variant",
        ),
        H::info(
            "SPLIT_READ_SUPPORT",
            C::Fixed(1),
            T::Integer,
            "Number of split read supplementary mappings supporting the variant",
        ),
        H::info(
            "EXTERNAL_CNV_CALLS",
            C::Fixed(1),
            T::String,
            "Comma-delimited list of external copy number calls that overlap with this variant in \
             format ID:CN:CNQ",
        ),
        H::info(
            "CTG_NAMES",
            C::Unbounded,
            T::String,
            "Name of local assembly contigs supporting this variant, formatted as \
             \"asm%06d:tig%05d\"",
        ),
        H::info(
            "TOTAL_MAPPINGS",
            C::Fixed(1),
            T::Integer,
            "Number of local assembly contigs supporting the variant, i.e. number of entries in \
             CTG_NAMES",
        ),
        H::info(
            "MAPPING_QUALITIES",
            C::Unbounded,
            T::Integer,
            "Mapping qualities of the contig alignments that support the variant",
        ),
        H::info(
            "HQ_MAPPINGS",
            C::Fixed(1),
            T::Integer,
            "Number of high-quality contig alignments that support the variant",
        ),
        H::info(
            "ALIGN_LENGTHS",
            C::Unbounded,
            T::Integer,
            "Minimum lengths of the flanking aligned region from each contig alignment",
        ),
        H::info(
            "MAX_ALIGN_LENGTH",
            C::Fixed(1),
            T::Integer,
            "Maximum of the values listed in ALIGN_LENGTHS",
        ),
        H::info(
            "CTG_GOOD_NONCANONICAL_MAPPING",
            C::Unbounded,
            T::String,
            "Good mapping of evidence contigs as listed in CTG_NAMES to non-canonical reference \
             chromosomes that could potentially offer better explanation of the assembly contig \
             without the SV record. One for each evidence assembly contig, if available; \
             otherwise a \".\". If no evidence contig has such mapping, this annotation is \
             omitted for the record.",
        ),
        H::info(
            "LINK",
            C::Unbounded,
            T::String,
            "ID(s) of other variants that are linked to this record, i.e. they constitute a \
             larger more complex variant",
        ),
        H::info(
            "MATEID",
            C::Unbounded,
            T::String,
            "ID(s) for mate(s) of a BND record",
        ),
        H::info(
            "IMPRECISE",
            C::Fixed(0),
            T::Flag,
            "Imprecise structural variation",
        ),
        H::info(
            "CIPOS",
            C::Fixed(2),
            T::Integer,
            "Confidence interval around POS for imprecise variants",
        ),
        H::info(
            "CIEND",
            C::Fixed(2),
            T::Integer,
            "Confidence interval around END for imprecise variants",
        ),
        H::info(
            "INV33",
            C::Fixed(0),
            T::Flag,
            "Whether the event represents a 3' to 5' breakpoint",
        ),
        H::info(
            "INV55",
            C::Fixed(0),
            T::Flag,
            "Whether the event represents a 5' to 3' breakpoint",
        ),
        H::info(
            "DUP_REPEAT_UNIT_REF_SPAN",
            C::Fixed(1),
            T::String,
            "Reference span of the suspected repeated unit in a tandem duplication",
        ),
        H::info(
            "DUP_SEQ_CIGARS",
            C::Unbounded,
            T::String,
            "CIGARs of the repeated sequence on the locally-assembled contigs when aligned to \
             DUP_REPEAT_UNIT_REF_SPAN (currently only available for repeats when \
             DUP_ANNOTATIONS_IMPRECISE is false)",
        ),
        H::info(
            "DUP_NUM",
            C::R,
            T::Integer,
            "Number of times the sequence is duplicated on reference and on the alternate alleles",
        ),
        H::info(
            "DUP_ANNOTATIONS_IMPRECISE",
            C::Fixed(0),
            T::Flag,
            "Whether the duplication annotations are from an experimental optimization procedure",
        ),
        H::info(
            "DUP_IMPRECISE_AFFECTED_RANGE",
            C::Unbounded,
            T::String,
            "Affected reference range for duplications annotated with the flag \
             DUP_ANNOTATIONS_IMPRECISE",
        ),
        H::info(
            "DUP_ORIENTATIONS",
            C::A,
            T::String,
            "Orientations of the duplicated sequence on alt allele relative to the copy on ref; \
             one group for each alt allele (currently only available for inverted duplication \
             variants)",
        ),
        H::info(
            "CONTRACTION",
            C::Fixed(0),
            T::Flag,
            "Tandem repeats contraction compared to reference",
        ),
        H::info(
            "EXPANSION",
            C::Fixed(0),
            T::Flag,
            "Tandem repeats expansion compared to reference",
        ),
        H::info(
            "ALT_ARRANGEMENT",
            C::Unbounded,
            T::String,
            "For CPX variants only; specifies how reference segments given in SEGMENTS are \
             re-arranged",
        ),
        H::info(
            "SEGMENTS",
            C::Unbounded,
            T::String,
            "For CPX variants only; segments of reference that are rearranged",
        ),
        H::info(
            "CPX_EVENT",
            C::Unbounded,
            T::String,
            "ID(s) of CPX(s) events from which current simple variant record is extracted",
        ),
        H::format(
            "CN",
            C::Fixed(1),
            T::Integer,
            "Copy number genotype for imprecise events",
        ),
        H::format(
            "CNQ",
            C::Fixed(1),
            T::Float,
            "Copy number genotype quality for imprecise events",
        ),
        H::filter(
            "LOW_MQ",
            "Assembly evidence based record that whose maximum value specified in \
             MAPPING_QUALITIES is lower than user specified threshold",
        ),
        H::filter(
            "SHORT_ALN",
            "Assembly evidence based record that whose MAPPING_QUALITIES value is lower than \
             user specified threshold",
        ),
        H::filter(
            "LOW_QS",
            "Depth-only copy number record whose QS value is lower than user specified threshold",
        ),
        H::filter(
            "FREQ",
            "Depth-only copy number record whose AF value is higher than user specified threshold",
        ),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repeat_structure_stops_at_the_first_exact_fit() {
        // `new TandemRepeatStructure(-10, -5)`: expansion, two copies to three of a unit of five.
        let repeat = tandem_repeat(-10, -5);
        assert_eq!(
            (
                repeat.lower,
                repeat.higher,
                repeat.unit_length,
                repeat.pseudo_homology
            ),
            (2, 3, 5, 0)
        );
    }

    #[test]
    fn a_group_keeps_its_slot_order_and_first_key() {
        let groups = append_only_group(
            vec![
                ("a".to_string(), 1),
                ("b".to_string(), 2),
                ("a".to_string(), 3),
            ],
            |key: &String| string_hash_code(key),
        );
        assert_eq!(groups.iter().map(|group| group.1.len()).sum::<usize>(), 3);
    }
}
