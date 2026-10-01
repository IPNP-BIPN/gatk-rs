//! The alignments of one locally assembled contig, and the configuration of them that
//! `StructuralVariantDiscoverer` believes.
//!
//! An assembled contig reaches the tool as several SAM lines: a primary, its supplementaries and,
//! when filters allow, secondaries. Each becomes an [`AlignmentInterval`] expressed along the
//! contig's 5' to 3' direction. The contig then keeps the best-scoring subset of those
//! alignments, splits any alignment carrying a gap of fifty bases or more into its pieces, and is
//! classified by what is left: one alignment, two (a simple chimera), or more (a complex event).
//!
//! Ported from `org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig`,
//! `org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments`
//! and the parts of `org.broadinstitute.hellbender.utils.read.CigarUtils` they call.

use std::cmp::Ordering;
use std::fmt;

/// `java.lang.IllegalArgumentException`, which `Utils.validateArg` throws.
pub const ILLEGAL_ARGUMENT: &str = "java.lang.IllegalArgumentException";
/// `java.lang.IllegalStateException`, which `Utils.validate` throws.
pub const ILLEGAL_STATE: &str = "java.lang.IllegalStateException";
/// `org.broadinstitute.hellbender.exceptions.GATKException`.
pub const GATK_EXCEPTION: &str = "org.broadinstitute.hellbender.exceptions.GATKException";
/// `GATKException.ShouldNeverReachHereException`, a nested class and printed as one.
pub const SHOULD_NEVER_REACH_HERE: &str =
    "org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException";
/// `java.lang.UnsupportedOperationException`.
pub const UNSUPPORTED_OPERATION: &str = "java.lang.UnsupportedOperationException";
/// `org.broadinstitute.hellbender.exceptions.UserException`, the one refusal that is the user's.
pub const USER_EXCEPTION: &str = "org.broadinstitute.hellbender.exceptions.UserException";

/// What the interpretation threw: the exception's class and its message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SvError {
    pub class: &'static str,
    pub message: String,
}

impl SvError {
    pub fn new(class: &'static str, message: impl Into<String>) -> SvError {
        SvError {
            class,
            message: message.into(),
        }
    }

    /// Whether `mainEntry` reports it through the user-error handler.
    pub fn is_user(&self) -> bool {
        self.class == USER_EXCEPTION
    }
}

/// `CigarOperator`, with the five predicates the SV code asks of it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    M,
    I,
    D,
    N,
    S,
    H,
    P,
    Eq,
    X,
}

impl CigarOp {
    pub fn from_char(character: char) -> Option<CigarOp> {
        Some(match character {
            'M' => CigarOp::M,
            'I' => CigarOp::I,
            'D' => CigarOp::D,
            'N' => CigarOp::N,
            'S' => CigarOp::S,
            'H' => CigarOp::H,
            'P' => CigarOp::P,
            '=' => CigarOp::Eq,
            'X' => CigarOp::X,
            _ => return None,
        })
    }

    pub fn to_char(self) -> char {
        match self {
            CigarOp::M => 'M',
            CigarOp::I => 'I',
            CigarOp::D => 'D',
            CigarOp::N => 'N',
            CigarOp::S => 'S',
            CigarOp::H => 'H',
            CigarOp::P => 'P',
            CigarOp::Eq => '=',
            CigarOp::X => 'X',
        }
    }

    pub fn consumes_read_bases(self) -> bool {
        matches!(
            self,
            CigarOp::M | CigarOp::I | CigarOp::S | CigarOp::Eq | CigarOp::X
        )
    }

    pub fn consumes_reference_bases(self) -> bool {
        matches!(
            self,
            CigarOp::M | CigarOp::D | CigarOp::N | CigarOp::Eq | CigarOp::X
        )
    }

    pub fn is_clipping(self) -> bool {
        matches!(self, CigarOp::S | CigarOp::H)
    }

    pub fn is_indel(self) -> bool {
        matches!(self, CigarOp::I | CigarOp::D)
    }

    pub fn is_alignment(self) -> bool {
        matches!(self, CigarOp::M | CigarOp::Eq | CigarOp::X)
    }
}

/// `htsjdk.samtools.Cigar`: its elements as (length, operator), in order.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Cigar(pub Vec<(i32, CigarOp)>);

impl Cigar {
    /// `TextCigarCodec.decode`, with `*` for the empty cigar.
    pub fn parse(text: &str) -> Option<Cigar> {
        if text == "*" {
            return Some(Cigar::default());
        }
        let mut elements = Vec::new();
        let mut length = String::new();
        for character in text.chars() {
            if character.is_ascii_digit() {
                length.push(character);
            } else {
                let op = CigarOp::from_char(character)?;
                elements.push((length.parse().ok()?, op));
                length.clear();
            }
        }
        length.is_empty().then_some(Cigar(elements))
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// `Cigar.getReadLength`: the bases the operators consume from the read.
    pub fn read_length(&self) -> i32 {
        self.0
            .iter()
            .filter(|(_, op)| op.consumes_read_bases())
            .map(|(length, _)| length)
            .sum()
    }

    /// `Cigar.getReferenceLength`.
    pub fn reference_length(&self) -> i32 {
        self.0
            .iter()
            .filter(|(_, op)| op.consumes_reference_bases())
            .map(|(length, _)| length)
            .sum()
    }

    /// `CigarUtils.invertCigar`: the elements in reverse order.
    pub fn inverted(&self) -> Cigar {
        Cigar(self.0.iter().rev().copied().collect())
    }

    /// `CigarUtils.countUnclippedReadBases`: read bases plus every clip, hard ones included.
    pub fn count_unclipped_read_bases(&self) -> i32 {
        self.0
            .iter()
            .filter(|(_, op)| op.is_clipping() || op.consumes_read_bases())
            .map(|(length, _)| length)
            .sum()
    }

    /// `CigarUtils.countClippedBases(cigar, tail, type)` for one clipping operator.
    fn count_clipped_of(&self, left: bool, kind: CigarOp) -> Result<i32, SvError> {
        let size = self.0.len();
        if size < 2 {
            if size == 1 && !self.0[0].1.is_clipping() {
                return Ok(0);
            }
            return Err(SvError::new(
                ILLEGAL_ARGUMENT,
                "cigar is empty or completely clipped.",
            ));
        }
        let mut result = 0;
        for n in 0..size {
            let (length, op) = self.0[if left { n } else { size - n - 1 }];
            if !op.is_clipping() {
                return Ok(result);
            } else if op == kind {
                result += length;
            }
        }
        Err(SvError::new(
            ILLEGAL_ARGUMENT,
            format!("Input cigar {self} is completely clipped."),
        ))
    }

    /// `CigarUtils.countClippedBases(cigar, tail)`: soft and hard clips at one end.
    pub fn count_clipped_bases(&self, left: bool) -> Result<i32, SvError> {
        Ok(self.count_clipped_of(left, CigarOp::S)? + self.count_clipped_of(left, CigarOp::H)?)
    }

    /// `CigarUtils.computeAssociatedDistOnRead`: how many read bases, walking from `start`,
    /// cover `ref_dist` reference bases.
    pub fn associated_distance_on_read(
        &self,
        start: i32,
        ref_dist: i32,
        backward: bool,
    ) -> Result<i32, SvError> {
        if !(ref_dist > 0 && start > 0) {
            return Err(SvError::new(
                ILLEGAL_ARGUMENT,
                format!("start {start} or distance {ref_dist} is non-positive."),
            ));
        }
        let elements: Vec<(i32, CigarOp)> = if backward {
            self.0.iter().rev().copied().collect()
        } else {
            self.0.clone()
        };
        let read_length: i32 = elements
            .iter()
            .filter(|(_, op)| op.consumes_read_bases())
            .map(|(length, _)| length)
            .sum();
        let skip = if backward {
            read_length - start
        } else {
            start - 1
        };
        let mut read_consumed = 0;
        let mut ref_consumed = 0;
        for (length, op) in elements {
            let before = read_consumed;
            if op.consumes_read_bases() {
                read_consumed += length;
            }
            if read_consumed <= skip {
                continue;
            }
            if op.consumes_reference_bases() {
                ref_consumed += length - (skip - before).max(0);
            }
            if ref_consumed >= ref_dist {
                let excess = (ref_consumed - ref_dist).max(0);
                return Ok(read_consumed
                    - skip
                    - if op.consumes_read_bases() { excess } else { 0 });
            }
        }
        Err(SvError::new(
            ILLEGAL_ARGUMENT,
            format!(
                "Cigar {self}does not contain at least {ref_dist} reference bases past red start {start}."
            ),
        ))
    }
}

impl fmt::Display for Cigar {
    /// `TextCigarCodec.encode`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            return write!(f, "*");
        }
        for (length, op) in &self.0 {
            write!(f, "{length}{}", op.to_char())?;
        }
        Ok(())
    }
}

/// `SimpleInterval`: a contig and a 1-based, closed span.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Interval {
    /// The constructor's own validation, which is what a window past a contig's end meets.
    pub fn new(contig: &str, start: i32, end: i32) -> Result<Interval, SvError> {
        if start > 0 && end >= start {
            Ok(Interval::of(contig, start, end))
        } else {
            Err(SvError::new(
                ILLEGAL_ARGUMENT,
                format!("Invalid interval. Contig:{contig} start:{start} end:{end}"),
            ))
        }
    }

    /// An interval the caller knows to be valid.
    pub fn of(contig: &str, start: i32, end: i32) -> Interval {
        Interval {
            contig: contig.to_string(),
            start,
            end,
        }
    }

    pub fn size(&self) -> i32 {
        self.end - self.start + 1
    }

    /// `overlaps`, which is `overlapsWithMargin(other, 0)`.
    pub fn overlaps(&self, other: &Interval) -> bool {
        self.contig == other.contig && self.start <= other.end && other.start <= self.end
    }

    pub fn contains(&self, other: &Interval) -> bool {
        self.contig == other.contig && self.start <= other.start && self.end >= other.end
    }
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start, self.end)
    }
}

/// The sequence dictionary the reads carry, which orders contigs and bounds breakpoints.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Dictionary {
    pub names: Vec<String>,
    pub lengths: Vec<i32>,
    /// Each record's `AS`, which a VCF contig line carries as `assembly`.
    pub assemblies: Vec<Option<String>>,
}

impl Dictionary {
    /// `getSequenceIndex`, -1 where the contig is absent.
    pub fn index(&self, contig: &str) -> i32 {
        self.names
            .iter()
            .position(|name| name == contig)
            .map_or(-1, |index| index as i32)
    }

    pub fn length(&self, contig: &str) -> Option<i32> {
        let index = self.index(contig);
        (index >= 0).then(|| self.lengths[index as usize])
    }

    /// `IntervalUtils.compareContigs`.
    pub fn compare_contigs(&self, first: &str, second: &str) -> Result<Ordering, SvError> {
        let (one, two) = (self.index(first), self.index(second));
        if one == -1 || two == -1 {
            return Err(SvError::new(
                ILLEGAL_ARGUMENT,
                "Can't do comparison because Locatables' contigs not found in sequence dictionary",
            ));
        }
        Ok(one.cmp(&two))
    }

    /// `IntervalUtils.compareLocatables`: contig, then start, then end.
    pub fn compare_locatables(
        &self,
        first: &Interval,
        second: &Interval,
    ) -> Result<Ordering, SvError> {
        Ok(self
            .compare_contigs(&first.contig, &second.contig)?
            .then(first.start.cmp(&second.start))
            .then(first.end.cmp(&second.end)))
    }
}

/// One SAM line of an assembled contig, as `GATKRead` exposes it to the tool.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigRead {
    pub name: String,
    /// `GATKRead.isUnmapped`: the flag, or no contig, or no start.
    pub unmapped: bool,
    pub reverse_strand: bool,
    pub supplementary: bool,
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub cigar: Cigar,
    pub mapping_quality: i32,
    pub bases: Vec<u8>,
    pub nm: Option<i32>,
    pub alignment_score: Option<i32>,
}

impl ContigRead {
    /// `GATKRead.commonToString`, which is what an exception about a read prints.
    pub fn common_to_string(&self) -> String {
        if self.unmapped || self.cigar.is_empty() {
            format!("{} UNMAPPED", self.name)
        } else {
            format!("{} {}:{}-{}", self.name, self.contig, self.start, self.end)
        }
    }
}

/// `ContigAlignmentsModifier.AlnModType`, printed by its one-letter code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlnModType {
    None,
    UndergoneOverlapRemoval,
    ExtractedFromLargerAlignment,
    FromSplitGappedAlignment,
}

impl AlnModType {
    fn code(self) -> &'static str {
        match self {
            AlnModType::None => "O",
            AlnModType::UndergoneOverlapRemoval => "H",
            AlnModType::ExtractedFromLargerAlignment => "E",
            AlnModType::FromSplitGappedAlignment => "S",
        }
    }
}

/// `AlignmentInterval.NO_NM` and `NO_AS`.
pub const NO_NM: i32 = -1;
pub const NO_AS: i32 = -1;

/// `AlignmentInterval`: one alignment of a contig, with its cigar along the contig's 5' to 3'
/// direction, so a reverse-strand line's cigar is the SAM one inverted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentInterval {
    pub reference_span: Interval,
    pub start_in_contig: i32,
    pub end_in_contig: i32,
    pub cigar: Cigar,
    pub forward_strand: bool,
    pub map_qual: i32,
    pub mismatches: i32,
    pub aln_score: i32,
    pub mod_type: AlnModType,
}

impl AlignmentInterval {
    /// `new AlignmentInterval(GATKRead)`, which refuses an unmapped line.
    pub fn from_read(read: &ContigRead) -> Result<AlignmentInterval, SvError> {
        if read.unmapped {
            return Err(SvError::new(
                ILLEGAL_ARGUMENT,
                format!(
                    "read being used to construct AlignmentInterval is unmapped: {}",
                    read.common_to_string()
                ),
            ));
        }
        let cigar = if read.reverse_strand {
            read.cigar.inverted()
        } else {
            read.cigar.clone()
        };
        let start_in_contig = 1 + cigar.count_clipped_bases(true)?;
        let end_in_contig =
            cigar.count_unclipped_read_bases() - cigar.count_clipped_bases(false)?;
        Ok(AlignmentInterval {
            reference_span: Interval::new(&read.contig, read.start, read.end)?,
            start_in_contig,
            end_in_contig,
            cigar,
            forward_strand: !read.reverse_strand,
            map_qual: read.mapping_quality,
            mismatches: read.nm.unwrap_or(NO_NM),
            aln_score: read.alignment_score.unwrap_or(NO_AS),
            mod_type: AlnModType::None,
        })
    }

    /// The full constructor, which validates the cigar against both spans.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        reference_span: Interval,
        start_in_contig: i32,
        end_in_contig: i32,
        cigar: Cigar,
        forward_strand: bool,
        map_qual: i32,
        mismatches: i32,
        aln_score: i32,
        mod_type: AlnModType,
    ) -> Result<AlignmentInterval, SvError> {
        // `checkValidArgument`: soft clips are counted after a terminal insertion becomes one.
        let count = cigar.0.len();
        let soft_clipped: i32 = cigar
            .0
            .iter()
            .enumerate()
            .filter(|(n, (_, op))| {
                *op == CigarOp::S
                    || (*op == CigarOp::I
                        && count >= 2
                        && (*n == 0
                            || *n == count - 1
                            || cigar.0[n - 1].1.is_clipping()
                            || cigar.0[n + 1].1.is_clipping()))
            })
            .map(|(_, (length, _))| length)
            .sum();
        let read_length = cigar.read_length() - soft_clipped;
        if cigar.reference_length() != reference_span.size()
            || read_length != end_in_contig - start_in_contig + 1
        {
            return Err(SvError::new(
                ILLEGAL_ARGUMENT,
                format!(
                    "Encountering invalid arguments for constructing alignment,\tcigar: {cigar} ref.span: {reference_span} read span: {start_in_contig}-{end_in_contig}"
                ),
            ));
        }
        Ok(AlignmentInterval {
            reference_span,
            start_in_contig,
            end_in_contig,
            cigar,
            forward_strand,
            map_qual,
            mismatches,
            aln_score,
            mod_type,
        })
    }

    pub fn contains_gap_of_equal_or_larger_size(&self, gap: i32) -> bool {
        self.cigar
            .0
            .iter()
            .any(|(length, op)| op.is_indel() && *length >= gap)
    }

    pub fn size_on_read(&self) -> i32 {
        self.end_in_contig - self.start_in_contig + 1
    }

    pub fn contains_on_ref(&self, other: &AlignmentInterval) -> bool {
        self.reference_span.contains(&other.reference_span)
    }

    pub fn contains_on_read(&self, other: &AlignmentInterval) -> bool {
        self.start_in_contig <= other.start_in_contig && self.end_in_contig >= other.end_in_contig
    }

    /// The cigar along the reference, which is the SAM one.
    pub fn cigar_along_reference(&self) -> Cigar {
        if self.forward_strand {
            self.cigar.clone()
        } else {
            self.cigar.inverted()
        }
    }

    /// `toSATagString`: `contig,start,strand,cigar,mq[,nm[,as]]`.
    pub fn to_sa_tag_string(&self) -> String {
        let mut text = format!(
            "{},{},{},{},{}",
            self.reference_span.contig,
            self.reference_span.start,
            if self.forward_strand { "+" } else { "-" },
            self.cigar_along_reference(),
            self.map_qual
        );
        if self.mismatches != NO_NM || self.aln_score != NO_AS {
            text.push(',');
            if self.mismatches == NO_NM {
                text.push_str(&format!(".,{}", self.aln_score));
            } else if self.aln_score == NO_AS {
                text.push_str(&self.mismatches.to_string());
            } else {
                text.push_str(&format!("{},{}", self.mismatches, self.aln_score));
            }
        }
        text
    }

    /// `toPackedString`, which is how an unused alignment is annotated.
    pub fn to_packed_string(&self) -> String {
        format!(
            "{}_{}_{}_{}_{}_{}_{}_{}_{}",
            self.start_in_contig,
            self.end_in_contig,
            self.reference_span,
            if self.forward_strand { "+" } else { "-" },
            self.cigar,
            self.map_qual,
            self.mismatches,
            self.aln_score,
            self.mod_type.code()
        )
    }
}

/// `AlignmentInterval.overlapOnContig`.
pub fn overlap_on_contig(one: &AlignmentInterval, two: &AlignmentInterval) -> i32 {
    ((one.end_in_contig + 1).min(two.end_in_contig + 1)
        - one.start_in_contig.max(two.start_in_contig))
    .max(0)
}

/// `AlignmentInterval.overlapOnRefSpan`, over half-open intervals.
pub fn overlap_on_ref_span(one: &AlignmentInterval, two: &AlignmentInterval) -> i32 {
    if one.reference_span.contig != two.reference_span.contig {
        return 0;
    }
    ((one.reference_span.end + 1).min(two.reference_span.end + 1)
        - one.reference_span.start.max(two.reference_span.start))
    .max(0)
}

/// `StructuralVariationDiscoveryArgumentCollection.STRUCTURAL_VARIANT_SIZE_LOWER_BOUND`, which is
/// also the gap an alignment is split at.
pub const STRUCTURAL_VARIANT_SIZE_LOWER_BOUND: i32 = 50;

/// `ContigAlignmentsModifier.splitGappedAlignment`: the pieces either side of every indel of at
/// least `sensitivity` bases, or the alignment itself when there are fewer than two.
pub fn split_gapped_alignment(
    one: &AlignmentInterval,
    sensitivity: i32,
    unclipped_contig_len: i32,
) -> Result<Vec<AlignmentInterval>, SvError> {
    let elements = &one.cigar.0;
    let mut count = elements.len();
    if count <= 1 {
        return Ok(vec![one.clone()]);
    }
    let mut contig_offset = 0;
    let mut index = 0;
    if elements[0].1 == CigarOp::H {
        contig_offset += elements[0].0;
        index += 1;
    }
    if elements[count - 1].1 == CigarOp::H {
        count -= 1;
    }
    // A reverse-strand alignment walks the reference backwards, which the negated offset keeps
    // increasing.
    let mut ref_offset = if one.forward_strand {
        one.reference_span.start
    } else {
        -one.reference_span.end
    };
    const NOT_SET: i32 = -1;
    let (mut start_contig, mut start_index, mut start_ref) = (NOT_SET, 0usize, 0);
    let (mut end_contig, mut end_index, mut end_ref) = (0, 0usize, 0);
    let mut result = Vec::new();
    let grab = |start_contig: i32,
                start_index: usize,
                start_ref: i32,
                end_contig: i32,
                end_index: usize,
                end_ref: i32|
     -> Result<AlignmentInterval, SvError> {
        let mut pieces = Vec::new();
        let mut initial_soft = start_contig;
        if elements[0].1 == CigarOp::H {
            pieces.push(elements[0]);
            initial_soft -= elements[0].0;
        }
        if initial_soft > 0 {
            pieces.push((initial_soft, CigarOp::S));
        }
        pieces.extend_from_slice(&elements[start_index..end_index]);
        let mut final_soft = unclipped_contig_len - end_contig;
        let last = elements[elements.len() - 1];
        if last.1 == CigarOp::H {
            final_soft -= last.0;
            if final_soft > 0 {
                pieces.push((final_soft, CigarOp::S));
            }
            pieces.push(last);
        } else if final_soft > 0 {
            pieces.push((final_soft, CigarOp::S));
        }
        let contig = &one.reference_span.contig;
        let span = if one.forward_strand {
            Interval::new(contig, start_ref, end_ref - 1)?
        } else {
            Interval::new(contig, -end_ref + 1, -start_ref)?
        };
        AlignmentInterval::build(
            span,
            start_contig + 1,
            end_contig,
            Cigar(pieces),
            one.forward_strand,
            one.map_qual,
            NO_NM,
            NO_AS,
            AlnModType::FromSplitGappedAlignment,
        )
    };
    while index < count {
        let (length, op) = elements[index];
        if op.is_alignment() {
            if start_contig == NOT_SET {
                start_contig = contig_offset;
                start_index = index;
                start_ref = ref_offset;
            }
            end_contig = contig_offset + length;
            end_index = index + 1;
            end_ref = ref_offset + length;
        } else if op.is_indel() && length >= sensitivity && start_contig != NOT_SET {
            result.push(grab(
                start_contig,
                start_index,
                start_ref,
                end_contig,
                end_index,
                end_ref,
            )?);
            start_contig = NOT_SET;
        }
        if op.consumes_read_bases() {
            contig_offset += length;
        }
        if op.consumes_reference_bases() {
            ref_offset += length;
        }
        index += 1;
    }
    if start_contig != NOT_SET {
        result.push(grab(
            start_contig,
            start_index,
            start_ref,
            end_contig,
            end_index,
            end_ref,
        )?);
    }
    if result.len() < 2 {
        return Ok(vec![one.clone()]);
    }
    Ok(result)
}

/// `AlignedContig.ALIGNMENT_MQ_THRESHOLD`: the mapping quality an alignment must beat to count.
const ALIGNMENT_MQ_THRESHOLD: i32 = 20;
/// `AlignedContig.ALIGNMENT_LOW_READ_UNIQUENESS_THRESHOLD`.
const ALIGNMENT_LOW_READ_UNIQUENESS_THRESHOLD: i32 = 10;
/// `AlignedContig.COVERAGE_MQ_NORMALIZATION_CONST`.
const COVERAGE_MQ_NORMALIZATION_CONST: f64 = 60.0;
/// `AlignedContig.ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST`.
const ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST: i32 = 10;
/// `AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME`.
pub const NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME: &str = "NONE";

/// `AlignedContig.getAlignmentIntervalComparator`: start on the contig, then reference contig
/// name, then reference start.
pub fn alignment_order(one: &AlignmentInterval, two: &AlignmentInterval) -> Ordering {
    one.start_in_contig
        .cmp(&two.start_in_contig)
        .then_with(|| one.reference_span.contig.cmp(&two.reference_span.contig))
        .then(one.reference_span.start.cmp(&two.reference_span.start))
}

/// `AlignedContig`: a contig's name, its sequence as assembled, and its alignments in order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignedContig {
    pub name: String,
    pub sequence: Vec<u8>,
    pub alignments: Vec<AlignmentInterval>,
}

/// `AlignedContig.GoodAndBadMappings`.
#[derive(Debug, Clone, PartialEq, Eq)]
struct GoodAndBadMappings {
    good: Vec<AlignmentInterval>,
    bad: Vec<AlignmentInterval>,
    non_canonical: Option<AlignmentInterval>,
}

/// `List.removeAll`: every element equal to one of `gone` leaves.
fn remove_all(list: &mut Vec<AlignmentInterval>, gone: &[AlignmentInterval]) {
    list.retain(|alignment| !gone.contains(alignment));
}

/// `Math.ulp(double)`: the distance to the next double of larger magnitude.
fn ulp(value: f64) -> f64 {
    let magnitude = value.abs();
    if magnitude == f64::MAX {
        return 2f64.powi(971);
    }
    f64::from_bits(magnitude.to_bits() + 1) - magnitude
}

impl AlignedContig {
    /// The constructor, which sorts the alignments.
    pub fn new(
        name: &str,
        sequence: Vec<u8>,
        mut alignments: Vec<AlignmentInterval>,
    ) -> AlignedContig {
        alignments.sort_by(alignment_order);
        AlignedContig {
            name: name.to_string(),
            sequence,
            alignments,
        }
    }

    /// `hasGoodMQ`: one alignment above MQ 20, or two, or one carrying a large gap.
    pub fn has_good_mq(&self) -> bool {
        if self.alignments.len() < 2 {
            return self
                .alignments
                .first()
                .is_some_and(|one| one.map_qual > ALIGNMENT_MQ_THRESHOLD);
        }
        let mut not_bad = 0;
        for alignment in &self.alignments {
            if alignment.map_qual > ALIGNMENT_MQ_THRESHOLD {
                if alignment
                    .contains_gap_of_equal_or_larger_size(STRUCTURAL_VARIANT_SIZE_LOWER_BOUND)
                {
                    return true;
                }
                not_bad += 1;
            }
        }
        not_bad > 1
    }

    /// `reconstructContigFromBestConfiguration`.
    pub fn reconstruct_from_best_configuration(
        &self,
        canonical: &[String],
        score_diff_tolerance: f64,
    ) -> Result<Vec<FineTunedContig>, SvError> {
        let best = self.pick_and_filter_configurations(canonical, score_diff_tolerance)?;
        if best.len() > 1 {
            // Every one of these is marked ambiguous, which no later step turns into a call.
            let mut result = Vec::new();
            for mappings in best {
                let mappings = remove_non_unique_mappings(split_gaps(mappings)?);
                if not_stitchable(&mappings.good)? {
                    result.push(self.fine_tuned(mappings, true));
                }
            }
            result.sort_by(|one, two| {
                one.contig
                    .alignments
                    .len()
                    .cmp(&two.contig.alignments.len())
                    .then_with(|| {
                        let sum = |tig: &FineTunedContig| -> i32 {
                            tig.contig.alignments.iter().map(|ai| ai.mismatches).sum()
                        };
                        sum(one).cmp(&sum(two))
                    })
            });
            Ok(result)
        } else {
            let first = best.into_iter().next().ok_or_else(|| {
                SvError::new(ILLEGAL_STATE, "no configuration picked for a contig")
            })?;
            let result = remove_non_unique_mappings(split_gaps(first)?);
            if not_stitchable(&result.good)? {
                Ok(vec![self.fine_tuned(result, false)])
            } else {
                Ok(Vec::new())
            }
        }
    }

    fn fine_tuned(&self, mappings: GoodAndBadMappings, ambiguous: bool) -> FineTunedContig {
        FineTunedContig {
            contig: AlignedContig::new(&self.name, self.sequence.clone(), mappings.good),
            insertion_mappings: mappings
                .bad
                .iter()
                .map(AlignmentInterval::to_packed_string)
                .collect(),
            ambiguous,
            non_canonical_sa_tag: mappings.non_canonical.as_ref().map_or_else(
                || NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME.to_string(),
                AlignmentInterval::to_sa_tag_string,
            ),
        }
    }

    fn pick_and_filter_configurations(
        &self,
        canonical: &[String],
        tolerance: f64,
    ) -> Result<Vec<GoodAndBadMappings>, SvError> {
        let picked = self.pick_best_configurations(canonical, tolerance)?;
        // `filterSecondaryConfigurationsByMappingQualityThreshold` at zero.
        if picked.len() == 1 {
            return Ok(picked);
        }
        let above: Vec<GoodAndBadMappings> = picked
            .iter()
            .filter(|mappings| {
                mappings
                    .good
                    .iter()
                    .map(|ai| ai.map_qual)
                    .min()
                    .unwrap_or(0)
                    > 0
            })
            .cloned()
            .collect();
        Ok(if above.len() != 1 { picked } else { above })
    }

    /// `pickBestConfigurations`: every subset of the good alignments scored, and the best kept.
    fn pick_best_configurations(
        &self,
        canonical: &[String],
        tolerance: f64,
    ) -> Result<Vec<GoodAndBadMappings>, SvError> {
        if self.alignments.len() == 1 {
            return Ok(vec![GoodAndBadMappings {
                good: vec![self.alignments[0].clone()],
                bad: Vec::new(),
                non_canonical: None,
            }]);
        }
        let is_canonical = |ai: &AlignmentInterval| canonical.contains(&ai.reference_span.contig);
        let max_canonical_score = |list: &[AlignmentInterval]| {
            list.iter()
                .filter(|ai| is_canonical(ai))
                .map(|ai| ai.aln_score)
                .max()
                .unwrap_or(0)
        };
        let first_max = max_canonical_score(&self.alignments);
        // `heuristicSpeedUpWhenFacingManyMappings`.
        let (mut good, bad): (Vec<AlignmentInterval>, Vec<AlignmentInterval>) =
            if self.alignments.len() > 10 {
                self.alignments.iter().cloned().partition(|ai| {
                    (!is_canonical(ai) && ai.aln_score > first_max)
                        || ai.map_qual > ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST
                })
            } else {
                (self.alignments.clone(), Vec::new())
            };
        let max_canonical = max_canonical_score(&good);
        let non_canonical = better_non_canonical_mapping(canonical, &good, max_canonical);
        if let Some(mapping) = &non_canonical {
            if let Some(position) = good.iter().position(|ai| ai == mapping) {
                good.remove(position);
            }
        }
        // `Sets.powerSet(new HashSet<>(goodMappings))`: equal alignments are one element.
        let mut unique: Vec<AlignmentInterval> = Vec::new();
        for alignment in &good {
            if !unique.contains(alignment) {
                unique.push(alignment.clone());
            }
        }
        let mut configurations = Vec::new();
        for mask in 0u64..(1u64 << unique.len()) {
            let mut configuration: Vec<AlignmentInterval> = unique
                .iter()
                .enumerate()
                .filter(|(bit, _)| mask & (1 << bit) != 0)
                .map(|(_, ai)| ai.clone())
                .collect();
            configuration.sort_by(alignment_order);
            configurations.push(configuration);
        }
        let scores: Vec<f64> = configurations
            .iter()
            .map(|configuration| score_of_configuration(configuration, canonical, max_canonical))
            .collect();
        let max_score = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if max_score == f64::NEG_INFINITY {
            return Err(SvError::new(
                GATK_EXCEPTION,
                format!(
                    "Cannot find best-scoring configuration on alignments of contig: {}",
                    self.name
                ),
            ));
        }
        let mut picked = Vec::new();
        for (configuration, score) in configurations.into_iter().zip(scores) {
            let tolerance = ulp(score).max(tolerance);
            if score >= max_score || max_score - score <= tolerance {
                let mut rest = good.clone();
                remove_all(&mut rest, &configuration);
                rest.extend(bad.iter().cloned());
                picked.push(GoodAndBadMappings {
                    good: configuration,
                    bad: rest,
                    non_canonical: non_canonical.clone(),
                });
            }
        }
        Ok(picked)
    }
}

/// `AlignedContig.getBetterNonCanonicalMapping`: a single non-canonical alignment that explains
/// the contig at least as well as the canonical ones, taken out of consideration.
fn better_non_canonical_mapping(
    canonical: &[String],
    good: &[AlignmentInterval],
    max_canonical_score: i32,
) -> Option<AlignmentInterval> {
    let (canonical_mappings, non_canonical): (Vec<AlignmentInterval>, Vec<AlignmentInterval>) =
        good.iter()
            .cloned()
            .partition(|ai| canonical.contains(&ai.reference_span.contig));
    if canonical_mappings.is_empty() {
        return None;
    }
    if non_canonical.len() == 1
        && (canonical_mappings.len() > 1
            || canonical_mappings[0]
                .contains_gap_of_equal_or_larger_size(STRUCTURAL_VARIANT_SIZE_LOWER_BOUND))
    {
        let canonical_score =
            score_of_configuration(&canonical_mappings, canonical, max_canonical_score);
        let non_canonical_score =
            score_of_configuration(&non_canonical, canonical, max_canonical_score);
        if canonical_score > non_canonical_score {
            None
        } else {
            Some(non_canonical[0].clone())
        }
    } else {
        None
    }
}

/// `computeScoreOfConfiguration`: the mapping-quality-weighted read coverage, less the bases two
/// alignments both claim.
fn score_of_configuration(
    configuration: &[AlignmentInterval],
    canonical: &[String],
    max_canonical_score: i32,
) -> f64 {
    let mut explained = 0.0f64;
    for alignment in configuration {
        let length = alignment.size_on_read();
        let weight = if canonical.contains(&alignment.reference_span.contig) {
            alignment.map_qual as f64 / COVERAGE_MQ_NORMALIZATION_CONST
        } else {
            (alignment.map_qual as f64 / COVERAGE_MQ_NORMALIZATION_CONST).max(
                if alignment.aln_score > max_canonical_score {
                    1.0
                } else {
                    0.0
                },
            )
        };
        explained += weight * length as f64;
    }
    let mut redundancy = 0;
    for i in 0..configuration.len().saturating_sub(1) {
        for j in i + 1..configuration.len() {
            redundancy += overlap_on_contig(&configuration[i], &configuration[j]);
        }
    }
    explained - redundancy as f64
}

/// `splitGapsAndDropAlignmentContainedByOtherOnRead`, the variant the tool uses.
fn split_gaps(configuration: GoodAndBadMappings) -> Result<GoodAndBadMappings, SvError> {
    let mut gap_split = Vec::new();
    for alignment in &configuration.good {
        if alignment.contains_gap_of_equal_or_larger_size(STRUCTURAL_VARIANT_SIZE_LOWER_BOUND) {
            gap_split.extend(split_gapped_alignment(
                alignment,
                STRUCTURAL_VARIANT_SIZE_LOWER_BOUND,
                alignment.cigar.count_unclipped_read_bases(),
            )?);
        } else {
            gap_split.push(alignment.clone());
        }
    }
    let mut bad = configuration.bad.clone();
    for i in 0..gap_split.len() {
        for j in i + 1..gap_split.len() {
            if gap_split[i].contains_on_read(&gap_split[j]) {
                bad.push(gap_split[j].clone());
            } else if gap_split[j].contains_on_read(&gap_split[i]) {
                bad.push(gap_split[i].clone());
            }
        }
    }
    remove_all(&mut gap_split, &bad);
    gap_split.sort_by(alignment_order);
    Ok(GoodAndBadMappings {
        good: gap_split,
        bad,
        non_canonical: configuration.non_canonical,
    })
}

/// `removeNonUniqueMappings` at MQ 20 and ten unique read bases, which only a configuration of
/// three or more alignments reaches.
fn remove_non_unique_mappings(mappings: GoodAndBadMappings) -> GoodAndBadMappings {
    if mappings.good.len() <= 2 {
        return mappings;
    }
    let mut selected = Vec::new();
    let mut low = mappings.bad.clone();
    for alignment in &mappings.good {
        if alignment.map_qual >= ALIGNMENT_MQ_THRESHOLD {
            selected.push(alignment.clone());
        } else {
            low.push(alignment.clone());
        }
    }
    // `removeDueToShortReadSpan`. The overlaps are held in a map keyed by the alignment, so two
    // equal alignments share the entry the later of them wrote.
    let overlaps = max_overlap_pairs(&selected);
    let lookup = |alignment: &AlignmentInterval| -> (i32, i32) {
        overlaps
            .iter()
            .rev()
            .find(|(key, _)| key == alignment)
            .map(|(_, value)| *value)
            .unwrap_or((-1, -1))
    };
    let mut kept = Vec::new();
    for alignment in selected {
        let (front, rear) = lookup(&alignment);
        let unique =
            alignment.end_in_contig - alignment.start_in_contig + 1 - front.max(0) - rear.max(0);
        if unique < ALIGNMENT_LOW_READ_UNIQUENESS_THRESHOLD {
            low.push(alignment);
        } else {
            kept.push(alignment);
        }
    }
    GoodAndBadMappings {
        good: kept,
        bad: low,
        non_canonical: mappings.non_canonical,
    }
}

/// `getMaxOverlapPairs`: for each alignment, the largest overlap with one before and one after.
fn max_overlap_pairs(configuration: &[AlignmentInterval]) -> Vec<(AlignmentInterval, (i32, i32))> {
    // (maxFront, maxRear), each an (index, bases) pair, as `TempMaxOverlapInfo` holds them.
    let mut info = vec![((-1, -1), (-1, -1)); configuration.len()];
    for i in 0..configuration.len().saturating_sub(1) {
        let mut max_rear_bases = -1;
        let mut max_rear_index: i32 = -1;
        for j in i + 1..configuration.len() {
            let overlap = overlap_on_contig(&configuration[i], &configuration[j]);
            if overlap > max_rear_bases {
                max_rear_bases = overlap;
                max_rear_index = j as i32;
            } else {
                break;
            }
        }
        if max_rear_bases > 0 {
            info[i] = (info[i].0, (max_rear_index, max_rear_bases));
            let rear = max_rear_index as usize;
            let old = info[rear];
            if old.0 .1 < max_rear_bases {
                info[rear] = ((i as i32, max_rear_bases), old.1);
            }
        }
    }
    configuration
        .iter()
        .zip(info)
        .map(|(alignment, (front, rear))| (alignment.clone(), (front.1, rear.1)))
        .collect()
}

/// `alignmentShouldNotBeStitchedTogether`.
fn not_stitchable(alignments: &[AlignmentInterval]) -> Result<bool, SvError> {
    Ok(alignments.len() != 2
        || !simple_chimera_with_stitchable_alignments(&alignments[0], &alignments[1])?)
}

/// `AlignedContig.simpleChimeraWithStichableAlignments`: two pieces that are really one
/// alignment, as in the signature of GATK ticket 4951.
pub fn simple_chimera_with_stitchable_alignments(
    one: &AlignmentInterval,
    two: &AlignmentInterval,
) -> Result<bool, SvError> {
    if one.start_in_contig > two.start_in_contig {
        return Err(SvError::new(
            ILLEGAL_ARGUMENT,
            format!(
                "Assumption that input intervals are sorted by their starts on read is violated.\tFirst: {}\tSecond: {}",
                one.to_packed_string(),
                two.to_packed_string()
            ),
        ));
    }
    if one.reference_span.contig != two.reference_span.contig
        || one.forward_strand != two.forward_strand
        || one.contains_on_read(two)
        || two.contains_on_read(one)
        || one.contains_on_ref(two)
        || two.contains_on_ref(one)
    {
        return Ok(false);
    }
    let ref_order_swap =
        one.forward_strand != (one.reference_span.start < two.reference_span.start);
    if ref_order_swap {
        return Ok(false);
    }
    let on_contig = overlap_on_contig(one, two);
    let on_ref = overlap_on_ref_span(one, two);
    if on_contig == 0 && on_ref == 0 {
        Ok(two.reference_span.start - one.reference_span.end == 1
            && two.start_in_contig - one.end_in_contig == 1)
    } else {
        Ok(on_contig == on_ref)
    }
}

/// `AssemblyContigWithFineTunedAlignments.AlignmentSignatureBasicType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentSignature {
    Normal,
    Unknown,
    SimpleChimera,
    Complex,
}

/// `AssemblyContigWithFineTunedAlignments`: a contig with the alignments it was left with, the
/// packed strings of those it lost, and whether another configuration explained it as well.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FineTunedContig {
    pub contig: AlignedContig,
    pub insertion_mappings: Vec<String>,
    pub ambiguous: bool,
    pub non_canonical_sa_tag: String,
}

impl FineTunedContig {
    /// `getAlignmentSignatureBasicType`.
    pub fn signature(&self) -> AlignmentSignature {
        let alignments = &self.contig.alignments;
        if self.ambiguous || self.has_incomplete_picture() {
            AlignmentSignature::Unknown
        } else if alignments.len() < 2 {
            if alignments.len() == 1 && self.insertion_mappings.is_empty() {
                AlignmentSignature::Normal
            } else {
                AlignmentSignature::Unknown
            }
        } else if alignments.len() == 2 {
            AlignmentSignature::SimpleChimera
        } else {
            AlignmentSignature::Complex
        }
    }

    /// `hasIncompletePicture`.
    fn has_incomplete_picture(&self) -> bool {
        let alignments = &self.contig.alignments;
        if alignments.len() <= 1 {
            false
        } else if alignments.len() == 2 {
            has_incomplete_picture_from_two_alignments(&alignments[0], &alignments[1])
        } else {
            self.has_incomplete_picture_from_multiple_alignments()
        }
    }

    fn has_incomplete_picture_from_multiple_alignments(&self) -> bool {
        let alignments = &self.contig.alignments;
        let head = &alignments[0];
        let tail = &alignments[alignments.len() - 1];
        if head.reference_span.contig != tail.reference_span.contig
            || head.forward_strand != tail.forward_strand
        {
            return true;
        }
        let (span_head, span_tail) = (&head.reference_span, &tail.reference_span);
        if span_head.contains(span_tail) || span_tail.contains(span_head) {
            return true;
        }
        let valid = Interval::of(
            &span_head.contig,
            span_head.start.min(span_tail.start),
            span_head.end.max(span_tail.end),
        );
        let not_complete_dup_region = alignments[1..alignments.len() - 1].iter().any(|middle| {
            middle.reference_span.overlaps(&valid) && !valid.contains(&middle.reference_span)
        });
        if not_complete_dup_region {
            return true;
        }
        if head.forward_strand {
            span_head.start >= span_tail.start
        } else {
            span_head.end <= span_tail.end
        }
    }
}

/// `AssemblyContigWithFineTunedAlignments.hasIncompletePictureFromTwoAlignments`.
pub fn has_incomplete_picture_from_two_alignments(
    head: &AlignmentInterval,
    tail: &AlignmentInterval,
) -> bool {
    let (one, two) = (&head.reference_span, &tail.reference_span);
    if one.contig != two.contig {
        return false;
    }
    if one.contains(two) || two.contains(one) {
        return true;
    }
    if head.forward_strand != tail.forward_strand {
        one.overlaps(two)
    } else if head.forward_strand {
        one.start > two.start && one.start <= two.end
    } else {
        two.start > one.start && two.start <= one.end
    }
}
