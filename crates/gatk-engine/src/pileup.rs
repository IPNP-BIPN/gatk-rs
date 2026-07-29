//! Ported from `org.broadinstitute.hellbender.utils.pileup.PileupElement` (GATK 4.6.2.0).
//!
//! One read's contribution to one locus. [`crate::alignment_state::AlignmentStateMachine`] decides
//! *where* a read stops; this decides what the caller sees there, and it is where the indel-aware
//! questions live: is this base immediately before an insertion, immediately after a deletion,
//! next to a soft clip.
//!
//! What makes those questions non-obvious is that the state machine never stops on an insertion or
//! a soft clip: it consumes them whole. So "before an insertion" is not a position the caller can
//! reach, it is a property of the *last* base of the preceding element, and every one of these
//! predicates is a walk through the cigar from the current element outwards.
//!
//! Three details are decisions rather than plumbing:
//!
//!  * **the two families of navigation differ.** `getAdjacentOperator` looks at exactly the next
//!    cigar element, whatever it is, while `getNearestOnGenomeCigarElement` skips everything that
//!    is not `M`, `=`, `X` or `D`. So `isBeforeInsertion` on `3M2S3I3M` is *false*, because the
//!    adjacent element is the soft clip and not the insertion, while `isBeforeDeletionStart`
//!    would still see a deletion through the clip;
//!  * **a deletion has no base and no quality of its own.** `getBase` answers `D` and every
//!    quality answers 16, whatever the read holds at that offset;
//!  * **the indel qualities are Q45 by default.** A read without `BI`/`BD` tags gets a flat
//!    array, so a pileup can report an insertion quality for a read that never had one.

use crate::read_utils;
use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;

/// `BaseUtils.Base.D.base`, the byte a deletion reports.
pub const DELETION_BASE: u8 = b'D';
/// `PileupElement.DELETION_QUAL`.
pub const DELETION_QUAL: u8 = 16;

/// `PileupElement.ON_GENOME_OPERATORS`.
fn is_on_genome(op: Op) -> bool {
    matches!(op, Op::M | Op::Eq | Op::X | Op::D)
}

/// The direction `getBetween` and the nearest-element walks run in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    Prev,
    Next,
}

/// `PileupElement`: a read, an offset into its bases, and where in its cigar that offset sits.
#[derive(Debug, Clone)]
pub struct PileupElement<'a> {
    pub read: &'a BamRecord,
    /// The offset into the read's bases. For a deletion this is the offset of the last `M`/`=`/`X`
    /// position, which is why the base and quality accessors have to special-case it.
    pub offset: i32,
    pub current_cigar_element: CigarElement,
    pub current_cigar_offset: i32,
    pub offset_in_current_cigar: i32,
}

impl<'a> PileupElement<'a> {
    pub fn new(
        read: &'a BamRecord,
        offset: i32,
        current_cigar_element: CigarElement,
        current_cigar_offset: i32,
        offset_in_current_cigar: i32,
    ) -> Self {
        // No bounds checking, as upstream: the reference documents dropping it because this class
        // is a HaplotypeCaller hotspot, and a port that added it would answer where the reference
        // panics.
        PileupElement {
            read,
            offset,
            current_cigar_element,
            current_cigar_offset,
            offset_in_current_cigar,
        }
    }

    /// `AlignmentStateMachine.makePileupElement`, which refuses both edges.
    pub fn from_state(
        read: &'a BamRecord,
        machine: &crate::alignment_state::AlignmentStateMachine,
    ) -> Option<Self> {
        if machine.is_left_edge() || machine.is_right_edge() {
            return None;
        }
        Some(PileupElement::new(
            read,
            machine.read_offset(),
            machine.current_cigar_element()?,
            machine.current_cigar_element_offset(),
            machine.offset_into_current_cigar_element(),
        ))
    }

    /// `createPileupForReadAndOffset`: run the machine until it reaches this read offset.
    ///
    /// Returns `None` where the reference throws `IllegalStateException`, which happens for an
    /// offset the alignment never visits: inside a soft clip, or inside an insertion.
    pub fn for_read_and_offset(read: &'a BamRecord, offset: i32) -> Option<Self> {
        let mut machine = crate::alignment_state::AlignmentStateMachine::new(read);
        while machine.step_forward_on_genome().ok()?.is_some() {
            if machine.read_offset() == offset {
                return PileupElement::from_state(read, &machine);
            }
        }
        None
    }

    pub fn is_deletion(&self) -> bool {
        self.current_cigar_element.op == Op::D
    }

    /// `getBase`: `D` for a deletion, whatever the read holds at the offset otherwise.
    pub fn base(&self) -> u8 {
        if self.is_deletion() {
            DELETION_BASE
        } else {
            self.read.read_bases[self.offset as usize]
        }
    }

    pub fn qual(&self) -> u8 {
        if self.is_deletion() {
            DELETION_QUAL
        } else {
            self.read.base_qualities[self.offset as usize]
        }
    }

    pub fn base_insertion_qual(&self) -> u8 {
        if self.is_deletion() {
            DELETION_QUAL
        } else {
            read_utils::base_insertion_qualities(self.read)[self.offset as usize]
        }
    }

    pub fn base_deletion_qual(&self) -> u8 {
        if self.is_deletion() {
            DELETION_QUAL
        } else {
            read_utils::base_deletion_qualities(self.read)[self.offset as usize]
        }
    }

    pub fn mapping_qual(&self) -> u8 {
        self.read.mapping_quality
    }

    pub fn at_start_of_current_cigar(&self) -> bool {
        self.offset_in_current_cigar == 0
    }

    pub fn at_end_of_current_cigar(&self) -> bool {
        self.offset_in_current_cigar == self.current_cigar_element.length as i32 - 1
    }

    /// `getAdjacentOperator`: exactly the neighbouring cigar element, on-genome or not.
    fn adjacent_operator(&self, direction: Direction) -> Option<Op> {
        let increment: i32 = match direction {
            Direction::Prev => -1,
            Direction::Next => 1,
        };
        let index = self.current_cigar_offset + increment;
        if index < 0 || index >= self.read.cigar.elements.len() as i32 {
            return None;
        }
        Some(self.read.cigar.elements[index as usize].op)
    }

    /// `getNearestOnGenomeCigarElement`: the first `M`, `=`, `X` or `D` in that direction, skipping
    /// everything else.
    fn nearest_on_genome(&self, direction: Direction) -> Option<CigarElement> {
        let increment: i32 = match direction {
            Direction::Prev => -1,
            Direction::Next => 1,
        };
        let count = self.read.cigar.elements.len() as i32;
        let mut index = self.current_cigar_offset + increment;
        while index >= 0 && index < count {
            let element = self.read.cigar.elements[index as usize];
            if is_on_genome(element.op) {
                return Some(element);
            }
            index += increment;
        }
        None
    }

    pub fn previous_on_genome_cigar_element(&self) -> Option<CigarElement> {
        self.nearest_on_genome(Direction::Prev)
    }

    pub fn next_on_genome_cigar_element(&self) -> Option<CigarElement> {
        self.nearest_on_genome(Direction::Next)
    }

    /// `getBetween`: the run of non-on-genome elements between this position and its neighbour,
    /// which stops at the first on-genome operator rather than at the neighbour itself.
    fn between(&self, direction: Direction) -> Vec<CigarElement> {
        let increment: i32 = match direction {
            Direction::Prev => -1,
            Direction::Next => 1,
        };
        let count = self.read.cigar.elements.len() as i32;
        let mut elements = Vec::new();
        let mut index = self.current_cigar_offset + increment;
        while index >= 0 && index < count {
            let element = self.read.cigar.elements[index as usize];
            if is_on_genome(element.op) {
                break;
            }
            elements.push(element);
            index += increment;
        }
        if increment < 0 {
            elements.reverse();
        }
        elements
    }

    /// `getBetweenPrevPosition`: empty unless this is the *first* position of its cigar element,
    /// because anywhere else there is nothing between here and the previous genomic position.
    pub fn between_prev_position(&self) -> Vec<CigarElement> {
        if self.at_start_of_current_cigar() {
            self.between(Direction::Prev)
        } else {
            Vec::new()
        }
    }

    pub fn between_next_position(&self) -> Vec<CigarElement> {
        if self.at_end_of_current_cigar() {
            self.between(Direction::Next)
        } else {
            Vec::new()
        }
    }

    fn is_immediately_after(&self, op: Op) -> bool {
        self.at_start_of_current_cigar() && self.adjacent_operator(Direction::Prev) == Some(op)
    }

    fn is_immediately_before(&self, op: Op) -> bool {
        self.at_end_of_current_cigar() && self.adjacent_operator(Direction::Next) == Some(op)
    }

    pub fn is_after_insertion(&self) -> bool {
        self.is_immediately_after(Op::I)
    }

    pub fn is_before_insertion(&self) -> bool {
        self.is_immediately_before(Op::I)
    }

    pub fn is_after_soft_clip(&self) -> bool {
        self.is_immediately_after(Op::S)
    }

    pub fn is_before_soft_clip(&self) -> bool {
        self.is_immediately_before(Op::S)
    }

    pub fn is_next_to_soft_clip(&self) -> bool {
        self.is_after_soft_clip() || self.is_before_soft_clip()
    }

    /// `isBeforeDeletionStart`: this is the last base of an element and the next *on-genome*
    /// element is a deletion. The on-genome walk is what makes this see a deletion through an
    /// intervening insertion, unlike `isBeforeInsertion`.
    pub fn is_before_deletion_start(&self) -> bool {
        !self.is_deletion()
            && self.at_end_of_current_cigar()
            && self.next_on_genome_cigar_element().map(|e| e.op) == Some(Op::D)
    }

    pub fn is_after_deletion_end(&self) -> bool {
        !self.is_deletion()
            && self.at_start_of_current_cigar()
            && self.previous_on_genome_cigar_element().map(|e| e.op) == Some(Op::D)
    }

    /// `getNextIndelCigarElement`: the deletion or insertion that follows this position, if any.
    fn next_indel_cigar_element(&self) -> Option<CigarElement> {
        if self.is_before_deletion_start() {
            self.next_on_genome_cigar_element()
        } else if self.is_before_insertion() {
            self.between_next_position().first().copied()
        } else {
            None
        }
    }

    pub fn length_of_immediately_following_indel(&self) -> u32 {
        self.next_indel_cigar_element()
            .map(|e| e.length)
            .unwrap_or(0)
    }

    /// `getBasesOfImmediatelyFollowingInsertion`: the inserted bases, read straight out of the
    /// read starting one past this offset.
    pub fn bases_of_immediately_following_insertion(&self) -> Option<Vec<u8>> {
        let element = self.next_indel_cigar_element()?;
        if element.op != Op::I {
            return None;
        }
        let from = (self.offset + 1) as usize;
        let to = from + element.length as usize;
        self.read.read_bases.get(from..to).map(|s| s.to_vec())
    }

    /// `isUsableBaseForAnnotation`. `MAPPING_QUALITY_UNAVAILABLE` is 255 and
    /// `MIN_USABLE_Q_SCORE` is 6.
    pub fn is_usable_base_for_annotation(&self) -> bool {
        !(self.is_deletion()
            || self.mapping_qual() == 0
            || self.mapping_qual() == 255
            || (self.qual() as i32) < 6)
    }
}
