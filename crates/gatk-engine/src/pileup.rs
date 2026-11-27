//! Ported from `org.broadinstitute.hellbender.utils.pileup.PileupElement` (GATK 4.7.0.0;
//! `PileupElement.java` is unchanged since 4.6.2.0).
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
//!
//! # Reading this file without reading Rust
//!
//! - `PileupElement<'a>` holds `&'a BamRecord`, a **borrowed** read: the element points at a read
//!   someone else owns, and the `'a` is the compiler's proof that the element cannot outlive it.
//!   Java's `final GATKRead read` field, with that guarantee added.
//! - `Option<T>` is the reference's nullable return: `Some(x)` is a value, `None` is `null`. A
//!   `?` after an `Option` returns `None` from the whole function at once, the way a null check
//!   followed by `return null` would.
//! - `x as usize` and `x as i32` are Java's casts. Every `as usize` below indexes an array with an
//!   `int` offset that the reference also uses unchecked.

use crate::read_utils;
use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;

/// `BaseUtils.Base.D.base`, the byte a deletion reports.
pub const DELETION_BASE: u8 = b'D';
/// `PileupElement.DELETION_QUAL`.
pub const DELETION_QUAL: u8 = 16;

/// `PileupElement.ON_GENOME_OPERATORS`, the `EnumSet` `M`, `=`, `X`, `D`.
///
/// `N` is not in it, so a skipped region is crossed by the nearest-element walks like an insertion
/// or a clip. `matches!` is a membership test against the listed variants.
fn is_on_genome(op: Op) -> bool {
    matches!(op, Op::M | Op::Eq | Op::X | Op::D)
}

/// The direction `getBetween` and the nearest-element walks run in.
///
/// The reference's `Direction` carries its increment (`PREV(-1), NEXT(1)`); here each walk maps
/// the variant to the increment with a `match` instead.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    Prev,
    Next,
}

/// `PileupElement`: a read, an offset into its bases, and where in its cigar that offset sits.
///
/// The fields are `pub` where the reference has getters (`getRead`, `getOffset`,
/// `getCurrentCigarElement`, `getCurrentCigarOffset`, `getOffsetInCurrentCigar`).
#[derive(Debug, Clone)]
pub struct PileupElement<'a> {
    pub read: &'a BamRecord,
    /// The offset into the read's bases. For a deletion this is the offset of the last `M`/`=`/`X`
    /// position, which is why the base and quality accessors have to special-case it.
    pub offset: i32,
    /// The cigar element this position is in.
    pub current_cigar_element: CigarElement,
    /// The index of `current_cigar_element` in the read's cigar; the walks below start from it.
    pub current_cigar_offset: i32,
    /// How far into `current_cigar_element` this position is, zero-based.
    pub offset_in_current_cigar: i32,
}

impl<'a> PileupElement<'a> {
    /// The five-argument `PileupElement` constructor.
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
    ///
    /// The reference validates `!(isLeftEdge() || isRightEdge())` and throws with the message
    /// "Cannot make a pileup element from an edge alignment state". This answers `None` instead,
    /// and the caller decides whether that is an error.
    pub fn from_state(
        read: &'a BamRecord,
        machine: &crate::alignment_state::AlignmentStateMachine,
    ) -> Option<Self> {
        if machine.is_left_edge() || machine.is_right_edge() {
            return None;
        }
        // The `?` after `current_cigar_element()` passes a missing element on as `None` too,
        // where the reference would construct the element with a null cigar element.
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
        // `step_forward_on_genome` answers `Result<Option<_>, _>`: `.ok()?` turns a malformed
        // read into `None` too, and `.is_some()` is the reference's `!= null`, false once the
        // machine steps off the right edge.
        while machine.step_forward_on_genome().ok()?.is_some() {
            if machine.read_offset() == offset {
                return PileupElement::from_state(read, &machine);
            }
        }
        None
    }

    /// `isDeletion`: the current element is `D`.
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

    /// `getQual`: [`DELETION_QUAL`] for a deletion, the base quality at the offset otherwise.
    pub fn qual(&self) -> u8 {
        if self.is_deletion() {
            DELETION_QUAL
        } else {
            self.read.base_qualities[self.offset as usize]
        }
    }

    /// `getBaseInsertionQual`: the `BI` quality at the offset, or the flat Q45 default when the
    /// read has no `BI` tag (`ReadUtils.getBaseInsertionQualities`).
    pub fn base_insertion_qual(&self) -> u8 {
        if self.is_deletion() {
            DELETION_QUAL
        } else {
            read_utils::base_insertion_qualities(self.read)[self.offset as usize]
        }
    }

    /// `getBaseDeletionQual`: as [`Self::base_insertion_qual`], for the `BD` tag.
    pub fn base_deletion_qual(&self) -> u8 {
        if self.is_deletion() {
            DELETION_QUAL
        } else {
            read_utils::base_deletion_qualities(self.read)[self.offset as usize]
        }
    }

    /// `getMappingQual`.
    pub fn mapping_qual(&self) -> u8 {
        self.read.mapping_quality
    }

    /// `atStartOfCurrentCigar`.
    pub fn at_start_of_current_cigar(&self) -> bool {
        self.offset_in_current_cigar == 0
    }

    /// `atEndOfCurrentCigar`: the last position of the current element.
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
        // The reference's `for` with `i += increment`, written as a `while` because the step is
        // signed and Rust's `for` ranges only count up.
        while index >= 0 && index < count {
            let element = self.read.cigar.elements[index as usize];
            if is_on_genome(element.op) {
                return Some(element);
            }
            index += increment;
        }
        None
    }

    /// `getPreviousOnGenomeCigarElement`.
    pub fn previous_on_genome_cigar_element(&self) -> Option<CigarElement> {
        self.nearest_on_genome(Direction::Prev)
    }

    /// `getNextOnGenomeCigarElement`.
    pub fn next_on_genome_cigar_element(&self) -> Option<CigarElement> {
        self.nearest_on_genome(Direction::Next)
    }

    /// `getBetween`: the run of non-on-genome elements between this position and its neighbour,
    /// which stops at the first on-genome operator rather than at the neighbour itself.
    ///
    /// Walking backwards collects the run in reverse, so it is reversed at the end, as upstream
    /// (`Collections.reverse(elements)`), to come back in cigar order.
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

    /// `getBetweenNextPosition`: the mirror of [`Self::between_prev_position`], non-empty only at
    /// the *last* position of the element.
    pub fn between_next_position(&self) -> Vec<CigarElement> {
        if self.at_end_of_current_cigar() {
            self.between(Direction::Next)
        } else {
            Vec::new()
        }
    }

    /// `isImmediatelyAfter(op)`: first position of the element, and the adjacent element (not the
    /// nearest on-genome one) is `op`.
    fn is_immediately_after(&self, op: Op) -> bool {
        self.at_start_of_current_cigar() && self.adjacent_operator(Direction::Prev) == Some(op)
    }

    /// `isImmediatelyBefore(op)`: last position of the element, and the adjacent element is `op`.
    fn is_immediately_before(&self, op: Op) -> bool {
        self.at_end_of_current_cigar() && self.adjacent_operator(Direction::Next) == Some(op)
    }

    /// `isAfterInsertion`.
    pub fn is_after_insertion(&self) -> bool {
        self.is_immediately_after(Op::I)
    }

    /// `isBeforeInsertion`. Adjacent element only: see the module note on `3M2S3I3M`.
    pub fn is_before_insertion(&self) -> bool {
        self.is_immediately_before(Op::I)
    }

    /// `isAfterSoftClip`.
    pub fn is_after_soft_clip(&self) -> bool {
        self.is_immediately_after(Op::S)
    }

    /// `isBeforeSoftClip`.
    pub fn is_before_soft_clip(&self) -> bool {
        self.is_immediately_before(Op::S)
    }

    /// `isNextToSoftClip`: either side.
    pub fn is_next_to_soft_clip(&self) -> bool {
        self.is_after_soft_clip() || self.is_before_soft_clip()
    }

    /// `isBeforeDeletionStart`: this is the last base of an element and the next *on-genome*
    /// element is a deletion. The on-genome walk is what makes this see a deletion through an
    /// intervening insertion, unlike `isBeforeInsertion`.
    ///
    /// `.map(|e| e.op) == Some(Op::D)` is the reference's `hasOperator`: a missing element maps to
    /// `None`, which is not equal to `Some(Op::D)`, so no separate null check is needed.
    pub fn is_before_deletion_start(&self) -> bool {
        !self.is_deletion()
            && self.at_end_of_current_cigar()
            && self.next_on_genome_cigar_element().map(|e| e.op) == Some(Op::D)
    }

    /// `isAfterDeletionEnd`: the mirror of [`Self::is_before_deletion_start`].
    pub fn is_after_deletion_end(&self) -> bool {
        !self.is_deletion()
            && self.at_start_of_current_cigar()
            && self.previous_on_genome_cigar_element().map(|e| e.op) == Some(Op::D)
    }

    /// `getNextIndelCigarElement`: the deletion or insertion that follows this position, if any.
    ///
    /// The reference checks each branch's answer with `Utils.validate` (the next on-genome element
    /// is a `D`; the first element between is an `I`). Both hold by construction of the predicate
    /// that chose the branch, so the port does not repeat them.
    fn next_indel_cigar_element(&self) -> Option<CigarElement> {
        if self.is_before_deletion_start() {
            self.next_on_genome_cigar_element()
        } else if self.is_before_insertion() {
            // `.first().copied()` is `get(0)`, answering `None` on an empty list rather than
            // throwing; `is_before_insertion` guarantees the list is not empty here.
            self.between_next_position().first().copied()
        } else {
            None
        }
    }

    /// `getLengthOfImmediatelyFollowingIndel`: that element's length, or 0 when there is none.
    pub fn length_of_immediately_following_indel(&self) -> u32 {
        self.next_indel_cigar_element()
            .map(|e| e.length)
            .unwrap_or(0)
    }

    /// `getBasesOfImmediatelyFollowingInsertion`: the inserted bases, read straight out of the
    /// read starting one past this offset.
    ///
    /// The reference returns a `String`; this returns the bytes. A following deletion answers
    /// `None`, as the reference's `null` does.
    pub fn bases_of_immediately_following_insertion(&self) -> Option<Vec<u8>> {
        let element = self.next_indel_cigar_element()?;
        if element.op != Op::I {
            return None;
        }
        let from = (self.offset + 1) as usize;
        let to = from + element.length as usize;
        // `.get(from..to)` is the slice `read_bases[from, to)`, or `None` if it runs past the
        // end, which a cigar whose `I` bases are all in the read cannot reach.
        self.read.read_bases.get(from..to).map(|s| s.to_vec())
    }

    /// `isUsableBaseForAnnotation`. `MAPPING_QUALITY_UNAVAILABLE` is 255 and
    /// `MIN_USABLE_Q_SCORE` is 6.
    ///
    /// ```java
    /// return !( p.isDeletion() ||
    ///         p.getMappingQual() == 0 ||
    ///         p.getMappingQual() == QualityUtils.MAPPING_QUALITY_UNAVAILABLE ||
    ///         ((int) p.getQual()) < QualityUtils.MIN_USABLE_Q_SCORE);
    /// ```
    ///
    /// The reference is a static method taking the element; this is a method on it.
    pub fn is_usable_base_for_annotation(&self) -> bool {
        !(self.is_deletion()
            || self.mapping_qual() == 0
            || self.mapping_qual() == 255
            || (self.qual() as i32) < 6)
    }
}
