//! Ported from `org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine`
//! (GATK 4.7.0.0; `AlignmentStateMachine.java` is unchanged since 4.6.2.0).
//!
//! One read, walked one *reference* base at a time. This is the bottom of every locus-based tool:
//! `LocusIteratorByState` runs one of these per read and merges their positions into a pileup, so
//! every depth, every allele count and every annotation that reads a pileup rests on the offsets
//! this produces.
//!
//! Four behaviours here are decisions rather than arithmetic:
//!
//!  * **a deletion is returned once per reference base it spans**, so a `10D` yields ten stops
//!    with the read offset frozen. `N` behaves identically and is a different logical claim;
//!  * **`I`, `S`, `H` and `P` are consumed whole inside one step.** The machine never stops on
//!    them, so a caller sees the base *after* an insertion and has to look backwards to find it,
//!    which is exactly what `PileupElement` does;
//!  * **stepping off the right edge advances the genome offset one past the read.** The
//!    reference's comment justifies it for a read that ends on an indel, so the final position
//!    models the next reference base after the indel, but the `genomeOffset++` under that comment
//!    (`AlignmentStateMachine.java:253`) is unconditional: a read ending in `M` gets it too;
//!  * **a cigar that starts or ends with a deletion is a malformed read**, and throws rather than
//!    being tolerated, even though the SAM specification permits it. Zero-length elements, by
//!    contrast, are skipped silently.
//!
//! # Reading this file without reading Rust
//!
//! - `AlignmentStateMachine<'a>` holds `&'a [CigarElement]`, a borrowed read-only view of the
//!   read's cigar. The `'a` says the machine cannot outlive the read it walks; Java's machine
//!   keeps the `GATKRead` itself alive instead.
//! - `Option<CigarElement>` is the reference's nullable `currentElement`: `None` is `null`, and
//!   `.map(|e| e.op)` is `currentElement == null ? null : currentElement.getOperator()`
//!   (`AlignmentStateMachine.java:189`).
//! - `Result<Option<Op>, MalformedRead>` is "an operator, or `null`, or a thrown
//!   `MalformedRead`". `Ok(None)` is the reference's `return null`, `Err(..)` its `throw`.
//! - `x as i32` and `x as usize` convert between the unsigned element length or index and Java's
//!   `int`. They never fail and would wrap on overflow; a cigar length past `i32::MAX` is not a
//!   read either side can hold.

use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;

/// `UserException.MalformedRead`, raised by the two cigars this machine refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MalformedRead {
    /// A deletion reached before any read base was consumed. The test is the reference's
    /// `if (readOffset < 0)` (`AlignmentStateMachine.java:273`), not "the first on-genome element":
    /// `2D5M` and `5H2D5M` refuse, while `5S2D5M` does not, because the soft clip has already
    /// moved the read offset to 4.
    StartsWithDeletion,
    /// The cigar's last element is a deletion, found when the walk runs out of elements.
    EndsWithDeletion,
}

/// `AlignmentStateMachine`: a read, and where along the genome we are in it.
pub struct AlignmentStateMachine<'a> {
    elements: &'a [CigarElement],
    /// `read.getLength()`, the number of read bases, which the right edge is tested against.
    read_length: i32,
    alignment_start: i32,
    /// Index of the current element in the cigar: -1 on the left edge, the element count once off
    /// the right edge.
    current_cigar_element_offset: i32,
    /// Index of the last read base processed; -1 before the first.
    read_offset: i32,
    /// Reference bases advanced from the alignment start; -1 on the left edge.
    genome_offset: i32,
    /// `None` on both edges, which is why the operator accessor is nullable upstream.
    current_element: Option<CigarElement>,
    offset_into_current_cigar_element: i32,
}

impl<'a> AlignmentStateMachine<'a> {
    /// `initializeAsLeftEdge`: one base before the alignment, so the first step lands on it.
    pub fn new(read: &'a BamRecord) -> Self {
        AlignmentStateMachine {
            elements: &read.cigar.elements,
            read_length: read.read_bases.len() as i32,
            alignment_start: read.alignment_start,
            current_cigar_element_offset: -1,
            read_offset: -1,
            genome_offset: -1,
            current_element: None,
            offset_into_current_cigar_element: -1,
        }
    }

    /// `isLeftEdge`: before the first step.
    pub fn is_left_edge(&self) -> bool {
        self.read_offset == -1
    }

    /// `isRightEdge`: after the step that returned no operator.
    pub fn is_right_edge(&self) -> bool {
        self.read_offset == self.read_length
    }

    /// `getReadOffset`: the read base the machine is on, or the base before a deletion.
    pub fn read_offset(&self) -> i32 {
        self.read_offset
    }

    /// `getGenomeOffset`.
    pub fn genome_offset(&self) -> i32 {
        self.genome_offset
    }

    /// `getGenomePosition`: 1-based, and meaningless on an edge, which the reference notes and
    /// does not guard against.
    ///
    /// The reference adds the offset to `read.getStart()` (`AlignmentStateMachine.java:124`), which
    /// is 0 for an unmapped read; this adds it to the stored alignment start. The two agree for a
    /// mapped read and differ for an unmapped one that carries a position.
    pub fn genome_position(&self) -> i32 {
        self.alignment_start + self.genome_offset
    }

    /// `getCurrentCigarElementOffset`.
    pub fn current_cigar_element_offset(&self) -> i32 {
        self.current_cigar_element_offset
    }

    /// `getOffsetIntoCurrentCigarElement`: 0 on the right edge, -1 on the left.
    pub fn offset_into_current_cigar_element(&self) -> i32 {
        self.offset_into_current_cigar_element
    }

    /// `getCurrentCigarElement`: `None` on either edge.
    pub fn current_cigar_element(&self) -> Option<CigarElement> {
        self.current_element
    }

    /// `getCigarOperator`: `None` on either edge.
    pub fn cigar_operator(&self) -> Option<Op> {
        self.current_element.map(|e| e.op)
    }

    /// `stepForwardOnGenome`: advance until the next on-genome element (`M`, `X`, `=`, `D`, `N`).
    ///
    /// Returns `Ok(None)` when the machine steps off the right edge, which is the reference's
    /// null and not an error.
    ///
    /// **How**: the reference's `while ( true )` loop. Each pass either moves to the next element
    /// (when the current one is used up), consumes an off-genome element whole and loops, or
    /// advances one base on an on-genome element and returns its operator.
    pub fn step_forward_on_genome(&mut self) -> Result<Option<Op>, MalformedRead> {
        // `loop` is `while (true)`; it is left only by a `return`.
        loop {
            // The reference's test (`AlignmentStateMachine.java:232`): no element yet, or the
            // current one used up.
            let exhausted = match self.current_element {
                None => true,
                Some(element) => {
                    self.offset_into_current_cigar_element + 1 >= element.length as i32
                }
            };
            if exhausted {
                self.current_cigar_element_offset += 1;
                if (self.current_cigar_element_offset as usize) < self.elements.len() {
                    self.current_element =
                        Some(self.elements[self.current_cigar_element_offset as usize]);
                    // Re-entered rather than stepped, so that a zero-length element is skipped
                    // by the same test that ends a normal one.
                    self.offset_into_current_cigar_element = -1;
                    continue;
                }
                if self.current_element.map(|e| e.op) == Some(Op::D) {
                    return Err(MalformedRead::EndsWithDeletion);
                }
                self.offset_into_current_cigar_element = 0;
                self.read_offset = self.read_length;
                self.current_element = None;
                // A read ending on an indel models the next reference base after it, so the
                // genome offset advances one past the read even though nothing was consumed.
                // The increment is unconditional, so a read ending on a match takes it as well.
                self.genome_offset += 1;
                return Ok(None);
            }

            self.offset_into_current_cigar_element += 1;
            // `expect` unwraps the `Option` and panics with that text if it is `None`. It cannot
            // be: the branch above either loaded an element or returned.
            let element = self.current_element.expect("checked above");
            // A `match` must name every operator, so the reference's unreachable `default:` throw
            // has no counterpart: an operator missing here would not compile.
            match element.op {
                // Hard clips and pads are consumed whole and consume no read bases.
                Op::H | Op::P => {
                    self.offset_into_current_cigar_element = element.length as i32;
                }
                // Insertions and soft clips are consumed whole and do consume read bases.
                Op::I | Op::S => {
                    self.offset_into_current_cigar_element = element.length as i32;
                    self.read_offset += element.length as i32;
                }
                Op::D => {
                    if self.read_offset < 0 {
                        return Err(MalformedRead::StartsWithDeletion);
                    }
                    self.genome_offset += 1;
                    return Ok(Some(element.op));
                }
                // A reference skip moves like a deletion and is not refused at either end.
                Op::N => {
                    self.genome_offset += 1;
                    return Ok(Some(element.op));
                }
                Op::M | Op::Eq | Op::X => {
                    self.read_offset += 1;
                    self.genome_offset += 1;
                    return Ok(Some(element.op));
                }
            }
        }
    }
}
