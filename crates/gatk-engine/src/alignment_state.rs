//! Ported from `org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine`
//! (GATK 4.6.2.0).
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
//!  * **a read that ends on an insertion still advances the genome offset one past its end.** The
//!    reference does this deliberately, so the final position models the next reference base after
//!    the indel rather than the last base of the read;
//!  * **a cigar that starts or ends with a deletion is a malformed read**, and throws rather than
//!    being tolerated, even though the SAM specification permits it. Zero-length elements, by
//!    contrast, are skipped silently.

use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;

/// `UserException.MalformedRead`, raised by the two cigars this machine refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MalformedRead {
    /// The cigar's first on-genome element is a deletion.
    StartsWithDeletion,
    /// The cigar's last element is a deletion.
    EndsWithDeletion,
}

/// `AlignmentStateMachine`: a read, and where along the genome we are in it.
pub struct AlignmentStateMachine<'a> {
    elements: &'a [CigarElement],
    read_length: i32,
    alignment_start: i32,
    current_cigar_element_offset: i32,
    read_offset: i32,
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

    pub fn is_left_edge(&self) -> bool {
        self.read_offset == -1
    }

    pub fn is_right_edge(&self) -> bool {
        self.read_offset == self.read_length
    }

    pub fn read_offset(&self) -> i32 {
        self.read_offset
    }

    pub fn genome_offset(&self) -> i32 {
        self.genome_offset
    }

    /// `getGenomePosition`: 1-based, and meaningless on an edge, which the reference notes and
    /// does not guard against.
    pub fn genome_position(&self) -> i32 {
        self.alignment_start + self.genome_offset
    }

    pub fn current_cigar_element_offset(&self) -> i32 {
        self.current_cigar_element_offset
    }

    pub fn offset_into_current_cigar_element(&self) -> i32 {
        self.offset_into_current_cigar_element
    }

    pub fn current_cigar_element(&self) -> Option<CigarElement> {
        self.current_element
    }

    pub fn cigar_operator(&self) -> Option<Op> {
        self.current_element.map(|e| e.op)
    }

    /// `stepForwardOnGenome`: advance until the next on-genome element (`M`, `X`, `=`, `D`, `N`).
    ///
    /// Returns `Ok(None)` when the machine steps off the right edge, which is the reference's
    /// null and not an error.
    pub fn step_forward_on_genome(&mut self) -> Result<Option<Op>, MalformedRead> {
        loop {
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
                self.genome_offset += 1;
                return Ok(None);
            }

            self.offset_into_current_cigar_element += 1;
            let element = self.current_element.expect("checked above");
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
