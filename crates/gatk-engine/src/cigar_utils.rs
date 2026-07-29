//! Ported from `org.broadinstitute.hellbender.utils.read.CigarUtils` (GATK 4.6.2.0), the clipping
//! half of it.
//!
//! These three functions are what `ReadClipper` runs on: they decide the cigar of every clipped
//! read, and therefore bytes in every output BAM a clipping tool writes.

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

use crate::cigar_builder::{CigarBuilder, CigarError};

/// `CigarUtils.countRefBasesAndClips`: reference-consuming elements plus **both** kinds of clip.
///
/// The distinction from `countRefBasesAndSoftClips`, which excludes hard clips, is what makes the
/// difference between the span of the read as sequenced and the span of what is left of it. This
/// one is used to place a read that has already lost bases to hard clipping.
pub fn count_ref_bases_and_clips(elements: &[CigarElement]) -> i32 {
    elements
        .iter()
        .filter(|e| e.op.consumes_reference_bases() || e.op == Op::S || e.op == Op::H)
        .map(|e| e.length as i32)
        .sum()
}

/// `CigarUtils.clipCigar(cigar, start, stop, clippingOperator)`.
///
/// `start` is inclusive and `stop` exclusive, both in **read** coordinates, and `start == 0` is
/// what makes this a left clip rather than a right one.
///
/// Two rules that look like details and are not:
///
///  * hard clips already in the cigar are copied through untouched, before anything else, so they
///    do not count towards the read coordinates the clip is expressed in;
///  * a deletion that sits exactly at the clip boundary is **dropped**, not clipped, because a
///    deletion at the edge of a clip describes nothing. That is the `elementStart != start &&
///    elementStart != stop` clause, and it is the reason `4M1I1D5M` clipped at 0..1 comes back as
///    `1S3M1D1I5M` rather than with the deletion where it was.
pub fn clip_cigar(
    cigar: &Cigar,
    start: i32,
    stop: i32,
    clipping_operator: Op,
) -> Result<Cigar, CigarError> {
    let clip_left = start == 0;
    let mut builder = CigarBuilder::default();

    let mut element_start = 0;
    for element in &cigar.elements {
        let operator = element.op;
        if operator == Op::H {
            builder.add(*element)?;
            continue;
        }
        let element_end = element_start
            + if operator.consumes_read_bases() {
                element.length as i32
            } else {
                0
            };

        if element_end <= start || element_start >= stop {
            // Outside the clipped span: copied, unless it is a deletion sitting on the boundary.
            if operator.consumes_read_bases() || (element_start != start && element_start != stop) {
                builder.add(*element)?;
            }
        } else {
            let unclipped_length = if clip_left {
                element_end - stop
            } else {
                start - element_start
            };
            let clipped_length = element.length as i32 - unclipped_length;

            if unclipped_length <= 0 {
                // Entirely inside the clip: an element that consumes read bases becomes clipping,
                // and one that does not simply disappears.
                if operator.consumes_read_bases() {
                    builder.add(CigarElement {
                        length: element.length,
                        op: clipping_operator,
                    })?;
                }
            } else if clip_left {
                builder.add(CigarElement {
                    length: clipped_length as u32,
                    op: clipping_operator,
                })?;
                builder.add(CigarElement {
                    length: unclipped_length as u32,
                    op: operator,
                })?;
            } else {
                builder.add(CigarElement {
                    length: unclipped_length as u32,
                    op: operator,
                })?;
                builder.add(CigarElement {
                    length: clipped_length as u32,
                    op: clipping_operator,
                })?;
            }
        }
        element_start = element_end;
    }

    builder.make(false)
}

/// `CigarUtils.alignmentStartShift`: how far the alignment start moves when the first `numClipped`
/// read bases are clipped away.
///
/// Hard clips are skipped outright, and a deletion immediately following the clipped span counts
/// as clipped too, which is what "this includes deletions immediately following clipping" in the
/// reference means.
pub fn alignment_start_shift(cigar: &Cigar, num_clipped: i32) -> i32 {
    let mut ref_bases_clipped = 0;
    let mut element_start = 0;
    for element in &cigar.elements {
        let operator = element.op;
        if operator == Op::H {
            continue;
        }
        let element_end = element_start
            + if operator.consumes_read_bases() {
                element.length as i32
            } else {
                0
            };
        if element_end <= num_clipped {
            if operator.consumes_reference_bases() {
                ref_bases_clipped += element.length as i32;
            }
        } else if element_start < num_clipped {
            // The clip lands inside this element, which therefore consumes read bases.
            let clipped_length = num_clipped - element_start;
            if operator.consumes_reference_bases() {
                ref_bases_clipped += clipped_length;
            }
            break;
        }
        element_start = element_end;
    }
    ref_bases_clipped
}

/// `CigarUtils.revertSoftClips`: every soft clip becomes an M, through the builder.
///
/// Through the builder is the whole point: `3S7M` does not become `3M7M`, it becomes `10M`.
pub fn revert_soft_clips(cigar: &Cigar) -> Result<Cigar, CigarError> {
    let mut builder = CigarBuilder::default();
    for element in &cigar.elements {
        builder.add(CigarElement {
            length: element.length,
            op: if element.op == Op::S {
                Op::M
            } else {
                element.op
            },
        })?;
    }
    builder.make(false)
}
