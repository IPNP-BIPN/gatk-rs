//! Ported from `org.broadinstitute.hellbender.utils.read.CigarUtils` (GATK 4.6.2.0): the clipping
//! half, and `calculateCigar`.
//!
//! The clipping functions are what `ReadClipper` runs on: they decide the cigar of every clipped
//! read, and therefore bytes in every output BAM a clipping tool writes. `calculateCigar` is what
//! the assembler runs on: it realigns a haplotype to the reference with Smith-Waterman, between two
//! pads of `N` so that the alignment cannot clip, then trims the pads and left-aligns the indels.

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

use crate::alignment_utils::{left_align_indels, trim_cigar_by_bases, AlignmentError};
use crate::cigar_builder::{CigarBuilder, CigarError};
use crate::smith_waterman::{SmithWatermanAligner, SwOverhangStrategy, SwParameters};

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

/// `CigarUtils.SW_PAD`: ten `N`, on both sides of both sequences.
const SW_PAD: &[u8] = b"NNNNNNNNNN";

/// `CigarUtils.calculateCigar`: the cigar of `alt_seq` against `ref_seq`, or `None` when the
/// aligner produced something unusable (an alignment starting past the first base, or a soft clip).
///
/// The two sequences are assumed to be anchored at both ends, as haplotypes of one assembly region
/// are, so there is no start or end coordinate. A deletion at either end is kept: the alt haplotype
/// must keep its reference span, as when the reference starts with N repeats of a long sequence and
/// the alt with N-1.
///
/// Trailing deletions that the padding trim removes are put back before left-alignment and may be
/// removed again by it; leading ones are handed to the left-aligner as a read start and put back as
/// a leading deletion afterwards.
pub fn calculate_cigar(
    ref_seq: &[u8],
    alt_seq: &[u8],
    aligner: &dyn SmithWatermanAligner,
    haplotype_to_reference_sw_parameters: &SwParameters,
    strategy: SwOverhangStrategy,
) -> Result<Option<Cigar>, AlignmentError> {
    if alt_seq.is_empty() {
        // Horrible edge case from the unit tests, where this path has no bases.
        return Ok(Some(Cigar::new(vec![CigarElement {
            length: ref_seq.len() as u32,
            op: Op::D,
        }])));
    }

    // Equal strings are trivial, and the check is O(n).
    if ref_seq == alt_seq {
        return Ok(Some(Cigar::new(vec![CigarElement {
            length: ref_seq.len() as u32,
            op: Op::M,
        }])));
    }

    let padded_ref = [SW_PAD, ref_seq, SW_PAD].concat();
    let padded_path = [SW_PAD, alt_seq, SW_PAD].concat();
    let alignment = aligner.align(
        &padded_ref,
        &padded_path,
        haplotype_to_reference_sw_parameters,
        strategy,
    )?;

    if is_sw_failure(&alignment.cigar, alignment.alignment_offset) {
        return Ok(None);
    }

    // Cut off the padding bases. The end is inclusive.
    let base_start = SW_PAD.len() as i32;
    let base_end = padded_path.len() as i32 - SW_PAD.len() as i32 - 1;
    let trimmed = trim_cigar_by_bases(&alignment.cigar, base_start, base_end)?;
    let mut non_standard = trimmed.cigar;

    // Leading deletions removed by the trim shift the alignment start to the right.
    let trimmed_leading_deletions = trimmed.leading_deletion_bases_removed;
    // Trailing ones are kept, in order to left-align.
    let trimmed_trailing_deletions = trimmed.trailing_deletion_bases_removed;
    if trimmed_trailing_deletions > 0 {
        non_standard.elements.push(CigarElement {
            length: trimmed_trailing_deletions,
            op: Op::D,
        });
    }

    let left_aligned = left_align_indels(
        &non_standard,
        ref_seq,
        alt_seq,
        trimmed_leading_deletions as i32,
    )?;

    // Leading deletions removed when trimming the padding and when left-aligning; trailing ones
    // were restored for the left-alignment, which may have removed them again.
    let total_leading = trimmed_leading_deletions + left_aligned.leading_deletion_bases_removed;
    let total_trailing = left_aligned.trailing_deletion_bases_removed;

    if total_leading == 0 && total_trailing == 0 {
        return Ok(Some(left_aligned.cigar));
    }
    let mut result = Vec::new();
    if total_leading > 0 {
        result.push(CigarElement {
            length: total_leading,
            op: Op::D,
        });
    }
    result.extend(left_aligned.cigar.elements.iter().copied());
    if total_trailing > 0 {
        result.push(CigarElement {
            length: total_trailing,
            op: Op::D,
        });
    }
    Ok(Some(Cigar::new(result)))
}

/// `isSWFailure`: the alignment must start at the first base, which the padding guarantees, and
/// must carry no soft clip, which downstream code cannot take.
fn is_sw_failure(cigar: &Cigar, alignment_offset: i32) -> bool {
    alignment_offset > 0 || cigar.elements.iter().any(|element| element.op == Op::S)
}
