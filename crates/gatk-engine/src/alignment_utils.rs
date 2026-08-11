//! Ported from `org.broadinstitute.hellbender.utils.read.AlignmentUtils` (GATK 4.6.2.0):
//! `leftAlignIndels`, `normalizeAlleles`, and the `IndexRange` the two of them push around.
//!
//! One indel can often sit at several positions and still mean the same haplotype. The convention
//! is the leftmost, and this is the function that enforces it: a cigar, a reference, a read and a
//! start in; a cigar and the count of leading deletion bases the rewrite removed out.
//!
//! # Two indels can cancel each other
//!
//! The cigar is walked **right to left**. An indel's reference and read ranges are accumulated and
//! not resolved until an alignment block or the start of the cigar is reached, so two indels with
//! too few matching bases between them are merged into one. Measured: `3M2I2M2D3M` over a
//! homopolymer comes back as `10M`, the insertion and the deletion having met and cancelled, and
//! `3M1D2M1D3M` comes back as `8M` with two reference bases removed.
//!
//! # A deletion that reaches the start is dropped, and the caller must move the read
//!
//! [`crate::cigar_builder::CigarBuilder`] drops a deletion that ends up leading and reports the
//! reference bases it removed. A plain homopolymer deletion is one of these: `4M1D3M` becomes `7M`
//! with one base removed, and the caller is expected to move the read right by one. Seven of the
//! twelve measured cases end this way, which makes it the ordinary outcome rather than the corner.
//!
//! # `normalize_alleles` can shift **right**
//!
//! It trims shared bases off the end of the alleles and then off the front, and the front trim
//! decrements the start shift. When the alleles end differently, the end trim stops at once, the
//! front trim still runs, and nothing can shift left afterwards: the function returns `-1`, which
//! [`left_align_indels`] reads as new matching bases on the left. The shift is signed for that
//! reason and not because signedness was inherited from Java.
//!
//! # Soft and hard clips are not alignment blocks
//!
//! An indel may not be shifted into one. The reference's own javadoc example is measured: `2S2M2I`
//! over `GGAA` and `TTAAAA` at read start 2 is `2S2I2M`, which also says that soft-clipped bases
//! are present in the read's byte array while hard-clipped ones are not.

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

use crate::cigar_builder::{CigarBuilder, CigarError};

/// `org.broadinstitute.hellbender.utils.IndexRange`, with the six shifts these two functions use.
///
/// A half-open range over a sequence, `[start, end)`. It is mutated in place because
/// `normalizeAlleles` adjusts what it is given and its caller reads the ranges afterwards rather
/// than the return value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexRange {
    pub start: i32,
    pub end: i32,
}

impl IndexRange {
    pub fn new(start: i32, end: i32) -> Self {
        IndexRange { start, end }
    }

    pub fn size(&self) -> i32 {
        self.end - self.start
    }

    /// `shift`: both ends.
    pub fn shift(&mut self, by: i32) {
        self.start += by;
        self.end += by;
    }

    /// `shiftLeft`.
    pub fn shift_left(&mut self, by: i32) {
        self.shift(-by);
    }

    /// `shiftStart`.
    pub fn shift_start(&mut self, by: i32) {
        self.start += by;
    }

    /// `shiftStartLeft`.
    pub fn shift_start_left(&mut self, by: i32) {
        self.start -= by;
    }

    /// `shiftEndLeft`.
    pub fn shift_end_left(&mut self, by: i32) {
        self.end -= by;
    }
}

/// Why the reference refused, rather than clamping.
///
/// Both messages are `IllegalArgumentException` out of a util, not `UserException` out of a tool,
/// so a tool that calls this does not get to dress them up.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlignmentError {
    /// `ParamUtils.isPositiveOrZero`.
    NegativeReadStart,
    /// `Utils.validateArg(necessaryRefLength <= ref.length, ...)`.
    PastTheReference,
    /// `Utils.validateArg(readIndelRange.getStart() == 0, ...)`.
    CigarMissesReadBases,
    /// `normalizeAlleles`' own preconditions.
    BadAlleleRanges(&'static str),
    /// A cigar the builder itself refuses, which this function does not produce on its own.
    Builder(CigarError),
}

impl AlignmentError {
    /// The reference's own message.
    pub fn message(&self) -> String {
        match self {
            AlignmentError::NegativeReadStart => {
                "read start within reference base array must be non-negative".to_string()
            }
            AlignmentError::PastTheReference => "read goes past end of reference".to_string(),
            AlignmentError::CigarMissesReadBases => {
                "Given cigar does not account for all bases of the read".to_string()
            }
            AlignmentError::BadAlleleRanges(message) => (*message).to_string(),
            AlignmentError::Builder(error) => format!("{error:?}"),
        }
    }

    pub fn class(&self) -> &'static str {
        match self {
            AlignmentError::Builder(_) => "java.lang.IllegalStateException",
            _ => "java.lang.IllegalArgumentException",
        }
    }
}

/// `CigarBuilder.Result`: the rewritten cigar and what the rewrite removed from its ends.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeftAlignment {
    pub cigar: Cigar,
    /// What the caller must move the read's start by.
    pub leading_deletion_bases_removed: u32,
    pub trailing_deletion_bases_removed: u32,
}

/// `CigarOperator.isIndel`.
fn is_indel(op: Op) -> bool {
    matches!(op, Op::I | Op::D)
}

/// `CigarOperator.isAlignment`: the three operators that line a read base up against a reference
/// base, which are the only ones an indel may be shifted into.
fn is_alignment(op: Op) -> bool {
    matches!(op, Op::M | Op::Eq | Op::X)
}

fn length_on_read(element: &CigarElement) -> i32 {
    if element.op.consumes_read_bases() {
        element.length as i32
    } else {
        0
    }
}

fn length_on_reference(element: &CigarElement) -> i32 {
    if element.op.consumes_reference_bases() {
        element.length as i32
    } else {
        0
    }
}

fn cigar_element(length: i32, op: Op) -> CigarElement {
    CigarElement {
        length: length.max(0) as u32,
        op,
    }
}

/// `leftAlignIndels`: the same alignment, with every indel at the leftmost position that still
/// represents it.
///
/// `read_start` is the **0-based** position on `ref` of the read's first aligned base. The tool
/// passes 0 because its reference context starts at the read; a haplotype caller passes something
/// else.
pub fn left_align_indels(
    cigar: &Cigar,
    reference: &[u8],
    read: &[u8],
    read_start: i32,
) -> Result<LeftAlignment, AlignmentError> {
    if read_start < 0 {
        return Err(AlignmentError::NegativeReadStart);
    }

    let elements = &cigar.elements;
    if !elements.iter().any(|element| is_indel(element.op)) {
        return Ok(LeftAlignment {
            cigar: cigar.clone(),
            leading_deletion_bases_removed: 0,
            trailing_deletion_bases_removed: 0,
        });
    }

    // We need reference bases from the start of the read to the rightmost indel, and no further.
    let last_indel = elements
        .iter()
        .rposition(|element| is_indel(element.op))
        .expect("an indel, because the scan above found one");
    let necessary_ref_length: i32 = read_start
        + elements[..=last_indel]
            .iter()
            .map(length_on_reference)
            .sum::<i32>();
    if necessary_ref_length > reference.len() as i32 {
        return Err(AlignmentError::PastTheReference);
    }

    // One base past the end of the read, and then right to left.
    let mut result_right_to_left: Vec<CigarElement> = Vec::new();
    let reference_length = cigar.reference_length() as i32;
    let mut ref_indel =
        IndexRange::new(read_start + reference_length, read_start + reference_length);
    let mut read_indel = IndexRange::new(read.len() as i32, read.len() as i32);

    for n in (0..elements.len()).rev() {
        let element = &elements[n];
        if is_indel(element.op) {
            // Accumulate. The indel is not shifted until an alignment block or the read start.
            ref_indel.shift_start_left(length_on_reference(element));
            read_indel.shift_start_left(length_on_read(element));
        } else if ref_indel.size() == 0 && read_indel.size() == 0 {
            result_right_to_left.push(*element);
            ref_indel.shift_left(length_on_reference(element));
            read_indel.shift_left(length_on_read(element));
        } else {
            // An indel to left-align, into this block if it is an alignment block and not at all
            // if it is a clip.
            let max_shift = if is_alignment(element.op) {
                element.length as i32
            } else {
                0
            };
            let mut bounds = [ref_indel, read_indel];
            let (start_shift, end_shift) =
                normalize_alleles(&[reference, read], &mut bounds, max_shift, true)?;
            ref_indel = bounds[0];
            read_indel = bounds[1];

            // New match alignments on the right, from having moved left.
            result_right_to_left.push(cigar_element(end_shift, Op::M));

            // `n == 0`: an indel at the very start of the cigar has nothing left to shift into.
            let emit_indel = n == 0 || start_shift < max_shift || !is_alignment(element.op);
            // We may have shifted RIGHT to make the alleles parsimonious.
            let new_match_on_left = if start_shift < 0 { -start_shift } else { 0 };
            let remaining_on_left = if start_shift < 0 {
                element.length as i32
            } else {
                element.length as i32 - start_shift
            };

            if emit_indel {
                result_right_to_left.push(cigar_element(ref_indel.size(), Op::D));
                result_right_to_left.push(cigar_element(read_indel.size(), Op::I));
                // Both ranges are now empty and point at the start of the left-aligned indel.
                ref_indel.shift_end_left(ref_indel.size());
                read_indel.shift_end_left(read_indel.size());

                ref_indel.shift_left(
                    new_match_on_left
                        + if element.op.consumes_reference_bases() {
                            remaining_on_left
                        } else {
                            0
                        },
                );
                read_indel.shift_left(
                    new_match_on_left
                        + if element.op.consumes_read_bases() {
                            remaining_on_left
                        } else {
                            0
                        },
                );
            }
            result_right_to_left.push(cigar_element(new_match_on_left, Op::M));
            result_right_to_left.push(cigar_element(remaining_on_left, element.op));
        }
    }

    // Indels at the start of the cigar, which had no non-indel element to their left.
    result_right_to_left.push(cigar_element(ref_indel.size(), Op::D));
    result_right_to_left.push(cigar_element(read_indel.size(), Op::I));

    if read_indel.start != 0 {
        return Err(AlignmentError::CigarMissesReadBases);
    }

    let mut builder = CigarBuilder::default();
    for element in result_right_to_left.iter().rev() {
        builder.add(*element).map_err(AlignmentError::Builder)?;
    }
    // `makeAndRecordDeletionsRemovedResult`, whose `make()` is `make(false)`.
    let cigar = builder.make(false).map_err(AlignmentError::Builder)?;
    Ok(LeftAlignment {
        cigar,
        leading_deletion_bases_removed: builder.leading_deletion_bases_removed(),
        trailing_deletion_bases_removed: builder.trailing_deletion_bases_removed(),
    })
}

/// `normalizeAlleles`: move a set of allele ranges as far left as they can go and still mean the
/// same event, trimming redundant shared bases first when asked.
///
/// Returns the start and end shifts, and **adjusts `bounds` in place**, which is how the caller
/// reads the new positions. The start shift is signed: trimming shared bases off the front moves
/// the range right, and a range that cannot then shift left keeps that negative shift.
pub fn normalize_alleles(
    sequences: &[&[u8]],
    bounds: &mut [IndexRange],
    max_shift: i32,
    trim: bool,
) -> Result<(i32, i32), AlignmentError> {
    if sequences.is_empty() {
        return Err(AlignmentError::BadAlleleRanges(
            "sequences must not be empty",
        ));
    }
    if sequences.len() != bounds.len() {
        return Err(AlignmentError::BadAlleleRanges(
            "Must have one initial allele range per sequence",
        ));
    }
    if bounds.iter().any(|bound| max_shift > bound.start) {
        return Err(AlignmentError::BadAlleleRanges(
            "maxShift goes past the start of a sequence",
        ));
    }

    let mut start_shift = 0;
    let mut end_shift = 0;

    let mut min_size = bounds
        .iter()
        .map(IndexRange::size)
        .min()
        .expect("a bound, because the lists are non-empty and the same length");

    // Redundant shared bases at the end of the alleles.
    while trim && min_size > 0 && last_base_on_right_is_same(sequences, bounds) {
        bounds.iter_mut().for_each(|bound| bound.shift_end_left(1));
        min_size -= 1;
        end_shift += 1;
    }

    // And at the front, which moves the range RIGHT.
    while trim && min_size > 0 && first_base_on_left_is_same(sequences, bounds) {
        bounds.iter_mut().for_each(|bound| bound.shift_start(1));
        min_size -= 1;
        start_shift -= 1;
    }

    // Then left as long as the bases on both sides agree across every sequence. An empty range's
    // last base on the right is its next base on the left, which is what makes an insertion's
    // reference range work.
    while start_shift < max_shift
        && next_base_on_left_is_same(sequences, bounds)
        && last_base_on_right_is_same(sequences, bounds)
    {
        bounds.iter_mut().for_each(|bound| bound.shift_left(1));
        start_shift += 1;
        end_shift += 1;
    }

    Ok((start_shift, end_shift))
}

/// The base at `at`, or `None` where Java would have thrown out of bounds. The callers' guards
/// make that unreachable: a trim loop runs only while the ranges are non-empty, and the shift loop
/// only while `start_shift < max_shift`, with `max_shift <= start` a precondition.
fn base(sequence: &[u8], at: i32) -> Option<u8> {
    usize::try_from(at)
        .ok()
        .and_then(|at| sequence.get(at))
        .copied()
}

fn same_base_at(
    sequences: &[&[u8]],
    bounds: &[IndexRange],
    at: impl Fn(&IndexRange) -> i32,
) -> bool {
    let Some(first) = base(sequences[0], at(&bounds[0])) else {
        return false;
    };
    sequences
        .iter()
        .zip(bounds)
        .all(|(sequence, bound)| base(sequence, at(bound)) == Some(first))
}

fn last_base_on_right_is_same(sequences: &[&[u8]], bounds: &[IndexRange]) -> bool {
    same_base_at(sequences, bounds, |bound| bound.end - 1)
}

fn first_base_on_left_is_same(sequences: &[&[u8]], bounds: &[IndexRange]) -> bool {
    same_base_at(sequences, bounds, |bound| bound.start)
}

fn next_base_on_left_is_same(sequences: &[&[u8]], bounds: &[IndexRange]) -> bool {
    same_base_at(sequences, bounds, |bound| bound.start - 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::text_parse::parse_cigar;

    fn aligned(cigar: &str, reference: &str, read: &str, start: i32) -> LeftAlignment {
        left_align_indels(
            &parse_cigar(cigar).unwrap(),
            reference.as_bytes(),
            read.as_bytes(),
            start,
        )
        .unwrap()
    }

    #[test]
    fn a_cigar_with_no_indel_comes_back_untouched() {
        let result = aligned("6M", "AAAAAA", "AAAAAA", 0);
        assert_eq!(result.cigar.to_text(), "6M");
        assert_eq!(result.leading_deletion_bases_removed, 0);
    }

    #[test]
    fn a_homopolymer_deletion_walks_off_the_start_and_moves_the_read() {
        let result = aligned("4M1D3M", "AAAAAAAT", "AAAAAAT", 0);
        assert_eq!(result.cigar.to_text(), "7M");
        assert_eq!(result.leading_deletion_bases_removed, 1);
    }

    #[test]
    fn an_insertion_and_a_deletion_that_meet_cancel() {
        let result = aligned("3M2I2M2D3M", "AAAAAAAAT", "AAAAAAAAAT", 0);
        assert_eq!(result.cigar.to_text(), "10M");
        assert_eq!(result.leading_deletion_bases_removed, 0);
    }

    #[test]
    fn an_indel_is_not_shifted_into_a_soft_clip() {
        // The reference's own javadoc example.
        let result = aligned("2S2M2I", "GGAA", "TTAAAA", 2);
        assert_eq!(result.cigar.to_text(), "2S2I2M");
    }

    #[test]
    fn the_two_refusals_are_the_reference_s() {
        let past =
            left_align_indels(&parse_cigar("4M1D3M").unwrap(), b"AAAA", b"AAAAAAT", 0).unwrap_err();
        assert_eq!(past.message(), "read goes past end of reference");
        assert_eq!(past.class(), "java.lang.IllegalArgumentException");

        let misses = left_align_indels(
            &parse_cigar("4M1D3M").unwrap(),
            b"AAAAAAAAAAAA",
            b"AAAAAAAAAAAA",
            0,
        )
        .unwrap_err();
        assert_eq!(
            misses.message(),
            "Given cigar does not account for all bases of the read"
        );
    }

    /// The shift that goes the other way, which is why it is signed.
    #[test]
    fn trimming_the_front_returns_a_negative_start_shift() {
        let mut bounds = [IndexRange::new(1, 3), IndexRange::new(1, 3)];
        let shifts = normalize_alleles(&[b"CAG", b"CAT"], &mut bounds, 1, true).unwrap();
        assert_eq!(shifts, (-1, 0));
        assert_eq!(bounds, [IndexRange::new(2, 3), IndexRange::new(2, 3)]);
    }

    #[test]
    fn max_shift_is_a_wall() {
        let mut bounds = [IndexRange::new(5, 5), IndexRange::new(5, 6)];
        let shifts = normalize_alleles(&[b"GAAAAT", b"GAAAAAT"], &mut bounds, 0, true).unwrap();
        assert_eq!(shifts, (0, 0), "nothing may move");
        assert_eq!(bounds, [IndexRange::new(5, 5), IndexRange::new(5, 6)]);
    }
}
