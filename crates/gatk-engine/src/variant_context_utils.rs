//! `GATKVariantContextUtils.leftAlignAndTrim`, ported from GATK 4.6.2.0.
//!
//! The workhorse under `LeftAlignAndTrimVariants`: it slides an indel as far left as the reference
//! lets it. [`crate::alignment_utils::normalize_alleles`] does the sliding; what is here is the
//! window that keeps widening around it, and what the record becomes afterwards.
//!
//! # The window widens, and it can stop short without saying so
//!
//! ```java
//! for (int leadingBases = Math.min(maxLeadingBases, 10); leadingBases <= maxLeadingBases; leadingBases = Math.min(2*leadingBases, maxLeadingBases)) {
//! ...
//! } else if (shifts.getLeft() == variantOffsetInRef && leadingBases < maxLeadingBases) {
//!     continue;
//! }
//! ```
//!
//! Ten bases first, then twenty, then forty, until the shift stops landing on the edge of the
//! slice. When the shift **does** land on the edge and the window is already at its maximum, the
//! loop falls through and the partly-aligned record is returned: the same deletion comes out at 15
//! with a window of 2, at 49 with 10, at 39 with 20, and at the start of its run with a wide
//! enough one. Nothing in the record says which of those happened.
//!
//! # What comes back unchanged
//!
//! A non-indel, a `maxLeadingBases` of zero or less, and a shift of zero all return the input
//! itself. The first is a **type** test, so a MIXED site is never aligned; the second is what the
//! caller produces for two adjacent variants, `min(maxLeadingBases, distanceToLastVariant - 1)`.

use crate::alignment_utils::{normalize_alleles, AlignmentError, IndexRange};

/// One allele: its bases and whether it is the reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allele {
    pub bases: Vec<u8>,
    pub is_reference: bool,
}

impl Allele {
    pub fn new(bases: &[u8], is_reference: bool) -> Allele {
        Allele {
            bases: bases.to_vec(),
            is_reference,
        }
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }
}

/// As much of a `VariantContext` as this function reads and writes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub contig: String,
    pub start: i32,
    pub stop: i32,
    /// The reference allele first, then the alternates.
    pub alleles: Vec<Allele>,
    /// One list of allele indices per sample, which the map rewrites.
    pub genotypes: Vec<Vec<usize>>,
}

impl Variant {
    /// `isIndel()`: the type is INDEL only when every alternate differs from the reference in
    /// length in the same way. Two alternates of different types are MIXED, which is not an indel.
    ///
    /// The same rule as [`crate::pileup_summary`]'s neighbours; it lives here because this is what
    /// decides whether the alignment runs at all.
    pub fn is_indel(&self) -> bool {
        if self.alleles.len() < 2 {
            return false;
        }
        let reference = &self.alleles[0];
        let mut kind: Option<bool> = None;
        for allele in &self.alleles[1..] {
            // An alternate is "indel-shaped" when its length differs from the reference's.
            let this = allele.len() != reference.len();
            match kind {
                None => kind = Some(this),
                Some(seen) if seen != this => return false,
                Some(_) => {}
            }
        }
        kind.unwrap_or(false)
    }
}

/// `leftAlignAndTrim(vc, ref, maxLeadingBases, trim)`.
///
/// `reference_bases` is the whole contig, one-based positions indexing from 1, which is what a
/// `ReferenceContext` hands back a slice of.
pub fn left_align_and_trim(
    variant: &Variant,
    reference_bases: &[u8],
    max_leading_bases: i32,
    trim: bool,
) -> Result<Variant, AlignmentError> {
    if !variant.is_indel() || max_leading_bases <= 0 {
        return Ok(variant.clone());
    }

    let mut leading_bases = max_leading_bases.min(10);
    loop {
        let ref_start = (variant.start - leading_bases).max(1);
        // The slice the reference context would return: from `ref_start` to the variant's end.
        let slice = &reference_bases[(ref_start - 1) as usize..variant.stop as usize];
        let variant_offset = variant.start - ref_start;

        // Each allele preceded by the reference bases before it, which is what gives the shift its
        // room. The trailing reference is NOT appended: the sequences end where the allele does.
        let sequences: Vec<Vec<u8>> = variant
            .alleles
            .iter()
            .map(|allele| {
                let mut sequence = slice[..variant_offset as usize].to_vec();
                sequence.extend_from_slice(&allele.bases);
                sequence
            })
            .collect();
        // `+1` to ignore the base shared by every allele.
        let mut ranges: Vec<IndexRange> = variant
            .alleles
            .iter()
            .map(|allele| IndexRange::new(variant_offset + 1, variant_offset + allele.len() as i32))
            .collect();

        let borrowed: Vec<&[u8]> = sequences
            .iter()
            .map(|sequence| sequence.as_slice())
            .collect();
        let (left, right) = normalize_alleles(&borrowed, &mut ranges, variant_offset, trim)?;

        if left == 0 && right == 0 {
            return Ok(variant.clone());
        }
        // The shift reached the edge of the slice and a wider one is allowed: try again.
        if left == variant_offset && leading_bases < max_leading_bases {
            leading_bases = (2 * leading_bases).min(max_leading_bases);
            continue;
        }

        let alleles: Vec<Allele> = variant
            .alleles
            .iter()
            .enumerate()
            .map(|(index, allele)| {
                let from = (variant_offset - left) as usize;
                let to = (variant_offset - right) as usize + allele.len();
                Allele::new(&sequences[index][from..to], index == 0)
            })
            .collect();

        return Ok(Variant {
            contig: variant.contig.clone(),
            start: variant.start - left,
            stop: variant.stop - right,
            // The reference builds this list from a HashMap's values, so the order is the map's.
            // Rebuilding it in record order is what the golden shows, and all it shows.
            alleles,
            genotypes: variant.genotypes.clone(),
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Ten G, eight A, twelve G, thirty T, ten G, a CA repeat, ten G: the dump's reference.
    fn reference() -> Vec<u8> {
        let mut bases = Vec::new();
        bases.extend(std::iter::repeat_n(b'G', 10));
        bases.extend(std::iter::repeat_n(b'A', 8));
        bases.extend(std::iter::repeat_n(b'G', 12));
        bases.extend(std::iter::repeat_n(b'T', 30));
        bases.extend(std::iter::repeat_n(b'G', 10));
        for _ in 0..10 {
            bases.extend_from_slice(b"CA");
        }
        bases.extend(std::iter::repeat_n(b'G', 10));
        bases
    }

    fn deletion(start: i32, reference: &str, alternate: &str) -> Variant {
        Variant {
            contig: "chr1".to_string(),
            start,
            stop: start + reference.len() as i32 - 1,
            alleles: vec![
                Allele::new(reference.as_bytes(), true),
                Allele::new(alternate.as_bytes(), false),
            ],
            genotypes: Vec::new(),
        }
    }

    #[test]
    fn a_narrow_window_stops_short_and_says_nothing() {
        let bases = reference();
        let variant = deletion(17, "AA", "A");
        let wide = left_align_and_trim(&variant, &bases, 1000, true).expect("aligned");
        assert_eq!((wide.start, wide.stop), (10, 11));
        let narrow = left_align_and_trim(&variant, &bases, 2, true).expect("aligned");
        assert_eq!((narrow.start, narrow.stop), (15, 16));
    }

    #[test]
    fn the_window_doubles_until_it_is_wide_enough() {
        let bases = reference();
        let variant = deletion(59, "TT", "T");
        for (window, start) in [(10, 49), (20, 39), (1000, 30)] {
            let result = left_align_and_trim(&variant, &bases, window, true).expect("aligned");
            assert_eq!(result.start, start, "window {window}");
        }
    }

    #[test]
    fn three_things_come_back_unchanged() {
        let bases = reference();
        let aligned = deletion(10, "GA", "G");
        assert_eq!(
            left_align_and_trim(&aligned, &bases, 1000, true).expect("same"),
            aligned
        );

        let variant = deletion(17, "AA", "A");
        assert_eq!(
            left_align_and_trim(&variant, &bases, 0, true).expect("same"),
            variant
        );

        // A snp is not an indel, and neither is a MIXED site.
        let snp = deletion(18, "A", "C");
        assert_eq!(
            left_align_and_trim(&snp, &bases, 1000, true).expect("same"),
            snp
        );
        let mut mixed = deletion(18, "AG", "A");
        mixed.alleles.push(Allele::new(b"CC", false));
        assert!(!mixed.is_indel());
        assert_eq!(
            left_align_and_trim(&mixed, &bases, 1000, true).expect("same"),
            mixed
        );
    }

    #[test]
    fn trimming_decides_whether_the_record_shrinks() {
        let bases = reference();
        // Alleles sharing a trailing base as well as the leading one.
        let variant = deletion(18, "AGG", "AG");
        let trimmed = left_align_and_trim(&variant, &bases, 1000, true).expect("aligned");
        assert_eq!((trimmed.start, trimmed.stop), (18, 19));
        let untrimmed = left_align_and_trim(&variant, &bases, 1000, false).expect("aligned");
        assert_eq!((untrimmed.start, untrimmed.stop), (17, 19));
    }
}
