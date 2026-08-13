//! `LeftAlignAndTrimVariants`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.variantutils.LeftAlignAndTrimVariants`
//! (GATK 4.6.2.0).
//!
//! The tool three engine bricks were built for. Its own logic is short: choose a window, call
//! [`gatk_engine::variant_context_utils::left_align_and_trim`], remember the record just written.
//!
//! # Aligning one record frees the next
//!
//! ```java
//! final int distanceToLastVariant = (lastVariant != null && splitVariant.contigsMatch(lastVariant)) ? splitVariant.getStart() - lastVariant.getEnd() : Integer.MAX_VALUE;
//! lastVariant = GATKVariantContextUtils.leftAlignAndTrim(splitVariant, ref, Math.min(maxLeadingBases, distanceToLastVariant - 1), !dontTrimAlleles);
//! ```
//!
//! `lastVariant` is the record **as written**, not as read. It has already moved left, so the
//! distance to the next one is larger than it was in the input: two indels a base apart become
//! eight apart by the time the second is measured, and the second moves after all. The bound is
//! there to stop two records crossing, and it loosens as records move out of the way.
//!
//! # A skipped record still bounds the next
//!
//! An indel longer than `--max-indel-length` is written untouched, and the assignment to
//! `lastVariant` happens on that branch too, so it bounds the record after it though it was never
//! aligned.

use gatk_engine::variant_context_utils::{
    left_align_and_trim, split_variant_context_to_biallelics, SplitError, Variant,
};

/// `DEFAULT_MAX_INDEL_SIZE`.
pub const DEFAULT_MAX_INDEL_SIZE: i32 = 200;
/// `DEFAULT_MAX_LEADING_BASES`.
pub const DEFAULT_MAX_LEADING_BASES: i32 = 1000;

/// The arguments this port reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Arguments {
    /// `--dont-trim-alleles`, which is passed down negated.
    pub dont_trim_alleles: bool,
    /// `--split-multi-allelics`.
    pub split_multiallelics: bool,
    /// `--max-indel-length`.
    pub max_indel_size: i32,
    /// `--max-leading-bases`.
    pub max_leading_bases: i32,
}

impl Default for Arguments {
    fn default() -> Arguments {
        Arguments {
            dont_trim_alleles: false,
            split_multiallelics: false,
            max_indel_size: DEFAULT_MAX_INDEL_SIZE,
            max_leading_bases: DEFAULT_MAX_LEADING_BASES,
        }
    }
}

/// `getIndelLengths` mapped through `abs` and maxed: the largest length change of any allele.
///
/// A record with a 300-base deletion and a one-base insertion is tested as 300, and a record with
/// no length change at all is tested as 0.
pub fn largest_indel_length(variant: &Variant) -> i32 {
    let reference = variant.alleles[0].len() as i32;
    variant.alleles[1..]
        .iter()
        .map(|allele| (allele.len() as i32 - reference).abs())
        .max()
        .unwrap_or(0)
}

/// The whole traversal: every record the tool writes, in order.
///
/// `reference_bases` is the whole of each contig, keyed by name, which is what a `ReferenceContext`
/// hands back a slice of.
pub fn left_align_and_trim_variants(
    variants: &[Variant],
    reference_bases: &dyn Fn(&str) -> Option<Vec<u8>>,
    arguments: Arguments,
) -> Result<Vec<Variant>, SplitError> {
    let mut written: Vec<Variant> = Vec::new();
    let mut last: Option<Variant> = None;

    for variant in variants {
        let pieces = if arguments.split_multiallelics {
            split_variant_context_to_biallelics(variant, false)?
        } else {
            vec![variant.clone()]
        };

        for piece in pieces {
            if largest_indel_length(&piece) > arguments.max_indel_size {
                // Written untouched, and it still becomes the record the next one is measured
                // against.
                last = Some(piece.clone());
                written.push(piece);
                continue;
            }

            // `Integer.MAX_VALUE` when there is no previous record on this contig, so the bound
            // only applies within a contig.
            let distance = match &last {
                Some(previous) if previous.contig == piece.contig => piece.start - previous.stop,
                _ => i32::MAX,
            };
            let window = arguments.max_leading_bases.min(distance - 1);
            let bases = reference_bases(&piece.contig).unwrap_or_default();
            let aligned = left_align_and_trim(&piece, &bases, window, !arguments.dont_trim_alleles)
                .map_err(SplitError::Alignment)?;
            last = Some(aligned.clone());
            written.push(aligned);
        }
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gatk_engine::variant_context_utils::Allele;

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

    fn variant(contig: &str, start: i32, reference: &str, alternates: &[&str]) -> Variant {
        let mut alleles = vec![Allele::new(reference.as_bytes(), true)];
        for alternate in alternates {
            alleles.push(Allele::new(alternate.as_bytes(), false));
        }
        Variant {
            contig: contig.to_string(),
            start,
            stop: start + reference.len() as i32 - 1,
            alleles,
            genotypes: Vec::new(),
            attributes: Vec::new(),
        }
    }

    fn bases(contig: &str) -> Option<Vec<u8>> {
        match contig {
            "chr1" => Some(reference()),
            _ => None,
        }
    }

    #[test]
    fn aligning_one_record_relaxes_the_bound_on_the_next() {
        // A base apart in the input; eight apart once the first has moved to 10-11.
        let variants = vec![
            variant("chr1", 17, "AA", &["A"]),
            variant("chr1", 19, "GG", &["G"]),
        ];
        let written =
            left_align_and_trim_variants(&variants, &bases, Arguments::default()).expect("written");
        assert_eq!(written[0].start, 10);
        assert_eq!(written[1].start, 18);
    }

    #[test]
    fn a_skipped_indel_still_bounds_the_next_record() {
        let long = variant("chr1", 34, &"T".repeat(26), &["T"]);
        let short = variant("chr1", 69, "GG", &["G"]);
        let variants = vec![long, short];

        let aligned =
            left_align_and_trim_variants(&variants, &bases, Arguments::default()).expect("written");
        assert_eq!(aligned[0].start, 30);

        let skipped = left_align_and_trim_variants(
            &variants,
            &bases,
            Arguments {
                max_indel_size: 5,
                ..Arguments::default()
            },
        )
        .expect("written");
        assert_eq!(skipped[0].start, 34, "written where it was");
        assert_eq!(skipped[0], variants[0], "and written untouched");
    }

    #[test]
    fn the_pieces_of_a_split_record_bound_each_other() {
        let variants = vec![variant("chr1", 17, "AA", &["A", "AAA"])];
        let written = left_align_and_trim_variants(
            &variants,
            &bases,
            Arguments {
                split_multiallelics: true,
                ..Arguments::default()
            },
        )
        .expect("written");
        assert_eq!(written.len(), 2);
        // The second piece is bounded by the first, which has already moved to 10-11.
        assert_eq!((written[0].start, written[1].start), (10, 12));
    }

    #[test]
    fn a_new_contig_is_not_bounded_at_all() {
        let variants = vec![
            variant("chr1", 59, "TT", &["T"]),
            variant("chr2", 19, "AA", &["A"]),
        ];
        let written = left_align_and_trim_variants(
            &variants,
            &|contig| match contig {
                "chr1" => Some(reference()),
                "chr2" => Some(b"TTTTTTTTTTAAAAAAAAAATTTTTTTTTT".to_vec()),
                _ => None,
            },
            Arguments::default(),
        )
        .expect("written");
        assert_eq!(written[0].start, 30);
        assert_eq!(written[1].start, 10);
    }
}
