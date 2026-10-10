//! `TandemRepeat` and the `GATKVariantContextUtils` repeat arithmetic it rests on (GATK 4.6.2.0).
//!
//! `STR`, `RU` and `RPA`: is this indel an expansion or contraction of a repeat, what is the repeat
//! unit, and how many copies does each allele carry?
//!
//! # The repeat arithmetic is the engine's
//!
//! `findRepeatedSubstring`, `findNumberOfRepetitions` and `getNumTandemRepeatUnits` are in
//! `gatk_engine::tandem_repeat_units`, with their notes; they are re-exported here.
//!
//! # `STR` is a boolean `true` in the map, so the key is written bare
//!
//! `map.put(GATKVCFConstants.STR_PRESENT_KEY, true)` puts a `Boolean`, and the VCF encoder writes a
//! flag as its key alone with no `=value`. `RU` is a `String` and `RPA` an `ArrayList<Integer>`.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
pub use gatk_engine::tandem_repeat_units::{
    find_number_of_repetitions, find_repeated_substring, num_tandem_repeat_units,
    num_tandem_repeat_units_for_bases,
};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `GATKVCFConstants.STR_PRESENT_KEY`.
pub const STR_PRESENT_KEY: &str = "STR";
/// `GATKVCFConstants.REPEAT_UNIT_KEY`.
pub const REPEAT_UNIT_KEY: &str = "RU";
/// `GATKVCFConstants.REPEATS_PER_ALLELE_KEY`.
pub const REPEATS_PER_ALLELE_KEY: &str = "RPA";

/// `VariantContext.Type`, as far as this annotation needs it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

/// `VariantContext.getType()`, which is a pairwise comparison of each alternate against the
/// reference and `MIXED` as soon as two alternates disagree.
pub fn variant_type(vc: &VariantContext) -> Option<VariantType> {
    match vc.alleles.len() {
        // `IllegalStateException`: a variant context with no alleles.
        0 => None,
        // Monomorphic independently of whether the one allele is the reference.
        1 => Some(VariantType::NoVariation),
        _ => {
            let reference = vc.reference();
            let mut kind: Option<VariantType> = None;
            for allele in vc.alleles.iter().skip(1) {
                let biallelic = if allele.is_symbolic() {
                    VariantType::Symbolic
                } else if reference.len() == allele.len() {
                    if allele.len() == 1 {
                        VariantType::Snp
                    } else {
                        VariantType::Mnp
                    }
                } else {
                    // Not a SNP, an MNP or symbolic, so necessarily an indel: the prefix check the
                    // reference used to make here was wrong for `CTTA -> C,CT,CA`.
                    VariantType::Indel
                };
                match kind {
                    None => kind = Some(biallelic),
                    Some(existing) if existing != biallelic => return Some(VariantType::Mixed),
                    Some(_) => {}
                }
            }
            kind
        }
    }
}

/// `VariantContext.isIndel()`.
pub fn is_indel(vc: &VariantContext) -> bool {
    variant_type(vc) == Some(VariantType::Indel)
}

/// `TandemRepeat`: `STR`, `RU` and `RPA`.
pub struct TandemRepeat;

impl TandemRepeat {
    /// `TandemRepeat.annotate`, over the reference window's bases as the caller already has them.
    ///
    /// The `+ 1` excludes the padding base the variant's reference and alternate alleles share, so
    /// the context starts one base after the variant's start.
    pub fn local_annotate(
        window_start: i64,
        window_bases: &[u8],
        vc: &VariantContext,
    ) -> Vec<(String, AnnotationValue)> {
        if !is_indel(vc) {
            return Vec::new();
        }
        let start_index = vc.start + 1 - window_start;
        if start_index < 0 || start_index as usize > window_bases.len() {
            // `Arrays.copyOfRange` with a negative or over-long start is an
            // `ArrayIndexOutOfBoundsException`, which no walker-built window can produce.
            return Vec::new();
        }
        let context = &window_bases[start_index as usize..];
        let Some((lengths, repeat_unit)) =
            num_tandem_repeat_units(vc.reference(), vc.alternate_alleles(), context)
        else {
            return Vec::new();
        };
        vec![
            (STR_PRESENT_KEY.to_string(), AnnotationValue::Flag(true)),
            (
                REPEAT_UNIT_KEY.to_string(),
                AnnotationValue::Str(String::from_utf8_lossy(&repeat_unit).into_owned()),
            ),
            (
                REPEATS_PER_ALLELE_KEY.to_string(),
                AnnotationValue::List(lengths.into_iter().map(AnnotationValue::Int).collect()),
            ),
        ]
    }
}

impl InfoFieldAnnotation for TandemRepeat {
    fn key_names(&self) -> Vec<&'static str> {
        vec![STR_PRESENT_KEY, REPEAT_UNIT_KEY, REPEATS_PER_ALLELE_KEY]
    }

    /// Without a reference window there is nothing to count against, and the reference
    /// dereferences the context straight away. Use [`TandemRepeat::local_annotate`].
    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        _vc: &VariantContext,
        _likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_partial_trailing_repeat_is_not_a_repeat() {
        assert_eq!(find_repeated_substring(b"ACTACT"), 3);
        // Eight bases, two and two thirds copies of ACT: reported as a unit of itself.
        assert_eq!(find_repeated_substring(b"ACTACTAC"), 8);
    }

    #[test]
    fn an_empty_string_gives_a_unit_of_one_zero_byte() {
        assert_eq!(find_repeated_substring(b""), 1);
    }

    #[test]
    fn leading_and_trailing_repeats_differ() {
        assert_eq!(find_number_of_repetitions(b"AT", b"GATAT", true), Some(0));
        assert_eq!(find_number_of_repetitions(b"AT", b"GATAT", false), Some(2));
        assert_eq!(
            find_number_of_repetitions(b"CCC", b"CCCCCCCC", true),
            Some(2)
        );
    }

    #[test]
    fn a_deletion_of_one_unit_is_counted_against_the_context() {
        // Reference GATCCACCACCAGTCGA, variant TCCA -> T, so the context that follows the padding
        // base is CCACCACCAGTCGA and still contains the deleted unit.
        let (counts, unit) =
            num_tandem_repeat_units_for_bases(b"CCA", b"", b"CCACCACCAGTCGA").expect("a repeat");
        assert_eq!(unit, b"CCA".to_vec());
        assert_eq!(counts, [3, 2]);
    }
}
