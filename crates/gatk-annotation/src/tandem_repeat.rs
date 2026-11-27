//! `TandemRepeat` and the `GATKVariantContextUtils` repeat arithmetic it rests on (GATK 4.7.0.0;
//! `TandemRepeat.java` and `GATKVariantContextUtils.java` are unchanged since 4.6.2.0).
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
//!
//! # Reading this file without reading Rust
//!
//! - `pub use` re-exports another crate's functions under this module's name, so a caller can
//!   write `tandem_repeat::find_repeated_substring` without knowing they live in the engine.
//! - `Option<T>` is "a `T` or nothing", which the compiler forces the caller to check. It stands in
//!   for Java's `null` return (`getNumTandemRepeatUnits` returns `null` for "not a repeat") and for
//!   an exception no caller here can reach.
//! - `let Some(x) = expr else { ... };` unwraps an `Option` or runs the block, which must leave the
//!   function. Java would write `if (x == null) return ...;` after the assignment.
//! - `&bases[i..]` is a borrowed view of the array from index `i` to the end, with no copy. Java's
//!   `Arrays.copyOfRange(bases, i, bases.length)` copies; the bytes seen are the same.

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
///
/// A Rust `enum` with no data on its variants is a Java `enum`. `#[derive(...)]` asks the compiler
/// to write `equals`, copying and debug printing for it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

/// `VariantContext.getType()` (htsjdk 5.0.0), which is `getType(false)`: `determineType` and
/// `determinePolymorphicType` with `ignoreNonRef` false.
///
/// **What**: the kind of variant, by a pairwise comparison of each alternate against the reference
/// (`typeOfBiallelicVariant`), and `MIXED` as soon as two alternates disagree.
///
/// **How the reference is skipped**: htsjdk tests object identity,
///
/// ```java
/// if (allele == REF) continue;
/// ```
///
/// and `VariantContext.makeAlleles` always moves the reference to index 0
/// (`alleleList.add(0, a)`), so skipping the first allele skips the same one.
///
/// **What is not reproduced**: `typeOfBiallelicVariant` throws `IllegalStateException` for a
/// symbolic reference allele. There is no such check here; a symbolic reference is classified as
/// if its text were bases.
pub fn variant_type(vc: &VariantContext) -> Option<VariantType> {
    // `match` on the allele count, with `_` as the `default:` arm. The compiler refuses a `match`
    // that leaves a value uncovered, which is why the catch-all is required.
    match vc.alleles.len() {
        // `IllegalStateException`: a variant context with no alleles.
        0 => None,
        // Monomorphic independently of whether the one allele is the reference.
        1 => Some(VariantType::NoVariation),
        _ => {
            let reference = vc.reference();
            // Java's `Type type = null;`, typed so that the "not yet seen" case must be handled.
            let mut kind: Option<VariantType> = None;
            for allele in vc.alleles.iter().skip(1) {
                // `if` is an expression here: each branch's last value is what `biallelic` becomes.
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
                // The three arms are the reference's `if (type == null)`, `else if (biallelicType
                // != type) return Type.MIXED;` and the silent fall-through. `Some(existing) if ...`
                // is a pattern with a guard: it matches only when the condition also holds.
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

/// `VariantContext.isIndel()`: `getType() == Type.INDEL`. A context with no alleles, where the
/// reference throws, answers false here.
pub fn is_indel(vc: &VariantContext) -> bool {
    variant_type(vc) == Some(VariantType::Indel)
}

/// `TandemRepeat`: `STR`, `RU` and `RPA`.
///
/// A struct with no fields: the annotation has no state, and the type exists only to carry the
/// `InfoFieldAnnotation` implementation below.
pub struct TandemRepeat;

impl TandemRepeat {
    /// `TandemRepeat.annotate` with `getNumTandemRepeatUnits(ref, vc)`, over the reference
    /// window's bases as the caller already has them.
    ///
    /// ```java
    /// final int startIndex = vc.getStart() + 1 - ref.getWindow().getStart();
    /// ```
    ///
    /// The Java comment on that line reads `+1 to exclude leading match base common to VC's ref
    /// and alt alleles`: the context starts one base after the variant's start.
    ///
    /// The engine's `num_tandem_repeat_units` answers `None` where
    /// `GATKVariantContextUtils.getNumTandemRepeatUnits` returns `null`: every alternate the same
    /// length as the reference, or any allele whose repetition count is zero. Either way the whole
    /// annotation is dropped, as the reference's `if (result == null)` branch does.
    pub fn local_annotate(
        window_start: i64,
        window_bases: &[u8],
        vc: &VariantContext,
    ) -> Vec<(String, AnnotationValue)> {
        if !is_indel(vc) {
            return Vec::new();
        }
        let start_index = vc.start + 1 - window_start;
        // `as usize` converts a signed 64-bit integer to the unsigned index type. A negative value
        // would wrap to a huge one, which is why the sign is tested first.
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
        // The three entries in the order of the reference's `LinkedHashMap` puts.
        vec![
            (STR_PRESENT_KEY.to_string(), AnnotationValue::Flag(true)),
            (
                REPEAT_UNIT_KEY.to_string(),
                // `new String(repeatUnit)`. The unit is made of reference bases, which are ASCII,
                // so the lossy UTF-8 decoding and Java's default charset give the same text.
                AnnotationValue::Str(String::from_utf8_lossy(&repeat_unit).into_owned()),
            ),
            (
                REPEATS_PER_ALLELE_KEY.to_string(),
                // The reference allele's count first, then one per alternate, as
                // `getNumTandemRepeatUnits` built the list. `.into_iter().map(..).collect()` is a
                // Java stream that wraps each `int` as an `AnnotationValue::Int`.
                AnnotationValue::List(lengths.into_iter().map(AnnotationValue::Int).collect()),
            ),
        ]
    }
}

// `impl Trait for Type` is Java's `implements`: the methods below are `TandemRepeat`'s
// `InfoFieldAnnotation` overrides.
impl InfoFieldAnnotation for TandemRepeat {
    /// `getKeyNames()`, in the reference's order.
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
