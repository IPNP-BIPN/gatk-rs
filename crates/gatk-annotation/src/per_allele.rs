//! `PerAlleleAnnotation` and its four members, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.6.2.0).
//!
//! `MBQ`, `MMQ`, `MPOS` and `MFRL`: the median, per allele, of one number taken from each read
//! that best supports it. One shared shape, four values, and four different answers to the same
//! question about an allele no read supports.
//!
//! # The empty case is a different number in each of the four
//!
//! ```java
//! BaseQuality       values.isEmpty() ? 0  : MathUtils.median(...)
//! MappingQuality    values.isEmpty() ? 60 : MathUtils.median(...)   // VALUE_FOR_NO_READS
//! ReadPosition      values.isEmpty() ? 50 : MathUtils.median(...)   // VALUE_FOR_NO_READS
//! FragmentLength    values.isEmpty() ? 0  : MathUtils.median(...)
//! ```
//!
//! The two that invent a value say why in a comment: *"we don't want a GGA mode allele with no
//! reads to prejudice us against a site"*. So an allele with no evidence is annotated as a good
//! mapping quality and a comfortable read position, and as a zero base quality. A port that
//! factored the empty case out would agree on every ordinary record and disagree on exactly the
//! sites the values were invented for.
//!
//! # Three filters run before a read contributes, and they are not the same filter
//!
//! ```java
//! .filter(ba -> ba.isInformative() && isUsableRead(ba.evidence))
//! .forEach(ba -> getValueForRead(ba.evidence, vc).ifPresent(v -> values.get(ba.allele).add(v)));
//! ```
//!
//! `isInformative` is the likelihood confidence against the log10 threshold; `isUsableRead` is
//! mapping quality neither `0` nor `255`; and `getValueForRead` can decline for its own reason.
//! `MappingQuality` therefore never sees a read of quality 0, so its median can never be 0 by way
//! of a read, only by way of the empty case, which is 60.
//!
//! # The key carries an `int[]`, and the alleles it covers are not the same set
//!
//! `ImmutableMap.of(getVcfKey(), statistics)` boxes a Java `int[]`, so the value has that class
//! and prints through `Arrays.toString` rather than as a list. And `includeRefAllele()` is
//! `false` by default and overridden to `true` in three of the four, so `MPOS` reports one number
//! per **alternate** allele while the other three report one per allele including the reference.
//!
//! # The median is commons-math3, not the obvious one
//!
//! `MathUtils.median` goes through commons-math3 `Percentile` and finishes with
//! `FastMath.round`, which is `(long) floor(x + 0.5)` and not `Math.round`. Both are in
//! [`jmath`]; this reaches the first. See htsjdk-rs decision 0023.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::read::mapping_quality;
use gatk_engine::read_utils::{self, BaseAt};
use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;
use jmath::percentile::{median_of_ints, EstimationType};

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `GATKVCFConstants.MEDIAN_BASE_QUALITY_KEY`.
pub const MEDIAN_BASE_QUALITY_KEY: &str = "MBQ";
/// `GATKVCFConstants.MEDIAN_MAPPING_QUALITY_KEY`.
pub const MEDIAN_MAPPING_QUALITY_KEY: &str = "MMQ";
/// `GATKVCFConstants.MEDIAN_READ_POSITON_KEY`, spelled as the reference spells it.
pub const MEDIAN_READ_POSITION_KEY: &str = "MPOS";
/// `GATKVCFConstants.MEDIAN_FRAGMENT_LENGTH_KEY`.
pub const MEDIAN_FRAGMENT_LENGTH_KEY: &str = "MFRL";

/// `QualityUtils.MAPPING_QUALITY_UNAVAILABLE`.
const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;

/// What one member of the family contributes: a value per read, a value for no reads, and whether
/// the reference allele is reported on.
pub trait PerAlleleAnnotation {
    fn vcf_key(&self) -> &'static str;

    /// `includeRefAllele()`, `false` in the parent and overridden in three of the four.
    fn include_ref_allele(&self) -> bool {
        false
    }

    /// `getValueForRead`. `None` is the reference's `OptionalInt.empty()`, which contributes
    /// nothing without stopping the traversal.
    fn value_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<i32>;

    /// `aggregate(values)`: the median, or the annotation's own value for no reads.
    fn value_for_no_reads(&self) -> i32;

    fn aggregate(&self, values: &[i32]) -> i32 {
        if values.is_empty() {
            self.value_for_no_reads()
        } else {
            median_of_ints(values, EstimationType::Legacy)
        }
    }
}

/// `PerAlleleAnnotation.isUsableRead`: a mapping quality that is neither absent nor zero.
fn is_usable_read(read: &BamRecord) -> bool {
    let quality = mapping_quality(read);
    quality != 0 && quality != MAPPING_QUALITY_UNAVAILABLE
}

/// `PerAlleleAnnotation.annotate`, shared by all four.
///
/// The traversal is over `bestAllelesBreakingTies`, so a read contributes to exactly one allele:
/// the one it best supports, ties broken by the reference's own rule.
pub fn annotate<A: PerAlleleAnnotation>(
    annotation: &A,
    _reference: Option<&ReferenceContext>,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Vec<(String, AnnotationValue)> {
    let Some(likelihoods) = likelihoods else {
        return Vec::new();
    };

    // `Collectors.toMap` over `likelihoods.alleles()`, so the buckets are the matrix's alleles and
    // not the variant's. An allele the matrix does not hold gets no bucket, which is why the
    // lookup below can be absent and is not an error.
    let mut values: Vec<(Allele, Vec<i32>)> = (0..likelihoods.number_of_alleles())
        .filter_map(|index| likelihoods.get_allele(index).cloned())
        .map(|allele| (allele, Vec::new()))
        .collect();

    for best in likelihoods.best_alleles_breaking_ties(None) {
        if !best.is_informative() {
            continue;
        }
        let Some(evidence) = likelihoods
            .sample_evidence(likelihoods.index_of_sample(&best.sample).unwrap_or(0))
            .and_then(|reads| reads.get(best.evidence_index))
        else {
            continue;
        };
        if !is_usable_read(evidence) {
            continue;
        }
        let Some(value) = annotation.value_for_read(evidence, vc) else {
            continue;
        };
        let Some(allele) = best.allele else {
            continue;
        };
        if let Some(bucket) = values.iter_mut().find(|(a, _)| *a == allele) {
            bucket.1.push(value);
        }
    }

    // `vc.getAlleles()`, filtered, in the variant's order: the matrix's order does not decide
    // this, and a variant allele the matrix never held aggregates an empty list.
    let statistics: Vec<AnnotationValue> = vc
        .alleles
        .iter()
        .filter(|allele| !allele.is_reference() || annotation.include_ref_allele())
        .map(|allele| {
            let collected = values
                .iter()
                .find(|(a, _)| a == allele)
                .map(|(_, v)| v.as_slice())
                .unwrap_or(&[]);
            AnnotationValue::Int(annotation.aggregate(collected))
        })
        .collect();

    vec![(
        annotation.vcf_key().to_string(),
        // `ImmutableMap.of(key, int[])`: an int array, whatever its length, and never a scalar.
        AnnotationValue::List(statistics),
    )]
}

/// `BaseQuality`: `MBQ`, the median base quality at the variant's start.
pub struct BaseQuality;

impl BaseQuality {
    /// `BaseQualityRankSumTest.getReadBaseQuality`, then `FastMath.round` and a narrowing cast.
    ///
    /// The guard is the reference's own: a read that starts after the variant, or ends before it,
    /// declines before the base is looked up at all.
    pub fn base_quality(read: &BamRecord, vc: &VariantContext) -> Option<i32> {
        let start = vc.start as i32;
        if start < read_utils::start(read) || read_utils::end(read) < start {
            return None;
        }
        match read_utils::read_base_quality_at_reference_coordinate(read, start) {
            BaseAt::Present(quality) => {
                // The reference boxes the byte into an OptionalDouble and rounds it back, which
                // is exact for every quality a BAM can hold.
                Some(jmath::fast_math::round(quality as f64) as i32)
            }
            _ => None,
        }
    }
}

impl PerAlleleAnnotation for BaseQuality {
    fn vcf_key(&self) -> &'static str {
        MEDIAN_BASE_QUALITY_KEY
    }

    fn include_ref_allele(&self) -> bool {
        true
    }

    fn value_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<i32> {
        Self::base_quality(read, vc)
    }

    fn value_for_no_reads(&self) -> i32 {
        0
    }
}

/// `MappingQuality`: `MMQ`, whose value for no reads is 60.
pub struct MappingQuality;

impl PerAlleleAnnotation for MappingQuality {
    fn vcf_key(&self) -> &'static str {
        MEDIAN_MAPPING_QUALITY_KEY
    }

    fn include_ref_allele(&self) -> bool {
        true
    }

    fn value_for_read(&self, read: &BamRecord, _vc: &VariantContext) -> Option<i32> {
        // Never declines: every read has a mapping quality, and the ones this annotation would
        // rather not see were already dropped by `isUsableRead`.
        Some(mapping_quality(read) as i32)
    }

    fn value_for_no_reads(&self) -> i32 {
        60
    }
}

/// `ReadPosition`: `MPOS`, the distance to the nearer end of the read.
///
/// The only member of the four that does **not** report on the reference allele.
pub struct ReadPosition;

impl ReadPosition {
    /// `ReadPosRankSumTest.getReadPosition`, guarded as `ReadPosition.getPosition` guards it.
    pub fn position(read: &BamRecord, vc: &VariantContext) -> Option<i32> {
        let start = vc.start as i32;
        if start < read_utils::start(read) || read_utils::end(read) < start {
            return None;
        }
        // A read that opens with an insertion sits one base past the variant's end, which looks
        // like no overlap and is answered with 0 rather than with nothing.
        if read_utils::start(read) == vc.stop as i32 + 1 && opens_with_insertion(read) {
            return Some(0);
        }
        let (index, _) = read_utils::read_index_for_reference_coordinate(
            read_utils::start(read),
            &read.cigar,
            start,
        );
        if index < 0 {
            return None;
        }
        // Hard clips are bases that were removed to fit an assembly region, so the distance is
        // measured as if they were still there.
        let (leading_hard, trailing_hard) = hard_clips(read);
        let left = leading_hard + index;
        let right = read.read_bases.len() as i32 - 1 - index + trailing_hard;
        Some(left.min(right))
    }
}

/// Whether the first non-clipping operator is an insertion.
fn opens_with_insertion(read: &BamRecord) -> bool {
    read.cigar
        .elements
        .iter()
        .find(|element| !matches!(element.op, Op::S | Op::H))
        .map(|element| element.op == Op::I)
        .unwrap_or(false)
}

/// The leading and trailing hard clip lengths, zero where the end is any other operator.
fn hard_clips(read: &BamRecord) -> (i32, i32) {
    let leading = read
        .cigar
        .elements
        .first()
        .filter(|e| e.op == Op::H)
        .map(|e| e.length as i32)
        .unwrap_or(0);
    let trailing = read
        .cigar
        .elements
        .last()
        .filter(|e| e.op == Op::H)
        .map(|e| e.length as i32)
        .unwrap_or(0);
    (leading, trailing)
}

impl PerAlleleAnnotation for ReadPosition {
    fn vcf_key(&self) -> &'static str {
        MEDIAN_READ_POSITION_KEY
    }

    fn value_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<i32> {
        Self::position(read, vc)
    }

    fn value_for_no_reads(&self) -> i32 {
        50
    }
}

/// `FragmentLength`: `MFRL`, the absolute template length.
pub struct FragmentLength;

impl PerAlleleAnnotation for FragmentLength {
    fn vcf_key(&self) -> &'static str {
        MEDIAN_FRAGMENT_LENGTH_KEY
    }

    fn include_ref_allele(&self) -> bool {
        true
    }

    fn value_for_read(&self, read: &BamRecord, _vc: &VariantContext) -> Option<i32> {
        // `Math.abs` of Integer.MIN_VALUE is itself, and the port keeps that rather than
        // saturating: a template length of Integer.MIN_VALUE is not reachable from a BAM, and an
        // absolute value that clamped would be a different function.
        Some(read.inferred_insert_size.wrapping_abs())
    }

    fn value_for_no_reads(&self) -> i32 {
        0
    }
}

/// The four wrapped as [`InfoFieldAnnotation`]s, so a caller can hold them beside the others.
macro_rules! info_annotation_for {
    ($name:ident, $key:expr) => {
        impl InfoFieldAnnotation for $name {
            fn key_names(&self) -> Vec<&'static str> {
                vec![$key]
            }

            fn annotate(
                &self,
                reference: Option<&ReferenceContext>,
                vc: &VariantContext,
                likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
            ) -> Vec<(String, AnnotationValue)> {
                annotate(self, reference, vc, likelihoods)
            }
        }
    };
}

info_annotation_for!(BaseQuality, MEDIAN_BASE_QUALITY_KEY);
info_annotation_for!(MappingQuality, MEDIAN_MAPPING_QUALITY_KEY);
info_annotation_for!(ReadPosition, MEDIAN_READ_POSITION_KEY);
info_annotation_for!(FragmentLength, MEDIAN_FRAGMENT_LENGTH_KEY);
