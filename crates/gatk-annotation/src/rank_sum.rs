//! `RankSumTest` and its four members, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.6.2.0).
//!
//! `BaseQRankSum`, `MQRankSum`, `ReadPosRankSum` and `ClippingRankSum`: a Mann-Whitney U test of
//! the alternate reads against the reference reads, reported as a Z score.
//!
//! # The alternate goes first, and that is the whole sign convention
//!
//! ```java
//! // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
//! mannWhitneyU.test(Doubles.toArray(altQuals), Doubles.toArray(refQuals), FIRST_DOMINATES);
//! ```
//!
//! Swapping the two arrays flips the sign of every value this family reports, and nothing
//! downstream would notice: a negative `MQRankSum` means the alternate reads have *lower* mapping
//! quality, and that reading depends entirely on the argument order here.
//!
//! # The value is a **string**, formatted to three decimals
//!
//! ```java
//! return Collections.singletonMap(getKeyNames().get(0), String.format("%.3f", zScore));
//! ```
//!
//! Not a Double: the annotation rounds to three decimals through `String.format` and puts the
//! text in the map. The rounding is `HALF_UP` on the decimal expansion, which is not what a
//! `Double`'s own rendering would do, and the VCF carries whatever this produced.
//!
//! # A NaN Z score is not written at all
//!
//! `MannWhitneyU` answers NaN when either series is empty, and this checks for it: the key is then
//! absent from the record rather than present with a placeholder. So a site with no alternate
//! reads at all has no `MQRankSum` field, which is different from having one that reads zero.
//!
//! # Three filters, and one of them is the annotation's own
//!
//! `isInformative` (the likelihood confidence), `isUsableRead` (mapping quality neither 0 nor 255),
//! and then a value that may be absent. `ReadPosRankSum` overrides `isUsableRead` to add a soft
//! clip check against `vc.getEnd() + 1`, so the four members do not see the same reads.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::mann_whitney::{self, TestType};
use gatk_engine::read::mapping_quality;
use gatk_engine::read_utils::{self, BaseAt};
use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `GATKVCFConstants.BASE_QUAL_RANK_SUM_KEY`.
pub const BASE_QUAL_RANK_SUM_KEY: &str = "BaseQRankSum";
/// `GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY`.
pub const MAP_QUAL_RANK_SUM_KEY: &str = "MQRankSum";
/// `GATKVCFConstants.READ_POS_RANK_SUM_KEY`.
pub const READ_POS_RANK_SUM_KEY: &str = "ReadPosRankSum";
/// `GATKVCFConstants.CLIPPING_RANK_SUM_KEY`.
pub const CLIPPING_RANK_SUM_KEY: &str = "ClippingRankSum";

/// `QualityUtils.MAPPING_QUALITY_UNAVAILABLE`.
const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;

/// `RankSumTest.INVALID_ELEMENT_FROM_READ`, which is a value and not an absence.
const INVALID_ELEMENT_FROM_READ: f64 = f64::NEG_INFINITY;

/// What one member contributes.
pub trait RankSumTest {
    fn vcf_key(&self) -> &'static str;

    /// `getElementForRead`. `None` is `OptionalDouble.empty()`.
    fn element_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<f64>;

    /// `isUsableRead`, which `ReadPosRankSumTest` overrides.
    fn is_usable_read(&self, read: &BamRecord, _vc: &VariantContext) -> bool {
        let quality = mapping_quality(read);
        quality != 0 && quality != MAPPING_QUALITY_UNAVAILABLE
    }
}

/// `RankSumTest.annotate`, shared by the four.
pub fn annotate<A: RankSumTest>(
    annotation: &A,
    _reference: Option<&ReferenceContext>,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Vec<(String, AnnotationValue)> {
    // `vc.getGenotypes()` empty is an early return, before any read is looked at.
    if vc.genotypes.is_empty() {
        return Vec::new();
    }

    let mut ref_quals: Vec<f64> = Vec::new();
    let mut alt_quals: Vec<f64> = Vec::new();

    if let Some(likelihoods) = likelihoods {
        for best in likelihoods.best_alleles_breaking_ties(None) {
            if !best.is_informative() {
                continue;
            }
            let Some(read) = likelihoods
                .sample_evidence(likelihoods.index_of_sample(&best.sample).unwrap_or(0))
                .and_then(|reads| reads.get(best.evidence_index))
            else {
                continue;
            };
            if !annotation.is_usable_read(read, vc) {
                continue;
            }
            let Some(value) = annotation.element_for_read(read, vc) else {
                continue;
            };
            // A read whose clipping goal was not reached, or whose position is inside a spanning
            // deletion, comes back as -infinity rather than as an absence, and is dropped here.
            if value == INVALID_ELEMENT_FROM_READ {
                continue;
            }
            let Some(allele) = best.allele else { continue };
            if allele.is_reference() {
                ref_quals.push(value);
            } else if vc.alleles.contains(&allele) {
                alt_quals.push(value);
            }
        }
    }

    if ref_quals.is_empty() && alt_quals.is_empty() {
        return Vec::new();
    }

    // The alternate is the first series. Swapping these two flips the sign of every reported
    // value in this family.
    let result = mann_whitney::test(&alt_quals, &ref_quals, TestType::FirstDominates);
    if result.z.is_nan() {
        // Absent from the record, which is not the same as present and zero.
        return Vec::new();
    }
    vec![(
        annotation.vcf_key().to_string(),
        // `String.format("%.3f", z)`: a String in the map, not a Double.
        AnnotationValue::Str(format_three_decimals(result.z)),
    )]
}

/// `String.format("%.3f", value)` for the values a Z score can take.
///
/// Java's `%f` rounds HALF_UP on the *decimal* expansion of the double, which is not Rust's
/// `{:.3}` (round-half-to-even on the same expansion). They differ on a value whose fourth decimal
/// is exactly 5 and whose expansion terminates there, which a Z score can be: `0.0625` prints
/// `0.063` in Java and `0.062` in Rust.
fn format_three_decimals(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    // The exact decimal expansion of the double, then a half-up rounding of it.
    let text = format!("{:.*}", 30, value.abs());
    let (whole, fraction) = text.split_once('.').expect("a decimal point");
    let mut digits: Vec<u8> = whole
        .bytes()
        .chain(fraction.bytes().take(3))
        .map(|b| b - b'0')
        .collect();
    let round_up = fraction.as_bytes()[3] >= b'5';
    if round_up {
        let mut index = digits.len();
        loop {
            if index == 0 {
                digits.insert(0, 1);
                break;
            }
            index -= 1;
            if digits[index] == 9 {
                digits[index] = 0;
            } else {
                digits[index] += 1;
                break;
            }
        }
    }
    let split = digits.len() - 3;
    let whole: String = digits[..split].iter().map(|d| (d + b'0') as char).collect();
    let fraction: String = digits[split..].iter().map(|d| (d + b'0') as char).collect();
    let sign = if value.is_sign_negative() { "-" } else { "" };
    format!("{sign}{whole}.{fraction}")
}

/// `BaseQualityRankSumTest`: `BaseQRankSum`.
pub struct BaseQualityRankSumTest;

impl RankSumTest for BaseQualityRankSumTest {
    fn vcf_key(&self) -> &'static str {
        BASE_QUAL_RANK_SUM_KEY
    }

    fn element_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<f64> {
        match read_utils::read_base_quality_at_reference_coordinate(read, vc.start as i32) {
            BaseAt::Present(quality) => Some(quality as f64),
            _ => None,
        }
    }
}

/// `MappingQualityRankSumTest`: `MQRankSum`.
pub struct MappingQualityRankSumTest;

impl RankSumTest for MappingQualityRankSumTest {
    fn vcf_key(&self) -> &'static str {
        MAP_QUAL_RANK_SUM_KEY
    }

    fn element_for_read(&self, read: &BamRecord, _vc: &VariantContext) -> Option<f64> {
        Some(mapping_quality(read) as f64)
    }
}

/// `ReadPosRankSumTest`: `ReadPosRankSum`, the one that overrides `isUsableRead`.
pub struct ReadPosRankSumTest;

impl RankSumTest for ReadPosRankSumTest {
    fn vcf_key(&self) -> &'static str {
        READ_POS_RANK_SUM_KEY
    }

    fn element_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<f64> {
        // The same computation `ReadPosition` uses, minus that annotation's own overlap guard:
        // this one guards through `is_usable_read` instead, and the two guards are not the same.
        read_position(read, vc)
    }

    fn is_usable_read(&self, read: &BamRecord, vc: &VariantContext) -> bool {
        let quality = mapping_quality(read);
        // `vc.getEnd() + 1` rather than `vc.getEnd()`, "in case of a leading indel".
        quality != 0
            && quality != MAPPING_QUALITY_UNAVAILABLE
            && read_utils::soft_start(read) <= vc.stop as i32 + 1
            && read_utils::soft_end(read) >= vc.start as i32
    }
}

/// `ReadPosRankSumTest.getReadPosition`, which `ReadPosition` also calls.
fn read_position(read: &BamRecord, vc: &VariantContext) -> Option<f64> {
    if read_utils::start(read) == vc.stop as i32 + 1 && opens_with_insertion(read) {
        return Some(0.0);
    }
    let (index, _) = read_utils::read_index_for_reference_coordinate(
        read_utils::start(read),
        &read.cigar,
        vc.start as i32,
    );
    if index < 0 {
        return None;
    }
    let (leading_hard, trailing_hard) = hard_clips(read);
    let left = leading_hard + index;
    let right = read.read_bases.len() as i32 - 1 - index + trailing_hard;
    Some(left.min(right) as f64)
}

fn opens_with_insertion(read: &BamRecord) -> bool {
    read.cigar
        .elements
        .iter()
        .find(|element| !matches!(element.op, Op::S | Op::H))
        .map(|element| element.op == Op::I)
        .unwrap_or(false)
}

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

/// `ClippingRankSumTest`: `ClippingRankSum`, over `AlignmentUtils.getNumHardClippedBases`.
pub struct ClippingRankSumTest;

impl RankSumTest for ClippingRankSumTest {
    fn vcf_key(&self) -> &'static str {
        CLIPPING_RANK_SUM_KEY
    }

    fn element_for_read(&self, read: &BamRecord, _vc: &VariantContext) -> Option<f64> {
        // `getNumHardClippedBases`: every H element, wherever it sits, not only the two ends.
        let clipped: u32 = read
            .cigar
            .elements
            .iter()
            .filter(|element| element.op == Op::H)
            .map(|element| element.length)
            .sum();
        Some(clipped as f64)
    }
}

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

info_annotation_for!(BaseQualityRankSumTest, BASE_QUAL_RANK_SUM_KEY);
info_annotation_for!(MappingQualityRankSumTest, MAP_QUAL_RANK_SUM_KEY);
info_annotation_for!(ReadPosRankSumTest, READ_POS_RANK_SUM_KEY);
info_annotation_for!(ClippingRankSumTest, CLIPPING_RANK_SUM_KEY);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_formatter_rounds_half_up_where_rust_would_round_to_even() {
        assert_eq!(format_three_decimals(0.0625), "0.063");
        assert_eq!(format_three_decimals(-0.0625), "-0.063");
        assert_eq!(format_three_decimals(1.2344), "1.234");
        assert_eq!(format_three_decimals(1.2345), "1.234");
        assert_eq!(format_three_decimals(0.0), "0.000");
        assert_eq!(format_three_decimals(-1.0), "-1.000");
        assert_eq!(format_three_decimals(9.9999), "10.000");
    }
}
