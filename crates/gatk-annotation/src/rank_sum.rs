//! `RankSumTest` and four of its members, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.7.0.0; `RankSumTest.java` and
//! the four members' files are unchanged since 4.6.2.0).
//!
//! `BaseQRankSum`, `MQRankSum`, `ReadPosRankSum` and `ClippingRankSum`: a Mann-Whitney U test of
//! the alternate reads against the reference reads, reported as a Z score. A fifth member,
//! `LikelihoodRankSumTest`, implements the same trait in `site_statistics.rs`.
//!
//! # The alternate goes first, and that is the whole sign convention
//!
//! `RankSumTest.java:48-49` (the call is one line in the source, wrapped here):
//!
//! ```java
//! // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
//! final MannWhitneyU.Result result = mannWhitneyU.test(Doubles.toArray(altQuals),
//! Doubles.toArray(refQuals), MannWhitneyU.TestType.FIRST_DOMINATES);
//! ```
//!
//! So a negative `MQRankSum` means the alternate reads have *lower* mapping quality, and that
//! reading depends entirely on the argument order here.
//!
//! Swapping the two arrays would **not** simply flip the sign. With either series at ten or more
//! (the normal approximation), `MannWhitneyU.calculateZ` keeps the same continuity correction
//! whichever series is first (`MannWhitneyU.java:334` and `:344`):
//!
//! ```java
//! correction = whichSide == TestType.FIRST_DOMINATES ? -.5 : .5;
//! ...
//! return (u - m - correction) / sigma;
//! ```
//!
//! U1 + U2 is n1 * n2, so in exact arithmetic the swapped Z is `-Z + 1/sigma` whenever the
//! correction applies, and exactly `-Z` only when `nties` is zero and the correction is dropped
//! (`MannWhitneyU.java:339-341`; see `gatk_engine::mann_whitney::test`). The golden shows it:
//! `tests/data/rank_sum.txt.gz` line 31, `identical-groups`, has twelve reads on each side with
//! the same twelve mapping qualities, so swapping the series changes nothing, and it reports
//! `MQRankSum=0.029` (0.5 / sigma), not zero.
//!
//! # The value is a **string**, formatted to three decimals
//!
//! ```java
//! return Collections.singletonMap(getKeyNames().get(0), String.format("%.3f", zScore));
//! ```
//!
//! Not a Double: the annotation rounds to three decimals through `String.format` and puts the
//! text in the map. The rounding is `HALF_UP` on the digits `Double.toString` would produce, not on
//! the exact decimal expansion (see `format_three_decimals`), and the VCF carries whatever this
//! produced.
//!
//! # A NaN Z score is not written at all
//!
//! `MannWhitneyU` answers NaN when either series is empty, and this checks for it: the key is then
//! absent from the record rather than present with a placeholder. So a site with no alternate
//! reads at all has no `MQRankSum` field, which is different from having one that reads zero
//! (`rank_sum.txt.gz` lines 34-37, `ref-only`, carry no value, nor do the `alt-only` rows after
//! them).
//!
//! # Three filters, and one of them is the annotation's own
//!
//! `isInformative` (the likelihood confidence), `isUsableRead` (mapping quality neither 0 nor 255),
//! and then a value that may be absent. `ReadPosRankSum` overrides `isUsableRead` to add a soft
//! clip check against `vc.getEnd() + 1`, so the four members do not see the same reads.
//!
//! # Reading this file without reading Rust
//!
//! - `trait RankSumTest` plays the part of the Java abstract class: each member is a `struct`
//!   with no fields, and `impl RankSumTest for X` is its `extends RankSumTest` with the overrides.
//!   A method with a body in the trait is a default the member may override, as in Java.
//! - `Option<f64>` is `OptionalDouble`: `Some(x)` is present, `None` is empty.
//! - `let Some(x) = expr else { continue; };` unwraps a present value into `x` or takes the
//!   `else` branch; there is no Java equivalent, and here it is the `if (!present) continue;`.
//! - `&BamRecord` is a borrowed, read-only reference: the function reads the record and the caller
//!   keeps it.
//! - `x as f64` and `x as i32` are Java's widening and narrowing casts.
//! - `|element| ...` is a lambda, and `.iter().filter(..).map(..).sum()` is Java's stream chain.
//! - A function's last expression, written without `;`, is its return value.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::java_format;
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
pub(crate) const INVALID_ELEMENT_FROM_READ: f64 = f64::NEG_INFINITY;

/// What one member contributes: the abstract and overridable methods of `RankSumTest`.
pub trait RankSumTest {
    /// `getKeyNames().get(0)`: the one INFO key this member writes.
    fn vcf_key(&self) -> &'static str;

    /// `getElementForRead`. `None` is `OptionalDouble.empty()`.
    fn element_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<f64>;

    /// `getElementForRead(read, vc, bestAllele)`, which only `LikelihoodRankSumTest` overrides.
    ///
    /// The parent's form delegates to the two-argument one and ignores the best allele entirely.
    /// That delegation is not here: [`annotate`] calls [`RankSumTest::element_for_read`] instead
    /// whenever [`RankSumTest::uses_best_allele`] is false, so this default is never reached and
    /// answers `None`. The override in `site_statistics.rs` returns the likelihood, as
    /// `LikelihoodRankSumTest.java:41` does:
    ///
    /// ```java
    /// return OptionalDouble.of(bestAllele.likelihood);
    /// ```
    fn element_for_best_allele(&self, _likelihood: f64) -> Option<f64> {
        None
    }

    /// Whether this member reads the best allele's likelihood rather than the read.
    ///
    /// The reference needs no such switch, because the override of the three-argument
    /// `getElementForRead` is itself the switch; here it picks which of the two methods above
    /// [`annotate`] calls.
    fn uses_best_allele(&self) -> bool {
        false
    }

    /// `isUsableRead`, which `ReadPosRankSumTest` overrides.
    fn is_usable_read(&self, read: &BamRecord, _vc: &VariantContext) -> bool {
        let quality = mapping_quality(read);
        quality != 0 && quality != MAPPING_QUALITY_UNAVAILABLE
    }
}

/// `RankSumTest.annotate` with `fillQualsFromLikelihood` inlined, shared by every member.
///
/// `<A: RankSumTest>` makes this a generic method over any member, Java's
/// `<A extends RankSumTest>`. The returned `Vec` of key and value pairs is the reference's
/// `Map<String, Object>`; an empty one is `Collections.emptyMap()`.
///
/// The filters run in the reference's order, `RankSumTest.java:63-66` (indentation reduced, and
/// the comment line wrapped):
///
/// ```java
/// if (bestAllele.isInformative() && isUsableRead(read, vc)) {
///     final OptionalDouble value = getElementForRead(read, vc, bestAllele);
///     // Bypass read if the clipping goal is not reached or the refloc is inside a
///     // spanning deletion
///     if (value.isPresent() && value.getAsDouble() != INVALID_ELEMENT_FROM_READ) {
/// ```
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

    // `if let Some(x) = y { .. }` runs the block only when `y` is present:
    // `if (likelihoods != null)`.
    if let Some(likelihoods) = likelihoods {
        for best in likelihoods.best_alleles_breaking_ties(None) {
            if !best.is_informative() {
                continue;
            }
            // The reference's `bestAllele.evidence` is the read itself; here it is an index into
            // the sample's reads, looked up again. `and_then` chains a second lookup that may also
            // be absent. The sample is one of this matrix's own, so `unwrap_or(0)` does not decide
            // anything.
            let Some(read) = likelihoods
                .sample_evidence(likelihoods.index_of_sample(&best.sample).unwrap_or(0))
                .and_then(|reads| reads.get(best.evidence_index))
            else {
                continue;
            };
            if !annotation.is_usable_read(read, vc) {
                continue;
            }
            let Some(value) = (if annotation.uses_best_allele() {
                annotation.element_for_best_allele(best.likelihood)
            } else {
                annotation.element_for_read(read, vc)
            }) else {
                continue;
            };
            // The reference's comment says -infinity marks a read whose clipping goal was not
            // reached, or whose position is inside a spanning deletion. No member in GATK 4.7.0.0
            // returns it (the constant appears only here and in `AS_RankSumTest`), so this check
            // is the reference's, kept for parity.
            if value == INVALID_ELEMENT_FROM_READ {
                continue;
            }
            let Some(allele) = best.allele else { continue };
            // `allele.isReference()`, else `vc.hasAllele(allele)`: a best allele the variant does
            // not carry joins neither series.
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

    // The alternate is the first series, as in the reference. Swapping the two negates Z exactly
    // only when `nties` is zero: whenever the continuity correction applies, the swapped Z is
    // `-Z + 1/sigma` (`MannWhitneyU.java:339-341` and `:344`; `rank_sum.txt.gz` line 31).
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
/// Java's `%f` rounds HALF_UP on the digits `Double.toString` would produce, which is neither
/// Rust's `{:.3}` (half-to-even on the *exact expansion*) nor a half-up rounding of that
/// expansion. Two values a Z score can take show both differences: `0.0625` prints `0.063` in
/// Java and `0.062` in Rust, and `1.2345` prints `1.235` in Java where a half-up rounding of the
/// expansion, which is `1.23449999999999993...`, prints `1.234`. The `string-format` golden
/// (`crates/gatk-engine/tests/data/string_format.txt.gz`) pins both.
pub(crate) fn format_three_decimals(value: f64) -> String {
    format_decimals(value, 3)
}

/// `String.format("%.Nf", value)`, half-up on the digits `Double.toString` would produce.
///
/// One implementation, in [`gatk_engine::java_format`], because `ClipReads` prints percentages
/// through the same rule and two copies of a rounding rule is how a rounding rule drifts.
pub(crate) fn format_decimals(value: f64, places: usize) -> String {
    java_format::format_decimals(value, places)
}

/// `BaseQualityRankSumTest`: `BaseQRankSum`.
pub struct BaseQualityRankSumTest;

impl RankSumTest for BaseQualityRankSumTest {
    fn vcf_key(&self) -> &'static str {
        BASE_QUAL_RANK_SUM_KEY
    }

    /// `BaseQualityRankSumTest.getReadBaseQuality`, `BaseQualityRankSumTest.java:42-43` (each
    /// line is one line in the source, wrapped here):
    ///
    /// ```java
    /// final Optional<Byte> readBaseQuality =
    ///         ReadUtils.getReadBaseQualityAtReferenceCoordinate(read, vc.getStart());
    /// return readBaseQuality.isPresent() ? OptionalDouble.of(readBaseQuality.get())
    ///         : OptionalDouble.empty();
    /// ```
    ///
    /// No overlap guard, unlike `BaseQuality.getBaseQuality` in `per_allele.rs`. `match` is a
    /// `switch` the compiler checks for completeness, and `_` is its `default`: both the absent
    /// base and the index the reference would throw on come back as empty here.
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

    /// `ReadPosRankSumTest.isUsableRead`, `ReadPosRankSumTest.java:49-50` (each line is one line
    /// in the source, wrapped here):
    ///
    /// ```java
    /// // we use vc.getEnd() + 1 in case of a leading indel -- if this isn't relevant
    /// // getReadPosition will return empty
    /// return super.isUsableRead(read, vc) && read.getSoftStart() <= vc.getEnd() + 1
    ///         && read.getSoftEnd() >= vc.getStart();
    /// ```
    fn is_usable_read(&self, read: &BamRecord, vc: &VariantContext) -> bool {
        let quality = mapping_quality(read);
        // `vc.getEnd() + 1` rather than `vc.getEnd()`, "in case of a leading indel".
        quality != 0
            && quality != MAPPING_QUALITY_UNAVAILABLE
            && read_utils::soft_start(read) <= vc.stop as i32 + 1
            && read_utils::soft_end(read) >= vc.start as i32
    }
}

/// `ReadPosRankSumTest.getReadPosition`, which `ReadPosition` also calls: the distance from the
/// variant's start to the nearer end of the read, counted in read bases.
///
/// Hard clips count as bases. `ReadPosRankSumTest.java:69` (one line in the source, wrapped here):
///
/// ```java
/// // hard clips at this point in the code are perfectly good bases that were clipped to make
/// // the read fit the assembly region
/// ```
///
/// A read that opens with an insertion right after the variant answers 0 (`:56-60`, "which
/// appears like no overlap"). `ReadPosRankSum` can reach that branch because its `isUsableRead`
/// allows `vc.getEnd() + 1`; `rank_sum.txt.gz` line 56, `leading-insertion`, reports
/// `ReadPosRankSum=-0.674`.
fn read_position(read: &BamRecord, vc: &VariantContext) -> Option<f64> {
    if read_utils::start(read) == vc.stop as i32 + 1 && opens_with_insertion(read) {
        return Some(0.0);
    }
    // `getReadIndexForReferenceCoordinate(read, ...)` walks from the SOFT start.
    // `let (index, _) = ..` takes the first of the returned pair (the reference's `getLeft()`)
    // and discards the operator. A negative index is `READ_INDEX_NOT_FOUND`.
    let (index, _) = read_utils::read_index_for_read(read, vc.start as i32);
    if index < 0 {
        return None;
    }
    let (leading_hard, trailing_hard) = hard_clips(read);
    let left = leading_hard + index;
    // `read.getLength()`: the bases stored in the record, which excludes hard clips.
    let right = read.read_bases.len() as i32 - 1 - index + trailing_hard;
    Some(left.min(right) as f64)
}

/// Whether the first CIGAR operator that is not a clip is an insertion: the reference's
/// `.filter(o -> !o.isClipping()).findFirst().orElse(null) == CigarOperator.INSERTION`.
///
/// `matches!(x, Op::S | Op::H)` is true when `x` is either operator. A CIGAR of clips only has no
/// such operator, and `unwrap_or(false)` answers `false`, as the reference's `null` compares.
fn opens_with_insertion(read: &BamRecord) -> bool {
    read.cigar
        .elements
        .iter()
        .find(|element| !matches!(element.op, Op::S | Op::H))
        .map(|element| element.op == Op::I)
        .unwrap_or(false)
}

/// The leading and trailing hard clip lengths, zero where that end is any other operator: the
/// reference's `leadingHardClips` and `trailingHardClips` (`ReadPosRankSumTest.java:73-74`).
///
/// `.first()` is an `Option`, so `.filter(..)` keeps it only for an `H`, `.map(..)` takes its
/// length, and `.unwrap_or(0)` supplies the zero.
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

// `macro_rules!` writes the same `impl InfoFieldAnnotation` block once per member, with `$name`
// and `$key` filled in: the code the four Java subclasses inherit from `RankSumTest.annotate`.
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
        // The double is `1.23449999999999993...`, and Java rounds the shortest digits
        // `1.2345`, not that expansion.
        assert_eq!(format_three_decimals(1.2345), "1.235");
        assert_eq!(format_three_decimals(0.0), "0.000");
        assert_eq!(format_three_decimals(-1.0), "-1.000");
        assert_eq!(format_three_decimals(9.9999), "10.000");
    }
}
