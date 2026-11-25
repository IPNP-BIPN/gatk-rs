//! `PerAlleleAnnotation` and its four members, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.7.0.0; the five Java files
//! are unchanged since 4.6.2.0).
//!
//! `MBQ`, `MMQ`, `MPOS` and `MFRL`: the median, per allele, of one number taken from each read
//! that best supports it. One shared shape, four values, and four different answers to the same
//! question about an allele no read supports.
//!
//! # The empty case is a different number in each of the four
//!
//! A summary of the four `aggregate` bodies (`BaseQuality.java:29`, `MappingQuality.java:28`,
//! `ReadPosition.java:36`, `FragmentLength.java:27`), not a quote:
//!
//! ```text
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
//! # An allele the matrix does not hold is a `NullPointerException`
//!
//! ```java
//! final Map<Allele, List<Integer>> values = likelihoods.alleles().stream()
//!         .collect(Collectors.toMap(a -> a, a -> new ArrayList<>()));
//! ...
//! vc.getAlleles().stream().filter(this::includeAllele).mapToInt(a -> aggregate(values.get(a)))
//! ```
//!
//! The buckets come from the **matrix's** alleles and the loop runs over the **variant's**, so a
//! variant allele the matrix never held looks up `null` and `aggregate` dereferences it. Not an
//! empty list, which would have been the value-for-no-reads: the annotation throws, and the
//! golden says so. Lines 62-65 of `tests/data/per_allele.txt.gz` (`allele-missing-from-matrix`,
//! one per annotation) are the suite's only exception rather than a value.
//!
//! # The median is commons-math3, not the obvious one
//!
//! `MathUtils.median` goes through commons-math3 `Percentile` and finishes with
//! `FastMath.round`, which is `(long) floor(x + 0.5)` and not `Math.round`. Both are in
//! [`jmath`]; this reaches the first. See htsjdk-rs decision 0023. `MathUtils.java:162` (one
//! line in the source, wrapped here):
//!
//! ```java
//! return (int) FastMath.round(new Median().evaluate(
//!         Arrays.stream(values).mapToDouble(n -> n).toArray()));
//! ```
//!
//! # Reading this file without reading Rust
//!
//! - `trait PerAlleleAnnotation` plays the part of the Java abstract class; each member is a
//!   field-less `struct`, and `impl PerAlleleAnnotation for X` holds its overrides. A trait
//!   method with a body is a default a member may override.
//! - `Option<i32>` is `OptionalInt`: `Some(v)` present, `None` empty.
//! - `Result<T, MissingBucket>` is how this port models a method that can throw: `Ok(value)` is a
//!   normal return, `Err(MissingBucket)` is the `NullPointerException`.
//! - `let Some(x) = expr else { continue; };` unwraps a present value into `x` or runs the
//!   `else` branch; there is no Java equivalent.
//! - `&BamRecord` and `&[i32]` are borrowed, read-only views: the function reads them and the
//!   caller keeps ownership. `&mut` (in `iter_mut`) is the writable form.
//! - `x as i32` is Java's `(int)` cast.
//! - `|a| ...` is a lambda, and `.iter().filter(..).map(..).collect()` is Java's stream chain.
//! - A function's last expression, written without `;`, is its return value.

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
    /// `getVcfKey()`: the one INFO key this member writes.
    fn vcf_key(&self) -> &'static str;

    /// `includeRefAllele()`, `false` in the parent and overridden in three of the four.
    fn include_ref_allele(&self) -> bool {
        false
    }

    /// `getValueForRead`. `None` is the reference's `OptionalInt.empty()`, which contributes
    /// nothing without stopping the traversal.
    fn value_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<i32>;

    /// What `aggregate` answers for an empty list: the constant in front of the `:` in each
    /// member's `aggregate` (0, 60, 50, 0; see the module note).
    fn value_for_no_reads(&self) -> i32;

    /// `aggregate(values)`: the median, or the annotation's own value for no reads.
    ///
    /// Each Java member writes the same body with its own constant, as `MappingQuality.java:28`:
    ///
    /// ```java
    /// return values.isEmpty() ? VALUE_FOR_NO_READS : MathUtils.median(Ints.toArray(values));
    /// ```
    ///
    /// so it is written once here and the constant is [`PerAlleleAnnotation::value_for_no_reads`].
    /// `MathUtils.median(int[])` uses `new Median()`, whose estimation type is `LEGACY`.
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
///
/// `<A: PerAlleleAnnotation>` makes this generic over any member, Java's
/// `<A extends PerAlleleAnnotation>`. `Ok(vec![])` is `Collections.emptyMap()`, and
/// `Err(MissingBucket)` is the exception the module note describes.
pub fn annotate<A: PerAlleleAnnotation>(
    annotation: &A,
    _reference: Option<&ReferenceContext>,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Result<Vec<(String, AnnotationValue)>, MissingBucket> {
    let Some(likelihoods) = likelihoods else {
        return Ok(Vec::new());
    };

    // `Collectors.toMap` over `likelihoods.alleles()`, so the buckets are the matrix's alleles and
    // not the variant's. An allele the matrix does not hold gets no bucket. A list of pairs
    // searched in order stands in for the `Map`; only lookups by allele are made, so the order
    // is not observable.
    let mut values: Vec<(Allele, Vec<i32>)> = (0..likelihoods.number_of_alleles())
        .filter_map(|index| likelihoods.get_allele(index).cloned())
        .map(|allele| (allele, Vec::new()))
        .collect();

    // `.filter(ba -> ba.isInformative() && isUsableRead(ba.evidence))`, then `getValueForRead`.
    for best in likelihoods.best_alleles_breaking_ties(None) {
        if !best.is_informative() {
            continue;
        }
        // The reference's `ba.evidence` is the read itself; here it is an index into the
        // sample's reads, looked up again (`and_then` chains a lookup that may be absent). The
        // sample is one of this matrix's own, so `unwrap_or(0)` does not decide anything.
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
        // `values.get(ba.allele).add(v)`. A best allele is chosen by index from the matrix's own
        // alleles, so its bucket exists; `iter_mut` hands out writable access to each pair, and
        // `bucket.1` is the pair's second element, the list.
        if let Some(bucket) = values.iter_mut().find(|(a, _)| *a == allele) {
            bucket.1.push(value);
        }
    }

    // `vc.getAlleles()`, filtered by `includeAllele`, in the variant's order: the matrix's order
    // does not decide this. A variant allele the matrix never held has no bucket at all, not an
    // empty one, so it is the exception below rather than the value for no reads
    // (`per_allele.txt.gz` lines 62-65).
    let mut statistics: Vec<AnnotationValue> = Vec::new();
    for allele in vc
        .alleles
        .iter()
        .filter(|allele| !allele.is_reference() || annotation.include_ref_allele())
    {
        // `values.get(a)` is null for an allele the matrix does not hold, and `aggregate`
        // dereferences it. The refusal is the reference's, not the port's.
        let Some((_, collected)) = values.iter().find(|(a, _)| a == allele) else {
            return Err(MissingBucket);
        };
        statistics.push(AnnotationValue::Int(annotation.aggregate(collected)));
    }

    Ok(vec![(
        annotation.vcf_key().to_string(),
        // `ImmutableMap.of(key, int[])`: an int array, whatever its length, and never a scalar.
        AnnotationValue::List(statistics),
    )])
}

/// The `NullPointerException` an allele outside the matrix produces.
///
/// It carries no message of its own: the reference's is the JVM's helpful-NPE text, which names
/// the local variable rather than the allele. `class` and `message` return the two strings the
/// golden records (`per_allele.txt.gz` lines 62-65).
///
/// `#[derive(..)]` has the compiler write the equality, copy and debug-print methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MissingBucket;

impl MissingBucket {
    pub fn class(&self) -> &'static str {
        "java.lang.NullPointerException"
    }

    pub fn message(&self) -> &'static str {
        "Cannot invoke \"java.util.List.isEmpty()\" because \"values\" is null"
    }
}

/// `BaseQuality`: `MBQ`, the median base quality at the variant's start.
pub struct BaseQuality;

impl BaseQuality {
    /// `BaseQuality.getBaseQuality`: `BaseQualityRankSumTest.getReadBaseQuality`, then
    /// `FastMath.round` and a narrowing cast. `BaseQualityHistogram` calls it too.
    ///
    /// The guard is the reference's own, `BaseQuality.java:44`: a read that starts after the
    /// variant, or ends before it, declines before the base is looked up at all.
    ///
    /// ```java
    /// if (vc.getStart() < read.getStart() || read.getEnd() < vc.getStart()) {
    /// ```
    pub fn base_quality(read: &BamRecord, vc: &VariantContext) -> Option<i32> {
        let start = vc.start as i32;
        if start < read_utils::start(read) || read_utils::end(read) < start {
            return None;
        }
        match read_utils::read_base_quality_at_reference_coordinate(read, start) {
            BaseAt::Present(quality) => {
                // The reference boxes the byte into an OptionalDouble and rounds it back, which
                // is exact for any whole number. The byte is Java's signed `byte` and this one is
                // an unsigned `u8`, so the two agree for qualities up to 127 only.
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
    /// `ReadPosition.getPosition`: `ReadPosRankSumTest.getReadPosition` behind the same guard as
    /// `BaseQuality.getBaseQuality` (`ReadPosition.java:48`), then an `(int)` cast of the double.
    ///
    /// The body repeats `read_position` in `rank_sum.rs` in `i32` rather than calling it; both
    /// follow `ReadPosRankSumTest.java:53-80`, where the hard clips are described.
    pub fn position(read: &BamRecord, vc: &VariantContext) -> Option<i32> {
        let start = vc.start as i32;
        if start < read_utils::start(read) || read_utils::end(read) < start {
            return None;
        }
        // `getReadPosition`'s leading-insertion case: a read that opens with an insertion and
        // starts one base past the variant's end is answered with 0. It cannot fire here, in the
        // reference either: a variant's end is never before its start, so such a read starts
        // after `vc.start` and the guard above has already returned. `ReadPosRankSum` is where
        // it is reached.
        if read_utils::start(read) == vc.stop as i32 + 1 && opens_with_insertion(read) {
            return Some(0);
        }
        // `getReadIndexForReferenceCoordinate(read, ...)` walks from the SOFT start. The pair's
        // first element is the index (`getLeft()`), and a negative one is not found.
        let (index, _) = read_utils::read_index_for_read(read, start);
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
///
/// `matches!(x, Op::S | Op::H)` is true when `x` is either clip; `find` gives an `Option`, and
/// `unwrap_or(false)` answers for a CIGAR of clips only, as the reference's `orElse(null)` does.
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
        // `OptionalInt.of(Math.abs(read.getFragmentLength()))` (`FragmentLength.java:40`).
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
                // The trait cannot carry the reference's exception, so a caller that needs to
                // distinguish it calls `annotate` directly. An empty map here is the absence of
                // keys, which is not what the reference does, and the suite tests the fallible
                // form rather than this one.
                annotate(self, reference, vc, likelihoods).unwrap_or_default()
            }
        }
    };
}

info_annotation_for!(BaseQuality, MEDIAN_BASE_QUALITY_KEY);
info_annotation_for!(MappingQuality, MEDIAN_MAPPING_QUALITY_KEY);
info_annotation_for!(ReadPosition, MEDIAN_READ_POSITION_KEY);
info_annotation_for!(FragmentLength, MEDIAN_FRAGMENT_LENGTH_KEY);
