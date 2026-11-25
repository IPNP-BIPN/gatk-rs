//! `OriginalAlignment`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator.OriginalAlignment`, with the `OA` tag
//! accessors of `AddOriginalAlignmentTags` (GATK 4.7.0.0).
//!
//! It counts the reads supporting the **best** alternate allele whose original alignment was on a
//! different contig, which is how a Mutect2 call on the mitochondrion is checked against a NuMT.
//! Four of its decisions are somewhere other than in this class.
//!
//! # The allele it counts for is chosen by `TLOD`, and a missing `TLOD` is a `-1`
//!
//! ```java
//! final double[] lods = Mutect2FilteringEngine.getTumorLogOdds(vc);
//! if (lods==null) {
//! ```
//!
//! whose branch logs a one-shot warning and returns `Collections.emptyMap()`; then
//!
//! ```java
//! final int indexOfMaxLod = MathUtils.maxElementIndex(lods);
//! final Allele altAlelle = vc.getAlternateAllele(indexOfMaxLod);
//! ```
//!
//! Absent, the annotation says nothing and logs once (`tests/data/original_alignment.txt.gz`, row
//! `no-tlod`, empty). This port does not log. **Present but `.`**, the getter turns the
//! missing element into `-1` and the conversion into `-ln(10)`, which is an ordinary number that
//! can win the maximum, so a site whose only stated `TLOD` is missing still gets an allele chosen
//! for it. That is settled in the `variant-getters` suite, and row `tlod-missing` here counts a
//! read (`OCM=1[java.lang.Long]`). The `-1` is `getAttributeAsDoubleArray`'s missing value,
//! `getAttributeAsDoubleArray(variantContext, attribute, () -> null, -1)`, and the conversion is
//! `log10ToLog`, `log10 * LOG_10`.
//!
//! `maxElementIndex` gives a tie to the first element, so the earliest alternate allele wins
//! (`if (array[i] > array[maxI])`, strictly greater; row `tlod-tie`).
//!
//! # The contig comparison is a string field, not a contig
//!
//! ```java
//! public static String getOAContig (final GATKRead read) {
//!     return(read.getAttributeAsString(OA_TAG_NAME).split(OA_SEPARATOR)[0]);
//! }
//! ```
//!
//! The `OA` value is `contig,start,strand,cigar,mapq,nm;`, and the writer replaces any comma
//! **inside the contig name** with an underscore before joining
//! (`read.getContig().replace(OA_SEPARATOR, "_"),`), so a round trip through the tag renames such
//! a contig and the comparison then fails against the reference's own name. An unmapped read is
//! written as `OAValue = "*,0,*,*,0,0;";`, so its original contig is the string `*`, which differs
//! from every real contig and therefore counts (row `unmapped-oa`, `OCM=1`).
//!
//! # `isInformative` is the log10 threshold whatever base the matrix is in
//!
//! The filter takes `ba.isInformative()`, which compares the confidence against `0.2` even after
//! `switchToNaturalLog`: its body is `return confidence > LOG_10_INFORMATIVE_THRESHOLD;`, and
//! `LOG_10_INFORMATIVE_THRESHOLD = 0.2`. That is measured in the `genotyper-allele-likelihoods`
//! suite; here it decides which reads are counted (row `uninformative`, `OCM=0`).
//!
//! # `Utils.nonNull(likelihoods)` rather than a guard
//!
//! Unlike `Coverage`, `MappingQualityZero` and `CountNs`, which all treat null likelihoods as
//! nothing to say, this one throws (row `null-likelihoods`,
//! `E:java.lang.IllegalArgumentException:Null object is not allowed here.`).
//! `OriginalAlignment::count` takes the matrix by reference, which cannot be null; the trait
//! method, whose argument is an `Option`, answers an empty map for `None` instead of throwing.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<Option<i64>, NonDoubleValue>` is three outcomes in one type: `Ok(Some(n))` a count,
//!   `Ok(None)` "nothing to say", `Err(e)` the exception the `TLOD` parse throws.
//! - `?` after a call that returns a `Result` hands an `Err` straight back to the caller, like an
//!   exception that is not caught; otherwise it unwraps the `Ok`.
//! - `let Some(x) = y else { return ...; };` unwraps an `Option` or leaves early; `continue`
//!   in its place skips to the next loop iteration.
//! - `match` with `Some(TagValue::Str(text)) => ...` and `_ => ...` is a type test and a cast in
//!   one pattern; `_` is every other case.
//! - `*b"OA"` is the two bytes `'O'`, `'A'` as a fixed-size array, a SAM tag name.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::variant_getters::{get_tumor_log_odds, max_element_index, NonDoubleValue};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_vcf::variant::VariantContext;

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY`.
pub const ORIGINAL_CONTIG_MISMATCH_KEY: &str = "OCM";

/// `AddOriginalAlignmentTags.OA_TAG_NAME`.
pub const OA_TAG: [u8; 2] = *b"OA";

/// `AddOriginalAlignmentTags.OA_SEPARATOR`.
pub const OA_SEPARATOR: char = ',';

/// `read.getAttributeAsString(OA_TAG_NAME)`.
///
/// Only a string-typed tag gives `Some`. The Java adapter's version ends with
/// `return attributeValue != null ? attributeValue.toString() : null;`, so an `OA` tag stored as
/// another type is still a string there; here it reads as absent.
pub fn oa_value(read: &BamRecord) -> Option<&str> {
    match read.tags.get(Tag::new(&OA_TAG)) {
        Some(TagValue::Str(text)) => Some(text.as_str()),
        _ => None,
    }
}

/// `AddOriginalAlignmentTags.getOAContig`: the first comma-separated field of the tag.
///
/// `None` only where the reference would have thrown on a null tag, which its caller avoids by
/// testing `hasAttribute` first.
///
/// `str::split` always yields at least one piece, the empty string for an empty input, so the
/// `unwrap_or` fallback is never taken. Java's `split` agrees on `""` and on a leading comma
/// (rows `oacontig empty` and `oacontig leading-comma` both give the empty string), but drops
/// trailing empty pieces, so a value of nothing but commas is an empty array there and `[0]`
/// throws, where this returns `Some("")`.
pub fn oa_contig(read: &BamRecord) -> Option<&str> {
    oa_value(read).map(|value| value.split(OA_SEPARATOR).next().unwrap_or(value))
}

pub struct OriginalAlignment;

impl OriginalAlignment {
    /// The body of `OriginalAlignment.annotate`: the count, with the failures the trait cannot
    /// express.
    ///
    /// `Err` is the `TLOD` parse failure; `Ok(None)` is the "no `TLOD`" case, which the reference
    /// reports by logging once and returning an empty map.
    pub fn count(
        vc: &VariantContext,
        likelihoods: &AlleleLikelihoods<BamRecord>,
        current_contig: &str,
    ) -> Result<Option<i64>, NonDoubleValue> {
        let Some(lods) = get_tumor_log_odds(vc)? else {
            return Ok(None);
        };
        // `None` for an empty `TLOD` list. The reference's `maxElementIndex` throws there:
        // `Utils.validateArg(array.length > 0, "array may not be empty");`.
        let Some(index) = max_element_index(&lods) else {
            return Ok(None);
        };
        // `getAlternateAllele(index)` past the end is an index error in the reference; a caller
        // reaching it has more TLOD entries than alternate alleles.
        let Some(alt_allele) = vc.alternate_alleles().get(index) else {
            return Ok(None);
        };

        // `likelihoods.bestAllelesBreakingTies()`, the no-argument form, hence `None` priorities.
        // The Java `BestAllele` holds the read itself (`ba.evidence`); this one holds the sample
        // and the read's index, so the read is looked up. The sample always comes from this
        // matrix, so the `unwrap_or(0)` fallback is not taken.
        let mut count = 0i64;
        for best in likelihoods.best_alleles_breaking_ties(None) {
            let evidence = likelihoods
                .sample_evidence(likelihoods.index_of_sample(&best.sample).unwrap_or(0))
                .and_then(|reads| reads.get(best.evidence_index));
            let Some(read) = evidence else { continue };
            // The four conditions in the reference's order: the tag is there, the call is
            // informative, the allele is the best alternate, and the original contig differs.
            // Java's `&&` stops at the first false one; each `continue` does the same.
            if oa_value(read).is_none() {
                continue;
            }
            if !best.is_informative() {
                continue;
            }
            // `ba.allele.equals(altAlelle)`. `allele` is an `Option`; `as_ref()` compares it
            // without taking it out of `best`.
            if best.allele.as_ref() != Some(alt_allele) {
                continue;
            }
            if oa_contig(read) != Some(current_contig) {
                count += 1;
            }
        }
        Ok(Some(count))
    }
}

impl InfoFieldAnnotation for OriginalAlignment {
    fn key_names(&self) -> Vec<&'static str> {
        vec![ORIGINAL_CONTIG_MISMATCH_KEY]
    }

    /// The reference reads the contig from `ref.getInterval().getContig()`, so a `None` reference
    /// context is the `NullPointerException` it would have thrown; here it is an empty map, which
    /// is the same absence of keys.
    ///
    /// The same holds for `None` likelihoods (the reference's `Utils.nonNull(likelihoods);`
    /// throws) and for an `Err` from [`OriginalAlignment::count`]: the `_` arm answers an empty map
    /// for both `Ok(None)` and the `TLOD` parse error, which `getTumorLogOdds` throws as a
    /// `GATKException` in the reference.
    fn annotate(
        &self,
        reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        let (Some(reference), Some(likelihoods)) = (reference, likelihoods) else {
            return Vec::new();
        };
        let Some(interval) = reference.interval() else {
            return Vec::new();
        };
        match OriginalAlignment::count(vc, likelihoods, &interval.contig) {
            Ok(Some(count)) => vec![(
                ORIGINAL_CONTIG_MISMATCH_KEY.to_string(),
                // ImmutableMap.of(key, nonChrMAlt), where nonChrMAlt is a long, so a boxed `Long`
                // (every count row of `original_alignment.txt.gz` is `[java.lang.Long]`).
                AnnotationValue::Long(count),
            )],
            _ => Vec::new(),
        }
    }
}
