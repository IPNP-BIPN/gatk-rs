//! `OriginalAlignment`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator.OriginalAlignment`, with the `OA` tag
//! accessors of `AddOriginalAlignmentTags`.
//!
//! It counts the reads supporting the **best** alternate allele whose original alignment was on a
//! different contig, which is how a Mutect2 call on the mitochondrion is checked against a NuMT.
//! Four of its decisions are somewhere other than in this class.
//!
//! # The allele it counts for is chosen by `TLOD`, and a missing `TLOD` is a `-1`
//!
//! ```java
//! final double[] lods = Mutect2FilteringEngine.getTumorLogOdds(vc);
//! if (lods == null) { warning.warn(...); return Collections.emptyMap(); }
//! final int indexOfMaxLod = MathUtils.maxElementIndex(lods);
//! final Allele altAlelle = vc.getAlternateAllele(indexOfMaxLod);
//! ```
//!
//! Absent, the annotation says nothing and logs once. **Present but `.`**, the getter turns the
//! missing element into `-1` and the conversion into `-ln(10)`, which is an ordinary number that
//! can win the maximum, so a site whose only stated `TLOD` is missing still gets an allele chosen
//! for it. That is settled in the `variant-getters` suite.
//!
//! `maxElementIndex` gives a tie to the first element, so the earliest alternate allele wins.
//!
//! # The contig comparison is a string field, not a contig
//!
//! ```java
//! public static String getOAContig(final GATKRead read) {
//!     return read.getAttributeAsString(OA_TAG_NAME).split(",")[0];
//! }
//! ```
//!
//! The `OA` value is `contig,start,strand,cigar,mapq,nm;`, and the writer replaces any comma
//! **inside the contig name** with an underscore before joining, so a round trip through the tag
//! renames such a contig and the comparison then fails against the reference's own name. An
//! unmapped read is written as `*,0,*,*,0,0;`, so its original contig is the string `*`, which
//! differs from every real contig and therefore counts.
//!
//! # `isInformative` is the log10 threshold whatever base the matrix is in
//!
//! The filter takes `ba.isInformative()`, which compares the confidence against `0.2` even after
//! `switchToNaturalLog`. That is measured in the `genotyper-allele-likelihoods` suite; here it
//! decides which reads are counted.
//!
//! # `Utils.nonNull(likelihoods)` rather than a guard
//!
//! Unlike `Coverage`, `MappingQualityZero` and `CountNs`, which all treat null likelihoods as
//! nothing to say, this one throws.

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
pub fn oa_contig(read: &BamRecord) -> Option<&str> {
    oa_value(read).map(|value| value.split(OA_SEPARATOR).next().unwrap_or(value))
}

pub struct OriginalAlignment;

impl OriginalAlignment {
    /// The count, with the failures the trait cannot express.
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
        let Some(index) = max_element_index(&lods) else {
            return Ok(None);
        };
        // `getAlternateAllele(index)` past the end is an index error in the reference; a caller
        // reaching it has more TLOD entries than alternate alleles.
        let Some(alt_allele) = vc.alternate_alleles().get(index) else {
            return Ok(None);
        };

        let mut count = 0i64;
        for best in likelihoods.best_alleles_breaking_ties(None) {
            let evidence = likelihoods
                .sample_evidence(likelihoods.index_of_sample(&best.sample).unwrap_or(0))
                .and_then(|reads| reads.get(best.evidence_index));
            let Some(read) = evidence else { continue };
            // The four conditions in the reference's order: the tag is there, the call is
            // informative, the allele is the best alternate, and the original contig differs.
            if oa_value(read).is_none() {
                continue;
            }
            if !best.is_informative() {
                continue;
            }
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
                // ImmutableMap.of(key, nonChrMAlt), where nonChrMAlt is a long.
                AnnotationValue::Long(count),
            )],
            _ => Vec::new(),
        }
    }
}
