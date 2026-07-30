//! `Coverage`, `MappingQualityZero` and `CountNs`, the three annotations that read the likelihood
//! matrix and nothing else, ported from `org.broadinstitute.hellbender.tools.walkers.annotator`.
//!
//! Their arithmetic is three counts. What separates them is the guard, the Java type, and what
//! they count *over*, and no two of the three agree on any of it.
//!
//! # Three guards, three different questions
//!
//! ```java
//! // Coverage
//! if (likelihoods == null || likelihoods.evidenceCount() == 0) { return Collections.emptyMap(); }
//! // MappingQualityZero
//! if (!vc.isVariant() || likelihoods == null) { return Collections.emptyMap(); }
//! // CountNs
//! if ( likelihoods == null ) { return Collections.emptyMap(); }
//! ```
//!
//! An empty matrix makes `Coverage` write nothing and `MappingQualityZero` write `0`. A
//! non-variant site makes `MappingQualityZero` write nothing and the other two write their count.
//! A port that factored the guard out would be wrong on both.
//!
//! # Three Java types for three counts
//!
//! `Coverage` and `MappingQualityZero` both go through `String.format("%d", ...)`, so their values
//! are **Strings**; `CountNs` puts the `long` straight into an `ImmutableMap`, so its value is a
//! **Long**. All three render identically and are three different objects to a consumer.
//!
//! # `CountNs` asks the read for a base at the variant's start
//!
//! ```java
//! final Optional<Byte> readBase = ReadUtils.getReadBaseAtReferenceCoordinate(read, vc.getStart());
//! return readBase.isPresent() && readBase.get() == 'N';
//! ```
//!
//! The comparison is against the byte `'N'`, upper case only, so a lower-case `n` is not counted.
//! And the lookup's own asymmetry applies: a base inside a soft clip cannot be reached even though
//! the index walk knows where it is, so a read whose `N` was clipped does not count.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::read::mapping_quality;
use gatk_engine::read_utils::{read_base_at_reference_coordinate, BaseAt};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `VCFConstants.DEPTH_KEY`.
pub const DEPTH_KEY: &str = "DP";
/// `VCFConstants.MAPPING_QUALITY_ZERO_KEY`.
pub const MAPPING_QUALITY_ZERO_KEY: &str = "MQ0";
/// `GATKVCFConstants.N_COUNT_KEY`.
pub const N_COUNT_KEY: &str = "NCount";

/// Every piece of evidence in the matrix, in sample then evidence order.
fn all_evidence(likelihoods: &AlleleLikelihoods<BamRecord>) -> Vec<&BamRecord> {
    (0..likelihoods.number_of_samples())
        .filter_map(|sample| likelihoods.sample_evidence(sample))
        .flatten()
        .collect()
}

/// `Coverage`: `DP`, the total evidence count, as a `String`.
pub struct Coverage;

impl InfoFieldAnnotation for Coverage {
    fn key_names(&self) -> Vec<&'static str> {
        vec![DEPTH_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        _vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        let Some(likelihoods) = likelihoods else {
            return Vec::new();
        };
        // The only one of the three that treats an empty matrix as nothing to say.
        if likelihoods.evidence_count() == 0 {
            return Vec::new();
        }
        vec![(
            DEPTH_KEY.to_string(),
            // String.format("%d", depth): a String, not an Integer.
            AnnotationValue::Str(likelihoods.evidence_count().to_string()),
        )]
    }
}

/// `MappingQualityZero`: `MQ0`, the reads with a mapping quality of zero, as a `String`.
pub struct MappingQualityZero;

impl InfoFieldAnnotation for MappingQualityZero {
    fn key_names(&self) -> Vec<&'static str> {
        vec![MAPPING_QUALITY_ZERO_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        // The only one of the three that asks about the site, and it does not ask about the
        // evidence count, so an empty matrix at a variant site writes a zero.
        if !vc.is_variant() {
            return Vec::new();
        }
        let Some(likelihoods) = likelihoods else {
            return Vec::new();
        };
        let mq0 = all_evidence(likelihoods)
            .into_iter()
            .filter(|read| mapping_quality(read) == 0)
            .count();
        vec![(
            MAPPING_QUALITY_ZERO_KEY.to_string(),
            AnnotationValue::Str(mq0.to_string()),
        )]
    }
}

/// `CountNs`: `NCount`, the reads whose base at the variant's start is `N`, as a `Long`.
pub struct CountNs;

impl CountNs {
    /// `doesReadHaveN`.
    pub fn does_read_have_n(read: &BamRecord, vc: &VariantContext) -> bool {
        matches!(
            read_base_at_reference_coordinate(read, vc.start as i32),
            BaseAt::Present(b'N')
        )
    }
}

impl InfoFieldAnnotation for CountNs {
    fn key_names(&self) -> Vec<&'static str> {
        vec![N_COUNT_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        let Some(likelihoods) = likelihoods else {
            return Vec::new();
        };
        let count = all_evidence(likelihoods)
            .into_iter()
            .filter(|read| CountNs::does_read_have_n(read, vc))
            .count();
        vec![(
            N_COUNT_KEY.to_string(),
            // ImmutableMap.of(N_COUNT_KEY, Count): a boxed Long, not a String and not an Integer.
            AnnotationValue::Long(count as i64),
        )]
    }
}
