//! `SampleList`, ported from `org.broadinstitute.hellbender.tools.walkers.annotator.SampleList`.
//!
//! ```java
//! if ( vc.isMonomorphicInSamples() || !vc.hasGenotypes() ) { return Collections.emptyMap(); }
//! for ( final Genotype genotype : vc.getGenotypesOrderedByName() ) {
//!     if ( genotype.isCalled() && !genotype.isHomRef() ) { ... append name ... }
//! }
//! ```
//!
//! Four decisions, none of which follows from the name "list of polymorphic samples".
//!
//! The guard is `isMonomorphicInSamples`, which is **not** "every genotype is hom-ref". A site with
//! alternate alleles and no genotypes at all is not monomorphic, so it survives the first disjunct
//! and is caught by the second; a site where every sample is a no-call *is* monomorphic and is
//! dropped. Both facts are settled in htsjdk-rs's `vcf-genotype-type` suite.
//!
//! The filter is `isCalled() && !isHomRef()`, so a partially called genotype, which is `MIXED` and
//! therefore called, is listed. A genotype that is `HET` only because it holds the reference bases
//! flagged as an alternate is listed too, and prints as `A/A` in the record beside it.
//!
//! The iteration order is `getGenotypesOrderedByName`, which is `String.compareTo` over the sample
//! names and not the order the genotypes are stored in.
//!
//! And the value is a single **String** with commas in it, under a header line declared
//! `UNBOUNDED` and `String`. It looks like a list and is not one, so a consumer splitting on the
//! comma is doing something the annotation never promised.

use htsjdk_vcf::genotype_type::{genotypes_ordered_by_name, is_monomorphic_in_samples};
use htsjdk_vcf::variant::VariantContext;

use gatk_engine::context::ReferenceContext;

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `GATKVCFConstants.SAMPLE_LIST_KEY`. The key is `Samples`, capitalised, which is unusual for a
/// VCF INFO key and is what the file carries.
pub const SAMPLE_LIST_KEY: &str = "Samples";

pub struct SampleList;

impl InfoFieldAnnotation for SampleList {
    fn key_names(&self) -> Vec<&'static str> {
        vec![SAMPLE_LIST_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
    ) -> Vec<(String, AnnotationValue)> {
        if is_monomorphic_in_samples(vc) || vc.genotypes.is_empty() {
            return Vec::new();
        }

        let mut samples = String::new();
        for genotype in genotypes_ordered_by_name(vc) {
            if genotype.is_called() && !genotype.is_hom_ref() {
                if !samples.is_empty() {
                    samples.push(',');
                }
                samples.push_str(&genotype.sample_name);
            }
        }

        // Reachable: a polymorphic site whose only non-hom-ref genotypes are filtered out of the
        // called counts still passes the guard, and then nothing matches the filter here.
        if samples.is_empty() {
            return Vec::new();
        }

        vec![(SAMPLE_LIST_KEY.to_string(), AnnotationValue::Str(samples))]
    }
}
