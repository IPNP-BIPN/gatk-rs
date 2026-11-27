//! `SampleList`, ported from `org.broadinstitute.hellbender.tools.walkers.annotator.SampleList`
//! (GATK 4.7.0.0), over htsjdk 5.0.0's genotype predicates. The body (indentation shortened):
//!
//! ```java
//! if ( vc.isMonomorphicInSamples() || !vc.hasGenotypes() ) {
//!     return Collections.emptyMap();
//! }
//!
//! final StringBuilder samples = new StringBuilder();
//! for ( final Genotype genotype : vc.getGenotypesOrderedByName() ) {
//!     if ( genotype.isCalled() && !genotype.isHomRef() ){
//!         if ( samples.length() > 0 ) {
//!             samples.append(",");
//!         }
//!         samples.append(genotype.getSampleName());
//!     }
//! }
//! ```
//!
//! Four decisions, none of which follows from the name "list of polymorphic samples".
//!
//! The guard is `isMonomorphicInSamples`, which is **not** "every genotype is hom-ref". A site with
//! alternate alleles and no genotypes at all is not monomorphic, so it survives the first disjunct
//! and is caught by the second; a site where every sample is a no-call *is* monomorphic and is
//! dropped. Both facts are settled in htsjdk-rs's `vcf-genotype-type` suite. htsjdk's test is
//! (line broken to fit)
//!
//! ```java
//! monomorphic = !isVariant() ||
//!         (hasGenotypes() && getCalledChrCount(getReference()) == getCalledChrCount());
//! ```
//!
//! and both counts skip filtered genotypes (`if (!g.isFiltered()) {`). The filter below does not
//! look at the filter field, so a filtered non-reference genotype is listed when an unfiltered one
//! makes the site polymorphic, and when only filtered genotypes carry an alternate the guard
//! fires: row `filtered-het` of `tests/data/annotations.txt.gz`, a filtered `0/1` beside a `0/0`,
//! is empty.
//!
//! The filter is `isCalled() && !isHomRef()`, so a partially called genotype, which is `MIXED` and
//! therefore called, is listed. A genotype that is `HET` only because it holds the reference bases
//! flagged as an alternate is listed too, and prints as `A/A` in the record beside it. Row
//! `mixed-and-hom-ref` gives `Samples=s1,s2` where `s1` is `0/.`: a genotype holding nothing but
//! the reference and a no-call is listed as non-reference.
//!
//! The iteration order is `getGenotypesOrderedByName`, which is `String.compareTo` over the sample
//! names and not the order the genotypes are stored in: htsjdk sorts the names with
//! `Collections.sort(sampleNamesInOrder);`. Row `name-order` stores `b`, `A`, `a`, `B` (the last
//! hom-ref) and writes `Samples=A,a,b`.
//!
//! And the value is a single **String** with commas in it, under a header line declared
//! `UNBOUNDED` and `String`. It looks like a list and is not one, so a consumer splitting on the
//! comma is doing something the annotation never promised.
//!
//! # Reading this file without reading Rust
//!
//! - `String::new()`, `push` and `push_str` are a `StringBuilder` and its `append`; the
//!   `String` grows in place.
//! - `&genotype.sample_name` lends the name to `push_str` without copying or giving it away.

use htsjdk_vcf::genotype_type::{genotypes_ordered_by_name, is_monomorphic_in_samples};
use htsjdk_vcf::variant::VariantContext;

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use htsjdk_bam::record::BamRecord;

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `GATKVCFConstants.SAMPLE_LIST_KEY`. The key is `Samples`, capitalised, which is unusual for a
/// VCF INFO key and is what the file carries.
pub const SAMPLE_LIST_KEY: &str = "Samples";

/// `SampleList`: `Samples`, the called non-hom-ref samples in name order. See the module note.
pub struct SampleList;

impl InfoFieldAnnotation for SampleList {
    fn key_names(&self) -> Vec<&'static str> {
        vec![SAMPLE_LIST_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        _likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
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

        // `if ( samples.length() == 0 ) {`. After the guard this is a corner case: a polymorphic
        // site has an unfiltered genotype carrying a called non-reference allele, and that
        // genotype is called and not `HOM_REF`, so it is listed. It is reached when the listed
        // names are all the empty string, or when a genotype holds a reference-flagged allele
        // with other bases than the site's reference: `determineType` calls it `HOM_REF`
        // (`firstCallAllele.isReference()`) while `countAllele` does not count it as the
        // reference, since `Allele.equals` compares bases.
        if samples.is_empty() {
            return Vec::new();
        }

        vec![(SAMPLE_LIST_KEY.to_string(), AnnotationValue::Str(samples))]
    }
}
