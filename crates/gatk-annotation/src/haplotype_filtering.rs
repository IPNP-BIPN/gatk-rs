//! `HaplotypeFilteringAnnotation`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator.HaplotypeFilteringAnnotation`
//! (GATK 4.6.2.0), with the `JumboInfoAnnotation` interface it is declared against.
//!
//! ```java
//! final Map<String, Object> result = new HashMap<>();
//! result.put(GATKVCFConstants.HAPLOTYPES_BEFORE_FILTERING_KEY, haplotypeLikelihoods.alleles().size());
//! result.put(GATKVCFConstants.HAPLOTYPES_FILTERED_KEY, haplotypeLikelihoods.getFilteredHaplotypeCount());
//! ```
//!
//! Two counts, and the whole annotation is those two lines. What is worth explaining is not the
//! arithmetic but why this one was the last portable annotation of the fifty-four rather than an
//! early one, and what its two inputs mean.
//!
//! # It reads the haplotype matrix, and only that one
//!
//! Every other ported annotation is an `InfoFieldAnnotation`, handed the read-by-allele matrix.
//! This one is a `JumboInfoAnnotation`, which takes three matrices, and it ignores the first two.
//! That is what made it wait: not jmath, and not the annotator, but the existence of a likelihood
//! matrix whose allele axis is a `Haplotype` (see [`gatk_engine::haplotype`]) and of the
//! `filteredHaplotypeCount` field on it.
//!
//! # `ASSEMBLED_HAPS` counts the alleles of a matrix, not the haplotypes the assembler built
//!
//! `alleles().size()` is the size of an `IndexedSet` under `Haplotype.equals`, which includes the
//! uniqueness value. So it is the number of haplotypes **the matrix has rows for** after any
//! duplicate the set swallowed, and the key's name ("before filtering") is relative to the
//! filtering step that wrote the second count, not to the assembler.
//!
//! # `FILTERED_HAPS` is read from a field, so a matrix nobody filtered reports zero
//!
//! `getFilteredHaplotypeCount()` is a getter over an `int` field that only
//! `AlleleFiltering.subsetHaplotypesByAlleles` ever writes:
//!
//! ```java
//! readLikelihoods.setFilteredHaplotypeCount(readLikelihoods.numberOfAlleles() - subsettedReadLikelihoodsFinal.numberOfAlleles());
//! ```
//!
//! That writer is the flow-based `HaplotypeCaller`'s allele filtering, which belongs to Milestone
//! G3. Until it is ported the field is whatever the caller set, and its default is `0` rather than
//! absent: this annotation writes `FILTERED_HAPS=0` where a reader might expect no key at all.
//! Both keys are always written, since there is no guard anywhere in the method.
//!
//! # The map is a `HashMap`, and this is the one place that is visible here
//!
//! Every other annotation in this crate builds a `LinkedHashMap` or a singleton. This one builds a
//! `HashMap`, so the reference's own iteration order over the two keys is decided by
//! `String.hashCode` rather than by the order of the two `put` calls. Nothing observable depends
//! on it, because `VariantAnnotatorEngine` copies the result into a `LinkedHashMap` that the
//! encoder then sorts; the returned vector below is in `put` order, which is the order the two
//! statements are written in and not the order a debugger would show.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::fragment::Fragment;
use gatk_engine::haplotype::Haplotype;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::info_annotation::AnnotationValue;

/// `GATKVCFConstants.HAPLOTYPES_BEFORE_FILTERING_KEY`. The constant's name and its value differ,
/// which is why the value is written here rather than paraphrased.
pub const HAPLOTYPES_BEFORE_FILTERING_KEY: &str = "ASSEMBLED_HAPS";

/// `GATKVCFConstants.HAPLOTYPES_FILTERED_KEY`.
pub const HAPLOTYPES_FILTERED_KEY: &str = "FILTERED_HAPS";

/// The haplotype-typed matrix a `JumboInfoAnnotation` is handed, which is one of **two** types.
///
/// ```java
/// haplotypeLikelihoods.isPresent() ? haplotypeLikelihoods.get() : readHaplotypeAlleleLikelihoods.get()
/// ```
///
/// `VariantAnnotatorEngine.addInfoAnnotations` holds two optionals, `Optional<AlleleLikelihoods
/// <Fragment, Haplotype>>` and `Optional<AlleleLikelihoods<GATKRead, Haplotype>>`, and passes
/// whichever is present into a parameter declared `AlleleLikelihoods<? extends Locatable,
/// Haplotype>`. The wildcard erases the difference in Java; here the two instantiations are
/// different types, so the branch the ternary takes is named instead of erased.
///
/// The annotation reads nothing that distinguishes them, which is the point: naming the branch
/// costs nothing today and stops a later annotation from silently assuming the evidence is reads.
pub enum HaplotypeLikelihoods<'a> {
    /// `haplotypeLikelihoods`, the fragment-grouped matrix, taken when it is present.
    ByFragment(&'a AlleleLikelihoods<Fragment, Haplotype>),
    /// `readHaplotypeAlleleLikelihoods`, taken only when the fragment-grouped one is absent.
    ByRead(&'a AlleleLikelihoods<BamRecord, Haplotype>),
}

impl HaplotypeLikelihoods<'_> {
    /// `alleles().size()`.
    pub fn number_of_alleles(&self) -> usize {
        match self {
            HaplotypeLikelihoods::ByFragment(likelihoods) => likelihoods.number_of_alleles(),
            HaplotypeLikelihoods::ByRead(likelihoods) => likelihoods.number_of_alleles(),
        }
    }

    /// `getFilteredHaplotypeCount()`.
    pub fn filtered_haplotype_count(&self) -> i32 {
        match self {
            HaplotypeLikelihoods::ByFragment(likelihoods) => likelihoods.filtered_haplotype_count(),
            HaplotypeLikelihoods::ByRead(likelihoods) => likelihoods.filtered_haplotype_count(),
        }
    }
}

/// `JumboInfoAnnotation`: an INFO annotation that is handed more than the read matrix.
///
/// `fragment_likelihoods` is `Option` because the call site passes an explicit `null` for it:
///
/// ```java
/// fragmentLikelihoods.isPresent() ? fragmentLikelihoods.get() : null,
/// ```
///
/// The haplotype matrix is not optional at the same call site, because the guard above it
/// (`(fragmentLikelihoods.isPresent() && haplotypeLikelihoods.isPresent()) ||
/// readHaplotypeAlleleLikelihoods.isPresent()`) guarantees one of the two is there. It is not
/// optional here either, so an implementation cannot be written to guard something the engine
/// never hands it.
pub trait JumboInfoAnnotation {
    /// `getKeyNames()`, in declaration order.
    fn key_names(&self) -> Vec<&'static str>;

    /// `annotate(ref, features, vc, likelihoods, fragmentLikelihoods, haplotypeLikelihoods)`.
    ///
    /// `features` is the `FeatureContext` the engine passes; no ported implementation reads it, so
    /// it is carried as a unit for now rather than typed against a context this crate would
    /// otherwise not depend on. When an implementation needs it, it becomes the real type here.
    fn annotate(
        &self,
        reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
        fragment_likelihoods: Option<&AlleleLikelihoods<Fragment>>,
        haplotype_likelihoods: &HaplotypeLikelihoods<'_>,
    ) -> Vec<(String, AnnotationValue)>;
}

/// `HaplotypeFilteringAnnotation`.
pub struct HaplotypeFilteringAnnotation;

impl JumboInfoAnnotation for HaplotypeFilteringAnnotation {
    fn key_names(&self) -> Vec<&'static str> {
        vec![HAPLOTYPES_BEFORE_FILTERING_KEY, HAPLOTYPES_FILTERED_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        _vc: &VariantContext,
        _likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
        _fragment_likelihoods: Option<&AlleleLikelihoods<Fragment>>,
        haplotype_likelihoods: &HaplotypeLikelihoods<'_>,
    ) -> Vec<(String, AnnotationValue)> {
        vec![
            (
                HAPLOTYPES_BEFORE_FILTERING_KEY.to_string(),
                // `alleles().size()` is an `int` in Java and is boxed as an `Integer`; the cast is
                // where a matrix with more than 2^31 alleles would overflow, which is also where
                // the reference's own `int` would.
                AnnotationValue::Int(haplotype_likelihoods.number_of_alleles() as i32),
            ),
            (
                HAPLOTYPES_FILTERED_KEY.to_string(),
                AnnotationValue::Int(haplotype_likelihoods.filtered_haplotype_count()),
            ),
        ]
    }
}
