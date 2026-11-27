//! `ChromosomeCounts`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator.ChromosomeCounts` (GATK 4.7.0.0) and,
//! through it, htsjdk 5.0.0's `VariantContextUtils.calculateChromosomeCounts`.
//!
//! The annotation itself is three statements, because everything it does happens in htsjdk
//! (`ChromosomeCounts.java` lines 47 to 52; the last line broken to fit):
//!
//! ```java
//! Utils.nonNull(vc);
//! if ( ! vc.hasGenotypes() ) {
//!     return Collections.emptyMap();
//! }
//!
//! return VariantContextUtils.calculateChromosomeCounts(vc, new LinkedHashMap<>(), true,
//!         Collections.emptySet());
//! ```
//!
//! Two of those three statements are decisions rather than plumbing.
//!
//! `removeStaleValues` is hard-coded `true`, so a site where nobody is called loses `AN`, `AC` and
//! `AF` entirely rather than reporting zeroes. The `false` branch of that function is unreachable
//! through this annotation and reachable through the library, which is why the library keeps it.
//! The golden row `anno ChromosomeCounts all-no-call` in `tests/data/annotations.txt.gz` is empty.
//!
//! `founderIds` is hard-coded **empty**, which the library reads as "the whole cohort". The
//! pedigree behaviour of `AF`, where the numerator and denominator come from the founders while the
//! `AC` beside it comes from everyone, therefore never fires here. It is measured in htsjdk-rs's
//! `vcf-chromosome-counts` suite, since a caller that passes a real founder set reaches it.
//!
//! And the `hasGenotypes()` guard gives the same output as leaving it out. Without it the library
//! is entered with no genotypes, `getCalledChrCount()` is 0, and with `removeStaleValues` on it
//! returns the `LinkedHashMap` it was handed, still empty:
//!
//! ```java
//! if (AN == 0 && removeStaleValues) {
//! ```
//!
//! The two routes return different objects (`Collections.emptyMap()` and that map) with no keys
//! in either. The htsjdk-rs port reports the second route as `remove_stale: true`, which this
//! file does not read, since a fresh map has nothing stale to remove.
//!
//! # `AC` and `AF` change Java type with the number of alternate alleles
//!
//! From `calculateChromosomeCounts` (indentation shortened):
//!
//! ```java
//! attributes.put(
//!   VCFConstants.ALLELE_COUNT_KEY, alleleCounts.size() == 1 ? alleleCounts.get(0) : alleleCounts);
//! ```
//!
//! and the same for `AF`. So one alternate allele puts an `Integer` (and a `Double`), two or more
//! put an `ArrayList`; the golden rows `one-het` and `two-alts` of `annotations.txt.gz` show
//! `AC=1[java.lang.Integer]` and `AC=(2,1)[java.util.ArrayList]`. A site with no alternate allele
//! keeps only `AN` (row `ref-only-site`, `AN=2[java.lang.Integer]`).
//!
//! Filtered genotypes are left out of every count, because `getCalledChrCount` skips them
//! (`if (!g.isFiltered()) {`): row `filtered-het` has a filtered `0/1` beside a `0/0` and reports
//! `AN=2`, `AC=0`.
//!
//! # Reading this file without reading Rust
//!
//! - `pub struct ChromosomeCounts;` is a type with no fields, a stateless class; the methods come
//!   from the `impl InfoFieldAnnotation for ChromosomeCounts` block.
//! - `if let Some(x) = y { ... }` runs the block only when the `Option` `y` holds a value, with
//!   that value bound to `x`. It is a null check and a cast in one.
//! - `Count::One(value) => ...` inside a `match` is one case of an enum, with its payload bound
//!   to `value`.

use htsjdk_vcf::chromosome_counts::{
    calculate_chromosome_counts, Count, Frequency, ALLELE_COUNT_KEY, ALLELE_FREQUENCY_KEY,
    ALLELE_NUMBER_KEY,
};
use htsjdk_vcf::variant::VariantContext;

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use htsjdk_bam::record::BamRecord;

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `ChromosomeCounts.keyNames`, in declaration order. The record's own order is the encoder's,
/// which sorts.
pub const KEY_NAMES: [&str; 3] = [ALLELE_NUMBER_KEY, ALLELE_COUNT_KEY, ALLELE_FREQUENCY_KEY];

/// `ChromosomeCounts`: `AN`, `AC` and `AF` from the called genotypes. See the module note.
pub struct ChromosomeCounts;

impl InfoFieldAnnotation for ChromosomeCounts {
    fn key_names(&self) -> Vec<&'static str> {
        vec![ALLELE_NUMBER_KEY, ALLELE_COUNT_KEY, ALLELE_FREQUENCY_KEY]
    }

    /// `ChromosomeCounts.annotate`. The reference context and the likelihoods are unused, as in
    /// the Java; a leading `_` on a parameter name tells the compiler that is deliberate.
    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        _likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        // `if ( ! vc.hasGenotypes() )`; in htsjdk `hasGenotypes()` is
        // `return !genotypes.isEmpty();`.
        if vc.genotypes.is_empty() {
            return Vec::new();
        }

        // removeStaleValues = true and founderIds = empty, both hard-coded by the annotation.
        let counts = calculate_chromosome_counts(vc, true, &[]);

        // `remove_stale` is not read: see the module note. Each `None` field is a key the library
        // did not put, so it is simply not pushed.
        let mut out = Vec::new();
        if let Some(an) = counts.allele_number {
            out.push((ALLELE_NUMBER_KEY.to_string(), AnnotationValue::Int(an)));
        }
        if let Some(ac) = counts.allele_count {
            out.push((
                ALLELE_COUNT_KEY.to_string(),
                match ac {
                    // One alternate allele gives an Integer, two or more give an ArrayList. The
                    // two render the same and are not the same object to a consumer.
                    Count::One(value) => AnnotationValue::Int(value),
                    // `.map(AnnotationValue::Int)` passes the case's constructor as the function,
                    // wrapping each count as an `Integer` element of the list.
                    Count::Many(values) => AnnotationValue::List(
                        values.into_iter().map(AnnotationValue::Int).collect(),
                    ),
                },
            ));
        }
        if let Some(af) = counts.allele_frequency {
            out.push((
                ALLELE_FREQUENCY_KEY.to_string(),
                match af {
                    Frequency::One(value) => AnnotationValue::Double(value),
                    Frequency::Many(values) => AnnotationValue::List(
                        values.into_iter().map(AnnotationValue::Double).collect(),
                    ),
                },
            ));
        }
        out
    }
}
