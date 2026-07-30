//! `ChromosomeCounts`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator.ChromosomeCounts`.
//!
//! The annotation itself is three lines, because everything it does happens in htsjdk:
//!
//! ```java
//! Utils.nonNull(vc);
//! if ( ! vc.hasGenotypes() ) { return Collections.emptyMap(); }
//! return VariantContextUtils.calculateChromosomeCounts(vc, new LinkedHashMap<>(), true, Collections.emptySet());
//! ```
//!
//! Two of those three lines are decisions rather than plumbing.
//!
//! `removeStaleValues` is hard-coded `true`, so a site where nobody is called loses `AN`, `AC` and
//! `AF` entirely rather than reporting zeroes. The `false` branch of that function is unreachable
//! through this annotation and reachable through the library, which is why the library keeps it.
//!
//! `founderIds` is hard-coded **empty**, which the library reads as "the whole cohort". The
//! pedigree behaviour of `AF`, where the numerator and denominator come from the founders while the
//! `AC` beside it comes from everyone, therefore never fires here. It is measured in htsjdk-rs's
//! `vcf-chromosome-counts` suite, since a caller that passes a real founder set reaches it.
//!
//! And the `hasGenotypes()` guard is not redundant with the library's own: without it the function
//! would be entered with an empty genotype list, and with `removeStaleValues` on that returns the
//! removal signal rather than an empty map. Both end with no keys written, by different routes.

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

pub struct ChromosomeCounts;

impl InfoFieldAnnotation for ChromosomeCounts {
    fn key_names(&self) -> Vec<&'static str> {
        vec![ALLELE_NUMBER_KEY, ALLELE_COUNT_KEY, ALLELE_FREQUENCY_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        _likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        if vc.genotypes.is_empty() {
            return Vec::new();
        }

        // removeStaleValues = true and founderIds = empty, both hard-coded by the annotation.
        let counts = calculate_chromosome_counts(vc, true, &[]);

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
