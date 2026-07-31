//! GATK's variant annotations, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.6.2.0).
//!
//! An annotation is a small function from a variant context to a handful of INFO or FORMAT keys.
//! What makes the archetype worth its own crate is not the arithmetic, which is usually a count,
//! but three things that are invisible in the field names:
//!
//! - **the Java type of the value put in the map is observable.** `Coverage` puts a `String` built
//!   with `String.format("%d", depth)`, `CountNs` puts a boxed `Long`, and `ChromosomeCounts` puts
//!   an `Integer` or an `ArrayList` depending on the alternate count. The encoder renders each
//!   differently in the edge cases, so [`AnnotationValue`] keeps the distinction rather than
//!   flattening everything to a string;
//! - **an annotation that has nothing to say returns an empty map, not a zero.** The key is then
//!   absent from the record. Every one of these annotations has at least one guard that reaches
//!   that branch, and the guards are not the obvious ones;
//! - **the key order in the returned map is not the order in the record.** The encoder sorts, so
//!   `getKeyNames()` is the declaration order and nothing else.

pub mod allele_specific_rank_sum;
pub mod allele_specific_strand_bias;
pub mod chromosome_counts;
pub mod coverage;
pub mod depth_per_allele;
pub mod flow;
pub mod heterozygosity;
pub mod info_annotation;
pub mod mapping_quality;
pub mod original_alignment;
pub mod per_allele;
pub mod rank_sum;
pub mod raw_gt_count;
pub mod read_grouping;
pub mod sample_list;
pub mod site_statistics;
pub mod strand_bias;
pub mod tandem_repeat;

pub use info_annotation::{AnnotationValue, InfoFieldAnnotation};
