//! The GATK engine, ported from `org.broadinstitute.hellbender.engine` and the utilities it rests
//! on (GATK 4.6.2.0).
//!
//! The first thing here is intervals, because they are what `-L` means and what every walker
//! restricts itself to. A tool given the wrong interval reads the wrong data and every number it
//! produces is wrong in a way no downstream comparison can attribute.

pub mod activity_profile;
pub mod alignment_state;
pub mod allele_likelihoods;
pub mod allele_list;
pub mod assembly_region;
pub mod assembly_region_iterator;
pub mod assembly_region_walker;
pub mod cigar_builder;
pub mod cigar_utils;
pub mod clipping;
pub mod context;
pub mod context_iterator;
pub mod downsampling;
pub mod features;
pub mod interval;
pub mod interval_args;
pub mod java_hash;
pub mod java_random;
pub mod jexl;
pub mod locus_iterator;
pub mod locus_shards;
pub mod permutation;
pub mod pileup;
pub mod read;
pub mod read_group;
pub mod read_pileup;
pub mod read_states;
pub mod read_utils;
pub mod reads;
pub mod reference;
pub mod variant_getters;
pub mod variant_source;
pub mod well19937c;
