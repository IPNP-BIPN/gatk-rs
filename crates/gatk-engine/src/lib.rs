//! The GATK engine, ported from `org.broadinstitute.hellbender.engine` and the utilities it rests
//! on (GATK 4.6.2.0).
//!
//! The first thing here is intervals, because they are what `-L` means and what every walker
//! restricts itself to. A tool given the wrong interval reads the wrong data and every number it
//! produces is wrong in a way no downstream comparison can attribute.

pub mod activity_profile;
pub mod alignment_state;
pub mod alignment_utils;
pub mod allele_likelihoods;
pub mod allele_list;
pub mod assembly_region;
pub mod assembly_region_iterator;
pub mod assembly_region_walker;
pub mod baq;
pub mod base_recalibration_engine;
pub mod base_utils;
pub mod bqsr_transformer;
pub mod cigar_builder;
pub mod cigar_utils;
pub mod clipping;
pub mod concordance_walker;
pub mod contamination_tables;
pub mod context;
pub mod context_iterator;
pub mod covariates;
pub mod downsampling;
pub mod feature_intervals;
pub mod features;
pub mod fisher_exact;
pub mod fragment;
pub mod gatk_report;
pub mod genotype_index;
pub mod haplotype;
pub mod histogram;
pub mod interval;
pub mod interval_args;
pub mod java_format;
pub mod java_hash;
pub mod java_random;
pub mod java_regex;
pub mod jexl;
pub mod locus_iterator;
pub mod locus_shards;
pub mod mann_whitney;
pub mod math_utils;
pub mod natural_log_utils;
pub mod overhang_fixing_manager;
pub mod permutation;
pub mod persistence_optimizer;
pub mod pileup;
pub mod pileup_summary;
pub mod qual_quantizer;
pub mod read;
pub mod read_group;
pub mod read_pileup;
pub mod read_states;
pub mod read_utils;
pub mod reads;
pub mod recal_datum;
pub mod recal_utils;
pub mod recalibration_report;
pub mod recalibration_tables;
pub mod reference;
pub mod sa_tag;
pub mod sam_pileup;
pub mod somatic_likelihoods;
pub mod subset_alleles;
pub mod tranches;
pub mod tsv_table;
pub mod variant_context_utils;
pub mod variant_getters;
pub mod variant_source;
pub mod well19937c;
