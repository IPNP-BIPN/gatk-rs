//! The walkers, ported from `org.broadinstitute.hellbender.engine` (GATK 4.6.2.0).
//!
//! A walker is the traversal a tool inherits: which records reach `apply`, in what order, and
//! what context each one arrives with. It lives above [`gatk_engine`] and [`gatk_readfilter`]
//! rather than inside either, because a traversal is a data source and a filter chain composed,
//! and neither half can see the other.

pub mod add_original_alignment_tags;
pub mod analyze_covariates;
pub mod annotate_vcf_with_bam_depth;
pub mod annotate_vcf_with_expected_allele_fraction;
pub mod apply_bqsr;
pub mod apply_vqsr;
pub mod base_recalibrator;
pub mod calculate_contamination;
pub mod calculate_mixing_fractions;
pub mod check_pileup;
pub mod clip_reads;
pub mod compare_base_qualities;
pub mod concordance;
pub mod convert_headerless_shard;
pub mod count_bases_in_reference;
pub mod count_false_positives;
pub mod count_reads;
pub mod count_variants;
pub mod counting_walkers;
pub mod dump_tabix_index;
pub mod evaluate_info_field_concordance;
pub mod filter_mutect_calls;
pub mod filter_variant_tranches;
pub mod fix_misencoded_base_quality_reads;
pub mod gather_vcfs;
pub mod get_sample_name;
pub mod interval_walker;
pub mod left_align_and_trim_variants;
pub mod left_align_indels;
pub mod locus_walker;
pub mod methylation_type_caller;
pub mod multi_pass;
pub mod post_process_reads_for_rsem;
pub mod print_distant_mates;
pub mod print_reads;
pub mod print_reads_header;
pub mod read_anonymizer;
pub mod read_walker;
pub mod reference_walker;
pub mod remove_nearby_indels;
pub mod revert_base_quality_scores;
pub mod sam_output;
pub mod select_variants;
pub mod split_n_cigar_reads;
pub mod split_reads;
pub mod transfer_read_tags;
pub mod unmark_duplicates;
pub mod update_vcf_sequence_dictionary;
pub mod validate_variants;
pub mod variant_filtration;
pub mod variants_to_table;
