//! The walkers, ported from `org.broadinstitute.hellbender.engine` (GATK 4.6.2.0).
//!
//! A walker is the traversal a tool inherits: which records reach `apply`, in what order, and
//! what context each one arrives with. It lives above [`gatk_engine`] and [`gatk_readfilter`]
//! rather than inside either, because a traversal is a data source and a filter chain composed,
//! and neither half can see the other.

pub mod add_original_alignment_tags;
pub mod clip_reads;
pub mod fix_misencoded_base_quality_reads;
pub mod interval_walker;
pub mod left_align_indels;
pub mod locus_walker;
pub mod multi_pass;
pub mod print_distant_mates;
pub mod print_reads;
pub mod read_walker;
pub mod revert_base_quality_scores;
pub mod sam_output;
pub mod split_reads;
pub mod unmark_duplicates;
