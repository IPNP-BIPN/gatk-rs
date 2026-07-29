//! The walkers, ported from `org.broadinstitute.hellbender.engine` (GATK 4.6.2.0).
//!
//! A walker is the traversal a tool inherits: which records reach `apply`, in what order, and
//! what context each one arrives with. It lives above [`gatk_engine`] and [`gatk_readfilter`]
//! rather than inside either, because a traversal is a data source and a filter chain composed,
//! and neither half can see the other.

pub mod interval_walker;
pub mod print_reads;
pub mod read_walker;
