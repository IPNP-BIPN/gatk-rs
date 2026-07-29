//! The GATK engine, ported from `org.broadinstitute.hellbender.engine` and the utilities it rests
//! on (GATK 4.6.2.0).
//!
//! The first thing here is intervals, because they are what `-L` means and what every walker
//! restricts itself to. A tool given the wrong interval reads the wrong data and every number it
//! produces is wrong in a way no downstream comparison can attribute.

pub mod alignment_state;
pub mod cigar_builder;
pub mod cigar_utils;
pub mod clipping;
pub mod context;
pub mod interval;
pub mod interval_args;
pub mod read;
pub mod read_group;
pub mod read_utils;
pub mod reads;
pub mod reference;
