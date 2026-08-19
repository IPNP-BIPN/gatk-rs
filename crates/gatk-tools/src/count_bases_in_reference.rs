//! Ported from `org.broadinstitute.hellbender.tools.walkers.fasta.CountBasesInReference`
//! (GATK 4.6.2.0).
//!
//! The first tool on [`crate::reference_walker`], and about as small as a GATK tool gets: one
//! `long[256]` indexed by the byte at each locus, printed in byte order.
//!
//! # What it counts is not what is in the file
//!
//! `referenceContext.getBase()` comes through `CachingIndexedFastaSequenceFile`, which upper-cases
//! and replaces every IUPAC ambiguity code with `N`. So a FASTA holding `acgtRYKMSWBDHV` counts as
//! four bases and ten `N`s, and the table has at most five rows however many symbols the file uses.
//! Soft-masking is invisible for the same reason.
//!
//! # The table is sparse and ordered by byte
//!
//! Only counts above zero are printed, and the loop runs 0..256, so the order is the ASCII order of
//! the bases rather than the order they were seen: `A`, `C`, `G`, `N`, `T`.
//!
//! The reference prints to standard output AND, when `-O` is given, writes the same text to the
//! file -- `print` rather than `println`, so there is no trailing newline beyond the one each row
//! already carries.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::interval_args::IntervalArguments;
use gatk_engine::reference::ReferenceFileSource;

use crate::reference_walker::{self, TraversalError};

/// `baseCounts`, one slot per byte value.
pub struct BaseCounts {
    counts: [i64; 256],
}

impl Default for BaseCounts {
    fn default() -> Self {
        BaseCounts { counts: [0; 256] }
    }
}

impl BaseCounts {
    /// The count for one byte, which is zero for every byte never seen.
    pub fn get(&self, base: u8) -> i64 {
        self.counts[base as usize]
    }

    /// `onTraversalSuccess`: `<char> : <count>\n` for every byte seen at least once, in byte order.
    pub fn report(&self) -> String {
        let mut text = String::new();
        for (byte, count) in self.counts.iter().enumerate() {
            if *count > 0 {
                text.push_str(&format!("{} : {}\n", byte as u8 as char, count));
            }
        }
        text
    }
}

/// `doWork`: traverse and count.
///
/// The window is the default one, the locus itself, so `getBase()` is the single base at it.
pub fn run(
    reference: &mut ReferenceFileSource,
    arguments: &IntervalArguments,
) -> Result<BaseCounts, TraversalError> {
    let applied =
        reference_walker::traverse(reference, arguments, |locus: &SimpleInterval| locus.clone())?;
    let mut counts = BaseCounts::default();
    for call in applied {
        // `getBase()` is the first byte of the window, which at the default window is its only one.
        counts.counts[call.bases[0] as usize] += 1;
    }
    Ok(counts)
}
