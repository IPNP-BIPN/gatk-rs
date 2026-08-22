//! `DownsampleByDuplicateSet`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.consensus.DownsampleByDuplicateSet` and
//! `org.broadinstitute.hellbender.engine.DuplicateSetWalker` (GATK 4.6.2.0).
//!
//! Whole molecules dropped rather than reads, so that a mixture keeps its family-size
//! distribution. The tool is a hundred and five lines over a hundred and fifty of walker, and the
//! walker is where the surprises are.
//!
//! # The last duplicate set of the file escapes every rejection rule
//!
//! ```java
//! private void processLastReadSet(){
//!     if (currentReadsWithSameUMI.getReads().size() > 0){
//!         apply(currentReadsWithSameUMI, ...);
//!     }
//! }
//! ```
//!
//! No `rejectSet` call. A trailing set that is too small, or has an odd number of reads, or is
//! short on one strand, is offered to the tool anyway: at `--min-reads 4` over ten two-read
//! molecules, exactly one molecule is written, and it is the last.
//!
//! # A rejected set does not consume a random draw
//!
//! The rejection happens before `apply`, so a rejectable molecule at the FRONT of a file leaves
//! every later molecule with the decision it had before. Adding one anywhere the walker accepts it
//! would have shifted the whole sequence instead.
//!
//! # A set with an odd number of reads is rejected at the defaults
//!
//! `size() % 2 == 1` is one of the three rules, described in the source as checking that the set is
//! paired, so a three-read molecule is dropped whatever the minimums say.

use gatk_engine::java_random::JavaRandom;

/// `RANDOM_SEED`.
pub const RANDOM_SEED: i64 = 142;

/// The arguments, with the walker's defaults beside the tool's own.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Arguments {
    pub fraction_to_keep: f64,
    /// `DEFAULT_MINIMUM_READS_PER_SET`.
    pub minimum_reads: usize,
    /// `DEFAULT_MINIMUM_READS_PER_STRAND`.
    pub minimum_reads_per_strand: usize,
}

impl Arguments {
    /// The defaults, for a given fraction. `--fraction-to-keep` is required, so it has none.
    pub fn keeping(fraction_to_keep: f64) -> Self {
        Arguments {
            fraction_to_keep,
            minimum_reads: 1,
            minimum_reads_per_strand: 0,
        }
    }
}

/// One read, reduced to what the walker reads from it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub name: String,
    /// The molecule number, the part of `MI:Z:<number>/<strand>` before the slash.
    pub molecule: i32,
    /// The strand suffix, which only feeds the per-strand minimum.
    pub strand: String,
}

/// What the walk refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DownsampleError {
    /// A molecule number went backwards.
    NotSortedByMoleculeId,
}

impl DownsampleError {
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    pub fn message(&self) -> String {
        "The input bam must be sorted by the molecule ID (MI) tag.".to_string()
    }
}

/// `rejectSet`: the three rules, in the order the reference applies them.
pub fn reject_set(set: &[Read], arguments: &Arguments) -> bool {
    // `MoleculeID.countStrands`, which counts the two suffixes separately.
    let a = set.iter().filter(|read| read.strand == "A").count();
    let b = set.iter().filter(|read| read.strand == "B").count();
    if a.min(b) < arguments.minimum_reads_per_strand {
        return true;
    }
    // "Check that the read set is paired", which an odd count fails whatever the minimums say.
    if set.len() % 2 == 1 {
        return true;
    }
    set.len() < arguments.minimum_reads
}

/// The whole traversal: reads grouped by molecule number, each group offered or rejected, and the
/// last group offered without ever being asked.
///
/// Returns the reads written, in order.
pub fn run(reads: &[Read], arguments: &Arguments) -> Result<Vec<Read>, DownsampleError> {
    let mut rng = JavaRandom::new(RANDOM_SEED);
    let mut written = Vec::new();
    let mut current: Vec<Read> = Vec::new();

    let offer = |set: &[Read], rng: &mut JavaRandom, written: &mut Vec<Read>| {
        // One draw per set that reaches `apply`, consumed whether the set is kept or not.
        if rng.next_double() < arguments.fraction_to_keep {
            written.extend(set.iter().cloned());
        }
    };

    for read in reads {
        let held = match current.first() {
            None => {
                current.push(read.clone());
                continue;
            }
            Some(held) => held.molecule,
        };
        if held > read.molecule {
            return Err(DownsampleError::NotSortedByMoleculeId);
        }
        if held < read.molecule {
            // The end of a set. A rejected one costs no draw, which is why the order matters.
            if !reject_set(&current, arguments) {
                offer(&current, &mut rng, &mut written);
            }
            current = vec![read.clone()];
            continue;
        }
        current.push(read.clone());
    }

    // `processLastReadSet`, which never asks `rejectSet`.
    if !current.is_empty() {
        offer(&current, &mut rng, &mut written);
    }
    Ok(written)
}
