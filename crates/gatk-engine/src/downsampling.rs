//! Ported from `org.broadinstitute.hellbender.utils.downsampling.ReservoirDownsampler`
//! (GATK 4.6.2.0), over [`crate::java_random`].
//!
//! Reservoir sampling, which decides *which* reads a deep pileup keeps. Four of its behaviours
//! change the answer and none is inherent to the algorithm:
//!
//!  * **the draw happens even when the slot is discarded.** `nextInt(totalReadsSeen)` is called for
//!    every read past the target, and only then is the slot compared against the target. So a read
//!    that is thrown away still advances the shared random stream, and every later draw anywhere in
//!    the run depends on how many reads were seen here;
//!  * **the bound grows with the count**, so the same read at the same position draws differently
//!    depending on how many reads preceded it;
//!  * **replacement is in place.** The winner takes the loser's index rather than being appended,
//!    so the reservoir's order is the order slots were *last written*, not the reads' order;
//!  * **`expectFewOverflows` chooses the backing list**, and the reference converts it to an
//!    `ArrayList` on the first overflow. That is invisible in the output, and it is the reason the
//!    class carries a flag that changes nothing observable, which is worth knowing before someone
//!    ports it as if it did.
//!
//! `setNonRandomReplacementMode` replaces the draw with `Math.abs(name.hashCode()) % totalReadsSeen`,
//! which is deterministic without the generator, and which `Math.abs` makes wrong for
//! `Integer.MIN_VALUE`: that hash stays negative and the modulo is negative too.

use crate::java_random::JavaRandom;

/// `ReservoirDownsampler`.
pub struct ReservoirDownsampler<'a, T> {
    target_sample_size: usize,
    reservoir: Vec<&'a T>,
    total_seen: i32,
    discarded: usize,
    end_of_input: bool,
    non_random_replacement: bool,
}

/// How a replacement slot is chosen.
pub enum SlotSource<'r> {
    /// `Utils.getRandomGenerator().nextInt(totalReadsSeen)`, from the shared stream.
    Random(&'r mut JavaRandom),
    /// `Math.abs(read.getName().hashCode()) % totalReadsSeen`, with `Math.abs` left as it is.
    NonRandom,
}

impl<'a, T> ReservoirDownsampler<'a, T> {
    /// The reference throws on a non-positive target, so this refuses too.
    pub fn new(target_sample_size: usize) -> Self {
        assert!(
            target_sample_size > 0,
            "Cannot do reservoir downsampling with a sample size <= 0"
        );
        ReservoirDownsampler {
            target_sample_size,
            reservoir: Vec::with_capacity(target_sample_size),
            total_seen: 0,
            discarded: 0,
            end_of_input: false,
            non_random_replacement: false,
        }
    }

    /// `setNonRandomReplacementMode`.
    pub fn set_non_random_replacement_mode(&mut self, on: bool) {
        self.non_random_replacement = on;
    }

    pub fn size(&self) -> usize {
        self.reservoir.len()
    }

    pub fn discarded(&self) -> usize {
        self.discarded
    }

    /// `submit`, with the name needed only by the non-random mode.
    ///
    /// The order of the two steps is the behaviour: the draw is taken first and unconditionally,
    /// then the slot decides whether anything is replaced.
    pub fn submit(&mut self, item: &'a T, name: &str, slots: &mut SlotSource<'_>) {
        assert!(
            !self.end_of_input,
            "attempt to submit read after end of input stream has been signaled"
        );
        self.total_seen += 1;
        if (self.total_seen as usize) <= self.target_sample_size {
            self.reservoir.push(item);
            return;
        }
        let slot = match slots {
            SlotSource::Random(random) => random.next_int_bound(self.total_seen),
            // `Math.abs(Integer.MIN_VALUE)` is `Integer.MIN_VALUE`, so a name hashing to it gives a
            // negative slot, which is less than the target and indexes out of bounds upstream.
            // Reproduced as an arithmetic fact rather than corrected.
            SlotSource::NonRandom => {
                crate::java_hash::string_hash_code(name).wrapping_abs() % self.total_seen
            }
        };
        if slot >= 0 && (slot as usize) < self.target_sample_size {
            self.reservoir[slot as usize] = item;
        }
        self.discarded += 1;
    }

    /// `signalEndOfInput`.
    pub fn signal_end_of_input(&mut self) {
        self.end_of_input = true;
    }

    /// `consumeFinalizedItems`: the reservoir, and the downsampler reset.
    pub fn consume_finalized_items(&mut self) -> Vec<&'a T> {
        let items = if self.end_of_input {
            std::mem::take(&mut self.reservoir)
        } else {
            Vec::new()
        };
        self.clear_items();
        items
    }

    /// `clearItems`, which also resets the count and the end-of-input flag but *not* the discarded
    /// statistic, which `resetStats` clears separately.
    pub fn clear_items(&mut self) {
        self.reservoir = Vec::with_capacity(self.target_sample_size);
        self.total_seen = 0;
        self.end_of_input = false;
    }

    /// `resetStats`.
    pub fn reset_stats(&mut self) {
        self.discarded = 0;
    }
}
