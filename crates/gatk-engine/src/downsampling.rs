//! Ported from `org.broadinstitute.hellbender.utils.downsampling.ReservoirDownsampler` and
//! `org.broadinstitute.hellbender.utils.downsampling.LevelingDownsampler` (GATK 4.6.2.0).
//!
//! The two draw from **different** static generators, which is the first thing to get right:
//! `ReservoirDownsampler` calls `Utils.getRandomGenerator()`, a `java.util.Random`
//! ([`crate::java_random`]), while `LevelingDownsampler` reaches
//! `Utils.getRandomDataGenerator()`, a `RandomDataGenerator` over a `Well19937c`
//! ([`crate::well19937c`]). Routing both through one stream would be wrong twice over: wrong
//! values, and wrong positions for every later consumer of either stream.
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
use crate::well19937c::Well19937c;

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

/// `LevelingDownsampler`: given several stacks and a total target, remove items evenly until the
/// sum fits.
///
/// Four behaviours decide the answer:
///
///  * **the plan is computed on sizes alone, before anything is removed.** A round-robin walk over
///    the stack sizes decrements one at a time, so the *distribution* of removals across stacks is
///    fixed by arithmetic and only which items go is random;
///  * **the walk stops when no stack can give**, counted by consecutive refusals rather than by a
///    scan, so a stack that becomes unmodifiable and then modifiable again is reconsidered;
///  * **each stack is reduced by a full `nextPermutation`**, which shuffles all its items and keeps
///    the head. That costs `size - 1` draws from the shared `Well19937c` per stack reduced, in
///    stack order, so the stacks' *order* changes which items every later stack keeps;
///  * **a stack that keeps everything takes no draw at all**, because `downsampleOneGroup` returns
///    before sampling when `numItemsToKeep >= group.size()`.
///
/// With `min_elements_per_stack` of zero a stack can be planned down to zero items, and then
/// `nextPermutation(n, 0)` throws rather than returning nothing. That refusal is reproduced, not
/// smoothed over: see [`crate::permutation::PermutationError::NotStrictlyPositive`].
pub struct LevelingDownsampler<T> {
    target_size: i64,
    min_elements_per_stack: usize,
    groups: Vec<Vec<T>>,
    groups_are_finalized: bool,
    discarded: usize,
}

impl<T> LevelingDownsampler<T> {
    /// `new LevelingDownsampler(targetSize)`, whose default `minElementsPerStack` is 1.
    pub fn new(target_size: i64) -> Self {
        LevelingDownsampler::with_minimum(target_size, 1)
    }

    /// `new LevelingDownsampler(targetSize, minElementsPerStack)`.
    pub fn with_minimum(target_size: i64, min_elements_per_stack: usize) -> Self {
        assert!(
            target_size >= 0,
            "targetSize must be >= 0 but got {target_size}"
        );
        LevelingDownsampler {
            target_size,
            min_elements_per_stack,
            groups: Vec::new(),
            groups_are_finalized: false,
            discarded: 0,
        }
    }

    /// `submit(T item)`.
    pub fn submit(&mut self, group: Vec<T>) {
        self.groups.push(group);
    }

    /// `size()`: the sum over the stacks, not the number of stacks.
    pub fn size(&self) -> usize {
        self.groups.iter().map(Vec::len).sum()
    }

    pub fn discarded(&self) -> usize {
        self.discarded
    }

    pub fn has_finalized_items(&self) -> bool {
        self.groups_are_finalized && !self.groups.is_empty()
    }

    pub fn has_pending_items(&self) -> bool {
        !self.groups_are_finalized && !self.groups.is_empty()
    }

    /// `signalEndOfInput()`, which is where the levelling actually happens: nothing is removed
    /// until every stack has been submitted, because the plan needs all the sizes.
    pub fn signal_end_of_input(
        &mut self,
        rng: &mut Well19937c,
    ) -> Result<(), crate::permutation::PermutationError> {
        let result = self.level_groups(rng);
        self.groups_are_finalized = true;
        result
    }

    /// `consumeFinalizedItems()`.
    pub fn consume_finalized_items(&mut self) -> Vec<Vec<T>> {
        if !self.has_finalized_items() {
            return Vec::new();
        }
        let groups = std::mem::take(&mut self.groups);
        self.groups_are_finalized = false;
        groups
    }

    /// `levelGroups()`.
    fn level_groups(
        &mut self,
        rng: &mut Well19937c,
    ) -> Result<(), crate::permutation::PermutationError> {
        let mut group_sizes: Vec<usize> = self.groups.iter().map(Vec::len).collect();
        let total_size: usize = group_sizes.iter().sum();

        if total_size as i64 <= self.target_size {
            return Ok(());
        }

        let mut to_remove = total_size as i64 - self.target_size;
        let mut current = 0usize;
        let mut consecutive_unmodifiable = 0usize;

        while to_remove > 0 && consecutive_unmodifiable < group_sizes.len() {
            if group_sizes[current] > self.min_elements_per_stack {
                group_sizes[current] -= 1;
                to_remove -= 1;
                consecutive_unmodifiable = 0;
            } else {
                consecutive_unmodifiable += 1;
            }
            current = (current + 1) % group_sizes.len();
        }

        // Reduced in submission order, which is also the order the draws are taken in.
        for (group, keep) in self.groups.iter_mut().zip(group_sizes) {
            self.discarded += downsample_one_group(group, keep, rng)?;
        }
        Ok(())
    }
}

/// `downsampleOneGroup`, returning what it counts through `incrementNumberOfDiscardedItems`.
///
/// The kept items stay in the stack's own order: the permutation is used as set membership, not as
/// an ordering, so shuffled indices produce an unshuffled result.
fn downsample_one_group<T>(
    group: &mut Vec<T>,
    keep: usize,
    rng: &mut Well19937c,
) -> Result<usize, crate::permutation::PermutationError> {
    if keep >= group.len() {
        // No draw is taken here at all.
        return Ok(0);
    }
    let indices = crate::permutation::sample_indices_without_replacement(
        group.len() as i32,
        keep as i32,
        rng,
    )?;
    let mut wanted = vec![false; group.len()];
    for index in indices {
        wanted[index as usize] = true;
    }
    let mut keeper = wanted.iter();
    let before = group.len();
    group.retain(|_| *keeper.next().unwrap_or(&false));
    Ok(before - group.len())
}
