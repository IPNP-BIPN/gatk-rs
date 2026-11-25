//! Ported from `org.broadinstitute.hellbender.utils.downsampling.ReservoirDownsampler` and
//! `org.broadinstitute.hellbender.utils.downsampling.LevelingDownsampler` (GATK 4.7.0.0; both Java
//! files are unchanged since 4.6.2.0), with `PositionalDownsampler` and `MutectDownsampler` from
//! the same package, each run to exhaustion over a sorted list rather than pushed one read at a
//! time.
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
//! `setNonRandomReplacementMode` replaces the draw with
//! `Math.abs(name.hashCode()) % totalReadsSeen`, which is deterministic without the generator, and
//! which `Math.abs` makes wrong for `Integer.MIN_VALUE`: that hash stays negative, and the modulo
//! is then negative or, when `totalReadsSeen` divides 2^31, zero.
//!
//! The reservoir goldens are `tests/data/reservoir.txt.gz` (`tests/reservoir.rs`), and the
//! levelling ones `tests/data/leveling.txt.gz` (`tests/leveling.rs`).
//!
//! # Reading this file without reading Rust
//!
//! - `ReservoirDownsampler<'a, T>` holds `&'a T`: **borrowed** references to reads the caller
//!   owns, where Java holds object references. `'a` is a lifetime, the compiler's proof that the
//!   reads outlive the reservoir; it changes nothing at run time. `.cloned()` turns the borrowed
//!   reads back into owned copies when they leave.
//! - `&mut T` is a writable borrow; only one may exist at a time, which is why the shared
//!   generator is passed down as `&mut JavaRandom` rather than held in a field.
//! - `enum SlotSource` is a type whose value is one of its variants, and a variant may carry data
//!   (`Random` carries the generator). `match` picks the arm for the variant it is handed, and the
//!   compiler refuses a `match` that leaves one out.
//! - `Option<T>` is "a value or nothing", Java's nullable reference with the check enforced;
//!   `Some(x)` is the value and `None` is `null`. `is_some_and(f)` is "non-null and `f` holds".
//! - `Result<(), E>` is "nothing, or the error `E`", how this codebase models a `void` method
//!   that throws; `?` after a call passes the error to the caller.
//! - `std::mem::take(&mut v)` moves the list out of the field and leaves an empty one behind,
//!   Java's "return the reference and allocate a fresh list".
//! - `|a, b| { ... }` is a closure, Java's lambda. `usize` is an unsigned index type Java lacks,
//!   and `x as usize` / `x as i32` convert between it and `int`.

use crate::java_random::JavaRandom;
use crate::well19937c::Well19937c;

/// `ReservoirDownsampler`.
///
/// `expectFewOverflows` and `isLinkedList` are not fields here: they choose between a `LinkedList`
/// and an `ArrayList` for the reservoir, which no output can tell apart (see the module note).
pub struct ReservoirDownsampler<'a, T> {
    /// `targetSampleSize`.
    target_sample_size: usize,
    /// `reservoir`.
    reservoir: Vec<&'a T>,
    /// `totalReadsSeen`, an `int` as upstream because it is the bound handed to `nextInt`.
    total_seen: i32,
    /// `Downsampler.numDiscardedItems`, the statistic `resetStats` clears.
    discarded: usize,
    /// `endOfInputStream`.
    end_of_input: bool,
    /// `nonRandomReplacementMode`.
    non_random_replacement: bool,
}

/// How a replacement slot is chosen: the two arms of the reference's
///
/// ```java
/// final int randomSlot = !nonRandomReplacementMode
///                             ? Utils.getRandomGenerator().nextInt(totalReadsSeen)
///                             : Math.abs(newRead.getName().hashCode()) % totalReadsSeen;
/// ```
///
/// The generator is a parameter rather than a static, so it is handed in with the variant that
/// uses it. `<'r>` is the lifetime of that borrow.
pub enum SlotSource<'r> {
    /// `Utils.getRandomGenerator().nextInt(totalReadsSeen)`, from the shared stream.
    Random(&'r mut JavaRandom),
    /// `Math.abs(read.getName().hashCode()) % totalReadsSeen`, with `Math.abs` left as it is.
    NonRandom,
}

impl<'a, T> ReservoirDownsampler<'a, T> {
    /// `new ReservoirDownsampler(targetSampleSize)`. The reference throws on a non-positive target,
    /// so this refuses too; a `usize` cannot be negative, so only zero reaches the check.
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

    /// `size()`: the reads held, pending or finalized.
    pub fn size(&self) -> usize {
        self.reservoir.len()
    }

    /// `getNumberOfDiscardedItems()`.
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
        // Counted before the comparison, so the first read past the target draws with bound
        // `target + 1`.
        self.total_seen += 1;
        if (self.total_seen as usize) <= self.target_sample_size {
            self.reservoir.push(item);
            return;
        }
        let slot = match slots {
            SlotSource::Random(random) => random.next_int_bound(self.total_seen),
            // `Math.abs(Integer.MIN_VALUE)` is `Integer.MIN_VALUE`, so a name hashing to it can
            // give a negative slot. `wrapping_abs` is that same `Math.abs`: plain `abs` would
            // panic on it in a debug build. Upstream the negative slot passes
            // `randomSlot < targetSampleSize` and `reservoir.set` then throws
            // `IndexOutOfBoundsException`; here the `slot >= 0` test below skips the replacement
            // instead. The arithmetic is the reference's, the outcome for that one hash is not.
            SlotSource::NonRandom => {
                crate::java_hash::string_hash_code(name).wrapping_abs() % self.total_seen
            }
        };
        // In place: the winner takes the loser's index, so the reservoir's order is the order its
        // slots were last written.
        if slot >= 0 && (slot as usize) < self.target_sample_size {
            self.reservoir[slot as usize] = item;
        }
        // Counted whether or not a slot was replaced: one read in, one read out, either way.
        self.discarded += 1;
    }

    /// `signalEndOfInput`.
    pub fn signal_end_of_input(&mut self) {
        self.end_of_input = true;
    }

    /// `consumeFinalizedItems`: the reservoir, and the downsampler reset.
    ///
    /// Before end of input the reference returns an empty list **without** clearing, so pending
    /// reads survive; this clears in every case. Every caller signals end of input first
    /// (`positional_downsample_with`, `mutect_downsample` and `read_states`), so the two never
    /// differ there.
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

/// `PositionalDownsampler` over a coordinate-sorted stream, run to exhaustion: the reads that
/// share an alignment start go through one reservoir of `target` slots, drawing from the shared
/// generator once a start holds more than `target` reads, and a read with no assigned position
/// passes through untouched. The positions are compared on the assigned contig and start, which
/// is what `ReadCoordinateComparator.compareCoordinates` asks.
///
/// The reference also throws when a read sorts before the previous one; the input here is taken
/// as sorted and not checked.
pub fn positional_downsample(
    reads: &[htsjdk_bam::record::BamRecord],
    target: usize,
    random: &mut JavaRandom,
) -> Vec<htsjdk_bam::record::BamRecord> {
    positional_downsample_with(reads, target, random, false)
}

/// The same with `nonRandomDownsamplingMode`, which `RampedHaplotypeCaller` asks for: the slot a
/// read takes is its name's hash modulo the count seen, and the generator is never drawn from.
pub fn positional_downsample_with(
    reads: &[htsjdk_bam::record::BamRecord],
    target: usize,
    random: &mut JavaRandom,
    non_random: bool,
) -> Vec<htsjdk_bam::record::BamRecord> {
    let mut out = Vec::with_capacity(reads.len());
    let mut reservoir: ReservoirDownsampler<'_, htsjdk_bam::record::BamRecord> =
        ReservoirDownsampler::new(target);
    let mut previous: Option<(i32, i32)> = None;
    // `finalizeReservoir`: signal end of input, take what survived, and reset the statistic so
    // the same reservoir serves the next position.
    let drain = |reservoir: &mut ReservoirDownsampler<'_, htsjdk_bam::record::BamRecord>,
                 out: &mut Vec<htsjdk_bam::record::BamRecord>| {
        reservoir.signal_end_of_input();
        out.extend(reservoir.consume_finalized_items().into_iter().cloned());
        reservoir.reset_stats();
    };
    for read in reads {
        let position = (read.reference_index, read.alignment_start);
        // `handlePositionalChange`, before the read is looked at. Two differences from upstream.
        // Every read with no contig compares equal there whatever its start, while these tuples
        // may differ; no such read enters the reservoir, so the extra drain finds it empty and
        // changes nothing. And upstream `finalizeReservoir(true)` throws
        // `ShouldNeverReachHereException` when the position moves on with an empty reservoir (a
        // read with a contig but start 0, then a placed read on that contig), where this drains
        // nothing and goes on: a divergence.
        if previous.is_some_and(|p| p != position) {
            drain(&mut reservoir, &mut out);
        }
        // `ReadUtils.readHasNoAssignedPosition`: no contig, or no start.
        if read.reference_index < 0 || read.alignment_start <= 0 {
            out.push(read.clone());
        } else {
            let mut slots = if non_random {
                SlotSource::NonRandom
            } else {
                SlotSource::Random(random)
            };
            reservoir.submit(read, &read.read_name, &mut slots);
        }
        previous = Some(position);
    }
    drain(&mut reservoir, &mut out);
    out
}

/// `MutectDownsampler` over a coordinate-sorted stream, run to exhaustion.
///
/// **Divergence**: the reference refuses a non-positive stride (`ParamUtils.isPositive(stride,
/// "stride must be > 0")`); this clamps it to one.
///
/// Reads are grouped in strides of `stride` bases from the first read of each; a stride holding
/// `max_suspicious_reads_per_alignment_start * stride` reads of mapping quality 50 or less is
/// dropped whole from the read that reaches that count on; and a stride with more than
/// `max_reads_per_alignment_start * stride` reads keeps a reservoir of that many among its reads of
/// mapping quality **above** 50 only, drawing from the shared generator. A zero or negative
/// maximum is no maximum. Each finished stride is sorted by start, as `consumeFinalizedItems`
/// sorts what it hands back.
///
/// **Why sorting per stride is the reference's sort**: `consumeFinalizedItems` sorts everything
/// finalized since the last call, and its caller, `PushToPullIterator.fillCache`, calls it as soon
/// as anything is finalized:
///
/// ```java
/// while ( inputElements.hasNext() && ! transformer.hasFinalizedItems() ) {
///     transformer.submit(inputElements.next());
/// }
/// ```
///
/// So a consumption is one finished stride or one read with no position, except at the end:
/// `fillCache` calls `signalEndOfInput()` before it consumes, so when the last read given opens
/// a stride of its own, that stride and the one it closed are sorted together. On one contig
/// that is the same order (every start in the later stride is larger). **Divergence**: when that
/// last read is on a new contig at a smaller start, the reference hands it over first and this
/// hands it over last.
pub fn mutect_downsample(
    reads: &[htsjdk_bam::record::BamRecord],
    max_reads_per_alignment_start: i32,
    max_suspicious_reads_per_alignment_start: i32,
    stride: i32,
    random: &mut JavaRandom,
) -> Vec<htsjdk_bam::record::BamRecord> {
    const SUSPICIOUS_MAPPING_QUALITY: u8 = 50;
    let stride = stride.max(1);
    let max_coverage = if max_reads_per_alignment_start <= 0 {
        usize::MAX
    } else {
        (max_reads_per_alignment_start * stride) as usize
    };
    let max_suspicious = if max_suspicious_reads_per_alignment_start <= 0 {
        usize::MAX
    } else {
        (stride * max_suspicious_reads_per_alignment_start) as usize
    };
    let mut out: Vec<htsjdk_bam::record::BamRecord> = Vec::with_capacity(reads.len());
    // `pendingReads`, `firstReadInStride` (as its contig and start), `rejectAllReadsInStride` and
    // `suspiciousReadCount`.
    let mut pending: Vec<&htsjdk_bam::record::BamRecord> = Vec::new();
    let mut first: Option<(i32, i32)> = None;
    let mut reject = false;
    let mut suspicious = 0usize;
    // `finalizePendingReads`. A closure over `max_coverage`; the state it changes is passed in,
    // because a closure that held those borrows itself would lock them for the whole loop.
    let finalize = |pending: &mut Vec<&htsjdk_bam::record::BamRecord>,
                    reject: bool,
                    out: &mut Vec<htsjdk_bam::record::BamRecord>,
                    random: &mut JavaRandom| {
        if !reject {
            let mut batch: Vec<htsjdk_bam::record::BamRecord> = if pending.len() <= max_coverage {
                pending.iter().map(|r| (*r).clone()).collect()
            } else {
                let mut reservoir: ReservoirDownsampler<'_, htsjdk_bam::record::BamRecord> =
                    ReservoirDownsampler::new(max_coverage);
                for read in pending
                    .iter()
                    .filter(|r| r.mapping_quality > SUSPICIOUS_MAPPING_QUALITY)
                {
                    reservoir.submit(read, &read.read_name, &mut SlotSource::Random(random));
                }
                reservoir.signal_end_of_input();
                reservoir
                    .consume_finalized_items()
                    .into_iter()
                    .cloned()
                    .collect()
            };
            // Stable, as `Collections.sort` is, so reads with one start keep their order.
            batch.sort_by_key(|r| r.alignment_start);
            out.extend(batch);
        }
        pending.clear();
    };
    for read in reads {
        // A read with no position is finalized at once and never joins a stride.
        if read.reference_index < 0 || read.alignment_start <= 0 {
            out.push(read.clone());
            continue;
        }
        // `handlePositionalChange`: a new contig, or a start a full stride past the stride's first
        // read, closes the stride.
        match first {
            None => first = Some((read.reference_index, read.alignment_start)),
            Some((contig, start)) => {
                if read.reference_index != contig || read.alignment_start >= start + stride {
                    finalize(&mut pending, reject, &mut out, random);
                    first = Some((read.reference_index, read.alignment_start));
                    reject = false;
                    suspicious = 0;
                }
            }
        }
        // Once rejected, the rest of the stride is dropped without being counted.
        if reject {
            continue;
        }
        if read.mapping_quality <= SUSPICIOUS_MAPPING_QUALITY {
            suspicious += 1;
        }
        // The read that reaches the count is still added, and dropped with the rest at finalize.
        reject |= suspicious >= max_suspicious;
        pending.push(read);
    }
    finalize(&mut pending, reject, &mut out, random);
    out
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
///    scan: the loop ends after `groupSizes.length` refusals in a row, one full lap in which no
///    stack is above `minElementsPerStack`. Sizes only fall, so a stack that refuses once refuses
///    for the rest of the walk;
///  * **each stack is reduced by a full `nextPermutation`**, which shuffles all its items and keeps
///    the head. That costs `size - 1` draws from the shared `Well19937c` per stack reduced, in
///    stack order, so the stacks' *order* changes which items every later stack keeps;
///  * **a stack that keeps everything takes no draw at all**, because `downsampleOneGroup` returns
///    before sampling when `numItemsToKeep >= group.size()`.
///
/// With `min_elements_per_stack` of zero a stack can be planned down to zero items, and then
/// `nextPermutation(n, 0)` throws rather than returning nothing. That refusal is reproduced, not
/// smoothed over: see [`crate::permutation::PermutationError::NotStrictlyPositive`].
///
/// Each stack is a `Vec<T>`, the reference's `T extends List<E>`.
pub struct LevelingDownsampler<T> {
    /// `targetSize`, a `long` upstream.
    target_size: i64,
    /// `minElementsPerStack`. Upstream it is an `int` checked `>= 0`; a `usize` cannot be negative.
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

    /// `submit(T item)`: one stack, appended. Nothing is removed until end of input.
    pub fn submit(&mut self, group: Vec<T>) {
        self.groups.push(group);
    }

    /// `size()`: the sum over the stacks, not the number of stacks.
    pub fn size(&self) -> usize {
        self.groups.iter().map(Vec::len).sum()
    }

    /// `getNumberOfDiscardedItems()`.
    pub fn discarded(&self) -> usize {
        self.discarded
    }

    /// `hasFinalizedItems()`.
    pub fn has_finalized_items(&self) -> bool {
        self.groups_are_finalized && !self.groups.is_empty()
    }

    /// `hasPendingItems()`.
    pub fn has_pending_items(&self) -> bool {
        !self.groups_are_finalized && !self.groups.is_empty()
    }

    /// `signalEndOfInput()`, which is where the levelling actually happens: nothing is removed
    /// until every stack has been submitted, because the plan needs all the sizes.
    ///
    /// On a refusal the reference throws before it sets `groupsAreFinalized`; this sets the flag
    /// either way and returns the error. Every caller stops on that error (`read_states` maps it,
    /// the suite moves on to the next case), so the flag is never read afterwards.
    pub fn signal_end_of_input(
        &mut self,
        rng: &mut Well19937c,
    ) -> Result<(), crate::permutation::PermutationError> {
        let result = self.level_groups(rng);
        self.groups_are_finalized = true;
        result
    }

    /// `consumeFinalizedItems()`: the stacks, then `clearItems()`. Before end of input, an empty
    /// list and no change, as upstream.
    pub fn consume_finalized_items(&mut self) -> Vec<Vec<T>> {
        if !self.has_finalized_items() {
            return Vec::new();
        }
        let groups = std::mem::take(&mut self.groups);
        self.groups_are_finalized = false;
        groups
    }

    /// `levelGroups()`: plan the per-stack sizes, then reduce each stack to its planned size.
    fn level_groups(
        &mut self,
        rng: &mut Well19937c,
    ) -> Result<(), crate::permutation::PermutationError> {
        let mut group_sizes: Vec<usize> = self.groups.iter().map(Vec::len).collect();
        let total_size: usize = group_sizes.iter().sum();

        if total_size as i64 <= self.target_size {
            return Ok(());
        }

        // The plan. `current` walks the stacks round-robin, taking one item from any stack still
        // above the minimum, until enough are gone or a full lap gives nothing.
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

        // Reduced in submission order, which is also the order the draws are taken in. `zip` pairs
        // each stack with its planned size; `iter_mut` lends each stack writably in turn.
        for (group, keep) in self.groups.iter_mut().zip(group_sizes) {
            self.discarded += downsample_one_group(group, keep, rng)?;
        }
        Ok(())
    }
}

/// `downsampleOneGroup`, returning what it counts through `incrementNumberOfDiscardedItems`.
///
/// The kept items stay in the stack's own order: the permutation is used as set membership, not as
/// an ordering, so shuffled indices produce an unshuffled result. The reference marks the kept
/// indices in a `BitSet` and copies the marked items back (or removes the unmarked ones through an
/// iterator, for a `LinkedList`); both keep the stack's order, and so does `retain`.
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
    // `retain` visits the items in order and keeps those its closure answers `true` for; the
    // closure ignores the item and reads the next flag, so item `i` is kept when `wanted[i]` is.
    let mut keeper = wanted.iter();
    let before = group.len();
    group.retain(|_| *keeper.next().unwrap_or(&false));
    Ok(before - group.len())
}
