//! The iterators between a shard of reads and an activity profile, ported from
//! `utils.iterators.ShardedIntervalIterator`, `utils.iterators.IntervalLocusIterator`,
//! `utils.iterators.ReadCachingIterator` and
//! `utils.locusiterator.IntervalAlignmentContextIterator` (GATK 4.6.2.0).
//!
//! # Why an empty pileup is a value and not an absence
//!
//! `LocusIteratorByState` yields nothing at all for a locus no read covers: a position with zero
//! depth is silently absent from its output. The assembly-region traversal cannot work that way,
//! because the activity profile has to see the gap in order to close a region over it. So
//! `AssemblyRegionIterator` wraps the locus iterator in `IntervalAlignmentContextIterator`, whose
//! whole job is to walk **every base of every interval** and manufacture an empty
//! `AlignmentContext` wherever the wrapped iterator has none. The comment upstream is unusually
//! direct about the stakes: "This is critical for reproducing GATK 3.x behavior!"
//!
//! A port that skipped the wrapper would produce regions that end at the last covered base instead
//! of where the activity actually dies, which moves every boundary downstream of a coverage gap.
//!
//! # The shard arithmetic is absolute, not relative
//!
//! `ShardedIntervalIterator` does not cut an interval into `shardSize` chunks starting at the
//! interval's own start. It computes shard *indices* over the interval's length with
//! `IntervalUtils.shardIndex`, then converts each index back to coordinates and adds the interval's
//! start. The consequence is visible with any shard size above 1: the first shard of `chr1:10-20`
//! at size 3 is `chr1:10-12`, and the last one is truncated against the interval's end rather than
//! running past it. At size 1, which is the only size `IntervalLocusIterator` ever uses, every
//! shard is one base and the arithmetic is invisible. It is ported anyway because the class is
//! shared.

use crate::interval::SimpleInterval;
use crate::variant_source::Located;
use htsjdk_bam::header::SamHeader;

/// `IntervalUtils.shardIndex`: which shard a 1-based offset falls in, numbered from zero.
pub fn shard_index(one_based_offset: i32, shard_size: i32) -> i32 {
    (one_based_offset - 1) / shard_size
}

/// `IntervalUtils.beginOfShard`, 1-based.
pub fn begin_of_shard(shard_index: i32, shard_size: i32) -> i32 {
    shard_index * shard_size + 1
}

/// `IntervalUtils.endOfShard`, 1-based and inclusive.
pub fn end_of_shard(shard_index: i32, shard_size: i32) -> i32 {
    begin_of_shard(shard_index + 1, shard_size) - 1
}

/// `IntervalUtils.compareContigs`.
///
/// `None` is the `IllegalArgumentException` the reference throws when either contig is absent from
/// the dictionary: the comparison is refused rather than answered with an arbitrary order.
pub fn compare_contigs(left: &str, right: &str, header: &SamHeader) -> Option<std::cmp::Ordering> {
    let index = |contig: &str| {
        header
            .sequences
            .iter()
            .position(|sequence| sequence.name == contig)
    };
    let left_index = index(left)?;
    let right_index = index(right)?;
    Some(left_index.cmp(&right_index))
}

/// `IntervalUtils.compareLocatables`: contig, then start, then **end**.
///
/// The end is the third key and it matters here: two loci at the same start but different ends are
/// ordered, and the interval iterator relies on that ordering to decide whether to advance.
pub fn compare_locatables(
    left: &dyn Located,
    right: &dyn Located,
    header: &SamHeader,
) -> Option<std::cmp::Ordering> {
    let contigs = compare_contigs(left.contig(), right.contig(), header)?;
    if contigs != std::cmp::Ordering::Equal {
        return Some(contigs);
    }
    Some(
        left.start()
            .cmp(&right.start())
            .then(left.stop().cmp(&right.stop())),
    )
}

/// `IntervalUtils.isAfter`: the first starts after the second **ends**.
pub fn is_after(left: &dyn Located, right: &dyn Located, header: &SamHeader) -> Option<bool> {
    let contigs = compare_contigs(left.contig(), right.contig(), header)?;
    Some(
        contigs == std::cmp::Ordering::Greater
            || (contigs == std::cmp::Ordering::Equal && left.start() > right.stop()),
    )
}

/// `IntervalUtils.isBefore`: the first **ends** before the second starts.
pub fn is_before(left: &dyn Located, right: &dyn Located, header: &SamHeader) -> Option<bool> {
    let contigs = compare_contigs(left.contig(), right.contig(), header)?;
    Some(
        contigs == std::cmp::Ordering::Less
            || (contigs == std::cmp::Ordering::Equal && left.stop() < right.start()),
    )
}

/// `ShardedIntervalIterator`, run to exhaustion.
///
/// `None` is the `Utils.validate` on a shard size of zero or less.
pub fn sharded_intervals(
    intervals: &[SimpleInterval],
    shard_size: i32,
) -> Option<Vec<SimpleInterval>> {
    if shard_size <= 0 {
        return None;
    }
    let mut out = Vec::new();
    for interval in intervals {
        let last = shard_index(interval.size(), shard_size);
        for index in shard_index(1, shard_size)..=last {
            let start = interval.start + begin_of_shard(index, shard_size) - 1;
            let end = (interval.start + end_of_shard(index, shard_size) - 1).min(interval.end);
            out.push(
                SimpleInterval::new(&interval.contig, start, end)
                    .expect("a shard of a valid interval is a valid interval"),
            );
        }
    }
    Some(out)
}

/// `IntervalLocusIterator`, run to exhaustion: one single-base interval per base of the input.
///
/// It is `ShardedIntervalIterator` at shard size 1, and nothing else. Worth writing out rather than
/// aliasing, because the constant is the entire behaviour.
pub fn interval_loci(intervals: &[SimpleInterval]) -> Vec<SimpleInterval> {
    sharded_intervals(intervals, 1).expect("a shard size of 1 is valid")
}

/// What `IntervalAlignmentContextIterator` emitted at one locus.
///
/// A covered locus carries the index of the wrapped iterator's context; an uncovered one carries
/// nothing, and upstream is a freshly constructed `AlignmentContext` with an empty `ReadPileup`.
/// Keeping it as an index rather than a clone is what lets the caller keep the pileup's borrows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EmittedLocus {
    pub interval: SimpleInterval,
    pub context: Option<usize>,
}

/// `IntervalAlignmentContextIterator`, run to exhaustion.
///
/// `contexts` is what `LocusIteratorByState` produced, in order; `intervals` is the interval list
/// the traversal was given, **before** being cut into loci, exactly as `AssemblyRegionIterator`
/// hands it over.
///
/// The state machine is the reference's, including the part that reads like an off-by-one and is
/// not: in the *uncovered* branch it advances the alignment context only when the new locus is
/// strictly after it. Advancing unconditionally would skip a context whenever two uncovered loci
/// precede a covered one.
///
/// The exhaustion rule is also the reference's, and it is what makes the loop terminate: once the
/// wrapped iterator runs dry, `advanceAlignmentContext` manufactures an empty context **at the
/// current locus**, so the comparison is an equality and the advance loop stops instead of
/// spinning.
pub fn interval_alignment_contexts<T: Located>(
    contexts: &[T],
    intervals: &[SimpleInterval],
    header: &SamHeader,
) -> Vec<EmittedLocus> {
    let mut out = Vec::new();
    let _ = for_each_interval_context(
        contexts.iter(),
        intervals,
        header,
        |locus, _: Option<&&T>| -> Result<(), ()> {
            out.push(locus);
            Ok(())
        },
    );
    out
}

/// `IntervalAlignmentContextIterator` as the iterator it is upstream: each locus is handed to `f`
/// with the context it carries, and neither the loci nor the contexts are held past their turn.
/// `EmittedLocus::context` counts the contexts drawn from `contexts`. The first error `f` answers
/// stops the walk.
pub fn for_each_interval_context<C: Located, E>(
    contexts: impl Iterator<Item = C>,
    intervals: &[SimpleInterval],
    header: &SamHeader,
    mut f: impl FnMut(EmittedLocus, Option<&C>) -> Result<(), E>,
) -> Result<(), E> {
    // `IntervalLocusIterator`: `interval_loci`, one base at a time.
    let mut loci = intervals.iter().flat_map(|interval| {
        (interval.start..=interval.end).map(move |position| {
            SimpleInterval::new(&interval.contig, position, position)
                .expect("a shard of a valid interval is a valid interval")
        })
    });
    let mut contexts = contexts.enumerate();

    // Position of the current context for comparison purposes: an exhausted iterator's context sits
    // at the current locus, so it compares equal to it.
    let context_position =
        |context: &Option<(usize, C)>, locus: &SimpleInterval| -> SimpleInterval {
            match context {
                Some((_, context)) => SimpleInterval {
                    contig: context.contig().to_string(),
                    start: context.start(),
                    end: context.stop(),
                },
                None => locus.clone(),
            }
        };

    // `currentInterval` and `currentAlignmentContext`. The context is `Some` when it is one of the
    // wrapped iterator's, and `None` when it is a manufactured empty one at the locus.
    let mut current_locus = loci.next();
    let mut current_context = contexts.next();

    // `advanceAlignmentContextToCurrentInterval`.
    let catch_up = |current_locus: &Option<SimpleInterval>,
                    current_context: &mut Option<(usize, C)>,
                    contexts: &mut dyn Iterator<Item = (usize, C)>| {
        let Some(locus) = current_locus else {
            *current_context = None;
            return;
        };
        loop {
            let position = context_position(current_context, locus);
            match compare_locatables(locus, &position, header) {
                Some(std::cmp::Ordering::Greater) => {
                    *current_context = contexts.next();
                }
                _ => break,
            }
        }
    };
    catch_up(&current_locus, &mut current_context, &mut contexts);

    while let Some(locus) = current_locus.take() {
        let position = context_position(&current_context, &locus);
        let overlaps = locus.overlaps(&position.contig, position.start, position.end);

        if overlaps {
            let context = current_context.as_ref();
            f(
                EmittedLocus {
                    interval: locus,
                    context: context.map(|(index, _)| *index),
                },
                context.map(|(_, context)| context),
            )?;
            current_locus = loci.next();
            catch_up(&current_locus, &mut current_context, &mut contexts);
        } else {
            f(
                EmittedLocus {
                    interval: locus,
                    context: None,
                },
                None,
            )?;
            current_locus = loci.next();
            if let Some(next) = &current_locus {
                let position = context_position(&current_context, next);
                if compare_locatables(next, &position, header) == Some(std::cmp::Ordering::Greater)
                {
                    catch_up(&current_locus, &mut current_context, &mut contexts);
                }
            }
        }
    }
    Ok(())
}

/// `ReadCachingIterator`: hand each item through and keep a copy until the client takes them.
///
/// A wrapper this thin is worth porting only because of what it guarantees to
/// `AssemblyRegionIterator`: the cache comes back in the order the wrapped iterator produced it,
/// which for a coordinate-sorted shard is coordinate order, and the region filler depends on that
/// to stop at the first read past the region rather than scanning the rest.
pub struct ReadCachingIterator<'a, T> {
    source: std::vec::IntoIter<&'a T>,
    cache: Vec<&'a T>,
}

impl<'a, T> ReadCachingIterator<'a, T> {
    pub fn new(items: Vec<&'a T>) -> ReadCachingIterator<'a, T> {
        ReadCachingIterator {
            source: items.into_iter(),
            cache: Vec::new(),
        }
    }

    /// `consumeCachedReads`: everything cached so far, and the cache is emptied by the call.
    pub fn consume_cached_reads(&mut self) -> Vec<&'a T> {
        std::mem::take(&mut self.cache)
    }
}

impl<'a, T> Iterator for ReadCachingIterator<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<&'a T> {
        let item = self.source.next()?;
        self.cache.push(item);
        Some(item)
    }
}

impl Located for SimpleInterval {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn stop(&self) -> i32 {
        self.end
    }
}
