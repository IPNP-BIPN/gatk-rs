//! Ported from `org.broadinstitute.hellbender.utils.locusiterator.ReadStateManager`,
//! `.PerSampleReadStateManager`, `.SamplePartitioner` and `.LIBSDownsamplingInfo`
//! (GATK 4.6.2.0).
//!
//! Between the per-read state machines and the pileup sits a bookkeeping layer that decides which
//! reads enter, under which sample, and when they leave. Four of its behaviours change what a
//! locus contains:
//!
//!  * **the sample map is a `LinkedHashMap`**, and the reference says so in capitals: the
//!    iteration order is the order the samples were given at construction, not sorted and not the
//!    header's. Every per-sample pileup a walker sees inherits that order;
//!  * **a read whose first step returns null is dropped, silently.** `addReadsToSample` builds a
//!    state machine per read and keeps it only if `stepForwardOnGenome()` is non-null, with a
//!    `todo` upstream saying this should have been an assertion. A read that is all insertions and
//!    soft clips therefore never reaches any pileup;
//!  * **the boundary is a genome position, not a read start.** `collectPendingReads` takes the
//!    left-most position among the states already in the system, which for a read part-way through
//!    a deletion is *not* where that read started, and admits only reads whose start equals it
//!    exactly on the same contig;
//!  * **a read with no read group has the null sample**, and a read whose sample was not declared
//!    at construction is a hard error rather than a new bucket.
//!
//! # What this refuses
//!
//! Downsampling is not ported. `ReservoirDownsampler` and `LevelingDownsampler` draw from
//! `Utils.getRandomGenerator`, so reproducing them means reproducing Java's `Random` and the exact
//! sequence of draws; until that exists, [`DownsamplingInfo::performing`] is refused rather than
//! approximated. The default is off (`--max-depth-per-sample 0` gives `NO_DOWNSAMPLING`), which is
//! the path every suite here exercises.

use crate::alignment_state::{AlignmentStateMachine, MalformedRead};
use crate::read_pileup::sample_name;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// `LIBSDownsamplingInfo`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DownsamplingInfo {
    pub performing: bool,
    pub to_coverage: i32,
}

impl DownsamplingInfo {
    /// `LocusIteratorByState.NO_DOWNSAMPLING`, which is `(false, -1)` and not `(false, 0)`.
    pub const NONE: DownsamplingInfo = DownsamplingInfo {
        performing: false,
        to_coverage: -1,
    };
}

/// What the bookkeeping layer can refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadStateError {
    /// A read whose sample was not among those declared at construction.
    UndeclaredSample(Option<String>),
    /// A cigar `AlignmentStateMachine` refuses.
    Malformed(MalformedRead),
    /// Downsampling, which this port does not reproduce. See the module doc.
    DownsamplingUnsupported,
}

/// One read, and where its state machine currently sits.
///
/// The machine borrows the read, so the reads outlive the manager; that is the reference's shape
/// too, where the `GATKRead` is held by the `AlignmentStateMachine`.
pub struct ReadState<'a> {
    pub read: &'a BamRecord,
    pub machine: AlignmentStateMachine<'a>,
}

/// `PerSampleReadStateManager` without the leveling downsampler.
#[derive(Default)]
pub struct PerSampleReadStateManager<'a> {
    /// `readStatesByAlignmentStart`, in insertion order.
    pub states: Vec<ReadState<'a>>,
}

impl<'a> PerSampleReadStateManager<'a> {
    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    pub fn size(&self) -> usize {
        self.states.len()
    }

    /// `addStatesAtNextAlignmentStart`: append, and report how many were added.
    ///
    /// With downsampling off this is exactly the count appended. With it on the reference subtracts
    /// what the leveling downsampler discarded, which is the branch this port refuses.
    pub fn add_states(&mut self, states: Vec<ReadState<'a>>) -> usize {
        let added = states.len();
        self.states.extend(states);
        added
    }

    /// `updateReadStates`: step every machine one reference base, and drop the ones that ran out.
    ///
    /// Returns how many were removed. A malformed cigar surfaces here rather than at construction,
    /// because the first step happened before the state was ever added.
    pub fn update_read_states(&mut self) -> Result<usize, ReadStateError> {
        let mut removed = 0;
        let mut kept = Vec::with_capacity(self.states.len());
        for mut state in std::mem::take(&mut self.states) {
            match state.machine.step_forward_on_genome() {
                Ok(None) => removed += 1,
                Ok(Some(_)) => kept.push(state),
                Err(error) => return Err(ReadStateError::Malformed(error)),
            }
        }
        self.states = kept;
        Ok(removed)
    }

    pub fn first(&self) -> Option<&ReadState<'a>> {
        self.states.first()
    }
}

/// `ReadStateManager`: one `PerSampleReadStateManager` per declared sample, in that order.
pub struct ReadStateManager<'a> {
    /// The declared samples, in the order given. `None` is the bucket for reads with no read group,
    /// which the reference keys under a null sample name.
    pub samples: Vec<Option<String>>,
    /// Parallel to `samples`, because the reference's `LinkedHashMap` iterates in that order and a
    /// hash map here would not.
    pub by_sample: Vec<PerSampleReadStateManager<'a>>,
    total_read_states: usize,
}

impl<'a> ReadStateManager<'a> {
    pub fn new(
        samples: Vec<Option<String>>,
        info: DownsamplingInfo,
    ) -> Result<Self, ReadStateError> {
        if info.performing {
            return Err(ReadStateError::DownsamplingUnsupported);
        }
        let by_sample = samples
            .iter()
            .map(|_| PerSampleReadStateManager::default())
            .collect();
        Ok(ReadStateManager {
            samples,
            by_sample,
            total_read_states: 0,
        })
    }

    pub fn size(&self) -> usize {
        self.total_read_states
    }

    pub fn is_empty(&self) -> bool {
        self.total_read_states == 0
    }

    /// `getFirst`: the first state of the first non-empty sample, in sample order.
    pub fn first(&self) -> Option<&ReadState<'a>> {
        self.by_sample.iter().find_map(|manager| manager.first())
    }

    /// `updateReadStates`, across every sample.
    pub fn update_read_states(&mut self) -> Result<(), ReadStateError> {
        for manager in &mut self.by_sample {
            self.total_read_states -= manager.update_read_states()?;
        }
        Ok(())
    }

    /// `collectPendingReads`: admit every read starting exactly at the current left-most position.
    ///
    /// `pending` is the peekable source, consumed from the front. Returns how many reads were
    /// admitted, so a caller can tell a no-op from a full one.
    ///
    /// The boundary is `getFirst().getGenomePosition()` when the system is non-empty, which is
    /// where the left-most read *currently is* rather than where it started. When the system is
    /// empty it is the next read's own start.
    pub fn collect_pending_reads(
        &mut self,
        pending: &mut std::collections::VecDeque<&'a BamRecord>,
        header: &SamHeader,
    ) -> Result<usize, ReadStateError> {
        let Some(next) = pending.front() else {
            return Ok(0);
        };
        let (contig_index, alignment_start) = if self.is_empty() {
            (next.reference_index, next.alignment_start)
        } else {
            let first = self.first().expect("non-empty");
            (first.read.reference_index, first.machine.genome_position())
        };

        // `SamplePartitioner`, with a pass-through downsampler: the reads are bucketed by sample
        // in submission order and handed back in that order.
        let mut buckets: Vec<Vec<&'a BamRecord>> =
            self.samples.iter().map(|_| Vec::new()).collect();
        while let Some(read) = pending.front() {
            if read.alignment_start != alignment_start || read.reference_index != contig_index {
                break;
            }
            let read = pending.pop_front().expect("peeked");
            let sample = sample_of(read, header);
            let index = self
                .samples
                .iter()
                .position(|s| *s == sample)
                .ok_or_else(|| ReadStateError::UndeclaredSample(sample.clone()))?;
            buckets[index].push(read);
        }

        let mut admitted = 0;
        for (index, reads) in buckets.into_iter().enumerate() {
            if reads.is_empty() {
                continue;
            }
            let mut states = Vec::with_capacity(reads.len());
            for read in reads {
                let mut machine = AlignmentStateMachine::new(read);
                match machine.step_forward_on_genome() {
                    // Dropped in silence: a read that is all insertions and soft clips never
                    // reaches any pileup, and upstream calls that a todo rather than an error.
                    Ok(None) => continue,
                    Ok(Some(_)) => states.push(ReadState { read, machine }),
                    Err(error) => return Err(ReadStateError::Malformed(error)),
                }
            }
            admitted += self.by_sample[index].add_states(states);
        }
        self.total_read_states += admitted;
        Ok(admitted)
    }
}

/// `SamplePartitioner.submitRead`'s sample lookup.
///
/// Two nulls collapse into one bucket here, as upstream: a read with no read group at all, and a
/// read whose read group declares no `SM`. The reference reaches the second through
/// `ReadUtils.getSampleName` only when `getReadGroup()` is non-null, and both end at the same key.
pub fn sample_of(read: &BamRecord, header: &SamHeader) -> Option<String> {
    let has_group = read.tags.get(htsjdk_bam::tag::Tag::new(b"RG")).is_some();
    if !has_group {
        return None;
    }
    sample_name(read, header)
}
