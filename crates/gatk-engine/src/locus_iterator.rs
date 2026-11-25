//! Ported from `org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState`
//! (GATK 4.7.0.0; `LocusIteratorByState.java` is unchanged since 4.6.2.0).
//!
//! One pileup per covered locus, assembled from the per-read state machines that
//! [`crate::read_states::ReadStateManager`] keeps in step. This is what a `LocusWalker` iterates,
//! so which elements land in which pileup is every locus-based tool's output.
//!
//! Four decisions sit in the twenty lines that build each pileup:
//!
//!  * **the two exclusions are not symmetric.** A read whose current operator is `N` is skipped
//!    *before* the adaptor test; a read whose operator is `D` is skipped *inside* it. So a read
//!    sitting in an adaptor with a deletion is excluded once, by the adaptor, and a read with an
//!    `N` is excluded whatever the adaptor says. All three tests are a `continue` with no side
//!    effect, so their order decides which test excludes a read, not whether it is excluded; it
//!    is kept as the reference's so that the two read side by side;
//!  * **the adaptor test is per base, not per read.** `isBaseInsideAdaptor` compares this locus
//!    against the adaptor boundary, so the same read contributes to some loci and not others;
//!  * **the pileup is monolithic and in sample order.** The elements of every sample are
//!    concatenated in the declaration order of the samples, and the reference explains it kept it
//!    that way for the HaplotypeCaller's benefit;
//!  * **an empty locus produces no context at all.** The iterator loops until it has at least one
//!    element, so a position covered only by `N`s or only by adaptor bases is skipped in silence.
//!    Emitting empty loci is a different class's job, and not this one's.
//!
//! # Reading this file without reading Rust
//!
//! - `'a` and `'h` are lifetimes: `Contexts<'a, 'h>` borrows the reads for `'a` and the header
//!   for `'h`, and the compiler proves neither is freed while the iterator is alive. A Java reader
//!   can read them as "this holds references to the caller's reads and header".
//! - `Result<Option<T>, E>` from `step` is three outcomes: a context, the end, or an error.
//!   `.transpose()` turns it into the `Option<Result<T, E>>` that Rust's `Iterator::next` returns,
//!   where `None` is Java's `hasNext() == false`.
//! - `let Some(x) = y else { return false; };` unwraps an `Option` or takes the early exit.
//! - `&dyn Fn(&BamRecord) -> bool` is a borrowed closure, Java's `Predicate<GATKRead>`.

use crate::alignment_state::AlignmentStateMachine;
use crate::pileup::PileupElement;
use crate::read_pileup::ReadPileup;
use crate::read_states::{ReadStateError, ReadStateManager};
use crate::read_utils;
use htsjdk_bam::cigar::Op;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// `ReadUtils.DEFAULT_ADAPTOR_SIZE`.
pub const DEFAULT_ADAPTOR_SIZE: i32 = 100;

/// `ReadUtils.isBaseInsideAdaptor`.
///
/// Two guards before the comparison, and both matter: a read whose adaptor boundary cannot be
/// computed is never inside one, and so is a read whose fragment is longer than 100 bases, which
/// is a *length* test rather than a position test. The length is the signed `TLEN`, so a negative
/// one never trips the second guard.
///
/// The comparison itself:
///
/// ```java
/// return read.isReverseStrand() ? basePos <= adaptorBoundary : basePos >= adaptorBoundary;
/// ```
pub fn is_base_inside_adaptor(read: &BamRecord, position: i32) -> bool {
    // `CANNOT_COMPUTE_ADAPTOR_BOUNDARY` is `None` here.
    let Some(boundary) = read_utils::adaptor_boundary(read) else {
        return false;
    };
    if read.inferred_insert_size > DEFAULT_ADAPTOR_SIZE {
        return false;
    }
    if crate::read::is_reverse_strand(read) {
        position <= boundary
    } else {
        position >= boundary
    }
}

/// `AlignmentContext`: a locus and its pileup.
pub struct AlignmentContext<'a> {
    pub contig: String,
    /// The 1-based position: the reference's `getLocation().getStart()`.
    pub position: i32,
    pub pileup: ReadPileup<'a>,
}

/// What the iterator was configured with: the two booleans of the six-argument
/// `LocusIteratorByState` constructor.
///
/// `Default` gives `LocusWalker`'s defaults. A caller porting a different construction (the
/// five-argument constructor passes `includeReadsWithNsAtLoci` as false and takes
/// `includeReadsWithDeletionAtLoci` from its caller) has to spell both out.
#[derive(Debug, Clone, Copy)]
pub struct LocusIteratorOptions {
    /// `includeReadsWithDeletionAtLoci`. `LocusWalker.includeDeletions()` defaults to true.
    pub include_deletions: bool,
    /// `includeReadsWithNsAtLoci`. `LocusWalker.includeNs()` defaults to false.
    pub include_ns: bool,
}

impl Default for LocusIteratorOptions {
    fn default() -> Self {
        LocusIteratorOptions {
            include_deletions: true,
            include_ns: false,
        }
    }
}

/// `LocusIteratorByState.lazyLoadNextAlignmentContext`, run to exhaustion.
///
/// Returns every context the iterator would yield, in order. The reference is lazy because the
/// HaplotypeCaller cares; the order and the contents are the same either way, and collecting them
/// is what makes the suite comparable line by line.
///
/// `samples` is not read: the samples are already in `states`, which was built with them.
pub fn contexts<'a>(
    reads: &'a [BamRecord],
    samples: Vec<Option<String>>,
    header: &SamHeader,
    options: LocusIteratorOptions,
    states: ReadStateManager<'a>,
) -> Result<Vec<AlignmentContext<'a>>, ReadStateError> {
    contexts_filtered(reads, samples, header, options, states, &|_| true)
}

/// The same, with the walker's read filter applied first.
///
/// The filter runs before the iterator sees anything, so a filtered read is *absent* from the
/// pileup rather than present and ignored: it changes the depth, not only the reported set.
pub fn contexts_filtered<'a>(
    reads: &'a [BamRecord],
    samples: Vec<Option<String>>,
    header: &SamHeader,
    options: LocusIteratorOptions,
    states: ReadStateManager<'a>,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<Vec<AlignmentContext<'a>>, ReadStateError> {
    // Collecting an iterator of `Result`s into one `Result<Vec<_>, _>` stops at the first `Err`
    // and returns it, or returns every context.
    Contexts::new(reads, samples, header, options, states, filter).collect()
}

/// `LocusIteratorByState` as the iterator it is upstream: one context at a time, so a traversal
/// holds one pileup rather than every pileup of a contig. The first error ends the iteration.
pub struct Contexts<'a, 'h> {
    /// The reads not yet handed to the state manager, in input order: the reference's
    /// `samIterator`, after the walker's filter.
    pending: std::collections::VecDeque<&'a BamRecord>,
    states: ReadStateManager<'a>,
    header: &'h SamHeader,
    options: LocusIteratorOptions,
    // The two streams the downsamplers draw from, and they are ONE EACH for the whole traversal.
    // `Utils.getRandomGenerator()` and `Utils.getRandomDataGenerator()` are static, seeded once at
    // class initialisation, and never reset between loci or samples: a generator built per locus
    // would restart the sequence and produce a different pileup at every depth.
    random: crate::java_random::JavaRandom,
    leveling: crate::well19937c::Well19937c,
    /// Set by the first error, after which the iterator yields nothing more.
    failed: bool,
}

impl<'a, 'h> Contexts<'a, 'h> {
    /// The constructor: filter the reads, keep them in order, and seed the two streams.
    pub fn new(
        reads: &'a [BamRecord],
        samples: Vec<Option<String>>,
        header: &'h SamHeader,
        options: LocusIteratorOptions,
        states: ReadStateManager<'a>,
        filter: &dyn Fn(&BamRecord) -> bool,
    ) -> Contexts<'a, 'h> {
        // `let _ =` discards the argument on purpose: the samples are already in `states`.
        let _ = samples;
        Contexts {
            pending: reads.iter().filter(|read| filter(read)).collect(),
            states,
            header,
            options,
            random: crate::java_random::JavaRandom::gatk(),
            leveling: crate::well19937c::Well19937c::gatk(),
            failed: false,
        }
    }

    /// `lazyLoadNextAlignmentContext`: advance until a locus has at least one element, or the
    /// reads run out (`Ok(None)`).
    fn step(&mut self) -> Result<Option<AlignmentContext<'a>>, ReadStateError> {
        // Borrow the fields once, so that `states` can be written while `self.pending` and the
        // generators are borrowed separately below.
        let (states, options, header) = (&mut self.states, self.options, self.header);
        // `readStates.hasNext()`: states in the system, or reads still to come.
        while !states.is_empty() || !self.pending.is_empty() {
            states.collect_pending_reads(
                &mut self.pending,
                header,
                &mut self.random,
                &mut self.leveling,
            )?;

            // `getLocation()`, which is null when nothing is in the system. The per-state loop
            // below is then empty, so the null never reaches a dereference. The location is the
            // first state's (contig index, 1-based position), as `readStates.getFirst()` is.
            let location = states
                .first()
                .map(|state| (state.read.reference_index, state.machine.genome_position()));

            let mut elements: Vec<PileupElement<'a>> = Vec::new();
            for per_sample in &states.by_sample {
                for state in &per_sample.states {
                    let op = state.machine.cigar_operator();
                    // First exclusion: an N is dropped before the adaptor is ever consulted.
                    if !options.include_ns && op == Some(Op::N) {
                        continue;
                    }
                    // A state exists, so `location` is `Some`; the `0` is never used.
                    let position = location.map(|(_, p)| p).unwrap_or(0);
                    if is_base_inside_adaptor(state.read, position) {
                        continue;
                    }
                    // Second exclusion: a D is dropped only for a read the adaptor test let
                    // through.
                    if !options.include_deletions && op == Some(Op::D) {
                        continue;
                    }
                    // `state.makePileupElement()`, which throws on an edge; `from_state` answers
                    // `None` there and the state is skipped instead.
                    if let Some(element) = PileupElement::from_state(state.read, &state.machine) {
                        elements.push(element);
                    }
                }
            }

            // Critical, and the reference says so: the states advance only after the current
            // offsets and location have been read.
            states.update_read_states()?;

            if !elements.is_empty() {
                // `.expect(..)` unwraps the `Option`, panicking with that text if it were `None`.
                let (reference_index, position) = location.expect("elements imply a location");
                let contig = header
                    .sequences
                    .get(reference_index as usize)
                    .map(|s| s.name.clone())
                    .unwrap_or_default();
                return Ok(Some(AlignmentContext {
                    contig: contig.clone(),
                    position,
                    pileup: ReadPileup::new(&contig, position, elements),
                }));
            }
        }
        Ok(None)
    }
}

/// The reference's `hasNext`/`next` pair, as one `next` that answers `None` at the end.
impl<'a> Iterator for Contexts<'a, '_> {
    type Item = Result<AlignmentContext<'a>, ReadStateError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        let step = self.step();
        self.failed = step.is_err();
        step.transpose()
    }
}

/// `AlignmentStateMachine.getLocation`, for a caller that wants the position without the context.
///
/// The reference answers a one-base `SimpleInterval`; this answers its start.
pub fn location_of(machine: &AlignmentStateMachine) -> i32 {
    machine.genome_position()
}
