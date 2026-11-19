//! `AssemblyRegionWalker`, ported from `org.broadinstitute.hellbender.engine.AssemblyRegionWalker`
//! (GATK 4.6.2.0).
//!
//! The last of the five walker archetypes, and the one `HaplotypeCaller` and `Mutect2` are built
//! on. Everything under it is already ported: [`crate::activity_profile`] decides where a region
//! starts and stops, [`crate::assembly_region`] is what a region is, [`crate::locus_shards`] turns
//! a coverage gap into a locus, and [`crate::assembly_region_iterator`] is the traversal. This is
//! the layer above them, and it is thin. What it adds is worth naming because three of the four
//! things are invisible from any tool's own output.
//!
//! # One shard per contig, not one per interval
//!
//! `makeReadShards` groups the traversal intervals by contig and builds **one**
//! `MultiIntervalLocalReadShard` per group. Two `-L` arguments on the same contig therefore share
//! a shard, share an activity profile, and can end up in the same assembly region if the padding
//! joins them; two on different contigs cannot, whatever their distance in the file.
//!
//! # `apply` gets the padded span, not the active one
//!
//! The `ReferenceContext` and `FeatureContext` handed to `apply` are built over
//! `region.getPaddedSpan()`. A tool therefore reads reference bases and features outside the
//! territory it is allowed to call variants in, which is the whole point of the padding: assembling
//! over more sequence improves the calls inside the primary span.
//!
//! # `--force-active` changes the flag and not one boundary
//!
//! `forceActive` is applied in the traversal loop, **after** the iterator has produced the region:
//! `if (assemblyRegionArgs.forceActive) { assemblyRegion.setIsActive(true); }`. The regions were
//! already cut by the real activity, so the argument marks every one of them active without moving
//! a single edge. A port that folded it into the evaluator instead would return one enormous active
//! region, which is a different traversal with the same name.
//!
//! # The default filters are a `LocusWalker`'s, not a `ReadWalker`'s
//!
//! `WellformedReadFilter` **and** `ReadFilterLibrary.MappedReadFilter`. An unmapped read parked at
//! its mate's coordinate reaches a `ReadWalker` and never reaches this traversal.

use crate::assembly_region::RegionError;
use crate::assembly_region_iterator::{
    assembly_regions, group_intervals_by_contig, ActivityAnswer, AssemblyRegionArgs, ReadShard,
    RegionCutter, RegionFiller, TraversedRegion,
};
use crate::interval::SimpleInterval;
use crate::java_random::JavaRandom;
use crate::locus_iterator::AlignmentContext;
use crate::locus_iterator::{self, LocusIteratorOptions};
use crate::locus_shards::for_each_interval_context;
use crate::read_states::{DownsamplingInfo, ReadStateError, ReadStateManager};
use crate::read_utils;
use crate::variant_source::Located;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// `getDefaultReadFilters()`, in the order the base class declares them.
///
/// Kept as names rather than as functions because what is being ported here is the *list*: which
/// filters a tool inherits before the user touches anything. Applying them is the `gatk-readfilter`
/// crate's job in the tools that use it.
pub const DEFAULT_READ_FILTERS: [&str; 2] = ["WellformedReadFilter", "MappedReadFilter"];

/// `createDownsampler()`: a `PositionalDownsampler` **only** when the argument is above zero.
///
/// Zero is not "downsample to zero" and not an error: it is the absence of a downsampler. A
/// negative value never reaches here, because `validate()` refuses it first.
pub fn needs_downsampler(args: &AssemblyRegionArgs) -> bool {
    args.max_reads_per_alignment_start > 0
}

/// `makeReadShards`: group the traversal intervals by contig, one shard per group.
///
/// `None` is the shard constructor's `Utils.validateArg(intervalPadding >= 0, ...)`, which cannot
/// fire through the command line because `validate()` has already refused a negative padding.
pub fn make_read_shards(
    intervals: &[SimpleInterval],
    padding: i32,
    header: &SamHeader,
) -> Option<Vec<ReadShard>> {
    group_intervals_by_contig(intervals)
        .iter()
        .map(|group| ReadShard::new(group, padding, header))
        .collect()
}

/// What the traversal refuses, from either of the two layers that can refuse.
///
/// Keeping them apart matters for one case in particular: disabling the default read filters lets
/// a read with no read group through, and the refusal comes from the **sample partition** inside
/// `LocusIteratorByState` rather than from anything in the walker. The reference reports it as an
/// `IllegalStateException` reading "Offered read with sample name null to SamplePartitioner but
/// this sample wasn't provided as one of possible samples at construction", which names the layer.
#[derive(Debug, Clone, PartialEq)]
pub enum WalkerError {
    Region(RegionError),
    ReadState(ReadStateError),
    /// The reads source failing under a shard's query, which a run only meets when the file
    /// changes under it: every read was decoded once before the traversal.
    Reads(String),
}

impl From<RegionError> for WalkerError {
    fn from(error: RegionError) -> WalkerError {
        WalkerError::Region(error)
    }
}

impl From<ReadStateError> for WalkerError {
    fn from(error: ReadStateError) -> WalkerError {
        WalkerError::ReadState(error)
    }
}

impl WalkerError {
    /// The Java class the reference throws, as the dump prints it.
    pub fn class(&self) -> &'static str {
        match self {
            WalkerError::ReadState(ReadStateError::UndeclaredSample(_)) => {
                "java.lang.IllegalStateException"
            }
            WalkerError::Region(error) => error.class(),
            WalkerError::Reads(_) => "UserException",
            _ => "java.lang.IllegalStateException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            WalkerError::ReadState(ReadStateError::UndeclaredSample(sample)) => format!(
                "Offered read with sample name {} to SamplePartitioner but this sample wasn't \
                 provided as one of possible samples at construction",
                sample.clone().unwrap_or_else(|| "null".to_string())
            ),
            WalkerError::Region(error) => error.message(),
            WalkerError::ReadState(other) => format!("{other:?}"),
            WalkerError::Reads(message) => message.clone(),
        }
    }
}

/// One locus as the traversal sees it: where it is and how deep it is.
///
/// The evaluator upstream receives the whole pileup, the reference context and the feature context.
/// A conformance probe needs the depth and nothing else, and passing the depth rather than the
/// pileup is what keeps this module free of the pileup's borrows.
pub struct LocusDepth {
    pub contig: String,
    pub position: i32,
    pub depth: usize,
}

impl Located for LocusDepth {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.position
    }
    fn stop(&self) -> i32 {
        self.position
    }
}

/// `AssemblyRegionWalker.traverse()`: every region every shard produces, in order.
///
/// `reads` must already have been through the read filters, as the shard's iterator applies them
/// before the traversal sees anything. `is_active` is the tool's `AssemblyRegionEvaluator`, given
/// the locus and its depth.
///
/// The shard loop is the reference's and its order matters: the shards come out of
/// `groupIntervalsByContig` in interval order, so the regions of chr1 all precede the regions of
/// chr2 even when the reads are interleaved.
pub fn traverse(
    reads: &[BamRecord],
    intervals: &[SimpleInterval],
    samples: &[Option<String>],
    args: &AssemblyRegionArgs,
    header: &SamHeader,
    is_active: &dyn Fn(&LocusDepth) -> f64,
) -> Result<Vec<TraversedRegion>, WalkerError> {
    let shards = make_read_shards(intervals, args.assembly_region_padding, header)
        .expect("a padding validate() has already accepted");

    let mut out = Vec::new();
    for shard in &shards {
        // The shard's iterator queries the reads over its **padded** intervals.
        let shard_reads: Vec<BamRecord> = reads
            .iter()
            .filter(|record| {
                let contig = contig_of(record, header);
                shard.padded_intervals.iter().any(|interval| {
                    interval.overlaps(&contig, read_utils::start(record), read_utils::end(record))
                })
            })
            .cloned()
            .collect();

        let states = ReadStateManager::new(samples.to_vec(), DownsamplingInfo::NONE)?;
        // `new LocusIteratorByState(reads, NONE, samples, header, true)`: the five-argument
        // constructor, deletions included and reads with an N at the locus left out.
        let contexts = locus_iterator::contexts(
            &shard_reads,
            samples.to_vec(),
            header,
            REGION_LOCUS_OPTIONS,
            states,
        )?;

        let loci: Vec<LocusDepth> = contexts
            .iter()
            .map(|context| LocusDepth {
                contig: context.contig.clone(),
                position: context.position,
                depth: context.pileup.size(),
            })
            .collect();

        let regions = assembly_regions(
            &loci,
            &shard_reads,
            shard,
            args,
            header,
            &|locus, context| match context {
                Some(context) => is_active(context),
                // A manufactured empty context. The evaluator is still called on it, with a pileup
                // of depth zero at that locus, and it is free to call the locus active: this is
                // where a port that skipped the empty contexts would silently stop asking.
                None => is_active(&LocusDepth {
                    contig: locus.interval.contig.clone(),
                    position: locus.interval.start,
                    depth: 0,
                }),
            },
        )?;
        out.extend(regions);
    }
    Ok(out)
}

/// What `AssemblyRegionIterator` builds its `LocusIteratorByState` with.
const REGION_LOCUS_OPTIONS: LocusIteratorOptions = LocusIteratorOptions {
    include_deletions: true,
    include_ns: false,
};

impl Located for AlignmentContext<'_> {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.position
    }
    fn stop(&self) -> i32 {
        self.position
    }
}

/// A shard's reads, one at a time, in the order the traversal reads them.
pub type ReadStream<'r> = Box<dyn Iterator<Item = Result<BamRecord, String>> + 'r>;

/// Where [`traverse_with_pileups`] reads a shard from: `query` answers the reads overlapping the
/// shard's padded intervals ([`overlaps_shard`]), through the tool's transformers and filters, in
/// file order, and is asked twice per shard.
pub struct WalkerReads<'r> {
    pub query: Box<dyn Fn(&ReadShard) -> ReadStream<'r> + 'r>,
    /// `Some(width)` when every stream is coordinate sorted: the loci are then walked `width` at a
    /// time and only the reads overlapping the window are held. `None` holds a whole shard.
    pub window: Option<i32>,
}

/// How wide [`WalkerReads::window`] walks: the reads held are the depth times this plus a read.
pub const WINDOW: i32 = 20_000;

impl<'r> WalkerReads<'r> {
    /// Every filtered read of the file, in hand: each shard picks out its own.
    pub fn held(reads: &'r [BamRecord], header: &'r SamHeader) -> WalkerReads<'r> {
        WalkerReads {
            query: Box::new(move |shard: &ReadShard| {
                let shard = shard.clone();
                Box::new(
                    reads
                        .iter()
                        .filter(move |read| overlaps_shard(read, &shard, header))
                        .cloned()
                        .map(Ok),
                )
            }),
            window: None,
        }
    }
}

/// The shard's iterator: a read is in it when it overlaps one of the **padded** intervals.
pub fn overlaps_shard(record: &BamRecord, shard: &ReadShard, header: &SamHeader) -> bool {
    let contig = contig_of(record, header);
    shard.padded_intervals.iter().any(|interval| {
        interval.overlaps(&contig, read_utils::start(record), read_utils::end(record))
    })
}

/// `createDownsampler()`, as the tool builds it, run a batch at a time: a batch ends where the
/// downsampler starts over, so running it batch by batch draws what running it over the whole
/// shard draws.
#[derive(Clone, Copy)]
pub struct Downsampler<'d> {
    /// Whether `next` starts a new batch after `previous`, in a batch that began with `first`.
    pub new_batch: &'d dyn Fn(&BamRecord, &BamRecord, &BamRecord) -> bool,
    /// The downsampler over one batch, drawing from the generator it is handed.
    pub run: &'d dyn Fn(&[BamRecord], &mut JavaRandom) -> Vec<BamRecord>,
}

/// A stream through a [`Downsampler`], a batch at a time.
struct Downsampled<'r, R> {
    input: std::iter::Peekable<ReadStream<'r>>,
    downsampler: Downsampler<'r>,
    random: R,
    out: std::collections::VecDeque<BamRecord>,
}

impl<'r, R: std::borrow::BorrowMut<JavaRandom>> Downsampled<'r, R> {
    fn new(input: ReadStream<'r>, downsampler: Downsampler<'r>, random: R) -> Self {
        Downsampled {
            input: input.peekable(),
            downsampler,
            random,
            out: Default::default(),
        }
    }
}

impl<R: std::borrow::BorrowMut<JavaRandom>> Iterator for Downsampled<'_, R> {
    type Item = Result<BamRecord, String>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(read) = self.out.pop_front() {
            return Some(Ok(read));
        }
        let mut batch = match self.input.next()? {
            Ok(first) => vec![first],
            Err(error) => return Some(Err(error)),
        };
        while let Some(Ok(next)) = self.input.peek() {
            let previous = batch.last().expect("a batch starts with a read");
            if (self.downsampler.new_batch)(&batch[0], previous, next) {
                break;
            }
            batch.push(
                self.input
                    .next()
                    .expect("just peeked")
                    .expect("just peeked"),
            );
        }
        self.out = (self.downsampler.run)(&batch, self.random.borrow_mut()).into();
        self.next()
    }
}

/// One shard's regions, filled as they are drawn from a second read of the shard.
pub struct ShardRegions<'h, 'r> {
    filler: RegionFiller<'h, Box<dyn Iterator<Item = BamRecord> + 'r>>,
    failure: std::rc::Rc<std::cell::RefCell<Option<String>>>,
}

impl Iterator for ShardRegions<'_, '_> {
    type Item = Result<TraversedRegion, WalkerError>;

    fn next(&mut self) -> Option<Self::Item> {
        let region = self.filler.next()?;
        // A region the stream failed under is short of reads: the failure is the answer.
        if let Some(error) = self.failure.borrow_mut().take() {
            return Some(Err(WalkerError::Reads(error)));
        }
        Some(region.map_err(WalkerError::from))
    }
}

/// `traverse()` as a tool with a real evaluator runs it: each shard's reads through the tool's
/// downsampler (`createDownsampler()`), then the evaluator handed the locus and its
/// `AlignmentContext`, `None` for a locus no read covers. The evaluator answers the probability
/// and, for a `HIGH_QUALITY_SOFT_CLIPS` state, its value.
///
/// Every locus of a shard is evaluated, and every draw its downsampler makes from `random` taken,
/// before the first of its regions is handed out, where the reference pulls reads lazily as the
/// loci advance and calls each region as it is ready: the draws are the same, and only their
/// interleaving with a draw the tool makes while calling a region could differ. That order is
/// kept here without holding the shard: the loci are walked over one read of it, and the regions
/// are filled from a second, through the downsampler again from a copy of `random` taken before
/// the first, which draws the same reads.
#[allow(clippy::too_many_arguments)]
pub fn traverse_with_pileups<'h, 's, 'r: 's, 'd: 's>(
    reads: &WalkerReads<'r>,
    intervals: &[SimpleInterval],
    samples: &[Option<String>],
    args: &AssemblyRegionArgs,
    header: &'h SamHeader,
    downsampler: Downsampler<'d>,
    random: &mut JavaRandom,
    is_active: &mut dyn FnMut(&SimpleInterval, Option<&AlignmentContext<'_>>) -> ActivityAnswer,
) -> Result<Vec<ShardRegions<'h, 's>>, WalkerError> {
    let shards = make_read_shards(intervals, args.assembly_region_padding, header)
        .expect("a padding validate() has already accepted");

    let mut out = Vec::new();
    for shard in &shards {
        let replay = random.clone();
        let regions = walk_loci(
            reads,
            shard,
            samples,
            args,
            header,
            downsampler,
            random,
            is_active,
        )?;
        let failure = std::rc::Rc::new(std::cell::RefCell::new(None));
        let slot = failure.clone();
        let input: ReadStream<'s> = (reads.query)(shard);
        let stream: Box<dyn Iterator<Item = BamRecord> + 's> = Box::new(
            Downsampled::new(input, downsampler, replay)
                .map_while(move |read| read.map_err(|error| *slot.borrow_mut() = Some(error)).ok()),
        );
        out.push(ShardRegions {
            filler: RegionFiller::new(stream, regions, args.force_active, header),
            failure,
        });
    }
    Ok(out)
}

/// The first read of a shard: every locus evaluated and the regions cut.
///
/// With a window the loci go `width` at a time, each window over a fresh `LocusIteratorByState`
/// fed the reads overlapping it and nothing else, whose contexts before the window are dropped.
/// A sorted stream makes that the same pileups: every read is taken in at its own start, so the
/// states at a locus are the reads overlapping it in stream order, whatever came before.
#[allow(clippy::too_many_arguments)]
fn walk_loci(
    reads: &WalkerReads<'_>,
    shard: &ReadShard,
    samples: &[Option<String>],
    args: &AssemblyRegionArgs,
    header: &SamHeader,
    downsampler: Downsampler<'_>,
    random: &mut JavaRandom,
    is_active: &mut dyn FnMut(&SimpleInterval, Option<&AlignmentContext<'_>>) -> ActivityAnswer,
) -> Result<Vec<crate::assembly_region::AssemblyRegion>, WalkerError> {
    // Each window is its loci and the bounds of its contexts, none for a whole shard.
    type Window = (Vec<SimpleInterval>, Option<(i32, i32)>);
    let windows: Vec<Window> = match reads.window {
        None => vec![(shard.intervals.clone(), None)],
        Some(width) => shard
            .intervals
            .iter()
            .flat_map(|interval| {
                (interval.start..=interval.end)
                    .step_by(width.max(1) as usize)
                    .map(move |start| {
                        let end = start.saturating_add(width.max(1) - 1).min(interval.end);
                        (
                            vec![SimpleInterval {
                                contig: interval.contig.clone(),
                                start,
                                end,
                            }],
                            Some((start, end)),
                        )
                    })
            })
            .collect(),
    };

    let mut stream = Downsampled::new((reads.query)(shard), downsampler, &mut *random).peekable();
    let mut held: Vec<BamRecord> = Vec::new();
    let mut cutter = RegionCutter::new(shard, args, header);
    for (loci, bounds) in &windows {
        if let Some((first, _)) = bounds {
            held.retain(|read| read_utils::end(read) >= *first);
        }
        loop {
            match stream.peek() {
                Some(Ok(read))
                    if bounds.is_none_or(|(_, last)| read_utils::start(read) <= last) => {}
                Some(Err(_)) => {
                    let error = stream
                        .next()
                        .expect("just peeked")
                        .expect_err("just peeked");
                    return Err(WalkerError::Reads(error));
                }
                _ => break,
            }
            held.push(stream.next().expect("just peeked").expect("just peeked"));
        }

        let states = ReadStateManager::new(samples.to_vec(), DownsamplingInfo::NONE)?;
        // Drawn one locus at a time: a contig's pileups held together are gigabytes.
        let mut failure = None;
        let contexts = locus_iterator::Contexts::new(
            &held,
            samples.to_vec(),
            header,
            REGION_LOCUS_OPTIONS,
            states,
            &|_| true,
        )
        .map_while(|context| context.map_err(|error| failure = Some(error)).ok())
        .skip_while(|context| bounds.is_some_and(|(first, _)| context.position < first))
        .take_while(|context| bounds.is_none_or(|(_, last)| context.position <= last));
        let walked = for_each_interval_context(contexts, loci, header, |locus, context| {
            cutter.add(&locus, || is_active(&locus.interval, context))
        });
        if let Some(error) = failure {
            return Err(error.into());
        }
        walked?;
    }
    drop(held);
    // The rest of the shard still goes through the downsampler: its draws come before any region's.
    for read in stream {
        read.map_err(WalkerError::Reads)?;
    }
    Ok(cutter.finish()?)
}

fn contig_of(record: &BamRecord, header: &SamHeader) -> String {
    usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.clone())
        .unwrap_or_else(|| "null".to_string())
}
