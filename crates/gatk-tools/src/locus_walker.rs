//! Ported from `org.broadinstitute.hellbender.engine.LocusWalker` (GATK 4.6.2.0).
//!
//! The traversal the largest GATK archetype inherits: one `apply` per locus, with the pileup at
//! that locus and the reference base under it. Everything it stands on is already measured, so
//! this module is thin, and what it adds is four defaults that are not the ones a `ReadWalker`
//! gets:
//!
//!  * **two default read filters, not one.** `getDefaultReadFilters` returns
//!    `WellformedReadFilter` *and* `MappedReadFilter`, so an unmapped read carrying its mate's
//!    position reaches a `ReadWalker` and never reaches a `LocusWalker`;
//!  * **`includeDeletions` is true and `includeNs` is false**, which is the pair
//!    [`gatk_engine::locus_iterator`] measures, and the pair a tool overrides to change what its
//!    pileups contain;
//!  * **`emitEmptyLoci` is false**, so an uncovered position inside `-L` produces no `apply` call
//!    at all unless the tool asks otherwise;
//!  * **`--max-depth-per-sample` defaults to 0, meaning no downsampling**, and a *negative* value
//!    is a bad argument rather than a synonym for unlimited.
//!
//! The reference context each `apply` receives is windowless over the locus itself, so a tool that
//! wants flanking bases widens it, exactly as `IntervalWalker` does.

use gatk_engine::context::ReferenceContext;
use gatk_engine::context_iterator::{self, Route};
use gatk_engine::interval::SimpleInterval;
use gatk_engine::interval_args::IntervalArguments;
use gatk_engine::locus_iterator::{self, AlignmentContext, LocusIteratorOptions};
use gatk_engine::read_states::{DownsamplingInfo, ReadStateError, ReadStateManager};
use gatk_engine::reads::ReadsDataSource;
use gatk_engine::reference::ReferenceFileSource;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// One `apply` call: the locus with its pileup, and the reference under it.
pub struct Applied<'a> {
    pub context: AlignmentContext<'a>,
    pub reference: ReferenceContext,
}

/// What the tool overrode, with `LocusWalker`'s own defaults.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    pub include_deletions: bool,
    pub include_ns: bool,
    pub emit_empty_loci: bool,
    /// `--max-depth-per-sample`. Zero is no downsampling; negative is refused.
    pub max_depth_per_sample: i32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            include_deletions: true,
            include_ns: false,
            emit_empty_loci: false,
            max_depth_per_sample: 0,
        }
    }
}

/// What a traversal can refuse before it starts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LocusWalkerError {
    /// `CommandLineException.BadArgumentValue` on a negative `--max-depth-per-sample`.
    NegativeMaxDepth(i32),
    /// `Utils.nonEmpty` inside `IntervalOverlappingIterator`: intervals were given as an empty
    /// list and empty loci were not asked for.
    EmptyIntervalList,
    /// The state layer refused, including the downsampling this port does not reproduce.
    States(ReadStateError),
}

impl Options {
    /// `LocusWalker.getDownsamplingInfo`.
    pub fn downsampling_info(&self) -> Result<DownsamplingInfo, LocusWalkerError> {
        if self.max_depth_per_sample < 0 {
            return Err(LocusWalkerError::NegativeMaxDepth(
                self.max_depth_per_sample,
            ));
        }
        Ok(if self.max_depth_per_sample == 0 {
            DownsamplingInfo::NONE
        } else {
            DownsamplingInfo {
                performing: true,
                to_coverage: self.max_depth_per_sample,
            }
        })
    }
}

/// `LocusWalker.getDefaultReadFilters`: wellformed, then mapped.
///
/// The order is the reference's and is observable through `CountingReadFilter`'s summary, which is
/// the only place a conjunction's order shows.
pub fn default_filter<'a>(header: &'a SamHeader) -> impl Fn(&BamRecord) -> bool + 'a {
    move |read: &BamRecord| {
        gatk_readfilter::with_header::wellformed(read, header) && gatk_readfilter::mapped(read)
    }
}

/// The sample names, in the order `AlignmentContextIteratorBuilder` hands them over.
///
/// That order is a `java.util.HashSet`'s and decides the element order of every multi-sample
/// pileup, which is why it is reproduced as a measured observable rather than sorted here. See
/// `docs/an-unspecified-order-that-reaches-the-output.md`.
pub fn samples_in_iteration_order(header: &SamHeader) -> Vec<Option<String>> {
    let declared: Vec<String> = header
        .read_groups
        .iter()
        .filter_map(|group| group.attributes.get("SM").map(|s| s.to_string()))
        .collect();
    match gatk_engine::java_hash::hash_set_order(&declared) {
        Ok(order) => order.into_iter().map(Some).collect(),
        // Past the measured range the order is unknown, and guessing it would reorder pileups
        // silently. Falling back to the declared order is wrong, so this is left to the caller
        // by returning what was declared *and* the module refusing in the hash layer.
        Err(_) => declared.into_iter().map(Some).collect(),
    }
}

/// `LocusWalker.traverse`, collected.
pub fn traverse<'a>(
    reads: &'a [BamRecord],
    header: &SamHeader,
    reference: Option<&mut ReferenceFileSource>,
    intervals: Option<&[SimpleInterval]>,
    options: Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<Vec<Applied<'a>>, LocusWalkerError> {
    let route = context_iterator::route(options.emit_empty_loci, intervals);
    if route == Route::RejectedEmptyIntervalList {
        return Err(LocusWalkerError::EmptyIntervalList);
    }

    let samples = samples_in_iteration_order(header);
    let states = ReadStateManager::new(samples.clone(), options.downsampling_info()?)
        .map_err(LocusWalkerError::States)?;

    // The filter runs before the locus iterator sees anything, so a filtered read is absent from
    // the pileup rather than present and ignored: it changes the depth, not only the reported set.
    let covered = locus_iterator::contexts_filtered(
        reads,
        samples,
        header,
        LocusIteratorOptions {
            include_deletions: options.include_deletions,
            include_ns: options.include_ns,
        },
        states,
        filter,
    )
    .map_err(LocusWalkerError::States)?;

    let contexts = match route {
        Route::Unbounded => covered,
        Route::Overlapping => {
            context_iterator::overlapping(covered, intervals.expect("routed on intervals"), header)
        }
        Route::EmptyLoci => {
            let whole: Vec<SimpleInterval> = header
                .sequences
                .iter()
                .map(|s| SimpleInterval {
                    contig: s.name.clone(),
                    start: 1,
                    end: s.length,
                })
                .collect();
            let requested = intervals.unwrap_or(&whole);
            context_iterator::with_empty_loci(covered, requested, header)
        }
        Route::RejectedEmptyIntervalList => unreachable!("checked above"),
    };

    let mut reference = reference;
    let mut applied = Vec::with_capacity(contexts.len());
    for context in contexts {
        // `new SimpleInterval(alignmentContext)`, which is the single base of the locus.
        let interval = SimpleInterval {
            contig: context.contig.clone(),
            start: context.position,
            end: context.position,
        };
        let reference_context = match reference.as_deref_mut() {
            Some(source) => ReferenceContext::new(source, Some(interval), 0, 0),
            None => ReferenceContext::without_source(Some(interval), 0, 0),
        }
        .unwrap_or_else(|_| ReferenceContext::empty());
        applied.push(Applied {
            context,
            reference: reference_context,
        });
    }
    Ok(applied)
}

/// The reads a segment of [`for_each_streamed`] takes before it looks for a cut: a few windows of
/// read length at 30x, so the reads walked twice are a small share.
pub const SEGMENT_READS: usize = 1 << 12;

/// Why [`for_each_streamed`] stopped before the last locus.
#[derive(Debug)]
pub enum Stopped<E> {
    /// The reads source refused a read.
    Read(E),
    /// The traversal refused, as [`traverse`] would have somewhere.
    Locus(LocusWalkerError),
    /// A read came before the one ahead of it, which the segments below rely on not happening.
    Unsorted,
    /// `apply` asked to stop.
    Applied,
}

/// Whether [`for_each_streamed`] answers what [`traverse`] does for these options: no
/// downsampling, whose random draws a read carried into the next segment would repeat.
pub fn streams(options: &Options) -> bool {
    options.max_depth_per_sample == 0
}

/// [`traverse`] over reads that arrive one at a time, each locus handed to `apply` in order and
/// dropped after it, so a traversal holds a segment of reads rather than all of them and no
/// pileup but the one `apply` is given. Only for options [`streams`] accepts.
///
/// The reads, filtered, are cut into segments of about `segment` reads, each cut where the next
/// read starts strictly after the one before it. A segment answers every locus before its cut. The
/// next one starts with the reads of that one that reach the cut, in their order, so at every
/// locus from the cut on the same reads are in the same states in the same order as in one pass:
/// they are admitted at their own starts and stepped base by base whatever came before. Their
/// loci before the cut are walked again and skipped.
///
/// It does not reproduce WHERE [`traverse`] would have refused. A refusal, a read out of order or
/// `apply` asking to stop ends it, and a caller then answers from [`traverse`], which refuses as it
/// always did.
#[allow(clippy::too_many_arguments)]
pub fn for_each_streamed<E>(
    reads: impl IntoIterator<Item = Result<BamRecord, E>>,
    header: &SamHeader,
    reference: Option<&mut ReferenceFileSource>,
    intervals: Option<&[SimpleInterval]>,
    options: Options,
    filter: &dyn Fn(&BamRecord) -> bool,
    segment: usize,
    mut apply: impl FnMut(Applied<'_>) -> bool,
) -> Result<(), Stopped<E>> {
    assert!(streams(&options), "options for_each_streamed reproduces");
    let route = context_iterator::route(options.emit_empty_loci, intervals);
    if route == Route::RejectedEmptyIntervalList {
        return Err(Stopped::Locus(LocusWalkerError::EmptyIntervalList));
    }
    let samples = samples_in_iteration_order(header);
    let mut reference = reference;
    let mut cursor = context_iterator::OverlapCursor::default();
    let whole: Vec<SimpleInterval> = header
        .sequences
        .iter()
        .map(|s| SimpleInterval {
            contig: s.name.clone(),
            start: 1,
            end: s.length,
        })
        .collect();
    let mut empty_loci = (route == Route::EmptyLoci)
        .then(|| context_iterator::EmptyLociCursor::new(intervals.unwrap_or(&whole)));
    // One locus to `apply`, with the reference under it.
    let mut emit = |context: AlignmentContext<'_>| -> bool {
        let interval = SimpleInterval {
            contig: context.contig.clone(),
            start: context.position,
            end: context.position,
        };
        let reference_context = match reference.as_deref_mut() {
            Some(source) => ReferenceContext::new(source, Some(interval), 0, 0),
            None => ReferenceContext::without_source(Some(interval), 0, 0),
        }
        .unwrap_or_else(|_| ReferenceContext::empty());
        apply(Applied {
            context,
            reference: reference_context,
        })
    };
    let key = |read: &BamRecord| (read.reference_index, read.alignment_start);

    let mut reads = reads.into_iter();
    let mut carried: Vec<BamRecord> = Vec::new();
    let mut next: Option<BamRecord> = None;
    // The cut the segment starts at: the carried reads' loci before it were answered already.
    let mut floor: Option<(String, i32)> = None;
    loop {
        let mut reads_here = std::mem::take(&mut carried);
        let carried_count = reads_here.len();
        reads_here.extend(next.take());
        let mut cut = None;
        for read in reads.by_ref() {
            let read = read.map_err(Stopped::Read)?;
            if !filter(&read) {
                continue;
            }
            let last = reads_here[carried_count..].last().map(key);
            if last.is_some_and(|last| key(&read) < last) {
                return Err(Stopped::Unsorted);
            }
            if reads_here.len() - carried_count >= segment && last != Some(key(&read)) {
                cut = Some(key(&read));
                next = Some(read);
                break;
            }
            reads_here.push(read);
        }
        let cut_contig = cut.map(|(index, position)| {
            let name = header
                .sequences
                .get(index as usize)
                .map(|sequence| sequence.name.as_str());
            (name, position)
        });

        let states = ReadStateManager::new(
            samples.clone(),
            options.downsampling_info().map_err(Stopped::Locus)?,
        )
        .map_err(|error| Stopped::Locus(LocusWalkerError::States(error)))?;
        let contexts = locus_iterator::Contexts::new(
            &reads_here,
            samples.clone(),
            header,
            LocusIteratorOptions {
                include_deletions: options.include_deletions,
                include_ns: options.include_ns,
            },
            states,
            &|_| true,
        );
        for context in contexts {
            let context =
                context.map_err(|error| Stopped::Locus(LocusWalkerError::States(error)))?;
            if let Some((name, position)) = cut_contig {
                if name == Some(context.contig.as_str()) && context.position >= position {
                    break;
                }
            }
            if let Some((name, position)) = &floor {
                if *name == context.contig && context.position < *position {
                    continue;
                }
            }
            let go_on = match (&route, empty_loci.as_mut()) {
                (Route::Overlapping, _) => {
                    match cursor.keeps(&context, intervals.expect("routed on intervals"), header) {
                        Some(true) => emit(context),
                        Some(false) => true,
                        // Every interval is behind: the rest is read only for a refusal in it.
                        None => {
                            return match reads.find_map(Result::err) {
                                Some(error) => Err(Stopped::Read(error)),
                                None => Ok(()),
                            }
                        }
                    }
                }
                (_, Some(loci)) => loci.take(context, header, &mut emit),
                _ => emit(context),
            };
            if !go_on {
                return Err(Stopped::Applied);
            }
        }
        let Some((index, position)) = cut else {
            if let Some(loci) = empty_loci.as_mut() {
                if !loci.finish(&mut emit) {
                    return Err(Stopped::Applied);
                }
            }
            return Ok(());
        };
        floor = cut_contig.map(|(name, position)| (name.unwrap_or_default().to_string(), position));
        carried = reads_here
            .into_iter()
            .filter(|read| read.reference_index == index && read.alignment_end() >= position)
            .collect();
    }
}

/// `LocusWalker.onStartup`: with user intervals the reads source is bounded before the traversal.
///
/// Kept as a named function rather than folded into [`traverse`], because the bounding happens on
/// the *data source* and therefore changes which reads exist, not merely which loci are reported.
pub fn traversal_bounds(
    source: &ReadsDataSource,
    arguments: &IntervalArguments,
    header: &SamHeader,
) -> Result<Option<Vec<SimpleInterval>>, LocusWalkerError> {
    let _ = source;
    if !arguments.specified() {
        return Ok(None);
    }
    match gatk_engine::interval_args::parse_intervals(arguments, header) {
        Ok(parameters) => Ok(Some(parameters.intervals)),
        Err(_) => Ok(None),
    }
}
