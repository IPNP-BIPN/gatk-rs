//! The assembly-region traversal, ported from `engine.AssemblyRegionIterator`,
//! `engine.MultiIntervalLocalReadShard` and `engine.spark.AssemblyRegionArgumentCollection`
//! (GATK 4.6.2.0).
//!
//! This is where [`crate::activity_profile`], [`crate::assembly_region`] and
//! [`crate::locus_shards`] meet: loci go in, assembly regions come out, and every region carries
//! the reads that overlap its padded span. `HaplotypeCaller` and `Mutect2` are this loop plus an
//! assembler.
//!
//! # A shard has two interval lists and they are not the same list
//!
//! `MultiIntervalLocalReadShard` builds both of them through `IntervalUtils.getIntervalsWithFlanks`,
//! the padded one with the padding and the unpadded one with **zero**. Passing zero is not a no-op:
//! the function sorts and merges with `IntervalMergingRule.ALL` whatever the padding is, so
//! `getIntervals()` is already sorted and already has its adjacent intervals joined. The reads are
//! queried over the padded list and the loci are walked over the unpadded one, so a read can be in
//! the shard without any of its bases being a locus.
//!
//! # The ordering that the comment insists on
//!
//! Inside the loop, the profile is popped **before** the current pileup is added, and the
//! `forceConversion` flag is `pileup.getStart() != profile.getEnd() + 1`. The upstream comment says
//! "Ordering matters here", and it does: adding first would make the profile contiguous with the
//! new pileup and the force would never fire, so a region would never be closed at a gap in the
//! loci. The gap is exactly what [`crate::locus_shards`] manufactures empty contexts to avoid
//! *inside* an interval; between two intervals the gap is real and this is what notices it.
//!
//! # A region is not ready when it is popped
//!
//! A region popped from the profile goes into a queue and only comes out once the locus iterator
//! has advanced **past the end of its padded span**, because until then the reads that belong in it
//! have not been read. That is the whole reason the read cache exists.
//!
//! # Why an eager port is faithful here
//!
//! Upstream the reads are pulled lazily and cached as the locus iterator consumes them; here the
//! contexts and the reads are both already in hand. The two agree because a region is only filled
//! after the locus iterator has passed the end of its padded span, so every read that could belong
//! to it has necessarily been consumed by then. The reference makes that explicit at the end of the
//! traversal, where it drains the underlying iterator by hand "to guarantee that the reads in the
//! final padded region end up in our read cache". A port that filled regions early would differ;
//! one that fills them at the same points cannot.

use crate::activity_profile::ActivityProfile;
use crate::assembly_region::{AssemblyRegion, RegionError};
use crate::interval::SimpleInterval;
use crate::interval_args::with_flanks;
use crate::locus_shards::{interval_alignment_contexts, is_after, EmittedLocus};
use crate::read_utils;
use crate::variant_source::Located;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// `AssemblyRegionArgumentCollection.DEFAULT_MIN_ASSEMBLY_REGION_SIZE`.
pub const DEFAULT_MIN_ASSEMBLY_REGION_SIZE: i32 = 50;
/// `DEFAULT_MAX_ASSEMBLY_REGION_SIZE`.
pub const DEFAULT_MAX_ASSEMBLY_REGION_SIZE: i32 = 300;
/// `DEFAULT_ASSEMBLY_REGION_PADDING`.
pub const DEFAULT_ASSEMBLY_REGION_PADDING: i32 = 100;
/// `DEFAULT_MAX_READS_PER_ALIGNMENT`.
pub const DEFAULT_MAX_READS_PER_ALIGNMENT: i32 = 50;
/// `DEFAULT_ACTIVE_PROB_THRESHOLD`.
pub const DEFAULT_ACTIVE_PROB_THRESHOLD: f64 = 0.002;
/// `DEFAULT_MAX_PROB_PROPAGATION_DISTANCE`.
pub const DEFAULT_MAX_PROB_PROPAGATION_DISTANCE: i32 = 50;

/// `AssemblyRegionArgumentCollection`, with the defaults the base class declares.
///
/// A tool may override any of the `defaultX()` methods, so these are the *base* defaults and not
/// necessarily `HaplotypeCaller`'s. Keeping them as constants rather than inlining them is what
/// makes an override visible as a difference.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyRegionArgs {
    pub min_assembly_region_size: i32,
    pub max_assembly_region_size: i32,
    pub active_prob_threshold: f64,
    pub max_prob_propagation_distance: i32,
    pub force_active: bool,
    pub assembly_region_padding: i32,
    pub max_reads_per_alignment_start: i32,
    pub indel_padding_for_genotyping: i32,
    pub snp_padding_for_genotyping: i32,
    pub str_padding_for_genotyping: i32,
    pub max_extension_into_region_padding: i32,
}

impl Default for AssemblyRegionArgs {
    fn default() -> Self {
        AssemblyRegionArgs {
            min_assembly_region_size: DEFAULT_MIN_ASSEMBLY_REGION_SIZE,
            max_assembly_region_size: DEFAULT_MAX_ASSEMBLY_REGION_SIZE,
            active_prob_threshold: DEFAULT_ACTIVE_PROB_THRESHOLD,
            max_prob_propagation_distance: DEFAULT_MAX_PROB_PROPAGATION_DISTANCE,
            force_active: false,
            assembly_region_padding: DEFAULT_ASSEMBLY_REGION_PADDING,
            max_reads_per_alignment_start: DEFAULT_MAX_READS_PER_ALIGNMENT,
            indel_padding_for_genotyping: 75,
            snp_padding_for_genotyping: 20,
            str_padding_for_genotyping: 75,
            max_extension_into_region_padding: 25,
        }
    }
}

/// What `AssemblyRegionArgumentCollection.validate` refuses, in the order it checks.
///
/// Every one of these is a `CommandLineException.BadArgumentValue`, but Barclay renders its two
/// constructors differently and the golden is what said so: the one-argument form prefixes
/// `Illegal argument value: `, while the two-argument form names the argument and its value. The
/// last two refusals use the second form, and they are the only two that do.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ArgsError {
    SizeNotPositive,
    MinAboveMax,
    NegativePadding,
    NegativeMaxReads,
    NegativeSnpPadding(i32),
    NegativeIndelPadding(i32),
}

impl ArgsError {
    /// The class the reference throws, as the dump prints it.
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue"
    }

    /// The message, verbatim. The `< 0` runs into the value without a space because upstream
    /// concatenates `"" + snpPaddingForGenotyping + "< 0"`.
    pub fn message(&self) -> String {
        match self {
            ArgsError::SizeNotPositive => {
                "Illegal argument value: min/max assembly region size must be > 0".to_string()
            }
            ArgsError::MinAboveMax => {
                "Illegal argument value: minAssemblyRegionSize must be <= maxAssemblyRegionSize"
                    .to_string()
            }
            ArgsError::NegativePadding => {
                "Illegal argument value: assemblyRegionPadding must be >= 0".to_string()
            }
            ArgsError::NegativeMaxReads => {
                "Illegal argument value: maxReadsPerAlignmentStart must be >= 0".to_string()
            }
            ArgsError::NegativeSnpPadding(value) => {
                format!("Argument paddingAroundSNPs has a bad value: {value}< 0")
            }
            ArgsError::NegativeIndelPadding(value) => {
                format!("Argument paddingAroundIndels has a bad value: {value}< 0")
            }
        }
    }
}

impl AssemblyRegionArgs {
    /// `AssemblyRegionArgumentCollection.validate()`.
    ///
    /// The order is the reference's and is observable: a collection with both a non-positive size
    /// and a negative padding reports the size, because that check comes first.
    pub fn validate(&self) -> Result<(), ArgsError> {
        if self.min_assembly_region_size <= 0 || self.max_assembly_region_size <= 0 {
            return Err(ArgsError::SizeNotPositive);
        }
        if self.min_assembly_region_size > self.max_assembly_region_size {
            return Err(ArgsError::MinAboveMax);
        }
        if self.assembly_region_padding < 0 {
            return Err(ArgsError::NegativePadding);
        }
        if self.max_reads_per_alignment_start < 0 {
            return Err(ArgsError::NegativeMaxReads);
        }
        if self.snp_padding_for_genotyping < 0 {
            return Err(ArgsError::NegativeSnpPadding(
                self.snp_padding_for_genotyping,
            ));
        }
        if self.indel_padding_for_genotyping < 0 {
            return Err(ArgsError::NegativeIndelPadding(
                self.indel_padding_for_genotyping,
            ));
        }
        Ok(())
    }
}

/// `IntervalUtils.groupIntervalsByContig`: one list per contig, over an **already sorted** input.
///
/// It groups on a change of contig rather than on the contig itself, so an unsorted list produces
/// two groups for the same contig instead of an error. The walker only ever hands it sorted
/// intervals, which is why nothing upstream notices.
pub fn group_intervals_by_contig(sorted: &[SimpleInterval]) -> Vec<Vec<SimpleInterval>> {
    let mut groups: Vec<Vec<SimpleInterval>> = Vec::new();
    let mut current: Vec<SimpleInterval> = Vec::new();
    let mut contig: Option<String> = None;

    for interval in sorted {
        if let Some(previous) = &contig {
            if previous != &interval.contig {
                groups.push(std::mem::take(&mut current));
            }
        }
        contig = Some(interval.contig.clone());
        current.push(interval.clone());
    }
    if !current.is_empty() {
        groups.push(current);
    }
    groups
}

/// `MultiIntervalLocalReadShard`: the two interval lists a shard carries.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadShard {
    /// `getIntervals()`, which is the input **sorted and merged** even though the padding is zero.
    pub intervals: Vec<SimpleInterval>,
    /// `getPaddedIntervals()`, padded then sorted and merged.
    pub padded_intervals: Vec<SimpleInterval>,
}

impl ReadShard {
    /// `new MultiIntervalLocalReadShard(intervals, intervalPadding, readsSource)`.
    ///
    /// `None` is the `Utils.validateArg(intervalPadding >= 0, ...)`.
    pub fn new(
        intervals: &[SimpleInterval],
        interval_padding: i32,
        header: &SamHeader,
    ) -> Option<ReadShard> {
        if interval_padding < 0 {
            return None;
        }
        Some(ReadShard {
            intervals: with_flanks(intervals.to_vec(), 0, header),
            padded_intervals: with_flanks(intervals.to_vec(), interval_padding, header),
        })
    }
}

/// One assembly region the traversal produced, with the reads it was filled with.
#[derive(Debug, Clone, PartialEq)]
pub struct TraversedRegion {
    pub region: AssemblyRegion,
    /// The activity the profile recorded, for the record: a region is active or not, and the
    /// per-locus probabilities that decided it are gone by the time it is emitted.
    pub is_active: bool,
}

/// A locus with its position, as the traversal sees it.
struct LocusPosition {
    contig: String,
    position: i32,
}

impl Located for LocusPosition {
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

/// `AssemblyRegionIterator`, run to exhaustion.
///
/// `contexts` is what `LocusIteratorByState` produced over the shard's reads, in order;
/// `is_active` is the tool's `AssemblyRegionEvaluator`, which upstream receives the pileup, the
/// reference context and the feature context and returns an `ActivityProfileState`. Here it
/// receives the index of the locus in the emitted sequence and the pileup size, which is all any
/// probe evaluator needs and all the conformance dump gives it.
///
/// `reads` must be the shard's reads in coordinate order, which is what the read cache guarantees.
pub fn assembly_regions<T: Located>(
    contexts: &[T],
    reads: &[BamRecord],
    shard: &ReadShard,
    args: &AssemblyRegionArgs,
    header: &SamHeader,
    is_active: &dyn Fn(&EmittedLocus, Option<&T>) -> f64,
) -> Result<Vec<TraversedRegion>, RegionError> {
    let loci = interval_alignment_contexts(contexts, &shard.intervals, header);

    let contig = shard
        .intervals
        .first()
        .map(|interval| interval.contig.clone())
        .unwrap_or_default();
    let contig_length = header
        .sequences
        .iter()
        .find(|sequence| sequence.name == contig)
        .map(|sequence| sequence.length)
        .unwrap_or(i32::MAX);

    let mut profile = ActivityProfile::band_pass(
        args.max_prob_propagation_distance,
        args.active_prob_threshold,
        // `BandPassActivityProfile.MAX_FILTER_SIZE` and `DEFAULT_SIGMA`, hard-coded at the call
        // site rather than taken from the arguments.
        crate::activity_profile::MAX_FILTER_SIZE,
        crate::activity_profile::DEFAULT_SIGMA,
        true,
        &contig,
        contig_length,
    );

    let mut pending: std::collections::VecDeque<AssemblyRegion> = Default::default();
    let mut ready: Vec<AssemblyRegion> = Vec::new();
    let mut previous_region_reads: Vec<BamRecord> = Vec::new();
    // The read cache: an index into `reads`, since the reads are already in coordinate order.
    let mut cache_index = 0usize;

    let push_popped = |profile: &mut ActivityProfile,
                       pending: &mut std::collections::VecDeque<AssemblyRegion>,
                       force: bool|
     -> Result<(), RegionError> {
        for popped in profile.pop_ready_regions(
            args.assembly_region_padding,
            args.min_assembly_region_size as usize,
            args.max_assembly_region_size as usize,
            force,
        ) {
            let region = AssemblyRegion::with_padding(
                popped.span.clone(),
                popped.is_active,
                args.assembly_region_padding,
                header,
            )?;
            pending.push_back(region);
        }
        Ok(())
    };

    for locus in loci.iter() {
        // Ordering matters: the profile is popped before the pileup is added, and the force flag
        // is "this locus does not continue the profile".
        if !profile.is_empty() {
            let force = locus.interval.start != profile.end() + 1;
            push_popped(&mut profile, &mut pending, force)?;
        }

        let context = locus.context.map(|position| &contexts[position]);
        profile.add(locus.interval.start, is_active(locus, context));

        // A pending region becomes ready only once the loci have advanced past the end of its
        // padded span, which is what guarantees its reads have been read.
        if let Some(front) = pending.front() {
            let here = LocusPosition {
                contig: locus.interval.contig.clone(),
                position: locus.interval.start,
            };
            let padded = front.padded_span().clone();
            if is_after(&here, &padded, header) == Some(true) {
                let region = pending.pop_front().expect("just peeked");
                ready.push(region);
            }
        }
    }

    // Out of loci: close the profile with forceConversion, then drain the queue.
    if !profile.is_empty() {
        push_popped(&mut profile, &mut pending, true)?;
    }
    while let Some(region) = pending.pop_front() {
        ready.push(region);
    }

    // Fill each region in turn, in the order they were emitted.
    let mut out = Vec::new();
    for mut region in ready {
        let padded = region.padded_span().clone();
        for read in &previous_region_reads {
            if padded.overlaps(
                &contig_of(read, header),
                read_utils::start(read),
                read_utils::end(read),
            ) {
                region.add(read.clone(), header)?;
            }
        }
        while cache_index < reads.len() {
            let read = &reads[cache_index];
            let located = LocusPosition {
                contig: contig_of(read, header),
                position: read_utils::start(read),
            };
            // `IntervalUtils.isAfter(read, paddedSpan, dict)` on the read, whose end matters: a
            // read is left in the cache only when its **start** is past the region's end.
            let read_interval = SimpleInterval {
                contig: located.contig.clone(),
                start: read_utils::start(read),
                end: read_utils::end(read),
            };
            if is_after(&read_interval, &padded, header) == Some(true) {
                break;
            }
            cache_index += 1;
            if padded.overlaps(
                &read_interval.contig,
                read_interval.start,
                read_interval.end,
            ) {
                region.add(read.clone(), header)?;
            }
        }

        previous_region_reads = region.reads().to_vec();
        let is_active = if args.force_active {
            region.set_is_active(true);
            true
        } else {
            region.is_active()
        };
        out.push(TraversedRegion { region, is_active });
    }

    Ok(out)
}

fn contig_of(record: &BamRecord, header: &SamHeader) -> String {
    usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.clone())
        .unwrap_or_else(|| "null".to_string())
}
