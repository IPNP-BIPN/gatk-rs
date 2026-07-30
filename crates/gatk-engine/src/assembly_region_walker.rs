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
    assembly_regions, group_intervals_by_contig, AssemblyRegionArgs, ReadShard, TraversedRegion,
};
use crate::interval::SimpleInterval;
use crate::locus_iterator::{self, LocusIteratorOptions};
use crate::read_states::{DownsamplingInfo, ReadStateManager};
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
) -> Result<Vec<TraversedRegion>, RegionError> {
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

        let states = ReadStateManager::new(samples.to_vec(), DownsamplingInfo::NONE)
            .expect("the traversal builds its LocusIteratorByState with DownsamplingMethod.NONE");
        // `new LocusIteratorByState(..., true)`: the keepUniqueReadList constructor, whose
        // deletion and N settings are both on, unlike a LocusWalker's.
        let contexts = locus_iterator::contexts(
            &shard_reads,
            samples.to_vec(),
            header,
            LocusIteratorOptions {
                include_deletions: true,
                include_ns: true,
            },
            states,
        )
        .expect("the shard's reads are what the iterator was built from");

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

fn contig_of(record: &BamRecord, header: &SamHeader) -> String {
    usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.clone())
        .unwrap_or_else(|| "null".to_string())
}
