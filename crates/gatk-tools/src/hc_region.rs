//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils`
//! (`assembleReads`, `getPaddedReferenceLoc`, `createReferenceHaplotype`) and HaplotypeCaller's
//! assembler and finalization settings (GATK 4.6.2.0): a region's reads assembled into haplotypes
//! the way `HaplotypeCallerEngine.callRegion` asks.
//!
//! The region is finalized with HaplotypeCaller's switches (soft clips reverted where the fragment
//! allows, tails under the minimum base quality less one clipped, overlapping mates corrected, the
//! hard-clipped pileup reads kept); the reference is the padded span widened by 500 on each side
//! within the contig; the reference haplotype is the padded span's bases, aligned at its offset
//! into that reference; and the assembler is `ReadThreadingAssemblerArgumentCollection`'s defaults
//! for HaplotypeCaller (k-mers 10 and 25, 128 paths, a prune factor of 2, dangling branches of 4
//! recovered, a minimum base quality of 10).

use gatk_engine::assembly_based_caller_utils::{finalize_region, FinalizeArguments, FinalizeError};
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::assembly_region_trimmer::AssemblyRegionTrimmer;
use gatk_engine::assembly_result_set::AssemblyResultSet;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::pair_hmm_likelihood_engine::PairHmmLikelihoodEngine;
use gatk_engine::read_pileup::sample_name;
use gatk_engine::read_realignment::realign_reads_to_their_best_haplotype;
use gatk_engine::read_threading_assembler::{
    AssemblerError, AssemblerSettings, ReadThreadingAssembler,
};
use gatk_engine::smith_waterman::{
    SmithWatermanAligner, ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, NEW_SW_PARAMETERS,
    STANDARD_NGS,
};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::genotyping_engine::GenotypingEngine;
use crate::hc_genotyping::{
    assign_genotype_likelihoods_full, CallAnnotator, HcGenotypingArguments,
};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::SamHeader;

/// `AssemblyBasedCallerUtils.REFERENCE_PADDING_FOR_ASSEMBLY`.
pub const REFERENCE_PADDING_FOR_ASSEMBLY: i32 = 500;

/// What the region assembly refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum RegionAssemblyError {
    Finalize(FinalizeError),
    Assembler(AssemblerError),
    IllegalState(String),
    /// A later stage of `callRegion`, with its message.
    Calling(String),
}

/// `HaplotypeCallerEngine.READ_LENGTH_FILTER_THRESHOLD`.
const READ_LENGTH_FILTER_THRESHOLD: i32 = 10;
/// `AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING`.
const MINIMUM_READ_LENGTH_AFTER_TRIMMING: i32 = 10;

/// What `callRegion` runs with, beside the region.
pub struct CallRegionContext<'a> {
    pub header: &'a SamHeader,
    pub samples: &'a [String],
    pub contig_bases: &'a [u8],
    pub assembly: &'a RegionAssemblyArguments,
    pub assembler: &'a ReadThreadingAssembler,
    pub aligner: &'a dyn SmithWatermanAligner,
    pub trimmer: &'a AssemblyRegionTrimmer,
    pub likelihood_engine: &'a PairHmmLikelihoodEngine,
    pub genotyping: &'a HcGenotypingArguments,
    /// `--mapping-quality-threshold-for-genotyping`, 20 by default.
    pub mapping_quality_threshold: u8,
}

/// `AlignmentUtils.unclippedReadLength`: the bases that are not soft-clipped.
fn unclipped_read_length(read: &BamRecord) -> i32 {
    let soft: u32 = read
        .cigar
        .elements
        .iter()
        .filter(|e| e.op == Op::S)
        .map(|e| e.length)
        .sum();
    read.read_bases.len() as i32 - soft as i32
}

/// `HaplotypeCallerEngine.callRegion` without reference confidence: the region's calls.
pub fn call_region(
    mut region: AssemblyRegion,
    context: &CallRegionContext<'_>,
    engine: &mut GenotypingEngine,
    annotator: Option<CallAnnotator<'_>>,
) -> Result<Vec<VariantContext>, RegionAssemblyError> {
    let calling = |e: String| RegionAssemblyError::Calling(e);
    // `referenceModelForNoVariation` answers nothing outside reference-confidence mode.
    if !region.is_active() || region.reads().is_empty() {
        return Ok(Vec::new());
    }
    let mut untrimmed = assemble_reads(
        &mut region,
        context.assembly,
        context.header,
        context.samples,
        context.contig_bases,
        context.assembler,
        context.aligner,
    )?;
    let events = untrimmed
        .variation_events(context.genotyping.max_mnp_distance)
        .map_err(|e| calling(format!("{e:?}")))?;
    let padded = region.padded_span().clone();
    let padded_bases = &context.contig_bases[padded.start as usize - 1..padded.end as usize];
    let trimmed = context
        .trimmer
        .trim(&region, &events, padded_bases, padded.start)
        .map_err(|e| calling(format!("{e:?}")))?;
    if trimmed.variant_span.is_none() {
        return Ok(Vec::new());
    }
    let variant_region = context
        .trimmer
        .variant_region(&trimmed, &region, context.header)
        .map_err(|e| calling(format!("{e:?}")))?;
    let mut assembly = untrimmed
        .trim_to(variant_region)
        .map_err(|e| calling(format!("{e:?}")))?;

    // The read stubs the trim left, then `filterNonPassingReads`, whose reads the annotations
    // still see.
    let genotyping_region = assembly
        .region_for_genotyping_mut()
        .expect("a trimmed set has its region");
    genotyping_region
        .retain_reads(|r| unclipped_read_length(r) >= MINIMUM_READ_LENGTH_AFTER_TRIMMING);
    let mut filtered: Vec<Vec<BamRecord>> = vec![Vec::new(); context.samples.len()];
    for read in genotyping_region.reads() {
        let failing = unclipped_read_length(read) < READ_LENGTH_FILTER_THRESHOLD
            || read.mapping_quality < context.mapping_quality_threshold
            || !gatk_readfilter::mate_on_same_contig_or_no_mapped_mate(read);
        if failing {
            let sample = sample_name(read, context.header);
            let index = context
                .samples
                .iter()
                .position(|s| Some(s) == sample.as_ref())
                .ok_or_else(|| calling("a filtered read of no known sample".to_string()))?;
            filtered[index].push(read.clone());
        }
    }
    let threshold = context.mapping_quality_threshold;
    genotyping_region.retain_reads(|read| {
        !(unclipped_read_length(read) < READ_LENGTH_FILTER_THRESHOLD
            || read.mapping_quality < threshold
            || !gatk_readfilter::mate_on_same_contig_or_no_mapped_mate(read))
    });
    if !assembly.is_variation_present() {
        return Ok(Vec::new());
    }
    let genotyping_region = assembly
        .region_for_genotyping()
        .expect("a trimmed set has its region")
        .clone();
    if genotyping_region.reads().is_empty() {
        return Ok(Vec::new());
    }

    // `splitReadsBySample`, then the likelihoods against the trimmed haplotypes.
    let mut reads: Vec<Vec<BamRecord>> = vec![Vec::new(); context.samples.len()];
    for read in genotyping_region.reads() {
        let sample = sample_name(read, context.header);
        let index = context
            .samples
            .iter()
            .position(|s| Some(s) == sample.as_ref())
            .ok_or_else(|| calling("a read of no known sample".to_string()))?;
        reads[index].push(read.clone());
    }
    let mut haplotypes: Vec<_> = assembly.haplotype_list().to_vec();
    let mut likelihoods = context
        .likelihood_engine
        .compute_read_likelihoods(&haplotypes, context.samples, &reads)
        .map_err(|e| calling(format!("{e:?}")))?;
    let ref_loc = assembly
        .padded_reference_loc()
        .expect("a trimmed set has its padded reference")
        .clone();
    let reference_haplotype = assembly
        .reference_haplotype()
        .expect("a trimmed set has its reference")
        .clone();
    realign_reads_to_their_best_haplotype(
        &mut likelihoods,
        &reference_haplotype,
        ref_loc.start,
        context.aligner,
        &ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
    )
    .map_err(|e| calling(format!("{e:?}")))?;
    let contig_length = context.contig_bases.len() as i32;
    assign_genotype_likelihoods_full(
        engine,
        context.genotyping,
        &mut haplotypes,
        &likelihoods,
        context.samples,
        assembly.full_reference_with_padding(),
        &ref_loc,
        genotyping_region.span(),
        contig_length,
        annotator,
        &filtered,
    )
    .map_err(|e| calling(format!("{e:?}")))
}

/// HaplotypeCaller's `minBaseQualityScore` and the switches `assembleReads` passes on.
#[derive(Debug, Clone)]
pub struct RegionAssemblyArguments {
    pub min_base_quality_score: u8,
    pub dont_use_soft_clipped_bases: bool,
    pub soft_clip_low_quality_ends: bool,
    pub override_softclip_fragment_check: bool,
    pub do_not_correct_overlapping_base_qualities: bool,
}

impl Default for RegionAssemblyArguments {
    fn default() -> Self {
        RegionAssemblyArguments {
            min_base_quality_score: 10,
            dont_use_soft_clipped_bases: false,
            soft_clip_low_quality_ends: false,
            override_softclip_fragment_check: false,
            do_not_correct_overlapping_base_qualities: false,
        }
    }
}

/// `hcArgs.createReadThreadingAssembler()` with HaplotypeCaller's defaults.
pub fn haplotype_caller_assembler(min_base_quality_score: u8) -> ReadThreadingAssembler {
    let settings = AssemblerSettings {
        max_allowed_paths: 128,
        kmer_sizes: vec![10, 25],
        dont_increase_kmer_sizes_for_cycles: false,
        allow_non_unique_kmers_in_ref: false,
        num_pruning_samples: 1,
        prune_factor: 2,
        use_adaptive_pruning: false,
        initial_error_rate_for_pruning: 0.001,
        pruning_log_odds_threshold: gatk_engine::allele_likelihoods::log10_to_log(1.0),
        pruning_seeding_log_odds_threshold: gatk_engine::allele_likelihoods::log10_to_log(4.0),
        max_unpruned_variants: 100,
        enable_legacy_graph_cycle_detection: false,
        min_matching_bases_to_dangling_end_recovery: -1,
    };
    let mut assembler = ReadThreadingAssembler::new(&settings).expect("valid defaults");
    assembler.set_recover_dangling_branches(true);
    assembler.set_recover_all_dangling_branches(false);
    assembler.set_min_dangling_branch_length(4);
    assembler.set_min_base_quality_to_use_in_assembly(min_base_quality_score);
    assembler
}

/// `getPaddedReferenceLoc(region, padding, reader)`: the padded span widened within the contig.
pub fn padded_reference_loc(
    region: &AssemblyRegion,
    padding: i32,
    contig_length: i32,
) -> SimpleInterval {
    let span = region.padded_span();
    SimpleInterval::new(
        &span.contig,
        (span.start - padding).max(1),
        (span.end + padding).min(contig_length),
    )
    .expect("a padded span inside the contig")
}

/// `assembleReads(region, ...)`: the region finalized in place, then assembled against
/// `contig_bases`, the whole of the region's contig.
pub fn assemble_reads(
    region: &mut AssemblyRegion,
    arguments: &RegionAssemblyArguments,
    header: &SamHeader,
    samples: &[String],
    contig_bases: &[u8],
    assembler: &ReadThreadingAssembler,
    aligner: &dyn SmithWatermanAligner,
) -> Result<AssemblyResultSet, RegionAssemblyError> {
    finalize_region(
        region,
        &FinalizeArguments {
            error_correct_reads: false,
            dont_use_soft_clipped_bases: arguments.dont_use_soft_clipped_bases,
            min_tail_quality: arguments.min_base_quality_score.wrapping_sub(1),
            correct_overlapping_base_qualities: !arguments
                .do_not_correct_overlapping_base_qualities,
            soft_clip_low_quality_ends: arguments.soft_clip_low_quality_ends,
            override_softclip_fragment_check: arguments.override_softclip_fragment_check,
            track_hardclipped_reads: true,
        },
        header,
        samples,
    )
    .map_err(RegionAssemblyError::Finalize)?;
    let contig_length = contig_bases.len() as i32;
    let ref_loc = padded_reference_loc(region, REFERENCE_PADDING_FOR_ASSEMBLY, contig_length);
    let full_reference = contig_bases[ref_loc.start as usize - 1..ref_loc.end as usize].to_vec();
    let ref_haplotype = reference_haplotype(region, contig_bases, &ref_loc)?;
    assembler
        .run_local_assembly(
            region,
            ref_haplotype,
            &full_reference,
            &ref_loc,
            header,
            aligner,
            &STANDARD_NGS,
            &NEW_SW_PARAMETERS,
        )
        .map_err(RegionAssemblyError::Assembler)
}

/// `createReferenceHaplotype(region, paddedReferenceLoc, reader)`: the padded span's bases, a
/// reference haplotype aligned at their offset into the padded reference.
pub fn reference_haplotype(
    region: &AssemblyRegion,
    contig_bases: &[u8],
    ref_loc: &SimpleInterval,
) -> Result<Haplotype, RegionAssemblyError> {
    let padded = region.padded_span();
    let alignment_start = padded.start - ref_loc.start;
    if alignment_start < 0 {
        return Err(RegionAssemblyError::IllegalState(format!(
            "Bad alignment start in createReferenceHaplotype {alignment_start}"
        )));
    }
    let bases = &contig_bases[padded.start as usize - 1..padded.end as usize];
    let mut haplotype = Haplotype::new(bases, true)
        .map_err(|e| RegionAssemblyError::IllegalState(e.to_string()))?;
    haplotype.set_genome_location(padded.clone());
    haplotype.set_alignment_start_hap_wrt_ref(alignment_start);
    haplotype
        .set_cigar(&Cigar::new(vec![CigarElement {
            length: bases.len() as u32,
            op: Op::M,
        }]))
        .map_err(|e| RegionAssemblyError::IllegalState(format!("{e:?}")))?;
    Ok(haplotype)
}
