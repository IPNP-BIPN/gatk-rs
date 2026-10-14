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

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::assembly_based_caller_utils::{finalize_region, FinalizeArguments, FinalizeError};
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::assembly_region_trimmer::{AssemblyRegionTrimmer, TrimResult};
use gatk_engine::assembly_result_set::{AssemblyResultSet, ResultSetError};
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::read_pileup::sample_name;
use gatk_engine::read_realignment::realign_reads_to_their_best_haplotype;
use gatk_engine::read_threading_assembler::{
    AssemblerError, AssemblerSettings, ReadThreadingAssembler,
};
use gatk_engine::smith_waterman::{
    SmithWatermanAligner, SwParameters, ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
    NEW_SW_PARAMETERS, STANDARD_NGS,
};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::ramp_files::{OnRamp, OnRampType};
use crate::reference_confidence_model::calculate_ref_confidence;

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
    /// `AssemblyRegionTrimmer.Result.getVariantRegion` on a region with no variation, which
    /// `--disable-optimizations` reaches: the reference's `IllegalStateException`.
    NoVariationPresent,
    /// A later stage of `callRegion`, with its message.
    Calling(String),
}

/// `HaplotypeCallerEngine.READ_LENGTH_FILTER_THRESHOLD`.
const READ_LENGTH_FILTER_THRESHOLD: i32 = 10;
/// `AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING`.
const MINIMUM_READ_LENGTH_AFTER_TRIMMING: i32 = 10;

/// `--assembly-region-out`'s header: `IGVUtils.printIGVFormatHeader(out, "line",
/// "AssemblyRegions")`.
pub const ASSEMBLY_REGION_OUT_HEADER: &str =
    "#track graphType=line\nChromosome\tStart\tEnd\tFeature\tAssemblyRegions\n";

/// `AssemblyRegionWalker.writeAssemblyRegion`: an end marker at the region's start, then the
/// region with its size, 1 when active and -1 when not, each `printIGVFormatRow`'s zero-based start
/// and `%.5f` values.
pub fn assembly_region_rows(region: &AssemblyRegion) -> String {
    let span = region.span();
    let value = if region.is_active() { 1.0 } else { -1.0 };
    format!(
        "{c}\t{s0}\t{s}\tend-marker\t{zero:.5}\n{c}\t{s0}\t{e}\tsize={size}\t{value:.5}\n",
        c = span.contig,
        s0 = span.start - 1,
        s = span.start,
        e = span.end,
        zero = 0.0,
        size = span.end - span.start + 1,
    )
}

/// What `callRegion` runs with, beside the region.
pub struct CallRegionContext<'a> {
    pub header: &'a SamHeader,
    pub samples: &'a [String],
    pub contig_bases: &'a [u8],
    pub assembly: &'a RegionAssemblyArguments,
    pub assembler: &'a ReadThreadingAssembler,
    pub aligner: &'a dyn SmithWatermanAligner,
    pub trimmer: &'a AssemblyRegionTrimmer,
    pub likelihood_engine: &'a crate::flow_likelihood_engine::ReadLikelihoodEngine,
    /// `filterStepLikelihoodCalculationEngine`, with stepwise filtering.
    pub stepwise_engine: Option<&'a crate::flow_likelihood_engine::ReadLikelihoodEngine>,
    pub genotyping: &'a HcGenotypingArguments,
    /// `--mapping-quality-threshold-for-genotyping`, 20 by default.
    pub mapping_quality_threshold: u8,
    /// RampedHaplotypeCaller's off ramp: the region's state is written to it at the ramp's point,
    /// and the region answers no calls.
    pub off_ramp: Option<&'a std::cell::RefCell<crate::ramp_files::OffRamp>>,
    /// RampedHaplotypeCaller's on ramp: the assembly or the likelihoods are read back from it
    /// instead of computed.
    pub on_ramp: Option<&'a crate::ramp_files::OnRamp>,
    /// `--disable-optimizations`: a region with no variation or no reads left is genotyped
    /// anyway. One the trimmer found no variation in then fails, as the reference's
    /// `getVariantRegion` does.
    pub disable_optimizations: bool,
    /// The haplotype collapsing settings.
    pub collapse: CollapseSettings,
    /// `--ramps-debug-post-assembler-on`: the post-assembler on ramp's assembly is checked against
    /// a fresh one, which then stands in its place.
    pub ramps_debug_post_assembler_on: bool,
    /// Pileup detection: its arguments and the loci walked, `region.getAlignmentData()` being
    /// those of a region's active span.
    pub pileup: Option<(
        &'a crate::pileup_based_alleles::PileupDetectionArguments,
        &'a [crate::pileup_based_alleles::LocusSummary],
    )>,
    /// `--use-pdhmm`: partially determined haplotypes and the engine that scores them.
    pub pd: Option<PdContext<'a>>,
    /// `--alleles`: the records, and `--force-call-filtered-alleles`.
    pub given_alleles: Option<(&'a [htsjdk_vcf::variant::VariantContext], bool)>,
    /// `--bam-output`.
    pub bamout: Option<&'a std::cell::RefCell<crate::haplotype_bam_writer::HaplotypeBamWriter>>,
}

/// What `--use-pdhmm` runs with.
pub struct PdContext<'a> {
    pub arguments: crate::pd_haplotypes::PdArguments,
    /// `--fallback-gga-if-pdhmm-fails`.
    pub gga_fallback: bool,
    /// `informativeReadOverlapMargin` with `--use-pdhmm-overlap-optimization`, else -1.
    pub overlap_range: i32,
    /// `PDPairHMMLikelihoodCalculationEngine`, built from the likelihood arguments.
    pub engine: &'a gatk_engine::pair_hmm_likelihood_engine::PairHmmLikelihoodEngine,
    /// `--pdhmm-results-file`'s text.
    pub results: Option<&'a std::cell::RefCell<String>>,
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
    if !region.is_active() {
        return reference_model_for_no_variation(region, true, context);
    }
    // The given events: the records over the padded span, split; HaplotypeCaller keeps those
    // starting in the region, the ramped engine every one.
    let given: Vec<gatk_engine::event_map::Event> = match context.given_alleles {
        Some((records, force_call_filtered)) => crate::given_alleles::given_events(
            records,
            &region.span().contig,
            region.padded_span(),
            if context.genotyping.ramped_engine {
                i32::MIN
            } else {
                region.span().start
            },
            force_call_filtered,
        )
        .map_err(calling)?,
        None => Vec::new(),
    };
    if given.is_empty() && region.reads().is_empty() {
        return reference_model_for_no_variation(region, true, context);
    }
    // The pre-assembler off ramp: the region's place, before anything is assembled.
    let ramp_kind = context.off_ramp.map(|ramp| ramp.borrow().kind);
    if ramp_kind == Some(crate::ramp_files::OffRampType::PreAssembler) {
        if let Some(ramp) = context.off_ramp {
            ramp.borrow_mut().add_assembly(
                region.span(),
                &crate::ramp_files::region_description(&region),
                "assemblyResult",
                None,
                None,
                None,
                context.header,
            );
        }
        return Ok(Vec::new());
    }
    let on_kind = context.on_ramp.map(|ramp| ramp.kind);
    // The post-assembler on ramp's flanks, which stand in for the trimmer's.
    let mut ramp_flanks: Option<(Option<AssemblyRegion>, Option<AssemblyRegion>)> = None;
    // With the debug switch, the ramp's assembly is only what the fresh one is checked against.
    let mut debug_assembly: Option<AssemblyResultSet> = None;
    let from_ramp = if let (Some(OnRampType::PostAssembler), Some(ramp)) =
        (on_kind, context.on_ramp)
    {
        // The assembly read back: its haplotypes, its padded reference and its flanks.
        let flank = |key: &str| {
            ramp.region(region.span(), key)
                .map(OnRamp::assembly_region)
                .transpose()
                .map_err(calling)
        };
        let left = flank("assemblyResult.nonVariantLeftFlankRegion")?;
        let right = flank("assemblyResult.nonVariantRightFlankRegion")?;
        let Some(object) = ramp.region(region.span(), "assemblyResult") else {
            return Ok(Vec::new());
        };
        let haplotypes = ramp.bam_haplotypes(region.span()).map_err(calling)?;
        let padded_loc = ramp
            .region_interval(region.span(), "assemblyResult", "paddedReferenceLoc")
            .map_err(calling)?;
        let genotyping =
            OnRamp::assembly_region(object.get("regionForGenotyping").ok_or_else(|| {
                calling("JSONObject[\"regionForGenotyping\"] not found.".to_string())
            })?)
            .map_err(calling)?;
        let mut set = AssemblyResultSet::new();
        set.set_padded_reference_loc(padded_loc.clone());
        set.set_region_for_genotyping(genotyping);
        for haplotype in haplotypes {
            set.add(haplotype).map_err(|e| calling(format!("{e:?}")))?;
        }
        set.set_full_reference_with_padding(
            context.contig_bases[padded_loc.start as usize - 1..padded_loc.end as usize].to_vec(),
        );
        // The reads are still assembled, for the region they leave to genotype.
        let untrimmed = assemble_reads(
            &mut region,
            context.assembly,
            context.header,
            context.samples,
            context.contig_bases,
            context.assembler,
            context.aligner,
        )?;
        set.set_region_for_genotyping(
            untrimmed
                .region_for_genotyping()
                .expect("an assembled set has its region")
                .clone(),
        );
        let events = set
            .variation_events(context.genotyping.max_mnp_distance)
            .map_err(|e| calling(format!("{e:?}")))?;
        let padded = region.padded_span().clone();
        let padded_bases = &context.contig_bases[padded.start as usize - 1..padded.end as usize];
        let trimmed = context
            .trimmer
            .trim(&region, &events, padded_bases, padded.start)
            .map_err(|e| calling(format!("{e:?}")))?;
        // `getVariantRegion()` validates rather than answers: no variation is an error here.
        if trimmed.variant_span.is_none() {
            return Err(calling(
                "java.lang.IllegalArgumentException: There is no variation present.".to_string(),
            ));
        }
        let variant_region = context
            .trimmer
            .variant_region(&trimmed, &region, context.header)
            .map_err(|e| calling(format!("{e:?}")))?;
        // "assembler bailing out": a trim the set refuses is a region with no calls.
        let assembly = match set.trim_to(variant_region) {
            Ok(assembly) => assembly,
            Err(ResultSetError::IllegalState(_)) => return Ok(Vec::new()),
            Err(e) => return Err(calling(format!("{e:?}"))),
        };
        if context.ramps_debug_post_assembler_on {
            debug_assembly = Some(assembly);
            None
        } else {
            ramp_flanks = Some((left, right));
            Some((assembly, trimmed))
        }
    } else {
        None
    };
    // The good pileup events, injected once the set is trimmed.
    let mut good_pileup_events: Vec<gatk_engine::event_map::Event> = Vec::new();
    let mut bad_pileup_events: Vec<gatk_engine::event_map::Event> = Vec::new();
    let (mut assembly, trimmed) = if let Some(assembled) = from_ramp {
        assembled
    } else {
        let post_filter = on_kind == Some(OnRampType::PostFilter);
        let mut untrimmed = if post_filter {
            // `bypassAssembly`: the region is finalized, and the set holds the reference alone.
            assemble_bypassed(&mut region, context)?
        } else {
            assemble_reads(
                &mut region,
                context.assembly,
                context.header,
                context.samples,
                context.contig_bases,
                context.assembler,
                context.aligner,
            )?
        };
        if let (true, Some(ramp)) = (post_filter, context.on_ramp) {
            if !ramp.has_region(region.span()) {
                return Ok(Vec::new());
            }
            for haplotype in ramp
                .haplotypes(region.span(), "readLikelihoods.haplotypes")
                .map_err(calling)?
            {
                if !haplotype.is_reference() {
                    untrimmed
                        .add(haplotype)
                        .map_err(|e| calling(format!("{e:?}")))?;
                }
            }
        }
        let mut events = untrimmed
            .variation_events(context.genotyping.max_mnp_distance)
            .map_err(|e| calling(format!("{e:?}")))?;
        // `goodAndBadPileupEvents`: the bad ones leave the trimmer's events, the good ones join.
        if let Some((detection, loci)) = context.pileup {
            let (good, bad) = crate::pileup_based_alleles::good_and_bad_pileup_events(
                crate::pileup_based_alleles::region_loci(loci, &region),
                detection,
                context.contig_bases,
            )
            .map_err(calling)?;
            events = crate::pileup_based_alleles::events_for_trimming(&events, &good, &bad);
            good_pileup_events = good;
            bad_pileup_events = bad;
        }
        // `givenAlleles.forEach(allVariationEvents::add)`, which the ramped engine does not do.
        if !given.is_empty() && !context.genotyping.ramped_engine {
            events = crate::pileup_based_alleles::events_for_trimming(&events, &given, &[]);
        }
        let padded = region.padded_span().clone();
        let padded_bases = &context.contig_bases[padded.start as usize - 1..padded.end as usize];
        let trimmed = context
            .trimmer
            .trim(&region, &events, padded_bases, padded.start)
            .map_err(|e| calling(format!("{e:?}")))?;
        let variant_region = match (post_filter, context.on_ramp) {
            // The genotyping region the likelihoods were computed over, from the ramp.
            (true, Some(ramp)) => {
                let span = ramp
                    .region_interval(region.span(), "readLikelihoods", "r4g")
                    .map_err(calling)?;
                let padded_span = ramp
                    .region_interval(region.span(), "readLikelihoods", "r4gp")
                    .map_err(calling)?;
                region
                    .trim(&span, &padded_span, context.header)
                    .map_err(|e| calling(format!("{e:?}")))?
            }
            _ => {
                if trimmed.variant_span.is_none() {
                    if !context.disable_optimizations {
                        return reference_model_for_no_variation(region, false, context);
                    }
                    return Err(RegionAssemblyError::NoVariationPresent);
                }
                context
                    .trimmer
                    .variant_region(&trimmed, &region, context.header)
                    .map_err(|e| calling(format!("{e:?}")))?
            }
        };
        let assembly = untrimmed
            .trim_to(variant_region)
            .map_err(|e| calling(format!("{e:?}")))?;
        (assembly, trimmed)
    };
    // `addGivenAlleles`, HaplotypeCaller's only.
    if !context.genotyping.ramped_engine {
        crate::given_alleles::add_given_alleles(
            &mut assembly,
            &given,
            context.genotyping.max_mnp_distance,
            context.aligner,
            &context.assembly.haplotype_to_reference_sw,
        )
        .map_err(calling)?;
    }
    // `generatePDHaplotypes`: the assembled haplotypes replaced by partially determined ones,
    // unless the construction gave up. The reference is kept aside: the set no longer has one.
    let assembled_reference = assembly.reference_haplotype().cloned();
    if let Some(pd) = &context.pd {
        crate::pd_haplotypes::generate_pd_haplotypes(
            &mut assembly,
            &bad_pileup_events,
            &good_pileup_events,
            context.aligner,
            &context.assembly.haplotype_to_reference_sw,
            &pd.arguments,
        )
        .map_err(calling)?;
    }
    // The pileup events go in the GGA way without PD haplotypes, or when they gave up and the
    // fallback is on: the haplotypes carrying a bad event out (only a set the PD engine built
    // event maps for has any to match), then the good events in.
    let fallback = context
        .pd
        .as_ref()
        .is_none_or(|pd| pd.gga_fallback && !assembly.has_overwritten_haps());
    if let (true, Some((detection, _))) = (fallback, context.pileup) {
        crate::pileup_based_alleles::remove_haplotypes_with_bad_alleles(
            &mut assembly,
            &bad_pileup_events,
        );
        crate::pileup_based_alleles::inject_pileup_events(
            &mut assembly,
            &region,
            detection,
            context.aligner,
            &context.assembly.haplotype_to_reference_sw,
            &good_pileup_events,
            context.genotyping.max_mnp_distance,
        )
        .map_err(calling)?;
    }
    if let Some(debug) = &debug_assembly {
        crate::ramp_files::compare_haplotypes(assembly.haplotype_list(), debug.haplotype_list())
            .map_err(calling)?;
        let reads = |set: &AssemblyResultSet| {
            set.region_for_genotyping()
                .expect("a trimmed set has its region")
                .reads()
                .to_vec()
        };
        crate::ramp_files::compare_reads(&reads(&assembly), &reads(debug), context.header)
            .map_err(calling)?;
    }
    // The post-assembler off ramp: the trimmed assembly, its spans and the flanks.
    if ramp_kind == Some(crate::ramp_files::OffRampType::PostAssembler) {
        let left = context
            .trimmer
            .non_variant_left_flank_region(&trimmed, &region, context.header)
            .map_err(|e| calling(format!("{e:?}")))?;
        let right = context
            .trimmer
            .non_variant_right_flank_region(&trimmed, &region, context.header)
            .map_err(|e| calling(format!("{e:?}")))?;
        let padded_loc = assembly
            .padded_reference_loc()
            .expect("a trimmed set has its padded reference")
            .clone();
        let genotyping = assembly
            .region_for_genotyping()
            .expect("a trimmed set has its region")
            .clone();
        if let Some(ramp) = context.off_ramp {
            ramp.borrow_mut().add_assembly(
                region.span(),
                &crate::ramp_files::region_description(&region),
                "assemblyResult",
                Some((&padded_loc, &genotyping, assembly.haplotype_list())),
                left.as_ref(),
                right.as_ref(),
                context.header,
            );
        }
        return Ok(Vec::new());
    }

    // The read stubs the trim left, then `filterNonPassingReads`, whose reads the annotations
    // still see.
    let genotyping_region = assembly
        .region_for_genotyping_mut()
        .expect("a trimmed set has its region");
    // The ramped engine measures a stub by its whole length, soft clips included.
    let ramped = context.genotyping.ramped_engine;
    if ramped {
        genotyping_region
            .retain_reads(|r| r.read_bases.len() as i32 >= MINIMUM_READ_LENGTH_AFTER_TRIMMING);
    } else {
        genotyping_region
            .retain_reads(|r| unclipped_read_length(r) >= MINIMUM_READ_LENGTH_AFTER_TRIMMING);
    }
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
    // HaplotypeCaller drops them before the likelihoods; the ramped engine only after, in
    // `genotype()`, so they are scored and realigned, then added again for the annotations.
    let threshold = context.mapping_quality_threshold;
    if !ramped {
        genotyping_region.retain_reads(|read| {
            !(unclipped_read_length(read) < READ_LENGTH_FILTER_THRESHOLD
                || read.mapping_quality < threshold
                || !gatk_readfilter::mate_on_same_contig_or_no_mapped_mate(read))
        });
    }
    if !assembly.is_variation_present() && !context.disable_optimizations {
        return reference_model_for_no_variation(region, false, context);
    }
    let genotyping_region = assembly
        .region_for_genotyping()
        .expect("a trimmed set has its region")
        .clone();
    if genotyping_region.reads().is_empty() && !context.disable_optimizations {
        return reference_model_for_no_variation(region, false, context);
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
    let ref_loc = assembly
        .padded_reference_loc()
        .expect("a trimmed set has its padded reference")
        .clone();
    let reference_haplotype = assembly
        .reference_haplotype()
        .or(assembled_reference.as_ref())
        .expect("a trimmed set has its reference")
        .clone();
    // The post-filter on ramp: the likelihoods are the file's, over the file's haplotypes and the
    // reads it names. What the engine would have computed in their place is never read.
    if let (Some(OnRampType::PostFilter), Some(ramp)) = (on_kind, context.on_ramp) {
        let mut all_reads: Vec<Vec<BamRecord>> = vec![Vec::new(); context.samples.len()];
        for read in region.reads() {
            let sample = sample_name(read, context.header);
            if let Some(index) = context
                .samples
                .iter()
                .position(|s| Some(s) == sample.as_ref())
            {
                all_reads[index].push(read.clone());
            }
        }
        let likelihoods = ramp
            .allele_likelihoods(
                region.span(),
                "readLikelihoods",
                context.samples,
                &reads,
                &all_reads,
            )
            .map_err(calling)?;
        haplotypes = (0..likelihoods.number_of_alleles())
            .filter_map(|index| likelihoods.get_allele(index).cloned())
            .collect();
        return finish_region(
            region,
            context,
            engine,
            annotator,
            &mut haplotypes,
            &likelihoods,
            &assembly,
            &ref_loc,
            &reference_haplotype,
            &genotyping_region,
            &filtered,
            &trimmed,
            ramp_flanks,
            &std::collections::BTreeSet::new(),
            None,
            &given,
        );
    }
    // Stepwise filtering scores with the flow aligner first; the ramped engine has no filtering
    // engine and scores with the main one, rescoring after the filter all the same.
    let first_engine = match (ramped, context.stepwise_engine) {
        (false, Some(engine)) => engine,
        _ => context.likelihood_engine,
    };
    let mut likelihoods = match (&context.pd, context.stepwise_engine) {
        // The PDHMM, over a set of partially determined haplotypes and without stepwise filtering.
        (Some(pd), None) if assembly.is_partially_determined_list() => {
            let mut results = pd.results.map(|cell| cell.borrow_mut());
            pd.engine
                .compute_read_likelihoods_partially_determined(
                    &haplotypes,
                    context.samples,
                    &reads,
                    pd.overlap_range,
                    results.as_deref_mut(),
                )
                .map_err(|e| match e {
                    gatk_engine::pair_hmm_likelihood_engine::LikelihoodEngineError::Thrown(m) => {
                        calling(m)
                    }
                    other => calling(format!("{other:?}")),
                })?
        }
        _ => first_engine
            .compute_read_likelihoods(&haplotypes, context.header, context.samples, &reads)
            .map_err(calling)?,
    };
    // Stepwise filtering: the main engine scores the reads again over the haplotypes the
    // filter kept, in their order.
    let rescore = |likelihoods: &mut AlleleLikelihoods<BamRecord, Haplotype>,
                   haplotypes: &[Haplotype]|
     -> Result<(), RegionAssemblyError> {
        if context.stepwise_engine.is_some() {
            *likelihoods = context
                .likelihood_engine
                .compute_read_likelihoods(haplotypes, context.header, context.samples, &reads)
                .map_err(calling)?;
        }
        Ok(())
    };
    // The collapsing engine `assembleReads` built for the region: a threshold above zero and a
    // run longer than it in the reference under the padded span. Read back from a post-assembler
    // on ramp, the set has none, and the ramped engine builds one only for a threshold given.
    let threshold = if on_kind == Some(OnRampType::PostAssembler) {
        context.collapse.raw
    } else {
        context.collapse.determined
    };
    let padded_span = region.padded_span().clone();
    let collapsing = (threshold > 0
        && crate::haplotype_collapsing::needs_collapsing(
            &context.contig_bases[padded_span.start as usize - 1..padded_span.end as usize],
            threshold,
        ))
    .then(|| crate::haplotype_collapsing::HaplotypeCollapsing {
        hmer_size_threshold: threshold,
        partial_mode: context.collapse.partial,
        full_ref: assembly.full_reference_with_padding().to_vec(),
        ref_loc: ref_loc.clone(),
        parameters: context.assembly.haplotype_to_reference_sw,
    });
    if let (false, Some(engine)) = (ramped, &collapsing) {
        uncollapse(engine, &mut likelihoods, &mut haplotypes, context.aligner).map_err(calling)?;
    }
    // `filterAlleles`: HaplotypeCaller's engine before realigning, the ramped one after its
    // pre-filter off ramp.
    let mut suspicious = std::collections::BTreeSet::new();
    // `preFilteringAlleleLikelihoods`: what the jumbo annotations see, its haplotypes given their
    // event maps and its count of filtered haplotypes set.
    let mut pre_filtering: Option<(AlleleLikelihoods<BamRecord, Haplotype>, Vec<Haplotype>)> = None;
    let filter =
        |likelihoods: &mut AlleleLikelihoods<BamRecord, Haplotype>,
         haplotypes: &mut Vec<Haplotype>,
         suspicious: &mut std::collections::BTreeSet<i32>,
         pre_filtering: &mut Option<(AlleleLikelihoods<BamRecord, Haplotype>, Vec<Haplotype>)>|
         -> Result<(), RegionAssemblyError> {
            if let Some(arguments) = &context.genotyping.filter_alleles {
                let filtered = crate::allele_filtering::filter_alleles(
                    likelihoods,
                    assembly.full_reference_with_padding(),
                    &ref_loc,
                    context.genotyping.max_mnp_distance,
                    arguments,
                    crate::allele_filtering::FilterModel::Germline {
                        ploidy: context.genotyping.sample_ploidy,
                    },
                )
                .map_err(calling)?;
                let mut before = likelihoods.clone();
                before.set_filtered_haplotype_count(filtered.filtered_haplotype_count);
                let mut before_haplotypes: Vec<Haplotype> = (0..before.number_of_alleles())
                    .filter_map(|index| before.get_allele(index).cloned())
                    .collect();
                gatk_engine::event_map::build_event_maps_for_haplotypes(
                    &mut before_haplotypes,
                    assembly.full_reference_with_padding(),
                    &ref_loc,
                    context.genotyping.max_mnp_distance,
                )
                .map_err(|e| calling(e.message()))?;
                *pre_filtering = Some((before, before_haplotypes));
                let kept = filtered.likelihoods;
                *haplotypes = (0..kept.number_of_alleles())
                    .filter_map(|index| kept.get_allele(index).cloned())
                    .collect();
                *likelihoods = kept;
                *suspicious = filtered.suspicious_locations;
            }
            Ok(())
        };
    if !ramped {
        filter(
            &mut likelihoods,
            &mut haplotypes,
            &mut suspicious,
            &mut pre_filtering,
        )?;
        rescore(&mut likelihoods, &haplotypes)?;
    }
    // No realignment with PD haplotypes, whose bases no longer say what they carry, whether or
    // not the construction gave up.
    let realign = context
        .pd
        .as_ref()
        .is_none_or(|pd| pd.arguments.use_determined_haplotypes);
    if realign {
        realign_reads_to_their_best_haplotype(
            &mut likelihoods,
            &reference_haplotype,
            ref_loc.start,
            context.aligner,
            // The ramped engine realigns with the haplotype-to-reference weights.
            if ramped {
                &context.assembly.haplotype_to_reference_sw
            } else {
                &context.assembly.read_to_haplotype_sw
            },
        )
        .map_err(|e| calling(format!("{e:?}")))?;
    }
    if let (true, Some(engine)) = (ramped, &collapsing) {
        uncollapse(engine, &mut likelihoods, &mut haplotypes, context.aligner).map_err(calling)?;
    }
    // The pre-filter off ramp: the realigned likelihoods, before any allele is filtered.
    if ramp_kind == Some(crate::ramp_files::OffRampType::PreFilter) {
        if let Some(ramp) = context.off_ramp {
            ramp.borrow_mut().add_likelihoods(
                region.span(),
                "readLikelihoods",
                &likelihoods,
                &genotyping_region,
            );
        }
        return Ok(Vec::new());
    }
    if ramped {
        filter(
            &mut likelihoods,
            &mut haplotypes,
            &mut suspicious,
            &mut pre_filtering,
        )?;
        rescore(&mut likelihoods, &haplotypes)?;
    }
    finish_region(
        region,
        context,
        engine,
        annotator,
        &mut haplotypes,
        &likelihoods,
        &assembly,
        &ref_loc,
        &reference_haplotype,
        &genotyping_region,
        &filtered,
        &trimmed,
        ramp_flanks,
        &suspicious,
        pre_filtering.as_ref().map(|(l, h)| (l, h.as_slice())),
        &given,
    )
}

/// `--flow-assembly-collapse-hmer-size` as given (`raw`) and as `assembleReads` resolves -1
/// (`determined`), and `--flow-assembly-collapse-partial-mode`.
#[derive(Debug, Clone, Copy, Default)]
pub struct CollapseSettings {
    pub raw: i32,
    pub determined: i32,
    pub partial: bool,
}

/// `possiblyUncollapseHaplotypesInReadLikelihoods`: the haplotypes put back to the reference's
/// hmer lengths and those made identical merged.
pub(crate) fn uncollapse(
    engine: &crate::haplotype_collapsing::HaplotypeCollapsing,
    likelihoods: &mut AlleleLikelihoods<BamRecord, Haplotype>,
    haplotypes: &mut Vec<Haplotype>,
    aligner: &dyn gatk_engine::smith_waterman::SmithWatermanAligner,
) -> Result<(), String> {
    let original: Vec<Haplotype> = (0..likelihoods.number_of_alleles())
        .filter_map(|index| likelihoods.get_allele(index).cloned())
        .collect();
    let uncollapsed = engine.uncollapse_hmers_in_haplotypes(&original, aligner)?;
    let groups = crate::haplotype_collapsing::identical_by_sequence(&uncollapsed)?;
    let new_to_old: Vec<(Haplotype, Vec<Haplotype>)> = groups
        .iter()
        .map(|(key, members)| {
            (
                uncollapsed[*key].clone(),
                members.iter().map(|&m| original[m].clone()).collect(),
            )
        })
        .collect();
    *likelihoods = likelihoods
        .marginalize(&new_to_old)
        .map_err(|e| format!("{e:?}"))?;
    *haplotypes = new_to_old.into_iter().map(|(h, _)| h).collect();
    Ok(())
}

/// The genotyping and, in reference-confidence mode, the flanks and the per-base records: what
/// `callRegion` does once the likelihoods are in hand, whether computed or read from an on ramp.
#[allow(clippy::too_many_arguments)]
fn finish_region(
    region: AssemblyRegion,
    context: &CallRegionContext<'_>,
    engine: &mut GenotypingEngine,
    annotator: Option<CallAnnotator<'_>>,
    haplotypes: &mut [Haplotype],
    likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>,
    assembly: &AssemblyResultSet,
    ref_loc: &SimpleInterval,
    reference_haplotype: &Haplotype,
    genotyping_region: &AssemblyRegion,
    filtered: &[Vec<BamRecord>],
    trimmed: &TrimResult,
    ramp_flanks: Option<(Option<AssemblyRegion>, Option<AssemblyRegion>)>,
    suspicious: &std::collections::BTreeSet<i32>,
    pre_filtering: Option<(&AlleleLikelihoods<BamRecord, Haplotype>, &[Haplotype])>,
    given: &[gatk_engine::event_map::Event],
) -> Result<Vec<VariantContext>, RegionAssemblyError> {
    let calling = |e: String| RegionAssemblyError::Calling(e);
    let contig_length = context.contig_bases.len() as i32;
    let mut bamout_tags = context
        .bamout
        .map(|_| crate::hc_genotyping::BamoutTags::default());
    let calls = assign_genotype_likelihoods_full(
        engine,
        context.genotyping,
        haplotypes,
        likelihoods,
        context.samples,
        assembly.full_reference_with_padding(),
        ref_loc,
        genotyping_region.span(),
        contig_length,
        annotator,
        filtered,
        suspicious,
        pre_filtering,
        given,
        bamout_tags.as_mut(),
    )
    .map_err(|e| calling(format!("{e:?}")))?;
    // `writeReadsAlignedToHaplotypes`, with the called haplotypes and the read tags genotyping
    // left.
    if let (Some(writer), Some(tags)) = (context.bamout, &bamout_tags) {
        let called: Vec<Haplotype> = tags.called.iter().map(|&h| haplotypes[h].clone()).collect();
        let reference = if context.disable_optimizations {
            assembly.reference_haplotype()
        } else {
            None
        };
        writer
            .borrow_mut()
            .write_region(
                haplotypes,
                ref_loc,
                &called,
                reference,
                likelihoods,
                genotyping_region.span(),
                &|read| tags.supported_for(read).map(str::to_string),
            )
            .map_err(calling)?;
    }
    if !context.genotyping.emit_reference_confidence {
        return Ok(calls);
    }

    // `containsCalls`: a genotype with any allele called.
    let contains_calls = calls.iter().any(|call| {
        call.genotypes
            .iter()
            .any(|g| g.alleles.iter().any(|a| !a.is_no_call()))
    });
    if !contains_calls {
        return reference_model_for_no_variation(region, false, context);
    }
    // The left flank, the genotyping region base by base with its calls, then the right flank.
    let mut result = Vec::new();
    let (left, right) = match ramp_flanks {
        Some(flanks) => flanks,
        None => (
            context
                .trimmer
                .non_variant_left_flank_region(trimmed, &region, context.header)
                .map_err(|e| calling(format!("{e:?}")))?,
            context
                .trimmer
                .non_variant_right_flank_region(trimmed, &region, context.header)
                .map_err(|e| calling(format!("{e:?}")))?,
        ),
    };
    if let Some(flank) = left {
        result.extend(reference_model_for_no_variation(flank, false, context)?);
    }
    let sample = &context.samples[0];
    // `Utils.nonNull(refHaplotype, message)`, an `IllegalArgumentException`: a set of partially
    // determined haplotypes has none.
    if assembly.reference_haplotype().is_none() {
        return Err(calling(
            "java.lang.IllegalArgumentException: refHaplotype cannot be null".to_string(),
        ));
    }
    result.extend(
        calculate_ref_confidence(
            &reference_haplotype.bases(),
            genotyping_region.span(),
            genotyping_region.padded_span().start,
            likelihoods.sample_evidence(0).unwrap_or(&[]),
            sample,
            context.genotyping.sample_ploidy,
            &calls,
            context.header,
            &context.genotyping.ref_confidence,
            None,
        )
        .map_err(engine_failure)?,
    );
    if let Some(flank) = right {
        result.extend(reference_model_for_no_variation(flank, false, context)?);
    }
    Ok(result)
}

/// `referenceModelForNoVariation(region, needsToBeFinalized)`: nothing outside reference-confidence
/// mode; otherwise the region finalized when it has not been, its non-passing reads removed, and
/// a reference-confidence record for each base of its span against its own reads, with no calls.
fn reference_model_for_no_variation(
    mut region: AssemblyRegion,
    needs_to_be_finalized: bool,
    context: &CallRegionContext<'_>,
) -> Result<Vec<VariantContext>, RegionAssemblyError> {
    if !context.genotyping.emit_reference_confidence {
        return Ok(Vec::new());
    }
    if needs_to_be_finalized {
        // `usePileupDetection`, off by default, is what tracks the hard-clipped reads here.
        let mut arguments = finalize_arguments(context.assembly);
        arguments.track_hardclipped_reads = false;
        finalize_region(&mut region, &arguments, context.header, context.samples)
            .map_err(RegionAssemblyError::Finalize)?;
    }
    let threshold = context.mapping_quality_threshold;
    region.retain_reads(|read| !fails_genotyping_filters(read, threshold));
    let padded = region.padded_span().clone();
    let ref_bases = &context.contig_bases[padded.start as usize - 1..padded.end as usize];
    let sample = &context.samples[0];
    let evidence: Vec<BamRecord> = region
        .reads()
        .iter()
        .filter(|read| sample_name(read, context.header).as_ref() == Some(sample))
        .cloned()
        .collect();
    calculate_ref_confidence(
        ref_bases,
        region.span(),
        padded.start,
        &evidence,
        sample,
        context.genotyping.sample_ploidy,
        &[],
        context.header,
        &context.genotyping.ref_confidence,
        None,
    )
    .map_err(engine_failure)
}

/// A reference exception as `class: message`, with the class fully qualified, which is how the
/// caller recognises one; a limitation as it was.
fn engine_failure(error: crate::genotyping_engine::EngineError) -> RegionAssemblyError {
    match error {
        crate::genotyping_engine::EngineError::Runtime { class, message } => {
            RegionAssemblyError::Calling(format!("{}: {message}", qualified_class(&class)))
        }
        other => RegionAssemblyError::Calling(format!("{other:?}")),
    }
}

/// A class the engine names short, fully qualified.
pub(crate) fn qualified_class(class: &str) -> String {
    match class {
        "GATKException" => "org.broadinstitute.hellbender.exceptions.GATKException".to_string(),
        short if !short.contains('.') => format!("java.lang.{short}"),
        full => full.to_string(),
    }
}

/// `filterNonPassingReads`' test: too short, mapped below the genotyping threshold, or with a
/// mate on another contig.
fn fails_genotyping_filters(read: &BamRecord, mapping_quality_threshold: u8) -> bool {
    unclipped_read_length(read) < READ_LENGTH_FILTER_THRESHOLD
        || read.mapping_quality < mapping_quality_threshold
        || !gatk_readfilter::mate_on_same_contig_or_no_mapped_mate(read)
}

/// The `finalizeRegion` switches HaplotypeCaller passes.
fn finalize_arguments(arguments: &RegionAssemblyArguments) -> FinalizeArguments {
    FinalizeArguments {
        error_correct_reads: false,
        dont_use_soft_clipped_bases: arguments.dont_use_soft_clipped_bases,
        min_tail_quality: arguments.min_base_quality_score.wrapping_sub(1),
        correct_overlapping_base_qualities: !arguments.do_not_correct_overlapping_base_qualities,
        soft_clip_low_quality_ends: arguments.soft_clip_low_quality_ends,
        override_softclip_fragment_check: arguments.override_softclip_fragment_check,
        track_hardclipped_reads: true,
    }
}

/// HaplotypeCaller's `minBaseQualityScore` and the switches `assembleReads` passes on.
#[derive(Debug, Clone)]
pub struct RegionAssemblyArguments {
    pub min_base_quality_score: u8,
    pub dont_use_soft_clipped_bases: bool,
    pub soft_clip_low_quality_ends: bool,
    pub override_softclip_fragment_check: bool,
    pub do_not_correct_overlapping_base_qualities: bool,
    /// `getDanglingEndSWParameters()`, `getHaplotypeToReferenceSWParameters()` and
    /// `getReadToHaplotypeSWParameters()`, the `--smith-waterman-*` values.
    pub dangling_end_sw: SwParameters,
    pub haplotype_to_reference_sw: SwParameters,
    pub read_to_haplotype_sw: SwParameters,
}

impl Default for RegionAssemblyArguments {
    fn default() -> Self {
        RegionAssemblyArguments {
            min_base_quality_score: 10,
            dont_use_soft_clipped_bases: false,
            soft_clip_low_quality_ends: false,
            override_softclip_fragment_check: false,
            do_not_correct_overlapping_base_qualities: false,
            dangling_end_sw: STANDARD_NGS,
            haplotype_to_reference_sw: NEW_SW_PARAMETERS,
            read_to_haplotype_sw: ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
        }
    }
}

/// `HaplotypeCallerReadThreadingAssemblerArgumentCollection`, as the command line sets it.
#[derive(Debug, Clone)]
pub struct AssemblerArguments {
    /// What the assembler is constructed with; `prune_factor` is `--min-pruning`, which adaptive
    /// pruning replaces with zero.
    pub settings: AssemblerSettings,
    /// `!--do-not-recover-dangling-branches`.
    pub recover_dangling_branches: bool,
    pub recover_all_dangling_branches: bool,
    pub min_dangling_branch_length: i32,
    /// `--disable-artificial-haplotype-recovery`.
    pub disable_artificial_haplotype_recovery: bool,
}

impl Default for AssemblerArguments {
    fn default() -> Self {
        AssemblerArguments {
            settings: AssemblerSettings {
                max_allowed_paths: 128,
                kmer_sizes: vec![10, 25],
                dont_increase_kmer_sizes_for_cycles: false,
                allow_non_unique_kmers_in_ref: false,
                num_pruning_samples: 1,
                prune_factor: 2,
                use_adaptive_pruning: false,
                initial_error_rate_for_pruning: 0.001,
                pruning_log_odds_threshold: gatk_engine::allele_likelihoods::log10_to_log(1.0),
                pruning_seeding_log_odds_threshold: gatk_engine::allele_likelihoods::log10_to_log(
                    4.0,
                ),
                max_unpruned_variants: 100,
                enable_legacy_graph_cycle_detection: false,
                min_matching_bases_to_dangling_end_recovery: -1,
                use_linked_debruijn_graphs: false,
            },
            recover_dangling_branches: true,
            recover_all_dangling_branches: false,
            min_dangling_branch_length: 4,
            disable_artificial_haplotype_recovery: false,
        }
    }
}

/// `hcArgs.createReadThreadingAssembler()`: the assembler the arguments describe, with the
/// minimum base quality `AssemblyBasedCallerArgumentCollection` sets on it.
pub fn assembler_from(
    arguments: &AssemblerArguments,
    min_base_quality_score: u8,
) -> Result<ReadThreadingAssembler, AssemblerError> {
    let mut settings = arguments.settings.clone();
    if settings.use_adaptive_pruning {
        settings.prune_factor = 0;
    }
    let mut assembler = ReadThreadingAssembler::new(&settings)?;
    assembler.set_recover_dangling_branches(arguments.recover_dangling_branches);
    assembler.set_recover_all_dangling_branches(arguments.recover_all_dangling_branches);
    assembler.set_min_dangling_branch_length(arguments.min_dangling_branch_length);
    assembler
        .set_artificial_haplotype_recovery_mode(arguments.disable_artificial_haplotype_recovery)?;
    assembler.set_min_base_quality_to_use_in_assembly(min_base_quality_score);
    Ok(assembler)
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
        use_linked_debruijn_graphs: false,
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
    finalize_region(region, &finalize_arguments(arguments), header, samples)
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
            &arguments.dangling_end_sw,
            &arguments.haplotype_to_reference_sw,
        )
        .map_err(RegionAssemblyError::Assembler)
}

/// `assembleReads` with `bypassAssembly`, as the post-filter on ramp runs it: the region finalized
/// as for an assembly, then `generateEmptyLLocalAssemblyResult`, a set holding the reference
/// haplotype alone with the whole region to genotype.
fn assemble_bypassed(
    region: &mut AssemblyRegion,
    context: &CallRegionContext<'_>,
) -> Result<AssemblyResultSet, RegionAssemblyError> {
    finalize_region(
        region,
        &finalize_arguments(context.assembly),
        context.header,
        context.samples,
    )
    .map_err(RegionAssemblyError::Finalize)?;
    let contig_length = context.contig_bases.len() as i32;
    let ref_loc = padded_reference_loc(region, REFERENCE_PADDING_FOR_ASSEMBLY, contig_length);
    let full_reference =
        context.contig_bases[ref_loc.start as usize - 1..ref_loc.end as usize].to_vec();
    let ref_haplotype = reference_haplotype(region, context.contig_bases, &ref_loc)?;
    let mut set = AssemblyResultSet::new();
    set.set_region_for_genotyping(region.clone());
    set.set_full_reference_with_padding(full_reference);
    set.set_padded_reference_loc(ref_loc);
    set.add(ref_haplotype)
        .map_err(|e| RegionAssemblyError::Calling(format!("{e:?}")))?;
    Ok(set)
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
