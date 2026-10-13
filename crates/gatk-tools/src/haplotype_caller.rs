//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller`
//! (`isActive`, `apply`) over `org.broadinstitute.hellbender.engine.AssemblyRegionWalker`'s
//! traversal (GATK 4.6.2.0): HaplotypeCaller in VCF mode, from filtered reads to calls.
//!
//! One contig at a time, as the walker's shards are: the reads overlapping the contig's padded
//! intervals go through the positional downsampler, the pileup at every locus through
//! [`ActiveRegionEvaluator`], the profile cuts the regions, and each region is called by
//! [`call_region`] as soon as its contig's regions are known. The calls come out in region order,
//! which is coordinate order since the regions do not overlap.
//!
//! What a reader would not guess: the samples are the read groups' `SM` values **sorted**,
//! because `ReadUtils.getSamplesFromHeader` collects them in a `TreeSet`.
//!
//! In reference-confidence mode each region's output is a record per base, which the command line
//! combines into blocks when it writes a GVCF.
//!
//! With `--sample-name` naming one of several samples, only that sample is genotyped and each
//! region keeps only its reads, but the activity profile reads every sample's reads as that one's.
//!
//! Not ported here: `--alleles`, pileup detection and the flow-based paths, all refused by the
//! command line.

use gatk_engine::allele_frequency_calculator::{AlleleFrequencyCalculator, Priors};
use gatk_engine::assembly_region_iterator::{group_intervals_by_contig, AssemblyRegionArgs};
use gatk_engine::assembly_region_trimmer::{AssemblyRegionTrimmer, TrimmerArguments};
use gatk_engine::assembly_region_walker::{needs_downsampler, traverse_with_pileups};
use gatk_engine::downsampling::positional_downsample_with;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_random::JavaRandom;
use gatk_engine::pair_hmm_likelihood_engine::{LikelihoodEngineArguments, PairHmmLikelihoodEngine};
use gatk_engine::smith_waterman::SmithWatermanJavaAligner;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::genotyping_engine::{Configuration, GenotypePriors, GenotypingEngine};
use crate::haplotype_caller_engine::{ActiveRegionEvaluator, ActivityType, IsActiveArguments};
use crate::hc_genotyping::{CallAnnotator, HcGenotypingArguments};
use crate::hc_region::{
    assembler_from, call_region, AssemblerArguments, CallRegionContext, RegionAssemblyArguments,
};
use crate::variant_annotator_engine::Engine as AnnotatorEngine;

/// Everything the command line decides for one run.
#[derive(Debug, Clone)]
pub struct HaplotypeCallerArguments {
    pub region: AssemblyRegionArgs,
    pub is_active: IsActiveArguments,
    /// The padding the trimmer reads, which is `--assembly-region-padding` and the
    /// genotyping paddings.
    pub trimmer: TrimmerArguments,
    pub assembly: RegionAssemblyArguments,
    /// The read-threading assembler's own arguments.
    pub assembler: AssemblerArguments,
    pub likelihoods: LikelihoodEngineArguments,
    pub genotyping: HcGenotypingArguments,
    pub calling: Configuration,
    pub priors: Priors,
    /// `--sample-name`, already cleared by the caller when it names the only sample: the samples
    /// genotyped, and the reads each region keeps.
    pub sample_name: Option<String>,
    /// `--mapping-quality-threshold-for-genotyping`.
    pub mapping_quality_threshold: u8,
    /// `--use-posteriors-to-calculate-qual`.
    pub use_posteriors_to_calculate_qual: bool,
    /// `nonRandomDownsamplingMode()`, which `RampedHaplotypeCaller` turns on.
    pub non_random_downsampling: bool,
    /// RampedHaplotypeCaller's `--off-ramp-type` and `--off-ramp-file`.
    pub off_ramp: Option<(crate::ramp_files::OffRampType, std::path::PathBuf)>,
    /// RampedHaplotypeCaller's `--on-ramp-type` and `--on-ramp-file`.
    pub on_ramp: Option<(crate::ramp_files::OnRampType, std::path::PathBuf)>,
    /// The BGZF level and deflater the run writes block-compressed files with.
    pub compression: (u32, htsjdk_bgzf::Deflater),
}

/// What a run refuses, with the Java class it would be thrown as.
#[derive(Debug, Clone, PartialEq)]
pub struct HaplotypeCallerError {
    pub class: String,
    pub message: String,
}

impl HaplotypeCallerError {
    fn new(class: &str, message: impl Into<String>) -> Self {
        HaplotypeCallerError {
            class: class.to_string(),
            message: message.into(),
        }
    }
}

/// `ReadUtils.getSamplesFromHeader`: the read groups' samples, unique and sorted.
pub fn samples_from_header(header: &SamHeader) -> Vec<String> {
    let mut samples: Vec<String> = header
        .read_groups
        .iter()
        .filter_map(|group| group.attributes.get("SM").map(str::to_string))
        .collect();
    samples.sort();
    samples.dedup();
    samples
}

/// The run: `reads` already filtered and in coordinate order, `intervals` the traversal
/// intervals, and `contig_bases` a contig's whole sequence by name.
pub fn call_variants(
    reads: &[BamRecord],
    header: &SamHeader,
    intervals: &[SimpleInterval],
    arguments: &HaplotypeCallerArguments,
    annotator: &AnnotatorEngine,
    contig_bases: &mut dyn FnMut(&str) -> Result<Vec<u8>, String>,
    random: &mut JavaRandom,
) -> Result<Vec<VariantContext>, HaplotypeCallerError> {
    let header_samples = samples_from_header(header);
    let partition: Vec<Option<String>> = header_samples.iter().cloned().map(Some).collect();
    // `initializeSamples`: one sample named, and the activity profile then reads every sample's
    // reads as that one's (`splitContextBySampleName(assumedSingleSample, header)`).
    let samples = match &arguments.sample_name {
        Some(name) => vec![name.clone()],
        None => header_samples,
    };
    let mut evaluator = ActiveRegionEvaluator::new(arguments.is_active.clone(), samples.clone());
    evaluator.configure_from_caller(
        arguments.calling.assignment_method,
        arguments.use_posteriors_to_calculate_qual,
        arguments.calling.annotate_all_sites_with_pls,
    );
    let mut engine = GenotypingEngine::new(
        arguments.calling,
        AlleleFrequencyCalculator::make_calculator(&arguments.priors),
    );
    // `resolveGenotypePriorCalculator` without DRAGstr: Hardy-Weinberg over the two
    // heterozygosities, which only `USE_POSTERIOR_PROBABILITIES` reads.
    engine.set_genotype_priors(GenotypePriors::assuming_hw(
        jmath::math::log10(arguments.priors.snp_heterozygosity),
        jmath::math::log10(arguments.priors.indel_heterozygosity),
    ));
    engine.set_use_posteriors_to_calculate_qual(arguments.use_posteriors_to_calculate_qual);
    let assembler = assembler_from(
        &arguments.assembler,
        arguments.assembly.min_base_quality_score,
    )
    .map_err(|e| HaplotypeCallerError::new("IllegalArgumentException", format!("{e:?}")))?;
    let likelihood_engine = PairHmmLikelihoodEngine::new(arguments.likelihoods.clone())
        .map_err(|e| HaplotypeCallerError::new("IllegalArgumentException", format!("{e:?}")))?;
    let genotyping = arguments.genotyping.clone();

    // `buildRamps`: the off ramp is opened before the traversal and written when it ends.
    let off_ramp = arguments.off_ramp.as_ref().map(|(kind, _)| {
        std::cell::RefCell::new(crate::ramp_files::OffRamp::new(
            *kind,
            arguments.compression,
        ))
    });
    // The on ramp is read whole when the engine is built.
    let on_ramp = arguments
        .on_ramp
        .as_ref()
        .map(|(kind, path)| crate::ramp_files::OnRamp::open(*kind, path))
        .transpose()
        .map_err(|message| HaplotypeCallerError::new("RuntimeException", message))?;
    let mut calls = Vec::new();
    for group in group_intervals_by_contig(intervals) {
        let contig = group[0].contig.clone();
        let bases = contig_bases(&contig)
            .map_err(|message| HaplotypeCallerError::new("UserException", message))?;
        let mut failure: Option<HaplotypeCallerError> = None;
        let regions = traverse_with_pileups(
            reads,
            &group,
            &partition,
            &arguments.region,
            header,
            &mut |reads| {
                // `createDownsampler()`: a `PositionalDownsampler` only above zero.
                if needs_downsampler(&arguments.region) {
                    positional_downsample_with(
                        reads,
                        arguments.region.max_reads_per_alignment_start as usize,
                        random,
                        arguments.non_random_downsampling,
                    )
                } else {
                    reads.to_vec()
                }
            },
            &mut |locus, context| {
                let Some(context) = context else {
                    return (0.0, None);
                };
                let ref_base = bases[locus.start as usize - 1];
                match evaluator.is_active(context, ref_base, header) {
                    Ok(state) => (
                        state.prob,
                        (state.kind == ActivityType::HighQualitySoftClips)
                            .then_some(state.result_value)
                            .flatten(),
                    ),
                    Err(error) => {
                        failure.get_or_insert(HaplotypeCallerError::new(
                            "IllegalStateException",
                            format!("{error:?}"),
                        ));
                        (0.0, None)
                    }
                }
            },
        )
        .map_err(|error| HaplotypeCallerError::new(error.class(), error.message()))?;
        if let Some(failure) = failure {
            return Err(failure);
        }

        let trimmer = AssemblyRegionTrimmer::new(arguments.trimmer.clone(), bases.len() as i32)
            .map_err(|e| HaplotypeCallerError::new("IllegalArgumentException", format!("{e:?}")))?;
        let context = CallRegionContext {
            header,
            samples: &samples,
            contig_bases: &bases,
            assembly: &arguments.assembly,
            assembler: &assembler,
            aligner: &SmithWatermanJavaAligner,
            trimmer: &trimmer,
            likelihood_engine: &likelihood_engine,
            genotyping: &genotyping,
            mapping_quality_threshold: arguments.mapping_quality_threshold,
            off_ramp: off_ramp.as_ref(),
            on_ramp: on_ramp.as_ref(),
        };
        for mut traversed in regions {
            // `removeReadsFromAllSamplesExcept`.
            if let Some(name) = &arguments.sample_name {
                traversed.region.retain_reads(|read| {
                    gatk_engine::read_pileup::sample_name(read, header).as_deref() == Some(name)
                });
            }
            let region_calls = call_region(
                traversed.region,
                &context,
                &mut engine,
                Some(CallAnnotator {
                    engine: annotator,
                    random,
                }),
            )
            .map_err(|e| HaplotypeCallerError::new("GATKException", format!("{e:?}")))?;
            calls.extend(region_calls);
        }
    }
    // `tearRamps`: the zip is written as the engine shuts down.
    if let (Some(ramp), Some((_, path))) = (&off_ramp, &arguments.off_ramp) {
        ramp.borrow()
            .write(path)
            .map_err(|message| HaplotypeCallerError::new("GATKException", message))?;
    }
    Ok(calls)
}

/// `IUPACReadTransformer(true)`, HaplotypeCaller's pre-filter transformer: every IUPAC ambiguity
/// code (and a lower-case `n`) becomes `N`, the four bases in either case are kept, and anything
/// else is `BaseUtils.convertIUPACtoN`'s strict refusal, which prints the base as the signed byte
/// Java concatenates.
pub fn iupac_to_n(read: &mut BamRecord) -> Result<(), String> {
    for base in read.read_bases.iter_mut() {
        match *base {
            b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => {}
            b'N' | b'n' | b'R' | b'r' | b'Y' | b'y' | b'M' | b'm' | b'K' | b'k' | b'W' | b'w'
            | b'S' | b's' | b'B' | b'b' | b'D' | b'd' | b'H' | b'h' | b'V' | b'v' => *base = b'N',
            other => {
                return Err(format!(
                    "Bad input: We encountered a non-standard non-IUPAC base in the provided \
                     input sequence: '{}'",
                    other as i8
                ))
            }
        }
    }
    Ok(())
}
