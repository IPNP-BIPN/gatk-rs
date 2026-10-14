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
use gatk_engine::pair_hmm_likelihood_engine::LikelihoodEngineArguments;
use gatk_engine::smith_waterman::SmithWatermanJavaAligner;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::genotyping_engine::{Configuration, GenotypePriors, GenotypingEngine};
use crate::haplotype_caller_engine::{ActiveRegionEvaluator, ActivityType, IsActiveArguments};
use crate::hc_genotyping::{CallAnnotator, HcGenotypingArguments};
use crate::hc_region::{
    assembler_from, call_region, AssemblerArguments, CallRegionContext, RegionAssemblyArguments,
    RegionAssemblyError,
};
use crate::variant_annotator_engine::Engine as AnnotatorEngine;

/// Everything the command line decides for one run.
#[derive(Debug, Clone)]
pub struct HaplotypeCallerArguments {
    pub region: AssemblyRegionArgs,
    pub is_active: IsActiveArguments,
    /// `PileupDetectionArgumentCollection`.
    pub pileup_detection: crate::pileup_based_alleles::PileupDetectionArguments,
    /// `--use-pdhmm` and its switches.
    pub pdhmm: Option<PdhmmArguments>,
    /// `--alleles`' records, and `--force-call-filtered-alleles`.
    pub given_alleles: Option<(Vec<htsjdk_vcf::variant::VariantContext>, bool)>,
    /// `--bam-output` and how it is written.
    pub bamout: Option<BamoutArguments>,
    /// `--assembly-region-out`.
    pub assembly_region_out: Option<std::path::PathBuf>,
    /// `--keep-rg`.
    pub keep_rg: Option<String>,
    /// `--just-determine-active-regions`.
    pub just_determine_active_regions: bool,
    /// `--dont-genotype`.
    pub dont_genotype: bool,
    /// `--haplotype-debug-histogram-output`.
    pub haplotype_histogram_output: Option<std::path::PathBuf>,
    /// `--alm-path` and `--alm-interval`.
    pub alm: Option<(std::path::PathBuf, Option<String>)>,
    /// `--debug-assembly-region-state`.
    pub assembly_state_output: Option<std::path::PathBuf>,
    /// `--debug-genotyper-output`.
    pub genotyper_debug_output: Option<std::path::PathBuf>,
    /// The padding the trimmer reads, which is `--assembly-region-padding` and the
    /// genotyping paddings.
    pub trimmer: TrimmerArguments,
    pub assembly: RegionAssemblyArguments,
    /// The read-threading assembler's own arguments.
    pub assembler: AssemblerArguments,
    pub likelihoods: LikelihoodEngineArguments,
    /// `--likelihood-calculation-engine FlowBased`: the flow engine's own arguments.
    pub flow_engine: Option<crate::flow_likelihood_engine::FlowEngineArguments>,
    /// `--likelihood-calculation-engine FlowBasedHMM`: its arguments.
    pub flow_hmm: Option<crate::flow_likelihood_engine::FlowHmmArguments>,
    /// `--use-flow-aligner-for-stepwise-hc-filtering`: the FlowBased engine that scores the reads
    /// for allele filtering, before the main engine scores them again over what is left.
    pub stepwise_flow: Option<crate::flow_likelihood_engine::FlowEngineArguments>,
    pub genotyping: HcGenotypingArguments,
    pub calling: Configuration,
    pub priors: Priors,
    /// `--sample-name`, already cleared by the caller when it names the only sample: the samples
    /// genotyped, and the reads each region keeps.
    pub sample_name: Option<String>,
    /// `--disable-optimizations`.
    pub disable_optimizations: bool,
    /// `--flow-assembly-collapse-hmer-size` (-1 asks the read groups) and
    /// `--flow-assembly-collapse-partial-mode`.
    pub flow_collapse_hmer_size: i32,
    pub flow_collapse_partial_mode: bool,
    /// `--ramps-debug-post-assembler-on`.
    pub ramps_debug_post_assembler_on: bool,
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
    /// The calls of the regions before the failing one. The reference has already handed them to
    /// its writer, and closes it in a `finally`, so an error the writer raises over them is the
    /// one the run reports.
    pub calls_before: Vec<VariantContext>,
}

impl HaplotypeCallerError {
    fn new(class: &str, message: impl Into<String>) -> Self {
        HaplotypeCallerError {
            class: class.to_string(),
            message: message.into(),
            calls_before: Vec::new(),
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

/// `--bam-output`, `--bam-writer-type`, `--create-output-bam-index`, `--create-output-bam-md5` and
/// the BGZF settings the writer inherits.
#[derive(Debug, Clone, PartialEq)]
pub struct BamoutArguments {
    pub kind: crate::haplotype_bam_writer::WriterType,
    pub path: std::path::PathBuf,
    pub create_index: bool,
    pub create_md5: bool,
    pub level: u32,
    pub deflater: htsjdk_bgzf::Deflater,
}

/// `HaplotypeBAMWriter.close()`: the BAM, its index (the output's extension replaced by `.bai`)
/// and its MD5.
pub fn close_bamout(
    writer: crate::haplotype_bam_writer::HaplotypeBamWriter,
    arguments: &BamoutArguments,
) -> Result<(), String> {
    let (bam, bai) = writer.finish(arguments.create_index, arguments.level, arguments.deflater)?;
    std::fs::write(&arguments.path, &bam).map_err(|e| e.to_string())?;
    if let Some(bai) = bai {
        std::fs::write(arguments.path.with_extension("bai"), bai).map_err(|e| e.to_string())?;
    }
    if arguments.create_md5 {
        let mut digest = arguments.path.clone().into_os_string();
        digest.push(".md5");
        std::fs::write(digest, crate::gather_bam_files::md5_file(&bam))
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// `--use-pdhmm`'s arguments.
#[derive(Debug, Clone, PartialEq)]
pub struct PdhmmArguments {
    /// `--make-determined-haps-from-pd-code`.
    pub use_determined_haplotypes: bool,
    /// `--fallback-gga-if-pdhmm-fails`.
    pub gga_fallback: bool,
    /// `--use-pdhmm-overlap-optimization`.
    pub overlap_optimization: bool,
    /// `--pdhmm-results-file`.
    pub results_file: Option<std::path::PathBuf>,
}

/// The fully qualified class a `class: message` names, `java.lang.` or GATK's own.
fn thrown_class(message: &str) -> Option<&str> {
    let (class, _) = message.split_once(": ")?;
    let qualified = class.starts_with("java.lang.") || class.starts_with("org.broadinstitute.");
    (qualified && !class.contains(' ')).then_some(class)
}

/// The run: `reads` already filtered and in coordinate order, `intervals` the traversal
/// intervals, and `contig_bases` a contig's whole sequence by name.
#[allow(clippy::too_many_arguments)]
pub fn call_variants(
    reads: &[BamRecord],
    header: &SamHeader,
    intervals: &[SimpleInterval],
    arguments: &HaplotypeCallerArguments,
    annotator: &AnnotatorEngine,
    contig_bases: &mut dyn FnMut(&str) -> Result<Vec<u8>, String>,
    random: &mut JavaRandom,
    assembly_variants: Option<
        &std::cell::RefCell<crate::assembly_variants_out::AssemblyVariantsOut>,
    >,
) -> Result<Vec<VariantContext>, HaplotypeCallerError> {
    // `HaplotypeCallerGenotypingDebugger.initialize`, as the engine starts: its stream is never
    // closed, but writes through, so the file holds what was printed before any failure.
    let _genotyper_debug = arguments
        .genotyper_debug_output
        .clone()
        .map(GenotyperDebugFile::open);
    let mut calls = Vec::new();
    // `PDPairHMMLikelihoodCalculationEngine`, built with the engine, and its results file, opened
    // then with its header line and closed as the engine shuts down, after a failure too.
    let pd_engine = arguments
        .pdhmm
        .as_ref()
        .map(|_| {
            // Its own results file is `--pdhmm-results-file`, not the PairHMM's.
            gatk_engine::pair_hmm_likelihood_engine::PairHmmLikelihoodEngine::new(
                gatk_engine::pair_hmm_likelihood_engine::LikelihoodEngineArguments {
                    results_file: None,
                    ..arguments.likelihoods.clone()
                },
            )
        })
        .transpose()
        .map_err(|e| HaplotypeCallerError::new("IllegalArgumentException", format!("{e:?}")))?;
    let results_file = arguments
        .pdhmm
        .as_ref()
        .and_then(|pdhmm| pdhmm.results_file.clone());
    let results = results_file.as_ref().map(|_| {
        std::cell::RefCell::new(
            "# hap-bases\thap-pd-bases\tread-bases\tread-qual\tread-ins-qual\tread-del-qual\tgcp\texpected-result"
                .to_string(),
        )
    });
    let bamout = arguments.bamout.as_ref().map(|b| {
        std::cell::RefCell::new(crate::haplotype_bam_writer::HaplotypeBamWriter::new(
            b.kind, header,
        ))
    });
    // `createAlleleLikelihoodWriter`: the interval is parsed as the engine starts, and the file
    // created then, written as regions go.
    let alm = match &arguments.alm {
        Some((_, interval)) => Some(std::cell::RefCell::new(
            crate::allele_likelihood_writer::AlleleLikelihoodWriter::new(
                crate::allele_likelihood_writer::parse_simple_interval(interval.as_deref())
                    .map_err(|message| {
                        let (class, rest) = message.split_once(": ").expect("class: message");
                        HaplotypeCallerError::new(class, rest.to_string())
                    })?,
            ),
        )),
        None => None,
    };
    let assembly_state = arguments
        .assembly_state_output
        .as_ref()
        .map(|_| std::cell::RefCell::new(String::new()));
    let region_out = arguments
        .assembly_region_out
        .as_ref()
        .map(|_| std::cell::RefCell::new(crate::hc_region::ASSEMBLY_REGION_OUT_HEADER.to_string()));
    let outcome = call_variants_into(
        reads,
        header,
        intervals,
        arguments,
        annotator,
        contig_bases,
        random,
        pd_engine.as_ref(),
        results.as_ref(),
        bamout.as_ref(),
        region_out.as_ref(),
        assembly_variants,
        alm.as_ref(),
        assembly_state.as_ref(),
        &mut calls,
    );
    if let (Some(path), Some(state)) = (&arguments.assembly_state_output, assembly_state) {
        std::fs::write(path, state.into_inner()).map_err(|e| {
            HaplotypeCallerError::new("UserException.CouldNotCreateOutputFile", e.to_string())
        })?;
    }
    if let (Some((path, _)), Some(alm)) = (&arguments.alm, alm) {
        std::fs::write(path, alm.into_inner().text).map_err(|e| {
            HaplotypeCallerError::new("UserException.CouldNotCreateOutputFile", e.to_string())
        })?;
    }
    if let (Some(path), Some(text)) = (&arguments.assembly_region_out, region_out) {
        std::fs::write(path, text.into_inner()).map_err(|e| {
            HaplotypeCallerError::new("UserException.CouldNotCreateOutputFile", e.to_string())
        })?;
    }
    // The writer closes with the engine, after a failure too.
    if let (Some(writer), Some(settings)) = (bamout, &arguments.bamout) {
        close_bamout(writer.into_inner(), settings)
            .map_err(|e| HaplotypeCallerError::new("UserException.CouldNotCreateOutputFile", e))?;
    }
    if let (Some(path), Some(text)) = (&results_file, &results) {
        std::fs::write(path, text.borrow().as_bytes()).map_err(|e| {
            HaplotypeCallerError::new("UserException.CouldNotCreateOutputFile", e.to_string())
        })?;
    }
    match outcome {
        Ok(()) => Ok(calls),
        Err(mut error) => {
            error.calls_before = calls;
            Err(error)
        }
    }
}

/// [`call_variants`], appending each region's calls to `calls` as it goes.
#[allow(clippy::too_many_arguments)]
fn call_variants_into(
    reads: &[BamRecord],
    header: &SamHeader,
    intervals: &[SimpleInterval],
    arguments: &HaplotypeCallerArguments,
    annotator: &AnnotatorEngine,
    contig_bases: &mut dyn FnMut(&str) -> Result<Vec<u8>, String>,
    random: &mut JavaRandom,
    pd_engine: Option<&gatk_engine::pair_hmm_likelihood_engine::PairHmmLikelihoodEngine>,
    pd_results: Option<&std::cell::RefCell<String>>,
    bamout: Option<&std::cell::RefCell<crate::haplotype_bam_writer::HaplotypeBamWriter>>,
    region_out: Option<&std::cell::RefCell<String>>,
    assembly_variants: Option<
        &std::cell::RefCell<crate::assembly_variants_out::AssemblyVariantsOut>,
    >,
    alm: Option<&std::cell::RefCell<crate::allele_likelihood_writer::AlleleLikelihoodWriter>>,
    assembly_state: Option<&std::cell::RefCell<String>>,
    calls: &mut Vec<VariantContext>,
) -> Result<(), HaplotypeCallerError> {
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
    let mut assembler = assembler_from(
        &arguments.assembler,
        arguments.assembly.min_base_quality_score,
    )
    .map_err(|e| HaplotypeCallerError::new("IllegalArgumentException", format!("{e:?}")))?;
    // `--haplotype-debug-histogram-output`: kept as the assembler runs, printed as the engine shuts
    // down, after a failure too.
    if arguments.haplotype_histogram_output.is_some() {
        assembler.keep_debug_histograms();
    }
    let assembler = assembler;
    struct HistogramsOnShutdown<'a> {
        assembler: &'a gatk_engine::read_threading_assembler::ReadThreadingAssembler,
        path: Option<&'a std::path::Path>,
    }
    impl Drop for HistogramsOnShutdown<'_> {
        fn drop(&mut self) {
            if let (Some(path), Some(histograms)) = (self.path, self.assembler.debug_histograms()) {
                let _ = std::fs::write(path, histograms.report());
            }
        }
    }
    let _histograms = HistogramsOnShutdown {
        assembler: &assembler,
        path: arguments.haplotype_histogram_output.as_deref(),
    };
    // `filterStepLikelihoodCalculationEngine`.
    let stepwise_engine = arguments.stepwise_flow.as_ref().map(|flow| {
        crate::flow_likelihood_engine::ReadLikelihoodEngine::FlowBased(
            crate::flow_likelihood_engine::FlowLikelihoodEngine {
                arguments: flow.clone(),
                log10_global_read_mismapping_rate: arguments
                    .likelihoods
                    .log10_global_read_mismapping_rate,
                expected_error_rate_per_base: arguments.likelihoods.expected_error_rate_per_base,
                dynamic_disqualification: arguments.likelihoods.dynamic_disqualification,
                read_disqualification_scale: arguments.likelihoods.read_disqualification_scale,
            },
        )
    });
    // `createLikelihoodCalculationEngine`.
    let likelihood_engine = crate::flow_likelihood_engine::ReadLikelihoodEngine::create(
        &arguments.likelihoods,
        &arguments.flow_engine,
        &arguments.flow_hmm,
    )
    .map_err(|(class, message)| HaplotypeCallerError::new(&class, message))?;
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
    for group in group_intervals_by_contig(intervals) {
        let contig = group[0].contig.clone();
        let bases = contig_bases(&contig)
            .map_err(|message| HaplotypeCallerError::new("UserException", message))?;
        let mut failure: Option<HaplotypeCallerError> = None;
        // `region.getAlignmentData()`: every locus walked, kept only with pileup detection.
        let detection = &arguments.pileup_detection;
        let mut loci: Vec<crate::pileup_based_alleles::LocusSummary> = Vec::new();
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
                let ref_base = bases[locus.start as usize - 1];
                // A locus a usable given record overlaps is active, before anything else.
                if let Some((records, force_call_filtered)) = &arguments.given_alleles {
                    if crate::given_alleles::forced_at(
                        records,
                        &locus.contig,
                        locus.start,
                        *force_call_filtered,
                    ) {
                        if detection.use_pileup_detection {
                            let summary = match context {
                                Some(context) => crate::pileup_based_alleles::LocusSummary::new(
                                    context,
                                    0.0,
                                    ref_base,
                                    &bases,
                                    i32::from(arguments.assembly.min_base_quality_score),
                                    detection,
                                ),
                                None => crate::pileup_based_alleles::LocusSummary::empty(
                                    &locus.contig,
                                    locus.start,
                                    0.0,
                                    ref_base,
                                ),
                            };
                            loci.push(summary);
                        }
                        return (1.0, None);
                    }
                }
                let Some(context) = context else {
                    if detection.use_pileup_detection {
                        loci.push(crate::pileup_based_alleles::LocusSummary::empty(
                            &locus.contig,
                            locus.start,
                            0.0,
                            ref_base,
                        ));
                    }
                    return (0.0, None);
                };
                let state = evaluator.is_active(context, ref_base, header);
                if detection.use_pileup_detection {
                    let score = state.as_ref().map_or(0.0, |s| s.original_active_prob);
                    loci.push(crate::pileup_based_alleles::LocusSummary::new(
                        context,
                        score,
                        ref_base,
                        &bases,
                        i32::from(arguments.assembly.min_base_quality_score),
                        detection,
                    ));
                }
                match state {
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

        // `assembleReads`: a threshold of -1 is the largest maximal class of the read groups.
        let collapse_threshold = if arguments.flow_collapse_hmer_size
            == crate::haplotype_collapsing::DETERMINE_COLLAPSE_THRESHOLD
        {
            crate::haplotype_collapsing::determine_collapse_threshold(header)
                .map_err(|(class, message)| HaplotypeCallerError::new(class, message))?
        } else {
            arguments.flow_collapse_hmer_size
        };
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
            stepwise_engine: stepwise_engine.as_ref(),
            genotyping: &genotyping,
            mapping_quality_threshold: arguments.mapping_quality_threshold,
            off_ramp: off_ramp.as_ref(),
            on_ramp: on_ramp.as_ref(),
            disable_optimizations: arguments.disable_optimizations,
            collapse: crate::hc_region::CollapseSettings {
                raw: arguments.flow_collapse_hmer_size,
                determined: collapse_threshold,
                partial: arguments.flow_collapse_partial_mode,
            },
            ramps_debug_post_assembler_on: arguments.ramps_debug_post_assembler_on,
            // RampedHaplotypeCallerEngine reads no pileup events.
            pileup: (detection.use_pileup_detection && !genotyping.ramped_engine)
                .then_some((detection, loci.as_slice())),
            // RampedHaplotypeCallerEngine makes no PD haplotypes.
            bamout,
            keep_rg: arguments.keep_rg.as_deref(),
            just_determine_active_regions: arguments.just_determine_active_regions,
            dont_genotype: arguments.dont_genotype,
            assembly_variants,
            alm,
            assembly_state,
            given_alleles: arguments
                .given_alleles
                .as_ref()
                .map(|(records, force)| (records.as_slice(), *force)),
            pd: match (&arguments.pdhmm, pd_engine) {
                (Some(pdhmm), Some(engine)) if !genotyping.ramped_engine => {
                    Some(crate::hc_region::PdContext {
                        arguments: crate::pd_haplotypes::PdArguments {
                            use_determined_haplotypes: pdhmm.use_determined_haplotypes,
                            snp_adjacent_to_assembly_indel: detection
                                .snp_adjacent_to_assembly_indel,
                            debug: detection.debug_stdout,
                        },
                        gga_fallback: pdhmm.gga_fallback,
                        overlap_range: if pdhmm.overlap_optimization {
                            genotyping.informative_read_overlap_margin
                        } else {
                            -1
                        },
                        engine,
                        results: pd_results,
                    })
                }
                _ => None,
            },
        };
        for mut traversed in regions {
            if let Some(out) = region_out {
                out.borrow_mut()
                    .push_str(&crate::hc_region::assembly_region_rows(&traversed.region));
            }
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
            .map_err(|e| match e {
                RegionAssemblyError::NoVariationPresent => HaplotypeCallerError::new(
                    "java.lang.IllegalStateException",
                    "There is no variation present.".to_string(),
                ),
                // A message the reference's own exception names, as `class: message`.
                RegionAssemblyError::Calling(message) if thrown_class(&message).is_some() => {
                    let class = thrown_class(&message).expect("just matched");
                    HaplotypeCallerError::new(class, message[class.len() + 2..].to_string())
                }
                other => HaplotypeCallerError::new("GATKException", format!("{other:?}")),
            })?;
            calls.extend(region_calls);
        }
    }
    // `tearRamps`: the zip is written as the engine shuts down.
    if let (Some(ramp), Some((_, path))) = (&off_ramp, &arguments.off_ramp) {
        ramp.borrow()
            .write(path)
            .map_err(|message| HaplotypeCallerError::new("GATKException", message))?;
    }
    Ok(())
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

/// The `--debug-genotyper-output` file, written with what was printed when the run ends.
struct GenotyperDebugFile(std::path::PathBuf);

impl GenotyperDebugFile {
    fn open(path: std::path::PathBuf) -> GenotyperDebugFile {
        gatk_engine::genotyper_debug::initialize();
        GenotyperDebugFile(path)
    }
}

impl Drop for GenotyperDebugFile {
    fn drop(&mut self) {
        let text = gatk_engine::genotyper_debug::take().unwrap_or_default();
        // The reference opens the file as the engine starts and reports nothing after.
        let _ = std::fs::write(&self.0, text);
    }
}
