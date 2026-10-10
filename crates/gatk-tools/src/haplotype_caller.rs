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
//! What a reader would not guess: physical phasing is **off** in VCF mode whatever
//! `--do-not-run-physical-phasing` says, since `validateAndInitializeArgs` turns it off unless
//! reference confidence is emitted; and the samples are the read groups' `SM` values **sorted**,
//! because `ReadUtils.getSamplesFromHeader` collects them in a `TreeSet`.
//!
//! Not ported here: reference-confidence (GVCF) mode, `--alleles`, pileup detection, DRAGEN mode
//! and the flow-based paths, all refused by the command line.

use gatk_engine::allele_frequency_calculator::{AlleleFrequencyCalculator, Priors};
use gatk_engine::assembly_region_iterator::{group_intervals_by_contig, AssemblyRegionArgs};
use gatk_engine::assembly_region_trimmer::{AssemblyRegionTrimmer, TrimmerArguments};
use gatk_engine::assembly_region_walker::traverse_with_pileups;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_random::JavaRandom;
use gatk_engine::pair_hmm_likelihood_engine::{LikelihoodEngineArguments, PairHmmLikelihoodEngine};
use gatk_engine::smith_waterman::SmithWatermanJavaAligner;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::genotyping_engine::{Configuration, GenotypingEngine};
use crate::haplotype_caller_engine::{ActiveRegionEvaluator, ActivityType, IsActiveArguments};
use crate::hc_genotyping::{CallAnnotator, HcGenotypingArguments};
use crate::hc_region::{
    call_region, haplotype_caller_assembler, CallRegionContext, RegionAssemblyArguments,
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
    pub likelihoods: LikelihoodEngineArguments,
    pub genotyping: HcGenotypingArguments,
    pub calling: Configuration,
    pub priors: Priors,
    /// `--mapping-quality-threshold-for-genotyping`.
    pub mapping_quality_threshold: u8,
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
    let samples = samples_from_header(header);
    let partition: Vec<Option<String>> = samples.iter().cloned().map(Some).collect();
    let mut evaluator = ActiveRegionEvaluator::new(arguments.is_active.clone(), samples.clone());
    let mut engine = GenotypingEngine::new(
        arguments.calling,
        AlleleFrequencyCalculator::make_calculator(&arguments.priors),
    );
    let assembler = haplotype_caller_assembler(arguments.assembly.min_base_quality_score);
    let likelihood_engine = PairHmmLikelihoodEngine::new(arguments.likelihoods.clone())
        .map_err(|e| HaplotypeCallerError::new("IllegalArgumentException", format!("{e:?}")))?;
    // Phasing needs reference confidence, which VCF mode never emits.
    let genotyping = HcGenotypingArguments {
        do_physical_phasing: false,
        ..arguments.genotyping.clone()
    };

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
            random,
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
        };
        for traversed in regions {
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
