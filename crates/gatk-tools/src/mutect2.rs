//! Ported from `org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine` (`isActive`,
//! `callRegion`) and `SomaticGenotypingEngine.callMutations` (GATK 4.6.2.0): Mutect2 from filtered
//! reads to somatic calls.
//!
//! The traversal is HaplotypeCaller's with `MutectDownsampler` in place of the positional one. A
//! locus is active when the tumor's likeliest alternate (a base, another substitution, or any
//! indel, the indels' qualities set by their length) has log odds of a somatic allele fraction
//! against error at least `--initial-tumor-lod`, unless a normal carries the same alternate in
//! more than three tenths of its reads with a quality sum above 100.
//!
//! A region is assembled (overlapping mates' base and indel qualities capped at half the PCR
//! qualities, unmarked duplicates of distant pairs dropped, adaptive pruning), trimmed, and its
//! reads' likelihoods taken in natural log and grouped into fragments. At each event start the
//! fragments are marginalized onto the merged alleles; an alternate is emitted when the tumor's
//! log odds (the Dirichlet evidence with and without it) clear `--tumor-lod-to-emit`, and
//! genotyped when a normal's diploid log odds clear `--normal-lod` as well. Each tumor genotype
//! carries every emitted allele and the posterior-mean allele fractions; a normal's is hom-ref.
//! The calls are trimmed, annotated with Mutect2's annotations over the read and fragment
//! matrices, phased, and given `ECNT` (the somatic events in the region) and `ECNTH` (those on the
//! best-supported haplotype carrying each alternate).
//!
//! In reference-confidence mode the activity and emission thresholds are minus infinity, every call
//! carries `<NON_REF>`, and each region is reported base by base by the somatic reference model
//! ([`crate::reference_confidence_model`]): its flanks and any region without a call entirely so.
//!
//! Not ported here: `--alleles`, `--germline-resource`,
//! `--panel-of-normals`, `--f1r2-tar-gz`, the bamout, a germline-sites fraction above zero (whose
//! generator is unseeded) and the Permutect datasets, all refused by the
//! command line.

use std::collections::HashMap;

use gatk_engine::alignment_utils::{normalize_alleles, IndexRange};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_mapping::{
    create_allele_mapper, make_merged_variant_context, variants_from_active_haplotypes,
};
use gatk_engine::assembly_based_caller_utils::{
    clean_overlapping_read_pairs_with_indels, finalize_region, FinalizeArguments,
};
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::assembly_region_iterator::{group_intervals_by_contig, AssemblyRegionArgs};
use gatk_engine::assembly_region_trimmer::{AssemblyRegionTrimmer, TrimmerArguments};
use gatk_engine::assembly_region_walker::traverse_with_pileups;
use gatk_engine::downsampling::mutect_downsample;
use gatk_engine::event_map::{build_event_maps_for_haplotypes, event_start_positions, Event};
use gatk_engine::fragment::Fragment;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_random::JavaRandom;
use gatk_engine::locus_iterator::AlignmentContext;
use gatk_engine::mutect_engine::log_likelihood_ratio;
use gatk_engine::natural_log_utils::{log_sum_exp, normalize_from_log_to_linear_space, posteriors};
use gatk_engine::pair_hmm_likelihood_engine::LikelihoodEngineArguments;
use gatk_engine::read_pileup::sample_name;
use gatk_engine::read_realignment::realign_reads_to_their_best_haplotype;
use gatk_engine::read_threading_assembler::ReadThreadingAssembler;
use gatk_engine::smith_waterman::SmithWatermanJavaAligner;
use gatk_engine::somatic_likelihoods::{
    allele_fractions_posterior, effective_log_multinomial_weights, log_dirichlet_normalization,
};
use gatk_engine::{read, read_utils};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

use crate::genotyping_engine::EngineError;
use crate::hc_genotyping::phase_calls;
use crate::hc_region::{
    assemble_reads, assembler_from, AssemblerArguments, RegionAssemblyArguments,
};
use crate::reference_confidence_model::{calculate_ref_confidence, RefConfidenceArguments};
use crate::variant_annotator_engine::{Engine as AnnotatorEngine, Site};

/// `Mutect2Engine.MAX_ALT_FRACTION_IN_NORMAL`.
const MAX_ALT_FRACTION_IN_NORMAL: f64 = 0.3;
/// `MAX_NORMAL_QUAL_SUM`.
const MAX_NORMAL_QUAL_SUM: i64 = 100;
/// `MINIMUM_BASE_QUALITY`, for active region determination.
const MINIMUM_BASE_QUALITY: u8 = 6;
/// `INDEL_START_QUAL` and `INDEL_CONTINUATION_QUAL`.
const INDEL_START_QUAL: i32 = 30;
const INDEL_CONTINUATION_QUAL: i32 = 10;
/// `HUGE_FRAGMENT_LENGTH`.
const HUGE_FRAGMENT_LENGTH: i32 = 1_000_000;
/// `AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING`.
const MINIMUM_READ_LENGTH_AFTER_TRIMMING: usize = 10;
/// `QualityUtils.MAX_QUAL`.
const MAX_QUAL: f64 = 93.0;

/// Everything the command line decides for one run.
#[derive(Debug, Clone)]
pub struct Mutect2Arguments {
    pub region: AssemblyRegionArgs,
    /// `PileupDetectionArgumentCollection`.
    pub pileup_detection: crate::pileup_based_alleles::PileupDetectionArguments,
    /// `--alleles`' records, and `--force-call-filtered-alleles`.
    pub given_alleles: Option<(Vec<VariantContext>, bool)>,
    /// `--bam-output` and how it is written.
    pub bamout: Option<crate::haplotype_caller::BamoutArguments>,
    /// `--assembly-region-out`.
    pub assembly_region_out: Option<std::path::PathBuf>,
    pub trimmer: TrimmerArguments,
    pub assembly: RegionAssemblyArguments,
    pub likelihoods: LikelihoodEngineArguments,
    /// `--likelihood-calculation-engine FlowBased`: the flow engine's own arguments.
    pub flow_engine: Option<crate::flow_likelihood_engine::FlowEngineArguments>,
    /// `--likelihood-calculation-engine FlowBasedHMM`: its arguments.
    pub flow_hmm: Option<crate::flow_likelihood_engine::FlowHmmArguments>,
    pub normal_samples: Vec<String>,
    /// `getInitialLogOdds()`, `getEmissionLogOdds()` and `--normal-lod`, natural log.
    pub initial_log_odds: f64,
    pub emission_log_odds: f64,
    pub normal_log_odds: f64,
    pub pcr_snv_qual: i32,
    pub pcr_indel_qual: i32,
    pub multiple_substitution_base_qual_correction: f64,
    pub callable_depth: usize,
    pub genotype_germline_sites: bool,
    /// `genotypeGermlineSites && rng.nextDouble() < genotypeGermlineSitesFraction`, which only a
    /// fraction of one or more makes deterministic: the generator is unseeded.
    pub consider_germline_active: bool,
    pub min_af: f64,
    /// `--flow-assembly-collapse-hmer-size`, -1 resolved from the read groups, and
    /// `--flow-assembly-collapse-partial-mode`.
    pub flow_collapse_hmer_size: i32,
    pub flow_collapse_partial_mode: bool,
    /// `--flow-filter-alleles`, with its thresholds.
    pub filter_alleles: Option<crate::allele_filtering::AlleleFilterArguments>,
    /// `getDefaultAlleleFrequency()`.
    pub default_af: f64,
    /// `--germline-resource`, read whole: `POPAF` from its `AF`, and in tumour-only calling the
    /// loci of its common alleles left inactive.
    pub germline_resource: Vec<VariantContext>,
    /// `--max-population-af`: the `AF` from which a germline allele makes its locus inactive.
    pub max_population_af: f64,
    /// `--panel-of-normals`, read whole: its loci left inactive, a call starting on one flagged
    /// `PON`.
    pub panel_of_normals: Vec<VariantContext>,
    /// `--genotype-pon-sites`: the panel's loci are not left inactive.
    pub genotype_pon_sites: bool,
    pub informative_read_overlap_margin: i32,
    pub max_mnp_distance: i32,
    pub phred_scaled_global_read_mismapping_rate: f64,
    pub max_suspicious_reads_per_alignment_start: i32,
    pub downsampling_stride: i32,
    /// The read-threading assembler's arguments.
    pub assembler: AssemblerArguments,
    /// `--independent-mates`: each read its own fragment.
    pub independent_mates: bool,
    /// `--emit-ref-confidence` other than `NONE`: every region is reported base by base, with
    /// `<NON_REF>` added to every call. The thresholds above are then minus infinity.
    pub emit_reference_confidence: bool,
    /// `SomaticReferenceConfidenceModel`'s prior, from `--minimum-allele-fraction`.
    pub somatic_ref_confidence: crate::reference_confidence_model::SomaticRefConfidence,
}

/// What a run produces: the calls, and the callable-site count the stats file carries.
#[derive(Debug, Clone)]
pub struct Mutect2Output {
    pub calls: Vec<VariantContext>,
    pub callable_sites: u64,
}

/// What a run refuses, with the Java class it would be thrown as.
#[derive(Debug, Clone, PartialEq)]
pub struct Mutect2Error {
    pub class: String,
    pub message: String,
    /// The calls made before the failure, which the reference's writer still writes as it closes.
    pub calls_before: Vec<VariantContext>,
}

fn failure(class: &str, message: impl Into<String>) -> Mutect2Error {
    Mutect2Error {
        class: class.to_string(),
        message: message.into(),
        calls_before: Vec::new(),
    }
}

/// The samples of the header, in `getSamplesFromHeader`'s sorted order.
fn samples_of(header: &SamHeader) -> Vec<String> {
    crate::haplotype_caller::samples_from_header(header)
}

/// The run: `reads` already filtered and in coordinate order.
#[allow(clippy::too_many_arguments)]
pub fn call_variants(
    reads: &[BamRecord],
    header: &SamHeader,
    intervals: &[SimpleInterval],
    arguments: &Mutect2Arguments,
    annotator: &AnnotatorEngine,
    contig_bases: &mut dyn FnMut(&str) -> Result<Vec<u8>, String>,
    random: &mut JavaRandom,
    assembly_variants: Option<
        &std::cell::RefCell<crate::assembly_variants_out::AssemblyVariantsOut>,
    >,
) -> Result<Mutect2Output, Mutect2Error> {
    // `createBamWriter`, closed with the engine, after a failure too.
    let bamout = arguments.bamout.as_ref().map(|b| {
        std::cell::RefCell::new(crate::haplotype_bam_writer::HaplotypeBamWriter::new(
            b.kind, header,
        ))
    });
    let region_out = arguments
        .assembly_region_out
        .as_ref()
        .map(|_| std::cell::RefCell::new(crate::hc_region::ASSEMBLY_REGION_OUT_HEADER.to_string()));
    let outcome = call_variants_with(
        reads,
        header,
        intervals,
        arguments,
        annotator,
        contig_bases,
        random,
        bamout.as_ref(),
        region_out.as_ref(),
        assembly_variants,
    );
    if let (Some(path), Some(text)) = (&arguments.assembly_region_out, region_out) {
        std::fs::write(path, text.into_inner())
            .map_err(|e| failure("UserException.CouldNotCreateOutputFile", e.to_string()))?;
    }
    if let (Some(writer), Some(settings)) = (bamout, &arguments.bamout) {
        crate::haplotype_caller::close_bamout(writer.into_inner(), settings)
            .map_err(|e| failure("UserException.CouldNotCreateOutputFile", e))?;
    }
    outcome
}

/// [`call_variants`], with the bamout writer.
#[allow(clippy::too_many_arguments)]
fn call_variants_with(
    reads: &[BamRecord],
    header: &SamHeader,
    intervals: &[SimpleInterval],
    arguments: &Mutect2Arguments,
    annotator: &AnnotatorEngine,
    contig_bases: &mut dyn FnMut(&str) -> Result<Vec<u8>, String>,
    random: &mut JavaRandom,
    bamout: Option<&std::cell::RefCell<crate::haplotype_bam_writer::HaplotypeBamWriter>>,
    region_out: Option<&std::cell::RefCell<String>>,
    assembly_variants: Option<
        &std::cell::RefCell<crate::assembly_variants_out::AssemblyVariantsOut>,
    >,
) -> Result<Mutect2Output, Mutect2Error> {
    let samples = samples_of(header);
    let partition: Vec<Option<String>> = samples.iter().cloned().map(Some).collect();
    // `MutectReadThreadingAssemblerArgumentCollection.makeReadThreadingAssembler`: adaptive
    // pruning unless disabled, and dangling branches always recovered.
    let assembler = assembler_from(
        &arguments.assembler,
        arguments.assembly.min_base_quality_score,
    )
    .map_err(|e| failure("IllegalArgumentException", format!("{e:?}")))?;
    // `createLikelihoodCalculationEngine`.
    let likelihood_engine = crate::flow_likelihood_engine::ReadLikelihoodEngine::create(
        &arguments.likelihoods,
        &arguments.flow_engine,
        &arguments.flow_hmm,
    )
    .map_err(|(class, message)| failure(&class, message))?;
    let mut tumor_buffer = PileupQualBuffer::default();
    let mut normal_buffer = PileupQualBuffer::default();
    let mut callable_sites = 0u64;
    let mut active_error: Option<Mutect2Error> = None;
    let mut calls = Vec::new();
    for group in group_intervals_by_contig(intervals) {
        let contig = group[0].contig.clone();
        let bases = contig_bases(&contig).map_err(|message| failure("UserException", message))?;
        // `region.getAlignmentData()`, kept only with pileup detection. Mutect2's `isActive` never
        // sets an original active probability, so every locus scores zero.
        let detection = &arguments.pileup_detection;
        let mut loci: Vec<crate::pileup_based_alleles::LocusSummary> = Vec::new();
        // The reads `isActive` gives `MZ`: every read of a non-empty pileup it examines.
        let mut seen: std::collections::HashSet<crate::pileup_based_alleles::ReadKey> =
            std::collections::HashSet::new();
        let regions = traverse_with_pileups(
            reads,
            &group,
            &partition,
            &arguments.region,
            header,
            &mut |reads| {
                mutect_downsample(
                    reads,
                    arguments.region.max_reads_per_alignment_start,
                    arguments.max_suspicious_reads_per_alignment_start,
                    arguments.downsampling_stride,
                    random,
                )
            },
            &mut |locus, context| {
                // A locus a usable given record overlaps is active, and nothing else is counted.
                let forced = arguments
                    .given_alleles
                    .as_ref()
                    .is_some_and(|(records, force)| {
                        crate::given_alleles::forced_at(records, &locus.contig, locus.start, *force)
                    });
                if detection.use_pileup_detection {
                    let ref_base = bases[locus.start as usize - 1];
                    loci.push(match context {
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
                    });
                }
                if forced {
                    return (1.0, None);
                }
                if let Some(context) = context.filter(|_| detection.use_pileup_detection) {
                    for element in &context.pileup.elements {
                        seen.insert(crate::pileup_based_alleles::read_key(element.read));
                    }
                }
                if active_error.is_some() {
                    return (0.0, None);
                }
                let prob = match is_active(
                    context,
                    bases[locus.start as usize - 1],
                    header,
                    arguments,
                    &mut tumor_buffer,
                    &mut normal_buffer,
                    &mut callable_sites,
                ) {
                    Ok(prob) => prob,
                    Err(error) => {
                        active_error = Some(error);
                        0.0
                    }
                };
                (prob, None)
            },
        )
        .map_err(|error| failure(error.class(), error.message()))?;
        // `isActive`'s throw ends the traversal there.
        if let Some(error) = active_error {
            return Err(error);
        }

        let trimmer = AssemblyRegionTrimmer::new(arguments.trimmer.clone(), bases.len() as i32)
            .map_err(|e| failure("IllegalArgumentException", format!("{e:?}")))?;
        let context = RegionContext {
            header,
            samples: &samples,
            contig_bases: &bases,
            arguments,
            assembler: &assembler,
            trimmer: &trimmer,
            likelihood_engine: &likelihood_engine,
            annotator,
            pileup_loci: &loci,
            bamout,
            assembly_variants,
        };
        for mut traversed in regions {
            if detection.use_pileup_detection {
                crate::pileup_based_alleles::tag_mismatch_percentage(
                    &mut traversed.region,
                    &seen,
                    &bases,
                );
            }
            if let Some(out) = region_out {
                out.borrow_mut()
                    .push_str(&crate::hc_region::assembly_region_rows(&traversed.region));
            }
            match call_region(traversed.region, &context, random) {
                Ok(region_calls) => calls.extend(region_calls),
                Err(mut error) => {
                    error.calls_before = calls;
                    return Err(error);
                }
            }
        }
    }
    Ok(Mutect2Output {
        calls,
        callable_sites,
    })
}

/// `PileupQualBuffer`: the qualities behind each kind of alternate at a locus, the four bases,
/// any other substitution, and indels.
#[derive(Default)]
struct PileupQualBuffer {
    buffers: [Vec<i8>; 6],
}

const OTHER_SUBSTITUTION: usize = 4;
const INDEL: usize = 5;

impl PileupQualBuffer {
    /// `accumulateQuals(pileup, refBase, pcrErrorQual)`.
    fn accumulate(
        &mut self,
        elements: &[&gatk_engine::pileup::PileupElement<'_>],
        position: i32,
        ref_base: u8,
        pcr_error_qual: i32,
        correction: f64,
    ) {
        for buffer in &mut self.buffers {
            buffer.clear();
        }
        for pe in elements {
            let indel_length = if pe.is_deletion() {
                pe.current_cigar_element.length as i32
            } else {
                pe.length_of_immediately_following_indel() as i32
            };
            if indel_length > 0 {
                self.buffers[INDEL].push(indel_qual(indel_length));
            } else if is_next_to_useful_soft_clip(pe) {
                self.buffers[INDEL].push(indel_qual(1));
            } else if pe.base() != ref_base && pe.qual() > MINIMUM_BASE_QUALITY {
                let read = pe.read;
                let mate_start = if !read::is_proper_pair(read) || read::mate_is_unmapped(read) {
                    i32::MAX
                } else {
                    read_utils::mate_start(read)
                };
                let overlaps_mate = mate_start <= position
                    && i64::from(position) < i64::from(mate_start) + read.read_bases.len() as i64;
                let qual = if overlaps_mate {
                    (pe.qual() as i32).min(pcr_error_qual / 2) as i8
                } else {
                    pe.qual() as i8
                };
                match simple_base_index(pe.base()) {
                    Some(index) => {
                        self.buffers[index].push((f64::from(qual) + correction).min(MAX_QUAL) as i8)
                    }
                    None => self.buffers[OTHER_SUBSTITUTION].push(qual),
                }
            }
        }
    }

    /// `likeliestIndexAndQuals`: the first kind with the largest positive quality sum, or the
    /// first kind when none has any.
    fn likeliest(&self) -> usize {
        let mut best = 0;
        let mut best_sum = 0i64;
        for n in 0..self.buffers.len() {
            let sum = self.qual_sum(n);
            if sum > best_sum {
                best_sum = sum;
                best = n;
            }
        }
        best
    }

    fn qual_sum(&self, index: usize) -> i64 {
        self.buffers[index].iter().map(|&q| i64::from(q)).sum()
    }
}

/// `BaseUtils.simpleBaseToBaseIndex`.
fn simple_base_index(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// `indelQual(indelLength)`.
fn indel_qual(length: i32) -> i8 {
    (INDEL_START_QUAL + (length - 1) * INDEL_CONTINUATION_QUAL).min(127) as i8
}

/// `isNextToUsefulSoftClip`: a soft clip beside the base whose own neighbour is of some quality.
fn is_next_to_useful_soft_clip(pe: &gatk_engine::pileup::PileupElement<'_>) -> bool {
    let offset = pe.offset;
    let quality_at = |i: i32| {
        usize::try_from(i)
            .ok()
            .and_then(|i| pe.read.base_qualities.get(i))
            .copied()
            .unwrap_or(0)
    };
    pe.qual() > MINIMUM_BASE_QUALITY
        && ((pe.is_before_soft_clip() && quality_at(offset + 1) > MINIMUM_BASE_QUALITY)
            || (pe.is_after_soft_clip() && quality_at(offset - 1) > MINIMUM_BASE_QUALITY))
}

/// `Mutect2Engine.isActive` for a locus, counting the callable sites as it goes.
fn is_active(
    context: Option<&AlignmentContext<'_>>,
    ref_base: u8,
    header: &SamHeader,
    arguments: &Mutect2Arguments,
    tumor: &mut PileupQualBuffer,
    normal: &mut PileupQualBuffer,
    callable_sites: &mut u64,
) -> Result<f64, Mutect2Error> {
    let Some(context) = context.filter(|c| !c.pileup.elements.is_empty()) else {
        return Ok(0.0);
    };
    let pileup = &context.pileup.elements;
    if pileup.len() >= arguments.callable_depth {
        *callable_sites += 1;
    }
    let is_normal = |pe: &&gatk_engine::pileup::PileupElement<'_>| {
        sample_name(pe.read, header).is_some_and(|s| arguments.normal_samples.contains(&s))
    };
    let tumor_pileup: Vec<&gatk_engine::pileup::PileupElement<'_>> =
        pileup.iter().filter(|pe| !is_normal(pe)).collect();
    tumor.accumulate(
        &tumor_pileup,
        context.position,
        ref_base,
        arguments.pcr_snv_qual,
        arguments.multiple_substitution_base_qual_correction,
    );
    let best_tumor = tumor.likeliest();
    let alt_quals: Vec<u8> = tumor.buffers[best_tumor].iter().map(|&q| q as u8).collect();
    let tumor_log_odds =
        log_likelihood_ratio((tumor_pileup.len() - alt_quals.len()) as i32, &alt_quals, 1);
    if tumor_log_odds < arguments.initial_log_odds {
        return Ok(0.0);
    }
    if !arguments.normal_samples.is_empty() && !arguments.genotype_germline_sites {
        let normal_pileup: Vec<&gatk_engine::pileup::PileupElement<'_>> =
            pileup.iter().filter(|pe| is_normal(pe)).collect();
        normal.accumulate(
            &normal_pileup,
            context.position,
            ref_base,
            arguments.pcr_snv_qual,
            arguments.multiple_substitution_base_qual_correction,
        );
        let best_normal = normal.likeliest();
        if best_normal == best_tumor {
            let normal_alt_count = normal.buffers[best_normal].len();
            let normal_qual_sum = normal.qual_sum(best_normal);
            if normal_alt_count as f64 > normal_pileup.len() as f64 * MAX_ALT_FRACTION_IN_NORMAL
                && normal_qual_sum > MAX_NORMAL_QUAL_SUM
            {
                return Ok(0.0);
            }
        }
    } else if !arguments.genotype_germline_sites {
        // Tumour-only: a germline allele at least `--max-population-af` common leaves the locus
        // inactive if it is an indel and the likeliest tumour allele is one, or if it is a
        // substitution and the likeliest tumour base is the germline REFERENCE's first base, which
        // is what the reference compares.
        for germline in overlapping(
            &arguments.germline_resource,
            &context.contig,
            context.position,
            context.position,
        ) {
            let frequencies = attribute_as_double_list(germline, "AF", 0.0)
                .map_err(|message| failure("java.lang.NumberFormatException", message))?;
            let germline_ref = &germline.alleles[0];
            for (index, frequency) in frequencies.iter().enumerate() {
                if *frequency < arguments.max_population_af {
                    continue;
                }
                // `getAlternateAllele(index)`, past the alternates when the AF has more values.
                let Some(germline_alt) = germline.alleles.get(index + 1) else {
                    return Err(failure(
                        "java.lang.IndexOutOfBoundsException",
                        format!(
                            "Index {} out of bounds for length {}",
                            index + 1,
                            germline.alleles.len()
                        ),
                    ));
                };
                let indel = best_tumor == INDEL && germline_alt.len() != germline_ref.len();
                let substitution = best_tumor < OTHER_SUBSTITUTION
                    && germline_alt.len() == germline_ref.len()
                    && Some(&b"ACGT"[best_tumor])
                        == germline_ref.display_string().as_bytes().first();
                if indel || substitution {
                    return Ok(0.0);
                }
            }
        }
    }
    if !arguments.genotype_pon_sites
        && overlapping(
            &arguments.panel_of_normals,
            &context.contig,
            context.position,
            context.position,
        )
        .next()
        .is_some()
    {
        return Ok(0.0);
    }
    if pileup.len() < arguments.callable_depth {
        *callable_sites += 1;
    }
    Ok(1.0)
}

/// `FeatureContext.getValues(input, interval)`: the records overlapping `start..=end`.
fn overlapping<'r>(
    records: &'r [VariantContext],
    contig: &'r str,
    start: i32,
    end: i32,
) -> impl Iterator<Item = &'r VariantContext> {
    records.iter().filter(move |vc| {
        vc.contig == contig && vc.start <= i64::from(end) && vc.stop >= i64::from(start)
    })
}

/// `Mutect2Engine.getAttributeAsDoubleList(vc, key, default)`: each value of the attribute, `.`
/// as the default; nothing when the attribute is absent.
fn attribute_as_double_list(
    vc: &VariantContext,
    key: &str,
    default: f64,
) -> Result<Vec<f64>, String> {
    let Some(text) = vc
        .attributes
        .iter()
        .find(|(name, _)| name == key)
        .and_then(|(_, value)| value.format())
    else {
        return Ok(Vec::new());
    };
    text.split(',')
        .map(|value| {
            if value == "." {
                Ok(default)
            } else {
                crate::variant_recalibrator_model::java_parse_double(value)
                    .ok_or_else(|| format!("For input string: \"{value}\""))
            }
        })
        .collect()
}

/// `Allele.toString()`: the display string, starred for the reference.
fn allele_to_string(allele: &Allele) -> String {
    format!(
        "{}{}",
        allele.display_string(),
        if allele.is_reference() { "*" } else { "" }
    )
}

/// `GATKVariantContextUtils.alleleIndices(alleles1, alleles2)`: for each of `alleles1`, the index
/// in `alleles2` of the allele it is once both are against the longer reference. The reference
/// maps through the `null` both mappings give it, which is index 0 (and so does an allele the
/// mapping drops, such as a symbolic one against a longer reference); two of `alleles2` that land
/// on one key are `Collectors.toMap`'s duplicate key.
fn allele_indices(
    alleles1: &[Allele],
    alleles2: &[Allele],
) -> Result<Vec<Option<usize>>, Mutect2Error> {
    let (ref1, ref2) = (&alleles1[0], &alleles2[0]);
    let common = if ref1.len() < ref2.len() {
        ref2
    } else if ref2.len() < ref1.len() {
        ref1
    } else if ref1 != ref2 {
        return Err(failure(
            "java.lang.IllegalStateException",
            format!(
                "The provided reference alleles do not appear to represent the same position, {} vs. {}",
                allele_to_string(ref1),
                allele_to_string(ref2)
            ),
        ));
    } else {
        ref1
    };
    let mapping = |reference: &Allele, alts: &[Allele]| {
        gatk_engine::allele_mapping::create_allele_mapping(common, reference, alts)
            .map_err(|e| failure("java.lang.IllegalStateException", format!("{e:?}")))
    };
    let to_common1 = mapping(ref1, &alleles1[1..])?;
    let to_common2 = mapping(ref2, &alleles2[1..])?;
    let key = |map: &[(Allele, Allele)], allele: &Allele| {
        map.iter()
            .find(|(from, _)| from == allele)
            .map(|(_, to)| to.clone())
    };
    let mut common_to_index2: Vec<(Option<Allele>, usize)> = Vec::new();
    for (n, allele) in alleles2.iter().enumerate() {
        let k = key(&to_common2, allele);
        if let Some((_, first)) = common_to_index2.iter().find(|(existing, _)| *existing == k) {
            return Err(failure(
                "java.lang.IllegalStateException",
                format!(
                    "Duplicate key {} (attempted merging values {first} and {n})",
                    k.as_ref().map_or("null".to_string(), allele_to_string)
                ),
            ));
        }
        common_to_index2.push((k, n));
    }
    Ok(alleles1
        .iter()
        .map(|allele| {
            let k = key(&to_common1, allele);
            common_to_index2
                .iter()
                .find(|(existing, _)| *existing == k)
                .map(|(_, n)| *n)
        })
        .collect())
}

/// `SomaticGenotypingEngine.getGermlineAltAlleleFrequencies(allAlleles, germlineVCs, default)`:
/// each alternate's frequency, the last germline record matching it winning; a record whose `AF`
/// is not one value per alternate is skipped.
fn germline_alt_allele_frequencies(
    all_alleles: &[Allele],
    germline: &[&VariantContext],
    default: f64,
) -> Result<Vec<f64>, Mutect2Error> {
    let mut frequencies = vec![default; all_alleles.len()];
    for vc in germline {
        let indices = allele_indices(all_alleles, &vc.alleles)?;
        let alt_frequencies = attribute_as_double_list(vc, "AF", default)
            .map_err(|message| failure("java.lang.NumberFormatException", message))?;
        if alt_frequencies.len() != vc.alleles.len() - 1 {
            continue;
        }
        for allele_index in 1..all_alleles.len() {
            if let Some(germline_index) = indices[allele_index] {
                frequencies[allele_index] = match germline_index.checked_sub(1) {
                    Some(i) => alt_frequencies[i],
                    None => {
                        return Err(failure(
                            "java.lang.IndexOutOfBoundsException",
                            format!(
                                "Index -1 out of bounds for length {}",
                                alt_frequencies.len()
                            ),
                        ))
                    }
                };
            }
        }
    }
    Ok(frequencies.split_off(1))
}

/// What every region of a contig is called with.
struct RegionContext<'a> {
    header: &'a SamHeader,
    samples: &'a [String],
    contig_bases: &'a [u8],
    arguments: &'a Mutect2Arguments,
    assembler: &'a ReadThreadingAssembler,
    trimmer: &'a AssemblyRegionTrimmer,
    likelihood_engine: &'a crate::flow_likelihood_engine::ReadLikelihoodEngine,
    annotator: &'a AnnotatorEngine,
    /// The loci walked, with pileup detection.
    pileup_loci: &'a [crate::pileup_based_alleles::LocusSummary],
    /// `--bam-output`.
    bamout: Option<&'a std::cell::RefCell<crate::haplotype_bam_writer::HaplotypeBamWriter>>,
    /// `--debug-assembly-variants-out`'s queue.
    assembly_variants:
        Option<&'a std::cell::RefCell<crate::assembly_variants_out::AssemblyVariantsOut>>,
}

/// A failure while calling: a message naming its own class as `class: message` keeps it.
fn calling(message: impl Into<String>) -> Mutect2Error {
    let message = message.into();
    if let Some((class, rest)) = message.split_once(": ") {
        let qualified = class.starts_with("java.lang.") || class.starts_with("org.broadinstitute.");
        if qualified && !class.contains(' ') {
            return failure(class, rest);
        }
    }
    failure("GATKException", message)
}

/// `Mutect2Engine.MIN_TAIL_QUALITY`.
const MIN_TAIL_QUALITY: u8 = 9;

/// `Mutect2Engine.referenceModelForNoVariation(region)`: nothing outside reference-confidence
/// mode; otherwise the region finalized as Mutect2 finalizes it for this (soft clips dropped, no
/// overlap correction) and a somatic reference-confidence record for each base of its span.
fn reference_model_for_no_variation(
    mut region: AssemblyRegion,
    context: &RegionContext<'_>,
) -> Result<Vec<VariantContext>, Mutect2Error> {
    if !context.arguments.emit_reference_confidence {
        return Ok(Vec::new());
    }
    finalize_region(
        &mut region,
        &FinalizeArguments {
            error_correct_reads: false,
            dont_use_soft_clipped_bases: true,
            min_tail_quality: MIN_TAIL_QUALITY,
            correct_overlapping_base_qualities: false,
            soft_clip_low_quality_ends: false,
            override_softclip_fragment_check: false,
            track_hardclipped_reads: false,
        },
        context.header,
        context.samples,
    )
    .map_err(|e| calling(format!("{e:?}")))?;
    let padded = region.padded_span().clone();
    let ref_bases = &context.contig_bases[padded.start as usize - 1..padded.end as usize];
    somatic_ref_confidence(
        ref_bases,
        region.span(),
        padded.start,
        region.reads(),
        &[],
        context,
    )
}

/// `calculateRefConfidence` under the somatic model, which takes exactly one sample.
fn somatic_ref_confidence(
    ref_bases: &[u8],
    span: &SimpleInterval,
    padded_start: i32,
    evidence: &[BamRecord],
    calls: &[VariantContext],
    context: &RegionContext<'_>,
) -> Result<Vec<VariantContext>, Mutect2Error> {
    if context.samples.len() != 1 {
        return Err(failure(
            "java.lang.IllegalArgumentException",
            format!(
                "readLikelihoods must contain exactly one sample but it contained {}",
                context.samples.len()
            ),
        ));
    }
    calculate_ref_confidence(
        ref_bases,
        span,
        padded_start,
        evidence,
        &context.samples[0],
        2,
        calls,
        context.header,
        &RefConfidenceArguments {
            flow_based_model: context.likelihood_engine.is_flow_based(),
            ..RefConfidenceArguments::default()
        },
        Some(&context.arguments.somatic_ref_confidence),
        None,
    )
    .map_err(|e| match e {
        EngineError::Limitation(what) => failure("PortLimitation", what),
        EngineError::Runtime { class, message } => failure(&class, message),
    })
}

/// `Mutect2Engine.callRegion`.
fn call_region(
    mut region: AssemblyRegion,
    context: &RegionContext<'_>,
    random: &mut JavaRandom,
) -> Result<Vec<VariantContext>, Mutect2Error> {
    let arguments = context.arguments;
    // `cleanOverlappingReadPairs(..., false, pcrSnvQual / 2, pcrIndelQual / 2)`.
    let mut reads = region.reads().to_vec();
    clean_overlapping_read_pairs_with_indels(
        &mut reads,
        context.samples,
        context.header,
        false,
        Some((arguments.pcr_snv_qual / 2) as u8),
        Some((arguments.pcr_indel_qual / 2) as u8),
    )
    .map_err(|e| calling(format!("{e:?}")))?;
    region.clear_reads();
    region
        .add_all(reads, context.header)
        .map_err(|e| calling(format!("{e:?}")))?;
    if !region.is_active() || region.reads().is_empty() {
        return reference_model_for_no_variation(region, context);
    }
    remove_unmarked_duplicates(&mut region, context.header)?;
    // The given events: the usable records over the padded span, split, wherever they start.
    let given = match &arguments.given_alleles {
        Some((records, force_call_filtered)) => crate::given_alleles::given_events(
            records,
            &region.span().contig,
            region.padded_span(),
            i32::MIN,
            *force_call_filtered,
        )
        .map_err(calling)?,
        None => Vec::new(),
    };

    let mut untrimmed = assemble_reads(
        &mut region,
        &arguments.assembly,
        context.header,
        context.samples,
        context.contig_bases,
        context.assembler,
        &SmithWatermanJavaAligner,
    )
    .map_err(|e| match e {
        // An exception the reference throws as it is, `class: message`.
        crate::hc_region::RegionAssemblyError::Calling(message)
            if message.split_once(": ").is_some_and(|(class, _)| {
                (class.starts_with("org.broadinstitute.") || class.starts_with("java.lang."))
                    && !class.contains(' ')
            }) =>
        {
            let (class, rest) = message.split_once(": ").expect("just checked");
            failure(class, rest)
        }
        e => calling(format!("{e:?}")),
    })?;
    // `addAssembledVariantsToEventMapOutput`.
    if let Some(queue) = context.assembly_variants {
        for event in untrimmed
            .variation_events(arguments.max_mnp_distance)
            .map_err(|e| calling(format!("{e:?}")))?
        {
            queue.borrow_mut().add_event(&event);
        }
    }
    let mut events = untrimmed
        .variation_events(arguments.max_mnp_distance)
        .map_err(|e| calling(format!("{e:?}")))?;
    // `goodAndBadPileupEvents`: the good events join the trimmer's, the bad ones are not taken
    // out of them.
    let (good_pileup_events, _) = crate::pileup_based_alleles::good_and_bad_pileup_events(
        crate::pileup_based_alleles::region_loci(context.pileup_loci, &region),
        &arguments.pileup_detection,
        context.contig_bases,
    )
    .map_err(calling)?;
    if !good_pileup_events.is_empty() {
        events =
            crate::pileup_based_alleles::events_for_trimming(&events, &good_pileup_events, &[]);
    }
    if !given.is_empty() {
        events = crate::pileup_based_alleles::events_for_trimming(&events, &given, &[]);
    }
    let padded = region.padded_span().clone();
    let padded_bases = &context.contig_bases[padded.start as usize - 1..padded.end as usize];
    let trimmed = context
        .trimmer
        .trim(&region, &events, padded_bases, padded.start)
        .map_err(|e| calling(format!("{e:?}")))?;
    if trimmed.variant_span.is_none() {
        return reference_model_for_no_variation(region, context);
    }
    let variant_region = context
        .trimmer
        .variant_region(&trimmed, &region, context.header)
        .map_err(|e| calling(format!("{e:?}")))?;
    let mut assembly = untrimmed
        .trim_to(variant_region)
        .map_err(|e| calling(format!("{e:?}")))?;
    // `addGivenAlleles`, then `injectPileupEvents`, after a `removeHaplotypesWithBadAlleles`
    // that the trimmed set's missing event maps leave without effect.
    crate::given_alleles::add_given_alleles(
        &mut assembly,
        &given,
        arguments.max_mnp_distance,
        &SmithWatermanJavaAligner,
        &arguments.assembly.haplotype_to_reference_sw,
    )
    .map_err(calling)?;
    crate::pileup_based_alleles::inject_pileup_events(
        &mut assembly,
        &region,
        &arguments.pileup_detection,
        &SmithWatermanJavaAligner,
        &arguments.assembly.haplotype_to_reference_sw,
        &good_pileup_events,
        arguments.max_mnp_distance,
    )
    .map_err(calling)?;
    if !assembly.is_variation_present() {
        return reference_model_for_no_variation(region, context);
    }
    // `removeReadStubs`: reads under ten bases long, soft clips included.
    let genotyping_region = assembly
        .region_for_genotyping_mut()
        .expect("a trimmed set has its region");
    genotyping_region.retain_reads(|r| r.read_bases.len() >= MINIMUM_READ_LENGTH_AFTER_TRIMMING);
    let genotyping_region = genotyping_region.clone();

    let mut by_sample: Vec<Vec<BamRecord>> = vec![Vec::new(); context.samples.len()];
    for read in genotyping_region.reads() {
        let sample = sample_name(read, context.header);
        let index = context
            .samples
            .iter()
            .position(|s| Some(s) == sample.as_ref())
            .ok_or_else(|| calling("a read of no known sample"))?;
        by_sample[index].push(read.clone());
    }
    let mut haplotypes: Vec<Haplotype> = assembly.haplotype_list().to_vec();
    let mut likelihoods = context
        .likelihood_engine
        .compute_read_likelihoods(&haplotypes, context.header, context.samples, &by_sample)
        .map_err(calling)?;
    likelihoods
        .switch_to_natural_log()
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
        &SmithWatermanJavaAligner,
        &arguments.assembly.read_to_haplotype_sw,
    )
    .map_err(|e| calling(format!("{e:?}")))?;
    // The collapsing engine `assembleReads` built: a threshold above zero and a run longer than it
    // under the padded span.
    let threshold = arguments.flow_collapse_hmer_size;
    if threshold > 0
        && crate::haplotype_collapsing::needs_collapsing(
            &context.contig_bases[padded.start as usize - 1..padded.end as usize],
            threshold,
        )
    {
        let engine = crate::haplotype_collapsing::HaplotypeCollapsing {
            hmer_size_threshold: threshold,
            partial_mode: arguments.flow_collapse_partial_mode,
            full_ref: assembly.full_reference_with_padding().to_vec(),
            ref_loc: ref_loc.clone(),
            parameters: arguments.assembly.haplotype_to_reference_sw,
        };
        crate::hc_region::uncollapse(
            &engine,
            &mut likelihoods,
            &mut haplotypes,
            &SmithWatermanJavaAligner,
        )
        .map_err(calling)?;
    }
    // `uncollapsedReadLikelihoods.alleles()`, which the bamout writes.
    let initial_haplotypes = haplotypes.clone();
    // `AlleleFilteringMutect`, on the natural-log matrix.
    if let Some(filter) = &arguments.filter_alleles {
        let alt_pseudocount = if arguments.min_af == 0.0 {
            1.0
        } else {
            1.0 - std::f64::consts::LN_2 / jmath::math::log(arguments.min_af)
        };
        let filtered = crate::allele_filtering::filter_alleles(
            &likelihoods,
            assembly.full_reference_with_padding(),
            &ref_loc,
            arguments.max_mnp_distance,
            filter,
            crate::allele_filtering::FilterModel::Somatic {
                ref_pseudocount: 1.0,
                alt_pseudocount,
            },
        )
        .map_err(calling)?;
        likelihoods = filtered.likelihoods;
        likelihoods.set_filtered_haplotype_count(filtered.filtered_haplotype_count);
        haplotypes = (0..likelihoods.number_of_alleles())
            .filter_map(|index| likelihoods.get_allele(index).cloned())
            .collect();
    }
    let window = (
        i64::from(padded.start),
        &context.contig_bases[padded.start as usize - 1..padded.end as usize],
    );
    let realigned: Vec<BamRecord> = likelihoods.sample_evidence(0).unwrap_or(&[]).to_vec();
    // The bamout reads the matrix as genotyping receives it, before any normalization.
    let mut bamout_tags = context.bamout.map(|_| {
        (
            crate::hc_genotyping::BamoutTags::default(),
            likelihoods.clone(),
        )
    });
    let calls = call_mutations(
        likelihoods,
        &mut haplotypes,
        assembly.full_reference_with_padding(),
        &ref_loc,
        genotyping_region.span(),
        window,
        context,
        random,
        &given,
        bamout_tags.as_mut().map(|(tags, _)| tags),
    )?;
    // `writeBamOutput(uncollapsedReadLikelihoods.alleles(), paddedReferenceLoc,
    // subsettedReadLikelihoodsFinal, calledHaplotypes, regionForGenotyping.getSpan())`.
    if let (Some(writer), Some((tags, scored))) = (context.bamout, &bamout_tags) {
        let called: Vec<Haplotype> = tags.called.iter().map(|&h| haplotypes[h].clone()).collect();
        writer
            .borrow_mut()
            .write_region(
                &initial_haplotypes,
                &ref_loc,
                &called,
                None,
                scored,
                genotyping_region.span(),
                &|read| tags.supported_for(read).map(str::to_string),
            )
            .map_err(calling)?;
    }
    if !arguments.emit_reference_confidence {
        return Ok(calls);
    }
    // `containsCalls`: a genotype with any allele called.
    let contains_calls = calls.iter().any(|call| {
        call.genotypes
            .iter()
            .any(|g| g.alleles.iter().any(|a| !a.is_no_call()))
    });
    if !contains_calls {
        return reference_model_for_no_variation(region, context);
    }
    // The left flank, the genotyping region base by base with its calls, then the right flank.
    let left = context
        .trimmer
        .non_variant_left_flank_region(&trimmed, &region, context.header)
        .map_err(|e| calling(format!("{e:?}")))?;
    let right = context
        .trimmer
        .non_variant_right_flank_region(&trimmed, &region, context.header)
        .map_err(|e| calling(format!("{e:?}")))?;
    let mut result = Vec::new();
    if let Some(flank) = left {
        result.extend(reference_model_for_no_variation(flank, context)?);
    }
    result.extend(somatic_ref_confidence(
        &reference_haplotype.bases(),
        genotyping_region.span(),
        // The trimmed reference haplotype spans the genotyping region's padded span.
        genotyping_region.padded_span().start,
        &realigned,
        &calls,
        context,
    )?);
    if let Some(flank) = right {
        result.extend(reference_model_for_no_variation(flank, context)?);
    }
    Ok(result)
}

/// `removeUnmarkedDuplicates`: among paired reads whose mate maps far away or to another contig,
/// grouped by sample and signed unclipped start, every read of a minority mate contig and all but
/// the first of a majority one.
fn remove_unmarked_duplicates(
    region: &mut AssemblyRegion,
    header: &SamHeader,
) -> Result<(), Mutect2Error> {
    let mut groups: Vec<(DuplicateKey, Vec<usize>)> = Vec::new();
    for (i, r) in region.reads().iter().enumerate() {
        let candidate = read::is_paired(r)
            && !read::mate_is_unmapped(r)
            && (r.mate_reference_index != r.reference_index
                || read::fragment_length(r).abs() > HUGE_FRAGMENT_LENGTH);
        if !candidate {
            continue;
        }
        let sign = if read::is_first_of_pair(r) { 1 } else { -1 };
        let key = (
            sample_name(r, header),
            sign * read_utils::unclipped_start(r),
        );
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some((_, list)) => list.push(i),
            None => groups.push((key, vec![i])),
        }
    }
    let mut duplicates: Vec<usize> = Vec::new();
    for (_, list) in &groups {
        let mut by_contig: Vec<(i32, Vec<usize>)> = Vec::new();
        for &i in list {
            let contig = region.reads()[i].mate_reference_index;
            match by_contig.iter_mut().find(|(c, _)| *c == contig) {
                Some((_, members)) => members.push(i),
                None => by_contig.push((contig, vec![i])),
            }
        }
        for (_, members) in by_contig {
            let skip = usize::from(members.len() > list.len() / 2);
            duplicates.extend(members.into_iter().skip(skip));
        }
    }
    if duplicates.is_empty() {
        return Ok(());
    }
    let kept: Vec<BamRecord> = region
        .reads()
        .iter()
        .enumerate()
        .filter(|(i, _)| !duplicates.contains(i))
        .map(|(_, r)| r.clone())
        .collect();
    let header_owned = header;
    region.clear_reads();
    region
        .add_all(kept, header_owned)
        .map_err(|e| calling(format!("{e:?}")))
}

/// A duplicate group's key: the sample and the signed unclipped start.
type DuplicateKey = (Option<String>, i32);

/// A likelihood matrix as `[allele][evidence]`.
pub(crate) type Matrix = Vec<Vec<f64>>;

/// `SomaticGenotypingEngine.callMutations`, with no given
/// alleles, germline resource or panel of normals.
#[allow(clippy::too_many_arguments)]
fn call_mutations(
    mut log_read_likelihoods: AlleleLikelihoods<BamRecord, Haplotype>,
    haplotypes: &mut [Haplotype],
    reference: &[u8],
    ref_loc: &SimpleInterval,
    window: &SimpleInterval,
    reference_window: (i64, &[u8]),
    context: &RegionContext<'_>,
    random: &mut JavaRandom,
    given: &[Event],
    mut bamout: Option<&mut crate::hc_genotyping::BamoutTags>,
) -> Result<Vec<VariantContext>, Mutect2Error> {
    let arguments = context.arguments;
    let illegal = |e: String| failure("IllegalArgumentException", e);
    build_event_maps_for_haplotypes(haplotypes, reference, ref_loc, arguments.max_mnp_distance)
        .map_err(|e| illegal(e.message()))?;
    let starts: Vec<i32> = event_start_positions(haplotypes)
        .map_err(|e| illegal(e.message()))?
        .into_iter()
        .filter(|&loc| window.start <= loc && loc <= window.end)
        .collect();
    if arguments.phred_scaled_global_read_mismapping_rate > 0.0 {
        let cap =
            -arguments.phred_scaled_global_read_mismapping_rate * std::f64::consts::LN_10 / 10.0;
        log_read_likelihoods
            .normalize_likelihoods(cap, true)
            .map_err(|e| illegal(format!("{e:?}")))?;
    }
    let fragments = log_read_likelihoods
        .group_by_fragment_with(arguments.independent_mates)
        .map_err(|e| illegal(format!("{e:?}")))?;
    let has_normal = !arguments.normal_samples.is_empty();
    let contig_length = context.contig_bases.len() as i32;
    let alt_pseudocount = if arguments.min_af == 0.0 {
        1.0
    } else {
        1.0 - std::f64::consts::LN_2 / jmath::math::log(arguments.min_af)
    };

    let mut potential: Vec<Event> = Vec::new();
    let mut called: Vec<usize> = Vec::new();
    let mut calls: Vec<VariantContext> = Vec::new();
    for loc in starts {
        let events = variants_from_active_haplotypes(loc, haplotypes, false)
            .map_err(|e| illegal(format!("{e:?}")))?;
        let Some(merged) =
            make_merged_variant_context(&events).map_err(|e| illegal(format!("{e:?}")))?
        else {
            continue;
        };
        // `ReferenceConfidenceUtils.addNonRefSymbolicAllele`: `<NON_REF>` joins the merged context.
        // The mapper skips symbolic alleles, so the marginalized matrices lack it until
        // `addNonReferenceAllele` adds its row, each read's median of the non-best alleles.
        let mut merged = merged;
        if arguments.emit_reference_confidence {
            merged
                .alleles
                .push(Allele::from_str("<NON_REF>", false).expect("a symbolic allele"));
        }
        let mapper = create_allele_mapper(&merged, loc, haplotypes, true)
            .map_err(|e| illegal(format!("{e:?}")))?;
        let new_to_old: Vec<(Allele, Vec<Haplotype>)> = mapper
            .iter()
            .map(|(allele, indices)| {
                (
                    allele.clone(),
                    indices.iter().map(|&h| haplotypes[h].clone()).collect(),
                )
            })
            .collect();
        let mut log_likelihoods = fragments
            .marginalize(&new_to_old)
            .map_err(|e| illegal(format!("{e:?}")))?;
        let overlap = SimpleInterval::new(&merged.contig, merged.start, merged.end)
            .and_then(|span| {
                span.expand_within_contig(arguments.informative_read_overlap_margin, contig_length)
            })
            .ok_or_else(|| illegal("no overlap interval".to_string()))?;
        let contig_index = context
            .header
            .sequences
            .iter()
            .position(|s| s.name == overlap.contig)
            .map(|i| i as i32);
        log_likelihoods.retain_evidence(|f: &Fragment| {
            Some(f.contig_index) == contig_index && f.start <= overlap.end && overlap.start <= f.end
        });
        if arguments.emit_reference_confidence {
            log_likelihoods.add_non_reference_allele();
        }

        let alleles: Vec<Allele> = (0..log_likelihoods.number_of_alleles())
            .map(|a| log_likelihoods.get_allele(a).expect("an allele").clone())
            .collect();
        let ref_index = alleles
            .iter()
            .position(Allele::is_reference)
            .ok_or_else(|| illegal("No ref allele found in likelihoods".to_string()))?;
        let is_normal = |s: usize| {
            let name = log_likelihoods.get_sample(s).expect("a sample");
            arguments.normal_samples.contains(name)
        };
        let samples = log_likelihoods.number_of_samples();
        let tumor_matrix =
            combined_matrix(&log_likelihoods, (0..samples).filter(|&s| !is_normal(s)));
        let normal_matrix =
            combined_matrix(&log_likelihoods, (0..samples).filter(|&s| is_normal(s)));
        let tumor_log_odds = somatic_log_odds(&tumor_matrix, ref_index, alt_pseudocount)?;
        let normal_log_odds = diploid_alt_log_odds(&normal_matrix, ref_index)?;
        let normal_artifact_log_odds =
            somatic_log_odds(&normal_matrix, ref_index, alt_pseudocount)?;

        let normal_threshold = arguments.normal_log_odds;
        // `allelesConsistentWithGivenAlleles(givenAlleles, mergedVC)`: kept whatever their odds.
        let forced = crate::given_alleles::alleles_consistent_with_given_alleles(
            given,
            &merged.alleles[0],
            &merged.alleles[1..],
        );
        let is_forced = |a: usize| forced.contains(&alleles[a]);
        let tumor_alts: Vec<usize> = merged_alternates(&merged.alleles, &alleles)
            .into_iter()
            .filter(|&a| is_forced(a) || tumor_log_odds[&a] > arguments.emission_log_odds)
            .collect();
        let to_genotype: Vec<usize> = tumor_alts
            .iter()
            .copied()
            .filter(|&a| {
                is_forced(a)
                    || !has_normal
                    || arguments.consider_germline_active
                    || normal_log_odds[&a] > normal_threshold
            })
            .collect();
        for &a in &merged_alternates(&merged.alleles, &alleles) {
            if tumor_log_odds[&a] > arguments.emission_log_odds
                && (!has_normal || normal_log_odds[&a] > normal_threshold)
            {
                let event = Event::new(
                    &merged.contig,
                    merged.start,
                    alleles[ref_index].clone(),
                    alleles[a].clone(),
                )
                .map_err(|e| illegal(e.message()))?;
                if !potential.contains(&event) {
                    potential.push(event);
                }
            }
        }
        if to_genotype.is_empty() {
            continue;
        }

        let emitted: Vec<usize> = std::iter::once(ref_index)
            .chain(tumor_alts.iter().copied())
            .collect();
        let emitted_alleles: Vec<Allele> = emitted.iter().map(|&a| alleles[a].clone()).collect();
        let mut vc = VariantContext::new(
            &merged.contig,
            i64::from(merged.start),
            emitted_alleles.clone(),
        );
        vc.stop = i64::from(merged.end);
        // `MathUtils.logToLog10`: `ln * LOG10_E`, which is not always `ln / ln(10)` in the last bit.
        let log10 = |x: f64| x * jmath::math::log10(std::f64::consts::E);
        let doubles =
            |values: Vec<f64>| Value::List(values.into_iter().map(Value::Double).collect());
        // `getNegativeLogPopulationAFAnnotation`, from the germline records overlapping the start.
        let germline: Vec<&VariantContext> = overlapping(
            &arguments.germline_resource,
            &merged.contig,
            merged.start,
            merged.start,
        )
        .collect();
        let frequencies =
            germline_alt_allele_frequencies(&emitted_alleles, &germline, arguments.default_af)?;
        vc.attributes.push((
            "POPAF".to_string(),
            doubles(
                frequencies
                    .into_iter()
                    .map(|f| -jmath::math::log10(f))
                    .collect(),
            ),
        ));
        vc.attributes.push((
            "TLOD".to_string(),
            doubles(
                tumor_alts
                    .iter()
                    .map(|a| log10(tumor_log_odds[a]))
                    .collect(),
            ),
        ));
        if has_normal {
            vc.attributes.push((
                "NALOD".to_string(),
                doubles(
                    tumor_alts
                        .iter()
                        .map(|a| log10(normal_artifact_log_odds[a]))
                        .collect(),
                ),
            ));
            vc.attributes.push((
                "NLOD".to_string(),
                doubles(
                    tumor_alts
                        .iter()
                        .map(|a| log10(normal_log_odds[a]))
                        .collect(),
                ),
            ));
        }
        if overlapping(
            &arguments.panel_of_normals,
            &merged.contig,
            merged.start,
            merged.start,
        )
        .next()
        .is_some()
        {
            vc.attributes.push(("PON".to_string(), Value::Bool(true)));
        }
        vc.genotypes = GenotypesContext::new(genotypes(
            &log_likelihoods,
            &emitted,
            &alleles,
            ref_index,
            &arguments.normal_samples,
        )?);

        let (trimmed_call, trim_map) = trim_alleles(&vc).map_err(illegal)?;
        // `Collectors.toMap`, a `HashMap`: the marginalised matrices list their alleles in
        // `Allele.hashCode` order, not the call's, and an annotation that reads the matrix by row
        // (`AllelePseudoDepth`, `BaseQualityHistogram`) writes them in that order.
        let order = gatk_annotation::depth_per_allele::marginalisation_order(
            &trim_map.iter().map(|(a, _)| a.clone()).collect::<Vec<_>>(),
        );
        let trim_map: Vec<(Allele, Vec<Allele>)> = order
            .iter()
            .filter_map(|allele| trim_map.iter().find(|(a, _)| a == allele).cloned())
            .collect();
        let trimmed_fragments = log_likelihoods
            .marginalize(&trim_map)
            .map_err(|e| illegal(format!("{e:?}")))?;
        let mut read_alleles = log_read_likelihoods
            .marginalize(&new_to_old)
            .map_err(|e| illegal(format!("{e:?}")))?;
        read_alleles.retain_evidence(|r: &BamRecord| {
            overlap.overlaps(&overlap.contig, read_utils::start(r), read_utils::end(r))
        });
        if arguments.emit_reference_confidence {
            read_alleles.add_non_reference_allele();
        }
        let trimmed_reads = read_alleles
            .marginalize(&trim_map)
            .map_err(|e| illegal(format!("{e:?}")))?;
        let site = Site {
            likelihoods: &trimmed_reads,
            window: reference_window,
            overlaps: Vec::new(),
            dbsnp: None,
            resources: Vec::new(),
            // `Optional.of(logFragmentLikelihoods)`, the region's fragment-by-haplotype matrix.
            haplotypes: Some(crate::variant_annotator_engine::HaplotypeMatrix::Fragments(
                &fragments,
                &*haplotypes,
            )),
        };
        let annotated = context
            .annotator
            .annotate_context_somatic(&trimmed_call, &site, &trimmed_fragments, random)
            .map_err(|e| match e {
                EngineError::Limitation(what) => failure("PortLimitation", what),
                EngineError::Runtime { class, message } => failure(&class, message),
            })?;
        // `annotateReadLikelihoodsWithSupportedAlleles(trimmedCall, trimmedLikelihoods,
        // Fragment::getReads)`, for the bamout.
        if let Some(tags) = bamout.as_deref_mut() {
            tags.add_supported_alleles_fragments(&trimmed_call, &trimmed_fragments)
                .map_err(|e| illegal(format!("{e:?}")))?;
        }
        for allele in &vc.alleles {
            if let Some((_, indices)) = mapper.iter().find(|(a, _)| a == allele) {
                for &h in indices {
                    if !called.iter().any(|&c| haplotypes[c] == haplotypes[h]) {
                        called.push(h);
                    }
                }
            }
        }
        calls.push(annotated);
    }

    if let Some(tags) = bamout {
        tags.called = called.clone();
    }
    let called_haplotypes: Vec<&Haplotype> = called.iter().map(|&h| &haplotypes[h]).collect();
    let output = phase_calls(calls, &called_haplotypes).map_err(|e| match e {
        EngineError::Limitation(what) => failure("PortLimitation", what),
        EngineError::Runtime { class, message } => failure(&class, message),
    })?;

    // The haplotype each read supports best, counted.
    let mut support = vec![0usize; log_read_likelihoods.number_of_alleles()];
    for best in log_read_likelihoods.best_alleles_breaking_ties(None) {
        if let Some(allele) = &best.allele {
            if let Some(index) = log_read_likelihoods.index_of_allele(allele) {
                support[index] += 1;
            }
        }
    }
    let likelihood_haplotypes: Vec<&Haplotype> = (0..log_read_likelihoods.number_of_alleles())
        .map(|a| {
            let allele = log_read_likelihoods.get_allele(a).expect("a haplotype");
            haplotypes
                .iter()
                .find(|h| *h == allele)
                .expect("the likelihoods' haplotypes are the assembly's")
        })
        .collect();
    let mut annotated = Vec::with_capacity(output.len());
    for mut call in output {
        let mut counts: Vec<i64> = Vec::new();
        for alt in call.alternate_alleles() {
            let event = Event::new(
                &call.contig,
                call.start as i32,
                call.alleles[0].clone(),
                alt.clone(),
            )
            .map_err(|e| illegal(e.message()))?;
            let mut carrying: Vec<usize> = (0..likelihood_haplotypes.len())
                .filter(|&h| {
                    likelihood_haplotypes[h]
                        .event_map()
                        .is_some_and(|map| map.events().any(|e| *e == event))
                })
                .collect();
            if carrying.is_empty() {
                continue;
            }
            carrying.sort_by(|a, b| support[*b].cmp(&support[*a]));
            let best = likelihood_haplotypes[carrying[0]];
            let count = best
                .event_map()
                .map(|map| map.events().filter(|e| potential.contains(e)).count())
                .unwrap_or(0);
            counts.push(count as i64);
        }
        call.attributes.push((
            "ECNTH".to_string(),
            Value::List(counts.into_iter().map(Value::Int).collect()),
        ));
        call.attributes
            .push(("ECNT".to_string(), Value::Int(potential.len() as i64)));
        annotated.push(call);
    }
    Ok(annotated)
}

/// The alternates of the merged context, as indices into the matrix's alleles.
fn merged_alternates(merged: &[Allele], alleles: &[Allele]) -> Vec<usize> {
    merged
        .iter()
        .filter(|a| !a.is_reference())
        .filter_map(|a| alleles.iter().position(|b| b == a))
        .collect()
}

/// `combinedLikelihoodMatrix`: the samples' evidence side by side.
fn combined_matrix(
    likelihoods: &AlleleLikelihoods<Fragment>,
    samples: impl Iterator<Item = usize>,
) -> Matrix {
    let alleles = likelihoods.number_of_alleles();
    let mut matrix: Matrix = vec![Vec::new(); alleles];
    for s in samples {
        for r in 0..likelihoods.sample_evidence_count(s) {
            for (a, row) in matrix.iter_mut().enumerate() {
                row.push(likelihoods.value(s, a, r));
            }
        }
    }
    matrix
}

fn column(matrix: &Matrix, r: usize) -> Vec<f64> {
    matrix.iter().map(|row| row[r]).collect()
}

fn evidence_count(matrix: &Matrix) -> usize {
    matrix.first().map_or(0, Vec::len)
}

fn numerical(e: impl std::fmt::Debug) -> Mutect2Error {
    failure("IllegalArgumentException", format!("{e:?}"))
}

/// `SomaticLikelihoodsEngine.logEvidence`.
pub(crate) fn log_evidence(matrix: &Matrix, prior: &[f64]) -> Result<f64, Mutect2Error> {
    let posterior = allele_fractions_posterior(matrix, prior, None)
        .map_err(numerical)?
        .values;
    let prior_contribution = log_dirichlet_normalization(prior);
    let posterior_contribution = -log_dirichlet_normalization(&posterior);
    let log_fractions =
        effective_log_multinomial_weights(&posterior).ok_or_else(|| numerical("digamma"))?;
    let mut rest = 0.0;
    for r in 0..evidence_count(matrix) {
        let read = column(matrix, r);
        let responsibilities =
            posteriors(&log_fractions, &read).ok_or_else(|| numerical("posteriors"))?;
        let entropy: f64 = responsibilities
            .iter()
            .map(|&x| {
                if x < 1e-8 {
                    0.0
                } else {
                    x * jmath::math::log(x)
                }
            })
            .fold(0.0, |acc, v| acc + v);
        let mut likelihood = 0.0;
        for (l, p) in read.iter().zip(&responsibilities) {
            likelihood += if *p < 1.0e-10 { 0.0 } else { l * p };
        }
        rest += likelihood - entropy;
    }
    Ok(prior_contribution + posterior_contribution + rest)
}

/// `somaticLogOdds`: per alternate, the log evidence with every allele less the log evidence
/// without that one, an empty matrix counting zero.
fn somatic_log_odds(
    matrix: &Matrix,
    ref_index: usize,
    alt_pseudocount: f64,
) -> Result<HashMap<usize, f64>, Mutect2Error> {
    let pseudocounts = |n: usize| -> Vec<f64> {
        (0..n)
            .map(|i| if i == 0 { 1.0 } else { alt_pseudocount })
            .collect()
    };
    let empty = evidence_count(matrix) == 0;
    let with_all = if empty {
        0.0
    } else {
        log_evidence(matrix, &pseudocounts(matrix.len()))?
    };
    let mut lods = HashMap::new();
    for a in (0..matrix.len()).filter(|&a| a != ref_index) {
        let without: Matrix = matrix
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != a)
            .map(|(_, row)| row.clone())
            .collect();
        let with_out = if empty {
            0.0
        } else {
            log_evidence(&without, &pseudocounts(without.len()))?
        };
        lods.insert(a, with_all - with_out);
    }
    Ok(lods)
}

/// `diploidAltLogOdds`: per alternate, the hom-ref log likelihood less the het one.
fn diploid_alt_log_odds(
    matrix: &Matrix,
    ref_index: usize,
) -> Result<HashMap<usize, f64>, Mutect2Error> {
    let reads = evidence_count(matrix);
    let mut hom_ref = 0.0;
    for value in &matrix[ref_index][..reads] {
        hom_ref += value;
    }
    let log_one_half = -std::f64::consts::LN_2;
    let mut result = HashMap::new();
    for a in (0..matrix.len()).filter(|&a| a != ref_index) {
        let mut het = 0.0;
        for (reference, alternate) in matrix[ref_index][..reads].iter().zip(&matrix[a][..reads]) {
            het += log_sum_exp(&[*reference, *alternate]).map_err(numerical)? + log_one_half;
        }
        result.insert(a, hom_ref - het);
    }
    Ok(result)
}

/// `addGenotypes`: per sample, the emitted alleles' effective counts as `AD` and the posterior
/// mean allele fractions as `AF`; a normal is hom-ref, a tumor carries every emitted allele.
fn genotypes(
    likelihoods: &AlleleLikelihoods<Fragment>,
    emitted: &[usize],
    alleles: &[Allele],
    ref_index: usize,
    normal_samples: &[String],
) -> Result<Vec<Genotype>, Mutect2Error> {
    let mut out = Vec::new();
    for s in 0..likelihoods.number_of_samples() {
        let sample = likelihoods.get_sample(s).expect("a sample").clone();
        let evidence = likelihoods.sample_evidence_count(s);
        let matrix: Matrix = emitted
            .iter()
            .map(|&a| (0..evidence).map(|r| likelihoods.value(s, a, r)).collect())
            .collect();
        let mut counts = vec![0.0; emitted.len()];
        for r in 0..evidence {
            let normalized =
                normalize_from_log_to_linear_space(&column(&matrix, r)).map_err(numerical)?;
            for (slot, value) in counts.iter_mut().zip(normalized) {
                *slot += value;
            }
        }
        let flat = vec![1.0; emitted.len()];
        let posterior = if evidence == 0 {
            flat.clone()
        } else {
            allele_fractions_posterior(&matrix, &flat, None)
                .map_err(numerical)?
                .values
        };
        let total: f64 = posterior.iter().sum();
        let means: Vec<f64> = posterior.iter().map(|p| p / total).collect();
        let reference = alleles[ref_index].clone();
        let genotype_alleles = if normal_samples.contains(&sample) {
            vec![reference.clone(), reference]
        } else {
            emitted.iter().map(|&a| alleles[a].clone()).collect()
        };
        let mut genotype = Genotype::new(&sample, genotype_alleles);
        genotype.ad = Some(
            counts
                .iter()
                .map(|&x| jmath::fast_math::round(x) as i32)
                .collect(),
        );
        genotype.extended.push((
            "AF".to_string(),
            Value::List(means[1..].iter().map(|&m| Value::Double(m)).collect()),
        ));
        out.push(genotype);
    }
    Ok(out)
}

/// `GATKVariantContextUtils.trimAlleles(vc, true, true)`, with the trimmed-to-untrimmed allele
/// map the likelihoods are marginalized through.
#[allow(clippy::type_complexity)]
fn trim_alleles(
    vc: &VariantContext,
) -> Result<(VariantContext, Vec<(Allele, Vec<Allele>)>), String> {
    let identity: Vec<(Allele, Vec<Allele>)> = vc
        .alleles
        .iter()
        .map(|a| (a.clone(), vec![a.clone()]))
        .collect();
    let span_del = |a: &Allele| a.display_string() == "*";
    if vc.alleles.len() <= 1 || vc.alleles.iter().any(|a| a.len() == 1 && !span_del(a)) {
        return Ok((vc.clone(), identity));
    }
    let comparable: Vec<&Allele> = vc
        .alleles
        .iter()
        .filter(|a| !a.is_symbolic() && !span_del(a))
        .collect();
    let owned: Vec<Vec<u8>> = comparable
        .iter()
        .map(|a| a.base_string().into_bytes())
        .collect();
    let sequences: Vec<&[u8]> = owned.iter().map(Vec::as_slice).collect();
    let mut ranges: Vec<IndexRange> = comparable
        .iter()
        .map(|a| IndexRange::new(0, a.len() as i32))
        .collect();
    let (left, right) =
        normalize_alleles(&sequences, &mut ranges, 0, true).map_err(|e| format!("{e:?}"))?;
    let end_trim = right;
    let start_trim = -left;
    let empty = ranges.iter().any(|r| r.size() == 0);
    let end_clip = if empty && start_trim == 0 {
        end_trim - 1
    } else {
        end_trim
    };
    let start_clip = if empty && start_trim > 0 {
        start_trim - 1
    } else {
        start_trim
    };
    let forward_end = start_clip - 1;
    if forward_end == -1 && end_clip == 0 {
        return Ok((vc.clone(), identity));
    }
    let mut map = Vec::new();
    let mut trimmed = Vec::new();
    for a in &vc.alleles {
        let new = if a.is_symbolic() || span_del(a) {
            a.clone()
        } else {
            let bases = a.base_string().into_bytes();
            Allele::create(
                &bases[(forward_end + 1) as usize..bases.len() - end_clip as usize],
                a.is_reference(),
            )
            .map_err(|e| e.to_string())?
        };
        map.push((new.clone(), vec![a.clone()]));
        trimmed.push((a.clone(), new));
    }
    let mut out = vc.clone();
    out.start = vc.start + i64::from(forward_end + 1);
    out.alleles = trimmed.iter().map(|(_, n)| n.clone()).collect();
    out.stop = out.start + out.alleles[0].len() as i64 - 1;
    let genotypes: Vec<Genotype> = vc
        .genotypes
        .iter()
        .map(|g| {
            let mut g = g.clone();
            g.alleles = g
                .alleles
                .iter()
                .map(|a| {
                    trimmed
                        .iter()
                        .find(|(old, _)| old == a)
                        .map(|(_, new)| new.clone())
                        .unwrap_or_else(|| a.clone())
                })
                .collect();
            g
        })
        .collect();
    out.genotypes = GenotypesContext::new(genotypes);
    Ok((out, map))
}

/// `Mutect2Engine.MIN_PALINDROME_SIZE`.
const MIN_PALINDROME_SIZE: i32 = 5;
/// `PalindromeArtifactClipReadTransformer.MIN_FRACTION_OF_MATCHING_BASES`.
const MIN_FRACTION_OF_MATCHING_BASES: f64 = 0.9;

/// `PalindromeArtifactClipReadTransformer.apply`, Mutect2's post-filter transformer unless
/// `--ignore-itr-artifacts`: a properly paired read whose clip or insertion on the side facing its
/// mate is, with five more bases, the reverse complement of the reference just inside the adaptor
/// boundary (nine tenths of the bases matching) has that clip hard-clipped away. `contig` is the
/// read's contig's whole sequence, `None` when the dictionary does not hold it.
pub fn palindrome_artifact_clip(
    read: &BamRecord,
    header: &SamHeader,
    contig: Option<&[u8]>,
) -> Result<BamRecord, String> {
    let Some(boundary) = read_utils::adaptor_boundary(read) else {
        return Ok(read.clone());
    };
    if !read::is_proper_pair(read) {
        return Ok(read.clone());
    }
    let elements = &read.cigar.elements;
    let (Some(first), Some(last)) = (elements.first(), elements.last()) else {
        return Ok(read.clone());
    };
    let upstream = read.inferred_insert_size > 0;
    let clip_like = |op: htsjdk_bam::cigar::Op| {
        matches!(op, htsjdk_bam::cigar::Op::S | htsjdk_bam::cigar::Op::I)
    };
    if (upstream && !clip_like(first.op)) || (!upstream && !clip_like(last.op)) {
        return Ok(read.clone());
    }
    let artifact = if upstream { first.length } else { last.length } as i32;
    let length = read.read_bases.len() as i32;
    let compared = (artifact + MIN_PALINDROME_SIZE).min(length);
    let ref_start = if upstream {
        boundary - compared
    } else {
        boundary + 1
    };
    let ref_end = if upstream {
        boundary - 1
    } else {
        boundary + compared
    };
    let Some(bases) = contig else {
        return Ok(read.clone());
    };
    if ref_start < 1 || ref_end > bases.len() as i32 {
        return Ok(read.clone());
    }
    if (upstream && ref_start < read_utils::start(read))
        || (!upstream && read_utils::end(read) < ref_end)
    {
        return Ok(read.clone());
    }
    let mut matches = 0i32;
    let mut index = if upstream { compared - 1 } else { length - 1 };
    for position in ref_start..=ref_end {
        let complement = match bases[position as usize - 1] {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            other => other,
        };
        if usize::try_from(index)
            .ok()
            .and_then(|i| read.read_bases.get(i))
            == Some(&complement)
        {
            matches += 1;
        }
        index -= 1;
    }
    if f64::from(matches) / f64::from(compared) >= MIN_FRACTION_OF_MATCHING_BASES {
        let mut clipper = gatk_engine::clipping::ReadClipper::new(read, Some(header));
        clipper.add_op(if upstream {
            gatk_engine::clipping::ClippingOp {
                start: 0,
                stop: artifact - 1,
            }
        } else {
            gatk_engine::clipping::ClippingOp {
                start: length - artifact,
                stop: length,
            }
        });
        clipper
            .clip_read(gatk_engine::clipping::ClippingRepresentation::HardclipBases)
            .map_err(|e| format!("{e:?}"))
    } else {
        Ok(read.clone())
    }
}
