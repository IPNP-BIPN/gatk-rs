//! `HaplotypeCaller` in VCF mode: the `AssemblyRegionWalker` startup around
//! [`gatk_tools::haplotype_caller::call_variants`], and the header `makeVCFHeader` writes.
//!
//! The reads are the input's, through `IUPACReadTransformer(true)` and then the resolved read
//! filters (`makeStandardHCReadFilters` by default); the traversal intervals are `-L`'s, or every
//! contig of the reads' dictionary. Every argument the port does not honour is refused when it is
//! set, rather than silently left at its default: GVCF output, `--alleles`, the assembler's own
//! knobs, pileup detection and the flow-based paths among them. DRAGEN mode is applied as the
//! reference applies it, while parsing ([`haplotype_caller_modes`]).

use super::*;
use gatk_annotation::catalogue;
use gatk_engine::allele_frequency_calculator::Priors;
use gatk_engine::assembly_region_iterator::AssemblyRegionArgs;
use gatk_engine::assembly_region_trimmer::TrimmerArguments;
use gatk_engine::pair_hmm_likelihood_engine::{LikelihoodEngineArguments, PcrErrorModel};
use gatk_engine::read_threading_assembler::AssemblerSettings;
use gatk_tools::genotyping_engine::{Configuration, SubsetMethod};
use gatk_tools::haplotype_caller::{self as hc, HaplotypeCallerArguments};
use gatk_tools::haplotype_caller_engine::IsActiveArguments;
use gatk_tools::hc_genotyping::HcGenotypingArguments;
use gatk_tools::hc_region::{AssemblerArguments, RegionAssemblyArguments};
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};

/// The tool's arguments this port reads. Any other argument the tool declares, apart from the
/// read-filter plugin's, is refused when the command line sets it.
const HONOURED: &[&str] = &[
    // The linked de Bruijn graph and its artificial haplotype recovery.
    "linked-de-bruijn-graph",
    "disable-artificial-haplotype-recovery",
    // The annotation plugin's pedigree arguments.
    "pedigree",
    "founder-id",
    "denovo-parent-gq-threshold",
    "denovo-depth-threshold",
    // The engine and the walker.
    "input",
    "output",
    "reference",
    "intervals",
    "exclude-intervals",
    "interval-set-rule",
    "interval-padding",
    "interval-exclusion-padding",
    "interval-merging-rule",
    "read-validation-stringency",
    "read-index",
    "sequence-dictionary",
    "seconds-between-progress-updates",
    "disable-sequence-dictionary-validation",
    "create-output-bam-index",
    "create-output-bam-md5",
    "create-output-variant-index",
    "create-output-variant-md5",
    "max-variants-per-shard",
    "lenient",
    "add-output-sam-program-record",
    "add-output-vcf-command-line",
    "cloud-prefetch-buffer",
    "cloud-index-prefetch-buffer",
    "disable-bam-index-caching",
    "sites-only-vcf-output",
    "tmp-dir",
    "help",
    "version",
    "arguments_file",
    "showHidden",
    "verbosity",
    "QUIET",
    "use-jdk-deflater",
    "use-jdk-inflater",
    "gcs-max-retries",
    "gcs-project-for-requester-pays",
    "gatk-config-file",
    "read-filter",
    "inverted-read-filter",
    "disable-read-filter",
    "disable-tool-default-read-filters",
    "annotation",
    "annotations-to-exclude",
    "annotation-group",
    "disable-tool-default-annotations",
    "enable-all-annotations",
    // The assembly region walker.
    "min-assembly-region-size",
    "max-assembly-region-size",
    "active-probability-threshold",
    "max-prob-propagation-distance",
    "force-active",
    "assembly-region-padding",
    "padding-around-indels",
    "padding-around-snps",
    "padding-around-strs",
    "max-extension-into-assembly-region-padding-legacy",
    "max-reads-per-alignment-start",
    "enable-legacy-assembly-region-trimming",
    // Genotyping.
    "standard-min-confidence-threshold-for-calling",
    "max-alternate-alleles",
    "max-genotype-count",
    "sample-ploidy",
    "heterozygosity",
    "indel-heterozygosity",
    "heterozygosity-stdev",
    "annotate-with-num-discovered-alleles",
    "do-not-run-physical-phasing",
    "allele-informative-reads-overlap-margin",
    "disable-spanning-event-genotyping",
    "max-mnp-distance",
    "mapping-quality-threshold-for-genotyping",
    "reference-model-deletion-quality",
    // Assembly and likelihoods.
    "min-base-quality-score",
    "dont-use-soft-clipped-bases",
    "override-fragment-softclip-check",
    "soft-clip-low-quality-ends",
    "do-not-correct-overlapping-quality",
    "pcr-indel-model",
    "phred-scaled-global-read-mismapping-rate",
    "base-quality-score-threshold",
    "pair-hmm-gap-continuation-penalty",
    "expected-mismatch-rate-for-read-disqualification",
    "enable-dynamic-read-disqualification-for-genotyping",
    "dynamic-read-disqualification-threshold",
    "disable-symmetric-hmm-normalizing",
    "disable-cap-base-qualities-to-map-quality",
    // Reference confidence.
    "emit-ref-confidence",
    "gvcf-gq-bands",
    "floor-blocks",
    "indel-size-to-eliminate-in-ref-model",
    // The read-threading assembler.
    "kmer-size",
    "dont-increase-kmer-sizes-for-cycles",
    "allow-non-unique-kmers-in-ref",
    "num-pruning-samples",
    "min-pruning",
    "adaptive-pruning",
    "adaptive-pruning-initial-error-rate",
    "pruning-lod-threshold",
    "pruning-seeding-lod-threshold",
    "max-unpruned-variants",
    "enable-legacy-graph-cycle-detection",
    "num-matching-bases-in-dangling-end-to-recover",
    "do-not-recover-dangling-branches",
    "recover-all-dangling-branches",
    "min-dangling-branch-length",
    "max-num-haplotypes-in-population",
    // Only the linked de Bruijn graph reads it, and that graph is not ported.
    "disable-artificial-haplotype-recovery",
    // Genotyping options.
    "all-site-pls",
    "output-mode",
    "genotype-assignment-method",
    "use-posteriors-to-calculate-qual",
    // Without contamination the filtered read map is the genotyping one either way.
    "use-filtered-reads-for-annotations",
    // Without --alleles, a bamout or a native PairHMM these change nothing.
    "force-call-filtered-alleles",
    "bam-writer-type",
    "native-pair-hmm-use-double-precision",
    "dont-use-dragstr-priors",
    "dont-use-dragstr-pair-hmm-scores",
    "sample-name",
    "disable-optimizations",
    // DRAGEN mode and DRAGstr.
    "dragen-mode",
    "apply-bqd",
    "apply-frd",
    "transform-dragen-mapping-quality",
    "max-effective-depth-adjustment-for-frd",
    "dragstr-params-path",
    "dragstr-het-hom-ratio",
    // Refused below beside --dragen-mode, and otherwise NONE.
    "flow-mode",
    // Logged only.
    "debug-assembly",
    // The implementations, refused below when they name a kernel this port does not carry.
    "pair-hmm-implementation",
    "smith-waterman",
    "likelihood-calculation-engine",
    // RampedHaplotypeCaller's own: the off ramp is written, the on ramp read, and the two debug
    // switches only log.
    "off-ramp-type",
    "off-ramp-file",
    "on-ramp-type",
    "on-ramp-file",
    "ramps-debug-reads",
    "ramps-debug-post-assembler-on",
    // The three Smith-Waterman weightings.
    "smith-waterman-dangling-end-match-value",
    "smith-waterman-dangling-end-mismatch-penalty",
    "smith-waterman-dangling-end-gap-open-penalty",
    "smith-waterman-dangling-end-gap-extend-penalty",
    "smith-waterman-haplotype-to-reference-match-value",
    "smith-waterman-haplotype-to-reference-mismatch-penalty",
    "smith-waterman-haplotype-to-reference-gap-open-penalty",
    "smith-waterman-haplotype-to-reference-gap-extend-penalty",
    "smith-waterman-read-to-haplotype-match-value",
    "smith-waterman-read-to-haplotype-mismatch-penalty",
    "smith-waterman-read-to-haplotype-gap-open-penalty",
    "smith-waterman-read-to-haplotype-gap-extend-penalty",
];

/// `getDragenVersion3412NameValuePairs`: what `--dragen-mode` sets, unless the command line set it.
const DRAGEN_3412_MODE: &[(&str, &str)] = &[
    ("apply-bqd", "true"),
    ("apply-frd", "true"),
    ("transform-dragen-mapping-quality", "true"),
    ("soft-clip-low-quality-ends", "true"),
    ("mapping-quality-threshold-for-genotyping", "1"),
    ("minimum-mapping-quality", "1"),
    ("allele-informative-reads-overlap-margin", "1"),
    ("disable-cap-base-qualities-to-map-quality", "true"),
    (
        "enable-dynamic-read-disqualification-for-genotyping",
        "true",
    ),
    ("expected-mismatch-rate-for-read-disqualification", "0.03"),
    ("genotype-assignment-method", "USE_POSTERIOR_PROBABILITIES"),
    ("padding-around-indels", "150"),
    ("standard-min-confidence-threshold-for-calling", "3.0"),
    ("use-posteriors-to-calculate-qual", "true"),
];

/// `HaplotypeCaller.customCommandLineValidation`, which runs as the last step of parsing: DRAGEN
/// mode beside flow mode is refused; otherwise DRAGEN mode lowers the resolved mapping-quality
/// filter to 1 **whatever the command line said**, and then sets its arguments where the command
/// line did not. `--dragen-378-concordance-mode` is refused before this, as an argument the port
/// does not carry.
pub fn haplotype_caller_modes(parser: &mut Parser, tool: &str) -> Result<(), Thrown> {
    let dragen = flag(parser, "dragen-mode");
    let flow = scalar(parser, "flow-mode").is_some_and(|mode| mode != "NONE");
    if dragen && flow {
        return Err(Thrown::user(
            "dragen mode and flow mode can't be both specified",
        ));
    }
    if dragen {
        let filter_resolved = resolve_read_filters(parser, tool)?
            .iter()
            .any(|filter| filter.name == "MappingQualityReadFilter");
        if filter_resolved {
            parser
                .override_value("minimum-mapping-quality", "1")
                .map_err(|error| Thrown::command_line(error.message))?;
        }
        parser
            .set_mode_values(DRAGEN_3412_MODE)
            .map_err(|error| Thrown::non_user(error.class, error.message))?;
    }
    Ok(())
}

/// `buildRamps`' off ramp: its type, and the file it is written to, which a type other than `NONE`
/// requires.
fn off_ramp(
    parser: &Parser,
) -> Result<Option<(gatk_tools::ramp_files::OffRampType, std::path::PathBuf)>, Thrown> {
    let Some(name) = scalar(parser, "off-ramp-type") else {
        return Ok(None);
    };
    let kind = gatk_tools::ramp_files::OffRampType::parse(&name)
        .ok_or_else(|| Thrown::command_line(format!("unknown off ramp type {name}")))?;
    let Some(kind) = kind else {
        return Ok(None);
    };
    let file = scalar(parser, "off-ramp-file").ok_or_else(|| {
        Thrown::non_user(
            "java.lang.RuntimeException",
            "rampFile must be specified".to_string(),
        )
    })?;
    Ok(Some((kind, std::path::PathBuf::from(file))))
}

/// `buildRamps`' on ramp, after the off ramp: its type, the file it is read from, and the rule
/// that an on ramp comes before any off ramp.
fn on_ramp(
    parser: &Parser,
    off: Option<gatk_tools::ramp_files::OffRampType>,
) -> Result<Option<(gatk_tools::ramp_files::OnRampType, std::path::PathBuf)>, Thrown> {
    use gatk_tools::ramp_files::{OffRampType, OnRampType};
    let Some(name) = scalar(parser, "on-ramp-type") else {
        return Ok(None);
    };
    let kind = OnRampType::parse(&name)
        .ok_or_else(|| Thrown::command_line(format!("unknown on ramp type {name}")))?;
    let Some(kind) = kind else {
        return Ok(None);
    };
    let file = scalar(parser, "on-ramp-file").ok_or_else(|| {
        Thrown::non_user(
            "java.lang.RuntimeException",
            "rampFile must be specified".to_string(),
        )
    })?;
    let out_of_order = match kind {
        OnRampType::PostFilter => off == Some(OffRampType::PreFilter),
        OnRampType::PostAssembler => off.is_some(),
    };
    if out_of_order {
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            "an on ramp, if present, must be before the off ramp".to_string(),
        ));
    }
    Ok(Some((kind, std::path::PathBuf::from(file))))
}

fn limitation(what: &str) -> Thrown {
    Thrown::non_user(
        PORT_LIMITATION,
        format!("{what} This message is the port's own and not GATK's."),
    )
}

/// The refusal of any declared argument the port does not read.
fn refuse_unported(parser: &Parser, tool: &str) -> Result<(), Thrown> {
    let declarations = gatk_tools::tool_declarations::declarations(tool).unwrap_or(&[]);
    for definition in parser.definitions() {
        let name = definition.long_name();
        if !definition.has_been_set() || HONOURED.contains(&name) {
            continue;
        }
        let read_filter_argument = declarations.iter().any(|declaration| {
            declaration.long_name == name
                && declaration.controlled_by == Some("GATKReadFilterPluginDescriptor")
        });
        if !read_filter_argument {
            return Err(limitation(&format!(
                "--{name} is a {tool} argument this port does not carry yet."
            )));
        }
    }
    Ok(())
}

fn pcr_model(parser: &Parser) -> PcrErrorModel {
    match scalar(parser, "pcr-indel-model").as_deref() {
        Some("NONE") => PcrErrorModel::None,
        Some("HOSTILE") => PcrErrorModel::Hostile,
        Some("AGGRESSIVE") => PcrErrorModel::Aggressive,
        _ => PcrErrorModel::Conservative,
    }
}

/// `HaplotypeCaller.doWork`.
pub fn haplotype_caller(parser: &Parser) -> Outcome {
    run(parser, "HaplotypeCaller")
}

/// `RampedHaplotypeCaller.doWork`: HaplotypeCaller with `RampedHaplotypeCallerEngine` and
/// non-random downsampling. With no ramp the engine still differs from HaplotypeCaller's where
/// `callRegion` measures read stubs, filters reads and realigns them.
pub fn ramped_haplotype_caller(parser: &Parser) -> Outcome {
    run(parser, "RampedHaplotypeCaller")
}

fn run(parser: &Parser, tool: &str) -> Outcome {
    let ramped = tool == "RampedHaplotypeCaller";
    refuse_unported(parser, tool)?;
    // The kernels this port carries: the Java LOGLESS_CACHING PairHMM behind the PairHMM engine.
    // FASTEST_AVAILABLE resolves to it wherever GKL's native kernel does not load, which is the
    // reference container; a kernel named outright is refused.
    match scalar(parser, "likelihood-calculation-engine").as_deref() {
        None | Some("PairHMM") => {}
        Some(other) => {
            return Err(limitation(&format!(
                "the {other} likelihood engine is not ported for HaplotypeCaller."
            )))
        }
    }
    match scalar(parser, "pair-hmm-implementation").as_deref() {
        None | Some("LOGLESS_CACHING") | Some("FASTEST_AVAILABLE") => {}
        Some(other) => {
            return Err(limitation(&format!(
                "the {other} PairHMM is not ported; only LOGLESS_CACHING is."
            )))
        }
    }
    let output_mode =
        scalar(parser, "output-mode").unwrap_or_else(|| "EMIT_VARIANTS_ONLY".to_string());
    let assignment_method = SubsetMethod::from_name(
        &scalar(parser, "genotype-assignment-method")
            .unwrap_or_else(|| "USE_PLS_TO_ASSIGN".to_string()),
    )
    .unwrap_or(SubsetMethod::UsePlsToAssign);
    let resolved = catalogue::resolve(
        &catalogue::AnnotationArguments {
            annotations: arguments(parser, "annotation"),
            groups: arguments(parser, "annotation-group"),
            excluded: arguments(parser, "annotations-to-exclude"),
            disable_tool_defaults: flag(parser, "disable-tool-default-annotations"),
            enable_all: flag(parser, "enable-all-annotations"),
        },
        &["StandardAnnotation", "StandardHCAnnotation"],
        &[],
    )
    .map_err(|error| Thrown::command_line(error.message()))?;
    pedigree_refusal(parser, &resolved)?;
    let emit_mode = scalar(parser, "emit-ref-confidence").unwrap_or_else(|| "NONE".to_string());
    let emit_reference_confidence = emit_mode != "NONE";
    // `filterReferenceConfidenceAnnotations`: StrandBiasBySample added, the four site
    // annotations a GVCF cannot carry removed.
    let mut resolved = resolved;
    if emit_reference_confidence {
        if !resolved
            .iter()
            .any(|entry| entry.name == "StrandBiasBySample")
        {
            resolved.push(catalogue::entry("StrandBiasBySample").expect("a known annotation"));
            resolved.sort_by(|a, b| a.name.cmp(b.name));
        }
        resolved.retain(|entry| {
            !matches!(
                entry.name,
                "ChromosomeCounts" | "FisherStrand" | "StrandOddsRatio" | "QualByDepth"
            )
        });
    }
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, tool)?;
    // `initializeSamples`.
    let header_samples = hc::samples_from_header(&header);
    let sample_name = match scalar(parser, "sample-name") {
        Some(name) if !header_samples.contains(&name) => {
            return Err(Thrown::command_line(
                "Argument --sample_name has a bad value: Specified name does not exist in input bam files",
            ))
        }
        // The only sample named: nothing to filter.
        Some(_) if header_samples.len() == 1 => None,
        other => other,
    };
    let sample_count = if sample_name.is_some() {
        1
    } else {
        header_samples.len()
    };
    if emit_reference_confidence && sample_count != 1 {
        return Err(Thrown::command_line(
            "Argument emit-ref-confidence has a bad value: Can only be used in single sample \
             mode currently. Use the --sample-name argument to run on a single sample out of a \
             multi-sample BAM file.",
        ));
    }
    if flag(parser, "floor-blocks") && !emit_reference_confidence {
        return Err(Thrown::user(
            "floor-blocks refers to GVCF blocks, so reference confidence mode \
             (emit-ref-confidence) must be specified.",
        ));
    }
    // `isAlleleSpecificExceptHmerLengthOrStrandBiasMode`, checked by the engine's constructor once
    // the arguments are validated: an `AS_` INFO or genotype annotation other than
    // `AS_StrandBiasMutectAnnotation`, in VCF mode.
    let allele_specific = resolved.iter().any(|entry| {
        matches!(
            entry.kind,
            gatk_annotation::catalogue::Kind::Info | gatk_annotation::catalogue::Kind::Genotype
        ) && entry.name.starts_with("AS_")
            && entry.name != "AS_StrandBiasMutectAnnotation"
    });
    if allele_specific && !emit_reference_confidence {
        return Err(Thrown::user(
            "Allele-specific annotations are not yet supported in the VCF mode",
        ));
    }
    // `validateAndInitializeArgs`: GVCF mode calls at a threshold of -0.0, which the
    // active-region engine then takes as its own minimum too.
    let calling_threshold = if emit_reference_confidence {
        -0.0
    } else {
        double_or(
            parser,
            "standard-min-confidence-threshold-for-calling",
            30.0,
        )
    };
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;

    // `getAllIntervalsForReference` over the reads' dictionary when no `-L` was given.
    let traversal: Vec<gatk_engine::interval::SimpleInterval> = if intervals.is_empty() {
        header
            .sequences
            .iter()
            .filter_map(|sequence| {
                gatk_engine::interval::SimpleInterval::new(&sequence.name, 1, sequence.length)
            })
            .collect()
    } else {
        intervals
    };

    // `makePreReadFilterTransformer`, then the filters.
    let filter = read_filter(parser, &filters, &header)?;
    let transform_mapping_quality = flag(parser, "transform-dragen-mapping-quality");
    let mut reads = Vec::new();
    for mut read in source
        .iter_all()
        .map_err(|error| Thrown::user(format!("{error:?}")))?
    {
        hc::iupac_to_n(&mut read).map_err(Thrown::user)?;
        if filter(&read) {
            // `makePostReadFilterTransformer`: DRAGEN's mapping quality swap, after the filters.
            if transform_mapping_quality {
                gatk_engine::dragen_genotypes::transform_mapping_quality(&mut read).map_err(
                    |message| Thrown::non_user("java.lang.IllegalArgumentException", message),
                )?;
            }
            reads.push(read);
        }
    }

    // `DragstrParamUtils.parse(hcArgs.likelihoodArgs.dragstrParams)`.
    let dragstr_params = match argument(parser, "dragstr-params-path") {
        None => None,
        Some(path) => {
            let text = std::fs::read_to_string(&path)
                .map_err(|_| Thrown::user(format!("Couldn't read file {path}")))?;
            Some(
                gatk_engine::dragstr::parse(&text, &path)
                    .map_err(|error| Thrown::user(format!("{error:?}")))?,
            )
        }
    };
    let ploidy = number_or(parser, "sample-ploidy", 2).max(0) as usize;
    let calling_confidence = calling_threshold;
    let max_alternate_alleles = number_or(parser, "max-alternate-alleles", 6).max(0) as usize;
    let priors = Priors {
        snp_heterozygosity: double_or(parser, "heterozygosity", 1e-3),
        indel_heterozygosity: double_or(parser, "indel-heterozygosity", 1.0 / 8000.0),
        heterozygosity_standard_deviation: double_or(parser, "heterozygosity-stdev", 0.01),
        sample_ploidy: ploidy,
    };
    let min_base_quality = number_or(parser, "min-base-quality-score", 10) as u8;
    let region = AssemblyRegionArgs {
        min_assembly_region_size: number_or(parser, "min-assembly-region-size", 50),
        max_assembly_region_size: number_or(parser, "max-assembly-region-size", 300),
        active_prob_threshold: double_or(parser, "active-probability-threshold", 0.002),
        max_prob_propagation_distance: number_or(parser, "max-prob-propagation-distance", 50),
        force_active: flag(parser, "force-active"),
        assembly_region_padding: number_or(parser, "assembly-region-padding", 100),
        max_reads_per_alignment_start: number_or(parser, "max-reads-per-alignment-start", 50),
        indel_padding_for_genotyping: number_or(parser, "padding-around-indels", 75),
        snp_padding_for_genotyping: number_or(parser, "padding-around-snps", 20),
        str_padding_for_genotyping: number_or(parser, "padding-around-strs", 75),
        max_extension_into_region_padding: number_or(
            parser,
            "max-extension-into-assembly-region-padding-legacy",
            25,
        ),
    };
    region
        .validate()
        .map_err(|error| Thrown::command_line(error.message()))?;
    let mismapping = double_or(parser, "phred-scaled-global-read-mismapping-rate", 45.0);
    let sw = smith_waterman_arguments(parser);
    // `buildRamps`: the off ramp first, then the on ramp.
    let off = if ramped { off_ramp(parser)? } else { None };
    artificial_haplotype_recovery_refusal(parser)?;
    let arguments = HaplotypeCallerArguments {
        trimmer: TrimmerArguments {
            assembly_region_padding: region.assembly_region_padding,
            indel_padding_for_genotyping: region.indel_padding_for_genotyping,
            snp_padding_for_genotyping: region.snp_padding_for_genotyping,
            str_padding_for_genotyping: region.str_padding_for_genotyping,
            max_extension_into_region_padding: region.max_extension_into_region_padding,
            enable_legacy_assembly_region_trimming: flag(
                parser,
                "enable-legacy-assembly-region-trimming",
            ),
        },
        region,
        assembler: assembler_arguments(parser, flag(parser, "adaptive-pruning")),
        is_active: IsActiveArguments {
            sample_ploidy: ploidy,
            snp_heterozygosity: priors.snp_heterozygosity,
            indel_heterozygosity: priors.indel_heterozygosity,
            heterozygosity_standard_deviation: priors.heterozygosity_standard_deviation,
            standard_confidence_for_calling: calling_confidence,
            max_alternate_alleles,
            min_base_quality_score: min_base_quality,
            ref_model_deletion_quality: number_or(parser, "reference-model-deletion-quality", 30)
                as u8,
        },
        assembly: RegionAssemblyArguments {
            min_base_quality_score: min_base_quality,
            dont_use_soft_clipped_bases: flag(parser, "dont-use-soft-clipped-bases"),
            soft_clip_low_quality_ends: flag(parser, "soft-clip-low-quality-ends"),
            override_softclip_fragment_check: flag(parser, "override-fragment-softclip-check"),
            do_not_correct_overlapping_base_qualities: flag(
                parser,
                "do-not-correct-overlapping-quality",
            ),
            dangling_end_sw: sw.0,
            haplotype_to_reference_sw: sw.1,
            read_to_haplotype_sw: sw.2,
        },
        likelihoods: LikelihoodEngineArguments {
            gap_continuation_penalty: number_or(parser, "pair-hmm-gap-continuation-penalty", 10)
                as i8,
            log10_global_read_mismapping_rate: if mismapping < 0.0 {
                -f64::MAX
            } else {
                mismapping / -10.0
            },
            pcr_error_model: pcr_model(parser),
            base_quality_score_threshold: number_or(parser, "base-quality-score-threshold", 18)
                as i8,
            dynamic_disqualification: flag(
                parser,
                "enable-dynamic-read-disqualification-for-genotyping",
            ),
            read_disqualification_scale: double_or(
                parser,
                "dynamic-read-disqualification-threshold",
                1.0,
            ),
            expected_error_rate_per_base: double_or(
                parser,
                "expected-mismatch-rate-for-read-disqualification",
                0.02,
            ),
            symmetrically_normalize_alleles_to_reference: !flag(
                parser,
                "disable-symmetric-hmm-normalizing",
            ),
            disable_cap_read_qualities_to_map_q: flag(
                parser,
                "disable-cap-base-qualities-to-map-quality",
            ),
            modify_soft_clipped_bases: !flag(parser, "soft-clip-low-quality-ends"),
            dragstr_params: if flag(parser, "dont-use-dragstr-pair-hmm-scores") {
                None
            } else {
                dragstr_params.clone()
            },
        },
        genotyping: HcGenotypingArguments {
            sample_ploidy: ploidy,
            informative_read_overlap_margin: number_or(
                parser,
                "allele-informative-reads-overlap-margin",
                2,
            ),
            disable_spanning_event_genotyping: flag(parser, "disable-spanning-event-genotyping"),
            max_genotype_count: number_or(parser, "max-genotype-count", 1024).max(0) as usize,
            max_mnp_distance: number_or(parser, "max-mnp-distance", 0),
            // Phasing is only for diploid samples, and only with reference confidence.
            do_physical_phasing: emit_reference_confidence
                && ploidy == 2
                && !flag(parser, "do-not-run-physical-phasing"),
            emit_reference_confidence,
            ref_confidence: gatk_tools::reference_confidence_model::RefConfidenceArguments {
                indel_informative_depth_indel_size: number_or(
                    parser,
                    "indel-size-to-eliminate-in-ref-model",
                    10,
                )
                .max(0) as usize,
                ref_model_deletion_quality: number_or(
                    parser,
                    "reference-model-deletion-quality",
                    30,
                ) as u8,
                use_soft_clipped_bases: !flag(parser, "override-fragment-softclip-check"),
            },
            ramped_engine: ramped,
            dragen_models: (flag(parser, "apply-bqd") || flag(parser, "apply-frd")).then(|| {
                gatk_tools::hc_genotyping::DragenModelArguments {
                    bqd: flag(parser, "apply-bqd"),
                    frd: flag(parser, "apply-frd"),
                    max_effective_depth_adjustment: number_or(
                        parser,
                        "max-effective-depth-adjustment-for-frd",
                        0,
                    ),
                }
            }),
            dragstr: dragstr_params.clone().map(|params| {
                gatk_tools::hc_genotyping::DragstrGenotyping {
                    params,
                    use_priors: !flag(parser, "dont-use-dragstr-priors"),
                    het_hom_ratio: number_or(parser, "dragstr-het-hom-ratio", 2),
                    snp_heterozygosity: priors.snp_heterozygosity,
                    indel_heterozygosity: priors.indel_heterozygosity,
                }
            }),
        },
        calling: Configuration {
            standard_confidence_for_calling: calling_confidence,
            max_alternate_alleles,
            sample_ploidy: ploidy,
            annotate_number_of_alleles_discovered: flag(
                parser,
                "annotate-with-num-discovered-alleles",
            ),
            emit_all_active_sites: output_mode == "EMIT_ALL_ACTIVE_SITES",
            allele_specific: false,
            emit_all_confident_sites: output_mode == "EMIT_ALL_CONFIDENT_SITES",
            annotate_all_sites_with_pls: emit_reference_confidence || flag(parser, "all-site-pls"),
            // `forceKeepAllele` is true for every allele with reference confidence.
            force_keep_all_alleles: emit_reference_confidence,
            assignment_method,
        },
        priors,
        sample_name: sample_name.clone(),
        disable_optimizations: flag(parser, "disable-optimizations"),
        ramps_debug_post_assembler_on: flag(parser, "ramps-debug-post-assembler-on"),
        mapping_quality_threshold: number_or(parser, "mapping-quality-threshold-for-genotyping", 20)
            as u8,
        use_posteriors_to_calculate_qual: flag(parser, "use-posteriors-to-calculate-qual"),
        non_random_downsampling: ramped,
        off_ramp: off.clone(),
        on_ramp: if ramped {
            on_ramp(parser, off.as_ref().map(|(kind, _)| *kind))?
        } else {
            None
        },
        compression: crate::runners::output_compression(parser),
    };
    let annotator = gatk_tools::variant_annotator_engine::Engine {
        resolved: resolved.clone(),
        overlap_names: Vec::new(),
        expressions: Vec::new(),
        allele_concordance: false,
        // The annotations' own arguments, which the plugin descriptor hands to their instances.
        pseudo_depth: std::cell::RefCell::new(
            gatk_annotation::allele_pseudo_depth::AllelePseudoDepth::new(
                double_or(parser, "dirichlet-prior-pseudo-count", 1.0),
                flag(parser, "dirichlet-keep-prior-in-count"),
                double_or(parser, "pseudo-count-weight-decay-rate", 1.0),
            ),
        ),
        assembly_complexity_reference_mode: flag(parser, "assembly-complexity-reference-mode"),
        pedigree: pedigree_arguments(parser)?,
    };

    let mut contig_bases = |contig: &str| -> Result<Vec<u8>, String> {
        let length = reference
            .sequences()
            .iter()
            .find(|(name, _)| name == contig)
            .map(|(_, length)| *length as i32)
            .ok_or_else(|| format!("contig {contig} is not in the reference"))?;
        reference
            .query(contig, 1, length)
            .map_err(|error| format!("{error:?}"))
    };
    let mut random = gatk_random();
    let called = hc::call_variants(
        &reads,
        &header,
        &traversal,
        &arguments,
        &annotator,
        &mut contig_bases,
        &mut random,
    );
    drop(random);
    // A failure mid-traversal still writes what came before it, and the writer's own refusal of
    // those records, raised as it closes, is what the run reports.
    let (mut written, failure) = match called {
        Ok(calls) => (calls, None),
        Err(error) => (
            error.calls_before,
            Some(Thrown::non_user(
                java_class_name(&error.class),
                error.message,
            )),
        ),
    };

    // `makeVCFHeader`.
    let compound =
        |key: &str, id: &str, number: Cardinality, line_type: LineType, description: &str| {
            HeaderLine::Compound {
                key: key.to_string(),
                id: id.to_string(),
                number,
                line_type,
                description: description.to_string(),
                extra: Vec::new(),
            }
        };
    let mut lines: Vec<HeaderLine> = default_tool_vcf_header_lines(parser, tool);
    if arguments.calling.annotate_number_of_alleles_discovered {
        lines.push(compound(
            "INFO",
            "NDA",
            Cardinality::Fixed(1),
            LineType::Integer,
            "Number of alternate alleles discovered (but not necessarily genotyped) at this site",
        ));
    }
    lines.extend(catalogue::descriptions(
        &resolved,
        emit_reference_confidence,
        false,
    ));
    lines.push(compound("INFO", "MLEAC", Cardinality::A, LineType::Integer, "Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed"));
    lines.push(compound("INFO", "MLEAF", Cardinality::A, LineType::Float, "Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed"));
    for id in ["GT", "GQ", "DP", "PL"] {
        lines.extend(htsjdk_vcf::standard_header_lines::standard_format_line(id));
    }
    if arguments.genotyping.do_physical_phasing {
        lines.push(compound("FORMAT", "PID", Cardinality::Fixed(1), LineType::String, "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group"));
        lines.push(compound("FORMAT", "PGT", Cardinality::Fixed(1), LineType::String, "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles"));
        lines.extend(htsjdk_vcf::standard_header_lines::standard_format_line(
            "PS",
        ));
    }
    lines.push(HeaderLine::Filter {
        id: "LowQual".to_string(),
        description: "Low quality".to_string(),
    });
    if assignment_method == SubsetMethod::UsePosteriorProbabilities {
        lines.push(compound(
            "FORMAT",
            "GP",
            Cardinality::G,
            LineType::Float,
            "genotype posterior in Phred Scale",
        ));
        lines.push(compound(
            "FORMAT",
            "PG",
            Cardinality::G,
            LineType::Float,
            "genotype priors in Phred Scale",
        ));
    }
    if argument(parser, "dragstr-params-path").is_some() {
        lines.push(compound(
            "INFO",
            "DRAGstrInfo",
            Cardinality::Fixed(2),
            LineType::Integer,
            "Indicates the period and repeat count",
        ));
        lines.push(compound(
            "INFO",
            "DRAGstrParams",
            Cardinality::Fixed(3),
            LineType::Float,
            "Parameters used (GOP, GCP, API)",
        ));
    }
    if emit_reference_confidence {
        lines.push(HeaderLine::Structured {
            key: "ALT".to_string(),
            fields: vec![
                ("ID".to_string(), "NON_REF".to_string()),
                (
                    "Description".to_string(),
                    "Represents any possible alternative allele not already represented at this \
                     location by REF and ALT"
                        .to_string(),
                ),
            ],
        });
    }
    for (index, sequence) in header.sequences.iter().enumerate() {
        lines.push(HeaderLine::Contig {
            index: index as i32,
            fields: vec![
                ("ID".to_string(), sequence.name.clone()),
                ("length".to_string(), sequence.length.to_string()),
            ],
        });
    }
    let mut unique: Vec<HeaderLine> = Vec::with_capacity(lines.len());
    for line in lines {
        if !unique.contains(&line) {
            unique.push(line);
        }
    }
    let mut vcf_header = VcfHeader {
        lines: unique,
        samples: match &sample_name {
            Some(name) => vec![name.clone()],
            None => hc::samples_from_header(&header),
        },
    };
    // `GVCFWriter`: the reference records combined into GQ-band blocks.
    if emit_mode == "GVCF" {
        let bands: Vec<i32> = arguments_or_default_bands(parser);
        let partitions = gatk_tools::gvcf_blocks::parse_partitions(&bands).map_err(|message| {
            Thrown::command_line(format!(
                "Argument GQBands has a bad value: are malformed: {message}"
            ))
        })?;
        gatk_tools::gvcf_blocks::add_ranges_to_header(&mut vcf_header, &partitions);
        let mut combiner = gatk_tools::gvcf_blocks::BlockCombiner::new(
            partitions,
            flag(parser, "floor-blocks"),
            None,
        );
        let mut combined = Vec::new();
        for vc in written {
            combiner
                .submit(vc)
                .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{error:?}")))?;
            combined.extend(combiner.drain());
        }
        combiner.end_of_input();
        combined.extend(combiner.drain());
        written = combined;
    }
    apply_sites_only(parser, &mut vcf_header, &mut written);
    // `HaplotypeCallerEngine.makeVCFWriter`, which never allows a key the header lacks.
    let text = write_vcf_strict_or_lenient(false, &vcf_header, &written)?;
    if let Some(failure) = failure {
        return Err(failure);
    }
    write_variant_output(parser, &output, &text)?;
    Ok(None)
}

/// `--gvcf-gq-bands`, which defaults to every GQ from 1 to 60 and then 70, 80, 90 and 99.
fn arguments_or_default_bands(parser: &Parser) -> Vec<i32> {
    let given: Vec<i32> = arguments(parser, "gvcf-gq-bands")
        .iter()
        .filter_map(|band| band.trim().parse().ok())
        .collect();
    if !given.is_empty() {
        return given;
    }
    (1..=60).chain([70, 80, 90, 99]).collect()
}

/// The three `--smith-waterman-*` groups, each defaulting to the constant the reference uses:
/// `(dangling end, haplotype to reference, read to haplotype)`.
pub(super) fn smith_waterman_arguments(
    parser: &Parser,
) -> (
    gatk_engine::smith_waterman::SwParameters,
    gatk_engine::smith_waterman::SwParameters,
    gatk_engine::smith_waterman::SwParameters,
) {
    use gatk_engine::smith_waterman::{
        SwParameters, ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, NEW_SW_PARAMETERS, STANDARD_NGS,
    };
    let group = |prefix: &str, default: SwParameters| SwParameters {
        match_value: number_or(
            parser,
            &format!("smith-waterman-{prefix}-match-value"),
            default.match_value,
        ),
        mismatch_penalty: number_or(
            parser,
            &format!("smith-waterman-{prefix}-mismatch-penalty"),
            default.mismatch_penalty,
        ),
        gap_open_penalty: number_or(
            parser,
            &format!("smith-waterman-{prefix}-gap-open-penalty"),
            default.gap_open_penalty,
        ),
        gap_extend_penalty: number_or(
            parser,
            &format!("smith-waterman-{prefix}-gap-extend-penalty"),
            default.gap_extend_penalty,
        ),
    };
    (
        group("dangling-end", STANDARD_NGS),
        group("haplotype-to-reference", NEW_SW_PARAMETERS),
        group(
            "read-to-haplotype",
            ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
        ),
    )
}

/// The read-threading assembler's arguments, which both callers declare; `adaptive` is
/// HaplotypeCaller's `--adaptive-pruning` or Mutect2's `!--disable-adaptive-pruning`.
pub(super) fn assembler_arguments(parser: &Parser, adaptive: bool) -> AssemblerArguments {
    AssemblerArguments {
        settings: AssemblerSettings {
            max_allowed_paths: number_or(parser, "max-num-haplotypes-in-population", 128).max(0)
                as usize,
            kmer_sizes: {
                let sizes: Vec<usize> = arguments(parser, "kmer-size")
                    .iter()
                    .filter_map(|k| k.trim().parse().ok())
                    .collect();
                if sizes.is_empty() {
                    vec![10, 25]
                } else {
                    sizes
                }
            },
            dont_increase_kmer_sizes_for_cycles: flag(
                parser,
                "dont-increase-kmer-sizes-for-cycles",
            ),
            allow_non_unique_kmers_in_ref: flag(parser, "allow-non-unique-kmers-in-ref"),
            num_pruning_samples: number_or(parser, "num-pruning-samples", 1).max(0) as usize,
            prune_factor: number_or(parser, "min-pruning", 2),
            use_adaptive_pruning: adaptive,
            initial_error_rate_for_pruning: double_or(
                parser,
                "adaptive-pruning-initial-error-rate",
                0.001,
            ),
            pruning_log_odds_threshold: double_or(
                parser,
                "pruning-lod-threshold",
                std::f64::consts::LN_10,
            ),
            pruning_seeding_log_odds_threshold: double_or(
                parser,
                "pruning-seeding-lod-threshold",
                4.0 * std::f64::consts::LN_10,
            ),
            max_unpruned_variants: number_or(parser, "max-unpruned-variants", 100),
            enable_legacy_graph_cycle_detection: flag(
                parser,
                "enable-legacy-graph-cycle-detection",
            ),
            min_matching_bases_to_dangling_end_recovery: number_or(
                parser,
                "num-matching-bases-in-dangling-end-to-recover",
                -1,
            ),
            use_linked_debruijn_graphs: flag(parser, "linked-de-bruijn-graph"),
        },
        recover_dangling_branches: !flag(parser, "do-not-recover-dangling-branches"),
        recover_all_dangling_branches: flag(parser, "recover-all-dangling-branches"),
        min_dangling_branch_length: number_or(parser, "min-dangling-branch-length", 4),
        disable_artificial_haplotype_recovery: flag(
            parser,
            "disable-artificial-haplotype-recovery",
        ),
    }
}

/// `ReadThreadingAssembler.setArtificialHaplotypeRecoveryMode`, which the engine's constructor
/// reaches: its test is the wrong way round, so disabling the recovery is refused WITH the linked
/// graph, and the refusal names an argument that does not exist.
pub(super) fn artificial_haplotype_recovery_refusal(parser: &Parser) -> Result<(), Thrown> {
    if flag(parser, "disable-artificial-haplotype-recovery")
        && flag(parser, "linked-de-bruijn-graph")
    {
        return Err(Thrown::user(
            "Argument '--experimental-dangling-branch-recover' requires '--linked-de-bruijn-graph' \
             to be set",
        ));
    }
    Ok(())
}
