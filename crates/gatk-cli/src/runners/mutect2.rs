//! `Mutect2`: the `AssemblyRegionWalker` startup around [`gatk_tools::mutect2::call_variants`],
//! the VCF it writes and the `.stats` file beside it.
//!
//! The reads are the input's through `makeStandardMutect2ReadFilters` and then the palindrome
//! artifact clipper; the samples named by `--normal-sample` are the normals and every other one a
//! tumor. Every argument the port does not honour is refused when it is set.

use super::*;
use gatk_annotation::catalogue;
use gatk_engine::assembly_region_iterator::AssemblyRegionArgs;
use gatk_engine::assembly_region_trimmer::TrimmerArguments;
use gatk_engine::pair_hmm_likelihood_engine::{LikelihoodEngineArguments, PcrErrorModel};
use gatk_tools::hc_region::RegionAssemblyArguments;
use gatk_tools::mutect2::{self as m2, Mutect2Arguments};
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};

/// `M2ArgumentCollection.DEFAULT_AF_FOR_TUMOR_ONLY_CALLING` and the tumor-normal one.
const DEFAULT_AF_FOR_TUMOR_ONLY_CALLING: f64 = 5e-8;
const DEFAULT_AF_FOR_TUMOR_NORMAL_CALLING: f64 = 1e-6;

/// The tool's arguments this port reads, beside the engine's and the read filters'.
const HONOURED: &[&str] = &[
    // The flow annotations' flow order.
    "flow-order-for-annotations",
    // The flow-based presets.
    "flow-mode",
    // Long hmer collapsing and allele filtering before genotyping.
    "flow-assembly-collapse-hmer-size",
    "flow-assembly-collapse-partial-mode",
    "flow-filter-alleles",
    "flow-filter-alleles-qual-threshold",
    "flow-filter-alleles-sor-threshold",
    "flow-filter-lone-alleles",
    // The linked de Bruijn graph and its artificial haplotype recovery.
    "linked-de-bruijn-graph",
    "disable-artificial-haplotype-recovery",
    // The annotation plugin's pedigree arguments.
    "pedigree",
    "founder-id",
    "denovo-parent-gq-threshold",
    "denovo-depth-threshold",
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
    "min-base-quality-score",
    "dont-use-soft-clipped-bases",
    "override-fragment-softclip-check",
    "soft-clip-low-quality-ends",
    "pcr-indel-model",
    "phred-scaled-global-read-mismapping-rate",
    "base-quality-score-threshold",
    "pair-hmm-gap-continuation-penalty",
    "expected-mismatch-rate-for-read-disqualification",
    "enable-dynamic-read-disqualification-for-genotyping",
    "dynamic-read-disqualification-threshold",
    "disable-symmetric-hmm-normalizing",
    "disable-cap-base-qualities-to-map-quality",
    "allele-informative-reads-overlap-margin",
    "max-mnp-distance",
    "recover-all-dangling-branches",
    "min-dangling-branch-length",
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
    // Mutect2's own.
    "normal-sample",
    "tumor-sample",
    "genotype-germline-sites",
    "genotype-germline-sites-fraction",
    "af-of-alleles-not-in-resource",
    "tumor-lod-to-emit",
    "initial-tumor-lod",
    "pcr-snv-qual",
    "pcr-indel-qual",
    "base-qual-correction-factor",
    "max-population-af",
    "downsampling-stride",
    "callable-depth",
    "max-suspicious-reads-per-alignment-start",
    "normal-lod",
    "ignore-itr-artifacts",
    "minimum-allele-fraction",
    "mitochondria-mode",
    // Somatic reference-confidence mode.
    "emit-ref-confidence",
    "gvcf-lod-band",
    // The read-threading assembler, shared with HaplotypeCaller.
    "allow-non-unique-kmers-in-ref",
    "dont-increase-kmer-sizes-for-cycles",
    "disable-adaptive-pruning",
    "kmer-size",
    "num-pruning-samples",
    "min-pruning",
    "adaptive-pruning-initial-error-rate",
    "pruning-lod-threshold",
    "pruning-seeding-lod-threshold",
    "max-unpruned-variants",
    "enable-legacy-graph-cycle-detection",
    "num-matching-bases-in-dangling-end-to-recover",
    "max-num-haplotypes-in-population",
    "disable-artificial-haplotype-recovery",
    "independent-mates",
    // The F1R2 counts, archived beside the stats.
    "f1r2-tar-gz",
    "f1r2-median-mq",
    "f1r2-min-bq",
    "f1r2-max-depth",
    // The germline resource and the panel of normals, read whole.
    "germline-resource",
    "panel-of-normals",
    // Without a Permutect dataset this changes nothing; the bamout is written; --debug-assembly
    // is logged only.
    "genotype-pon-sites",
    "permutect-dataset-mode",
    "bam-output",
    "bam-writer-type",
    "assembly-region-out",
    "capture-assembly-failure-bam",
    // Kept by the assembler and never printed: only HaplotypeCaller's engine prints it.
    "haplotype-debug-histogram-output",
    "debug-assembly-variants-out",
    // Declared for Mutect2 and read only by HaplotypeCaller's engine.
    "alm-path",
    "alm-interval",
    "pair-hmm-results-file",
    "debug-assembly",
    "graph-output",
    "debug-graph-transformations",
    "error-correct-reads",
    "kmer-length-for-read-error-correction",
    // Read only by the k-mer corrector, which is never given any k-mers.
    "min-observations-for-kmer-to-be-solid",
    "error-correction-log-odds",
    // Force-calling.
    "alleles",
    "force-call-filtered-alleles",
    // Without a native PairHMM or a DRAGstr model these change nothing.
    "native-pair-hmm-use-double-precision",
    "dont-use-dragstr-pair-hmm-scores",
    // The implementations, refused below when they name a kernel this port does not carry.
    "pair-hmm-implementation",
    "smith-waterman",
    "likelihood-calculation-engine",
];

fn limitation(what: &str) -> Thrown {
    Thrown::non_user(
        PORT_LIMITATION,
        format!("{what} This message is the port's own and not GATK's."),
    )
}

fn refuse_unported(parser: &Parser) -> Result<(), Thrown> {
    let declarations = gatk_tools::tool_declarations::declarations("Mutect2").unwrap_or(&[]);
    for definition in parser.definitions() {
        let name = definition.long_name();
        if !definition.has_been_set()
            || HONOURED.contains(&name)
            || super::FLOW_MATRIX_ARGUMENTS.contains(&name)
            || super::PILEUP_DETECTION_ARGUMENTS.contains(&name)
        {
            continue;
        }
        let read_filter_argument = declarations.iter().any(|declaration| {
            declaration.long_name == name
                && declaration.controlled_by == Some("GATKReadFilterPluginDescriptor")
        });
        if !read_filter_argument {
            return Err(limitation(&format!(
                "--{name} is a Mutect2 argument this port does not carry yet."
            )));
        }
    }
    let fraction = double_or(parser, "genotype-germline-sites-fraction", 1.0);
    if flag(parser, "genotype-germline-sites") && fraction > 0.0 && fraction < 1.0 {
        return Err(limitation(
            "--genotype-germline-sites-fraction below one draws from an unseeded generator.",
        ));
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

/// `Mutect2.doWork`.
/// `M2ArgumentCollection.getMitochondriaModeNameValuePairs`, each value as `String.valueOf` writes
/// it: the default allele frequency 4e-3, the emission and initial LODs 0, all dangling branches
/// recovered, a pruning threshold of `log10ToLog(-4)`, and `OriginalAlignment`.
fn mitochondria_mode() -> Vec<(&'static str, String)> {
    vec![
        ("af-of-alleles-not-in-resource", "0.004".to_string()),
        ("tumor-lod-to-emit", "0.0".to_string()),
        ("initial-tumor-lod", "0.0".to_string()),
        ("recover-all-dangling-branches", "true".to_string()),
        // `MathUtils.log10ToLog(-4)`, minus four times the natural logarithm of ten, as
        // `Double.toString` writes it.
        ("pruning-lod-threshold", "-9.210340371976184".to_string()),
        ("annotation", "OriginalAlignment".to_string()),
    ]
}

/// `Mutect2.customCommandLineValidation`: `--mitochondria-mode` sets its arguments where the
/// command line did not.
pub fn mutect2_modes(parser: &mut Parser) -> Result<(), Thrown> {
    if flag(parser, "mitochondria-mode") {
        let pairs = mitochondria_mode();
        let borrowed: Vec<(&str, &str)> = pairs.iter().map(|(k, v)| (*k, v.as_str())).collect();
        parser
            .set_mode_values(&borrowed)
            .map_err(|error| Thrown::non_user(error.class, error.message))?;
    }
    // `FlowMode`: STANDARD, or ADVANCED's pairs before STANDARD's.
    if let Some(mode) = scalar(parser, "flow-mode").filter(|mode| mode != "NONE") {
        let pairs = flow_mode_pairs(
            &mode,
            &[
                ("pruning-lod-threshold", "3.0"),
                (
                    "enable-dynamic-read-disqualification-for-genotyping",
                    "true",
                ),
                ("dynamic-read-disqualification-threshold", "10"),
                ("minimum-mapping-quality", "1"),
            ],
        );
        parser
            .set_mode_values(&pairs)
            .map_err(|error| Thrown::non_user(error.class, error.message))?;
    }
    Ok(())
}

pub fn mutect2(parser: &Parser) -> Outcome {
    refuse_unported(parser)?;
    match scalar(parser, "likelihood-calculation-engine").as_deref() {
        None | Some("PairHMM") | Some("FlowBased") | Some("FlowBasedHMM") => {}
        Some(other) => {
            return Err(limitation(&format!(
                "the {other} likelihood engine is not ported for Mutect2."
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
    let resolved = catalogue::resolve(
        &catalogue::AnnotationArguments {
            annotations: arguments(parser, "annotation"),
            groups: arguments(parser, "annotation-group"),
            excluded: arguments(parser, "annotations-to-exclude"),
            disable_tool_defaults: flag(parser, "disable-tool-default-annotations"),
            enable_all: flag(parser, "enable-all-annotations"),
        },
        &["StandardMutectAnnotation"],
        &[],
    )
    .map_err(|error| Thrown::command_line(error.message()))?;
    pedigree_refusal(parser, &resolved)?;
    let annotation_flow_order = flow_order_for_annotations(parser, &resolved)?;
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "Mutect2")?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let reference_path = argument(parser, "reference").ok_or_else(|| {
        Thrown::command_line("Argument reference was missing: Argument 'reference' is required")
    })?;
    let mut reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let samples = gatk_tools::haplotype_caller::samples_from_header(&header);
    let normal_samples = arguments(parser, "normal-sample");
    for normal in &normal_samples {
        if !samples.contains(normal) {
            return Err(Thrown::user(format!(
                "Bad input: Sample {normal} is not in BAM header: [{}]",
                samples.join(", ")
            )));
        }
    }

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

    let mut bases_of: std::collections::HashMap<String, Vec<u8>> = Default::default();
    let mut contig_bases = |contig: &str| -> Result<Vec<u8>, String> {
        if let Some(bases) = bases_of.get(contig) {
            return Ok(bases.clone());
        }
        let length = reference
            .sequences()
            .iter()
            .find(|(name, _)| name == contig)
            .map(|(_, length)| *length as i32)
            .ok_or_else(|| format!("contig {contig} is not in the reference"))?;
        let bases = reference
            .query(contig, 1, length)
            .map_err(|error| format!("{error:?}"))?;
        bases_of.insert(contig.to_string(), bases.clone());
        Ok(bases)
    };

    // The filters, then `makePostReadFilterTransformer`'s palindrome clipper.
    let filter = read_filter(parser, &filters, &header)?;
    let clip_itr = !flag(parser, "ignore-itr-artifacts");
    let mut reads = Vec::new();
    for read in source
        .iter_all()
        .map_err(|error| Thrown::user(format!("{error:?}")))?
    {
        if !filter(&read) {
            continue;
        }
        let read = if clip_itr {
            let contig = usize::try_from(read.reference_index)
                .ok()
                .and_then(|i| header.sequences.get(i))
                .map(|s| s.name.clone());
            let bases = match &contig {
                Some(name) => contig_bases(name).ok(),
                None => None,
            };
            m2::palindrome_artifact_clip(&read, &header, bases.as_deref())
                .map_err(|e| Thrown::non_user(PORT_FAILURE, e))?
        } else {
            read
        };
        reads.push(read);
    }

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
    let sw = super::haplotype_caller::smith_waterman_arguments(parser);
    let log10_to_log = gatk_engine::allele_likelihoods::log10_to_log;
    let default_af = {
        let given = double_or(parser, "af-of-alleles-not-in-resource", -1.0);
        if given >= 0.0 {
            given
        } else if normal_samples.is_empty() {
            DEFAULT_AF_FOR_TUMOR_ONLY_CALLING
        } else {
            DEFAULT_AF_FOR_TUMOR_NORMAL_CALLING
        }
    };
    let genotype_germline_sites = flag(parser, "genotype-germline-sites");
    // `getInitialLogOdds()` and `getEmissionLogOdds()`: minus infinity in reference-confidence mode.
    let emit_mode = scalar(parser, "emit-ref-confidence").unwrap_or_else(|| "NONE".to_string());
    let emit_reference_confidence = emit_mode != "NONE";
    let min_af = double_or(parser, "minimum-allele-fraction", 0.0);
    super::haplotype_caller::artificial_haplotype_recovery_refusal(parser)?;
    let arguments = Mutect2Arguments {
        given_alleles: super::given_alleles(parser)?,
        assembly_region_out: argument(parser, "assembly-region-out").map(std::path::PathBuf::from),
        bamout: argument(parser, "bam-output").map(|path| {
            let (level, deflater) = super::output_compression(parser);
            gatk_tools::haplotype_caller::BamoutArguments {
                kind: scalar(parser, "bam-writer-type")
                    .and_then(|name| gatk_tools::haplotype_bam_writer::WriterType::parse(&name))
                    .unwrap_or(gatk_tools::haplotype_bam_writer::WriterType::CalledHaplotypes),
                path: std::path::PathBuf::from(path),
                create_index: flag(parser, "create-output-bam-index"),
                create_md5: flag(parser, "create-output-bam-md5"),
                level,
                deflater,
            }
        }),
        pileup_detection: super::pileup_detection_arguments(parser),
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
        // `assembleReads(..., correctOverlappingBaseQualities = false)`.
        assembly: RegionAssemblyArguments {
            min_base_quality_score: number_or(parser, "min-base-quality-score", 10) as u8,
            dont_use_soft_clipped_bases: flag(parser, "dont-use-soft-clipped-bases"),
            soft_clip_low_quality_ends: flag(parser, "soft-clip-low-quality-ends"),
            override_softclip_fragment_check: flag(parser, "override-fragment-softclip-check"),
            do_not_correct_overlapping_base_qualities: true,
            dangling_end_sw: sw.0,
            haplotype_to_reference_sw: sw.1,
            read_to_haplotype_sw: sw.2,
            capture_assembly_failure_bam: flag(parser, "capture-assembly-failure-bam")
                .then(|| super::output_compression(parser)),
            error_correct_reads: flag(parser, "error-correct-reads"),
            kmer_length_for_read_error_correction: number_or(
                parser,
                "kmer-length-for-read-error-correction",
                25,
            ),
        },
        flow_engine: flow_engine_arguments(parser),
        flow_hmm: flow_hmm_arguments(parser, pcr_model(parser), None),
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
            modify_soft_clipped_bases: true,
            dragstr_params: None,
            results_file: argument(parser, "pair-hmm-results-file").map(std::path::PathBuf::from),
        },
        normal_samples: normal_samples.clone(),
        initial_log_odds: if emit_reference_confidence {
            f64::NEG_INFINITY
        } else {
            log10_to_log(double_or(parser, "initial-tumor-lod", 2.0))
        },
        emission_log_odds: if emit_reference_confidence {
            f64::NEG_INFINITY
        } else {
            log10_to_log(double_or(parser, "tumor-lod-to-emit", 3.0))
        },
        normal_log_odds: log10_to_log(double_or(parser, "normal-lod", 2.2)),
        pcr_snv_qual: number_or(parser, "pcr-snv-qual", 40),
        pcr_indel_qual: number_or(parser, "pcr-indel-qual", 40),
        multiple_substitution_base_qual_correction: double_or(
            parser,
            "base-qual-correction-factor",
            5.0,
        ),
        callable_depth: number_or(parser, "callable-depth", 10).max(0) as usize,
        genotype_germline_sites,
        consider_germline_active: genotype_germline_sites
            && double_or(parser, "genotype-germline-sites-fraction", 1.0) >= 1.0,
        min_af,
        flow_collapse_hmer_size: match number_or(parser, "flow-assembly-collapse-hmer-size", 0) {
            -1 => gatk_tools::haplotype_collapsing::determine_collapse_threshold(&header)
                .map_err(|(class, message)| Thrown::non_user(class, message))?,
            size => size,
        },
        flow_collapse_partial_mode: flag(parser, "flow-assembly-collapse-partial-mode"),
        filter_alleles: super::haplotype_caller::allele_filter_arguments(parser),
        default_af,
        germline_resource: match argument(parser, "germline-resource") {
            Some(path) => resource_records(&path)?,
            None => Vec::new(),
        },
        max_population_af: double_or(parser, "max-population-af", 0.01),
        panel_of_normals: match argument(parser, "panel-of-normals") {
            Some(path) => resource_records(&path)?,
            None => Vec::new(),
        },
        genotype_pon_sites: flag(parser, "genotype-pon-sites"),
        f1r2: argument(parser, "f1r2-tar-gz").map(|_| gatk_tools::collect_f1r2_counts::Args {
            min_median_map_qual: number_or(parser, "f1r2-median-mq", 50),
            min_base_quality: number_or(parser, "f1r2-min-bq", 20),
            max_depth: number_or(
                parser,
                "f1r2-max-depth",
                gatk_tools::collect_f1r2_counts::DEFAULT_MAX_DEPTH,
            ),
        }),
        informative_read_overlap_margin: number_or(
            parser,
            "allele-informative-reads-overlap-margin",
            2,
        ),
        max_mnp_distance: number_or(parser, "max-mnp-distance", 1),
        phred_scaled_global_read_mismapping_rate: mismapping,
        max_suspicious_reads_per_alignment_start: number_or(
            parser,
            "max-suspicious-reads-per-alignment-start",
            0,
        ),
        downsampling_stride: number_or(parser, "downsampling-stride", 1),
        assembler: {
            let mut assembler = super::haplotype_caller::assembler_arguments(
                parser,
                !flag(parser, "disable-adaptive-pruning"),
            );
            assembler.recover_dangling_branches = true;
            assembler
        },
        independent_mates: flag(parser, "independent-mates"),
        emit_reference_confidence,
        somatic_ref_confidence:
            gatk_tools::reference_confidence_model::SomaticRefConfidence::from_min_af(min_af),
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
        flow_order_for_annotations: annotation_flow_order.clone(),
    };
    // `onTraversalStart`, once the VCF writer exists.
    let f1r2_tar_gz = argument(parser, "f1r2-tar-gz");
    if let Some(path) = &f1r2_tar_gz {
        let absolute = java_absolute_path(path);
        if !absolute.ends_with(".tar.gz") {
            return Err(Thrown::user(format!(
                "Couldn't write file {absolute} because f1r2-tar-gz file must end in .tar.gz"
            )));
        }
    }
    let mut random = gatk_random();
    let assembly_variants = argument(parser, "debug-assembly-variants-out").map(|path| {
        (
            path,
            std::cell::RefCell::new(gatk_tools::assembly_variants_out::AssemblyVariantsOut::new(
                header.sequences.iter().map(|s| s.name.clone()).collect(),
            )),
        )
    });
    let called = m2::call_variants(
        &reads,
        &header,
        &traversal,
        &arguments,
        &annotator,
        &mut contig_bases,
        &mut random,
        assembly_variants.as_ref().map(|(_, queue)| queue),
    );
    // `--debug-assembly-variants-out`: an empty header, the queue, indexed on the fly.
    if let Some((path, queue)) = &assembly_variants {
        let empty = htsjdk_vcf::header::VcfHeader {
            lines: Vec::new(),
            samples: Vec::new(),
        };
        let text = super::write_vcf_strict_or_lenient(false, &empty, &queue.borrow().records())?;
        super::write_variant_file(
            parser,
            path,
            &text,
            true,
            false,
            Some(
                header
                    .sequences
                    .iter()
                    .map(|s| (s.name.clone(), s.length))
                    .collect(),
            ),
        )?;
    }
    drop(random);
    // A failure mid-traversal still writes the calls before it, as the writer closes.
    let mut f1r2_files = None;
    let (called_calls, callable_sites, failure) = match called {
        Ok(output) => {
            f1r2_files = output.f1r2_files;
            (output.calls, output.callable_sites, None)
        }
        Err(error) if error.class == "PortLimitation" => return Err(limitation(&error.message)),
        Err(error) => (
            error.calls_before,
            0,
            // A `UserException` is reported as a user error, by its message alone.
            Some(if error.class.contains("UserException") {
                Thrown::user(error.message)
            } else {
                Thrown::non_user(java_class_name(&error.class), error.message)
            }),
        ),
    };

    // `Mutect2Engine.writeHeader`.
    let compound = |key: &str, id: &str, number: Cardinality, line_type: LineType, text: &str| {
        HeaderLine::Compound {
            key: key.to_string(),
            id: id.to_string(),
            number,
            line_type,
            description: text.to_string(),
            extra: Vec::new(),
        }
    };
    let mut lines: Vec<HeaderLine> = vec![
        HeaderLine::Unstructured {
            key: "MutectVersion".to_string(),
            value: "2.2".to_string(),
        },
        HeaderLine::Unstructured {
            key: "filtering_status".to_string(),
            value: "Warning: unfiltered Mutect 2 calls.  Please run FilterMutectCalls to remove \
                    false positives."
                .to_string(),
        },
    ];
    lines.extend(catalogue::descriptions(&resolved, false, false));
    lines.extend(default_tool_vcf_header_lines(parser, "Mutect2"));
    use Cardinality::{Fixed, A};
    use LineType::{Flag, Float, Integer, String as Text};
    for (id, number, line_type, text) in [
        ("NLOD", A, Float, "Normal log 10 likelihood ratio of diploid het or hom alt genotypes"),
        ("TLOD", A, Float, "Log 10 likelihood ratio score of variant existing versus not existing"),
        ("NALOD", A, Float, "Log 10 odds of artifact in normal with same allele fraction as tumor"),
        ("ECNT", Fixed(1), Integer, "Number of potential somatic events in the assembly region"),
        ("ECNTH", A, Integer, "Number of somatic events in best supporting haplotype for each alt allele"),
        ("PON", Fixed(0), Flag, "site found in panel of normals"),
        ("POPAF", A, Float, "negative log 10 population allele frequencies of alt alleles"),
        ("GERMQ", Fixed(1), Integer, "Phred-scaled quality that alt alleles are not germline variants"),
        ("CONTQ", Fixed(1), Float, "Phred-scaled qualities that alt allele are not due to contamination"),
        ("SEQQ", Fixed(1), Integer, "Phred-scaled quality that alt alleles are not sequencing errors"),
        ("STRQ", Fixed(1), Integer, "Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors"),
        ("ROQ", Fixed(1), Float, "Phred-scaled qualities that alt allele are not due to read orientation artifact"),
        ("STRANDQ", Fixed(1), Integer, "Phred-scaled quality of strand bias artifact"),
        ("OCM", Fixed(1), Integer, "Number of alt reads whose original alignment doesn't match the current contig."),
        ("NCount", Fixed(1), Integer, "Count of N bases in the pileup"),
        ("AS_UNIQ_ALT_READ_COUNT", A, Integer, "Number of reads with unique start and mate end positions for each alt at a variant site"),
    ] {
        lines.push(compound("INFO", id, number, line_type, text));
    }
    for id in ["GT", "AD", "GQ", "DP", "PL", "PS"] {
        lines.extend(htsjdk_vcf::standard_header_lines::standard_format_line(id));
    }
    lines.push(compound(
        "FORMAT",
        "AF",
        A,
        Float,
        "Allele fractions of alternate alleles in the tumor",
    ));
    lines.push(compound("FORMAT", "PID", Fixed(1), Text, "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group"));
    lines.push(compound("FORMAT", "PGT", Fixed(1), Text, "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles"));
    for sample in samples.iter().filter(|s| !normal_samples.contains(s)) {
        lines.push(HeaderLine::Unstructured {
            key: "tumor_sample".to_string(),
            value: sample.clone(),
        });
    }
    for sample in &normal_samples {
        lines.push(HeaderLine::Unstructured {
            key: "normal_sample".to_string(),
            value: sample.clone(),
        });
    }
    // `referenceConfidenceModel.getVCFHeaderLines()` and the genotypes' `TLOD`.
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
        lines.push(compound(
            "FORMAT",
            "TLOD",
            A,
            Float,
            "Log 10 likelihood ratio score of variant existing versus not existing",
        ));
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
    let mut vcf_header = VcfHeader { lines, samples };
    let mut written = called_calls;
    // `SomaticGVCFWriter`: the reference records combined into blocks by TLOD band.
    if emit_mode == "GVCF" {
        let bands: Vec<f64> = {
            let given: Vec<f64> = super::arguments(parser, "gvcf-lod-band")
                .iter()
                .filter_map(|band| band.trim().parse().ok())
                .collect();
            if given.is_empty() {
                // `for (double i = -2.5; i <= 1; i = i + 0.5)`.
                let mut defaults = Vec::new();
                let mut i = -2.5f64;
                while i <= 1.0 {
                    defaults.push(i);
                    i += 0.5;
                }
                defaults
            } else {
                given
            }
        };
        let (partitions, precision) = gatk_tools::gvcf_blocks::parse_lod_partitions(&bands)
            .map_err(|message| Thrown::non_user("java.lang.IllegalArgumentException", message))?;
        gatk_tools::gvcf_blocks::add_ranges_to_header(&mut vcf_header, &partitions);
        let mut combiner =
            gatk_tools::gvcf_blocks::SomaticBlockCombiner::new(partitions, precision);
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
    let text = write_vcf_honouring_lenient(parser, &vcf_header, &written)?;
    write_variant_output(parser, &output, &text)?;
    if let Some(failure) = failure {
        return Err(failure);
    }
    // `MutectStats.writeToFile`, beside the output.
    let stats = format!(
        "statistic\tvalue\ncallable\t{}\n",
        gatk_engine::tsv_table::java_double_to_string(callable_sites as f64)
    );
    std::fs::write(format!("{output}.stats"), stats)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}.stats: {error}")))?;
    // `writeHistograms` and `closeAndArchiveFiles`: the members `./<sample><extension>`, which the
    // reference lists in the file system's order and the port sorts.
    if let (Some(path), Some(mut files)) = (&f1r2_tar_gz, f1r2_files) {
        files.sort_by(|left, right| left.0.cmp(&right.0));
        let entries: Vec<(String, Vec<u8>)> = files
            .into_iter()
            .map(|(name, text)| (format!("./{name}"), text.into_bytes()))
            .collect();
        write_file(path, &tar_gz(&entries))?;
    }
    // `onTraversalSuccess` returns "SUCCESS", which `handleResult` prints.
    Ok(Some("SUCCESS".to_string()))
}

/// A feature input the engine queries, read whole: a path that is not there is refused as
/// `GATKTool` refuses it, by its URI.
fn resource_records(path: &str) -> Result<Vec<htsjdk_vcf::variant::VariantContext>, Thrown> {
    if std::fs::metadata(path).is_err() {
        return Err(Thrown {
            failure: Failure::User,
            exception:
                "org.broadinstitute.hellbender.exceptions.UserException$CouldNotReadInputFile",
            message: Some(format!(
                "Couldn't read file file://{}. Error was: It doesn't exist.",
                java_absolute_path(path)
            )),
        });
    }
    Ok(feature_variants(path)?.0)
}
