//! `HaplotypeCaller` in VCF mode: the `AssemblyRegionWalker` startup around
//! [`gatk_tools::haplotype_caller::call_variants`], and the header `makeVCFHeader` writes.
//!
//! The reads are the input's, through `IUPACReadTransformer(true)` and then the resolved read
//! filters (`makeStandardHCReadFilters` by default); the traversal intervals are `-L`'s, or every
//! contig of the reads' dictionary. Every argument the port does not honour is refused when it is
//! set, rather than silently left at its default: GVCF output, `--alleles`, the assembler's own
//! knobs, pileup detection, DRAGEN mode and the flow-based paths among them.

use super::*;
use gatk_annotation::catalogue;
use gatk_engine::allele_frequency_calculator::Priors;
use gatk_engine::assembly_region_iterator::AssemblyRegionArgs;
use gatk_engine::assembly_region_trimmer::TrimmerArguments;
use gatk_engine::pair_hmm_likelihood_engine::{LikelihoodEngineArguments, PcrErrorModel};
use gatk_tools::genotyping_engine::{Configuration, SubsetMethod};
use gatk_tools::haplotype_caller::{self as hc, HaplotypeCallerArguments};
use gatk_tools::haplotype_caller_engine::IsActiveArguments;
use gatk_tools::hc_genotyping::HcGenotypingArguments;
use gatk_tools::hc_region::RegionAssemblyArguments;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};

/// The tool's arguments this port reads. Any other argument the tool declares, apart from the
/// read-filter plugin's, is refused when the command line sets it.
const HONOURED: &[&str] = &[
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
];

fn limitation(what: &str) -> Thrown {
    Thrown::non_user(
        PORT_LIMITATION,
        format!("{what} This message is the port's own and not GATK's."),
    )
}

/// The refusal of any declared argument the port does not read.
fn refuse_unported(parser: &Parser) -> Result<(), Thrown> {
    let declarations =
        gatk_tools::tool_declarations::declarations("HaplotypeCaller").unwrap_or(&[]);
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
                "--{name} is a HaplotypeCaller argument this port does not carry yet."
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

/// `HaplotypeCaller.doWork` in VCF mode.
pub fn haplotype_caller(parser: &Parser) -> Outcome {
    refuse_unported(parser)?;
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
    let ReadWalkerStart {
        source,
        header,
        intervals,
        filters,
    } = read_walker_startup(parser, "HaplotypeCaller")?;
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
    let mut reads = Vec::new();
    for mut read in source
        .iter_all()
        .map_err(|error| Thrown::user(format!("{error:?}")))?
    {
        hc::iupac_to_n(&mut read).map_err(Thrown::user)?;
        if filter(&read) {
            reads.push(read);
        }
    }

    let ploidy = number_or(parser, "sample-ploidy", 2).max(0) as usize;
    let calling_confidence = double_or(
        parser,
        "standard-min-confidence-threshold-for-calling",
        30.0,
    );
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
            dragstr_params: None,
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
            do_physical_phasing: false,
        },
        calling: Configuration {
            standard_confidence_for_calling: calling_confidence,
            max_alternate_alleles,
            sample_ploidy: ploidy,
            annotate_number_of_alleles_discovered: flag(
                parser,
                "annotate-with-num-discovered-alleles",
            ),
            emit_all_active_sites: false,
            allele_specific: false,
            emit_all_confident_sites: false,
            annotate_all_sites_with_pls: false,
            force_keep_all_alleles: false,
            assignment_method: SubsetMethod::UsePlsToAssign,
        },
        priors,
        mapping_quality_threshold: number_or(parser, "mapping-quality-threshold-for-genotyping", 20)
            as u8,
    };
    let annotator = gatk_tools::variant_annotator_engine::Engine {
        resolved: resolved.clone(),
        overlap_names: Vec::new(),
        expressions: Vec::new(),
        allele_concordance: false,
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
    let mut written =
        called.map_err(|error| Thrown::non_user(java_class_name(&error.class), error.message))?;

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
    let mut lines: Vec<HeaderLine> = default_tool_vcf_header_lines(parser, "HaplotypeCaller");
    if arguments.calling.annotate_number_of_alleles_discovered {
        lines.push(compound(
            "INFO",
            "NDA",
            Cardinality::Fixed(1),
            LineType::Integer,
            "Number of alternate alleles discovered (but not necessarily genotyped) at this site",
        ));
    }
    lines.extend(catalogue::descriptions(&resolved, false, false));
    lines.push(compound("INFO", "MLEAC", Cardinality::A, LineType::Integer, "Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed"));
    lines.push(compound("INFO", "MLEAF", Cardinality::A, LineType::Float, "Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed"));
    for id in ["GT", "GQ", "DP", "PL"] {
        lines.extend(htsjdk_vcf::standard_header_lines::standard_format_line(id));
    }
    lines.push(HeaderLine::Filter {
        id: "LowQual".to_string(),
        description: "Low quality".to_string(),
    });
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
        samples: hc::samples_from_header(&header),
    };
    apply_sites_only(parser, &mut vcf_header, &mut written);
    let text = write_vcf_honouring_lenient(parser, &vcf_header, &written)?;
    write_variant_output(parser, &output, &text)?;
    Ok(None)
}
