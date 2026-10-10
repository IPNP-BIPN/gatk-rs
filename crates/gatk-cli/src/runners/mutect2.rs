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
        if !definition.has_been_set() || HONOURED.contains(&name) {
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
pub fn mutect2(parser: &Parser) -> Outcome {
    refuse_unported(parser)?;
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
    let arguments = Mutect2Arguments {
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
            modify_soft_clipped_bases: true,
            dragstr_params: None,
        },
        normal_samples: normal_samples.clone(),
        initial_log_odds: log10_to_log(double_or(parser, "initial-tumor-lod", 2.0)),
        emission_log_odds: log10_to_log(double_or(parser, "tumor-lod-to-emit", 3.0)),
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
        min_af: double_or(parser, "minimum-allele-fraction", 0.0),
        default_af,
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
        recover_all_dangling_branches: flag(parser, "recover-all-dangling-branches"),
        min_dangling_branch_length: number_or(parser, "min-dangling-branch-length", 4),
    };
    let annotator = gatk_tools::variant_annotator_engine::Engine {
        resolved: resolved.clone(),
        overlap_names: Vec::new(),
        expressions: Vec::new(),
        allele_concordance: false,
    };
    let mut random = gatk_random();
    let called = m2::call_variants(
        &reads,
        &header,
        &traversal,
        &arguments,
        &annotator,
        &mut contig_bases,
        &mut random,
    );
    drop(random);
    let called = called.map_err(|error| {
        if error.class == "PortLimitation" {
            limitation(&error.message)
        } else {
            Thrown::non_user(java_class_name(&error.class), error.message)
        }
    })?;

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
    let mut written = called.calls;
    apply_sites_only(parser, &mut vcf_header, &mut written);
    let text = write_vcf_honouring_lenient(parser, &vcf_header, &written)?;
    write_variant_output(parser, &output, &text)?;
    // `MutectStats.writeToFile`, beside the output.
    let stats = format!(
        "statistic\tvalue\ncallable\t{}\n",
        gatk_engine::tsv_table::java_double_to_string(called.callable_sites as f64)
    );
    std::fs::write(format!("{output}.stats"), stats)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}.stats: {error}")))?;
    Ok(None)
}
