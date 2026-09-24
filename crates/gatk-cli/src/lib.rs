//! The dispatcher: what `gatk <Tool> <args>` does before a tool ever sees an argument.
//!
//! `Main.mainEntry` is three steps and an exit status: route the command line, run whatever the
//! route named, and turn whatever came back into a status. The routing and the statuses are ported
//! in [`gatk_tools::main_entry`]; what is here is the entry point that puts them together, writes
//! to the two streams the reference writes to, and returns the status the reference would have
//! exited with.
//!
//! # What this does not do yet
//!
//! One thing, and it is a deliberate boundary rather than an omission: **a WALKER's own usage**
//! carries a conditional block per read filter the plugin descriptor discovered, and neither the
//! blocks nor their order is measured, so [`tool_usage`] answers for a tool with no controlled
//! argument and nothing for the rest. Its command line does reach the parser: the plugin trim runs
//! over the ownership table measured in `plugin-argument-ownership`.
//!
//! The main usage listing WAS the other one. It is three hundred and seventy-three lines of tool
//! names, program groups and one-line summaries, and none of them was in a golden; they are now
//! (`main-usage`), and [`usage`] reproduces the listing line for line.
//!
//! Ported from `org.broadinstitute.hellbender.Main`.

pub mod command_line;
pub mod definitions;
pub mod runners;
pub mod variant_bridge;

use gatk_tools::main_entry::{self, Failure, Route, Stream, Thrown};

/// The reference's own pins, which `printVersionInfo` reads off the jar's manifest.
pub const TOOLKIT_NAME: &str = "The Genome Analysis Toolkit (GATK)";
pub const TOOLKIT_VERSION: &str = "4.6.2.0";
pub const HTSJDK_VERSION: &str = "4.2.0";
pub const PICARD_VERSION: &str = "3.4.0";

/// `getCommandLineName()`, which GATK leaves empty: the first line of the usage therefore reads
/// `USAGE:  <program name>`, with two spaces where a toolkit's name would have been.
pub const COMMAND_LINE_NAME: &str = "";

/// The first line of the main usage, which Barclay writes in colour.
///
/// The escapes are the reference's: bold red for the label, green for the program name. They are
/// bytes of the output like any other, which is why they are written out rather than described.
pub const USAGE_FIRST_LINE: &str =
    "\u{1b}[1m\u{1b}[31mUSAGE:  \u{1b}[32m<program name>\u{1b}[1m\u{1b}[31m [-h]";

/// What one run wrote and what it returned, which is what `mainEntry` turns into an exit status.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Outcome {
    pub stdout: String,
    pub stderr: String,
    pub status: i32,
}

/// A tool the port can actually run, and what it does with the arguments after its name.
///
/// The list is empty of walkers on purpose: a tool becomes reachable here the moment its name
/// resolves, and runnable only once its arguments have a declaration to parse against. The gap
/// between the two is Milestone C's, and it is visible rather than papered over.
pub type Runner = fn(&[String]) -> Result<Option<String>, Thrown>;

/// `Main.instanceMain`, with the two streams returned rather than written.
pub fn run(args: &[String]) -> Outcome {
    let mut stdout = String::new();
    let mut stderr = String::new();
    match main_entry::route(args) {
        Route::Usage(stream) => {
            // `handleResult(null)` prints nothing, and the run is a success: asking for help is
            // not an error, which is the one path that reaches the runner with no program.
            let target = match stream {
                Stream::Stdout => &mut stdout,
                Stream::Stderr => &mut stderr,
            };
            target.push_str(&usage());
            Outcome {
                stdout,
                stderr,
                status: 0,
            }
        }
        Route::Version => {
            let (first, rest) = main_entry::version_lines(
                TOOLKIT_NAME,
                TOOLKIT_VERSION,
                HTSJDK_VERSION,
                PICARD_VERSION,
            );
            // The first line goes to `System.out` whatever stream the method was handed, and the
            // two under it go to the stream. Here both are stdout, which is what the command line
            // asks for, and the split is only visible when a caller hands it something else.
            stdout.push_str(&first);
            stdout.push_str(&rest);
            Outcome {
                stdout,
                stderr,
                status: 0,
            }
        }
        Route::Unknown { message } => {
            // The refusal is the usage AND the message: the usage goes to stderr first, and the
            // UserException follows it through `handleUserException`.
            stderr.push_str(&usage());
            stderr.push_str(&main_entry::user_exception_report(&message));
            Outcome {
                stdout,
                stderr,
                status: main_entry::exit_status(Failure::User),
            }
        }
        Route::Tool { name } => {
            // The parse comes before the run, and it only happens for a tool whose whole argument
            // surface converts: a parser missing a third of its arguments would refuse a command
            // line the reference accepts, which is a worse answer than refusing to parse at all.
            // A tool asked for its own help answers with its usage and returns nothing, which is
            // a success: `-h` after a tool name is the tool's argument and not the dispatcher's.
            let tool_args = main_entry::tool_arguments(args);
            let help = tool_args.iter().any(|arg| arg == "-h" || arg == "--help");
            if help {
                if let Some(usage) = tool_usage(&name) {
                    stderr.push_str(&usage);
                    return Outcome {
                        stdout,
                        stderr,
                        status: 0,
                    };
                }
                // A tool whose usage this port cannot lay out is not refused for asking: the
                // reference answers `-h` with the usage and a zero status, and a parse refusal
                // about an unrecognised `-h` would be a message the reference never writes.
            }
            if let Some(error) = (!help).then(|| parse_failure(&name, tool_args)).flatten() {
                // `mainEntry` prints the PROGRAM's usage before the message on this path, which
                // the port now does for a tool whose usage it can lay out.
                if let Some(usage) = tool_usage(&name) {
                    stderr.push_str(&usage);
                }
                stderr.push_str(&main_entry::user_exception_report(&error));
                return Outcome {
                    stdout,
                    stderr,
                    status: main_entry::exit_status(Failure::CommandLine),
                };
            }
            match runner(&name) {
                None => {
                    stderr.push_str(&main_entry::user_exception_report(&not_ported(&name)));
                    Outcome {
                        stdout,
                        stderr,
                        status: main_entry::exit_status(Failure::User),
                    }
                }
                Some(runner) => match runner(main_entry::tool_arguments(args)) {
                    Ok(result) => {
                        if let Some(printed) = main_entry::tool_returned(result.as_deref()) {
                            stdout.push_str(&printed);
                            stdout.push('\n');
                        }
                        Outcome {
                            stdout,
                            stderr,
                            status: 0,
                        }
                    }
                    Err(thrown) => {
                        // Which handler this reaches is the throwable's own answer: a
                        // `UserException` is decorated with the banner, and anything else is
                        // printed as `handleNonUserException` prints it, which is the exception's
                        // CLASS in front of its message and no banner at all (#1020).
                        stderr.push_str(&thrown.report());
                        Outcome {
                            stdout,
                            stderr,
                            status: thrown.status(),
                        }
                    }
                },
            }
        }
    }
}

/// Whether this port can hand a tool's command line to the ported Barclay parser at all.
///
/// Two conditions. Every argument the tool declares has to convert to a definition, which since
/// the value classes were measured is true of all of them. And every argument a PLUGIN DESCRIPTOR
/// controls has to be one the ownership table names, because the trim that removes an unselected
/// plugin's arguments before the required check needs to know which plugin declared each one. The
/// table is measured in `plugin-argument-ownership`, and it covers all twenty-eight of the read
/// filter arguments the declared tools carry, so a walker is now parseable.
///
/// What is still missing is narrower than the tool: a tool hands its descriptor a list of DEFAULT
/// filters, whose arguments are allowed with no `--read-filter` on the command line, and the port
/// has no per-tool list to hand [`gatk_barclay::Parser::with_default_plugins`]. It refuses a
/// command line that sets a default filter's own argument without naming the filter, which the
/// reference accepts. None of the default filters of the tools declared here carries an argument,
/// so no command line of theirs is affected. See gatk-rs#987.
pub fn parseable(tool: &str) -> bool {
    gatk_tools::tool_declarations::declarations(tool)
        .map(|list| {
            definitions::missing(list).is_empty()
                && list
                    .iter()
                    .filter(|declaration| declaration.controlled_by.is_some())
                    .all(|declaration| {
                        gatk_tools::plugin_ownership::owner(declaration.long_name).is_some()
                    })
        })
        .unwrap_or(false)
}

/// Whether a tool's usage can be laid out from its declarations alone.
///
/// It can, for every declared tool. A walker's conditional blocks are one per read filter that
/// declares an argument, in the ownership table's order; the two arguments the descriptor answers
/// for print the catalogue and the tool's own defaults; and the mutex sentence names the target
/// definition's FIELD, which `mutex-target-names` measured. `CountReads`'s two hundred and
/// ninety-seven lines are the reference's, and the test beside this file compares them as one
/// string rather than counting how many agree.
pub fn usage_composable(tool: &str) -> bool {
    parseable(tool)
}

/// A tool's usage as this port composes it, whatever the dispatcher does with it.
///
/// [`tool_usage`] is what the dispatcher answers `-h` with and is gated on [`parseable`]; this is
/// the composition itself, which a test compares against the golden for every declared tool.
pub fn composed_usage(tool: &str) -> Option<String> {
    let list = gatk_tools::tool_declarations::declarations(tool)?;
    let summary = gatk_tools::tool_declarations::summary(tool)?;
    let (required, optional, advanced) = gatk_tools::usage_text::sections_for(Some(tool), list);
    let conditional = gatk_tools::usage_text::conditionals(list);
    Some(gatk_tools::usage_text::render(
        tool,
        summary,
        TOOLKIT_VERSION,
        gatk_tools::tool_declarations::maturity(tool),
        &required,
        &optional,
        &advanced,
        &conditional,
    ))
}

/// A tool's own usage, for a tool whose whole argument surface the port can lay out.
pub fn tool_usage(tool: &str) -> Option<String> {
    if !usage_composable(tool) {
        return None;
    }
    composed_usage(tool)
}

/// The refusal the ported parser makes of a command line, if the tool is parseable at all.
pub fn parse_failure(tool: &str, args: &[String]) -> Option<String> {
    if !parseable(tool) {
        return None;
    }
    let list = gatk_tools::tool_declarations::declarations(tool)?;
    let mut parser = parser_for(tool, list);
    let borrowed: Vec<&str> = args.iter().map(String::as_str).collect();
    parser
        .parse_arguments(&borrowed)
        .err()
        .map(|error| error.message)
}

/// The tools this port can run, which is three of them.
///
/// A name that resolves and has no runner is not a name that does not resolve: the reference would
/// have run it, and saying so is the honest answer. What it is NOT is a refusal the reference
/// makes, which is why [`not_ported`] says whose refusal it is.
pub fn runner(name: &str) -> Option<Runner> {
    match name {
        "CountReads" => Some(run_count_reads),
        "CountVariants" => Some(run_count_variants),
        "CreateHadoopBamSplittingIndex" => Some(run_create_hadoop_bam_splitting_index),
        "PrintReads" => Some(run_print_reads),
        "GatherVcfsCloud" => Some(run_gather_vcfs_cloud),
        "ApplyBQSR" => Some(run_apply_bqsr),
        "IndexFeatureFile" => Some(run_index_feature_file),
        "PrintBGZFBlockInformation" => Some(run_print_bgzf_block_information),
        "CountBases" => Some(run_count_bases),
        "FlagStat" => Some(run_flag_stat),
        "CountBasesInReference" => Some(run_count_bases_in_reference),
        "SplitIntervals" => Some(run_split_intervals),
        "PreprocessIntervals" => Some(run_preprocess_intervals),
        "Pileup" => Some(run_pileup),
        "CheckPileup" => Some(run_check_pileup),
        "CalculateContamination" => Some(run_calculate_contamination),
        "PrintReadsHeader" => Some(run_print_reads_header),
        "GetPileupSummaries" => Some(run_get_pileup_summaries),
        "UnmarkDuplicates" => Some(run_unmark_duplicates),
        "FastaReferenceMaker" => Some(run_fasta_reference_maker),
        "CheckReferenceCompatibility" => Some(run_check_reference_compatibility),
        "AnnotateVcfWithExpectedAlleleFraction" => {
            Some(run_annotate_vcf_with_expected_allele_fraction)
        }
        "CalculateMixingFractions" => Some(run_calculate_mixing_fractions),
        "AnnotateVcfWithBamDepth" => Some(run_annotate_vcf_with_bam_depth),
        "CountFalsePositives" => Some(run_count_false_positives),
        "EvaluateInfoFieldConcordance" => Some(run_evaluate_info_field_concordance),
        "CallCopyRatioSegments" => Some(run_call_copy_ratio_segments),
        "VariantFiltration" => Some(run_variant_filtration),
        "CollectAllelicCounts" => Some(run_collect_allelic_counts),
        "FilterIntervals" => Some(run_filter_intervals),
        "PathSeqBuildReferenceTaxonomy" => Some(run_path_seq_build_reference_taxonomy),
        "PathSeqBuildKmers" => Some(run_path_seq_build_kmers),
        "CondenseDepthEvidence" => Some(run_condense_depth_evidence),
        "PrintSVEvidence" => Some(run_print_sv_evidence),
        "SiteDepthtoBAF" => Some(run_site_depth_to_baf),
        "MergeMutectStats" => Some(run_merge_mutect_stats),
        "GatherPileupSummaries" => Some(run_gather_pileup_summaries),
        "GatherNormalArtifactData" => Some(run_gather_normal_artifact_data),
        "MergeAnnotatedRegions" => Some(run_merge_annotated_regions),
        "MergeAnnotatedRegionsByAnnotation" => Some(run_merge_annotated_regions_by_annotation),
        "TagGermlineEvents" => Some(run_tag_germline_events),
        "GatherTranches" => Some(run_gather_tranches),
        "NuMTFilterTool" => Some(run_numt_filter_tool),
        "MTLowHeteroplasmyFilterTool" => Some(run_mt_low_heteroplasmy_filter_tool),
        "ExampleMultiFeatureWalker" => Some(run_example_multi_feature_walker),
        "CalculateAverageCombinedAnnotations" => Some(run_calculate_average_combined_annotations),
        "FilterVariantTranches" => Some(run_filter_variant_tranches),
        "ApplyVQSR" => Some(run_apply_vqsr),
        "VCFComparator" => Some(run_vcf_comparator),
        "SVStratify" => Some(run_sv_stratify),
        "SVConcordance" => Some(run_sv_concordance),
        "SVCluster" => Some(run_sv_cluster),
        "GroupedSVCluster" => Some(run_grouped_sv_cluster),
        "SVAnnotate" => Some(run_sv_annotate),
        "ReferenceBlockConcordance" => Some(run_reference_block_concordance),
        "CombineSegmentBreakpoints" => Some(run_combine_segment_breakpoints),
        "MergeMutect2CallsWithMC3" => Some(run_merge_mutect2_calls_with_mc3),
        "FilterFuncotations" => Some(run_filter_funcotations),
        "AnalyzeCovariates" => Some(run_analyze_covariates),
        "ASEReadCounter" => Some(run_ase_read_counter),
        "ConvertHeaderlessHadoopBamShardToBam" => {
            Some(run_convert_headerless_hadoop_bam_shard_to_bam)
        }
        "DownsampleByDuplicateSet" => Some(run_downsample_by_duplicate_set),
        "GetNormalArtifactData" => Some(run_get_normal_artifact_data),
        "CollectF1R2Counts" => Some(run_collect_f1r2_counts),
        "CreateSomaticPanelOfNormals" => Some(run_create_somatic_panel_of_normals),
        "SplitCRAM" => Some(run_split_cram),
        "FilterMutectCalls" => Some(run_filter_mutect_calls),
        "ComposeSTRTableFile" => Some(run_compose_str_table_file),
        "CalibrateDragstrModel" => Some(run_calibrate_dragstr_model),
        "LearnReadOrientationModel" => Some(run_learn_read_orientation_model),
        "GeneExpressionEvaluation" => Some(run_gene_expression_evaluation),
        "CRAMIssue8768Detector" => Some(run_cram_issue_8768_detector),
        "ValidateVariants" => Some(run_validate_variants),
        "LeftAlignAndTrimVariants" => Some(run_left_align_and_trim_variants),
        "GatherBQSRReports" => Some(run_gather_bqsr_reports),
        "Concordance" => Some(run_concordance),
        "DepthOfCoverage" => Some(run_depth_of_coverage),
        "CalculateGenotypePosteriors" => Some(run_calculate_genotype_posteriors),
        "DenoiseReadCounts" => Some(run_denoise_read_counts),
        "ValidateBasicSomaticShortMutations" => Some(run_validate_basic_somatic_short_mutations),
        "PrintReadCounts" => Some(run_print_read_counts),
        "CollectSVEvidence" => Some(run_collect_sv_evidence),
        "FastaAlternateReferenceMaker" => Some(run_fasta_alternate_reference_maker),
        "CompareReferences" => Some(run_compare_references),
        "CollectReadCounts" => Some(run_collect_read_counts),
        "GetSampleName" => Some(run_get_sample_name),
        "PrintDistantMates" => Some(run_print_distant_mates),
        "CompareIntervalLists" => Some(run_compare_interval_lists),
        "FixMisencodedBaseQualityReads" => Some(run_fix_misencoded_base_quality_reads),
        "AnnotateIntervals" => Some(run_annotate_intervals),
        "RevertBaseQualityScores" => Some(run_revert_base_quality_scores),
        "AddOriginalAlignmentTags" => Some(run_add_original_alignment_tags),
        "LeftAlignIndels" => Some(run_left_align_indels),
        "DumpTabixIndex" => Some(run_dump_tabix_index),
        "ReadAnonymizer" => Some(run_read_anonymizer),
        "PrintFileDiagnostics" => Some(run_print_file_diagnostics),
        "SplitReads" => Some(run_split_reads),
        "ClipReads" => Some(run_clip_reads),
        "SplitNCigarReads" => Some(run_split_n_cigar_reads),
        "MethylationTypeCaller" => Some(run_methylation_type_caller),
        "BaseRecalibrator" => Some(run_base_recalibrator),
        "GtfToBed" => Some(run_gtf_to_bed),
        "CallableLoci" => Some(run_callable_loci),
        "ShiftFasta" => Some(run_shift_fasta),
        "TransferReadTags" => Some(run_transfer_read_tags),
        "PostProcessReadsForRSEM" => Some(run_post_process_reads_for_rsem),
        "VariantsToTable" => Some(run_variants_to_table),
        "CompareBaseQualities" => Some(run_compare_base_qualities),
        "RemoveNearbyIndels" => Some(run_remove_nearby_indels),
        "UpdateVCFSequenceDictionary" => Some(run_update_vcf_sequence_dictionary),
        "SelectVariants" => Some(run_select_variants),
        _ => None,
    }
}

/// The runners, each of which needs the parsed command line rather than the raw one.
fn run_apply_bqsr(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::apply_bqsr(&parsed("ApplyBQSR", args)?)
}

fn run_annotate_intervals(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::annotate_intervals(&parsed("AnnotateIntervals", args)?)
}

fn run_revert_base_quality_scores(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::revert_base_quality_scores(&parsed("RevertBaseQualityScores", args)?)
}

fn run_add_original_alignment_tags(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::add_original_alignment_tags(&parsed("AddOriginalAlignmentTags", args)?)
}

fn run_left_align_indels(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::left_align_indels(&parsed("LeftAlignIndels", args)?)
}

fn run_dump_tabix_index(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::dump_tabix_index(&parsed("DumpTabixIndex", args)?)
}

fn run_read_anonymizer(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::read_anonymizer(&parsed("ReadAnonymizer", args)?)
}

fn run_print_file_diagnostics(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::print_file_diagnostics(&parsed("PrintFileDiagnostics", args)?)
}

fn run_split_reads(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::split_reads(&parsed("SplitReads", args)?)
}

fn run_clip_reads(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::clip_reads(&parsed("ClipReads", args)?)
}

fn run_split_n_cigar_reads(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::split_n_cigar_reads(&parsed("SplitNCigarReads", args)?)
}

fn run_methylation_type_caller(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::methylation_type_caller(&parsed("MethylationTypeCaller", args)?)
}

fn run_base_recalibrator(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::base_recalibrator(&parsed("BaseRecalibrator", args)?)
}

fn run_gtf_to_bed(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::gtf_to_bed(&parsed("GtfToBed", args)?)
}

fn run_callable_loci(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::callable_loci(&parsed("CallableLoci", args)?)
}

fn run_shift_fasta(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::shift_fasta(&parsed("ShiftFasta", args)?)
}

fn run_transfer_read_tags(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::transfer_read_tags(&parsed("TransferReadTags", args)?)
}

fn run_post_process_reads_for_rsem(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::post_process_reads_for_rsem(&parsed("PostProcessReadsForRSEM", args)?)
}

fn run_variants_to_table(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::variants_to_table(&parsed("VariantsToTable", args)?)
}

fn run_compare_base_qualities(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::compare_base_qualities(&parsed("CompareBaseQualities", args)?)
}

fn run_remove_nearby_indels(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::remove_nearby_indels(&parsed("RemoveNearbyIndels", args)?)
}

fn run_update_vcf_sequence_dictionary(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::update_vcf_sequence_dictionary(&parsed("UpdateVCFSequenceDictionary", args)?)
}

fn run_select_variants(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::select_variants(&parsed("SelectVariants", args)?)
}

fn run_compare_interval_lists(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::compare_interval_lists(&parsed("CompareIntervalLists", args)?)
}

fn run_fix_misencoded_base_quality_reads(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::fix_misencoded_base_quality_reads(&parsed("FixMisencodedBaseQualityReads", args)?)
}

fn run_get_sample_name(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::get_sample_name(&parsed("GetSampleName", args)?)
}

fn run_print_distant_mates(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::print_distant_mates(&parsed("PrintDistantMates", args)?)
}

fn run_collect_read_counts(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::collect_read_counts(&parsed("CollectReadCounts", args)?)
}

fn run_fasta_reference_maker(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::fasta_reference_maker(&parsed("FastaReferenceMaker", args)?)
}

fn run_annotate_vcf_with_expected_allele_fraction(
    args: &[String],
) -> Result<Option<String>, Thrown> {
    runners::annotate_vcf_with_expected_allele_fraction(&parsed(
        "AnnotateVcfWithExpectedAlleleFraction",
        args,
    )?)
}

fn run_calculate_mixing_fractions(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::calculate_mixing_fractions(&parsed("CalculateMixingFractions", args)?)
}

fn run_annotate_vcf_with_bam_depth(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::annotate_vcf_with_bam_depth(&parsed("AnnotateVcfWithBamDepth", args)?)
}

fn run_count_false_positives(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::count_false_positives(&parsed("CountFalsePositives", args)?)
}

fn run_evaluate_info_field_concordance(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::evaluate_info_field_concordance(&parsed("EvaluateInfoFieldConcordance", args)?)
}

fn run_call_copy_ratio_segments(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::call_copy_ratio_segments(&parsed("CallCopyRatioSegments", args)?)
}

fn run_variant_filtration(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::variant_filtration(&parsed("VariantFiltration", args)?)
}

fn run_collect_allelic_counts(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::collect_allelic_counts(&parsed("CollectAllelicCounts", args)?)
}

fn run_filter_intervals(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::filter_intervals(&parsed("FilterIntervals", args)?)
}

fn run_downsample_by_duplicate_set(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::downsample_by_duplicate_set(&parsed("DownsampleByDuplicateSet", args)?)
}

fn run_ase_read_counter(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::ase_read_counter(&parsed("ASEReadCounter", args)?)
}

fn run_convert_headerless_hadoop_bam_shard_to_bam(
    args: &[String],
) -> Result<Option<String>, Thrown> {
    runners::convert_headerless_hadoop_bam_shard_to_bam(&parsed(
        "ConvertHeaderlessHadoopBamShardToBam",
        args,
    )?)
}

fn run_filter_variant_tranches(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::filter_variant_tranches(&parsed("FilterVariantTranches", args)?)
}

fn run_apply_vqsr(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::apply_vqsr(&parsed("ApplyVQSR", args)?)
}

fn run_vcf_comparator(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::vcf_comparator(&parsed("VCFComparator", args)?)
}

fn run_sv_stratify(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::sv_stratify(&parsed("SVStratify", args)?)
}

fn run_sv_concordance(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::sv_concordance(&parsed("SVConcordance", args)?)
}

fn run_sv_cluster(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::sv_cluster(&parsed("SVCluster", args)?)
}

fn run_grouped_sv_cluster(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::grouped_sv_cluster(&parsed("GroupedSVCluster", args)?)
}

fn run_sv_annotate(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::sv_annotate(&parsed("SVAnnotate", args)?)
}

fn run_reference_block_concordance(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::reference_block_concordance(&parsed("ReferenceBlockConcordance", args)?)
}

fn run_combine_segment_breakpoints(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::combine_segment_breakpoints(&parsed("CombineSegmentBreakpoints", args)?)
}

fn run_merge_mutect2_calls_with_mc3(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::merge_mutect2_calls_with_mc3(&parsed("MergeMutect2CallsWithMC3", args)?)
}

fn run_filter_funcotations(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::filter_funcotations(&parsed("FilterFuncotations", args)?)
}

fn run_analyze_covariates(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::analyze_covariates(&parsed("AnalyzeCovariates", args)?)
}

fn run_calculate_average_combined_annotations(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::calculate_average_combined_annotations(&parsed(
        "CalculateAverageCombinedAnnotations",
        args,
    )?)
}

fn run_example_multi_feature_walker(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::example_multi_feature_walker(&parsed("ExampleMultiFeatureWalker", args)?)
}

fn run_mt_low_heteroplasmy_filter_tool(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::mt_low_heteroplasmy_filter_tool(&parsed("MTLowHeteroplasmyFilterTool", args)?)
}

fn run_numt_filter_tool(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::numt_filter_tool(&parsed("NuMTFilterTool", args)?)
}

fn run_gather_tranches(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::gather_tranches(&parsed("GatherTranches", args)?)
}

fn run_merge_annotated_regions_by_annotation(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::merge_annotated_regions_by_annotation(&parsed(
        "MergeAnnotatedRegionsByAnnotation",
        args,
    )?)
}

fn run_tag_germline_events(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::tag_germline_events(&parsed("TagGermlineEvents", args)?)
}

fn run_merge_annotated_regions(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::merge_annotated_regions(&parsed("MergeAnnotatedRegions", args)?)
}

fn run_merge_mutect_stats(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::merge_mutect_stats(&parsed("MergeMutectStats", args)?)
}

fn run_gather_pileup_summaries(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::gather_pileup_summaries(&parsed("GatherPileupSummaries", args)?)
}

fn run_gather_normal_artifact_data(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::gather_normal_artifact_data(&parsed("GatherNormalArtifactData", args)?)
}

fn run_site_depth_to_baf(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::site_depth_to_baf(&parsed("SiteDepthtoBAF", args)?)
}

fn run_print_sv_evidence(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::print_sv_evidence(&parsed("PrintSVEvidence", args)?)
}

fn run_condense_depth_evidence(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::condense_depth_evidence(&parsed("CondenseDepthEvidence", args)?)
}

fn run_path_seq_build_kmers(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::path_seq_build_kmers(&parsed("PathSeqBuildKmers", args)?)
}

fn run_path_seq_build_reference_taxonomy(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::path_seq_build_reference_taxonomy(&parsed("PathSeqBuildReferenceTaxonomy", args)?)
}

fn run_get_normal_artifact_data(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::get_normal_artifact_data(&parsed("GetNormalArtifactData", args)?)
}

fn run_collect_f1r2_counts(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::collect_f1r2_counts(&parsed("CollectF1R2Counts", args)?)
}

fn run_create_somatic_panel_of_normals(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::create_somatic_panel_of_normals(&parsed("CreateSomaticPanelOfNormals", args)?)
}

fn run_split_cram(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::split_cram(&parsed("SplitCRAM", args)?)
}

fn run_filter_mutect_calls(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::filter_mutect_calls(&parsed("FilterMutectCalls", args)?)
}

fn run_compose_str_table_file(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::compose_str_table_file(&parsed("ComposeSTRTableFile", args)?)
}

fn run_calibrate_dragstr_model(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::calibrate_dragstr_model(&parsed("CalibrateDragstrModel", args)?)
}

fn run_learn_read_orientation_model(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::learn_read_orientation_model(&parsed("LearnReadOrientationModel", args)?)
}

fn run_gene_expression_evaluation(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::gene_expression_evaluation(&parsed("GeneExpressionEvaluation", args)?)
}

fn run_cram_issue_8768_detector(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::cram_issue_8768_detector(&parsed("CRAMIssue8768Detector", args)?)
}

fn run_validate_variants(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::validate_variants(&parsed("ValidateVariants", args)?)
}

fn run_left_align_and_trim_variants(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::left_align_and_trim_variants(&parsed("LeftAlignAndTrimVariants", args)?)
}

fn run_gather_bqsr_reports(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::gather_bqsr_reports(&parsed("GatherBQSRReports", args)?)
}

fn run_concordance(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::concordance(&parsed("Concordance", args)?)
}

fn run_depth_of_coverage(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::depth_of_coverage(&parsed("DepthOfCoverage", args)?)
}

fn run_calculate_genotype_posteriors(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::calculate_genotype_posteriors(&parsed("CalculateGenotypePosteriors", args)?)
}

fn run_denoise_read_counts(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::denoise_read_counts(&parsed("DenoiseReadCounts", args)?)
}

fn run_print_read_counts(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::print_read_counts(&parsed("PrintReadCounts", args)?)
}

fn run_collect_sv_evidence(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::collect_sv_evidence(&parsed("CollectSVEvidence", args)?)
}

fn run_validate_basic_somatic_short_mutations(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::validate_basic_somatic_short_mutations(&parsed(
        "ValidateBasicSomaticShortMutations",
        args,
    )?)
}

fn run_fasta_alternate_reference_maker(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::fasta_alternate_reference_maker(&parsed("FastaAlternateReferenceMaker", args)?)
}

fn run_compare_references(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::compare_references(&parsed("CompareReferences", args)?)
}

fn run_check_reference_compatibility(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::check_reference_compatibility(&parsed("CheckReferenceCompatibility", args)?)
}

fn run_get_pileup_summaries(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::get_pileup_summaries(&parsed("GetPileupSummaries", args)?)
}

fn run_unmark_duplicates(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::unmark_duplicates(&parsed("UnmarkDuplicates", args)?)
}

fn run_calculate_contamination(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::calculate_contamination(&parsed("CalculateContamination", args)?)
}

fn run_print_reads_header(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::print_reads_header(&parsed("PrintReadsHeader", args)?)
}

fn run_check_pileup(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::check_pileup(&parsed("CheckPileup", args)?)
}

fn run_pileup(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::pileup(&parsed("Pileup", args)?)
}

fn run_preprocess_intervals(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::preprocess_intervals(&parsed("PreprocessIntervals", args)?)
}

fn run_split_intervals(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::split_intervals(&parsed("SplitIntervals", args)?)
}

fn run_count_bases_in_reference(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::count_bases_in_reference(&parsed("CountBasesInReference", args)?)
}

fn run_count_bases(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::count_bases(&parsed("CountBases", args)?)
}

fn run_flag_stat(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::flag_stat(&parsed("FlagStat", args)?)
}

fn run_gather_vcfs_cloud(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::gather_vcfs_cloud(&parsed("GatherVcfsCloud", args)?)
}

fn run_print_reads(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::print_reads(&parsed("PrintReads", args)?)
}

fn run_create_hadoop_bam_splitting_index(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::create_hadoop_bam_splitting_index(&parsed("CreateHadoopBamSplittingIndex", args)?)
}

fn run_count_variants(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::count_variants(&parsed("CountVariants", args)?)
}

fn run_count_reads(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::count_reads(&parsed("CountReads", args)?)
}

fn run_index_feature_file(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::index_feature_file(&parsed("IndexFeatureFile", args)?)
}

fn run_print_bgzf_block_information(args: &[String]) -> Result<Option<String>, Thrown> {
    runners::print_bgzf_block_information(&parsed("PrintBGZFBlockInformation", args)?)
}

/// The tool's own parser over a command line, for a caller that wants the parse and not the run.
///
/// The dispatcher reaches this through [`parsed`]; the suite that compares the expanded command
/// line reaches it directly, because what it measures is the parser's state and not a tool's work.
pub fn parse_for(tool: &str, args: &[String]) -> Result<gatk_barclay::Parser, String> {
    parsed(tool, args).map_err(|thrown| thrown.message.unwrap_or_default())
}

/// The tool's own parser, over the command line the dispatcher was handed.
fn parsed(tool: &str, args: &[String]) -> Result<gatk_barclay::Parser, Thrown> {
    let list = gatk_tools::tool_declarations::declarations(tool)
        .expect("the declarations of a tool with a runner");
    let mut parser = parser_for(tool, list);
    let borrowed: Vec<&str> = args.iter().map(String::as_str).collect();
    parser
        .parse_arguments(&borrowed)
        .map_err(|error| Thrown::command_line(error.message))?;
    Ok(parser)
}

/// The tool's parser, with the filters its descriptor was handed as defaults.
///
/// The defaults are what closes the last gap the trim had: a default counts as selected, so a
/// default filter's own arguments are accepted with no `--read-filter` on the command line, and
/// the list is the tool's rather than the descriptor's.
fn parser_for(
    tool: &str,
    list: &'static [gatk_tools::tool_declarations::Declaration],
) -> gatk_barclay::Parser {
    let parser = gatk_barclay::Parser::new(definitions::definitions(list));
    // `@PositionalArguments`, which one declared tool has: `CompareBaseQualities` takes exactly two
    // SAM files that way. It is not part of the declarations golden, because
    // `getNamedArgumentDefinitions` does not carry a positional argument at all.
    let parser = match gatk_tools::plugin_ownership::positional_arguments(tool) {
        None => parser,
        Some((minimum, maximum)) => parser.with_positional_arguments(minimum, maximum),
    };
    let parser = match gatk_tools::plugin_ownership::default_filters(tool) {
        None => parser,
        Some(defaults) => {
            parser.with_default_plugins(defaults.iter().map(|name| (*name).to_string()).collect())
        }
    };
    // `GATKReadFilterPluginDescriptor.validateAndResolvePlugins()`, which the parser runs before
    // it walks the definitions. The resolution itself is the same one the runner asks for; what
    // this places is WHEN it refuses, and the difference is visible on a command line that breaks
    // two rules at once (#1070).
    let owned = tool.to_string();
    parser.with_plugin_validation(move |parser| {
        runners::resolve_read_filters_in(parser, &owned).map(|_| ())
    })
}

/// The port's own refusal, which the reference has no equivalent of.
///
/// It is worded so that nobody mistakes it for GATK's: a tool that GATK runs and this does not is
/// a gap in the port, and a golden will never contain this string.
pub fn not_ported(name: &str) -> String {
    format!("{name} is a GATK tool that this port does not carry yet. This message is the port's own and not GATK's.")
}

/// The main usage: all three hundred and seventy-three lines of it.
///
/// `getCommandLineName()` is EMPTY for GATK, which is why the first line reads `USAGE:  <program
/// name>` with two spaces. The listing itself is [`gatk_tools::main_usage`], compared against the
/// reference's own line for line by the `main-usage` suite.
pub fn usage() -> String {
    gatk_tools::main_usage::usage(COMMAND_LINE_NAME)
}
