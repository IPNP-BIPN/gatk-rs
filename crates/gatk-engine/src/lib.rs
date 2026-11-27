//! The GATK engine, ported from `org.broadinstitute.hellbender.engine` and the utilities it rests
//! on. The target is GATK 4.7.0.0; each module's own header names the Java files it ports and the
//! version it was checked against.
//!
//! The first thing here is intervals, because they are what `-L` means and what every walker
//! restricts itself to. A tool given the wrong interval reads the wrong data and every number it
//! produces is wrong in a way no downstream comparison can attribute.
//!
//! # A map of the modules
//!
//! Every module below is listed once, under what it ports. The tools themselves live in the
//! `gatk-tools` crate; this crate is what they share.
//!
//! **Intervals and arguments** (`utils.IntervalUtils`, `GenomeLoc`, the argument collections):
//! [`interval`] resolves one `-L` string; [`interval_args`] and [`interval_arguments`] are the
//! collection that combines several with padding and set rules; [`feature_intervals`] is `-L`
//! given a feature file; [`interval_list_scatter`] is Picard's `IntervalList` scatter modes.
//!
//! **Data sources and traversal** (`org.broadinstitute.hellbender.engine`): [`reads`],
//! [`mod@reference`], [`features`] and [`variant_source`] are the read, reference and feature data
//! sources; [`context`] is the `ReferenceContext` and `ReadsContext` a walker is handed;
//! [`context_iterator`], [`locus_iterator`], [`read_states`] and [`alignment_state`] are the locus
//! traversal under `LocusWalker`; [`locus_shards`], [`activity_profile`], [`assembly_region`],
//! [`assembly_region_iterator`] and [`assembly_region_walker`] are the assembly-region traversal;
//! [`multi_feature_walker`] and [`concordance_walker`] are the two walker bases that merge several
//! inputs.
//!
//! **Downsampling and the shared random streams** (`utils.downsampling`, and the two static
//! generators in `utils.Utils`): [`downsampling`] and [`allele_biased_downsampling`];
//! [`java_random`] (`java.util.Random`, from its specification), [`well19937c`] and
//! [`permutation`] (commons-math3, Apache 2.0).
//!
//! **Reads, cigars and pileups** (`utils.read`, `utils.clipping`, `utils.pileup`): [`read`],
//! [`read_utils`] and [`read_group`] are the `GATKRead` accessors; [`cigar_builder`],
//! [`cigar_utils`], [`clipping`] and [`alignment_utils`] build, clip and left-align cigars;
//! [`read_realignment`] realigns reads through their best haplotype; [`sa_tag`] is the `SA` tag;
//! [`base_utils`] is what a base is; [`pileup`], [`read_pileup`] and [`sam_pileup`] are pileups;
//! [`fragment`] groups a read pair's evidence; [`overhang_fixing_manager`] is the overhang
//! clipping of `SplitNCigarReads`.
//!
//! **Local assembly** (`tools.walkers.haplotypecaller`): [`assembly_based_caller_utils`],
//! [`assembly_region_trimmer`], [`assembly_result_set`], [`read_threading_assembler`],
//! [`read_threading_graph`], [`seq_graph`], [`kbest_haplotype`] and [`graph_dot`] take a region's
//! reads to candidate haplotypes; [`haplotype`], [`haplotype_alignment`],
//! [`partially_determined`] and [`event_map`] are the haplotype and the variants it carries;
//! [`smith_waterman`] is the aligner; [`pileup_read_error_corrector`] is
//! `--error-correction-log-odds`; [`allele_mapping`] turns haplotypes into the alleles genotyped;
//! [`genotyper_debug`] is `--debug-genotyper-output`.
//!
//! **Likelihoods and genotyping** (`utils.pairhmm`, `utils.genotyper`, `tools.walkers.genotyper`):
//! [`pair_hmm`], [`pd_pair_hmm`], [`pair_hmm_likelihood_engine`] and [`dragstr`] score reads
//! against haplotypes; [`baq`] is base alignment quality (`utils.baq`); [`allele_likelihoods`] and
//! [`allele_list`] are the likelihood matrix and its axes; [`genotype_index`],
//! [`genotype_likelihood_calculator`], [`dragen_genotypes`], [`allele_frequency_calculator`] and
//! [`subset_alleles`] are genotyping and QUAL; [`variant_context_utils`], [`tandem_repeat_units`]
//! and [`variant_getters`] are the `GATKVariantContextUtils` and `VariantContextGetters` pieces the
//! callers and annotations use.
//!
//! **Mutect2 and its filters** (`tools.walkers.mutect`): [`somatic_likelihoods`];
//! [`filtering_engine`], [`mutect_engine`], [`accumulate_data`] and [`apply_filters`] are
//! `Mutect2FilteringEngine` and its passes; [`mutect_filter_list`], [`mutect_hard_filters`],
//! [`allele_filter`], [`contamination_filter`], [`germline_filter`], [`haplotype_filter`],
//! [`normal_artifact_filter`], [`slippage_filter`] and [`strand_artifact_filter`] are the
//! filters; [`error_probabilities`], [`filtering_stats`] and [`threshold_calculator`] are what
//! the filters produce; [`somatic_clustering_model`] and [`allele_fraction_cluster`] are the
//! clustering model.
//!
//! **Contamination and validation** (`tools.walkers.contamination`, `tools.walkers.validation`):
//! [`contamination_model`], [`contamination_segmenter`], [`contamination_tables`],
//! [`pileup_summary`] and [`kernel_segmenter`] are `CalculateContamination`;
//! [`allele_pileup_counter`], [`basic_somatic_short_mutation_validator`], [`beta_binomial`] and
//! [`somatic_validation_power`] are the basic short-mutation validator.
//!
//! **Base recalibration** (`utils.recalibration`, `transformers`): [`base_recalibration_engine`],
//! [`covariates`], [`recal_datum`], [`recalibration_tables`], [`qual_quantizer`],
//! [`recal_utils`], [`recalibration_report`] and [`bqsr_transformer`]; [`gatk_report`] is the
//! report format they are written in.
//!
//! **Numerics** (`utils.MathUtils` and its neighbours, and the commons-math3 classes GATK calls):
//! [`math_utils`], [`natural_log_utils`], [`fisher_exact`], [`mann_whitney`] and [`histogram`];
//! [`brent_optimizer`], [`singular_value_decomposition`], [`persistence_optimizer`] and
//! [`copy_number_mcmc`].
//!
//! **Behaviour of libraries outside GATK that shows in the output**: [`java_format`] (`%f`),
//! [`java_hash`] (the iteration order a `HashSet` gives), [`java_regex`] (a measured
//! `java.util.regex` subset), [`jexl`] (Commons JEXL 2.1.1) and [`kryo`] (Kryo 4's `Output`).
//! Each header says how that behaviour was established.
//!
//! **PathSeq's k-mer file**: [`hopscotch`] is GATK's own `tools.spark.utils.LongHopscotchSet`,
//! whose iteration order the file written through `kryo` carries.
//!
//! **Tables and files**: [`tsv_table`] (`utils.tsv`), [`tranches`] (VQSR's tranches file) and
//! [`sample_db`] (pedigrees).
//!
//! **[`gpu`]**: a second implementation of paths that are already byte-identical, shipped only
//! while it produces the same bytes.

pub mod accumulate_data;
pub mod activity_profile;
pub mod alignment_state;
pub mod alignment_utils;
pub mod allele_biased_downsampling;
pub mod allele_filter;
pub mod allele_fraction_cluster;
pub mod allele_frequency_calculator;
pub mod allele_likelihoods;
pub mod allele_list;
pub mod allele_mapping;
pub mod allele_pileup_counter;
pub mod apply_filters;
pub mod assembly_based_caller_utils;
pub mod assembly_region;
pub mod assembly_region_iterator;
pub mod assembly_region_trimmer;
pub mod assembly_region_walker;
pub mod assembly_result_set;
pub mod baq;
pub mod base_recalibration_engine;
pub mod base_utils;
pub mod basic_somatic_short_mutation_validator;
pub mod beta_binomial;
pub mod bqsr_transformer;
pub mod brent_optimizer;
pub mod cigar_builder;
pub mod cigar_utils;
pub mod clipping;
pub mod concordance_walker;
pub mod contamination_filter;
pub mod contamination_model;
pub mod contamination_segmenter;
pub mod contamination_tables;
pub mod context;
pub mod context_iterator;
pub mod copy_number_mcmc;
pub mod covariates;
pub mod downsampling;
pub mod dragen_genotypes;
pub mod dragstr;
pub mod error_probabilities;
pub mod event_map;
pub mod feature_intervals;
pub mod features;
pub mod filtering_engine;
pub mod filtering_stats;
pub mod fisher_exact;
pub mod fragment;
pub mod gatk_report;
pub mod genotype_index;
pub mod genotype_likelihood_calculator;
pub mod genotyper_debug;
pub mod germline_filter;
pub mod gpu;
pub mod graph_dot;
pub mod haplotype;
pub mod haplotype_alignment;
pub mod haplotype_filter;
pub mod histogram;
pub mod hopscotch;
pub mod interval;
pub mod interval_args;
pub mod interval_arguments;
pub mod interval_list_scatter;
pub mod java_format;
pub mod java_hash;
pub mod java_random;
pub mod java_regex;
pub mod jexl;
pub mod kbest_haplotype;
pub mod kernel_segmenter;
pub mod kryo;
pub mod locus_iterator;
pub mod locus_shards;
pub mod mann_whitney;
pub mod math_utils;
pub mod multi_feature_walker;
pub mod mutect_engine;
pub mod mutect_filter_list;
pub mod mutect_hard_filters;
pub mod natural_log_utils;
pub mod normal_artifact_filter;
pub mod overhang_fixing_manager;
pub mod pair_hmm;
pub mod pair_hmm_likelihood_engine;
pub mod partially_determined;
pub mod pd_pair_hmm;
pub mod permutation;
pub mod persistence_optimizer;
pub mod pileup;
pub mod pileup_read_error_corrector;
pub mod pileup_summary;
pub mod qual_quantizer;
pub mod read;
pub mod read_group;
pub mod read_pileup;
pub mod read_realignment;
pub mod read_states;
pub mod read_threading_assembler;
pub mod read_threading_graph;
pub mod read_utils;
pub mod reads;
pub mod recal_datum;
pub mod recal_utils;
pub mod recalibration_report;
pub mod recalibration_tables;
pub mod reference;
pub mod sa_tag;
pub mod sam_pileup;
pub mod sample_db;
pub mod seq_graph;
pub mod singular_value_decomposition;
pub mod slippage_filter;
pub mod smith_waterman;
pub mod somatic_clustering_model;
pub mod somatic_likelihoods;
pub mod somatic_validation_power;
pub mod strand_artifact_filter;
pub mod subset_alleles;
pub mod tandem_repeat_units;
pub mod threshold_calculator;
pub mod tranches;
pub mod tsv_table;
pub mod variant_context_utils;
pub mod variant_getters;
pub mod variant_source;
pub mod well19937c;
