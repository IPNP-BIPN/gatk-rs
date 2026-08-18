//! `FilterMutectCalls`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls`
//! (GATK 4.6.2.0).
//!
//! The tool the whole Mutect filtering stack exists for: unfiltered calls and their `.stats` table
//! in, a filtered VCF and a filtering-stats file out. Every computation under it is ported and
//! oracle-backed; this is the driver that sequences them.
//!
//! # Four passes, and the first record is filtered by a model that has seen the last
//!
//! ```java
//! private static final int NUMBER_OF_LEARNING_PASSES = 2;
//! protected int numberOfPasses() { return NUMBER_OF_LEARNING_PASSES + 2; }
//! ```
//!
//! Passes 0, 1 and 2 accumulate; 0 and 1 learn the parameters afterwards; 2 learns the **threshold
//! alone**, so that "the final threshold used corresponds exactly to the filters"; 3 applies and
//! writes. A one-record input and a many-record input therefore filter the same record differently,
//! and the golden runs exactly that pair.
//!
//! # The header is rewritten, not appended to
//!
//! Mutect2's `filtering_status` line is dropped and replaced under the same key, every `##FILTER`
//! line in `MUTECT_FILTER_NAMES` is added whether or not its filter runs, and `AS_FilterStatus` and
//! `STRQ` arrive as `##INFO`.

use gatk_engine::accumulate_data::{accumulate_data, action_after_pass, AfterPassAction};
use gatk_engine::apply_filters::{apply_filters, AppliedRecord};
use gatk_engine::error_probabilities::{by_type, combined, kept, ErrorType};
use gatk_engine::filtering_engine::{error_probabilities_by_filter, EngineArguments, Record};
use gatk_engine::filtering_stats::{write_summary, FilterStats};
use gatk_engine::haplotype_filter::{
    FilterAnswer as HaplotypeAnswer, FilterIdentity, FilteredHaplotypeFilter,
};
use gatk_engine::mutect_filter_list::{
    filter_line, FILTERED_FILTERING_STATUS, FILTERING_STATUS_VCF_KEY, INFO_LINES,
    MUTECT_FILTER_NAMES,
};
use gatk_engine::somatic_clustering_model::{PriorArguments, SomaticClusteringModel};
use gatk_engine::strand_artifact_filter::{self as strand, EStep, LearnedParameters};
use gatk_engine::threshold_calculator::{Strategy, ThresholdCalculator};

/// `FilterMutectCalls.NUMBER_OF_LEARNING_PASSES`.
pub const NUMBER_OF_LEARNING_PASSES: i32 = 2;

/// `M2FiltersArgumentCollection`'s threshold defaults.
pub const DEFAULT_INITIAL_POSTERIOR_THRESHOLD: f64 = 0.1;
pub const DEFAULT_MAX_FALSE_DISCOVERY_RATE: f64 = 0.05;
pub const DEFAULT_F_SCORE_BETA: f64 = 1.0;

/// The tool's arguments, as far as this driver reads them.
#[derive(Debug, Clone, PartialEq)]
pub struct ToolArguments {
    pub engine: EngineArguments,
    pub strategy: Strategy,
    pub initial_posterior_threshold: f64,
    pub max_false_discovery_rate: f64,
    pub f_score_beta: f64,
    /// The `callable` value of the Mutect stats table, `None` when it says fewer than one.
    pub callable_sites: Option<f64>,
}

impl Default for ToolArguments {
    fn default() -> Self {
        Self {
            engine: EngineArguments::default(),
            strategy: Strategy::OptimalFScore,
            initial_posterior_threshold: DEFAULT_INITIAL_POSTERIOR_THRESHOLD,
            max_false_discovery_rate: DEFAULT_MAX_FALSE_DISCOVERY_RATE,
            f_score_beta: DEFAULT_F_SCORE_BETA,
            callable_sites: None,
        }
    }
}

/// What `onTraversalStart` refuses when the Mutect stats table is not beside the VCF.
///
/// A missing table is a `UserException` naming the file, not a silent default: the empirical priors
/// would otherwise be learned from nothing without anyone being told.
#[derive(Debug, Clone, PartialEq)]
pub struct MissingStatsTable {
    pub path: String,
}

impl MissingStatsTable {
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException$CouldNotReadInputFile"
    }

    pub fn message(&self) -> String {
        format!(
            "Mutect stats table {} not found.  When Mutect2 outputs a file calls.vcf it also \
             creates a calls.vcf.stats file.  Perhaps this file was not moved along with the vcf, \
             or perhaps it was not delocalized from a virtual machine while running in the cloud.",
            self.path
        )
    }
}

/// What the run produced.
#[derive(Debug, Clone, PartialEq)]
pub struct Output {
    /// One applied result per record, in input order.
    pub records: Vec<AppliedRecord>,
    /// The filtering-stats file, as written.
    pub filtering_stats: String,
}

/// `onTraversalStart`'s header lines: every `##FILTER`, the two `##INFO`, and the rewritten
/// `##filtering_status`.
///
/// The lines carry no leading `##`: `VCFHeaderLine.toString` is the key and the value, and the
/// writer adds the hashes.
pub fn output_header_lines(input_filter_lines: &[String]) -> Vec<String> {
    // `headerLines` is a `Set`, and `VCFHeader` writes its metadata in sorted order, so the output's
    // order is the strings' and not the order they were added in.
    let mut lines: Vec<String> = Vec::new();
    lines.push(format!(
        "##{FILTERING_STATUS_VCF_KEY}={FILTERED_FILTERING_STATUS}"
    ));
    for name in MUTECT_FILTER_NAMES {
        lines.push(format!(
            "##{}",
            filter_line(name).unwrap_or_else(|| panic!("no line for {name}"))
        ));
    }
    for (_, line) in INFO_LINES {
        lines.push(format!("##{line}"));
    }
    // The input's own `##FILTER` lines survive: only the `filtering_status` line is dropped.
    for line in input_filter_lines {
        if !lines.contains(line) {
            lines.push(line.clone());
        }
    }
    lines.sort();
    lines
}

/// One record's answers, kept so the accumulation and the application read the same thing.
struct Answers {
    engine: Vec<gatk_engine::filtering_engine::EngineAnswer>,
}

impl Answers {
    fn for_haplotype(&self) -> Vec<HaplotypeAnswer> {
        self.engine
            .iter()
            .map(|answer| HaplotypeAnswer {
                name: answer.name.to_string(),
                identity: match answer.class {
                    "GermlineFilter" => FilterIdentity::Germline,
                    "NormalArtifactFilter" => FilterIdentity::NormalArtifact,
                    _ => FilterIdentity::Other,
                },
                error_type: answer.error_type,
                probabilities: answer.probabilities.clone(),
            })
            .collect()
    }
}

/// `FilterMutectCalls` end to end: four passes over the records, then the two outputs.
pub fn run(records: &[Record], arguments: &ToolArguments) -> Output {
    let priors = PriorArguments {
        mitochondria: arguments.engine.list.mitochondria,
        ..PriorArguments::new()
    };
    let mut model = SomaticClusteringModel::new(priors, arguments.callable_sites);
    let mut haplotype = FilteredHaplotypeFilter::new(
        arguments
            .engine
            .max_distance_to_filtered_call_on_same_haplotype,
    );
    let mut strand_steps: Vec<EStep> = Vec::new();
    let mut strand_learned = LearnedParameters::default();
    let mut threshold = ThresholdCalculator::new(
        arguments.strategy,
        arguments.initial_posterior_threshold,
        arguments.max_false_discovery_rate,
        arguments.f_score_beta,
    );

    let mut applied = Vec::new();
    let mut statistics = Statistics::default();

    for pass in 0..(NUMBER_OF_LEARNING_PASSES + 2) {
        statistics = Statistics::default();
        applied = Vec::new();
        for record in records {
            let answers = Answers {
                engine: error_probabilities_by_filter(
                    &mut model,
                    &haplotype,
                    &strand_learned,
                    &arguments.engine,
                    record,
                )
                .expect("the golden's records are answerable"),
            };
            if pass <= NUMBER_OF_LEARNING_PASSES {
                accumulate(
                    &answers,
                    record,
                    &mut model,
                    &mut haplotype,
                    &mut strand_steps,
                    &mut threshold,
                    arguments,
                );
            } else {
                let result = apply(&answers, record, threshold.threshold());
                statistics.record(&answers, &result, threshold.threshold());
                applied.push(result);
            }
        }
        match action_after_pass(pass) {
            Some(AfterPassAction::LearnParameters) => {
                // `filters.forEach(learnParametersAndClearAccumulatedData)`, then the model, then
                // the threshold, then the statistics.
                haplotype.learn_parameters_and_clear_accumulated_data();
                strand_learned = strand::learn_parameters(&strand_steps);
                strand_steps.clear();
                model
                    .learn_and_clear_accumulated_data()
                    .expect("the golden's data is learnable");
                threshold
                    .relearn()
                    .expect("the golden's posteriors are learnable");
            }
            // The filters are frozen here on purpose: only the threshold moves.
            Some(AfterPassAction::LearnThresholdOnly) => {
                threshold
                    .relearn()
                    .expect("the golden's posteriors are learnable");
            }
            _ => {}
        }
    }

    Output {
        filtering_stats: statistics.write(&applied, &threshold, &model),
        records: applied,
    }
}

/// `accumulateData` plus the two filters that keep state of their own.
fn accumulate(
    answers: &Answers,
    record: &Record,
    model: &mut SomaticClusteringModel,
    haplotype: &mut FilteredHaplotypeFilter,
    strand_steps: &mut Vec<EStep>,
    threshold: &mut ThresholdCalculator,
    arguments: &ToolArguments,
) {
    let all: Vec<_> = answers
        .engine
        .iter()
        .map(|answer| answer.as_error_probability())
        .collect();
    let applied_answers = kept(&all);
    let artifact =
        by_type(&applied_answers, ErrorType::Artifact).expect("the lists are rectangular");
    let non_somatic =
        by_type(&applied_answers, ErrorType::NonSomatic).expect("the lists are rectangular");
    let combined_probabilities = combined(&all).expect("the lists are rectangular");

    // `filters.forEach(f -> f.accumulateDataForLearning(vc, errorProbabilities, this))`.
    haplotype.accumulate_data_for_learning(
        &answers.for_haplotype(),
        &record.phased_genotypes(),
        record.start,
    );
    if let Some(table) = &record.strand_bias_table {
        if let Ok(parsed) = strand::parse_strand_bias_table(table) {
            let sizes: Vec<i32> = record
                .alternates
                .iter()
                .map(|alternate| strand::indel_size(record.reference_length, alternate.allele))
                .collect();
            if let Ok(steps) = strand::calculate_artifact_probabilities(
                &parsed,
                &sizes,
                LearnedParameters::default().strand_artifact_prior,
                LearnedParameters::default().alpha_strand,
                LearnedParameters::default().beta_strand,
            ) {
                strand_steps.extend(steps);
            }
        }
    }

    let mut depths = record
        .tumour_allele_depths()
        .expect("the golden's records carry tumour depths");
    let mut accumulated = Vec::new();
    let _ = accumulate_data(
        model,
        &mut accumulated,
        &record.alternates,
        &mut depths,
        record.tumor_log_odds().as_deref(),
        &artifact,
        &non_somatic,
        &combined_probabilities,
        record.reference_length,
    );
    threshold.add(&accumulated);
    let _ = arguments;
}

/// `applyFiltersAndAccumulateOutputStats`, without the statistics.
fn apply(answers: &Answers, record: &Record, threshold: f64) -> AppliedRecord {
    let applied: Vec<_> = answers
        .engine
        .iter()
        .map(|answer| {
            let annotation = annotation_for(answer.class).map(str::to_string);
            // `applyFilters` checks the required annotations a SECOND time, for the annotation
            // alone: a per-site filter whose annotations are missing has already been zeroed by its
            // base class, and this decides whether its posterior is written beside the zero.
            answer.as_applied(
                annotation,
                required_annotations_present(answer.class, record),
            )
        })
        .collect();
    let alleles: Vec<_> = record.alternates.iter().map(|a| a.allele).collect();
    apply_filters(&applied, &alleles, threshold).expect("the golden's records are applicable")
}

/// `requiredInfoAnnotations().stream().allMatch(vc::hasAttribute)`, for the filters that annotate.
fn required_annotations_present(class: &str, record: &Record) -> bool {
    match class {
        "TumorEvidenceFilter" => record.tumor_log_10_odds.is_some(),
        "ContaminationFilter" => record.population_af.is_some(),
        "GermlineFilter" => record.tumor_log_10_odds.is_some() && record.population_af.is_some(),
        "StrandArtifactFilter" => record.strand_bias_table.is_some(),
        "PolymeraseSlippageFilter" => {
            record.repeats_per_allele.is_some() && record.repeat_unit.is_some()
        }
        _ => true,
    }
}

/// `phredScaledPosteriorAnnotationName`, by class.
fn annotation_for(class: &str) -> Option<&'static str> {
    match class {
        "TumorEvidenceFilter" => Some("SEQQ"),
        "ContaminationFilter" => Some("CONTQ"),
        "GermlineFilter" => Some("GERMQ"),
        "StrandArtifactFilter" => Some("STRANDQ"),
        "PolymeraseSlippageFilter" => Some("STRQ"),
        _ => None,
    }
}

/// `FilteringOutputStats`.
#[derive(Debug, Default)]
struct Statistics {
    pass: f64,
    true_positives: f64,
    false_positives: f64,
    false_negatives: f64,
    /// Per filter name, in the order the filters answered.
    per_filter: Vec<(String, f64, f64)>,
}

impl Statistics {
    fn record(&mut self, answers: &Answers, _applied: &AppliedRecord, threshold: f64) {
        let all: Vec<_> = answers
            .engine
            .iter()
            .map(|answer| answer.as_error_probability())
            .collect();
        let per_allele = combined(&all).expect("the lists are rectangular");
        let threshold = threshold - gatk_engine::apply_filters::EPSILON;
        let is_filtered: Vec<bool> = per_allele.iter().map(|p| *p > threshold).collect();
        for probability in &per_allele {
            if *probability > threshold {
                self.false_negatives += 1.0 - probability;
            } else {
                self.pass += 1.0;
                self.false_positives += probability;
                self.true_positives += 1.0 - probability;
            }
        }
        for (index, combined_probability) in per_allele.iter().enumerate() {
            for answer in &answers.engine {
                let allele_probability = answer.probabilities[index];
                let entry = self.entry(answer.name);
                if allele_probability > gatk_engine::apply_filters::EPSILON
                    && allele_probability > threshold - gatk_engine::apply_filters::EPSILON
                {
                    entry.2 += 1.0 - combined_probability;
                } else if !is_filtered[index] {
                    entry.1 += allele_probability;
                }
            }
        }
    }

    fn entry(&mut self, name: &str) -> &mut (String, f64, f64) {
        if let Some(index) = self.per_filter.iter().position(|(key, _, _)| key == name) {
            return &mut self.per_filter[index];
        }
        self.per_filter.push((name.to_string(), 0.0, 0.0));
        self.per_filter.last_mut().expect("just pushed")
    }

    /// `writeFilteringStats`.
    fn write(
        &self,
        _applied: &[AppliedRecord],
        threshold: &ThresholdCalculator,
        model: &SomaticClusteringModel,
    ) -> String {
        let total_true_variants = self.true_positives + self.false_negatives;
        let stats: Vec<FilterStats> = self
            .per_filter
            .iter()
            .filter(|(_, false_positives, false_negatives)| {
                *false_positives > 0.0 || *false_negatives > 0.0
            })
            .map(|(name, false_positives, false_negatives)| FilterStats {
                filter_name: name.clone(),
                false_positive_count: *false_positives,
                false_discovery_rate: false_positives / self.pass,
                false_negative_count: *false_negatives,
                false_negative_rate: false_negatives / total_true_variants,
            })
            .collect();
        write_summary(
            &stats,
            &model.clustering_metadata(),
            threshold.threshold(),
            self.pass,
            self.true_positives,
            self.false_positives,
            self.false_negatives,
        )
    }
}
