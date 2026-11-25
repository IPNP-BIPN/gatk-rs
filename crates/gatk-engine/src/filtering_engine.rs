//! The filtering engine assembled: the eighteen filters over one record, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine`
//! (GATK 4.7.0.0; `Mutect2FilteringEngine.java` and `ErrorProbabilities.java` are unchanged since
//! 4.6.2.0).
//!
//! Every filter body is ported and measured on its own. This is the thing that holds them: it
//! builds the list a mode calls for, asks all of them about one record in construction order, and
//! hands the answers to [`crate::error_probabilities`] and [`crate::apply_filters`].
//!
//! # Three different ways of saying nothing
//!
//! A per-allele filter whose required annotation is missing answers an **empty list**, which
//! `ErrorProbabilities` drops. A per-site filter in the same position answers **`0.0` per allele**,
//! which it counts. And [`crate::contamination_filter`] can answer **`NaN`**, which it combines.
//! The engine has to keep the three apart, so the private `FilterOutcome` does.
//!
//! # Seventeen of the eighteen answer a fully annotated record
//!
//! `StrictStrandBiasFilter` is switched off by its default of zero reads required on each strand
//! and answers an empty list, which is dropped. A record carrying every annotation the filters read
//! therefore has seventeen entries, not eighteen.
//!
//! # The mode reaches the arithmetic, not only the list
//!
//! Mitochondrial mode drops six filters *and* changes the somatic priors, so the tumour-evidence
//! probability of one record is `2.970295405296338E-14` by default and `9.383000947252313E-18` in
//! mitochondrial mode.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing" (`Some(x)` or `None`). Every `Option` field of [`Record`]
//!   is an INFO annotation, and `None` is `!vc.hasAttribute(key)`, which is what each filter's
//!   `requiredInfoAnnotations` check reads. `.as_deref()` turns an owned `Option<Vec<T>>` into an
//!   `Option<&[T]>` borrowed view without copying.
//! - `Result<T, E>` is "a value or an error", this codebase's Java exception. `.map_err(f)`
//!   rewrites the error, and the `?` after it returns that error to the caller at once.
//! - `match` is Java's `switch` that must cover every case; the dispatcher's last arm,
//!   `other => panic!(...)`, is the case no built filter reaches.
//! - `|x| expr` is a lambda, Java's `x -> expr`.

use crate::accumulate_data::AccumulationAllele;
use crate::allele_filter::{alt_data_by_allele, sum_ads_over_samples, GenotypeData};
use crate::apply_filters::{FilterAnswer as AppliedAnswer, FilterKind};
use crate::contamination_filter::contamination_error_probabilities;
use crate::error_probabilities::{ErrorType, FilterAnswer};
use crate::germline_filter::germline_error_probabilities;
use crate::haplotype_filter::{FilteredHaplotypeFilter, PhasedGenotype};
use crate::mutect_engine::round_finite_precision_errors;
use crate::mutect_filter_list::{
    build_filters_list, FilterArguments, FILTERS, READ_ORIENTATION_FILTER,
};
use crate::mutect_hard_filters as hard;
use crate::normal_artifact_filter::normal_artifact_error_probabilities;
use crate::slippage_filter::slippage_error_probabilities;
use crate::somatic_clustering_model::{AlternateAllele, SomaticClusteringModel};
use crate::strand_artifact_filter as strand;

/// The defaults of `M2FiltersArgumentCollection` the eighteen filters read.
///
/// `list` decides which filters exist; each other field is the `MTFAC` value the reference passes
/// to a filter's constructor in `buildFiltersList`, named after it in snake case. [`Default`]
/// holds the reference's declared defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineArguments {
    pub list: FilterArguments,
    pub min_median_base_quality: i32,
    /// `getMinMedianMappingQuality()`, resolved by the caller because the getter memoises.
    pub min_median_mapping_quality: i32,
    pub long_indel_length: i32,
    pub unique_alt_read_count: i32,
    pub contamination_estimate: f64,
    pub min_reads_on_each_strand: i32,
    pub min_median_read_position: i32,
    pub min_af: f64,
    pub normal_pileup_p_value_threshold: f64,
    pub n_ratio: f64,
    pub max_events_in_region: i32,
    pub max_events_in_haplotype: i32,
    pub num_alt_alleles_threshold: usize,
    pub max_median_fragment_length_difference: i32,
    pub min_slippage_length: i32,
    pub slippage_rate: f64,
    pub max_distance_to_filtered_call_on_same_haplotype: i32,
}

impl Default for EngineArguments {
    fn default() -> Self {
        Self {
            list: FilterArguments::default(),
            min_median_base_quality: 20,
            min_median_mapping_quality: 30,
            long_indel_length: 5,
            unique_alt_read_count: hard::DEFAULT_MIN_UNIQUE_ALT_READS,
            contamination_estimate: crate::contamination_filter::DEFAULT_CONTAMINATION,
            min_reads_on_each_strand: 0,
            min_median_read_position: 1,
            min_af: hard::DEFAULT_MIN_AF,
            normal_pileup_p_value_threshold:
                crate::normal_artifact_filter::DEFAULT_NORMAL_P_VALUE_THRESHOLD,
            n_ratio: hard::DEFAULT_MAX_N_RATIO,
            max_events_in_region: 3,
            max_events_in_haplotype: 2,
            num_alt_alleles_threshold: 1,
            max_median_fragment_length_difference: 10000,
            min_slippage_length: crate::slippage_filter::DEFAULT_MIN_SLIPPAGE_LENGTH,
            slippage_rate: crate::slippage_filter::DEFAULT_SLIPPAGE_RATE,
            max_distance_to_filtered_call_on_same_haplotype:
                crate::haplotype_filter::DEFAULT_MAX_INTRA_HAPLOTYPE_DISTANCE,
        }
    }
}

/// One record, as much of it as the eighteen filters read.
///
/// The reference hands every filter the whole `VariantContext`; here the caller extracts the
/// annotations once. An `Option` field is an INFO annotation, `None` when the record lacks it.
/// The per-genotype vectors (`genotypes`, `allele_fractions`, `phasing`, and the three resolved
/// lookups) are indexed alike, one entry per sample in the record's genotype order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Record {
    pub start: i32,
    pub reference_length: i32,
    pub alternates: Vec<AccumulationAllele>,
    pub genotypes: Vec<GenotypeData<i32>>,
    /// `AF` per genotype, empty where the genotype has none.
    pub allele_fractions: Vec<Vec<f64>>,
    /// `PGT`/`PID` per genotype.
    pub phasing: Vec<(Option<String>, Option<String>)>,
    pub tumor_log_10_odds: Option<Vec<f64>>,
    pub normal_artifact_log_10_odds: Option<Vec<f64>>,
    pub normal_log_10_odds: Option<Vec<f64>>,
    pub population_af: Option<Vec<f64>>,
    pub median_base_quality: Option<Vec<i32>>,
    pub median_mapping_quality: Option<Vec<i32>>,
    pub median_fragment_length: Option<Vec<i32>>,
    pub median_read_position: Option<Vec<i32>>,
    pub unique_alt_read_count: Option<Vec<i32>>,
    pub strand_bias_table: Option<String>,
    pub n_count: Option<i32>,
    pub event_count_in_region: Option<i32>,
    pub event_count_in_haplotype: Option<i32>,
    pub repeats_per_allele: Option<Vec<String>>,
    pub repeat_unit: Option<String>,
    pub in_panel_of_normals: bool,
    /// `getIndelLengths()`, `None` for a record that is not an indel.
    pub indel_lengths: Option<Vec<i32>>,
    /// `ContaminationFilter`'s `contaminationBySample.getOrDefault(sample, estimate)` per genotype,
    /// empty when no table names any sample, which is the estimate for all.
    pub contaminations: Vec<f64>,
    /// `GermlineFilter`'s minor allele fraction of the first segment over the record per genotype,
    /// `0.5` without one; empty for all `0.5`.
    pub minor_allele_fractions: Vec<f64>,
    /// `ReadOrientationFilter.artifactProbability` per genotype, `None` for a record that is
    /// neither a SNP nor an MNP, which the filter answers with zero.
    pub read_orientation_posteriors: Option<Vec<f64>>,
}

impl Record {
    /// The alternate alleles as the clustering model reads them.
    pub fn alleles(&self) -> Vec<AlternateAllele> {
        self.alternates.iter().map(|a| a.allele).collect()
    }

    /// `vc.getNAlleles() - 1`.
    pub fn alternate_count(&self) -> usize {
        self.alternates.len()
    }

    /// `sumADsOverSamples(vc, true, false)`, which the clustering model is handed and writes back
    /// through.
    pub fn tumour_allele_depths(
        &self,
    ) -> Result<Vec<i32>, crate::allele_filter::AlleleDepthTooShort> {
        sum_ads_over_samples(self.alternates.len() + 1, &self.genotypes, true, false)
    }

    /// `Mutect2FilteringEngine.getTumorLogOdds(vc)`:
    ///
    /// ```java
    /// return tumorLog10Odds == null ? null
    ///         : MathUtils.applyToArrayInPlace(tumorLog10Odds, MathUtils::log10ToLog);
    /// ```
    ///
    /// (One line in the reference, wrapped here.)
    pub fn tumor_log_odds(&self) -> Option<Vec<f64>> {
        self.tumor_log_10_odds.as_ref().map(|odds| {
            odds.iter()
                .map(|value| crate::allele_likelihoods::log10_to_log(*value))
                .collect()
        })
    }

    /// Each genotype with its `AF`, `PGT` and `PID`, the shape `FilteredHaplotypeFilter` reads.
    pub fn phased_genotypes(&self) -> Vec<PhasedGenotype> {
        // `.enumerate()` pairs each element with its index, as `IntStream.range` would.
        self.genotypes
            .iter()
            .enumerate()
            .map(|(index, genotype)| PhasedGenotype {
                tumor: genotype.tumor,
                allele_fractions: self.allele_fractions[index].clone(),
                phasing_gt: self.phasing[index].0.clone(),
                phasing_id: self.phasing[index].1.clone(),
            })
            .collect()
    }
}

/// What one filter answered, kept apart from what it means.
///
/// `Nothing` is the empty list a per-allele filter answers when its annotations are missing, which
/// `ErrorProbabilities` drops; a per-site filter in the same position answers `PerSite(0.0)`, which
/// it counts. A `PerAllele` list may hold `NaN` (`ContaminationFilter`'s uncovered allele), which
/// is an answer and is kept.
#[derive(Debug, Clone, PartialEq)]
enum FilterOutcome {
    Nothing,
    PerAllele(Vec<f64>),
    PerSite(f64),
}

/// One filter's answer, with what the engine needs to combine and to apply it.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineAnswer {
    pub class: &'static str,
    pub name: &'static str,
    pub error_type: ErrorType,
    pub kind: FilterKind,
    /// One probability per alternate allele, empty when the filter was not evaluated.
    pub probabilities: Vec<f64>,
}

impl EngineAnswer {
    /// The shape [`crate::error_probabilities`] combines.
    pub fn as_error_probability(&self) -> FilterAnswer {
        FilterAnswer {
            error_type: self.error_type,
            probabilities: self.probabilities.clone(),
        }
    }

    /// The shape [`crate::apply_filters`] applies.
    pub fn as_applied(&self, annotation: Option<String>, required_present: bool) -> AppliedAnswer {
        AppliedAnswer {
            name: self.name.to_string(),
            kind: self.kind,
            probabilities: self.probabilities.clone(),
            annotation,
            required_annotations_present: required_present,
        }
    }
}

/// What the assembly refuses: whatever the filter under it refused.
///
/// `class` is the Java exception the reference throws at that point; `message` is the filter's own
/// error, rendered for diagnosis rather than as the JVM's message.
#[derive(Debug, Clone, PartialEq)]
pub enum EngineError {
    Filter {
        class: &'static str,
        message: String,
    },
}

/// Wraps a filter's refusal with the Java exception class the caller names. `<E: Debug>` is a
/// generic parameter: any error type that can be printed, as Java's `<E>` with a bound.
fn refuse<E: std::fmt::Debug>(class: &'static str, error: E) -> EngineError {
    EngineError::Filter {
        class,
        message: format!("{error:?}"),
    }
}

/// `new ErrorProbabilities(filters, vc, engine, referenceContext)`, as far as the per-filter
/// answers.
///
/// The answers come back in construction order, and only the ones that were evaluated: an empty
/// list is dropped here as `ErrorProbabilities` drops it (see below for when), which is why a fully
/// annotated record has seventeen entries and not eighteen.
///
/// ```java
/// // remove filters that were not applied. i.e. returned empty list
/// .entrySet().stream().filter(entry -> !entry.getValue().isEmpty())
/// .collect(toMap(Map.Entry::getKey, Map.Entry::getValue, (a, b) -> a, LinkedHashMap::new));
///
/// // if vc has symbolic alleles, remove them from each filter list
/// alleleProbabilitiesByFilter.replaceAll((k, v) ->
///         GATKVariantContextUtils.removeDataForSymbolicAltAlleles(vc, v));
/// ```
///
/// (`ErrorProbabilities.java` lines 29 to 31 and 33 to 34; the last statement is one line there.)
///
/// The reference drops the empty lists **before** the symbolic removal; this drops them after it.
/// A list that only the removal empties (a filter's answer for a record whose every alternate is
/// symbolic) is therefore kept, empty, by the reference and absent here.
///
/// The rounding is the base classes': `Mutect2AlleleFilter.errorProbabilities` rounds each
/// probability with `roundFinitePrecisionErrors`, and `Mutect2VariantFilter.errorProbabilities`
/// rounds its one probability and then copies it with `Collections.nCopies(numAltAlleles, result)`.
pub fn error_probabilities_by_filter(
    model: &mut SomaticClusteringModel,
    haplotype: &FilteredHaplotypeFilter,
    strand_learned: &strand::LearnedParameters,
    arguments: &EngineArguments,
    record: &Record,
) -> Result<Vec<EngineAnswer>, EngineError> {
    let built = build_filters_list(&arguments.list);
    let mut answers = Vec::new();
    for class in built {
        // The identity table holds the eighteen filters of `FILTERS` and, chained after them, the
        // optional `ReadOrientationFilter`; `find` answers the first whose class matches.
        let identity = FILTERS
            .iter()
            .chain([&READ_ORIENTATION_FILTER])
            .find(|filter| filter.class == class)
            .unwrap_or_else(|| panic!("no identity for {class}"));
        let outcome = answer_for(model, haplotype, strand_learned, arguments, record, class)?;
        // A per-site answer is rounded once and copied to every alternate, as `nCopies` copies it;
        // the clamp is idempotent, so a filter that already rounded its own answer is unchanged.
        let (kind, probabilities) = match outcome {
            FilterOutcome::Nothing => (FilterKind::PerAllele, Vec::new()),
            FilterOutcome::PerAllele(values) => (
                FilterKind::PerAllele,
                values
                    .into_iter()
                    .map(round_finite_precision_errors)
                    .collect(),
            ),
            FilterOutcome::PerSite(value) => (
                FilterKind::PerSite,
                vec![round_finite_precision_errors(value); record.alternate_count()],
            ),
        };
        // `alleleProbabilitiesByFilter.replaceAll((k, v) -> removeDataForSymbolicAltAlleles(vc,
        // v))`: the constructor drops every symbolic allele's entry from every filter's list,
        // whatever that list's length. `apply_filters` removes them a SECOND time, which is the
        // defect that module documents; the engine must remove exactly once here.
        let probabilities = remove_symbolic(&probabilities, record);
        if probabilities.is_empty() {
            continue;
        }
        answers.push(EngineAnswer {
            class: identity.class,
            name: identity.filter_name,
            error_type: identity.error_type,
            kind,
            probabilities,
        });
    }
    Ok(answers)
}

/// `GATKVariantContextUtils.removeDataForSymbolicAltAlleles(vc, data)`: the entries at the symbolic
/// alternates' indices, dropped.
///
/// The data holds no reference entry, so position `i` is alternate `i`; a position past the last
/// alternate is kept, as `removeItemsByIndex` keeps an index nothing asked it to remove.
fn remove_symbolic(probabilities: &[f64], record: &Record) -> Vec<f64> {
    probabilities
        .iter()
        .enumerate()
        .filter(|(index, _)| {
            !record
                .alternates
                .get(*index)
                .map(|alternate| alternate.allele.symbolic)
                .unwrap_or(false)
        })
        .map(|(_, value)| *value)
        .collect()
}

/// One filter's answer, chosen by its class.
///
/// Each arm is one filter's `errorProbabilities`, with its `requiredInfoAnnotations` check made
/// explicit: a per-allele filter whose annotation is `None` answers `Nothing`, a per-site one
/// `PerSite(0.0)`. The exception class each `refuse` names is what the reference throws there.
fn answer_for(
    model: &mut SomaticClusteringModel,
    haplotype: &FilteredHaplotypeFilter,
    strand_learned: &strand::LearnedParameters,
    arguments: &EngineArguments,
    record: &Record,
    class: &str,
) -> Result<FilterOutcome, EngineError> {
    let alternates = record.alleles();
    let alternate_count = record.alternate_count();
    let outcome = match class {
        "TumorEvidenceFilter" => {
            let odds = record.tumor_log_odds();
            let probabilities = crate::allele_filter::tumor_evidence_error_probabilities(
                model,
                odds.as_deref(),
                &record.genotypes,
                &alternates,
                record.reference_length,
            );
            if probabilities.is_empty() {
                FilterOutcome::Nothing
            } else {
                FilterOutcome::PerAllele(probabilities)
            }
        }
        // `vc.getAttributeAsIntList(MEDIAN_BASE_QUALITY_KEY, 0)`. `f64::from` widens the `int`
        // threshold to the `double` the filter's constructor takes, exactly.
        "BaseQualityFilter" => match &record.median_base_quality {
            None => FilterOutcome::Nothing,
            Some(qualities) => FilterOutcome::PerAllele(
                hard::base_quality_artifacts(
                    qualities,
                    f64::from(arguments.min_median_base_quality),
                )
                .into_iter()
                .map(hard::error_probability)
                .collect(),
            ),
        },
        "MappingQualityFilter" => match &record.median_mapping_quality {
            None => FilterOutcome::Nothing,
            Some(qualities) => FilterOutcome::PerAllele(
                hard::mapping_quality_artifacts(
                    qualities,
                    record.indel_lengths.as_deref(),
                    f64::from(arguments.min_median_mapping_quality),
                    arguments.long_indel_length,
                )
                .map_err(|error| refuse("java.lang.IndexOutOfBoundsException", error))?
                .into_iter()
                .map(hard::error_probability)
                .collect(),
            ),
        },
        "DuplicatedAltReadFilter" => match &record.unique_alt_read_count {
            None => FilterOutcome::Nothing,
            Some(counts) => FilterOutcome::PerAllele(
                hard::duplicated_alt_read_artifacts(counts, arguments.unique_alt_read_count)
                    .into_iter()
                    .map(hard::error_probability)
                    .collect(),
            ),
        },
        "StrandArtifactFilter" => match &record.strand_bias_table {
            None => FilterOutcome::Nothing,
            Some(table) => {
                let parsed = strand::parse_strand_bias_table(table)
                    .map_err(|error| refuse("java.lang.NumberFormatException", error))?;
                // `StrandArtifactFilter` removes the symbolic alleles from the TABLE itself, before
                // anything is computed, so the totals it sums exclude them. The indel sizes are NOT
                // filtered with it, which is the pairing hazard that module records.
                let table_entries: Vec<Vec<i32>> = parsed
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| {
                        *index == 0
                            || !record
                                .alternates
                                .get(index - 1)
                                .map(|alternate| alternate.allele.symbolic)
                                .unwrap_or(false)
                    })
                    .map(|(_, entry)| entry.clone())
                    .collect();
                let sizes: Vec<i32> = alternates
                    .iter()
                    .map(|alternate| strand::indel_size(record.reference_length, *alternate))
                    .collect();
                let steps = strand::calculate_artifact_probabilities(
                    &table_entries,
                    &sizes,
                    strand_learned.strand_artifact_prior,
                    strand_learned.alpha_strand,
                    strand_learned.beta_strand,
                )
                .map_err(|error| refuse("java.lang.IllegalArgumentException", error))?;
                let probabilities = strand::error_probabilities(&steps);
                if probabilities.is_empty() {
                    FilterOutcome::Nothing
                } else {
                    FilterOutcome::PerAllele(probabilities)
                }
            }
        },
        "ContaminationFilter" => match &record.population_af {
            None => FilterOutcome::Nothing,
            Some(frequencies) => {
                let contaminations = if record.contaminations.is_empty() {
                    vec![arguments.contamination_estimate; record.genotypes.len()]
                } else {
                    record.contaminations.clone()
                };
                let probabilities = contamination_error_probabilities(
                    model,
                    Some(frequencies),
                    &record.genotypes,
                    &contaminations,
                    &alternates,
                    record.reference_length,
                )
                .map_err(|error| refuse("java.lang.ArrayIndexOutOfBoundsException", error))?;
                if probabilities.is_empty() {
                    FilterOutcome::Nothing
                } else {
                    FilterOutcome::PerAllele(probabilities)
                }
            }
        },
        "StrictStrandBiasFilter" => match &record.strand_bias_table {
            None => FilterOutcome::Nothing,
            Some(table) => {
                let parsed = strand::parse_strand_bias_table(table)
                    .map_err(|error| refuse("java.lang.NumberFormatException", error))?;
                // Unlike the `StrandArtifactFilter` arm, the table is NOT stripped of symbolic
                // entries here, while the reference strips it before mapping over the alternates:
                // `sbs = GATKVariantContextUtils.removeDataForSymbolicAlleles(vc, sbs);`. The two
                // agree when every symbolic alternate comes after every other one.
                let artifacts =
                    hard::strict_strand_artifacts(&parsed, arguments.min_reads_on_each_strand);
                if artifacts.is_empty() {
                    FilterOutcome::Nothing
                } else {
                    FilterOutcome::PerAllele(
                        artifacts.into_iter().map(hard::error_probability).collect(),
                    )
                }
            }
        },
        "ReadPositionFilter" => match &record.median_read_position {
            None => FilterOutcome::Nothing,
            Some(positions) => FilterOutcome::PerAllele(
                hard::read_position_artifacts(
                    positions,
                    f64::from(arguments.min_median_read_position),
                )
                .into_iter()
                .map(hard::error_probability)
                .collect(),
            ),
        },
        // No required annotation at all, so this one always answers.
        "MinAlleleFractionFilter" => {
            // `getAltDataByAllele(vc, g -> g.hasExtendedAttribute(AF) && isTumor(g), ...)`. The
            // allele keys only need to be distinct, so their indices stand in for the alleles.
            let alleles: Vec<String> = (0..=alternate_count).map(|i| i.to_string()).collect();
            let genotypes: Vec<GenotypeData<f64>> = record
                .genotypes
                .iter()
                .enumerate()
                .map(|(index, genotype)| GenotypeData {
                    tumor: genotype.tumor,
                    allele_depths: genotype.allele_depths.clone(),
                    values: record.allele_fractions[index].clone(),
                })
                .collect();
            let gathered: Vec<Vec<f64>> =
                alt_data_by_allele(&alleles, &genotypes, |g| g.tumor && !g.values.is_empty())
                    .into_iter()
                    .map(|(_, values)| values)
                    .collect();
            FilterOutcome::PerAllele(
                hard::min_allele_fraction_artifacts(&gathered, arguments.min_af)
                    .into_iter()
                    .map(hard::error_probability)
                    .collect(),
            )
        }
        "NormalArtifactFilter" => {
            let probabilities = normal_artifact_error_probabilities(
                model,
                record.tumor_log_10_odds.as_deref(),
                record.normal_artifact_log_10_odds.as_deref(),
                record.median_base_quality.as_deref().unwrap_or(&[]),
                &record.genotypes,
                alternate_count + 1,
                arguments.normal_pileup_p_value_threshold,
            )
            .map_err(|error| refuse("java.lang.IndexOutOfBoundsException", error))?;
            // The filter already answered one copy per alternate; take the one value back so the
            // per-site path copies it. `unwrap_or(&0.0)` covers a record with no alternate.
            FilterOutcome::PerSite(*probabilities.first().unwrap_or(&0.0))
        }
        // `requiredInfoAnnotations` is `N_COUNT_KEY`, so a missing `NCount` is the 0.0 of the base
        // class before the filter is ever asked.
        "NRatioFilter" => match record.n_count {
            None => FilterOutcome::PerSite(0.0),
            Some(count) => {
                let depths =
                    sum_ads_over_samples(alternate_count + 1, &record.genotypes, true, true)
                        .map_err(|error| {
                            refuse("java.lang.ArrayIndexOutOfBoundsException", error)
                        })?;
                FilterOutcome::PerSite(hard::error_probability(hard::n_ratio_is_artifact(
                    &depths,
                    count,
                    arguments.n_ratio,
                )))
            }
        },
        "PanelOfNormalsFilter" => FilterOutcome::PerSite(hard::error_probability(
            hard::panel_of_normals_is_artifact(record.in_panel_of_normals),
        )),
        // Both `ECNT` and `ECNTH` are required. `ECNTH` is declared `VCFHeaderLineCount.A`, one
        // per alternate, and the reference takes the maximum over the whole list; this record
        // field holds one value, passed as a one-element list, so only that value reaches it.
        "ClusteredEventsFilter" => {
            match (
                record.event_count_in_haplotype,
                record.event_count_in_region,
            ) {
                (Some(haplotype_count), Some(region_count)) => {
                    FilterOutcome::PerSite(hard::error_probability(
                        hard::clustered_events_is_artifact(
                            &[haplotype_count],
                            region_count,
                            arguments.max_events_in_region,
                            arguments.max_events_in_haplotype,
                        )
                        .map_err(|error| refuse("java.util.NoSuchElementException", error))?,
                    ))
                }
                _ => FilterOutcome::PerSite(0.0),
            }
        }
        "MultiallelicFilter" => match &record.tumor_log_10_odds {
            None => FilterOutcome::PerSite(0.0),
            Some(odds) => FilterOutcome::PerSite(hard::error_probability(
                hard::multiallelic_is_artifact(Some(odds), arguments.num_alt_alleles_threshold)
                    .map_err(|error| refuse("java.lang.NullPointerException", error))?,
            )),
        },
        "FragmentLengthFilter" => match &record.median_fragment_length {
            None => FilterOutcome::PerSite(0.0),
            Some(lengths) => FilterOutcome::PerSite(hard::error_probability(
                hard::fragment_length_is_artifact(
                    lengths,
                    f64::from(arguments.max_median_fragment_length_difference),
                )
                .map_err(|error| refuse("java.lang.IndexOutOfBoundsException", error))?,
            )),
        },
        "PolymeraseSlippageFilter" => {
            let probabilities = slippage_error_probabilities(
                model,
                record.repeats_per_allele.as_deref(),
                record.repeat_unit.as_deref(),
                &record.genotypes,
                &alternates,
                record.reference_length,
                arguments.min_slippage_length,
                arguments.slippage_rate,
            )
            .map_err(|error| refuse("java.lang.NumberFormatException", error))?;
            FilterOutcome::PerSite(*probabilities.first().unwrap_or(&0.0))
        }
        // `calculateErrorProbability`: the alt-depth-weighted median of the tumour genotypes'
        // posteriors, with no required annotation.
        "ReadOrientationFilter" => match &record.read_orientation_posteriors {
            None => FilterOutcome::PerSite(0.0),
            Some(posteriors) => {
                // `(int) MathUtils.sum(ADs) - ADs[0]` per tumour genotype, paired with its
                // posterior. `.zip` walks the genotypes and the posteriors in step.
                let mut pairs: Vec<(i32, f64)> = record
                    .genotypes
                    .iter()
                    .zip(posteriors)
                    .filter(|(genotype, _)| genotype.tumor)
                    .map(|(genotype, posterior)| {
                        let depths = &genotype.allele_depths;
                        (
                            depths.iter().sum::<i32>() - depths.first().copied().unwrap_or(0),
                            *posterior,
                        )
                    })
                    .collect();
                FilterOutcome::PerSite(crate::allele_filter::weighted_median_posterior_probability(
                    &mut pairs,
                ))
            }
        },
        "FilteredHaplotypeFilter" => {
            let probabilities = haplotype
                .error_probabilities(&record.phased_genotypes(), record.start, alternate_count)
                .map_err(|error| refuse("java.util.NoSuchElementException", error))?;
            FilterOutcome::PerSite(*probabilities.first().unwrap_or(&0.0))
        }
        "GermlineFilter" => {
            // An empty vector is the caller's shorthand for 0.5 everywhere (no segment table).
            let minor = if record.minor_allele_fractions.is_empty() {
                vec![0.5; record.genotypes.len()]
            } else {
                record.minor_allele_fractions.clone()
            };
            let probabilities = germline_error_probabilities(
                model,
                record.tumor_log_10_odds.as_deref(),
                record.population_af.as_deref(),
                record.normal_log_10_odds.as_deref(),
                &record.genotypes,
                &record.allele_fractions,
                &minor,
                &alternates,
                record.reference_length,
            )
            .map_err(|error| refuse("java.lang.IllegalArgumentException", error))?;
            FilterOutcome::PerSite(*probabilities.first().unwrap_or(&0.0))
        }
        // `build_filters_list` names only the classes above, so this arm is a programming error.
        other => panic!("no filter named {other}"),
    };
    Ok(outcome)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::somatic_clustering_model::PriorArguments;

    fn real() -> AccumulationAllele {
        AccumulationAllele {
            allele: AlternateAllele {
                length: 1,
                symbolic: false,
            },
            non_ref: false,
        }
    }

    /// A record with every annotation the filters read, biallelic.
    fn annotated() -> Record {
        Record {
            start: 100,
            reference_length: 1,
            alternates: vec![real()],
            genotypes: vec![
                GenotypeData {
                    tumor: true,
                    allele_depths: vec![80, 20],
                    values: Vec::new(),
                },
                GenotypeData {
                    tumor: false,
                    allele_depths: vec![90, 1],
                    values: Vec::new(),
                },
            ],
            allele_fractions: vec![vec![0.2], Vec::new()],
            phasing: vec![
                (Some("0|1".to_string()), Some("100_A_C".to_string())),
                (None, None),
            ],
            tumor_log_10_odds: Some(vec![20.0]),
            normal_artifact_log_10_odds: Some(vec![2.0]),
            normal_log_10_odds: None,
            population_af: Some(vec![6.0]),
            median_base_quality: Some(vec![30, 30]),
            median_mapping_quality: Some(vec![60, 60]),
            median_fragment_length: Some(vec![300, 300]),
            median_read_position: Some(vec![25]),
            unique_alt_read_count: Some(vec![8]),
            strand_bias_table: Some("40,40|10,10".to_string()),
            n_count: Some(0),
            event_count_in_region: Some(1),
            event_count_in_haplotype: Some(1),
            repeats_per_allele: Some(vec!["10".to_string(), "10".to_string()]),
            repeat_unit: Some("A".to_string()),
            in_panel_of_normals: false,
            indel_lengths: None,
            ..Default::default()
        }
    }

    fn answers(record: &Record, mitochondria: bool) -> Vec<EngineAnswer> {
        let arguments = EngineArguments {
            list: FilterArguments {
                mitochondria,
                ..FilterArguments::default()
            },
            ..EngineArguments::default()
        };
        let mut model = SomaticClusteringModel::new(
            PriorArguments {
                mitochondria,
                ..PriorArguments::new()
            },
            None,
        );
        let haplotype = FilteredHaplotypeFilter::new(100);
        let strand = strand::LearnedParameters::default();
        error_probabilities_by_filter(&mut model, &haplotype, &strand, &arguments, record)
            .expect("answered")
    }

    /// Eighteen filters are built and seventeen answer: the strict-strand one is switched off by
    /// its default and an empty list is dropped rather than counted.
    #[test]
    fn seventeen_of_the_eighteen_answer_a_fully_annotated_record() {
        let answered = answers(&annotated(), false);
        assert_eq!(answered.len(), 17);
        assert!(!answered
            .iter()
            .any(|answer| answer.class == "StrictStrandBiasFilter"));
        // In mitochondrial mode twelve are built and eleven answer.
        assert_eq!(answers(&annotated(), true).len(), 11);
    }

    /// A bare record leaves every per-allele filter unevaluated and every per-site one at zero.
    #[test]
    fn a_bare_record_is_answered_only_by_the_per_site_filters() {
        let bare = Record {
            start: 100,
            reference_length: 1,
            alternates: vec![real()],
            genotypes: annotated().genotypes,
            allele_fractions: annotated().allele_fractions,
            phasing: annotated().phasing,
            ..Record::default()
        };
        let answered = answers(&bare, false);
        // `MinAlleleFractionFilter` requires no annotation, so it answers per allele; the rest that
        // answer are the per-site ones, whose missing annotations are a zero rather than a silence.
        assert!(answered
            .iter()
            .all(|answer| answer.probabilities.iter().all(|value| *value == 0.0)));
        assert!(answered
            .iter()
            .any(|answer| answer.class == "MinAlleleFractionFilter"));
    }

    /// The mode reaches the arithmetic: the same record's tumour-evidence probability differs.
    #[test]
    fn the_mode_changes_the_priors_and_not_only_the_list() {
        let record = annotated();
        let default = answers(&record, false);
        let mitochondria = answers(&record, true);
        let evidence = |answered: &[EngineAnswer]| {
            answered
                .iter()
                .find(|answer| answer.class == "TumorEvidenceFilter")
                .expect("the first filter")
                .probabilities[0]
        };
        assert_ne!(evidence(&default), evidence(&mitochondria));
    }
}
