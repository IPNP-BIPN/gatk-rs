//! `SelectVariants`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants` (GATK 4.6.2.0).
//!
//! The first piece of it: the sample selection, which `createSampleNameInclusionList` decides once
//! in `onTraversalStart` before any record is read. Four arguments feed it, `-sn`, `-se`, `-xl-sn`
//! and `-xl-se`, and what they do is not what their names suggest.
//!
//! # Matching nothing selects everything, twice over
//!
//! The accumulated set is empty in two different situations that the code cannot tell apart: no
//! sample was asked for, and every sample asked for turned out not to exist. Both reach
//! `if (samples.isEmpty()) { samples.addAll(vcfSamples); noSamplesSpecified = true; }`, so
//!
//!  * `-se zzz`, an expression matching nothing, outputs **the whole cohort**;
//!  * and `-sn ghost --allow-nonoverlapping-command-line-samples`, whose only name is removed as
//!    missing, does the same.
//!
//! Neither is an empty output and neither is a refusal. A pipeline that selects one sample by
//! pattern and gets the pattern slightly wrong silently carries every sample forward.
//!
//! # The expressions are `find()`, not `matches()`
//!
//! `Utils.filterCollectionByExpressions` compiles each expression and searches with it, so `-se s1`
//! selects a sample named `xs10`. See [`gatk_engine::java_regex`], which is where the search lives
//! and where the reference's compile failure is reproduced: an uncompilable expression comes out as
//! the regex engine's own `PatternSyntaxException`, not wrapped into a `UserException`.
//!
//! # The missing-name check looks at `-sn` only
//!
//! `samplesNotInHeader` is seeded from `sampleNames` alone, never from the expression results, so a
//! name that matches nothing is a four-paragraph refusal while an expression that matches nothing
//! is silent. The names in that message keep the order they were given, comma-separated, while
//! everything else here is sorted: both the header's sample set and the accumulated set are
//! `TreeSet`s, so `-sn tumor -sn s0` writes `s0` first.
//!
//! # Exclusion beats inclusion, and empties by two different routes
//!
//! Exclusion is applied after the empty-means-all rule, so it can empty the set again, and that
//! emptiness **is** a refusal. `noSamplesSpecified` is then `false` either because something was
//! included or because something was excluded, which is why excluding every sample with nothing
//! included reaches the same `UserException` as excluding exactly what was included.
//!
//! # `no_samples_specified` is what decides whether the record is touched at all
//!
//! `subsetGenotypesBySampleNames` returns the record unchanged when it is set, so a run that
//! selects everything keeps the record's INFO as it was, without the `AC`, `AF` and `AN` that a
//! genuinely subset record gains. An exclusion naming a sample that is not in the file clears the
//! flag without removing a column, and the record still comes out unrewritten, but by the later
//! check that the subset has as many samples and alleles as the original.
//!
//! **The genotype columns move regardless.** The writer emits calls in the header's sample order,
//! and the header's order is this sorted list, so a file whose samples are declared `tumor, s1,
//! NA12891, xs10, s0, NA12878` comes out `NA12878, NA12891, s0, s1, tumor, xs10` even when no
//! record was touched. The values follow their sample names rather than their places, which is why
//! nothing is lost, and it is also why a line-by-line comparison of input and output shows a
//! difference that is not a rewriting.

use gatk_engine::genotype_index::GenotypeIndexError;
use gatk_engine::java_hash::compare_strings;
use gatk_engine::java_regex::{self, PatternSyntaxError};
use gatk_engine::jexl::{create_expression, Value as JexlValue};
use gatk_engine::subset_alleles::{subset_alleles, AssignmentMethod, Genotype};
use gatk_engine::variant_context_utils::{trim_alleles, Allele, Variant};

/// The four sample arguments and the flag that forgives a missing name.
#[derive(Debug, Clone, Default)]
pub struct SampleArguments {
    /// `-sn`, in the order given: a `LinkedHashSet` in the reference.
    pub sample_names: Vec<String>,
    /// `-se`.
    pub sample_expressions: Vec<String>,
    /// `-xl-sn`.
    pub exclude_sample_names: Vec<String>,
    /// `-xl-se`.
    pub exclude_sample_expressions: Vec<String>,
    pub allow_nonoverlapping_command_line_samples: bool,
}

/// What the selection decided.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleSelection {
    /// The samples to keep, sorted, which is the order the output header carries.
    pub samples: Vec<String>,
    /// Whether nothing was asked for, which is what keeps a record from being rewritten.
    pub no_samples_specified: bool,
}

/// The two refusals, plus the regex engine's own.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SampleSelectionError {
    /// Names given to `-sn` that the header does not have, in the order they were given.
    SamplesNotInHeader(Vec<String>),
    /// The accumulated set emptied by exclusion.
    AllExcluded,
    /// An expression the pattern compiler refused.
    BadPattern(PatternSyntaxError),
}

impl SampleSelectionError {
    /// The message the reference's exception carries, without the prefix its class adds.
    pub fn message(&self) -> String {
        match self {
            // `String.format("%s%n%n%s%n%n%s%n%n%s", ...)`: `%n` is the platform separator, and the
            // oracle runs on Linux.
            SampleSelectionError::SamplesNotInHeader(names) => format!(
                "Samples entered on command line (through -sf or -sn) that are not present in the \
                 VCF.\n\nA list of these samples:\n\n{}\n\nTo ignore these samples, run with \
                 --allow-nonoverlapping-command-line-samples",
                names.join(",")
            ),
            SampleSelectionError::AllExcluded => {
                "All samples requested to be included were also requested to be excluded."
                    .to_string()
            }
            SampleSelectionError::BadPattern(error) => error.message(),
        }
    }

    pub fn java_class(&self) -> &'static str {
        match self {
            SampleSelectionError::SamplesNotInHeader(_) => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
            SampleSelectionError::AllExcluded => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            SampleSelectionError::BadPattern(error) => error.java_class(),
        }
    }
}

/// `createSampleNameInclusionList`, given the header's samples in whatever order the file has them.
pub fn create_sample_name_inclusion_list(
    header_samples: &[String],
    arguments: &SampleArguments,
) -> Result<SampleSelection, SampleSelectionError> {
    // `VcfUtils.getSortedSampleSet` is a TreeSet, so the order the file declares is already lost.
    // Its order is `String.compareTo`'s, which is UTF-16 code units and not UTF-8 bytes: measured,
    // the two disagree above the BMP, where a sample named with a supplementary character sorts
    // before one named `Ａ` in Java and after it under Rust's own `Ord`.
    let mut vcf_samples: Vec<String> = header_samples.to_vec();
    vcf_samples.sort_by(|left, right| compare_strings(left, right));
    vcf_samples.dedup();

    let from_expressions = filter(&vcf_samples, &arguments.sample_expressions)?;

    // Seeded from `-sn` alone: the expressions are never checked against the header.
    let not_in_header: Vec<String> = arguments
        .sample_names
        .iter()
        .filter(|name| !vcf_samples.contains(name))
        .cloned()
        .collect();

    let mut samples = SortedSet::new();
    samples.extend(&arguments.sample_names);
    samples.extend(&from_expressions);

    if !not_in_header.is_empty() {
        if arguments.allow_nonoverlapping_command_line_samples {
            samples.remove_all(&not_in_header);
        } else {
            return Err(SampleSelectionError::SamplesNotInHeader(not_in_header));
        }
    }

    // The empty set that means "all", reached either by asking for nothing or by asking only for
    // what does not exist.
    let mut no_samples_specified = false;
    if samples.is_empty() {
        samples.extend(&vcf_samples);
        no_samples_specified = true;
    }

    let excluded_from_expressions = filter(&vcf_samples, &arguments.exclude_sample_expressions)?;
    samples.remove_all(&arguments.exclude_sample_names);
    samples.remove_all(&excluded_from_expressions);
    no_samples_specified = no_samples_specified
        && arguments.exclude_sample_names.is_empty()
        && excluded_from_expressions.is_empty();

    if samples.is_empty() && !no_samples_specified {
        return Err(SampleSelectionError::AllExcluded);
    }

    Ok(SampleSelection {
        samples: samples.into_vec(),
        no_samples_specified,
    })
}

fn filter(values: &[String], expressions: &[String]) -> Result<Vec<String>, SampleSelectionError> {
    java_regex::filter_collection_by_expressions(values, expressions, false)
        .map_err(SampleSelectionError::BadPattern)
}

/// The `TreeSet<String>` the reference accumulates into: sorted, and without duplicates.
///
/// Sorted by `String.compareTo`, which is UTF-16 code-unit order. Rust's own `Ord` for `str` is
/// UTF-8 byte order, and the two disagree for every supplementary character: `"\u{1f600}"` sorts
/// before `"\u{ff21}"` in Java and after it here. Sample names are ASCII in practice, and a sort
/// order that is right in practice is exactly the kind of thing that is wrong once.
struct SortedSet {
    values: Vec<String>,
}

impl SortedSet {
    fn new() -> SortedSet {
        SortedSet { values: Vec::new() }
    }

    fn extend(&mut self, values: &[String]) {
        for value in values {
            if let Err(at) = self
                .values
                .binary_search_by(|held| compare_strings(held, value))
            {
                self.values.insert(at, value.clone());
            }
        }
    }

    fn remove_all(&mut self, values: &[String]) {
        self.values.retain(|value| !values.contains(value));
    }

    fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    fn into_vec(self) -> Vec<String> {
        self.values
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    /// `TreeSet<String>` is `String.compareTo`, which is UTF-16 code units. Measured against the
    /// reference: Java orders a supplementary character before `Ａ`, Rust's own `Ord` after it.
    #[test]
    fn the_sorted_order_is_javas_and_not_rusts() {
        let samples = names(&["\u{ff21}", "\u{1f600}", "s0"]);
        let selection =
            create_sample_name_inclusion_list(&samples, &SampleArguments::default()).expect("all");
        assert_eq!(
            selection.samples,
            names(&["s0", "\u{1f600}", "\u{ff21}"]),
            "the emoji sorts before Ａ, as its leading surrogate 0xd83d is below 0xff21"
        );

        let mut rusts_own = samples.clone();
        rusts_own.sort();
        assert_ne!(selection.samples, rusts_own);
    }
}

/// The four arguments that change what is written once the samples are known.
#[derive(Debug, Clone, Default)]
pub struct SubsetArguments {
    /// `--remove-unused-alternates`, which also forces the subsetting path for a whole-cohort run.
    pub remove_unused_alternates: bool,
    /// `--preserve-alleles`, which skips the trimming the subset would otherwise end with.
    pub preserve_alleles: bool,
    /// `--keep-original-ac`.
    pub keep_original_chr_counts: bool,
    /// `--keep-original-dp`.
    pub keep_original_depth: bool,
}

/// A record and the names its genotype columns carry, which the engine's `Variant` does not hold.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub variant: Variant,
    pub samples: Vec<String>,
}

/// `subsetGenotypesBySampleNames` followed by `addAnnotations`, which between them decide the whole
/// of the output record.
///
/// Returns the record unchanged where the reference returns `vc` itself, which is what keeps a
/// whole-cohort run from gaining annotations. The two early exits are both the reference's: the
/// first is `noSamplesSpecified && !removeUnusedAlternates`, before anything is decoded, and the
/// second is the subset having as many samples and as many alleles as the original.
pub fn subset_record(
    record: &Record,
    selection: &SampleSelection,
    arguments: &SubsetArguments,
) -> Result<Record, GenotypeIndexError> {
    if selection.no_samples_specified && !arguments.remove_unused_alternates {
        return Ok(record.clone());
    }

    // `subContextFromSamples`: the kept columns, in the record's own order.
    let kept_samples: Vec<usize> = record
        .samples
        .iter()
        .enumerate()
        .filter(|(_, name)| selection.samples.contains(name))
        .map(|(index, _)| index)
        .collect();
    let kept_genotypes: Vec<Genotype> = kept_samples
        .iter()
        .map(|index| record.variant.genotypes[*index].clone())
        .collect();

    // `rederiveAllelesFromGenotypes`: the alleles some kept genotype actually calls, in the
    // record's order, with the reference put back when no genotype called it.
    let kept_alleles: Vec<usize> = if arguments.remove_unused_alternates {
        let mut called = vec![false; record.variant.alleles.len()];
        let mut added_reference = false;
        for genotype in &kept_genotypes {
            for allele in genotype.alleles.iter().flatten() {
                added_reference = added_reference || *allele == 0;
                called[*allele] = true;
            }
        }
        if !added_reference {
            called[0] = true;
        }
        (0..record.variant.alleles.len())
            .filter(|index| called[*index])
            .collect()
    } else {
        (0..record.variant.alleles.len()).collect()
    };

    // The second early exit, which is why an exclusion naming nobody leaves the record alone.
    if kept_samples.len() == record.variant.genotypes.len()
        && kept_alleles.len() == record.variant.alleles.len()
    {
        return Ok(record.clone());
    }

    let genotypes = if kept_alleles.len() == record.variant.alleles.len() {
        kept_genotypes
    } else {
        subset_alleles(
            &kept_genotypes,
            2,
            record.variant.alleles.len(),
            &kept_alleles,
            AssignmentMethod::DoNotAssignGenotypes,
        )?
    };

    let mut subset = Variant {
        contig: record.variant.contig.clone(),
        start: record.variant.start,
        stop: record.variant.stop,
        alleles: kept_alleles
            .iter()
            .map(|index| record.variant.alleles[*index].clone())
            .collect(),
        genotypes,
        // The MLE tags describe a calling that no longer applies, so the reference strips them
        // rather than recomputing them. AC and AF are recomputed below, not stripped.
        attributes: record
            .variant
            .attributes
            .iter()
            .filter(|(key, _)| key != "MLEAC" && key != "MLEAF")
            .cloned()
            .collect(),
    };
    let names: Vec<String> = kept_samples
        .iter()
        .map(|index| record.samples[*index].clone())
        .collect();

    add_annotations(&mut subset, &record.variant, &kept_alleles, arguments);

    let variant = if arguments.preserve_alleles {
        subset
    } else {
        // The trimmer refuses nothing this path can produce: the alleles come from the record's
        // own list, and a subset of them is still a set of plain alleles.
        trim_alleles(&subset, true, true).expect("a subset of the record's own alleles")
    };
    Ok(Record {
        variant,
        samples: names,
    })
}

/// `addAnnotations`, in the reference's order: the originals first, then the recount, then the
/// depth.
fn add_annotations(
    subset: &mut Variant,
    original: &Variant,
    kept_alleles: &[usize],
    arguments: &SubsetArguments,
) {
    if arguments.keep_original_chr_counts {
        // The per-allele originals are reordered to the new allele list; `.` stands where an
        // allele the original counted no longer exists. AN is a single number and is copied.
        let reorder = |value: &str| -> String {
            if kept_alleles.len() == original.alleles.len() {
                return value.to_string();
            }
            let parts: Vec<&str> = value.split(',').collect();
            let mapped: Vec<String> = kept_alleles
                .iter()
                .skip(1)
                .map(|index| {
                    parts
                        .get(index - 1)
                        .map(|part| part.to_string())
                        .unwrap_or_else(|| ".".to_string())
                })
                .collect();
            if mapped.is_empty() {
                ".".to_string()
            } else {
                mapped.join(",")
            }
        };
        for (key, original_key) in [("AC", "AC_Orig"), ("AF", "AF_Orig")] {
            if let Some((_, value)) = original.attributes.iter().find(|(name, _)| name == key) {
                let reordered = reorder(value);
                set_attribute(subset, original_key, &reordered);
            }
        }
        if let Some((_, value)) = original.attributes.iter().find(|(name, _)| name == "AN") {
            set_attribute(subset, "AN_Orig", &value.clone());
        }
    }

    calculate_chromosome_counts(subset);

    if arguments.keep_original_depth {
        if let Some((_, value)) = original.attributes.iter().find(|(name, _)| name == "DP") {
            set_attribute(subset, "DP_Orig", &value.clone());
        }
    }

    // The depth is summed over the KEPT genotypes, skipping the filtered ones, and written only
    // where at least one of them had a DP at all: where none does, the record keeps its own.
    let mut saw_depth = false;
    let mut depth = 0;
    for genotype in &subset.genotypes {
        if is_filtered(genotype) {
            continue;
        }
        if let Some(value) = genotype.dp {
            depth += value;
            saw_depth = true;
        }
    }
    if saw_depth {
        set_attribute(subset, "DP", &depth.to_string());
    }
}

/// `VariantContextUtils.calculateChromosomeCounts(builder, false)`.
///
/// AN is the called chromosome count, AC the count per alternate and AF each of those over AN,
/// and all three skip a FILTERED genotype: a genotype carrying FT is not a called one. AC and AF
/// are removed outright where no alternate is left, which is what an ALT column of `.` means.
fn calculate_chromosome_counts(variant: &mut Variant) {
    if variant.genotypes.is_empty() {
        return;
    }
    let called: Vec<Genotype> = variant
        .genotypes
        .iter()
        .filter(|genotype| !is_filtered(genotype))
        .cloned()
        .collect();
    let allele_number: i32 = called
        .iter()
        .map(|genotype| genotype.alleles.iter().flatten().count() as i32)
        .sum();
    set_attribute(variant, "AN", &allele_number.to_string());

    if variant.alleles.len() < 2 {
        variant
            .attributes
            .retain(|(key, _)| key != "AC" && key != "AF");
        return;
    }
    let mut counts: Vec<String> = Vec::new();
    let mut frequencies: Vec<String> = Vec::new();
    for index in 1..variant.alleles.len() {
        let count: i32 = called
            .iter()
            .map(|genotype| {
                genotype
                    .alleles
                    .iter()
                    .filter(|allele| **allele == Some(index))
                    .count() as i32
            })
            .sum();
        counts.push(count.to_string());
        // `(double) count / AN`, formatted by the writer rather than here, which is why 0.5 is
        // `0.500` and 0.0 is `0.00`.
        let frequency = if allele_number == 0 {
            0.0
        } else {
            f64::from(count) / f64::from(allele_number)
        };
        frequencies.push(htsjdk_vcf::variant::format_vcf_double(frequency));
    }
    set_attribute(variant, "AC", &counts.join(","));
    set_attribute(variant, "AF", &frequencies.join(","));
}

/// `Genotype.isFiltered()`, which `PASS` and an absent field are not.
fn is_filtered(genotype: &Genotype) -> bool {
    genotype
        .attributes
        .iter()
        .any(|(key, value)| key == "FT" && value != "PASS" && value != ".")
}

fn set_attribute(variant: &mut Variant, key: &str, value: &str) {
    match variant.attributes.iter_mut().find(|(name, _)| name == key) {
        Some((_, held)) => *held = value.to_string(),
        None => variant
            .attributes
            .push((key.to_string(), value.to_string())),
    }
}

/// `VariantContext.Type`, as `determineType` decides it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

/// `determineType`: one alternate decides alone, several must agree or the record is MIXED.
///
/// The spanning deletion is not symbolic here: `*` is one base against a one-base reference, so a
/// record whose only alternate is `*` is a SNP, and `--select-type-to-include SNP` keeps it.
pub fn variant_type(alleles: &[Allele]) -> VariantType {
    if alleles.len() < 2 {
        return VariantType::NoVariation;
    }
    let reference = &alleles[0];
    let mut kind: Option<VariantType> = None;
    for alternate in &alleles[1..] {
        let this = if alternate.is_symbolic() {
            VariantType::Symbolic
        } else if alternate.len() == reference.len() {
            if reference.len() == 1 {
                VariantType::Snp
            } else {
                VariantType::Mnp
            }
        } else {
            VariantType::Indel
        };
        match kind {
            None => kind = Some(this),
            Some(seen) if seen != this => return VariantType::Mixed,
            Some(_) => {}
        }
    }
    kind.unwrap_or(VariantType::NoVariation)
}

/// `getIndelLengths`, which is `null` for anything but an INDEL or a MIXED record.
fn indel_lengths(alleles: &[Allele]) -> Option<Vec<i32>> {
    match variant_type(alleles) {
        VariantType::Indel | VariantType::Mixed => Some(
            alleles[1..]
                .iter()
                .map(|alternate| alternate.len() as i32 - alleles[0].len() as i32)
                .collect(),
        ),
        _ => None,
    }
}

/// `--restrict-alleles-to`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlleleRestriction {
    All,
    Biallelic,
    Multiallelic,
}

/// The arguments that decide which records survive.
#[derive(Debug, Clone)]
pub struct FilterArguments {
    pub types_to_include: Vec<VariantType>,
    pub types_to_exclude: Vec<VariantType>,
    pub allele_restriction: AlleleRestriction,
    pub max_indel_size: i32,
    pub min_indel_size: i32,
    pub keep_ids: Vec<String>,
    pub exclude_ids: Vec<String>,
    pub exclude_filtered: bool,
    pub exclude_non_variants: bool,
    pub max_filtered_genotypes: i32,
    pub min_filtered_genotypes: i32,
    pub max_fraction_filtered_genotypes: f64,
    pub min_fraction_filtered_genotypes: f64,
    pub max_nocall_number: i32,
    pub max_nocall_fraction: f64,
    pub select_expressions: Vec<String>,
    pub select_genotype_expressions: Vec<String>,
    pub invert_select: bool,
    pub apply_jexl_filters_first: bool,
}

impl Default for FilterArguments {
    fn default() -> FilterArguments {
        FilterArguments {
            types_to_include: Vec::new(),
            types_to_exclude: Vec::new(),
            allele_restriction: AlleleRestriction::All,
            max_indel_size: i32::MAX,
            min_indel_size: 0,
            keep_ids: Vec::new(),
            exclude_ids: Vec::new(),
            exclude_filtered: false,
            exclude_non_variants: false,
            max_filtered_genotypes: i32::MAX,
            min_filtered_genotypes: 0,
            max_fraction_filtered_genotypes: 1.0,
            min_fraction_filtered_genotypes: 0.0,
            max_nocall_number: i32::MAX,
            max_nocall_fraction: 1.0,
            select_expressions: Vec::new(),
            select_genotype_expressions: Vec::new(),
            invert_select: false,
            apply_jexl_filters_first: false,
        }
    }
}

impl FilterArguments {
    /// `considerFilteredGenotypes`: the gate runs only where an argument moved a default.
    fn consider_filtered_genotypes(&self) -> bool {
        self.max_filtered_genotypes != i32::MAX
            || self.min_filtered_genotypes != 0
            || self.max_fraction_filtered_genotypes != 1.0
            || self.min_fraction_filtered_genotypes != 0.0
    }

    fn consider_no_call_genotypes(&self) -> bool {
        self.max_nocall_number != i32::MAX || self.max_nocall_fraction != 1.0
    }
}

/// What an expression can refuse with, both of which reach the user.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SelectError {
    /// The expression compiled and then failed on a record, which a per-allele annotation does.
    Invalid { index: usize, genotype: bool },
    /// The expression did not compile, refused by the argument parser before any record was read.
    Unparseable { index: usize, text: String },
}

impl SelectError {
    pub fn message(&self) -> String {
        match self {
            SelectError::Invalid { index, genotype } => format!(
                "Invalid JEXL expression detected for {}-{index}\nSee \
                 https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions \
                 for documentation on using JEXL in GATK",
                if *genotype { "select-genotype" } else { "select" }
            ),
            // The reference's own string, missing the space after the argument name.
            SelectError::Unparseable { index, text } => format!(
                "Argument select-{index}has a bad value. Invalid expression used ({text}). Please \
                 see the JEXL docs for correct syntax."
            ),
        }
    }

    pub fn java_class(&self) -> &'static str {
        match self {
            SelectError::Invalid { .. } => "org.broadinstitute.hellbender.exceptions.UserException",
            SelectError::Unparseable { .. } => "java.lang.IllegalArgumentException",
        }
    }
}

/// As much of a record as the filtering reads, beside the [`Record`] the subsetting works on.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterRecord {
    pub id: String,
    /// The FILTER column: empty for `.` or `PASS`, the names otherwise.
    pub filters: Vec<String>,
    /// The INFO fields as an expression reads them, which is the decoded value's `toString`. A
    /// `Number=A` field decodes to a list, so its value here is `[2]` or `[1, 1]`, and nothing
    /// numeric compares against that: an expression over a per-allele annotation is a refusal
    /// rather than a false, which is what the golden holds.
    pub info: std::collections::HashMap<String, String>,
    /// Per sample, the fields a genotype expression reads.
    pub genotype_fields: Vec<std::collections::HashMap<String, String>>,
}

/// `applyFirstRoundOfFiltering` plus the two genotype-count gates, which run before the subset.
///
/// Returns whether the record survives. The Mendelian and random-fraction gates are not here:
/// one needs a pedigree and the other a random generator, and neither is measured.
pub fn keeps_before_subset(
    record: &Record,
    filter_record: &FilterRecord,
    arguments: &FilterArguments,
    selection: &SampleSelection,
) -> Result<bool, SelectError> {
    if arguments.exclude_filtered && !filter_record.filters.is_empty() {
        return Ok(false);
    }

    // `makeVariantFilter` runs before `apply` and is the same decision, so it is here.
    let selected_types = selected_types(arguments);
    if !selected_types.contains(&variant_type(&record.variant.alleles)) {
        return Ok(false);
    }
    if !arguments.keep_ids.is_empty() && !arguments.keep_ids.contains(&filter_record.id) {
        return Ok(false);
    }
    if arguments.exclude_ids.contains(&filter_record.id) {
        return Ok(false);
    }

    let biallelic = record.variant.alleles.len() == 2;
    match arguments.allele_restriction {
        AlleleRestriction::Biallelic if !biallelic => return Ok(false),
        AlleleRestriction::Multiallelic if biallelic => return Ok(false),
        _ => {}
    }

    // Both gates are about ABSOLUTE length and both reject the RECORD: one alternate out of range
    // takes the whole record with it.
    if let Some(lengths) = indel_lengths(&record.variant.alleles) {
        if lengths.iter().any(|length| {
            length.abs() > arguments.max_indel_size || length.abs() < arguments.min_indel_size
        }) {
            return Ok(false);
        }
    }

    if arguments.apply_jexl_filters_first && !passes_jexl_filters(filter_record, arguments)? {
        return Ok(false);
    }

    // The two counting gates, over the SELECTED samples.
    let kept: Vec<usize> = record
        .samples
        .iter()
        .enumerate()
        .filter(|(_, name)| selection.samples.contains(name))
        .map(|(index, _)| index)
        .collect();
    let sample_count = selection.samples.len();

    if arguments.consider_filtered_genotypes() {
        let filtered = kept
            .iter()
            .filter(|index| is_filtered(&record.variant.genotypes[**index]))
            .count() as i32;
        // `numFilteredSamples / samples.size()` is INT OVER INT in the reference, assigned to a
        // double afterwards: the fraction is 0 unless every sample is filtered. Reproduced, since
        // a port that divided properly would keep different records.
        let fraction = if sample_count == 0 {
            0.0
        } else {
            f64::from(filtered / sample_count as i32)
        };
        if filtered > arguments.max_filtered_genotypes
            || filtered < arguments.min_filtered_genotypes
            || fraction > arguments.max_fraction_filtered_genotypes
            || fraction < arguments.min_fraction_filtered_genotypes
        {
            return Ok(false);
        }
    }

    if arguments.consider_no_call_genotypes() {
        let no_calls = kept
            .iter()
            .filter(|index| {
                let genotype = &record.variant.genotypes[**index];
                !genotype.alleles.is_empty() && genotype.alleles.iter().all(Option::is_none)
            })
            .count() as i32;
        // One line below the other in the reference, and this one casts.
        let fraction = if sample_count == 0 {
            0.0
        } else {
            f64::from(no_calls) / sample_count as f64
        };
        if no_calls > arguments.max_nocall_number || fraction > arguments.max_nocall_fraction {
            return Ok(false);
        }
    }

    Ok(true)
}

/// `--exclude-non-variants` and the JEXL expressions that did not run first, both of which see the
/// record as the subset left it.
pub fn keeps_after_subset(
    subset: &Record,
    filter_record: &FilterRecord,
    arguments: &FilterArguments,
) -> Result<bool, SelectError> {
    if arguments.exclude_non_variants {
        // `isPolymorphicInSamples`: some genotype calls an alternate. `isSpanningDeletionOnly`
        // takes the other half: a record whose only alternate is `*` is not a variant either.
        let polymorphic = subset
            .variant
            .genotypes
            .iter()
            .any(|genotype| genotype.alleles.iter().flatten().any(|allele| *allele > 0));
        let spanning_only =
            subset.variant.alleles.len() == 2 && subset.variant.alleles[1].is_span_del();
        if !polymorphic || spanning_only {
            return Ok(false);
        }
    }

    if !arguments.apply_jexl_filters_first && !passes_jexl_filters(filter_record, arguments)? {
        return Ok(false);
    }
    Ok(true)
}

/// `passesJexlFilters`: the expressions are OR-ed, and `--invert-select` inverts EACH of them
/// before the or, which is not the complement of the whole.
fn passes_jexl_filters(
    record: &FilterRecord,
    arguments: &FilterArguments,
) -> Result<bool, SelectError> {
    if arguments.select_expressions.is_empty() && arguments.select_genotype_expressions.is_empty() {
        return Ok(true);
    }
    for (index, text) in arguments.select_expressions.iter().enumerate() {
        let expression = create_expression(text).map_err(|_| SelectError::Unparseable {
            index,
            text: text.clone(),
        })?;
        let matched = match expression.evaluate(&record.info) {
            Ok(JexlValue::Bool(value)) => value,
            Ok(_) => false,
            // A per-allele annotation reaches here: the comparison is not defined over a list, and
            // the reference turns the engine's complaint into a UserException naming the index.
            Err(_) => {
                return Err(SelectError::Invalid {
                    index,
                    genotype: false,
                })
            }
        };
        if matched != arguments.invert_select {
            return Ok(true);
        }
    }
    for (index, text) in arguments.select_genotype_expressions.iter().enumerate() {
        let expression = create_expression(text).map_err(|_| SelectError::Unparseable {
            index,
            text: text.clone(),
        })?;
        // Any genotype matching is enough, and the result joins the or above rather than an and.
        for fields in &record.genotype_fields {
            let matched = match expression.evaluate(fields) {
                Ok(JexlValue::Bool(value)) => value,
                Ok(_) => false,
                Err(_) => {
                    return Err(SelectError::Invalid {
                        index,
                        genotype: true,
                    })
                }
            };
            if matched != arguments.invert_select {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// `createSampleTypeInclusionList`: every type when none is asked for, minus the exclusions, which
/// is why an exclusion beats an inclusion of the same type.
fn selected_types(arguments: &FilterArguments) -> Vec<VariantType> {
    let all = [
        VariantType::NoVariation,
        VariantType::Snp,
        VariantType::Mnp,
        VariantType::Indel,
        VariantType::Symbolic,
        VariantType::Mixed,
    ];
    let included: Vec<VariantType> = if arguments.types_to_include.is_empty() {
        all.to_vec()
    } else {
        arguments.types_to_include.clone()
    };
    included
        .into_iter()
        .filter(|kind| !arguments.types_to_exclude.contains(kind))
        .collect()
}

/// `--set-filtered-gt-to-nocall`, `--drop-info-annotation` and `--drop-genotype-annotation`.
#[derive(Debug, Clone, Default)]
pub struct OutputArguments {
    pub set_filtered_genotypes_to_no_call: bool,
    pub info_annotations_to_drop: Vec<String>,
    pub genotype_annotations_to_drop: Vec<String>,
}

/// `setFilteredGenotypeToNocall`, which replaces the call and keeps the filter that caused it.
///
/// The counts are recomputed **here** rather than left stale, and only when something was actually
/// replaced: `calculateChromosomeCounts` is called with `removeStaleValues`, so a whole-cohort run
/// that rewrites nothing else still rewrites this record's AC and AN, and adds the AF an input
/// without one never had.
pub fn set_filtered_genotypes_to_no_call(record: &mut Record) {
    let mut replaced = false;
    for genotype in &mut record.variant.genotypes {
        let called = !genotype.alleles.is_empty() && genotype.alleles.iter().any(Option::is_some);
        if called && is_filtered(genotype) {
            replaced = true;
            genotype.alleles = vec![None; genotype.alleles.len()];
        }
    }
    if replaced {
        calculate_chromosome_counts(&mut record.variant);
    }
}

/// `dropAnnotations`, which is a no-op when nothing is named.
///
/// The genotype keys it can reach are the EXTENDED attributes alone: GT, GQ, DP, AD and PL live in
/// their own fields and survive whatever is asked for.
pub fn drop_annotations(record: &mut Record, arguments: &OutputArguments) {
    if arguments.info_annotations_to_drop.is_empty()
        && arguments.genotype_annotations_to_drop.is_empty()
    {
        return;
    }
    record
        .variant
        .attributes
        .retain(|(key, _)| !arguments.info_annotations_to_drop.contains(key));
    if arguments.genotype_annotations_to_drop.is_empty() {
        return;
    }
    for genotype in &mut record.variant.genotypes {
        genotype
            .attributes
            .retain(|(key, _)| !arguments.genotype_annotations_to_drop.contains(key));
    }
}

/// The writer's queue: records are held until the record being read is at or past them.
///
/// `apply` drains with `<=` against the current record's start, or entirely when the contig
/// changes, then adds the record it just finished. `onTraversalSuccess` drains the rest. The queue
/// exists because trimming moves a record RIGHT, so a file written in the order it was read would
/// not be sorted; it is the tool repairing an order it broke itself.
#[derive(Debug, Default)]
pub struct PendingWriter {
    pending: Vec<Record>,
}

impl PendingWriter {
    pub fn new() -> PendingWriter {
        PendingWriter {
            pending: Vec::new(),
        }
    }

    /// What is written before `record` is read, in order.
    pub fn drain_before(&mut self, contig: &str, start: i32) -> Vec<Record> {
        let mut written = Vec::new();
        // `PriorityQueue.peek` is the smallest start; ties keep insertion order here, which is the
        // order the reference's heap gives for equal keys of a two-element comparison.
        while let Some(head) = self.head() {
            let same_contig = self.pending[head].variant.contig == contig;
            if same_contig && self.pending[head].variant.start > start {
                break;
            }
            written.push(self.pending.remove(head));
        }
        written
    }

    /// The record joins the queue rather than the file.
    pub fn add(&mut self, record: Record) {
        self.pending.push(record);
    }

    /// `onTraversalSuccess`: whatever is left, in start order.
    pub fn drain(&mut self) -> Vec<Record> {
        let mut written = Vec::new();
        while let Some(head) = self.head() {
            written.push(self.pending.remove(head));
        }
        written
    }

    fn head(&self) -> Option<usize> {
        self.pending
            .iter()
            .enumerate()
            .min_by_key(|(index, record)| (record.variant.start, *index))
            .map(|(index, _)| index)
    }
}

/// `--discordance` and `--concordance`, which need the other file's records at the same position.
#[derive(Debug, Clone, Default)]
pub struct ComparisonArguments {
    pub discordance_only: bool,
    pub concordance_only: bool,
    /// `--exclude-filtered`, which both comparisons read and neither is changed by. See
    /// [`have_same_genotypes`].
    pub exclude_filtered: bool,
}

/// `haveSameGenotypes`, which compares ALLELE SETS and refuses anything filtered.
///
/// Two things about it are worth stating, because both are surprising and both are measured:
///
///  * the comparison is `a1.containsAll(a2) && a2.containsAll(a1)`, so multiplicity is invisible:
///    `0/1` matches `1/0`, and `1/1` matches a haploid `1`;
///  * **a filtered genotype never matches anything, including another filtered genotype.** The
///    first clause is `g1.isCalled() && g2.isFiltered()`, and `isCalled()` is about alleles rather
///    than filters, so a genotype carrying an FT is still called. Two identically filtered
///    genotypes therefore take that clause and are declared different, which makes the third
///    clause, `g1.isFiltered() && g2.isFiltered() && excludeFiltered`, unreachable. It is written
///    out here anyway, in the reference's order, because a reader who deletes it will conclude the
///    behaviour is a bug in the port rather than in what it ports.
pub fn have_same_genotypes(
    left: Option<(&Genotype, &[Allele])>,
    right: Option<(&Genotype, &[Allele])>,
    exclude_filtered: bool,
) -> bool {
    let (Some((left, left_alleles)), Some((right, right_alleles))) = (left, right) else {
        return false;
    };
    let called = |genotype: &Genotype| {
        !genotype.alleles.is_empty() && genotype.alleles.iter().any(Option::is_some)
    };
    if (called(left) && is_filtered(right))
        || (called(right) && is_filtered(left))
        || (is_filtered(left) && is_filtered(right) && exclude_filtered)
    {
        return false;
    }
    let bases = |genotype: &Genotype, alleles: &[Allele]| -> Vec<Vec<u8>> {
        genotype
            .alleles
            .iter()
            .map(|allele| match allele {
                Some(index) => alleles[*index].bases.clone(),
                None => b".".to_vec(),
            })
            .collect()
    };
    let left_bases = bases(left, left_alleles);
    let right_bases = bases(right, right_alleles);
    left_bases.iter().all(|allele| right_bases.contains(allele))
        && right_bases.iter().all(|allele| left_bases.contains(allele))
}

/// `isDiscordant`: without a sample it is only "the other file has nothing here".
pub fn is_discordant(
    record: &Record,
    others: &[Record],
    selection: &SampleSelection,
    arguments: &ComparisonArguments,
) -> bool {
    if selection.no_samples_specified {
        return others.is_empty();
    }
    for (index, name) in record.samples.iter().enumerate() {
        if !selection.samples.contains(name) {
            continue;
        }
        let genotype = &record.variant.genotypes[index];
        // `sampleHasVariant`, whose `isFiltered && !excludeFiltered` half is unreachable for the
        // same reason as above: a filtered genotype with alleles is still called.
        let hom_ref =
            !genotype.alleles.is_empty() && genotype.alleles.iter().all(|a| *a == Some(0));
        let called = !genotype.alleles.is_empty() && genotype.alleles.iter().any(Option::is_some);
        if hom_ref || !(called || (is_filtered(genotype) && !arguments.exclude_filtered)) {
            continue;
        }
        let found = others.iter().any(|other| {
            have_same_genotypes(
                Some((genotype, &record.variant.alleles)),
                genotype_of(other, name),
                arguments.exclude_filtered,
            )
        });
        if !found {
            return true;
        }
    }
    false
}

/// `isConcordant`: every selected sample must match somewhere, which is not the negation of
/// discordance.
pub fn is_concordant(
    record: &Record,
    others: &[Record],
    selection: &SampleSelection,
    arguments: &ComparisonArguments,
) -> bool {
    if others.is_empty() {
        return false;
    }
    if selection.no_samples_specified {
        return true;
    }
    for (index, name) in record.samples.iter().enumerate() {
        if !selection.samples.contains(name) {
            continue;
        }
        let genotype = &record.variant.genotypes[index];
        let found = others.iter().any(|other| {
            have_same_genotypes(
                Some((genotype, &record.variant.alleles)),
                genotype_of(other, name),
                arguments.exclude_filtered,
            )
        });
        if !found {
            return false;
        }
    }
    true
}

fn genotype_of<'a>(record: &'a Record, sample: &str) -> Option<(&'a Genotype, &'a [Allele])> {
    record
        .samples
        .iter()
        .position(|name| name == sample)
        .map(|index| {
            (
                &record.variant.genotypes[index],
                &record.variant.alleles[..],
            )
        })
}
