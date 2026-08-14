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
use gatk_engine::subset_alleles::{subset_alleles, AssignmentMethod, Genotype};
use gatk_engine::variant_context_utils::{trim_alleles, Variant};

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
