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

use gatk_engine::java_regex::{self, PatternSyntaxError};

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
    let mut vcf_samples: Vec<String> = header_samples.to_vec();
    vcf_samples.sort();
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
struct SortedSet {
    values: Vec<String>,
}

impl SortedSet {
    fn new() -> SortedSet {
        SortedSet { values: Vec::new() }
    }

    fn extend(&mut self, values: &[String]) {
        for value in values {
            if let Err(at) = self.values.binary_search(value) {
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
