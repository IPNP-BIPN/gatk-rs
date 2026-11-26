//! `NuMTFilterTool`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NuMTFilterTool` (GATK 4.7.0.0;
//! `NuMTFilterTool.java` is unchanged since 4.6.2.0).
//!
//! A mitochondrial call filtered when its alternate depth is low enough to be a nuclear insertion
//! of mitochondrial DNA rather than a real heteroplasmy.
//!
//! # It does nothing at its defaults
//!
//! ```java
//! if (maxNuMTAutosomalCopies > 0 && medianAutosomalCoverage > 0) {
//!     maxAltDepthCutoff = getMaxAltDepthCutoff(maxNuMTAutosomalCopies, medianAutosomalCoverage);
//! }
//! ```
//!
//! `medianAutosomalCoverage` defaults to **zero**, so the cutoff stays at its initial `0`,
//! `max(AD) < 0` is never true, and the output is the input with one more header line. The tool
//! has to be told the coverage to do anything at all.
//!
//! # An ordinary VCF makes it throw
//!
//! `getMergedASFilterString` validates that the decoded `AS_FilterStatus` has one entry per
//! alternate allele. An absent attribute decodes to an **empty** list, so a VCF that has never
//! been through Mutect2's filtering fails with `lists are not the same size` as soon as anything
//! is filtered. The same file passes untouched while the cutoff is zero: the failure depends on
//! the data and on an argument, not on the file's shape alone.
//!
//! # A record with no alternate allele is filtered
//!
//! ```java
//! if (!appliedFilter.contains(Boolean.FALSE)) { vcb.filter(filterName()); }
//! ```
//!
//! The site filter is applied when no alternate escapes it, and an **empty** list of decisions
//! escapes nothing, so `A -> .` comes out carrying `possible_numt`.
//!
//! # The depth compared is the maximum across samples
//!
//! Not the sum: two samples of fifty each, under a cutoff of seventy-nine, are filtered while
//! their sum is a hundred. And a genotype without `AD` is skipped by the precondition, so an
//! allele nobody reports leaves an empty list whose maximum is taken as zero and is therefore
//! filtered.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing" (`Some(x)` or `None`), Java's nullable reference checked
//!   by the compiler. `Result<T, E>` is "a value or an error" (`Ok(x)` or `Err(e)`); a reference
//!   exception is modelled as an `Err`.
//! - `?` after a call returns its `Err` to the caller at once, as an uncaught Java exception would.
//! - `.iter().map(..).collect()` is Java's `stream().map(..).collect(..)`; `|x| ...` is a lambda.
//!   `a.iter().zip(b)` walks two sequences in step and stops at the shorter.
//! - `&[T]` and `&str` are borrowed, read-only views of a list and of a string.

use jmath::poisson::{inverse_cumulative_probability, PoissonError};

/// `LOWER_BOUND_PROB`: the cutoff is the Poisson quantile at `1 - LOWER_BOUND_PROB`.
const LOWER_BOUND_PROB: f64 = 0.01;

/// `GATKVCFConstants.POSSIBLE_NUMT_FILTER_NAME`.
pub const FILTER_NAME: &str = "possible_numt";

/// `GATKVCFConstants.SITE_LEVEL_FILTERS`, the placeholder an unfiltered allele carries.
pub const SITE_LEVEL_FILTERS: &str = "SITE";

/// The arguments, with the tool's own defaults: `--autosomal-coverage`
/// (`medianAutosomalCoverage`, default `DEFAULT_MEDIAN_AUTOSOMAL_COVERAGE = 0`) and
/// `--max-numt-autosomal-copies` (`maxNuMTAutosomalCopies`, default
/// `DEFAULT_MAX_NUMT_AUTOSOMAL_COPIES = 4`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Arguments {
    pub median_autosomal_coverage: f64,
    pub max_numt_autosomal_copies: f64,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            // Zero, which is why the tool does nothing until it is told otherwise.
            median_autosomal_coverage: 0.0,
            max_numt_autosomal_copies: 4.0,
        }
    }
}

/// What the run refuses, each with the Java exception class and message the reference would
/// report.
#[derive(Debug, Clone, PartialEq)]
pub enum NuMTError {
    /// `Utils.validateArg(isFiltered.size() == alleleFilters.size(), ...)`.
    ListsNotTheSameSize,
    /// The Poisson underneath refused. Its message is the Rust error's debug rendering, not a
    /// commons-math message.
    Poisson(PoissonError),
}

impl NuMTError {
    /// The Java exception class. `match` must name every variant, so a new one cannot be added
    /// without deciding its class here.
    pub fn java_class(&self) -> &str {
        match self {
            NuMTError::ListsNotTheSameSize => "java.lang.IllegalArgumentException",
            NuMTError::Poisson(_) => "java.lang.IllegalArgumentException",
        }
    }

    /// The exception's message.
    pub fn message(&self) -> String {
        match self {
            NuMTError::ListsNotTheSameSize => "lists are not the same size".to_string(),
            NuMTError::Poisson(error) => format!("{error:?}"),
        }
    }
}

/// `getMaxAltDepthCutoff`: the 99th percentile of a Poisson whose mean is half the expected NuMT
/// coverage.
///
/// ```java
/// final PoissonDistribution autosomalCoverage =
///         new PoissonDistribution(medianAutosomalCoverage * maxNuMTAutosomalCopies / 2.0);
/// return autosomalCoverage.inverseCumulativeProbability(1 - LOWER_BOUND_PROB);
/// ```
///
/// The mean is computed in the reference's operand order, coverage times copies, then halved.
/// The quantile is `jmath`'s port of commons-math's `inverseCumulativeProbability`.
pub fn max_alt_depth_cutoff(
    max_numt_autosomal_copies: f64,
    median_autosomal_coverage: f64,
) -> Result<i32, PoissonError> {
    inverse_cumulative_probability(
        median_autosomal_coverage * max_numt_autosomal_copies / 2.0,
        1.0 - LOWER_BOUND_PROB,
    )
}

/// `onTraversalStart`'s one computation: the cutoff, or zero when either argument is not positive.
///
/// Zero is the reference's field initialiser, `private int maxAltDepthCutoff = 0;`, which the
/// guarded assignment quoted in the module header leaves in place.
pub fn cutoff_for(arguments: &Arguments) -> Result<i32, PoissonError> {
    if arguments.max_numt_autosomal_copies > 0.0 && arguments.median_autosomal_coverage > 0.0 {
        return max_alt_depth_cutoff(
            arguments.max_numt_autosomal_copies,
            arguments.median_autosomal_coverage,
        );
    }
    Ok(0)
}

/// One record, reduced to what this filter reads and writes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// One entry per allele, reference first, as the VCF declares them.
    pub alleles: Vec<String>,
    /// Per genotype, the `AD` list if it has one. A genotype without `AD` is `None` and is skipped.
    pub allele_depths: Vec<Option<Vec<i32>>>,
    /// The `FILTER` column, empty when the record is unfiltered.
    pub filters: Vec<String>,
    /// The `AS_FilterStatus` attribute, absent when the record has none.
    pub as_filter_status: Option<String>,
}

/// `decodeASFilters`: `|` between alleles, `,` within, and an absent attribute decodes to an EMPTY
/// list rather than to one empty entry per allele.
///
/// ```java
/// return AnnotationUtils.decodeAnyASListWithRawDelim(vc.getCommonInfo()
///                 .getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, "")).stream()
///         .map(filters -> AnnotationUtils.decodeAnyASList(filters).stream()
///                 .map(String::trim).collect(Collectors.toList()))
///         .collect(Collectors.toList());
/// ```
///
/// Both decoders strip `[` and `]` first (`replaceAll(BRACKET_REGEX, "")`), which is the
/// `replace(['[', ']'], "")` here; the brackets are what a list-valued attribute's string form
/// carries.
pub fn decode_as_filters(attribute: Option<&str>) -> Vec<Vec<String>> {
    // `unwrap_or("")` is `getAttributeAsString(key, "")`: an absent attribute reads as "".
    let text = attribute.unwrap_or("").replace(['[', ']'], "");
    if text.is_empty() {
        // `StringUtils.splitByWholeSeparatorPreserveAllTokens("", "|")` is an EMPTY array, which is
        // what makes an ordinary VCF fail the length check.
        return Vec::new();
    }
    text.split('|')
        .map(|allele| {
            allele
                .split(',')
                .map(|filter| filter.trim().to_string())
                .collect()
        })
        .collect()
}

/// `encodeASFilters`: `,` within an allele, `|` between alleles.
///
/// The reference's `encodeAnyASListWithRawDelim` also strips brackets after joining; nothing
/// here can contain one, since decoding removed them and the filter name has none.
pub fn encode_as_filters(filters: &[Vec<String>]) -> String {
    filters
        .iter()
        .map(|allele| allele.join(","))
        .collect::<Vec<String>>()
        .join("|")
}

/// `addAlleleFilters`: the SITE placeholder is REPLACED, and anything else is unioned through a
/// `LinkedHashSet`, so insertion order is kept and a repeat is dropped.
///
/// ```java
/// } else if (currentAlleleFilters.isEmpty() || (currentAlleleFilters.size() == 1
///         && currentAlleleFilters.contains(GATKVCFConstants.SITE_LEVEL_FILTERS))) {
/// ```
///
/// The reference's first branch, an empty `newFilters`, cannot arise: this tool always adds
/// exactly one name.
fn add_allele_filters(current: &[String], new_filter: &str) -> Vec<String> {
    if current.is_empty() || (current.len() == 1 && current[0] == SITE_LEVEL_FILTERS) {
        return vec![new_filter.to_string()];
    }
    let mut updated: Vec<String> = current.to_vec();
    if !updated.iter().any(|filter| filter == new_filter) {
        updated.push(new_filter.to_string());
    }
    updated
}

/// `getMergedASFilterString`, whose length check is what an ordinary VCF fails.
///
/// ```java
/// Utils.validateArg(isFiltered.size() == alleleFilters.size(), "lists are not the same size");
/// ```
///
/// `is_filtered` has one entry per alternate allele, and each `true` adds `filter_name` to that
/// allele's list; a `false` leaves the list exactly as decoded.
pub fn merged_as_filter_string(
    record: &Record,
    is_filtered: &[bool],
    filter_name: &str,
) -> Result<String, NuMTError> {
    // `as_deref()` turns the record's `Option<String>` into the `Option<&str>` view the decoder
    // borrows, without copying the text.
    let allele_filters = decode_as_filters(record.as_filter_status.as_deref());
    if is_filtered.len() != allele_filters.len() {
        return Err(NuMTError::ListsNotTheSameSize);
    }
    let updated: Vec<Vec<String>> = allele_filters
        .iter()
        .zip(is_filtered)
        .map(|(filters, filtered)| {
            if *filtered {
                add_allele_filters(filters, filter_name)
            } else {
                filters.clone()
            }
        })
        .collect();
    Ok(encode_as_filters(&updated))
}

/// `getDataByAllele` reduced to what this tool asks for: per alternate allele, the depths reported
/// by the genotypes that have `AD` at all.
///
/// ```java
/// while(alleleDataIterator.hasNext() && dataByAlleleIterator.hasNext())
///     dataByAlleleIterator.next().add(alleleDataIterator.next());
/// ```
///
/// One list per allele in VCF order (the reference's `LinkedHashMap` keyed by allele), then the
/// reference allele's list removed.
fn depths_by_alternate(record: &Record) -> Vec<Vec<i32>> {
    let mut by_allele: Vec<Vec<i32>> = vec![Vec::new(); record.alleles.len()];
    // `.flatten()` over `Option`s yields only the `Some` values: Java's
    // `.filter(Genotype::hasAD)` precondition.
    for depths in record.allele_depths.iter().flatten() {
        // The two iterators are walked together and stop at the shorter, so a genotype whose AD is
        // too short fills only the alleles it reaches.
        for (slot, depth) in by_allele.iter_mut().zip(depths) {
            slot.push(*depth);
        }
    }
    // The reference is the first allele, and the tool drops it by equality rather than by index:
    // `.filter(entry -> !variant.getReference().equals(entry.getKey()))`. Skipping the first
    // slot is the same thing because the reference allele is always first, and no alternate
    // equals it (an `Allele`'s equality includes its reference flag).
    by_allele.into_iter().skip(1).collect()
}

/// `apply`: the record as the writer receives it.
///
/// ```java
/// .map(entry -> entry.getValue().stream().max(Integer::compare).orElse(0) < maxAltDepthCutoff)
///         .collect(Collectors.toList());
/// ```
///
/// Per alternate allele, `true` when the largest depth any sample reports is under the cutoff.
/// `max()` yields an `Option`, `None` for an empty list, and `unwrap_or(0)` is `orElse(0)`.
pub fn apply(record: &Record, cutoff: i32) -> Result<Record, NuMTError> {
    let applied: Vec<bool> = depths_by_alternate(record)
        .iter()
        .map(|depths| depths.iter().max().copied().unwrap_or(0) < cutoff)
        .collect();
    let mut out = record.clone();
    // An EMPTY list contains no `false` either, so a record with no alternate allele is filtered.
    if !applied.contains(&false) {
        out.filters.push(FILTER_NAME.to_string());
    }
    // `AS_FilterStatus` is rewritten only when some allele is filtered, which is why the length
    // check, and its failure, depends on the cutoff.
    if applied.contains(&true) {
        out.as_filter_status = Some(merged_as_filter_string(record, &applied, FILTER_NAME)?);
    }
    Ok(out)
}
