//! `NuMTFilterTool`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NuMTFilterTool` (GATK 4.6.2.0).
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

use jmath::poisson::{inverse_cumulative_probability, PoissonError};

/// `LOWER_BOUND_PROB`.
const LOWER_BOUND_PROB: f64 = 0.01;

/// `GATKVCFConstants.POSSIBLE_NUMT_FILTER_NAME`.
pub const FILTER_NAME: &str = "possible_numt";

/// `GATKVCFConstants.SITE_LEVEL_FILTERS`, the placeholder an unfiltered allele carries.
pub const SITE_LEVEL_FILTERS: &str = "SITE";

/// The arguments, with the tool's own defaults.
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

/// What the run refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum NuMTError {
    /// `Utils.validateArg(isFiltered.size() == alleleFilters.size(), ...)`.
    ListsNotTheSameSize,
    /// The Poisson underneath refused.
    Poisson(PoissonError),
}

impl NuMTError {
    pub fn java_class(&self) -> &str {
        match self {
            NuMTError::ListsNotTheSameSize => "java.lang.IllegalArgumentException",
            NuMTError::Poisson(_) => "java.lang.IllegalArgumentException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            NuMTError::ListsNotTheSameSize => "lists are not the same size".to_string(),
            NuMTError::Poisson(error) => format!("{error:?}"),
        }
    }
}

/// `getMaxAltDepthCutoff`: the 99th percentile of a Poisson whose mean is half the expected NuMT
/// coverage.
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
pub fn decode_as_filters(attribute: Option<&str>) -> Vec<Vec<String>> {
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

/// `encodeASFilters`.
pub fn encode_as_filters(filters: &[Vec<String>]) -> String {
    filters
        .iter()
        .map(|allele| allele.join(","))
        .collect::<Vec<String>>()
        .join("|")
}

/// `addAlleleFilters`: the SITE placeholder is REPLACED, and anything else is unioned through a
/// `LinkedHashSet`, so insertion order is kept and a repeat is dropped.
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
pub fn merged_as_filter_string(
    record: &Record,
    is_filtered: &[bool],
    filter_name: &str,
) -> Result<String, NuMTError> {
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
fn depths_by_alternate(record: &Record) -> Vec<Vec<i32>> {
    let mut by_allele: Vec<Vec<i32>> = vec![Vec::new(); record.alleles.len()];
    for depths in record.allele_depths.iter().flatten() {
        // The two iterators are walked together and stop at the shorter, so a genotype whose AD is
        // too short fills only the alleles it reaches.
        for (slot, depth) in by_allele.iter_mut().zip(depths) {
            slot.push(*depth);
        }
    }
    // The reference is the first allele, and the tool drops it by equality rather than by index.
    by_allele.into_iter().skip(1).collect()
}

/// `apply`: the record as the writer receives it.
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
    if applied.contains(&true) {
        out.as_filter_status = Some(merged_as_filter_string(record, &applied, FILTER_NAME)?);
    }
    Ok(out)
}
