//! Ported from `org.broadinstitute.hellbender.tools.walkers.qc.Pileup` (GATK 4.6.2.0).
//!
//! One samtools-style line per locus, with three optional columns. Everything here is a format
//! rather than a computation, which is why the golden compares whole files.
//!
//! # Three filters on top of the walker's
//!
//! `getDefaultReadFilters` calls `super` and appends `NotDuplicate`, `PassesVendorQualityCheck` and
//! `NotSecondaryAlignment`. A duplicate that `CountReads` counts never reaches this tool, and a
//! window covered only by one produces an EMPTY FILE rather than a line with no bases, because a
//! locus with no reads is not emitted at all.
//!
//! # The deletion filter runs before everything, including the deletion count
//!
//! `makeFilteredPileup(pe -> !pe.isDeletion())` happens first, and the verbose column then counts
//! deletions in THAT pileup. So the count is always zero, and a locus covered only by a deletion
//! prints an empty base string with the reference base still in place.
//!
//! # The features column is always there
//!
//! `String.format("%s %s", pileupString, features)` with no metadata gives a line ending in a
//! space. Every line in the golden's plain case does.

use gatk_engine::read_pileup::ReadPileup;
use gatk_readfilter::{self as filters, with_header};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// `VERBOSE_DELIMITER`, which the reference's own comment calls ugly.
pub const VERBOSE_DELIMITER: char = '@';

/// The filters this tool adds to `LocusWalker`'s, in its order.
pub const ADDITIONAL_READ_FILTERS: [&str; 3] = [
    "NotDuplicateReadFilter",
    "PassesVendorQualityCheckReadFilter",
    "NotSecondaryAlignmentReadFilter",
];

/// `getDefaultReadFilters`: the walker's two, then this tool's three.
pub fn default_read_filter(read: &BamRecord, header: &SamHeader) -> bool {
    with_header::wellformed(read, header)
        && filters::mapped(read)
        && filters::not_duplicate(read)
        && filters::passes_vendor_quality_check(read)
        && filters::not_secondary_alignment(read)
}

/// `insertLengthOutput`: one fragment length per read, comma joined, in pileup order.
pub fn insert_length_output(pileup: &ReadPileup) -> String {
    pileup
        .elements
        .iter()
        .map(|element| element.read.inferred_insert_size.to_string())
        .collect::<Vec<String>>()
        .join(",")
}

/// `createVerboseOutput`: the deletion count, a space, then one entry per read.
///
/// The count is of the pileup it is handed, which the caller has already stripped of deletions, so
/// it is zero on every line this tool writes. Kept as a count rather than a constant because the
/// function is `@VisibleForTesting` and the reference tests it on unfiltered pileups.
pub fn create_verbose_output(pileup: &ReadPileup) -> String {
    let mut text = String::new();
    text.push_str(
        &pileup
            .number_of_elements(|element| element.is_deletion())
            .to_string(),
    );
    text.push(' ');
    let entries: Vec<String> = pileup
        .elements
        .iter()
        .map(|element| {
            format!(
                "{}{VERBOSE_DELIMITER}{}{VERBOSE_DELIMITER}{}{VERBOSE_DELIMITER}{}",
                element.read.read_name,
                element.offset,
                element.read.read_bases.len(),
                element.read.mapping_quality
            )
        })
        .collect();
    text.push_str(&entries.join(","));
    text
}

/// `getFeaturesString`: the overlapping features, or an empty string.
///
/// The brackets are only added when there is something inside them, and the caller prints the
/// result either way, which is what leaves a trailing space on a line with no metadata.
pub fn features_string(features: &[String]) -> String {
    if features.is_empty() {
        return String::new();
    }
    format!("[Feature(s): {}]", features.join(", "))
}

/// `apply`: the line for one locus, newline included.
///
/// `reference_base` is `'N'` when the run has no reference, which the reference decides with
/// `hasReference()` rather than by looking at the context.
pub fn line(
    pileup: &ReadPileup,
    reference_base: char,
    features: &[String],
    output_insert_length: bool,
    show_verbose: bool,
) -> String {
    // The deletions go before anything is printed, and everything below sees the filtered pileup.
    let filtered = pileup.filtered(|element| !element.is_deletion());
    let mut text = format!(
        "{} {}",
        filtered.pileup_string(reference_base),
        features_string(features)
    );
    if output_insert_length {
        text.push(' ');
        text.push_str(&insert_length_output(&filtered));
    }
    if show_verbose {
        text.push(' ');
        text.push_str(&create_verbose_output(&filtered));
    }
    text.push('\n');
    text
}
