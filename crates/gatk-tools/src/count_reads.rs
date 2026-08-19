//! Ported from `org.broadinstitute.hellbender.tools.CountReads` and
//! `org.broadinstitute.hellbender.tools.CountBases` (GATK 4.6.2.0).
//!
//! Two tools in one module because they are one tool with two increments: both are a `ReadWalker`
//! whose `apply` adds to a counter and whose `onTraversalSuccess` prints it. Everything that
//! decides their answer happens before `apply`, in [`crate::read_walker`] and in the filters.
//!
//! # The default filters are not nothing
//!
//! `ReadWalker.getDefaultReadFilters` is `WellformedReadFilter`, and neither tool adds to it or
//! takes from it. So a read with no read group and a read with an `N` cigar operator are absent
//! from both counts, and the count of a file is not the number of records in it: eight of eleven in
//! the golden's fixture, and eleven only under `--disable-tool-default-read-filters`.
//!
//! # `CountBases` counts the sequence, not the span
//!
//! `count += read.getLength()`, which is the number of bases the record carries. A read with a
//! ten-base deletion spans twenty reference bases and contributes ten; a read with no sequence at
//! all contributes zero while still being one read to `CountReads`.
//!
//! # Both print twice
//!
//! `OptionalTextOutputArgumentCollection.print` writes the number to `-O` when one was given, with
//! no trailing newline, and `CountReads` additionally logs it. The log line is not output and is
//! not reproduced here; the file is.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::read;
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_readfilter::with_header;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// `ReadWalker.getDefaultReadFilters`, which neither tool overrides.
pub const DEFAULT_READ_FILTERS: [&str; 1] = ["WellformedReadFilter"];

/// The conjunction [`DEFAULT_READ_FILTERS`] names, which is one filter.
pub fn default_read_filter(read: &BamRecord, header: &SamHeader) -> bool {
    with_header::wellformed(read, header)
}

/// `CountReads.doWork`: one per read that reaches `apply`.
pub fn count_reads(
    source: &ReadsDataSource,
    intervals: &[SimpleInterval],
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<i64, ReadsError> {
    Ok(crate::read_walker::traverse(source, intervals, filter)?.len() as i64)
}

/// `CountBases.doWork`: `read.getLength()` per read that reaches `apply`.
pub fn count_bases(
    source: &ReadsDataSource,
    intervals: &[SimpleInterval],
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<i64, ReadsError> {
    let reads = crate::read_walker::traverse(source, intervals, filter)?;
    Ok(reads.iter().map(|read| read::length(read) as i64).sum())
}

/// What either tool writes to `-O`: the number, and nothing else.
///
/// `print` rather than `println`, so there is no trailing newline. A port that used `println` would
/// write a file one byte longer than the reference's.
pub fn output(count: i64) -> String {
    count.to_string()
}
