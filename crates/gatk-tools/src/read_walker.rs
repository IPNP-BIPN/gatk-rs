//! Ported from `org.broadinstitute.hellbender.engine.ReadWalker` and the parts of `GATKTool` it
//! runs on (GATK 4.6.2.0).
//!
//! Every read-based tool is this traversal, so what it hands to `apply` is what the tool sees.
//! Three things in it are decisions rather than plumbing:
//!
//!  * **the filter sits between the two transformers.** `getTransformedReadStream` is
//!    `map(preTransformer).filter(filter).map(postTransformer)`, so the filter judges a read the
//!    walker never sees, and the walker sees a read the filter never judged;
//!  * **`getReadInterval` returns null for a mapped read too.** It is null when the read is
//!    unmapped *or* when its coordinates do not form a valid `SimpleInterval`. A mapped read with
//!    an empty cigar has an alignment end one before its start, so it reaches `apply` with an
//!    empty `ReferenceContext`: a read with no reference under it, and no error anywhere;
//!  * **intervals change which reads exist, not just which are reported.** With `-L` the data
//!    source is bounded before the traversal starts, and unplaced unmapped reads are excluded
//!    unless the user asked for them by name.

use gatk_engine::context::ReferenceContext;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::read;
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_engine::reference::ReferenceFileSource;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// One `apply` call: the read, and the reference context built from its interval.
pub struct Applied {
    pub read: BamRecord,
    pub context: ReferenceContext,
}

/// `ReadWalker.getReadInterval`.
///
/// Both halves matter. `isUnmapped` is the adapter's three criteria, not the flag; and
/// `SimpleInterval.isValid` is what excludes a mapped read whose end precedes its start.
pub fn read_interval(record: &BamRecord, header: &SamHeader) -> Option<SimpleInterval> {
    if read::is_unmapped(record) {
        return None;
    }
    let contig = header.sequences.get(record.reference_index as usize)?;
    SimpleInterval::new(&contig.name, record.alignment_start, record.alignment_end())
}

/// The reads a traversal presents, in order.
///
/// `intervals` empty is an unbounded traversal, which includes the unplaced unmapped reads at the
/// tail; with intervals the data source is bounded and they are excluded, because
/// `getTraversalParameters` only sets `traverseUnmappedReads` when the user names them.
pub fn traverse(
    source: &ReadsDataSource,
    intervals: &[SimpleInterval],
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<Vec<BamRecord>, ReadsError> {
    traverse_with_bounds(source, intervals, false, filter)
}

/// The same, with `-L unmapped` honoured.
///
/// `ReadsPathDataSource.setTraversalBounds` makes a traversal bounded when it has intervals **or**
/// when unmapped reads were asked for, and `SamReaderQueryingIterator.loadNextIterator` runs the
/// interval query first and the unmapped query second. So the order is not a choice: the unplaced
/// reads arrive as a tail after every interval, and `-L unmapped` on its own is a bounded
/// traversal of nothing but that tail rather than an unbounded one.
///
/// The distinction the reference draws twice, in two comments, is between an unmapped read with no
/// position and an unmapped read carrying its mate's: only the first is in this tail, the second is
/// returned by an interval query overlapping that position.
pub fn traverse_with_bounds(
    source: &ReadsDataSource,
    intervals: &[SimpleInterval],
    traverse_unmapped: bool,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<Vec<BamRecord>, ReadsError> {
    traverse_with_bounds_mut(source, intervals, traverse_unmapped, &mut |read| {
        filter(read)
    })
}

/// The same traversal, with a filter that may keep state.
///
/// The reference's filter always does: `getTransformedReadStream` is given a `CountingReadFilter`,
/// whose counters are what the summary line is made of, so "the filter is a pure predicate" is an
/// assumption of the callers here rather than of the engine. The multi-pass read walker is the
/// first caller that needs the counts, and it goes through this form; the one above stays for the
/// tools that pass a predicate and read nothing back.
pub fn traverse_with_bounds_mut(
    source: &ReadsDataSource,
    intervals: &[SimpleInterval],
    traverse_unmapped: bool,
    filter: &mut dyn FnMut(&BamRecord) -> bool,
) -> Result<Vec<BamRecord>, ReadsError> {
    let bounded = !intervals.is_empty() || traverse_unmapped;
    let records = if !bounded {
        source.iter_all()?
    } else {
        let mut records = if intervals.is_empty() {
            Vec::new()
        } else {
            source.query(intervals)?
        };
        if traverse_unmapped {
            records.extend(source.query_unmapped()?);
        }
        records
    };
    // map(pre).filter(f).map(post): both transformers are the identity for a tool that declares
    // none, which is every tool measured here. The order is kept because a tool that declares one
    // changes what the filter sees, not only what apply sees.
    Ok(records.into_iter().filter(|read| filter(read)).collect())
}

/// The traversal with its reference contexts, which is what `apply` actually receives.
pub fn traverse_with_reference(
    source: &ReadsDataSource,
    reference: Option<&mut ReferenceFileSource>,
    intervals: &[SimpleInterval],
    traverse_unmapped: bool,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<Vec<Applied>, ReadsError> {
    let header = source.header().clone();
    let records = traverse_with_bounds(source, intervals, traverse_unmapped, filter)?;

    let mut applied = Vec::with_capacity(records.len());
    match reference {
        // No reference: the context still carries the read's interval as its window, and only its
        // *bases* are empty. Measured, not assumed: the golden shows the reference reporting
        // `chr1:10-19` with no bases for a read a tool ran over without `-R`.
        None => {
            for read in records {
                let interval = read_interval(&read, &header);
                let context = ReferenceContext::without_source(interval, 0, 0)
                    .unwrap_or_else(|_| ReferenceContext::empty());
                applied.push(Applied { read, context });
            }
        }
        Some(reference) => {
            for read in records {
                let interval = read_interval(&read, &header);
                // `new ReferenceContext(reference, readInterval)` is windowless: the window is the
                // read's own span, with no extra bases either side.
                let context = ReferenceContext::new(reference, interval, 0, 0)
                    .unwrap_or_else(|_| ReferenceContext::empty());
                applied.push(Applied { read, context });
            }
        }
    }
    Ok(applied)
}
