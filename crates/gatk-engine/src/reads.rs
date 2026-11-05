//! Ported from `org.broadinstitute.hellbender.engine.ReadsPathDataSource` and
//! `utils.iterators.SamReaderQueryingIterator` (GATK 4.6.2.0), over htsjdk's
//! `BAMFileReader.createIndexIterator`, `BAMQueryMultipleIntervalsIteratorFilter`,
//! `QueryInterval.optimizeIntervals` and `Chunk.optimizeChunkList` (htsjdk 4.2.0).
//!
//! # Why an interval query is not "the reads that overlap the interval"
//!
//! The index says which *blocks* to read, and a filter says which records to keep. The filter is
//! the part with the decisions, and three of them are invisible from the outside:
//!
//!  * **It is stateful and single-pass.** `intervalIndex` only ever advances. An interval the
//!    scan has walked past can never match again, and once every interval is behind the current
//!    record the iteration *stops* rather than continuing to the end of the file. On a
//!    coordinate-sorted BAM that is an optimisation; on any other order it silently truncates.
//!  * **Unmapped reads that carry their mate's position are kept.** Their end would otherwise be
//!    `getAlignmentEnd()`, which is `0` for an unmapped record, so every one of them would sort
//!    before every interval and be dropped. The filter special-cases them to `end = start`.
//!  * **A mapped record with an empty cigar ends before it starts.** `getAlignmentEnd()` is
//!    `start + referenceLength - 1`, and a `*` cigar has reference length zero, so `end` is
//!    `start - 1` and the record fails `alignmentEnd < interval.start` for an interval that
//!    begins exactly at its start.
//!
//! Above that, GATK converts each `SimpleInterval` through the *reads* sequence dictionary and
//! runs `QueryInterval.optimizeIntervals`, which sorts and merges not only overlapping intervals
//! but **abutting** ones. Two adjacent intervals become one, so a record spanning the join is
//! returned once rather than twice, and the count of returned reads is not the sum of the
//! per-interval counts.
//!
//! # What is a dependency and what is a port
//!
//! The `.bai` bytes are parsed by [`htsjdk_bam::index::read_bai`], under the rule in
//! `docs/when-a-dependency-is-cheaper-than-a-port.md`: a bin's chunk list is what the format
//! says it is. Everything that decides *which records come back* is ported here: `regionToBins`,
//! the linear-index minimum offset, `optimizeChunkList`, the interval optimisation and the
//! record filter. Records are decompressed by `htsjdk-bgzf` and decoded by `htsjdk-bam`, which
//! are this programme's own ports, because what a record *is* is htsjdk's decision.

use std::path::Path;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bgzf::read::BgzfReader;
use htsjdk_bgzf::vfp;

use crate::interval::SimpleInterval;

/// `BAMFileConstants.FIXED_BLOCK_SIZE`: the part of a record before its name.
const FIXED_BLOCK_SIZE: usize = 32;

/// What a query refuses to answer, rather than answering wrongly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadsError {
    /// `IntervalUtils.convertSimpleIntervalToQueryInterval`: the contig is not in the reads
    /// sequence dictionary. A `UserException` there, not a silently empty result.
    ContigNotInDictionary(String),
    /// The BAM, or its index, could not be read at all.
    Io(String),
    /// The bytes are not a BAM, or a record does not decode.
    Malformed(String),
    /// A record's length was read and the bytes it promises are not there.
    ///
    /// `BinaryCodec.readBytes` reads until it has the length it was asked for and raises
    /// `RuntimeEOFException` when the stream ends first, counting both numbers. It is reached by
    /// seeking to an offset that is not a record boundary, which is what an index built for a
    /// DIFFERENT file does: the four bytes there read as a length of over a gigabyte and the file
    /// has a few hundred left. The reference names the file it was reading, so this carries it.
    PrematureEof {
        expected: i64,
        received: usize,
        file: String,
    },
}

/// `htsjdk.samtools.QueryInterval`: a reference *index*, not a name, and a 1-based closed span.
///
/// `end <= 0` means "to the end of the contig", which is why the comparison operators below
/// cannot simply compare the numbers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QueryInterval {
    pub reference_index: i32,
    pub start: i32,
    pub end: i32,
}

impl QueryInterval {
    /// `QueryInterval.compareTo`, including its treatment of `end == 0` as the largest end.
    ///
    /// Written out rather than derived because `end == 0` sorts *after* every other end on the
    /// same start, which no ordinary numeric ordering gives.
    pub fn compare(&self, other: &QueryInterval) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        match self.reference_index.cmp(&other.reference_index) {
            Ordering::Equal => {}
            other => return other,
        }
        match self.start.cmp(&other.start) {
            Ordering::Equal => {}
            other => return other,
        }
        if self.end == other.end {
            Ordering::Equal
        } else if self.end == 0 {
            Ordering::Greater
        } else if other.end == 0 {
            Ordering::Less
        } else {
            self.end.cmp(&other.end)
        }
    }

    /// `QueryInterval.endsAtStartOf`: abutting, which counts as mergeable.
    pub fn ends_at_start_of(&self, other: &QueryInterval) -> bool {
        self.reference_index == other.reference_index && self.end + 1 == other.start
    }

    /// `QueryInterval.overlaps`, with `end == 0` standing for `Integer.MAX_VALUE`.
    pub fn overlaps(&self, other: &QueryInterval) -> bool {
        if self.reference_index != other.reference_index {
            return false;
        }
        let this_end = if self.end == 0 { i32::MAX } else { self.end };
        let other_end = if other.end == 0 { i32::MAX } else { other.end };
        // CoordMath.overlaps
        (self.start >= other.start && self.start <= other_end)
            || (other.start >= self.start && other.start <= this_end)
    }
}

/// `QueryInterval.optimizeIntervals`: sort, then merge overlapping **and abutting** intervals.
///
/// The abutting case is the one worth naming: `1:100-200` and `1:201-300` come back as
/// `1:100-300`, so a read crossing 200/201 is returned once. htsjdk's query API requires this to
/// have been done and asserts it, so it is not optional.
pub fn optimize_intervals(intervals: &[QueryInterval]) -> Vec<QueryInterval> {
    if intervals.is_empty() {
        return Vec::new();
    }
    let mut sorted = intervals.to_vec();
    sorted.sort_by(QueryInterval::compare);

    let mut unique: Vec<QueryInterval> = Vec::new();
    let mut previous = sorted[0];
    for next in &sorted[1..] {
        if previous.ends_at_start_of(next) || previous.overlaps(next) {
            // A zero end on either side swallows the other: it already runs to the contig's end.
            let new_end = if previous.end == 0 || next.end == 0 {
                0
            } else {
                previous.end.max(next.end)
            };
            previous = QueryInterval {
                reference_index: previous.reference_index,
                start: previous.start,
                end: new_end,
            };
        } else {
            unique.push(previous);
            previous = *next;
        }
    }
    unique.push(previous);
    unique
}

/// `IntervalUtils.convertSimpleIntervalToQueryInterval`.
///
/// The dictionary is the **reads'** dictionary, so an interval naming a contig that exists in the
/// reference but not in the BAM header is a user error rather than an empty result.
pub fn convert_simple_interval_to_query_interval(
    interval: &SimpleInterval,
    header: &SamHeader,
) -> Result<QueryInterval, ReadsError> {
    let index = header
        .sequences
        .iter()
        .position(|s| s.name == interval.contig)
        .ok_or_else(|| ReadsError::ContigNotInDictionary(interval.contig.clone()))?;
    Ok(QueryInterval {
        reference_index: index as i32,
        start: interval.start,
        end: interval.end,
    })
}

/// `htsjdk.samtools.BAMIteratorFilter.IntervalComparison`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntervalComparison {
    /// The interval lies before the record: advance to the next interval.
    Before,
    /// The interval lies after the record: advance to the next record.
    After,
    Contained,
    Overlapping,
}

/// `BAMQueryMultipleIntervalsIteratorFilter.compareIntervalToRecord`.
///
/// The unmapped-with-a-position case is the reason this cannot be `interval.overlaps(read)`:
/// `getAlignmentEnd()` returns 0 for any record with the unmapped flag, so without the
/// special case every mate-placed unmapped read would compare `After` and never be returned.
pub fn compare_interval_to_record(
    interval: &QueryInterval,
    record: &BamRecord,
) -> IntervalComparison {
    let interval_end = if interval.end <= 0 {
        i32::MAX
    } else {
        interval.end
    };
    let alignment_end = if record.read_unmapped() && record.alignment_start != 0 {
        // Unmapped read at the coordinate of its mate.
        record.alignment_start
    } else {
        record.alignment_end()
    };

    if interval.reference_index < record.reference_index {
        IntervalComparison::Before
    } else if interval.reference_index > record.reference_index {
        IntervalComparison::After
    } else if interval_end < record.alignment_start {
        IntervalComparison::Before
    } else if alignment_end < interval.start {
        IntervalComparison::After
    } else if interval.start <= record.alignment_start && alignment_end <= interval_end {
        // CoordMath.encloses
        IntervalComparison::Contained
    } else {
        IntervalComparison::Overlapping
    }
}

/// `BAMIteratorFilter.FilteringIteratorState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterState {
    MatchesFilter,
    ContinueIteration,
    StopIteration,
}

/// `BAMQueryMultipleIntervalsIteratorFilter`: stateful, and that is the point.
///
/// `interval_index` never goes backwards, and `StopIteration` ends the traversal outright.
pub struct IntervalFilter<'a> {
    intervals: &'a [QueryInterval],
    contained: bool,
    interval_index: usize,
}

impl<'a> IntervalFilter<'a> {
    pub fn new(intervals: &'a [QueryInterval], contained: bool) -> Self {
        IntervalFilter {
            intervals,
            contained,
            interval_index: 0,
        }
    }

    pub fn compare_to_filter(&mut self, record: &BamRecord) -> FilterState {
        while self.interval_index < self.intervals.len() {
            match compare_interval_to_record(&self.intervals[self.interval_index], record) {
                IntervalComparison::Before => self.interval_index += 1,
                IntervalComparison::After => return FilterState::ContinueIteration,
                IntervalComparison::Contained => return FilterState::MatchesFilter,
                IntervalComparison::Overlapping => {
                    return if self.contained {
                        FilterState::ContinueIteration
                    } else {
                        FilterState::MatchesFilter
                    };
                }
            }
        }
        FilterState::StopIteration
    }
}

/// `htsjdk.samtools.Chunk`: a half-open span of BGZF virtual file pointers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chunk {
    pub start: u64,
    pub end: u64,
}

/// `Chunk.overlaps`, which is not the interval test it sounds like.
///
/// It sorts the pair first and then asks whether the earlier one's *end* is past the later
/// one's *start*, comparing block address then block offset. A virtual pointer is
/// `address << 16 | offset`, so that comparison is the pointers' numeric order.
fn chunks_overlap(a: &Chunk, b: &Chunk) -> bool {
    let order = (a.start, a.end).cmp(&(b.start, b.end));
    match order {
        std::cmp::Ordering::Equal => true,
        std::cmp::Ordering::Less => a.end > b.start,
        std::cmp::Ordering::Greater => b.end > a.start,
    }
}

/// `Chunk.isAdjacentTo`: the pointers *touch* exactly, end to start. Not "the same BGZF block",
/// which is what the comment about pointing directly at a block might suggest.
fn chunks_adjacent(a: &Chunk, b: &Chunk) -> bool {
    a.end == b.start || a.start == b.end
}

/// `Chunk.optimizeChunkList`.
///
/// Two rules, both htsjdk's: chunks ending at or before the linear index's minimum offset are
/// dropped, and chunks that overlap or exactly touch are coalesced into the running one.
pub fn optimize_chunk_list(chunks: &[Chunk], minimum_offset: u64) -> Vec<Chunk> {
    let mut sorted = chunks.to_vec();
    sorted.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));

    let mut result: Vec<Chunk> = Vec::new();
    for chunk in sorted {
        if chunk.end <= minimum_offset {
            continue; // linear index optimisation
        }
        match result.last_mut() {
            None => result.push(chunk),
            Some(last) => {
                if !chunks_overlap(last, &chunk) && !chunks_adjacent(last, &chunk) {
                    result.push(chunk);
                } else if chunk.end > last.end {
                    last.end = chunk.end;
                }
            }
        }
    }
    result
}

/// `GenomicIndexUtil.regionToBins`: the bins that may hold a record overlapping `[start, stop]`.
///
/// 1-based inclusive in, and the ends are decremented before shifting, which htsjdk's own
/// comment calls "suspicious" and then keeps.
pub fn region_to_bins(start_pos: i32, end_pos: i32) -> Option<Vec<usize>> {
    const MAX_POS: i32 = 0x1FFF_FFFF;
    let start = if start_pos <= 0 {
        0
    } else {
        (start_pos - 1) & MAX_POS
    };
    let end = if end_pos <= 0 {
        MAX_POS
    } else {
        (end_pos - 1) & MAX_POS
    };
    if start > end {
        return None;
    }
    let mut bins = vec![0usize];
    for (offset, shift) in [(1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)] {
        for k in (offset + (start >> shift))..=(offset + (end >> shift)) {
            bins.push(k as usize);
        }
    }
    Some(bins)
}

/// `LinearIndex.getMinimumOffset`: the earliest virtual offset any record covering `start_pos`
/// can have. Out of range means no constraint, which is `0`.
pub fn minimum_offset(linear_index: &[u64], start_pos: i32) -> u64 {
    const BAM_LIDX_SHIFT: i32 = 14;
    let start = if start_pos <= 0 { 0 } else { start_pos - 1 };
    let bin = (start >> BAM_LIDX_SHIFT) as usize;
    linear_index.get(bin).copied().unwrap_or(0)
}

/// `ReadsPathDataSource` over a single indexed BAM.
///
/// The whole compressed file is held in memory. That is right for the conformance corpus and
/// wrong for a genome, and it is a property of this slice rather than of the port: the seeking
/// below only needs a byte offset, so a file-backed reader replaces it without touching the
/// query semantics.
pub struct ReadsDataSource {
    compressed: Vec<u8>,
    header: SamHeader,
    /// The path the BAM was opened from, which a `RuntimeEOFException` names.
    path: String,
    /// One entry per reference sequence: its bins' chunks and its linear index.
    index: Vec<ReferenceIndex>,
    /// The virtual pointer no record before which is read: Disq's guessed first record for a
    /// Spark tool, and 0, which clips nothing, for everything else.
    split_start: u64,
}

struct ReferenceIndex {
    bins: Vec<(usize, Vec<Chunk>)>,
    linear: Vec<u64>,
}

impl ReadsDataSource {
    /// Open a BAM and its `.bai`, which must exist: GATK raises a `UserException` naming the
    /// `samtools index` command rather than falling back to a linear scan.
    pub fn open(bam: &Path, bai: &Path) -> Result<ReadsDataSource, ReadsError> {
        let compressed = std::fs::read(bam).map_err(|e| ReadsError::Io(e.to_string()))?;
        let decompressed = htsjdk_bgzf::read::decompress_all(&compressed)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        let header = BamReader::new(&decompressed)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?
            .header
            .text;

        let index = Self::read_reference_indexes(bai)?;
        Ok(ReadsDataSource {
            compressed,
            header,
            path: bam.display().to_string(),
            index,
            split_start: 0,
        })
    }

    fn read_reference_indexes(bai: &Path) -> Result<Vec<ReferenceIndex>, ReadsError> {
        // The `.bai` is parsed by the port's own reader rather than by a dependency. It is the
        // same file either way; what differs is who answers for it, and htsjdk-bam's `read_bai`
        // is measured by the `textual-index` suite, which parses six `.bai` files and reprints
        // each through `TextualBAMIndexWriter`'s format.
        let bai_bytes = std::fs::read(bai).map_err(|e| ReadsError::Io(e.to_string()))?;
        let parsed = htsjdk_bam::index::read_bai(&bai_bytes)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        Ok(parsed
            .references
            .iter()
            .map(|reference| ReferenceIndex {
                // The pseudo-bin is not a bin: `read_bai` keeps it as `metadata`, so the map here
                // carries exactly the bins a query walks.
                bins: reference
                    .bins
                    .iter()
                    .map(|bin| {
                        (
                            bin.bin_number as usize,
                            bin.chunks
                                .iter()
                                .map(|c| Chunk {
                                    start: c.start,
                                    end: c.end,
                                })
                                .collect(),
                        )
                    })
                    .collect(),
                linear: reference.linear_index.clone(),
            })
            .collect())
    }

    /// Open a BAM with no index, which is what `new ReadsPathDataSource(path)` does.
    ///
    /// The index is not a property of the file for the reference either: `ReadsPathDataSource`
    /// takes one when it is given one and refuses only when a query needs it. A queryname-sorted
    /// BAM cannot be indexed at all, so a tool that walks one whole, like
    /// `crate::transfer_read_tags`'s reference does, must be able to open it. Every interval query
    /// on a source opened this way finds no bins and returns nothing, which is the shape a caller
    /// with no intervals never reaches.
    pub fn open_unindexed(bam: &Path) -> Result<ReadsDataSource, ReadsError> {
        let compressed = std::fs::read(bam).map_err(|e| ReadsError::Io(e.to_string()))?;
        let decompressed = htsjdk_bgzf::read::decompress_all(&compressed)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        let header = BamReader::new(&decompressed)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?
            .header
            .text;
        Ok(ReadsDataSource {
            compressed,
            header,
            path: bam.display().to_string(),
            index: Vec::new(),
            split_start: 0,
        })
    }

    pub fn header(&self) -> &SamHeader {
        &self.header
    }

    /// Fill in the file a `RuntimeEOFException` names, which the cursor cannot know.
    fn name_the_file(&self, error: ReadsError) -> ReadsError {
        match error {
            ReadsError::PrematureEof {
                expected, received, ..
            } => ReadsError::PrematureEof {
                expected,
                received,
                file: self.path.clone(),
            },
            other => other,
        }
    }

    /// `BinningIndexContent.getChunksOverlapping` then `Chunk.optimizeChunkList`.
    fn chunks_overlapping(&self, interval: &QueryInterval) -> Vec<Chunk> {
        let Some(reference) = self.index.get(interval.reference_index as usize) else {
            return Vec::new();
        };
        let Some(wanted) = region_to_bins(interval.start, interval.end) else {
            return Vec::new();
        };
        let mut chunks = Vec::new();
        // htsjdk walks the bit set in ascending bin order, so the pre-sort order is that and not
        // the index's own bin order. optimizeChunkList sorts anyway, but the two differ before it.
        let mut wanted = wanted;
        wanted.sort_unstable();
        for bin_id in wanted {
            if let Some((_, bin)) = reference.bins.iter().find(|(id, _)| *id == bin_id) {
                chunks.extend(bin.iter().copied());
            }
        }
        if chunks.is_empty() {
            return Vec::new();
        }
        optimize_chunk_list(&chunks, minimum_offset(&reference.linear, interval.start))
    }

    /// `ReadsPathDataSource.query`: the records overlapping the intervals, in file order.
    ///
    /// The intervals go through `convertSimpleIntervalToQueryInterval` and
    /// `optimizeIntervals` first, exactly as `SamReaderQueryingIterator` does, so abutting
    /// intervals merge and a read is returned once rather than once per interval.
    pub fn query(&self, intervals: &[SimpleInterval]) -> Result<Vec<BamRecord>, ReadsError> {
        self.query_within(intervals, self.split_start, u64::MAX)
    }

    /// `query` over what `span.removeContentsBefore(splitSpan)` and `removeContentsAfter` leave of
    /// the index's chunks, with a filter of its own: each Disq partition opens a new iterator.
    pub fn query_within(
        &self,
        intervals: &[SimpleInterval],
        from: u64,
        to: u64,
    ) -> Result<Vec<BamRecord>, ReadsError> {
        let converted: Result<Vec<_>, _> = intervals
            .iter()
            .map(|i| convert_simple_interval_to_query_interval(i, &self.header))
            .collect();
        let optimized = optimize_intervals(&converted?);
        if optimized.is_empty() {
            return Ok(Vec::new());
        }

        // BAMFileReader.getFileSpan: one span per interval, merged with minimumOffset 0.
        let mut spans = Vec::new();
        for interval in &optimized {
            spans.extend(self.chunks_overlapping(interval));
        }
        let span = optimize_chunk_list(&spans, 0);

        let mut filter = IntervalFilter::new(&optimized, false);
        let mut kept = Vec::new();
        'chunks: for chunk in span {
            // `removeContentsBefore(splitSpan)`: Disq's split clips the index's chunks.
            if chunk.end <= from || chunk.start >= to {
                continue;
            }
            let start = chunk.start.max(from);
            let end = chunk.end.min(to);
            let mut cursor = BlockCursor::seek(&self.compressed, start)?;
            while cursor.virtual_pos() < end {
                let Some(record) = cursor.next_record().map_err(|e| self.name_the_file(e))? else {
                    break;
                };
                match filter.compare_to_filter(&record) {
                    FilterState::MatchesFilter => kept.push(record),
                    FilterState::ContinueIteration => {}
                    FilterState::StopIteration => break 'chunks,
                }
            }
        }
        Ok(kept)
    }

    /// `ReadsPathDataSource.queryUnmapped`: the unplaced unmapped reads at the tail of the file.
    ///
    /// htsjdk seeks to the last linear-index offset and then skips forward to the first record
    /// with no reference, which is why a file whose unmapped reads are not at the end returns
    /// nothing here.
    pub fn query_unmapped(&self) -> Result<Vec<BamRecord>, ReadsError> {
        let start = match self.start_of_last_linear_bin() {
            Some(pointer) => pointer,
            // No mapped reads in the file, so start at the first record.
            None => self.first_record_pointer()?,
        };
        let mut cursor = BlockCursor::seek(&self.compressed, start)?;
        let mut kept = Vec::new();
        let mut reached = false;
        while let Some(record) = cursor.next_record().map_err(|e| self.name_the_file(e))? {
            if !reached {
                if record.reference_index != -1 {
                    continue;
                }
                reached = true;
            }
            kept.push(record);
        }
        Ok(kept)
    }

    /// `AbstractBAMFileIndex.getStartOfLastLinearBin`: the **last** entry of the last reference
    /// that has a linear index at all, which is not the same as the largest entry.
    ///
    /// The loop overwrites rather than maximising, so a reference whose linear index is shorter
    /// than an earlier one's still wins by being later. htsjdk's own comment explains the shape:
    /// no read may align to the last sequences in the dictionary.
    fn start_of_last_linear_bin(&self) -> Option<u64> {
        let mut last = None;
        for reference in &self.index {
            if let Some(&offset) = reference.linear.last() {
                last = Some(offset);
            }
        }
        last
    }

    /// The virtual pointer of the first record, which htsjdk keeps as `mFirstRecordPointer`.
    ///
    /// The header can span BGZF blocks, so this walks blocks counting decompressed bytes rather
    /// than assuming the records start in the first one.
    fn first_record_pointer(&self) -> Result<u64, ReadsError> {
        let decompressed = htsjdk_bgzf::read::decompress_all(&self.compressed)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        let header_end = BamReader::new(&decompressed)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?
            .records_start();

        let mut consumed = 0usize;
        let mut address = 0u64;
        while (address as usize) < self.compressed.len() {
            let mut reader = BgzfReader::new(&self.compressed[address as usize..]);
            let Some(block) = reader
                .next_block()
                .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?
            else {
                break;
            };
            if consumed + block.data.len() > header_end {
                let offset = (header_end - consumed) as u32;
                return vfp::make_file_pointer(address, offset)
                    .map_err(|e| ReadsError::Malformed(format!("{e:?}")));
            }
            consumed += block.data.len();
            address += block.block_compressed_size as u64;
        }
        Ok(0)
    }

    /// A full traversal, which is what a walker with no `-L` does.
    pub fn iter_all(&self) -> Result<Vec<BamRecord>, ReadsError> {
        if self.split_start != 0 {
            let mut cursor = BlockCursor::seek(&self.compressed, self.split_start)?;
            let mut records = Vec::new();
            while let Some(record) = cursor.next_record().map_err(|e| self.name_the_file(e))? {
                records.push(record);
            }
            return Ok(records);
        }
        let decompressed = htsjdk_bgzf::read::decompress_all(&self.compressed)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        let reader =
            BamReader::new(&decompressed).map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        reader
            .map(|r| r.map_err(|e| ReadsError::Malformed(format!("{e:?}"))))
            .collect()
    }
}

impl ReadsDataSource {
    /// A BAM and, under intervals, its `.bai`, inflating the header's blocks and no others: the
    /// whole file inflated at once is what a split exists to avoid.
    ///
    /// shortcut: the compressed file is still held whole, as [`ReadsDataSource`] says.
    pub fn open_for_splits(bam: &Path, bai: Option<&Path>) -> Result<ReadsDataSource, ReadsError> {
        let compressed = std::fs::read(bam).map_err(|e| ReadsError::Io(e.to_string()))?;
        let mut inflated = Vec::new();
        let mut address = 0usize;
        let header = loop {
            match BamReader::new(&inflated) {
                Ok(reader) => break reader.header.text,
                Err(htsjdk_bam::reader::ReadError::Truncated { .. }) => {}
                Err(error) => return Err(ReadsError::Malformed(format!("{error:?}"))),
            }
            let mut reader = BgzfReader::new(&compressed[address.min(compressed.len())..]);
            let Some(block) = reader
                .next_block()
                .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?
            else {
                let error = BamReader::new(&inflated).err();
                return Err(ReadsError::Malformed(format!("{error:?}")));
            };
            address += block.block_compressed_size;
            inflated.extend_from_slice(&block.data);
        };
        let index = match bai {
            Some(bai) => Self::read_reference_indexes(bai)?,
            None => Vec::new(),
        };
        Ok(ReadsDataSource {
            compressed,
            header,
            path: bam.display().to_string(),
            index,
            split_start: 0,
        })
    }

    /// `BamSource.getPathChunks`: a chunk for each split, or none for a split no record starts in,
    /// which is still a partition, with no reads.
    ///
    /// With a `.sbi` beside the BAM the index says where a split's records run. Without one, the
    /// split is its blocks, the start is the first position `BamRecordGuesser` accepts in them, and
    /// the end is the split's own, so a record belongs to the split its block starts in, however
    /// far it runs. A record the guesser rejects at the head of a split is in no partition at all.
    pub fn disq_partitions(
        &self,
        split_size: i32,
        use_nio: bool,
    ) -> Result<Vec<Option<Chunk>>, DisqSplitError> {
        let splits = disq_path_splits(self.compressed.len() as i64, split_size, use_nio);
        let sbi_path = format!("{}.sbi", self.path);
        if Path::new(&sbi_path).exists() {
            let bytes = std::fs::read(&sbi_path).map_err(|e| ReadsError::Io(e.to_string()))?;
            let offsets = load_sbi(&bytes)?;
            return Ok(splits
                .into_iter()
                .map(|(start, end)| sbi_chunk(&offsets, start, end))
                .collect());
        }
        let blocks = block_table(&self.compressed)?;
        let total = blocks
            .last()
            .map_or(0, |block| block.linear as usize + block.size as usize);
        let mut stream = InflatedStream {
            data: &self.compressed,
            blocks: &blocks,
            total,
            window: Vec::new(),
            window_start: 0,
            next_block: 0,
        };
        let lengths: Vec<i32> = self.header.sequences.iter().map(|s| s.length).collect();
        let mut partitions = Vec::with_capacity(splits.len());
        for (start, end) in splits {
            partitions.push(first_read_in_split(
                &self.compressed,
                &blocks,
                &mut stream,
                &lengths,
                start,
                end,
            )?);
        }
        Ok(partitions)
    }

    /// `BAMFileReader.getIterator(chunks)`: a record belongs to the chunk it starts in.
    pub fn read_chunk(&self, chunk: &Chunk) -> Result<Vec<BamRecord>, ReadsError> {
        let mut cursor = BlockCursor::seek(&self.compressed, chunk.start)?;
        let mut records = Vec::new();
        while cursor.virtual_pos() < chunk.end {
            let Some(record) = cursor.next_record().map_err(|e| self.name_the_file(e))? else {
                break;
            };
            records.push(record);
        }
        Ok(records)
    }

    /// Disq's `BamSource.getFirstReadInPartition` with no `.sbi`: the reads start where its
    /// `BamRecordGuesser` first accepts a record, not where the header ends.
    ///
    /// The guess scans every uncompressed byte from the start of the first block, the header's
    /// own, and accepts a position from which ten records chain, or the file ends after at least
    /// one. A mapped record with no cigar is never a start, so every record before one that lies
    /// within ten of it is never read: the reads are lost, not refused. No accepted position in
    /// the first ten million is a partition with no reads.
    ///
    /// shortcut: one partition, which is what a file under the split size gets; a larger file is
    /// split and each split guesses its own start, which [`Self::disq_partitions`] does.
    pub fn with_disq_split_start(mut self) -> Result<ReadsDataSource, ReadsError> {
        let mut blocks = Vec::new();
        let mut address = 0usize;
        while address < self.compressed.len() {
            let mut reader = BgzfReader::new(&self.compressed[address..]);
            let Some(block) = reader
                .next_block()
                .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?
            else {
                break;
            };
            let size = block.block_compressed_size;
            blocks.push((address as u64, block.data));
            address += size;
        }
        let flat: Vec<u8> = blocks
            .iter()
            .flat_map(|(_, data)| data.iter().copied())
            .collect();
        let lengths: Vec<i32> = self.header.sequences.iter().map(|s| s.length).collect();
        let mut linear = 0usize;
        let mut index = 0usize;
        'blocks: for (address, data) in &blocks {
            for up in 0..data.len() {
                index += 1;
                if index > 10_000_000 {
                    break 'blocks;
                }
                if guess_chain(&mut flat.as_slice(), &lengths, linear + up)? {
                    self.split_start = vfp::make_file_pointer(*address, up as u32)
                        .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
                    return Ok(self);
                }
            }
            linear += data.len();
        }
        // Nothing accepted: past every block, so nothing is read.
        self.split_start = vfp::make_file_pointer(self.compressed.len() as u64, 0)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        Ok(self)
    }
}

/// What the guesser reads: the decompressed stream, addressed from its first byte.
trait GuessedStream {
    fn total(&self) -> usize;
    fn bytes(&mut self, at: usize, length: usize) -> Result<&[u8], ReadsError>;
}

impl GuessedStream for &[u8] {
    fn total(&self) -> usize {
        self.len()
    }
    fn bytes(&mut self, at: usize, length: usize) -> Result<&[u8], ReadsError> {
        Ok(&self[at..at + length])
    }
}

/// `BamRecordGuesser.checkRecordStart`: ten records chaining from `at`, or the end of the
/// stream after at least one.
fn guess_chain(
    flat: &mut impl GuessedStream,
    lengths: &[i32],
    mut at: usize,
) -> Result<bool, ReadsError> {
    for successful in 0..10 {
        match guess_record(flat, lengths, at)? {
            Guess::Start(next) => at = next,
            Guess::NoStart => return Ok(false),
            Guess::Eof => return Ok(successful > 0),
        }
    }
    Ok(true)
}

enum Guess {
    Start(usize),
    NoStart,
    Eof,
}

/// `BamRecordGuesser.checkRecordStartInternal`, over the decompressed stream.
fn guess_record(
    flat: &mut impl GuessedStream,
    lengths: &[i32],
    at: usize,
) -> Result<Guess, ReadsError> {
    let total = flat.total();
    if at + 36 > total {
        return Ok(Guess::Eof);
    }
    let fixed: [u8; 36] = flat.bytes(at, 36)?.try_into().expect("36 bytes");
    let int = |offset: usize| {
        i32::from_le_bytes([
            fixed[offset],
            fixed[offset + 1],
            fixed[offset + 2],
            fixed[offset + 3],
        ])
    };
    let count = lengths.len() as i32;
    let locus_ok = |id: i32, pos: i32| {
        !(id < -1 || id >= count || pos < -1 || (id >= 0 && pos > lengths[id as usize]))
    };
    let remaining = int(0);
    if !locus_ok(int(4), int(8)) || !locus_ok(int(24), int(28)) {
        return Ok(Guess::NoStart);
    }
    let name_length = (int(12) & 0xff) as usize;
    if name_length < 2 {
        return Ok(Guess::NoStart);
    }
    let flags = (int(16) as u32 >> 16) as u16;
    let cigar_ops = (int(16) & 0xffff) as usize;
    let l_seq = int(20);
    let seq_length = l_seq.wrapping_add(l_seq.wrapping_add(1) / 2);
    let mapped = flags & 0x4 == 0;
    let primary = flags & 0x100 == 0 && flags & 0x800 == 0;
    if mapped && ((seq_length == 0 && primary) || cigar_ops == 0) {
        return Ok(Guess::NoStart);
    }
    let name = at + 36;
    if name + name_length > total {
        return Ok(Guess::Eof);
    }
    let bytes = flat.bytes(name, name_length)?;
    if bytes[name_length - 1] != 0 {
        return Ok(Guess::NoStart);
    }
    let valid = |b: u8| (b'!'..=b'?').contains(&b) || (b'A'..=b'~').contains(&b);
    if !bytes[..name_length - 1].iter().all(|&b| valid(b)) {
        return Ok(Guess::NoStart);
    }
    let mut cigar = name + name_length;
    for _ in 0..cigar_ops {
        if cigar + 4 > total {
            return Ok(Guess::Eof);
        }
        let held = flat.bytes(cigar, 4)?;
        let op = i32::from_le_bytes([held[0], held[1], held[2], held[3]]);
        if op == -1 {
            return Ok(Guess::Eof);
        }
        if op & 0xf > 8 {
            return Ok(Guess::NoStart);
        }
        cigar += 4;
    }
    let zero_min = (32 + name_length as i64 + 4 * cigar_ops as i64) + seq_length as i64;
    if (remaining as i64) < zero_min {
        return Ok(Guess::NoStart);
    }
    let next = at + 4 + remaining as usize;
    if next > total {
        return Ok(Guess::Eof);
    }
    Ok(Guess::Start(next))
}

/// `FileInputFormat.SPLIT_SLOP`: a last split may run a tenth over.
const SPLIT_SLOP: f64 = 1.1;
/// `fs.local.block.size`, which caps a split on a local file whatever was asked for.
const LOCAL_BLOCK_SIZE: i64 = 32 * 1024 * 1024;
/// `PathSplitSource`'s default when `--use-nio` is given no size.
const NIO_SPLIT_SIZE: i64 = 128 * 1024 * 1024;
/// `BamSource.MAX_READ_SIZE`.
const MAX_READ_SIZE: usize = 10_000_000;

/// `PathSplitSource.getPathSplits` over a file of `length` bytes, as (start, end).
///
/// `splitSize` reached Disq through `(int) splitSize`, and only a positive one is set on the
/// configuration, so `2^32` asks for the default and `--bam-partition-size` above the block size
/// changes nothing on a local file.
pub fn disq_path_splits(length: i64, split_size: i32, use_nio: bool) -> Vec<(i64, i64)> {
    let mut splits = Vec::new();
    if use_nio {
        let actual = if split_size <= 0 {
            NIO_SPLIT_SIZE
        } else {
            i64::from(split_size)
        };
        let count = (length as f64 / actual as f64).ceil() as i64;
        for index in 0..count {
            let start = index * actual;
            splits.push((start, (start + actual).min(length)));
        }
        return splits;
    }
    let maximum = if split_size > 0 {
        i64::from(split_size)
    } else {
        i64::MAX
    };
    let size = maximum.clamp(1, LOCAL_BLOCK_SIZE);
    let mut remaining = length;
    while remaining as f64 / size as f64 > SPLIT_SLOP {
        splits.push((length - remaining, length - remaining + size));
        remaining -= size;
    }
    if remaining != 0 {
        splits.push((length - remaining, length));
    }
    splits
}

/// What Disq's split computation can stop on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DisqSplitError {
    Reads(ReadsError),
    /// `BgzfBlockGuesser` accepted bytes no block starts at, which its reader would then inflate
    /// as a block: the port's own refusal, carrying the address.
    NotABlock(i64),
    /// `SBIIndex.load` said no.
    Sbi(String),
}

impl From<ReadsError> for DisqSplitError {
    fn from(error: ReadsError) -> DisqSplitError {
        DisqSplitError::Reads(error)
    }
}

/// One real block of the file, and where its bytes fall in the decompressed stream.
#[derive(Debug, Clone, Copy)]
struct BlockEntry {
    address: u64,
    linear: u64,
    size: u32,
}

/// The file's blocks, read from their headers and trailers with nothing inflated.
fn block_table(data: &[u8]) -> Result<Vec<BlockEntry>, ReadsError> {
    let mut blocks = Vec::new();
    let mut address = 0usize;
    let mut linear = 0u64;
    let malformed = |address: usize| {
        ReadsError::Malformed(format!("no BGZF block where one was due at {address}"))
    };
    while address < data.len() {
        let head = data.get(address..address + 12).ok_or(malformed(address))?;
        if head[..4] != [0x1f, 0x8b, 8, 4] {
            return Err(malformed(address));
        }
        let extra = usize::from(u16::from_le_bytes([head[10], head[11]]));
        let fields = data
            .get(address + 12..address + 12 + extra)
            .ok_or(malformed(address))?;
        let mut at = 0;
        let mut size = None;
        while at + 4 <= fields.len() {
            let length = usize::from(u16::from_le_bytes([fields[at + 2], fields[at + 3]]));
            if fields[at..at + 2] == *b"BC" && length == 2 && at + 6 <= fields.len() {
                size = Some(usize::from(u16::from_le_bytes([
                    fields[at + 4],
                    fields[at + 5],
                ])));
            }
            at += 4 + length;
        }
        let compressed = size.ok_or(malformed(address))? + 1;
        let trailer = data
            .get(address + compressed - 4..address + compressed)
            .ok_or(malformed(address))?;
        let inflated = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
        blocks.push(BlockEntry {
            address: address as u64,
            linear,
            size: inflated,
        });
        linear += u64::from(inflated);
        address += compressed;
    }
    Ok(blocks)
}

/// How far behind the guesser's position its window is allowed to fall.
const WINDOW_KEPT: usize = 16 * 1024 * 1024;

/// The decompressed stream, inflated as the guesser walks it, keeping a window behind.
struct InflatedStream<'a> {
    data: &'a [u8],
    blocks: &'a [BlockEntry],
    total: usize,
    window: Vec<u8>,
    window_start: usize,
    next_block: usize,
}

impl GuessedStream for InflatedStream<'_> {
    fn total(&self) -> usize {
        self.total
    }

    fn bytes(&mut self, at: usize, length: usize) -> Result<&[u8], ReadsError> {
        if at < self.window_start || at > self.window_start + self.window.len() + WINDOW_KEPT {
            let block = self
                .blocks
                .partition_point(|block| block.linear as usize + block.size as usize <= at);
            self.next_block = block;
            self.window.clear();
            self.window_start = self.blocks.get(block).map_or(at, |b| b.linear as usize);
        } else if at > self.window_start + 2 * WINDOW_KEPT {
            let dropped = at - WINDOW_KEPT - self.window_start;
            self.window.drain(..dropped);
            self.window_start += dropped;
        }
        while self.window_start + self.window.len() < at + length {
            let Some(block) = self.blocks.get(self.next_block) else {
                break;
            };
            let mut reader = BgzfReader::new(&self.data[block.address as usize..]);
            let inflated = reader
                .next_block()
                .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?
                .ok_or_else(|| ReadsError::Malformed("a block that is not there".to_string()))?;
            self.window.extend_from_slice(&inflated.data);
            self.next_block += 1;
        }
        let from = at - self.window_start;
        Ok(&self.window[from..from + length])
    }
}

/// `BgzfBlockGuesser.guessNextBGZFPos`: (pos, cSize, uSize), or null at the end or on a short read.
fn guess_next_block(data: &[u8], mut p: i64, end: i64) -> Option<(i64, i32, i32)> {
    const MAGIC: u32 = 0x0408_8b1f;
    const MAGIC_SUB: u32 = 0x0002_4342;
    let at = |p: i64, n: usize| -> Option<&[u8]> {
        usize::try_from(p)
            .ok()
            .and_then(|p| data.get(p..p.checked_add(n)?))
    };
    let int = |p: i64| at(p, 4).map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]));
    let short = |p: i64| at(p, 2).map(|b| i64::from(u16::from_le_bytes([b[0], b[1]])));
    loop {
        loop {
            let n = int(p)?;
            if n == MAGIC {
                break;
            }
            if n >> 8 == MAGIC & 0x00ff_ffff {
                p += 1;
            } else if n >> 16 == MAGIC & 0x0000_ffff {
                p += 2;
            } else {
                p += 3;
            }
            if p >= end {
                return None;
            }
        }
        let start = p;
        p += 10;
        let xlen = short(p)?;
        p += 2;
        let sub_end = p + xlen;
        while p < sub_end {
            let found = int(p)?;
            if found != MAGIC_SUB {
                p += 4 + i64::from(found >> 16);
                continue;
            }
            let bsize = short(p + 4)?;
            p += 6;
            while p < sub_end {
                p += 4 + i64::from(int(p)? >> 16);
            }
            if p != sub_end {
                break;
            }
            p += bsize - xlen - 19 + 4;
            let size = int(p)?;
            return Some((start, (p + 4 - start) as i32, size as i32));
        }
        p = start + 4;
    }
}

/// `SBIIndex.load`: the virtual offsets.
fn load_sbi(bytes: &[u8]) -> Result<Vec<u64>, DisqSplitError> {
    let eof = || DisqSplitError::Sbi("Premature EOF".to_string());
    let long = |at: usize| -> Result<u64, DisqSplitError> {
        let held = bytes.get(at..at + 8).ok_or_else(eof)?;
        Ok(u64::from_le_bytes(held.try_into().expect("eight bytes")))
    };
    if bytes.get(..4) != Some(b"SBI\x01".as_slice()) {
        return Err(DisqSplitError::Sbi(
            "Invalid file header in SBI".to_string(),
        ));
    }
    // The magic, the length, the MD5, the UUID, the record count, the granularity, then how many.
    let count = long(4 + 8 + 16 + 16 + 8 + 8)? as i64;
    let mut offsets = Vec::new();
    let mut previous = -1i64;
    for index in 0..count.max(0) as usize {
        let current = long(68 + 8 * index)? as i64;
        if previous > current {
            return Err(DisqSplitError::Sbi(format!(
                "Invalid SBI; offsets not in order: {previous:#x} > {current:#x}"
            )));
        }
        offsets.push(current as u64);
        previous = current;
    }
    if offsets.is_empty() {
        return Err(DisqSplitError::Sbi(
            "Invalid SBI format: should contain at least one offset".to_string(),
        ));
    }
    Ok(offsets)
}

/// `SBIIndex.getChunk`.
fn sbi_chunk(offsets: &[u64], start: i64, end: i64) -> Option<Chunk> {
    let last = *offsets.last().expect("at least one offset");
    let maximum = vfp::block_address(last) as i64;
    let ceiling = |address: i64| {
        let pointer = (address.min(maximum) as u64) << 16;
        offsets
            .get(offsets.partition_point(|offset| *offset < pointer))
            .copied()
    };
    let (from, to) = (ceiling(start)?, ceiling(end)?);
    (from != to).then_some(Chunk {
        start: from,
        end: to,
    })
}

/// `getFirstReadInPartition`, over the blocks `BgzfBlockSource` guesses in one split.
fn first_read_in_split(
    data: &[u8],
    blocks: &[BlockEntry],
    stream: &mut InflatedStream<'_>,
    lengths: &[i32],
    split_start: i64,
    split_end: i64,
) -> Result<Option<Chunk>, DisqSplitError> {
    let mut index = 0usize;
    let mut start = split_start;
    while start < split_end {
        let Some((position, compressed, size)) = guess_next_block(data, start, split_end) else {
            break;
        };
        start = position + i64::from(compressed);
        let entry = blocks
            .binary_search_by_key(&(position as u64), |block| block.address)
            .map_err(|_| DisqSplitError::NotABlock(position))?;
        let linear = blocks[entry].linear as usize;
        for up in 0..size.max(0) as usize {
            index += 1;
            if index > MAX_READ_SIZE {
                return Ok(None);
            }
            if guess_chain(stream, lengths, linear + up)? {
                let pointer = |address: i64, offset: u32| {
                    vfp::make_file_pointer(address as u64, offset)
                        .map_err(|e| ReadsError::Malformed(format!("{e:?}")))
                };
                return Ok(Some(Chunk {
                    start: pointer(position, up as u32)?,
                    end: pointer(split_end, 0)?,
                }));
            }
        }
    }
    Ok(None)
}

/// Reads records from an arbitrary virtual file pointer, crossing BGZF block boundaries.
///
/// `BAMFileIndexIterator` checks the position *before* pulling a record, so a record that starts
/// inside a chunk is read whole even when it runs past the chunk's end.
struct BlockCursor<'a> {
    data: &'a [u8],
    /// Byte offset of the current block in the compressed file.
    block_address: u64,
    /// Byte offset of the block after the current one.
    next_address: u64,
    block: Vec<u8>,
    offset: usize,
}

impl<'a> BlockCursor<'a> {
    fn seek(data: &'a [u8], pointer: u64) -> Result<BlockCursor<'a>, ReadsError> {
        let address = vfp::block_address(pointer);
        let mut cursor = BlockCursor {
            data,
            block_address: address,
            next_address: address,
            block: Vec::new(),
            offset: 0,
        };
        cursor.load_block(address)?;
        cursor.offset = vfp::block_offset(pointer) as usize;
        Ok(cursor)
    }

    fn load_block(&mut self, address: u64) -> Result<bool, ReadsError> {
        if address as usize >= self.data.len() {
            return Ok(false);
        }
        let mut reader = BgzfReader::new(&self.data[address as usize..]);
        match reader
            .next_block()
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?
        {
            None => Ok(false),
            Some(block) => {
                self.block_address = address;
                self.next_address = address + block.block_compressed_size as u64;
                self.block = block.data;
                self.offset = 0;
                Ok(true)
            }
        }
    }

    /// The virtual pointer of the next byte, matching `getFilePointer`: a pointer never dangles
    /// past the end of a block, it becomes the next block's zero offset.
    fn virtual_pos(&self) -> u64 {
        if self.offset >= self.block.len() {
            vfp::make_file_pointer(self.next_address, 0).unwrap_or(0)
        } else {
            vfp::make_file_pointer(self.block_address, self.offset as u32).unwrap_or(0)
        }
    }

    /// Reads up to `n` bytes, crossing into following blocks as needed.
    ///
    /// A short result is the stream ending, and the caller decides what that means: the length
    /// word ending short is the end of the records, and a record BODY ending short is
    /// `RuntimeEOFException`. Both need the count, which is why this returns what it read rather
    /// than `None`.
    fn read_up_to(&mut self, n: usize) -> Result<Vec<u8>, ReadsError> {
        let mut out = Vec::with_capacity(n.min(self.data.len()));
        while out.len() < n {
            if self.offset >= self.block.len() {
                let next = self.next_address;
                if !self.load_block(next)? {
                    return Ok(out);
                }
                // An empty block (the BGZF terminator) is not the end of the stream by itself.
                if self.block.is_empty() && self.next_address as usize >= self.data.len() {
                    return Ok(out);
                }
                continue;
            }
            let take = (n - out.len()).min(self.block.len() - self.offset);
            out.extend_from_slice(&self.block[self.offset..self.offset + take]);
            self.offset += take;
        }
        Ok(out)
    }

    fn next_record(&mut self) -> Result<Option<BamRecord>, ReadsError> {
        let size = self.read_up_to(4)?;
        if size.len() < 4 {
            // `BAMRecordCodec.decode` catches the EOF on the LENGTH and answers null, which is
            // how a traversal ends rather than how it fails.
            return Ok(None);
        }
        let block_size = i32::from_le_bytes([size[0], size[1], size[2], size[3]]);
        // A length under the fixed block size, negative ones included, never reaches a read at
        // all: `BAMRecordCodec.decode` compares it first and raises `SAMFormatException`, which is
        // what `BamRecord::decode` answers below with the bytes it was given.
        let body = if (block_size as usize) < FIXED_BLOCK_SIZE || block_size < 0 {
            self.read_up_to(block_size.max(0) as usize)?
        } else {
            // The body is read in TWO calls, and which one runs out decides the numbers the
            // exception carries: the fixed block first, then the variable part after it. So a
            // stream with a hundred bytes left and a length of a gigabyte reports the SECOND
            // call, expecting the length less the fixed block and having received what was left
            // after it.
            let fixed = self.read_up_to(FIXED_BLOCK_SIZE)?;
            if fixed.len() < FIXED_BLOCK_SIZE {
                return Err(ReadsError::PrematureEof {
                    expected: FIXED_BLOCK_SIZE as i64,
                    received: fixed.len(),
                    file: String::new(),
                });
            }
            let rest_length = block_size as usize - FIXED_BLOCK_SIZE;
            let rest = self.read_up_to(rest_length)?;
            if rest.len() < rest_length {
                // The EOF here is not caught: `BinaryCodec.readBytes` raises, and the file it
                // names is filled in by the caller, which is the one that knows the path.
                return Err(ReadsError::PrematureEof {
                    expected: rest_length as i64,
                    received: rest.len(),
                    file: String::new(),
                });
            }
            let mut body = fixed;
            body.extend_from_slice(&rest);
            body
        };
        let mut bytes = size;
        bytes.extend_from_slice(&body);
        match BamRecord::decode(&bytes) {
            Ok(Some((record, _))) => Ok(Some(record)),
            Ok(None) => Ok(None),
            Err(e) => Err(ReadsError::Malformed(format!("{e:?}"))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(reference_index: i32, start: i32, end: i32) -> QueryInterval {
        QueryInterval {
            reference_index,
            start,
            end,
        }
    }

    #[test]
    fn a_local_split_stops_at_the_block_size_and_the_slop_keeps_a_tenth() {
        let mib = 1024 * 1024;
        // 35 MiB over 32 is under the slop, and 36 is over it; asking for more changes nothing.
        assert_eq!(disq_path_splits(35 * mib, 0, false), vec![(0, 35 * mib)]);
        assert_eq!(
            disq_path_splits(36 * mib, 64 * mib as i32, false),
            vec![(0, 32 * mib), (32 * mib, 36 * mib)]
        );
        assert_eq!(disq_path_splits(0, 0, false), Vec::new());
        // NIO has no slop and no block size.
        assert_eq!(
            disq_path_splits(11, 5, true),
            vec![(0, 5), (5, 10), (10, 11)]
        );
    }

    #[test]
    fn an_sbi_holds_its_count_then_its_offsets() {
        let mut bytes = b"SBI\x01".to_vec();
        bytes.extend_from_slice(&[0; 56]);
        bytes.extend_from_slice(&2i64.to_le_bytes());
        bytes.extend_from_slice(&7i64.to_le_bytes());
        bytes.extend_from_slice(&9i64.to_le_bytes());
        assert_eq!(load_sbi(&bytes).unwrap(), vec![7, 9]);
    }

    #[test]
    fn an_sbi_chunk_is_the_offsets_between_its_ceilings() {
        let offsets = [1 << 16, 100 << 16, (100 << 16) + 7, 300 << 16];
        let chunk = sbi_chunk(&offsets, 0, 200).unwrap();
        assert_eq!((chunk.start, chunk.end), (1 << 16, 300 << 16));
        // No offset starts between them, so the split is empty.
        assert_eq!(sbi_chunk(&offsets, 101, 200), None);
        // An end past the last offset's block stops at that offset.
        assert_eq!(sbi_chunk(&offsets, 301, 400), None);
    }

    #[test]
    fn abutting_intervals_merge_so_a_read_is_returned_once() {
        let merged = optimize_intervals(&[interval(0, 201, 300), interval(0, 100, 200)]);
        assert_eq!(merged, vec![interval(0, 100, 300)]);
    }

    #[test]
    fn a_zero_end_swallows_what_follows_it_on_the_same_contig() {
        let merged = optimize_intervals(&[interval(0, 100, 0), interval(0, 150, 200)]);
        assert_eq!(merged, vec![interval(0, 100, 0)]);
    }

    #[test]
    fn intervals_on_different_contigs_never_merge() {
        let merged = optimize_intervals(&[interval(1, 1, 10), interval(0, 1, 10)]);
        assert_eq!(merged, vec![interval(0, 1, 10), interval(1, 1, 10)]);
    }

    #[test]
    fn region_to_bins_always_includes_bin_zero() {
        let bins = region_to_bins(1, 100).unwrap();
        assert!(bins.contains(&0));
        // 1-based [1,100] falls in the first 16 kb window, bin 4681.
        assert!(bins.contains(&4681));
        assert!(region_to_bins(100, 1).is_none());
    }

    #[test]
    fn the_linear_index_bounds_a_query_and_missing_entries_do_not() {
        let linear = vec![100u64, 200, 300];
        assert_eq!(minimum_offset(&linear, 1), 100);
        assert_eq!(minimum_offset(&linear, 16_385), 200);
        assert_eq!(minimum_offset(&linear, 1_000_000), 0);
    }
}
