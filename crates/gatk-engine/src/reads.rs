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
//! The `.bai` bytes are parsed by [`noodles_bam::bai`], under the rule in
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
    /// One entry per reference sequence: its bins' chunks and its linear index.
    index: Vec<ReferenceIndex>,
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

        let bai = noodles_bam::bai::fs::read(bai).map_err(|e| ReadsError::Io(e.to_string()))?;
        let index = bai
            .reference_sequences()
            .iter()
            .map(|reference| ReferenceIndex {
                bins: reference
                    .bins()
                    .iter()
                    .map(|(id, bin)| {
                        (
                            *id,
                            bin.chunks()
                                .iter()
                                .map(|c| Chunk {
                                    start: u64::from(c.start()),
                                    end: u64::from(c.end()),
                                })
                                .collect(),
                        )
                    })
                    .collect(),
                linear: reference.index().iter().map(|p| u64::from(*p)).collect(),
            })
            .collect();

        Ok(ReadsDataSource {
            compressed,
            header,
            index,
        })
    }

    pub fn header(&self) -> &SamHeader {
        &self.header
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
            let mut cursor = BlockCursor::seek(&self.compressed, chunk.start)?;
            while cursor.virtual_pos() < chunk.end {
                let Some(record) = cursor.next_record()? else {
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
        while let Some(record) = cursor.next_record()? {
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
        let decompressed = htsjdk_bgzf::read::decompress_all(&self.compressed)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        let reader =
            BamReader::new(&decompressed).map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        reader
            .map(|r| r.map_err(|e| ReadsError::Malformed(format!("{e:?}"))))
            .collect()
    }
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

    /// Reads `n` bytes, crossing into following blocks as needed.
    fn read_exact(&mut self, n: usize) -> Result<Option<Vec<u8>>, ReadsError> {
        let mut out = Vec::with_capacity(n);
        while out.len() < n {
            if self.offset >= self.block.len() {
                let next = self.next_address;
                if !self.load_block(next)? {
                    return Ok(None);
                }
                // An empty block (the BGZF terminator) is not the end of the stream by itself.
                if self.block.is_empty() && self.next_address as usize >= self.data.len() {
                    return Ok(None);
                }
                continue;
            }
            let take = (n - out.len()).min(self.block.len() - self.offset);
            out.extend_from_slice(&self.block[self.offset..self.offset + take]);
            self.offset += take;
        }
        Ok(Some(out))
    }

    fn next_record(&mut self) -> Result<Option<BamRecord>, ReadsError> {
        let Some(size) = self.read_exact(4)? else {
            return Ok(None);
        };
        let block_size = i32::from_le_bytes([size[0], size[1], size[2], size[3]]);
        if block_size < 0 {
            return Err(ReadsError::Malformed(format!(
                "negative record length {block_size}"
            )));
        }
        let Some(body) = self.read_exact(block_size as usize)? else {
            return Ok(None);
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
