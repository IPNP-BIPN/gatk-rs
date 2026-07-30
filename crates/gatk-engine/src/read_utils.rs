//! Ported from `org.broadinstitute.hellbender.utils.read.ReadUtils` and the part of `BaseUtils`
//! it reaches (GATK 4.6.2.0).
//!
//! This is the arithmetic every walker, annotation and clipping operation stands on: a reference
//! position in, a read offset out. An off-by-one here does not fail loudly, it reads the
//! neighbouring base, and every number computed from that base is wrong by a plausible amount.
//!
//! # The fiction that makes it work
//!
//! `getReadIndexForReferenceCoordinate` advances a cigar treating **S as consuming reference**:
//!
//! ```java
//! lastRefPosOfElement += operator.consumesReferenceBases() || operator == CigarOperator.S
//!         ? element.getLength() : 0;
//! ```
//!
//! That is false of SAM's soft clip and true of this function, because it is called with the
//! *soft* start rather than the alignment start, so the clipped bases occupy the reference
//! positions they would have occupied had they aligned. Porting the honest rule, or calling it
//! with `getStart()`, shifts every answer on a soft-clipped read.

use htsjdk_bam::cigar::{Cigar, Op};
use htsjdk_bam::record::BamRecord;

use crate::read;

/// `ReadUtils.READ_INDEX_NOT_FOUND`.
pub const READ_INDEX_NOT_FOUND: i32 = -1;

/// `ReadConstants.UNSET_POSITION`.
pub const UNSET_POSITION: i32 = 0;

/// `GATKRead.getStart()`: the sentinel for an unmapped read, not its stored start.
pub fn start(record: &BamRecord) -> i32 {
    if read::is_unmapped(record) {
        UNSET_POSITION
    } else {
        record.alignment_start
    }
}

/// `GATKRead.getEnd()`.
pub fn end(record: &BamRecord) -> i32 {
    if read::is_unmapped(record) {
        UNSET_POSITION
    } else {
        record.alignment_end()
    }
}

/// `ReadUtils.getSoftStart`: the start the read would have if its soft clips had aligned.
///
/// The loop walks from the front subtracting soft clips and **skipping** hard clips, and stops at
/// the first element that is neither.
pub fn soft_start(record: &BamRecord) -> i32 {
    let mut soft_start = start(record);
    for element in &record.cigar.elements {
        match element.op {
            Op::S => soft_start -= element.length as i32,
            Op::H => {}
            _ => break,
        }
    }
    soft_start
}

/// `ReadUtils.getSoftEnd`.
///
/// The asymmetry with [`soft_start`] is the reference's: if the walk from the back finds no
/// aligned base at all, the soft end is reset to the alignment end rather than kept. `64H14S` is
/// the case the reference's own comment names, and it is why this cannot be written as the mirror
/// of the other function.
pub fn soft_end(record: &BamRecord) -> i32 {
    let mut found_aligned_base = false;
    let mut soft_end = end(record);
    for element in record.cigar.elements.iter().rev() {
        match element.op {
            Op::S => soft_end += element.length as i32,
            Op::H => {}
            _ => {
                found_aligned_base = true;
                break;
            }
        }
    }
    if !found_aligned_base {
        soft_end = end(record);
    }
    soft_end
}

/// `GATKRead.getUnclippedStart`, which is `SAMUtils.getUnclippedStart` under the adapter.
///
/// Unlike [`soft_start`], hard clips count: the walk subtracts `S` **and** `H` and stops at the
/// first element that is neither. An unmapped read has no start to unclip, so the sentinel wins.
pub fn unclipped_start(record: &BamRecord) -> i32 {
    if read::is_unmapped(record) {
        return UNSET_POSITION;
    }
    let mut unclipped = record.alignment_start;
    for element in &record.cigar.elements {
        match element.op {
            Op::S | Op::H => unclipped -= element.length as i32,
            _ => break,
        }
    }
    unclipped
}

/// `ReadUtils.hasWellDefinedFragmentSize`: can this read's adaptor be found from its mate?
///
/// Five refusals, and the order is the reference's. The commented-out `isProperlyPaired` check is
/// still commented out upstream, with the note that the flag is not always set correctly in BAMs;
/// reproducing the *live* code means not restoring it.
pub fn has_well_defined_fragment_size(record: &BamRecord) -> bool {
    if record.inferred_insert_size == 0 {
        return false; // mates on another contig, or an unmapped pair
    }
    if !read::is_paired(record) {
        return false;
    }
    if read::is_unmapped(record) || read::mate_is_unmapped(record) {
        return false;
    }
    if read::is_reverse_strand(record) == read::mate_is_reverse_strand(record) {
        return false;
    }
    if read::is_reverse_strand(record) {
        end(record) > record.mate_alignment_start
    } else {
        start(record) <= record.mate_alignment_start + record.inferred_insert_size
    }
}

/// `ReadUtils.getAdaptorBoundary`, or `None` for `CANNOT_COMPUTE_ADAPTOR_BOUNDARY`.
///
/// The two branches are not symmetric: on the reverse strand the boundary is the mate's start
/// minus one, a *measured* coordinate; on the forward strand it is this read's start plus the
/// absolute insert size, which is an inference and can land outside the read.
pub fn adaptor_boundary(record: &BamRecord) -> Option<i32> {
    if !has_well_defined_fragment_size(record) {
        return None;
    }
    if read::is_reverse_strand(record) {
        Some(record.mate_alignment_start - 1)
    } else {
        Some(start(record) + record.inferred_insert_size.abs())
    }
}

/// `ReadUtils.getLastInsertionOffset`, or `None` where the reference throws.
///
/// It indexes the last cigar element without checking that there is one, so a read with no cigar
/// raises an IndexOutOfBoundsException rather than answering zero.
pub fn last_insertion_offset(record: &BamRecord) -> Option<i32> {
    let last = record.cigar.elements.last()?;
    Some(if last.op == Op::I {
        last.length as i32
    } else {
        0
    })
}

/// `BaseUtils.getComplement`, which **throws** on anything that is not ACGTN.
///
/// Not a silent identity and not an N: a read of `*` bases, which `BaseUtils.isRegularBase`
/// happily calls regular, cannot be complemented at all.
pub fn complement(base: u8) -> Option<u8> {
    Some(match base {
        b'a' | b'A' => b'T',
        b'c' | b'C' => b'G',
        b'g' | b'G' => b'C',
        b't' | b'T' => b'A',
        b'n' | b'N' => b'N',
        _ => return None,
    })
}

/// `ReadUtils.getBasesReverseComplement`, or `None` where a base has no complement.
pub fn bases_reverse_complement(record: &BamRecord) -> Option<String> {
    let mut out = String::with_capacity(record.read_bases.len());
    for base in record.read_bases.iter().rev() {
        out.push(complement(*base)? as char);
    }
    Some(out)
}

/// `ReadUtils.getReadIndexForReferenceCoordinate(alignmentStart, cigar, refCoord)`.
///
/// Returns the read offset and the cigar operator it landed in, or `READ_INDEX_NOT_FOUND` with no
/// operator when the coordinate is before the start or past every element.
pub fn read_index_for_reference_coordinate(
    alignment_start: i32,
    cigar: &Cigar,
    ref_coord: i32,
) -> (i32, Option<Op>) {
    if ref_coord < alignment_start {
        return (READ_INDEX_NOT_FOUND, None);
    }
    let mut last_read_pos = 0;
    let mut last_ref_pos = alignment_start;
    for element in &cigar.elements {
        let first_read_pos = last_read_pos;
        let first_ref_pos = last_ref_pos;
        if element.op.consumes_read_bases() {
            last_read_pos += element.length as i32;
        }
        // S counts as consuming reference here; see the module comment.
        if element.op.consumes_reference_bases() || element.op == Op::S {
            last_ref_pos += element.length as i32;
        }
        if first_ref_pos <= ref_coord && ref_coord < last_ref_pos {
            let offset = if element.op.consumes_read_bases() {
                ref_coord - first_ref_pos
            } else {
                0
            };
            return (first_read_pos + offset, Some(element.op));
        }
    }
    (READ_INDEX_NOT_FOUND, None)
}

/// The same, for a read: the walk starts at the **soft** start.
pub fn read_index_for_read(record: &BamRecord, ref_coord: i32) -> (i32, Option<Op>) {
    read_index_for_reference_coordinate(soft_start(record), &record.cigar, ref_coord)
}

/// The outcome of asking a read for the base at a reference position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseAt {
    /// The reference answers `Optional.empty()`.
    Absent,
    Present(u8),
    /// The reference throws: the index it computed is outside the array it then indexes.
    Threw,
}

/// `ReadUtils.getReadBaseAtReferenceCoordinate`.
///
/// The bounds test uses `getStart()`/`getEnd()`, the alignment span, while the index comes from
/// the *soft* start. A soft-clipped base therefore cannot be reached through this function even
/// though the index walk knows where it is.
pub fn read_base_at_reference_coordinate(record: &BamRecord, ref_coord: i32) -> BaseAt {
    if ref_coord < start(record) || end(record) < ref_coord {
        return BaseAt::Absent;
    }
    let (index, op) = read_index_for_read(record, ref_coord);
    if index == READ_INDEX_NOT_FOUND || !op.is_some_and(|op| op.consumes_read_bases()) {
        return BaseAt::Absent;
    }
    match record.read_bases.get(index as usize) {
        Some(base) => BaseAt::Present(*base),
        None => BaseAt::Threw,
    }
}

/// `ReadUtils.getReadBaseQualityAtReferenceCoordinate`.
///
/// Not the mirror of the base accessor: this one tests only that the operator is non-null, where
/// the other also tests the index, and it indexes the quality array, which a read may not have.
/// A read with no qualities throws here and answers absent there.
pub fn read_base_quality_at_reference_coordinate(record: &BamRecord, ref_coord: i32) -> BaseAt {
    if ref_coord < start(record) || end(record) < ref_coord {
        return BaseAt::Absent;
    }
    let (index, op) = read_index_for_read(record, ref_coord);
    if !op.is_some_and(|op| op.consumes_read_bases()) {
        return BaseAt::Absent;
    }
    match record.base_qualities.get(index as usize) {
        Some(quality) => BaseAt::Present(*quality),
        None => BaseAt::Threw,
    }
}

/// `ReadUtils.isInsideRead`.
pub fn is_inside_read(record: &BamRecord, ref_coord: i32) -> bool {
    ref_coord >= start(record) && ref_coord <= end(record)
}

/// `ReadUtils.DEFAULT_INSERTION_DELETION_QUAL`, the flat Q45 the reference assumes when a read
/// carries no recalibrated indel qualities.
pub const DEFAULT_INSERTION_DELETION_QUAL: u8 = 45;

/// `ReadUtils.BQSR_BASE_INSERTION_QUALITIES`.
pub const BQSR_BASE_INSERTION_QUALITIES: [u8; 2] = *b"BI";
/// `ReadUtils.BQSR_BASE_DELETION_QUALITIES`.
pub const BQSR_BASE_DELETION_QUALITIES: [u8; 2] = *b"BD";

/// `SAMUtils.fastqToPhred` over a tag: the tag is FASTQ text, not raw phred bytes.
///
/// Returns `None` where the reference returns null, which is what makes the caller fall back to
/// the flat default rather than to an empty array.
fn indel_qualities(record: &BamRecord, tag: [u8; 2]) -> Option<Vec<u8>> {
    match record.tags.get(htsjdk_bam::tag::Tag::new(&tag)) {
        Some(htsjdk_bam::tag::TagValue::Str(text)) => {
            Some(text.bytes().map(|b| b.wrapping_sub(33)).collect())
        }
        _ => None,
    }
}

/// `ReadUtils.getBaseInsertionQualities`.
///
/// The fallback array is as long as the read's **quality** count, not its base count, which is a
/// difference for a read whose qualities are absent: the array is then empty and indexing it is
/// an error rather than a Q45 answer.
pub fn base_insertion_qualities(record: &BamRecord) -> Vec<u8> {
    indel_qualities(record, BQSR_BASE_INSERTION_QUALITIES)
        .unwrap_or_else(|| vec![DEFAULT_INSERTION_DELETION_QUAL; record.base_qualities.len()])
}

/// `ReadUtils.getBaseDeletionQualities`.
pub fn base_deletion_qualities(record: &BamRecord) -> Vec<u8> {
    indel_qualities(record, BQSR_BASE_DELETION_QUALITIES)
        .unwrap_or_else(|| vec![DEFAULT_INSERTION_DELETION_QUAL; record.base_qualities.len()])
}

/// `ReadUtils.getSAMFlagsForRead`: the flags **recomputed** from the accessors, not the stored
/// word.
///
/// This is not `record.flags`, and the difference is measurable. Every bit is rebuilt from a
/// `GATKRead` accessor, and three of those accessors are not the flag test they look like:
/// `isFirstOfPair` is `isPaired() && 0x40`, `isSecondOfPair` is `isPaired() && 0x80`, and
/// `isProperlyPaired` is `isPaired() && 0x2`. So a record carrying 0x40 without 0x1 keeps that bit
/// in its own header and loses it here. `SAM_READ_STRAND_FLAG` is likewise conditioned on the read
/// being mapped, so an unmapped reverse-strand read loses 0x10.
///
/// It matters because `ReadCoordinateComparator` breaks ties on this value, so two reads at the
/// same position can order differently under the stored flags and under these.
pub fn sam_flags_for_read(read: &BamRecord) -> i32 {
    let mut flags = 0i32;
    if read::is_paired(read) {
        flags |= read::flags::READ_PAIRED as i32;
    }
    if read::is_proper_pair(read) {
        flags |= read::flags::PROPER_PAIR as i32;
    }
    if read::is_unmapped(read) {
        flags |= read::flags::READ_UNMAPPED as i32;
    }
    if read::is_paired(read) && read::mate_is_unmapped(read) {
        flags |= read::flags::MATE_UNMAPPED as i32;
    }
    if !read::is_unmapped(read) && read::is_reverse_strand(read) {
        flags |= read::flags::READ_REVERSE_STRAND as i32;
    }
    if read::is_paired(read) && !read::mate_is_unmapped(read) && read::mate_is_reverse_strand(read)
    {
        flags |= read::flags::MATE_REVERSE_STRAND as i32;
    }
    if read::is_first_of_pair(read) {
        flags |= read::flags::FIRST_OF_PAIR as i32;
    }
    if read::is_second_of_pair(read) {
        flags |= read::flags::SECOND_OF_PAIR as i32;
    }
    if read::is_secondary_alignment(read) {
        flags |= read::flags::NOT_PRIMARY_ALIGNMENT as i32;
    }
    if read::fails_vendor_quality_check(read) {
        flags |= read::flags::READ_FAILS_VENDOR_QUALITY_CHECK as i32;
    }
    if read::is_duplicate(read) {
        flags |= read::flags::DUPLICATE_READ as i32;
    }
    if read::is_supplementary_alignment(read) {
        flags |= read::flags::SUPPLEMENTARY_ALIGNMENT as i32;
    }
    flags
}

/// `ReadUtils.getMateReferenceIndex`: `-1` when the mate is unmapped, whatever the record stores.
pub fn mate_reference_index(read: &BamRecord) -> i32 {
    if read::mate_is_unmapped(read) {
        read::NO_ALIGNMENT_REFERENCE_INDEX
    } else {
        read.mate_reference_index
    }
}

/// `ReadCoordinateComparator.compareCoordinates`.
///
/// The *assigned* position, not the one `getStart()` reports: an unmapped read placed at its
/// mate's coordinate sorts there rather than at the end, which is what interleaves unmapped reads
/// with mapped ones in a coordinate-sorted BAM. A read with no assigned reference sorts **after**
/// everything, and two of them compare equal.
pub fn compare_coordinates(first: &BamRecord, second: &BamRecord) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let first_ref = first.reference_index;
    let second_ref = second.reference_index;

    if first_ref == read::NO_ALIGNMENT_REFERENCE_INDEX {
        return if second_ref == read::NO_ALIGNMENT_REFERENCE_INDEX {
            Ordering::Equal
        } else {
            Ordering::Greater
        };
    }
    if second_ref == read::NO_ALIGNMENT_REFERENCE_INDEX {
        return Ordering::Less;
    }
    first_ref
        .cmp(&second_ref)
        .then(first.alignment_start.cmp(&second.alignment_start))
}

/// `ReadCoordinateComparator.compare`: coordinates first, then six tie-breakers in this order.
///
/// The strand tie-breaker is inverted on purpose and is the reference's: a read on the reverse
/// strand sorts **after** one on the forward strand at the same position. The comment upstream
/// says it mimics `SAMRecordCoordinateComparator`.
///
/// The name is compared with Java's `String.compareTo`, which is UTF-16 code-unit order, and the
/// flags are the recomputed ones from [`sam_flags_for_read`] rather than the stored word.
pub fn compare_read_coordinate(first: &BamRecord, second: &BamRecord) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let result = compare_coordinates(first, second);
    if result != Ordering::Equal {
        return result;
    }

    if read::is_reverse_strand(first) != read::is_reverse_strand(second) {
        return if read::is_reverse_strand(first) {
            Ordering::Greater
        } else {
            Ordering::Less
        };
    }

    let result = crate::java_hash::compare_strings(&first.read_name, &second.read_name);
    if result != Ordering::Equal {
        return result;
    }
    let result = sam_flags_for_read(first).cmp(&sam_flags_for_read(second));
    if result != Ordering::Equal {
        return result;
    }
    let result = (first.mapping_quality as i32).cmp(&(second.mapping_quality as i32));
    if result != Ordering::Equal {
        return result;
    }
    if read::is_paired(first) && read::is_paired(second) {
        let result = mate_reference_index(first).cmp(&mate_reference_index(second));
        if result != Ordering::Equal {
            return result;
        }
        let result = mate_start(first).cmp(&mate_start(second));
        if result != Ordering::Equal {
            return result;
        }
    }
    first.inferred_insert_size.cmp(&second.inferred_insert_size)
}

/// `GATKRead.getMateStart()`: the sentinel for an unmapped mate, not the stored value.
pub fn mate_start(read: &BamRecord) -> i32 {
    if read::mate_is_unmapped(read) {
        UNSET_POSITION
    } else {
        read.mate_alignment_start
    }
}
