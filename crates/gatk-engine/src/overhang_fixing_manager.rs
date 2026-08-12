//! `OverhangFixingManager`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.rnaseq.OverhangFixingManager` (GATK 4.6.2.0).
//!
//! The second piece of `SplitNCigarReads`, and the one that decides what the output reads look
//! like. The tool splits a read at every `N` and hands the resulting family here; this class holds
//! the families in a priority queue, remembers where the splices were, and soft-clips the overhang a
//! read leaves on the far side of a splice when the bases there disagree with the reference.
//!
//! # Three ways a span is refused before a single base is compared
//!
//! ```java
//! if ( spanToTest < 1 || spanToTest > maxBasesInOverhang || spanToTest > readLength / 2 ) {
//!     return false;
//! }
//! ```
//!
//! The last one is integer division against a **strict** comparison, so on a ten-base read a span of
//! five is still tested and six is refused, while on a nine-base read five is already refused. The
//! same span, two different answers, one base of read length apart.
//!
//! And there are two ways to say yes: more than `maxMismatchesInOverhang` mismatches returns early,
//! and a span where at least `(span+1)/2` bases mismatch returns true at the end even when the
//! tolerance was never exceeded. One mismatch out of two is a mismatch by the second rule while the
//! tolerance is one.
//!
//! # The queue is Java's, heap and all
//!
//! `waitingReadGroups` is a `PriorityQueue`, whose poll order for **equal** elements is decided by
//! the shape of its binary heap rather than by anything in the comparator. The class's own comment
//! says it does not guarantee coordinate order, so a port that used a stable sort would produce a
//! different, equally valid file, and no golden could hold. [`JavaPriorityQueue`] is therefore the
//! reference's `siftUp`/`siftDown` transcribed, so the two poll the same reads in the same order.
//!
//! # The mate key is deliberately asymmetric
//!
//! `makeKey(name, !isFirstOfPair, oldStart)` stores, `makeKey(name, isFirstOfPair, mateStart)` looks
//! up. A read stores the key its **mate** will search for, and the position in it is the read's
//! **old** start, so a mate already pointing at the clipped read's new position misses.

use std::collections::HashMap;

use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use crate::clipping::{self, ClipError};
use crate::interval::SimpleInterval;
use crate::read;
use crate::read_utils;
use crate::sa_tag;

/// The arguments `SplitNCigarReads` hands the manager.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OverhangArguments {
    /// `--max-reads-in-memory`, 150000 in the tool.
    pub max_records_in_memory: usize,
    /// `--max-mismatches-in-overhang`.
    pub max_mismatches_in_overhang: i32,
    /// `--max-bases-in-overhang`.
    pub max_bases_in_overhang: i32,
    /// `--do-not-fix-overhangs`.
    pub do_not_fix_overhangs: bool,
    /// `--process-secondary-alignments`.
    pub process_secondary_reads: bool,
}

impl Default for OverhangArguments {
    fn default() -> Self {
        OverhangArguments {
            max_records_in_memory: 150_000,
            max_mismatches_in_overhang: 1,
            max_bases_in_overhang: 40,
            do_not_fix_overhangs: false,
            process_secondary_reads: false,
        }
    }
}

/// The reference source a splice asks for its bases: contig, 1-based inclusive start and stop.
///
/// A callback rather than a reader, because the manager's only need of a reference is this one
/// query and the tool already holds an open one.
pub type ReferenceQuery<'a> = &'a mut dyn FnMut(&str, i32, i32) -> Result<Vec<u8>, String>;

/// `MAX_SPLICES_TO_KEEP`.
pub const MAX_SPLICES_TO_KEEP: usize = 1000;

/// What the manager refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OverhangError {
    /// `activateWriting` called twice.
    AlreadyWriting,
    /// An empty family, which `Utils.nonEmpty` refuses.
    EmptyReadGroup,
    /// The reference could not be read for a splice.
    Reference(String),
    /// A clip that the clipper itself refused.
    Clip(ClipError),
    /// A tag on a read of a family that could not be parsed.
    SaTag(sa_tag::SaTagError),
}

impl OverhangError {
    /// The message the reference carries, for the two it words itself.
    pub fn message(&self) -> String {
        match self {
            OverhangError::AlreadyWriting => {
                "Cannot activate writing for OverhangClippingManager multiple times".to_string()
            }
            OverhangError::EmptyReadGroup => {
                "readGroup added to manager is empty, which is not allowed".to_string()
            }
            OverhangError::Reference(text) => text.clone(),
            OverhangError::Clip(error) => format!("{error:?}"),
            OverhangError::SaTag(error) => error.message(),
        }
    }
}

/// `OverhangFixingManager.Splice`: a location and the reference bases under it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Splice {
    pub loc: SimpleInterval,
    pub reference: Vec<u8>,
}

/// `OverhangFixingManager.SplitRead`: a read with its soft-clip-included location, and what it
/// looked like when it arrived.
#[derive(Debug, Clone)]
pub struct SplitRead {
    pub read: BamRecord,
    /// None when the read is unmapped or consumes no reference bases.
    pub unclipped_loc: Option<SimpleInterval>,
    old_cigar: Cigar,
    old_start: i32,
}

impl SplitRead {
    /// `new SplitRead(read)`.
    pub fn new(record: &BamRecord, header: &SamHeader) -> SplitRead {
        let mut split = SplitRead {
            old_cigar: record.cigar.clone(),
            old_start: record.alignment_start,
            read: record.clone(),
            unclipped_loc: None,
        };
        split.set_read(record.clone(), header);
        split
    }

    /// `setRead`: the location is only assigned when the read is mapped **and** its soft start is
    /// strictly before its soft end, so a read consuming no reference is left without one.
    pub fn set_read(&mut self, record: BamRecord, header: &SamHeader) {
        let soft_start = read_utils::soft_start(&record);
        let soft_end = read_utils::soft_end(&record);
        self.unclipped_loc = if !read::is_unmapped(&record) && soft_start < soft_end {
            contig_of(&record, header)
                .and_then(|contig| SimpleInterval::new(contig, soft_start, soft_end))
        } else {
            None
        };
        self.read = record;
    }

    /// `hasBeenOverhangClipped`: either the cigar or the start changed.
    pub fn has_been_overhang_clipped(&self) -> bool {
        self.old_cigar != self.read.cigar || self.old_start != self.read.alignment_start
    }
}

/// `OverhangFixingManager`.
///
/// The writer is a `Vec` here rather than a trait object: the only caller collects the reads and
/// writes one BAM at the end, and a collected list is what a suite can compare.
pub struct OverhangFixingManager {
    header: SamHeader,
    arguments: OverhangArguments,
    splices: Vec<Splice>,
    waiting_read_groups: JavaPriorityQueue,
    waiting_reads: usize,
    output_to_file: bool,
    mate_changed_reads: HashMap<String, (i32, String)>,
    /// What the reference hands to `writer.addRead`, in the order it hands it over.
    pub written: Vec<BamRecord>,
}

impl OverhangFixingManager {
    pub fn new(header: &SamHeader, arguments: OverhangArguments) -> OverhangFixingManager {
        OverhangFixingManager {
            header: header.clone(),
            arguments,
            splices: Vec::new(),
            waiting_read_groups: JavaPriorityQueue::new(),
            waiting_reads: 0,
            output_to_file: false,
            mate_changed_reads: HashMap::new(),
            written: Vec::new(),
        }
    }

    /// `getNReadsInQueue`.
    pub fn reads_in_queue(&self) -> usize {
        self.waiting_reads
    }

    /// `getSplicesForTesting`, which the suite compares against.
    pub fn splices(&self) -> &[Splice] {
        &self.splices
    }

    /// `addSplicePosition`.
    ///
    /// Returns `None` for a splice already held **and** for a manager with overhang fixing turned
    /// off, which is the reference's own conflation: the two are told apart by whether anything is
    /// held afterwards.
    ///
    /// `reference` is asked for the bases only once the splice is known to be new, and the waiting
    /// reads are run against it there and then.
    pub fn add_splice_position(
        &mut self,
        contig: &str,
        start: i32,
        end: i32,
        reference: ReferenceQuery<'_>,
    ) -> Result<Option<Splice>, OverhangError> {
        if self.arguments.do_not_fix_overhangs {
            return Ok(None);
        }

        let Some(loc) = SimpleInterval::new(contig, start, end) else {
            return Ok(None);
        };
        if self.splices.iter().any(|splice| splice.loc == loc) {
            return Ok(None);
        }

        let bases = reference(contig, start, end).map_err(OverhangError::Reference)?;
        let splice = Splice {
            loc: loc.clone(),
            reference: bases,
        };

        // The contig is compared against the FIRST splice of the sorted set, not the last one
        // added, which is the same thing only because the set holds one contig at a time.
        let same_contig = self
            .splices
            .first()
            .map(|first| first.loc.contig == contig)
            .unwrap_or(true);
        if !same_contig {
            self.splices.clear();
        }

        let arguments = self.arguments.clone();
        let header = self.header.clone();
        for group in self.waiting_read_groups.iter_mut() {
            for split in group.iter_mut() {
                fix_split(split, &splice, &arguments, &header)?;
            }
        }

        self.insert_splice(splice.clone());

        if self.splices.len() > MAX_SPLICES_TO_KEEP {
            self.clean_splices();
        }
        Ok(Some(splice))
    }

    /// The `TreeSet` insertion: ordered by `GenomeLoc.compareTo`, which is contig index, then start,
    /// then stop.
    fn insert_splice(&mut self, splice: Splice) {
        let key = self.loc_key(&splice.loc);
        let position = self
            .splices
            .iter()
            .position(|held| self.loc_key(&held.loc) > key)
            .unwrap_or(self.splices.len());
        self.splices.insert(position, splice);
    }

    fn loc_key(&self, loc: &SimpleInterval) -> (i32, i32, i32) {
        let index = self
            .header
            .sequences
            .iter()
            .position(|sequence| sequence.name == loc.contig)
            .map(|index| index as i32)
            .unwrap_or(-1);
        (index, loc.start, loc.end)
    }

    /// `cleanSplices`: the lowest half of the set, dropped.
    fn clean_splices(&mut self) {
        let target = self.splices.len() / 2;
        self.splices.drain(0..target);
    }

    /// `addReadGroup`.
    pub fn add_read_group(&mut self, group: &[BamRecord]) -> Result<(), OverhangError> {
        if group.is_empty() {
            return Err(OverhangError::EmptyReadGroup);
        }

        // The flush decision is taken from the queue BEFORE the new family joins it.
        let too_many_reads = self.waiting_reads >= self.arguments.max_records_in_memory;
        let top_contig = self
            .waiting_read_groups
            .peek()
            .map(|group| group[0].read.clone());
        let first_new = &group[0];
        let encountered_new_contig = match &top_contig {
            Some(top) => {
                self.waiting_reads > 0
                    && !read::is_unmapped(top)
                    && !read::is_unmapped(first_new)
                    && contig_of(top, &self.header) != contig_of(first_new, &self.header)
            }
            None => false,
        };

        if too_many_reads || encountered_new_contig {
            // A new contig empties the queue; pressure alone empties half of it.
            let target = if encountered_new_contig {
                0
            } else {
                self.arguments.max_records_in_memory / 2
            };
            self.write_reads(target)?;
        }

        let mut new_group: Vec<SplitRead> = group
            .iter()
            .map(|record| SplitRead::new(record, &self.header))
            .collect();

        for splice in &self.splices {
            for split in new_group.iter_mut() {
                fix_split(split, splice, &self.arguments, &self.header)?;
            }
        }

        self.waiting_reads += new_group.len();
        self.waiting_read_groups.push(new_group, &self.header);
        Ok(())
    }

    /// `flush`.
    pub fn flush(&mut self) -> Result<(), OverhangError> {
        self.write_reads(0)
    }

    /// `writeReads(targetQueueSize)`.
    ///
    /// Before writing is active this records nothing but the mate repairs, and only for a family
    /// whose **first** read is not secondary and has actually been clipped.
    fn write_reads(&mut self, target: usize) -> Result<(), OverhangError> {
        while self.waiting_reads > target {
            let Some(group) = self.waiting_read_groups.poll(&self.header) else {
                break;
            };
            self.waiting_reads -= group.len();

            if self.output_to_file {
                let mut family: Vec<BamRecord> =
                    group.iter().map(|split| split.read.clone()).collect();
                repair_supplementary_tags(&mut family, &self.header)?;
                self.written.extend(family);
            } else if !read::is_secondary_alignment(&group[0].read)
                && group[0].has_been_overhang_clipped()
            {
                self.set_mate_changed(&group[0]);
            }
        }
        Ok(())
    }

    /// `SplitRead.setMateChanged`: the key a read stores is the one its **mate** will look up.
    fn set_mate_changed(&mut self, split: &SplitRead) {
        if read::is_unmapped(&split.read) {
            return;
        }
        let key = make_key(
            &split.read.read_name,
            !read::is_first_of_pair(&split.read),
            split.old_start,
        );
        self.mate_changed_reads.insert(
            key,
            (split.read.alignment_start, split.read.cigar.to_text()),
        );
    }

    /// `activateWriting`: flushes what is waiting, forgets the splices, and refuses a second call.
    pub fn activate_writing(&mut self) -> Result<(), OverhangError> {
        if self.output_to_file {
            return Err(OverhangError::AlreadyWriting);
        }
        self.flush()?;
        self.splices.clear();
        self.output_to_file = true;
        Ok(())
    }

    /// `setPredictedMateInformation`: move a mate to where its partner was clipped to.
    ///
    /// Returns whether the read was edited. Before writing is active this is always false, which is
    /// what makes the first pass a recording pass.
    pub fn set_predicted_mate_information(&self, record: &mut BamRecord) -> bool {
        if !self.output_to_file {
            return false;
        }
        if record.read_bases.is_empty() || !read::is_paired(record) {
            return false;
        }
        let key = make_key(
            &record.read_name,
            read::is_first_of_pair(record),
            record.mate_alignment_start,
        );
        let Some((start, cigar)) = self.mate_changed_reads.get(&key) else {
            return false;
        };
        record.mate_alignment_start = *start;
        if record.tags.get(Tag::new(b"MC")).is_some() {
            record
                .tags
                .insert(Tag::new(b"MC"), TagValue::Str(cigar.clone()));
        }
        true
    }
}

/// `SplitNCigarReads.repairSupplementaryTags`.
///
/// It lives here because it is what the manager calls on the way out, and it is the only place the
/// three tags are cleared: NM, MD and NH break when a read is split, and nothing recomputes them.
pub fn repair_supplementary_tags(
    family: &mut Vec<BamRecord>,
    header: &SamHeader,
) -> Result<(), OverhangError> {
    for record in family.iter_mut() {
        for name in [b"NM", b"MD", b"NH"] {
            record.tags.remove(Tag::new(name));
        }
    }
    if family.len() > 1 {
        let mut primary = family.remove(0);
        sa_tag::set_reads_as_supplemental(&mut primary, family, header)
            .map_err(OverhangError::SaTag)?;
        family.insert(0, primary);
    }
    Ok(())
}

/// `fixSplit`: clip the overhang this read leaves across this splice, if its bases disagree.
pub fn fix_split(
    split: &mut SplitRead,
    splice: &Splice,
    arguments: &OverhangArguments,
    header: &SamHeader,
) -> Result<(), OverhangError> {
    let Some(read_loc) = split.unclipped_loc.clone() else {
        return Ok(());
    };
    if !overlaps(&splice.loc, &read_loc) {
        return Ok(());
    }
    if !arguments.process_secondary_reads && read::is_secondary_alignment(&split.read) {
        return Ok(());
    }

    let record = split.read.clone();
    // The bases the cigar consumes that are not clipping, which is what stops a soft-clipped read
    // from being clipped again on the strength of bases it no longer aligns.
    let read_bases_length: i32 = record
        .cigar
        .elements
        .iter()
        .filter(|element| element.op.consumes_read_bases() && !is_clipping(element.op))
        .map(|element| element.length as i32)
        .sum();

    if is_left_overhang(&read_loc, &splice.loc) {
        let overhang = splice.loc.end - record.alignment_start + 1;
        if overhanging_bases_mismatch(
            &record.read_bases,
            record.alignment_start - read_loc.start,
            read_bases_length,
            &splice.reference,
            splice.reference.len() as i32 - overhang,
            overhang,
            arguments,
        ) {
            let clipped = clipping::soft_clip_by_read_coordinates(
                &record,
                Some(header),
                0,
                splice.loc.end - read_loc.start,
            )
            .map_err(OverhangError::Clip)?;
            split.set_read(clipped, header);
        }
    } else if is_right_overhang(&read_loc, &splice.loc) {
        let overhang = read_loc.end - splice.loc.start + 1;
        let length = record.read_bases.len() as i32;
        if overhanging_bases_mismatch(
            &record.read_bases,
            length - overhang,
            read_bases_length,
            &splice.reference,
            0,
            read_utils::end(&record) - splice.loc.start + 1,
            arguments,
        ) {
            let clipped = clipping::soft_clip_by_read_coordinates(
                &record,
                Some(header),
                length - overhang,
                length - 1,
            )
            .map_err(OverhangError::Clip)?;
            split.set_read(clipped, header);
        }
    }
    Ok(())
}

/// `isLeftOverhang`: starts inside the splice, past its start, and runs beyond its end.
pub fn is_left_overhang(read_loc: &SimpleInterval, splice_loc: &SimpleInterval) -> bool {
    read_loc.start <= splice_loc.end
        && read_loc.start > splice_loc.start
        && read_loc.end > splice_loc.end
}

/// `isRightOverhang`: ends inside the splice, before its end, and starts before its start.
pub fn is_right_overhang(read_loc: &SimpleInterval, splice_loc: &SimpleInterval) -> bool {
    read_loc.end >= splice_loc.start
        && read_loc.end < splice_loc.end
        && read_loc.start < splice_loc.start
}

/// `overhangingBasesMismatch`.
///
/// Three refusals before anything is compared, then two ways to say yes: the tolerance, and the
/// half rule at the end. The `spanToTest > readLength / 2` gate is integer division against a
/// strict comparison, so a span of five is tested on a ten-base read and refused on a nine-base one.
pub fn overhanging_bases_mismatch(
    read_bases: &[u8],
    read_start_index: i32,
    read_length: i32,
    reference: &[u8],
    reference_start_index: i32,
    span_to_test: i32,
    arguments: &OverhangArguments,
) -> bool {
    if span_to_test < 1
        || span_to_test > arguments.max_bases_in_overhang
        || span_to_test > read_length / 2
    {
        return false;
    }

    let mut mismatches = 0;
    for i in 0..span_to_test {
        let read_index = read_start_index + i;
        let reference_index = reference_start_index + i;
        // Java would throw ArrayIndexOutOfBounds here; a port that indexed out of range would
        // panic, so an out-of-range span is treated as no mismatch at that base. Nothing in the
        // reference reaches it: the caller's arithmetic keeps both indices inside their arrays.
        let (Some(&base), Some(&expected)) = (
            usize::try_from(read_index)
                .ok()
                .and_then(|i| read_bases.get(i)),
            usize::try_from(reference_index)
                .ok()
                .and_then(|i| reference.get(i)),
        ) else {
            continue;
        };
        if base != expected {
            mismatches += 1;
            if mismatches > arguments.max_mismatches_in_overhang {
                return true;
            }
        }
    }

    mismatches >= (span_to_test + 1) / 2
}

/// `makeKey(name, firstOfPair, mateStart)`, whose `@` is forbidden in read names by the SAM spec.
pub fn make_key(name: &str, first_of_pair: bool, mate_start: i32) -> String {
    format!("{name}@{}@{mate_start}", if first_of_pair { 1 } else { 0 })
}

/// `GenomeLoc.overlapsP`, for two locations already known to be on some contig.
fn overlaps(first: &SimpleInterval, second: &SimpleInterval) -> bool {
    first.contig == second.contig && first.start <= second.end && second.start <= first.end
}

/// `CigarOperator.isClipping`.
fn is_clipping(op: htsjdk_bam::cigar::Op) -> bool {
    use htsjdk_bam::cigar::Op;
    matches!(op, Op::S | Op::H)
}

fn contig_of<'a>(record: &BamRecord, header: &'a SamHeader) -> Option<&'a str> {
    usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.as_str())
}

/// `java.util.PriorityQueue`, transcribed for the one element type this holds.
///
/// The order two equal elements come out in is a property of the heap's shape and not of the
/// comparator, and the manager's own documentation says its output is not guaranteed to be sorted.
/// Reproducing `siftUp` and `siftDown` is therefore the only way the two implementations write the
/// same file.
pub struct JavaPriorityQueue {
    queue: Vec<Vec<SplitRead>>,
}

impl JavaPriorityQueue {
    pub fn new() -> JavaPriorityQueue {
        JavaPriorityQueue { queue: Vec::new() }
    }

    pub fn len(&self) -> usize {
        self.queue.len()
    }

    pub fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }

    pub fn peek(&self) -> Option<&Vec<SplitRead>> {
        self.queue.first()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Vec<SplitRead>> {
        self.queue.iter_mut()
    }

    /// `offer(e)`: append, then sift it up.
    pub fn push(&mut self, element: Vec<SplitRead>, header: &SamHeader) {
        let mut k = self.queue.len();
        self.queue.push(element);
        while k > 0 {
            let parent = (k - 1) >> 1;
            if compare(&self.queue[k], &self.queue[parent], header) >= std::cmp::Ordering::Equal {
                break;
            }
            self.queue.swap(k, parent);
            k = parent;
        }
    }

    /// `poll()`: take the head, move the last element to the root and sift it down.
    pub fn poll(&mut self, header: &SamHeader) -> Option<Vec<SplitRead>> {
        if self.queue.is_empty() {
            return None;
        }
        let result = self.queue.swap_remove(0);
        if self.queue.len() > 1 {
            self.sift_down(0, header);
        }
        Some(result)
    }

    fn sift_down(&mut self, mut k: usize, header: &SamHeader) {
        let size = self.queue.len();
        let half = size >> 1;
        while k < half {
            let mut child = 2 * k + 1;
            let right = child + 1;
            if right < size
                && compare(&self.queue[child], &self.queue[right], header)
                    == std::cmp::Ordering::Greater
            {
                child = right;
            }
            if compare(&self.queue[k], &self.queue[child], header) <= std::cmp::Ordering::Equal {
                break;
            }
            self.queue.swap(k, child);
            k = child;
        }
    }
}

impl Default for JavaPriorityQueue {
    fn default() -> Self {
        JavaPriorityQueue::new()
    }
}

/// `SplitReadComparator`: the first read of each family, through `ReadCoordinateComparator`.
fn compare(first: &[SplitRead], second: &[SplitRead], _header: &SamHeader) -> std::cmp::Ordering {
    read_utils::compare_read_coordinate(&first[0].read, &second[0].read)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn arguments() -> OverhangArguments {
        OverhangArguments::default()
    }

    #[test]
    fn the_span_gate_is_strict_and_halves_downwards() {
        let read = b"ACGTACGTAC";
        let two_off = b"TTGTACGTAC";
        // Six of ten is refused, five of ten is tested, five of nine is refused.
        assert!(!overhanging_bases_mismatch(
            read,
            0,
            10,
            two_off,
            0,
            6,
            &arguments()
        ));
        assert!(overhanging_bases_mismatch(
            read,
            0,
            10,
            two_off,
            0,
            5,
            &arguments()
        ));
        assert!(!overhanging_bases_mismatch(
            read,
            0,
            9,
            two_off,
            0,
            5,
            &arguments()
        ));
    }

    #[test]
    fn one_mismatch_of_two_is_a_mismatch_although_the_tolerance_is_one() {
        let read = b"ACGTACGTAC";
        let one_off = b"TCGTACGTAC";
        // The tolerance is never exceeded; the half rule is what returns true.
        assert!(overhanging_bases_mismatch(
            read,
            0,
            10,
            one_off,
            0,
            2,
            &arguments()
        ));
        assert!(!overhanging_bases_mismatch(
            read,
            0,
            10,
            one_off,
            0,
            4,
            &arguments()
        ));
    }

    #[test]
    fn the_overhang_predicates_are_mixed_strict_and_non_strict() {
        let splice = SimpleInterval::new("chr1", 50, 60).expect("a valid splice");
        let starts_at_splice_start = SimpleInterval::new("chr1", 50, 70).expect("valid");
        let starts_one_later = SimpleInterval::new("chr1", 51, 70).expect("valid");
        let spans_the_splice = SimpleInterval::new("chr1", 40, 70).expect("valid");
        let ends_inside = SimpleInterval::new("chr1", 40, 55).expect("valid");

        assert!(!is_left_overhang(&starts_at_splice_start, &splice));
        assert!(is_left_overhang(&starts_one_later, &splice));
        assert!(!is_left_overhang(&spans_the_splice, &splice));
        assert!(!is_right_overhang(&spans_the_splice, &splice));
        assert!(is_right_overhang(&ends_inside, &splice));
    }

    #[test]
    fn the_mate_key_carries_the_other_end_of_the_pair() {
        assert_eq!(make_key("read1", true, 100), "read1@1@100");
        assert_eq!(make_key("read1", false, 100), "read1@0@100");
    }
}
