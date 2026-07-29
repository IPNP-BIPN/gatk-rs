//! Ported from `org.broadinstitute.hellbender.utils.clipping.ClippingOp`, `ClippingRepresentation`
//! and `ReadClipper` (GATK 4.6.2.0).
//!
//! Clipping rewrites the read itself: its bases, its qualities, its cigar, its alignment start and,
//! for a flow read group, some of its tags. Everything downstream reads the clipped read, so a
//! wrong clip is not a wrong number, it is a different read.
//!
//! Three things here are decisions rather than mechanics:
//!
//!  * a hard clip that would remove every base returns an **empty read**, which is unmapped with a
//!    mapping quality of zero and no attributes except its read group. It is not a read of length
//!    zero at the old position;
//!  * a soft clip cannot take the whole read. `applySoftClipBases` caps the clip two bases short
//!    (GATK issue #2022), and the cap uses the *capped* stop for the cigar and the **uncapped**
//!    stop for the alignment-start shift;
//!  * reverting soft clips can move a read's start before the contig, which SAM cannot represent.
//!    The reference sets the position to 1, hard-clips the overhang away, and sets it to 1 again,
//!    because the first assignment is what makes the read's cached length right for the clip.

use htsjdk_bam::cigar::{Cigar, Op};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use crate::cigar_builder::CigarError;
use crate::cigar_utils::{alignment_start_shift, clip_cigar, revert_soft_clips};
use crate::read;
use crate::read_utils;

/// `ClippingRepresentation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClippingRepresentation {
    /// Clipped bases become N.
    WriteNs,
    /// Clipped bases get a quality of zero.
    WriteQ0s,
    /// Both.
    WriteNsQ0s,
    SoftclipBases,
    HardclipBases,
    RevertSoftclippedBases,
}

/// What the reference throws rather than returning a read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ClipError {
    /// `applySoftClipBases`: soft-clipping an unmapped read, or the middle of a read.
    CannotSoftClip,
    /// A cigar the builder refuses; see [`CigarError`].
    Cigar(CigarError),
    /// `clipByReferenceCoordinates`: the arguments do not describe one tail.
    BadCoordinates,
    /// An `arraycopy` past the end of the array the reference copies from.
    IndexOutOfBounds,
}

impl From<CigarError> for ClipError {
    fn from(error: CigarError) -> Self {
        ClipError::Cigar(error)
    }
}

/// `FlowBasedRead.HARD_CLIPPED_TAGS`: the per-base tags a hard clip trims with the read.
const HARD_CLIPPED_TAGS: [&[u8; 2]; 7] = [b"tp", b"kr", b"ti", b"kh", b"kf", b"kd", b"t0"];

/// `ReadUtils.emptyRead`.
///
/// Unmapped, mapping quality zero, no cigar, no bases, no qualities, and **no attributes except
/// the read group**, which is copied back after the rest are cleared.
pub fn empty_read(record: &BamRecord) -> BamRecord {
    let read_group = record
        .tags
        .iter()
        .find(|(tag, _)| tag.name() == *b"RG")
        .map(|(tag, value)| (*tag, value.clone()));
    let mut empty = record.clone();
    empty.flags |= read::flags::READ_UNMAPPED;
    empty.mapping_quality = 0;
    empty.cigar = Cigar::default();
    empty.read_bases.clear();
    empty.base_qualities.clear();
    empty.tags = htsjdk_bam::tag::Tags::new();
    if let Some((tag, value)) = read_group {
        empty.tags.insert(tag, value);
    }
    empty
}

/// `SAMRecordToGATKReadAdapter.hardClipAttributes`.
///
/// Applies **only** to reads in a read group that declares a flow order, and only to a fixed
/// whitelist of tags whose length matches the original read's. Everything else is left alone, so a
/// hard clip on an ordinary read touches no tag at all.
fn hard_clip_attributes(
    record: &mut BamRecord,
    header: Option<&SamHeader>,
    new_start: usize,
    new_length: usize,
    original_length: usize,
) {
    let Some(header) = header else { return };
    let has_flow_order = crate::read_group::flow_order(record, header).is_some();
    if !has_flow_order {
        return;
    }
    for name in HARD_CLIPPED_TAGS {
        let tag = htsjdk_bam::Tag::new(name);
        let Some(value) = record.tags.get(tag).cloned() else {
            continue;
        };
        let clipped = match value {
            htsjdk_bam::tag::TagValue::ByteArray { values, unsigned } => {
                if values.len() != original_length {
                    continue;
                }
                htsjdk_bam::tag::TagValue::ByteArray {
                    values: values[new_start..new_start + new_length].to_vec(),
                    unsigned,
                }
            }
            htsjdk_bam::tag::TagValue::Str(text) => {
                if text.len() != original_length {
                    continue;
                }
                htsjdk_bam::tag::TagValue::Str(text[new_start..new_start + new_length].to_string())
            }
            // Any other type is left alone: the reference skips a value that is neither an array
            // nor a string.
            _ => continue,
        };
        record.tags.insert(tag, clipped);
    }
}

/// `ClippingOp`: a half-open-looking but **inclusive** span of read offsets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClippingOp {
    pub start: i32,
    pub stop: i32,
}

impl ClippingOp {
    /// `ClippingOp.apply`.
    pub fn apply(
        &self,
        algorithm: ClippingRepresentation,
        record: &BamRecord,
        header: Option<&SamHeader>,
    ) -> Result<BamRecord, ClipError> {
        match algorithm {
            ClippingRepresentation::WriteNs => {
                let mut copy = record.clone();
                self.write_ns(&mut copy);
                Ok(copy)
            }
            ClippingRepresentation::WriteQ0s => {
                let mut copy = record.clone();
                self.write_q0s(&mut copy);
                Ok(copy)
            }
            ClippingRepresentation::WriteNsQ0s => {
                let mut copy = record.clone();
                self.write_ns(&mut copy);
                self.write_q0s(&mut copy);
                Ok(copy)
            }
            ClippingRepresentation::HardclipBases => {
                hard_clip_bases(record, self.start, self.stop, header)
            }
            ClippingRepresentation::SoftclipBases => self.soft_clip_bases(record),
            ClippingRepresentation::RevertSoftclippedBases => {
                revert_softclipped_bases(record, header)
            }
        }
    }

    /// `overwriteFromStartToStop`, which stops at the array's end rather than failing.
    fn overwrite(&self, values: &mut [u8], new_value: u8) {
        let start = self.start as usize;
        let stop = (self.stop as usize + 1).min(values.len());
        if start < stop {
            values[start..stop].fill(new_value);
        }
    }

    fn write_ns(&self, record: &mut BamRecord) {
        // Nucleotide.N.encodeAsByte() is upper-case N.
        let mut bases = record.read_bases.clone();
        self.overwrite(&mut bases, b'N');
        record.read_bases = bases;
    }

    fn write_q0s(&self, record: &mut BamRecord) {
        let mut quals = record.base_qualities.clone();
        self.overwrite(&mut quals, 0);
        record.base_qualities = quals;
    }

    /// `applySoftClipBases`.
    fn soft_clip_bases(&self, record: &BamRecord) -> Result<BamRecord, ClipError> {
        if read::is_unmapped(record) {
            return Err(ClipError::CannotSoftClip);
        }
        let length = record.read_bases.len() as i32;
        // GATK issue #2022: a read cannot be entirely soft-clipped.
        if length <= 2 {
            return Ok(record.clone());
        }
        let my_stop = self.stop.min(self.start + length - 2);
        if !(self.start <= 0 || my_stop == length - 1) {
            return Err(ClipError::CannotSoftClip);
        }
        let mut copy = record.clone();
        let old_cigar = record.cigar.clone();
        copy.cigar = clip_cigar(&old_cigar, self.start, my_stop + 1, Op::S)?;
        // The shift uses the *uncapped* stop while the cigar used the capped one. That asymmetry
        // is the reference's, and it only shows on a read the cap actually changed.
        let shift = if self.start == 0 {
            alignment_start_shift(&old_cigar, self.stop + 1)
        } else {
            0
        };
        copy.alignment_start = read_utils::start(record) + shift;
        Ok(copy)
    }
}

/// `ClippingOp.applyHardClipBases`.
pub fn hard_clip_bases(
    record: &BamRecord,
    start: i32,
    stop: i32,
    header: Option<&SamHeader>,
) -> Result<BamRecord, ClipError> {
    let original_length = record.read_bases.len();
    let new_length = original_length as i32 - (stop - start + 1);
    // A hard clip that removes everything is an empty read, not a zero-length one: the difference
    // is the mapping quality, the flags and the attributes.
    if new_length == 0 {
        return Ok(empty_read(record));
    }
    let new_cigar = if read::is_unmapped(record) {
        Cigar::default()
    } else {
        clip_cigar(&record.cigar, start, stop + 1, Op::H)?
    };
    let copy_start = if start == 0 { (stop + 1) as usize } else { 0 };
    let new_length = new_length as usize;

    let mut clipped = record.clone();
    hard_clip_attributes(
        &mut clipped,
        header,
        copy_start,
        new_length,
        original_length,
    );
    // `System.arraycopy` from the read's bases and qualities. A read whose qualities are absent
    // has an empty array, and the copy then runs past its end: the reference raises
    // ArrayIndexOutOfBoundsException rather than producing a read with no qualities. Measured on
    // the corpus record that carries no qualities, which this port first got wrong by treating
    // absent qualities as a special case worth surviving.
    let end = copy_start + new_length;
    if record.read_bases.len() < end || record.base_qualities.len() < end {
        return Err(ClipError::IndexOutOfBounds);
    }
    clipped.read_bases = record.read_bases[copy_start..end].to_vec();
    clipped.base_qualities = record.base_qualities[copy_start..end].to_vec();
    clipped.cigar = new_cigar;
    if start == 0 && !read::is_unmapped(record) {
        clipped.alignment_start =
            read_utils::start(record) + alignment_start_shift(&record.cigar, stop + 1);
    }
    Ok(clipped)
}

/// `ClippingOp.applyRevertSoftClippedBases`.
pub fn revert_softclipped_bases(
    record: &BamRecord,
    header: Option<&SamHeader>,
) -> Result<BamRecord, ClipError> {
    let elements = &record.cigar.elements;
    let clipped_at_an_end = elements
        .first()
        .is_some_and(|e| matches!(e.op, Op::S | Op::H))
        || elements
            .last()
            .is_some_and(|e| matches!(e.op, Op::S | Op::H));
    if elements.is_empty() || !clipped_at_an_end {
        return Ok(record.clone());
    }
    let mut unclipped = record.clone();
    unclipped.cigar = revert_soft_clips(&record.cigar)?;
    let new_start = read_utils::soft_start(record);
    if new_start <= 0 {
        // The unclipped read would start before the contig, which SAM cannot hold: the reference
        // sets the position to 1, hard-clips the overhang, and sets it to 1 again, because the
        // first assignment is what makes the read's length right for the clip.
        unclipped.alignment_start = 1;
        let mut unclipped = hard_clip_bases(&unclipped, 0, -new_start, header)?;
        if !read::is_unmapped(&unclipped) {
            unclipped.alignment_start = 1;
        }
        Ok(unclipped)
    } else {
        unclipped.alignment_start = new_start;
        Ok(unclipped)
    }
}

/// `ReadClipper`: a read plus the clips to apply to it.
pub struct ReadClipper<'a> {
    read: BamRecord,
    header: Option<&'a SamHeader>,
    ops: Vec<ClippingOp>,
}

impl<'a> ReadClipper<'a> {
    pub fn new(read: &BamRecord, header: Option<&'a SamHeader>) -> ReadClipper<'a> {
        ReadClipper {
            read: read.clone(),
            header,
            ops: Vec::new(),
        }
    }

    pub fn add_op(&mut self, op: ClippingOp) {
        self.ops.push(op);
    }

    /// `ReadClipper.clipRead`.
    ///
    /// An op whose start is past the current read is skipped entirely, and one whose stop is past
    /// it is shortened: the read shrinks as ops are applied, so the second op of a pair is
    /// interpreted against a read the first one already cut.
    pub fn clip_read(&mut self, algorithm: ClippingRepresentation) -> Result<BamRecord, ClipError> {
        if self.ops.is_empty() {
            return Ok(self.read.clone());
        }
        let mut clipped = self.read.clone();
        for op in &self.ops {
            let read_length = clipped.read_bases.len() as i32;
            if op.start < read_length {
                let fixed = if op.stop >= read_length {
                    ClippingOp {
                        start: op.start,
                        stop: read_length - 1,
                    }
                } else {
                    *op
                };
                clipped = fixed.apply(algorithm, &clipped, self.header)?;
            }
        }
        self.ops.clear();
        if clipped.read_bases.is_empty() {
            return Ok(empty_read(&clipped));
        }
        Ok(clipped)
    }

    /// `ReadClipper.clipByReferenceCoordinates`: exactly one of the two bounds is `< 0`.
    pub fn clip_by_reference_coordinates(
        &mut self,
        ref_start: i32,
        ref_stop: i32,
        algorithm: ClippingRepresentation,
    ) -> Result<BamRecord, ClipError> {
        if self.read.read_bases.is_empty() {
            return Ok(self.read.clone());
        }
        if algorithm == ClippingRepresentation::SoftclipBases && read::is_unmapped(&self.read) {
            return Err(ClipError::CannotSoftClip);
        }
        let (start, stop) = if ref_start < 0 {
            if ref_stop < 0 {
                return Err(ClipError::BadCoordinates);
            }
            let (index, operator) = read_utils::read_index_for_read(&self.read, ref_stop);
            // A reference coordinate inside a deletion returns the position *after* it, and the
            // stop here is inclusive, so it is decremented: the deletion is left unclipped rather
            // than one base too much being taken.
            let stop = index
                - if operator.is_some_and(|op| op.consumes_read_bases()) {
                    0
                } else {
                    1
                };
            (0, stop)
        } else {
            if ref_stop >= 0 {
                return Err(ClipError::BadCoordinates);
            }
            let (index, _) = read_utils::read_index_for_read(&self.read, ref_start);
            (index, self.read.read_bases.len() as i32 - 1)
        };

        if start == read_utils::READ_INDEX_NOT_FOUND || stop == read_utils::READ_INDEX_NOT_FOUND {
            return Ok(self.read.clone());
        }
        if start < 0 || stop > self.read.read_bases.len() as i32 - 1 {
            return Err(ClipError::BadCoordinates);
        }
        if start > stop {
            return Err(ClipError::BadCoordinates);
        }
        if start > 0 && stop < self.read.read_bases.len() as i32 - 1 {
            // Clipping the middle of a read is refused, whatever the representation.
            return Err(ClipError::BadCoordinates);
        }
        self.add_op(ClippingOp { start, stop });
        self.clip_read(algorithm)
    }
}

/// `ReadClipper.hardClipByReferenceCoordinatesLeftTail`.
pub fn hard_clip_by_reference_coordinates_left_tail(
    read: &BamRecord,
    header: Option<&SamHeader>,
    ref_stop: i32,
) -> Result<BamRecord, ClipError> {
    ReadClipper::new(read, header).clip_by_reference_coordinates(
        -1,
        ref_stop,
        ClippingRepresentation::HardclipBases,
    )
}

/// `ReadClipper.hardClipByReferenceCoordinatesRightTail`.
pub fn hard_clip_by_reference_coordinates_right_tail(
    read: &BamRecord,
    header: Option<&SamHeader>,
    ref_start: i32,
) -> Result<BamRecord, ClipError> {
    ReadClipper::new(read, header).clip_by_reference_coordinates(
        ref_start,
        -1,
        ClippingRepresentation::HardclipBases,
    )
}

/// `ReadClipper.hardClipBothEndsByReferenceCoordinates`.
///
/// The right tail goes first, and the left bound is then checked against the **clipped** read:
/// hard-clipping one end can remove adjacent deletions, which moves the other end's coordinate.
pub fn hard_clip_both_ends_by_reference_coordinates(
    read: &BamRecord,
    header: Option<&SamHeader>,
    left: i32,
    right: i32,
) -> Result<BamRecord, ClipError> {
    if read.read_bases.is_empty() || left == right {
        return Ok(empty_read(read));
    }
    let left_tail = hard_clip_by_reference_coordinates_right_tail(read, header, right)?;
    if left > read_utils::end(&left_tail) {
        return Ok(empty_read(read));
    }
    hard_clip_by_reference_coordinates_left_tail(&left_tail, header, left)
}

/// `ReadClipper.hardClipToRegion`.
pub fn hard_clip_to_region(
    read: &BamRecord,
    header: Option<&SamHeader>,
    ref_start: i32,
    ref_stop: i32,
) -> Result<BamRecord, ClipError> {
    let alignment_start = read_utils::start(read);
    let alignment_stop = read_utils::end(read);
    if alignment_start <= ref_stop && alignment_stop >= ref_start {
        if alignment_start < ref_start && alignment_stop > ref_stop {
            hard_clip_both_ends_by_reference_coordinates(read, header, ref_start - 1, ref_stop + 1)
        } else if alignment_start < ref_start {
            hard_clip_by_reference_coordinates_left_tail(read, header, ref_start - 1)
        } else if alignment_stop > ref_stop {
            hard_clip_by_reference_coordinates_right_tail(read, header, ref_stop + 1)
        } else {
            Ok(read.clone())
        }
    } else {
        Ok(empty_read(read))
    }
}

/// `ReadClipper.revertSoftClippedBases`.
pub fn revert_soft_clipped_bases(
    read: &BamRecord,
    header: Option<&SamHeader>,
) -> Result<BamRecord, ClipError> {
    if read.read_bases.is_empty() {
        return Ok(read.clone());
    }
    let mut clipper = ReadClipper::new(read, header);
    // The op is (0, 0) and is ignored by the representation, which rewrites the whole cigar.
    clipper.add_op(ClippingOp { start: 0, stop: 0 });
    clipper.clip_read(ClippingRepresentation::RevertSoftclippedBases)
}

/// `ReadClipper.clipLowQualEnds`.
///
/// The right end is found first and the ops are added right-then-left, which matters: `clipRead`
/// applies them in order against a read that the previous op has already shortened.
pub fn clip_low_qual_ends(
    read: &BamRecord,
    header: Option<&SamHeader>,
    low_qual: u8,
    algorithm: ClippingRepresentation,
) -> Result<BamRecord, ClipError> {
    if read.read_bases.is_empty() {
        return Ok(read.clone());
    }
    // The length is the read's, and the qualities are indexed with it: a read whose qualities are
    // absent has an empty array, and `getBaseQuality` runs past its end rather than answering.
    let read_length = read.read_bases.len() as i32;
    let quality = |index: i32| -> Option<u8> { read.base_qualities.get(index as usize).copied() };
    let mut right = read_length - 1;
    while right >= 0 {
        match quality(right) {
            None => return Err(ClipError::IndexOutOfBounds),
            Some(q) if q <= low_qual => right -= 1,
            Some(_) => break,
        }
    }
    let mut left = 0;
    while left < read_length {
        match quality(left) {
            None => return Err(ClipError::IndexOutOfBounds),
            Some(q) if q <= low_qual => left += 1,
            Some(_) => break,
        }
    }
    if left > right {
        return Ok(empty_read(read));
    }
    let mut clipper = ReadClipper::new(read, header);
    if right < read_length - 1 {
        clipper.add_op(ClippingOp {
            start: right + 1,
            stop: read_length - 1,
        });
    }
    if left > 0 {
        clipper.add_op(ClippingOp {
            start: 0,
            stop: left - 1,
        });
    }
    clipper.clip_read(algorithm)
}
