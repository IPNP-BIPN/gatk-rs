//! The `GATKRead` accessors, ported from `SAMRecordToGATKReadAdapter` (GATK 4.6.2.0).
//!
//! Every walker, filter and annotation sees a read through this adapter, and three of its
//! accessors are not the flag test they look like. `isUnmapped` is:
//!
//! ```java
//! samRecord.getReadUnmappedFlag()
//!     || samRecord.getReferenceName() == null
//!     || samRecord.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)
//!     || samRecord.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START
//! ```
//!
//! Three criteria, of which the 0x4 flag is one. A record with the flag clear, a reference index
//! of -1 and a start of 0 is mapped by the flag and unmapped by GATK. The same shape applies to
//! `mateIsUnmapped`, and `isFirstOfPair` is `isPaired() && firstOfPairFlag`, not the 0x40 flag on
//! its own.
//!
//! They live in one place because a second implementation of them is a second definition of what
//! a read *is*.

use htsjdk_bam::record::BamRecord;

/// The SAM flag bits, as `SAMFlag` defines them.
pub mod flags {
    pub const READ_PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const READ_UNMAPPED: u16 = 0x4;
    pub const MATE_UNMAPPED: u16 = 0x8;
    pub const READ_REVERSE_STRAND: u16 = 0x10;
    pub const MATE_REVERSE_STRAND: u16 = 0x20;
    pub const FIRST_OF_PAIR: u16 = 0x40;
    pub const SECOND_OF_PAIR: u16 = 0x80;
    pub const NOT_PRIMARY_ALIGNMENT: u16 = 0x100;
    pub const READ_FAILS_VENDOR_QUALITY_CHECK: u16 = 0x200;
    pub const DUPLICATE_READ: u16 = 0x400;
    pub const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;
}

/// `SAMRecord.NO_ALIGNMENT_START`.
pub const NO_ALIGNMENT_START: i32 = 0;
/// The index htsjdk uses for `SAMRecord.NO_ALIGNMENT_REFERENCE_NAME` ("*").
pub const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;
/// `QualityUtils.MAPPING_QUALITY_UNAVAILABLE`.
pub const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;

/// The `GATKRead` accessors the filters are written against.
///
/// They live here rather than inside each filter because their definitions are where the
/// divergences hide: three of them are not the flag test they look like.
pub fn is_paired(read: &BamRecord) -> bool {
    read.flags & flags::READ_PAIRED != 0
}

/// `SAMRecordToGATKReadAdapter.isUnmapped`: the flag, an absent reference, or a zero start.
pub fn is_unmapped(read: &BamRecord) -> bool {
    read.flags & flags::READ_UNMAPPED != 0
        || read.reference_index == NO_ALIGNMENT_REFERENCE_INDEX
        || read.alignment_start == NO_ALIGNMENT_START
}

/// `SAMRecordToGATKReadAdapter.mateIsUnmapped`, same three criteria applied to the mate.
///
/// The Java asserts `isPaired()` first and throws otherwise; here that is the caller's job, and
/// every caller in this crate is a filter that has already tested pairing.
pub fn mate_is_unmapped(read: &BamRecord) -> bool {
    read.flags & flags::MATE_UNMAPPED != 0
        || read.mate_reference_index == NO_ALIGNMENT_REFERENCE_INDEX
        || read.mate_alignment_start == NO_ALIGNMENT_START
}

/// `isPaired() && firstOfPairFlag`, not the 0x40 flag alone.
pub fn is_first_of_pair(read: &BamRecord) -> bool {
    is_paired(read) && read.flags & flags::FIRST_OF_PAIR != 0
}

/// `isPaired() && secondOfPairFlag`.
pub fn is_second_of_pair(read: &BamRecord) -> bool {
    is_paired(read) && read.flags & flags::SECOND_OF_PAIR != 0
}

pub fn is_proper_pair(read: &BamRecord) -> bool {
    is_paired(read) && read.flags & flags::PROPER_PAIR != 0
}

pub fn is_duplicate(read: &BamRecord) -> bool {
    read.flags & flags::DUPLICATE_READ != 0
}

pub fn is_secondary_alignment(read: &BamRecord) -> bool {
    read.flags & flags::NOT_PRIMARY_ALIGNMENT != 0
}

pub fn is_supplementary_alignment(read: &BamRecord) -> bool {
    read.flags & flags::SUPPLEMENTARY_ALIGNMENT != 0
}

pub fn fails_vendor_quality_check(read: &BamRecord) -> bool {
    read.flags & flags::READ_FAILS_VENDOR_QUALITY_CHECK != 0
}

pub fn is_reverse_strand(read: &BamRecord) -> bool {
    read.flags & flags::READ_REVERSE_STRAND != 0
}

pub fn mate_is_reverse_strand(read: &BamRecord) -> bool {
    read.flags & flags::MATE_REVERSE_STRAND != 0
}

pub fn fragment_length(read: &BamRecord) -> i32 {
    read.inferred_insert_size
}

/// `GATKRead.getLength()`: the number of bases, which is what the record carries.
pub fn length(read: &BamRecord) -> usize {
    read.read_bases.len()
}

/// `GATKRead.getBaseQualityCount()`.
///
/// htsjdk encodes "qualities absent" as a run of 0xFF of the read's length, and htsjdk-rs
/// represents that as an empty vector, so an absent quality array counts as zero here exactly
/// as `SAMRecord.getBaseQualities()` returns an empty array for it.
pub fn base_quality_count(read: &BamRecord) -> usize {
    read.base_qualities.len()
}

pub fn mapping_quality(read: &BamRecord) -> u8 {
    read.mapping_quality
}

/// Whether the record carries an `RG` tag at all.
///
/// Presence is the whole semantics, and this was worth checking rather than assuming: an
/// earlier version of this comment claimed the reference resolves the tag against the header's
/// `@RG` lines and drops a read naming an undeclared group. It does not.
/// `SAMRecordToGATKReadAdapter.getReadGroup` returns the raw string attribute, so a read whose
/// `RG` is `rg_absent` passes `HasReadGroupReadFilter` on a header that has never heard of it.
/// The conformance corpus now carries exactly that record, and the reference keeps it.
pub fn has_read_group(read: &BamRecord) -> bool {
    read.tags.iter().any(|(tag, _)| tag.name() == *b"RG")
}
