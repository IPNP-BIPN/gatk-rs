//! Ported from `org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary` (GATK 4.6.2.0).
//!
//! The read filters are the first thing ported in this repository, on purpose. They are stateless,
//! they touch no floating point, and every tool that reads reads runs a chain of them, so the 55
//! of them are both the cheapest thing to get right and the widest-reaching. A wrong filter does
//! not produce a wrong number, it produces a different set of reads, and every number downstream
//! inherits that.
//!
//! # The part that is easy to get wrong
//!
//! A filter reads its predicate off `GATKRead`, not off the SAM flags, and the two are not the
//! same thing. `SAMRecordToGATKReadAdapter.isUnmapped` is:
//!
//! ```java
//! samRecord.getReadUnmappedFlag()
//!     || samRecord.getReferenceName() == null
//!     || samRecord.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)
//!     || samRecord.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START
//! ```
//!
//! Three criteria, of which the 0x4 flag is one. A record with the flag clear, a reference index of
//! -1 and a start of 0 is mapped by the flag and unmapped by GATK, and a port that tested the flag
//! alone would keep it. The same shape applies to `mateIsUnmapped`, and `isFirstOfPair` is
//! `isPaired() && firstOfPairFlag`, not the 0x40 flag on its own.
//!
//! Each filter below therefore records which GATKRead accessor it goes through, and the accessors
//! are implemented once, in [`read`], rather than inlined per filter.

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
pub mod read {
    use super::*;

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

    /// Whether the record carries an `RG` tag at all, which is what `getReadGroup() != null` tests.
    ///
    /// Only presence, deliberately. `SAMRecord.getReadGroup()` resolves the tag against the
    /// header's `@RG` lines and returns null when the header has no such group, so a record whose
    /// `RG` names a group the header does not declare is filtered out by GATK and kept here. That
    /// gap closes when the ported filters take a header; it is recorded rather than left implicit,
    /// because it is the kind of difference that only shows on a malformed file.
    pub fn has_read_group(read: &BamRecord) -> bool {
        read.tags.iter().any(|(tag, _)| tag.name() == *b"RG")
    }
}

/// A read filter: `true` keeps the read, exactly as `ReadFilter.test` does.
pub type ReadFilter = fn(&BamRecord) -> bool;

/// `ReadFilterLibrary.AllowAllReadsReadFilter`: "Do not filter out any read."
pub fn allow_all_reads(_read: &BamRecord) -> bool {
    true
}

/// `ReadFilterLibrary.MappedReadFilter`: filter out unmapped reads.
pub fn mapped(read: &BamRecord) -> bool {
    !read::is_unmapped(read)
}

/// `ReadFilterLibrary.MappingQualityAvailableReadFilter`: filter out MAPQ 255.
pub fn mapping_quality_available(read: &BamRecord) -> bool {
    read::mapping_quality(read) != MAPPING_QUALITY_UNAVAILABLE
}

/// `ReadFilterLibrary.MappingQualityNotZeroReadFilter`.
pub fn mapping_quality_not_zero(read: &BamRecord) -> bool {
    read::mapping_quality(read) != 0
}

/// `ReadFilterLibrary.NotDuplicateReadFilter`.
pub fn not_duplicate(read: &BamRecord) -> bool {
    !read::is_duplicate(read)
}

/// `ReadFilterLibrary.NotSecondaryAlignmentReadFilter`.
pub fn not_secondary_alignment(read: &BamRecord) -> bool {
    !read::is_secondary_alignment(read)
}

/// `ReadFilterLibrary.NotSupplementaryAlignmentReadFilter`.
pub fn not_supplementary_alignment(read: &BamRecord) -> bool {
    !read::is_supplementary_alignment(read)
}

/// `ReadFilterLibrary.PassesVendorQualityCheckReadFilter`.
pub fn passes_vendor_quality_check(read: &BamRecord) -> bool {
    !read::fails_vendor_quality_check(read)
}

/// `ReadFilterLibrary.PairedReadFilter`.
pub fn paired(read: &BamRecord) -> bool {
    read::is_paired(read)
}

/// `ReadFilterLibrary.NotProperlyPairedReadFilter`: keep reads that are *not* properly paired.
pub fn not_properly_paired(read: &BamRecord) -> bool {
    !read::is_proper_pair(read)
}

/// `ReadFilterLibrary.FirstOfPairReadFilter`.
pub fn first_of_pair(read: &BamRecord) -> bool {
    read::is_first_of_pair(read)
}

/// `ReadFilterLibrary.SecondOfPairReadFilter`.
pub fn second_of_pair(read: &BamRecord) -> bool {
    read::is_second_of_pair(read)
}

/// `ReadFilterLibrary.NonZeroFragmentLengthReadFilter`.
pub fn non_zero_fragment_length(read: &BamRecord) -> bool {
    read::fragment_length(read) != 0
}

/// `ReadFilterLibrary.MatchingBasesAndQualsReadFilter`.
pub fn matching_bases_and_quals(read: &BamRecord) -> bool {
    read::length(read) == read::base_quality_count(read)
}

/// `ReadFilterLibrary.HasReadGroupReadFilter`.
pub fn has_read_group(read: &BamRecord) -> bool {
    read::has_read_group(read)
}

/// `ReadFilterLibrary.MateDifferentStrandReadFilter`.
///
/// Keep only paired reads whose mate maps to the opposite strand, both ends mapped.
pub fn mate_different_strand(read: &BamRecord) -> bool {
    read::is_paired(read)
        && !read::is_unmapped(read)
        && !read::mate_is_unmapped(read)
        && read::mate_is_reverse_strand(read) != read::is_reverse_strand(read)
}

/// `ReadFilterLibrary.MateOnSameContigOrNoMappedMateReadFilter`.
///
/// Keep a read that is unpaired, whose mate is unmapped, or whose mate is on this read's contig.
pub fn mate_on_same_contig_or_no_mapped_mate(read: &BamRecord) -> bool {
    !read::is_paired(read)
        || read::mate_is_unmapped(read)
        || read.mate_reference_index == read.reference_index
}

/// The filters by the name GATK exposes on the command line (`--read-filter <Name>`).
pub fn by_name(name: &str) -> Option<ReadFilter> {
    Some(match name {
        "AllowAllReadsReadFilter" => allow_all_reads as ReadFilter,
        "MappedReadFilter" => mapped,
        "MappingQualityAvailableReadFilter" => mapping_quality_available,
        "MappingQualityNotZeroReadFilter" => mapping_quality_not_zero,
        "NotDuplicateReadFilter" => not_duplicate,
        "NotSecondaryAlignmentReadFilter" => not_secondary_alignment,
        "NotSupplementaryAlignmentReadFilter" => not_supplementary_alignment,
        "PassesVendorQualityCheckReadFilter" => passes_vendor_quality_check,
        "PairedReadFilter" => paired,
        "NotProperlyPairedReadFilter" => not_properly_paired,
        "FirstOfPairReadFilter" => first_of_pair,
        "SecondOfPairReadFilter" => second_of_pair,
        "NonZeroFragmentLengthReadFilter" => non_zero_fragment_length,
        "MatchingBasesAndQualsReadFilter" => matching_bases_and_quals,
        "HasReadGroupReadFilter" => has_read_group,
        "MateDifferentStrandReadFilter" => mate_different_strand,
        "MateOnSameContigOrNoMappedMateReadFilter" => mate_on_same_contig_or_no_mapped_mate,
        _ => return None,
    })
}

/// Every filter ported so far, by name. The conformance harness iterates this, so a filter that is
/// added here is exercised against the oracle without touching the harness.
pub const PORTED: &[&str] = &[
    "AllowAllReadsReadFilter",
    "FirstOfPairReadFilter",
    "HasReadGroupReadFilter",
    "MappedReadFilter",
    "MappingQualityAvailableReadFilter",
    "MappingQualityNotZeroReadFilter",
    "MatchingBasesAndQualsReadFilter",
    "MateDifferentStrandReadFilter",
    "MateOnSameContigOrNoMappedMateReadFilter",
    "NonZeroFragmentLengthReadFilter",
    "NotDuplicateReadFilter",
    "NotProperlyPairedReadFilter",
    "NotSecondaryAlignmentReadFilter",
    "NotSupplementaryAlignmentReadFilter",
    "PairedReadFilter",
    "PassesVendorQualityCheckReadFilter",
    "SecondOfPairReadFilter",
];

#[cfg(test)]
mod tests {
    use super::*;

    fn mapped_read() -> BamRecord {
        BamRecord {
            reference_index: 0,
            alignment_start: 100,
            mapping_quality: 60,
            read_bases: vec![b'A'; 10],
            base_qualities: vec![30; 10],
            ..BamRecord::default()
        }
    }

    #[test]
    fn unmapped_is_three_criteria_not_one() {
        // The flag is clear, so a port that tested the flag alone would call this mapped.
        let mut read = mapped_read();
        read.reference_index = NO_ALIGNMENT_REFERENCE_INDEX;
        assert!(read::is_unmapped(&read));
        assert!(!mapped(&read));

        let mut read = mapped_read();
        read.alignment_start = NO_ALIGNMENT_START;
        assert!(read::is_unmapped(&read));

        let mut read = mapped_read();
        read.flags |= flags::READ_UNMAPPED;
        assert!(read::is_unmapped(&read));

        assert!(mapped(&mapped_read()));
    }

    #[test]
    fn first_of_pair_requires_pairing() {
        // 0x40 without 0x1 is not first-of-pair to GATK, whatever the flag says.
        let mut read = mapped_read();
        read.flags |= flags::FIRST_OF_PAIR;
        assert!(!first_of_pair(&read));

        read.flags |= flags::READ_PAIRED;
        assert!(first_of_pair(&read));
    }

    #[test]
    fn matching_bases_and_quals_counts_absent_qualities_as_zero() {
        let mut read = mapped_read();
        assert!(matching_bases_and_quals(&read));
        read.base_qualities.clear();
        assert!(!matching_bases_and_quals(&read));
    }

    #[test]
    fn mate_on_same_contig_or_no_mapped_mate() {
        let mut read = mapped_read();
        // Unpaired: kept.
        assert!(super::mate_on_same_contig_or_no_mapped_mate(&read));

        read.flags |= flags::READ_PAIRED;
        read.mate_reference_index = 1;
        read.mate_alignment_start = 50;
        assert!(!super::mate_on_same_contig_or_no_mapped_mate(&read));

        read.mate_reference_index = 0;
        assert!(super::mate_on_same_contig_or_no_mapped_mate(&read));
    }

    #[test]
    fn every_ported_name_resolves() {
        for name in PORTED {
            assert!(by_name(name).is_some(), "{name} is listed but not wired");
        }
    }
}
