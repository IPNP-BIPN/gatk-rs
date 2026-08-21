//! `SetNmMdAndUqTags`, ported from `picard.sam.SetNmMdAndUqTags` and the two `fix` methods of
//! `picard.sam.AbstractAlignmentMerger` (Picard 3.4.0).
//!
//! NM, MD and UQ recalculated against the reference for every mapped record. `SetNmAndUqTags` is a
//! deprecated subclass that adds nothing, so this is both tools.
//!
//! # Two NM tags, computed by two different functions
//!
//! ```java
//! SequenceUtil.calculateMdAndNmTags(record, referenceBases, true, !isBisulfiteSequence);
//! if (isBisulfiteSequence) {
//!     record.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(record, referenceBases, 0, isBisulfiteSequence));
//! }
//! ```
//!
//! For an ordinary read one walk writes both MD and NM. For a bisulfite read the NM is written
//! again by `calculateSamNmTag`, which forgives a C read as a T, while the MD from the first walk
//! is left as it was. The golden shows the consequence in a single record: eight Ts over eight
//! reference Cs come out `MD:Z:0C0C0C0C0C0C0C0C0` with `NM:i:0`. A port that computed one number
//! and wrote it twice would disagree with the reference on exactly this record.
//!
//! # UQ needs qualities, and it is set even when it is zero
//!
//! `fixUq` returns without touching anything when the read has no qualities, so such a record
//! keeps whatever UQ it arrived with. Otherwise the tag is always written, `UQ:i:0` included.
//!
//! # An unmapped read is not touched at all
//!
//! `fixRecord` checks `getReadUnmappedFlag()` first, so a wrong NM, MD or UQ on an unmapped read
//! survives the run.

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::md_nm::calculate_md_and_nm;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::{calculate_sam_nm_tag, sum_qualities_of_mismatches};
use htsjdk_bam::tag::{Tag, TagValue};

/// `SAMFlag.READ_UNMAPPED`.
const READ_UNMAPPED: u16 = 0x4;
/// `SAMFlag.READ_REVERSE_STRAND`, which only the bisulfite comparison reads.
const READ_REVERSE_STRAND: u16 = 0x10;

/// The arguments the run reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Arguments {
    pub is_bisulfite_sequence: bool,
    pub set_only_uq: bool,
}

/// What the tool refuses before it reads a record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SetTagsError {
    /// The header's `SO` is anything but `coordinate`. The message names what was found, and an
    /// absent `SO` is `unsorted`.
    NotCoordinateSorted { found: String },
}

impl SetTagsError {
    pub fn java_class(&self) -> &str {
        "htsjdk.samtools.SAMException"
    }

    pub fn message(&self) -> String {
        match self {
            SetTagsError::NotCoordinateSorted { found } => {
                format!("Input must be coordinate-sorted for this program to run. Found: {found}")
            }
        }
    }
}

/// `fixRecord`, for one record and the bases of the contig it sits on.
///
/// The reference bases are the whole contig, which is what `ReferenceSequenceFileWalker.get`
/// hands over, so the offset is always zero.
pub fn fix_record(record: &mut BamRecord, reference_bases: &[u8], arguments: &Arguments) {
    if record.flags & READ_UNMAPPED != 0 {
        return;
    }
    if !arguments.set_only_uq {
        fix_nm_md(record, reference_bases, arguments.is_bisulfite_sequence);
    }
    fix_uq(record, reference_bases, arguments.is_bisulfite_sequence);
}

/// `AbstractAlignmentMerger.fixNmMdAndUq`, minus the `fixUq` it ends with.
fn fix_nm_md(record: &mut BamRecord, reference_bases: &[u8], is_bisulfite_sequence: bool) {
    // `calculateMdAndNmTags(record, ref, true, !isBisulfiteSequence)`: the NM this walk produces
    // is only kept when the read is not bisulfite treated.
    let (md, nm) = calculate_md_and_nm(
        record.alignment_start,
        &record.cigar,
        &record.read_bases,
        reference_bases,
    );
    record.tags.insert(Tag::new(b"MD"), TagValue::Str(md));
    if is_bisulfite_sequence {
        // The second function, which forgives the conversion and disagrees with the MD above.
        let blocks = alignment_blocks(&record.cigar, record.alignment_start);
        let nm = calculate_sam_nm_tag(
            &record.read_bases,
            &blocks,
            &record.cigar,
            reference_bases,
            0,
            record.flags & READ_REVERSE_STRAND != 0,
            true,
        );
        record
            .tags
            .insert(Tag::new(b"NM"), TagValue::Int(i64::from(nm)));
    } else {
        record
            .tags
            .insert(Tag::new(b"NM"), TagValue::Int(i64::from(nm)));
    }
}

/// `AbstractAlignmentMerger.fixUq`.
fn fix_uq(record: &mut BamRecord, reference_bases: &[u8], is_bisulfite_sequence: bool) {
    // `record.getBaseQualities() != SAMRecord.NULL_QUALS`, which here is an empty vector.
    if record.base_qualities.is_empty() {
        return;
    }
    let blocks = alignment_blocks(&record.cigar, record.alignment_start);
    let qualities = sum_qualities_of_mismatches(
        &record.read_bases,
        &record.base_qualities,
        &blocks,
        record.alignment_start,
        reference_bases,
        0,
        record.flags & READ_REVERSE_STRAND != 0,
        is_bisulfite_sequence,
    )
    .expect("the offset is zero and an alignment start is at least one");
    record
        .tags
        .insert(Tag::new(b"UQ"), TagValue::Int(i64::from(qualities)));
}

/// `doWork()`'s sort-order check, which reads the header's claim and nothing else.
pub fn check_sort_order(sort_order: Option<&str>) -> Result<(), SetTagsError> {
    let found = sort_order.unwrap_or("unsorted");
    if found == "coordinate" {
        Ok(())
    } else {
        Err(SetTagsError::NotCoordinateSorted {
            found: found.to_string(),
        })
    }
}
