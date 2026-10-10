//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils`
//! (`finalizeRegion`, `cleanOverlappingReadPairs`), with `FragmentCollection.create` and
//! `FragmentUtils.adjustQualsOfOverlappingPairedFragments` (GATK 4.6.2.0): the reads of an assembly
//! region as the assembler and the genotyper see them.
//!
//! Each read is first rid of its soft clips (hard-clipped, or reverted when soft-clipped bases
//! are used and the fragment size is well defined, the original soft span then kept in the `os`
//! and `oe` tags), then of its low-quality ends (soft- or hard-clipped below the minimum tail
//! quality, 6 when reads are error-corrected), then of its adaptor, then of whatever lies outside
//! the padded span. A read left empty, unmapped or with no aligned base is dropped. The survivors
//! are sorted by coordinate, and the two mates of a fragment that overlap have the qualities of
//! their shared bases capped at half the PCR error quality (20), or zeroed where the bases differ.
//!
//! The second list, of reads hard-clipped whatever the soft-clip setting, is kept only when asked
//! for (`trackHardclippedReads`) and becomes the region's hard-clipped pileup reads.

use htsjdk_bam::cigar::Op;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use crate::assembly_region::{AssemblyRegion, RegionError};
use crate::clipping::{
    clip_low_qual_ends, hard_clip_adaptor_sequence, hard_clip_soft_clipped_bases,
    hard_clip_to_region, revert_soft_clipped_bases, ClipError, ClippingRepresentation,
};
use crate::read;
use crate::read_pileup::sample_name;
use crate::read_utils::{self, READ_INDEX_NOT_FOUND};

/// `HaplotypeCallerEngine.MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION`.
pub const MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION: u8 = 6;
/// `FragmentUtils.HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL`: the Phred score of a 1e-4 PCR error rate,
/// 40, halved.
pub const HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL: u8 = 20;
/// `ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG`.
pub const ORIGINAL_SOFTCLIP_START_TAG: &[u8; 2] = b"os";
/// `ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG`.
pub const ORIGINAL_SOFTCLIP_END_TAG: &[u8; 2] = b"oe";

/// `finalizeRegion`'s switches.
#[derive(Debug, Clone)]
pub struct FinalizeArguments {
    pub error_correct_reads: bool,
    pub dont_use_soft_clipped_bases: bool,
    pub min_tail_quality: u8,
    pub correct_overlapping_base_qualities: bool,
    pub soft_clip_low_quality_ends: bool,
    pub override_softclip_fragment_check: bool,
    pub track_hardclipped_reads: bool,
}

/// What `finalizeRegion` refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum FinalizeError {
    Clip(ClipError),
    Region(RegionError),
    /// `NullPointerException`: a read whose sample is not in the sample list.
    UnknownSample(Option<String>),
    /// `IllegalArgumentException` from `FragmentCollection.create`: reads out of start order.
    OutOfOrder(String),
    /// `ArrayIndexOutOfBoundsException` in the overlap adjustment.
    IndexOutOfBounds,
}

impl From<ClipError> for FinalizeError {
    fn from(e: ClipError) -> Self {
        FinalizeError::Clip(e)
    }
}

/// `finalizeRegion`: a region already finalized is left alone.
pub fn finalize_region(
    region: &mut AssemblyRegion,
    arguments: &FinalizeArguments,
    header: &SamHeader,
    samples: &[String],
) -> Result<(), FinalizeError> {
    if region.is_finalized() {
        return Ok(());
    }
    let min_tail_quality = if arguments.error_correct_reads {
        MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION
    } else {
        arguments.min_tail_quality
    };
    let tail = if arguments.soft_clip_low_quality_ends {
        ClippingRepresentation::SoftclipBases
    } else {
        ClippingRepresentation::HardclipBases
    };
    let mut reads_to_use = Vec::new();
    let mut hard_clipped_reads_to_use = Vec::new();
    for original in region.reads() {
        let unclipped = if arguments.dont_use_soft_clipped_bases
            || !(arguments.override_softclip_fragment_check
                || read_utils::has_well_defined_fragment_size(original))
        {
            hard_clip_soft_clipped_bases(original, Some(header), 0)?
        } else {
            revert_soft_clipped_bases_keeping_span(original, header)?
        };
        let read = clip_low_qual_ends(&unclipped, Some(header), min_tail_quality, tail)?;
        if let Some(kept) = hard_clip_and_possibly_keep(region, &read, header)? {
            reads_to_use.push(kept);
        }
        if arguments.track_hardclipped_reads {
            let hard = clip_low_qual_ends(
                &hard_clip_soft_clipped_bases(original, Some(header), 0)?,
                Some(header),
                min_tail_quality,
                ClippingRepresentation::HardclipBases,
            )?;
            if let Some(kept) = hard_clip_and_possibly_keep(region, &hard, header)? {
                hard_clipped_reads_to_use.push(kept);
            }
        }
    }
    reads_to_use.sort_by(read_utils::compare_read_coordinate);
    hard_clipped_reads_to_use.sort_by(read_utils::compare_read_coordinate);
    if arguments.correct_overlapping_base_qualities {
        clean_overlapping_read_pairs(&mut reads_to_use, samples, header, true, None)?;
        clean_overlapping_read_pairs(&mut hard_clipped_reads_to_use, samples, header, true, None)?;
    }
    region.clear_reads();
    region
        .add_all(reads_to_use, header)
        .map_err(FinalizeError::Region)?;
    region
        .add_hard_clipped_pileup_reads(hard_clipped_reads_to_use, header)
        .map_err(FinalizeError::Region)?;
    region.set_finalized(true);
    Ok(())
}

/// `AssemblyBasedCallerUtils.revertSoftClippedBases`: the reverted read, with the soft span it had
/// in `os` and `oe`.
fn revert_soft_clipped_bases_keeping_span(
    read: &BamRecord,
    header: &SamHeader,
) -> Result<BamRecord, FinalizeError> {
    let soft_start = read_utils::start(read);
    let soft_end = read_utils::end(read);
    let mut result = revert_soft_clipped_bases(read, Some(header))?;
    result.tags.insert(
        Tag::new(ORIGINAL_SOFTCLIP_START_TAG),
        TagValue::Int(soft_start.into()),
    );
    result.tags.insert(
        Tag::new(ORIGINAL_SOFTCLIP_END_TAG),
        TagValue::Int(soft_end.into()),
    );
    Ok(result)
}

/// `HardClipAndPossiblyAddToCollection`: the read rid of its adaptor and cut to the padded span,
/// if anything of it is left there.
fn hard_clip_and_possibly_keep(
    region: &AssemblyRegion,
    read: &BamRecord,
    header: &SamHeader,
) -> Result<Option<BamRecord>, FinalizeError> {
    if read_utils::start(read) > read_utils::end(read) || read::is_unmapped(read) {
        return Ok(None);
    }
    let adaptor_clipped = hard_clip_adaptor_sequence(read, Some(header))?;
    if adaptor_clipped.read_bases.is_empty() || adaptor_clipped.cigar.read_length() == 0 {
        return Ok(None);
    }
    let padded = region.padded_span();
    let clipped = hard_clip_to_region(&adaptor_clipped, Some(header), padded.start, padded.end)?;
    let contig = header
        .sequences
        .get(adaptor_clipped.reference_index as usize)
        .map(|s| s.name.as_str())
        .unwrap_or("*");
    let overlaps = padded.overlaps(
        contig,
        read_utils::start(&adaptor_clipped),
        read_utils::end(&adaptor_clipped),
    );
    if read_utils::start(&clipped) <= read_utils::end(&clipped)
        && !clipped.read_bases.is_empty()
        && overlaps
    {
        Ok(Some(clipped))
    } else {
        Ok(None)
    }
}

/// `cleanOverlappingReadPairs`: per sample, every overlapping pair `FragmentCollection.create`
/// finds is adjusted. The qualities are changed in place.
pub fn clean_overlapping_read_pairs(
    reads: &mut [BamRecord],
    samples: &[String],
    header: &SamHeader,
    set_conflicting_to_zero: bool,
    half_of_pcr_snv_qual: Option<u8>,
) -> Result<(), FinalizeError> {
    clean_overlapping_read_pairs_with_indels(
        reads,
        samples,
        header,
        set_conflicting_to_zero,
        half_of_pcr_snv_qual,
        None,
    )
}

/// The same with `halfOfPcrIndelQual`, which Mutect2 passes: the overlapping bases' insertion and
/// deletion qualities (`BI` and `BD`, 45 where a read carries none) are capped at it too.
pub fn clean_overlapping_read_pairs_with_indels(
    reads: &mut [BamRecord],
    samples: &[String],
    header: &SamHeader,
    set_conflicting_to_zero: bool,
    half_of_pcr_snv_qual: Option<u8>,
    half_of_pcr_indel_qual: Option<u8>,
) -> Result<(), FinalizeError> {
    let mut by_sample: Vec<Vec<usize>> = vec![Vec::new(); samples.len()];
    for (index, read) in reads.iter().enumerate() {
        let sample = sample_name(read, header);
        let slot = samples
            .iter()
            .position(|s| Some(s) == sample.as_ref())
            .ok_or(FinalizeError::UnknownSample(sample))?;
        by_sample[slot].push(index);
    }
    for indices in by_sample {
        for (first, second) in overlapping_pairs(reads, &indices)? {
            adjust_quals_of_overlapping_paired_fragments(
                reads,
                first,
                second,
                set_conflicting_to_zero,
                half_of_pcr_snv_qual,
                half_of_pcr_indel_qual,
            )?;
        }
    }
    Ok(())
}

/// `FragmentCollection.create(...).getOverlappingPairs()`: a paired read whose mate starts within
/// it waits for that mate, and the two make a pair in the order they arrived.
fn overlapping_pairs(
    reads: &[BamRecord],
    indices: &[usize],
) -> Result<Vec<(usize, usize)>, FinalizeError> {
    let mut pairs = Vec::new();
    let mut waiting: Vec<(String, usize)> = Vec::new();
    let mut last_start = -1;
    for &i in indices {
        let read = &reads[i];
        let start = read_utils::start(read);
        if start < last_start {
            return Err(FinalizeError::OutOfOrder(read.read_name.clone()));
        }
        last_start = start;
        let singleton = !read::is_paired(read)
            || read::mate_is_unmapped(read)
            || read_utils::mate_start(read) == 0
            || read_utils::mate_start(read) > read_utils::end(read);
        if singleton {
            continue;
        }
        match waiting.iter().position(|(name, _)| name == &read.read_name) {
            Some(position) => {
                let (_, first) = waiting.remove(position);
                pairs.push((first, i));
            }
            None => waiting.push((read.read_name.clone(), i)),
        }
    }
    Ok(pairs)
}

/// `FragmentUtils.adjustQualsOfOverlappingPairedFragments`.
fn adjust_quals_of_overlapping_paired_fragments(
    reads: &mut [BamRecord],
    left: usize,
    right: usize,
    set_conflicting_to_zero: bool,
    half_of_pcr_snv_qual: Option<u8>,
    half_of_pcr_indel_qual: Option<u8>,
) -> Result<(), FinalizeError> {
    let in_order = read_utils::soft_start(&reads[left]) < read_utils::soft_start(&reads[right]);
    let (first, second) = if in_order {
        (left, right)
    } else {
        (right, left)
    };
    let (a, b) = (&reads[first], &reads[second]);
    if read_utils::end(a) < read_utils::start(b) || a.reference_index != b.reference_index {
        return Ok(());
    }
    let (offset, operator) = read_utils::read_index_for_read(a, read_utils::start(b));
    if offset == READ_INDEX_NOT_FOUND || matches!(operator, Some(Op::S | Op::H)) {
        return Ok(());
    }
    let first_end_base = read_utils::read_index_for_read(a, read_utils::end(a)).0;
    let second_end_base = read_utils::read_index_for_read(b, read_utils::end(b)).0;
    let second_offset = read_utils::read_index_for_read(b, read_utils::start(b)).0;
    let overlapping = (first_end_base - offset).min(second_end_base - second_offset) + 1;
    let half = half_of_pcr_snv_qual.unwrap_or(HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL);
    let mut first_quals = a.base_qualities.clone();
    let mut second_quals = b.base_qualities.clone();
    for i in 0..overlapping.max(0) {
        let fi = (offset + i) as usize;
        let si = (second_offset + i) as usize;
        let (Some(&fb), Some(&sb)) = (a.read_bases.get(fi), b.read_bases.get(si)) else {
            return Err(FinalizeError::IndexOutOfBounds);
        };
        if fi >= first_quals.len() || si >= second_quals.len() {
            return Err(FinalizeError::IndexOutOfBounds);
        }
        if fb == sb {
            first_quals[fi] = first_quals[fi].min(half);
            second_quals[si] = second_quals[si].min(half);
        } else if set_conflicting_to_zero {
            first_quals[fi] = 0;
            second_quals[si] = 0;
        }
    }
    reads[first].base_qualities = first_quals;
    reads[second].base_qualities = second_quals;
    if let Some(max_indel) = half_of_pcr_indel_qual {
        for tag in [b"BD", b"BI"] {
            let mut first_indel = indel_qualities(&reads[first], tag);
            let mut second_indel = indel_qualities(&reads[second], tag);
            for i in 0..overlapping.max(0) {
                let fi = (offset + i) as usize;
                let si = (second_offset + i) as usize;
                if fi >= first_indel.len() || si >= second_indel.len() {
                    return Err(FinalizeError::IndexOutOfBounds);
                }
                first_indel[fi] = first_indel[fi].min(max_indel);
                second_indel[si] = second_indel[si].min(max_indel);
            }
            set_indel_qualities(&mut reads[first], tag, &first_indel);
            set_indel_qualities(&mut reads[second], tag, &second_indel);
        }
    }
    Ok(())
}

/// `ReadUtils.getBaseInsertionQualities` and `getBaseDeletionQualities`: the tag's phred+33
/// string, or `DEFAULT_INSERTION_DELETION_QUAL` (45) at every base when the read has none.
fn indel_qualities(read: &BamRecord, tag: &[u8; 2]) -> Vec<u8> {
    match read.tags.get(htsjdk_bam::tag::Tag::new(tag)) {
        Some(htsjdk_bam::tag::TagValue::Str(text)) => {
            text.bytes().map(|b| b.wrapping_sub(33)).collect()
        }
        _ => vec![45; read.read_bases.len()],
    }
}

/// `ReadUtils.setInsertionBaseQualities` and `setDeletionBaseQualities`: the qualities written
/// back as a phred+33 string.
fn set_indel_qualities(read: &mut BamRecord, tag: &[u8; 2], quals: &[u8]) {
    let text: String = quals.iter().map(|&q| (q + 33) as char).collect();
    read.tags.insert(
        htsjdk_bam::tag::Tag::new(tag),
        htsjdk_bam::tag::TagValue::Str(text),
    );
}
