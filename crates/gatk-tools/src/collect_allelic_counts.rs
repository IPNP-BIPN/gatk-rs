//! Ported from `org.broadinstitute.hellbender.tools.copynumber.CollectAllelicCounts` and
//! `datacollection.AllelicCountCollector` (GATK 4.6.2.0).
//!
//! Reference and alternate counts at every locus of the requested intervals, which is what the
//! copy-number tools segment on.
//!
//! # Every locus is a row, except the ones that are not
//!
//! `emitEmptyLoci` is true, so a position with no reads is `0 0 <ref> N` rather than a gap. But the
//! collector returns early when the REFERENCE base is not one of `ACGT`, and that early return
//! happens after the walker has already emitted the locus -- so an `N` run in the reference is a
//! HOLE in the table, and it is the only hole there is.
//!
//! # The alternate count is not the alternate base's count
//!
//! `altReadCount = totalBaseCount - refReadCount`, and the reference's own comment says so: "we take
//! alt = total - ref instead of the actual alt count". A locus carrying three different
//! non-reference bases has all three in the alternate count, while the alternate BASE names only
//! the most common of them.
//!
//! # The base quality threshold is the collector's, not a filter's
//!
//! It defaults to 20 and is applied per PILEUP ELEMENT, so a read whose one base at this locus is
//! at quality 19 still passes every read filter and still contributes nothing here.

use gatk_engine::read_pileup::ReadPileup;
use gatk_readfilter::{self as filters, Parameterized};
use htsjdk_bam::record::BamRecord;

/// `DEFAULT_MINIMUM_MAPPING_QUALITY`, which is this tool's own.
pub const DEFAULT_MINIMUM_MAPPING_QUALITY: i32 = 30;
/// `DEFAULT_MINIMUM_BASE_QUALITY`, applied inside the collector.
pub const DEFAULT_MINIMUM_BASE_QUALITY: u8 = 20;

/// The four filters `DEFAULT_ADDITIONAL_READ_FILTERS` adds to the walker's, in its order.
pub const ADDITIONAL_READ_FILTERS: [&str; 4] = [
    "MappedReadFilter",
    "NonZeroReferenceLengthAlignmentReadFilter",
    "NotDuplicateReadFilter",
    "MappingQualityReadFilter",
];

/// `getDefaultReadFilters`: the walker's, then these four.
pub fn additional_read_filter(read: &BamRecord) -> bool {
    filters::mapped(read)
        && filters::non_zero_reference_length_alignment(read)
        && filters::not_duplicate(read)
        && Parameterized::MappingQuality {
            min: DEFAULT_MINIMUM_MAPPING_QUALITY,
            max: None,
        }
        .test(read)
}

/// `AllelicCountCollector.BASES`: the four the counter knows.
const BASES: [u8; 4] = *b"ACGT";

/// One locus's row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AllelicCount {
    /// The contig.
    pub contig: String,
    /// The position, one-based.
    pub position: i32,
    /// Reads carrying the reference base.
    pub reference_count: i32,
    /// Every other `ACGT` read, whatever base it carried.
    pub alternate_count: i32,
    /// The reference base.
    pub reference_nucleotide: u8,
    /// The most common non-reference base, or `N` when there is none.
    pub alternate_nucleotide: u8,
}

/// `collectAtLocus`, or nothing when the locus is skipped.
///
/// `None` is the reference base not being one of `ACGT`, which is a warning and a return upstream,
/// and therefore a missing row rather than an empty one.
pub fn collect_at_locus(
    reference_base: u8,
    pileup: &ReadPileup,
    contig: &str,
    position: i32,
    minimum_base_quality: u8,
) -> Option<AllelicCount> {
    if !BASES.contains(&reference_base) {
        return None;
    }

    // Deletions first, then the quality threshold: both are per element, and an `N` survives both
    // and is then not counted by either total.
    let mut counts = [0i32; 4];
    for element in &pileup.elements {
        if element.is_deletion() || element.qual() < minimum_base_quality {
            continue;
        }
        if let Some(index) = BASES.iter().position(|base| *base == element.base()) {
            counts[index] += 1;
        }
    }

    let total: i32 = counts.iter().sum();
    let reference_index = BASES
        .iter()
        .position(|base| *base == reference_base)
        .expect("the reference base is one of the four");
    let reference_count = counts[reference_index];
    // "we take alt = total - ref instead of the actual alt count".
    let alternate_count = total - reference_count;
    let alternate_nucleotide = if alternate_count == 0 {
        b'N'
    } else {
        infer_alternate(&counts, reference_index)
    };

    Some(AllelicCount {
        contig: contig.to_string(),
        position,
        reference_count,
        alternate_count,
        reference_nucleotide: reference_base,
        alternate_nucleotide,
    })
}

/// `inferAltFromPileupBaseCounts`: the most common base that is not the reference.
///
/// The reference sorts descending by count and takes the first. Java's sort is stable, so a tie is
/// broken by the order of `BASES`, which is `A C G T`.
fn infer_alternate(counts: &[i32; 4], reference_index: usize) -> u8 {
    let mut candidates: Vec<usize> = (0..4).filter(|index| *index != reference_index).collect();
    candidates.sort_by(|left, right| counts[*right].cmp(&counts[*left]));
    BASES[candidates[0]]
}

/// The table's header: a SAM header, then the column names.
///
/// The read group is the tool's own -- `ID:GATKCopyNumber` -- carrying the sample the reads named,
/// so the table is self-describing without the BAM beside it.
pub fn header(sequences: &[(String, i32)], sample: &str) -> String {
    let mut text = String::from("@HD\tVN:1.6\n");
    for (name, length) in sequences {
        text.push_str(&format!("@SQ\tSN:{name}\tLN:{length}\n"));
    }
    text.push_str(&format!("@RG\tID:GATKCopyNumber\tSM:{sample}\n"));
    text.push_str("CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n");
    text
}

/// The whole file: the header, then one row per locus that was not skipped.
pub fn write(sequences: &[(String, i32)], sample: &str, counts: &[AllelicCount]) -> String {
    let mut text = header(sequences, sample);
    for count in counts {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            count.contig,
            count.position,
            count.reference_count,
            count.alternate_count,
            count.reference_nucleotide as char,
            count.alternate_nucleotide as char
        ));
    }
    text
}
