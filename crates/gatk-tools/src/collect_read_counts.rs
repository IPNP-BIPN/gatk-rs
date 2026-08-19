//! Ported from `org.broadinstitute.hellbender.tools.copynumber.CollectReadCounts`
//! (GATK 4.6.2.0).
//!
//! Read counts per interval. The tool's whole body is one lookup per read, and the lookup is the
//! behaviour worth porting carefully.
//!
//! # A read is counted by its start
//!
//! `intervalCachedOverlapDetector.getOverlap(new SimpleInterval(contig, read.getStart(),
//! read.getStart()))`: the query is ONE BASE WIDE, at the read's start. So a read spanning three
//! intervals is counted in exactly one of them, and a read that starts before an interval and
//! covers the whole of it is counted in NEITHER -- that interval's row is a zero even though every
//! base of it was read.
//!
//! # The detector holds one contig at a time
//!
//! It is rebuilt whenever the contig changes, from the intervals on THAT contig, which is why a
//! read can never match an interval elsewhere. With a coordinate-sorted BAM that is one rebuild per
//! contig; with an unsorted one it would be a rebuild per switch, and the counts would be the same.
//!
//! # Every requested interval is a row
//!
//! The output is built from the interval list rather than from the counts, so an interval nothing
//! started in is a zero and the row order is the interval argument's -- which is sorted, whatever
//! order the user gave.

use gatk_engine::interval::SimpleInterval;
use gatk_readfilter::{self as filters, with_header, Parameterized};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// `DEFAULT_MINIMUM_MAPPING_QUALITY`, which is this tool's own and not the library's.
pub const DEFAULT_MINIMUM_MAPPING_QUALITY: i32 = 30;

/// The four filters this tool adds to the walker's two, in its order.
pub const ADDITIONAL_READ_FILTERS: [&str; 4] = [
    "MappedReadFilter",
    "NonZeroReferenceLengthAlignmentReadFilter",
    "NotDuplicateReadFilter",
    "MappingQualityReadFilter",
];

/// `getDefaultReadFilters`: the walker's wellformed and mapped, then these four.
pub fn default_read_filter(read: &BamRecord, header: &SamHeader) -> bool {
    with_header::wellformed(read, header)
        && filters::mapped(read)
        && filters::non_zero_reference_length_alignment(read)
        && filters::not_duplicate(read)
        && Parameterized::MappingQuality {
            min: DEFAULT_MINIMUM_MAPPING_QUALITY,
            max: None,
        }
        .test(read)
}

/// `apply`, for one read: the index of the interval its START falls in, if any.
///
/// `intervals` is the list for the read's contig. The reference's detector answers ONE interval; the
/// intervals it is built from do not overlap, because the tool requires `OVERLAPPING_ONLY` merging,
/// so the first match is the only match.
pub fn interval_of_start(
    read_contig: &str,
    read_start: i32,
    intervals: &[SimpleInterval],
) -> Option<usize> {
    intervals.iter().position(|interval| {
        interval.contig == read_contig && interval.start <= read_start && read_start <= interval.end
    })
}

/// `onTraversalSuccess`: one count per interval, in the interval list's order.
pub fn count(reads: &[(&str, i32)], intervals: &[SimpleInterval]) -> Vec<i32> {
    let mut counts = vec![0; intervals.len()];
    for (contig, start) in reads {
        if let Some(index) = interval_of_start(contig, *start, intervals) {
            counts[index] += 1;
        }
    }
    counts
}

/// The TSV's header: a SAM header, the tool's own read group, then the three columns.
pub fn header(sequences: &[(String, i32)], sample: &str) -> String {
    let mut text = String::from("@HD\tVN:1.6\n");
    for (name, length) in sequences {
        text.push_str(&format!("@SQ\tSN:{name}\tLN:{length}\n"));
    }
    text.push_str(&format!("@RG\tID:GATKCopyNumber\tSM:{sample}\n"));
    text.push_str("CONTIG\tSTART\tEND\tCOUNT\n");
    text
}

/// The whole TSV.
pub fn write(
    sequences: &[(String, i32)],
    sample: &str,
    intervals: &[SimpleInterval],
    counts: &[i32],
) -> String {
    let mut text = header(sequences, sample);
    for (interval, count) in intervals.iter().zip(counts) {
        text.push_str(&format!(
            "{}\t{}\t{}\t{count}\n",
            interval.contig, interval.start, interval.end
        ));
    }
    text
}
