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

/// [`interval_of_start`] by binary search: the reference's `CachedOverlapDetector` answers in
/// `O(log n)`, and a linear scan per read was `O(reads x bins)`, 35 s on a 2 Mb slice of chr20.
///
/// A contig whose intervals overlap (which `OVERLAPPING_ONLY` merging rules out) keeps the linear
/// first match, so the answer is [`interval_of_start`]'s on every input.
/// Per contig, the intervals' `(start, end, index)` sorted by start, or `None` where two overlap.
type ByContig<'a> = std::collections::HashMap<&'a str, Option<Vec<(i32, i32, usize)>>>;

pub struct Bins<'a> {
    intervals: &'a [SimpleInterval],
    by_contig: ByContig<'a>,
}

impl<'a> Bins<'a> {
    pub fn new(intervals: &'a [SimpleInterval]) -> Bins<'a> {
        let mut by_contig = ByContig::new();
        for (index, interval) in intervals.iter().enumerate() {
            if let Some(bins) = by_contig
                .entry(&interval.contig)
                .or_insert(Some(Vec::new()))
            {
                bins.push((interval.start, interval.end, index));
            }
        }
        for bins in by_contig.values_mut() {
            let sorted = bins.as_mut().map(|bins| {
                bins.sort_by_key(|&(start, _, index)| (start, index));
                bins.windows(2).all(|pair| pair[0].1 < pair[1].0)
            });
            if sorted == Some(false) {
                *bins = None;
            }
        }
        Bins {
            intervals,
            by_contig,
        }
    }

    pub fn index_of(&self, contig: &str, start: i32) -> Option<usize> {
        match self.by_contig.get(contig)? {
            Some(bins) => {
                let after = bins.partition_point(|&(begin, _, _)| begin <= start);
                let &(_, end, index) = bins.get(after.checked_sub(1)?)?;
                (start <= end).then_some(index)
            }
            None => interval_of_start(contig, start, self.intervals),
        }
    }
}

/// `onTraversalSuccess`: one count per interval, in the interval list's order.
pub fn count(reads: &[(&str, i32)], intervals: &[SimpleInterval]) -> Vec<i32> {
    let bins = Bins::new(intervals);
    let mut counts = vec![0; intervals.len()];
    for (contig, start) in reads {
        if let Some(index) = bins.index_of(contig, *start) {
            counts[index] += 1;
        }
    }
    counts
}

/// The TSV's header: a SAM header, the tool's own read group, then the three columns.
pub fn header(sequences: &[htsjdk_bam::header::SequenceRecord], sample: &str) -> String {
    let mut text = crate::model_segments_tool::header(sample, Some(sequences));
    text.push_str("CONTIG\tSTART\tEND\tCOUNT\n");
    text
}

/// The whole TSV.
pub fn write(
    sequences: &[htsjdk_bam::header::SequenceRecord],
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

/// `HDF5SimpleCountCollection.write`: the sample name, the dictionary as
/// `SAMTextHeaderCodec` encodes a header holding only it, the intervals as `HDF5Utils` writes
/// them, and the counts as a one-row matrix of doubles.
pub fn write_hdf5(
    sequences: &[htsjdk_bam::header::SequenceRecord],
    sample: &str,
    intervals: &[SimpleInterval],
    counts: &[i32],
) -> crate::hdf5_writer::Hdf5File {
    let mut file = crate::hdf5_writer::Hdf5File::new();
    file.string_array("/sample_metadata/sample_name", &[sample.to_string()]);
    let dictionary = htsjdk_bam::header::SamHeader {
        sequences: sequences.to_vec(),
        ..Default::default()
    }
    .encode();
    file.string_array("/locatable_metadata/sequence_dictionary", &[dictionary]);
    let located: Vec<crate::create_read_count_panel_of_normals::Interval> = intervals
        .iter()
        .map(
            |interval| crate::create_read_count_panel_of_normals::Interval {
                contig: interval.contig.clone(),
                start: interval.start,
                end: interval.end,
            },
        )
        .collect();
    crate::create_read_count_panel_of_normals::write_intervals(&mut file, "/intervals", &located);
    file.double_matrix(
        "/counts/values",
        &[counts.iter().map(|&count| f64::from(count)).collect()],
    );
    file
}
