//! Ported from `org.broadinstitute.hellbender.tools.PrintDistantMates` (GATK 4.6.2.0).
//!
//! The sixth whole tool of the record-transform archetype, and the first that **moves** a read
//! rather than editing it where it lies: the read is rewritten into its mate's position, unmapped,
//! with its old alignment kept in an `OA` tag and a `DM` tag saying it was moved.
//!
//! ```java
//! public void apply(GATKRead read, ReferenceContext r, FeatureContext f) {
//!     final GATKRead copy = doDistantMateAlterations(read);
//!     outputWriter.addRead(copy);
//! }
//! ```
//!
//! # A third pattern for the default read filters
//!
//! ```java
//! final List<ReadFilter> readFilters = new ArrayList<>(super.getDefaultReadFilters());
//! readFilters.add(ReadFilterLibrary.PAIRED);
//! ...
//! ```
//!
//! `PrintReads` takes `GATKTool`'s default. `UnmarkDuplicates` and `RevertBaseQualityScores`
//! replace it. This one **extends** it, so the measured list is `WellformedReadFilter`, then
//! `PairedReadFilter`, `PrimaryLineReadFilter`, `NotDuplicateReadFilter` and
//! `MateDistantReadFilter`, in that order. Three patterns inside one archetype, which is the kind
//! of difference an archetype hides: see [`default_read_filter`].
//!
//! # The same tag, spelled two ways
//!
//! [`crate::add_original_alignment_tags`] builds the same six-field `OA` string, and the two tools
//! disagree about a missing `NM`: this one writes **nothing** between the last comma and the
//! semicolon, the other writes the four characters `null`. Measured on one read of one fixture,
//! run through both tools: `chr1,200,+,10M,60,;` here against `chr1,200,+,10M,60,null;` there. A
//! reader of the tag cannot tell which tool wrote it, and only one of the two spellings parses
//! back.
//!
//! # The writer is not told the reads are sorted
//!
//! `createSAMWriter(output, false)` where every other tool of the archetype passes `true`. That is
//! not a formality here: the transform moves each read to its mate, so the traversal order is not
//! the output order. Measured, three reads leaving the traversal for `chr2:600`, `chr1:2500` and
//! `chr2:150` are written `chr1:2500`, `chr2:150`, `chr2:600`: htsjdk re-sorts them by
//! [`htsjdk_bam::coordinate::compare`], and the index is the index of that order.
//!
//! # `undoDistantMateAlterations` is the inverse, and one thing escapes it
//!
//! Measured as a round trip against each fixture record's own SAM text, tag block included, it
//! comes back identical: the `NM` cleared into the `OA` and set again lands back in its old place,
//! because htsjdk keeps the tag list sorted by the packed tag rather than by insertion order.
//!
//! Its own guard is narrower than it looks. Three malformed `OA` values raise
//! `UserException: can't recover alignment from OA tag: ...` as the `catch` promises, but a value
//! naming a contig the dictionary does not hold does **not**: `setReferenceName` stores the name
//! and resolves it lazily, so the failure is an `IllegalArgumentException` thrown out of somewhere
//! else entirely, `Reference index for 'chr9' not found in sequence dictionary.` This port cannot
//! defer it (a record here holds an index, not a name), so it reports it as
//! [`UndoError::UnknownContig`], with the reference's own message.

use std::cmp::Ordering;

use htsjdk_bam::coordinate;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text_parse::parse_cigar;

use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_readfilter::{not_duplicate, paired, primary_line, with_header, Parameterized};

use crate::sam_output::{header_for_sam_writer, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK PrintDistantMates";

/// `DISTANT_MATE_TAG`, whose value is always the empty string: the tag's presence is the message.
pub const DISTANT_MATE_TAG: &[u8; 2] = b"DM";
/// `SAMTag.OA`.
pub const ORIGINAL_ALIGNMENT_TAG: &[u8; 2] = b"OA";
/// `SAMTag.NM`, which this tool moves out of the record and into the `OA`.
pub const EDIT_DISTANCE_TAG: &[u8; 2] = b"NM";

/// `MateDistantReadFilter.DEFAULT_MATE_TOO_DISTANT_THRESHOLD`.
pub const DEFAULT_MATE_TOO_DISTANT_THRESHOLD: i32 = 1000;

/// `SAMFlag.READ_UNMAPPED`.
const READ_UNMAPPED: u16 = 0x4;
/// `SAMFlag.READ_REVERSE_STRAND`.
const READ_REVERSE_STRAND: u16 = 0x10;

/// `getDefaultReadFilters()`, in the order the list is built.
///
/// The first entry is `super.getDefaultReadFilters()`; the other four are this tool's additions.
pub const DEFAULT_READ_FILTERS: [&str; 5] = [
    "WellformedReadFilter",
    "PairedReadFilter",
    "PrimaryLineReadFilter",
    "NotDuplicateReadFilter",
    "MateDistantReadFilter",
];

/// The conjunction [`DEFAULT_READ_FILTERS`] names, which is what a run with no `--read-filter`
/// applies.
///
/// `mate_too_distant_length` is `MateDistantReadFilter`'s one argument; the default is
/// [`DEFAULT_MATE_TOO_DISTANT_THRESHOLD`].
pub fn default_read_filter(
    read: &BamRecord,
    header: &SamHeader,
    mate_too_distant_length: i32,
) -> bool {
    with_header::wellformed(read, header)
        && paired(read)
        && primary_line(read)
        && not_duplicate(read)
        && Parameterized::MateDistant {
            threshold: mate_too_distant_length,
        }
        .test(read)
}

/// `isDistantMate`: the `DM` tag is present, whatever it holds.
pub fn is_distant_mate(read: &BamRecord) -> bool {
    read.tags.get(Tag::new(DISTANT_MATE_TAG)).is_some()
}

/// The `OA` value this tool writes: `contig,start,strand,cigar,mapq,NM;`.
///
/// A missing `NM` leaves the field **empty**, which is where this tool and
/// [`crate::add_original_alignment_tags::original_alignment_value`] part company.
pub fn original_alignment_value(read: &BamRecord, header: &SamHeader) -> String {
    let contig = sequence_name(header, read.reference_index);
    let strand = if read.flags & READ_REVERSE_STRAND != 0 {
        "-"
    } else {
        "+"
    };
    let edit_distance = match read.tags.get(Tag::new(EDIT_DISTANCE_TAG)) {
        Some(TagValue::Int(value)) => value.to_string(),
        Some(other) => format!("{other:?}"),
        None => String::new(),
    };
    format!(
        "{contig},{},{strand},{},{},{edit_distance};",
        read.alignment_start,
        read.cigar.to_text(),
        read.mapping_quality
    )
}

/// The name of a reference by its index, which is what `getContig` resolves through the header.
fn sequence_name(header: &SamHeader, index: i32) -> String {
    usize::try_from(index)
        .ok()
        .and_then(|at| header.sequences.get(at))
        .map(|sequence| sequence.name.clone())
        .unwrap_or_default()
}

/// The index of a reference by name, which is what `setReferenceName` defers and htsjdk later
/// resolves.
fn sequence_index(header: &SamHeader, name: &str) -> Option<i32> {
    header
        .sequences
        .iter()
        .position(|sequence| sequence.name == name)
        .and_then(|at| i32::try_from(at).ok())
}

/// `doDistantMateAlterations`: the read, moved onto its mate and unmapped.
///
/// The `OA` is built from the read as it stands, before anything is cleared, and the copy is the
/// thing that changes: the reference returns a copy and leaves its argument alone.
pub fn do_distant_mate_alterations(read: &BamRecord, header: &SamHeader) -> BamRecord {
    let mut copy = read.clone();
    let original = original_alignment_value(read, header);
    copy.tags.remove(Tag::new(EDIT_DISTANCE_TAG));
    copy.tags
        .insert(Tag::new(ORIGINAL_ALIGNMENT_TAG), TagValue::Str(original));
    // The tag's value is the empty string: `setAttribute(DISTANT_MATE_TAG, "")`.
    copy.tags
        .insert(Tag::new(DISTANT_MATE_TAG), TagValue::Str(String::new()));
    // `setPosition` clears the unmapped flag and `setIsUnmapped` sets it again, in that order.
    copy.reference_index = read.mate_reference_index;
    copy.alignment_start = read.mate_alignment_start;
    copy.flags &= !READ_UNMAPPED;
    // `SAMRecord.NO_ALIGNMENT_CIGAR`.
    copy.cigar = htsjdk_bam::cigar::Cigar::default();
    copy.mapping_quality = 0;
    copy.flags |= READ_UNMAPPED;
    copy
}

/// Why `undoDistantMateAlterations` could not put a read back.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UndoError {
    /// The `catch` the reference wraps everything in.
    Unrecoverable(String),
    /// What the reference does **not** raise here: `setReferenceName` keeps the name and resolves
    /// it lazily, so this surfaces later, out of whatever first asks for the index.
    UnknownContig(String),
}

impl UndoError {
    /// The reference's own message.
    pub fn message(&self) -> String {
        match self {
            UndoError::Unrecoverable(oa) => format!("can't recover alignment from OA tag: {oa}"),
            UndoError::UnknownContig(name) => {
                format!("Reference index for '{name}' not found in sequence dictionary.")
            }
        }
    }

    pub fn class(&self) -> &'static str {
        match self {
            UndoError::Unrecoverable(_) => "org.broadinstitute.hellbender.exceptions.UserException",
            UndoError::UnknownContig(_) => "java.lang.IllegalArgumentException",
        }
    }
}

/// `undoDistantMateAlterations`: the read as it was before [`do_distant_mate_alterations`].
///
/// A read with no `OA` is returned unchanged, and the reference returns *the same object* rather
/// than a copy, which a caller that then mutates the result can observe.
pub fn undo_distant_mate_alterations(
    read: &BamRecord,
    header: &SamHeader,
) -> Result<BamRecord, UndoError> {
    let Some(TagValue::Str(oa)) = read.tags.get(Tag::new(ORIGINAL_ALIGNMENT_TAG)) else {
        return Ok(read.clone());
    };
    let oa = oa.clone();
    let mut copy = read.clone();
    copy.tags.remove(Tag::new(DISTANT_MATE_TAG));
    copy.tags.remove(Tag::new(ORIGINAL_ALIGNMENT_TAG));

    // `oaValue.split(",")`, then five indexed reads and two parseInts, all inside one try.
    let tokens: Vec<&str> = oa.split(',').collect();
    let unrecoverable = || UndoError::Unrecoverable(oa.clone());
    let contig = *tokens.first().ok_or_else(unrecoverable)?;
    let start: i32 = tokens
        .get(1)
        .ok_or_else(unrecoverable)?
        .parse()
        .map_err(|_| unrecoverable())?;
    let strand = *tokens.get(2).ok_or_else(unrecoverable)?;
    let cigar = *tokens.get(3).ok_or_else(unrecoverable)?;
    let mapping_quality: i32 = tokens
        .get(4)
        .ok_or_else(unrecoverable)?
        .parse()
        .map_err(|_| unrecoverable())?;
    let edit_distance = *tokens.get(5).ok_or_else(unrecoverable)?;

    let cigar = parse_cigar(cigar).map_err(|_| unrecoverable())?;
    let mapping_quality =
        u8::try_from(mapping_quality).map_err(|_| UndoError::Unrecoverable(oa.clone()))?;

    // The name is resolved here rather than kept, which is the one place this port cannot be lazy.
    let index = sequence_index(header, contig)
        .ok_or_else(|| UndoError::UnknownContig(contig.to_string()))?;
    copy.reference_index = index;
    copy.alignment_start = start;
    copy.flags &= !READ_UNMAPPED;
    if strand == "-" {
        copy.flags |= READ_REVERSE_STRAND;
    } else {
        copy.flags &= !READ_REVERSE_STRAND;
    }
    copy.cigar = cigar;
    copy.mapping_quality = mapping_quality;
    // `if (tokens[5].length() > 1)`: the field carries the trailing semicolon, so a length of one
    // is an empty NM and the tag stays absent.
    if edit_distance.len() > 1 {
        let value: i64 = edit_distance[..edit_distance.len() - 1]
            .parse()
            .map_err(|_| unrecoverable())?;
        copy.tags
            .insert(Tag::new(EDIT_DISTANCE_TAG), TagValue::Int(value));
    }
    Ok(copy)
}

/// `PrintDistantMates`: every read the traversal reaches, moved onto its mate, coordinate-sorted.
///
/// The sort is `preSorted = false`, not a tidying step: the transform makes the traversal order
/// wrong, and htsjdk's sorting collection sorts what it is handed with
/// [`htsjdk_bam::coordinate::compare`]. A stable sort, because two records that compare equal keep
/// the order they arrived in.
pub fn print_distant_mates(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ReadsError> {
    print_distant_mates_with(
        source,
        options,
        filter,
        htsjdk_bgzf::DEFAULT_COMPRESSION_LEVEL,
        htsjdk_bgzf::Deflater::Jdk,
    )
}

/// The same, with the BGZF compression named: see [`crate::sam_output::write_records_with`].
///
/// A real `gatk` run writes at level TWO through GKL rather than at htsjdk's default of five, so a
/// runner that wants the reference's BYTES has to say so.
pub fn print_distant_mates_with(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ReadsError> {
    let records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    let input_header = source.header().clone();
    let mut altered: Vec<BamRecord> = records
        .iter()
        .map(|record| do_distant_mate_alterations(record, &input_header))
        .collect();
    altered.sort_by(coordinate::compare);
    let header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    crate::sam_output::write_records_with(
        &header,
        &altered,
        options.create_output_bam_index,
        level,
        deflater,
    )
}

/// The output order of a set of altered records, for a caller that wants the order without the
/// bytes. Exposed because the order is the measurement.
pub fn output_order(records: &mut [BamRecord]) {
    records.sort_by(coordinate::compare);
}

/// `SAMRecordCoordinateComparator.compare` under its own name, so a reader of this module does not
/// have to know which crate the writer's comparator lives in.
pub fn coordinate_compare(a: &BamRecord, b: &BamRecord) -> Ordering {
    coordinate::compare(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::SequenceRecord;

    fn header() -> SamHeader {
        let mut header = SamHeader::default();
        header.sequences.push(SequenceRecord::new("chr1", 3000));
        header.sequences.push(SequenceRecord::new("chr2", 3000));
        header
    }

    fn read(name: &str) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            flags: 0x1,
            reference_index: 0,
            alignment_start: 100,
            mapping_quality: 60,
            cigar: parse_cigar("10M").unwrap(),
            mate_reference_index: 1,
            mate_alignment_start: 600,
            ..BamRecord::default()
        }
    }

    #[test]
    fn a_missing_nm_leaves_the_field_empty() {
        let read = read("r0");
        assert_eq!(
            original_alignment_value(&read, &header()),
            "chr1,100,+,10M,60,;"
        );
    }

    #[test]
    fn an_nm_that_is_there_is_written_out() {
        let mut read = read("r0");
        read.tags.insert(Tag::new(b"NM"), TagValue::Int(3));
        assert_eq!(
            original_alignment_value(&read, &header()),
            "chr1,100,+,10M,60,3;"
        );
    }

    #[test]
    fn the_read_lands_on_its_mate_unmapped_and_uncigared() {
        let read = read("r0");
        let moved = do_distant_mate_alterations(&read, &header());
        assert_eq!(moved.reference_index, 1);
        assert_eq!(moved.alignment_start, 600);
        assert_eq!(moved.mapping_quality, 0);
        assert!(moved.cigar.is_empty());
        assert_eq!(moved.flags & READ_UNMAPPED, READ_UNMAPPED);
        assert!(is_distant_mate(&moved));
        // The mate fields are not touched, so the moved read still points at where it now is.
        assert_eq!(moved.mate_reference_index, 1);
    }

    #[test]
    fn the_round_trip_puts_the_tags_back_where_they_were() {
        let header = header();
        let mut read = read("r0");
        read.tags
            .insert(Tag::new(b"RG"), TagValue::Str("rg1".into()));
        read.tags.insert(Tag::new(b"NM"), TagValue::Int(3));
        let moved = do_distant_mate_alterations(&read, &header);
        let back = undo_distant_mate_alterations(&moved, &header).unwrap();
        assert_eq!(back, read);
    }

    #[test]
    fn a_read_with_no_nm_comes_back_with_no_nm() {
        let header = header();
        let read = read("r0");
        let moved = do_distant_mate_alterations(&read, &header);
        let back = undo_distant_mate_alterations(&moved, &header).unwrap();
        assert_eq!(back, read);
        assert!(back.tags.get(Tag::new(b"NM")).is_none());
    }

    #[test]
    fn a_read_with_no_oa_is_returned_as_it_was() {
        let header = header();
        let read = read("r0");
        assert_eq!(undo_distant_mate_alterations(&read, &header).unwrap(), read);
    }

    #[test]
    fn an_oa_it_cannot_parse_is_the_reference_s_user_exception() {
        let header = header();
        let mut read = read("r0");
        for oa in ["garbage", "chr1,100,+,10M,60", "chr1,x,+,10M,60,3;"] {
            read.tags
                .insert(Tag::new(b"OA"), TagValue::Str(oa.to_string()));
            let error = undo_distant_mate_alterations(&read, &header).unwrap_err();
            assert_eq!(error, UndoError::Unrecoverable(oa.to_string()));
            assert_eq!(
                error.message(),
                format!("can't recover alignment from OA tag: {oa}")
            );
        }
    }

    #[test]
    fn a_contig_the_dictionary_does_not_hold_is_a_different_refusal() {
        let header = header();
        let mut read = read("r0");
        read.tags.insert(
            Tag::new(b"OA"),
            TagValue::Str("chr9,100,+,10M,60,3;".to_string()),
        );
        let error = undo_distant_mate_alterations(&read, &header).unwrap_err();
        assert_eq!(error.class(), "java.lang.IllegalArgumentException");
        assert_eq!(
            error.message(),
            "Reference index for 'chr9' not found in sequence dictionary."
        );
    }

    /// The traversal order and the output order disagree, which is the whole point of
    /// `preSorted = false`.
    #[test]
    fn the_output_is_coordinate_ordered_not_traversal_ordered() {
        let header = header();
        let mut traversal = Vec::new();
        for (name, mate_ref, mate_start) in [("r0", 1, 600), ("r4", 0, 2500), ("r5", 1, 150)] {
            let mut record = read(name);
            record.mate_reference_index = mate_ref;
            record.mate_alignment_start = mate_start;
            traversal.push(do_distant_mate_alterations(&record, &header));
        }
        output_order(&mut traversal);
        let names: Vec<&str> = traversal.iter().map(|r| r.read_name.as_str()).collect();
        assert_eq!(names, ["r4", "r5", "r0"]);
    }
}
