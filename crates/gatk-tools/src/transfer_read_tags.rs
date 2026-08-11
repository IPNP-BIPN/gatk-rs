//! Ported from `org.broadinstitute.hellbender.tools.walkers.qc.TransferReadTags` (GATK 4.6.2.0).
//!
//! The tenth whole tool of the record-transform archetype, and the first that is not a walker.
//!
//! # The traversal is the tool's own
//!
//! It extends `GATKTool` and overrides `traverse()` rather than inheriting a walker's, so **none of
//! the archetype's usual machinery runs**: no `CountingReadFilter`, no read transformer, no interval
//! bound. `directlyAccessEngineReadsDataSource().iterator()` is the whole traversal on the aligned
//! side, and a second `ReadsPathDataSource` is opened by hand for the unmapped side.
//!
//! A port that reached for [`crate::read_walker`] would apply `WellformedReadFilter` and drop reads
//! this tool writes. Measured: the fixture's reads carry no read group, which the wellformed filter
//! rejects, and every one of them comes out.
//!
//! # Every tag is transferred as a string
//!
//! ```java
//! final String tagValue = originRead.getAttributeAsString(tagName);
//! Utils.nonNull(tagValue, "The attribute is empty: read " + currentUnmappedRead.getName());
//! updatedRead.setAttribute(tagName, tagValue);
//! ```
//!
//! `getAttributeAsString` renders whatever was there and `setAttribute(name, String)` writes a `Z`,
//! so `XI:i:42` in the unmapped file arrives as `XI:Z:42` and `XN:f:1.5` as `XN:Z:1.5`. The
//! characters are the same and the type is not, which is why the golden prints the Java class of
//! every tag value beside it: nothing else in the row would show it.
//!
//! [`attribute_as_string`] is local rather than shared with
//! `gatk_readfilter::jexl_filter::attribute_as_string`, which claims the same reference method: that
//! one renders an array tag as the Rust enum's debug form, and the adapter decodes a `byte[]` as
//! UTF-8. The two agree on every scalar type. The divergence is on the array branch, which no
//! golden here or there covers, and it is recorded rather than propagated.
//!
//! # An aligned read past the end of the unmapped file is silently dropped
//!
//! ```java
//! } else if (diff > 0){
//!     while (unmappedSamIterator.hasNext()){
//!         ...
//!     }
//! }
//! ```
//!
//! The catch-up loop is bounded by `hasNext()`, so when the unmapped side runs out the loop simply
//! ends: nothing is written, nothing is logged, and the outer loop moves to the next aligned read.
//! Measured: an aligned file of `a1, a9` against an unmapped file ending at `a6` produces a one-read
//! output and no error at all.
//!
//! This is the one behaviour here a reasonable port would get wrong by being careful, because the
//! obvious reading of "the aligned file must be a subset" is that a missing read is refused. It is
//! refused in one direction and dropped in the other.
//!
//! # The writer is not told the reads are sorted
//!
//! `createSAMWriter(..., false)` on a queryname header, so htsjdk sorts what it is handed with
//! [`htsjdk_bam::query_name::compare`], which has six tie-breaks after the name. Measured: a name
//! whose second-of-pair record is written before its first-of-pair record comes out the other way
//! round. The comparison the *traversal* makes is
//! [`htsjdk_bam::query_name::compare_read_names`], the name alone with no tie-break, which is a
//! different function of the same class.
//!
//! No output carries an index, whatever `--create-output-bam-index` says: a queryname-sorted file
//! has nothing to index.
//!
//! # The refusals come from four layers
//!
//! `--read-tags` omitted is refused by Barclay before the tool is built at all, so the
//! `Utils.nonEmpty(readTags, ...)` written for it in `onTraversalStart` is unreachable.
//! [`TransferError::NoReadTags`] is that unreachable check, kept because the reference keeps it.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::query_name;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use gatk_engine::reads::{ReadsDataSource, ReadsError};

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK TransferReadTags";

/// What the tool refuses on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TransferError {
    /// `Utils.nonEmpty(readTags, "read tags may not be empty")`, which Barclay never lets the tool
    /// reach: a `List` argument with no `optional = true` is required, so the observed refusal is
    /// `CommandLineException$MissingArgument`.
    NoReadTags,
    /// `Utils.validate(sortOrderAlignedReads == queryname, "aligned sam must be sorted by
    /// queryname")`. The unmapped side is not checked, and the reference says why: the field "is
    /// often not populated".
    AlignedNotQueryNameSorted,
    /// `UserException("Unmapped sam iterator is empty and aligned sam iterator is not.")`.
    UnmappedEmptyAndAlignedIsNot,
    /// `IllegalStateException`: an aligned read that is behind the unmapped cursor. Raised from
    /// either of two sites with the same message, and never for a read that is *ahead* of the last
    /// unmapped read, which is dropped instead.
    NotInUnmapped { aligned: String, unmapped: String },
    /// `Utils.nonNull(tagValue, "The attribute is empty: read " + currentUnmappedRead.getName())`,
    /// which names the **unmapped** read.
    AttributeEmpty { unmapped: String },
}

impl TransferError {
    /// The message the reference raises, for the two that carry one worth comparing.
    pub fn message(&self) -> String {
        match self {
            TransferError::NoReadTags => "read tags may not be empty".to_string(),
            TransferError::AlignedNotQueryNameSorted => {
                "aligned sam must be sorted by queryname".to_string()
            }
            TransferError::UnmappedEmptyAndAlignedIsNot => {
                "Unmapped sam iterator is empty and aligned sam iterator is not.".to_string()
            }
            TransferError::NotInUnmapped { aligned, unmapped } => format!(
                "A read found in the aligned bam is not found in the unmapped bam. This tool \
                 assumes reads in both input files are query-name sorted lexicographically (i.e. \
                 by Picard SortSam but not by samtools sort): aligned read = {aligned}, unmapped \
                 read = {unmapped}"
            ),
            TransferError::AttributeEmpty { unmapped } => {
                format!("The attribute is empty: read {unmapped}")
            }
        }
    }
}

/// `SAMFileHeader.getSortOrder() == queryname`.
fn is_query_name_sorted(header: &SamHeader) -> bool {
    header.attributes.get("SO") == Some("queryname")
}

/// `GATKRead.getAttributeAsString`: the value's `toString`, except a `byte[]`, which is decoded.
///
/// The reference's own branch order, including the empty-array case that returns `""` rather than
/// decoding nothing.
pub fn attribute_as_string(read: &BamRecord, name: &str) -> Option<String> {
    let bytes = name.as_bytes();
    if bytes.len() != 2 {
        return None;
    }
    match read.tags.get(Tag::new(&[bytes[0], bytes[1]]))? {
        TagValue::Str(text) => Some(text.clone()),
        TagValue::Char(c) => Some((*c as char).to_string()),
        // Every integral Java box type reaches `toString` as a decimal, and htsjdk narrows them all
        // to one variant.
        TagValue::Int(value) => Some(value.to_string()),
        TagValue::Float(value) => Some(format_java_float(*value)),
        // `instanceof byte[]`: decoded as UTF-8, and `""` when empty.
        TagValue::ByteArray { values, .. } => Some(if values.is_empty() {
            String::new()
        } else {
            String::from_utf8_lossy(&values.iter().map(|v| *v as u8).collect::<Vec<u8>>())
                .into_owned()
        }),
        // A `short[]`, `int[]` or `float[]` is not a `byte[]`, so the adapter falls through to
        // `Object.toString`, which is `[S@1b6d3586`: an identity hash, not reproducible and not
        // meaningful. The reference produces it; this port refuses to invent one.
        _ => None,
    }
}

/// `Float.toString`, which always prints a decimal point.
fn format_java_float(value: f32) -> String {
    if value == value.trunc() && value.is_finite() && value.abs() < 1e7 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

/// `updateReadTags`: the aligned read, with each named tag taken from the unmapped one as a string.
pub fn update_read_tags(
    target: &BamRecord,
    origin: &BamRecord,
    read_tags: &[String],
) -> Result<BamRecord, TransferError> {
    let mut updated = target.clone();
    for name in read_tags {
        let value =
            attribute_as_string(origin, name).ok_or_else(|| TransferError::AttributeEmpty {
                unmapped: origin.read_name.clone(),
            })?;
        let bytes = name.as_bytes();
        updated
            .tags
            .insert(Tag::new(&[bytes[0], bytes[1]]), TagValue::Str(value));
    }
    Ok(updated)
}

/// What a run produces: the output BAM and its index, which is always absent here.
pub type RunResult = Result<Result<(Vec<u8>, Option<Vec<u8>>), TransferError>, ReadsError>;

/// `TransferReadTags`: the aligned reads, each carrying the named tags of its unmapped counterpart.
pub fn transfer_read_tags(
    aligned: &ReadsDataSource,
    unmapped: &ReadsDataSource,
    read_tags: &[String],
    options: &Options,
) -> RunResult {
    if read_tags.is_empty() {
        return Ok(Err(TransferError::NoReadTags));
    }
    if !is_query_name_sorted(aligned.header()) {
        return Ok(Err(TransferError::AlignedNotQueryNameSorted));
    }

    let aligned_reads = aligned.iter_all()?;
    let unmapped_reads = unmapped.iter_all()?;

    // `onTraversalStart` pulls the first unmapped read before the traversal begins, and refuses
    // only when there is an aligned read to match it against.
    let mut unmapped_cursor = 0usize;
    if unmapped_reads.is_empty() {
        if !aligned_reads.is_empty() {
            return Ok(Err(TransferError::UnmappedEmptyAndAlignedIsNot));
        }
        // "Input data contains no reads. Output will also contain no reads." is a warning, and the
        // empty output is written anyway.
        let header = header_for_sam_writer(aligned.header(), TOOL_NAME, options);
        return Ok(Ok(write_records(&header, &[], false)?));
    }

    let mut records = Vec::with_capacity(aligned_reads.len());
    for target in &aligned_reads {
        let current = &unmapped_reads[unmapped_cursor];
        let diff = query_name::compare_read_names(&target.read_name, &current.read_name);
        match diff {
            std::cmp::Ordering::Equal => match update_read_tags(target, current, read_tags) {
                Ok(read) => records.push(read),
                Err(error) => return Ok(Err(error)),
            },
            std::cmp::Ordering::Greater => {
                // Play the unmapped reads forward until they catch up. If they run out first, this
                // aligned read is written nowhere and nothing says so.
                let mut matched = None;
                while unmapped_cursor + 1 < unmapped_reads.len() {
                    unmapped_cursor += 1;
                    let current = &unmapped_reads[unmapped_cursor];
                    match query_name::compare_read_names(&target.read_name, &current.read_name) {
                        std::cmp::Ordering::Greater => continue,
                        std::cmp::Ordering::Equal => {
                            matched = Some(unmapped_cursor);
                            break;
                        }
                        std::cmp::Ordering::Less => {
                            return Ok(Err(TransferError::NotInUnmapped {
                                aligned: target.read_name.clone(),
                                unmapped: current.read_name.clone(),
                            }))
                        }
                    }
                }
                if let Some(index) = matched {
                    match update_read_tags(target, &unmapped_reads[index], read_tags) {
                        Ok(read) => records.push(read),
                        Err(error) => return Ok(Err(error)),
                    }
                }
            }
            std::cmp::Ordering::Less => {
                return Ok(Err(TransferError::NotInUnmapped {
                    aligned: target.read_name.clone(),
                    unmapped: current.read_name.clone(),
                }))
            }
        }
    }

    // `preSorted = false` on a queryname header: htsjdk's sorting collection, whose comparator is
    // the full queryname one with its six tie-breaks, not the name-only comparison the traversal
    // made. Stable, because two records that compare equal keep the order they were written in.
    records.sort_by(query_name::compare);

    let header = header_for_sam_writer(aligned.header(), TOOL_NAME, options);
    // No index: the output is queryname sorted and there is nothing to index, whatever
    // `--create-output-bam-index` says.
    Ok(Ok(write_records(&header, &records, false)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(name: &str, flags: u16) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            flags,
            ..BamRecord::default()
        }
    }

    fn tagged(name: &str, tag: &[u8; 2], value: TagValue) -> BamRecord {
        let mut record = read(name, 4);
        record.tags.insert(Tag::new(tag), value);
        record
    }

    #[test]
    fn an_integer_tag_arrives_as_a_string() {
        let origin = tagged("a1", b"XI", TagValue::Int(42));
        let updated =
            update_read_tags(&read("a1", 0), &origin, &["XI".to_string()]).expect("transferred");
        assert_eq!(
            updated.tags.get(Tag::new(b"XI")),
            Some(&TagValue::Str("42".to_string())),
            "the characters are the same and the type is not"
        );
    }

    #[test]
    fn and_a_float_through_javas_tostring() {
        let origin = tagged("a1", b"XN", TagValue::Float(1.5));
        let updated =
            update_read_tags(&read("a1", 0), &origin, &["XN".to_string()]).expect("transferred");
        assert_eq!(
            updated.tags.get(Tag::new(b"XN")),
            Some(&TagValue::Str("1.5".to_string()))
        );
        // `Float.toString` always prints a decimal point, which Rust's `{}` does not.
        assert_eq!(format_java_float(2.0), "2.0");
    }

    #[test]
    fn a_missing_tag_names_the_unmapped_read() {
        let error = update_read_tags(&read("a1", 0), &read("a3", 4), &["RX".to_string()])
            .expect_err("refused");
        assert_eq!(
            error,
            TransferError::AttributeEmpty {
                unmapped: "a3".to_string()
            }
        );
        assert_eq!(error.message(), "The attribute is empty: read a3");
    }

    #[test]
    fn a_byte_array_tag_is_decoded_rather_than_printed() {
        let origin = tagged(
            "a1",
            b"XB",
            TagValue::ByteArray {
                values: vec![72, 105],
                unsigned: false,
            },
        );
        assert_eq!(attribute_as_string(&origin, "XB").as_deref(), Some("Hi"));
    }

    #[test]
    fn the_refusal_messages_are_the_references() {
        assert_eq!(
            TransferError::AlignedNotQueryNameSorted.message(),
            "aligned sam must be sorted by queryname"
        );
        assert_eq!(
            TransferError::UnmappedEmptyAndAlignedIsNot.message(),
            "Unmapped sam iterator is empty and aligned sam iterator is not."
        );
        let not_found = TransferError::NotInUnmapped {
            aligned: "a1".to_string(),
            unmapped: "a4".to_string(),
        };
        assert!(not_found
            .message()
            .ends_with("aligned read = a1, unmapped read = a4"));
    }
}
