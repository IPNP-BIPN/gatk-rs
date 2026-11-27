//! `ReadUtils.getSAMReadGroupRecord`: a read's `RG` resolved against the header's `@RG` lines.
//!
//! Ported from `org.broadinstitute.hellbender.utils.read.ReadUtils` (GATK 4.6.2.0).
//!
//! The resolution is where "the record mentions a group" and "the header knows that group" stop
//! being the same statement. A read whose `RG` names a group the header does not declare resolves
//! to nothing here, while `HasReadGroupReadFilter`, which tests the raw attribute, keeps it.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<&'a ReadGroup>` is "an `@RG` line or nothing", which is the reference's nullable
//!   `SAMReadGroupRecord`. `None` is Java's `null`.
//! - The `'a` ties the answer to the header: what comes back is a view into `header`, not a copy,
//!   so it lives exactly as long as the header does. Java needs no such mark because its garbage
//!   collector keeps the header alive for it.
//! - `?` after an `Option` returns `None` from the whole function at once when there is nothing
//!   there, which is the reference's `readGroupName != null ? ... : null`.

use htsjdk_bam::header::{ReadGroup, SamHeader};
use htsjdk_bam::record::BamRecord;

/// The `@RG` line the read's `RG` attribute names, if the header declares it.
///
/// The reference, `ReadUtils.java:230-231`:
///
/// ```java
/// final String readGroupName = read.getReadGroup();
/// return readGroupName != null ? header.getReadGroup(readGroupName) : null;
/// ```
///
/// `getReadGroup` is a cast, `return (String)samRecord.getAttribute(SAMTag.RG.getBinaryTag());`
/// (`SAMRecordToGATKReadAdapter.java:384`), so an `RG` stored as anything but a string throws a
/// `ClassCastException` there. Here it is read as the empty id instead, and resolves only if the
/// header declares a group with that id.
pub fn resolve<'a>(record: &BamRecord, header: &'a SamHeader) -> Option<&'a ReadGroup> {
    // `find_map` walks the tags and stops at the first one the closure answers `Some` for, which
    // is the first `RG`. `cond.then_some(x)` is `cond ? x : null`.
    let id = record.tags.iter().find_map(|(tag, value)| {
        (tag.name() == *b"RG").then_some(match value {
            htsjdk_bam::tag::TagValue::Str(text) => text.as_str(),
            _ => "",
        })
    })?;
    // `header.getReadGroup(readGroupName)`: the `@RG` line with that ID, or nothing.
    header.read_groups.iter().find(|group| group.id == id)
}

/// One field of the resolved `@RG` line, such as `SM` or `FO`; nothing when the read has no
/// group, the header does not declare it, or the line lacks the field.
pub fn attribute<'a>(record: &BamRecord, header: &'a SamHeader, key: &str) -> Option<&'a str> {
    resolve(record, header)?.attributes.get(key)
}

/// `SAMReadGroupRecord.getFlowOrder`, which is what marks a read group as flow-based.
pub fn flow_order<'a>(record: &BamRecord, header: &'a SamHeader) -> Option<&'a str> {
    attribute(record, header, "FO")
}
