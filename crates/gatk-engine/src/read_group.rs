//! `ReadUtils.getSAMReadGroupRecord`: a read's `RG` resolved against the header's `@RG` lines.
//!
//! Ported from `org.broadinstitute.hellbender.utils.read.ReadUtils` (GATK 4.6.2.0).
//!
//! The resolution is where "the record mentions a group" and "the header knows that group" stop
//! being the same statement. A read whose `RG` names a group the header does not declare resolves
//! to nothing here, while `HasReadGroupReadFilter`, which tests the raw attribute, keeps it.

use htsjdk_bam::header::{ReadGroup, SamHeader};
use htsjdk_bam::record::BamRecord;

/// The `@RG` line the read's `RG` attribute names, if the header declares it.
pub fn resolve<'a>(record: &BamRecord, header: &'a SamHeader) -> Option<&'a ReadGroup> {
    let id = record.tags.iter().find_map(|(tag, value)| {
        (tag.name() == *b"RG").then_some(match value {
            htsjdk_bam::tag::TagValue::Str(text) => text.as_str(),
            _ => "",
        })
    })?;
    header.read_groups.iter().find(|group| group.id == id)
}

/// One field of the resolved `@RG` line.
pub fn attribute<'a>(record: &BamRecord, header: &'a SamHeader, key: &str) -> Option<&'a str> {
    resolve(record, header)?.attributes.get(key)
}

/// `SAMReadGroupRecord.getFlowOrder`, which is what marks a read group as flow-based.
pub fn flow_order<'a>(record: &BamRecord, header: &'a SamHeader) -> Option<&'a str> {
    attribute(record, header, "FO")
}
