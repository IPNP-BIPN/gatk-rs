//! `SATagBuilder`, ported from `org.broadinstitute.hellbender.utils.SATagBuilder` (GATK 4.6.2.0).
//!
//! The SA tag is how a split read says where its other pieces went. `SplitNCigarReads` is what
//! needs it: `repairSupplementaryTags` marks every piece after the first as supplementary and hands
//! the family to [`set_reads_as_supplemental`].
//!
//! # A unit is six fields, split with a limit of -1
//!
//! ```java
//! String[] values = SATag.split(",", -1);
//! if (values.length != 6) { throw new GATKException("Could not parse SATag: " + SATag); }
//! ```
//!
//! The `-1` is the whole difference between accepting `chr1,100,+,10M,60,;` and refusing it: Java's
//! zero-limit `split` drops trailing empty strings, so the same tag would arrive with five fields
//! and be refused. An NM that is present but empty is a real thing a writer produces, which is why
//! the limit is there.
//!
//! Three fields are validated and each accepts `*`: the position and the mapping quality must not
//! parse as negative, and the cigar must match `\*|([0-9]+[MIDNSHPX=])+`. That pattern accepts
//! `1M1M1M` and refuses the empty string. **NM is not validated at all**, so any text survives.
//!
//! # The order the units come out in is the point
//!
//! `addTag` puts a **non-supplementary** read at the front of the list and a supplementary one at
//! the back, so the primary alignment is the first unit of everyone's tag.
//! [`set_reads_as_supplemental`] marks every read but the first as supplementary **before** any
//! builder is constructed, so the marking is what decides the order rather than the argument list.
//!
//! Existing tags are preserved: the builder parses the read's own SA at construction, so a primary
//! that already claimed a piece keeps that claim ahead of its new siblings. And the primary is
//! never un-marked: a read that arrives already carrying 0x800 keeps it.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use crate::read::flags;

/// The `SA` tag.
pub const SA_TAG: Tag = Tag(((b'A' as i16) << 8) | b'S' as i16);

/// `GATKException` from `SARead`'s parsing constructor.
///
/// The message is the reference's, because a port that refuses the same input with different words
/// is a port a caller cannot switch to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SaTagError {
    /// Not six comma-separated fields.
    Unit(String),
    /// A position that is neither `*` nor a non-negative number.
    Position(String),
    /// A cigar that does not match the pattern.
    Cigar(String),
    /// A mapping quality that is neither `*` nor a non-negative number.
    MappingQuality(String),
}

impl SaTagError {
    /// The exact text `GATKException` carries.
    pub fn message(&self) -> String {
        match self {
            SaTagError::Unit(tag) => format!("Could not parse SATag: {tag}"),
            SaTagError::Position(tag) => format!("Could not parse POS in SATag: {tag}"),
            SaTagError::Cigar(tag) => format!("Could not parse cigar in SATag: {tag}"),
            SaTagError::MappingQuality(tag) => format!("Could not parse MapQ in SATag: {tag}"),
        }
    }
}

/// `SATagBuilder.SARead`: one unit of an SA tag, kept as the six strings it is written as.
///
/// The fields stay strings because that is what the reference stores: a position of `*` and an NM
/// of `not-a-number` both round trip, which they could not if this parsed them into numbers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SaRead {
    /// `null` for a read with no contig, which prints as `*`.
    pub contig: Option<String>,
    pub pos: String,
    pub strand: String,
    pub cigar: String,
    pub mapq: String,
    pub nm: String,
}

impl SaRead {
    /// The unit a read becomes.
    ///
    /// An absent NM is `*`. Every other field comes from the read itself, so an unmapped read gives
    /// `*,0,+,*,0,*;`: the position and the mapping quality are the read's own zeroes, not the
    /// defaults in `toString`, which are unreachable from a read.
    pub fn from_record(record: &BamRecord, header: &SamHeader) -> SaRead {
        SaRead {
            contig: contig_of(record, header).map(str::to_string),
            pos: record.alignment_start.to_string(),
            strand: if record.flags & flags::READ_REVERSE_STRAND != 0 {
                "-".to_string()
            } else {
                "+".to_string()
            },
            cigar: record.cigar.to_text(),
            mapq: record.mapping_quality.to_string(),
            nm: match record.tags.get(Tag::new(b"NM")) {
                Some(value) => attribute_as_string(value),
                None => "*".to_string(),
            },
        }
    }

    /// `new SARead(String)`: the parsing constructor, with its three validations.
    pub fn parse(unit: &str) -> Result<SaRead, SaTagError> {
        let values: Vec<&str> = split_keeping_trailing(unit);
        if values.len() != 6 {
            return Err(SaTagError::Unit(unit.to_string()));
        }
        if values[1] != "*" && parses_negative(values[1]) {
            return Err(SaTagError::Position(unit.to_string()));
        }
        if !matches_cigar_pattern(values[3]) {
            return Err(SaTagError::Cigar(unit.to_string()));
        }
        if values[4] != "*" && parses_negative(values[4]) {
            return Err(SaTagError::MappingQuality(unit.to_string()));
        }
        Ok(SaRead {
            contig: Some(values[0].to_string()),
            pos: values[1].to_string(),
            strand: values[2].to_string(),
            cigar: values[3].to_string(),
            mapq: values[4].to_string(),
            nm: values[5].to_string(),
        })
    }

    /// `SARead.toString`, trailing `;` included.
    ///
    /// The strand is normalised here and nowhere else: anything that is not exactly `-` prints as
    /// `+`, so a unit parsed with a strand of `x` does not survive the round trip.
    pub fn to_text(&self) -> String {
        format!(
            "{},{},{},{},{},{};",
            self.contig.as_deref().unwrap_or("*"),
            self.pos,
            if self.strand == "-" { "-" } else { "+" },
            self.cigar,
            self.mapq,
            self.nm,
        )
    }
}

/// `SATagBuilder`: the units a read's SA tag will be written from.
///
/// The reference holds the read and reaches back into it for `isPrimary` and for its own unit. Here
/// both are captured when the builder is constructed, which is the same thing for the only caller
/// that exists: `setReadsAsSupplemental` sets every supplementary flag before it builds anything.
#[derive(Debug, Clone)]
pub struct SaTagBuilder {
    units: Vec<SaRead>,
    this_read: SaRead,
    primary: bool,
}

impl SaTagBuilder {
    /// The builder for a read, with its existing SA tag already parsed.
    pub fn new(record: &BamRecord, header: &SamHeader) -> Result<SaTagBuilder, SaTagError> {
        let units = match record.tags.get(SA_TAG) {
            Some(value) => parse_tag(&attribute_as_string(value))?,
            None => Vec::new(),
        };
        Ok(SaTagBuilder {
            units,
            this_read: SaRead::from_record(record, header),
            primary: record.flags & flags::SUPPLEMENTARY_ALIGNMENT == 0,
        })
    }

    /// `clear()`.
    pub fn clear(&mut self) -> &mut Self {
        self.units.clear();
        self
    }

    /// `addTag(SATagBuilder)`: a primary other read goes to the front, a supplementary one to the
    /// back.
    pub fn add_tag(&mut self, other: &SaTagBuilder) -> &mut Self {
        if other.primary {
            self.units.insert(0, other.this_read.clone());
        } else {
            self.units.push(other.this_read.clone());
        }
        self
    }

    /// `addTag(GATKRead)`, which asks the read's own flag rather than a builder's.
    pub fn add_tag_record(&mut self, record: &BamRecord, header: &SamHeader) -> &mut Self {
        let unit = SaRead::from_record(record, header);
        if record.flags & flags::SUPPLEMENTARY_ALIGNMENT == 0 {
            self.units.insert(0, unit);
        } else {
            self.units.push(unit);
        }
        self
    }

    /// `removeTag(contig, start)`: every unit whose contig and position match, by string.
    pub fn remove_tag(&mut self, contig: Option<&str>, start: i32) -> &mut Self {
        let start = start.to_string();
        self.units
            .retain(|unit| !(unit.contig.as_deref() == contig && unit.pos == start));
        self
    }

    /// The units, in the order they will be written.
    pub fn units(&self) -> &[SaRead] {
        &self.units
    }

    /// `getTag()`: the units concatenated, each already carrying its `;`.
    pub fn tag(&self) -> String {
        self.units.iter().map(SaRead::to_text).collect()
    }

    /// `setSATag()`: write the tag onto the read, **unless there is nothing to write**.
    ///
    /// A builder with no units leaves the read alone, so a read with no SA tag and no additions
    /// does not gain an empty one.
    pub fn set_sa_tag(&self, record: &mut BamRecord) {
        if self.units.is_empty() {
            return;
        }
        record.tags.insert(SA_TAG, TagValue::Str(self.tag()));
    }
}

/// `parseSATag`: a whole tag split on `;`.
pub fn parse_tag(tag: &str) -> Result<Vec<SaRead>, SaTagError> {
    // Java's `split(";")` drops trailing empty strings, so the `;` every unit ends with does not
    // produce a seventh, empty unit.
    let mut units = Vec::new();
    for piece in split_dropping_trailing(tag, ';') {
        units.push(SaRead::parse(piece)?);
    }
    Ok(units)
}

/// `SATagBuilder.setReadsAsSupplemental(primary, supplementalReads)`.
///
/// Every read but the primary is marked supplementary first, then every read gets a unit for every
/// other read and never for itself. The primary is not un-marked: one that arrives supplementary
/// stays that way, and then no read in the family is the primary line.
pub fn set_reads_as_supplemental(
    primary: &mut BamRecord,
    supplemental: &mut [BamRecord],
    header: &SamHeader,
) -> Result<(), SaTagError> {
    let mut builders = Vec::with_capacity(supplemental.len() + 1);
    builders.push(SaTagBuilder::new(primary, header)?);

    for record in supplemental.iter_mut() {
        record.flags |= flags::SUPPLEMENTARY_ALIGNMENT;
        builders.push(SaTagBuilder::new(record, header)?);
    }

    let snapshot = builders.clone();
    for (i, builder) in builders.iter_mut().enumerate() {
        for (j, other) in snapshot.iter().enumerate() {
            if i != j {
                builder.add_tag(other);
            }
        }
    }

    let mut iter = builders.iter();
    if let Some(builder) = iter.next() {
        builder.set_sa_tag(primary);
    }
    for (builder, record) in iter.zip(supplemental.iter_mut()) {
        builder.set_sa_tag(record);
    }
    Ok(())
}

/// `GATKRead.getAttributeAsString`, for the two value kinds an SA or NM tag can hold.
fn attribute_as_string(value: &TagValue) -> String {
    match value {
        TagValue::Str(text) => text.clone(),
        TagValue::Int(number) => number.to_string(),
        TagValue::Char(byte) => (*byte as char).to_string(),
        TagValue::Float(number) => number.to_string(),
        // A byte array is decoded as text rather than printed as an object, which is the one case
        // the adapter special-cases.
        TagValue::ByteArray { values, .. } => {
            String::from_utf8_lossy(&values.iter().map(|&v| v as u8).collect::<Vec<_>>())
                .into_owned()
        }
        other => format!("{other:?}"),
    }
}

/// The contig `read.getContig()` resolves to, or none.
fn contig_of<'a>(record: &BamRecord, header: &'a SamHeader) -> Option<&'a str> {
    usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map(|sequence| sequence.name.as_str())
}

/// `String.split(",", -1)`: every field, trailing empties included.
fn split_keeping_trailing(text: &str) -> Vec<&str> {
    text.split(',').collect()
}

/// `String.split(separator)` with no limit: trailing empty fields are dropped.
fn split_dropping_trailing(text: &str, separator: char) -> Vec<&str> {
    let mut pieces: Vec<&str> = text.split(separator).collect();
    while pieces.last() == Some(&"") {
        pieces.pop();
    }
    pieces
}

/// `Integer.parseInt(value) < 0`, where text that is not a number is not negative.
///
/// The reference lets a `NumberFormatException` out of here rather than catching it, so a position
/// of `abc` is a different failure from a position of `-1`. Nothing in GATK reaches that path with
/// text, and this reproduces the comparison only.
fn parses_negative(value: &str) -> bool {
    value.starts_with('-') && value[1..].chars().all(|c| c.is_ascii_digit())
}

/// `\*|([0-9]+[MIDNSHPX=])+`, matched against the whole string.
///
/// It counts nothing: `1M1M1M` matches, and so does a cigar whose lengths do not add up to the
/// read. What it refuses is an empty string and any operator outside the set.
fn matches_cigar_pattern(cigar: &str) -> bool {
    if cigar == "*" {
        return true;
    }
    if cigar.is_empty() {
        return false;
    }
    let bytes = cigar.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let digits = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if i == digits || i == bytes.len() {
            return false;
        }
        if !matches!(
            bytes[i],
            b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'='
        ) {
            return false;
        }
        i += 1;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_split_keeps_a_trailing_empty_field() {
        // Six fields with an empty NM parses; five does not.
        assert!(SaRead::parse("chr1,100,+,10M,60,").is_ok());
        assert_eq!(
            SaRead::parse("chr1,100,+,10M,60").unwrap_err(),
            SaTagError::Unit("chr1,100,+,10M,60".to_string())
        );
    }

    #[test]
    fn the_strand_is_the_one_field_that_does_not_round_trip() {
        let unit = SaRead::parse("chr1,100,x,10M,60,2").expect("parses");
        assert_eq!(unit.to_text(), "chr1,100,+,10M,60,2;");
    }

    #[test]
    fn nm_is_never_validated() {
        let unit = SaRead::parse("chr1,100,+,10M,60,not-a-number").expect("parses");
        assert_eq!(unit.to_text(), "chr1,100,+,10M,60,not-a-number;");
    }

    #[test]
    fn the_cigar_pattern_counts_nothing_and_refuses_the_empty_string() {
        assert!(matches_cigar_pattern("*"));
        assert!(matches_cigar_pattern("1M1M1M"));
        assert!(matches_cigar_pattern("10M2I5S"));
        assert!(!matches_cigar_pattern(""));
        assert!(!matches_cigar_pattern("10Z"));
        assert!(!matches_cigar_pattern("M10"));
        assert!(!matches_cigar_pattern("10"));
    }

    #[test]
    fn a_negative_position_and_a_negative_mapping_quality_are_refused() {
        assert_eq!(
            SaRead::parse("chr1,-1,+,10M,60,2").unwrap_err(),
            SaTagError::Position("chr1,-1,+,10M,60,2".to_string())
        );
        assert_eq!(
            SaRead::parse("chr1,100,+,10M,-1,2").unwrap_err(),
            SaTagError::MappingQuality("chr1,100,+,10M,-1,2".to_string())
        );
        // Both accept `*`.
        assert!(SaRead::parse("chr1,*,+,10M,*,2").is_ok());
    }
}
