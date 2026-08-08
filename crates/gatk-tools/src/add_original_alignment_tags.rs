//! Ported from `org.broadinstitute.hellbender.tools.AddOriginalAlignmentTags` (GATK 4.6.2.0).
//!
//! The fifth whole tool of the record-transform archetype, and the first that writes tags rather
//! than changing the read.
//!
//! ```java
//! public void apply(GATKRead read, ReferenceContext r, FeatureContext f) {
//!     addOATag(read);
//!     addMateContigTag(read);
//!     outputWriter.addRead(read);
//! }
//! ```
//!
//! # It aborts on an unpaired read
//!
//! `addMateContigTag` calls `mateIsUnmapped`, which htsjdk defines only for a **paired** read. On a
//! read that is not paired it throws, and the whole run stops. An unpaired read is ordinary in an
//! ordinary file, so this is the tool's behaviour on ordinary input rather than an edge case:
//! measured, three of five runs over a fixture holding one unpaired read among six abort, and the
//! one that succeeds does so because its interval excludes it.
//!
//! # The comma escaping cannot fire
//!
//! The OA format replaces a comma in the contig name with an underscore. A comma is not a legal
//! sequence name character at all, so no valid file reaches that branch: htsjdk refuses `chr,1`
//! against `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`. The replacement is
//! reproduced anyway, because reproducing it costs one call and knowing it is dead costs a
//! measurement.
//!
//! # A missing NM prints as the word `null`
//!
//! `getAttributeAsString` returns null and the format string prints it, so a read with no `NM`
//! gets `chr1,120,+,10M,60,null;`. Four characters, in the file, forever.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use gatk_engine::reads::{ReadsDataSource, ReadsError};

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK AddOriginalAlignmentTags";

/// `MATE_CONTIG_TAG_NAME`.
pub const MATE_CONTIG_TAG: &[u8; 2] = b"XM";
/// `OA_TAG_NAME`.
pub const ORIGINAL_ALIGNMENT_TAG: &[u8; 2] = b"OA";
/// `OA_SEPARATOR`, which is also the character the contig name may not contain.
pub const OA_SEPARATOR: char = ',';

/// `SAMFlag.READ_PAIRED`.
const READ_PAIRED: u16 = 0x1;
/// `SAMFlag.READ_UNMAPPED`.
const READ_UNMAPPED: u16 = 0x4;
/// `SAMFlag.MATE_UNMAPPED`.
const MATE_UNMAPPED: u16 = 0x8;
/// `SAMFlag.READ_REVERSE_STRAND`.
const READ_REVERSE_STRAND: u16 = 0x10;

/// Why the tool stopped: the one thing it refuses, and it refuses it on ordinary input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnpairedRead;

impl UnpairedRead {
    /// htsjdk's own message, out of `SAMRecord.requireReadPaired`.
    pub fn message(&self) -> String {
        "Cannot get mate information for an unpaired read".to_string()
    }

    pub fn class(&self) -> &'static str {
        "java.lang.IllegalStateException"
    }
}

/// `addOATag`: the read's current alignment, formatted, before anything changes it.
///
/// Six fields and a trailing semicolon. The `NM` field is whatever `getAttributeAsString` returned,
/// which for an absent tag is the string `null`.
pub fn original_alignment_value(read: &BamRecord, header: &SamHeader) -> String {
    if read.flags & READ_UNMAPPED != 0 {
        return "*,0,*,*,0,0;".to_string();
    }
    let contig = sequence_name(header, read.reference_index).replace(OA_SEPARATOR, "_");
    let strand = if read.flags & READ_REVERSE_STRAND != 0 {
        "-"
    } else {
        "+"
    };
    let edit_distance = match read.tags.get(Tag::new(b"NM")) {
        Some(TagValue::Int(value)) => value.to_string(),
        Some(other) => format!("{other:?}"),
        // The reference formats a null into the string rather than leaving the field out.
        None => "null".to_string(),
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

/// `addMateContigTag`: the mate's contig, or `*` when the mate is unmapped.
///
/// Refuses an unpaired read, because `mateIsUnmapped` is only defined for a paired one.
pub fn mate_contig_value(read: &BamRecord, header: &SamHeader) -> Result<String, UnpairedRead> {
    if read.flags & READ_PAIRED == 0 {
        return Err(UnpairedRead);
    }
    if read.flags & MATE_UNMAPPED != 0 {
        return Ok("*".to_string());
    }
    Ok(sequence_name(header, read.mate_reference_index))
}

/// Both tags, in the reference's order: the OA is written before the mate contig is asked for, so
/// a read that aborts the run has already had its OA set.
pub fn add_tags(read: &mut BamRecord, header: &SamHeader) -> Result<(), UnpairedRead> {
    let original = original_alignment_value(read, header);
    read.tags
        .insert(Tag::new(ORIGINAL_ALIGNMENT_TAG), TagValue::Str(original));
    let mate = mate_contig_value(read, header)?;
    read.tags
        .insert(Tag::new(MATE_CONTIG_TAG), TagValue::Str(mate));
    Ok(())
}

/// What a run produces: the output BAM and its index, or the refusal, or a failure to read.
pub type RunResult = Result<Result<(Vec<u8>, Option<Vec<u8>>), UnpairedRead>, ReadsError>;

/// `AddOriginalAlignmentTags`: every read the traversal reaches, with two tags added.
pub fn add_original_alignment_tags(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> RunResult {
    let mut records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    let input_header = source.header().clone();
    for record in &mut records {
        if let Err(error) = add_tags(record, &input_header) {
            return Ok(Err(error));
        }
    }
    let header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    Ok(Ok(write_records(
        &header,
        &records,
        options.create_output_bam_index,
    )?))
}
