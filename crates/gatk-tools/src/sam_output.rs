//! What every read-writing tool pays before it does anything of its own.
//!
//! Ported from `org.broadinstitute.hellbender.engine.GATKTool`: `getHeaderForSAMWriter`,
//! `createProgramGroupID` and `createSAMWriter` (GATK 4.6.2.0).
//!
//! This module exists because of a measurement rather than a taste for structure. `PrintReads` was
//! the first whole tool in this crate, and the question G2's calibration gate asks is what the
//! *second* one costs. Splitting the shared part out is the only way to answer it honestly: what
//! stays here is what the first tool paid for once, and what a new tool file contains is its
//! marginal cost.
//!
//! # The `@PG` record is where the bytes are decided
//!
//! Three things about it are observable and easy to get wrong:
//!
//!  * **the tool name and the command line are different strings.** `getToolName()` returns
//!    `GATK PrintReads`, with a space; the command line the tool records begins `PrintReads`. Both
//!    land in the same `@PG` record, in different fields;
//!  * **the ID is made unique by suffixing.** `.1`, `.2`, and so on until one is free, so running a
//!    tool over its own output appends a second record rather than replacing or duplicating the
//!    first;
//!  * **the attribute order is the order the setters run**, not alphabetical and not the order the
//!    SAM specification lists: `VN`, then `CL`, then `PN`, after the `ID` every `@PG` line leads
//!    with.
//!
//! `getHeaderForSAMWriter` mutates the reads header in place rather than copying it. That is not
//! observable in a single run's output, and is reproduced anyway.

use gatk_engine::reads::ReadsError;
use htsjdk_bam::header::{ProgramRecord, SamHeader};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::writer::BamWriter;

use gatk_engine::interval::SimpleInterval;

/// The arguments a read-writing tool has in common. Only those that change the output bytes.
///
/// Every tool in this archetype declares the same 42 common arguments through Barclay, and all but
/// these are either inputs or have no effect on what is written.
pub struct Options<'a> {
    /// `-L`. Empty is an unbounded traversal.
    pub intervals: Vec<SimpleInterval>,
    /// `--create-output-bam-index`, whose default is true.
    pub create_output_bam_index: bool,
    /// `--add-output-sam-program-record`, whose default is true.
    pub add_output_sam_program_record: bool,
    /// The expanded command line the tool records in `CL`. An input, not something a port can
    /// invent: Barclay builds it from every argument including the defaults.
    pub command_line: &'a str,
    /// `--version`, which lands in `VN`.
    pub version: &'a str,
}

impl Default for Options<'_> {
    fn default() -> Self {
        Options {
            intervals: Vec::new(),
            create_output_bam_index: true,
            add_output_sam_program_record: true,
            command_line: "",
            version: "4.6.2.0",
        }
    }
}

/// `GATKTool.createProgramGroupID`: the tool name, then the tool name with `.1`, `.2`, until one
/// is free.
pub fn create_program_group_id(header: &SamHeader, tool_name: &str) -> String {
    let taken = |id: &str| header.programs.iter().any(|record| record.id == id);
    if !taken(tool_name) {
        return tool_name.to_string();
    }
    let mut count = 1;
    loop {
        let candidate = format!("{tool_name}.{count}");
        if !taken(&candidate) {
            return candidate;
        }
        count += 1;
    }
}

/// `GATKTool.getHeaderForSAMWriter`: the reads header with the tool's `@PG` appended.
pub fn header_for_sam_writer(header: &SamHeader, tool_name: &str, options: &Options) -> SamHeader {
    let mut header = header.clone();
    if !options.add_output_sam_program_record {
        return header;
    }
    let mut record = ProgramRecord::new(&create_program_group_id(&header, tool_name));
    record.attributes.set("VN", options.version);
    record.attributes.set("CL", options.command_line);
    record.attributes.set("PN", tool_name);
    header.programs.push(record);
    header
}

/// `createSAMWriter(OUTPUT, true)` and everything written through it, as bytes.
///
/// The `true` is `preSorted`, which every tool in this archetype passes: the traversal emits reads
/// in the order the source held them, so the writer is not asked to sort and the index it builds
/// is the index of that order.
pub fn write_records(
    header: &SamHeader,
    records: &[BamRecord],
    create_index: bool,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ReadsError> {
    let writer = BamWriter::new(Vec::new(), header).map_err(|e| ReadsError::Io(e.to_string()))?;
    let mut writer = if create_index {
        writer.with_index()
    } else {
        writer
    };
    for record in records {
        writer
            .write(record)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
    }
    if create_index {
        let (bam, bai) = writer
            .finish_with_index()
            .map_err(|e| ReadsError::Io(e.to_string()))?;
        Ok((bam, Some(bai)))
    } else {
        let bam = writer.finish().map_err(|e| ReadsError::Io(e.to_string()))?;
        Ok((bam, None))
    }
}
