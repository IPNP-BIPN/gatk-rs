//! Ported from `org.broadinstitute.hellbender.tools.PrintReads` and the header handling in
//! `GATKTool.getHeaderForSAMWriter` (GATK 4.6.2.0).
//!
//! The first whole tool here, and the first whose claim is bytes. `PrintReads` is a `ReadWalker`
//! whose `apply` writes the read out, so the interesting part is not the traversal, which is
//! [`crate::read_walker`], but what happens to the header on the way.
//!
//! Two things there decide the output's bytes:
//!
//!  * a `@PG` record is appended, with `ID` = the **tool name**, `VN` = the GATK version, `CL` =
//!    the whole expanded command line and `PN` = the tool name again. The tool name is
//!    `GATK PrintReads`, with a space, while the command line begins `PrintReads`: they are not
//!    the same string, which is measurable and easy to get wrong;
//!  * the ID is made unique by appending `.1`, `.2` and so on until it is free, so running the
//!    tool over its own output produces `GATK PrintReads.1` rather than replacing the first
//!    record or emitting a duplicate.
//!
//! `getHeaderForSAMWriter` mutates the reads header in place rather than copying it, which is not
//! observable in one run's output and is reproduced anyway.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use htsjdk_bam::header::{ProgramRecord, SamHeader};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::writer::BamWriter;

/// `GATKTool.getToolName()` for this tool, which is not the string the command line starts with.
pub const TOOL_NAME: &str = "GATK PrintReads";

/// What the tool was asked to do. Only the arguments that change the output bytes are here.
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
///
/// The attribute order is the order the setters run in, which is what the writer emits: `VN`,
/// then `CL`, then `PN`, after the `ID` every `@PG` line leads with.
pub fn header_for_sam_writer(header: &SamHeader, options: &Options) -> SamHeader {
    let mut header = header.clone();
    if !options.add_output_sam_program_record {
        return header;
    }
    let mut record = ProgramRecord::new(&create_program_group_id(&header, TOOL_NAME));
    record.attributes.set("VN", options.version);
    record.attributes.set("CL", options.command_line);
    record.attributes.set("PN", TOOL_NAME);
    header.programs.push(record);
    header
}

/// `PrintReads`: the reads that survive the traversal, written back out.
///
/// Returns the BAM's bytes and, when asked for, its `.bai`.
pub fn print_reads(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ReadsError> {
    let records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    let header = header_for_sam_writer(source.header(), options);

    let writer = BamWriter::new(Vec::new(), &header).map_err(|e| ReadsError::Io(e.to_string()))?;
    let mut writer = if options.create_output_bam_index {
        writer.with_index()
    } else {
        writer
    };
    for record in &records {
        writer
            .write(record)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
    }
    if options.create_output_bam_index {
        let (bam, bai) = writer
            .finish_with_index()
            .map_err(|e| ReadsError::Io(e.to_string()))?;
        Ok((bam, Some(bai)))
    } else {
        let bam = writer.finish().map_err(|e| ReadsError::Io(e.to_string()))?;
        Ok((bam, None))
    }
}
