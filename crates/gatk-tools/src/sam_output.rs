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
            version: "4.7.0.0",
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
    write_records_with(
        header,
        records,
        create_index,
        htsjdk_bgzf::DEFAULT_COMPRESSION_LEVEL,
        htsjdk_bgzf::Deflater::Jdk,
    )
}

/// The same, with the BGZF compression named.
///
/// htsjdk writes with the JDK deflater at level five and GATK does not: `GATKConfig` sets
/// `samjdk.compression_level` to TWO and GATK replaces the static factory with Intel's GKL, so a
/// tool run from a command line writes different bytes from the same records
/// (`gatk-config`, #1032). The suites here call the first form because a dump reaches the tool
/// through `instanceMain` and never installs the config; a runner calls this one.
pub fn write_records_with(
    header: &SamHeader,
    records: &[BamRecord],
    create_index: bool,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
) -> Result<(Vec<u8>, Option<Vec<u8>>), ReadsError> {
    write_stream(
        Vec::new(),
        header,
        records.iter().cloned().map(Ok),
        create_index,
        level,
        deflater,
    )
}

/// The same, a record at a time into `out`: the BAM is written as the reads arrive and the `.bai`
/// is built beside it, so neither the reads nor the output are held.
///
/// A refusal from the traversal stops the writing where it came, and `out` holds what came
/// before it; a caller that must leave no file behind writes to a scratch path, as
/// [`write_bam_file`] does.
pub fn write_stream<W: std::io::Write>(
    out: W,
    header: &SamHeader,
    records: impl IntoIterator<Item = Result<BamRecord, ReadsError>>,
    create_index: bool,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
) -> Result<(W, Option<Vec<u8>>), ReadsError> {
    let writer = BamWriter::with_compression(out, header, level, deflater)
        .map_err(|e| ReadsError::Io(e.to_string()))?;
    let mut writer = if create_index {
        writer.with_index()
    } else {
        writer
    };
    for record in records {
        writer
            .write(&record?)
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

/// What a BAM written to disk leaves: its `.bai` when one was asked for, and the MD5 of its bytes
/// when one was, taken as they were written rather than read back.
pub struct Written {
    pub index: Option<Vec<u8>>,
    pub md5: Option<String>,
}

/// [`write_stream`] to `path`, through a scratch file beside it that is renamed over `path` only
/// once everything is written: a refusal partway leaves `path` as it was, which is what writing
/// the bytes only after the traversal had finished did.
pub fn write_bam_file(
    path: &std::path::Path,
    header: &SamHeader,
    records: impl IntoIterator<Item = Result<BamRecord, ReadsError>>,
    create_index: bool,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
    md5: bool,
) -> Result<Result<Written, std::io::Error>, ReadsError> {
    // A device, a pipe or a directory is written, or refused, where it is: there is nothing to
    // rename over it.
    let in_place = std::fs::metadata(path).is_ok_and(|metadata| !metadata.is_file());
    let scratch = if in_place {
        path.to_path_buf()
    } else {
        scratch_beside(path)
    };
    let file = match std::fs::File::create(&scratch) {
        Ok(file) => file,
        Err(error) => return Ok(Err(error)),
    };
    let out = Md5Writer {
        inner: std::io::BufWriter::with_capacity(1 << 20, file),
        // Hashing every byte costs as much as a tenth of the deflate, so only when asked.
        md5: md5.then(<md5::Md5 as md5::Digest>::new),
    };
    let written = write_stream(out, header, records, create_index, level, deflater);
    let (out, index) = match written {
        Ok(written) => written,
        Err(error) => {
            if !in_place {
                let _ = std::fs::remove_file(&scratch);
            }
            return Err(error);
        }
    };
    let md5 = out.md5.map(|md5| {
        format!(
            "{:032x}",
            u128::from_be_bytes(md5::Digest::finalize(md5).into())
        )
    });
    let finished = out
        .inner
        .into_inner()
        .map_err(|error| error.into_error())
        .and_then(|_| {
            if in_place {
                Ok(())
            } else {
                std::fs::rename(&scratch, path)
            }
        });
    if let Err(error) = finished {
        if !in_place {
            let _ = std::fs::remove_file(&scratch);
        }
        return Ok(Err(error));
    }
    Ok(Ok(Written { index, md5 }))
}

/// A name in `path`'s directory that no other writer uses.
fn scratch_beside(path: &std::path::Path) -> std::path::PathBuf {
    let name = path
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    path.with_file_name(format!(".{name}.{}.partial", std::process::id()))
}

/// The bytes on their way to the file, hashed as they pass.
struct Md5Writer<W> {
    inner: W,
    md5: Option<md5::Md5>,
}

impl<W: std::io::Write> std::io::Write for Md5Writer<W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(bytes)?;
        if let Some(md5) = &mut self.md5 {
            md5::Digest::update(md5, &bytes[..written]);
        }
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}
