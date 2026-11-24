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
    if create_index {
        return write_indexed(out, header, records, level, deflater);
    }
    let mut writer = BamWriter::with_compression(out, header, level, deflater)
        .map_err(|e| ReadsError::Io(e.to_string()))?;
    for record in records {
        writer
            .write(&record?)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
    }
    let bam = writer.finish().map_err(|e| ReadsError::Io(e.to_string()))?;
    Ok((bam, None))
}

/// [`write_stream`] with the `.bai` that `BamWriter::with_index` builds, byte for byte, without
/// its cost: that writer asks for the file pointer around every record, and a pointer waits for
/// every block being deflated, so the whole file was deflated on this thread. Here the blocks are
/// deflated in parallel and each record's pointers are worked out once its blocks have come out:
/// a record's place in the uncompressed stream is known as it is written, and each block's address
/// and uncompressed size are read off its BGZF frame as it reaches `out`.
fn write_indexed<W: std::io::Write>(
    out: W,
    header: &SamHeader,
    records: impl IntoIterator<Item = Result<BamRecord, ReadsError>>,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
) -> Result<(W, Option<Vec<u8>>), ReadsError> {
    let io = |e: std::io::Error| ReadsError::Io(e.to_string());
    // Where the first record starts: the header's uncompressed length, from the same writer's
    // header alone.
    let probe = BamWriter::with_compression(Vec::new(), header, 0, deflater)
        .and_then(BamWriter::finish)
        .map_err(io)?;
    let header_length = htsjdk_bgzf::decompress_all(&probe)
        .map_err(|e| ReadsError::Io(format!("{e:?}")))?
        .len() as u64;

    let blocks = std::rc::Rc::new(std::cell::RefCell::new(std::collections::VecDeque::new()));
    let spy = BlockSpy {
        inner: out,
        partial: Vec::new(),
        address: 0,
        uncompressed: 0,
        blocks: blocks.clone(),
    };
    let mut writer = BamWriter::with_compression(spy, header, level, deflater).map_err(io)?;
    let lengths: Vec<i32> = header.sequences.iter().map(|s| s.length).collect();
    let mut pending = Pointers {
        indexer: htsjdk_bam::index::BamIndexer::new(lengths.clone()),
        pending: std::collections::VecDeque::new(),
        blocks,
    };
    let mut position = header_length;
    let mut encoded = Vec::new();
    for record in records {
        let record = record?;
        writer
            .write(&record)
            .map_err(|e| ReadsError::Malformed(format!("{e:?}")))?;
        // `BamWriter::write` encoded it already; its length is all that is wanted here.
        encoded.clear();
        let _ = record.encode_into(&mut encoded);
        // `BAMRecordCodec.warnIfReferenceIsTooLargeForBinField`: bin 0 past the bin field's span.
        let too_large = record.reference_index >= 0
            && lengths
                .get(record.reference_index as usize)
                .is_some_and(|length| *length > htsjdk_bam::bin::BIN_GENOMIC_SPAN);
        let bin = if too_large || record.alignment_start == htsjdk_bam::bin::NO_ALIGNMENT_START {
            0
        } else {
            htsjdk_bam::compute_indexing_bin(record.alignment_start, record.alignment_end())
                .unwrap_or(0)
        };
        let start = position;
        position += encoded.len() as u64;
        pending.pending.push_back(Pending {
            reference_index: record.reference_index,
            alignment_start: record.alignment_start,
            alignment_end: record.alignment_end(),
            bin,
            unmapped: record.read_unmapped(),
            start,
            stop: position,
        });
        pending.resolve(false);
    }
    let spy = writer.finish().map_err(io)?;
    pending.resolve(true);
    debug_assert!(pending.pending.is_empty());
    Ok((spy.inner, Some(pending.indexer.finish())))
}

/// One BGZF block as it went out: where it starts, compressed and uncompressed, and its sizes.
struct Block {
    address: u64,
    uncompressed_start: u64,
    uncompressed_size: u64,
    compressed_size: u64,
}

/// The writer's output, passed through and cut into its BGZF blocks.
struct BlockSpy<W> {
    inner: W,
    /// The bytes of a block not yet whole.
    partial: Vec<u8>,
    address: u64,
    uncompressed: u64,
    blocks: std::rc::Rc<std::cell::RefCell<std::collections::VecDeque<Block>>>,
}

impl<W: std::io::Write> std::io::Write for BlockSpy<W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.inner.write_all(bytes)?;
        self.partial.extend_from_slice(bytes);
        let mut used = 0;
        // A BGZF frame: BSIZE (the block's length less one) at bytes 16-17, ISIZE in the last four.
        while self.partial.len() - used >= 18 {
            let frame = &self.partial[used..];
            let length = usize::from(u16::from_le_bytes([frame[16], frame[17]])) + 1;
            if frame.len() < length {
                break;
            }
            let size = u32::from_le_bytes(frame[length - 4..length].try_into().unwrap());
            self.blocks.borrow_mut().push_back(Block {
                address: self.address,
                uncompressed_start: self.uncompressed,
                uncompressed_size: u64::from(size),
                compressed_size: length as u64,
            });
            self.address += length as u64;
            self.uncompressed += u64::from(size);
            used += length;
        }
        self.partial.drain(..used);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// A written record waiting for the blocks its pointers fall in.
struct Pending {
    reference_index: i32,
    alignment_start: i32,
    alignment_end: i32,
    bin: i32,
    unmapped: bool,
    start: u64,
    stop: u64,
}

struct Pointers {
    indexer: htsjdk_bam::index::BamIndexer,
    pending: std::collections::VecDeque<Pending>,
    blocks: std::rc::Rc<std::cell::RefCell<std::collections::VecDeque<Block>>>,
}

impl Pointers {
    /// Hands the indexer every pending record whose pointers are known, in order. `finished`
    /// says the stream is closed, so its last block, which is the only one not full, is final.
    fn resolve(&mut self, finished: bool) {
        while let Some(record) = self.pending.front() {
            let (Some(start), Some(stop)) = (
                self.pointer(record.start, finished),
                self.pointer(record.stop, finished),
            ) else {
                return;
            };
            let record = self.pending.pop_front().unwrap();
            self.indexer.process(
                record.reference_index,
                record.alignment_start,
                record.alignment_end,
                record.bin,
                record.unmapped,
                htsjdk_bam::index::Chunk { start, end: stop },
            );
        }
    }

    /// `BlockCompressedOutputStream.getFilePointer()` when the stream stood at `position`. A full
    /// block is deflated the moment it fills, so the pointer at its end is the next block's start;
    /// the last block is deflated at the close, after its last pointer was taken, which is why
    /// that one points at its own end. Blocks no pointer can reach again are let go.
    fn pointer(&self, position: u64, finished: bool) -> Option<u64> {
        let mut blocks = self.blocks.borrow_mut();
        while blocks.len() > 1
            && blocks[0].uncompressed_start + blocks[0].uncompressed_size < position
        {
            blocks.pop_front();
        }
        let full = htsjdk_bgzf::DEFAULT_UNCOMPRESSED_BLOCK_SIZE as u64;
        for block in blocks.iter() {
            let end = block.uncompressed_start + block.uncompressed_size;
            let place = if position < end {
                Some((block.address, position - block.uncompressed_start))
            } else if position == end && block.uncompressed_size == full {
                Some((block.address + block.compressed_size, 0))
            } else if position == end && finished && block.uncompressed_size > 0 {
                Some((block.address, block.uncompressed_size))
            } else {
                None
            };
            if let Some((address, offset)) = place {
                return htsjdk_bgzf::vfp::make_file_pointer(address, offset as u32).ok();
            }
        }
        None
    }
}

/// What a BAM written to disk leaves: its `.bai` when one was asked for, and the MD5 of its bytes,
/// taken as they were written rather than read back.
pub struct Written {
    pub index: Option<Vec<u8>>,
    pub md5: String,
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
        md5: <md5::Md5 as md5::Digest>::new(),
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
    let md5 = format!(
        "{:032x}",
        u128::from_be_bytes(md5::Digest::finalize(out.md5).into())
    );
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
    md5: md5::Md5,
}

impl<W: std::io::Write> std::io::Write for Md5Writer<W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(bytes)?;
        md5::Digest::update(&mut self.md5, &bytes[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

#[cfg(test)]
mod indexed_writer_tests {
    use super::*;

    /// The `.bai` and the BAM are `BamWriter::with_index`'s byte for byte, over record counts that
    /// put the last record's end at every place in a block, an unmapped tail included.
    #[test]
    fn the_index_is_the_one_the_indexing_writer_builds() {
        let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n");
        for index in 0..900 {
            let length = 40 + (index * 37) % 211;
            let bases: String = (0..length)
                .map(|i| ['A', 'C', 'G', 'T'][(i + index) % 4])
                .collect();
            let qualities: String = (0..length)
                .map(|i| (b'#' + (i % 40) as u8) as char)
                .collect();
            if index < 880 {
                text.push_str(&format!(
                    "r{index}\t0\tchr1\t{}\t60\t{length}M\t*\t0\t0\t{bases}\t{qualities}\n",
                    1 + index * 97
                ));
            } else {
                text.push_str(&format!(
                    "u{index}\t4\t*\t0\t0\t*\t*\t0\t0\t{bases}\t{qualities}\n"
                ));
            }
        }
        let (header, records) = htsjdk_bam::sam_file::read_sam(&text).unwrap();
        for count in (0..=records.len()).step_by(7).chain([1, 2, records.len()]) {
            let mut expected =
                BamWriter::with_compression(Vec::new(), &header, 5, htsjdk_bgzf::Deflater::Jdk)
                    .unwrap()
                    .with_index();
            for record in &records[..count] {
                expected.write(record).unwrap();
            }
            let (expected_bam, expected_bai) = expected.finish_with_index().unwrap();
            let (bam, bai) = write_stream(
                Vec::new(),
                &header,
                records[..count].iter().cloned().map(Ok),
                true,
                5,
                htsjdk_bgzf::Deflater::Jdk,
            )
            .unwrap();
            assert_eq!(bam, expected_bam, "{count} records");
            assert_eq!(bai.unwrap(), expected_bai, "{count} records");
        }
    }
}
