//! `CreateHadoopBamSplittingIndex`: the `.sbi` of a BAM, and the tool around htsjdk's writer.
//!
//! Reading the BAM is not ported. What is ported is where the index lands, which records it holds
//! an offset for, what its last entry is, and the bytes the writer lays those out in.
//!
//! Ported from `org.broadinstitute.hellbender.tools.spark.CreateHadoopBamSplittingIndex` and
//! `htsjdk.samtools.SBIIndexWriter` (GATK 4.7.0.0, where `CreateHadoopBamSplittingIndex.java` is
//! unchanged since 4.6.2.0, and htsjdk 5.0.0).
//!
//! The functions follow `createBaiAndSplittingIndex`, the `--create-bai` path, which drives
//! `SBIIndexWriter` one record at a time:
//!
//! ```java
//! BAMFileSpan filePointer = (BAMFileSpan) read.getFileSource().getFilePointer();
//! indexer.processRecord(filePointer.getFirstOffset());
//! ```
//!
//! # Reading this file without reading Rust
//!
//! - `u64` is an unsigned 64-bit integer; the index stores Java `long`s, and every value written
//!   here is non-negative, so the bit patterns agree.
//! - `Result<(), String>` is "nothing, or an error message": `Ok(())` passes, `Err(..)` is the
//!   exception the reference throws.
//! - `Option` is "a value or nothing"; `rfind` returns one (the index of the last match, or
//!   `None`), and `.map(..).unwrap_or(..)` transforms the value or falls back.
//! - `&name[a..b]` is a borrowed slice of text by byte index, Java's `substring(a, b)`.

/// `SBIIndexWriter.DEFAULT_GRANULARITY`: one offset every this many records.
pub const DEFAULT_GRANULARITY: u64 = 4096;

/// `SBIIndex.SBI_MAGIC`, whose fourth byte is a version rather than a letter:
///
/// ```java
/// static final byte[] SBI_MAGIC = "SBI\1".getBytes();
/// ```
///
/// `\1` in a Java string literal is the octal escape for the byte 1.
pub const SBI_MAGIC: [u8; 4] = [b'S', b'B', b'I', 1];

/// `FileExtensions.SBI`.
pub const SBI_EXTENSION: &str = ".sbi";
/// `FileExtensions.BAI_INDEX`.
pub const BAI_EXTENSION: &str = ".bai";

/// The message a granularity of nought or less is refused with, the last argument of
/// `new CommandLineException.BadArgumentValue(SPLITTING_INDEX_GRANULARITY_LONG_NAME, ..)`.
pub const GRANULARITY_MESSAGE: &str = "Granularity must be > 0";

/// `BlockCompressedFilePointerUtil.makeFilePointer`, on an offset with no block offset: the block
/// address occupies the top forty-eight bits and the low sixteen (the offset within the
/// uncompressed block) are zero. `<<` is the same shift as Java's.
pub fn make_file_pointer(block_address: u64) -> u64 {
    block_address << 16
}

/// `getOutputFile`, which APPENDS rather than replacing.
///
/// `reads.bam` becomes `reads.bam.sbi` and not `reads.sbi`, which is the opposite of what
/// `BuildBamIndex` does with an output argument of the same shape.
///
/// ```java
/// return new File(input.getPath() + FileExtensions.SBI);
/// ```
///
/// This is the branch with no `-O`; an explicit `-O` is used as given, and the reference only logs
/// a warning when it does not end in `bam.sbi`.
pub fn default_output(input: &str) -> String {
    format!("{input}{SBI_EXTENSION}")
}

/// `IOUtils.replaceExtension`, which is how the `.bai` companion is named from the index.
///
/// ```java
/// return FilenameUtils.removeExtension(path) + '.' + extensionNoLeadingDot;
/// ```
///
/// `FilenameUtils` is Apache Commons IO, which has no local source here. This cuts at the last `.`
/// anywhere in the text, including one in a directory name.
pub fn bai_companion(index: &str) -> String {
    match index.rfind('.') {
        Some(dot) => format!("{}{BAI_EXTENSION}", &index[..dot]),
        None => format!("{index}{BAI_EXTENSION}"),
    }
}

/// `assertIsBam`, whose message names the extension it found and not the file.
///
/// `BamFileIoUtils.isBamFile` is `SamReader.Type.BAM_TYPE.hasValidFileExtension(file.getName())`,
/// which is `fileName.endsWith("." + fileExtension())`, so the test is case-sensitive and `x.BAM`
/// is refused. The extension in the message is Guava's `Files.getFileExtension` in the reference
/// (no local source); here it is the text after the last `.`, or nothing.
pub fn assert_is_bam(name: &str) -> Result<(), String> {
    if name.ends_with(".bam") {
        return Ok(());
    }
    let extension = name.rfind('.').map(|dot| &name[dot + 1..]).unwrap_or("");
    Err(format!(
        "A splitting index is only relevant for a bam file, but a file with extension {extension} \
         was specified."
    ))
}

/// `doWork`'s first line, which happens before anything is opened.
///
/// ```java
/// if( granularity <= 0) {
/// ```
pub fn assert_granularity(granularity: i64) -> Result<(), String> {
    if granularity <= 0 {
        return Err(GRANULARITY_MESSAGE.to_string());
    }
    Ok(())
}

/// The entries the index holds, given every record's virtual offset.
///
/// One entry every `granularity` records, counting from the first, and then ONE MORE: the offset
/// the next record would have been written at. That last entry comes from the last record's chunk
/// end when there is a last record, and from the file's length when there is not, which is why an
/// empty BAM's last entry points past the file rather than inside it.
///
/// `processRecord` is `if (recordCount++ % granularity == 0)`, so the first record always has an
/// entry and the test reads the count BEFORE it is incremented: that is the zero-based `index`
/// below. `next_start` is the caller's `nextStart`, already defaulted to
/// [`make_file_pointer`] of the file length.
pub fn offsets(record_offsets: &[u64], granularity: u64, next_start: u64) -> Vec<u64> {
    // An iterator chain, Java's stream: number each offset, keep those whose index is a multiple
    // of the granularity, drop the index again, collect into a list.
    let mut entries: Vec<u64> = record_offsets
        .iter()
        .enumerate()
        .filter(|(index, _)| (*index as u64).is_multiple_of(granularity))
        .map(|(_, offset)| *offset)
        .collect();
    entries.push(next_start);
    entries
}

/// `SBIIndexWriter.finish`: the bytes of the index itself.
///
/// The md5 and the uuid the header has room for are written as zeroes. The writer has nowhere to
/// get either from, and nothing reads them back.
///
/// `finish(finalVirtualOffset, dataFileLength)` passes `null` for both, which become `EMPTY_MD5`
/// and `EMPTY_UUID`, sixteen zero bytes each. The layout follows the writer's calls in order:
///
/// ```java
/// binaryCodec.writeBytes(SBIIndex.SBI_MAGIC);
/// binaryCodec.writeLong(header.getFileLength());
/// binaryCodec.writeBytes(header.getMd5());
/// binaryCodec.writeBytes(header.getUuid());
/// binaryCodec.writeLong(header.getTotalNumberOfRecords());
/// binaryCodec.writeLong(header.getGranularity());
/// binaryCodec.writeLong(virtualOffsetCount);
/// ```
///
/// `BinaryCodec` writes little-endian, hence `to_le_bytes`. `entries` already holds the final
/// offset, so its length is `virtualOffsetCount`.
pub fn write(file_length: u64, records: u64, granularity: u64, entries: &[u64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(4 + 8 + 16 + 16 + 24 + entries.len() * 8);
    out.extend_from_slice(&SBI_MAGIC);
    out.extend_from_slice(&file_length.to_le_bytes());
    out.extend_from_slice(&[0u8; 16]);
    out.extend_from_slice(&[0u8; 16]);
    out.extend_from_slice(&records.to_le_bytes());
    out.extend_from_slice(&granularity.to_le_bytes());
    out.extend_from_slice(&(entries.len() as u64).to_le_bytes());
    for entry in entries {
        out.extend_from_slice(&entry.to_le_bytes());
    }
    out
}
