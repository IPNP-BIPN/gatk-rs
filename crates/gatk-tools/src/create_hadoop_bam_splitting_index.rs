//! `CreateHadoopBamSplittingIndex`: the `.sbi` of a BAM, and the tool around htsjdk's writer.
//!
//! Reading the BAM is not ported. What is ported is where the index lands, which records it holds
//! an offset for, what its last entry is, and the bytes the writer lays those out in.
//!
//! Ported from `org.broadinstitute.hellbender.tools.spark.CreateHadoopBamSplittingIndex` and
//! `htsjdk.samtools.SBIIndexWriter` in GATK 4.6.2.0.

/// `SBIIndexWriter.DEFAULT_GRANULARITY`: one offset every this many records.
pub const DEFAULT_GRANULARITY: u64 = 4096;

/// `SBIIndexWriter.SBI_MAGIC`, whose fourth byte is a version rather than a letter.
pub const SBI_MAGIC: [u8; 4] = [b'S', b'B', b'I', 1];

/// `FileExtensions.SBI`.
pub const SBI_EXTENSION: &str = ".sbi";
/// `FileExtensions.BAI_INDEX`.
pub const BAI_EXTENSION: &str = ".bai";

/// The message a granularity of nought or less is refused with.
pub const GRANULARITY_MESSAGE: &str = "Granularity must be > 0";

/// `BlockCompressedFilePointerUtil.makeFilePointer`, on an offset with no block offset: the block
/// address occupies the top forty-eight bits.
pub fn make_file_pointer(block_address: u64) -> u64 {
    block_address << 16
}

/// `getOutputFile`, which APPENDS rather than replacing.
///
/// `reads.bam` becomes `reads.bam.sbi` and not `reads.sbi`, which is the opposite of what
/// `BuildBamIndex` does with an output argument of the same shape.
pub fn default_output(input: &str) -> String {
    format!("{input}{SBI_EXTENSION}")
}

/// `IOUtils.replaceExtension`, which is how the `.bai` companion is named from the index.
pub fn bai_companion(index: &str) -> String {
    match index.rfind('.') {
        Some(dot) => format!("{}{BAI_EXTENSION}", &index[..dot]),
        None => format!("{index}{BAI_EXTENSION}"),
    }
}

/// `assertIsBam`, whose message names the extension it found and not the file.
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
pub fn offsets(record_offsets: &[u64], granularity: u64, next_start: u64) -> Vec<u64> {
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
