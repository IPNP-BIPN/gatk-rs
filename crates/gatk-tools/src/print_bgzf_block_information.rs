//! `PrintBGZFBlockInformation`, ported from
//! `org.broadinstitute.hellbender.tools.PrintBGZFBlockInformation` (GATK 4.6.2.0).
//!
//! A bgzf file walked block by block. The tool reads the framing itself, adapted from
//! `BlockCompressedInputStream` rather than using it, so nothing here decompresses anything: the
//! report is offsets, sizes and the terminator check.
//!
//! # A premature terminator is reported twice, and always one block late
//!
//! ```java
//! if ( previousBlockInfo != null && previousBlockInfo.uncompressedSize == 0 ) {
//!     nonFinalTerminatorBlockIndices.add(blockNumber - 1);
//! ```
//!
//! The check fires while printing the block AFTER the terminator, so the banner is printed above
//! that block and names `blockNumber - 1`, which is the terminator's own number. The same numbers
//! come back at the end joined with a comma and no space.
//!
//! # The first line names the file and every refusal names the path
//!
//! `bgzfPath.getFileName()` reaches the report, while the `UserException`s carry the whole path,
//! which is why [`report`] takes the name and [`Refusal`] takes the path.

use htsjdk_bgzf::{BLOCK_HEADER_LENGTH, EMPTY_GZIP_BLOCK, MAX_COMPRESSED_BLOCK_SIZE};

/// One block, as the tool's own parser reads it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Block {
    pub offset: u64,
    /// The framed length, `BSIZE + 1`.
    pub compressed_size: usize,
    /// The `ISIZE` trailer, which a terminator block leaves at 0.
    pub uncompressed_size: i32,
}

/// What the run refuses, all of them `UserException.CouldNotReadInputFile`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// `onStartup`, before anything is opened.
    DoesNotExist { path: String },
    /// `IOUtil.isBlockCompressed` said no, which a regular gzip and a plain file both earn.
    NotBlockCompressed { path: String },
    /// The parser ran out of file inside a block.
    PrematureEndOfFile { path: String },
    /// A header that is present but too short.
    IncorrectHeaderSize { path: String },
    /// A BSIZE that cannot be a block.
    UnexpectedBlockLength { length: usize, path: String },
}

impl Refusal {
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException$CouldNotReadInputFile"
    }

    pub fn message(&self) -> String {
        match self {
            Refusal::DoesNotExist { path } => format!("File {path} does not exist"),
            Refusal::NotBlockCompressed { path } => format!(
                "Couldn't read file file://{path}. Error was: File is not a valid BGZF file. Could \
                 be a regular GZIP file, or some other non-BGZF format."
            ),
            Refusal::PrematureEndOfFile { path } => format!(
                "Couldn't read file. Error was: Error while parsing BGZF file. with exception: \
                 Premature end of file: {path}"
            ),
            Refusal::IncorrectHeaderSize { path } => format!(
                "Couldn't read file. Error was: Error while parsing BGZF file. with exception: \
                 Incorrect header size for file: {path}"
            ),
            Refusal::UnexpectedBlockLength { length, path } => format!(
                "Couldn't read file. Error was: Error while parsing BGZF file. with exception: \
                 Unexpected compressed block length: {length} for {path}"
            ),
        }
    }
}

/// `IOUtil.isBlockCompressed`, which is the startup check: a gzip whose extra field carries the
/// `BC` subfield of length 2.
pub fn is_block_compressed(file: &[u8]) -> bool {
    if file.len() < 18 {
        return false;
    }
    if file[0] != 31 || file[1] != 139 || file[2] != 8 || file[3] & 4 == 0 {
        return false;
    }
    let extra_length = u16::from_le_bytes([file[10], file[11]]) as usize;
    if extra_length < 6 || 12 + extra_length > file.len() {
        return false;
    }
    // Walk the extra subfields looking for `BC` with a payload of two bytes.
    let mut position = 12;
    let end = 12 + extra_length;
    while position + 4 <= end {
        let id = (file[position], file[position + 1]);
        let payload = u16::from_le_bytes([file[position + 2], file[position + 3]]) as usize;
        if id == (b'B', b'C') && payload == 2 {
            return true;
        }
        position += 4 + payload;
    }
    false
}

/// `processNextBlock` run to exhaustion: every block's frame, or the refusal the walk earned.
pub fn blocks(file: &[u8], path: &str) -> Result<Vec<Block>, Refusal> {
    let mut blocks = Vec::new();
    let mut offset = 0usize;
    while offset < file.len() {
        let available = file.len() - offset;
        if available < BLOCK_HEADER_LENGTH {
            // `readBytes` returned something short of a header, which is not EOF.
            return Err(Refusal::IncorrectHeaderSize {
                path: path.to_string(),
            });
        }
        let block_length = u16::from_le_bytes([file[offset + 16], file[offset + 17]]) as usize + 1;
        if !(BLOCK_HEADER_LENGTH..=MAX_COMPRESSED_BLOCK_SIZE).contains(&block_length) {
            return Err(Refusal::UnexpectedBlockLength {
                length: block_length,
                path: path.to_string(),
            });
        }
        if available < block_length {
            return Err(Refusal::PrematureEndOfFile {
                path: path.to_string(),
            });
        }
        let trailer = offset + block_length - 4;
        let uncompressed_size = i32::from_le_bytes([
            file[trailer],
            file[trailer + 1],
            file[trailer + 2],
            file[trailer + 3],
        ]);
        blocks.push(Block {
            offset: offset as u64,
            compressed_size: block_length,
            uncompressed_size,
        });
        offset += block_length;
    }
    Ok(blocks)
}

/// `doWork()`: the whole report, given the file's own name as the first line quotes it.
pub fn report(file: &[u8], file_name: &str, path: &str) -> (String, Option<Refusal>) {
    let mut text = format!("BGZF block information for file: {file_name}\n\n");
    // `doWork` asks `BlockCompressedInputStream.isValidFile` BEFORE it walks anything, so a file
    // that is not block compressed is refused for what it is rather than for how its first bytes
    // frame: a long plain file has enough bytes for a header and frames as a premature end, which
    // is the refusal the walk earns and not the one the tool gives.
    //
    // The check was in this module and nothing called it: the suite asked it directly and the
    // walk did not, so the two disagreed the moment a caller ran the tool rather than the suite.
    // A covering array run against the binary is what found that.
    if !is_block_compressed(file) {
        return (
            text,
            Some(Refusal::NotBlockCompressed {
                path: path.to_string(),
            }),
        );
    }
    let (found, refusal) = match blocks(file, path) {
        Ok(found) => (found, None),
        // The refusal comes after the blocks already printed, and those stay on disk.
        Err(refusal) => (blocks_before_failure(file), Some(refusal)),
    };

    let mut premature: Vec<usize> = Vec::new();
    let mut previous: Option<Block> = None;
    for (index, block) in found.iter().enumerate() {
        let number = index + 1;
        // The check runs while printing the block AFTER the terminator, and names that one's
        // predecessor.
        if let Some(previous) = previous {
            if previous.uncompressed_size == 0 {
                premature.push(number - 1);
                text.push_str("*******************************************************\n");
                text.push_str("ERROR: Premature BGZF 0-byte terminator block was found\n");
                text.push_str(&format!("at block number: {}\n", number - 1));
                text.push_str("*******************************************************\n");
                text.push('\n');
            }
        }
        text.push_str(&format!(
            "Block #{number} at file offset {}\n",
            block.offset
        ));
        text.push_str(&format!("\t- compressed size: {}\n", block.compressed_size));
        text.push_str(&format!(
            "\t- uncompressed size: {}\n",
            block.uncompressed_size
        ));
        text.push('\n');
        previous = Some(*block);
    }

    // A run that failed mid-file never reaches the banners: the exception leaves doWork.
    if refusal.is_some() {
        return (text, refusal);
    }

    match previous {
        Some(block) if block.uncompressed_size == 0 => {
            text.push_str(
                "***************************************************************************\n",
            );
            text.push_str(&format!(
                "Final BGZF 0-byte terminator block FOUND as expected at block number {}\n",
                found.len()
            ));
            text.push_str(
                "***************************************************************************\n",
            );
            text.push('\n');
        }
        _ => {
            text.push_str("******************************************************\n");
            text.push_str("ERROR: Final BGZF 0-byte terminator block was MISSING!\n");
            text.push_str("******************************************************\n");
            text.push('\n');
        }
    }

    if !premature.is_empty() {
        text.push_str("***********************************************************\n");
        text.push_str("ERROR: Premature BGZF 0-byte terminator block(s) were found\n");
        text.push_str(&format!(
            "at block number(s): {}\n",
            premature
                .iter()
                .map(|number| number.to_string())
                .collect::<Vec<String>>()
                .join(",")
        ));
        text.push_str("***********************************************************\n");
        text.push('\n');
    }

    (text, refusal)
}

/// The blocks the walk managed to print before it ran out of file.
fn blocks_before_failure(file: &[u8]) -> Vec<Block> {
    let mut blocks = Vec::new();
    let mut offset = 0usize;
    while offset + BLOCK_HEADER_LENGTH <= file.len() {
        let block_length = u16::from_le_bytes([file[offset + 16], file[offset + 17]]) as usize + 1;
        if !(BLOCK_HEADER_LENGTH..=MAX_COMPRESSED_BLOCK_SIZE).contains(&block_length)
            || offset + block_length > file.len()
        {
            break;
        }
        let trailer = offset + block_length - 4;
        blocks.push(Block {
            offset: offset as u64,
            compressed_size: block_length,
            uncompressed_size: i32::from_le_bytes([
                file[trailer],
                file[trailer + 1],
                file[trailer + 2],
                file[trailer + 3],
            ]),
        });
        offset += block_length;
    }
    blocks
}

/// The terminator block's own length, which the report prints as a compressed size of 28.
pub const TERMINATOR_LENGTH: usize = EMPTY_GZIP_BLOCK.len();
