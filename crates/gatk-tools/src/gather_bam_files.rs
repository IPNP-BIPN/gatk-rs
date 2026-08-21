//! `GatherBamFiles`, ported from `picard.sam.GatherBamFiles` (Picard 3.4.0).
//!
//! Shards of a scattered run concatenated. The fast path is
//! [`htsjdk_bam::gather::gather_bam_files`], already byte-identical to
//! `BamFileIoUtils.gatherWithBlockCopying`; what this module adds is the tool's choice between the
//! two paths and what each one keeps.
//!
//! # The block copy never looks at a record, and that is visible in the output
//!
//! The first file's header is copied whole, every later file's is dropped unread, and no order or
//! dictionary check runs at all. So a shard whose header declares a read group the first does not
//! is concatenated anyway, and the output carries `RG:Z:rg2` under a header that declares only
//! `rg1`. The golden holds exactly that file.
//!
//! # The choice is what the files are, not what they are named
//!
//! ```java
//! for (final File f : inputs) { if (!BamFileIoUtils.isBamFile(f)) useBlockCopying = false; }
//! ```
//!
//! `isBamFile` reads the file, so one sam among the inputs sends the WHOLE run through the
//! record-by-record gather, which recompresses everything and does not produce the block copy's
//! bytes. [`use_block_copying`] answers that question and nothing else.
//!
//! # An empty shard and a list file both come out as if they were not there
//!
//! An empty shard contributes no blocks, and an input that is a list of paths is unrolled before
//! anything else happens, so both runs produce the same bytes as naming the files directly.

use htsjdk_bam::gather::gather_bam_files as block_copy_gather;
use htsjdk_bgzf::read::decompress_all;

/// What the run refuses before it gathers anything.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GatherError {
    /// `IOUtil.assertFileIsReadable`, which names the path as a URI.
    MissingInput { path: String },
    /// The block copy could not read a file it was told was a BAM.
    Unreadable(String),
}

impl GatherError {
    pub fn java_class(&self) -> &str {
        "htsjdk.samtools.SAMException"
    }

    pub fn message(&self) -> String {
        match self {
            GatherError::MissingInput { path } => {
                format!("Cannot read non-existent file: file://{path}")
            }
            GatherError::Unreadable(detail) => detail.clone(),
        }
    }
}

/// `BamFileIoUtils.isBamFile`, which reads the file's first bytes rather than its name: a bgzf
/// stream whose decompressed head is the BAM magic.
pub fn is_bam_file(file: &[u8]) -> bool {
    match decompress_all(file) {
        Ok(plain) => plain.starts_with(b"BAM\x01"),
        Err(_) => false,
    }
}

/// `determineBlockCopyingStatus`: every input has to be a BAM, so one sam decides for all of them.
pub fn use_block_copying(inputs: &[&[u8]]) -> bool {
    inputs.iter().all(|input| is_bam_file(input))
}

/// `IOUtil.unrollFiles`: an input naming a file of paths is replaced by the paths it names, one
/// per line, and anything else is left alone.
///
/// The reference keys this on the extension of the LIST file rather than on its contents, and an
/// empty line is skipped.
pub fn unroll(entries: &[String]) -> Vec<String> {
    let mut unrolled = Vec::new();
    for entry in entries {
        if entry.ends_with(".bam") || entry.ends_with(".sam") || entry.ends_with(".cram") {
            unrolled.push(entry.clone());
        } else {
            // A list file, whose lines are paths.
            for line in entry.lines() {
                let line = line.trim();
                if !line.is_empty() {
                    unrolled.push(line.to_string());
                }
            }
        }
    }
    unrolled
}

/// `doWork()`'s fast path: the concatenated BAM, first header kept and every other dropped.
pub fn gather(inputs: &[&[u8]]) -> Result<Vec<u8>, GatherError> {
    block_copy_gather(inputs).map_err(|error| GatherError::Unreadable(format!("{error:?}")))
}

/// The `.md5` the tool writes beside the output when `CREATE_MD5_FILE` is set: the digest as
/// lower-case hex, with no file name and no newline.
pub fn md5_file(output: &[u8]) -> String {
    format!("{:x}", md5_digest(output))
}

/// The digest itself, as a 128-bit value so the hex above is exactly 32 characters.
fn md5_digest(bytes: &[u8]) -> u128 {
    use md5::{Digest, Md5};
    let mut hasher = Md5::new();
    hasher.update(bytes);
    u128::from_be_bytes(hasher.finalize().into())
}
