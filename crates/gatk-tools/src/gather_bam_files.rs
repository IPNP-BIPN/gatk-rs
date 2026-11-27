//! `GatherBamFiles`, ported from `picard.sam.GatherBamFiles` (Picard 3.5.0; `GatherBamFiles.java`
//! is unchanged since 3.4.0).
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
//! # The choice is made on the file NAMES, and this module's helper reads the bytes instead
//!
//! ```java
//! for (final File f : inputs) {
//!     if (!BamFileIoUtils.isBamFile(f)) {
//!         useBlockCopying = false;
//!     }
//! }
//! ```
//!
//! One input that is not a BAM sends the WHOLE run through the record-by-record gather, which
//! recompresses everything and does not produce the block copy's bytes. In htsjdk 5.0.0
//! `isBamFile` looks only at the name:
//!
//! ```java
//! return ((file != null) && SamReader.Type.BAM_TYPE.hasValidFileExtension(file.getName()));
//! ```
//!
//! [`is_bam_file`] and [`use_block_copying`] below read the bytes instead, so they agree with the
//! reference only when every name tells the truth (the golden's fixtures do). The command line
//! (`crates/gatk-cli/src/picard/gather_bam_files.rs`) makes the choice by the `.bam` suffix, as
//! the reference does.
//!
//! # An empty shard and a list file both come out as if they were not there
//!
//! An empty shard contributes no blocks, and an input that is a list of paths is unrolled before
//! anything else happens, so both runs produce the same bytes as naming the files directly.
//!
//! # Reading this file without reading Rust
//!
//! - `&[u8]` is a borrowed, read-only view of a byte array (Java's `byte[]`); `&[&[u8]]` is a
//!   list of such views, one per input file. The caller keeps ownership of the bytes.
//! - `Result<T, E>` is "a value or an error", how this codebase models a thrown exception, and
//!   `.map_err(..)` converts one error type into another on the way out.
//! - A function's last expression is its return value; `return` is written only for early exits.

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
    /// The class the reference throws, which is the same for both refusals.
    pub fn java_class(&self) -> &str {
        "htsjdk.samtools.SAMException"
    }

    /// The exception's message, as `getMessage()` would return it.
    pub fn message(&self) -> String {
        // `match` must name every variant; the compiler refuses one that forgets a case.
        match self {
            GatherError::MissingInput { path } => {
                format!("Cannot read non-existent file: file://{path}")
            }
            GatherError::Unreadable(detail) => detail.clone(),
        }
    }
}

/// Whether `file` is a BAM by its CONTENT: a bgzf stream whose decompressed head is the BAM magic.
///
/// This is not what `BamFileIoUtils.isBamFile` does: in htsjdk 5.0.0 that answers by the `.bam`
/// extension of the file's name (see the module note). The two agree whenever the names are
/// honest.
pub fn is_bam_file(file: &[u8]) -> bool {
    match decompress_all(file) {
        Ok(plain) => plain.starts_with(b"BAM\x01"),
        Err(_) => false,
    }
}

/// `determineBlockCopyingStatus`: every input has to be a BAM, so one sam decides for all of them.
///
/// Each input is judged by [`is_bam_file`], so by content rather than by the reference's name test.
pub fn use_block_copying(inputs: &[&[u8]]) -> bool {
    inputs.iter().all(|input| is_bam_file(input))
}

/// `IOUtil.unrollFiles`: an input naming a file of paths is replaced by the paths it names, one
/// per line, and anything else is left alone.
///
/// The reference keys this on the extension of the LIST file rather than on its contents, and an
/// empty line is skipped.
///
/// Each entry here is either a path ending in `.bam`, `.sam` or `.cram`, or the TEXT of a list
/// file (the caller has already read it). Unlike `IOUtil.unrollPaths`, which pushes every listed
/// path back on its stack and unrolls a list named inside a list, this goes one level deep only.
/// The command line has its own recursive `unroll` and does not call this one.
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
    // ZERO PADDED to thirty-two. `Md5CalculatingOutputStream.makeHash` prints the digest with
    // `new BigInteger(1, md5.digest()).toString(16)` and then left-pads it with zeros when
    // `hash.length() != 32`, and a digest whose first nibble is zero is thirty-one characters
    // without that padding: a covering-array row over ApplyBQSR caught one (`0f6124c6...` against
    // `f6124c6...`), and no golden had happened to hold one. `{:032x}` is that padding.
    format!("{:032x}", md5_digest(output))
}

/// [`md5_file`] of what `reader` yields, read in pieces rather than held whole.
///
/// `impl std::io::Read` accepts any readable source (a file, a byte slice); Java would take an
/// `InputStream`.
pub fn md5_of_reader(mut reader: impl std::io::Read) -> std::io::Result<String> {
    use md5::{Digest, Md5};
    let mut hasher = Md5::new();
    let mut buffer = vec![0u8; 1 << 16];
    loop {
        // A read of zero bytes is end of stream. An `Interrupted` error is retried rather than
        // returned, which is the Rust convention for a read the OS cut short.
        match reader.read(&mut buffer) {
            Ok(0) => break,
            Ok(read) => hasher.update(&buffer[..read]),
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(format!(
        "{:032x}",
        u128::from_be_bytes(hasher.finalize().into())
    ))
}

/// The digest itself, as a 128-bit value so the hex above is exactly 32 characters.
///
/// `from_be_bytes` reads the sixteen digest bytes most significant first, which is how
/// `new BigInteger(1, digest)` reads them.
fn md5_digest(bytes: &[u8]) -> u128 {
    use md5::{Digest, Md5};
    let mut hasher = Md5::new();
    hasher.update(bytes);
    u128::from_be_bytes(hasher.finalize().into())
}
