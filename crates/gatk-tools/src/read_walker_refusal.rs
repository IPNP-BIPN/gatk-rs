//! What a read walker does with a file that is not a BAM.
//!
//! `ReadsDataSource` opens whatever it is given and only then finds out what it is, so the
//! refusals a walker makes are the READER'S rather than the tool's, and there are three of them
//! with different statuses.
//!
//! # The three, and the two that are not refusals
//!
//!  * **a path that does not exist, and a path that is a directory**, are
//!    `UserException$CouldNotReadInputFile`, which is status two. The message wraps htsjdk's inside
//!    GATK's, and htsjdk's names the file as a URI with a trailing slash for a directory;
//!  * **a file that is not a BAM at all** is read as a TEXT SAM stream, and the first line that is
//!    not a header runs out of fields: `SAMFormatException`, which is no `UserException` at all,
//!    so `Main` reports it as a non-user failure and the status is THREE. A block-compressed file
//!    is the same refusal without the file's name in it, because the reader is looking at a
//!    decompressed stream by then;
//!  * **the same file with an INTERVAL** is a different refusal again: an interval needs a
//!    sequence dictionary, the dictionary is asked for before any record is read, and an empty one
//!    is an `IllegalArgumentException`. One file, two refusals, decided by an argument that has
//!    nothing to do with the file;
//!  * **an empty file is not refused**: it is a BAM with no records, and the tool returns zero;
//!  * **and a BAM is a BAM.**
//!
//! Ported from `htsjdk.samtools.SamReaderFactory`, `htsjdk.samtools.SAMTextReader`,
//! `org.broadinstitute.hellbender.engine.ReadsPathDataSource` and
//! `org.broadinstitute.hellbender.exceptions.UserException` as measured in
//! `read-walker-refusals`.

/// What the reference throws, and whether `Main` calls it the user's fault.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// `UserException$CouldNotReadInputFile`, status two.
    CouldNotRead { path: String, reason: String },
    /// `SAMFormatException`, status three: the stream is not SAM text either.
    NotSamText {
        /// The file's name where the reader still had it, which a decompressed stream does not.
        file: Option<String>,
        line: String,
    },
    /// `IllegalArgumentException`, status three: an interval was given and the dictionary is empty.
    EmptyDictionary,
}

pub const COULD_NOT_READ: &str =
    "org.broadinstitute.hellbender.exceptions.UserException$CouldNotReadInputFile";
pub const SAM_FORMAT: &str = "htsjdk.samtools.SAMFormatException";
pub const ILLEGAL_ARGUMENT: &str = "java.lang.IllegalArgumentException";

impl Refusal {
    /// The exception class the reference throws.
    pub fn exception(&self) -> &'static str {
        match self {
            Refusal::CouldNotRead { .. } => COULD_NOT_READ,
            Refusal::NotSamText { .. } => SAM_FORMAT,
            Refusal::EmptyDictionary => ILLEGAL_ARGUMENT,
        }
    }

    /// Whether it is a `UserException`, which is what decides the status: two, or three.
    pub fn is_user(&self) -> bool {
        matches!(self, Refusal::CouldNotRead { .. })
    }

    /// The message, in the reference's own words.
    pub fn message(&self) -> String {
        match self {
            Refusal::CouldNotRead { path, reason } => {
                format!("Couldn't read file. Error was: {path} with exception: {reason}")
            }
            Refusal::NotSamText { file, line } => {
                let named = match file {
                    Some(name) => format!(" File {name};"),
                    None => String::new(),
                };
                format!(
                    "Error parsing text SAM file. Not enough fields;{named} Line 1\nLine: {line}"
                )
            }
            Refusal::EmptyDictionary => "Dictionary cannot have size zero".to_string(),
        }
    }
}

/// htsjdk's own wording for a path it cannot open, which GATK wraps rather than replaces.
pub fn cannot_read(path: &str, is_directory: bool) -> String {
    if is_directory {
        format!("Cannot read file because it is a directory: file://{path}/")
    } else {
        format!("Cannot read non-existent file: file://{path}")
    }
}

/// The BAM magic, which is what a reader looks at first once the stream is decompressed.
pub const BAM_MAGIC: [u8; 4] = *b"BAM\x01";

/// Whether the bytes are a BGZF member, which the reader decompresses before deciding anything.
pub fn is_block_compressed(bytes: &[u8]) -> bool {
    bytes.len() >= 4 && bytes[0] == 0x1f && bytes[1] == 0x8b && bytes[2] == 0x08
}

/// What a read walker makes of a file, given its bytes and whether an interval was asked for.
///
/// `decompressed` is what the reader would be looking at: the file's own bytes for a plain file,
/// and the inflated stream for a block-compressed one. `None` where the caller could not inflate
/// it, which is a file that claims to be BGZF and is not.
pub fn refusal(
    path: &str,
    exists: bool,
    is_directory: bool,
    decompressed: Option<&[u8]>,
    was_compressed: bool,
    has_intervals: bool,
) -> Option<Refusal> {
    if !exists || is_directory {
        return Some(Refusal::CouldNotRead {
            path: path.to_string(),
            reason: cannot_read(path, is_directory),
        });
    }
    let bytes = decompressed?;
    if bytes.starts_with(&BAM_MAGIC) {
        return None;
    }
    if has_intervals {
        // The dictionary is asked for before the first record, and a stream that is not a BAM has
        // none. That includes an EMPTY file, which is a stream of no records with no header
        // either: it is accepted where no interval was given and refused where one was.
        return Some(Refusal::EmptyDictionary);
    }
    if bytes.is_empty() {
        // Nothing to parse and nothing to refuse: the tool returns zero.
        return None;
    }
    let line = bytes
        .split(|byte| *byte == b'\n')
        .next()
        .map(|line| {
            String::from_utf8_lossy(line)
                .trim_end_matches('\r')
                .to_string()
        })
        .unwrap_or_default();
    Some(Refusal::NotSamText {
        // A decompressed stream has lost the file's name by the time the reader complains.
        file: if was_compressed {
            None
        } else {
            Some(path.to_string())
        },
        line,
    })
}
