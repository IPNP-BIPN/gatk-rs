//! `AddCommentsToBam`, ported from `picard.sam.AddCommentsToBam` (Picard 3.5.0;
//! `AddCommentsToBam.java` is unchanged since Picard 3.4.0).
//!
//! A BAM whose header gains `@CO` lines and whose records are copied block for block. The copy is
//! [`htsjdk_bam::reheader::reheader_bam`], already byte-identical to
//! `BamFileIoUtils.reheaderBamFile`.
//!
//! # The prefix is added once, and the tab is part of it
//!
//! ```java
//! if (!comment.startsWith(SAMTextHeaderCodec.COMMENT_PREFIX)) {
//!     comment = SAMTextHeaderCodec.COMMENT_PREFIX + comment;
//! }
//! ```
//!
//! `SAMFileHeader.addComment` (htsjdk 5.0.0) prefixes `@CO\t` only when the comment does not
//! already start with it, so a comment that carries the prefix comes out exactly once. htsjdk-rs
//! prefixed unconditionally until this port measured the reference and the fix landed upstream.
//!
//! # The sam refusal reads the name, and the file is checked much later
//!
//! ```java
//! if (INPUT.getAbsolutePath().endsWith(".sam")) {
//!     throw new PicardException("SAM files are not supported");
//! }
//! ```
//!
//! That is the whole check: a BAM named `.sam` is refused, and a sam named `.bam` is not. The
//! second one fails afterwards, in the block copy, for having no valid GZIP block at its end.
//!
//! # The newline check the tool carries cannot be reached
//!
//! `doWork` refuses a comment holding a newline with a `PicardException`, but Picard's parser
//! refuses the argument first. [`CommentError::ContainsNewline`] is kept because it is the tool's
//! own behaviour for a caller that is not the command line, and the golden records the parser's
//! refusal instead.
//!
//! # Reading this file without reading Rust
//!
//! - `&[u8]` is a borrowed, read-only view of a byte array (Java's `byte[]`); `&str` is the same
//!   for a string. The caller keeps ownership.
//! - `Result<T, E>` is "a value or an error"; it is how this codebase models a Java exception.
//!   `Err(..)` is the throw, and `?` after a call passes that call's error straight up to the
//!   caller, as an uncaught exception would.
//! - `.map_err(|e| ..)` rewrites the error a call returned before `?` passes it on: here, into
//!   the variant naming the Java exception class.
//! - `match` must name every variant of an enum; the compiler refuses one that misses a case.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::reheader::reheader_bam;
use htsjdk_bgzf::read::decompress_all;

/// What the run refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CommentError {
    /// The input's path ends `.sam`, whatever the file holds.
    SamNotSupported,
    /// A comment holding a newline, which the command line never delivers.
    ContainsNewline,
    /// The block copy could not read the input.
    NotBlockCompressed { path: String },
    /// Anything else the copy refused.
    Unreadable(String),
}

impl CommentError {
    /// The fully qualified class of the exception the reference throws for this refusal: the
    /// tool's own two are `PicardException`, and the two the block copy raises are htsjdk's.
    pub fn java_class(&self) -> &str {
        match self {
            CommentError::SamNotSupported | CommentError::ContainsNewline => {
                "picard.PicardException"
            }
            CommentError::NotBlockCompressed { .. } | CommentError::Unreadable(_) => {
                "htsjdk.samtools.SAMException"
            }
        }
    }

    /// The exception's message, verbatim, as the reference's `getMessage()` would print it.
    pub fn message(&self) -> String {
        match self {
            CommentError::SamNotSupported => "SAM files are not supported".to_string(),
            CommentError::ContainsNewline => "Comments can not contain a new line".to_string(),
            CommentError::NotBlockCompressed { path } => {
                format!("file://{path} does not have a valid GZIP block at the end of the file.")
            }
            CommentError::Unreadable(detail) => detail.clone(),
        }
    }
}

/// The tool's own check, which reads the path and not the file:
/// `INPUT.getAbsolutePath().endsWith(".sam")`, case-sensitive, so `x.SAM` passes.
pub fn is_refused_by_name(path: &str) -> bool {
    path.ends_with(".sam")
}

/// `doWork()`: the rewritten BAM.
///
/// `path` is the input's, because both refusals quote it: the first by its suffix and the second
/// by its whole name.
///
/// The reference reads the header before the newline check (`getFileHeader(INPUT)` comes first
/// in `doWork`); here the newline check comes first. The two orders differ only for a caller that
/// passes a newline and an unreadable file, which the command line cannot (see the module note).
pub fn add_comments(bam: &[u8], path: &str, comments: &[String]) -> Result<Vec<u8>, CommentError> {
    if is_refused_by_name(path) {
        return Err(CommentError::SamNotSupported);
    }
    for comment in comments {
        if comment.contains('\n') {
            return Err(CommentError::ContainsNewline);
        }
    }
    // A file that does not decompress as bgzf is refused with the message
    // `reheaderBamFile`'s block copy gives, which quotes the path as a `file://` URL.
    let plain = decompress_all(bam).map_err(|_| CommentError::NotBlockCompressed {
        path: path.to_string(),
    })?;
    let reader =
        BamReader::new(&plain).map_err(|error| CommentError::Unreadable(format!("{error:?}")))?;
    // `clone()` copies the parsed header so it can be changed; the reader's own stays untouched.
    let mut header: SamHeader = reader.header.text.clone();
    for comment in comments {
        header.add_comment(comment);
    }
    // The original compressed bytes go to the copy, not `plain`: the records' blocks are copied
    // as they are, as `reheaderBamFile` does.
    reheader_bam(&header, bam).map_err(|error| CommentError::Unreadable(format!("{error:?}")))
}
