//! `AddCommentsToBam`, ported from `picard.sam.AddCommentsToBam` (Picard 3.4.0).
//!
//! A BAM whose header gains `@CO` lines and whose records are copied block for block. The copy is
//! [`htsjdk_bam::reheader::reheader_bam`], already byte-identical to `BamFileIoUtils.reheaderBamFile`.
//!
//! # The prefix is added once, and the tab is part of it
//!
//! `SAMFileHeader.addComment` prefixes `@CO\t` only when the comment does not already start with
//! it, so a comment that carries the prefix comes out exactly once. htsjdk-rs prefixed
//! unconditionally until this port measured the reference and the fix landed upstream.
//!
//! # The sam refusal reads the name, and the file is checked much later
//!
//! ```java
//! if (INPUT.getAbsolutePath().endsWith(".sam")) { throw new PicardException("SAM files are not supported"); }
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

/// The tool's own check, which reads the path and not the file.
pub fn is_refused_by_name(path: &str) -> bool {
    path.ends_with(".sam")
}

/// `doWork()`: the rewritten BAM.
///
/// `path` is the input's, because both refusals quote it: the first by its suffix and the second
/// by its whole name.
pub fn add_comments(bam: &[u8], path: &str, comments: &[String]) -> Result<Vec<u8>, CommentError> {
    if is_refused_by_name(path) {
        return Err(CommentError::SamNotSupported);
    }
    for comment in comments {
        if comment.contains('\n') {
            return Err(CommentError::ContainsNewline);
        }
    }
    let plain = decompress_all(bam).map_err(|_| CommentError::NotBlockCompressed {
        path: path.to_string(),
    })?;
    let reader =
        BamReader::new(&plain).map_err(|error| CommentError::Unreadable(format!("{error:?}")))?;
    let mut header: SamHeader = reader.header.text.clone();
    for comment in comments {
        header.add_comment(comment);
    }
    reheader_bam(&header, bam).map_err(|error| CommentError::Unreadable(format!("{error:?}")))
}
