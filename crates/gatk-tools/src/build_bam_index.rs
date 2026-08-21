//! `BuildBamIndex`, ported from `picard.sam.BuildBamIndex` (Picard 3.4.0).
//!
//! The `.bai` of a coordinate-sorted BAM. The index itself is
//! [`htsjdk_bam::build_index::build_bam_index`], already byte-identical to htsjdk's `BAMIndexer`;
//! what this module carries is the tool's own behaviour around it, which is two refusals and a
//! naming rule with a trap in it.
//!
//! # The default output lands beside the process, not beside the input
//!
//! ```java
//! final String baseFileName = inputPath.getFileName().toString();
//! OUTPUT = new File(baseFileName.substring(0, index) + FileExtensions.BAI_INDEX);
//! ```
//!
//! `getFileName()` drops the directory and `new File(name)` resolves against the working
//! directory, so indexing `shards/sorted.bam` with no OUTPUT writes `./sorted.bai` and leaves the
//! shard directory untouched. [`default_output`] answers the file name alone, and the caller is
//! the one that resolves it, which is where the reference resolves it too.
//!
//! # And the extension is replaced only for a name that ends `.bam`
//!
//! Anything else keeps its whole name and gains `.bai`, so `sorted.bam.copy` indexes to
//! `sorted.bam.copy.bai`.
//!
//! # The sort order is read from the header's claim and nothing else
//!
//! A header saying `queryname` and one saying `unsorted` are refused by the same message, and no
//! record is looked at. A file whose header says `coordinate` over records that are not is
//! indexed without complaint, which htsjdk's own writer will not produce and this port does not
//! guard against either.

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::reader::BamReader;
use htsjdk_bgzf::read::decompress_all;

/// What the tool refuses, both as the `SAMException` the reference raises.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IndexError {
    /// The reader's type is not BAM, which a `.sam` reaches whatever its name is.
    NotBam,
    /// The header's `SO` is not `coordinate`.
    NotCoordinateSorted,
    /// The file is a BAM the reader could not finish. No run in the golden reaches it: what the
    /// reference does with a truncated BAM is not measured here.
    Unreadable(String),
}

impl IndexError {
    pub fn java_class(&self) -> &str {
        match self {
            IndexError::NotBam | IndexError::NotCoordinateSorted => "htsjdk.samtools.SAMException",
            IndexError::Unreadable(_) => "htsjdk.samtools.SAMFormatException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            IndexError::NotBam => "Input file must be bam file, not sam file.".to_string(),
            IndexError::NotCoordinateSorted => {
                "Input bam file must be sorted by coordinate".to_string()
            }
            IndexError::Unreadable(detail) => detail.clone(),
        }
    }
}

/// `FileExtensions.BAI_INDEX`.
pub const BAI_INDEX: &str = ".bai";

/// The file name the tool writes when `OUTPUT` is not given, which is resolved against the
/// working directory and not against the input's.
pub fn default_output(input_file_name: &str) -> String {
    if let Some(stem) = input_file_name.strip_suffix(".bam") {
        format!("{stem}{BAI_INDEX}")
    } else {
        format!("{input_file_name}{BAI_INDEX}")
    }
}

/// Whether the bytes are a BAM at all, which is what `bam.type()` answers. A sam file, however it
/// is named, is not one.
pub fn is_bam(file: &[u8]) -> bool {
    match decompress_all(file) {
        Ok(plain) => plain.starts_with(b"BAM\x01"),
        Err(_) => false,
    }
}

/// `doWork()`: the index of one file, or the refusal it earns.
pub fn build(file: &[u8]) -> Result<Vec<u8>, IndexError> {
    // `bam.type()` is decided by what the reader could open, so a file that is not bgzf at all is
    // not a BAM either: a sam file fails to decompress and earns the same refusal as one that
    // decompresses without the magic.
    let plain = match decompress_all(file) {
        Ok(plain) => plain,
        Err(_) => return Err(IndexError::NotBam),
    };
    if !plain.starts_with(b"BAM\x01") {
        return Err(IndexError::NotBam);
    }
    let reader =
        BamReader::new(&plain).map_err(|error| IndexError::Unreadable(format!("{error:?}")))?;
    // `getSortOrder()` reads `SO` and defaults to `unsorted` when the header has none.
    let order = reader
        .header
        .text
        .attributes
        .get("SO")
        .unwrap_or("unsorted")
        .to_string();
    if order != "coordinate" {
        return Err(IndexError::NotCoordinateSorted);
    }
    build_bam_index(file).map_err(|error| IndexError::Unreadable(format!("{error:?}")))
}
