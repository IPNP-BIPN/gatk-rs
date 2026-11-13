//! `ConvertCountsToDepthFile`, ported from
//! `org.broadinstitute.hellbender.tools.sv.ConvertCountsToDepthFile`, `SimpleCountCodec` and the
//! block-compressed branch of `FeatureOutputStream` (GATK 4.7.0.0, which adds the tool).
//!
//! A `.counts.tsv` rewritten as a one-sample `.rd.txt`: every `SimpleCount` becomes a
//! `DepthEvidence` with one count, written by `DepthEvidenceCodec`, so a count `1 100` comes out
//! `0 100`.
//!
//! # The tool's own sample refusal cannot be reached
//!
//! ```java
//! if (sampleName == null) { sampleName = header.getSampleName();
//!     if (sampleName == null) throw new UserException.BadInput("Sample name is not provided ...
//! ```
//!
//! The header is `SimpleCountCodec.readActualHeader`'s `SimpleSampleLocatableMetadata`, whose
//! constructor already refuses a null sample (`Utils.nonEmpty`) while the feature reader is being
//! opened. So a read group with no `SM` is refused before `onTraversalStart`, wrapped twice, and
//! `--sample-name` cannot rescue it; [`header_error_chain`] is that wrapping.
//!
//! # The name decides the format
//!
//! `FeatureOutputCodecFinder.find` picks the codec by the output's name and `FeatureOutputStream`
//! compresses when the name has a block-compressed extension, at `--compression-level`, with a
//! tabix index beside it built over the counts header's dictionary. A plain name ignores the level.

use crate::condense_depth_evidence::{self as depth, DepthEvidence};
use crate::print_read_counts::{read_sample_name, PrintError, SimpleCount};
use crate::sv_feature_codecs::{self as codecs, Encoding};
use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_tribble::tabix::{FeatureRef, TabixError, TabixFormat, TabixIndexCreator, ZERO_BASED};

/// `SimpleCountCollection`'s column line, which ends the header.
pub const COLUMN_HEADER: &str = "CONTIG\tSTART\tEND\tCOUNT";

/// `DepthEvidenceCodec.getTabixFormat()`.
pub const DEPTH_TABIX: TabixFormat = TabixFormat::new(ZERO_BASED, 1, 2, 3, b'#', 1);

/// What `SimpleCountCodec.readActualHeader` refuses, before the tool runs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HeaderError {
    /// `UserException.MalformedFile` with no file named, so the message says `Unknown file`.
    Malformed(&'static str),
    /// `MetadataUtils.fromHeader`'s sample refusals.
    Sample(PrintError),
}

impl HeaderError {
    pub fn java_class(&self) -> &'static str {
        match self {
            HeaderError::Malformed(_) => {
                "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile"
            }
            HeaderError::Sample(_) => "java.lang.IllegalArgumentException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            HeaderError::Malformed(problem) => format!("Unknown file is malformed: {problem}"),
            HeaderError::Sample(error) => error.message(),
        }
    }
}

/// The two wrappers a header refusal picks up on its way out of `FeatureDataSource`, outermost
/// first, the codec's own exception last.
pub fn header_error_chain(path: &str, inner: &HeaderError) -> Vec<(String, String)> {
    vec![
        (
            "org.broadinstitute.hellbender.exceptions.GATKException".to_string(),
            format!("Error initializing feature reader for path {path}"),
        ),
        (
            "htsjdk.tribble.TribbleException$MalformedFeatureFile".to_string(),
            format!(
                "Unable to parse header with error: {}, for input source: {path}",
                inner.message()
            ),
        ),
        (inner.java_class().to_string(), inner.message()),
    ]
}

/// A `.counts.tsv` as the codec reads it: the SAM header, its one sample, and the records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Counts {
    pub header: SamHeader,
    pub sample: String,
    pub records: Vec<SimpleCount>,
}

/// `SimpleCountCodec.readActualHeader`, then `decode` over every other line.
///
/// `Err(Ok(..))` is a header the codec refuses; `Err(Err(line))` is a record `decode` refuses,
/// which Tribble wraps in a message naming an iterator by its identity hash, so only the line is
/// kept.
pub fn read_counts(text: &str) -> Result<Counts, Result<HeaderError, String>> {
    let mut lines = text.lines();
    let mut sam = String::new();
    for line in lines.by_ref() {
        if line.starts_with('@') {
            sam.push_str(line);
            sam.push('\n');
            continue;
        }
        if sam.is_empty() {
            return Err(Ok(HeaderError::Malformed(
                "SAM header lines must be at the beginning of the file.",
            )));
        }
        if !line.starts_with(COLUMN_HEADER) {
            return Err(Ok(HeaderError::Malformed(
                "File does not have a column header.",
            )));
        }
        break;
    }
    let header = htsjdk_bam::reader::parse_header_text(&sam);
    // `new SimpleSampleMetadata(readSampleName(header))`, whose `Utils.nonEmpty` is the refusal a
    // read group with no SM reaches.
    let sample = match read_sample_name(&header) {
        Ok(Some(sample)) if !sample.is_empty() => sample,
        Ok(_) => return Err(Ok(HeaderError::Sample(PrintError::NullSampleName))),
        Err(error) => return Err(Ok(HeaderError::Sample(error))),
    };
    let mut records = Vec::new();
    for line in lines {
        if line.starts_with('@') || line.starts_with(COLUMN_HEADER) {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let number = |index: usize| {
            fields
                .get(index)
                .and_then(|field| field.parse::<i32>().ok())
        };
        match (fields.first(), number(1), number(2), number(3)) {
            (Some(contig), Some(start), Some(end), Some(count)) => records.push(SimpleCount {
                contig: contig.to_string(),
                start,
                end,
                count,
            }),
            _ => return Err(Err(line.to_string())),
        }
    }
    Ok(Counts {
        header,
        sample,
        records,
    })
}

/// What `onTraversalStart` refuses about the output's name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutputError {
    /// `FeatureOutputCodecFinder.find`, naming the path as given.
    NoOutputCodec { path: String },
    /// The tool's own `BadInput`, naming the path by `getURIString()`.
    WrongOutputType { uri: String },
}

impl OutputError {
    pub fn java_class(&self) -> &'static str {
        match self {
            OutputError::NoOutputCodec { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            OutputError::WrongOutputType { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
        }
    }

    pub fn message(&self) -> String {
        match self {
            OutputError::NoOutputCodec { path } => codecs::no_output_codec(path),
            OutputError::WrongOutputType { uri } => {
                format!("Bad input: {uri} is not a valid format for DepthEvidence")
            }
        }
    }
}

/// The finder's refusal and then the tool's; `Ok` is how the depth evidence will be laid out.
pub fn check_output(path: &str, uri: &str) -> Result<Encoding, OutputError> {
    match codecs::find(path) {
        Some(codec) if codec.feature_type == "DepthEvidence" => Ok(codec.encoding),
        Some(_) => Err(OutputError::WrongOutputType {
            uri: uri.to_string(),
        }),
        None => Err(OutputError::NoOutputCodec {
            path: path.to_string(),
        }),
    }
}

/// `apply`: one count, one sample, the same closed interval.
pub fn to_depth(count: &SimpleCount) -> DepthEvidence {
    DepthEvidence {
        contig: count.contig.clone(),
        start: count.start,
        end: count.end,
        counts: vec![count.count],
    }
}

/// The text output: `makeSink`'s header line, then `DepthEvidenceCodec.encode` per record.
pub fn write_text(sample: &str, records: &[DepthEvidence]) -> String {
    depth::write(&[sample.to_string()], records)
}

/// Why [`write_block_compressed`] stopped.
#[derive(Debug)]
pub enum BlockError {
    /// The index creator's refusal, which the reference raises from `write`.
    Tabix(TabixError),
    /// The deflater, which fails only where GKL's ISA-L levels cannot be reached on this host.
    Deflate(std::io::Error),
}

/// The block-compressed output and its `.tbi`, as `FeatureOutputStream` writes them.
///
/// Each feature is indexed at the virtual offset BEFORE it is written, and the index is closed at
/// the offset after `close()`, which is the end of the last data block: the terminator is written
/// by the codec underneath and never moves `mBlockAddress`. The data is deflated at the tool's
/// `--compression-level`, the index at the static default (`index_level`).
pub fn write_block_compressed(
    sample: &str,
    records: &[DepthEvidence],
    dictionary: &[SequenceRecord],
    level: u32,
    index_level: u32,
    deflater: htsjdk_bgzf::Deflater,
) -> Result<(Vec<u8>, Vec<u8>), BlockError> {
    use std::io::Write;
    let text = write_text(sample, records);
    let mut lines = text.split_inclusive('\n');
    let mut writer = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), level, deflater);
    let mut creator = TabixIndexCreator::new(DEPTH_TABIX);
    writer
        .write_all(lines.next().unwrap_or_default().as_bytes())
        .map_err(BlockError::Deflate)?;
    for (record, line) in records.iter().zip(lines) {
        let description = format!("{}:{}-{}", record.contig, record.start, record.end);
        // shortcut: a contig the dictionary lacks is a NullPointerException in the reference;
        // here it indexes with no length, upgrade if a golden ever reaches it.
        let sequence_length = dictionary
            .iter()
            .find(|sequence| sequence.name == record.contig)
            .map_or(0, |sequence| sequence.length);
        let offset = writer.file_pointer().map_err(BlockError::Deflate)? as i64;
        creator
            .add_feature(
                FeatureRef {
                    contig: &record.contig,
                    start: record.start,
                    end: record.end,
                    description: &description,
                    sequence_length,
                },
                offset,
            )
            .map_err(BlockError::Tabix)?;
        writer
            .write_all(line.as_bytes())
            .map_err(BlockError::Deflate)?;
    }
    let data = writer.into_inner().map_err(BlockError::Deflate)?;
    let end = ((data.len() - htsjdk_bgzf::EMPTY_GZIP_BLOCK.len()) as i64) << 16;
    let index = creator.finish(end).map_err(BlockError::Tabix)?;
    // `TabixIndex.write` is the body through a `BlockCompressedOutputStream` at the static level.
    let mut framed = htsjdk_bgzf::BgzfWriter::with_deflater(Vec::new(), index_level, deflater);
    framed
        .write_all(&index.write_body())
        .map_err(BlockError::Deflate)?;
    Ok((data, framed.into_inner().map_err(BlockError::Deflate)?))
}
