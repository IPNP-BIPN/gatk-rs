//! `PrintReadCounts`, ported from `org.broadinstitute.hellbender.tools.sv.PrintReadCounts`
//! and the two codecs it accepts, `DepthEvidenceCodec` and `SimpleCountCodec` (GATK 4.6.2.0).
//!
//! A depth-evidence file, or a counts file, rewritten as one counts file per sample for the CNV
//! callers. The walk is a plain `FeatureWalker`; what is here is the two headers it can be handed,
//! the header it builds in return, and what the files on disk look like when the run does not
//! finish.
//!
//! # The two inputs disagree about what a header is
//!
//! An `.rd.txt` header is one line of column names, so `DepthEvidenceCodec.readActualHeader` can
//! name the samples but hands back a null dictionary:
//!
//! ```java
//! return new SVFeaturesHeader(DepthEvidence.class.getSimpleName(), "unknown", null,
//!                             headerCols.subList(3, headerCols.size()));
//! ```
//!
//! so the run needs `--sequence-dictionary` and refuses without one through a message that
//! misspells the argument it is asking for. A `.counts.tsv` carries a whole SAM header instead,
//! which names both its dictionary and, through its read groups, its one sample, and there
//! `--sequence-dictionary` is never consulted.
//!
//! # The output header is built, not copied
//!
//! Either way the output is a fresh `SAMFileHeader` over the dictionary plus one read group whose
//! ID is always `GATKCopyNumber`, then the column line. So a counts file fed back in loses its
//! `@PG`, its `@CO` and its read group's own ID, and keeps only its dictionary, its sample name and
//! its records. The two paths add the read group and the dictionary in opposite orders, which
//! makes no difference: `SAMTextHeaderCodec` writes `@HD`, then `@SQ`, then `@RG`, whatever order
//! they arrived in.
//!
//! # The coordinates change base and the records do not
//!
//! `DepthEvidenceCodec.decode` adds one to the start it reads, and the tool writes `getStart()`
//! rather than re-encoding, so a bin written `0 100` comes out `1 100`. A `.counts.tsv` is already
//! one-based, and `SimpleCountCodec.encode` gives back exactly what was read.
//!
//! # A run that does not finish leaves half a header
//!
//! Every writer is a `BufferedWriter` over a `FileWriter`, and `SAMTextHeaderCodec.encode` flushes
//! when it is done. The tool's own writes, the column line and every record, do not. So a crash
//! partway through leaves the `@HD`, `@SQ` and `@RG` lines on disk and nothing else, the column
//! line included. [`JavaWriter`] carries the 8192-character buffer that decides this. The
//! `FileWriter` underneath adds a second 8192-byte stage inside its `StreamEncoder`, which no case
//! measured so far reaches and which is not modelled here.
//!
//! # Nothing checks that the sample names are distinct
//!
//! Two identical column names in an `.rd.txt` header produce two writers over one path, and the
//! output file list names that path twice. Both `FileWriter`s truncate on open and each keeps its
//! own offset, so the second writer's header lands over the first's and, the two being the same
//! length, the surviving file is the second writer's throughout. [`Disk`] carries those offsets,
//! because a model that simply overwrites by name would agree here by accident and diverge as soon
//! as two samples wrote different amounts.

use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use std::collections::BTreeMap;

/// `MetadataUtils.GATK_CNV_READ_GROUP_ID`, which every output header carries whatever the input
/// read group was called.
pub const CNV_READ_GROUP_ID: &str = "GATKCopyNumber";

/// The line the tool writes itself, rather than through any codec.
pub const COLUMN_HEADER: &str = "CONTIG\tSTART\tEND\tCOUNT";

/// `BufferedWriter.defaultCharBufferSize`.
pub const BUFFER_CHARS: usize = 8192;

/// One count, one-based and closed, as a `.counts.tsv` holds it and as `SimpleCount` carries it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SimpleCount {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub count: i32,
}

/// `SimpleCountCodec.encode`.
pub fn encode_count(count: &SimpleCount) -> String {
    format!(
        "{}\t{}\t{}\t{}",
        count.contig, count.start, count.end, count.count
    )
}

/// One depth bin, one-based and closed as `DepthEvidence` carries it, which is NOT how the
/// `.rd.txt` line it came from was written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DepthEvidence {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub counts: Vec<i32>,
}

/// `DepthEvidenceCodec.decode`, which returns null for the header line and adds one to the start.
pub fn decode_depth(line: &str) -> Option<DepthEvidence> {
    if line.starts_with("#Chr") {
        return None;
    }
    let tokens: Vec<&str> = line.split('\t').collect();
    let counts = tokens[3..]
        .iter()
        .map(|token| token.parse::<i32>().unwrap_or(0))
        .collect();
    Some(DepthEvidence {
        contig: tokens[0].to_string(),
        // Adjust for 0-based indexing, as the reference puts it.
        start: tokens[1].parse::<i32>().unwrap_or(0) + 1,
        end: tokens[2].parse::<i32>().unwrap_or(0),
        counts,
    })
}

/// `DepthEvidenceCodec.readActualHeader`: the sample names are every column past the third, and
/// the dictionary is null.
pub fn depth_header_samples(line: &str) -> Vec<String> {
    line.split('\t').skip(3).map(str::to_string).collect()
}

/// An `.rd.txt` reduced to what the tool reads from it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DepthFile {
    pub samples: Vec<String>,
    pub records: Vec<DepthEvidence>,
}

/// A `.counts.tsv` reduced to what the tool reads from it: its own SAM header, and its records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CountsFile {
    pub header: SamHeader,
    pub records: Vec<SimpleCount>,
}

/// The driving feature file, whichever of the two acceptable types it holds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Input {
    Depth(DepthFile),
    Counts(CountsFile),
}

/// One `-L`, one-based and closed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Interval {
    fn overlaps(&self, contig: &str, start: i32, end: i32) -> bool {
        self.contig == contig && self.start <= end && start <= self.end
    }
}

/// What the run refuses, and what refuses it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrintError {
    /// `onTraversalStart` with an `.rd.txt` header and no `--sequence-dictionary`. The message
    /// misspells the argument.
    NoDictionary,
    /// `onTraversalStart` with a header that is neither kind.
    NoHeader,
    /// `MetadataUtils.readSampleName` with no read groups at all.
    NoReadGroups,
    /// `MetadataUtils.readSampleName` with read groups naming more than one sample.
    ManySampleNames(Vec<String>),
    /// `Utils.nonEmpty` on the null `readSampleName` gave back, which is what a read group with no
    /// `SM` produces. `readSampleName`'s own "does not contain a sample name" is unreachable from
    /// a header that has read groups: a lone null survives `distinct` and the emptiness test.
    NullSampleName,
    /// `FeatureWalker.initializeDrivingFeatures` on any other feature type.
    WrongFeatureType { path: String },
    /// A record with fewer counts than the header has samples, which is not a checked failure.
    IndexOutOfBounds { index: usize, length: usize },
}

impl PrintError {
    pub fn java_class(&self) -> &str {
        match self {
            PrintError::NoDictionary
            | PrintError::NoHeader
            | PrintError::WrongFeatureType { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            PrintError::NoReadGroups
            | PrintError::ManySampleNames(_)
            | PrintError::NullSampleName => "java.lang.IllegalArgumentException",
            PrintError::IndexOutOfBounds { .. } => "java.lang.ArrayIndexOutOfBoundsException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            // The reference misspells "dictionary" here, and the two spaces are its own.
            PrintError::NoDictionary => {
                "No dictionary available.  Supply one with --sequence-dictonary.".to_string()
            }
            PrintError::NoHeader => "Input file has no header.".to_string(),
            PrintError::NoReadGroups => {
                "The input header does not contain any read groups.  Cannot determine a sample name."
                    .to_string()
            }
            // StringUtils.join renders a null element as nothing at all.
            PrintError::ManySampleNames(names) => format!(
                "The input header contains more than one unique sample name: {}",
                names.join(", ")
            ),
            PrintError::NullSampleName => {
                "The string is null: string must not be null or empty".to_string()
            }
            PrintError::WrongFeatureType { path } => {
                format!("File {path} contains features of the wrong type.")
            }
            PrintError::IndexOutOfBounds { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
        }
    }
}

/// The two wrappers a header refusal picks up on its way out.
///
/// The header is parsed by the feature reader rather than by the tool, so `readSampleName`'s
/// `IllegalArgumentException` never reaches the caller as itself: tribble catches it and
/// `FeatureDataSource` catches that. The chain is what the run prints, outermost first.
pub fn reader_error_chain(path: &str, inner: &PrintError) -> Vec<(String, String)> {
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

/// `MetadataUtils.readSampleName`.
///
/// The distinct list is built from `SAMReadGroupRecord::getSample`, which is nullable, so a read
/// group with no `SM` contributes a null that survives both `distinct` and the emptiness test and
/// is handed back as the sample name.
pub fn read_sample_name(header: &SamHeader) -> Result<Option<String>, PrintError> {
    if header.read_groups.is_empty() {
        return Err(PrintError::NoReadGroups);
    }
    let mut distinct: Vec<Option<String>> = Vec::new();
    for group in &header.read_groups {
        let sample = group.attributes.get("SM").map(str::to_string);
        if !distinct.contains(&sample) {
            distinct.push(sample);
        }
    }
    if distinct.len() > 1 {
        return Err(PrintError::ManySampleNames(
            distinct
                .into_iter()
                .map(|name| name.unwrap_or_default())
                .collect(),
        ));
    }
    Ok(distinct.into_iter().next().flatten())
}

/// `createWriter`'s filename: raw concatenation, so a prefix that does not end in a separator
/// glues onto the sample name.
pub fn output_name(prefix: &str, sample: &str) -> String {
    format!("{prefix}{sample}.counts.tsv")
}

/// The header every output carries: the dictionary, and one read group whose ID is fixed.
pub fn output_header(sample: &str, sequences: &[SequenceRecord]) -> SamHeader {
    let mut header = SamHeader::new();
    header.sequences = sequences.to_vec();
    let mut group = ReadGroup::new(CNV_READ_GROUP_ID);
    group.attributes.set("SM", sample);
    header.read_groups.push(group);
    header
}

/// The files, each as a byte sequence with a length, so that two writers over one path behave the
/// way two file descriptors do rather than the way two map entries do.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Disk {
    files: BTreeMap<String, Vec<u8>>,
}

impl Disk {
    pub fn new() -> Self {
        Self::default()
    }

    /// `new FileWriter(path)`, which truncates.
    fn truncate(&mut self, path: &str) {
        self.files.insert(path.to_string(), Vec::new());
    }

    /// A write at this descriptor's own offset, which is not the file's length.
    fn write_at(&mut self, path: &str, offset: usize, bytes: &[u8]) {
        let file = self.files.entry(path.to_string()).or_default();
        if file.len() < offset {
            file.resize(offset, 0);
        }
        let end = offset + bytes.len();
        if file.len() < end {
            file.resize(end, 0);
        }
        file[offset..end].copy_from_slice(bytes);
    }

    pub fn read(&self, path: &str) -> Option<String> {
        self.files
            .get(path)
            .map(|bytes| String::from_utf8_lossy(bytes).into_owned())
    }

    /// Every file, in name order.
    pub fn files(&self) -> Vec<(String, String)> {
        self.files
            .iter()
            .map(|(name, bytes)| (name.clone(), String::from_utf8_lossy(bytes).into_owned()))
            .collect()
    }
}

/// A `BufferedWriter` over a `FileWriter`: 8192 characters held back, an offset of its own, and a
/// flush that only `SAMTextHeaderCodec` and `close` ever call.
#[derive(Debug, Clone)]
pub struct JavaWriter {
    path: String,
    offset: usize,
    buffer: String,
}

impl JavaWriter {
    /// `new FileWriter(path)`, truncating, wrapped in a `BufferedWriter`.
    pub fn create(disk: &mut Disk, path: &str) -> Self {
        disk.truncate(path);
        JavaWriter {
            path: path.to_string(),
            offset: 0,
            buffer: String::new(),
        }
    }

    pub fn write(&mut self, disk: &mut Disk, text: &str) {
        self.buffer.push_str(text);
        while self.buffer.chars().count() >= BUFFER_CHARS {
            let split = self
                .buffer
                .char_indices()
                .nth(BUFFER_CHARS)
                .map_or(self.buffer.len(), |(index, _)| index);
            let chunk: String = self.buffer.drain(..split).collect();
            self.push(disk, &chunk);
        }
    }

    pub fn new_line(&mut self, disk: &mut Disk) {
        self.write(disk, "\n");
    }

    pub fn flush(&mut self, disk: &mut Disk) {
        let chunk = std::mem::take(&mut self.buffer);
        self.push(disk, &chunk);
    }

    pub fn close(&mut self, disk: &mut Disk) {
        self.flush(disk);
    }

    fn push(&mut self, disk: &mut Disk, chunk: &str) {
        if chunk.is_empty() {
            return;
        }
        disk.write_at(&self.path, self.offset, chunk.as_bytes());
        self.offset += chunk.len();
    }
}

/// What one run leaves behind.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Run {
    /// Every file written, in name order, whether or not the run finished.
    pub disk: Disk,
    /// The names the tool chose, in sample order, which may repeat.
    pub names: Vec<String>,
    /// The refusal, if the run did not finish.
    pub error: Option<PrintError>,
}

/// The whole tool: `onTraversalStart`, then `apply` per record, then `onTraversalSuccess`.
///
/// `dictionary` is `--sequence-dictionary`, consulted only for an `.rd.txt` input. `list_path` is
/// `--output-file-list`, written and closed before a single record. `intervals` is `-L`, which
/// subsets the records and never the files.
pub fn run(
    input: &Input,
    dictionary: Option<&[SequenceRecord]>,
    output_prefix: &str,
    list_path: Option<&str>,
    intervals: &[Interval],
) -> Run {
    let mut disk = Disk::new();
    let (samples, sequences) = match input {
        Input::Depth(file) => match dictionary {
            Some(sequences) => (file.samples.clone(), sequences.to_vec()),
            None => {
                return Run {
                    disk,
                    names: Vec::new(),
                    error: Some(PrintError::NoDictionary),
                }
            }
        },
        Input::Counts(file) => match read_sample_name(&file.header) {
            Ok(Some(sample)) => (vec![sample], file.header.sequences.clone()),
            Ok(None) => {
                return Run {
                    disk,
                    names: Vec::new(),
                    error: Some(PrintError::NullSampleName),
                }
            }
            Err(error) => {
                return Run {
                    disk,
                    names: Vec::new(),
                    error: Some(error),
                }
            }
        },
    };

    // The list is opened, written and closed inside the try-with-resources that creates the
    // writers, so it is complete before any record is.
    let mut names = Vec::new();
    let mut writers = Vec::new();
    let mut list = list_path.map(|path| JavaWriter::create(&mut disk, path));
    for sample in &samples {
        let name = output_name(output_prefix, sample);
        if let Some(writer) = list.as_mut() {
            writer.write(&mut disk, &format!("{sample}\t{name}"));
            writer.new_line(&mut disk);
        }
        let mut writer = JavaWriter::create(&mut disk, &name);
        // SAMTextHeaderCodec.encode flushes when it is done, which is why a crash still leaves
        // the header on disk.
        writer.write(&mut disk, &output_header(sample, &sequences).encode());
        writer.flush(&mut disk);
        writer.write(&mut disk, COLUMN_HEADER);
        writer.new_line(&mut disk);
        names.push(name);
        writers.push(writer);
    }
    if let Some(mut writer) = list {
        writer.close(&mut disk);
    }

    let mut error = None;
    match input {
        Input::Counts(file) => {
            for record in &file.records {
                if !kept(intervals, &record.contig, record.start, record.end) {
                    continue;
                }
                writers[0].write(&mut disk, &encode_count(record));
                writers[0].new_line(&mut disk);
            }
        }
        Input::Depth(file) => {
            'records: for record in &file.records {
                if !kept(intervals, &record.contig, record.start, record.end) {
                    continue;
                }
                let interval_fields =
                    format!("{}\t{}\t{}\t", record.contig, record.start, record.end);
                for (index, writer) in writers.iter_mut().enumerate() {
                    let Some(count) = record.counts.get(index) else {
                        error = Some(PrintError::IndexOutOfBounds {
                            index,
                            length: record.counts.len(),
                        });
                        break 'records;
                    };
                    writer.write(&mut disk, &format!("{interval_fields}{count}"));
                    writer.new_line(&mut disk);
                }
            }
        }
    }

    // onTraversalSuccess, which a refusal never reaches, so the buffers die with it.
    if error.is_none() {
        for writer in writers.iter_mut() {
            writer.close(&mut disk);
        }
    }

    Run { disk, names, error }
}

/// An empty `-L` is no filtering at all, which is the traversal's own rule.
fn kept(intervals: &[Interval], contig: &str, start: i32, end: i32) -> bool {
    intervals.is_empty()
        || intervals
            .iter()
            .any(|interval| interval.overlaps(contig, start, end))
}
