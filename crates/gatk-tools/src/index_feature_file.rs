//! `IndexFeatureFile`, ported from `org.broadinstitute.hellbender.tools.IndexFeatureFile`
//! (GATK 4.6.2.0).
//!
//! The tool every dump in this repository's harness already calls to make its fixtures. What it
//! builds is decided by the file's NAME, three ways, and the index writing itself is htsjdk-rs's
//! [`htsjdk_tribble::index_write`], already measured against htsjdk.
//!
//! # The extension chooses the index, and the same records give different bytes
//!
//! ```java
//! if (IOUtil.hasBlockCompressedExtension(featurePath.toPath())) { ... TABIX ... }
//! else if (featurePath.getURIString().endsWith(GVCF_FILE_EXTENSION)) { ... createLinearIndex(128000) }
//! else { ... createDynamicIndex(FOR_SEEK_TIME) }
//! ```
//!
//! So one file copied to `reads.vcf` and to `reads.g.vcf` indexes to two different files: a dynamic
//! index that chose its own layout, and a linear index with a bin width of 128000.
//!
//! # The index carries the source file's identity, mtime included
//!
//! The header holds the file's URI, its size and its `lastModified` in milliseconds, so indexing
//! the same bytes twice never gives the same file. The golden zeroes those eight bytes, and so does
//! [`Source::timestamp`] by default: a caller that wants the reference's own bytes has to supply
//! the real mtime.
//!
//! # Tabix is measured but not built here
//!
//! A block compressed input takes the tabix branch, which needs a tabix writer htsjdk-rs does not
//! have yet. [`index_kind`] answers [`IndexKind::Tabix`] for those, [`default_output`] names the
//! `.tbi`, and the refusal that guards its extension is ported; the bytes are not.

use htsjdk_tribble::index::{IntervalChrIndex, TribbleIndex, INTERVAL_TREE, LINEAR, VERSION};
use htsjdk_tribble::index_write::{
    BalanceApproach, BuiltIndex, DynamicIndexCreator, Feature, LinearIndexCreator,
};

/// `OPTIMAL_GVCF_INDEX_BIN_SIZE`.
pub const OPTIMAL_GVCF_INDEX_BIN_SIZE: i32 = 128_000;
/// `GVCF_FILE_EXTENSION`.
pub const GVCF_FILE_EXTENSION: &str = ".g.vcf";

/// Which index the file's name asks for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    /// `createDynamicIndex(FOR_SEEK_TIME)`, which is everything not named otherwise.
    Dynamic,
    /// `createLinearIndex(OPTIMAL_GVCF_INDEX_BIN_SIZE)`.
    Linear,
    /// `IndexType.TABIX`, whose bytes this port does not write.
    Tabix,
}

/// The codecs this port knows, which is what decides whether a file can be indexed at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Codec {
    Vcf,
    Bed,
}

/// What the run refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// The input is not readable, which for the reference means it does not exist.
    CouldNotReadInputFile { path: String },
    /// `FeatureManager.getCodecForFile` found nothing.
    NoSuitableCodecs { path: String },
    /// A block compressed input with an output that is not a `.tbi`.
    WrongIndexExtension { path: String },
    /// The features are not in order, which Tribble raises and the tool wraps.
    CouldNotIndexFile { path: String, detail: String },
    /// THE PORT'S OWN refusal, which the reference never makes: a block-compressed input needs a
    /// tabix index, and this port has no writer for one.
    ///
    /// It is a separate variant because the alternative was borrowing
    /// [`Refusal::WrongIndexExtension`], and that put the reference's words on the port's gap: a
    /// covering array run against the binary reported the two as the same answer, and the row
    /// where the reference writes a tabix index and the port cannot read as a refusal the
    /// reference agrees with. A gap has to say it is one.
    TabixIsNotWritten { path: String },
}

impl Refusal {
    pub fn java_class(&self) -> &str {
        match self {
            Refusal::CouldNotReadInputFile { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$CouldNotReadInputFile"
            }
            Refusal::NoSuitableCodecs { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$NoSuitableCodecs"
            }
            Refusal::WrongIndexExtension { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            Refusal::CouldNotIndexFile { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$CouldNotIndexFile"
            }
            // No Java class: the reference does not make this refusal, and naming one of its
            // exceptions here would be a claim about the reference.
            Refusal::TabixIsNotWritten { .. } => "",
        }
    }

    pub fn message(&self) -> String {
        match self {
            Refusal::CouldNotReadInputFile { path } => {
                format!("Couldn't read file file://{path}")
            }
            Refusal::NoSuitableCodecs { path } => {
                format!("Cannot read file://{path} because no suitable codecs found")
            }
            Refusal::WrongIndexExtension { path } => {
                format!("The index for {path} must be written to a file with a \".tbi\" extension")
            }
            Refusal::CouldNotIndexFile { path, detail } => {
                format!("Error while trying to create index for {path}. Error was: {detail}")
            }
            Refusal::TabixIsNotWritten { path } => format!(
                "{path} needs a tabix index, which this port does not write yet. This message is \
                 the port's own and not GATK's."
            ),
        }
    }
}

/// How the tool renders a wrapped Tribble exception: `UserException.CouldNotIndexFile` builds its
/// message from the cause's class NAME WITH DOTS, so the nested class comes out
/// `htsjdk.tribble.TribbleException.MalformedFeatureFile` and not with the `$` a JVM would print.
fn wrapped(error: &htsjdk_tribble::index_write::CreateError) -> String {
    format!("{}: {}", error.class().replace('$', "."), error.message())
}

/// The file whose identity the index header records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Source {
    /// The absolute path, which becomes `file://<path>` in the header.
    pub path: String,
    /// `lastModified()` in milliseconds. The golden masks it, so this defaults to zero and a
    /// caller that wants the reference's bytes supplies the real one.
    pub timestamp: i64,
}

impl Source {
    pub fn new(path: &str) -> Self {
        Source {
            path: path.to_string(),
            timestamp: 0,
        }
    }
}

/// `createAppropriateIndexInMemory`'s choice, made from the name alone.
pub fn index_kind(file_name: &str) -> IndexKind {
    if file_name.ends_with(".gz") || file_name.ends_with(".bgz") {
        IndexKind::Tabix
    } else if file_name.ends_with(GVCF_FILE_EXTENSION) {
        IndexKind::Linear
    } else {
        IndexKind::Dynamic
    }
}

/// `determineFileName` with no OUTPUT given: `Tribble.indexPath` or `Tribble.tabixIndexPath`, both
/// of which APPEND to the whole name rather than replacing anything.
pub fn default_output(file_name: &str) -> String {
    match index_kind(file_name) {
        IndexKind::Tabix => format!("{file_name}.tbi"),
        _ => format!("{file_name}.idx"),
    }
}

/// `FeatureManager.getCodecForFile`, for the two formats this port reads.
pub fn codec_for(file_name: &str) -> Option<Codec> {
    let stripped = file_name
        .strip_suffix(".gz")
        .or_else(|| file_name.strip_suffix(".bgz"))
        .unwrap_or(file_name);
    if stripped.ends_with(".vcf") {
        Some(Codec::Vcf)
    } else if stripped.ends_with(".bed") {
        Some(Codec::Bed)
    } else {
        None
    }
}

/// Every feature of the file with the byte offset its line starts at, which is what the creators
/// index against.
pub fn features(text: &str, codec: Codec) -> Vec<(Feature, i64)> {
    let mut found = Vec::new();
    let mut offset: i64 = 0;
    for line in text.split_inclusive('\n') {
        let start_of_line = offset;
        offset += line.len() as i64;
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("track") {
            continue;
        }
        let columns: Vec<&str> = trimmed.split('\t').collect();
        let feature = match codec {
            Codec::Vcf => {
                if columns.len() < 8 {
                    continue;
                }
                let start: i32 = match columns[1].parse() {
                    Ok(start) => start,
                    Err(_) => continue,
                };
                // `VariantContext.getEnd()`: the reference allele's span, or END when it is given.
                let end = end_attribute(columns[7]).unwrap_or(start + columns[3].len() as i32 - 1);
                Feature {
                    contig: columns[0].to_string(),
                    start,
                    end,
                }
            }
            Codec::Bed => {
                if columns.len() < 3 {
                    continue;
                }
                let start: i32 = match columns[1].parse() {
                    Ok(start) => start,
                    Err(_) => continue,
                };
                let end: i32 = match columns[2].parse() {
                    Ok(end) => end,
                    Err(_) => continue,
                };
                // BEDCodec turns a half-open zero-based interval into a closed one-based one.
                Feature {
                    contig: columns[0].to_string(),
                    start: start + 1,
                    end,
                }
            }
        };
        found.push((feature, start_of_line));
    }
    found
}

fn end_attribute(info: &str) -> Option<i32> {
    info.split(';')
        .find_map(|field| field.strip_prefix("END="))
        .and_then(|value| value.parse().ok())
}

/// `IndexFactory.createIndex` for the two branches that write a `.idx`, then `index.write(path)`.
pub fn build(text: &str, source: &Source, file_name: &str) -> Result<Vec<u8>, Refusal> {
    let codec = codec_for(file_name).ok_or_else(|| Refusal::NoSuitableCodecs {
        path: source.path.clone(),
    })?;
    let found = features(text, codec);
    let ordered: Vec<Feature> = found.iter().map(|(feature, _)| feature.clone()).collect();
    // Tribble's own complaint names the plain path here, where the index header holds a URI.
    htsjdk_tribble::index_write::check_ordering(&ordered, &source.path).map_err(|error| {
        Refusal::CouldNotIndexFile {
            path: source.path.clone(),
            detail: wrapped(&error),
        }
    })?;

    let indexed_path = format!("file://{}", source.path);
    let file_size = text.len() as i64;
    let (index_type, properties, contigs, interval_contigs) = match index_kind(file_name) {
        IndexKind::Tabix => {
            // The caller is expected to have taken the tabix branch itself: this port has no
            // writer for it, and the golden records those bytes for a later brick. The refusal is
            // the PORT'S own, so that nothing downstream can read it as the reference's.
            return Err(Refusal::TabixIsNotWritten {
                path: source.path.clone(),
            });
        }
        IndexKind::Linear => {
            let mut creator = LinearIndexCreator::new(OPTIMAL_GVCF_INDEX_BIN_SIZE);
            for (feature, position) in &found {
                creator.add_feature(feature, *position);
            }
            let contigs = creator.finalize(file_size, Vec::new()).map_err(|error| {
                Refusal::CouldNotIndexFile {
                    path: source.path.clone(),
                    detail: wrapped(&error),
                }
            })?;
            (LINEAR, Vec::new(), contigs, Vec::new())
        }
        IndexKind::Dynamic => {
            let mut creator = DynamicIndexCreator::new(BalanceApproach::ForSeekTime);
            for (feature, position) in &found {
                creator.add_feature(feature, *position);
            }
            let properties = creator.properties();
            let built =
                creator
                    .finalize(file_size)
                    .map_err(|error| Refusal::CouldNotIndexFile {
                        path: source.path.clone(),
                        detail: wrapped(&error),
                    })?;
            match built {
                BuiltIndex::Linear(contigs) => (LINEAR, properties, contigs, Vec::new()),
                BuiltIndex::IntervalTree(contigs) => {
                    (INTERVAL_TREE, properties, Vec::new(), contigs)
                }
            }
        }
    };

    let index = TribbleIndex {
        index_type,
        version: VERSION,
        indexed_path,
        indexed_file_size: file_size,
        indexed_file_timestamp: source.timestamp,
        indexed_file_md5: String::new(),
        flags: 0,
        properties,
        contigs,
        interval_contigs,
    };
    index.write().map_err(|error| Refusal::CouldNotIndexFile {
        path: source.path.clone(),
        detail: format!("{error:?}"),
    })
}

/// The two lists a caller needs when it holds an explicit output: a block compressed input refuses
/// anything that does not end `.tbi`, and a plain one takes whatever it is given.
pub fn check_output(file_name: &str, output_name: &str, path: &str) -> Result<(), Refusal> {
    if index_kind(file_name) == IndexKind::Tabix && !output_name.ends_with(".tbi") {
        return Err(Refusal::WrongIndexExtension {
            path: path.to_string(),
        });
    }
    Ok(())
}

/// Kept so a caller can name the types the golden's rows carry.
pub fn is_interval_tree(index: &[u8]) -> bool {
    index.len() > 8 && i32::from_le_bytes([index[4], index[5], index[6], index[7]]) == INTERVAL_TREE
}

/// The interval-tree contigs of an index this port built, for a caller that wants to look inside
/// without re-reading the file.
pub fn interval_contigs(index: &TribbleIndex) -> &[IntervalChrIndex] {
    &index.interval_contigs
}
