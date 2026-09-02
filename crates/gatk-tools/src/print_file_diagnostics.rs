//! `PrintFileDiagnostics`, ported from
//! `org.broadinstitute.hellbender.tools.PrintFileDiagnostics` and the analyzers under
//! `org.broadinstitute.hellbender.tools.filediagnostics` (GATK 4.6.2.0).
//!
//! An index printed as text. The BAI branch is `htsjdk.samtools.TextualBAMIndexWriter`, which no
//! other tool in this repository reaches, and it is the only place a `.bai` is written as anything
//! but bytes.
//!
//! # The analyzer is chosen by the name
//!
//! ```java
//! if (inputPath.isCram()) { ... } else if (hasExtension(CRAM_INDEX)) { ... }
//! else if (hasExtension(BAI_INDEX)) { ... } else { throw new RuntimeException(...); }
//! ```
//!
//! Nothing reads the file to decide, and the refusal quotes the RAW argument.
//!
//! # An empty reference is printed with different spacing
//!
//! `writeNullContent` writes `n_bin=0` and `n_intv=0` with no space after the `=`, where every
//! other line writes `n_bin= 4` with one. A port that formatted both the same way would differ
//! from the reference on any file with an unused contig.
//!
//! # The metadata bin is counted with the others and printed apart
//!
//! `n_bin` is the real bins plus one whenever the pseudo-bin is there, and the pseudo-bin is then
//! printed after every real bin, out of numeric order, always claiming two chunks: the first pair
//! is offsets, the second is the aligned and unaligned counts printed in the same hexadecimal.

use htsjdk_bam::index::{read_bai, BaiParseError, BamIndex};

/// Which analyzer the name selects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Analyzer {
    Cram,
    Crai,
    Bai,
}

/// What the run refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiagnosticsError {
    /// An extension no analyzer claims, which quotes the raw argument.
    Unsupported { raw: String },
    /// The `.bai` could not be read.
    Unreadable { detail: String },
}

impl DiagnosticsError {
    pub fn java_class(&self) -> &str {
        match self {
            DiagnosticsError::Unsupported { .. } => "java.lang.RuntimeException",
            DiagnosticsError::Unreadable { .. } => "htsjdk.samtools.SAMException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            DiagnosticsError::Unsupported { raw } => {
                format!("Unsupported diagnostic file type: {raw}")
            }
            DiagnosticsError::Unreadable { detail } => detail.clone(),
        }
    }
}

/// `HTSAnalyzerFactory.getFileAnalyzer`, which reads the name and nothing else.
pub fn analyzer_for(raw: &str) -> Result<Analyzer, DiagnosticsError> {
    if raw.ends_with(".cram") {
        Ok(Analyzer::Cram)
    } else if raw.ends_with(".crai") {
        Ok(Analyzer::Crai)
    } else if raw.ends_with(".bai") {
        Ok(Analyzer::Bai)
    } else {
        Err(DiagnosticsError::Unsupported {
            raw: raw.to_string(),
        })
    }
}

/// `BAIAnalyzer.doAnalysis`: the whole report of one `.bai`.
pub fn bai_report(bai: &[u8]) -> Result<String, DiagnosticsError> {
    let index = read_bai(bai).map_err(|error| DiagnosticsError::Unreadable {
        detail: match error {
            BaiParseError::NotABai => "Not a BAM index file".to_string(),
            BaiParseError::Truncated => "Truncated BAM index file".to_string(),
        },
    })?;
    Ok(render(&index))
}

/// `TextualBAMIndexWriter` over an index already read.
///
/// The renderer moved to `htsjdk-bam::textual_index`: the format is htsjdk's, not this tool's, and
/// it is measured there against the reference over the same `.bai` files the index suite measures
/// as bytes. This tool keeps the analyzer that chooses it and the refusal that names the argument.
pub fn render(index: &BamIndex) -> String {
    htsjdk_bam::textual_index::render(index)
}
