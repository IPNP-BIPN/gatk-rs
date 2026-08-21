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

use htsjdk_bam::bin::{bin_summary_string, MAX_BINS};
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

/// `Long.toString(value, 16)`, which is lower case and unpadded and prints a negative value with a
/// leading minus rather than as two's complement.
fn hex(value: u64) -> String {
    let signed = value as i64;
    if signed < 0 {
        format!("-{:x}", signed.unsigned_abs())
    } else {
        format!("{signed:x}")
    }
}

/// `Chunk.toString`, which is the two virtual pointers as `block:offset`.
fn chunk_string(pointer: u64) -> String {
    format!("{}:{}", pointer >> 16, pointer & 0xFFFF)
}

/// `BlockCompressedFilePointerUtil.asAddressOffsetString`.
fn address_offset(pointer: u64) -> String {
    format!("{}:{}", pointer >> 16, pointer & 0xFFFF)
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
pub fn render(index: &BamIndex) -> String {
    let mut out = String::new();
    out.push_str(&format!("n_ref={}\n", index.references.len()));
    for (reference, content) in index.references.iter().enumerate() {
        let bins: Vec<_> = content
            .bins
            .iter()
            .filter(|bin| bin.bin_number != MAX_BINS)
            .collect();
        if bins.is_empty() {
            // `writeNullContent`, whose two lines have no space after the `=`.
            out.push_str(&format!("Reference {reference} has n_bin=0\n"));
            out.push_str(&format!("Reference {reference} has n_intv=0\n"));
            continue;
        }
        let counted = bins.len() + usize::from(content.metadata.is_some());
        out.push_str(&format!("Reference {reference} has n_bin= {counted}\n"));
        for bin in &bins {
            out.push_str(&format!(
                "  Ref {reference} bin {} ({}) has n_chunk= {}\n",
                bin.bin_number,
                bin_summary_string(bin.bin_number),
                bin.chunks.len()
            ));
            if bin.chunks.is_empty() {
                out.push('\n');
            }
            for chunk in &bin.chunks {
                out.push_str(&format!(
                    "     Chunk: {}-{} start: {} end: {}\n",
                    chunk_string(chunk.start),
                    chunk_string(chunk.end),
                    hex(chunk.start),
                    hex(chunk.end)
                ));
            }
        }

        // `writeChunkMetaData`, always bin 37450 and always two chunks when metadata is there.
        match content.metadata {
            None => {
                out.push_str(&format!("  Ref {reference} bin 37450 has n_chunk= 0\n"));
                out.push('\n');
            }
            Some(metadata) => {
                out.push_str(&format!("  Ref {reference} bin 37450 has n_chunk= 2\n"));
                out.push_str(&format!(
                    "     Chunk:  start: {} end: {}\n",
                    hex(metadata.first_offset),
                    hex(metadata.last_offset)
                ));
                out.push_str(&format!(
                    "     Chunk:  start: {} end: {}\n",
                    hex(metadata.aligned as u64),
                    hex(metadata.unaligned as u64)
                ));
            }
        }

        if content.linear_index.is_empty() {
            out.push_str(&format!("Reference {reference} has n_intv= 0\n"));
            continue;
        }
        out.push_str(&format!(
            "Reference {reference} has n_intv= {}\n",
            content.linear_index.len()
        ));
        for (window, entry) in content.linear_index.iter().enumerate() {
            if *entry != 0 {
                out.push_str(&format!(
                    "  Ref {reference} ioffset for {window} is {}\n",
                    address_offset(*entry)
                ));
            }
        }
    }
    out.push_str(&format!(
        "No Coordinate Count={}\n",
        index.no_coordinate_records
    ));
    out
}
