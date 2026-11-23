//! `VCFCodec.decodeLoc` and `BEDCodec.decodeLoc`: where a feature is, and nothing else about it.
//!
//! Tribble reads a feature file twice over: once for the whole record, and once for its LOCUS
//! alone. The second is what an index is built from and what a walker's traversal is decided by,
//! and it is much less than the first: a contig, a start and a stop, taken off the line without
//! parsing an allele or a genotype.
//!
//! # The stop is not the position
//!
//! For a VCF it is `END` where the INFO field carries one, and `POS + len(REF) - 1` otherwise.
//! That is the span a query matches, which is why `-L chr1:605-606` reaches a record whose
//! position it does not contain. For a BED it is the interval's own end, and the START is
//! one-based where the file's is zero-based.
//!
//! # A line that cannot be read is skipped, not refused
//!
//! A header, a blank line, a `track` line, a row with too few columns and a position that is not a
//! number all produce no feature and no error. The codec's own `decodeLoc` returns null for them
//! and Tribble asks for the next line.
//!
//! # Where this belongs
//!
//! In htsjdk, and it is here. `VCFCodec` and `BEDCodec` are `htsjdk.tribble`'s, and htsjdk-rs's
//! own list of classes a ported call site reaches would have this as a row if anything had named
//! it. Two tools reach it now -- `IndexFeatureFile`, which indexes the features, and
//! `CountVariants`, which counts them -- which is the shape that row exists to catch.
//!
//! Ported from `htsjdk.tribble.readers.LineIteratorImpl` as driven by
//! `htsjdk.variant.vcf.VCFCodec.decodeLoc` and `htsjdk.tribble.bed.BEDCodec.decodeLoc`.

use htsjdk_tribble::index_write::Feature;

/// The codecs this port knows, which is what decides whether a file can be indexed at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Codec {
    Vcf,
    Bed,
}

/// `FeatureManager.getCodecForFile`, for the two formats this port reads.
///
/// `VCFCodec.canDecode` reads the FILE, not its name: the first bytes, raw or through a gzip
/// stream, must be `##fileformat=VCFv4`, so `x.blockgz.gz` is a VCF and an empty `x.vcf` is not.
/// `BEDCodec.canDecode` reads the name alone, case-insensitively, past a block compressed
/// extension. A path that cannot be opened falls back to its name, which leaves the caller's own
/// "could not read" refusal where it was.
pub fn codec_for(file_name: &str) -> Option<Codec> {
    match std::fs::File::open(file_name) {
        Ok(file) => {
            if is_vcf_stream(file) {
                Some(Codec::Vcf)
            } else {
                is_bed_name(file_name).then_some(Codec::Bed)
            }
        }
        Err(_) => codec_for_name(file_name),
    }
}

/// `isVCFStream` over the raw bytes, then over a gzip stream (which also reads BGZF).
fn is_vcf_stream(mut file: std::fs::File) -> bool {
    use std::io::{Read, Seek};
    const MAGIC: &[u8] = b"##fileformat=VCFv4";
    let starts = |reader: &mut dyn Read| {
        let mut buffer = [0u8; MAGIC.len()];
        reader.read_exact(&mut buffer).is_ok() && buffer == MAGIC
    };
    starts(&mut file)
        || (file.rewind().is_ok() && starts(&mut flate2::read::MultiGzDecoder::new(file)))
}

/// `BEDCodec.canDecode`: `IOUtil.hasBlockCompressedExtension` stripped, then `.bed` ignoring case.
fn is_bed_name(file_name: &str) -> bool {
    let lower = file_name.to_lowercase();
    let stripped = [".gz", ".gzip", ".bgz", ".bgzf"]
        .iter()
        .find_map(|extension| lower.strip_suffix(extension))
        .unwrap_or(&lower);
    stripped.ends_with(".bed")
}

/// The codec a NAME suggests, for a path that is not (yet) on disk.
pub fn codec_for_name(file_name: &str) -> Option<Codec> {
    let stripped = file_name
        .strip_suffix(".gz")
        .or_else(|| file_name.strip_suffix(".bgz"))
        .unwrap_or(file_name);
    if stripped.ends_with(".vcf") {
        Some(Codec::Vcf)
    } else if is_bed_name(file_name) {
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
