//! `FeatureOutputCodecFinder` and the `canDecode` of the ten codecs it knows, ported from
//! `org.broadinstitute.hellbender.utils.codecs` (GATK 4.6.2.0).
//!
//! The SV evidence tools choose what they write by the output's NAME, and so does the engine when
//! it picks a codec to read one: five text codecs and five binary (`.bci`) ones, each answering
//! for one suffix and nothing else.
//!
//! Two details are the codecs' rather than the finder's. A text codec lower-cases the name and
//! strips ONE block-compressed extension (`IOUtil.hasBlockCompressedExtension`: `.gz`, `.gzip`,
//! `.bgz` or `.bgzf`) before it tests the suffix, so `counts.RD.TXT.BGZ` is depth evidence. A
//! binary codec lower-cases and strips nothing, so `counts.rd.bci.gz` is no codec's at all.
//!
//! The finder asks every codec and refuses a name two of them claim. No suffix here is the tail
//! of another, so at most one codec ever answers, and this port asks the five text codecs first
//! and the five binary ones after without changing which one that is.
//!
//! # Reading this file without reading Rust
//!
//! - `&str` is a borrowed, read-only view of a string; `&'static str` is one that lives for the
//!   whole program (a literal), Java's `static final String`.
//! - `Option<Codec>` is "a codec or nothing", Java's nullable return checked by the compiler:
//!   `Some(codec)` is the value, `None` the `null`.
//! - `|x| ...` is a closure, Java's lambda `x -> ...`. `.iter().any(..)` and `.iter().find(..)`
//!   are Java's `stream().anyMatch(..)` and `stream().filter(..).findFirst()`.
//! - `&lower[..n]` is the slice of `lower` before byte `n`, Java's `substring(0, n)`.

/// How a codec lays the features out on disk.
///
/// A Rust `enum` variant may carry its own fields, so `Text` holds whether the name asked for
/// block compression and `Bci` holds nothing. Java would need a small class hierarchy for this.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Encoding {
    /// Tab-separated text, block compressed when the name says so.
    Text { block_compressed: bool },
    /// The binary container, `.bci`.
    Bci,
}

/// The codec a name selects: the feature type it produces and how it is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Codec {
    /// `getFeatureType().getSimpleName()`.
    pub feature_type: &'static str,
    pub encoding: Encoding,
}

/// `FORMAT_SUFFIX` of each text codec, in `FeatureOutputCodecFinder`'s own order, paired with
/// the simple name of the feature type it produces. Each is one line of its codec, for example
/// `DepthEvidenceCodec.java`:
///
/// ```java
/// public static final String FORMAT_SUFFIX = ".rd.txt";
/// ```
const TEXT: [(&str, &str); 5] = [
    (".baf.txt", "BafEvidence"),
    (".rd.txt", "DepthEvidence"),
    (".pe.txt", "DiscordantPairEvidence"),
    (".sd.txt", "SiteDepth"),
    (".sr.txt", "SplitReadEvidence"),
];

/// The `*_BCI_FILE_EXTENSION` of each binary codec, in the same order (`BafEvidenceBCICodec`
/// through `SplitReadEvidenceBCICodec`, the finder's last five entries).
const BCI: [(&str, &str); 5] = [
    (".baf.bci", "BafEvidence"),
    (".rd.bci", "DepthEvidence"),
    (".pe.bci", "DiscordantPairEvidence"),
    (".sd.bci", "SiteDepth"),
    (".sr.bci", "SplitReadEvidence"),
];

/// `FileExtensions.BLOCK_COMPRESSED`.
const BLOCK_COMPRESSED: [&str; 4] = [".gz", ".gzip", ".bgz", ".bgzf"];

/// `IOUtil.hasBlockCompressedExtension(String)` (htsjdk 5.0.0):
///
/// ```java
/// String cleanedPath = stripQueryStringIfPathIsAnHttpUrl(fileName);
/// for (final String extension : FileExtensions.BLOCK_COMPRESSED) {
///     if (cleanedPath.toLowerCase().endsWith(extension)) return true;
/// }
/// return false;
/// ```
///
/// The first line is not ported: a name starting `http://` or `https://` keeps its query string
/// here, so `https://h/x.rd.txt.gz?y` is not block compressed in this port and is in
/// htsjdk.
pub fn has_block_compressed_extension(name: &str) -> bool {
    let lower = name.to_lowercase();
    BLOCK_COMPRESSED
        .iter()
        .any(|extension| lower.ends_with(extension))
}

/// `FeatureOutputCodecFinder.find`, short of its refusal: `None` is the
/// `No feature output codec found for ...` the caller raises. No two codecs claim one name, so the
/// finder's `Found multiple output codecs` cannot be reached.
///
/// Every text codec's `canDecode` is the same four lines over its own suffix, here
/// `DepthEvidenceCodec`'s:
///
/// ```java
/// String toDecode = path.toLowerCase();
/// if ( IOUtil.hasBlockCompressedExtension(toDecode) ) {
///     toDecode = toDecode.substring(0, toDecode.lastIndexOf('.'));
/// }
/// return toDecode.endsWith(FORMAT_SUFFIX);
/// ```
///
/// and every binary codec's is one, here `DepthEvidenceBCICodec`'s:
///
/// ```java
/// return path.toLowerCase().endsWith(RD_BCI_FILE_EXTENSION);
/// ```
pub fn find(name: &str) -> Option<Codec> {
    let lower = name.to_lowercase();
    let block_compressed = has_block_compressed_extension(&lower);
    // `if` is an expression in Rust, so `stripped` takes the value of whichever branch runs.
    // A block-compressed name always has a `.`, so the `unwrap_or` fallback (no dot: keep the
    // whole name) is never taken; it stands in for Java's `lastIndexOf` without a panic.
    let stripped = if block_compressed {
        &lower[..lower.rfind('.').unwrap_or(lower.len())]
    } else {
        &lower[..]
    };
    // `if let Some((_, feature_type)) = ..` runs the block only when `find` found an entry, and
    // names the entry's second half; `_` discards the suffix.
    if let Some((_, feature_type)) = TEXT.iter().find(|(suffix, _)| stripped.ends_with(suffix)) {
        return Some(Codec {
            feature_type,
            encoding: Encoding::Text { block_compressed },
        });
    }
    // The binary suffixes are tested against the unstripped name, as their `canDecode` does.
    // `.map(..)` on an `Option` turns `Some(entry)` into `Some(codec)` and leaves `None` alone.
    BCI.iter()
        .find(|(suffix, _)| lower.ends_with(suffix))
        .map(|(_, feature_type)| Codec {
            feature_type,
            encoding: Encoding::Bci,
        })
}

/// `find`'s refusal, the message of its `UserException`:
///
/// ```java
/// throw new UserException("No feature output codec found for " + outputFileName);
/// ```
pub fn no_output_codec(name: &str) -> String {
    format!("No feature output codec found for {name}")
}
