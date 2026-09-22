//! `FeatureOutputCodecFinder` and the `canDecode` of the ten codecs it knows, ported from
//! `org.broadinstitute.hellbender.utils.codecs` (GATK 4.6.2.0).
//!
//! The SV evidence tools choose what they write by the output's NAME, and so does the engine when it
//! picks a codec to read one: five text codecs and five binary (`.bci`) ones, each answering for
//! one suffix and nothing else.
//!
//! Two details are the codecs' rather than the finder's. A text codec lower-cases the name and
//! strips ONE block-compressed extension (`IOUtil.hasBlockCompressedExtension`: `.gz`, `.gzip`,
//! `.bgz` or `.bgzf`) before it tests the suffix, so `counts.RD.TXT.BGZ` is depth evidence. A binary
//! codec lower-cases and strips nothing, so `counts.rd.bci.gz` is no codec's at all.

/// How a codec lays the features out on disk.
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

/// `FORMAT_SUFFIX` of each text codec, in `FeatureOutputCodecFinder`'s own order.
const TEXT: [(&str, &str); 5] = [
    (".baf.txt", "BafEvidence"),
    (".rd.txt", "DepthEvidence"),
    (".pe.txt", "DiscordantPairEvidence"),
    (".sd.txt", "SiteDepth"),
    (".sr.txt", "SplitReadEvidence"),
];

/// The `*_BCI_FILE_EXTENSION` of each binary codec, in the same order.
const BCI: [(&str, &str); 5] = [
    (".baf.bci", "BafEvidence"),
    (".rd.bci", "DepthEvidence"),
    (".pe.bci", "DiscordantPairEvidence"),
    (".sd.bci", "SiteDepth"),
    (".sr.bci", "SplitReadEvidence"),
];

/// `FileExtensions.BLOCK_COMPRESSED`.
const BLOCK_COMPRESSED: [&str; 4] = [".gz", ".gzip", ".bgz", ".bgzf"];

/// `IOUtil.hasBlockCompressedExtension`.
pub fn has_block_compressed_extension(name: &str) -> bool {
    let lower = name.to_lowercase();
    BLOCK_COMPRESSED
        .iter()
        .any(|extension| lower.ends_with(extension))
}

/// `FeatureOutputCodecFinder.find`, short of its refusal: `None` is the
/// `No feature output codec found for ...` the caller raises. No two codecs claim one name, so the
/// finder's `Found multiple output codecs` cannot be reached.
pub fn find(name: &str) -> Option<Codec> {
    let lower = name.to_lowercase();
    let block_compressed = has_block_compressed_extension(&lower);
    let stripped = if block_compressed {
        &lower[..lower.rfind('.').unwrap_or(lower.len())]
    } else {
        &lower[..]
    };
    if let Some((_, feature_type)) = TEXT.iter().find(|(suffix, _)| stripped.ends_with(suffix)) {
        return Some(Codec {
            feature_type,
            encoding: Encoding::Text { block_compressed },
        });
    }
    BCI.iter()
        .find(|(suffix, _)| lower.ends_with(suffix))
        .map(|(_, feature_type)| Codec {
            feature_type,
            encoding: Encoding::Bci,
        })
}

/// `find`'s refusal.
pub fn no_output_codec(name: &str) -> String {
    format!("No feature output codec found for {name}")
}
