//! The Feature-file branch of `-L`, ported from `IntervalUtils.featureFileToIntervals` and the
//! codec resolution `FeatureManager.isFeatureFile` performs (GATK 4.7.0.0; `IntervalUtils.java`
//! and `FeatureManager.java` are unchanged since 4.6.2.0; the codecs are htsjdk 5.0.0's).
//!
//! This fills the seam [`crate::interval_args::NoFeatureSources`] was standing in for. The
//! `interval-file` suite has measured the answer since it was written; what was missing was a
//! codec to produce it.
//!
//! # A Feature file is recognised by codec, not by extension
//!
//! `FeatureManager` asks every registered codec's `canDecode`, and a codec answers on the path.
//! `BEDCodec.canDecode` strips one block-compressed extension and then tests for `.bed`, so
//! `regions.BED.gz` is a Feature file and `regions.bed.gz.gz` is not. `IntervalListCodec.canDecode`
//! strips nothing, so `.interval_list.bgz` is not one either. Nothing looks inside the file, which
//! is why a `.list` whose contents are BED is **not** a Feature file and reaches the interval
//! reader instead, where it fails as a malformed locus. That row is in the golden.
//!
//! ```java
//! if (IOUtil.hasBlockCompressedExtension(path)) {
//!     toDecode = path.substring(0, path.lastIndexOf("."));
//! } else {
//!     toDecode = path;
//! }
//! return toDecode.toLowerCase().endsWith(FileExtensions.BED);
//! ```
//!
//! `IntervalListCodec.canDecode` is `s.endsWith(FileExtensions.INTERVAL_LIST)` or
//! `s.endsWith(FileExtensions.COMPRESSED_INTERVAL_LIST)`. The BED test lower-cases the path and
//! the interval-list test does not, so `regions.BED` is a Feature file and
//! `regions.INTERVAL_LIST` is not.
//!
//! # The coordinates come from the codec, not from a conversion here
//!
//! ```text
//! chr1  0  10        the BED line
//! chr1:1-10          the interval -L produces
//! ```
//!
//! The +1 is `BEDCodec`'s `StartOffset.ONE`, applied while decoding. This module does not shift
//! anything: `createGenomeLoc(feature, true)` takes the feature's own start and end. A port that
//! converted here as well would be off by one, and a port that converted here *instead* would
//! agree until someone constructed the codec at `StartOffset.ZERO`. `IntervalListCodec` needs no
//! shift at all: its file is 1-based already.
//!
//! ```java
//! // The BED format uses a first-base-is-zero convention,  Tribble features use 1 => add 1.
//! int start = Integer.parseInt(tokens[1]) + startOffsetValue;
//!
//! int end = start;
//! if (tokenCount > 2) {
//!     end = Integer.parseInt(tokens[2]);
//! }
//! ```
//!
//! Only the start moves: BED's end is exclusive and 0-based, which is the same number as a
//! 1-based inclusive end. A two-column line has no end at all and gets the shifted start, so
//! `chr1  9` is the one-base interval `chr1:10-10`.
//!
//! # One of the three codecs opens the file
//!
//! `BEDCodec` and `IntervalListCodec` answer `canDecode` on the path; `VCFCodec` reads the first
//! eighteen bytes and looks for `##fileformat=VCFv4`. So the extension decides for two of the
//! three and the contents decide for the other, and a `.list` holding a VCF body is a Feature file
//! while a `.list` holding a BED body is not. Both rows are in the golden.
//!
//! `AbstractVCFCodec.canDecodeFile` tries the bytes as they are, then through a
//! `GZIPInputStream`, then through a `BlockCompressedInputStream`, and any `IOException` answers
//! `false`, so a gzipped VCF is claimed too and a file that cannot be opened is simply not one.
//!
//! # An interval list validates against two dictionaries, not one
//!
//! `IntervalListCodec.readActualHeader` reads the file's **own** `@SQ` lines and decodes every
//! record against those, and only then does `createGenomeLoc` check the result against the
//! **reference** dictionary. The two are different checks with different outcomes: a contig the
//! file does not declare is dropped and the file still loads, while a contig the file declares and
//! the reference does not kills the argument. A port with one dictionary cannot produce both.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, E>` is "a `T` or an error `E`", how a reference exception is modelled, and `?`
//!   after a call returns its error to the caller at once;
//! - a decoder answers `Ok(Some(record))`, `Ok(None)` (Java's `null`: a line with no record) or
//!   `Err(error)`, and a `match` must name all three;
//! - `Option<Result<..>>` from [`FeatureIntervals::intervals_from_feature_file`] is three
//!   answers: `None` "not a Feature file", `Some(Ok(..))` its intervals, `Some(Err(..))` a
//!   Feature file that failed;
//! - `impl FeatureIntervals for RegisteredCodecs` is Java's `implements`.

use std::path::Path;

use htsjdk_bam::header::SamHeader;
use htsjdk_tribble::bed::{self, StartOffset};
use htsjdk_tribble::{interval_list, vcf};

use crate::interval::SimpleInterval;
use crate::interval_args::{FeatureIntervals, IntervalArgumentError};

/// The codecs this port registers, which is the subset of GATK's that reads intervals.
///
/// BED, interval list and VCF. GATK registers more (BCF, gVCF variants, table codecs), and none of
/// them produce intervals for `-L` in a way this suite reaches.
pub struct RegisteredCodecs;

impl RegisteredCodecs {
    /// `IntervalUtils.featureFileToIntervals` over a BED body: one interval per feature, in file
    /// order.
    ///
    /// A line the codec answers `null` for (blank, header, or fewer than two fields) contributes
    /// nothing, because the reference's iterator skips nulls rather than stopping.
    ///
    /// ```java
    /// for ( final Feature feature : dataSource ) {
    ///     featureIntervals.add(parser.createGenomeLoc(feature, true));
    /// }
    /// ```
    pub fn bed_intervals(
        &self,
        text: &str,
        header: &SamHeader,
        path_hint: &str,
    ) -> Result<Vec<SimpleInterval>, IntervalArgumentError> {
        let mut intervals = Vec::new();
        for line in text.lines() {
            let feature = match bed::decode(line, StartOffset::One) {
                Ok(Some(feature)) => feature,
                Ok(None) => continue,
                // A malformed BED line reaches the caller as a file that could not be read, which
                // is the shape the reference's Tribble layer wraps a codec failure in. The exact
                // exception is not measured yet: the golden has no malformed-BED case, so this is
                // named as the nearest existing refusal rather than invented.
                Err(_) => {
                    return Err(IntervalArgumentError::IntervalFileMissing(
                        path_hint.to_string(),
                    ))
                }
            };
            // `createGenomeLoc(feature, true)` validates against the dictionary, so a contig the
            // header does not hold is refused here and not silently kept. The port reaches that
            // check by rendering `contig:start-end` and parsing it as an `-L` locus.
            intervals.push(crate::interval::parse_interval(
                &format!("{}:{}-{}", feature.contig, feature.start, feature.end),
                header,
            )?);
        }
        Ok(intervals)
    }

    /// The same over an `.interval_list`, whose codec carries a dictionary of its own.
    ///
    /// The dictionary is the file's `@SQ` lines, read by `readActualHeader` before the first
    /// record, so a file with no `@SQ` at all decodes every record to `null` and loads as no
    /// intervals rather than as an error.
    pub fn interval_list_intervals(
        &self,
        text: &str,
        reference: &SamHeader,
    ) -> Result<Vec<SimpleInterval>, IntervalArgumentError> {
        let file_dictionary = htsjdk_bam::reader::parse_header_text(text);
        let mut intervals = Vec::new();
        for line in text.lines() {
            let record = match interval_list::decode(line, Some(&file_dictionary)) {
                Ok(Some(record)) => record,
                Ok(None) => continue,
                Err(error) => return Err(codec_failure(error)),
            };
            intervals.push(crate::interval::parse_interval(
                &format!("{}:{}-{}", record.contig, record.start, record.end),
                reference,
            )?);
        }
        Ok(intervals)
    }
}

/// Which of the codec's refusals `parseIntervalArguments` converts around its call to
/// `featureFileToIntervals`, and which it lets through.
///
/// ```java
/// } catch (final IllegalArgumentException e){
/// ```
///
/// and the handler throws `UserException.MalformedFile` with the message
/// `"Failure while loading intervals from file."`. Only `IllegalArgumentException` is caught. A
/// `TribbleException`, which is what a wrong field count raises, is neither an
/// `IllegalArgumentException` nor a `UserException`, so it leaves the engine as itself: the same
/// malformed file produces two different exception classes depending on which of its lines is
/// wrong.
fn codec_failure(error: interval_list::IntervalListError) -> IntervalArgumentError {
    match error {
        interval_list::IntervalListError::NoDictionary
        | interval_list::IntervalListError::FieldCount { .. } => {
            IntervalArgumentError::FeatureCodecRefused(error.message())
        }
        other => IntervalArgumentError::FeatureFileMalformed(other.message()),
    }
}

/// Which codec claims a path.
///
/// `FeatureManager.getCandidateCodecsForFile` asks EVERY registered codec and collects all that
/// answer; it does not stop at the first:
///
/// ```java
/// if ( codec.canDecode(featureFile.toAbsolutePath().toUri().toString()) ) {
///     candidateCodecs.add(codec);
/// }
/// ```
///
/// `isFeatureFile` needs only a non-empty list, but opening the source goes through
/// `getCodecForFile`, which refuses more than one candidate with a `GATKException` ("Multiple
/// codecs found able to decode file"). The port asks the two path-based codecs first and the VCF
/// one only if neither claims the path. Their two extensions are disjoint, but the VCF codec
/// reads the bytes whatever the name, so a `.bed` or `.interval_list` file whose body starts
/// `##fileformat=VCFv4` is claimed by two codecs in the reference and refused, where the port
/// takes the path codec. No golden row has such a file.
enum Claimed {
    Bed,
    IntervalList,
    Vcf,
}

impl RegisteredCodecs {
    /// `IntervalUtils.featureFileToIntervals` over a VCF: one interval per record, `getStart()` to
    /// `getEnd()`.
    ///
    /// `getEnd()` is the `END` INFO attribute when the record carries one, and the reference
    /// allele's span otherwise, so a symbolic deletion contributes the interval it declares rather
    /// than one base. The record decoder is htsjdk's and already oracle-backed; this only turns
    /// each decoded record into a locus.
    pub fn vcf_intervals(
        &self,
        text: &str,
        reference: &SamHeader,
        path_hint: &str,
    ) -> Result<Vec<SimpleInterval>, IntervalArgumentError> {
        let frame = match htsjdk_vcf::header_parse::read_header_frame(text) {
            Ok(frame) => frame,
            // A file whose magic matched and whose header the codec refuses. The engine never
            // gets a reader, so this is a GATKException rather than the codec's own exception:
            // a malformed header and a malformed record in the same file are two classes.
            Err(error) => return Err(IntervalArgumentError::FeatureSourceFailed(error.0)),
        };
        // The version the file declared, which the codec's text transformer needs: the percent
        // escapes are read back only from 4.3 onwards.
        let version = frame.version;
        let header = htsjdk_vcf::header::VcfHeader {
            lines: Vec::new(),
            samples: frame.samples,
        };
        let mut intervals = Vec::new();
        for (number, line) in text.lines().enumerate() {
            match htsjdk_vcf::record_parse::decode_line(line, &header, number, version) {
                Ok(Some(decoded)) => {
                    let variant = decoded.variant;
                    intervals.push(crate::interval::parse_interval(
                        &format!("{}:{}-{}", variant.contig, variant.start, variant.stop),
                        reference,
                    )?);
                }
                Ok(None) => continue,
                // A record the codec cannot decode is a TribbleException, which
                // `featureFileToIntervals` does not catch, so it leaves the engine as itself.
                Err(error) => {
                    // Marks the parameter as used; this message does not name the file.
                    let _ = path_hint;
                    return Err(IntervalArgumentError::FeatureCodecRefused(format!(
                        "{error:?}"
                    )));
                }
            }
        }
        Ok(intervals)
    }
}

impl FeatureIntervals for RegisteredCodecs {
    fn intervals_from_feature_file(
        &self,
        path: &Path,
        header: &SamHeader,
    ) -> Option<Result<Vec<SimpleInterval>, IntervalArgumentError>> {
        // `canDecode` is asked before the file is opened, exactly as `FeatureManager` asks it, so
        // a path no codec claims never becomes a read attempt. The codecs are asked in turn and
        // the first to answer wins; see [`Claimed`] for where that differs from the reference,
        // which collects every codec that answers.
        //
        // The reference hands `canDecode` the absolute file URI; the port hands it the path as
        // given (`to_string_lossy` replaces any bytes that are not UTF-8). The suffix tests read
        // the same on both forms unless URI escaping changes the end of the name.
        let name = path.to_string_lossy().to_string();
        let claimed = if bed::can_decode(&name) {
            Claimed::Bed
        } else if interval_list::can_decode(&name) {
            Claimed::IntervalList
        } else {
            // The VCF codec reads the file to answer, and a path it cannot open answers `false`
            // through a caught IOException rather than raising, so an unreadable file is simply
            // not a Feature file and falls through to the branch that reports it.
            let Ok(bytes) = std::fs::read(path) else {
                return None;
            };
            if !vcf::can_decode(&bytes) {
                return None;
            }
            Claimed::Vcf
        };
        // The file is read as UTF-8 text with no decompression, so a block-compressed file that a
        // codec claimed (`regions.bed.gz`, a gzipped VCF) fails here and is reported as
        // `IntervalFileMissing`. No golden row has a compressed Feature file.
        let text = match std::fs::read_to_string(path) {
            Ok(text) => text,
            Err(_) => return Some(Err(IntervalArgumentError::IntervalFileMissing(name))),
        };
        Some(match claimed {
            Claimed::Bed => self.bed_intervals(&text, header, &name),
            Claimed::IntervalList => self.interval_list_intervals(&text, header),
            Claimed::Vcf => self.vcf_intervals(&text, header, &name),
        })
    }
}
