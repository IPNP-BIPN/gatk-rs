//! The Feature-file branch of `-L`, ported from `IntervalUtils.featureFileToIntervals` and the
//! codec resolution `FeatureManager.isFeatureFile` performs (GATK 4.6.2.0).
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
//! # One of the three codecs opens the file
//!
//! `BEDCodec` and `IntervalListCodec` answer `canDecode` on the path; `VCFCodec` reads the first
//! eighteen bytes and looks for `##fileformat=VCFv4`. So the extension decides for two of the
//! three and the contents decide for the other, and a `.list` holding a VCF body is a Feature file
//! while a `.list` holding a BED body is not. Both rows are in the golden.
//!
//! # An interval list validates against two dictionaries, not one
//!
//! `IntervalListCodec.readActualHeader` reads the file's **own** `@SQ` lines and decodes every
//! record against those, and only then does `createGenomeLoc` check the result against the
//! **reference** dictionary. The two are different checks with different outcomes: a contig the
//! file does not declare is dropped and the file still loads, while a contig the file declares and
//! the reference does not kills the argument. A port with one dictionary cannot produce both.

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
            // header does not hold is refused here and not silently kept.
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

/// Which of the codec's refusals `featureFileToIntervals` converts and which it lets through.
///
/// ```java
/// catch (final IllegalArgumentException e) { throw new UserException.MalformedFile(...); }
/// ```
///
/// Only `IllegalArgumentException` is caught. A `TribbleException`, which is what a wrong field
/// count raises, is neither an `IllegalArgumentException` nor a `UserException`, so it leaves the
/// engine as itself: the same malformed file produces two different exception classes depending on
/// which of its lines is wrong.
fn codec_failure(error: interval_list::IntervalListError) -> IntervalArgumentError {
    match error {
        interval_list::IntervalListError::NoDictionary
        | interval_list::IntervalListError::FieldCount { .. } => {
            IntervalArgumentError::FeatureCodecRefused(error.message())
        }
        other => IntervalArgumentError::FeatureFileMalformed(other.message()),
    }
}

/// Which codec claims a path, asked in the reference's own order.
///
/// `FeatureManager` asks every registered codec and takes the first that answers, and the two
/// path-based codecs are asked without opening anything. The VCF one needs the bytes, so it is
/// asked last here, which changes nothing: a `.bed` path and a VCF body cannot both be true of a
/// file the other two claim, since neither of them looks inside.
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
        let header = htsjdk_vcf::header::VcfHeader {
            lines: Vec::new(),
            samples: frame.samples,
        };
        let mut intervals = Vec::new();
        for (number, line) in text.lines().enumerate() {
            match htsjdk_vcf::record_parse::decode_line(line, &header, number) {
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
        // the first to answer wins, which is safe here because their extensions are disjoint.
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
