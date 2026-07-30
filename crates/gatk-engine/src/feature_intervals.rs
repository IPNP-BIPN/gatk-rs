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
//! `regions.BED.gz` is a Feature file and `regions.bed.gz.gz` is not. Nothing looks inside the
//! file, which is why a `.list` whose contents are BED is **not** a Feature file and reaches the
//! interval reader instead, where it fails as a malformed locus. That row is in the golden.
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
//! agree until someone constructed the codec at `StartOffset.ZERO`.

use std::path::Path;

use htsjdk_bam::header::SamHeader;
use htsjdk_tribble::bed::{self, StartOffset};

use crate::interval::SimpleInterval;
use crate::interval_args::{FeatureIntervals, IntervalArgumentError};

/// The codecs this port registers, which is the subset of GATK's that reads intervals.
///
/// `.interval_list` is **not** here: htsjdk's `IntervalList` reader is its own thing and lands with
/// its own slice. Until then this source answers `None` for one, which sends it down the
/// interval-file branch exactly as an unregistered codec would, and the suite says so.
pub struct BedFeatureSource;

impl BedFeatureSource {
    /// `IntervalUtils.featureFileToIntervals`: one interval per feature, in file order.
    ///
    /// A line the codec answers `null` for (blank, header, or fewer than two fields) contributes
    /// nothing, because the reference's iterator skips nulls rather than stopping.
    pub fn intervals(
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
}

impl FeatureIntervals for BedFeatureSource {
    fn intervals_from_feature_file(
        &self,
        path: &Path,
        header: &SamHeader,
    ) -> Option<Result<Vec<SimpleInterval>, IntervalArgumentError>> {
        // `canDecode` is asked before the file is opened, exactly as `FeatureManager` asks it, so
        // a path no codec claims never becomes a read attempt.
        if !bed::can_decode(&path.to_string_lossy()) {
            return None;
        }
        let text = match std::fs::read_to_string(path) {
            Ok(text) => text,
            Err(_) => {
                return Some(Err(IntervalArgumentError::IntervalFileMissing(
                    path.to_string_lossy().to_string(),
                )))
            }
        };
        Some(self.intervals(&text, header, &path.to_string_lossy()))
    }
}
