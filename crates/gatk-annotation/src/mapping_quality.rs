//! `RMSMappingQuality` and `MappingQualityZero`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.6.2.0).
//!
//! `MQ` and `MQ0`: two counts over the same evidence, disagreeing about which reads exist.
//!
//! # `MQ` is a root mean square over a **long** sum, and the two counts are not the same set
//!
//! ```java
//! long mq = read.getMappingQuality();
//! if (mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE) { squareSum += mq * mq; numReadsUsed++; }
//! ```
//!
//! Reads with a mapping quality of 255 ("unavailable") are dropped from **both** the numerator and
//! the denominator, so they do not depress `MQ`. Reads with a mapping quality of zero are kept in
//! both, so they do. The sum is a `long`, which is exact for any read count a genome can hold, and
//! only the final division reaches floating point.
//!
//! With every read unavailable the divisor is zero, `0 / 0.0` is `NaN`, and the annotation writes
//! the four characters `NaN` into the INFO field rather than declining. That is a real record a
//! downstream parser has to survive.
//!
//! # The raw form is a tuple in a string, and it is not the annotation
//!
//! `RAW_MQandDP` is `String.format("%d,%d", sumOfSquares, depth)`, written by the gVCF path and
//! summed across samples before `finalizeRawData` turns it into `MQ`. The finalised value from the
//! raw tuple and the direct value from the likelihoods take the same square root of the same two
//! numbers, so they agree, and the port shares [`rms_from_tuple`] between them to make that
//! structural rather than a coincidence.
//!
//! # `MQ0` returns zero where every other annotation returns nothing
//!
//! ```java
//! //NOTE: unlike other annotations, this one returns 0 if likelihoods are empty
//! ```
//!
//! An empty matrix is a written `MQ0=0`, and a **null** matrix is an absent key. The two states are
//! distinguishable in the output, which for every other annotation in this crate they are not.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::read::mapping_quality;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::VariantContext;

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use crate::rank_sum::format_decimals;

/// `VCFConstants.RMS_MAPPING_QUALITY_KEY`.
pub const RMS_MAPPING_QUALITY_KEY: &str = "MQ";
/// `VCFConstants.MAPPING_QUALITY_ZERO_KEY`.
pub const MAPPING_QUALITY_ZERO_KEY: &str = "MQ0";
/// `GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY`.
pub const RAW_MAPPING_QUALITY_WITH_DEPTH_KEY: &str = "RAW_MQandDP";
/// `GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_DEPRECATED`.
pub const RAW_RMS_MAPPING_QUALITY_DEPRECATED: &str = "RAW_MQ";

/// `QualityUtils.MAPPING_QUALITY_UNAVAILABLE`.
const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;
/// `RMSMappingQuality.NUM_LIST_ENTRIES`.
const NUM_LIST_ENTRIES: usize = 2;

/// What the raw-data parser refuses, which the reference raises as `UserException.BadInput`.
#[derive(Debug, Clone, PartialEq)]
pub enum RawMqError {
    /// Not two comma-separated values after the brackets were stripped.
    WrongNumberOfValues { found: usize },
    /// A value that is not a `long`.
    Malformed { raw: String },
}

/// `RMSMappingQuality.calculateRawData`: the `(sum of squares, reads used)` tuple.
pub fn raw_mapping_quality_data(likelihoods: &AlleleLikelihoods<BamRecord>) -> (i64, i64) {
    let mut square_sum = 0i64;
    let mut reads_used = 0i64;
    for sample in 0..likelihoods.number_of_samples() {
        let Some(reads) = likelihoods.sample_evidence(sample) else {
            continue;
        };
        for read in reads {
            let mq = mapping_quality(read);
            if mq != MAPPING_QUALITY_UNAVAILABLE {
                let mq = mq as i64;
                square_sum += mq * mq;
                reads_used += 1;
            }
        }
    }
    (square_sum, reads_used)
}

/// `RMSMappingQuality.makeFinalizedAnnotationString`.
///
/// The division is `long / (double)` and the square root is `Math.sqrt`, which IEEE-754 requires
/// to be correctly rounded, so it is the one function in this file that needs no oracle.
pub fn rms_from_tuple(sum_of_squares: i64, num_of_reads: i64) -> String {
    format_decimals(
        jmath::math::sqrt(sum_of_squares as f64 / num_of_reads as f64),
        2,
    )
}

/// `RMSMappingQuality.makeRawAnnotationString`.
pub fn raw_annotation_string(sum_of_squares: i64, num_of_reads: i64) -> String {
    format!("{sum_of_squares},{num_of_reads}")
}

/// `RMSMappingQuality.parseRawDataString`.
///
/// Strips every square bracket anywhere in the string, not only at the ends, then splits on a comma
/// followed by any number of spaces, which is what `AbstractCollection.toString` produces. So
/// `"[10, 2]"` and `"10,2"` parse the same, and `"1[0],2"` parses as ten.
pub fn parse_raw_data_string(raw: &str) -> Result<(i64, i64), RawMqError> {
    let stripped: String = raw
        .trim()
        .chars()
        .filter(|c| *c != '[' && *c != ']')
        .collect();
    let parsed = split_on_comma_spaces(&stripped);
    if parsed.len() != NUM_LIST_ENTRIES {
        return Err(RawMqError::WrongNumberOfValues {
            found: parsed.len(),
        });
    }
    let square_sum = parsed[0]
        .parse::<i64>()
        .map_err(|_| RawMqError::Malformed {
            raw: raw.to_string(),
        })?;
    let total_dp = parsed[1]
        .parse::<i64>()
        .map_err(|_| RawMqError::Malformed {
            raw: raw.to_string(),
        })?;
    Ok((square_sum, total_dp))
}

/// `String.split(", *")`, which drops trailing empty tokens but keeps leading and interior ones.
fn split_on_comma_spaces(text: &str) -> Vec<String> {
    let mut parts: Vec<String> = Vec::new();
    let bytes = text.as_bytes();
    let mut start = 0usize;
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b',' {
            parts.push(text[start..i].to_string());
            i += 1;
            while i < bytes.len() && bytes[i] == b' ' {
                i += 1;
            }
            start = i;
        } else {
            i += 1;
        }
    }
    parts.push(text[start..].to_string());
    while parts.len() > 1 && parts.last().is_some_and(|part| part.is_empty()) {
        parts.pop();
    }
    parts
}

/// `RMSMappingQuality`: `MQ`, computed straight from the likelihoods.
pub struct RmsMappingQuality;

impl InfoFieldAnnotation for RmsMappingQuality {
    fn key_names(&self) -> Vec<&'static str> {
        vec![RMS_MAPPING_QUALITY_KEY, RAW_MAPPING_QUALITY_WITH_DEPTH_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        _vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        // `evidenceCount() < 1`, not "no informative reads": every read in the matrix counts,
        // whatever its likelihoods say.
        let Some(likelihoods) = likelihoods else {
            return Vec::new();
        };
        if likelihoods.evidence_count() < 1 {
            return Vec::new();
        }
        let (square_sum, reads_used) = raw_mapping_quality_data(likelihoods);
        vec![(
            RMS_MAPPING_QUALITY_KEY.to_string(),
            AnnotationValue::Str(rms_from_tuple(square_sum, reads_used)),
        )]
    }
}

impl RmsMappingQuality {
    /// `RMSMappingQuality.annotateRawData`: the gVCF form, `RAW_MQandDP`.
    ///
    /// Guarded on `evidenceCount() == 0` rather than `< 1`, which is the same test written twice
    /// in the same class.
    pub fn annotate_raw_data(
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        let Some(likelihoods) = likelihoods else {
            return Vec::new();
        };
        if likelihoods.evidence_count() == 0 {
            return Vec::new();
        }
        let (square_sum, reads_used) = raw_mapping_quality_data(likelihoods);
        vec![(
            RAW_MAPPING_QUALITY_WITH_DEPTH_KEY.to_string(),
            AnnotationValue::Str(raw_annotation_string(square_sum, reads_used)),
        )]
    }

    /// `RMSMappingQuality.combineRawData`: two tuples added componentwise.
    pub fn combine_raw_data(tuples: &[(i64, i64)]) -> Option<(i64, i64)> {
        let mut combined: Option<(i64, i64)> = None;
        for (square_sum, depth) in tuples {
            combined = Some(match combined {
                None => (*square_sum, *depth),
                Some((s, d)) => (s + square_sum, d + depth),
            });
        }
        combined
    }

    /// `RMSMappingQuality.finalizeRawData` on the modern key alone.
    ///
    /// The deprecated `RAW_MQ` path is refused rather than ported: reaching it requires the
    /// `allow-old-rms-mapping-quality-annotation-data` argument, and without it the reference
    /// throws.
    pub fn finalize_raw_data(raw: &str) -> Result<String, RawMqError> {
        let (square_sum, depth) = parse_raw_data_string(raw)?;
        Ok(rms_from_tuple(square_sum, depth))
    }
}

/// `MappingQualityZero`: `MQ0`, a count of reads whose mapping quality is exactly zero.
pub struct MappingQualityZero;

impl InfoFieldAnnotation for MappingQualityZero {
    fn key_names(&self) -> Vec<&'static str> {
        vec![MAPPING_QUALITY_ZERO_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        let Some(likelihoods) = likelihoods else {
            return Vec::new();
        };
        if !vc.is_variant() {
            return Vec::new();
        }
        let mut count = 0i64;
        for sample in 0..likelihoods.number_of_samples() {
            let Some(reads) = likelihoods.sample_evidence(sample) else {
                continue;
            };
            count += reads
                .iter()
                .filter(|read| mapping_quality(read) == 0)
                .count() as i64;
        }
        vec![(
            MAPPING_QUALITY_ZERO_KEY.to_string(),
            // `String.format("%d", mq0)`: a String, not a Long, so the encoder does not reformat
            // it.
            AnnotationValue::Str(count.to_string()),
        )]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_empty_denominator_is_written_out_as_nan() {
        assert_eq!(rms_from_tuple(0, 0), "NaN");
    }

    #[test]
    fn the_raw_tuple_parses_with_or_without_its_brackets() {
        assert_eq!(parse_raw_data_string("[10, 2]"), Ok((10, 2)));
        assert_eq!(parse_raw_data_string("10,2"), Ok((10, 2)));
        assert_eq!(
            parse_raw_data_string("10"),
            Err(RawMqError::WrongNumberOfValues { found: 1 })
        );
    }

    #[test]
    fn the_root_mean_square_is_not_the_mean() {
        // Two reads at 60 and one at 0: the mean is 40, the root mean square is about 49.
        assert_eq!(rms_from_tuple(60 * 60 * 2, 3), "48.99");
    }
}
