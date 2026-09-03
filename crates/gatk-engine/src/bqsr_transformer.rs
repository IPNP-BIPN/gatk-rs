//! `BQSRReadTransformer`, ported from `org.broadinstitute.hellbender.transformers` (GATK 4.6.2.0).
//!
//! A read in, the same read out with every base quality replaced. Everything else in the
//! recalibration machinery is an input to this one function: [`crate::recal_datum`] holds the
//! counts, [`crate::covariates`] turns a read into keys, [`crate::recalibration_tables`] stores the
//! datums those keys reach, and [`crate::qual_quantizer`] collapses the answer.
//!
//! # The estimate is `y3 + y4 - y2`
//!
//! ```java
//! final double empiricalQualityForReadGroup = readGroupDatum == null ? priorQualityScore : readGroupDatum.getEmpiricalQuality(priorQualityScore);
//! final double posteriorEmpiricalQualityForReportedQuality = qualityScoreDatum == null ? empiricalQualityForReadGroup :
//!         qualityScoreDatum.getEmpiricalQuality(empiricalQualityForReadGroup);
//! for( final RecalDatum specialCovariateDatum : specialCovariateDatums ) {
//!     if (specialCovariateDatum != null) {
//!         deltaSpecialCovariates += specialCovariateDatum.getEmpiricalQuality(posteriorEmpiricalQualityForReportedQuality) - posteriorEmpiricalQualityForReportedQuality;
//!     }
//! }
//! ```
//!
//! The first two covariates chain: the read group's empirical quality is the prior for the reported
//! quality's. The special covariates do not chain; each contributes its own empirical quality minus
//! that posterior. **A null datum contributes nothing at all**, which is not the same as
//! contributing a zero delta.
//!
//! # And it depends on which reads came before
//!
//! `getEmpiricalQuality(prior)` computes on the first call and returns the cached value ever after,
//! whatever prior is asked for, and the datums live in the table and are shared across every read.
//! So the first read to reach a datum fixes its empirical quality for the whole run. The reference's
//! own comment says so:
//!
//! ```text
//! // TODO: the prior is ignored if the empirical quality for the datum is already cached.
//! ```
//!
//! The golden measures it directly: the same datum asked with a prior of 25 and then of 45 answers
//! 25.0 both times, where a fresh datum asked with 45 answers 45.0. This port takes
//! `&mut RecalibrationTables` for that reason: the mutation is real and hiding it would make the
//! second run of the same table give a different answer for no visible cause.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use crate::covariates::{
    CovariateError, PerReadCovariateMatrix, ReadGroupCovariate, StandardCovariateList,
    BASE_QUALITY_COVARIATE_DEFAULT_INDEX, NUM_REQUIRED_COVARITES,
};
use crate::math_utils::{fast_round, qual_to_prob};
use crate::qual_quantizer::{bound_qual, QuantizationInfo, MAX_SAM_QUAL_SCORE, MIN_USABLE_Q_SCORE};
use crate::read_group;
use crate::recal_datum::{EventType, MAX_RECALIBRATED_Q_SCORE};
use crate::recalibration_tables::{RecalibrationTables, SharedDatum};

/// `QualityUtils.MAX_QUAL`, the length of a static quantization mapping.
pub const MAX_QUAL: usize = 254;

/// The reference's `READ_GROUP_MISSING_IN_RECAL_TABLE_CODE`.
const READ_GROUP_MISSING_IN_RECAL_TABLE_CODE: i32 = -1;

/// `ApplyBQSRArgumentCollection`, with the defaults the reference declares.
#[derive(Debug, Clone, PartialEq)]
pub struct ApplyBqsrArguments {
    /// `--preserve-qscores-less-than`. A quality below this is left alone entirely, not even
    /// quantized.
    pub preserve_qscores_less_than: i32,
    /// `--quantize-quals`. Zero means no quantization at all; a positive value re-quantizes the
    /// report's map; a negative one keeps the report's.
    pub quantization_levels: i32,
    /// `--static-quantized-quals`, which is a **separate** mechanism from the one above and does not
    /// consult the quantization info.
    pub static_quantization_quals: Vec<i32>,
    /// `--round-down-quantized`.
    pub round_down: bool,
    /// `--emit-original-quals`.
    pub emit_original_quals: bool,
    /// `--use-original-qualities`.
    pub use_original_base_qualities: bool,
    /// `--global-qscore-prior`. Used only when strictly positive.
    pub global_qscore_prior: f64,
    /// `--allow-missing-read-group`. See [`BqsrError::ReadGroupNotInTable`]: it covers a narrower
    /// case than its name suggests.
    pub allow_missing_read_groups: bool,
}

impl Default for ApplyBqsrArguments {
    fn default() -> Self {
        ApplyBqsrArguments {
            preserve_qscores_less_than: MIN_USABLE_Q_SCORE,
            quantization_levels: 0,
            static_quantization_quals: Vec::new(),
            round_down: false,
            emit_original_quals: false,
            use_original_base_qualities: false,
            global_qscore_prior: -1.0,
            allow_missing_read_groups: false,
        }
    }
}

/// The two ways the transformer stops.
#[derive(Debug, Clone, PartialEq)]
pub enum BqsrError {
    /// `GATKException`: the **covariate** does not know this read group, and
    /// `--allow-missing-read-group` was not set. This is the case the flag is for, and it happens
    /// when the recalibration report was written from a different set of read groups than the BAM
    /// has.
    ReadGroupNotInTable(String),
    /// `GATKException`: the covariate knows the read group but the **table** holds no datum for it.
    /// The flag does not cover this, because it is tested before this lookup.
    NoReadGroupDatum(String),
    /// Anything the covariates refuse while computing the keys.
    Covariate(CovariateError),
}

impl BqsrError {
    /// The Java class the reference throws.
    ///
    /// Both of the transformer's own refusals are `GATKException`s rather than `UserException`s,
    /// which is what the documentation above already said and what the caller did not read: a
    /// recalibration table written from a different BAM leaves exit 3, not 2.
    pub fn java_class(&self) -> &'static str {
        match self {
            BqsrError::ReadGroupNotInTable(_) | BqsrError::NoReadGroupDatum(_) => {
                "org.broadinstitute.hellbender.exceptions.GATKException"
            }
            BqsrError::Covariate(error) => error.java_class(),
        }
    }

    pub fn message(&self) -> String {
        match self {
            // No space after the full stop, which is the reference's concatenation.
            BqsrError::ReadGroupNotInTable(group) => format!(
                "Read group {group} not found in the recalibration table.Set the \
                 allow-missing-read-group command line argument to ignore this error."
            ),
            BqsrError::NoReadGroupDatum(identifier) => {
                format!("readGroupDatum for {identifier} is null")
            }
            BqsrError::Covariate(error) => error.message(),
        }
    }
}

/// `hierarchicalBayesianQualityEstimate`.
///
/// See the module note. The datums are mutated, because reading an empirical quality caches it.
pub fn hierarchical_bayesian_quality_estimate(
    prior_quality_score: f64,
    read_group_datum: Option<&SharedDatum>,
    quality_score_datum: Option<&SharedDatum>,
    special_covariate_datums: &[Option<SharedDatum>],
) -> f64 {
    let empirical_for_read_group = match read_group_datum {
        None => prior_quality_score,
        Some(datum) => datum
            .borrow_mut()
            .empirical_quality_with_prior(prior_quality_score),
    };
    let posterior = match quality_score_datum {
        None => empirical_for_read_group,
        Some(datum) => datum
            .borrow_mut()
            .empirical_quality_with_prior(empirical_for_read_group),
    };

    let mut delta = 0.0;
    for datum in special_covariate_datums.iter().flatten() {
        delta += datum.borrow_mut().empirical_quality_with_prior(posterior) - posterior;
    }
    posterior + delta
}

/// `getBoundedIntegerQual`: `fastRound` and then the `[1, 93]` clamp.
///
/// `fastRound` is `(int)(x + 0.5)`, not `Math.round`, so the double below a half rounds twice and
/// answers 1 where `Math.round` answers 0.
pub fn bounded_integer_qual(recalibrated: f64) -> u8 {
    bound_qual(fast_round(recalibrated), MAX_RECALIBRATED_Q_SCORE)
}

/// `constructStaticQuantizedMapping(staticQuantizedQuals, roundDown)`.
///
/// **It sorts the caller's own list**, which this port reproduces by taking `&mut Vec<i32>`.
///
/// Three things decide the mapping:
///
///  * every quality **below** `MIN_USABLE_Q_SCORE` maps to itself, so the special codes survive;
///  * with one static quality, everything from there up takes it;
///  * otherwise the rounding is in **probability space**, not in Phred space, so the midpoint
///    between two static qualities is not their arithmetic mean. `roundDown` skips the comparison
///    and always takes the lower one.
#[allow(clippy::ptr_arg)]
pub fn construct_static_quantized_mapping(
    // `&mut Vec` and not `&mut [i32]` on purpose: the reference sorts the caller's own list, and a
    // slice would hide that this call has an effect the caller can see.
    static_quantized_quals: &mut Vec<i32>,
    round_down: bool,
) -> Vec<u8> {
    if static_quantized_quals.is_empty() {
        // `createIdentityMatrix(MAX_QUAL)`, whose entries wrap past 127 exactly as a Java byte does.
        return (0..MAX_QUAL).map(|i| i as u8).collect();
    }
    let mut mapping = vec![0u8; MAX_QUAL];
    static_quantized_quals.sort();

    for (i, slot) in mapping
        .iter_mut()
        .enumerate()
        .take(MIN_USABLE_Q_SCORE as usize)
    {
        *slot = i as u8;
    }

    if static_quantized_quals.len() == 1 {
        let only = static_quantized_quals[0] as u8;
        for slot in mapping
            .iter_mut()
            .take(MAX_QUAL)
            .skip(MIN_USABLE_Q_SCORE as usize)
        {
            *slot = only;
        }
        return mapping;
    }

    let mut previous_qual = MIN_USABLE_Q_SCORE;
    let mut previous_prob = qual_to_prob(previous_qual as f64);
    for next_qual in static_quantized_quals.iter().copied() {
        let next_prob = qual_to_prob(next_qual as f64);
        for i in previous_qual..next_qual {
            if let Some(slot) = mapping.get_mut(i as usize) {
                *slot = if round_down {
                    previous_qual as u8
                } else {
                    let this_prob = qual_to_prob(i as f64);
                    if this_prob - previous_prob > next_prob - this_prob {
                        next_qual as u8
                    } else {
                        previous_qual as u8
                    }
                };
            }
        }
        previous_qual = next_qual;
        previous_prob = next_prob;
    }
    for slot in mapping
        .iter_mut()
        .take(MAX_QUAL)
        .skip(previous_qual.max(0) as usize)
    {
        *slot = previous_qual as u8;
    }
    mapping
}

/// `BQSRReadTransformer`.
///
/// It holds the tables mutably, because applying it to a read caches empirical qualities in them.
pub struct BqsrReadTransformer<'a> {
    tables: &'a mut RecalibrationTables,
    covariates: &'a StandardCovariateList,
    header: &'a SamHeader,
    preserve_q_less_than: i32,
    constant_quality_score_prior: f64,
    emit_original_quals: bool,
    use_original_base_qualities: bool,
    total_covariate_count: usize,
    static_quantized_mapping: Option<Vec<u8>>,
    quantized_quals: Vec<u8>,
    allow_missing_read_groups: bool,
}

impl<'a> BqsrReadTransformer<'a> {
    /// The constructor, whose first act is to decide the quantization.
    ///
    /// Zero levels means the identity map, a positive value **different from the report's** means a
    /// re-quantization, and anything else keeps the report's map. The static quantization is a
    /// separate mechanism that does not consult the quantization info at all.
    pub fn new(
        header: &'a SamHeader,
        tables: &'a mut RecalibrationTables,
        quantization_info: &mut QuantizationInfo,
        covariates: &'a StandardCovariateList,
        arguments: &ApplyBqsrArguments,
    ) -> Result<BqsrReadTransformer<'a>, BqsrError> {
        if arguments.quantization_levels == 0 {
            quantization_info.no_quantization();
        } else if arguments.quantization_levels > 0
            && arguments.quantization_levels != quantization_info.quantization_levels
        {
            // The reference lets this failure propagate; the histogram always has 94 bins here, so
            // only a level count of zero could fail, and that took the branch above.
            let _ = quantization_info.quantize_quality_scores(arguments.quantization_levels);
        }

        let static_quantized_mapping = if arguments.static_quantization_quals.is_empty() {
            None
        } else {
            let mut quals = arguments.static_quantization_quals.clone();
            Some(construct_static_quantized_mapping(
                &mut quals,
                arguments.round_down,
            ))
        };

        Ok(BqsrReadTransformer {
            tables,
            covariates,
            header,
            preserve_q_less_than: arguments.preserve_qscores_less_than,
            constant_quality_score_prior: arguments.global_qscore_prior,
            emit_original_quals: arguments.emit_original_quals,
            use_original_base_qualities: arguments.use_original_base_qualities,
            total_covariate_count: covariates.size(),
            static_quantized_mapping,
            quantized_quals: quantization_info.quantized_quals.clone(),
            allow_missing_read_groups: arguments.allow_missing_read_groups,
        })
    }

    /// `apply(read)`: the read with every usable base quality replaced.
    pub fn apply(&mut self, original: &BamRecord) -> Result<BamRecord, BqsrError> {
        let mut read = original.clone();
        if self.use_original_base_qualities {
            // `ReadUtils.resetOriginalBaseQualities`, which puts the OQ tag back into the qualities
            // and leaves the read alone when there is no OQ.
            if let Some(original_qualities) = read_original_qualities(&read) {
                read.base_qualities = original_qualities;
            }
        }

        if self.emit_original_quals && !has_tag(&read, b"OQ") {
            let fastq: String = read
                .base_qualities
                .iter()
                .map(|quality| (quality + 33) as char)
                .collect();
            set_tag(&mut read, b"OQ", &fastq);
        }

        let mut matrix = PerReadCovariateMatrix::new(read.read_bases.len(), self.covariates.size());
        self.covariates
            .populate_per_read_covariate_matrix(&read, self.header, &mut matrix, false)
            .map_err(BqsrError::Covariate)?;

        // The reference clears the indel qualities here. They are tags, and this port drops them
        // the same way.
        clear_tag(&mut read, b"BI");
        clear_tag(&mut read, b"BD");

        let identifier = read_group::resolve(&read, self.header)
            .map(ReadGroupCovariate::read_group_identifier)
            .ok_or(BqsrError::Covariate(CovariateError::NoReadGroupInHeader))?;
        let rg_key = self.covariates.read_group.key_from_value(&identifier);

        let mut recalibrated = read.base_qualities.clone();
        let pre_update = read.base_qualities.clone();

        if rg_key == READ_GROUP_MISSING_IN_RECAL_TABLE_CODE {
            if !self.allow_missing_read_groups {
                // The reference names the read's own `RG` attribute here, not the identifier.
                return Err(BqsrError::ReadGroupNotInTable(
                    read_group_attribute(&read).unwrap_or_default(),
                ));
            }
            // Quantized, but not recalibrated: the table has nothing to recalibrate it with.
            for (index, quality) in recalibrated.iter_mut().enumerate() {
                *quality = match &self.static_quantized_mapping {
                    Some(mapping) => mapping[pre_update[index] as usize],
                    None => self.quantized_quals[pre_update[index] as usize],
                };
            }
            read.base_qualities = recalibrated;
            return Ok(read);
        }

        let substitution = EventType::BaseSubstitution.ordinal() as i32;
        let read_group_datum = self
            .tables
            .read_group_table()
            .get2_keys(rg_key, substitution)
            .ok()
            .flatten()
            .ok_or_else(|| BqsrError::NoReadGroupDatum(identifier.clone()))?;

        let mut specials: Vec<Option<SharedDatum>> =
            vec![None; self.total_covariate_count - NUM_REQUIRED_COVARITES];

        // Indexed rather than iterated because the loop both reads and writes the slot, and the
        // reference's own `continue` leaves it untouched.
        #[allow(clippy::needless_range_loop)]
        for offset in 0..recalibrated.len() {
            if (recalibrated[offset] as i32) < self.preserve_q_less_than {
                continue;
            }
            // The array is cleared and reused, so a covariate that is skipped this time does not
            // keep the datum it had last time.
            specials.iter_mut().for_each(|slot| *slot = None);

            let keys = matrix.covariates_at_offset(offset, EventType::BaseSubstitution);
            let reported = keys[BASE_QUALITY_COVARIATE_DEFAULT_INDEX];

            let quality_score_datum = self
                .tables
                .quality_score_table()
                .get3_keys(rg_key, reported, substitution)
                .ok()
                .flatten();

            for j in NUM_REQUIRED_COVARITES..self.total_covariate_count {
                // A key of -1 is skipped and its datum stays absent, which is why the first bases
                // of a read, which have no context, are recalibrated from fewer covariates.
                if keys[j] >= 0 {
                    specials[j - NUM_REQUIRED_COVARITES] = self.tables.all_tables[j]
                        .get4_keys(rg_key, reported, keys[j], substitution)
                        .ok()
                        .flatten();
                }
            }

            let prior = if self.constant_quality_score_prior > 0.0 {
                self.constant_quality_score_prior
            } else {
                read_group_datum.borrow().reported_quality()
            };
            let raw = hierarchical_bayesian_quality_estimate(
                prior,
                Some(&read_group_datum),
                quality_score_datum.as_ref(),
                &specials,
            );
            let quantized = self.quantized_quals[bounded_integer_qual(raw) as usize];
            // Quantized twice when the static mapping is on, which the reference's own TODO calls
            // out and this port keeps.
            recalibrated[offset] = match &self.static_quantized_mapping {
                None => quantized,
                Some(mapping) => mapping[quantized as usize],
            };
        }

        read.base_qualities = recalibrated;
        Ok(read)
    }
}

fn has_tag(read: &BamRecord, name: &[u8; 2]) -> bool {
    read.tags.iter().any(|(tag, _)| tag.name() == *name)
}

fn read_group_attribute(read: &BamRecord) -> Option<String> {
    read.tags.iter().find_map(|(tag, value)| {
        (tag.name() == *b"RG").then(|| match value {
            htsjdk_bam::tag::TagValue::Str(text) => text.clone(),
            other => format!("{other:?}"),
        })
    })
}

/// `ReadUtils.resetOriginalBaseQualities`: the `OQ` tag decoded back into qualities.
fn read_original_qualities(read: &BamRecord) -> Option<Vec<u8>> {
    read.tags.iter().find_map(|(tag, value)| {
        (tag.name() == *b"OQ").then(|| match value {
            htsjdk_bam::tag::TagValue::Str(text) => {
                text.bytes().map(|character| character - 33).collect()
            }
            _ => read.base_qualities.clone(),
        })
    })
}

fn set_tag(read: &mut BamRecord, name: &[u8; 2], value: &str) {
    read.tags.insert(
        htsjdk_bam::tag::Tag::new(name),
        htsjdk_bam::tag::TagValue::Str(value.to_string()),
    );
}

fn clear_tag(read: &mut BamRecord, name: &[u8; 2]) {
    read.tags.remove(htsjdk_bam::tag::Tag::new(name));
}

/// `QualityUtils.MAX_SAM_QUAL_SCORE`, re-exported so a reader of this module sees the ceiling the
/// recalibrated qualities are clamped to.
pub const RECALIBRATED_CEILING: i32 = MAX_SAM_QUAL_SCORE;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recal_datum::RecalDatum;
    use std::cell::RefCell;
    use std::rc::Rc;

    fn datum(observations: i64, mismatches: f64, quality: i8) -> SharedDatum {
        Rc::new(RefCell::new(
            RecalDatum::new(observations, mismatches, quality).unwrap(),
        ))
    }

    #[test]
    fn a_null_datum_contributes_nothing_rather_than_a_zero_delta() {
        // Every datum absent: the prior comes straight back.
        assert_eq!(
            hierarchical_bayesian_quality_estimate(25.0, None, None, &[None, None]),
            25.0
        );
        // A special covariate with no read group or quality datum takes its delta against the
        // prior itself.
        let special = datum(1000, 1.0, 30);
        let estimate =
            hierarchical_bayesian_quality_estimate(25.0, None, None, &[Some(special), None]);
        assert_eq!(estimate, 25.0);
    }

    #[test]
    fn the_datum_cache_makes_the_estimate_order_dependent() {
        let shared = datum(1000, 1.0, 30);
        let first = hierarchical_bayesian_quality_estimate(25.0, None, Some(&shared), &[]);
        let second = hierarchical_bayesian_quality_estimate(45.0, None, Some(&shared), &[]);
        assert_eq!(first, second, "the second prior is ignored");
        let fresh = datum(1000, 1.0, 30);
        let alone = hierarchical_bayesian_quality_estimate(45.0, None, Some(&fresh), &[]);
        assert_ne!(
            alone, second,
            "a fresh datum answers the prior it was given"
        );
    }

    #[test]
    fn the_rounding_is_fast_round_and_then_the_clamp() {
        // Rounded twice, where Math.round answers 0.
        assert_eq!(bounded_integer_qual(0.499_999_999_999_999_94), 1);
        assert_eq!(bounded_integer_qual(0.0), 1);
        assert_eq!(bounded_integer_qual(-1.0), 1);
        assert_eq!(bounded_integer_qual(93.5), 93);
        assert_eq!(bounded_integer_qual(200.0), 93);
        assert_eq!(bounded_integer_qual(1.5), 2);
    }

    #[test]
    fn the_static_mapping_sorts_the_callers_list() {
        let mut quals = vec![40, 10, 30, 20];
        let mapping = construct_static_quantized_mapping(&mut quals, false);
        assert_eq!(quals, vec![10, 20, 30, 40]);
        // Everything below the usable minimum is itself.
        for (i, value) in mapping.iter().enumerate().take(MIN_USABLE_Q_SCORE as usize) {
            assert_eq!(*value, i as u8);
        }
        // And everything above the largest static quality is that quality.
        assert_eq!(mapping[93], 40);
    }

    #[test]
    fn one_static_quality_swallows_everything_above_the_minimum() {
        let mut quals = vec![30];
        let mapping = construct_static_quantized_mapping(&mut quals, false);
        assert_eq!(mapping[5], 5);
        assert_eq!(mapping[6], 30);
        assert_eq!(mapping[93], 30);
    }
}
