//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller`
//! (`PairHMMLikelihoodCalculationEngine`, `ReadLikelihoodCalculationEngine`'s statics and
//! defaults, `StandardPairHMMInputScoreImputator`) and `utils.pairhmm.PairHMM.computeLog10Likelihoods`,
//! GATK 4.6.2.0: every read against every haplotype, the way HaplotypeCaller and Mutect2 score
//! their reads.
//!
//! Each read is prepared, scored with the Java `LoglessPairHMM` ([`crate::pair_hmm`]), and the
//! matrix is then normalized and filtered:
//!
//! 1. the soft clips are hard-clipped away unless soft-clipped bases are kept;
//! 2. the PCR indel error model caps each base's insertion and deletion quality by the length of
//!    the tandem repeat that ends there;
//! 3. base qualities are capped by the mapping quality and squashed to Q6 under the threshold
//!    (Q18 by default), indel qualities squashed to Q6 under Q6;
//! 4. the PairHMM runs with those qualities and a constant gap continuation penalty;
//! 5. `normalizeLikelihoods` keeps every likelihood within the mismapping rate of the best;
//! 6. `filterPoorlyModeledEvidence` removes the reads no haplotype explains.
//!
//! What a reader would not guess:
//!
//! * the PairHMM keeps its matrices between haplotypes and only re-seeds the deletion row when the
//!   haplotype length changes, but `computeLog10Likelihoods` always passes `recacheReadValues`, so
//!   no state crosses a call and the stateless recursion gives the same bits;
//! * the filtering threshold reads the length of the PREPARED qualities, so a hard-clipped read is
//!   judged by its clipped length; the quality comparisons are on signed bytes, as in Java;
//! * `findTandemRepeatUnits` counts the repeat that ENDS at the base backwards and the one that
//!   starts after it forwards, adds them only when the units agree, and otherwise re-counts the
//!   forward unit backwards over the whole prefix.

use htsjdk_bam::record::BamRecord;

use crate::allele_likelihoods::{AlleleLikelihoods, LikelihoodsError};
use crate::allele_list::{AlleleList, SampleList};
use crate::clipping::{hard_clip_soft_clipped_bases, ClipError};
use crate::haplotype::Haplotype;
use crate::pair_hmm::{fast_round, read_likelihood_given_haplotype_log10};
use crate::read_utils::{base_deletion_qualities, base_insertion_qualities};

/// `ReadLikelihoodCalculationEngine.MAX_STR_UNIT_LENGTH`.
pub const MAX_STR_UNIT_LENGTH: usize = 8;
/// `ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH`.
pub const MAX_REPEAT_LENGTH: i32 = 20;
/// `ReadLikelihoodCalculationEngine.DEFAULT_EXPECTED_ERROR_RATE_PER_BASE`.
pub const DEFAULT_EXPECTED_ERROR_RATE_PER_BASE: f64 = 0.02;
/// `PairHMM.BASE_QUALITY_SCORE_THRESHOLD`.
pub const BASE_QUALITY_SCORE_THRESHOLD: u8 = 18;
/// `QualityUtils.MIN_USABLE_Q_SCORE`.
pub const MIN_USABLE_Q_SCORE: u8 = 6;
/// `MIN_ADJUSTED_QSCORE`.
const MIN_ADJUSTED_QSCORE: i32 = 10;
/// `INITIAL_QSCORE`.
const INITIAL_QSCORE: f64 = 40.0;

/// `PCRErrorModel`, with its rate factor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PcrErrorModel {
    None,
    Hostile,
    Aggressive,
    Conservative,
}

impl PcrErrorModel {
    pub fn rate_factor(self) -> f64 {
        match self {
            PcrErrorModel::None => 0.0,
            PcrErrorModel::Hostile => 1.0,
            PcrErrorModel::Aggressive => 2.0,
            PcrErrorModel::Conservative => 3.0,
        }
    }

    pub fn has_rate_factor(self) -> bool {
        self.rate_factor() != 0.0
    }
}

/// What the engine refuses, with the reference's exception.
#[derive(Debug, Clone, PartialEq)]
pub enum LikelihoodEngineError {
    /// `IllegalArgumentException("gap continuation penalty must be non-negative")`.
    NegativeGapContinuationPenalty,
    /// `IllegalArgumentException("log10globalReadMismappingRate must be negative")`.
    PositiveMismappingRate,
    /// `IllegalArgumentException("baseQualityScoreThreshold must be greater than or equal to 6 ...")`.
    BaseQualityThresholdTooLow,
    /// `PairHMM.initialize`: `haplotypeMaxLength must be > 0`.
    NoHaplotypeBases,
    /// `IllegalStateException`: a PairHMM answer above zero, or not a number.
    InvalidLikelihood(f64),
    Clip(ClipError),
    Likelihoods(LikelihoodsError),
}

/// A read as the PairHMM sees it: `createQualityModifiedRead`'s bases and three quality arrays.
#[derive(Debug, Clone, PartialEq)]
pub struct ProcessedRead {
    pub bases: Vec<u8>,
    pub qualities: Vec<u8>,
    pub insertion_qualities: Vec<u8>,
    pub deletion_qualities: Vec<u8>,
}

/// `PairHMMLikelihoodCalculationEngine`'s arguments, defaulted as HaplotypeCaller defaults them.
#[derive(Debug, Clone, PartialEq)]
pub struct LikelihoodEngineArguments {
    pub gap_continuation_penalty: i8,
    pub log10_global_read_mismapping_rate: f64,
    pub pcr_error_model: PcrErrorModel,
    pub base_quality_score_threshold: i8,
    pub dynamic_disqualification: bool,
    pub read_disqualification_scale: f64,
    pub expected_error_rate_per_base: f64,
    pub symmetrically_normalize_alleles_to_reference: bool,
    pub disable_cap_read_qualities_to_map_q: bool,
    pub modify_soft_clipped_bases: bool,
}

impl Default for LikelihoodEngineArguments {
    fn default() -> Self {
        LikelihoodEngineArguments {
            gap_continuation_penalty: 10,
            log10_global_read_mismapping_rate: -4.5,
            pcr_error_model: PcrErrorModel::Conservative,
            base_quality_score_threshold: BASE_QUALITY_SCORE_THRESHOLD as i8,
            dynamic_disqualification: false,
            read_disqualification_scale: 1.0,
            expected_error_rate_per_base: DEFAULT_EXPECTED_ERROR_RATE_PER_BASE,
            symmetrically_normalize_alleles_to_reference: true,
            disable_cap_read_qualities_to_map_q: false,
            modify_soft_clipped_bases: false,
        }
    }
}

/// `PairHMMLikelihoodCalculationEngine` with the Java `LOGLESS_CACHING` PairHMM and no DRAGstr
/// parameters.
#[derive(Debug, Clone)]
pub struct PairHmmLikelihoodEngine {
    arguments: LikelihoodEngineArguments,
    /// `pcrIndelErrorModelCache`, indexed by repeat length; empty without a rate factor.
    pcr_indel_error_model_cache: Vec<u8>,
}

/// `getErrorModelAdjustedQual(repeatLength, rateFactor)`.
pub fn error_model_adjusted_qual(repeat_length: i32, rate_factor: f64) -> u8 {
    // `Math.exp`; the answer is rounded to an integer, and the `pcr` rows measure every value.
    let penalty =
        jmath::strict_exp::exp(f64::from(repeat_length) / (rate_factor * std::f64::consts::PI));
    MIN_ADJUSTED_QSCORE.max(fast_round(INITIAL_QSCORE - penalty + 1.0)) as u8
}

impl PairHmmLikelihoodEngine {
    /// The fourteen-argument constructor, `resultsFile` absent.
    pub fn new(arguments: LikelihoodEngineArguments) -> Result<Self, LikelihoodEngineError> {
        if arguments.gap_continuation_penalty < 0 {
            return Err(LikelihoodEngineError::NegativeGapContinuationPenalty);
        }
        if arguments.log10_global_read_mismapping_rate > 0.0 {
            return Err(LikelihoodEngineError::PositiveMismappingRate);
        }
        let mut cache = Vec::new();
        if arguments.pcr_error_model.has_rate_factor() {
            let rate = arguments.pcr_error_model.rate_factor();
            cache = (0..=MAX_REPEAT_LENGTH)
                .map(|i| error_model_adjusted_qual(i, rate))
                .collect();
        }
        if arguments.base_quality_score_threshold < MIN_USABLE_Q_SCORE as i8 {
            return Err(LikelihoodEngineError::BaseQualityThresholdTooLow);
        }
        Ok(PairHmmLikelihoodEngine {
            arguments,
            pcr_indel_error_model_cache: cache,
        })
    }

    /// `modifyReadQualities`: the prepared copies, in order. The base qualities of each are what
    /// the reference stores on the original read as `HMMQuals`.
    pub fn modify_read_qualities(
        &self,
        reads: &[BamRecord],
    ) -> Result<Vec<ProcessedRead>, LikelihoodEngineError> {
        let mut result = Vec::with_capacity(reads.len());
        for read in reads {
            let clipped;
            let maybe_unclipped = if self.arguments.modify_soft_clipped_bases {
                read
            } else {
                clipped = hard_clip_soft_clipped_bases(read, None, 0)
                    .map_err(LikelihoodEngineError::Clip)?;
                &clipped
            };
            let bases = maybe_unclipped.read_bases.clone();
            let mut qualities = maybe_unclipped.base_qualities.clone();
            let mut insertion = base_insertion_qualities(maybe_unclipped);
            let mut deletion = base_deletion_qualities(maybe_unclipped);
            self.apply_pcr_error_model(&bases, &mut insertion, &mut deletion);
            self.cap_minimum_read_qualities(
                maybe_unclipped.mapping_quality,
                &mut qualities,
                &mut insertion,
                &mut deletion,
            );
            result.push(ProcessedRead {
                bases,
                qualities,
                insertion_qualities: insertion,
                deletion_qualities: deletion,
            });
        }
        Ok(result)
    }

    /// `applyPCRErrorModel`: the base before each position takes the cache's quality for the
    /// tandem repeat ending there, when lower.
    fn apply_pcr_error_model(&self, bases: &[u8], insertion: &mut [u8], deletion: &mut [u8]) {
        if self.arguments.pcr_error_model == PcrErrorModel::None {
            return;
        }
        for i in 1..bases.len() {
            let (_, repeat_length) = find_tandem_repeat_units(bases, i - 1);
            let cap = self.pcr_indel_error_model_cache[repeat_length as usize];
            insertion[i - 1] = insertion[i - 1].min(cap);
            deletion[i - 1] = deletion[i - 1].min(cap);
        }
    }

    /// `capMinimumReadQualities`; `setToFixedValueIfTooLow` compares signed bytes.
    fn cap_minimum_read_qualities(
        &self,
        mapping_quality: u8,
        qualities: &mut [u8],
        insertion: &mut [u8],
        deletion: &mut [u8],
    ) {
        let fixed = |value: u8, min: i8| {
            if (value as i8) < min {
                MIN_USABLE_Q_SCORE
            } else {
                value
            }
        };
        for i in 0..qualities.len() {
            if !self.arguments.disable_cap_read_qualities_to_map_q {
                qualities[i] = qualities[i].min(mapping_quality);
            }
            qualities[i] = fixed(qualities[i], self.arguments.base_quality_score_threshold);
            insertion[i] = fixed(insertion[i], MIN_USABLE_Q_SCORE as i8);
            deletion[i] = fixed(deletion[i], MIN_USABLE_Q_SCORE as i8);
        }
    }

    /// `computeReadLikelihoods(haplotypeList, header, samples, perSampleReadList, filterPoorly)`:
    /// `reads` holds each sample's reads in sample order.
    pub fn compute_read_likelihoods(
        &self,
        haplotypes: &[Haplotype],
        samples: &[String],
        reads: &[Vec<BamRecord>],
    ) -> Result<AlleleLikelihoods<BamRecord, Haplotype>, LikelihoodEngineError> {
        // `initializePairHMM`: refused for haplotypes with no bases, before any read is looked at.
        let haplotype_max_length = haplotypes
            .iter()
            .map(|h| h.bases().len())
            .max()
            .unwrap_or(0);
        if haplotype_max_length == 0 {
            return Err(LikelihoodEngineError::NoHaplotypeBases);
        }
        let haplotype_bases: Vec<Vec<u8>> = haplotypes.iter().map(Haplotype::bases).collect();

        let mut values = Vec::with_capacity(samples.len());
        // The length of each read's `HMMQuals`, which the fixed threshold reads, and the qualities,
        // which the dynamic one reads.
        let mut hmm_qualities: Vec<Vec<Vec<u8>>> = Vec::with_capacity(samples.len());
        for sample_reads in reads {
            let processed = self.modify_read_qualities(sample_reads)?;
            let mut sample_values = vec![vec![0.0; processed.len()]; haplotypes.len()];
            for (r, read) in processed.iter().enumerate() {
                let gap_continuation =
                    vec![self.arguments.gap_continuation_penalty as u8; read.bases.len()];
                for (a, haplotype) in haplotype_bases.iter().enumerate() {
                    let likelihood = read_likelihood_given_haplotype_log10(
                        haplotype,
                        &read.bases,
                        &read.qualities,
                        &read.insertion_qualities,
                        &read.deletion_qualities,
                        &gap_continuation,
                    );
                    if likelihood > 0.0 || likelihood.is_nan() {
                        return Err(LikelihoodEngineError::InvalidLikelihood(likelihood));
                    }
                    sample_values[a][r] = likelihood;
                }
            }
            hmm_qualities.push(processed.into_iter().map(|p| p.qualities).collect());
            values.push(sample_values);
        }

        let mut result = AlleleLikelihoods::new(
            SampleList::new(samples),
            AlleleList::new(haplotypes),
            reads.to_vec(),
            values,
        )
        .map_err(LikelihoodEngineError::Likelihoods)?;
        result
            .normalize_likelihoods(
                self.arguments.log10_global_read_mismapping_rate,
                self.arguments.symmetrically_normalize_alleles_to_reference,
            )
            .map_err(LikelihoodEngineError::Likelihoods)?;

        let fixed = |length: usize, cap: bool| {
            log10_min_true_likelihood(length, self.arguments.expected_error_rate_per_base, cap)
        };
        let dynamic = self.arguments.dynamic_disqualification;
        let scale = self.arguments.read_disqualification_scale;
        result
            .filter_poorly_modeled_evidence(|sample, evidence, _| {
                let qualities = &hmm_qualities[sample][evidence];
                if dynamic {
                    let threshold = log10_dynamic_read_qual_threshold(qualities, scale);
                    let maximum = fixed(qualities.len(), false);
                    if threshold < maximum {
                        threshold
                    } else {
                        maximum
                    }
                } else {
                    fixed(qualities.len(), true)
                }
            })
            .map_err(LikelihoodEngineError::Likelihoods)?;
        // `filterPoorlyModeledEvidence` removed entries; the HMM qualities it read are dropped
        // with them by the reference too, being transient attributes of the removed reads.
        Ok(result)
    }
}

/// `log10MinTrueLikelihood(maximumErrorPerBase, capLikelihoods)` for a read whose prepared
/// qualities are `length` long: Q40 for each tolerated error, at most two when capped.
pub fn log10_min_true_likelihood(length: usize, maximum_error_per_base: f64, cap: bool) -> f64 {
    let errors = (length as f64 * maximum_error_per_base).ceil();
    let errors = if cap { 2.0f64.min(errors) } else { errors };
    errors * -4.0
}

/// `calculateLog10DynamicReadQualThreshold`: the table's means summed, plus the scale times the
/// square root of its variances summed, as a phred score.
pub fn log10_dynamic_read_qual_threshold(qualities: &[u8], scale: f64) -> f64 {
    let mut sum_mean = 0.0;
    let mut sum_variance = 0.0;
    for &quality in qualities {
        let bq = quality as usize;
        let entry = if bq <= 1 {
            0
        } else {
            bq.min(MAXIMUM_DYNAMIC_QUAL_THRESHOLD_ENTRY_BASEQ) - 1
        };
        let mean = entry * 3 + 1;
        sum_mean += DYNAMIC_READ_QUAL_THRESH_LOOKUP_TABLE[mean];
        sum_variance += DYNAMIC_READ_QUAL_THRESH_LOOKUP_TABLE[mean + 1];
    }
    let threshold = sum_mean + scale * jmath::math::sqrt(sum_variance);
    // `QualityUtils.qualToErrorProbLog10`.
    threshold * -0.1
}

const MAXIMUM_DYNAMIC_QUAL_THRESHOLD_ENTRY_BASEQ: usize = 40;

/// `dynamicReadQualThreshLookupTable`: base quality, mean, variance.
#[rustfmt::skip]
const DYNAMIC_READ_QUAL_THRESH_LOOKUP_TABLE: [f64; 120] = [
    1.0,  5.996842844, 0.196616587, 2.0,  5.870018422, 1.388545569, 3.0,  5.401558531, 5.641990128,
    4.0,  4.818940919, 10.33176216, 5.0,  4.218758304, 14.25799688, 6.0,  3.646319832, 17.02880749,
    7.0,  3.122346753, 18.64537883, 8.0,  2.654731979, 19.27521677, 9.0,  2.244479156, 19.13584613,
    10.0, 1.88893867,  18.43922003, 11.0, 1.583645342, 17.36842261, 12.0, 1.3233807, 16.07088712,
    13.0, 1.102785365, 14.65952563, 14.0, 0.916703025, 13.21718577, 15.0, 0.760361881, 11.80207947,
    16.0, 0.629457387, 10.45304833, 17.0, 0.520175654, 9.194183767, 18.0, 0.42918208,  8.038657241,
    19.0, 0.353590663, 6.991779595, 20.0, 0.290923699, 6.053379213, 21.0, 0.23906788,  5.219610436,
    22.0, 0.196230431, 4.484302033, 23.0, 0.160897421, 3.839943445, 24.0, 0.131795374, 3.27839108,
    25.0, 0.1078567,   2.791361596, 26.0, 0.088189063, 2.370765375, 27.0, 0.072048567, 2.008921719,
    28.0, 0.058816518, 1.698687797, 29.0, 0.047979438, 1.433525748, 30.0, 0.039111985, 1.207526336,
    31.0, 0.031862437, 1.015402928, 32.0, 0.025940415, 0.852465956, 33.0, 0.021106532, 0.714585285,
    34.0, 0.017163711, 0.598145851, 35.0, 0.013949904, 0.500000349, 36.0, 0.011332027, 0.41742159,
    37.0, 0.009200898, 0.348056286, 38.0, 0.007467036, 0.289881373, 39.0, 0.006057179, 0.241163527,
    40.0, 0.004911394, 0.200422214,
];

/// `GATKVariantContextUtils.findNumberOfRepetitions(repeatUnit, testString, leadingRepeats)`; the
/// offset form the engine calls is this one over sub-slices.
fn find_number_of_repetitions(
    repeat_unit: &[u8],
    test_string: &[u8],
    leading_repeats: bool,
) -> i32 {
    if test_string.is_empty() {
        return 0;
    }
    let unit = repeat_unit.len() as i64;
    let length_difference = test_string.len() as i64 - unit;
    let mut repeats = 0;
    let mut start = if leading_repeats {
        0
    } else {
        length_difference
    };
    while (leading_repeats && start <= length_difference) || (!leading_repeats && start >= 0) {
        if test_string[start as usize..(start + unit) as usize] != *repeat_unit {
            return repeats;
        }
        repeats += 1;
        start += if leading_repeats { unit } else { -unit };
    }
    repeats
}

/// `ReadLikelihoodCalculationEngine.findTandemRepeatUnits(readBases, offset)`: the repeat unit
/// around `offset` and its length in units, capped at `MAX_REPEAT_LENGTH`.
pub fn find_tandem_repeat_units(bases: &[u8], offset: usize) -> (Vec<u8>, i32) {
    let mut max_bw = 0;
    let mut best_bw_unit = vec![bases[offset]];
    for str_len in 1..=MAX_STR_UNIT_LENGTH {
        if offset + 1 < str_len {
            break;
        }
        let unit = &bases[offset + 1 - str_len..offset + 1];
        max_bw = find_number_of_repetitions(unit, &bases[..offset + 1], false);
        if max_bw > 1 {
            best_bw_unit = unit.to_vec();
            break;
        }
    }
    let mut best_unit = best_bw_unit.clone();
    let mut max_rl = max_bw;
    if offset < bases.len() - 1 {
        let mut best_fw_unit = vec![bases[offset + 1]];
        let mut max_fw = 0;
        for str_len in 1..=MAX_STR_UNIT_LENGTH {
            if offset + str_len + 1 > bases.len() {
                break;
            }
            let unit = &bases[offset + 1..offset + str_len + 1];
            max_fw = find_number_of_repetitions(unit, &bases[offset + 1..], true);
            if max_fw > 1 {
                best_fw_unit = unit.to_vec();
                break;
            }
        }
        if best_fw_unit == best_bw_unit {
            max_rl = max_bw + max_fw;
        } else {
            // The backward unit may still be part of the forward one: `TTCTT(C) CCC` at (C) is
            // `(TTC)2` backwards and `(C)3` forwards, and is `(C)4`.
            let backwards = find_number_of_repetitions(&best_fw_unit, &bases[..offset + 1], false);
            max_rl = max_fw + backwards;
        }
        best_unit = best_fw_unit;
    }
    (best_unit, max_rl.min(MAX_REPEAT_LENGTH))
}
