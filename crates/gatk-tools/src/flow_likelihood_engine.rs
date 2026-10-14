//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.
//! FlowBasedAlignmentLikelihoodEngine` and `org.broadinstitute.hellbender.utils.haplotype.
//! FlowBasedHaplotype` (GATK 4.6.2.0): the callers' `FlowBased` read likelihoods.
//!
//! Each sample's reads become flow reads under the flow order and maximal hmer of the sample's
//! FIRST read's read group, whichever group the others belong to; each haplotype becomes a key in
//! that order. A read is trimmed, in bases, to the first haplotype's genome span, which marks it
//! invalid (negative infinity against every haplotype) when fewer than ten bases would remain.
//! The haplotype is clipped to the flows the read's trimmed span covers on it (through the
//! haplotype's CIGAR), widened by four flows each way, and the read's matrix is read along it from
//! the first flow that reads the read's first base, stepping four flows, the best sum kept.
//!
//! The matrix is then normalized to the reference within the mismapping rate, and the poorly
//! modelled reads dropped against a threshold of their own: at most three errors at the expected
//! rate per base and two at the filling value, both base-ten.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use gatk_engine::haplotype::Haplotype;
use gatk_engine::pair_hmm_likelihood_engine::log10_dynamic_read_qual_threshold;

use crate::flow_based_read::{find_clipping, FlowArguments, FlowRead, FlowReadError};
use crate::flow_pairhmm_align_reads_to_haplotypes::{clipped_score, FlowHaplotype};

/// `MAX_ERRORS_FOR_READ_CAP`.
const MAX_ERRORS_FOR_READ_CAP: f64 = 3.0;
/// `MAX_CATASTROPHIC_ERRORS_FOR_READ_CAP`.
const MAX_CATASTROPHIC_ERRORS_FOR_READ_CAP: f64 = 2.0;
/// `ALIGNMENT_UNCERTAINTY`.
const ALIGNMENT_UNCERTAINTY: i32 = 4;
/// `FlowBasedRead.MINIMAL_READ_LENGTH`.
const MINIMAL_READ_LENGTH: i32 = 10;

/// `FlowBasedAlignmentArgumentCollection`: the matrix's arguments and the engine's own.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowEngineArguments {
    pub flow: FlowArguments,
    /// `--flow-order-cycle-length`.
    pub flow_order_cycle_length: usize,
    /// `--flow-likelihood-optimized-comp`.
    pub optimized: bool,
    /// `--trim-to-haplotype`.
    pub trim_to_haplotype: bool,
    /// `--exact-matching`.
    pub exact_matching: bool,
}

impl Default for FlowEngineArguments {
    fn default() -> Self {
        FlowEngineArguments {
            flow: FlowArguments::default(),
            flow_order_cycle_length: 4,
            optimized: false,
            trim_to_haplotype: true,
            exact_matching: false,
        }
    }
}

/// `FlowBasedAlignmentLikelihoodEngine`.
#[derive(Debug, Clone)]
pub struct FlowLikelihoodEngine {
    pub arguments: FlowEngineArguments,
    pub log10_global_read_mismapping_rate: f64,
    pub expected_error_rate_per_base: f64,
    pub dynamic_disqualification: bool,
    pub read_disqualification_scale: f64,
}

/// Why the engine stopped: the reference's exception, or a layout the port does not read.
#[derive(Debug, Clone, PartialEq)]
pub enum FlowEngineError {
    Read(FlowReadError),
    /// `GATKException.ShouldNeverReachHereException` and `IllegalStateException`.
    Thrown(&'static str, String),
    Likelihoods(String),
}

impl FlowEngineError {
    /// A reference exception as `class: message`; anything else as the port describes it.
    pub fn render(&self) -> String {
        match self {
            FlowEngineError::Read(error) if !error.port_limitation && !error.class.is_empty() => {
                format!("{}: {}", error.class, error.message)
            }
            FlowEngineError::Thrown(class, message) => {
                format!("{}: {message}", crate::hc_region::qualified_class(class))
            }
            other => format!("{other:?}"),
        }
    }
}

/// A read as the engine scores it: its flow matrix, and its alignment start and end, length and
/// trimming in bases.
struct ScoredRead {
    flow: FlowRead,
    start: i32,
    end: i32,
    length: i32,
    trim_left: i32,
    trim_right: i32,
    base_clipped: bool,
}

/// A haplotype as the engine scores it.
struct ScoredHaplotype {
    flow: FlowHaplotype,
    start: i32,
    end: i32,
    length: i32,
    cigar: Option<htsjdk_bam::cigar::Cigar>,
}

impl FlowLikelihoodEngine {
    /// `computeReadLikelihoods(haplotypes, header, samples, perSampleReadList, filterPoorly)`,
    /// normalized and filtered.
    pub fn compute_read_likelihoods(
        &self,
        haplotypes: &[Haplotype],
        header: &SamHeader,
        samples: &[String],
        reads: &[Vec<BamRecord>],
    ) -> Result<AlleleLikelihoods<BamRecord, Haplotype>, FlowEngineError> {
        let mut values = Vec::with_capacity(samples.len());
        for sample_reads in reads {
            values.push(self.sample_matrix(haplotypes, header, sample_reads)?);
        }
        let mut result = AlleleLikelihoods::new(
            SampleList::new(samples),
            AlleleList::new(haplotypes),
            reads.to_vec(),
            values,
        )
        .map_err(|e| FlowEngineError::Likelihoods(format!("{e:?}")))?;
        result
            .normalize_likelihoods(self.log10_global_read_mismapping_rate, true)
            .map_err(|e| FlowEngineError::Likelihoods(format!("{e:?}")))?;
        let expected = self.expected_error_rate_per_base;
        let filling = self.arguments.flow.filling_value;
        let dynamic = self.dynamic_disqualification;
        let scale = self.read_disqualification_scale;
        result
            .filter_poorly_modeled_evidence(|_, _, read: &BamRecord| {
                let fixed = |cap: bool| {
                    log10_min_true_likelihood(read.read_bases.len(), expected, filling, cap)
                };
                if dynamic {
                    let threshold = log10_dynamic_read_qual_threshold(&read.base_qualities, scale);
                    let maximum = fixed(false);
                    if threshold < maximum {
                        threshold
                    } else {
                        maximum
                    }
                } else {
                    fixed(true)
                }
            })
            .map_err(|e| FlowEngineError::Likelihoods(format!("{e:?}")))?;
        Ok(result)
    }

    /// `computeReadLikelihoods(likelihoodMatrix, header)`: one sample's values, by haplotype.
    fn sample_matrix(
        &self,
        haplotypes: &[Haplotype],
        header: &SamHeader,
        reads: &[BamRecord],
    ) -> Result<Vec<Vec<f64>>, FlowEngineError> {
        let cycle = self.arguments.flow_order_cycle_length;
        let (flow_order, max_class) = match reads.first() {
            Some(first) => {
                let info = crate::flow_based_read::read_group_info(first, header)
                    .map_err(FlowEngineError::Read)?;
                (substring(&info.flow_order, cycle)?, info.max_class)
            }
            None => (first_usable_flow_order(header, cycle)?, 0),
        };
        let mut processed = Vec::with_capacity(reads.len());
        for read in reads {
            let mut flow = FlowRead::new(read, &flow_order, max_class, &self.arguments.flow)
                .map_err(FlowEngineError::Read)?;
            flow.apply_alignment().map_err(FlowEngineError::Read)?;
            processed.push(ScoredRead {
                flow,
                start: gatk_engine::read_utils::start(read),
                end: gatk_engine::read_utils::end(read),
                length: read.read_bases.len() as i32,
                trim_left: 0,
                trim_right: 0,
                base_clipped: false,
            });
        }
        let mut scored_haplotypes = Vec::with_capacity(haplotypes.len());
        for haplotype in haplotypes {
            let bases = haplotype.bases();
            let flow = FlowHaplotype::new(&bases, &flow_order).ok_or_else(|| {
                FlowEngineError::Thrown(
                    "org.broadinstitute.hellbender.exceptions.GATKException",
                    "baseArrayToKey".to_string(),
                )
            })?;
            let location = haplotype.genome_location();
            scored_haplotypes.push(ScoredHaplotype {
                flow,
                start: location.map(|l| l.start).unwrap_or(0),
                end: location.map(|l| l.end).unwrap_or(0),
                length: bases.len() as i32,
                cigar: haplotype.cigar().cloned(),
            });
        }
        if self.arguments.trim_to_haplotype {
            if let Some(first) = scored_haplotypes.first() {
                let (haplotype_start, haplotype_end) = (first.start, first.end);
                for read in &mut processed {
                    let left = (haplotype_start - read.start).max(0);
                    let right = (read.end - haplotype_end).max(0);
                    read.flow
                        .apply_base_clipping(left, right, true, read.length, MINIMAL_READ_LENGTH)
                        .map_err(FlowEngineError::Read)?;
                    read.base_clipped = true;
                    read.trim_left = left;
                    read.trim_right = right;
                }
            }
        }
        let mut values = vec![vec![0.0; processed.len()]; haplotypes.len()];
        for (h, haplotype) in scored_haplotypes.iter().enumerate() {
            for (r, read) in processed.iter().enumerate() {
                values[h][r] = if self.arguments.exact_matching {
                    self.exact_matching(haplotype, read)?
                } else {
                    self.haplotype_read_matching(haplotype, read)?
                };
            }
        }
        Ok(values)
    }

    fn check_read(&self, read: &ScoredRead) -> Result<(), FlowEngineError> {
        if !read.base_clipped && self.arguments.trim_to_haplotype {
            return Err(FlowEngineError::Thrown(
                "org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException",
                "Reads should be trimmed to the haplotype".to_string(),
            ));
        }
        Ok(())
    }

    /// `haplotypeReadMatchingExactLength`.
    fn exact_matching(
        &self,
        haplotype: &ScoredHaplotype,
        read: &ScoredRead,
    ) -> Result<f64, FlowEngineError> {
        self.check_read(read)?;
        if !read.flow.valid {
            return Ok(f64::NEG_INFINITY);
        }
        Ok(clipped_score(
            &haplotype.flow,
            &read.flow,
            0,
            haplotype.flow.key.len(),
            self.arguments.optimized,
        ))
    }

    /// `haplotypeReadMatching`.
    fn haplotype_read_matching(
        &self,
        haplotype: &ScoredHaplotype,
        read: &ScoredRead,
    ) -> Result<f64, FlowEngineError> {
        self.check_read(read)?;
        if !read.flow.valid {
            return Ok(f64::NEG_INFINITY);
        }
        let cigar = haplotype.cigar.clone().unwrap_or_default();
        let index = |coordinate: i32| {
            gatk_engine::read_utils::read_index_for_reference_coordinate(
                haplotype.start,
                &cigar,
                coordinate,
            )
            .0
        };
        let haplotype_start = index(read.start + read.trim_left);
        let haplotype_end = index(read.end - read.trim_right);
        let haplotype_length = haplotype_end - haplotype_start;
        let uncertainty = (read.length - haplotype_length).max(0);
        let left_clip = (haplotype_start - uncertainty).max(0);
        let right_clip = (haplotype.length - haplotype_end - 1 - uncertainty).max(0);
        if left_clip >= haplotype.length || right_clip >= haplotype.length {
            return Ok(f64::NEG_INFINITY);
        }
        let (clip_left, left_hmer_clip) = find_clipping(left_clip, &haplotype.flow.key);
        let reversed: Vec<i32> = haplotype.flow.key.iter().rev().copied().collect();
        let (clip_right, right_hmer_clip) = find_clipping(right_clip, &reversed);
        let key_length = haplotype.flow.key.len() as i32;
        if clip_left >= key_length || clip_right >= key_length {
            return Ok(f64::NEG_INFINITY);
        }
        if left_hmer_clip < 0 || right_hmer_clip < 0 {
            return Err(FlowEngineError::Thrown(
                "org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException",
                "Negative hmer clips found. Check".to_string(),
            ));
        }
        let clip_left = if clip_left - ALIGNMENT_UNCERTAINTY > 0 {
            clip_left - ALIGNMENT_UNCERTAINTY
        } else {
            0
        };
        let clip_right = if key_length - clip_right + ALIGNMENT_UNCERTAINTY < key_length {
            key_length - clip_right + ALIGNMENT_UNCERTAINTY
        } else {
            key_length
        };
        Ok(clipped_score(
            &haplotype.flow,
            &read.flow,
            clip_left as usize,
            clip_right as usize,
            self.arguments.optimized,
        ))
    }
}

/// The engine's `log10MinTrueLikelihood`: the expected errors, at least three when capped, at
/// the expected rate, and the catastrophic ones, at least two when capped, at the filling value.
pub fn log10_min_true_likelihood(
    length: usize,
    expected_error_rate: f64,
    filling_value: f64,
    cap: bool,
) -> f64 {
    let log10_error_rate = jmath::math::log10(expected_error_rate);
    let large_event_error_rate = filling_value.max(0.000001);
    let log10_catastrophic = jmath::math::log10(large_event_error_rate);
    let errors = (length as f64 * expected_error_rate).ceil();
    let catastrophic = (length as f64 * large_event_error_rate).ceil();
    let (errors, catastrophic) = if cap {
        (
            MAX_ERRORS_FOR_READ_CAP.max(errors),
            MAX_CATASTROPHIC_ERRORS_FOR_READ_CAP.max(catastrophic),
        )
    } else {
        (errors, catastrophic)
    };
    errors * log10_error_rate + catastrophic * log10_catastrophic
}

/// `flowOrder.substring(0, cycle)`.
fn substring(order: &str, cycle: usize) -> Result<String, FlowEngineError> {
    order.get(..cycle).map(str::to_string).ok_or_else(|| {
        FlowEngineError::Thrown(
            "java.lang.StringIndexOutOfBoundsException",
            format!("begin 0, end {cycle}, length {}", order.len()),
        )
    })
}

/// `FlowBasedReadUtils.findFirstUsableFlowOrder`.
fn first_usable_flow_order(header: &SamHeader, cycle: usize) -> Result<String, FlowEngineError> {
    for group in &header.read_groups {
        if let Some(order) = group.attributes.get("FO") {
            if order.len() >= cycle {
                return Ok(order[..cycle].to_string());
            }
        }
    }
    Err(FlowEngineError::Thrown(
        "org.broadinstitute.hellbender.exceptions.GATKException",
        "Unable to perform flow based operations without the flow order".to_string(),
    ))
}

/// `ReadLikelihoodCalculationEngine`, as `createLikelihoodCalculationEngine` picks it.
#[derive(Debug, Clone)]
pub enum ReadLikelihoodEngine {
    PairHmm(gatk_engine::pair_hmm_likelihood_engine::PairHmmLikelihoodEngine),
    FlowBased(FlowLikelihoodEngine),
    FlowBasedHmm(FlowHmmEngine),
}

impl ReadLikelihoodEngine {
    /// `computeReadLikelihoods(haplotypes, header, samples, perSampleReadList, true)`, the error
    /// as the reference's class and message would read.
    pub fn compute_read_likelihoods(
        &self,
        haplotypes: &[Haplotype],
        header: &SamHeader,
        samples: &[String],
        reads: &[Vec<BamRecord>],
    ) -> Result<AlleleLikelihoods<BamRecord, Haplotype>, String> {
        match self {
            ReadLikelihoodEngine::PairHmm(engine) => engine
                .compute_read_likelihoods(haplotypes, samples, reads)
                .map_err(|e| format!("{e:?}")),
            ReadLikelihoodEngine::FlowBased(engine) => engine
                .compute_read_likelihoods(haplotypes, header, samples, reads)
                .map_err(|e| e.render()),
            ReadLikelihoodEngine::FlowBasedHmm(engine) => engine
                .compute_read_likelihoods(haplotypes, header, samples, reads)
                .map_err(|e| e.render()),
        }
    }

    /// Whether it is the `FlowBased` engine, which Mutect2's reference model asks.
    pub fn is_flow_based(&self) -> bool {
        matches!(self, ReadLikelihoodEngine::FlowBased(_))
    }

    /// `createLikelihoodCalculationEngine`: FlowBasedHMM, FlowBased or the PairHMM, as the
    /// arguments ask.
    pub fn create(
        likelihoods: &gatk_engine::pair_hmm_likelihood_engine::LikelihoodEngineArguments,
        flow: &Option<FlowEngineArguments>,
        flow_hmm: &Option<FlowHmmArguments>,
    ) -> Result<ReadLikelihoodEngine, (String, String)> {
        if let Some(hmm) = flow_hmm {
            let engine = FlowHmmEngine {
                arguments: hmm.clone(),
                log10_global_read_mismapping_rate: likelihoods.log10_global_read_mismapping_rate,
                expected_error_rate_per_base: likelihoods.expected_error_rate_per_base,
                dynamic_disqualification: likelihoods.dynamic_disqualification,
                read_disqualification_scale: likelihoods.read_disqualification_scale,
            };
            engine.check().map_err(|e| match e {
                FlowEngineError::Thrown(class, message) => (class.to_string(), message),
                other => ("GATKException".to_string(), format!("{other:?}")),
            })?;
            return Ok(ReadLikelihoodEngine::FlowBasedHmm(engine));
        }
        if let Some(flow) = flow {
            return Ok(ReadLikelihoodEngine::FlowBased(FlowLikelihoodEngine {
                arguments: flow.clone(),
                log10_global_read_mismapping_rate: likelihoods.log10_global_read_mismapping_rate,
                expected_error_rate_per_base: likelihoods.expected_error_rate_per_base,
                dynamic_disqualification: likelihoods.dynamic_disqualification,
                read_disqualification_scale: likelihoods.read_disqualification_scale,
            }));
        }
        gatk_engine::pair_hmm_likelihood_engine::PairHmmLikelihoodEngine::new(likelihoods.clone())
            .map(ReadLikelihoodEngine::PairHmm)
            .map_err(|e| ("IllegalArgumentException".to_string(), format!("{e:?}")))
    }
}

/// `FlowBasedHMMEngine`'s own arguments: the gap continuation penalty, the PCR model, DRAGstr's
/// table, and the three `--flow-hmm-engine-*` values.
#[derive(Debug, Clone)]
pub struct FlowHmmArguments {
    pub flow: FlowEngineArguments,
    pub constant_gcp: u8,
    pub pcr_error_model: gatk_engine::pair_hmm_likelihood_engine::PcrErrorModel,
    pub dragstr_params: Option<gatk_engine::dragstr::DragstrParams>,
    pub min_usable_indel_score: i32,
    pub flat_insertion_penalty: u8,
    pub flat_deletion_penalty: u8,
}

/// `FlowBasedHMMEngine`.
#[derive(Debug, Clone)]
pub struct FlowHmmEngine {
    pub arguments: FlowHmmArguments,
    pub log10_global_read_mismapping_rate: f64,
    pub expected_error_rate_per_base: f64,
    pub dynamic_disqualification: bool,
    pub read_disqualification_scale: f64,
}

/// `FlowBasedPairHMM.FLOW_SIZE`.
const FLOW_SIZE: usize = 4;

impl FlowHmmEngine {
    /// The constructor's check: the HMM steps four flows at a time.
    pub fn check(&self) -> Result<(), FlowEngineError> {
        let cycle = self.arguments.flow.flow_order_cycle_length;
        if cycle != 4 {
            return Err(FlowEngineError::Thrown(
                "org.broadinstitute.hellbender.exceptions.GATKException",
                format!(
                    "FlowBasedHMMEngine requires flow order of 4 elements but cycle length is \
                     specified as {cycle}"
                ),
            ));
        }
        Ok(())
    }

    /// `computeReadLikelihoods(haplotypes, header, samples, perSampleReadList, filterPoorly)`.
    pub fn compute_read_likelihoods(
        &self,
        haplotypes: &[Haplotype],
        header: &SamHeader,
        samples: &[String],
        reads: &[Vec<BamRecord>],
    ) -> Result<AlleleLikelihoods<BamRecord, Haplotype>, FlowEngineError> {
        let mut values = Vec::with_capacity(samples.len());
        for sample_reads in reads {
            values.push(self.sample_matrix(haplotypes, header, sample_reads)?);
        }
        let mut result = AlleleLikelihoods::new(
            SampleList::new(samples),
            AlleleList::new(haplotypes),
            reads.to_vec(),
            values,
        )
        .map_err(|e| FlowEngineError::Likelihoods(format!("{e:?}")))?;
        result
            .normalize_likelihoods(self.log10_global_read_mismapping_rate, true)
            .map_err(|e| FlowEngineError::Likelihoods(format!("{e:?}")))?;
        let expected = self.expected_error_rate_per_base;
        let dynamic = self.dynamic_disqualification;
        let scale = self.read_disqualification_scale;
        // Its own `log10MinTrueLikelihood`: three expected errors and two catastrophic ones at
        // least, at a fixed 0.001, whether capped or not.
        let floor = |length: usize| {
            let errors = 3.0f64.max((length as f64 * expected).ceil());
            let catastrophic = 2.0f64.max((length as f64 * 0.001).ceil());
            errors * jmath::math::log10(expected) + catastrophic * jmath::math::log10(0.001)
        };
        result
            .filter_poorly_modeled_evidence(|_, _, read: &BamRecord| {
                let maximum = floor(read.read_bases.len());
                if dynamic {
                    let threshold = log10_dynamic_read_qual_threshold(&read.base_qualities, scale);
                    if threshold < maximum {
                        threshold
                    } else {
                        maximum
                    }
                } else {
                    maximum
                }
            })
            .map_err(|e| FlowEngineError::Likelihoods(format!("{e:?}")))?;
        Ok(result)
    }

    /// `getErrorModelAdjustedQual` with the engine's floor, as the `byte` it is cast to.
    fn pcr_cache(&self) -> Option<Vec<i8>> {
        let model = self.arguments.pcr_error_model;
        if !model.has_rate_factor() {
            return None;
        }
        let max = gatk_engine::pair_hmm_likelihood_engine::MAX_REPEAT_LENGTH;
        Some(
            (0..=max)
                .map(|i| {
                    let penalty = jmath::strict_exp::exp(
                        f64::from(i) / (model.rate_factor() * std::f64::consts::PI),
                    );
                    let rounded = gatk_engine::pair_hmm::fast_round(40.0 - penalty + 1.0);
                    i32::from(self.arguments.min_usable_indel_score as i8).max(rounded) as i8
                })
                .collect(),
        )
    }

    fn sample_matrix(
        &self,
        haplotypes: &[Haplotype],
        header: &SamHeader,
        reads: &[BamRecord],
    ) -> Result<Vec<Vec<f64>>, FlowEngineError> {
        let cycle = self.arguments.flow.flow_order_cycle_length;
        let (flow_order, max_class) = match reads.first() {
            Some(first) => {
                let info = crate::flow_based_read::read_group_info(first, header)
                    .map_err(FlowEngineError::Read)?;
                (substring(&info.flow_order, cycle)?, info.max_class)
            }
            None => (first_usable_flow_order(header, cycle)?, 0),
        };
        let cache = self.pcr_cache();
        let min_usable = self.arguments.min_usable_indel_score as i8;
        let mut processed = Vec::with_capacity(reads.len());
        for read in reads {
            let mut flow = FlowRead::new(read, &flow_order, max_class, &self.arguments.flow.flow)
                .map_err(FlowEngineError::Read)?;
            flow.apply_alignment().map_err(FlowEngineError::Read)?;
            let bases = &read.read_bases;
            let (mut insertion, mut deletion, gcp): (Vec<i8>, Vec<i8>, Vec<i8>) = match &self
                .arguments
                .dragstr_params
            {
                Some(params) => {
                    let (gop, gcp) = gatk_engine::dragstr::dragstr_impute(params, bases)
                        .map_err(|e| FlowEngineError::Likelihoods(format!("{e:?}")))?;
                    let signed = |v: Vec<u8>| v.into_iter().map(|q| q as i8).collect::<Vec<_>>();
                    (signed(gop.clone()), signed(gop), signed(gcp))
                }
                None => (
                    vec![self.arguments.flat_insertion_penalty as i8; read.base_qualities.len()],
                    vec![self.arguments.flat_deletion_penalty as i8; read.base_qualities.len()],
                    vec![self.arguments.constant_gcp as i8; bases.len()],
                ),
            };
            // `applyPCRErrorModel`, compared as unsigned, then `capMinimumReadIndelQualities`,
            // compared as signed.
            if let Some(cache) = &cache {
                for i in 1..bases.len() {
                    let (_, repeat) =
                        gatk_engine::pair_hmm_likelihood_engine::find_tandem_repeat_units(
                            bases,
                            i - 1,
                        );
                    let cap = cache[repeat as usize] as u8;
                    insertion[i - 1] = (insertion[i - 1] as u8).min(cap) as i8;
                    deletion[i - 1] = (deletion[i - 1] as u8).min(cap) as i8;
                }
            }
            for i in 0..insertion.len() {
                insertion[i] = insertion[i].max(min_usable);
                deletion[i] = deletion[i].max(min_usable);
            }
            let to_key = |quals: &[i8], fill: i8| -> Result<Vec<u8>, FlowEngineError> {
                base_array_to_key_space(bases, flow.key.len(), quals, fill, &flow.flow_order)
            };
            let flow_insertion = to_key(&insertion, 40)?;
            let flow_deletion = to_key(&deletion, 40)?;
            let flow_gcp = to_key(&gcp, 10)?;
            processed.push((flow, flow_insertion, flow_deletion, flow_gcp));
        }
        let scored: Vec<FlowHaplotype> = haplotypes
            .iter()
            .map(|h| {
                FlowHaplotype::new(&h.bases(), &flow_order).ok_or_else(|| {
                    FlowEngineError::Thrown(
                        "org.broadinstitute.hellbender.exceptions.GATKException",
                        "baseArrayToKey".to_string(),
                    )
                })
            })
            .collect::<Result<_, _>>()?;
        let mut values = vec![vec![0.0; processed.len()]; haplotypes.len()];
        for (r, (read, ins, del, gcp)) in processed.iter().enumerate() {
            for (h, haplotype) in scored.iter().enumerate() {
                let read_first = read.flow_order.first().copied().unwrap_or(0);
                let start = haplotype
                    .flow_order
                    .iter()
                    .position(|b| *b == read_first)
                    .unwrap_or(0);
                values[h][r] = flow_pair_hmm(start, haplotype, read, ins, del, gcp);
            }
        }
        Ok(values)
    }
}

/// `FlowBasedKeyCodec.baseArrayToKeySpace`: each flow's smallest value over the bases it reads,
/// a flow reading none the last flow's value (the fill before any).
fn base_array_to_key_space(
    bases: &[u8],
    key_length: usize,
    values: &[i8],
    fill: i8,
    flow_order: &[u8],
) -> Result<Vec<u8>, FlowEngineError> {
    if bases.len() != values.len() {
        return Err(FlowEngineError::Thrown(
            "java.lang.IllegalArgumentException",
            "Read and qual arrays do not match".to_string(),
        ));
    }
    // `getFlowOrder()`: the first cycle of the read's flow order.
    let order = &flow_order[..flow_order.len().min(4)];
    let mut result = vec![0u8; key_length];
    let mut loc = 0;
    let mut flow_number = 0;
    let mut last = fill;
    while loc < bases.len() {
        let flow_base = order[flow_number % order.len()];
        let slot = result.get_mut(flow_number).ok_or_else(|| {
            FlowEngineError::Thrown(
                "java.lang.ArrayIndexOutOfBoundsException",
                format!("Index {flow_number} out of bounds for length {key_length}"),
            )
        })?;
        if bases[loc] != flow_base && bases[loc] != b'N' {
            *slot = last as u8;
        } else {
            let mut qual = i8::MAX;
            while loc < bases.len() && (bases[loc] == flow_base || bases[loc] == b'N') {
                qual = values[loc].min(qual);
                loc += 1;
            }
            *slot = qual as u8;
            last = qual;
        }
        flow_number += 1;
    }
    Ok(result)
}

/// `subComputekeayLikelihoodGivenHaplotypeKeysLog10`: the PairHMM over flows, a step of four
/// flows on either axis, a cell's prior the read's probability of the haplotype's hmer where the
/// two flows read the same base and zero elsewhere.
fn flow_pair_hmm(
    hap_start: usize,
    haplotype: &FlowHaplotype,
    read: &FlowRead,
    insertion: &[u8],
    deletion: &[u8],
    gcp: &[u8],
) -> f64 {
    use gatk_engine::pair_hmm::{
        initial_condition, initial_condition_log10, qual_to_trans_probs, DELETION_TO_DELETION,
        INDEL_TO_MATCH, INSERTION_TO_INSERTION, MATCH_TO_DELETION, MATCH_TO_INSERTION,
        MATCH_TO_MATCH,
    };
    let read_len = read.key.len() + 1 + FLOW_SIZE;
    let hap_len = haplotype.key.len() + 1 + FLOW_SIZE;
    let mut matches = vec![vec![0.0f64; hap_len]; read_len];
    let mut insertions = vec![vec![0.0f64; hap_len]; read_len];
    let mut deletions = vec![vec![0.0f64; hap_len]; read_len];
    let initial = initial_condition() / haplotype.key.len() as f64;
    for row in deletions.iter_mut().take(FLOW_SIZE + 1) {
        for cell in row.iter_mut() {
            *cell = initial;
        }
    }
    let mut transition = vec![[0.0f64; 6]; read_len];
    for i in 0..insertion.len() {
        transition[i + 1] = qual_to_trans_probs(insertion[i], deletion[i], gcp[i]);
    }
    let mut prior = vec![vec![0.0f64; hap_len]; read_len];
    for i in 0..read.key.len() {
        for j in 0..haplotype.key.len() {
            if haplotype.flow_order[j] == read.flow_order[i] {
                prior[i + 1 + FLOW_SIZE][j + 1 + FLOW_SIZE] = read.prob(i, haplotype.key[j]);
            }
        }
    }
    for i in 1 + FLOW_SIZE..read_len {
        let t = transition[i - FLOW_SIZE];
        let mut j = ((hap_start + i) % FLOW_SIZE) + FLOW_SIZE;
        while j < hap_len {
            matches[i][j] = prior[i][j]
                * (matches[i - 1][j - 1] * t[MATCH_TO_MATCH]
                    + insertions[i - 1][j - 1] * t[INDEL_TO_MATCH]
                    + deletions[i - 1][j - 1] * t[INDEL_TO_MATCH]);
            insertions[i][j] = matches[i - FLOW_SIZE][j] * t[MATCH_TO_INSERTION]
                + insertions[i - FLOW_SIZE][j] * t[INSERTION_TO_INSERTION];
            deletions[i][j] = matches[i][j - FLOW_SIZE] * t[MATCH_TO_DELETION]
                + deletions[i][j - FLOW_SIZE] * t[DELETION_TO_DELETION];
            j += FLOW_SIZE;
        }
    }
    let end = read_len - 1;
    let mut total = 0.0;
    for j in 1..hap_len {
        total += matches[end][j] + insertions[end][j];
    }
    jmath::math::log10(total) - initial_condition_log10()
}
