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
                .map_err(|e| format!("{e:?}")),
        }
    }

    /// Whether it is one of the flow-based engines, which the reference model asks.
    pub fn is_flow_based(&self) -> bool {
        matches!(self, ReadLikelihoodEngine::FlowBased(_))
    }
}
