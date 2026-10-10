//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine`
//! (`isActive`) and `ReferenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny` (GATK 4.6.2.0):
//! how likely a locus is to hold a variant, which the activity profile turns into assembly
//! regions.
//!
//! Each sample's pileup is scored as reference against anything else: a base at or under the
//! minimum base quality is skipped (a deletion never is, and is scored at the reference-model
//! deletion quality), and a base is "alternate" if it differs from the reference or sits in, next
//! to, or beside an indel or a soft clip. The likelihoods are then stored as PLs, so what the rest
//! of the computation reads is **the PL-rounded likelihoods**, not the ones just computed. One
//! sample's state is the Dirichlet posterior that it is not hom-ref; several samples' is the QUAL
//! of a genotyping engine run on a fake `N`/`<FAKE_ALT>` site at a calling confidence of 4, as a
//! probability. Alternate bases next to a soft clip add the read's soft-clipped bases of quality
//! above 28 to a running average, and an average above 6 marks the state as high-quality soft
//! clips with that average as its value.
//!
//! The force-calling and pileup-detection branches are not ported.

use gatk_engine::allele_frequency_calculator::{AlleleFrequencyCalculator, Priors};
use gatk_engine::locus_iterator::AlignmentContext;
use gatk_engine::math_utils::{max_element_index, qual_to_prob};
use gatk_engine::pair_hmm::approximate_log10_sum_log10;
use gatk_engine::pileup::PileupElement;
use gatk_engine::qual_quantizer::qual_to_prob_log10;
use gatk_engine::read_pileup::ReadPileup;
use htsjdk_bam::cigar::Op;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, VariantContext};

use crate::calculate_genotype_posteriors::gls_to_pls;
use crate::genotyping_engine::{Configuration, EngineError, GenotypingEngine, SubsetMethod};

/// `HaplotypeCallerEngine.AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD`.
const AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD: f64 = 6.0;
/// `HaplotypeCallerEngine.MAXMIN_CONFIDENCE_FOR_CONSIDERING_A_SITE_AS_POSSIBLE_VARIANT_IN_ACTIVE_REGION_DISCOVERY`.
const MAXMIN_CONFIDENCE_FOR_ACTIVE_REGION_DISCOVERY: f64 = 4.0;
/// `HaplotypeCallerEngine.MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY`.
const MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY: usize = 2;
/// `ReferenceConfidenceModel.HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD`.
const HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD: u8 = 28;
/// `ReferenceConfidenceModel.REF_MODEL_DELETION_QUAL`.
pub const REF_MODEL_DELETION_QUAL: u8 = 30;

/// `ActivityProfileState.Type`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActivityType {
    None,
    HighQualitySoftClips,
}

impl ActivityType {
    pub fn name(self) -> &'static str {
        match self {
            ActivityType::None => "NONE",
            ActivityType::HighQualitySoftClips => "HIGH_QUALITY_SOFT_CLIPS",
        }
    }
}

/// `ActivityProfileState`, with the original active probability `isActive` sets on it.
#[derive(Debug, Clone, PartialEq)]
pub struct ActivityState {
    pub prob: f64,
    pub kind: ActivityType,
    /// `getResultValue()`, null for a locus with no pileup.
    pub result_value: Option<f64>,
    pub original_active_prob: f64,
}

/// The `HaplotypeCallerArgumentCollection` fields `isActive` reads.
#[derive(Debug, Clone)]
pub struct IsActiveArguments {
    pub sample_ploidy: usize,
    pub snp_heterozygosity: f64,
    pub indel_heterozygosity: f64,
    pub heterozygosity_standard_deviation: f64,
    pub standard_confidence_for_calling: f64,
    pub max_alternate_alleles: usize,
    pub min_base_quality_score: u8,
    pub ref_model_deletion_quality: u8,
}

impl Default for IsActiveArguments {
    fn default() -> Self {
        IsActiveArguments {
            sample_ploidy: 2,
            snp_heterozygosity: 0.001,
            indel_heterozygosity: 1.25e-4,
            heterozygosity_standard_deviation: 0.01,
            standard_confidence_for_calling: 30.0,
            max_alternate_alleles: 6,
            min_base_quality_score: 10,
            ref_model_deletion_quality: REF_MODEL_DELETION_QUAL,
        }
    }
}

/// `MathUtils.RunningAverage`: Welford's running mean.
#[derive(Debug, Clone, Default)]
pub struct RunningAverage {
    mean: f64,
    count: u64,
}

impl RunningAverage {
    pub fn add(&mut self, observation: f64) {
        self.count += 1;
        self.mean += (observation - self.mean) / self.count as f64;
    }

    pub fn mean(&self) -> f64 {
        self.mean
    }
}

/// The active-region evaluator: `isActive` with the engine it builds for it.
pub struct ActiveRegionEvaluator {
    arguments: IsActiveArguments,
    ploidy: usize,
    calculator: AlleleFrequencyCalculator,
    engine: GenotypingEngine,
    samples: Vec<String>,
}

impl ActiveRegionEvaluator {
    /// `initializeActiveRegionEvaluationGenotyperEngine`: the caller's arguments with a ploidy of
    /// at least two, a calling confidence of at most 4, and variant sites only.
    pub fn new(arguments: IsActiveArguments, samples: Vec<String>) -> Self {
        let ploidy = arguments
            .sample_ploidy
            .max(MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY);
        let priors = Priors {
            snp_heterozygosity: arguments.snp_heterozygosity,
            indel_heterozygosity: arguments.indel_heterozygosity,
            heterozygosity_standard_deviation: arguments.heterozygosity_standard_deviation,
            sample_ploidy: ploidy,
        };
        let calculator = AlleleFrequencyCalculator::make_calculator(&priors);
        let configuration = Configuration {
            standard_confidence_for_calling: MAXMIN_CONFIDENCE_FOR_ACTIVE_REGION_DISCOVERY
                .min(arguments.standard_confidence_for_calling),
            max_alternate_alleles: arguments.max_alternate_alleles,
            sample_ploidy: ploidy,
            annotate_number_of_alleles_discovered: false,
            emit_all_active_sites: false,
            allele_specific: false,
            emit_all_confident_sites: false,
            annotate_all_sites_with_pls: false,
            force_keep_all_alleles: false,
            assignment_method: SubsetMethod::UsePlsToAssign,
        };
        let engine = GenotypingEngine::new(
            configuration,
            AlleleFrequencyCalculator::make_calculator(&priors),
        );
        ActiveRegionEvaluator {
            arguments,
            ploidy,
            calculator,
            engine,
            samples,
        }
    }

    /// `isActive(context, ref, features)`, `ref_base` being the reference base at the locus.
    pub fn is_active(
        &mut self,
        context: &AlignmentContext<'_>,
        ref_base: u8,
        header: &SamHeader,
    ) -> Result<ActivityState, EngineError> {
        if context.pileup.is_empty() {
            return Ok(ActivityState {
                prob: 0.0,
                kind: ActivityType::None,
                result_value: None,
                original_active_prob: 0.0,
            });
        }
        // `splitContextBySampleName`: the one sample takes the whole pileup; otherwise each sample
        // with reads, in the pileup's sample order.
        let split: Vec<(String, ReadPileup<'_>)> = if self.samples.len() == 1 {
            vec![(
                self.samples[0].clone(),
                ReadPileup::new(
                    &context.contig,
                    context.position,
                    context.pileup.elements.clone(),
                ),
            )]
        } else {
            context
                .pileup
                .split_by_sample(header, None)
                .map_err(|message| EngineError::Runtime {
                    class: "UserException.ReadMissingReadGroup".to_string(),
                    message,
                })?
                .into_iter()
                .filter(|(_, pileup)| !pileup.is_empty())
                .collect()
        };

        let mut soft_clips = RunningAverage::default();
        let mut genotypes = Vec::with_capacity(split.len());
        for (sample, pileup) in &split {
            let likelihoods = self.ref_vs_any_likelihoods(pileup, ref_base, &mut soft_clips);
            let mut genotype = Genotype::new(sample, vec![Allele::no_call(); self.ploidy]);
            genotype.pl = Some(gls_to_pls(&likelihoods));
            genotypes.push(genotype);
        }

        // What `genotypes.get(0).getLikelihoods().getAsVector()` gives back: the PLs over -10.
        let first_gls = gls_of(&genotypes[0]);
        let prob = if genotypes.len() == 1 {
            self.calculator
                .calculate_single_sample_biallelic_non_ref_posterior(&first_gls, true)
        } else {
            let alleles = vec![
                Allele::create(b"N", true).expect("the fake reference"),
                Allele::create(b"<FAKE_ALT>", false).expect("the fake alternate"),
            ];
            let mut vc = VariantContext::new(&context.contig, i64::from(context.position), alleles);
            vc.genotypes = GenotypesContext::new(genotypes.clone());
            match self.engine.calculate_genotypes(&vc)? {
                None => 0.0,
                Some(out) => qual_to_prob(-10.0 * out.log10_p_error),
            }
        };
        let mean = soft_clips.mean();
        let kind = if mean > AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD {
            ActivityType::HighQualitySoftClips
        } else {
            ActivityType::None
        };
        let max = max_element_index(&first_gls, 0, first_gls.len());
        Ok(ActivityState {
            prob,
            kind,
            result_value: Some(mean),
            original_active_prob: first_gls[max] - first_gls[0],
        })
    }

    /// `calcGenotypeLikelihoodsOfRefVsAny(ploidy, pileup, refBase, minBaseQual, hqSoftClips,
    /// false)`.
    fn ref_vs_any_likelihoods(
        &self,
        pileup: &ReadPileup<'_>,
        ref_base: u8,
        soft_clips: &mut RunningAverage,
    ) -> Vec<f64> {
        let count = self.ploidy + 1;
        let log10_ploidy = jmath::math::log10(self.ploidy as f64);
        let log10_one_third = -jmath::math::log10(3.0);
        let mut likelihoods = vec![0.0; count];
        let mut read_count = 0usize;
        for element in &pileup.elements {
            let qual = if element.is_deletion() {
                self.arguments.ref_model_deletion_quality
            } else {
                element.qual()
            };
            if (qual as i8) <= (self.arguments.min_base_quality_score as i8)
                && !element.is_deletion()
            {
                continue;
            }
            read_count += 1;
            let is_alt = is_alt_before_assembly(element, ref_base);
            let qual_to_error_log10 = f64::from(qual) * -0.1;
            let (reference, non_ref) = if is_alt {
                (
                    qual_to_error_log10 + log10_one_third,
                    qual_to_prob_log10(qual),
                )
            } else {
                (
                    qual_to_prob_log10(qual),
                    qual_to_error_log10 + log10_one_third,
                )
            };
            likelihoods[0] += reference + log10_ploidy;
            likelihoods[count - 1] += non_ref + log10_ploidy;
            let mut j = count as i64 - 2;
            for (i, slot) in likelihoods.iter_mut().enumerate().take(count - 1).skip(1) {
                *slot += approximate_log10_sum_log10(
                    reference + jmath::math::log10(j as f64),
                    non_ref + jmath::math::log10(i as f64),
                );
                j -= 1;
            }
            if is_alt && element.is_next_to_soft_clip() {
                soft_clips.add(f64::from(count_high_quality_soft_clips(
                    element.read,
                    HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD,
                )));
            }
        }
        let denominator = read_count as f64 * log10_ploidy;
        for value in &mut likelihoods {
            *value -= denominator;
        }
        likelihoods
    }
}

/// `GenotypeLikelihoods.fromPLs(pls).getAsVector()`.
fn gls_of(genotype: &Genotype) -> Vec<f64> {
    genotype
        .pl
        .as_ref()
        .expect("a PL")
        .iter()
        .map(|&pl| f64::from(pl) / -10.0)
        .collect()
}

/// `ReferenceConfidenceModel.isAltBeforeAssembly`.
fn is_alt_before_assembly(element: &PileupElement<'_>, ref_base: u8) -> bool {
    element.base() != ref_base
        || element.is_deletion()
        || element.is_before_deletion_start()
        || element.is_after_deletion_end()
        || element.is_before_insertion()
        || element.is_after_insertion()
        || element.is_next_to_soft_clip()
}

/// `AlignmentUtils.countHighQualitySoftClips`: soft-clipped bases of quality above the threshold.
pub fn count_high_quality_soft_clips(read: &BamRecord, threshold: u8) -> i32 {
    let mut count = 0;
    let mut position = 0usize;
    for element in &read.cigar.elements {
        let length = element.length as usize;
        if element.op == Op::S {
            for _ in 0..length {
                if read.base_qualities.get(position).copied().unwrap_or(0) > threshold {
                    count += 1;
                }
                position += 1;
            }
        } else if element.op.consumes_read_bases() {
            position += length;
        }
    }
    count
}
