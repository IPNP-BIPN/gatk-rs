//! Ported from `org.broadinstitute.hellbender.tools.walkers.genotyper.DRAGENGenotypesModel`,
//! `GenotypeLikelihoodCalculatorDRAGEN` and `FRDBQDUtils`, and from
//! `org.broadinstitute.hellbender.transformers.DRAGENMappingQualityReadTransformer` (GATK 4.6.2.0):
//! the genotype models HaplotypeCaller switches on with `--apply-bqd` and `--apply-frd`.
//!
//! Each sample's genotype likelihoods start as the standard model's
//! ([`compute_log10_genotype_likelihoods`]). Each enabled model then computes a likelihood for every
//! **homozygous** genotype and the sample keeps, genotype by genotype, the larger of the two:
//!
//!  * base quality dropout (BQD) asks whether a homozygous call's dissent is an error that set in
//!    along the reads, so it sorts each strand's reads by how far into the read the site falls and
//!    looks for the best split: before it the reads are a mixture of the error allele and the call,
//!    after it they are the call alone, and the split pays the mean base quality before it, less a
//!    homopolymer adjustment. Indels and the error allele equal to the call are skipped;
//!  * foreign read detection (FRD) asks whether the dissent comes from reads placed here by
//!    mistake, so it tries every distinct mapping-quality threshold the reads give, lets only the
//!    reads at or below it support the foreign allele, and keeps the best of the forward, the
//!    reverse and both strands.
//!
//! The reads are each sample's reads in the likelihood matrix plus those the HMM disqualified,
//! which BQD counts for the split and the mean base quality but FRD only for its thresholds.
//!
//! What a reader would not guess:
//!
//!  * a read's position at the site is taken from its **soft** start, so a soft clip covering the
//!    site gives the clipped base's quality, and a deletion gives the quality of the base after it;
//!  * the "feather ends" the reads are sorted by count the hard clips before (forward) or after
//!    (reverse) the read, plus the offset, and a read not covering the site sorts with offset -1;
//!  * FRD's thresholds are the reads' mapping qualities mapped through DRAGEN's piecewise-linear
//!    table (0, 30, 40, 45, 50, 50 at 0, 30, 60, 100, 200 and 256), times -0.1, plus the allele's
//!    prior: none for the reference, the flat SNP prior 34.77 for a SNP, and the DRAGstr API (or
//!    the same flat prior without DRAGstr) for an indel, also times -0.1;
//!  * the depth cap (`--max-effective-depth-adjustment-for-frd`) blends the FRD score towards the
//!    standard model's likelihood for the het of the call and the foreign allele, read at the
//!    **diploid** index of that pair whatever the ploidy, less `log10(2)`;
//!  * FRD sees the likelihoods BQD has already raised.

use crate::genotyper_debug as debug;
use crate::tsv_table::java_double_to_string;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_vcf::allele::Allele;

use crate::dragstr::{DragstrError, DragstrParams, DragstrReferenceAnalyzer, StrLookupError};
use crate::genotype_index::{allele_counts_to_index, genotype_count};
use crate::genotype_likelihood_calculator::compute_log10_genotype_likelihoods;
use crate::math_utils::{java_max, java_min, pow10};
use crate::pair_hmm::approximate_log10_sum_log10;
use crate::{read, read_utils};

/// `DRAGENGenotypesModel.FLAT_SNP_HET_PRIOR`.
pub const FLAT_SNP_HET_PRIOR: f64 = 34.77;
/// `DRAGENGenotypesModel.BQD_HOMOPOLYMER_PHRED_ADJUSTMENT_FACTOR`.
pub const BQD_HOMOPOLYMER_PHRED_ADJUSTMENT_FACTOR: f64 = 5.0;
/// `GenotypeLikelihoodCalculatorDRAGEN.BQD_FIXED_ERROR_RATE`.
const BQD_FIXED_ERROR_RATE: f64 = 0.5;
/// `GenotypeLikelihoodCalculatorDRAGEN.PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE`.
const PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE: f64 = 2.5;

/// `DRAGENMappingQualityReadTransformer.EXTENDED_MAPPING_QUALITY_READ_TAG`.
pub const EXTENDED_MAPPING_QUALITY_READ_TAG: [u8; 2] = *b"XQ";

const MQ_TABLE_X: [i32; 6] = [0, 30, 60, 100, 200, 256];
const MQ_TABLE_Y: [i32; 6] = [0, 30, 40, 45, 50, 50];

/// `mapMappingQualityToPhredLikelihoodScore`: DRAGEN's piecewise-linear map from a mapping quality
/// to the phred score of a mismapping. Past 256 it is the reference's `GATKException`, whose
/// message is returned.
pub fn map_mapping_quality_to_phred_likelihood_score(value: i32) -> Result<f64, String> {
    for i in 1..MQ_TABLE_X.len() {
        if value <= MQ_TABLE_X[i] {
            let xfactor =
                f64::from(value - MQ_TABLE_X[i - 1]) / f64::from(MQ_TABLE_X[i] - MQ_TABLE_X[i - 1]);
            return Ok(f64::from(MQ_TABLE_Y[i - 1])
                + (xfactor * f64::from(MQ_TABLE_Y[i] - MQ_TABLE_Y[i - 1])));
        }
    }
    Err(format!(
        "Something went wrong trying to map an an invalid XQ tag '{value}' val must fall between 0 and 250"
    ))
}

/// `DRAGENMappingQualityReadTransformer.apply`: a read carrying `XQ` takes it as its mapping
/// quality. A value outside 0 to 255 is `setMappingQuality`'s `IllegalArgumentException`, and a
/// value that is not an integer `getAttributeAsInteger`'s type mismatch; both are returned as
/// their messages.
pub fn transform_mapping_quality(record: &mut BamRecord) -> Result<(), String> {
    let extended = match record
        .tags
        .get(Tag::new(&EXTENDED_MAPPING_QUALITY_READ_TAG))
    {
        None => return Ok(()),
        Some(TagValue::Int(value)) => *value,
        // `Integer.parseInt(attributeValue.toString())`, which does not trim.
        Some(TagValue::Str(text)) => i64::from(
            text.parse::<i32>()
                .map_err(|_| "Attribute XQ not of (or convertible to) type integer".to_string())?,
        ),
        Some(TagValue::Char(c)) => {
            i64::from((*c as char).to_digit(10).ok_or_else(|| {
                "Attribute XQ not of (or convertible to) type integer".to_string()
            })?)
        }
        Some(_) => return Err("Attribute XQ not of (or convertible to) type integer".to_string()),
    };
    if !(0..=255).contains(&extended) {
        return Err(format!(
            "mapping quality must be >= 0 and <= 255 \nMappingQuality: {extended}"
        ));
    }
    record.mapping_quality = extended as u8;
    Ok(())
}

/// `FRDBQDUtils.computeForwardHomopolymerAdjustment`: five per base of the error base's run, up to
/// four, ending just before the site.
pub fn forward_homopolymer_adjustment(padded_reference: &[u8], offset: i32, error_base: u8) -> f64 {
    let mut length = 0;
    let mut position = offset - 1;
    while position >= 0 && length < 4 {
        if error_base != padded_reference[position as usize] {
            break;
        }
        length += 1;
        position -= 1;
    }
    BQD_HOMOPOLYMER_PHRED_ADJUSTMENT_FACTOR * f64::from(length)
}

/// `FRDBQDUtils.computeReverseHomopolymerAdjustment`: the same, for the run starting just after
/// the site.
pub fn reverse_homopolymer_adjustment(padded_reference: &[u8], offset: i32, error_base: u8) -> f64 {
    let mut length = 0;
    let mut position = offset + 1;
    while (position as usize) < padded_reference.len() && length < 4 {
        if error_base != padded_reference[position as usize] {
            break;
        }
        length += 1;
        position += 1;
    }
    BQD_HOMOPOLYMER_PHRED_ADJUSTMENT_FACTOR * f64::from(length)
}

/// The prior the reference's model reads for an indel alternate allele: DRAGstr's API at the
/// period and repeat count one base past the event's offset, or the flat SNP prior without DRAGstr.
pub fn indel_prior(
    dragstrs: Option<&DragstrReferenceAnalyzer>,
    params: &DragstrParams,
    offset_for_ref_into_event: i32,
) -> Result<f64, IndelPriorError> {
    match dragstrs {
        None => {
            debug::println("No API from DRAGStrs found, falling back on snp het prior for indels");
            Ok(FLAT_SNP_HET_PRIOR)
        }
        Some(strs) => {
            let period = strs
                .period(offset_for_ref_into_event + 1)
                .map_err(IndelPriorError::Lookup)?;
            let repeats = strs
                .repeat_length(offset_for_ref_into_event + 1)
                .map_err(IndelPriorError::Lookup)?;
            let api = params
                .api(period, repeats)
                .map_err(IndelPriorError::Params)?;
            debug::println(&format!(
                "API found: {} with period used: {period}  and repeats: {repeats}",
                java_double_to_string(api)
            ));
            Ok(api)
        }
    }
}

/// Why [`indel_prior`] failed.
#[derive(Debug)]
pub enum IndelPriorError {
    Lookup(StrLookupError),
    Params(DragstrError),
}

/// `DRAGENGenotypesModel.DragenReadContainer`: what the models read of one read at the site.
#[derive(Debug, Clone)]
struct DragenRead {
    reverse: bool,
    /// The read's quality at the site, `None` where `getReadIndexForReferenceCoordinate` found no
    /// offset.
    base_quality: Option<i32>,
    forward_feather_end: i32,
    reverse_feather_end: i32,
    /// The read's column in the sample's likelihood matrix, `None` for a read the HMM disqualified.
    index: Option<usize>,
    phred_mapping_quality: Result<f64, String>,
    /// `phredPFValue`, set by FRD for each foreign allele in turn.
    phred_pf: f64,
}

impl DragenRead {
    fn new(record: &BamRecord, variant_offset: i32, index: Option<usize>) -> DragenRead {
        let offset = read_utils::read_index_for_read(record, variant_offset).0;
        let base_quality = (offset != read_utils::READ_INDEX_NOT_FOUND)
            .then(|| i32::from(record.base_qualities[offset as usize]));
        let length = read::length(record) as i32;
        DragenRead {
            reverse: read::is_reverse_strand(record),
            base_quality,
            forward_feather_end: (read_utils::soft_start(record)
                - read_utils::unclipped_start(record))
                + offset,
            reverse_feather_end: (read_utils::unclipped_end(record) - read_utils::soft_end(record))
                + (length - offset),
            index,
            phred_mapping_quality: map_mapping_quality_to_phred_likelihood_score(i32::from(
                record.mapping_quality,
            )),
            phred_pf: 0.0,
        }
    }

    fn quality_or_zero(&self) -> i32 {
        self.base_quality.unwrap_or(0)
    }
}

/// One sample's input to the model.
pub struct DragenSample<'a> {
    /// The sample's read-by-allele log10 likelihoods, `[allele][read]`, in the genotyping alleles'
    /// order.
    pub likelihoods: &'a [Vec<f64>],
    /// The reads of the matrix, in its order.
    pub reads: &'a [BamRecord],
    /// The reads the HMM disqualified (`filteredSampleEvidence`).
    pub filtered: &'a [BamRecord],
}

/// `DRAGENGenotypesModel`.
#[derive(Debug, Clone)]
pub struct DragenGenotypesModel {
    compute_bqd: bool,
    compute_frd: bool,
    allele_padding: i32,
    max_effective_depth_adjustment: i32,
}

impl DragenGenotypesModel {
    /// The constructor, which refuses a model with both models off.
    pub fn new(
        use_bqd: bool,
        use_frd: bool,
        allele_padding: i32,
        max_effective_depth_adjustment: i32,
    ) -> Result<Self, String> {
        if !(use_bqd || use_frd) {
            return Err("DRAGENGenotypesModel is intended for computing BQD/FRD adjustments but both BQD and FRD are disabled, use IndependentSamplesGenotypesModel instead".to_string());
        }
        Ok(DragenGenotypesModel {
            compute_bqd: use_bqd,
            compute_frd: use_frd,
            allele_padding,
            max_effective_depth_adjustment,
        })
    }

    /// `calculateLikelihoods` for one sample: its log10 genotype likelihoods, in canonical order.
    ///
    /// `variant_calling_subset_start` is the start of the interval `setVariantCallingSubsetUsed`
    /// recorded, `api` the indel prior from [`indel_prior`].
    #[allow(clippy::too_many_arguments)]
    pub fn sample_log10_likelihoods(
        &self,
        ploidy: usize,
        alleles: &[Allele],
        sample: &DragenSample,
        variant_calling_subset_start: i32,
        padded_reference: &[u8],
        offset_for_ref_into_event: i32,
        api: f64,
    ) -> Result<Vec<f64>, String> {
        let variant_offset = variant_calling_subset_start + self.allele_padding;
        let mut forward = Vec::new();
        let mut reverse = Vec::new();
        for (index, record) in sample.reads.iter().enumerate() {
            let container = DragenRead::new(record, variant_offset, Some(index));
            if container.reverse {
                &mut reverse
            } else {
                &mut forward
            }
            .push(container);
        }
        for record in sample.filtered {
            let container = DragenRead::new(record, variant_offset, None);
            if container.reverse {
                &mut reverse
            } else {
                &mut forward
            }
            .push(container);
        }
        forward.sort_by(|a, b| {
            b.forward_feather_end
                .cmp(&a.forward_feather_end)
                .then(a.quality_or_zero().cmp(&b.quality_or_zero()))
        });
        reverse.sort_by(|a, b| {
            b.reverse_feather_end
                .cmp(&a.reverse_feather_end)
                .then(a.quality_or_zero().cmp(&b.quality_or_zero()))
        });

        let mut likelihoods = compute_log10_genotype_likelihoods(ploidy, sample.likelihoods);
        debug::println("\n Standard Genotyping Likelihoods Results:");
        debug::println(&java_doubles(&likelihoods));
        if self.compute_bqd {
            let bqd = bqd_likelihoods(
                ploidy,
                alleles,
                sample.likelihoods,
                &forward,
                &reverse,
                padded_reference,
                offset_for_ref_into_event,
            )?;
            apply_adjustment(&mut likelihoods, "BQD", &bqd);
        }
        if self.compute_frd {
            let mut all = forward;
            all.extend(reverse);
            let frd = frd_likelihoods(
                ploidy,
                alleles,
                sample.likelihoods,
                &likelihoods,
                &mut all,
                FLAT_SNP_HET_PRIOR,
                api,
                self.max_effective_depth_adjustment,
            )?;
            apply_adjustment(&mut likelihoods, "FRD", &frd);
        }
        debug::println("merged matrix:");
        debug::println(&java_doubles(&likelihoods));
        Ok(likelihoods)
    }
}

/// `applyLikelihoodsAdjusmentToBaseline`: the larger of the two, genotype by genotype.
fn apply_adjustment(initial: &mut [f64], name: &str, adjustment: &[f64]) {
    debug::println(&format!("{name} adjusted likelihoods:"));
    debug::println(&java_doubles(adjustment));
    for (value, adjusted) in initial.iter_mut().zip(adjustment) {
        *value = java_max(*value, *adjusted);
    }
}

/// `Allele.length()`, which is zero for a symbolic allele.
fn java_length(allele: &Allele) -> usize {
    if allele.is_symbolic() {
        0
    } else {
        allele.len()
    }
}

/// The index of the homozygous genotype of an allele.
fn homozygous_index(allele: usize, ploidy: usize) -> usize {
    allele_counts_to_index(&[allele, ploidy]).expect("an even-length allele count array")
}

/// `calculateBQDLikelihoods`.
fn bqd_likelihoods(
    ploidy: usize,
    alleles: &[Allele],
    likelihoods: &[Vec<f64>],
    forward: &[DragenRead],
    reverse: &[DragenRead],
    padded_reference: &[u8],
    offset_for_ref_into_event: i32,
) -> Result<Vec<f64>, String> {
    let count = genotype_count(ploidy, alleles.len()).map_err(|e| e.message())?;
    let mut output = vec![f64::NEG_INFINITY; count];
    let ref_length = java_length(&alleles[0]);
    for gt in 0..alleles.len() {
        let index = homozygous_index(gt, ploidy);
        for error in 0..alleles.len() {
            if gt == error
                || java_length(&alleles[gt]) != ref_length
                || java_length(&alleles[error]) != ref_length
            {
                continue;
            }
            let error_base = alleles[error].display_string().as_bytes()[0];
            let forward_adjustment = forward_homopolymer_adjustment(
                padded_reference,
                offset_for_ref_into_event,
                error_base,
            );
            let reverse_adjustment = reverse_homopolymer_adjustment(
                padded_reference,
                offset_for_ref_into_event,
                error_base,
            );
            let forward_score =
                bqd_for_strand(likelihoods, forward, forward_adjustment, true, gt, error);
            let reverse_score =
                bqd_for_strand(likelihoods, reverse, reverse_adjustment, false, gt, error);
            let score = (forward_score + reverse_score) * -0.1;
            output[index] = java_max(output[index], score);
        }
    }
    Ok(output)
}

/// `computeBQDModelForStrandData`: the phred score of the best split of one strand's sorted reads.
fn bqd_for_strand(
    likelihoods: &[Vec<f64>],
    reads: &[DragenRead],
    homopolymer_adjustment: f64,
    forwards: bool,
    homozygous: usize,
    error: usize,
) -> f64 {
    if reads.is_empty() {
        return 0.0;
    }
    debug::println(&format!(
        "errorAllele index: {error} theta: {} homopolymerAdjustment: {}",
        if forwards { "1" } else { "2" },
        java_double_to_string(homopolymer_adjustment)
    ));
    let log10_error = jmath::math::log10(BQD_FIXED_ERROR_RATE);
    let log10_non_error = jmath::math::log10(1.0 - BQD_FIXED_ERROR_RATE);
    let size = reads.len();
    let mut error_allele = vec![0.0; size + 1];
    let mut mean_quality = vec![0.0; size + 1];
    let mut genotype = vec![0.0; size + 1];
    let mut total_quality = 0.0;
    let mut quality_denominator = 0i32;
    for i in 1..=size {
        let container = &reads[i - 1];
        let (homozygous_contribution, error_contribution) = match container.index {
            Some(read) => (likelihoods[homozygous][read], likelihoods[error][read]),
            None => (0.0, 0.0),
        };
        let phred = if homozygous_contribution == 0.0 && error_contribution == 0.0 {
            0.0
        } else {
            -10.0
                * approximate_log10_sum_log10(
                    error_contribution + log10_error,
                    homozygous_contribution + log10_non_error,
                )
        };
        error_allele[i] = error_allele[i - 1] + phred;
        genotype[i] = genotype[i - 1] + -10.0 * homozygous_contribution;
        if let Some(quality) = container.base_quality {
            total_quality += f64::from(quality);
            quality_denominator += 1;
        }
        let denominator = if quality_denominator == 0 {
            1
        } else {
            quality_denominator
        };
        mean_quality[i] = java_max(
            0.0,
            ((total_quality / f64::from(denominator)) * PHRED_SCALED_ADJUSTMENT_FOR_BQ_SCORE)
                - homopolymer_adjustment,
        );
    }
    let mut min_score = f64::INFINITY;
    let mut n_used = 0;
    for n in 0..mean_quality.len() {
        let remaining = genotype[genotype.len() - 1] - genotype[n];
        let score = mean_quality[n] + error_allele[n] + remaining;
        if debug::is_enabled() {
            debug::println(&format!(
                "n={n}: {}, cum_phred_bq={}, cum_prob_r_Error={}, prob_G_remaining={}",
                fixed(score, 0, 2),
                fixed(mean_quality[n], 0, 2),
                fixed(error_allele[n], 0, 2),
                fixed(remaining, 0, 2)
            ));
        }
        if min_score > score {
            min_score = score;
            n_used = n;
        }
    }
    // The reference prints theta as 1 or 0 and the n's suffix as 1 or 2.
    debug::println(&format!(
        "theta={} n{}={n_used:2}, best_phred_score ={} q_mean={}, alpha={}, Ph(E)={};  ",
        i32::from(forwards),
        if forwards { 1 } else { 2 },
        fixed(min_score, 5, 2),
        fixed(mean_quality[n_used], 5, 2),
        fixed(BQD_FIXED_ERROR_RATE, 4, 2),
        fixed(error_allele[n_used], 4, 2)
    ));
    min_score
}

/// `calculateFRDLikelihoods`.
#[allow(clippy::too_many_arguments)]
fn frd_likelihoods(
    ploidy: usize,
    alleles: &[Allele],
    likelihoods: &[Vec<f64>],
    ploidy_model_likelihoods: &[f64],
    reads: &mut [DragenRead],
    snp_apriori_het: f64,
    indel_apriori_het: f64,
    max_effective_depth: i32,
) -> Result<Vec<f64>, String> {
    let count = genotype_count(ploidy, alleles.len()).map_err(|e| e.message())?;
    let mut output = vec![f64::NEG_INFINITY; count];
    let ref_length = java_length(&alleles[0]);
    let depth = likelihoods.first().map_or(0, Vec::len) as i32;
    let log10_one_half = jmath::math::log10(0.5);
    for foreign in 0..alleles.len() {
        let is_indel = java_length(&alleles[foreign]) != ref_length;
        let prior = if foreign == 0 {
            0.0
        } else {
            (if is_indel {
                indel_apriori_het
            } else {
                snp_apriori_het
            }) * -0.1
        };
        let thresholds = critical_values(reads, prior)?;
        if debug::is_enabled() {
            let values: Vec<String> = thresholds
                .iter()
                .map(|t| java_double_to_string(*t))
                .collect();
            debug::println(&format!(
                "fIndex: {foreign} criticalValues: \n{}",
                values.join("\n")
            ));
        }
        for gt in 0..alleles.len() {
            if gt == foreign {
                continue;
            }
            let index = homozygous_index(gt, ploidy);
            debug::println(&format!("indexForGT {index}"));
            debug::println("\nForwards Strands: ");
            let forward =
                frd_for_strand(likelihoods, gt, foreign, reads, |c| !c.reverse, &thresholds);
            debug::println("\nReverse Strands: ");
            let reverse =
                frd_for_strand(likelihoods, gt, foreign, reads, |c| c.reverse, &thresholds);
            debug::println("\nBoth Strands: ");
            let both = frd_for_strand(likelihoods, gt, foreign, reads, |_| true, &thresholds);
            // The reference prints the three arrays' identities, which no run reproduces.
            debug::println(&format!(
                "gtAlleleIndex : {gt} fAlleleIndex: {foreign} forwards: [D@0 reverse: [D@0 both: [D@0"
            ));
            let mut best = forward;
            if best[0] < reverse[0] {
                best = reverse;
            }
            if best[0] < both[0] {
                best = both;
            }
            if max_effective_depth > 0 {
                let best_score = best[0] - best[1];
                let closest = allele_counts_to_index(&[gt, 1, foreign, 1])
                    .expect("an even-length allele count array");
                let ploidy_model =
                    ploidy_model_likelihoods
                        .get(closest)
                        .copied()
                        .ok_or_else(|| {
                            format!(
                                "Index {closest} out of bounds for length {}",
                                ploidy_model_likelihoods.len()
                            )
                        })?
                        - -log10_one_half;
                let adjusted = ploidy_model
                    + ((best_score - ploidy_model)
                        * ((f64::from(depth.min(max_effective_depth)) * 1.0) / f64::from(depth)));
                output[index] = java_max(output[index], adjusted + best[1]);
                if debug::is_enabled() {
                    debug::println(&format!(
                        "best FRD likelihoods: {} P(F) score used: {}  use MaxEffectiveDepth: {max_effective_depth}",
                        java_double_to_string(best_score),
                        java_double_to_string(best[1])
                    ));
                    debug::println(&format!(
                        "Using array index {closest} for mixture gt with likelihood of {} adjusted based on depth: {depth}",
                        java_double_to_string(ploidy_model)
                    ));
                    debug::println(&format!("p_rG_adj : {}", java_double_to_string(adjusted)));
                }
            } else {
                output[index] = java_max(output[index], best[0]);
            }
        }
    }
    Ok(output)
}

/// `computeCriticalValues`: each read's `log10 P(F)` threshold, stored on the read, and the
/// distinct thresholds in ascending order (`HashSet<Double>` sorted by `Double.compareTo`).
fn critical_values(
    reads: &mut [DragenRead],
    log10_mapq_prior_adjustment: f64,
) -> Result<Vec<f64>, String> {
    let mut thresholds = Vec::with_capacity(reads.len());
    for container in reads.iter_mut() {
        let value = container.phred_mapping_quality.clone()? * -0.1 + log10_mapq_prior_adjustment;
        container.phred_pf = value;
        thresholds.push(value);
    }
    thresholds.sort_by(|a, b| a.total_cmp(b));
    // `Double.equals` compares `doubleToLongBits`, under which every NaN is one value.
    thresholds.dedup_by(|a, b| (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits());
    Ok(thresholds)
}

/// `computeFRDModelForStrandData`: the best `log10 P(F) + log10 P(reads | foreign mixture)` over the
/// thresholds, and the threshold that gave it.
fn frd_for_strand(
    likelihoods: &[Vec<f64>],
    homozygous: usize,
    foreign: usize,
    reads: &[DragenRead],
    selected: impl Fn(&DragenRead) -> bool,
    thresholds: &[f64],
) -> [f64; 2] {
    if reads.is_empty() {
        return [f64::NEG_INFINITY, 0.0];
    }
    let mut max_lpspi = f64::NEG_INFINITY;
    let mut lpf_applied = 0.0;
    let mut counter = 0;
    let foreign_given_cutoff = |container: &DragenRead, read: usize, log_prob_f: f64| {
        if container.phred_pf + 0.0000001 <= log_prob_f {
            f64::NEG_INFINITY
        } else {
            likelihoods[foreign][read]
        }
    };
    for &log_prob_f in thresholds {
        let mut ratio = 0.0;
        let mut denominator = 0.0;
        let mut local_max = f64::NEG_INFINITY;
        for container in reads {
            let Some(read) = container.index else {
                continue;
            };
            if selected(container) {
                let foreign_likelihood = foreign_given_cutoff(container, read, log_prob_f);
                let genotype_likelihood = likelihoods[homozygous][read];
                ratio += pow10(
                    foreign_likelihood
                        - approximate_log10_sum_log10(foreign_likelihood, genotype_likelihood),
                );
                denominator += 1.0;
            }
        }
        let foreign_fraction = java_min(ratio / denominator, 0.5);
        let log10_foreign = jmath::math::log10(foreign_fraction);
        let log10_not_foreign = jmath::math::log10(1.0 - foreign_fraction);
        let mut cumulative = 0.0;
        for container in reads {
            let Some(read) = container.index else {
                continue;
            };
            let genotype_likelihood = likelihoods[homozygous][read];
            if selected(container) {
                let foreign_likelihood = foreign_given_cutoff(container, read, log_prob_f);
                cumulative += approximate_log10_sum_log10(
                    log10_foreign + foreign_likelihood,
                    log10_not_foreign + genotype_likelihood,
                );
            } else {
                cumulative += genotype_likelihood;
            }
        }
        let lpsi = log_prob_f + cumulative;
        local_max = java_max(local_max, lpsi);
        if debug::is_enabled() {
            debug::println(&format!(
                "beta: {} localMaxLpspi: {} for lpf: {} with LP_R_GF: {} index: {counter}",
                java_double_to_string(foreign_fraction),
                java_double_to_string(local_max),
                java_double_to_string(log_prob_f),
                java_double_to_string(cumulative)
            ));
            counter += 1;
        }
        if local_max > max_lpspi {
            max_lpspi = java_max(max_lpspi, local_max);
            lpf_applied = log_prob_f;
        }
    }
    [max_lpspi, lpf_applied]
}

/// `Arrays.toString(double[])`.
fn java_doubles(values: &[f64]) -> String {
    let values: Vec<String> = values.iter().map(|v| java_double_to_string(*v)).collect();
    format!("[{}]", values.join(", "))
}

/// `String.format("%W.Pf", value)`: right-aligned in `width`.
fn fixed(value: f64, width: usize, places: usize) -> String {
    let text = if value.is_nan() {
        "NaN".to_string()
    } else if value.is_infinite() {
        if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string()
    } else {
        crate::java_format::format_decimals(value, places)
    };
    format!("{text:>width$}")
}
