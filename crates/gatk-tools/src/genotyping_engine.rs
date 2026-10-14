//! `GenotypingEngine.calculateGenotypes`, as `MinimalGenotypingEngine` runs it for GenotypeGVCFs.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine`,
//! `org.broadinstitute.hellbender.tools.walkers.genotyper.MinimalGenotypingEngine` and
//! `org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils`
//! (GATK 4.6.2.0).
//!
//! Given a merged record whose samples carry likelihoods, the engine decides the site's QUAL, which
//! alternates survive, and every sample's call against the alleles kept.
//!
//! # QUAL is the no-variant posterior, and only a variant site is written
//!
//! The Dirichlet fit of [`gatk_engine::allele_frequency_calculator`] gives the posterior that no
//! alternate is present; QUAL is minus ten times its log. An alternate is kept when its own
//! absent-posterior passes the calling threshold, and a site whose alternates all fail is
//! monomorphic. A monomorphic site is dropped unless output is forced, which is what
//! `--include-non-variant-sites` and `--force-output-intervals` ask for; a variant site below the
//! threshold is written with the `LowQual` filter.
//!
//! # A spanning deletion is kept only when a deletion upstream still covers the site
//!
//! The engine remembers every deletion it has emitted. A `*` whose site no earlier deletion covers
//! is spurious and dropped, and a site whose only surviving alternate is `*` is dropped too unless
//! output is forced.
//!
//! # The calls are PREFER_PLS
//!
//! GenotypeGVCFs asks for `PREFER_PLS`: a sample whose likelihoods are informative is called from
//! them, and one whose likelihoods are flat or absent keeps its original call, matched to the
//! alleles kept. A hom-ref or no-call with a GQ of 0 is made a no-call first, and one with a DP of
//! 0 as well loses everything but its name.

use std::collections::VecDeque;

use gatk_engine::allele_frequency_calculator::{AfCalculationResult, AlleleFrequencyCalculator};
use gatk_engine::genotype_index::{
    genotype_count, genotypes_in_canonical_order, subsetted_pl_indices,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

use crate::reference_confidence_merger::{
    best_match_to_original, gq_log10_from_likelihoods, java_round, max_element_index, non_ref,
    SUM_GL_THRESH_NOCALL,
};

/// `GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED`.
const MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED: usize = 50;

/// `GATKVCFConstants.MLE_ALLELE_COUNT_KEY` and its siblings.
pub const MLE_ALLELE_COUNT_KEY: &str = "MLEAC";
pub const MLE_ALLELE_FREQUENCY_KEY: &str = "MLEAF";
pub const AS_QUAL_KEY: &str = "AS_QUAL";
pub const NUMBER_OF_DISCOVERED_ALLELES_KEY: &str = "NDA";
pub const LOW_QUAL_FILTER_NAME: &str = "LowQual";

/// What the engine is configured with: `StandardCallerArgumentCollection` as
/// `GenotypeGVCFsEngine.createMinimalArgs` fills it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Configuration {
    /// `--standard-min-confidence-threshold-for-calling`.
    pub standard_confidence_for_calling: f64,
    /// `--max-alternate-alleles`.
    pub max_alternate_alleles: usize,
    /// `--sample-ploidy`.
    pub sample_ploidy: usize,
    /// `--annotate-with-num-discovered-alleles`.
    pub annotate_number_of_alleles_discovered: bool,
    /// `OutputMode.EMIT_ALL_ACTIVE_SITES`, which forced output sets.
    pub emit_all_active_sites: bool,
    /// `doAlleleSpecificCalcs`: whether any allele-specific info annotation was resolved.
    pub allele_specific: bool,
    /// `OutputMode.EMIT_ALL_CONFIDENT_SITES`, which ReblockGVCF sets: the emit threshold is the
    /// call threshold, whatever the best guess.
    pub emit_all_confident_sites: bool,
    /// `annotateAllSitesWithPLs`, which ReblockGVCF sets: a monomorphic site's QUAL is the
    /// probability that only the reference exists, as a polymorphic one's is.
    pub annotate_all_sites_with_pls: bool,
    /// `HaplotypeCallerGenotypingEngine.forceKeepAllele` in GVCF mode: every alternate is output.
    pub force_keep_all_alleles: bool,
    /// `genotypeArgs.genotypeAssignmentMethod`, which assigns the emitted genotypes.
    pub assignment_method: SubsetMethod,
}

/// What the engine refuses, in the reference's classes, or the port's own limitation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EngineError {
    /// A reference exception, by class and message.
    Runtime { class: String, message: String },
    /// A path this port does not carry yet.
    Limitation(String),
}

impl From<gatk_engine::allele_frequency_calculator::AfCalcError> for EngineError {
    fn from(error: gatk_engine::allele_frequency_calculator::AfCalcError) -> Self {
        EngineError::Runtime {
            class: error.java_class().to_string(),
            message: error.message().to_string(),
        }
    }
}

/// The engine, with the upstream deletions it has emitted.
pub struct GenotypingEngine {
    pub configuration: Configuration,
    calculator: AlleleFrequencyCalculator,
    /// `upstreamDeletionsLoc`: contig, start and end of each deletion emitted so far.
    upstream_deletions: VecDeque<(String, i64, i64)>,
    /// The genotype prior calculator HaplotypeCaller passes to `calculateGenotypes`; `None` for
    /// the tools that pass null.
    genotype_priors: Option<GenotypePriors>,
    /// `--use-posteriors-to-calculate-qual`.
    use_posteriors_to_calculate_qual: bool,
}

fn is_span_del(allele: &Allele) -> bool {
    !allele.is_reference() && allele.display_string() == htsjdk_vcf::allele::SPAN_DEL_STRING
}

fn is_hom_ref(genotype: &Genotype) -> bool {
    crate::reference_confidence_merger::is_hom_ref(genotype)
}

fn is_no_call(genotype: &Genotype) -> bool {
    crate::reference_confidence_merger::is_no_call(genotype)
}

/// `GenotypeLikelihoods.GLsToPLs`: each likelihood less the maximum, times minus ten, rounded.
pub fn pls_of(likelihoods: &[f64]) -> Vec<i32> {
    let adjust = likelihoods
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    likelihoods
        .iter()
        .map(|gl| java_round((-10.0 * (gl - adjust)).min(i32::MAX as f64)) as i32)
        .collect()
}

/// `MathUtils.scaleLogSpaceArrayForNumericalStability`: subtract the maximum.
fn scale_log_space(values: &[f64]) -> Vec<f64> {
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    values.iter().map(|value| value - max).collect()
}

/// `GenotypeBuilder.log10PError(e)` in htsjdk 4.2.0: the GQ it implies, NOT capped.
///
/// ```java
/// return GQ((int)Math.round(pLog10Error * -10));
/// ```
///
/// The cap at `MAX_GENOTYPE_QUAL` is the VCF encoder's (`GQAccessor`), so a GQ over 99 is written
/// as 99 while an annotation such as `GQ_MEAN` still reads the genotype's real value.
pub(crate) fn gq_of_log10(log10: f64) -> i32 {
    java_round(log10 * -10.0) as i32
}

/// `GenotypeAssignmentMethod`. GenotypeGVCFs reaches `PREFER_PLS` for the calls it emits and
/// `BEST_MATCH_TO_ORIGINAL` for the reduction to `--max-alternate-alleles`; ReblockGVCF's
/// `--genotype-assignment-method` reaches the rest.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubsetMethod {
    SetToNoCall,
    SetToNoCallNoAnnotations,
    UsePlsToAssign,
    PreferPls,
    BestMatchToOriginal,
    /// Throws, since no genotype prior calculator is ever passed where the tools reach it.
    UsePosteriorProbabilities,
    /// Neither branch of `makeGenotypeCall` touches the alleles.
    DoNotAssignGenotypes,
    /// As `DO_NOT_ASSIGN_GENOTYPES`: `makeGenotypeCall` has no branch for it.
    UsePosteriorsAnnotation,
}

impl SubsetMethod {
    /// The enum constant's name, as `--genotype-assignment-method` takes it.
    pub fn from_name(name: &str) -> Option<SubsetMethod> {
        Some(match name {
            "SET_TO_NO_CALL" => SubsetMethod::SetToNoCall,
            "SET_TO_NO_CALL_NO_ANNOTATIONS" => SubsetMethod::SetToNoCallNoAnnotations,
            "USE_PLS_TO_ASSIGN" => SubsetMethod::UsePlsToAssign,
            "PREFER_PLS" => SubsetMethod::PreferPls,
            "BEST_MATCH_TO_ORIGINAL" => SubsetMethod::BestMatchToOriginal,
            "USE_POSTERIOR_PROBABILITIES" => SubsetMethod::UsePosteriorProbabilities,
            "DO_NOT_ASSIGN_GENOTYPES" => SubsetMethod::DoNotAssignGenotypes,
            "USE_POSTERIORS_ANNOTATION" => SubsetMethod::UsePosteriorsAnnotation,
            _ => return None,
        })
    }
}

/// `GenotypePriorCalculator`: log10 genotype priors from a het and a hom value per allele type
/// (reference, SNP, indel, anything else), the SNP ones spread over the three other bases.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenotypePriors {
    het: [f64; 4],
    hom: [f64; 4],
    diff: [f64; 4],
}

const PRIOR_REF: usize = 0;
const PRIOR_SNP: usize = 1;
const PRIOR_INDEL: usize = 2;
const PRIOR_OTHER: usize = 3;

impl GenotypePriors {
    fn new(snp: (f64, f64), indel: (f64, f64), other: (f64, f64)) -> Self {
        // `LOG10_SNP_NORMALIZATION_CONSTANT`: log10 of the three bases a SNP can turn into.
        let normalization = jmath::math::log10(3.0);
        let mut het = [0.0; 4];
        let mut hom = [0.0; 4];
        het[PRIOR_SNP] = snp.0 - normalization;
        hom[PRIOR_SNP] = snp.1 - normalization;
        het[PRIOR_INDEL] = indel.0;
        hom[PRIOR_INDEL] = indel.1;
        het[PRIOR_OTHER] = other.0;
        hom[PRIOR_OTHER] = other.1;
        let mut diff = [0.0; 4];
        for i in 0..4 {
            diff[i] = hom[i] - het[i];
        }
        GenotypePriors { het, hom, diff }
    }

    /// `assumingHW(log10SnpHet, log10IndelHet)`: a hom is twice the het, and anything that is
    /// neither a SNP nor an indel takes the larger of the two.
    pub fn assuming_hw(log10_snp_het: f64, log10_indel_het: f64) -> Self {
        let other = gatk_engine::math_utils::java_max(log10_snp_het, log10_indel_het);
        GenotypePriors::new(
            (log10_snp_het, log10_snp_het * 2.0),
            (log10_indel_het, log10_indel_het * 2.0),
            (other, other * 2.0),
        )
    }

    /// `givenHetToHomRatio(log10SnpHet, log10IndelHet, log10OtherHet, hetHomRatio)`: each hom is
    /// its het less `log10(hetHomRatio)`.
    pub fn given_het_to_hom_ratio(
        log10_snp_het: f64,
        log10_indel_het: f64,
        log10_other_het: f64,
        het_hom_ratio: f64,
    ) -> Self {
        let log10_ratio = jmath::math::log10(het_hom_ratio);
        GenotypePriors::new(
            (log10_snp_het, log10_snp_het - log10_ratio),
            (log10_indel_het, log10_indel_het - log10_ratio),
            (log10_other_het, log10_other_het - log10_ratio),
        )
    }

    /// `givenDragstrParams(params, period, repeats, snpHeterozygosity, hetHomRatio)`: the indel het
    /// is DRAGstr's API at the site's STR, times -0.1, and anything else takes the larger of it
    /// and the SNP het. `snp_heterozygosity` is already a log10, as the reference passes it.
    pub fn given_dragstr_params(
        params: &gatk_engine::dragstr::DragstrParams,
        period: i32,
        repeats: i32,
        snp_heterozygosity: f64,
        het_hom_ratio: f64,
    ) -> Result<Self, gatk_engine::dragstr::DragstrError> {
        let indel_het = -0.1 * params.api(period, repeats)?;
        let other_het = gatk_engine::math_utils::java_max(snp_heterozygosity, indel_het);
        Ok(GenotypePriors::given_het_to_hom_ratio(
            snp_heterozygosity,
            indel_het,
            other_het,
            het_hom_ratio,
        ))
    }

    /// `getLog10Priors(ploidy, alleles)`: zero for the hom-ref genotype, and for every other the
    /// sum over its alleles of the het value (raised by the difference per extra copy) or, at
    /// exactly two copies, the hom value.
    pub fn log10_priors(&self, ploidy: usize, alleles: &[Allele]) -> Result<Vec<f64>, EngineError> {
        let reference_length = alleles.first().map_or(0, Allele::len);
        let mut types = Vec::with_capacity(alleles.len());
        for allele in alleles {
            types.push(if allele.is_reference() {
                PRIOR_REF
            } else if !allele.is_no_call() && !allele.is_symbolic() {
                if allele.len() == reference_length {
                    PRIOR_SNP
                } else {
                    PRIOR_INDEL
                }
            } else if matches!(allele.display_string().as_str(), "<INS>" | "<DEL>") {
                return Err(EngineError::Runtime {
                    class: "java.lang.IllegalArgumentException".to_string(),
                    message: format!("cannot handle symbolic indels: {}", allele.display_string()),
                });
            } else {
                PRIOR_OTHER
            });
        }
        let genotypes = genotypes_in_canonical_order(ploidy, alleles.len());
        let mut result = vec![0.0; genotypes.len()];
        for (index, genotype) in genotypes.iter().enumerate().skip(1) {
            let mut sum = 0.0;
            let mut a = 0;
            while a < genotype.len() {
                let allele = genotype[a];
                let count = genotype[a..].iter().take_while(|&&x| x == allele).count();
                let t = types[allele];
                sum += if count == 2 {
                    self.hom[t]
                } else {
                    self.het[t] + self.diff[t] * (count as f64 - 1.0)
                };
                a += count;
            }
            result[index] = sum;
        }
        Ok(result)
    }
}

/// `GATKVariantContextUtils.formatGP`: two decimals, with trailing zeros (and a bare point) cut.
fn format_gp(value: f64) -> String {
    let formatted = gatk_engine::java_format::format_decimals(value, 2);
    if let Some(stripped) = formatted.strip_suffix("00") {
        stripped.trim_end_matches('.').to_string()
    } else if let Some(stripped) = formatted.strip_suffix('0') {
        stripped.to_string()
    } else {
        formatted
    }
}

/// `MathUtils.log10SumLog10(a, b)`, the two-argument form.
fn log10_sum_log10_two(a: f64, b: f64) -> f64 {
    if a > b {
        a + jmath::math::log10(1.0 + gatk_engine::math_utils::pow10(b - a))
    } else {
        b + jmath::math::log10(1.0 + gatk_engine::math_utils::pow10(a - b))
    }
}

/// `MathUtils.log10SumLog10(a, b, c)`, the three-argument form.
fn log10_sum_log10_three(a: f64, b: f64, c: f64) -> f64 {
    use gatk_engine::math_utils::pow10;
    if a >= b && a >= c {
        a + jmath::math::log10(1.0 + pow10(b - a) + pow10(c - a))
    } else if b >= c {
        b + jmath::math::log10(1.0 + pow10(a - b) + pow10(c - b))
    } else {
        c + jmath::math::log10(1.0 + pow10(a - c) + pow10(b - c))
    }
}

/// `MathUtils.log10sumLog10(array, start, finish)`, the lower-case form: one element is returned
/// as it is, and the sum is `1.0 + (the other terms)`.
fn log10sum_log10(values: &[f64], start: usize, finish: usize) -> f64 {
    if finish - start < 2 {
        return if finish == start {
            f64::NEG_INFINITY
        } else {
            values[start]
        };
    }
    let max_index = gatk_engine::math_utils::max_element_index(values, start, finish);
    let max_value = values[max_index];
    if max_value == f64::NEG_INFINITY {
        return max_value;
    }
    let mut sum = 0.0;
    for (i, value) in values.iter().enumerate().take(finish).skip(start) {
        if i != max_index {
            sum += gatk_engine::math_utils::pow10(value - max_value);
        }
    }
    max_value + jmath::math::log10(1.0 + sum)
}

/// `GATKVariantContextUtils.getGQLog10FromPosteriors`.
fn gq_log10_from_posteriors(best: usize, posteriors: &[f64]) -> Option<f64> {
    match posteriors.len() {
        0 | 1 => None,
        2 => Some(if best == 0 {
            posteriors[1]
        } else {
            posteriors[0]
        }),
        3 => Some(gatk_engine::math_utils::java_min(
            0.0,
            log10_sum_log10_two(
                posteriors[if best == 0 { 2 } else { best - 1 }],
                posteriors[if best == 2 { 0 } else { best + 1 }],
            ),
        )),
        n => Some(gatk_engine::math_utils::java_min(
            0.0,
            if best == 0 {
                gatk_engine::math_utils::log10_sum_log10_range(posteriors, 1, n)
            } else if best == n - 1 {
                gatk_engine::math_utils::log10_sum_log10_range(posteriors, 0, best)
            } else {
                log10_sum_log10_two(
                    log10sum_log10(posteriors, 0, best),
                    log10sum_log10(posteriors, best + 1, n),
                )
            },
        )),
    }
}

/// `makeGenotypeCall(ploidy, gb, method, likelihoods, allelesToUse, originalGT, null)`.
pub(crate) fn make_genotype_call(
    ploidy: usize,
    genotype: &mut Genotype,
    method: SubsetMethod,
    likelihoods: Option<&[f64]>,
    targets: &[Allele],
    original: &Genotype,
) -> Result<(), EngineError> {
    make_genotype_call_with_priors(
        ploidy,
        genotype,
        method,
        likelihoods,
        targets,
        original,
        None,
    )
}

/// `makeGenotypeCall` with the genotype prior calculator `USE_POSTERIOR_PROBABILITIES` reads.
pub(crate) fn make_genotype_call_with_priors(
    ploidy: usize,
    genotype: &mut Genotype,
    method: SubsetMethod,
    likelihoods: Option<&[f64]>,
    targets: &[Allele],
    original: &Genotype,
    priors: Option<&GenotypePriors>,
) -> Result<(), EngineError> {
    // Before the method: a hom-ref or no-call with GQ 0 is a no-call, and with DP 0 as well it
    // loses everything but its name. `SET_TO_NO_CALL` alone skips it, "so CombineGVCFs and
    // GenomicsDB keep attributes".
    if method != SubsetMethod::SetToNoCall
        && (is_hom_ref(original) || is_no_call(original))
        && original.gq == Some(0)
    {
        genotype.alleles = vec![Allele::no_call(); ploidy];
        if original.dp == Some(0) {
            genotype.pl = None;
            genotype.dp = None;
            genotype.ad = None;
            genotype.gq = None;
            genotype.extended.clear();
            return Ok(());
        }
    }
    match method {
        SubsetMethod::SetToNoCall => genotype.alleles = vec![Allele::no_call(); ploidy],
        SubsetMethod::UsePlsToAssign | SubsetMethod::PreferPls => {
            let informative =
                likelihoods.filter(|gls| gls.iter().sum::<f64>() < SUM_GL_THRESH_NOCALL);
            let usable = if method == SubsetMethod::PreferPls {
                informative
            } else {
                likelihoods
            };
            let Some(gls) = usable else {
                genotype.alleles = best_match_to_original(targets, &original.alleles);
                return Ok(());
            };
            let best = max_element_index(gls);
            let called: Vec<Allele> = genotypes_in_canonical_order(ploidy, targets.len())
                .into_iter()
                .nth(best)
                .unwrap_or_default()
                .into_iter()
                .map(|index| targets[index].clone())
                .collect();
            let gq = gq_log10_from_likelihoods(best, gls);
            if called.contains(&non_ref()) {
                genotype.alleles = vec![targets[0].clone(); ploidy];
                genotype.pl = Some(vec![0; gls.len()]);
                genotype.gq = Some(0);
            } else if best == 0 && gq > SUM_GL_THRESH_NOCALL {
                genotype.alleles = vec![Allele::no_call(); ploidy];
            } else {
                genotype.alleles = called;
            }
            if targets.len() > 1 {
                genotype.gq = Some(gq_of_log10(gq));
            }
        }
        SubsetMethod::SetToNoCallNoAnnotations => {
            genotype.alleles = vec![Allele::no_call(); ploidy];
            genotype.gq = None;
            genotype.ad = None;
            genotype.pl = None;
            genotype.extended.clear();
        }
        SubsetMethod::BestMatchToOriginal => {
            // "no-calls are just going to cause problems": a GQ-0 call whose best likelihood was
            // already the reference stays uncalled.
            let uninformative = original.gq == Some(0)
                && original.pl.as_ref().is_none_or(|pl| pl.first() == Some(&0));
            genotype.alleles = if uninformative {
                vec![Allele::no_call(); ploidy]
            } else {
                best_match_to_original(targets, &original.alleles)
            };
        }
        SubsetMethod::UsePosteriorProbabilities => {
            let Some(priors) = priors else {
                return Err(EngineError::Runtime {
                    class: "org.broadinstitute.hellbender.exceptions.GATKException".to_string(),
                    message: "cannot uses posteriors without an genotype prior calculator present"
                        .to_string(),
                });
            };
            let Some(gls) = likelihoods else {
                return Err(EngineError::Runtime {
                    class: "java.lang.NullPointerException".to_string(),
                    message: "genotype likelihoods are null".to_string(),
                });
            };
            let log10_priors = priors.log10_priors(ploidy, targets)?;
            let posteriors: Vec<f64> = log10_priors.iter().zip(gls).map(|(p, l)| p + l).collect();
            let max = posteriors.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let normalized: Vec<f64> = posteriors.iter().map(|v| v - max).collect();
            let phred = |values: &[f64]| -> Value {
                Value::List(
                    values
                        .iter()
                        .map(|&v| Value::Str(format_gp(if v == 0.0 { 0.0 } else { v * -10.0 })))
                        .collect(),
                )
            };
            genotype
                .extended
                .push(("GP".to_string(), phred(&normalized)));
            genotype
                .extended
                .push(("PG".to_string(), phred(&log10_priors)));
            let best = max_element_index(&posteriors);
            if !targets.is_empty() {
                // `log10PError`; the encoder caps what it writes at 99.
                genotype.gq = gq_log10_from_posteriors(best, &normalized).map(gq_of_log10);
            }
            genotype.alleles = genotypes_in_canonical_order(ploidy, targets.len())
                .into_iter()
                .nth(best)
                .unwrap_or_default()
                .into_iter()
                .map(|index| targets[index].clone())
                .collect();
        }
        SubsetMethod::DoNotAssignGenotypes | SubsetMethod::UsePosteriorsAnnotation => {}
    }
    Ok(())
}

/// `AlleleSubsettingUtils.subsetAlleles(genotypes, defaultPloidy, originalAlleles, allelesToKeep,
/// null, PREFER_PLS)`.
pub fn subset_alleles_prefer_pls(
    genotypes: &[Genotype],
    default_ploidy: usize,
    original_alleles: &[Allele],
    alleles_to_keep: &[Allele],
) -> Result<Vec<Genotype>, EngineError> {
    subset_alleles(
        genotypes,
        default_ploidy,
        original_alleles,
        alleles_to_keep,
        SubsetMethod::PreferPls,
    )
}

/// `AlleleSubsettingUtils.subsetAlleles(genotypes, defaultPloidy, originalAlleles, allelesToKeep,
/// null, method)`.
pub fn subset_alleles(
    genotypes: &[Genotype],
    default_ploidy: usize,
    original_alleles: &[Allele],
    alleles_to_keep: &[Allele],
    method: SubsetMethod,
) -> Result<Vec<Genotype>, EngineError> {
    subset_alleles_with_length_annotations(
        genotypes,
        default_ploidy,
        original_alleles,
        alleles_to_keep,
        method,
        &[],
    )
}

/// `AlleleSubsettingUtils.subsetAlleles(genotypes, defaultPloidy, originalAlleles, allelesToKeep,
/// null, method, alleleBasedLengthAnnots)`: the FORMAT attributes whose header count is A, R or G
/// are removed, except F1R2 and F2R1, which are remapped to the alleles kept.
pub fn subset_alleles_with_length_annotations(
    genotypes: &[Genotype],
    default_ploidy: usize,
    original_alleles: &[Allele],
    alleles_to_keep: &[Allele],
    method: SubsetMethod,
    allele_based_length_annotations: &[String],
) -> Result<Vec<Genotype>, EngineError> {
    subset_alleles_with_priors(
        genotypes,
        default_ploidy,
        original_alleles,
        alleles_to_keep,
        method,
        allele_based_length_annotations,
        None,
    )
}

/// `AlleleSubsettingUtils.subsetAlleles` with a genotype prior calculator, which
/// `USE_POSTERIOR_PROBABILITIES` assigns from.
pub fn subset_alleles_with_priors(
    genotypes: &[Genotype],
    default_ploidy: usize,
    original_alleles: &[Allele],
    alleles_to_keep: &[Allele],
    method: SubsetMethod,
    allele_based_length_annotations: &[String],
    priors: Option<&GenotypePriors>,
) -> Result<Vec<Genotype>, EngineError> {
    let kept: Vec<usize> = alleles_to_keep
        .iter()
        .map(|allele| {
            original_alleles
                .iter()
                .position(|original| original == allele)
                .ok_or_else(|| EngineError::Runtime {
                    class: "java.lang.IllegalArgumentException".to_string(),
                    message: format!(
                        "allele {} is not in the original list",
                        allele.display_string()
                    ),
                })
        })
        .collect::<Result<_, _>>()?;
    let index_error =
        |error: gatk_engine::genotype_index::GenotypeIndexError| EngineError::Runtime {
            class: error.java_class().to_string(),
            message: error.message(),
        };
    let mut out = Vec::with_capacity(genotypes.len());
    for g in genotypes {
        let ploidy = if g.ploidy() > 0 {
            g.ploidy()
        } else {
            default_ploidy
        };
        let indices = subsetted_pl_indices(ploidy, &kept).map_err(index_error)?;
        let expected = genotype_count(ploidy, original_alleles.len()).map_err(index_error)?;

        let mut new_likelihoods: Option<Vec<f64>> = None;
        let mut new_log10_gq = f64::NEG_INFINITY;
        if let Some(pl) = &g.pl {
            if pl.len() == expected {
                let original: Vec<f64> = pl.iter().map(|p| *p as f64 / -10.0).collect();
                let subset: Vec<f64> = indices.iter().map(|index| original[*index]).collect();
                let scaled = scale_log_space(&subset);
                new_log10_gq = if scaled.len() > 1 {
                    gq_log10_from_likelihoods(max_element_index(&scaled), &scaled)
                } else {
                    // `g.getGQ()/-10.0`, which is -1 over -10 when there is no GQ.
                    g.gq.unwrap_or(-1) as f64 / -10.0
                };
                new_likelihoods = Some(scaled);
            }
        } else if let Some(gq) = g.gq {
            new_log10_gq = -0.1 * gq as f64;
        }

        let mut built = g.clone();
        built
            .extended
            .retain(|(key, _)| !matches!(key.as_str(), "PP" | "GP" | "PG"));
        let mut removed_allele_fraction = false;
        let mut remapped: Vec<(String, Value)> = Vec::new();
        for (key, value) in &built.extended {
            if !allele_based_length_annotations.contains(key) {
                remapped.push((key.clone(), value.clone()));
            } else if key == "F1R2" || key == "F2R1" {
                let counts: Vec<i32> = crate::reference_confidence_merger::value_to_string(value)
                    .split(',')
                    .map(|count| count.trim().parse().unwrap_or(0))
                    .collect();
                remapped.push((
                    key.clone(),
                    Value::List(
                        kept.iter()
                            .map(|index| {
                                Value::Int(i64::from(counts.get(*index).copied().unwrap_or(0)))
                            })
                            .collect(),
                    ),
                ));
            } else if key == "AF" {
                removed_allele_fraction = true;
            }
        }
        built.extended = remapped;
        built.pl = None;
        built.gq = None;
        if new_log10_gq != f64::NEG_INFINITY && g.gq.is_some() {
            built.gq = Some(gq_of_log10(new_log10_gq));
        }
        built.pl = new_likelihoods.as_deref().map(pls_of);

        make_genotype_call_with_priors(
            g.ploidy(),
            &mut built,
            method,
            new_likelihoods.as_deref(),
            alleles_to_keep,
            g,
            priors,
        )?;

        if g.extended.iter().any(|(key, _)| key == "SAC") {
            return Err(EngineError::Limitation(
                "subsetting SAC to the alleles kept is not ported yet.".to_string(),
            ));
        }
        if let (Some(ad), true) = (&g.ad, built.ad.is_some()) {
            let new_ad: Vec<i32> = kept.iter().map(|index| ad[*index]).collect();
            // "if we have recalculated AD and the original genotype had AF but was then removed,
            // then recalculate AF based on AD counts", without the reference's entry.
            if removed_allele_fraction {
                let total: f64 = new_ad.iter().map(|count| f64::from(*count)).sum();
                built.extended.push((
                    "AF".to_string(),
                    Value::List(
                        new_ad[1..]
                            .iter()
                            .map(|count| Value::Double(f64::from(*count) / total))
                            .collect(),
                    ),
                ));
            }
            built.ad = Some(new_ad);
        }
        out.push(built);
    }
    Ok(out)
}

/// `AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, defaultPloidy, numAltAllelesToKeep,
/// false)`: the reference, `<NON_REF>` if present, and the alternates whose summed likelihood
/// margins are largest.
///
/// Each sample adds, to every alternate of its most likely genotype, how far that genotype's
/// likelihood is from the hom-ref's. The samples are visited in name order, which is the order
/// the sums are accumulated in. Ties keep the lower index: the sort is stable and descending.
pub fn most_likely_alleles(
    vc: &VariantContext,
    default_ploidy: usize,
    keep: usize,
) -> Result<Vec<Allele>, EngineError> {
    most_likely_alleles_ensuring_alt(vc, default_ploidy, keep, false)
}

/// `calculateMostLikelyAlleles(vc, defaultPloidy, numAltAllelesToKeep, ensureReturnContainsAlt)`:
/// with the flag set, a sample's best genotype is looked for past the hom-ref, so a confident
/// reference call still votes for its best alternate.
pub fn most_likely_alleles_ensuring_alt(
    vc: &VariantContext,
    default_ploidy: usize,
    keep: usize,
    ensure_return_contains_alt: bool,
) -> Result<Vec<Allele>, EngineError> {
    let has_non_ref = vc.alleles.contains(&non_ref());
    let proper_alternates = vc.alleles.len() - if has_non_ref { 2 } else { 1 };
    if keep >= proper_alternates {
        return Ok(vc.alleles.clone());
    }
    let mut sums = vec![0.0f64; vc.alleles.len()];
    let mut ordered: Vec<&Genotype> = vc.genotypes.iter().collect();
    ordered.sort_by(|a, b| a.sample_name.cmp(&b.sample_name));
    for genotype in ordered {
        let Some(pl) = &genotype.pl else { continue };
        let gls: Vec<f64> = pl.iter().map(|p| *p as f64 / -10.0).collect();
        // `maxElementIndex(glsVector, start, length)`: the first maximum past `start`.
        let start = usize::from(ensure_return_contains_alt).min(gls.len().saturating_sub(1));
        let best = start + max_element_index(&gls[start..]);
        let margin = (gls[best] - gls[0]).abs();
        let ploidy = if genotype.ploidy() > 0 {
            genotype.ploidy()
        } else {
            default_ploidy
        };
        let called = genotypes_in_canonical_order(ploidy, vc.alleles.len())
            .into_iter()
            .nth(best)
            .unwrap_or_default();
        for (allele, sum) in sums.iter_mut().enumerate().skip(1) {
            if called.contains(&allele) {
                *sum += margin;
            }
        }
    }
    let non_ref_index = vc.alleles.iter().position(|a| *a == non_ref());
    let mut candidates: Vec<usize> = (1..vc.alleles.len())
        .filter(|index| Some(*index) != non_ref_index)
        .collect();
    candidates.sort_by(|a, b| sums[*b].total_cmp(&sums[*a]));
    candidates.truncate(keep);
    Ok(vc
        .alleles
        .iter()
        .enumerate()
        .filter(|(index, _)| {
            *index == 0 || Some(*index) == non_ref_index || candidates.contains(index)
        })
        .map(|(_, allele)| allele.clone())
        .collect())
}

/// `GATKVariantContextUtils.subsetToRefOnly`: every sample the reference, keeping only DP and GQ.
pub fn subset_to_ref_only(vc: &VariantContext, default_ploidy: usize) -> Vec<Genotype> {
    let reference = vc.reference().clone();
    vc.genotypes
        .iter()
        .map(|g| {
            let ploidy = if g.ploidy() == 0 {
                default_ploidy
            } else {
                g.ploidy()
            };
            let mut made = Genotype::new(&g.sample_name, vec![reference.clone(); ploidy]);
            made.dp = g.dp;
            made.gq = g.gq;
            made
        })
        .collect()
}

impl GenotypingEngine {
    pub fn new(configuration: Configuration, calculator: AlleleFrequencyCalculator) -> Self {
        GenotypingEngine {
            configuration,
            calculator,
            upstream_deletions: VecDeque::new(),
            genotype_priors: None,
            use_posteriors_to_calculate_qual: false,
        }
    }

    /// The genotype prior calculator `calculateGenotypes` is handed, as HaplotypeCaller's
    /// `resolveGenotypePriorCalculator` builds it.
    pub fn set_genotype_priors(&mut self, priors: GenotypePriors) {
        self.genotype_priors = Some(priors);
    }

    /// `--use-posteriors-to-calculate-qual`.
    pub fn set_use_posteriors_to_calculate_qual(&mut self, value: bool) {
        self.use_posteriors_to_calculate_qual = value;
    }

    /// `cannotBeGenotyped`.
    fn cannot_be_genotyped(vc: &VariantContext) -> bool {
        !(vc.alleles.len() <= MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED
            && vc
                .genotypes
                .iter()
                .filter(|g| !(is_no_call(g) || is_hom_ref(g)))
                .all(|g| g.pl.is_some())
            && vc.genotypes.iter().any(|g| g.pl.is_some()))
    }

    /// `recordDeletions`.
    fn record_deletions(&mut self, vc: &VariantContext, emitted: &[Allele]) {
        while let Some((contig, _, end)) = self.upstream_deletions.front() {
            if *contig != vc.contig || *end < vc.start {
                self.upstream_deletions.pop_front();
            } else {
                break;
            }
        }
        for allele in emitted {
            let deletion_size = vc.reference().len() as i64 - allele.len() as i64;
            if deletion_size > 0 {
                self.upstream_deletions.push_back((
                    vc.contig.clone(),
                    vc.start,
                    vc.start + deletion_size,
                ));
            }
        }
    }

    /// `isVcCoveredByDeletion`: an upstream deletion starting before the site and reaching it.
    fn is_covered_by_deletion(&self, vc: &VariantContext) -> bool {
        self.upstream_deletions.iter().any(|(contig, start, end)| {
            *contig == vc.contig && *start < vc.start && vc.start <= *end
        })
    }

    /// `calculateOutputAlleleSubset`: the alternates kept, their MLE counts, and whether the site
    /// is monomorphic.
    fn output_allele_subset(
        &self,
        result: &AfCalculationResult,
        alleles: &[Allele],
        vc: &VariantContext,
        forced: &[Allele],
    ) -> (Vec<Allele>, Vec<i32>, bool) {
        let mut output = Vec::new();
        let mut mle_counts = Vec::new();
        let mut monomorphic = true;
        let alternative_count = alleles.len() - 1;
        for (index, allele) in alleles.iter().enumerate().skip(1) {
            let alt = index - 1;
            let lone_non_ref = alternative_count == 1 && *allele == non_ref();
            let plausible =
                result.passes_threshold(alt, self.configuration.standard_confidence_for_calling);
            let spurious_span_del = is_span_del(allele) && !self.is_covered_by_deletion(vc);
            let to_output = (plausible
                || self.configuration.force_keep_all_alleles
                || lone_non_ref
                || forced.contains(allele))
                && !spurious_span_del;
            monomorphic &= !(plausible && !spurious_span_del);
            if to_output {
                output.push(allele.clone());
                mle_counts.push(result.allele_counts_of_mle[alt]);
            }
        }
        (output, mle_counts, monomorphic)
    }

    /// `calculateGenotypes(vc, null, emptyList())`.
    pub fn calculate_genotypes(
        &mut self,
        vc: &VariantContext,
    ) -> Result<Option<VariantContext>, EngineError> {
        self.calculate_genotypes_forcing(vc, &[])
    }

    /// `calculateGenotypes(vc, gpc, givenAlleles)`, `forced` being the alternates
    /// `allelesConsistentWithGivenAlleles` found: each is output whatever its confidence, and a
    /// site with one is emitted whatever its quality.
    pub fn calculate_genotypes_forcing(
        &mut self,
        vc: &VariantContext,
        forced: &[Allele],
    ) -> Result<Option<VariantContext>, EngineError> {
        if Self::cannot_be_genotyped(vc) || vc.genotypes.is_empty() {
            return Ok(None);
        }
        let default_ploidy = self.configuration.sample_ploidy;
        // The reduction to `--max-alternate-alleles`: the AF calculation sees only the most likely
        // alternates, with the genotypes matched to them, while everything after it still reads
        // the ORIGINAL record.
        let (reduced_alleles, reduced_genotypes) = if self.configuration.max_alternate_alleles
            < vc.alleles.len() - 1
        {
            let keep =
                most_likely_alleles(vc, default_ploidy, self.configuration.max_alternate_alleles)?;
            let genotypes = if keep.len() == 1 {
                subset_to_ref_only(vc, default_ploidy)
            } else {
                subset_alleles(
                    &vc.genotypes,
                    default_ploidy,
                    &vc.alleles,
                    &keep,
                    SubsetMethod::BestMatchToOriginal,
                )?
            };
            (keep, genotypes)
        } else {
            (vc.alleles.clone(), vc.genotypes.to_vec())
        };
        let result = self.calculator.calculate_with_ploidy(
            &vc.contig,
            vc.start,
            &reduced_alleles,
            &reduced_genotypes,
            default_ploidy,
        )?;
        let (alternatives, mle_counts, monomorphic) =
            self.output_allele_subset(&result, &reduced_alleles, vc, forced);

        // `+ 0.0` turns a -0.0 into 0.0, twice, as the reference writes it.
        let log10_confidence = if !monomorphic || self.configuration.annotate_all_sites_with_pls {
            result.log10_posterior_of_no_variant + 0.0
        } else {
            result.log10_prob_variant_present() + 0.0
        };
        let phred = (-10.0 * log10_confidence) + 0.0;

        let passes_call = phred >= self.configuration.standard_confidence_for_calling;
        // `passesEmitThreshold`: under EMIT_VARIANTS_ONLY a variant site that passes the call,
        // under EMIT_ALL_CONFIDENT_SITES any site that does.
        let passes_emit =
            (self.configuration.emit_all_confident_sites || !monomorphic) && passes_call;
        let first_is_non_ref = alternatives.first().is_some_and(|a| *a == non_ref());
        if !passes_emit
            && !self.configuration.emit_all_active_sites
            && !first_is_non_ref
            && forced.is_empty()
        {
            return Ok(None);
        }
        if !self.configuration.emit_all_active_sites
            && alternatives.len() == 1
            && is_span_del(&alternatives[0])
        {
            return Ok(None);
        }

        let mut output_alleles = vec![vc.reference().clone()];
        output_alleles.extend(alternatives.iter().cloned());
        self.record_deletions(vc, &output_alleles);

        let mut out = VariantContext::new(&vc.contig, vc.start, output_alleles.clone());
        out.stop = vc.stop;
        out.log10_p_error = log10_confidence;
        if !passes_call {
            out.filters = Some(vec![LOW_QUAL_FILTER_NAME.to_string()]);
        }
        let genotypes = if output_alleles.len() == 1 {
            subset_to_ref_only(vc, default_ploidy)
        } else {
            subset_alleles_with_priors(
                &vc.genotypes,
                default_ploidy,
                &vc.alleles,
                &output_alleles,
                self.configuration.assignment_method,
                &[],
                self.genotype_priors.as_ref(),
            )?
        };
        // `usePosteriorProbabilitiesToCalculateQual`: QUAL again, from the posteriors the
        // genotypes now carry, read back from their two-decimal GP strings.
        if self.use_posteriors_to_calculate_qual
            && genotypes
                .iter()
                .any(|g| g.extended.iter().any(|(key, _)| key == "GP"))
        {
            let no_variant = phred_no_variant_posterior(&output_alleles, &genotypes);
            let log10_no_variant = no_variant * -0.1;
            let update = if !monomorphic || self.configuration.annotate_all_sites_with_pls {
                log10_no_variant + 0.0
            } else {
                gatk_engine::allele_frequency_calculator::log10_one_minus_pow10(log10_no_variant)
                    + 0.0
            };
            if !update.is_nan() {
                out.log10_p_error = update;
            }
        }
        out.attributes = self.compose_call_attributes(
            vc,
            &reduced_alleles,
            &mle_counts,
            &result,
            &output_alleles,
            &genotypes,
        );
        out.genotypes = GenotypesContext::new(genotypes);
        Ok(Some(out))
    }

    /// `composeCallAttributes`: MLEAC and MLEAF, AS_QUAL for allele-specific runs, and NDA.
    fn compose_call_attributes(
        &self,
        vc: &VariantContext,
        genotyped_alleles: &[Allele],
        mle_counts: &[i32],
        result: &AfCalculationResult,
        output_alleles: &[Allele],
        genotypes: &[Genotype],
    ) -> Vec<(String, Value)> {
        let mut attributes = Vec::new();
        if !mle_counts.is_empty() {
            attributes.push((
                MLE_ALLELE_COUNT_KEY.to_string(),
                Value::List(mle_counts.iter().map(|c| Value::Int(*c as i64)).collect()),
            ));
            // `calculateMLEAlleleFrequencies`: over the called alleles of the subset genotypes.
            let an = genotypes
                .iter()
                .flat_map(|g| g.alleles.iter())
                .filter(|a| !a.is_no_call())
                .count() as f64;
            attributes.push((
                MLE_ALLELE_FREQUENCY_KEY.to_string(),
                Value::List(
                    mle_counts
                        .iter()
                        .map(|c| {
                            // `Math.min(1.0, count / AN)`, which keeps the NaN of no chromosomes.
                            let frequency = *c as f64 / an;
                            Value::Double(if frequency.is_nan() {
                                frequency
                            } else {
                                frequency.min(1.0)
                            })
                        })
                        .collect(),
                ),
            ));
        }
        if self.configuration.allele_specific {
            let mut quals: Vec<i32> = Vec::new();
            if result.log10_p_ref_by_allele.len() > 1 {
                for allele in output_alleles.iter().skip(1) {
                    let index = genotyped_alleles
                        .iter()
                        .position(|a| a == allele)
                        .expect("an allele the AF calculation saw")
                        - 1;
                    quals.push(java_round(result.log10_p_ref_by_allele[index] * -10.0) as i32);
                }
            } else {
                quals.push(java_round(result.log10_posterior_of_no_variant * -10.0) as i32);
            }
            attributes.push((
                AS_QUAL_KEY.to_string(),
                Value::List(quals.into_iter().map(|q| Value::Int(q as i64)).collect()),
            ));
        }
        if self.configuration.annotate_number_of_alleles_discovered {
            attributes.push((
                NUMBER_OF_DISCOVERED_ALLELES_KEY.to_string(),
                Value::Int((vc.alleles.len() - 1) as i64),
            ));
        }
        attributes
    }
}

/// `phredNoVariantPosteriorProbability`: over the genotypes, the sum of each one's
/// `extractPNoAlt`, the NaN ones skipped and NaN when every one is.
fn phred_no_variant_posterior(alleles: &[Allele], genotypes: &[Genotype]) -> f64 {
    let mut total = f64::NAN;
    for genotype in genotypes {
        let value = extract_p_no_alt(alleles, genotype);
        if value.is_nan() {
            continue;
        }
        total = if total.is_nan() { value } else { total + value };
    }
    total
}

/// `extractPNoAlt`: the hom-ref posterior less the phred sum of all of them (or, with a spanning
/// deletion, the phred sum of the genotypes made of the reference and `*` alone), from the GP
/// strings, NaN when a genotype carries none.
fn extract_p_no_alt(alleles: &[Allele], genotype: &Genotype) -> f64 {
    let posteriors: Vec<f64> = match genotype.extended.iter().find(|(key, _)| key == "GP") {
        Some((_, Value::List(values))) => values
            .iter()
            .map(|v| match v {
                Value::Str(text) => text.parse().unwrap_or(f64::NAN),
                Value::Double(d) => *d,
                Value::Int(i) => *i as f64,
                _ => f64::NAN,
            })
            .collect(),
        _ => vec![f64::NAN],
    };
    let phred_sum = |phreds: &[f64]| -> f64 {
        match phreds.len() {
            0 => f64::MAX,
            1 => phreds[0],
            2 => -10.0 * log10_sum_log10_two(phreds[0] * -0.1, phreds[1] * -0.1),
            3 => {
                -10.0 * log10_sum_log10_three(phreds[0] * -0.1, phreds[1] * -0.1, phreds[2] * -0.1)
            }
            _ => {
                let log10: Vec<f64> = phreds.iter().map(|p| p * -0.1).collect();
                -10.0 * gatk_engine::math_utils::log10_sum_log10(&log10)
            }
        }
    };
    match alleles.iter().position(is_span_del) {
        None => posteriors[0] - gatk_engine::math_utils::java_max(0.0, phred_sum(&posteriors)),
        Some(span_del) => {
            let ploidy = genotype.ploidy();
            let genotypes = genotypes_in_canonical_order(ploidy, alleles.len());
            let non_variant: Vec<f64> = (0..=ploidy)
                .map(|n| {
                    let mut wanted = vec![0usize; ploidy - n];
                    wanted.extend(std::iter::repeat_n(span_del, n));
                    let index = genotypes.iter().position(|g| *g == wanted).unwrap_or(0);
                    posteriors.get(index).copied().unwrap_or(f64::NAN)
                })
                .collect();
            gatk_engine::math_utils::java_max(0.0, phred_sum(&non_variant))
                - gatk_engine::math_utils::java_max(0.0, phred_sum(&posteriors))
        }
    }
}
