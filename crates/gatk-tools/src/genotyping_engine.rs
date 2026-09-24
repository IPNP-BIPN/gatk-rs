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

/// `GenotypeBuilder.log10PError(e)`: the GQ it implies, capped at `VCFConstants.MAX_GENOTYPE_QUAL`.
///
/// ```java
/// return GQ((int)Math.round(Math.min(pLog10Error * -10, VCFConstants.MAX_GENOTYPE_QUAL)));
/// ```
pub(crate) fn gq_of_log10(log10: f64) -> i32 {
    java_round((log10 * -10.0).min(MAX_GENOTYPE_QUAL)) as i32
}

/// `VCFConstants.MAX_GENOTYPE_QUAL`.
const MAX_GENOTYPE_QUAL: f64 = 99.0;

/// The two `GenotypeAssignmentMethod`s the genotyper reaches: `PREFER_PLS` for the calls it
/// emits, `BEST_MATCH_TO_ORIGINAL` for the reduction to `--max-alternate-alleles`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubsetMethod {
    PreferPls,
    BestMatchToOriginal,
}

/// `makeGenotypeCall(ploidy, gb, method, likelihoods, allelesToUse, originalGT, gpc)`.
fn make_genotype_call(
    ploidy: usize,
    genotype: &mut Genotype,
    method: SubsetMethod,
    likelihoods: Option<&[f64]>,
    targets: &[Allele],
    original: &Genotype,
) {
    // Before the method: a hom-ref or no-call with GQ 0 is a no-call, and with DP 0 as well it
    // loses everything but its name.
    if (is_hom_ref(original) || is_no_call(original)) && original.gq == Some(0) {
        genotype.alleles = vec![Allele::no_call(); ploidy];
        if original.dp == Some(0) {
            genotype.pl = None;
            genotype.dp = None;
            genotype.ad = None;
            genotype.gq = None;
            genotype.extended.clear();
            return;
        }
    }
    if method == SubsetMethod::BestMatchToOriginal {
        // "no-calls are just going to cause problems": a GQ-0 call whose best likelihood was
        // already the reference stays uncalled.
        let uninformative =
            original.gq == Some(0) && original.pl.as_ref().is_none_or(|pl| pl.first() == Some(&0));
        genotype.alleles = if uninformative {
            vec![Allele::no_call(); ploidy]
        } else {
            best_match_to_original(targets, &original.alleles)
        };
        return;
    }
    let informative = likelihoods.filter(|gls| gls.iter().sum::<f64>() < SUM_GL_THRESH_NOCALL);
    let Some(gls) = informative else {
        genotype.alleles = best_match_to_original(targets, &original.alleles);
        return;
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
        built.pl = None;
        built.gq = None;
        if new_log10_gq != f64::NEG_INFINITY && g.gq.is_some() {
            built.gq = Some(gq_of_log10(new_log10_gq));
        }
        built.pl = new_likelihoods.as_deref().map(pls_of);

        make_genotype_call(
            g.ploidy(),
            &mut built,
            method,
            new_likelihoods.as_deref(),
            alleles_to_keep,
            g,
        );

        if g.extended.iter().any(|(key, _)| key == "SAC") {
            return Err(EngineError::Limitation(
                "subsetting SAC to the alleles kept is not ported yet.".to_string(),
            ));
        }
        if let (Some(ad), true) = (&g.ad, built.ad.is_some()) {
            built.ad = Some(kept.iter().map(|index| ad[*index]).collect());
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
        let best = max_element_index(&gls);
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
        }
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
            // `forceKeepAllele` is `annotateAllSitesWithPLs`, which GenotypeGVCFs never sets.
            let to_output = (plausible || lone_non_ref) && !spurious_span_del;
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
            self.output_allele_subset(&result, &reduced_alleles, vc);

        // `+ 0.0` turns a -0.0 into 0.0, twice, as the reference writes it.
        let log10_confidence = if !monomorphic {
            result.log10_posterior_of_no_variant + 0.0
        } else {
            result.log10_prob_variant_present() + 0.0
        };
        let phred = (-10.0 * log10_confidence) + 0.0;

        let passes_call = phred >= self.configuration.standard_confidence_for_calling;
        // `passesEmitThreshold` under EMIT_VARIANTS_ONLY: a variant site that passes the call.
        let passes_emit = !monomorphic && passes_call;
        let first_is_non_ref = alternatives.first().is_some_and(|a| *a == non_ref());
        if !passes_emit && !self.configuration.emit_all_active_sites && !first_is_non_ref {
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
            subset_alleles_prefer_pls(&vc.genotypes, default_ploidy, &vc.alleles, &output_alleles)?
        };
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
                        .map(|c| Value::Double((*c as f64 / an).min(1.0)))
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
