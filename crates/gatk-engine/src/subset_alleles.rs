//! `AlleleSubsettingUtils.subsetAlleles`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.genotyper` (GATK 4.6.2.0).
//!
//! What a split does to a genotype that carries likelihoods. The index machinery is
//! [`crate::genotype_index`]; this is what the reference does with it.
//!
//! # Keeping every allele still rewrites the GQ
//!
//! The PLs are permuted through [`crate::genotype_index::subsetted_pl_indices`] and then rescaled
//! by subtracting the maximum, so the smallest phred of the subset is 0 whatever it was. The GQ is
//! then **recomputed from the new PLs**, which changes it even when the permutation was the
//! identity: a genotype with `50,0,60,40,30,70` and `GQ 50` keeps its PLs and comes out with
//! `GQ 30`, because 30 is what the second-best genotype was all along.
//!
//! It is kept unchanged in exactly one case, the subset to the reference alone:
//!
//! ```java
//! } else {  //if we subset to just ref allele, keep the GQ
//!     newLog10GQ = g.getGQ()/-10.0;
//! ```
//!
//! # The no-data cases come before the assignment method
//!
//! A hom-ref or no-call genotype with `GQ 0` is turned into a no-call whatever method was asked
//! for, and if its DP is 0 as well everything is cleared: no PL, no DP, no AD, no GQ, no
//! attributes. That is a different test from the one `BEST_MATCH_TO_ORIGINAL` applies afterwards.

use crate::genotype_index::{genotype_count, subsetted_pl_indices, GenotypeIndexError};

/// `GenotypeAssignmentMethod`, as far as a split reaches it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssignmentMethod {
    /// `DO_NOT_ASSIGN_GENOTYPES`, which falls through every branch of `makeGenotypeCall`: the
    /// calls stay as they were while the likelihoods and the GQ are rewritten for the new allele
    /// list. SelectVariants' `--remove-unused-alternates` is what asks for it.
    DoNotAssignGenotypes,
    BestMatchToOriginal,
    UsePlsToAssign,
    SetToNoCall,
    SetToNoCallNoAnnotations,
}

/// `GATKVariantContextUtils.SUM_GL_THRESH_NOCALL`, which is a GQ near zero and not a confident one.
pub const SUM_GL_THRESH_NOCALL: f64 = -0.1;

/// One sample's call, as much of it as subsetting touches.
#[derive(Debug, Clone, PartialEq)]
pub struct Genotype {
    /// Allele indices into the record's list; `None` is a no-call.
    pub alleles: Vec<Option<usize>>,
    /// Phred-scaled likelihoods, in the canonical genotype order.
    pub pl: Option<Vec<i32>>,
    pub gq: Option<i32>,
    pub ad: Option<Vec<i32>>,
    pub dp: Option<i32>,
    /// Everything else, in the order it was set.
    pub attributes: Vec<(String, String)>,
}

impl Genotype {
    /// `isHomRef()`, which needs every call to be the reference.
    fn is_hom_ref(&self) -> bool {
        !self.alleles.is_empty() && self.alleles.iter().all(|allele| *allele == Some(0))
    }

    fn is_no_call(&self) -> bool {
        !self.alleles.is_empty() && self.alleles.iter().all(Option::is_none)
    }
}

/// `subsetAlleles(genotypes, defaultPloidy, originalAlleles, allelesToKeep, null, method)`.
///
/// `kept` is the index in the original allele list of each allele the new list keeps, in the new
/// list's order.
pub fn subset_alleles(
    genotypes: &[Genotype],
    default_ploidy: usize,
    original_allele_count: usize,
    kept: &[usize],
    method: AssignmentMethod,
) -> Result<Vec<Genotype>, GenotypeIndexError> {
    let mut out = Vec::with_capacity(genotypes.len());
    for genotype in genotypes {
        out.push(subset_one(
            genotype,
            default_ploidy,
            original_allele_count,
            kept,
            method,
        )?);
    }
    Ok(out)
}

fn subset_one(
    genotype: &Genotype,
    default_ploidy: usize,
    original_allele_count: usize,
    kept: &[usize],
    method: AssignmentMethod,
) -> Result<Genotype, GenotypeIndexError> {
    let ploidy = if genotype.alleles.is_empty() {
        default_ploidy
    } else {
        genotype.alleles.len()
    };
    let expected = genotype_count(ploidy, original_allele_count)?;
    let indices = subsetted_pl_indices(ploidy, kept)?;

    // The PLs, permuted and rescaled. An array of the wrong length is dropped entirely.
    let new_likelihoods: Option<Vec<f64>> = match &genotype.pl {
        Some(pl) if pl.len() == expected => {
            let original: Vec<f64> = pl.iter().map(|phred| f64::from(*phred) / -10.0).collect();
            let subset: Vec<f64> = indices.iter().map(|index| original[*index]).collect();
            Some(scale_log_space(&subset))
        }
        _ => None,
    };

    // The GQ, recomputed from the new likelihoods, kept when the subset is the reference alone,
    // and carried through a different branch when there were no likelihoods to begin with.
    let mut new_log10_gq = f64::NEG_INFINITY;
    match (&genotype.pl, &new_likelihoods) {
        (Some(_), Some(likelihoods)) if likelihoods.len() > 1 => {
            let best = max_element_index(likelihoods);
            new_log10_gq = gq_log10_from_likelihoods(best, likelihoods);
        }
        (Some(_), Some(_)) => {
            // Subset to the reference alone: the old GQ stands.
            if let Some(gq) = genotype.gq {
                new_log10_gq = f64::from(gq) / -10.0;
            }
        }
        (Some(_), None) => {}
        (None, _) => {
            if let Some(gq) = genotype.gq {
                new_log10_gq = -0.1 * f64::from(gq);
            }
        }
    }

    let mut result = Genotype {
        alleles: genotype.alleles.clone(),
        pl: new_likelihoods
            .as_ref()
            .map(|likelihoods| pls_of(likelihoods)),
        // `gb.noPL().noGQ()`: the old GQ is invalid, and only the new one is put back, and only
        // when the genotype had one to begin with.
        gq: if new_log10_gq != f64::NEG_INFINITY && genotype.gq.is_some() {
            Some((new_log10_gq * -10.0).round() as i32)
        } else {
            None
        },
        // The reference builds from a COPY of the genotype, so the builder starts with the old
        // AD and `gb.noPL().noGQ().noAttributes()` clears only those three.
        ad: genotype.ad.clone(),
        dp: genotype.dp,
        // Posteriors and priors are invalid for a new allele list, whatever else happens.
        attributes: genotype
            .attributes
            .iter()
            .filter(|(key, _)| !matches!(key.as_str(), "GP" | "PG" | "PP"))
            .cloned()
            .collect(),
    };

    make_genotype_call(
        &mut result,
        ploidy,
        method,
        new_likelihoods.as_deref(),
        genotype,
        kept,
    );

    // ```java
    // if (g.hasAD() && gb.makeWithShallowCopy().hasAD()) {
    // ```
    //
    // The builder must STILL have an AD for the new one to be set, so the no-data branch that
    // cleared it wins. Otherwise the AD is PERMUTED, not summed, and a dropped allele takes its
    // depth with it.
    if let (Some(ad), true) = (&genotype.ad, result.ad.is_some()) {
        result.ad = Some(kept.iter().map(|index| ad[*index]).collect());
    }
    Ok(result)
}

/// `makeGenotypeCall`, as far as a split reaches it.
fn make_genotype_call(
    result: &mut Genotype,
    ploidy: usize,
    method: AssignmentMethod,
    likelihoods: Option<&[f64]>,
    original: &Genotype,
    kept: &[usize],
) {
    // Before the method is consulted at all: a hom-ref or no-call with GQ 0.
    if method != AssignmentMethod::SetToNoCall
        && (original.is_hom_ref() || original.is_no_call())
        && original.gq == Some(0)
    {
        result.alleles = vec![None; ploidy];
        if original.dp == Some(0) {
            result.pl = None;
            result.dp = None;
            result.ad = None;
            result.gq = None;
            result.attributes.clear();
            return;
        }
    }

    match method {
        // `DO_NOT_ASSIGN_GENOTYPES` matches none of the reference's branches, so the calls stay as
        // the builder copied them. The likelihoods and the GQ above were rewritten regardless,
        // which is what makes it different from doing nothing at all.
        AssignmentMethod::DoNotAssignGenotypes => {}
        AssignmentMethod::SetToNoCall => result.alleles = vec![None; ploidy],
        AssignmentMethod::SetToNoCallNoAnnotations => {
            result.alleles = vec![None; ploidy];
            result.gq = None;
            result.ad = None;
            result.pl = None;
            result.attributes.clear();
        }
        AssignmentMethod::UsePlsToAssign => {
            if let Some(likelihoods) = likelihoods {
                let best = max_element_index(likelihoods);
                let gq = gq_log10_from_likelihoods(best, likelihoods);
                let called = genotype_of_index(ploidy, likelihoods.len(), best);
                if best == 0 && gq > SUM_GL_THRESH_NOCALL {
                    // The reference is most likely AND the call is not informative.
                    result.alleles = vec![None; ploidy];
                } else {
                    result.alleles = called.into_iter().map(Some).collect();
                }
            } else {
                result.alleles = best_match_to_original(&original.alleles, kept, ploidy);
            }
        }
        AssignmentMethod::BestMatchToOriginal => {
            let uninformative =
                original.gq == Some(0) && original.pl.as_ref().is_none_or(|pl| pl[0] == 0);
            if uninformative {
                result.alleles = vec![None; ploidy];
            } else {
                result.alleles = best_match_to_original(&original.alleles, kept, ploidy);
            }
        }
    }
}

/// `bestMatchToOriginalGT`: a call the new list still has, or the reference.
///
/// The alleles are indices into the ORIGINAL list on the way in and into the NEW one on the way
/// out, which is the whole of the translation: an allele the new list keeps moves to its new
/// position, and one it does not becomes the reference.
fn best_match_to_original(
    alleles: &[Option<usize>],
    kept: &[usize],
    ploidy: usize,
) -> Vec<Option<usize>> {
    if alleles.is_empty() {
        return vec![None; ploidy];
    }
    alleles
        .iter()
        .map(|allele| allele.map(|index| kept.iter().position(|old| *old == index).unwrap_or(0)))
        .collect()
}

/// The genotype at one index of a likelihood vector, as allele indices.
fn genotype_of_index(ploidy: usize, count: usize, index: usize) -> Vec<usize> {
    let alleles = (1..=count).find(|alleles| {
        genotype_count(ploidy, *alleles)
            .map(|size| size >= count)
            .unwrap_or(false)
    });
    let alleles = alleles.unwrap_or(1);
    crate::genotype_index::genotypes_in_canonical_order(ploidy, alleles)
        .into_iter()
        .nth(index)
        .unwrap_or_default()
}

/// `scaleLogSpaceArrayForNumericalStability`: subtract the maximum.
fn scale_log_space(values: &[f64]) -> Vec<f64> {
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    values.iter().map(|value| value - max).collect()
}

/// `MathUtils.maxElementIndex`: the FIRST index of the maximum.
fn max_element_index(values: &[f64]) -> usize {
    let mut best = 0;
    for (index, value) in values.iter().enumerate() {
        if *value > values[best] {
            best = index;
        }
    }
    best
}

/// `GenotypeLikelihoods.getGQLog10FromLikelihoods`.
///
/// The best likelihood minus the best of the others, negated; and when the chosen genotype is not
/// the most likely one, the log10 of one minus its normalised probability instead.
fn gq_log10_from_likelihoods(chosen: usize, likelihoods: &[f64]) -> f64 {
    let mut other = f64::NEG_INFINITY;
    for (index, value) in likelihoods.iter().enumerate() {
        if index != chosen {
            // `>=`, so the LAST of equal maxima wins, which changes nothing about the value.
            if *value >= other {
                other = *value;
            }
        }
    }
    let quality = likelihoods[chosen] - other;
    if quality < 0.0 {
        let total: f64 = likelihoods.iter().map(|value| 10f64.powf(*value)).sum();
        let normalised = 10f64.powf(likelihoods[chosen]) / total;
        (1.0 - normalised).log10()
    } else {
        -quality
    }
}

/// `GLsToPLs`: rescale to the maximum and round, which is what a PL array is.
fn pls_of(likelihoods: &[f64]) -> Vec<i32> {
    let max = likelihoods
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    likelihoods
        .iter()
        .map(|value| (-10.0 * (value - max)).round() as i32)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn het() -> Genotype {
        Genotype {
            alleles: vec![Some(0), Some(1)],
            pl: Some(vec![50, 0, 60, 40, 30, 70]),
            gq: Some(50),
            ad: Some(vec![10, 12, 8]),
            dp: Some(30),
            attributes: Vec::new(),
        }
    }

    #[test]
    fn keeping_every_allele_still_rewrites_the_gq() {
        let out = subset_alleles(
            &[het()],
            2,
            3,
            &[0, 1, 2],
            AssignmentMethod::BestMatchToOriginal,
        )
        .expect("subset");
        assert_eq!(out[0].pl, Some(vec![50, 0, 60, 40, 30, 70]));
        assert_eq!(out[0].gq, Some(30));
    }

    #[test]
    fn the_pls_are_permuted_and_rescaled() {
        let out = subset_alleles(
            &[het()],
            2,
            3,
            &[0, 2],
            AssignmentMethod::BestMatchToOriginal,
        )
        .expect("subset");
        assert_eq!(out[0].pl, Some(vec![10, 0, 30]));
        assert_eq!(out[0].gq, Some(10));
        assert_eq!(out[0].ad, Some(vec![10, 8]));
    }
}
