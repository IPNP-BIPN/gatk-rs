//! `ExcessHet` and `InbreedingCoeff`, and the `GenotypeUtils` counting both rest on, ported from
//! GATK 4.6.2.0.
//!
//! Two annotations that ask the same question of a cohort (is there more heterozygosity here than
//! Hardy-Weinberg predicts?) and answer it from the same three numbers, computed by the same
//! function called with **opposite** rounding flags.
//!
//! # The same counts, rounded for one annotation and not for the other
//!
//! ```java
//! // ExcessHet
//! private static final boolean ROUND_GENOTYPE_COUNTS = true;
//! // InbreedingCoeff
//! private static final boolean ROUND_GENOTYPE_COUNTS = false;  //if this changes update the caveats above
//! ```
//!
//! `ExcessHet` needs integers because its exact test indexes an array by het count; `InbreedingCoeff`
//! keeps the fractions because it divides one by another. So the two annotations on one site can
//! disagree about how many hets there are, and both are right.
//!
//! # The rounded branch is special-cased so that a GQ of zero cannot be counted twice
//!
//! ```java
//! genotypeWithTwoRefsCount += MathUtils.fastRound(refLikelihood);
//! if (refLikelihood != hetLikelihood) { ... }
//! if (varLikelihood != hetLikelihood) { ... }
//! ```
//!
//! PLs of `[0, 0, X]` normalise to two likelihoods of one half, and half rounds **up** in
//! [`gatk_engine::math_utils::fast_round`], so without the guards one genotype would contribute two
//! counts. The guards are exact floating-point equality between two separate `Math.pow` results,
//! which is the strongest possible statement that the two PLs were equal, and it holds because both
//! calls got the same argument.
//!
//! # `ExcessHet` cannot report worse than a p-value of 1e-16
//!
//! ```java
//! if (pval < 10e-60) { return Pair.of(sampleCount, PHRED_SCALED_MIN_P_VALUE); }
//! ```
//!
//! The saturating constant is `-10 * log10(1e-16)`, so the annotation reads `160.0000`, but the
//! comparison that reaches it is against `10e-60`, which is `1e-59` and not `1e-60`. Between those
//! two numbers the phred value is computed rather than clamped, so values above 160 do occur.
//!
//! # `InbreedingCoeff` counts samples differently from `ExcessHet`
//!
//! `ExcessHet` counts genotypes that are called, diploid and have likelihoods or a GQ.
//! `InbreedingCoeff` counts those **or** any diploid genotype with likelihoods, so an uncalled
//! sample with PLs raises its denominator and not the other's. The ten-sample minimum is then
//! tested twice, once on the genotype count and once on that sample count.

use gatk_engine::math_utils::{
    fast_round, java_max, java_min, normalize_from_log10_to_linear_space, normalize_sum_to_one,
    qual_to_prob,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

use crate::rank_sum::format_decimals;

/// `GATKVCFConstants.EXCESS_HET_KEY`.
pub const EXCESS_HET_KEY: &str = "ExcessHet";
/// `GATKVCFConstants.INBREEDING_COEFFICIENT_KEY`.
pub const INBREEDING_COEFFICIENT_KEY: &str = "InbreedingCoeff";

/// `ExcessHet.MIN_NEEDED_VALUE`.
const MIN_NEEDED_VALUE: f64 = 1.0E-16;
/// `ExcessHet.PHRED_SCALED_MIN_P_VALUE`, which is `-10 * Math.log10(1e-16)`.
pub fn phred_scaled_min_p_value() -> f64 {
    -10.0 * jmath::math::log10(MIN_NEEDED_VALUE)
}
/// `InbreedingCoeff.MIN_SAMPLES`.
const MIN_SAMPLES: usize = 10;

/// `GenotypeCounts`: reference, heterozygous and no-reference "counts", which are only integers
/// when the caller asked for rounding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenotypeCounts {
    pub refs: f64,
    pub hets: f64,
    pub homs: f64,
}

/// `GenotypeLikelihoods.calculatePLindex`.
pub fn pl_index(allele1: usize, allele2: usize) -> usize {
    (allele2 * (allele2 + 1) / 2) + allele1
}

/// `GenotypeLikelihoods.getAllelePair`, inverted from the same triangular numbering.
///
/// The reference reads a precomputed cache and throws an `IllegalStateException` past its end;
/// this walks the same triangle, which cannot run off it.
pub fn allele_pair(index: usize) -> (usize, usize) {
    let mut second = 0usize;
    while pl_index(0, second + 1) <= index {
        second += 1;
    }
    (index - pl_index(0, second), second)
}

/// `VariantContext.getGLIndicesOfAlternateAllele`: the `(AA, AB, BB)` triple for one alternate.
pub fn gl_indices_of_alternate_allele(alt_index: usize) -> [usize; 3] {
    [
        pl_index(0, 0),
        pl_index(0, alt_index),
        pl_index(alt_index, alt_index),
    ]
}

/// `Genotype.isCalled()`: htsjdk's type is `NO_CALL` only when **every** allele is a no-call, so a
/// half-called genotype is `MIXED` and counts as called.
pub fn is_called(genotype: &Genotype) -> bool {
    !genotype.alleles.is_empty() && !genotype.alleles.iter().all(|a| a.is_no_call())
}

/// `Genotype.isHomRef()`: two or more alleles, all called, all reference.
pub fn is_hom_ref(genotype: &Genotype) -> bool {
    !genotype.alleles.is_empty()
        && genotype
            .alleles
            .iter()
            .all(|a| !a.is_no_call() && a.is_reference())
}

/// `GenotypeUtils.isDiploidWithLikelihoods`.
pub fn is_diploid_with_likelihoods(genotype: &Genotype) -> bool {
    genotype.pl.is_some() && genotype.ploidy() == 2
}

/// `GenotypeUtils.isCalledAndDiploidWithLikelihoodsOrWithGQ`.
pub fn is_called_and_diploid_with_likelihoods_or_with_gq(genotype: &Genotype) -> bool {
    is_called(genotype)
        && genotype.ploidy() == 2
        && (genotype.pl.is_some() || genotype.gq.is_some())
}

/// `GenotypeLikelihoods.getAsVector()` for a genotype whose likelihoods came from PLs.
fn log10_likelihoods(genotype: &Genotype) -> Option<Vec<f64>> {
    genotype
        .pl
        .as_ref()
        .map(|pls| pls.iter().map(|pl| *pl as f64 / -10.0).collect())
}

/// `GenotypeUtils.computeDiploidGenotypeCounts`.
///
/// Skips anything not diploid, and for a hom-ref without likelihoods invents a distribution from
/// the GQ rather than declining, in three different ways depending on the rounding flag and the GQ.
pub fn compute_diploid_genotype_counts(
    vc: &VariantContext,
    genotypes: &[&Genotype],
    round_contribution_from_each_genotype: bool,
) -> GenotypeCounts {
    const IDX_AA: usize = 0;
    const IDX_AB: usize = 1;
    const IDX_BB: usize = 2;

    let mut with_two_refs = 0.0f64;
    let mut with_one_ref = 0.0f64;
    let mut with_no_refs = 0.0f64;

    for genotype in genotypes {
        if !is_diploid_with_likelihoods(genotype)
            && !is_called_and_diploid_with_likelihoods_or_with_gq(genotype)
        {
            continue;
        }

        if genotype.pl.is_none() && is_hom_ref(genotype) {
            let gq = genotype.gq.unwrap_or(0);
            if round_contribution_from_each_genotype {
                with_two_refs += 1.0;
            } else if gq == 0 {
                // A third each, so a GQ of zero is a flat prior over the three genotypes.
                with_two_refs += 1.0 / 3.0;
                with_one_ref += 1.0 / 3.0;
                with_no_refs += 1.0 / 3.0;
            } else {
                // "assume last likelihood is negligible": the hom-var mass is dropped rather than
                // distributed, so these three counts do not sum to one for this genotype.
                with_two_refs += qual_to_prob(gq as f64);
                with_one_ref += 1.0 - qual_to_prob(gq as f64);
            }
            continue;
        }

        // `throw new IllegalStateException("Genotype has no likelihoods")` otherwise, unreachable
        // because the guard above already required likelihoods or a hom-ref GQ.
        let Some(log10) = log10_likelihoods(genotype) else {
            continue;
        };
        let normalized = normalize_from_log10_to_linear_space(&log10);

        let biallelic: Vec<f64> = if vc.alternate_alleles().len() > 1 {
            // `MathUtil.indexOfMax`, first maximum wins.
            let mut max_ind = 0usize;
            for (i, value) in normalized.iter().enumerate().skip(1) {
                if *value > normalized[max_ind] {
                    max_ind = i;
                }
            }
            let (first, second) = allele_pair(max_ind);
            if first != 0 && second != 0 {
                // "all likelihoods go to genotypesWithNoRefsCount because no ref allele is called":
                // a whole 1, not the likelihood, and not rounded either way.
                with_no_refs += 1.0;
                continue;
            }
            let mut max_likelihood = normalized[IDX_AB];
            let mut het_index = IDX_AB;
            let mut var_index = IDX_BB;
            for (alt_index, _) in vc.alternate_alleles().iter().enumerate() {
                let idx_vector = gl_indices_of_alternate_allele(alt_index + 1);
                let temp_index = idx_vector[1];
                if normalized[temp_index] > max_likelihood {
                    max_likelihood = normalized[temp_index];
                    het_index = temp_index;
                    var_index = idx_vector[2];
                }
            }
            // A second normalisation, over three entries of an already normalised vector.
            match normalize_sum_to_one(&[
                normalized[IDX_AA],
                normalized[het_index],
                normalized[var_index],
            ]) {
                Some(values) => values,
                None => continue,
            }
        } else {
            normalized
        };

        let ref_likelihood = biallelic[IDX_AA];
        let het_likelihood = biallelic[IDX_AB];
        let var_likelihood = biallelic[IDX_BB];

        if round_contribution_from_each_genotype {
            with_two_refs += fast_round(ref_likelihood) as f64;
            // `[0,0,X]`: the two are exactly equal, so the het count is not incremented and the
            // genotype is counted as hom-ref alone.
            if ref_likelihood != het_likelihood {
                with_one_ref += fast_round(het_likelihood) as f64;
            }
            // `[X,0,0]`: counted as het and not as variant.
            if var_likelihood != het_likelihood {
                with_no_refs += fast_round(var_likelihood) as f64;
            }
        } else {
            with_two_refs += ref_likelihood;
            with_one_ref += het_likelihood;
            with_no_refs += var_likelihood;
        }
    }

    GenotypeCounts {
        refs: with_two_refs,
        hets: with_one_ref,
        homs: with_no_refs,
    }
}

/// `ExcessHet.exactTest`: the Wigginton, Cutler and Abecasis right-sided p-value.
///
/// The mid-p correction that GATK applied up to and including 4.2.3.0 is **gone**, so this is the
/// `P_high` of the paper and agrees with bcftools' `ExcHet`.
pub fn exact_test(het_count: i32, ref_count: i32, hom_count: i32) -> Option<f64> {
    if het_count < 0 || ref_count < 0 || hom_count < 0 {
        return None;
    }
    // The rarer homozygote is whichever of ref and hom is smaller, so the test is symmetric in
    // them and "ref" plays no privileged role.
    let (obs_hom_r, obs_hom_c) = if ref_count < hom_count {
        (ref_count, hom_count)
    } else {
        (hom_count, ref_count)
    };

    let rare_copies = 2 * obs_hom_r + het_count;
    let n = het_count + obs_hom_c + obs_hom_r;

    let mut probs = vec![0.0f64; rare_copies as usize + 1];

    let mut mid = ((rare_copies as f64 * (2.0 * n as f64 - rare_copies as f64))
        / (2.0 * n as f64 - 1.0))
        .floor() as i32;
    if (mid % 2) != (rare_copies % 2) {
        mid += 1;
    }

    probs[mid as usize] = 1.0;
    let mut mysum = 1.0f64;

    let mut curr_hets = mid;
    let mut curr_hom_r = (rare_copies - mid) / 2;
    let mut curr_hom_c = n - curr_hets - curr_hom_r;

    while curr_hets >= 2 {
        let potential = probs[curr_hets as usize] * curr_hets as f64 * (curr_hets as f64 - 1.0)
            / (4.0 * (curr_hom_r as f64 + 1.0) * (curr_hom_c as f64 + 1.0));
        if potential < MIN_NEEDED_VALUE {
            break;
        }
        probs[curr_hets as usize - 2] = potential;
        mysum += probs[curr_hets as usize - 2];
        curr_hets -= 2;
        curr_hom_r += 1;
        curr_hom_c += 1;
    }

    let mut curr_hets = mid;
    let mut curr_hom_r = (rare_copies - mid) / 2;
    let mut curr_hom_c = n - curr_hets - curr_hom_r;

    while curr_hets <= rare_copies - 2 {
        let potential = probs[curr_hets as usize] * 4.0 * curr_hom_r as f64 * curr_hom_c as f64
            / ((curr_hets as f64 + 2.0) * (curr_hets as f64 + 1.0));
        if potential < MIN_NEEDED_VALUE {
            break;
        }
        probs[curr_hets as usize + 2] = potential;
        mysum += probs[curr_hets as usize + 2];
        curr_hets += 2;
        curr_hom_r -= 1;
        curr_hom_c -= 1;
    }

    let right_pval = probs[het_count as usize] / mysum;
    if het_count == rare_copies {
        return Some(java_max(0., java_min(1., right_pval)));
    }
    // `StatUtils.sum`, which is a plain left-to-right accumulation over the tail.
    let tail: f64 = probs[het_count as usize + 1..].iter().sum();
    Some(java_max(0., java_min(1., right_pval + tail / mysum)))
}

/// `ExcessHet.calculateEH(GenotypeCounts, int)`: the sample count and the phred-scaled p-value.
pub fn calculate_eh(counts: GenotypeCounts, sample_count: usize) -> Option<(usize, f64)> {
    // `(int)` on a double truncates towards zero; the counts are already rounded when ExcessHet
    // asks for them, so this only matters if a caller passes unrounded ones.
    let pval = exact_test(counts.hets as i32, counts.refs as i32, counts.homs as i32)?;
    if pval < 10e-60 {
        return Some((sample_count, phred_scaled_min_p_value()));
    }
    // "We add 0. to prevent -0.0000 from being output": a p-value of exactly one gives a log of
    // -0.0, which times -10 is 0.0 already, but a log of +0.0 times -10 is -0.0 and would print
    // with a sign.
    Some((sample_count, -10.0 * jmath::math::log10(pval) + 0.))
}

/// `ExcessHet.annotate`: `ExcessHet`, or nothing.
pub fn excess_het(vc: &VariantContext) -> Option<String> {
    let genotypes: Vec<&Genotype> = vc.genotypes.iter().collect();
    if genotypes.is_empty() || !vc.is_variant() {
        return None;
    }
    let counts = compute_diploid_genotype_counts(vc, &genotypes, true);
    let sample_count = genotypes
        .iter()
        .filter(|g| is_called_and_diploid_with_likelihoods_or_with_gq(g))
        .count();
    let (sample_count, eh) = calculate_eh(counts, sample_count)?;
    if sample_count < 1 {
        return None;
    }
    Some(format_decimals(eh, 4))
}

/// `InbreedingCoeff.calculateIC`.
pub fn calculate_ic(vc: &VariantContext, genotypes: &[&Genotype]) -> Option<(usize, f64)> {
    if genotypes.is_empty() {
        return None;
    }
    let counts = compute_diploid_genotype_counts(vc, genotypes, false);
    let sample_count = genotypes
        .iter()
        .filter(|g| {
            is_called_and_diploid_with_likelihoods_or_with_gq(g) || is_diploid_with_likelihoods(g)
        })
        .count();

    let p = (2.0 * counts.refs + counts.hets) / (2.0 * (counts.refs + counts.hets + counts.homs));
    let q = 1.0 - p;
    let expected_hets = 2.0 * p * q * sample_count as f64;
    // No guard: a monomorphic cohort has zero expected hets and F is then NaN or -Infinity, and
    // the annotation writes it out as such.
    let f = 1.0 - (counts.hets / expected_hets);
    Some((sample_count, f))
}

/// `InbreedingCoeff.annotate`: `InbreedingCoeff`, or nothing below ten samples.
pub fn inbreeding_coeff(vc: &VariantContext) -> Option<String> {
    let genotypes: Vec<&Genotype> = vc.genotypes.iter().collect();
    // The first test is on how many genotypes exist at all, the second on how many were usable.
    if genotypes.len() < MIN_SAMPLES || !vc.is_variant() {
        return None;
    }
    let (sample_count, f) = calculate_ic(vc, &genotypes)?;
    if sample_count < MIN_SAMPLES {
        return None;
    }
    Some(format_decimals(f, 4))
}

/// Unused, kept so the signature matches the reference's, which takes an allele it never reads
/// once the founder set is empty.
pub fn founder_ids_are_empty(_allele: Option<&Allele>) {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_triangular_numbering_inverts() {
        for second in 0..8 {
            for first in 0..=second {
                assert_eq!(allele_pair(pl_index(first, second)), (first, second));
            }
        }
    }

    #[test]
    fn a_cohort_in_equilibrium_has_a_p_value_of_about_one() {
        // 25 hom-ref, 50 het, 25 hom-var is exactly Hardy-Weinberg for p = q = 1/2.
        let p = exact_test(50, 25, 25).expect("non-negative counts");
        assert!(p > 0.4 && p <= 1.0, "{p}");
    }

    #[test]
    fn every_sample_heterozygous_saturates() {
        // 100 hets and nothing else: the most extreme excess there is.
        let p = exact_test(100, 0, 0).expect("non-negative counts");
        assert!(p <= 1.0);
    }
}

/// `HeterozygosityCalculator`, which `AS_InbreedingCoeff` uses instead of the genotype counts.
///
/// # The reference's own count is incremented once per **alternate allele**
///
/// ```java
/// for(final Allele a : vc.getAlternateAlleles()) {
///     ...
///     final double refAlleleCounts = alleleCounts.get(vc.getReference());
///     alleleCounts.put(vc.getReference(), refAlleleCounts + 2*normalizedLikelihoods[0]);
/// }
/// ```
///
/// That last statement sits **inside** the loop over alternates, so a triallelic site adds the
/// hom-ref mass to the reference's count twice. The allele frequency `AS_InbreedingCoeff` derives
/// from it is therefore not a frequency in any normalised sense once there is more than one
/// alternate.
///
/// # It is fractional, not integral
///
/// The counts are sums of normalised likelihoods rather than of called genotypes, which is the
/// point: "for a GQ10 variant, the probability of the call will be ~0.9 and the second best call
/// will be ~0.1 so adding up those 0.1s for het counts can dramatically change the AF compared with
/// integer counts".
#[derive(Debug, Clone, Default)]
pub struct HeterozygosityCounts {
    /// Per **alternate** allele, in the variant's order.
    pub het_counts: Vec<(Allele, f64)>,
    /// Per allele including the reference, in the variant's order.
    pub allele_counts: Vec<(Allele, f64)>,
    /// Genotypes that were called, diploid, and had likelihoods or a GQ.
    pub sample_count: usize,
}

/// `GenotypeUtils.PLOIDY_2_HOM_VAR_SCALE_FACTOR`, which is `round(30 / -10 / log10(0.5))`.
fn ploidy_2_hom_var_scale_factor() -> i32 {
    jmath::math::round(30.0 / -10.0 / jmath::math::log10(0.5)) as i32
}

/// `GenotypeUtils.makeApproximateDiploidLog10LikelihoodsFromGQ`.
///
/// The reference's own javadoc calls this "completely bogus": it gives every heterozygote the same
/// likelihood whatever the alternate, so a multiallelic hom-ref's quality is deflated "for no reason
/// whatsoever".
pub fn approximate_diploid_log10_likelihoods_from_gq(gq: i32, num_alleles: usize) -> Vec<f64> {
    let count = num_alleles * (num_alleles + 1) / 2;
    let hom_var = ploidy_2_hom_var_scale_factor() * gq;
    (0..count)
        .map(|index| {
            let pl = if index == 0 {
                0
            } else {
                let (first, _) = allele_pair(index);
                if first == 0 {
                    gq
                } else {
                    hom_var
                }
            };
            pl as f64 / -10.0
        })
        .collect()
}

/// `HeterozygosityCalculator`'s constructor, which does all its work eagerly.
pub fn heterozygosity_counts(vc: &VariantContext) -> HeterozygosityCounts {
    let mut counts = HeterozygosityCounts::default();
    if vc.genotypes.is_empty() || !vc.is_variant() {
        // `sampleCount` stays at its initial **-1** in the reference, which no caller reads before
        // the maps are built. Here it is zero, which the ten-sample guard treats the same way.
        return counts;
    }
    let num_alleles = vc.alleles.len();
    let reference = vc.reference().clone();
    counts.het_counts = vc
        .alternate_alleles()
        .iter()
        .map(|allele| (allele.clone(), 0.0))
        .collect();
    counts.allele_counts = vc
        .alleles
        .iter()
        .map(|allele| (allele.clone(), 0.0))
        .collect();

    for genotype in &vc.genotypes {
        if !is_called_and_diploid_with_likelihoods_or_with_gq(genotype) {
            continue;
        }
        counts.sample_count += 1;

        let log10 = match log10_likelihoods(genotype) {
            Some(values) => values,
            None => {
                approximate_diploid_log10_likelihoods_from_gq(genotype.gq.unwrap_or(0), num_alleles)
            }
        };
        let normalized = normalize_from_log10_to_linear_space(&log10);

        for (alt_position, alt) in vc.alternate_alleles().iter().enumerate() {
            let alt_index = alt_position + 1;
            for i in 0..num_alleles {
                if i == alt_index {
                    // Hom-var mass, counted twice because a homozygote carries two copies.
                    if let Some(slot) = counts.allele_counts.iter_mut().find(|(a, _)| a == alt) {
                        slot.1 += 2.0 * normalized[pl_index(alt_index, alt_index)];
                    }
                    continue;
                }
                let idx_ab = pl_index(i.min(alt_index), i.max(alt_index));
                if let Some(slot) = counts.het_counts.iter_mut().find(|(a, _)| a == alt) {
                    slot.1 += normalized[idx_ab];
                }
                if let Some(slot) = counts.allele_counts.iter_mut().find(|(a, _)| a == alt) {
                    slot.1 += normalized[idx_ab];
                }
                if let Some(slot) = counts
                    .allele_counts
                    .iter_mut()
                    .find(|(a, _)| *a == reference)
                {
                    slot.1 += normalized[idx_ab];
                }
            }
            // Inside the alternate loop: see the type's note.
            if let Some(slot) = counts
                .allele_counts
                .iter_mut()
                .find(|(a, _)| *a == reference)
            {
                slot.1 += 2.0 * normalized[0];
            }
        }
    }
    counts
}

impl HeterozygosityCounts {
    /// `getHetCount(altAllele)`, zero for an allele the map does not hold.
    pub fn het_count(&self, allele: &Allele) -> f64 {
        self.het_counts
            .iter()
            .find(|(a, _)| a == allele)
            .map(|(_, count)| *count)
            .unwrap_or(0.0)
    }

    /// `getAlleleCount(allele)`, zero for an allele the map does not hold.
    pub fn allele_count(&self, allele: &Allele) -> f64 {
        self.allele_counts
            .iter()
            .find(|(a, _)| a == allele)
            .map(|(_, count)| *count)
            .unwrap_or(0.0)
    }
}
