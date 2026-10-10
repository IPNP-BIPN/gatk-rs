//! Ported from `org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator`
//! and `IndependentSampleGenotypesModel` (GATK 4.6.2.0): one sample's genotype likelihoods from
//! its read-by-allele likelihoods.
//!
//! Each genotype, in canonical order, is the product over reads of the average, over the genotype's
//! allele copies, of the read's likelihood for the allele, computed three ways by how many distinct
//! alleles the genotype holds:
//!
//!  * one: the sum of that allele's log10 likelihoods;
//!  * two: per read `approximateLog10SumLog10` of each allele's log10 likelihood plus the log10 of
//!    its copy count, summed, less `reads * log10(ploidy)`. The approximation is the table-based
//!    one, so these genotypes carry its rounding;
//!  * three or more (possible only with three alleles and a ploidy above two): the likelihoods are
//!    rescaled per read by the read's maximum and taken out of log space once, each read's
//!    copy-weighted sum logged, and the rescaling added back.
//!
//! The likelihoods are then stored as PLs (`GenotypeLikelihoods.fromLog10Likelihoods`), which is
//! all a caller of the genotyping model sees.

use crate::genotype_index::{allele_counts_of, genotypes_in_canonical_order};
use crate::math_utils::pow10;
use crate::pair_hmm::approximate_log10_sum_log10;

/// `GenotypeLikelihoods.MAX_PL`.
const MAX_PL: f64 = i32::MAX as f64;

/// `computeLog10GenotypeLikelihoods(ploidy, log10AlleleLikelihoods)`, the matrix being indexed
/// `[allele][read]`.
pub fn compute_log10_genotype_likelihoods(
    ploidy: usize,
    by_allele_and_read: &[Vec<f64>],
) -> Vec<f64> {
    let allele_count = by_allele_and_read.len();
    let read_count = by_allele_and_read.first().map_or(0, Vec::len);
    let log10_ploidy = jmath::math::log10(ploidy as f64);
    let triallelic_possible = allele_count > 2 && ploidy > 2;
    let rescaled = if triallelic_possible {
        Some(rescaled_non_log_likelihoods(by_allele_and_read))
    } else {
        None
    };
    let genotypes = genotypes_in_canonical_order(ploidy, allele_count);
    let mut result = vec![0.0; genotypes.len()];
    for (index, genotype) in genotypes.iter().enumerate() {
        let counts = allele_counts_of(genotype);
        result[index] = match counts.len() {
            // `MathUtils.sum`, from zero.
            1 => by_allele_and_read[counts[0].0]
                .iter()
                .fold(0.0, |a, b| a + b),
            2 => {
                let (allele1, count1) = counts[0];
                let (allele2, _) = counts[1];
                let log10_count1 = jmath::math::log10(count1 as f64);
                let log10_count2 = jmath::math::log10((ploidy - count1) as f64);
                let lks1 = &by_allele_and_read[allele1];
                let lks2 = &by_allele_and_read[allele2];
                let mut sum = 0.0;
                for r in 0..read_count {
                    sum +=
                        approximate_log10_sum_log10(lks1[r] + log10_count1, lks2[r] + log10_count2);
                }
                sum - read_count as f64 * log10_ploidy
            }
            _ => {
                let (values, log10_rescaling) = rescaled
                    .as_ref()
                    .expect("three alleles over a ploidy above two");
                let mut per_read = vec![0.0; read_count];
                for &(allele, count) in &counts {
                    for (r, slot) in per_read.iter_mut().enumerate() {
                        *slot += count as f64 * values[allele][r];
                    }
                }
                let mut sum = 0.0;
                for value in &per_read {
                    sum += jmath::math::log10(*value);
                }
                sum - read_count as f64 * log10_ploidy + log10_rescaling
            }
        };
    }
    result
}

/// `rescaledNonLogLikelihoods`: each read's likelihoods less its maximum, out of log space, and the
/// sum of the maxima.
fn rescaled_non_log_likelihoods(by_allele_and_read: &[Vec<f64>]) -> (Vec<Vec<f64>>, f64) {
    let read_count = by_allele_and_read.first().map_or(0, Vec::len);
    let mut maxima = vec![f64::NEG_INFINITY; read_count];
    for row in by_allele_and_read {
        for (r, value) in row.iter().enumerate() {
            maxima[r] = maxima[r].max(*value);
        }
    }
    let rescaled: Vec<Vec<f64>> = by_allele_and_read
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .map(|(r, value)| pow10(value - maxima[r]))
                .collect()
        })
        .collect();
    let scale_factor = maxima.iter().fold(0.0, |a, b| a + b);
    (rescaled, scale_factor)
}

/// `GenotypeLikelihoods.GLsToPLs`: the distances from the most likely genotype, times ten,
/// rounded as `Math.round` does and capped at `Integer.MAX_VALUE`.
pub fn gls_to_pls(log10_likelihoods: &[f64]) -> Vec<i32> {
    let adjust = log10_likelihoods
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    log10_likelihoods
        .iter()
        .map(|value| jmath::math::round((-10.0 * (value - adjust)).min(MAX_PL)) as i32)
        .collect()
}
