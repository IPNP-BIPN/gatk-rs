//! `org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils` (GATK 4.6.2.0):
//! the callers' contamination downsampling, which removes a fraction of a sample's reads so the
//! allele counts that remain look most like a clean het or hom site.
//!
//! The draws come from `Utils.getRandomDataGenerator()`, the static `Well19937c` the leveling
//! downsampler also uses, and nothing resets it between runs in one JVM: a run starts from
//! wherever the previous one left it.

use crate::permutation::{sample_indices_without_replacement, PermutationError};
use crate::well19937c::Well19937c;

/// `targetAlleleCounts`: of the counts with `remove` reads taken from one allele, or half of them
/// from each of two, the first that scores lowest; the original counts unless one beats them.
pub fn target_allele_counts(counts: &[i32], remove: i32) -> Vec<i32> {
    let mut best_score = score_allele_counts(counts);
    let mut best = counts.to_vec();
    let per_allele = remove / 2;
    for i in 0..counts.len() {
        for j in i..counts.len() {
            let mut candidate = counts.to_vec();
            if i == j {
                candidate[i] = (candidate[i] - remove).max(0);
            } else {
                candidate[i] = (candidate[i] - per_allele).max(0);
                candidate[j] = (candidate[j] - per_allele).max(0);
            }
            let score = score_allele_counts(&candidate);
            if score < best_score {
                best_score = score;
                best = candidate;
            }
        }
    }
    best
}

/// `scoreAlleleCounts`: how far the counts are from a het (the two largest equal, nothing else) or
/// a hom (everything but the largest zero).
fn score_allele_counts(counts: &[i32]) -> i32 {
    if counts.len() < 2 {
        return 0;
    }
    let mut sorted = counts.to_vec();
    sorted.sort_unstable();
    let max = sorted[sorted.len() - 1];
    let next = sorted[sorted.len() - 2];
    let remainder = sorted.iter().sum::<i32>() - max - next;
    (max - next + remainder).min((next + remainder).abs())
}

/// `selectAlleleBiasedEvidence`: the evidence to remove, given each allele's evidence in allele
/// order. The count removed is `(int) (total * fraction)`; an allele above its target loses a
/// random subset (`sampleIndicesWithoutReplacement`), or everything when the subset is all of it.
pub fn select_allele_biased_evidence<T: Clone>(
    by_allele: &[Vec<T>],
    fraction: f64,
    random: &mut Well19937c,
) -> Result<Vec<T>, PermutationError> {
    let total: usize = by_allele.iter().map(Vec::len).sum();
    let remove = (total as f64 * fraction) as i32;
    let counts: Vec<i32> = by_allele.iter().map(|e| e.len() as i32).collect();
    let target = target_allele_counts(&counts, remove);
    let mut out = Vec::new();
    for (i, evidence) in by_allele.iter().enumerate() {
        if counts[i] <= target[i] {
            continue;
        }
        let k = counts[i] - target[i];
        if k >= evidence.len() as i32 {
            out.extend(evidence.iter().cloned());
        } else {
            for index in sample_indices_without_replacement(evidence.len() as i32, k, random)? {
                out.push(evidence[index as usize].clone());
            }
        }
    }
    Ok(out)
}
