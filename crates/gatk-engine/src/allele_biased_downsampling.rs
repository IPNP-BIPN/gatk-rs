//! `org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils` (GATK 4.7.0.0;
//! `AlleleBiasedDownsamplingUtils.java` is unchanged since 4.6.2.0): the callers' contamination
//! downsampling, which removes a fraction of a sample's reads so the allele counts that remain
//! look most like a clean het or hom site.
//!
//! The draws come from `Utils.getRandomDataGenerator()`, the static `Well19937c` the leveling
//! downsampler also uses, and nothing resets it between runs in one JVM: a run starts from
//! wherever the previous one left it.
//!
//! `loadContaminationFile`, the other half of the Java class, is not here.
//!
//! # Reading this file without reading Rust
//!
//! - `&[i32]` is a borrowed, read-only view of an `int[]`; `.to_vec()` copies it into a new
//!   `Vec<i32>`, which is Java's `alleleCounts.clone()`.
//! - `&mut Well19937c` is a writable borrow of the caller's generator, so every draw advances the
//!   one shared stream rather than a copy of it.
//! - `Result<Vec<T>, PermutationError>` is "the list, or the refusal `nextPermutation` throws";
//!   `?` after a call hands that refusal straight to the caller, Java's `throws` at one call site.
//! - `<T: Clone>` makes the function generic over the evidence type, as Java's `<EVIDENCE>` is;
//!   `Clone` is what lets it copy an item into the output, where Java copies the reference.

use crate::permutation::{sample_indices_without_replacement, PermutationError};
use crate::well19937c::Well19937c;

/// `targetAlleleCounts`: of the counts with `remove` reads taken from one allele, or half of them
/// from each of two, the first that scores lowest; the original counts unless one beats them.
///
/// **How**: every pair `i <= j` of alleles, in index order. `i == j` takes all of `remove` from
/// one allele; `i != j` takes `remove / 2` from each, integer division, so an odd `remove` loses
/// one read on the floor in the pair case. Each count is floored at zero.
///
/// **Why strictly less**: the reference keeps a candidate only on `score < maxScore`, so the first
/// candidate to reach the lowest score wins a tie, and the untouched counts win against every
/// candidate that merely equals them. Which allele loses reads at a tied site is decided here.
pub fn target_allele_counts(counts: &[i32], remove: i32) -> Vec<i32> {
    let mut best_score = score_allele_counts(counts);
    let mut best = counts.to_vec();
    // `i32` division truncates towards zero, as Java's `int` division does.
    let per_allele = remove / 2;
    for i in 0..counts.len() {
        for j in i..counts.len() {
            let mut candidate = counts.to_vec();
            // The reference splits these "so we don't lose on the floor (since we divided by 2)".
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
///
/// ```java
/// return Math.min(maxCount - nextBestCount + remainderCount,
///         Math.abs(nextBestCount + remainderCount));
/// ```
///
/// (One line in the reference, wrapped here.)
///
/// Fewer than two alleles score zero, which no candidate can beat, so such a site is never
/// downsampled by [`target_allele_counts`].
fn score_allele_counts(counts: &[i32]) -> i32 {
    if counts.len() < 2 {
        return 0;
    }
    let mut sorted = counts.to_vec();
    // An unstable sort is safe here: equal integers cannot be told apart, so any sort of them
    // gives the array `Arrays.sort` gives.
    sorted.sort_unstable();
    let max = sorted[sorted.len() - 1];
    let next = sorted[sorted.len() - 2];
    let remainder = sorted.iter().sum::<i32>() - max - next;
    (max - next + remainder).min((next + remainder).abs())
}

/// `selectAlleleBiasedEvidence`: the evidence to remove, given each allele's evidence in allele
/// order. The count removed is `(int) (total * fraction)`; an allele above its target loses a
/// random subset (`sampleIndicesWithoutReplacement`), or everything when the subset is all of it.
///
/// The reference takes a `Map<A, List<EVIDENCE>>` and first drops the no-call bin:
///
/// ```java
/// alleles.remove(Allele.NO_CALL);    // ignore the no-call bin
/// ```
///
/// **Possible divergence**: nothing here drops it. The reference's map has one bin per matrix
/// allele (`evidenceByBestAlleleMap` fills a `LinkedHashMap` in allele order), and so does
/// `by_allele` from its one caller, `AlleleLikelihoods::contamination_downsampling`; a matrix
/// allele equal to `Allele.NO_CALL` would be skipped upstream and counted here. Whether any
/// caller's matrix can hold one is not established.
///
/// **Why the draws are in allele order**: each allele that loses a subset spends `size - 1` draws
/// of the shared stream (see [`crate::permutation`]), so the order the alleles are visited in is
/// the order the stream is consumed, and every later draw in the run depends on it.
pub fn select_allele_biased_evidence<T: Clone>(
    by_allele: &[Vec<T>],
    fraction: f64,
    random: &mut Well19937c,
) -> Result<Vec<T>, PermutationError> {
    let total: usize = by_allele.iter().map(Vec::len).sum();
    // `as f64` widens the count to a `double`; `as i32` then truncates towards zero, which is
    // Java's `(int)` cast on a finite value.
    let remove = (total as f64 * fraction) as i32;
    let counts: Vec<i32> = by_allele.iter().map(|e| e.len() as i32).collect();
    let target = target_allele_counts(&counts, remove);
    let mut out = Vec::new();
    for (i, evidence) in by_allele.iter().enumerate() {
        // The reference's `filter(i -> alleleCounts[i] > targetAlleleCounts[i])`, which also makes
        // its `numElementsToRemove == 0` branch unreachable.
        if counts[i] <= target[i] {
            continue;
        }
        let k = counts[i] - target[i];
        // `downsampleElements`: all of it when the subset would be all of it, which takes no draw.
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
