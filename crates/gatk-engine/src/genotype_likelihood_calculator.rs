//! Ported from `org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator`
//! and `IndependentSampleGenotypesModel` (GATK 4.7.0.0; both `.java` files are unchanged since
//! 4.6.2.0): one sample's genotype likelihoods from its read-by-allele likelihoods.
//!
//! Each genotype, in canonical order, is the product over reads of the average, over the genotype's
//! allele copies, of the read's likelihood for the allele, computed three ways by how many distinct
//! alleles the genotype holds:
//!
//!  * one: the sum of that allele's log10 likelihoods;
//!  * two: per read `approximateLog10SumLog10` of each allele's log10 likelihood plus the log10 of
//!    its copy count, summed, less `reads * log10(ploidy)`. The approximation is the table-based
//!    one, so these genotypes carry its rounding;
//!  * three or more (possible only with more than two alleles and a ploidy above two): the
//!    likelihoods are rescaled per read by the read's maximum and taken out of log space once,
//!    each read's copy-weighted sum logged, and the rescaling added back.
//!
//! The Java hands the vector to `GenotypeLikelihoods.fromLog10Likelihoods`, which keeps it as log10
//! values (`log10Likelihoods = asVector;`); PLs are derived from it when asked, by `GLsToPLs`,
//! which is [`gls_to_pls`] below.
//!
//! `IndependentSampleGenotypesModel.calculateLikelihoods` is only the per-sample loop around that
//! call, after remapping the alleles with `AlleleLikelihoodMatrixMapper`. Its callers in
//! `gatk-tools` run that loop themselves (`hc_genotyping`), so no function here carries its name.
//!
//! # The order of the floating-point operations
//!
//! A PL is a rounded difference of these sums, so a last-bit change can move one across a rounding
//! boundary. Each case adds in the order the Java does:
//!
//!  * one allele: `MathUtils.sum(row)`, which is `double s = 0.0;` then `s += v;` left to right;
//!  * two alleles: the per-read terms are accumulated by `IndexRange.sum` (`double result = 0;`,
//!    then `result += lambda.applyAsDouble(i);`), and the average over ploidy is taken **once, at
//!    the end**, as one subtraction of `readCount` times `Math.log10` of the ploidy, not as a
//!    per-read subtraction. The second count is `ploidy - count1`, not a second lookup;
//!  * three or more: a per-read buffer starts at zero and each distinct allele, in ascending order,
//!    adds `perReadBuffer[r] += f * rescaledNonLogLikelihoods[a][r]` (an `int` count times a
//!    `double`), so each read's sum is accumulated allele by allele. Then the `Math.log10` of each
//!    entry is summed from zero, and the expression is, left to right, that sum, minus `readCount`
//!    times `Math.log10` of the ploidy, plus the rescaling
//!    (`GenotypeLikelihoodCalculator.java:106`).
//!
//! The constants and functions behind them: `Math.log10` is `jmath::math::log10` (correctly
//! rounded, as `jmath` measured it), `Math.pow(10.0, x)` is [`crate::math_utils::pow10`] (the
//! deferred intrinsic, see that module), and `MathUtils.approximateLog10SumLog10` with its
//! `JacobianLogTable` is `pair_hmm::approximate_log10_sum_log10`.
//!
//! # Reading this file without reading Rust
//!
//! - `&[Vec<f64>]` is a borrowed, read-only view of a list of rows, Java's `double[][]`, indexed
//!   `[allele][read]`;
//! - `Option<T>` is "a value or nothing", `Some(v)` or `None`: here it stands for the Java's
//!   `null` when the rescaled matrix is not needed;
//! - `match x { 1 => .., 2 => .., _ => .. }` is a `switch` whose `_` arm is the `default`; the
//!   value of the chosen arm is the value of the whole `match`;
//! - `.fold(0.0, |a, b| a + b)` adds left to right from `0.0`, which states the Java's starting
//!   value rather than leaving it to the library's `sum()`;
//! - `x as f64` converts an integer to `double`, exactly here, as Java's promotion does.

use crate::genotype_index::{allele_counts_of, genotypes_in_canonical_order};
use crate::math_utils::pow10;
use crate::pair_hmm::approximate_log10_sum_log10;

/// `GenotypeLikelihoods.MAX_PL` (`public static final int MAX_PL = Integer.MAX_VALUE;`), held as a
/// `double` because the Java compares it with one in `Math.min`, which promotes it.
const MAX_PL: f64 = i32::MAX as f64;

/// `computeLog10GenotypeLikelihoods(ploidy, log10AlleleLikelihoods)`, the matrix being indexed
/// `[allele][read]`.
///
/// **What**: one log10 likelihood per genotype, in canonical order, from the class javadoc's model
/// `product_{all reads} [[sum_{alleles in genotype} Prob(read | allele)]/ploidy]`.
///
/// **How**: the three cases and their order of operations are in the module note. The two-allele
/// case stays in log space with the table approximation; only the three-allele case leaves it.
pub fn compute_log10_genotype_likelihoods(
    ploidy: usize,
    by_allele_and_read: &[Vec<f64>],
) -> Vec<f64> {
    let allele_count = by_allele_and_read.len();
    // The Java's `evidenceCount()`. `.first()` is the first row or `None`, and `.map_or(0, ..)`
    // answers 0 reads when there are no rows at all.
    let read_count = by_allele_and_read.first().map_or(0, Vec::len);
    // The Java evaluates `Math.log10` of the ploidy afresh in each genotype's expression; the
    // value is the same every time, so it is computed once here.
    let log10_ploidy = jmath::math::log10(ploidy as f64);
    // The Java's `triallelicGenotypesPossible = alleleCount > 2 && ploidy > 2`: only then can a
    // genotype hold three distinct alleles, and only then is the rescaled matrix built.
    let triallelic_possible = allele_count > 2 && ploidy > 2;
    let rescaled = if triallelic_possible {
        Some(rescaled_non_log_likelihoods(by_allele_and_read))
    } else {
        None
    };
    let genotypes = genotypes_in_canonical_order(ploidy, allele_count);
    let mut result = vec![0.0; genotypes.len()];
    for (index, genotype) in genotypes.iter().enumerate() {
        // The Java's `gac.distinctAlleleCount()`, `alleleIndexAt(rank)` and `alleleCountAt(rank)`
        // are this list's length and its pairs' two halves.
        let counts = allele_counts_of(genotype);
        result[index] = match counts.len() {
            // `MathUtils.sum`, from zero.
            1 => by_allele_and_read[counts[0].0]
                .iter()
                .fold(0.0, |a, b| a + b),
            2 => {
                // `let (a, b) = pair` unpacks a pair into two names; `_` discards the half not
                // used, because the Java takes the second count as `ploidy - count1`.
                let (allele1, count1) = counts[0];
                let (allele2, _) = counts[1];
                let log10_count1 = jmath::math::log10(count1 as f64);
                let log10_count2 = jmath::math::log10((ploidy - count1) as f64);
                let lks1 = &by_allele_and_read[allele1];
                let lks2 = &by_allele_and_read[allele2];
                // `IndexRange.sum`: from zero, one read at a time, in read order.
                let mut sum = 0.0;
                for r in 0..read_count {
                    sum +=
                        approximate_log10_sum_log10(lks1[r] + log10_count1, lks2[r] + log10_count2);
                }
                // The average over ploidy, once for all reads.
                sum - read_count as f64 * log10_ploidy
            }
            // Three or more distinct alleles. Ploidy zero also lands here, with no alleles at all;
            // the Java then reaches `Arrays.fill(perReadBuffer,0, readCount, 0);` with
            // `perReadBuffer` null, and this `expect` stops the program instead.
            _ => {
                // `.as_ref()` borrows the pair inside the `Option` without taking it, and
                // `.expect(..)` unwraps it or panics with that message.
                let (values, log10_rescaling) = rescaled
                    .as_ref()
                    .expect("three alleles over a ploidy above two");
                // The Java's `perReadBuffer`, zero-filled for each genotype.
                let mut per_read = vec![0.0; read_count];
                // Allele outer, read inner, as `forEachAlleleIndexAndCount` with an
                // `IndexRange.forEach` inside. `for &(allele, count) in &counts` walks the pairs
                // without consuming the list, unpacking each; `iter_mut()` gives a writable
                // reference per entry, and `*slot +=` adds through it.
                for &(allele, count) in &counts {
                    for (r, slot) in per_read.iter_mut().enumerate() {
                        *slot += count as f64 * values[allele][r];
                    }
                }
                let mut sum = 0.0;
                for value in &per_read {
                    sum += jmath::math::log10(*value);
                }
                // Left to right, as the Java's one expression: `(sum - reads * log10 ploidy)`,
                // then plus the rescaling.
                sum - read_count as f64 * log10_ploidy + log10_rescaling
            }
        };
    }
    result
}

/// `rescaledNonLogLikelihoods`: each read's likelihoods less its maximum, out of log space, and the
/// sum of the maxima.
///
/// **How**, in the Java's order: the per-read maximum over all alleles, starting from negative
/// infinity (`FastMath.max`, allele outer, read inner); each entry less its read's maximum, then
/// `Math.pow(10.0, x)` of the difference; and the scale factor `MathUtils.sum(perReadMaxima, 0,
/// readCount)`, from zero in read order.
///
/// **Why a new matrix**: the Java subtracts in place (`log10LikelihoodsByAlleleAndRead[a][r] -=
/// perReadMaxima[r];`) in the array it got from `asRealMatrix().getData()`. This returns a fresh
/// one and leaves the caller's rows untouched; the value of each entry is the same subtraction
/// followed by the same `pow`. The tuple `(matrix, scale)` is the Java's `Pair` of the two.
fn rescaled_non_log_likelihoods(by_allele_and_read: &[Vec<f64>]) -> (Vec<Vec<f64>>, f64) {
    let read_count = by_allele_and_read.first().map_or(0, Vec::len);
    let mut maxima = vec![f64::NEG_INFINITY; read_count];
    for row in by_allele_and_read {
        for (r, value) in row.iter().enumerate() {
            maxima[r] = maxima[r].max(*value);
        }
    }
    // `.iter().map(..).collect()` is Java's stream, map, collect: one new row per allele, each
    // entry `pow10(value - maximum)`.
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
///
/// ```java
/// pls[i] = (int) Math.round(Math.min(-10 * (GLs[i] - adjust), MAX_PL));
/// ```
///
/// The cap comes before the rounding. `adjust` is `maxPL(GLs)`, a `Math.max` fold from negative
/// infinity. `jmath::math::round` returns an `i64`, as `Math.round` returns a `long`, and `as i32`
/// keeps its low 32 bits, as Java's `(int)` does; after the cap the value always fits anyway.
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
