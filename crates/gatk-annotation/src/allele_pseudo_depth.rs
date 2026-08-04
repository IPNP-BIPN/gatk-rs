//! Ported from `org.broadinstitute.hellbender.tools.walkers.annotator.AllelePseudoDepth`
//! (GATK 4.6.2.0).
//!
//! `DD` and `DF`: an allele depth and an allele fraction, both read off the Dirichlet posterior
//! that [`gatk_engine::somatic_likelihoods`] computes, and both emitted as **strings** rather than
//! as numbers.
//!
//! This is the annotation G1.9 was opened for. The chain under it calls `Math.exp`, whose only
//! exact port would be a transcription of GPL2-only HotSpot source (htsjdk-rs decision 0014, and
//! htsjdk-rs #71 for why no route round it is open), so the annotation was held back as
//! licence-blocked. Two measurements since have taken that argument apart:
//!
//!  * `NaturalLogUtils` on FDLIBM is within **1 ulp** of the reference (G1.9.1, #97);
//!  * the fixed point does not amplify it. 36 values, every iteration count matched, worst
//!    divergence **zero ulp** (G1.9.2, #98).
//!
//! And then the output goes through `DecimalFormat` at two and four decimals, which is about
//! twelve orders of magnitude coarser than a last bit. The rounding argument was meant to be the
//! load-bearing one; after G1.9.2 it is belt and braces.
//!
//! # `Math.pow` enters here, and it carries no bound
//!
//! [`Self::calculate_weights`] ends on `Math.pow(10, secondAdjusted)`. That is a **second**
//! unported intrinsic. It now has a bound — htsjdk-rs decision 0027 measured fdlibm's `pow` at
//! **1 ulp** from `Math.pow` over 404,964 points, the same figure 0025 got for `exp` — but this
//! port does not use fdlibm here, and the reason is measured rather than assumed.
//!
//! Switching all three of gatk-rs's `pow` call sites to [`jmath::strict_math::pow`] was tried. It
//! **broke a passing byte-identity claim**: `heterozygosity_and_mq` moved from `3.3333333333333335`
//! to `3.3333333333333344`, two ulp on the emitted value, from one ulp in `pow` amplified by the
//! arithmetic above it. The host's `powf` agrees with `Math.pow` on 99.9378% of the corpus where
//! fdlibm agrees on 98.5317%, so on the points these suites reach, the libm is simply closer.
//!
//! What that trade buys and costs is worth stating plainly, because it is the opposite of the one
//! `NaturalLogUtils` makes for `exp`. fdlibm is fixed and the host libm is whatever the machine
//! ships, so fdlibm would make the port host-independent — at the price of a suite that passes
//! today. The evidence that the libm is safe here is not an assumption: these suites pass on
//! aarch64 locally and on x86-64 in CI, which is two platforms agreeing. If a third ever disagrees,
//! the bounded alternative now exists and this is the note that says so.
//!
//! Three settings of `weightDecay` reach three different sets of intrinsics, so the conformance
//! suite has to carry all three:
//!
//! | `weightDecay` | what runs |
//! |---|---|
//! | `1.0`, the default | one `pow` per read; the second is skipped by `!= 1.0` |
//! | `0.0` | no weights at all, `calculateWeights` returns null before any `pow` |
//! | anything else | two `pow` per read |
//!
//! Worth naming in passing: `secondAdjusted` is a difference of **natural** logarithms, and it is
//! fed to `Math.pow(10, …)`. Base ten on a natural log. Transcribed, not corrected.
//!
//! # The prior array is cached, shared, and then written through
//!
//! ```java
//! final double[] prior = composePriorPseudoCounts(allelesToEmit.size());   // cached, per size
//! …
//! if (sampleMatrix.evidence().isEmpty()) {
//!     posteriors = prior;                                                  // the same array
//! }
//! …
//! if (!keepPriorInCount) {
//!     for (int i = 0; i < posteriors.length; i++) { posteriors[i] -= prior[i]; }
//! }
//! ```
//!
//! `composePriorPseudoCounts` memoises one array per allele count and hands out **that array**, not
//! a copy. On the empty-evidence branch `posteriors` is that array. So the subtraction at the end
//! is `prior[i] -= prior[i]`, which zeroes the cache — and every later genotype with the same
//! allele count then gets a prior of zeros, a different annotation entirely.
//!
//! It is a real, order-dependent state bug in the reference, it is reachable from the default
//! settings (`keepPriorInCount` is false by default), and a port that quietly copied the array
//! would be more correct and less faithful. [`AllelePseudoDepth::annotate`] reproduces it, and the
//! conformance suite calls the annotation twice to catch it.
//!
//! # `evidence().get(row)` is indexed by allele
//!
//! ```java
//! final RealMatrixChangingVisitor log10ToLnTransformer = new DefaultRealMatrixChangingVisitor() {
//!     public double visit(int row, int column, double value) {
//!         return Math.max(value, -.1 * sampleMatrixForAlleles.evidence().get(row).getMappingQuality())
//!                 * MathUtils.LOG_10;
//!     }
//! };
//! ```
//!
//! `visit` receives `(row, column)` as `(allele, read)`, and the mapping quality is looked up at
//! `evidence().get(row)` — the evidence list indexed by the **allele** number. So every entry in
//! an allele's row is floored using one read's mapping quality, the read that happens to share the
//! allele's index, and a site with more alleles than reads throws `IndexOutOfBoundsException`
//! before any of it matters.
//!
//! This only runs when the likelihoods are still in log10, which is why it has survived: the
//! HaplotypeCaller hands over natural-log likelihoods and takes the other branch.

use std::collections::HashMap;

use gatk_engine::natural_log_utils::NonFiniteSum;
use gatk_engine::somatic_likelihoods::{allele_fractions_posterior, sum};
use jmath::math::log;

use crate::decimal_format::{DEPTH_FORMAT, FRACTION_FORMAT};

/// `GATKVCFConstants.PSEUDO_DEPTH_KEY`.
pub const PSEUDO_DEPTH_KEY: &str = "DD";

/// `GATKVCFConstants.PSEUDO_FRACTION_KEY`.
pub const PSEUDO_FRACTION_KEY: &str = "DF";

/// `VCFConstants.INFO_FIELD_ARRAY_SEPARATOR`.
const ARRAY_SEPARATOR: &str = ",";

/// `MathUtils.LOG_10`, which is `Math.log(10)` and not a literal.
fn log_10() -> f64 {
    log(10.0)
}

/// One sample's slice of the likelihoods, in the shape the annotation reads it.
///
/// `log_likelihoods` is `[allele][read]`, matching the reference's `RealMatrix` whose rows are
/// alleles. `mapping_qualities` is the evidence list, which the log10 branch indexes by allele
/// rather than by read; see the module docs.
#[derive(Debug, Clone, Copy)]
pub struct SampleMatrix<'a> {
    pub log_likelihoods: &'a [Vec<f64>],
    pub mapping_qualities: &'a [i32],
    /// `AlleleLikelihoods.isNaturalLog()`. False means the values are log10 and get converted.
    pub is_natural_log: bool,
}

impl SampleMatrix<'_> {
    fn evidence_count(&self) -> usize {
        self.mapping_qualities.len()
    }
}

/// What one call put in the genotype builder, as the two strings it actually writes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PseudoDepths {
    /// `DD`, each posterior through `DecimalFormat("#.##")`, comma joined.
    pub depth: String,
    /// `DF`, each frequency through `DecimalFormat("#.####")`, comma joined.
    pub fraction: String,
}

/// The reasons a call produces nothing or fails, kept apart because they are not the same event.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PseudoDepthError {
    /// `evidence().get(row)` with `row` an allele index, on a site with more alleles than reads.
    /// The reference throws `IndexOutOfBoundsException` out of the visitor.
    EvidenceIndexOutOfBounds {
        allele: usize,
        evidence_count: usize,
    },
    /// `logSumExp` refused a non-finite accumulator somewhere inside the fixed point.
    NonFiniteSum,
}

/// The annotation, with the argument defaults and the memo that the reference keeps.
#[derive(Debug, Clone)]
pub struct AllelePseudoDepth {
    /// `--dirichlet-prior-pseudo-count`, default 1.0, a flat prior.
    pub prior: f64,
    /// `--dirichlet-keep-prior-in-count`, default false.
    pub keep_prior_in_count: bool,
    /// `--pseudo-count-weight-decay-rate`, default 1.0, `minValue = 0.0`.
    pub weight_decay: f64,
    /// `priorPseudoCounts`, the `Int2ObjectOpenHashMap` keyed on allele count. Not a cache in the
    /// harmless sense: callers get the stored array itself, and one branch writes through it.
    prior_pseudo_counts: HashMap<usize, Vec<f64>>,
}

impl Default for AllelePseudoDepth {
    fn default() -> Self {
        Self {
            prior: 1.0,
            keep_prior_in_count: false,
            weight_decay: 1.0,
            prior_pseudo_counts: HashMap::new(),
        }
    }
}

impl AllelePseudoDepth {
    /// The three arguments, with an empty memo.
    ///
    /// The memo is not a constructor parameter because it is not configuration: it is state the
    /// annotation accumulates across genotypes, and one branch writes through it. Two calls on the
    /// same object are not two independent calls, which is why this exists rather than a struct
    /// literal.
    pub fn new(prior: f64, keep_prior_in_count: bool, weight_decay: f64) -> Self {
        Self {
            prior,
            keep_prior_in_count,
            weight_decay,
            prior_pseudo_counts: HashMap::new(),
        }
    }

    /// `getKeyNames()`, in declaration order.
    pub fn key_names() -> [&'static str; 2] {
        [PSEUDO_DEPTH_KEY, PSEUDO_FRACTION_KEY]
    }

    /// `annotate(ref, vc, g, gb, likelihoods)`.
    ///
    /// `allele_indices` is `vc.getAlleles()` expressed against the matrix: the emitted alleles, in
    /// emit order, by their row in `matrix`. That is what `SubsettedLikelihoodMatrix` does, and
    /// when it is the identity the reference skips the wrapper — which changes nothing observable,
    /// because the wrapper only renumbers.
    ///
    /// `None` is the reference's silent return: a null likelihoods object, or fewer than two
    /// alleles to emit. Neither puts a key in the genotype, and an absent key is not an empty one.
    pub fn annotate(
        &mut self,
        allele_indices: &[usize],
        matrix: Option<&SampleMatrix>,
    ) -> Result<Option<PseudoDepths>, PseudoDepthError> {
        let Some(matrix) = matrix else {
            return Ok(None);
        };
        if allele_indices.len() <= 1 {
            return Ok(None);
        }

        let allele_count = allele_indices.len();
        let prior = self.compose_prior_pseudo_counts(allele_count).to_vec();

        // The empty-evidence branch aliases the cached prior array; every other branch allocates.
        // Which one ran is what decides whether the subtraction below reaches the memo.
        let aliases_the_prior = matrix.evidence_count() == 0;
        let mut posteriors = if aliases_the_prior {
            prior.clone()
        } else {
            let likelihoods = self.compose_input_likelihood_matrix(allele_indices, matrix)?;
            let weights = self.calculate_weights(&likelihoods);
            allele_fractions_posterior(&likelihoods, &prior, weights.as_deref())
                .map_err(|NonFiniteSum| PseudoDepthError::NonFiniteSum)?
                .values
        };

        // `MathUtils.normalizeSumToOne`, which allocates, so the frequencies are taken from the
        // posteriors *before* the prior is subtracted out of them.
        let frequencies = normalize_sum_to_one(&posteriors);

        if !self.keep_prior_in_count {
            for (value, prior) in posteriors.iter_mut().zip(&prior) {
                *value -= prior;
            }
            if aliases_the_prior {
                // `posteriors` was the stored array, so this subtraction just zeroed the memo for
                // every later genotype with this allele count. Reproduced, not repaired.
                if let Some(cached) = self.prior_pseudo_counts.get_mut(&allele_count) {
                    for value in cached.iter_mut() {
                        *value -= *value;
                    }
                }
            }
        }

        Ok(Some(PseudoDepths {
            depth: join(&posteriors, |value| DEPTH_FORMAT.format(value)),
            fraction: join(&frequencies, |value| FRACTION_FORMAT.format(value)),
        }))
    }

    /// `composePriorPseudoCounts(numberOfAlleles)`, memoised per allele count.
    fn compose_prior_pseudo_counts(&mut self, allele_count: usize) -> &[f64] {
        self.prior_pseudo_counts
            .entry(allele_count)
            .or_insert_with(|| vec![self.prior; allele_count])
    }

    /// `composeInputLikelihoodMatrix(likelihoods, sampleMatrixForAlleles)`.
    ///
    /// The natural-log branch returns the matrix untouched. The log10 branch floors each entry at
    /// `-0.1 * mappingQuality` before scaling, and looks that mapping quality up by **allele**
    /// index; see the module docs.
    fn compose_input_likelihood_matrix(
        &self,
        allele_indices: &[usize],
        matrix: &SampleMatrix,
    ) -> Result<Vec<Vec<f64>>, PseudoDepthError> {
        let rows: Vec<Vec<f64>> = allele_indices
            .iter()
            .map(|index| matrix.log_likelihoods[*index].clone())
            .collect();
        if matrix.is_natural_log {
            return Ok(rows);
        }

        let log_10 = log_10();
        let mut converted = Vec::with_capacity(rows.len());
        for (allele, row) in rows.into_iter().enumerate() {
            let Some(quality) = matrix.mapping_qualities.get(allele) else {
                return Err(PseudoDepthError::EvidenceIndexOutOfBounds {
                    allele,
                    evidence_count: matrix.evidence_count(),
                });
            };
            let floor = -0.1 * f64::from(*quality);
            converted.push(row.into_iter().map(|v| v.max(floor) * log_10).collect());
        }
        Ok(converted)
    }

    /// `calculateWeights(lkMatrix)`, one weight per read.
    ///
    /// `None` is the reference's `null`, which `alleleFractionsPosterior` reads as "no weights".
    /// A negative decay throws; the argument is declared `minValue = 0.0`, so the parser refuses it
    /// first and this is unreachable through the command line.
    fn calculate_weights(&self, likelihoods: &[Vec<f64>]) -> Option<Vec<f64>> {
        if self.weight_decay == 0.0 {
            return None;
        }
        assert!(
            self.weight_decay >= 0.0,
            "the weight decay must be 0 or greater"
        );
        let read_count = likelihoods.first().map_or(0, Vec::len);
        let mut weights = Vec::with_capacity(read_count);
        for read in 0..read_count {
            // `best` starts at the first allele's entry and `secondBest` at negative infinity, so
            // the loop runs from allele 1 and a single-allele matrix would leave secondBest there.
            let mut best = likelihoods[0][read];
            let mut second_best = f64::NEG_INFINITY;
            for row in likelihoods.iter().skip(1) {
                let value = row[read];
                if value > best {
                    second_best = best;
                    best = value;
                } else if value > second_best {
                    second_best = value;
                }
            }
            let second_adjusted = second_best - best;
            // Base ten, on a difference of natural logs. Transcribed as written.
            let mut weight = 1.0 - 10f64.powf(second_adjusted);
            if self.weight_decay != 1.0 {
                weight = weight.powf(self.weight_decay);
            }
            weights.push(weight);
        }
        Some(weights)
    }
}

/// `MathUtils.normalizeSumToOne(array)`. The sum accumulates in index order, and the division is
/// `applyToArray`, which allocates.
fn normalize_sum_to_one(values: &[f64]) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let total = sum(values);
    values.iter().map(|value| value / total).collect()
}

/// `Arrays.stream(…).mapToObj(format).collect(joining(","))`.
fn join(values: &[f64], format: impl Fn(f64) -> String) -> String {
    values
        .iter()
        .map(|value| format(*value))
        .collect::<Vec<_>>()
        .join(ARRAY_SEPARATOR)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix<'a>(
        rows: &'a [Vec<f64>],
        qualities: &'a [i32],
        natural_log: bool,
    ) -> SampleMatrix<'a> {
        SampleMatrix {
            log_likelihoods: rows,
            mapping_qualities: qualities,
            is_natural_log: natural_log,
        }
    }

    #[test]
    fn a_null_matrix_and_a_single_allele_both_emit_nothing() {
        let mut annotation = AllelePseudoDepth::default();
        assert_eq!(annotation.annotate(&[0, 1], None), Ok(None));
        let rows = vec![vec![-0.1], vec![-3.0]];
        let matrix = matrix(&rows, &[60], true);
        assert_eq!(annotation.annotate(&[0], Some(&matrix)), Ok(None));
        assert_eq!(annotation.annotate(&[], Some(&matrix)), Ok(None));
    }

    /// The state bug, which is the reason the suite calls the annotation twice.
    #[test]
    fn an_empty_evidence_call_zeroes_the_memo_for_every_later_call() {
        let mut annotation = AllelePseudoDepth::default();
        let empty: Vec<Vec<f64>> = vec![Vec::new(), Vec::new()];
        let empty_matrix = matrix(&empty, &[], true);

        // With no evidence the posteriors are the prior, so the depths are the prior minus itself.
        let first = annotation
            .annotate(&[0, 1], Some(&empty_matrix))
            .expect("no failure")
            .expect("two alleles");
        assert_eq!(first.depth, "0,0");
        assert_eq!(first.fraction, "0.5,0.5");

        // And the memo is now zeros, which the next call inherits as its prior.
        assert_eq!(
            annotation.prior_pseudo_counts.get(&2),
            Some(&vec![0.0, 0.0]),
            "the cached prior was written through"
        );
    }

    /// Keeping the prior in the count is the switch that stops the write-through, because the
    /// subtraction is what does the damage.
    #[test]
    fn keeping_the_prior_leaves_the_memo_alone() {
        let mut annotation = AllelePseudoDepth::new(1.0, true, 1.0);
        let empty: Vec<Vec<f64>> = vec![Vec::new(), Vec::new()];
        let empty_matrix = matrix(&empty, &[], true);
        let first = annotation
            .annotate(&[0, 1], Some(&empty_matrix))
            .expect("no failure")
            .expect("two alleles");
        assert_eq!(first.depth, "1,1");
        assert_eq!(
            annotation.prior_pseudo_counts.get(&2),
            Some(&vec![1.0, 1.0])
        );
    }

    /// More alleles than reads, on the log10 branch, reaches the reference's exception.
    #[test]
    fn the_log10_branch_indexes_the_evidence_by_allele() {
        let mut annotation = AllelePseudoDepth::default();
        let rows = vec![vec![-0.1, -0.2], vec![-3.0, -0.1], vec![-6.0, -5.0]];
        let qualities = [60, 20];
        let log10 = matrix(&rows, &qualities, false);
        assert_eq!(
            annotation.annotate(&[0, 1, 2], Some(&log10)),
            Err(PseudoDepthError::EvidenceIndexOutOfBounds {
                allele: 2,
                evidence_count: 2
            })
        );
        // Two alleles and two reads is within bounds, and quietly uses read 0's and read 1's
        // mapping qualities as the floors for the two allele rows.
        let mut annotation = AllelePseudoDepth::default();
        assert!(annotation.annotate(&[0, 1], Some(&log10)).is_ok());
    }

    /// The three weight regimes, which reach three different sets of intrinsics.
    #[test]
    fn the_decay_rate_decides_how_many_powers_run() {
        let likelihoods = vec![vec![-0.01, -4.0], vec![-4.0, -0.01]];
        let default = AllelePseudoDepth::default();
        let weights = default.calculate_weights(&likelihoods).expect("weights");
        assert_eq!(weights.len(), 2);

        let none = AllelePseudoDepth::new(1.0, false, 0.0);
        assert_eq!(none.calculate_weights(&likelihoods), None);

        let squared = AllelePseudoDepth::new(1.0, false, 2.0);
        let squared = squared.calculate_weights(&likelihoods).expect("weights");
        for (plain, squared) in weights.iter().zip(&squared) {
            assert_eq!(*squared, plain.powf(2.0));
        }
    }
}
