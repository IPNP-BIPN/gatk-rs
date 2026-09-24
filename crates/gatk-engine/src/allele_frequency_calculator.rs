//! `AlleleFrequencyCalculator` and `AFCalculationResult`: the Dirichlet fit behind QUAL.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AlleleFrequencyCalculator`,
//! `org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculationResult`,
//! `org.broadinstitute.hellbender.utils.Dirichlet` and
//! `org.broadinstitute.hellbender.utils.GenotypeUtils` (GATK 4.6.2.0).
//!
//! What the genotyping engine asks before it calls anything: the probability that the site holds
//! no variant, which becomes QUAL, the probability that each alternate is absent, which decides
//! the alleles kept, and the integer allele counts that become MLEAC.
//!
//! # The fit
//!
//! Allele frequencies start flat, one over the allele count. Each round turns them into genotype
//! posteriors per sample, sums those into effective allele counts, adds the prior pseudocounts and
//! takes the Dirichlet mean as the next frequencies. It stops once no count moves by more than
//! 0.1. The prior only enters from the second round, so one strong het is not drowned by the
//! reference pseudocount before it has any counts of its own.
//!
//! # Which arithmetic
//!
//! `Math.log10` and `Math.pow` are HotSpot intrinsics: [`jmath::math::log10`] and
//! [`crate::math_utils::pow10`] stand in for them, as everywhere else in this crate. The
//! combination count goes through commons-math's `CombinatoricsUtils.factorialLog`, which is
//! `FastMath.log` of an exact factorial below 21, so it is [`jmath::fast_math::log`] here.

use crate::genotype_index;
use crate::math_utils::{log10_sum_log10, pow10};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::Genotype;

/// `THRESHOLD_FOR_ALLELE_COUNT_CONVERGENCE`.
const CONVERGENCE: f64 = 0.1;

/// `AFCalculationResult.EPSILON`.
const EPSILON: f64 = 1.0e-10;

/// `GenotypeUtils.PLOIDY_2_HOM_VAR_SCALE_FACTOR`.
const PLOIDY_2_HOM_VAR_SCALE_FACTOR: i32 = 20;

/// `GenotypeCalculationArgumentCollection`'s defaults: `HomoSapiensConstants.SNP_HETEROZYGOSITY`,
/// `INDEL_HETEROZYGOSITY`, a standard deviation of 0.01 and `DEFAULT_PLOIDY`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Priors {
    pub snp_heterozygosity: f64,
    pub indel_heterozygosity: f64,
    pub heterozygosity_standard_deviation: f64,
    pub sample_ploidy: usize,
}

impl Default for Priors {
    fn default() -> Self {
        Priors {
            snp_heterozygosity: 1e-3,
            indel_heterozygosity: 1.0 / 8000.0,
            heterozygosity_standard_deviation: 0.01,
            sample_ploidy: 2,
        }
    }
}

/// What `calculate` refuses, in the reference's classes and words.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AfCalcError {
    /// `IllegalStateException`.
    IllegalState(String),
    /// `IllegalArgumentException`.
    IllegalArgument(String),
}

impl AfCalcError {
    /// The exception's fully qualified class.
    pub fn java_class(&self) -> &'static str {
        match self {
            AfCalcError::IllegalState(_) => "java.lang.IllegalStateException",
            AfCalcError::IllegalArgument(_) => "java.lang.IllegalArgumentException",
        }
    }

    /// The exception's message.
    pub fn message(&self) -> &str {
        match self {
            AfCalcError::IllegalState(message) | AfCalcError::IllegalArgument(message) => message,
        }
    }
}

/// `AlleleFrequencyCalculator`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlleleFrequencyCalculator {
    ref_pseudocount: f64,
    snp_pseudocount: f64,
    indel_pseudocount: f64,
    default_ploidy: usize,
}

/// `AFCalculationResult`.
#[derive(Debug, Clone, PartialEq)]
pub struct AfCalculationResult {
    /// `getAlleleCountsOfMLE`: one per alternate, in the site's order.
    pub allele_counts_of_mle: Vec<i32>,
    /// `log10ProbOnlyRefAlleleExists`.
    pub log10_posterior_of_no_variant: f64,
    /// The log10 posterior that each alternate is absent, in the site's order.
    pub log10_p_ref_by_allele: Vec<f64>,
}

impl AfCalculationResult {
    /// `log10ProbVariantPresent`: `MathUtils.log10OneMinusPow10` of the no-variant posterior.
    pub fn log10_prob_variant_present(&self) -> f64 {
        log10_one_minus_pow10(self.log10_posterior_of_no_variant)
    }

    /// `passesThreshold(allele, phredScaleQualThreshold)`, for the alternate at `alt_index`.
    pub fn passes_threshold(&self, alt_index: usize, phred_scale_qual_threshold: f64) -> bool {
        self.log10_p_ref_by_allele[alt_index] + EPSILON < phred_scale_qual_threshold * -0.1
    }
}

/// `MathUtils.log10OneMinusPow10`.
pub fn log10_one_minus_pow10(a: f64) -> f64 {
    if a > 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return f64::NEG_INFINITY;
    }
    let ln10 = jmath::math::log(10.0);
    log1mexp_with_host_exp(a * ln10) * (1.0 / ln10)
}

/// `NaturalLogUtils.log1mexp`, with the host's `exp` standing in for `Math.exp`.
///
/// The shared [`crate::natural_log_utils::log1mexp`] calls fdlibm, which is what the somatic
/// likelihood suites were measured against. Here the golden says otherwise: a no-variant posterior
/// of about 1e-19 goes through `log1p(-exp(a))` with `log1p` returning its tiny argument
/// unchanged, so `exp(a)` reaches the output as it is, and fdlibm's answer is two ulp from the
/// reference's where the host libm's is identical. It is the same trade [`pow10`] makes for
/// `Math.pow`: the intrinsic cannot be transcribed (htsjdk-rs decision 0014), and the host is the
/// closer stand-in on the points measured.
fn log1mexp_with_host_exp(a: f64) -> f64 {
    if a > 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return f64::NEG_INFINITY;
    }
    if a < jmath::math::log(0.5) {
        (-std::hint::black_box(a).exp()).ln_1p()
    } else {
        jmath::math::log(-(a.exp_m1()))
    }
}

/// `MathUtils.log10SumLog10(a, b)`, the two-argument form, which is not the array one.
fn log10_sum_log10_pair(a: f64, b: f64) -> f64 {
    if a > b {
        a + jmath::math::log10(1.0 + pow10(b - a))
    } else {
        b + jmath::math::log10(1.0 + pow10(a - b))
    }
}

/// `CombinatoricsUtils.factorialLog`: `FastMath.log` of the exact factorial below 21, and the sum
/// of the logarithms of `2..=n` from there.
fn factorial_log(n: usize) -> f64 {
    if n < 21 {
        let factorial: u64 = (2..=n as u64).product();
        return jmath::fast_math::log(factorial as f64);
    }
    (2..=n).map(|i| jmath::fast_math::log(i as f64)).sum()
}

/// `GenotypeAlleleCounts.log10CombinationCount`: the log10 of the number of orderings of the
/// genotype's alleles, `ploidy! / prod(count!)` taken in natural logs and converted.
fn log10_combination_count(ploidy: usize, counts: &[(usize, usize)]) -> f64 {
    let ln = factorial_log(ploidy) - counts.iter().map(|&(_, c)| factorial_log(c)).sum::<f64>();
    ln * jmath::math::log10(std::f64::consts::E)
}

/// `GenotypeUtils.genotypeIsUsableForAFCalculation`.
fn is_usable(g: &Genotype) -> bool {
    g.pl.is_some() || (g.is_hom_ref() && g.gq.is_some() && g.ploidy() == 2)
}

/// `GenotypeUtils.makeApproximateDiploidLog10LikelihoodsFromGQ`: 0 for the hom-ref, GQ for every
/// genotype holding the reference, twenty GQ for every other, and back from phred.
fn approximate_diploid_likelihoods(gq: i32, allele_count: usize) -> Vec<f64> {
    genotype_index::genotypes_in_canonical_order(2, allele_count)
        .iter()
        .enumerate()
        .map(|(index, genotype)| {
            let pl = if index == 0 {
                0
            } else if genotype.contains(&0) {
                gq
            } else {
                PLOIDY_2_HOM_VAR_SCALE_FACTOR * gq
            };
            pl as f64 / -10.0
        })
        .collect()
}

/// `MathUtils.normalizeLog10(array)`, the log-space normalisation.
fn normalize_log10(array: &[f64]) -> Vec<f64> {
    let sum = log10_sum_log10(array);
    array.iter().map(|x| x - sum).collect()
}

/// `Dirichlet.log10MeanWeights`.
fn log10_mean_weights(alpha: &[f64]) -> Vec<f64> {
    let sum: f64 = alpha.iter().sum();
    alpha.iter().map(|x| jmath::math::log10(x / sum)).collect()
}

fn is_span_del(allele: &Allele) -> bool {
    !allele.is_reference() && allele.display_string() == htsjdk_vcf::allele::SPAN_DEL_STRING
}

/// One genotype shape, with the combination counts computed once.
struct Shape {
    genotypes: Vec<Vec<(usize, usize)>>,
    log10_combinations: Vec<f64>,
}

impl Shape {
    fn new(ploidy: usize, allele_count: usize) -> Shape {
        let genotypes: Vec<Vec<(usize, usize)>> =
            genotype_index::genotypes_in_canonical_order(ploidy, allele_count)
                .iter()
                .map(|g| genotype_index::allele_counts_of(g))
                .collect();
        let log10_combinations = genotypes
            .iter()
            .map(|counts| log10_combination_count(ploidy, counts))
            .collect();
        Shape {
            genotypes,
            log10_combinations,
        }
    }
}

impl AlleleFrequencyCalculator {
    /// The constructor, from the three pseudocounts and the default ploidy.
    pub fn new(
        ref_pseudocount: f64,
        snp_pseudocount: f64,
        indel_pseudocount: f64,
        default_ploidy: usize,
    ) -> Self {
        AlleleFrequencyCalculator {
            ref_pseudocount,
            snp_pseudocount,
            indel_pseudocount,
            default_ploidy,
        }
    }

    /// `makeCalculator(GenotypeCalculationArgumentCollection)`.
    pub fn make_calculator(priors: &Priors) -> Self {
        let sd = priors.heterozygosity_standard_deviation;
        let ref_pseudocount = priors.snp_heterozygosity / (sd * sd);
        AlleleFrequencyCalculator::new(
            ref_pseudocount,
            priors.snp_heterozygosity * ref_pseudocount,
            priors.indel_heterozygosity * ref_pseudocount,
            priors.sample_ploidy,
        )
    }

    /// `getPloidy`.
    pub fn ploidy(&self) -> usize {
        self.default_ploidy
    }

    /// `calculate(vc)`: the site's alleles, reference first, and its genotypes.
    pub fn calculate(
        &self,
        contig: &str,
        start: i64,
        alleles: &[Allele],
        genotypes: &[Genotype],
    ) -> Result<AfCalculationResult, AfCalcError> {
        self.calculate_with_ploidy(contig, start, alleles, genotypes, self.default_ploidy)
    }

    /// `calculate(vc, defaultPloidy)`.
    pub fn calculate_with_ploidy(
        &self,
        contig: &str,
        start: i64,
        alleles: &[Allele],
        genotypes: &[Genotype],
        default_ploidy: usize,
    ) -> Result<AfCalculationResult, AfCalcError> {
        if !genotypes.iter().any(|g| g.pl.is_some()) {
            return Err(AfCalcError::IllegalState(format!(
                "VariantContext  at {contig}:{start}must contain at least one genotype with \
                 likelihoods -- did this VC exceed the max number of alt alleles?"
            )));
        }
        if alleles.len() <= 1 {
            return Err(AfCalcError::IllegalArgument(format!(
                "VariantContext  at {contig}:{start}has only a single reference allele, but \
                 getLog10PNonRef requires at least alternate allele"
            )));
        }
        self.fit(alleles, genotypes, default_ploidy)
    }

    /// `log10NormalizedGenotypePosteriors`.
    fn posteriors(
        &self,
        g: &Genotype,
        shape: &Shape,
        log10_allele_frequencies: &[f64],
    ) -> Result<Vec<f64>, AfCalcError> {
        let likelihoods: Vec<f64> = match &g.pl {
            Some(pl) => pl.iter().map(|&p| p as f64 / -10.0).collect(),
            None => {
                // Only a hom-ref or a no-call reaches here: `is_usable` let nothing else through.
                approximate_diploid_likelihoods(g.gq.unwrap_or(0), log10_allele_frequencies.len())
            }
        };
        if likelihoods.len() != shape.genotypes.len() {
            return Err(AfCalcError::IllegalState(
                "Ploidy, allele count, and genotype likelihoods are inconsistent".to_string(),
            ));
        }
        let unnormalized: Vec<f64> = shape
            .genotypes
            .iter()
            .enumerate()
            .map(|(index, counts)| {
                let prior: f64 = counts
                    .iter()
                    .map(|&(allele, count)| count as f64 * log10_allele_frequencies[allele])
                    .sum();
                shape.log10_combinations[index] + likelihoods[index] + prior
            })
            .collect();
        Ok(normalize_log10(&unnormalized))
    }

    /// `effectiveAlleleCounts`.
    fn effective_allele_counts(
        &self,
        genotypes: &[Genotype],
        shapes: &mut Shapes,
        log10_allele_frequencies: &[f64],
    ) -> Result<Vec<f64>, AfCalcError> {
        let allele_count = log10_allele_frequencies.len();
        let mut log10_result = vec![f64::NEG_INFINITY; allele_count];
        for g in genotypes.iter().filter(|g| is_usable(g)) {
            let shape = shapes.get(g.ploidy(), allele_count);
            let posteriors = self.posteriors(g, shape, log10_allele_frequencies)?;
            for (index, counts) in shape.genotypes.iter().enumerate() {
                for &(allele, count) in counts {
                    log10_result[allele] = log10_sum_log10_pair(
                        log10_result[allele],
                        posteriors[index] + jmath::math::log10(count as f64),
                    );
                }
            }
        }
        Ok(log10_result.into_iter().map(pow10).collect())
    }

    /// The private `calculate`, which does the fit.
    fn fit(
        &self,
        alleles: &[Allele],
        genotypes: &[Genotype],
        default_ploidy: usize,
    ) -> Result<AfCalculationResult, AfCalcError> {
        let allele_count = alleles.len();
        let ref_length = alleles[0].len();
        let prior_pseudocounts: Vec<f64> = alleles
            .iter()
            .map(|a| {
                if a.is_reference() {
                    self.ref_pseudocount
                } else if a.len() == ref_length {
                    self.snp_pseudocount
                } else {
                    self.indel_pseudocount
                }
            })
            .collect();

        let mut shapes = Shapes::default();
        let mut allele_counts = vec![0.0; allele_count];
        let flat = -jmath::math::log10(allele_count as f64);
        let mut log10_allele_frequencies = vec![flat; allele_count];
        let mut maximum_difference = f64::INFINITY;
        while maximum_difference > CONVERGENCE {
            let new_counts =
                self.effective_allele_counts(genotypes, &mut shapes, &log10_allele_frequencies)?;
            // `Arrays.stream(...).map(Math::abs).max()`, which is NaN-propagating only through
            // `Math.max`'s own rules; every count here is finite.
            maximum_difference = allele_counts
                .iter()
                .zip(&new_counts)
                .map(|(old, new)| (old - new).abs())
                .fold(f64::NEG_INFINITY, f64::max);
            allele_counts = new_counts;
            let posterior_pseudocounts: Vec<f64> = prior_pseudocounts
                .iter()
                .zip(&allele_counts)
                .map(|(p, c)| p + c)
                .collect();
            log10_allele_frequencies = log10_mean_weights(&posterior_pseudocounts);
        }

        let span_del_index = alleles.iter().position(is_span_del);
        let mut log10_p_of_zero_counts_by_allele = vec![0.0; allele_count];
        let mut log10_p_no_variant = 0.0;
        for g in genotypes.iter().filter(|g| is_usable(g)) {
            let ploidy = if g.ploidy() == 0 {
                default_ploidy
            } else {
                g.ploidy()
            };
            let shape = shapes.get(g.ploidy(), allele_count);
            let posteriors = self.posteriors(g, shape, &log10_allele_frequencies)?;
            match span_del_index {
                None => log10_p_no_variant += posteriors[0],
                Some(span_del) => {
                    let non_variant: Vec<f64> = (0..=ploidy)
                        .map(|n| {
                            let index = genotype_index::allele_counts_to_index(&[
                                0,
                                ploidy - n,
                                span_del,
                                n,
                            ])
                            .expect("an even-length count array");
                            posteriors[index]
                        })
                        .collect();
                    log10_p_no_variant += java_min_zero(log10_sum_log10(&non_variant));
                }
            }
            if allele_count == 2 && span_del_index.is_none() {
                continue;
            }
            let shape = shapes.get(ploidy, allele_count);
            let mut absent: Vec<Vec<f64>> = vec![Vec::new(); allele_count];
            for (index, counts) in shape.genotypes.iter().enumerate() {
                for (allele, buffer) in absent.iter_mut().enumerate() {
                    if !counts.iter().any(|&(present, _)| present == allele) {
                        buffer.push(posteriors[index]);
                    }
                }
            }
            for (allele, buffer) in absent.iter().enumerate() {
                log10_p_of_zero_counts_by_allele[allele] += java_min_zero(log10_sum_log10(buffer));
            }
        }
        if allele_count == 2 && span_del_index.is_none() {
            log10_p_of_zero_counts_by_allele[1] = log10_p_no_variant;
        }

        let integer_counts: Vec<i32> = allele_counts
            .iter()
            .map(|&x| jmath::math::round(x) as i32)
            .collect();
        Ok(AfCalculationResult {
            allele_counts_of_mle: integer_counts[1..].to_vec(),
            log10_posterior_of_no_variant: log10_p_no_variant,
            log10_p_ref_by_allele: log10_p_of_zero_counts_by_allele[1..].to_vec(),
        })
    }
}

/// `Math.min(0, x)`, which keeps a NaN and prefers -0.0 over 0.0.
fn java_min_zero(x: f64) -> f64 {
    crate::math_utils::java_min(0.0, x)
}

/// The shapes met so far, keyed by ploidy, for one allele count.
#[derive(Default)]
struct Shapes {
    by_ploidy: Vec<(usize, usize, Shape)>,
}

impl Shapes {
    fn get(&mut self, ploidy: usize, allele_count: usize) -> &Shape {
        let position = self
            .by_ploidy
            .iter()
            .position(|(p, a, _)| *p == ploidy && *a == allele_count);
        let position = match position {
            Some(position) => position,
            None => {
                self.by_ploidy
                    .push((ploidy, allele_count, Shape::new(ploidy, allele_count)));
                self.by_ploidy.len() - 1
            }
        };
        &self.by_ploidy[position].2
    }
}
