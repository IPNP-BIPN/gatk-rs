//! `CreateSomaticPanelOfNormals`: which sites recur across normals often enough to be artefacts.
//!
//! A site is kept when enough samples carry an alternate that is more likely to be an artefact than
//! to be germline, and it is written with the fraction of samples that carried it and a beta shape
//! fitted to their counts.
//!
//! Reading the input, which is normally a GenomicsDB, is not ported. The site rules, the germline
//! test and the beta fit are.
//!
//! The reference is
//! `org.broadinstitute.hellbender.tools.walkers.mutect.CreateSomaticPanelOfNormals` (GATK
//! 4.7.0.0; `CreateSomaticPanelOfNormals.java` is unchanged since 4.6.2.0).
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing" (`Some(x)` or `None`), here an AD array that may be
//!   absent; a `match` on it must handle both cases.
//! - `&[T]` is a borrowed, read-only view of a list; `Vec<&Genotype>` is a list of references to
//!   genotypes the site still owns, so nothing is copied.
//! - `.iter().filter(..).collect()` is Java's `stream().filter(..).collect(Collectors.toList())`;
//!   `filter_map` filters and maps in one step, keeping only the `Some` results.
//! - `|x| ..` is a closure, Java's lambda. `<F: Fn(&Site) -> f64>` makes a function generic over
//!   any callable from a site to a `double`.
//! - `.expect("..")` unwraps an `Option` or `Result` and aborts with that text if there is none;
//!   it marks a case the reference does not check because it cannot happen.
//! - `i32` is Java's `int`, `usize` an unsigned index type, `f64` a `double`; `x as f64` is the
//!   cast `(double) x`.

use gatk_engine::beta_binomial::BetaBinomialDistribution;
use gatk_engine::brent_optimizer::maximize;
use gatk_engine::math_utils::normalize_sum_to_one;
use gatk_engine::somatic_clustering_model::binomial_probability;

/// The tool's own hyperparameters. `ARTIFACT_PRIOR`: "the prior probability that any given site
/// has an artifact", as the reference's comment puts it.
pub const ARTIFACT_PRIOR: f64 = 0.001;
/// `ARTIFACT_ALPHA`: with `ARTIFACT_BETA`, the beta prior on an artefact's allele fraction.
pub const ARTIFACT_ALPHA: f64 = 1.0;
/// `ARTIFACT_BETA`.
pub const ARTIFACT_BETA: f64 = 7.0;
/// Below this, an allele frequency is treated as no frequency at all.
pub const NEGLIGIBLE_ALLELE_FREQUENCY: f64 = 1.0e-8;
/// `DEFAULT_MIN_SAMPLE_COUNT`, the default of `--min-sample-count`.
pub const DEFAULT_MIN_SAMPLE_COUNT: usize = 2;
/// `DEFAULT_MAX_GERMLINE_PROBABILITY`, the default of `--max-germline-probability`.
pub const DEFAULT_MAX_GERMLINE_PROBABILITY: f64 = 0.5;

/// One sample's call at one site. `allele_depths` is AD, reference first.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    /// The sample name, which no decision here reads.
    pub sample: String,
    /// `None` when the genotype carries no AD at all.
    pub allele_depths: Option<Vec<i32>>,
}

impl Genotype {
    /// `altCount`: the total depth less the reference depth, and zero without AD.
    ///
    /// ```java
    /// return g.hasAD() ? (int) MathUtils.sum(g.getAD()) - g.getAD()[0] : 0;
    /// ```
    ///
    /// So every alternate's depth counts, including a spanning deletion's.
    pub fn alt_count(&self) -> i32 {
        match &self.allele_depths {
            Some(depths) => depths.iter().sum::<i32>() - depths[0],
            None => 0,
        }
    }

    /// `(int) MathUtils.sum(g.getAD())`, as `hasArtifact` computes it; zero without AD.
    pub fn total_count(&self) -> i32 {
        match &self.allele_depths {
            Some(depths) => depths.iter().sum(),
            None => 0,
        }
    }
}

/// One site, as the walker sees it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    /// `vc.getContig()`.
    pub contig: String,
    /// `vc.getStart()`.
    pub position: i32,
    /// The reference allele's bases.
    pub reference: String,
    /// Empty for a site with no alternate at all.
    pub alternates: Vec<String>,
    /// `vc.getGenotypes()`, in header order.
    pub genotypes: Vec<Genotype>,
}

/// The spanning deletion, which is not an alternate worth a panel entry.
pub const SPAN_DEL: &str = "*";

impl Site {
    /// The two site shapes that are skipped before anything is counted: no alternate at all, and
    /// one whose ONLY alternate is the spanning deletion.
    ///
    /// ```java
    /// if (vc.getNAlleles() == 1) {
    ///     return;
    /// } else if (vc.getNAlleles() == 2 && vc.getAlternateAllele(0).basesMatch(Allele.SPAN_DEL)) {
    ///     return;
    /// }
    /// ```
    pub fn is_skipped(&self) -> bool {
        self.alternates.is_empty() || (self.alternates.len() == 1 && self.alternates[0] == SPAN_DEL)
    }

    /// `vc.getNAlleles() > 2`: more than one alternate.
    pub fn is_multiallelic(&self) -> bool {
        self.alternates.len() > 1
    }
}

/// `germlineProbability`.
///
/// A frequency below the negligible threshold, or above one, returns EXACTLY zero, which is why no
/// threshold can drop a site the germline resource does not mention: zero is below every one.
///
/// ```java
/// if (alleleFrequency < NEGLIGIBLE_ALLELE_FREQUENCY || alleleFrequency > 1) {
/// ```
///
/// The Rust test is the negation of "inside the closed range", which agrees with that line for
/// every number but not for `NaN`: the reference lets a `NaN` frequency through, and this returns
/// zero for it.
pub fn germline_probability(allele_frequency: f64, alt_count: i32, total_count: i32) -> f64 {
    if !(NEGLIGIBLE_ALLELE_FREQUENCY..=1.0).contains(&allele_frequency) {
        return 0.0;
    }
    // `alleleFrequency * (1 - alleleFrequency) * 2`, left to right as written.
    let het_prior = allele_frequency * (1.0 - allele_frequency) * 2.0;
    // `MathUtils.square(alleleFrequency)`, one multiplication.
    let hom_prior = allele_frequency * allele_frequency;
    let het_likelihood = binomial_probability(total_count, alt_count, 0.5);
    let hom_likelihood = binomial_probability(total_count, alt_count, 0.98);
    let artifact_likelihood =
        BetaBinomialDistribution::new(ARTIFACT_ALPHA, ARTIFACT_BETA, total_count)
            .expect("a valid shape")
            .probability(alt_count)
            .expect("a valid count");
    // Two entries, not three, whatever the comment above them says: the germline half is already
    // summed and the artefact half is the other. The reference's comment:
    // `the indices of this array are 0 -- germline het, 1 -- germline hom alt, 2 -- not germline`.
    let relative = [
        het_prior * het_likelihood + hom_prior * hom_likelihood,
        ARTIFACT_PRIOR * artifact_likelihood,
    ];
    // `MathUtils.sum(..) < 0 ? 0 : MathUtils.normalizeSumToOne(..)[0]`: a plain left-to-right
    // sum; `::<f64>` names the type `sum` accumulates into.
    if relative.iter().sum::<f64>() < 0.0 {
        return 0.0;
    }
    normalize_sum_to_one(&relative).expect("a positive sum")[0]
}

/// `hasArtifact`: a genotype with no alternate read never counts, whatever the frequency.
///
/// Strictly less than the threshold, as `germlineProbability(..) < maxGermlineProbability`.
pub fn has_artifact(
    genotype: &Genotype,
    population_allele_frequency: f64,
    max_germline_probability: f64,
) -> bool {
    let alt_count = genotype.alt_count();
    if alt_count == 0 {
        return false;
    }
    germline_probability(
        population_allele_frequency,
        alt_count,
        genotype.total_count(),
    ) < max_germline_probability
}

/// The genotypes a site contributes.
///
/// A MULTIALLELIC site skips the germline test entirely and counts every genotype, which the source
/// marks as a TODO rather than as a rule:
///
/// ```java
/// final List<Genotype> variantGenotypes = vc.getNAlleles() > 2 ? vc.getGenotypes() :
///         vc.getGenotypes().stream()
///         .filter(g -> hasArtifact(g, germlineAF)).collect(Collectors.toList());
/// ```
///
/// (The first line is wrapped here after the colon.)
pub fn variant_genotypes(
    site: &Site,
    germline_allele_frequency: f64,
    max_germline_probability: f64,
) -> Vec<&Genotype> {
    if site.is_multiallelic() {
        return site.genotypes.iter().collect();
    }
    site.genotypes
        .iter()
        .filter(|genotype| {
            has_artifact(
                genotype,
                germline_allele_frequency,
                max_germline_probability,
            )
        })
        .collect()
}

/// A fitted beta.
#[derive(Debug, Clone, Copy, PartialEq)]
///
/// `BetaDistributionShape`, written to the `BETA` INFO field as `[alpha, beta]`.
pub struct BetaShape {
    /// `getAlpha()`.
    pub alpha: f64,
    /// `getBeta()`.
    pub beta: f64,
}

/// `fitBeta`.
///
/// The base shape is the empirical mean, scaled so that alpha and beta keep their ratio, and the
/// scale is found by a Brent search over the summed beta-binomial log likelihood. The search runs
/// with the tool's own tolerances, so its answer carries them: see the `brent-optimizer` suite.
///
/// `alt_and_ref_counts` is the reference's `List<int[]>` of `{altCount, refCount}`, here a list
/// of pairs. The integer sums are `mapToInt(..).sum()`, exact. The log-likelihood is summed with
/// a plain `.sum()`, where the reference's is `DoubleStream.sum`
/// (`altAndRefCounts.stream().mapToDouble(..).sum()`), which is compensated: see
/// [`gatk_engine::allele_fraction_cluster::double_stream_sum`]. The two are different sums.
pub fn fit_beta(alt_and_ref_counts: &[(i32, i32)]) -> BetaShape {
    let total_alt: i32 = alt_and_ref_counts.iter().map(|(alt, _)| *alt).sum();
    let total_ref: i32 = alt_and_ref_counts
        .iter()
        .map(|(_, reference)| *reference)
        .sum();
    let min = total_alt.min(total_ref);
    // `(totalAltCount + 1.0) / (min + 1)`: "keeping the ratio of alpha and beta equal to the
    // ratio of baseAlpha and baseBeta gives the empirical mean", in the reference's words.
    let base_alpha = (total_alt as f64 + 1.0) / (min as f64 + 1.0);
    let base_beta = (total_ref as f64 + 1.0) / (min as f64 + 1.0);

    let log_likelihood = |scale: f64| -> f64 {
        let alpha = base_alpha * scale;
        let beta = base_beta * scale;
        alt_and_ref_counts
            .iter()
            .map(|(alt, reference)| {
                BetaBinomialDistribution::new(alpha, beta, alt + reference)
                    .expect("a valid shape")
                    .log_probability(*alt)
                    .expect("a valid count")
            })
            .sum()
    };

    // `OptimizationUtils.max(logLikelihood, 0.01, 100, 1, 0.01, 0.1, 100).getPoint()`: search
    // [0.01, 100] from a guess of 1, relative tolerance 0.01, absolute 0.1, 100 evaluations.
    let scale = maximize(log_likelihood, 0.01, 100.0, 1.0, 0.01, 0.1, 100)
        .expect("a search inside its budget")
        .point;
    BetaShape {
        alpha: base_alpha * scale,
        beta: base_beta * scale,
    }
}

/// One panel entry.
#[derive(Debug, Clone, PartialEq)]
pub struct PanelEntry {
    /// The site's contig.
    pub contig: String,
    /// The site's start.
    pub position: i32,
    /// The site's reference allele.
    pub reference: String,
    /// The site's alternates, all of them, the spanning deletion included.
    pub alternates: Vec<String>,
    /// The survivors over ALL samples in the header, not over the survivors.
    pub fraction: f64,
    /// The `BETA` INFO field.
    pub beta: BetaShape,
}

/// `apply` over a whole input.
///
/// `germline_allele_frequency` answers the resource for one site, summed over its alternates, and
/// gives zero where the resource says nothing. `sample_count` is
/// `getHeaderForVariants().getNGenotypeSamples()`, the `FRACTION` denominator.
pub fn build_panel<F: Fn(&Site) -> f64>(
    sites: &[Site],
    sample_count: usize,
    germline_allele_frequency: F,
    min_sample_count: usize,
    max_germline_probability: f64,
) -> Vec<PanelEntry> {
    let mut out = Vec::new();
    for site in sites {
        if site.is_skipped() {
            continue;
        }
        let frequency = germline_allele_frequency(site);
        let survivors = variant_genotypes(site, frequency, max_germline_probability);
        if survivors.len() < min_sample_count {
            continue;
        }
        // `(double) variantGenotypes.size() / numSamples`.
        let fraction = survivors.len() as f64 / sample_count as f64;
        // A genotype without AD counted for the site and contributes nothing to the fit:
        // `.filter(Genotype::hasAD).map(g -> new int[] {altCount(g), g.getAD()[0]})`.
        // `as_ref()` borrows the AD inside the `Option` rather than moving it out.
        let counts: Vec<(i32, i32)> = survivors
            .iter()
            .filter_map(|genotype| {
                genotype
                    .allele_depths
                    .as_ref()
                    .map(|depths| (genotype.alt_count(), depths[0]))
            })
            .collect();
        out.push(PanelEntry {
            contig: site.contig.clone(),
            position: site.position,
            reference: site.reference.clone(),
            alternates: site.alternates.clone(),
            fraction,
            beta: fit_beta(&counts),
        });
    }
    out
}
