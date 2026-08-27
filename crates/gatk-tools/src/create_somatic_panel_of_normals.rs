//! `CreateSomaticPanelOfNormals`: which sites recur across normals often enough to be artefacts.
//!
//! A site is kept when enough samples carry an alternate that is more likely to be an artefact than
//! to be germline, and it is written with the fraction of samples that carried it and a beta shape
//! fitted to their counts.
//!
//! Reading the input, which is normally a GenomicsDB, is not ported. The site rules, the germline
//! test and the beta fit are.

use gatk_engine::beta_binomial::BetaBinomialDistribution;
use gatk_engine::brent_optimizer::maximize;
use gatk_engine::math_utils::normalize_sum_to_one;
use gatk_engine::somatic_clustering_model::binomial_probability;

/// The tool's own hyperparameters.
pub const ARTIFACT_PRIOR: f64 = 0.001;
pub const ARTIFACT_ALPHA: f64 = 1.0;
pub const ARTIFACT_BETA: f64 = 7.0;
/// Below this, an allele frequency is treated as no frequency at all.
pub const NEGLIGIBLE_ALLELE_FREQUENCY: f64 = 1.0e-8;
pub const DEFAULT_MIN_SAMPLE_COUNT: usize = 2;
pub const DEFAULT_MAX_GERMLINE_PROBABILITY: f64 = 0.5;

/// One sample's call at one site. `allele_depths` is AD, reference first.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    pub sample: String,
    /// `None` when the genotype carries no AD at all.
    pub allele_depths: Option<Vec<i32>>,
}

impl Genotype {
    /// `altCount`: the total depth less the reference depth, and zero without AD.
    pub fn alt_count(&self) -> i32 {
        match &self.allele_depths {
            Some(depths) => depths.iter().sum::<i32>() - depths[0],
            None => 0,
        }
    }

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
    pub contig: String,
    pub position: i32,
    pub reference: String,
    /// Empty for a site with no alternate at all.
    pub alternates: Vec<String>,
    pub genotypes: Vec<Genotype>,
}

/// The spanning deletion, which is not an alternate worth a panel entry.
pub const SPAN_DEL: &str = "*";

impl Site {
    /// The two site shapes that are skipped before anything is counted: no alternate at all, and
    /// one whose ONLY alternate is the spanning deletion.
    pub fn is_skipped(&self) -> bool {
        self.alternates.is_empty() || (self.alternates.len() == 1 && self.alternates[0] == SPAN_DEL)
    }

    pub fn is_multiallelic(&self) -> bool {
        self.alternates.len() > 1
    }
}

/// `germlineProbability`.
///
/// A frequency below the negligible threshold, or above one, returns EXACTLY zero, which is why no
/// threshold can drop a site the germline resource does not mention: zero is below every one.
pub fn germline_probability(allele_frequency: f64, alt_count: i32, total_count: i32) -> f64 {
    if !(NEGLIGIBLE_ALLELE_FREQUENCY..=1.0).contains(&allele_frequency) {
        return 0.0;
    }
    let het_prior = allele_frequency * (1.0 - allele_frequency) * 2.0;
    let hom_prior = allele_frequency * allele_frequency;
    let het_likelihood = binomial_probability(total_count, alt_count, 0.5);
    let hom_likelihood = binomial_probability(total_count, alt_count, 0.98);
    let artifact_likelihood =
        BetaBinomialDistribution::new(ARTIFACT_ALPHA, ARTIFACT_BETA, total_count)
            .expect("a valid shape")
            .probability(alt_count)
            .expect("a valid count");
    // Two entries, not three, whatever the comment above them says: the germline half is already
    // summed and the artefact half is the other.
    let relative = [
        het_prior * het_likelihood + hom_prior * hom_likelihood,
        ARTIFACT_PRIOR * artifact_likelihood,
    ];
    if relative.iter().sum::<f64>() < 0.0 {
        return 0.0;
    }
    normalize_sum_to_one(&relative).expect("a positive sum")[0]
}

/// `hasArtifact`: a genotype with no alternate read never counts, whatever the frequency.
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
/// marks as a TODO rather than as a rule.
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
pub struct BetaShape {
    pub alpha: f64,
    pub beta: f64,
}

/// `fitBeta`.
///
/// The base shape is the empirical mean, scaled so that alpha and beta keep their ratio, and the
/// scale is found by a Brent search over the summed beta-binomial log likelihood. The search runs
/// with the tool's own tolerances, so its answer carries them: see the `brent-optimizer` suite.
pub fn fit_beta(alt_and_ref_counts: &[(i32, i32)]) -> BetaShape {
    let total_alt: i32 = alt_and_ref_counts.iter().map(|(alt, _)| *alt).sum();
    let total_ref: i32 = alt_and_ref_counts
        .iter()
        .map(|(_, reference)| *reference)
        .sum();
    let min = total_alt.min(total_ref);
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
    pub contig: String,
    pub position: i32,
    pub reference: String,
    pub alternates: Vec<String>,
    /// The survivors over ALL samples in the header, not over the survivors.
    pub fraction: f64,
    pub beta: BetaShape,
}

/// `apply` over a whole input.
///
/// `germline_allele_frequency` answers the resource for one site, summed over its alternates, and
/// gives zero where the resource says nothing.
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
        let fraction = survivors.len() as f64 / sample_count as f64;
        // A genotype without AD counted for the site and contributes nothing to the fit.
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
