//! `CalculateGenotypePosteriors` and `PosteriorProbabilitiesUtils`, ported from GATK 4.6.2.0.
//!
//! Genotype likelihoods turned into posteriors under a Dirichlet-multinomial prior built from
//! allele counts. Where those counts come from is the whole of the tool. Reading the VCF is not
//! ported, and neither is the family-prior half: `FamilyLikelihoods` is a second algorithm, and a
//! run with no pedigree skips it, which is the configuration measured here.
//!
//! # A site nothing supplies counts for gets a flat prior
//!
//! ```java
//! final boolean useFlatPriors = (!vc1.isSNP() && opts.useFlatPriorsForIndels)
//!         || (resources.isEmpty() && !useDiscoveredACForMissing && numRefSamplesFromMissingResources == 0);
//! ```
//!
//! A flat prior is 1.0 for every genotype, so `PG` comes out all zeros and `PP` equals `PL`.
//!
//! # And the input's own samples only count when there are ten of them
//!
//! ```java
//! final boolean useDiscoveredACForMissing = !opts.ignoreInputSamplesForMissingResources
//!         && (vc1.getNSamples() >= minSamplesToUseInputs || numRefSamplesFromMissingResources != 0);
//! ```
//!
//! `minSamplesToUseInputs` is 10. Below that, a site the panel does not carry falls back to flat
//! priors unless `--num-reference-samples-if-no-call` was given.
//!
//! # MLEAC is preferred, and the reference count is AN minus the alternates
//!
//! ```java
//! if ( context.hasAttribute(MLEAC) && !useAC ) { ac = getAlleleCounts(MLEAC, context); }
//! else if ( context.hasAttribute(AC) ) { ac = getAlleleCounts(AC, context); }
//! ```
//!
//! A panel with `AN=2000;MLEAC=20;AC=200` contributes 1980 reference chromosomes and 20 alternate
//! ones, and `--default-to-allele-count` makes it 1800 and 200. A panel carrying no MLEAC at all
//! falls through to AC and gives the same answer as the flag. The reference count is floored at
//! zero, "because occasionally an MLEAC value will sneak in that's greater than the AN".
//!
//! # The reference samples enter as chromosomes, and only where the panel is silent
//!
//! `--num-reference-samples-if-no-call` is doubled before it enters the prior, and it reaches
//! `calculatePosteriorProbs` as zero at every site the panel does carry.
//!
//! # The SNP pseudocount is chosen by allele LENGTH
//!
//! An allele the same length as the reference takes the SNP pseudocount, a symbolic one takes the
//! larger of the two, and everything else takes the indel one. So a single site can mix them.
//!
//! # And the genotypes are recalled from the posteriors
//!
//! `makeGenotypeCall` with `USE_PLS_TO_ASSIGN` picks the most likely posterior, so a `GT` can
//! change and the `GQ` is recomputed from the `PP` rather than kept from the `PL`.

use jmath::gamma::log_gamma;
use std::collections::BTreeMap;

/// `PosteriorProbabilitiesUtils.minSamplesToUseInputs`.
pub const MIN_SAMPLES_TO_USE_INPUTS: usize = 10;
/// `HomoSapiensConstants.SNP_HETEROZYGOSITY`, which is both defaults.
pub const DEFAULT_PRIOR: f64 = 1e-3;
/// `HomoSapiensConstants.DEFAULT_PLOIDY`.
pub const DEFAULT_PLOIDY: usize = 2;

/// One allele, which is the reference one or not.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Allele {
    pub bases: String,
    pub is_ref: bool,
}

impl Allele {
    pub fn is_symbolic(&self) -> bool {
        self.bases.starts_with('<') || self.bases.starts_with('*')
    }

    pub fn is_non_ref(&self) -> bool {
        self.bases == "<NON_REF>"
    }
}

/// One sample's genotype, reduced to what the posteriors read and write.
#[derive(Debug, Clone, PartialEq)]
pub struct Genotype {
    pub sample: String,
    /// Allele indices into the record's allele list.
    pub alleles: Vec<usize>,
    pub depth: Option<i32>,
    pub likelihoods: Option<Vec<i32>>,
    /// `PP`, which a second pass would read in preference to `PL`.
    pub posteriors: Option<Vec<i32>>,
}

/// One record, reduced the same way.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub id: String,
    pub start: i32,
    pub alleles: Vec<Allele>,
    /// The INFO attributes this reads: `AC`, `AN` and `MLEAC`, as they appear in the file.
    pub attributes: BTreeMap<String, String>,
    pub genotypes: Vec<Genotype>,
}

impl Record {
    pub fn reference(&self) -> &Allele {
        self.alleles.first().expect("a reference allele")
    }

    pub fn alternates(&self) -> Vec<&Allele> {
        self.alleles.iter().skip(1).collect()
    }

    /// `VariantContext.isSNP`, which does NOT require the site to be biallelic: a reference of one
    /// base and every alternate of one base is a SNP however many alternates there are.
    ///
    /// Getting this wrong makes `--use-flat-priors-for-indels` flatten a triallelic SNP, which is
    /// what the first version of this port did and what the golden caught.
    pub fn is_snp(&self) -> bool {
        self.reference().bases.len() == 1
            && self
                .alternates()
                .iter()
                .all(|allele| allele.bases.len() == 1 && !allele.is_symbolic())
    }

    fn attribute_ints(&self, key: &str) -> Option<Vec<i32>> {
        self.attributes.get(key).map(|value| {
            value
                .split(',')
                .map(|part| part.parse().expect("an integer attribute"))
                .collect()
        })
    }

    /// `getCalledChrCount(allele)` over the genotypes.
    fn called_chr_count(&self, allele: usize) -> i32 {
        self.genotypes
            .iter()
            .flat_map(|genotype| genotype.alleles.iter())
            .filter(|index| **index == allele)
            .count() as i32
    }

    fn called_chr_count_total(&self) -> i32 {
        self.genotypes
            .iter()
            .map(|genotype| genotype.alleles.len() as i32)
            .sum()
    }
}

/// `PosteriorProbabilitiesOptions`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    pub snp_prior_dirichlet: f64,
    pub indel_prior_dirichlet: f64,
    pub use_input_samples_allele_counts: bool,
    pub use_mleac: bool,
    pub ignore_input_samples_for_missing_resources: bool,
    pub use_flat_priors_for_indels: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            snp_prior_dirichlet: DEFAULT_PRIOR,
            indel_prior_dirichlet: DEFAULT_PRIOR,
            use_input_samples_allele_counts: true,
            use_mleac: true,
            ignore_input_samples_for_missing_resources: false,
            use_flat_priors_for_indels: false,
        }
    }
}

/// `MathUtils.dirichletMultinomial`, over `Gamma.logGamma`.
pub fn dirichlet_multinomial(params: &[f64], counts: &[i32]) -> f64 {
    let dirichlet_sum: f64 = params.iter().sum();
    let count_sum: f64 = counts.iter().map(|count| f64::from(*count)).sum();
    let mut value = log_gamma(count_sum + 1.0) + log_gamma(dirichlet_sum)
        - log_gamma(dirichlet_sum + count_sum);
    for (param, count) in params.iter().zip(counts) {
        let count = f64::from(*count);
        value += log_gamma(count + param) - log_gamma(*param) - log_gamma(count + 1.0);
    }
    // `logToLog10`.
    value / std::f64::consts::LN_10
}

/// `getDirichletPrior`, over the diploid genotype order AA AB BB AC BC CC.
pub fn get_dirichlet_prior(known_counts: &[f64], ploidy: usize, flat: bool) -> Vec<f64> {
    assert_eq!(
        ploidy, 2,
        "genotype priors are only implemented for ploidy 2"
    );
    let mut priors = Vec::with_capacity(known_counts.len() * (known_counts.len() + 1) / 2);
    for allele2 in 0..known_counts.len() {
        for allele1 in 0..=allele2 {
            if flat {
                priors.push(1.0);
            } else {
                let mut counts = vec![0; known_counts.len()];
                counts[allele1] += 1;
                counts[allele2] += 1;
                priors.push(dirichlet_multinomial(known_counts, &counts));
            }
        }
    }
    priors
}

/// `GenotypeLikelihoods.GLsToPLs`: the maximum becomes zero and the rest are rounded distances.
pub fn gls_to_pls(log10: &[f64]) -> Vec<i32> {
    let adjust = log10.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    log10
        .iter()
        .map(|value| ((-10.0 * (value - adjust)) + 0.5).floor() as i32)
        .collect()
}

/// `GenotypeLikelihoods.getGQLog10FromLikelihoods`, for the chosen genotype.
pub fn gq_log10_from_likelihoods(chosen: usize, likelihoods: &[f64]) -> f64 {
    let mut qual = f64::NEG_INFINITY;
    for (index, value) in likelihoods.iter().enumerate() {
        if index != chosen && *value >= qual {
            qual = *value;
        }
    }
    let qual = likelihoods[chosen] - qual;
    if qual < 0.0 {
        let maximum = likelihoods
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let total: f64 = likelihoods
            .iter()
            .map(|value| 10f64.powf(value - maximum))
            .sum();
        let chosen_probability = 10f64.powf(likelihoods[chosen] - maximum) / total;
        (1.0 - chosen_probability).log10()
    } else {
        -qual
    }
}

/// `hasRealLikelihoods`: a genotype at depth zero whose PLs are all zero carries no data.
pub fn has_real_likelihoods(genotype: &Genotype) -> bool {
    let Some(likelihoods) = &genotype.likelihoods else {
        return false;
    };
    if genotype.depth == Some(0) {
        return likelihoods.iter().max().copied().unwrap_or(0) > 0;
    }
    true
}

/// `parsePosteriorsIntoProbSpace`: the `PP` if there is one, otherwise the `PL`.
pub fn parse_posteriors_into_prob_space(genotype: &Genotype) -> Option<Vec<f64>> {
    if let Some(posteriors) = &genotype.posteriors {
        return Some(
            posteriors
                .iter()
                .map(|value| f64::from(*value) / -10.0)
                .collect(),
        );
    }
    if !has_real_likelihoods(genotype) {
        return None;
    }
    genotype.likelihoods.as_ref().map(|likelihoods| {
        likelihoods
            .iter()
            .map(|value| f64::from(*value) / -10.0)
            .collect()
    })
}

/// `addAlleleCounts`: MLEAC first unless `use_ac`, then AC, then the genotypes themselves.
///
/// The alleles are keyed by their remapped selves; this fixture never needs a remapping, since
/// every resource shares the input's reference allele, and the padding path is not ported.
pub fn add_allele_counts(counts: &mut BTreeMap<Allele, i32>, record: &Record, use_ac: bool) {
    let allele_counts: Vec<i32> = if record.attributes.contains_key("MLEAC") && !use_ac {
        record.attribute_ints("MLEAC").expect("an MLEAC")
    } else if record.attributes.contains_key("AC") {
        record.attribute_ints("AC").expect("an AC")
    } else {
        (1..record.alleles.len())
            .map(|index| record.called_chr_count(index))
            .collect()
    };
    let alternate_sum: i32 = allele_counts.iter().sum();

    for (index, allele) in record.alleles.iter().enumerate() {
        let count = if allele.is_ref {
            // The reference count is never written down, so it is AN minus the alternates, and
            // floored at zero because an MLEAC can exceed the AN.
            let total = match record.attributes.get("AN") {
                Some(value) => value.parse::<i32>().expect("an AN"),
                None => record.called_chr_count_total(),
            };
            (total - alternate_sum).max(0)
        } else {
            allele_counts[index - 1]
        };
        *counts.entry(allele.clone()).or_insert(0) += count;
    }
}

/// What one record comes out as.
#[derive(Debug, Clone, PartialEq)]
pub struct Posteriors {
    /// The `PG` INFO field, absent for a hom-ref block.
    pub prior: Option<Vec<i32>>,
    /// Per sample: the recalled allele indices, the `GQ` and the `PP`.
    pub genotypes: Vec<CalledGenotype>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CalledGenotype {
    pub sample: String,
    pub alleles: Vec<usize>,
    pub gq: Option<i32>,
    pub posteriors: Option<Vec<i32>>,
}

/// `calculatePosteriorProbs`, over one record and the resources that start where it does.
pub fn calculate_posterior_probs(
    record: &Record,
    resources: &[Record],
    num_ref_samples_from_missing_resources: i32,
    options: &Options,
) -> Posteriors {
    let use_discovered_ac_for_missing = !options.ignore_input_samples_for_missing_resources
        && (record.genotypes.len() >= MIN_SAMPLES_TO_USE_INPUTS
            || num_ref_samples_from_missing_resources != 0);

    let mut total_allele_counts: BTreeMap<Allele, i32> = BTreeMap::new();
    let reference_allele_count_for_missing = if resources.is_empty() {
        DEFAULT_PLOIDY as i32 * num_ref_samples_from_missing_resources
    } else {
        0
    };

    for resource in resources {
        if resource.start == record.start {
            add_allele_counts(&mut total_allele_counts, resource, !options.use_mleac);
        }
    }
    if (options.use_input_samples_allele_counts && !resources.is_empty())
        || (resources.is_empty() && use_discovered_ac_for_missing)
    {
        add_allele_counts(&mut total_allele_counts, record, !options.use_mleac);
    }

    let common_ref = record.reference().clone();
    *total_allele_counts.entry(common_ref.clone()).or_insert(0) +=
        reference_allele_count_for_missing;

    // The pseudocount is chosen by LENGTH, not by the site's type.
    let mut allele_counts: Vec<f64> = record
        .alleles
        .iter()
        .map(|allele| {
            let observed = f64::from(*total_allele_counts.get(allele).unwrap_or(&0));
            if allele.bases.len() == common_ref.bases.len() && !allele.is_symbolic() {
                options.snp_prior_dirichlet + observed
            } else if allele.is_symbolic() {
                options
                    .snp_prior_dirichlet
                    .max(options.indel_prior_dirichlet)
                    + observed
            } else {
                options.indel_prior_dirichlet + observed
            }
        })
        .collect();

    // Every count belonging to an allele the input did not carry goes to <NON_REF>, which also
    // takes the larger of the two pseudocounts rather than one of them.
    if let Some(non_ref_index) = record.alleles.iter().position(Allele::is_non_ref) {
        let resource_only: i32 = total_allele_counts
            .iter()
            .filter(|(allele, _)| !record.alleles.contains(allele))
            .map(|(_, count)| *count)
            .sum();
        allele_counts[non_ref_index] = options
            .snp_prior_dirichlet
            .max(options.indel_prior_dirichlet)
            + f64::from(resource_only);
    }

    let use_flat_priors = (!record.is_snp() && options.use_flat_priors_for_indels)
        || (resources.is_empty()
            && !use_discovered_ac_for_missing
            && num_ref_samples_from_missing_resources == 0);

    let prior = get_dirichlet_prior(&allele_counts, DEFAULT_PLOIDY, use_flat_priors);
    let mut called = Vec::new();
    for genotype in &record.genotypes {
        let likelihoods = parse_posteriors_into_prob_space(genotype);
        match likelihoods {
            None => called.push(CalledGenotype {
                sample: genotype.sample.clone(),
                alleles: genotype.alleles.clone(),
                gq: None,
                posteriors: None,
            }),
            Some(likelihoods) => {
                assert_eq!(
                    likelihoods.len(),
                    prior.len(),
                    "likelihoods not of correct size"
                );
                let posterior: Vec<f64> = likelihoods
                    .iter()
                    .zip(&prior)
                    .map(|(likelihood, prior)| likelihood + prior)
                    .collect();
                let best = max_element_index(&posterior);
                let alleles = genotype_alleles(DEFAULT_PLOIDY, best);
                let gq = gq_log10_from_likelihoods(best, &posterior);
                called.push(CalledGenotype {
                    sample: genotype.sample.clone(),
                    // A recalled genotype whose alleles include <NON_REF> is written as hom-ref.
                    alleles: if alleles
                        .iter()
                        .any(|index| record.alleles[*index].is_non_ref())
                    {
                        vec![0; DEFAULT_PLOIDY]
                    } else {
                        alleles
                    },
                    gq: Some(gq_to_phred(gq)),
                    posteriors: Some(gls_to_pls(&posterior)),
                });
            }
        }
    }

    // A hom-ref block keeps its counts and gains no prior.
    let is_hom_ref_block =
        record.alternates().len() == 1 && record.alleles.iter().any(Allele::is_non_ref);
    Posteriors {
        prior: if is_hom_ref_block {
            None
        } else {
            Some(gls_to_pls(&prior))
        },
        genotypes: called,
    }
}

/// `MathUtils.maxElementIndex`, which keeps the FIRST maximum.
fn max_element_index(values: &[f64]) -> usize {
    let mut best = 0;
    for (index, value) in values.iter().enumerate() {
        if *value > values[best] {
            best = index;
        }
    }
    best
}

/// `GenotypesCache.get(ploidy, index).asAlleleList`, for the diploid order AA AB BB AC BC CC.
pub fn genotype_alleles(ploidy: usize, index: usize) -> Vec<usize> {
    assert_eq!(ploidy, 2, "only the diploid order is ported");
    let mut position = 0;
    for allele2 in 0.. {
        for allele1 in 0..=allele2 {
            if position == index {
                return vec![allele1, allele2];
            }
            position += 1;
        }
    }
    unreachable!()
}

/// `Genotype.getGQ`, which is the log10 error rounded and capped at 99.
fn gq_to_phred(log10_perror: f64) -> i32 {
    let value = ((-10.0 * log10_perror) + 0.5).floor() as i64;
    value.clamp(0, 99) as i32
}

/// `apply`: the resources that start where the record does, and the reference-sample count that
/// only applies when there are none.
pub fn apply(
    record: &Record,
    resources: &[Record],
    num_ref_if_missing: i32,
    options: &Options,
    skip_population_priors: bool,
) -> Option<Posteriors> {
    if skip_population_priors {
        return None;
    }
    let matching: Vec<Record> = resources
        .iter()
        .filter(|resource| resource.start == record.start)
        .cloned()
        .collect();
    let num_ref = if matching.is_empty() {
        num_ref_if_missing
    } else {
        0
    };
    Some(calculate_posterior_probs(
        record, &matching, num_ref, options,
    ))
}
