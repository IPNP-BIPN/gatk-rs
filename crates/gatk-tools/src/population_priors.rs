//! `--population-callset` in HaplotypeCaller's reference-confidence mode, ported from
//! `ReferenceConfidenceModel.calculateRefConfidence(..., applyPriors, VCpriors)` and
//! `PosteriorProbabilitiesUtils.calculatePosteriorProbs` (GATK 4.6.2.0).
//!
//! With a supporting callset every record the model emits goes through the posteriors, whether a
//! prior matches it or not: a reference block, and a call that starts at the position. The priors
//! are the callset's records starting there; without one, the counts come from
//! `--num-reference-samples-if-no-call` alone, and a flat prior otherwise. The model's options are
//! fixed: the SNP and indel heterozygosities as pseudocounts, the input's own counts ignored, MLEAC
//! preferred. A call gains `PG` and recomputed `AC`, `AF` and `AN`; a hom-ref block gains neither.
//! Every genotype gains `PP`, and its call and `GQ` come from the posteriors.

use htsjdk_vcf::variant::{format_vcf_double, Value, VariantContext};

use crate::calculate_genotype_posteriors::{self as posteriors, record_from_vc, Options, Record};

/// The supporting callset and the reference samples assumed where it is silent.
#[derive(Debug, Clone)]
pub struct PopulationPriors {
    /// Each record with its contig.
    pub records: Vec<(String, Record)>,
    /// `--num-reference-samples-if-no-call`.
    pub num_ref_if_missing: i32,
}

/// `ReferenceConfidenceModel`'s `PosteriorProbabilitiesOptions`.
const OPTIONS: Options = Options {
    snp_prior_dirichlet: 0.001,
    indel_prior_dirichlet: 1.25e-4,
    use_input_samples_allele_counts: false,
    use_mleac: true,
    ignore_input_samples_for_missing_resources: true,
    use_flat_priors_for_indels: false,
};

impl PopulationPriors {
    /// `calculatePosteriorProbs(vc, getMatchingPriors(position, ...), numRefSamplesForPrior, options)`,
    /// which refuses a record whose largest ploidy (two without genotypes) is not two: a haploid
    /// or triploid `--ploidy-regions` interval meets that refusal at its first record.
    pub fn apply(&self, vc: &VariantContext) -> Result<VariantContext, String> {
        let max_ploidy = vc
            .genotypes
            .iter()
            .map(|genotype| genotype.alleles.len())
            .max()
            .filter(|&ploidy| ploidy > 0)
            .unwrap_or(2);
        if max_ploidy != 2 {
            return Err("Genotype posteriors not yet implemented for ploidy != 2".to_string());
        }
        let matching: Vec<Record> = self
            .records
            .iter()
            .filter(|(contig, record)| *contig == vc.contig && i64::from(record.start) == vc.start)
            .map(|(_, record)| record.clone())
            .collect();
        let answer = posteriors::calculate_posterior_probs(
            &record_from_vc(vc),
            &matching,
            self.num_ref_if_missing,
            &OPTIONS,
        );
        let mut out = vc.clone();
        for called in &answer.genotypes {
            let Some(genotype) = out
                .genotypes
                .iter_mut()
                .find(|genotype| genotype.sample_name == called.sample)
            else {
                continue;
            };
            let Some(values) = &called.posteriors else {
                continue;
            };
            genotype.alleles = called
                .alleles
                .iter()
                .filter_map(|index| vc.alleles.get(*index).cloned())
                .collect();
            if let Some(gq) = called.gq {
                genotype.gq = Some(gq);
            }
            genotype.extended.retain(|(key, _)| key != "PP");
            genotype
                .extended
                .push(("PP".to_string(), Value::Str(join(values))));
        }
        if let Some(prior) = &answer.prior {
            set(&mut out, "PG", join(prior));
            chromosome_counts(&mut out);
        }
        Ok(out)
    }
}

fn join(values: &[i32]) -> String {
    values
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn set(vc: &mut VariantContext, key: &str, value: String) {
    vc.attributes.retain(|(name, _)| name != key);
    vc.attributes.push((key.to_string(), Value::Str(value)));
}

/// `VariantContextUtils.calculateChromosomeCounts(builder, true)`: AN over the called alleles of
/// the unfiltered genotypes, and AC and AF for each alternate.
fn chromosome_counts(vc: &mut VariantContext) {
    if vc.genotypes.is_empty() {
        return;
    }
    let called = || {
        vc.genotypes
            .iter()
            .filter(|g| g.filters.as_deref().is_none_or(|f| f == "PASS" || f == "."))
    };
    let allele_number: i32 = called()
        .map(|g| g.alleles.iter().filter(|a| !a.is_no_call()).count() as i32)
        .sum();
    let mut counts = Vec::new();
    let mut frequencies = Vec::new();
    for allele in &vc.alleles[1..] {
        let count: i32 = called()
            .map(|g| g.alleles.iter().filter(|a| *a == allele).count() as i32)
            .sum();
        counts.push(count.to_string());
        let frequency = if allele_number == 0 {
            0.0
        } else {
            f64::from(count) / f64::from(allele_number)
        };
        frequencies.push(format_vcf_double(frequency));
    }
    set(vc, "AN", allele_number.to_string());
    if vc.alleles.len() > 1 {
        set(vc, "AC", counts.join(","));
        set(vc, "AF", frequencies.join(","));
    }
}
