//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine`
//! (`assignGenotypeLikelihoods`, `calculateGLsForThisEvent`) (GATK 4.6.2.0): the calls
//! HaplotypeCaller makes from a region's haplotypes and read likelihoods.
//!
//! At each event start inside the active window, the events at the start become one merged context
//! and an allele mapper ([`gatk_engine::allele_mapping`]); the read-by-haplotype matrix is
//! marginalized onto the merged alleles and keeps the reads overlapping the context widened by the
//! informative-read margin (2 by default, within the contig); each sample's genotype likelihoods
//! come from that matrix ([`gatk_engine::genotype_likelihood_calculator`]) as PLs; the genotyping
//! engine then decides QUAL, the alternates kept and the genotypes; and a call that lost alleles is
//! reverse-trimmed.
//!
//! What a reader would not guess: the genotyping alleles are the **matrix's** when it has as many
//! alleles as the merged context, and the context's otherwise, so the PL order follows the allele
//! mapper whenever the two agree in size.
//!
//! Not ported here: the annotations and physical phasing (their own bricks), reference-confidence
//! mode, contamination downsampling, DRAGstr priors and the BQD/FRD genotyping models, and the
//! reduction of sites with more alleles than the genotype count allows (refused).

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_mapping::{
    create_allele_mapper, make_merged_variant_context, replace_span_dels,
    variants_from_active_haplotypes,
};
use gatk_engine::event_map::{build_event_maps_for_haplotypes, event_start_positions};
use gatk_engine::genotype_index::genotype_count;
use gatk_engine::genotype_likelihood_calculator::{compute_log10_genotype_likelihoods, gls_to_pls};
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::read_utils;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, VariantContext};

use crate::genotyping_engine::{EngineError, GenotypingEngine};
use crate::variant_trim::reverse_trim_alleles;

/// The `HaplotypeCallerArgumentCollection` fields the genotyping engine reads.
#[derive(Debug, Clone)]
pub struct HcGenotypingArguments {
    pub sample_ploidy: usize,
    /// `--informative-read-overlap-margin`.
    pub informative_read_overlap_margin: i32,
    /// `--disable-spanning-event-genotyping`.
    pub disable_spanning_event_genotyping: bool,
    /// `--max-genotype-count`.
    pub max_genotype_count: usize,
    /// The maximum MNP distance the event maps are built with.
    pub max_mnp_distance: i32,
}

impl Default for HcGenotypingArguments {
    fn default() -> Self {
        HcGenotypingArguments {
            sample_ploidy: 2,
            informative_read_overlap_margin: 2,
            disable_spanning_event_genotyping: false,
            max_genotype_count: 1024,
            max_mnp_distance: 0,
        }
    }
}

/// `assignGenotypeLikelihoods`, without annotations or phasing: the calls, in event-start order.
#[allow(clippy::too_many_arguments)]
pub fn assign_genotype_likelihoods(
    engine: &mut GenotypingEngine,
    arguments: &HcGenotypingArguments,
    haplotypes: &mut [Haplotype],
    read_likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>,
    samples: &[String],
    reference: &[u8],
    ref_loc: &SimpleInterval,
    active_region_window: &SimpleInterval,
    contig_length: i32,
) -> Result<Vec<VariantContext>, EngineError> {
    let refused = |message: String| EngineError::Runtime {
        class: "IllegalArgumentException".to_string(),
        message,
    };
    build_event_maps_for_haplotypes(haplotypes, reference, ref_loc, arguments.max_mnp_distance)
        .map_err(|e| refused(e.message()))?;
    let starts = event_start_positions(haplotypes).map_err(|e| refused(e.message()))?;
    let ploidy = arguments.sample_ploidy;
    let spanning = !arguments.disable_spanning_event_genotyping;
    let mut calls = Vec::new();
    for loc in starts {
        if loc < active_region_window.start || loc > active_region_window.end {
            continue;
        }
        let events = variants_from_active_haplotypes(loc, haplotypes, spanning)
            .map_err(|e| refused(format!("{e:?}")))?;
        let ref_base = Allele::create(&reference[(loc - ref_loc.start) as usize..][..1], true)
            .map_err(|e| refused(e.to_string()))?;
        let replaced = replace_span_dels(events, &ref_base, loc);
        let Some(merged) =
            make_merged_variant_context(&replaced).map_err(|e| refused(format!("{e:?}")))?
        else {
            continue;
        };
        let merged_allele_count = merged.alleles.len();
        let mapper = create_allele_mapper(&merged, loc, haplotypes, spanning)
            .map_err(|e| refused(format!("{e:?}")))?;
        // `removeAltAllelesIfTooManyGenotypes`.
        if mapper.len() > max_acceptable_allele_count(ploidy, arguments.max_genotype_count) {
            return Err(EngineError::Limitation(
                "a site with more alleles than --max-genotype-count allows".to_string(),
            ));
        }
        let new_to_old: Vec<(Allele, Vec<Haplotype>)> = mapper
            .iter()
            .map(|(allele, indices)| {
                (
                    allele.clone(),
                    indices.iter().map(|&h| haplotypes[h].clone()).collect(),
                )
            })
            .collect();
        let mut marginal = read_likelihoods
            .marginalize(&new_to_old)
            .map_err(|e| refused(format!("{e:?}")))?;
        let overlap = SimpleInterval::new(&merged.contig, merged.start, merged.end)
            .and_then(|span| {
                span.expand_within_contig(arguments.informative_read_overlap_margin, contig_length)
            })
            .ok_or_else(|| refused("no overlap interval".to_string()))?;
        marginal.retain_evidence(|read| {
            overlap.overlaps(
                &overlap.contig,
                read_utils::start(read),
                read_utils::end(read),
            )
        });

        let genotypes = genotypes_for_event(&marginal, &merged.alleles, samples, ploidy)?;
        let mut vc = VariantContext::new(
            &merged.contig,
            i64::from(merged.start),
            merged.alleles.clone(),
        );
        vc.stop = i64::from(merged.end);
        vc.genotypes = GenotypesContext::new(genotypes);
        if let Some(call) = engine.calculate_genotypes(&vc)? {
            let call = if call.alleles.len() == merged_allele_count {
                call
            } else {
                reverse_trim_alleles(&call).map_err(refused)?
            };
            calls.push(call);
        }
    }
    Ok(calls)
}

/// `calculateGLsForThisEvent`: each sample, no-called, with the PLs of its genotype likelihoods
/// over the genotyping alleles.
fn genotypes_for_event(
    marginal: &AlleleLikelihoods<BamRecord, Allele>,
    merged_alleles: &[Allele],
    samples: &[String],
    ploidy: usize,
) -> Result<Vec<Genotype>, EngineError> {
    let genotyping_alleles: Vec<Allele> = if marginal.number_of_alleles() == merged_alleles.len() {
        (0..marginal.number_of_alleles())
            .map(|a| marginal.get_allele(a).expect("an allele").clone())
            .collect()
    } else {
        merged_alleles.to_vec()
    };
    let permutation: Vec<usize> = genotyping_alleles
        .iter()
        .map(|allele| {
            marginal
                .index_of_allele(allele)
                .ok_or_else(|| EngineError::Runtime {
                    class: "IllegalArgumentException".to_string(),
                    message: format!(
                        "allele {} missing from the likelihoods",
                        allele.display_string()
                    ),
                })
        })
        .collect::<Result<_, _>>()?;
    let mut genotypes = Vec::with_capacity(samples.len());
    for (s, sample) in samples.iter().enumerate() {
        let evidence = marginal.sample_evidence_count(s);
        let matrix: Vec<Vec<f64>> = permutation
            .iter()
            .map(|&a| (0..evidence).map(|r| marginal.value(s, a, r)).collect())
            .collect();
        let likelihoods = compute_log10_genotype_likelihoods(ploidy, &matrix);
        let mut genotype = Genotype::new(sample, vec![Allele::no_call(); ploidy]);
        genotype.pl = Some(gls_to_pls(&likelihoods));
        genotypes.push(genotype);
    }
    Ok(genotypes)
}

/// `GenotypeIndexCalculator.computeMaxAcceptableAlleleCount(ploidy, maxGenotypeCount)`: the most
/// alleles whose genotype count stays within the limit.
fn max_acceptable_allele_count(ploidy: usize, max_genotype_count: usize) -> usize {
    let mut alleles = 1;
    while genotype_count(ploidy, alleles + 1).is_ok_and(|count| count <= max_genotype_count) {
        alleles += 1;
    }
    alleles
}
