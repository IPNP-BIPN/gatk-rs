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
//! The annotations are the given annotator engine's, run before the trim as `makeAnnotatedCall`
//! runs them; with none, the call is left as genotyped.
//!
//! Physical phasing (`AssemblyBasedCallerUtils.phaseCalls`) runs over the trimmed calls when
//! asked for: calls on the same called haplotypes are phased `0|1` together, calls on
//! complementary ones `0|1` against `1|0`, and a call whose partner already sits in another group
//! empties the whole mapping.
//!
//! Not ported here: reference-confidence
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
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

use crate::genotyping_engine::{EngineError, GenotypingEngine};
use crate::variant_annotator_engine::{Engine as AnnotatorEngine, Site};
use crate::variant_trim::reverse_trim_alleles;
use gatk_engine::java_random::JavaRandom;

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
    /// `doPhysicalPhasing`, which the tool sets from `--do-not-run-physical-phasing`. Off here
    /// unless asked for, so that the genotyping can be measured on its own.
    pub do_physical_phasing: bool,
}

impl Default for HcGenotypingArguments {
    fn default() -> Self {
        HcGenotypingArguments {
            sample_ploidy: 2,
            informative_read_overlap_margin: 2,
            disable_spanning_event_genotyping: false,
            max_genotype_count: 1024,
            max_mnp_distance: 0,
            do_physical_phasing: false,
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
    assign_genotype_likelihoods_annotated(
        engine,
        arguments,
        haplotypes,
        read_likelihoods,
        samples,
        reference,
        ref_loc,
        active_region_window,
        contig_length,
        None,
    )
}

/// The annotation engine `makeAnnotatedCall` runs, and the tool's random generator, which
/// `QualByDepth` draws from.
pub struct CallAnnotator<'a> {
    pub engine: &'a AnnotatorEngine,
    pub random: &'a mut JavaRandom,
}

/// `assignGenotypeLikelihoods` with the annotations: each call is annotated from the
/// read-by-allele matrix it was genotyped with, then reverse-trimmed if it lost alleles.
#[allow(clippy::too_many_arguments)]
pub fn assign_genotype_likelihoods_annotated(
    engine: &mut GenotypingEngine,
    arguments: &HcGenotypingArguments,
    haplotypes: &mut [Haplotype],
    read_likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>,
    samples: &[String],
    reference: &[u8],
    ref_loc: &SimpleInterval,
    active_region_window: &SimpleInterval,
    contig_length: i32,
    annotator: Option<CallAnnotator<'_>>,
) -> Result<Vec<VariantContext>, EngineError> {
    assign_genotype_likelihoods_full(
        engine,
        arguments,
        haplotypes,
        read_likelihoods,
        samples,
        reference,
        ref_loc,
        active_region_window,
        contig_length,
        annotator,
        &[],
    )
}

/// `assignGenotypeLikelihoods` as `callRegion` runs it: with the annotations, and with each
/// sample's reads filtered before genotyping, which the annotations still see
/// (`addEvidence(overlappingFilteredReads, 0)`) wherever they overlap a call.
#[allow(clippy::too_many_arguments)]
pub fn assign_genotype_likelihoods_full(
    engine: &mut GenotypingEngine,
    arguments: &HcGenotypingArguments,
    haplotypes: &mut [Haplotype],
    read_likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>,
    samples: &[String],
    reference: &[u8],
    ref_loc: &SimpleInterval,
    active_region_window: &SimpleInterval,
    contig_length: i32,
    mut annotator: Option<CallAnnotator<'_>>,
    filtered_reads: &[Vec<BamRecord>],
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
    // `calledHaplotypes`: the haplotypes behind any allele of any call, once each.
    let mut called: Vec<usize> = Vec::new();
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
            for allele in &call.alleles {
                if let Some((_, indices)) = mapper.iter().find(|(a, _)| a == allele) {
                    for &h in indices {
                        if !called.iter().any(|&c| haplotypes[c] == haplotypes[h]) {
                            called.push(h);
                        }
                    }
                }
            }
            // `makeAnnotatedCall`: the annotations read the matrix the call was genotyped with,
            // over a reference context whose window is the padded reference.
            let call = match annotator.as_mut() {
                Some(annotator) => {
                    // `prepareReadAlleleLikelihoodsForAnnotation`: the genotyping matrix, with the
                    // filtered reads overlapping the call added at no likelihood.
                    let mut for_annotation = marginal.clone();
                    for (s, reads) in filtered_reads.iter().enumerate() {
                        let overlapping: Vec<BamRecord> = reads
                            .iter()
                            .filter(|read| {
                                overlap.overlaps(
                                    &overlap.contig,
                                    read_utils::start(read),
                                    read_utils::end(read),
                                )
                            })
                            .cloned()
                            .collect();
                        for_annotation.add_evidence(s, &overlapping, 0.0);
                    }
                    let site = Site {
                        likelihoods: &for_annotation,
                        window: (i64::from(ref_loc.start), reference),
                        overlaps: Vec::new(),
                        dbsnp: None,
                        resources: Vec::new(),
                    };
                    annotator
                        .engine
                        .annotate_context(&call, &site, annotator.random)?
                }
                None => call,
            };
            let call = if call.alleles.len() == merged_allele_count {
                call
            } else {
                reverse_trim_alleles(&call).map_err(refused)?
            };
            calls.push(call);
        }
    }
    if arguments.do_physical_phasing {
        let called: Vec<&Haplotype> = called.iter().map(|&h| &haplotypes[h]).collect();
        return phase_calls(calls, &called);
    }
    Ok(calls)
}

/// `PhaseGroup`: which alternate index of a het genotype holds the phased alternate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PhaseGroup {
    Phase01,
    Phase10,
}

impl PhaseGroup {
    fn description(self) -> &'static str {
        match self {
            PhaseGroup::Phase01 => "0|1",
            PhaseGroup::Phase10 => "1|0",
        }
    }

    fn alt_allele_index(self) -> usize {
        match self {
            PhaseGroup::Phase01 => 1,
            PhaseGroup::Phase10 => 0,
        }
    }

    fn other(self) -> PhaseGroup {
        match self {
            PhaseGroup::Phase01 => PhaseGroup::Phase10,
            PhaseGroup::Phase10 => PhaseGroup::Phase01,
        }
    }
}

/// `isSiteSpecificAltAllele`: neither the reference, `<NON_REF>` nor `*`.
fn is_site_specific_alt(allele: &Allele) -> bool {
    !(allele.is_reference()
        || allele.display_string() == "<NON_REF>"
        || allele.display_string() == htsjdk_vcf::allele::SPAN_DEL_STRING)
}

/// `AssemblyBasedCallerUtils.phaseCalls(calls, calledHaplotypes)`.
pub fn phase_calls(
    calls: Vec<VariantContext>,
    called: &[&Haplotype],
) -> Result<Vec<VariantContext>, EngineError> {
    // `constructHaplotypeMapping`: a call with exactly one site-specific alternate, tied to the
    // called haplotypes whose event map holds it at the call's start.
    let haplotype_map: Vec<Vec<usize>> = calls
        .iter()
        .map(|call| {
            let alts: Vec<&Allele> = call
                .alternate_alleles()
                .iter()
                .filter(|a| is_site_specific_alt(a))
                .collect();
            if alts.len() != 1 {
                return Vec::new();
            }
            let alt = alts[0];
            (0..called.len())
                .filter(|&h| {
                    called[h].event_map().is_some_and(|map| {
                        map.events()
                            .any(|e| i64::from(e.start()) == call.start && e.alt_allele() == alt)
                    })
                })
                .collect()
        })
        .collect();

    // `constructPhaseSetMapping`.
    let mut available: Vec<usize> = Vec::new();
    for set in &haplotype_map {
        for &h in set {
            if !available.contains(&h) {
                available.push(h);
            }
        }
    }
    let total = available.len();
    let mut mapping: Vec<Option<(usize, PhaseGroup)>> = vec![None; calls.len()];
    let mut counter = 0usize;
    let contains_all = |a: &[usize], b: &[usize]| b.iter().all(|x| a.contains(x));
    'outer: for i in 0..calls.len().saturating_sub(1) {
        let with_call = &haplotype_map[i];
        if with_call.is_empty() {
            continue;
        }
        let call_on_all = with_call.len() == total;
        let mut call_available = with_call.clone();
        for j in i + 1..calls.len() {
            let with_comp = &haplotype_map[j];
            if with_comp.is_empty() {
                continue;
            }
            let comp_on_all = with_comp.len() == total;
            if (with_call.len() == with_comp.len() && contains_all(with_call, with_comp))
                || (call_on_all && contains_all(&call_available, with_comp))
                || comp_on_all
            {
                if mapping[i].is_none() {
                    if mapping[j].is_some() {
                        mapping = vec![None; calls.len()];
                        break 'outer;
                    }
                    mapping[i] = Some((counter, PhaseGroup::Phase01));
                    mapping[j] = Some((counter, PhaseGroup::Phase01));
                    call_available.retain(|h| with_comp.contains(h));
                    counter += 1;
                } else if mapping[j].is_none() {
                    mapping[j] = mapping[i];
                }
            } else if with_call.len() + with_comp.len() == total
                && !with_call.iter().any(|h| with_comp.contains(h))
            {
                if mapping[i].is_none() {
                    if mapping[j].is_some() {
                        mapping = vec![None; calls.len()];
                        break 'outer;
                    }
                    mapping[i] = Some((counter, PhaseGroup::Phase01));
                    mapping[j] = Some((counter, PhaseGroup::Phase10));
                    counter += 1;
                } else if mapping[j].is_none() {
                    let (group, phase) = mapping[i].expect("mapped");
                    mapping[j] = Some((group, phase.other()));
                }
            }
        }
    }

    // `constructPhaseGroups`, up to the number of distinct groups.
    let mut groups: Vec<usize> = mapping.iter().flatten().map(|(g, _)| *g).collect();
    groups.sort_unstable();
    groups.dedup();
    let mut phased = calls.clone();
    for count in 0..groups.len() {
        let indexes: Vec<usize> = (0..calls.len())
            .filter(|&i| mapping[i].is_some_and(|(g, _)| g == count))
            .collect();
        if indexes.len() < 2 {
            return Err(EngineError::Runtime {
                class: "IllegalStateException".to_string(),
                message: "Somehow we have a group of phased variants that has fewer than 2 members"
                    .to_string(),
            });
        }
        let first = &calls[indexes[0]];
        let unique_id = format!(
            "{}_{}_{}",
            first.start,
            first.reference().display_string(),
            first.alternate_alleles()[0].display_string()
        );
        let phase_set = first.start;
        for &i in &indexes {
            let (_, phase) = mapping[i].expect("mapped");
            phased[i] = phase_vc(&calls[i], &unique_id, phase, phase_set);
        }
    }
    Ok(phased)
}

/// `htsjdk Genotype.isHet()`: every allele called, and not all the same.
fn is_het(genotype: &Genotype) -> bool {
    if genotype.alleles.is_empty() || genotype.alleles.iter().any(Allele::is_no_call) {
        return false;
    }
    genotype.alleles.iter().any(|a| *a != genotype.alleles[0])
}

/// `phaseVC`: every genotype phased with PID, PGT and PS, a het one reversed when its phased
/// alternate index holds no site-specific alternate.
fn phase_vc(vc: &VariantContext, id: &str, phase: PhaseGroup, phase_set: i64) -> VariantContext {
    let mut out = vc.clone();
    let genotypes: Vec<Genotype> = vc
        .genotypes
        .iter()
        .map(|g| {
            let mut genotype = g.clone();
            if is_het(g) && !is_site_specific_alt(&g.alleles[phase.alt_allele_index()]) {
                genotype.alleles.reverse();
            }
            genotype.phased = true;
            set_extended(&mut genotype, "PID", Value::Str(id.to_string()));
            set_extended(
                &mut genotype,
                "PGT",
                Value::Str(phase.description().to_string()),
            );
            set_extended(&mut genotype, "PS", Value::Int(phase_set));
            genotype
        })
        .collect();
    out.genotypes = GenotypesContext::new(genotypes);
    out
}

/// `GenotypeBuilder.attribute(key, value)`: replaced in place, or appended.
fn set_extended(genotype: &mut Genotype, key: &str, value: Value) {
    match genotype.extended.iter_mut().find(|(k, _)| k == key) {
        Some((_, slot)) => *slot = value,
        None => genotype.extended.push((key.to_string(), value)),
    }
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
