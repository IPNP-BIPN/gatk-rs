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
//! With `--apply-bqd` or `--apply-frd` each sample is genotyped by the DRAGEN models
//! ([`gatk_engine::dragen_genotypes`]), over the reads whose SOFT-unclipped span overlaps the
//! context and the reads the HMM disqualified. Those disqualified reads are one list shared by
//! every event of the region and narrowed in place by each event's overlap, as `marginalize` and
//! `retainEvidence` share and narrow it in the reference. With a DRAGstr table the priors are
//! resolved per event from the STR at the site, and an indel call is annotated with the period,
//! the repeat count and the table's values there.
//!
//! Not ported here: contamination downsampling, and the reduction of sites with more alleles than
//! the genotype count allows (refused).

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_mapping::{
    create_allele_mapper, make_merged_variant_context, replace_span_dels,
    variants_from_active_haplotypes,
};
use gatk_engine::dragen_genotypes::{indel_prior, DragenGenotypesModel, DragenSample};
use gatk_engine::dragstr::{DragstrParams, DragstrReferenceAnalyzer};
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

use crate::genotyping_engine::{EngineError, GenotypePriors, GenotypingEngine};
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
    /// Reference-confidence mode: `<NON_REF>` added to every event and to its likelihoods, and
    /// the reducible annotations written raw.
    pub emit_reference_confidence: bool,
    /// The reference-confidence model's switches, read only in that mode.
    pub ref_confidence: crate::reference_confidence_model::RefConfidenceArguments,
    /// `RampedHaplotypeCallerEngine` in place of HaplotypeCaller's, which `callRegion` follows
    /// where the two part.
    pub ramped_engine: bool,
    /// `--apply-bqd` and `--apply-frd`: the DRAGEN genotype models in place of the independent
    /// samples model, and the soft-unclipped read overlap that comes with them.
    pub dragen_models: Option<DragenModelArguments>,
    /// `--dragstr-params-path`: DRAGstr's genotype priors and the STR annotations.
    pub dragstr: Option<DragstrGenotyping>,
    /// `--flow-filter-alleles`, with its thresholds.
    pub filter_alleles: Option<crate::allele_filtering::AlleleFilterArguments>,
}

/// The DRAGEN genotype models' switches.
#[derive(Debug, Clone)]
pub struct DragenModelArguments {
    pub bqd: bool,
    pub frd: bool,
    /// `--max-effective-depth-adjustment-for-frd`.
    pub max_effective_depth_adjustment: i32,
}

/// What `HaplotypeCallerGenotypingEngine` reads of DRAGstr.
#[derive(Debug, Clone)]
pub struct DragstrGenotyping {
    pub params: DragstrParams,
    /// Not `--dont-use-dragstr-priors`.
    pub use_priors: bool,
    /// `--dragstr-het-hom-ratio`.
    pub het_hom_ratio: i32,
    /// `--heterozygosity` and `--indel-heterozygosity`, not yet in log10.
    pub snp_heterozygosity: f64,
    pub indel_heterozygosity: f64,
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
            emit_reference_confidence: false,
            ref_confidence: Default::default(),
            ramped_engine: false,
            dragen_models: None,
            dragstr: None,
            filter_alleles: None,
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
        &std::collections::BTreeSet::new(),
        None,
        &[],
        None,
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
    suspicious_locations: &std::collections::BTreeSet<i32>,
    pre_filtering: Option<(&AlleleLikelihoods<BamRecord, Haplotype>, &[Haplotype])>,
    given: &[gatk_engine::event_map::Event],
    mut bamout: Option<&mut BamoutTags>,
) -> Result<Vec<VariantContext>, EngineError> {
    let refused = |message: String| EngineError::Runtime {
        class: "IllegalArgumentException".to_string(),
        message,
    };
    build_event_maps_for_haplotypes(haplotypes, reference, ref_loc, arguments.max_mnp_distance)
        .map_err(|e| refused(e.message()))?;
    let starts = event_start_positions(haplotypes).map_err(|e| refused(e.message()))?;
    // `constructDragstrReferenceSTRAnalyzerIfNecessary`: only with DRAGstr's priors in use and an
    // indel among the haplotypes' events.
    let dragstrs = match &arguments.dragstr {
        Some(dragstr)
            if dragstr.use_priors
                && !starts.is_empty()
                && haplotypes.iter().any(|h| {
                    h.event_map().is_some_and(|events| {
                        events.events().any(|e| {
                            contains_inline_indel(
                                e.ref_allele(),
                                std::slice::from_ref(e.alt_allele()),
                            )
                        })
                    })
                }) =>
        {
            let first = *starts.first().expect("a start");
            let last = *starts.last().expect("a start");
            Some(
                DragstrReferenceAnalyzer::of(
                    reference,
                    first - ref_loc.start,
                    last - ref_loc.start + 2,
                    dragstr.params.maximum_period(),
                )
                .map_err(|m| EngineError::Runtime {
                    class: "IndexOutOfBoundsException".to_string(),
                    message: m,
                })?,
            )
        }
        _ => None,
    };
    let dragen_model = match &arguments.dragen_models {
        Some(models) => Some(
            DragenGenotypesModel::new(
                models.bqd,
                models.frd,
                arguments.informative_read_overlap_margin,
                models.max_effective_depth_adjustment,
            )
            .map_err(|message| EngineError::Runtime {
                class: "GATKException".to_string(),
                message,
            })?,
        ),
        None => None,
    };
    // `marginalize` hands every event's matrix the SAME list of disqualified reads, and
    // `retainEvidence` narrows that list in place, so each event sees the disqualified reads that
    // overlapped every earlier event's interval as well as its own.
    let mut shared_filtered: Vec<Vec<BamRecord>> = (0..samples.len())
        .map(|s| read_likelihoods.filtered_evidence(s).to_vec())
        .collect();
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
        let collapsed = replaced.iter().any(|c| c.collapsed);
        let Some(mut merged) =
            make_merged_variant_context(&replaced).map_err(|e| refused(format!("{e:?}")))?
        else {
            continue;
        };
        let mut merged_allele_count = merged.alleles.len();
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
        let qualifies = |read: &BamRecord| {
            if dragen_model.is_some() {
                // `softUnclippedReadOverlapsInterval`.
                let soft_start = read_utils::soft_start(read);
                let soft_end = read_utils::soft_end(read);
                soft_start <= overlap.end && soft_end >= overlap.start && soft_start <= soft_end
            } else {
                overlap.overlaps(
                    &overlap.contig,
                    read_utils::start(read),
                    read_utils::end(read),
                )
            }
        };
        marginal.retain_evidence(qualifies);
        for reads in shared_filtered.iter_mut() {
            reads.retain(qualifies);
        }
        if arguments.emit_reference_confidence {
            // `addNonRefSymbolicAllele` on the merged context, `addNonReferenceAllele` on the
            // matrix, which fills its row from the concrete alleles.
            merged
                .alleles
                .push(Allele::from_str("<NON_REF>", false).expect("a symbolic allele"));
            marginal.add_non_reference_allele();
            merged_allele_count += 1;
        }

        let offset = loc - ref_loc.start;
        let dragen = match &dragen_model {
            Some(model) => {
                let params = arguments
                    .dragstr
                    .as_ref()
                    .map(|d| d.params.clone())
                    .unwrap_or_else(DragstrParams::default_params);
                let api = indel_prior(dragstrs.as_ref(), &params, offset)
                    .map_err(|e| refused(format!("{e:?}")))?;
                Some(DragenEvent {
                    model,
                    filtered: &shared_filtered,
                    subset_start: overlap.start,
                    reference,
                    offset,
                    api,
                })
            }
            None => None,
        };
        let genotypes =
            genotypes_for_event(&marginal, &merged.alleles, samples, ploidy, dragen.as_ref())?;
        // `resolveGenotypePriorCalculator`.
        if let Some(dragstr) = &arguments.dragstr {
            let log10_snp = jmath::math::log10(dragstr.snp_heterozygosity);
            let log10_indel = jmath::math::log10(dragstr.indel_heterozygosity);
            let priors = if !dragstr.use_priors {
                GenotypePriors::assuming_hw(log10_snp, log10_indel)
            } else {
                match &dragstrs {
                    None => GenotypePriors::given_het_to_hom_ratio(
                        log10_snp,
                        log10_indel,
                        jmath::math::log10(gatk_engine::math_utils::java_min(
                            dragstr.snp_heterozygosity,
                            dragstr.indel_heterozygosity,
                        )),
                        f64::from(dragstr.het_hom_ratio),
                    ),
                    Some(strs) => {
                        let period = strs.period(offset + 1).map_err(|e| refused(e.message))?;
                        let repeats = strs
                            .repeat_length(offset + 1)
                            .map_err(|e| refused(e.message))?;
                        GenotypePriors::given_dragstr_params(
                            &dragstr.params,
                            period,
                            repeats,
                            log10_snp,
                            f64::from(dragstr.het_hom_ratio),
                        )
                        .map_err(|e| refused(format!("{e:?}")))?
                    }
                }
            };
            engine.set_genotype_priors(priors);
        }
        let mut vc = VariantContext::new(
            &merged.contig,
            i64::from(merged.start),
            merged.alleles.clone(),
        );
        vc.stop = i64::from(merged.end);
        vc.genotypes = GenotypesContext::new(genotypes);
        // `allelesConsistentWithGivenAlleles(givenAlleles, mergedVC)`.
        let forced = crate::given_alleles::alleles_consistent_with_given_alleles(
            given,
            &vc.alleles[0],
            &vc.alleles[1..],
        );
        if let Some(call) = engine.calculate_genotypes_forcing(&vc, &forced)? {
            for allele in &call.alleles {
                if let Some((_, indices)) = mapper.iter().find(|(a, _)| a == allele) {
                    for &h in indices {
                        if !called.iter().any(|&c| haplotypes[c] == haplotypes[h]) {
                            called.push(h);
                        }
                    }
                }
            }
            // `annotateReadLikelihoodsWithSupportedAlleles(call, readAlleleLikelihoods)`, for the
            // bamout: the matrix the annotations read, marginalized to the call's alleles.
            if let Some(tags) = bamout.as_deref_mut() {
                let mut for_tags = marginal.clone();
                if call.alleles.len() != for_tags.number_of_alleles() {
                    for_tags.update_non_ref_allele_likelihoods(&call.alleles);
                }
                tags.add_supported_alleles(&call, &for_tags)
                    .map_err(|e| refused(format!("{e:?}")))?;
            }
            // `makeAnnotatedCall`: the annotations read the matrix the call was genotyped with,
            // over a reference context whose window is the padded reference.
            let call = match annotator.as_mut() {
                Some(annotator) => {
                    // `prepareReadAlleleLikelihoodsForAnnotation`: the genotyping matrix, with the
                    // filtered reads overlapping the call added at no likelihood.
                    let mut for_annotation = marginal.clone();
                    if call.alleles.len() != for_annotation.number_of_alleles() {
                        for_annotation.update_non_ref_allele_likelihoods(&call.alleles);
                    }
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
                        // `preFilteringAlleleLikelihoods`, which without `--flow-filter-alleles`
                        // is the matrix genotyped here, its haplotypes the ones given event maps.
                        haplotypes: Some(match pre_filtering {
                            Some((likelihoods, haplotypes)) => {
                                crate::variant_annotator_engine::HaplotypeMatrix::Reads(
                                    likelihoods,
                                    haplotypes,
                                )
                            }
                            None => crate::variant_annotator_engine::HaplotypeMatrix::Reads(
                                read_likelihoods,
                                &*haplotypes,
                            ),
                        }),
                    };
                    annotator.engine.annotate_context_raw(
                        &call,
                        &site,
                        annotator.random,
                        arguments.emit_reference_confidence,
                    )?
                }
                None => call,
            };
            // `EXT_COLLAPSED_TAG`, from the merged context, after the annotations.
            let mut call = call;
            if collapsed {
                set_attribute(&mut call, "XC", Value::Str("1".to_string()));
            }
            let mut call = if call.alleles.len() == merged_allele_count {
                call
            } else {
                reverse_trim_alleles(&call).map_err(refused)?
            };
            if suspicious_locations.contains(&loc) {
                set_attribute(
                    &mut call,
                    "SUSP_NOISY_ADJACENT_TP_VARIANT",
                    Value::Bool(true),
                );
            }
            // `annotateVariantContextWithDragstrParametersUsed`.
            if let (Some(strs), Some(dragstr)) = (&dragstrs, &arguments.dragstr) {
                if contains_inline_indel(&call.alleles[0], &call.alleles[1..]) {
                    let period = strs.period(offset + 1).map_err(|e| refused(e.message))?;
                    let repeats = strs
                        .repeat_length(offset + 1)
                        .map_err(|e| refused(e.message))?;
                    let table = |v: Result<f64, gatk_engine::dragstr::DragstrError>| {
                        v.map(|v| Value::Str(gatk_engine::java_format::format_decimals(v, 1)))
                            .map_err(|e| refused(format!("{e:?}")))
                    };
                    let values = vec![
                        table(dragstr.params.gop(period, repeats))?,
                        table(dragstr.params.gcp(period, repeats))?,
                        table(dragstr.params.api(period, repeats))?,
                    ];
                    set_attribute(
                        &mut call,
                        "DRAGstrInfo",
                        Value::List(vec![
                            Value::Int(i64::from(period)),
                            Value::Int(i64::from(repeats)),
                        ]),
                    );
                    set_attribute(&mut call, "DRAGstrParams", Value::List(values));
                }
            }
            calls.push(call);
        }
    }
    if let Some(tags) = bamout {
        tags.called = called.clone();
    }
    if arguments.do_physical_phasing {
        let called: Vec<&Haplotype> = called.iter().map(|&h| &haplotypes[h]).collect();
        return phase_calls(calls, &called);
    }
    Ok(calls)
}

/// What a bamout needs of genotyping: the haplotypes called, in the order the reference's
/// `HashSet` was filled, and the `XA` tag each read collects, by read.
#[derive(Debug, Clone, Default)]
pub struct BamoutTags {
    /// Indices into the genotyped haplotypes.
    pub called: Vec<usize>,
    /// `XA` by read name and flags.
    pub supported: Vec<((String, u16), String)>,
}

impl BamoutTags {
    /// `annotateReadLikelihoodsWithSupportedAlleles(vc, likelihoods)`: each read informative for
    /// one of the call's alleles gains `contig:start=index`, after a comma when it has one.
    pub fn add_supported_alleles(
        &mut self,
        call: &VariantContext,
        likelihoods: &AlleleLikelihoods<BamRecord>,
    ) -> Result<(), gatk_engine::allele_likelihoods::LikelihoodsError> {
        let new_to_old: Vec<(Allele, Vec<Allele>)> =
            gatk_annotation::depth_per_allele::marginalisation_order(&call.alleles)
                .into_iter()
                .map(|a| (a.clone(), vec![a]))
                .collect();
        let marginal = likelihoods.marginalize(&new_to_old)?;
        for sample in 0..marginal.number_of_samples() {
            let evidence = marginal.sample_evidence(sample).unwrap_or(&[]);
            for best in marginal.best_alleles_breaking_ties_for_sample(sample, None) {
                if !best.is_informative() {
                    continue;
                }
                let Some(allele) = best.allele else { continue };
                let index = call
                    .alleles
                    .iter()
                    .position(|a| *a == allele)
                    .map_or(-1, |i| i as i64);
                let read = &evidence[best.evidence_index];
                let text = format!("{}:{}={index}", call.contig, call.start);
                self.append((read.read_name.clone(), read.flags), &text);
            }
        }
        Ok(())
    }

    /// The same over fragments, `annotateReadLikelihoodsWithSupportedAlleles(vc, likelihoods,
    /// Fragment::getReads)`: every read of an informative fragment gains the tag.
    pub fn add_supported_alleles_fragments(
        &mut self,
        call: &VariantContext,
        likelihoods: &AlleleLikelihoods<gatk_engine::fragment::Fragment>,
    ) -> Result<(), gatk_engine::allele_likelihoods::LikelihoodsError> {
        let new_to_old: Vec<(Allele, Vec<Allele>)> =
            gatk_annotation::depth_per_allele::marginalisation_order(&call.alleles)
                .into_iter()
                .map(|a| (a.clone(), vec![a]))
                .collect();
        let marginal = likelihoods.marginalize(&new_to_old)?;
        for sample in 0..marginal.number_of_samples() {
            let evidence = marginal.sample_evidence(sample).unwrap_or(&[]);
            for best in marginal.best_alleles_breaking_ties_for_sample(sample, None) {
                if !best.is_informative() {
                    continue;
                }
                let Some(allele) = best.allele else { continue };
                let index = call
                    .alleles
                    .iter()
                    .position(|a| *a == allele)
                    .map_or(-1, |i| i as i64);
                let text = format!("{}:{}={index}", call.contig, call.start);
                for read in &evidence[best.evidence_index].reads {
                    self.append((read.read_name.clone(), read.flags), &text);
                }
            }
        }
        Ok(())
    }

    fn append(&mut self, key: (String, u16), text: &str) {
        match self.supported.iter_mut().find(|(k, _)| *k == key) {
            Some((_, value)) => {
                value.push_str(", ");
                value.push_str(text);
            }
            None => self.supported.push((key, text.to_string())),
        }
    }

    /// The `XA` a read collected.
    pub fn supported_for(&self, read: &BamRecord) -> Option<&str> {
        self.supported
            .iter()
            .find(|((name, flags), _)| *name == read.read_name && *flags == read.flags)
            .map(|(_, value)| value.as_str())
    }
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
    dragen: Option<&DragenEvent<'_>>,
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
        let likelihoods = match dragen {
            None => compute_log10_genotype_likelihoods(ploidy, &matrix),
            Some(event) => event
                .model
                .sample_log10_likelihoods(
                    ploidy,
                    &genotyping_alleles,
                    &DragenSample {
                        likelihoods: &matrix,
                        reads: marginal.sample_evidence(s).unwrap_or(&[]),
                        filtered: &event.filtered[s],
                    },
                    event.subset_start,
                    event.reference,
                    event.offset,
                    event.api,
                )
                .map_err(|message| EngineError::Runtime {
                    class: "GATKException".to_string(),
                    message,
                })?,
        };
        let mut genotype = Genotype::new(sample, vec![Allele::no_call(); ploidy]);
        genotype.pl = Some(gls_to_pls(&likelihoods));
        genotypes.push(genotype);
    }
    Ok(genotypes)
}

/// What the DRAGEN models read of one event beyond its matrix.
struct DragenEvent<'a> {
    model: &'a DragenGenotypesModel,
    filtered: &'a [Vec<BamRecord>],
    /// The start of the interval `setVariantCallingSubsetUsed` records.
    subset_start: i32,
    reference: &'a [u8],
    offset: i32,
    api: f64,
}

/// `GATKVariantContextUtils.containsInlineIndel(ref, alts)`: an alternate that is neither
/// symbolic nor the spanning deletion and differs from the reference in length.
fn contains_inline_indel(reference: &Allele, alternates: &[Allele]) -> bool {
    let length = |a: &Allele| if a.is_symbolic() { 0 } else { a.len() };
    alternates.iter().any(|alt| {
        !alt.is_symbolic() && alt.display_string() != "*" && length(alt) != length(reference)
    })
}

/// `VariantContextBuilder.attribute(key, value)`: replaced in place, or appended.
fn set_attribute(vc: &mut VariantContext, key: &str, value: Value) {
    match vc.attributes.iter_mut().find(|(k, _)| k == key) {
        Some((_, slot)) => *slot = value,
        None => vc.attributes.push((key.to_string(), value)),
    }
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
