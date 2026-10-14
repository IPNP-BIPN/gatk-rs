//! Force-calling, `--alleles`: ported from the given-allele paths of `HaplotypeCallerEngine`,
//! `Mutect2Engine`, `AssemblyResultSet.addGivenAlleles` and
//! `AssemblyBasedCallerUtils.allelesConsistentWithGivenAlleles` (GATK 4.6.2.0).
//!
//! A locus a given record overlaps is active whatever its pileup says. A region's given events are
//! the records over its padded span, filtered ones dropped unless `--force-call-filtered-alleles`,
//! each split into one event per alternate and kept when it starts at or after the region's
//! start: an event in the right padding stays. The trimmer sees them beside the assembled ones;
//! each one the assembly lacks is inserted into the five best haplotypes (the reference first,
//! then the best scores); and an alternate consistent with a given event (the same length change,
//! one allele a prefix of the other) is genotyped and emitted whatever its confidence.

use gatk_engine::assembly_result_set::AssemblyResultSet;
use gatk_engine::event_map::Event;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::smith_waterman::{SmithWatermanAligner, SwParameters};
use gatk_engine::variant_context_utils::{trim_alleles, Allele as TrimAllele, Variant};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

/// `AssemblyBasedCallerUtils.NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO`.
const NUM_HAPLOTYPES_TO_INJECT_INTO: usize = 5;

/// `vc.isNotFiltered()`, or `--force-call-filtered-alleles`.
fn usable(vc: &VariantContext, force_call_filtered: bool) -> bool {
    force_call_filtered || !vc.is_filtered()
}

fn overlaps(vc: &VariantContext, contig: &str, start: i32, end: i32) -> bool {
    vc.contig == contig && vc.start <= i64::from(end) && i64::from(start) <= vc.stop
}

/// `isActive`'s test: a usable given record over the locus.
pub fn forced_at(
    alleles: &[VariantContext],
    contig: &str,
    position: i32,
    force_call_filtered: bool,
) -> bool {
    alleles
        .iter()
        .any(|vc| overlaps(vc, contig, position, position) && usable(vc, force_call_filtered))
}

/// `splitVariantContextToEvents(vc, false, BEST_MATCH_TO_ORIGINAL, false)`: one event per
/// alternate, each right trimmed when the record had to be split.
fn split_to_events(vc: &VariantContext) -> Result<Vec<Event>, String> {
    if vc.alleles.len() <= 1 {
        return Ok(Vec::new());
    }
    let local = |a: &Allele| TrimAllele::new(a.display_string().as_bytes(), a.is_reference());
    let mut out = Vec::new();
    for alt in &vc.alleles[1..] {
        let pair = if vc.alleles.len() == 2 {
            vec![local(&vc.alleles[0]), local(alt)]
        } else {
            let variant = Variant {
                contig: vc.contig.clone(),
                start: vc.start as i32,
                stop: vc.stop as i32,
                alleles: vec![local(&vc.alleles[0]), local(alt)],
                genotypes: Vec::new(),
                attributes: Vec::new(),
            };
            trim_alleles(&variant, false, true)
                .map_err(|e| format!("{e:?}"))?
                .alleles
        };
        let make =
            |a: &TrimAllele| Allele::create(&a.bases, a.is_reference).map_err(|e| format!("{e:?}"));
        out.push(
            Event::new(
                &vc.contig,
                vc.start as i32,
                make(&pair[0])?,
                make(&pair[1])?,
            )
            .map_err(|e| e.message())?,
        );
    }
    Ok(out)
}

/// The region's given events: the usable records over `padded`, split, and those starting at or
/// after `span_start`, in the order the file holds them.
pub fn given_events(
    alleles: &[VariantContext],
    contig: &str,
    padded: &SimpleInterval,
    span_start: i32,
    force_call_filtered: bool,
) -> Result<Vec<Event>, String> {
    let mut out = Vec::new();
    for vc in alleles {
        if !overlaps(vc, contig, padded.start, padded.end) || !usable(vc, force_call_filtered) {
            continue;
        }
        for event in split_to_events(vc)? {
            if event.start() >= span_start {
                out.push(event);
            }
        }
    }
    Ok(out)
}

/// `allelesAreConsistent(ref1, alt1, ref2, alt2)`.
fn alleles_are_consistent(ref1: &Allele, alt1: &Allele, ref2: &Allele, alt2: &Allele) -> bool {
    if alt1.is_symbolic() || alt2.is_symbolic() {
        return false;
    }
    let diff1 = alt1.len() as i64 - ref1.len() as i64;
    let diff2 = alt2.len() as i64 - ref2.len() as i64;
    let (a1, a2) = (alt1.display_string(), alt2.display_string());
    diff1 == diff2
        && if a1.len() < a2.len() {
            a2.as_bytes()[..a1.len()] == *a1.as_bytes()
        } else {
            a1.as_bytes()[..a2.len()] == *a2.as_bytes()
        }
}

/// `allelesConsistentWithGivenAlleles(givenAlleles, mergedVC)`.
pub fn alleles_consistent_with_given_alleles(
    given: &[Event],
    reference: &Allele,
    alternates: &[Allele],
) -> Vec<Allele> {
    if given.is_empty() {
        return Vec::new();
    }
    alternates
        .iter()
        .filter(|allele| {
            given
                .iter()
                .any(|g| alleles_are_consistent(g.ref_allele(), g.alt_allele(), reference, allele))
        })
        .cloned()
        .collect()
}

/// `addGivenAlleles(givenAlleles, maxMnpDistance, aligner, swParameters)` on a trimmed set.
pub fn add_given_alleles(
    assembly: &mut AssemblyResultSet,
    given: &[Event],
    max_mnp_distance: i32,
    aligner: &dyn SmithWatermanAligner,
    haplotype_to_reference_sw: &SwParameters,
) -> Result<(), String> {
    if given.is_empty() {
        return Ok(());
    }
    let error = |e: gatk_engine::assembly_result_set::ResultSetError| format!("{e:?}");
    let reference = assembly
        .reference_haplotype()
        .cloned()
        .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
    let assembled = assembly.variation_events(max_mnp_distance).map_err(error)?;
    // The reference first, then the best scores, five of them, as a list.
    let mut ranked: Vec<Haplotype> = assembly.haplotype_list().to_vec();
    ranked.sort_by(|a, b| {
        b.is_reference()
            .cmp(&a.is_reference())
            .then_with(|| b.score().total_cmp(&a.score()))
    });
    ranked.truncate(NUM_HAPLOTYPES_TO_INJECT_INTO);
    for event in given {
        if assembled.iter().any(|e| e == event)
            || crate::pileup_based_alleles::is_symbolic(event.alt_allele())
        {
            continue;
        }
        for base in &ranked {
            if let Some(inserted) = crate::pileup_based_alleles::make_haplotype_with_inserted_event(
                base,
                &reference,
                event,
                aligner,
                haplotype_to_reference_sw,
            )? {
                assembly.add(inserted).map_err(error)?;
            }
        }
    }
    assembly.variation_events(max_mnp_distance).map_err(error)?;
    Ok(())
}
