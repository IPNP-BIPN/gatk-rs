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
//!
//! # Reading this file without reading Rust
//!
//! - `&[VariantContext]` is a borrowed, read-only view of a list; the caller keeps it;
//! - `Result<T, String>` is "a value or an error message", how this port carries a Java
//!   exception; `?` after a call hands its error straight back to the caller, Java's `throws`
//!   at that one call site, and `.map_err(..)` converts the error type first;
//! - `|a| ...` is a closure, Java's lambda `a -> ...`;
//! - `.iter().any(..)` is `stream().anyMatch(..)`, `.filter(..).cloned().collect()` is
//!   `stream().filter(..).collect(..)` with each kept element copied;
//! - `as i32` / `as i64` converts between integer widths; narrowing keeps the low bits;
//! - `if let Some(x) = f()` runs the block only when `f` returned a value, not `None` (null).

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

/// `vc.isNotFiltered()`, or `--force-call-filtered-alleles`:
/// `.filter(vc -> hcArgs.forceCallFiltered || vc.isNotFiltered())`.
fn usable(vc: &VariantContext, force_call_filtered: bool) -> bool {
    force_call_filtered || !vc.is_filtered()
}

/// Whether the record's `[start, stop]` meets `[start, end]` on the same contig, both closed.
///
/// Stands in for the feature query `features.getValues(hcArgs.alleles, ref)` (one locus) and
/// `features.getValues(hcArgs.alleles)` (the padded region). The `i64::from` widens the `i32`
/// bounds to the record's `i64` positions without loss.
fn overlaps(vc: &VariantContext, contig: &str, start: i32, end: i32) -> bool {
    vc.contig == contig && vc.start <= i64::from(end) && i64::from(start) <= vc.stop
}

/// `isActive`'s test: a usable given record over the locus.
///
/// When it holds, `HaplotypeCallerEngine.isActive` returns
/// `new ActivityProfileState(ref.getInterval(), 1.0)` before looking at the pileup.
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
///
/// `splitVariantContextToEvents` is `splitVariantContextToBiallelics(..)` mapped through
/// `Event::ofWithoutAttributes`. A biallelic record passes through as it is; a multiallelic one
/// becomes one (reference, alternate) pair per alternate, then `trim_alleles(.., false, true)`
/// trims the pair's common suffix only (`trimLeft` is `false`). Genotypes and attributes are left
/// empty because an `Event` keeps neither.
fn split_to_events(vc: &VariantContext) -> Result<Vec<Event>, String> {
    if vc.alleles.len() <= 1 {
        return Ok(Vec::new());
    }
    // Converts a VCF allele to the trimmer's own allele type, by display string.
    let local = |a: &Allele| TrimAllele::new(a.display_string().as_bytes(), a.is_reference());
    let mut out = Vec::new();
    // `&vc.alleles[1..]` is a view of the alternates: every allele after the reference.
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
///
/// ```java
/// .filter(event -> event.getStart() >= region.getSpan().getStart())
/// ```
///
/// The quote is cut before its trailing comment. An event in the left padding is dropped, one
/// in the right padding is kept.
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
///
/// Two (ref, alt) pairs from different records "could be the same" when they change the length
/// by the same amount and one alternate's bases are a prefix of the other's:
///
/// ```java
/// return (sizeDiff1 == sizeDiff2) && (alt1.length() < alt2.length() ?
///         alt1.basesMatch(Arrays.copyOf(alt2.getBases(), alt1.length())) :
///         alt2.basesMatch(Arrays.copyOf(alt1.getBases(), alt2.length())));
/// ```
///
/// The lengths are widened to `i64` before subtracting, since `usize` (an unsigned index) would
/// wrap below zero where Java's `int` goes negative. Both alternates are non-symbolic by then,
/// so the display string is the bases.
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

/// `allelesConsistentWithGivenAlleles(givenAlleles, mergedVC)`: the merged record's alternates
/// that match at least one given event.
///
/// Java collects into a `Set` (`Collectors.toSet()`); this returns the alternates in record
/// order, without repeats only because a record's alternates are already distinct.
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

/// `AssemblyResultSet.addGivenAlleles(givenAlleles, maxMnpDistance, aligner, swParameters)` on a
/// trimmed set.
///
/// Each given event the assembly did not find, and that is not symbolic, is inserted into each of
/// the five base haplotypes, and every insertion that succeeds is added to the set. The variation
/// events are then regenerated, as `regenerateVariationEvents(maxMnpDistance)` does.
///
/// The base haplotypes:
///
/// `Comparator.comparing(Haplotype::isReference).thenComparingDouble(hap -> hap.getScore())`
/// then `.reversed()`, then
/// `.limit(AssemblyBasedCallerUtils.NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO)`.
/// `reversed()` applies to the whole comparator, so the reference comes first and then the scores
/// descend. The Rust comparator swaps `a` and `b` for the same effect. Both sorts are stable, so
/// equal scores keep the haplotype list's order.
///
/// The Java looks the given event up only among the assembled events with the same start
/// (`assembledEventsByStart`, keyed by `Event::getStart`); here it is looked up among all of
/// them, which is the same answer because equal events have equal starts.
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
    // One closure that formats any result-set error, reused by every `map_err` below.
    let error = |e: gatk_engine::assembly_result_set::ResultSetError| format!("{e:?}");
    // `getReferenceHaplotype()` can return null in Java, and the NullPointerException would come
    // from the first `makeHaplotypeWithInsertedEvent`; here a missing reference is refused before
    // the loop. `.cloned()` copies the haplotype so `assembly` can be modified below.
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
    // `regenerateVariationEvents(maxMnpDistance)`. The list is discarded: the call is kept for
    // building the event map of every haplotype, the inserted ones included, and for its errors.
    assembly.variation_events(max_mnp_distance).map_err(error)?;
    Ok(())
}
