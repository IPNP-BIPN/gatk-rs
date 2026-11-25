//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils`
//! (`getVariantsFromActiveHaplotypes`, `makeMergedVariantContext`, `createAlleleMapper`), with
//! `HaplotypeCallerGenotypingEngine.replaceSpanDels` and the allele side of
//! `GATKVariantContextUtils.simpleMerge` (GATK
//! 4.7.0.0; the three Java files are unchanged since
//! 4.6.2.0): from assembled haplotypes to the alleles
//! genotyped at one event start, and which haplotypes support each.
//!
//! At a start, the distinct events overlapping it (by `Event.equals`, the first haplotype to carry
//! one naming it `HC<index>`) become biallelic contexts; one that starts earlier, a spanning
//! deletion, becomes the reference base against `*`. The merge takes the longest reference allele
//! (two of one length must be equal), extends every alternate of a shorter one by the reference's
//! extra bases (`*` is left as it is), and keeps the alleles in arrival order, the reference first.
//! Each haplotype then goes to the reference if no event overlaps the start, to the merged allele
//! its event at the start maps to, or, for an event that starts earlier, to `*` (or the reference
//! when spanning events are not genotyped).
//!
//! What a reader would not guess:
//!
//!  * `simpleMerge` takes the merged allele list from a `LinkedHashSet` in which a context with a
//!    shorter reference adds its extended alternates **before** the reference; htsjdk's
//!    `VariantContext` then moves the reference to the front, so the order is the reference
//!    followed by the alternates as they arrived;
//!  * the merged context ends where its longest input ends, and among inputs of one length the
//!    first wins;
//!  * a haplotype with several events overlapping the start is credited once per event that starts
//!    there and maps, but stops at the first that starts earlier (the reference's own comment asks
//!    why that is a `break`: `//Why is there a break here? Shouldn't this be a continue?`).
//!
//! The contexts built from events carry no genotypes, filters, IDs or attributes (the flow-based
//! collapsed tag is not ported), so the rest of `simpleMerge`, which merges those, has nothing to
//! do here and is not ported.
//!
//! The `allele-mapping` suite (`allele_mapping.txt.gz`, 255 rows) measures all four steps; its
//! `mapper spanning-deletion 1020 * 2` and `mapper spanning-deletion-off 1020 A* 0,2` rows are the
//! two branches of the spanning event.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, MappingError>` is "a value or an error", and `?` after a call returns its error
//!   from the enclosing function at once (Java's `throws` at one call site). `impl
//!   From<AlleleError> for MappingError` is what lets `?` convert an allele error into this file's
//!   error type.
//! - `let Some(first) = contexts.first() else { return ... };` unwraps an `Option` or takes the
//!   early exit; there is no Java equivalent.
//! - `match` with `Some(r) if ...` arms tries each guard in order and takes the first that holds,
//!   like an `if`/`else if` chain over the cases.
//! - `EventContext { start: loc, ..c }` builds a copy of `c` with the named fields replaced.
//! - The `LinkedHashMap`s of the reference are `Vec`s of pairs here: insertion order is kept and a
//!   key is found by a linear search with `==`, which is `Allele.equals`.

use htsjdk_vcf::allele::{Allele, AlleleError, SPAN_DEL_STRING};

use crate::event_map::Event;
use crate::haplotype::Haplotype;

/// A biallelic context built from an event, as `Event.convertToVariantContext(source)` makes it.
#[derive(Debug, Clone, PartialEq)]
pub struct EventContext {
    pub source: String,
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub reference: Allele,
    pub alternate: Allele,
    /// The `XC` attribute of the event it was made from, which a spanning deletion keeps.
    pub collapsed: bool,
}

/// The merged context: its source, span, and alleles with the reference first.
#[derive(Debug, Clone, PartialEq)]
pub struct MergedContext {
    pub source: String,
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub alleles: Vec<Allele>,
}

impl MergedContext {
    pub fn reference(&self) -> &Allele {
        &self.alleles[0]
    }

    pub fn alternates(&self) -> &[Allele] {
        &self.alleles[1..]
    }
}

/// What the mapping refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MappingError {
    /// `IllegalStateException`: two reference alleles of one length that differ.
    InconsistentReferences(String),
    /// `NullPointerException`: a haplotype without an event map.
    NoEventMap,
    /// `IllegalStateException`: contexts that do not start together, or a reference longer than
    /// the one it is mapped to.
    IllegalState(String),
    Allele(AlleleError),
}

impl From<AlleleError> for MappingError {
    fn from(e: AlleleError) -> Self {
        MappingError::Allele(e)
    }
}

/// `Allele.SPAN_DEL`.
pub fn span_del() -> Allele {
    Allele::create(SPAN_DEL_STRING.as_bytes(), false).expect("the spanning deletion allele")
}

/// `getVariantsFromActiveHaplotypes(loc, haplotypes, includeSpanningEvents)`.
pub fn variants_from_active_haplotypes(
    loc: i32,
    haplotypes: &[Haplotype],
    include_spanning_events: bool,
) -> Result<Vec<EventContext>, MappingError> {
    // `final Set<Event> events = new HashSet<>();`, used only for `events.add(event)`'s answer, so
    // a list searched with `Event.equals` gives the same first-seen order and the same answers.
    let mut seen: Vec<Event> = Vec::new();
    let mut result = Vec::new();
    for (number, haplotype) in haplotypes.iter().enumerate() {
        let map = haplotype.event_map().ok_or(MappingError::NoEventMap)?;
        for event in map.overlapping_events(loc) {
            if !include_spanning_events && event.start() != loc {
                continue;
            }
            if seen.contains(event) {
                continue;
            }
            seen.push(event.clone());
            // `"HC" + hapNumber++` counts every haplotype, so the name is the index of the first
            // haplotype that carries the event, not a count of the events.
            result.push(EventContext {
                source: format!("HC{number}"),
                contig: event.contig().to_string(),
                start: event.start(),
                end: event.end(),
                reference: event.ref_allele().clone(),
                alternate: event.alt_allele().clone(),
                collapsed: event.is_collapsed(),
            });
        }
    }
    Ok(result)
}

/// `HaplotypeCallerGenotypingEngine.replaceSpanDels`: a context starting before `loc` becomes the
/// reference base at `loc` against `*`.
pub fn replace_span_dels(
    contexts: Vec<EventContext>,
    ref_allele: &Allele,
    loc: i32,
) -> Vec<EventContext> {
    contexts
        .into_iter()
        .map(|c| {
            if c.start == loc {
                c
            } else {
                EventContext {
                    start: loc,
                    end: loc,
                    reference: ref_allele.clone(),
                    alternate: span_del(),
                    ..c
                }
            }
        })
        .collect()
}

/// `makeMergedVariantContext`: `None` for no contexts.
pub fn make_merged_variant_context(
    contexts: &[EventContext],
) -> Result<Option<MergedContext>, MappingError> {
    let Some(first) = contexts.first() else {
        return Ok(None);
    };
    // `determineReferenceAllele`: the longest, and two of one length must be equal.
    let mut reference: Option<&Allele> = None;
    for c in contexts {
        reference = Some(match reference {
            None => &c.reference,
            Some(r) if r.len() < c.reference.len() => &c.reference,
            Some(r) if c.reference.len() < r.len() => r,
            // The reference's message goes on with `, %s vs. %s` and the two alleles, which this
            // one leaves out; at one start of one reference the branch is not reached.
            Some(r) if *r != c.reference => {
                return Err(MappingError::InconsistentReferences(format!(
                    "The provided variant file(s) have inconsistent references for the same \
                     position(s) at {}:{}",
                    c.contig, c.start
                )))
            }
            Some(r) => r,
        });
    }
    let reference = reference.expect("a context").clone();

    let mut alleles: Vec<Allele> = Vec::new();
    let mut add = |allele: Allele| {
        if !alleles.contains(&allele) {
            alleles.push(allele);
        }
    };
    let mut longest = first;
    for c in contexts {
        if longest.start != c.start {
            return Err(MappingError::IllegalState(
                "attempting to merge VariantContexts with different start sites".to_string(),
            ));
        }
        // `VariantContextUtils.getSize(vc) > getSize(longestVC)`, strictly, so of two equal lengths
        // the first stays.
        if c.end - c.start > longest.end - longest.start {
            longest = c;
        }
        // `resolveIncompatibleAlleles`: the context's own alleles, or its alternates extended and
        // then the reference.
        if c.reference == reference {
            add(c.reference.clone());
            add(c.alternate.clone());
        } else {
            for (_, mapped) in
                create_allele_mapping(&reference, &c.reference, std::slice::from_ref(&c.alternate))?
            {
                add(mapped);
            }
            add(reference.clone());
        }
    }
    // htsjdk's `VariantContext` puts the reference first. `remove` then `insert(0, ..)` keeps
    // the alternates in the order the `LinkedHashSet` received them.
    let position = alleles
        .iter()
        .position(|a| a.is_reference())
        .expect("the reference allele");
    let reference_allele = alleles.remove(position);
    alleles.insert(0, reference_allele);
    Ok(Some(MergedContext {
        source: first.source.clone(),
        contig: longest.contig.clone(),
        start: longest.start,
        end: longest.end,
        alleles,
    }))
}

/// `GATKVariantContextUtils.createAlleleMapping(refAllele, inputRef, inputAlts)`: each alternate,
/// and what it becomes against the longer reference.
pub fn create_allele_mapping(
    ref_allele: &Allele,
    input_ref: &Allele,
    input_alts: &[Allele],
) -> Result<Vec<(Allele, Allele)>, MappingError> {
    if ref_allele.len() < input_ref.len() {
        return Err(MappingError::IllegalState(format!(
            "BUG: inputRef={} is longer than refAllele={}",
            input_ref.display_string(),
            ref_allele.display_string()
        )));
    }
    if ref_allele.len() == input_ref.len() {
        return Ok(input_alts.iter().map(|a| (a.clone(), a.clone())).collect());
    }
    // `Arrays.copyOfRange(refAllele.getBases(), inputRef.length(), refAllele.length())`.
    let extra = &ref_allele.display_string().into_bytes()[input_ref.len()..];
    let star = span_del();
    let mut map = Vec::new();
    for a in input_alts {
        // `isNonSymbolicExtendableAllele(a)`, then `Allele.extend`; a reference or symbolic
        // alternate other than `*` is left out of the map altogether.
        if !(a.is_reference() || a.is_symbolic() || *a == star) {
            let mut bases = a.display_string().into_bytes();
            bases.extend_from_slice(extra);
            map.push((a.clone(), Allele::create(&bases, false)?));
        } else if *a == star {
            map.push((a.clone(), a.clone()));
        }
    }
    Ok(map)
}

/// `createAlleleMapper(mergedVC, loc, haplotypes, emitSpanningDels)`: each merged allele (the
/// symbolic ones left out) with the indices of the haplotypes supporting it, `*` appended when a
/// spanning event first needs it.
pub fn create_allele_mapper(
    merged: &MergedContext,
    loc: i32,
    haplotypes: &[Haplotype],
    emit_spanning_dels: bool,
) -> Result<Vec<(Allele, Vec<usize>)>, MappingError> {
    let reference = merged.reference().clone();
    let mut result: Vec<(Allele, Vec<usize>)> = vec![(reference.clone(), Vec::new())];
    for a in merged.alternates().iter().filter(|a| !a.is_symbolic()) {
        result.push((a.clone(), Vec::new()));
    }
    let slot = |result: &mut Vec<(Allele, Vec<usize>)>, allele: &Allele| {
        result.iter().position(|(a, _)| a == allele)
    };
    for (h, haplotype) in haplotypes.iter().enumerate() {
        // A partially determined haplotype supports nothing away from its determined position.
        if haplotype
            .partially_determined()
            .is_some_and(|pd| pd.determined_position != loc)
        {
            continue;
        }
        let map = haplotype.event_map().ok_or(MappingError::NoEventMap)?;
        let overlapping = map.overlapping_events(loc);
        if overlapping.is_empty() {
            result[0].1.push(h);
            continue;
        }
        for event in overlapping {
            if event.start() == loc {
                let ref_len = event.ref_allele().len();
                if ref_len == reference.len() {
                    if let Some(i) = slot(&mut result, event.alt_allele()) {
                        result[i].1.push(h);
                    }
                } else if ref_len < reference.len() {
                    let mapping = create_allele_mapping(
                        &reference,
                        event.ref_allele(),
                        std::slice::from_ref(event.alt_allele()),
                    )?;
                    // `spanningEventAlleleMappingToMergedVc.get(altAllele)`: the map has at most
                    // the one entry, and none when the alternate was left out of it.
                    if let Some((_, remapped)) = mapping.first() {
                        if let Some(i) = slot(&mut result, remapped) {
                            result[i].1.push(h);
                        }
                    }
                }
            } else if emit_spanning_dels {
                let star = span_del();
                let i = match slot(&mut result, &star) {
                    Some(i) => i,
                    None => {
                        result.push((star, Vec::new()));
                        result.len() - 1
                    }
                };
                result[i].1.push(h);
                break;
            } else {
                result[0].1.push(h);
                break;
            }
        }
    }
    Ok(result)
}
