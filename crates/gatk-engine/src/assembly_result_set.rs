//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResult` and
//! `AssemblyResultSet` (GATK 4.7.0.0; both files are unchanged since 4.6.2.0): what one assembly
//! graph gave, and the haplotypes of a region the assembler hands to genotyping.
//!
//! The result set is a `LinkedHashSet<Haplotype>`: insertion-ordered, compared by
//! `Haplotype.equals` (bases, reference flag, uniqueness value), and an element added again keeps
//! its first position. In the sequence-graph mode the assembler adds haplotypes without their
//! assembly result, so the per-k-mer map the reference also keeps stays empty and is not ported
//! here; trimming, event maps and given alleles come with the bricks that need them.
//!
//! # Reading this file without reading Rust
//!
//! - the set is a `Vec<Haplotype>` (a growable list) searched with `==`, which is
//!   `Haplotype.equals`; that keeps the insertion order and the first-position rule of a
//!   `LinkedHashSet` without a hash;
//! - the reference haplotype is held as its index in that list (`Option<usize>`, "an index or
//!   nothing"), not as a second pointer to the same object, so it has to be kept in step when the
//!   list changes;
//! - `Result<T, ResultSetError>` is "a value or an error", how the reference's exceptions are
//!   modelled; a trailing `?` returns the error to the caller at once, and `.map_err(f)` wraps
//!   another module's error in this one's;
//! - `match x { Some(r) if r == index => ..., other => other }` tests a pattern and a condition
//!   together; the final `other` arm catches the rest, unchanged;
//! - `.collect::<Result<_, _>>()?` turns a list of results into one: the list of values, or the
//!   first error.
//!
//! A long Java line quoted here is broken at a token boundary to fit the file's width.

use crate::assembly_region::AssemblyRegion;
use crate::event_map::{build_event_maps_for_haplotypes, Event};
use crate::haplotype::Haplotype;
use crate::interval::SimpleInterval;
use crate::seq_graph::SeqGraph;

/// `AssemblyResult.Status`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyStatus {
    /// The reference was shorter than the k-mer size.
    Failed,
    /// The graph lost its reference source or sink while being cleaned up.
    JustAssembledReference,
    /// `/** Assembly succeeded, and the graph has some meaningful structure */`.
    AssembledSomeVariation,
}

impl AssemblyStatus {
    /// The constant's name, as `Status.toString()` prints it.
    pub fn name(self) -> &'static str {
        match self {
            AssemblyStatus::Failed => "FAILED",
            AssemblyStatus::JustAssembledReference => "JUST_ASSEMBLED_REFERENCE",
            AssemblyStatus::AssembledSomeVariation => "ASSEMBLED_SOME_VARIATION",
        }
    }
}

/// `AssemblyResult`: a status, the sequence graph it came with, and the k-mer size of the graph.
///
/// The reference also keeps the read threading graph, which only the linked de Bruijn mode reads
/// back; its k-mer size is all the sequence-graph mode asks of it.
#[derive(Debug, Clone)]
pub struct AssemblyResult {
    /// `getStatus()`.
    pub status: AssemblyStatus,
    /// `getSeqGraph()`, `null` for a `FAILED` result.
    pub seq_graph: Option<SeqGraph>,
    /// `getKmerSize()`, which throws on a `FAILED` result that has no graph.
    pub kmer_size: Option<usize>,
    /// `getDiscoveredHaplotypes()`.
    pub discovered_haplotypes: Vec<Haplotype>,
}

/// What the result set refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResultSetError {
    /// `IllegalStateException`: no reference haplotype to trim, or none left after trimming.
    IllegalState(&'static str),
    /// A haplotype the trim or the event maps refuse, carried as the refusal's text.
    Haplotype(String),
    /// `Utils.nonNull(h.getGenomeLocation(), "haplotype genomeLocation cannot be null")`.
    NoGenomeLocation,
    /// `IllegalStateException("the assembly-result-set already have a reference haplotype that is
    /// different")`.
    SecondReference,
}

/// `AssemblyResultSet`, as far as the assembler fills it.
///
/// Not kept: `assemblyResultByKmerSize`, `assemblyResultByHaplotype` and `kmerSizes` (see the
/// module note), the collapsing engine, and the `variationEvents` cache with its
/// `lastMaxMnpDistanceUsed` (see [`AssemblyResultSet::variation_events`]).
#[derive(Debug, Clone, Default)]
pub struct AssemblyResultSet {
    /// `haplotypes`, the `LinkedHashSet<Haplotype>`.
    haplotypes: Vec<Haplotype>,
    /// The index of `refHaplotype` in `haplotypes`.
    reference: Option<usize>,
    /// `variationPresent`, the raw flag; [`AssemblyResultSet::is_variation_present`] also asks
    /// for more than one haplotype.
    variation_present: bool,
    region_for_genotyping: Option<AssemblyRegion>,
    /// `fullReferenceWithPadding`: the reference bases of the padded span, which every
    /// haplotype's `alignmentStartHapwrtRef` indexes.
    full_reference_with_padding: Vec<u8>,
    /// `paddedReferenceLoc`: where `full_reference_with_padding` lies on the contig.
    padded_reference_loc: Option<SimpleInterval>,
    /// `isPartiallyDeterminedList`.
    partially_determined: bool,
    /// `!originalAssemblyHaps.isEmpty()`: the PD engine replaced the assembled haplotypes.
    overwritten: bool,
}

impl AssemblyResultSet {
    /// `new AssemblyResultSet()`: empty. `Default` fills every field with its zero value
    /// (empty lists, `None`, `false`).
    pub fn new() -> AssemblyResultSet {
        AssemblyResultSet::default()
    }

    /// `add(Haplotype)`: false if an equal haplotype is already there.
    ///
    /// ```java
    /// Utils.nonNull(h.getGenomeLocation(), "haplotype genomeLocation cannot be null");
    /// if (haplotypes.contains(h)) {
    ///     return false;
    /// }
    /// haplotypes.add(h);
    /// updateReferenceHaplotype(h);
    /// return true;
    /// ```
    ///
    /// This `add` leaves the variation flag alone; only `add(Haplotype, AssemblyResult)`, not
    /// ported, raises it. A refused second reference is already in the list when the error
    /// leaves, as it is in the reference, which adds before it checks.
    pub fn add(&mut self, haplotype: Haplotype) -> Result<bool, ResultSetError> {
        if haplotype.genome_location().is_none() {
            return Err(ResultSetError::NoGenomeLocation);
        }
        if self.haplotypes.contains(&haplotype) {
            return Ok(false);
        }
        let is_reference = haplotype.is_reference();
        self.haplotypes.push(haplotype);
        self.update_reference_haplotype(is_reference)?;
        Ok(true)
    }

    /// `updateReferenceHaplotype`, for the haplotype just pushed: it becomes the reference if it
    /// is one, and a second reference is refused. The two cannot be equal, because `add` has
    /// already turned an equal haplotype away.
    fn update_reference_haplotype(&mut self, is_reference: bool) -> Result<(), ResultSetError> {
        if !is_reference {
            return Ok(());
        }
        if self.reference.is_some() {
            return Err(ResultSetError::SecondReference);
        }
        self.reference = Some(self.haplotypes.len() - 1);
        Ok(())
    }

    /// `replaceAllHaplotypes`: cleared, then each added, the variation flag raised by any
    /// non-reference haplotype and never lowered. The flag is raised even for a haplotype `add`
    /// turned away as a duplicate, as in the reference, whose loop does not look at `add`'s
    /// answer.
    pub fn replace_all_haplotypes(&mut self, list: &[Haplotype]) -> Result<(), ResultSetError> {
        self.haplotypes.clear();
        self.reference = None;
        for haplotype in list {
            self.add(haplotype.clone())?;
            if !haplotype.is_reference() {
                self.variation_present = true;
            }
        }
        Ok(())
    }

    /// `setPartiallyDeterminedMode`.
    pub fn set_partially_determined_mode(&mut self, value: bool) {
        self.partially_determined = value;
    }

    /// `isPartiallyDeterminedList()`.
    pub fn is_partially_determined_list(&self) -> bool {
        self.partially_determined
    }

    /// `storeAssemblyHaplotypes()`, as far as anything reads it: whether there were any.
    ///
    /// The reference copies the set into `originalAssemblyHaps`, and its only reader is
    /// `hasOverwrittenHaps()`, `return !originalAssemblyHaps.isEmpty();`. Only that answer is
    /// kept.
    pub fn store_assembly_haplotypes(&mut self) {
        self.overwritten = !self.haplotypes.is_empty();
    }

    /// `hasOverwrittenHaps()`.
    pub fn has_overwritten_haps(&self) -> bool {
        self.overwritten
    }

    /// `removeHaplotype`: out of the list, the reference's place kept in step.
    ///
    /// `position` finds the first equal haplotype (the set holds at most one). The index of the
    /// reference haplotype shifts down by one when an earlier haplotype is removed. One difference
    /// in kind: the reference's `removeHaplotype` never touches `refHaplotype`, so removing the
    /// reference haplotype there leaves `getReferenceHaplotype()` answering it; here it answers
    /// `None`.
    pub fn remove_haplotype(&mut self, haplotype: &Haplotype) {
        if let Some(index) = self.haplotypes.iter().position(|h| h == haplotype) {
            self.haplotypes.remove(index);
            self.reference = match self.reference {
                Some(r) if r == index => None,
                Some(r) if r > index => Some(r - 1),
                other => other,
            };
        }
    }

    /// `getHaplotypeList()`, in insertion order. The reference returns a fresh copy; a read-only
    /// borrow gives the caller the same guarantee.
    pub fn haplotype_list(&self) -> &[Haplotype] {
        &self.haplotypes
    }

    /// `getHaplotypeCount()`.
    pub fn haplotype_count(&self) -> usize {
        self.haplotypes.len()
    }

    /// `getReferenceHaplotype()`.
    pub fn reference_haplotype(&self) -> Option<&Haplotype> {
        self.reference.map(|i| &self.haplotypes[i])
    }

    /// The position of `getReferenceHaplotype()` in `getHaplotypeList()`.
    pub fn reference_index(&self) -> Option<usize> {
        self.reference
    }

    /// `isVariationPresent()`: the flag, and more than one haplotype.
    /// `return variationPresent && haplotypes.size() > 1;`
    pub fn is_variation_present(&self) -> bool {
        self.variation_present && self.haplotypes.len() > 1
    }

    /// `getRegionForGenotyping()`.
    pub fn region_for_genotyping(&self) -> Option<&AssemblyRegion> {
        self.region_for_genotyping.as_ref()
    }

    /// The region for genotyping, for a caller that removes reads from it as `callRegion` does.
    pub fn region_for_genotyping_mut(&mut self) -> Option<&mut AssemblyRegion> {
        self.region_for_genotyping.as_mut()
    }

    /// `setRegionForGenotyping`.
    pub fn set_region_for_genotyping(&mut self, region: AssemblyRegion) {
        self.region_for_genotyping = Some(region);
    }

    /// `getFullReferenceWithPadding()`.
    pub fn full_reference_with_padding(&self) -> &[u8] {
        &self.full_reference_with_padding
    }

    /// `setFullReferenceWithPadding`.
    pub fn set_full_reference_with_padding(&mut self, bases: Vec<u8>) {
        self.full_reference_with_padding = bases;
    }

    /// `getPaddedReferenceLoc()`.
    pub fn padded_reference_loc(&self) -> Option<&SimpleInterval> {
        self.padded_reference_loc.as_ref()
    }

    /// `setPaddedReferenceLoc`.
    pub fn set_padded_reference_loc(&mut self, location: SimpleInterval) {
        self.padded_reference_loc = Some(location);
    }

    /// `getVariationEvents(maxMnpDistance)`: every haplotype's event map rebuilt against the padded
    /// reference, and their events in `HAPLOTYPE_EVENT_COMPARATOR` order (start, reference
    /// length, alternate bases), once each. The variation flag becomes whether any haplotype is not
    /// the reference.
    ///
    /// ```java
    /// public static final Comparator<Event> HAPLOTYPE_EVENT_COMPARATOR =
    ///         Comparator.comparingInt(Event::getStart)
    ///         .thenComparingInt(e -> e.refAllele().length())
    ///         .thenComparing(e -> e.altAllele());
    /// ```
    ///
    /// The reference collects into a `TreeSet` with that comparator, which keeps the first of
    /// the events it calls equal. A stable sort followed by `dedup_by` (drop each element equal to
    /// the one before it) keeps the same one. The alternate alleles are compared by display
    /// string: every event's alternate is a non-reference allele of plain bases, and htsjdk's
    /// `SimpleAllele.compareTo` orders two of those byte by byte, then by length, which for
    /// ASCII bases is the string order.
    ///
    /// The reference caches the events and rebuilds them only when the MNP distance changes, the
    /// cache is empty, or a non-reference haplotype has no event map; this rebuilds on every call.
    pub fn variation_events(
        &mut self,
        max_mnp_distance: i32,
    ) -> Result<Vec<Event>, ResultSetError> {
        let ref_loc = self
            .padded_reference_loc
            .clone()
            .ok_or(ResultSetError::IllegalState("no padded reference location"))?;
        build_event_maps_for_haplotypes(
            &mut self.haplotypes,
            &self.full_reference_with_padding,
            &ref_loc,
            max_mnp_distance,
        )
        .map_err(|e| ResultSetError::Haplotype(e.message()))?;
        let mut events: Vec<Event> = Vec::new();
        for haplotype in &self.haplotypes {
            for event in haplotype.event_map().expect("just built").events() {
                events.push(event.clone());
            }
        }
        events.sort_by(|a, b| {
            a.start()
                .cmp(&b.start())
                .then(a.ref_allele().len().cmp(&b.ref_allele().len()))
                .then(
                    a.alt_allele()
                        .display_string()
                        .cmp(&b.alt_allele().display_string()),
                )
        });
        events.dedup_by(|a, b| {
            a.start() == b.start()
                && a.ref_allele().len() == b.ref_allele().len()
                && a.alt_allele().display_string() == b.alt_allele().display_string()
        });
        self.variation_present = self.haplotypes.iter().any(|h| !h.is_reference());
        Ok(events)
    }

    /// `trimTo(trimmedAssemblyRegion)`: every haplotype cut to the region's padded span (one
    /// dropped when it cannot be, the reference refused), duplicates collapsed with the
    /// reference's original winning, the reference rebuilt as one, and the lot sorted by length
    /// and then bases. The variation flag is whether any **untrimmed** haplotype is not the
    /// reference.
    ///
    /// The reference trims every haplotype before it tests for a missing reference haplotype
    /// (`calculateOriginalByTrimmedHaplotypes` comes first, then `if (refHaplotype == null)`), so
    /// there a refused trim wins over `"refHaplotype is null"`; here the test comes first.
    pub fn trim_to(&self, region: AssemblyRegion) -> Result<AssemblyResultSet, ResultSetError> {
        if self.reference.is_none() {
            return Err(ResultSetError::IllegalState("refHaplotype is null"));
        }
        let span = region.padded_span().clone();
        // `trimDownHaplotypes`: a map from trimmed to original, keyed by `Haplotype.equals`.
        // The trimmed haplotypes are made with `ignoreRefState` true, so the reference's trim
        // equals any other that trims to the same bases, and the reference then takes the entry.
        let mut by_trimmed: Vec<(Haplotype, usize)> = Vec::new();
        for (index, haplotype) in self.haplotypes.iter().enumerate() {
            let trimmed = haplotype
                .trim(&span, true)
                .map_err(|e| ResultSetError::Haplotype(format!("{e:?}")))?;
            match trimmed {
                Some(trimmed) => match by_trimmed.iter().position(|(t, _)| *t == trimmed) {
                    Some(position) => {
                        if haplotype.is_reference() {
                            by_trimmed.remove(position);
                            by_trimmed.push((trimmed, index));
                        }
                    }
                    None => by_trimmed.push((trimmed, index)),
                },
                None if haplotype.is_reference() => {
                    return Err(ResultSetError::IllegalState(
                        "trimming eliminates the reference haplotype",
                    ))
                }
                None => {}
            }
        }
        let mut fixed: Vec<Haplotype> = by_trimmed
            .into_iter()
            .map(|(trimmed, original)| {
                if self.haplotypes[original].is_reference() {
                    // A new reference haplotype, which does not carry the k-mer size over.
                    let mut reference = Haplotype::new(&trimmed.bases(), true)
                        .map_err(|e| ResultSetError::Haplotype(e.to_string()))?;
                    if let Some(cigar) = trimmed.cigar() {
                        reference
                            .set_cigar(cigar)
                            .map_err(|e| ResultSetError::Haplotype(format!("{e:?}")))?;
                    }
                    if let Some(location) = trimmed.genome_location() {
                        reference.set_genome_location(location.clone());
                    }
                    reference.set_score(trimmed.score());
                    reference
                        .set_alignment_start_hap_wrt_ref(trimmed.alignment_start_hap_wrt_ref());
                    Ok(reference)
                } else {
                    Ok(trimmed)
                }
            })
            .collect::<Result<_, ResultSetError>>()?;
        // `SIZE_AND_BASE_ORDER`:
        // `Comparator.comparingInt((Haplotype hap) -> hap.getBases().length)
        //  .thenComparing(Allele::getBaseString)`. The reference sorts a `HashMap`'s keys, whose
        // order is a hash order; no two keys tie here (equal bases were merged above), so where
        // they start does not change where they end.
        fixed.sort_by(|a, b| {
            a.bases()
                .len()
                .cmp(&b.bases().len())
                .then_with(|| a.bases().cmp(&b.bases()))
        });
        let mut result = AssemblyResultSet::new();
        for haplotype in fixed {
            result.add(haplotype)?;
        }
        result.set_region_for_genotyping(region);
        result.set_full_reference_with_padding(self.full_reference_with_padding.clone());
        if let Some(location) = &self.padded_reference_loc {
            result.set_padded_reference_loc(location.clone());
        }
        result.variation_present = self.haplotypes.iter().any(|h| !h.is_reference());
        if result.reference.is_none() {
            return Err(ResultSetError::IllegalState(
                "missing reference haplotype in the trimmed set",
            ));
        }
        Ok(result)
    }

    /// The haplotypes, for a caller that builds their event maps.
    pub fn haplotypes_mut(&mut self) -> &mut [Haplotype] {
        &mut self.haplotypes
    }
}
