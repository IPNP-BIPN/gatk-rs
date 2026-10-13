//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResult` and
//! `AssemblyResultSet` (GATK 4.6.2.0): what one assembly graph gave, and the haplotypes of a
//! region the assembler hands to genotyping.
//!
//! The result set is a `LinkedHashSet<Haplotype>`: insertion-ordered, compared by
//! `Haplotype.equals` (bases, reference flag, uniqueness value), and an element added again keeps
//! its first position. In the sequence-graph mode the assembler adds haplotypes without their
//! assembly result, so the per-k-mer map the reference also keeps stays empty and is not ported
//! here; trimming, event maps and given alleles come with the bricks that need them.

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
    pub status: AssemblyStatus,
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
    /// A haplotype the trim or the event maps refuse.
    Haplotype(String),
    /// `Utils.nonNull(h.getGenomeLocation(), "haplotype genomeLocation cannot be null")`.
    NoGenomeLocation,
    /// `IllegalStateException("the assembly-result-set already have a reference haplotype that is
    /// different")`.
    SecondReference,
}

/// `AssemblyResultSet`, as far as the assembler fills it.
#[derive(Debug, Clone, Default)]
pub struct AssemblyResultSet {
    haplotypes: Vec<Haplotype>,
    /// The index of `refHaplotype` in `haplotypes`.
    reference: Option<usize>,
    variation_present: bool,
    region_for_genotyping: Option<AssemblyRegion>,
    full_reference_with_padding: Vec<u8>,
    padded_reference_loc: Option<SimpleInterval>,
}

impl AssemblyResultSet {
    pub fn new() -> AssemblyResultSet {
        AssemblyResultSet::default()
    }

    /// `add(Haplotype)`: false if an equal haplotype is already there.
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

    /// `updateReferenceHaplotype`, for the haplotype just pushed.
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
    /// non-reference haplotype and never lowered.
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

    /// `getHaplotypeList()`.
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
    pub fn trim_to(&self, region: AssemblyRegion) -> Result<AssemblyResultSet, ResultSetError> {
        if self.reference.is_none() {
            return Err(ResultSetError::IllegalState("refHaplotype is null"));
        }
        let span = region.padded_span().clone();
        // `trimDownHaplotypes`: a map from trimmed to original, keyed by `Haplotype.equals`.
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
        // `SIZE_AND_BASE_ORDER`.
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
