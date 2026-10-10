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
}
