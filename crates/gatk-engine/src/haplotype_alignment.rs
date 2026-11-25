//! Ported from `org.broadinstitute.hellbender.utils.haplotype.Haplotype` (GATK 4.7.0.0;
//! `Haplotype.java` is unchanged since 4.6.2.0): the fields the assembler sets on a haplotype
//! after it is built, and their accessors.
//!
//! `Haplotype` itself keeps its identity (bases, reference flag, uniqueness value) in
//! `haplotype.rs`; this is the rest of its state, which `equals` never reads: the genome location,
//! the CIGAR against the reference haplotype, the event map read off that CIGAR, the offset of the
//! haplotype in the padded reference (`alignmentStartHapwrtRef`, an index into bases held in
//! memory and not a contig position), the score of the path it came from, and the k-mer size of
//! that path's graph.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing", Java's nullable field checked by the compiler: the
//!   reference's `genomeLocation = null` and `eventMap = null` are `None` here;
//! - `Result<T, E>` is "a value or an error", how a reference exception is modelled;
//!   `.map_err(Variant)` wraps the inner error in this file's own error type, `.ok_or(e)` turns a
//!   `None` into the error `e`, and a trailing `?` returns the error to the caller at once;
//! - `let Some(x) = expr else { return ...; };` unwraps the value or takes the early exit;
//! - `impl Haplotype { ... }` adds methods to `Haplotype`, which is declared in another file:
//!   Rust lets a type's methods be split across files of its crate, which Java does not;
//! - `&self` is a read-only `this`, `&mut self` a writable one.

use htsjdk_bam::cigar::Cigar;

use crate::alignment_utils::{
    bases_covering_ref_interval, trim_cigar_by_reference, AlignmentError,
};
use crate::cigar_builder::{CigarBuilder, CigarError};
use crate::event_map::EventMap;
use crate::haplotype::Haplotype;
use crate::interval::SimpleInterval;

/// The assembler's state on a `Haplotype`: the `Haplotype` fields that are not its identity,
/// gathered in one struct that `Haplotype` holds, so that `haplotype.rs` stays the identity alone.
#[derive(Debug, Clone)]
pub struct HaplotypeAlignment {
    /// `Locatable genomeLocation = null`, held as the `SimpleInterval` every caller passes.
    genome_location: Option<SimpleInterval>,
    /// `Cigar cigar`, always as `setCigar` consolidated it.
    cigar: Option<Cigar>,
    /// `alignmentStartHapwrtRef`:
    /// `//NOTE: this is the offset to a supposed array of reference bases held in memory and has
    /// nothing to do with start positions`.
    alignment_start_hap_wrt_ref: i32,
    score: f64,
    /// `kmerSize`, debug information about the graph the haplotype came from.
    kmer_size: i32,
    event_map: Option<EventMap>,
    /// `isCollapsed`: the flow-based collapsing engine lengthened one of its hmers.
    collapsed: bool,
}

impl Default for HaplotypeAlignment {
    /// The field initialisers: `score` starts at `Double.NaN`, not 0.
    ///
    /// ```java
    /// private double score = Double.NaN;
    /// ```
    fn default() -> Self {
        HaplotypeAlignment {
            genome_location: None,
            cigar: None,
            alignment_start_hap_wrt_ref: 0,
            score: f64::NAN,
            kmer_size: 0,
            event_map: None,
            collapsed: false,
        }
    }
}

/// What `trim` refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HaplotypeTrimError {
    /// `NullPointerException`: no genome location or no CIGAR. The string names which field was
    /// missing (`&'static str` is a string constant built into the program).
    Missing(&'static str),
    /// `IllegalArgumentException`: a span the haplotype does not contain.
    NotContained,
    /// A refusal from `AlignmentUtils.getBasesCoveringRefInterval` or `trimCigarByReference`.
    Alignment(AlignmentError),
    /// A refusal from `setCigar` on the trimmed haplotype.
    Cigar(HaplotypeCigarError),
    /// A refusal from `new Haplotype(newBases, ...)`, the `Allele` constructor underneath it.
    Allele(htsjdk_vcf::allele::AlleleError),
}

/// What `setCigar` refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HaplotypeCigarError {
    /// The `CigarBuilder`'s own refusal.
    Builder(CigarError),
    /// `IllegalArgumentException("Read length <length> not equal to the read length of the cigar
    /// <cigar read length> <consolidated cigar>")`.
    ReadLength {
        length: usize,
        cigar_read_length: u32,
        cigar: String,
    },
}

impl Haplotype {
    /// `getGenomeLocation()`.
    pub fn genome_location(&self) -> Option<&SimpleInterval> {
        self.alignment.genome_location.as_ref()
    }

    /// `setGenomeLocation(Locatable)`.
    pub fn set_genome_location(&mut self, location: SimpleInterval) {
        self.alignment.genome_location = Some(location);
    }

    /// `getCigar()`.
    pub fn cigar(&self) -> Option<&Cigar> {
        self.alignment.cigar.as_ref()
    }

    /// `setCigar(Cigar)`: the cigar is consolidated through `new CigarBuilder(false)`, so adjacent
    /// elements of one operator merge, and its read length must be the haplotype's length.
    ///
    /// ```java
    /// this.cigar = new CigarBuilder(false).addAll(cigar).make();
    /// Utils.validateArg( this.cigar.getReadLength() == length(), () -> "Read length " + length()
    ///         + " not equal to the read length of the cigar " + cigar.getReadLength() + " "
    ///         + this.cigar);
    /// ```
    ///
    /// The message mixes the two cigars: the read length is the argument's, the text is the
    /// consolidated one's, and the error below carries them the same way. One difference of
    /// order: the reference assigns the field before it validates, so a cigar of the wrong read
    /// length is still stored when the exception leaves; here it leaves the field as it was.
    /// (The `validateArg` line is broken at token boundaries to fit this file's width.)
    pub fn set_cigar(&mut self, cigar: &Cigar) -> Result<(), HaplotypeCigarError> {
        let mut builder = CigarBuilder::new(false);
        for element in &cigar.elements {
            builder
                .add(*element)
                .map_err(HaplotypeCigarError::Builder)?;
        }
        // `make()` is `make(false)`: an empty cigar is refused.
        let made = builder.make(false).map_err(HaplotypeCigarError::Builder)?;
        let length = self.as_allele().len();
        // `as usize` widens the `u32` read length to the index type `len()` returns; no value
        // is lost.
        if made.read_length() as usize != length {
            return Err(HaplotypeCigarError::ReadLength {
                length,
                cigar_read_length: cigar.read_length(),
                cigar: made.to_text(),
            });
        }
        self.alignment.cigar = Some(made);
        Ok(())
    }

    /// `getAlignmentStartHapwrtRef()`.
    pub fn alignment_start_hap_wrt_ref(&self) -> i32 {
        self.alignment.alignment_start_hap_wrt_ref
    }

    /// `setAlignmentStartHapwrtRef(int)`.
    pub fn set_alignment_start_hap_wrt_ref(&mut self, value: i32) {
        self.alignment.alignment_start_hap_wrt_ref = value;
    }

    /// `getScore()`.
    pub fn score(&self) -> f64 {
        self.alignment.score
    }

    /// `setScore(double)`.
    pub fn set_score(&mut self, score: f64) {
        self.alignment.score = score;
    }

    /// `getKmerSize()`.
    pub fn kmer_size(&self) -> i32 {
        self.alignment.kmer_size
    }

    /// `getEventMap()`.
    pub fn event_map(&self) -> Option<&EventMap> {
        self.alignment.event_map.as_ref()
    }

    /// `isCollapsed()`.
    pub fn is_collapsed(&self) -> bool {
        self.alignment.collapsed
    }

    /// `setCollapsed(collapsed)`.
    pub fn set_collapsed(&mut self, collapsed: bool) {
        self.alignment.collapsed = collapsed;
    }

    /// `setEventMap(EventMap)`.
    pub fn set_event_map(&mut self, map: EventMap) {
        self.alignment.event_map = Some(map);
    }

    /// `trim(loc, ignoreRefState)`: the part of the haplotype aligned to `loc`, which it must
    /// contain, or `None` when an end of `loc` falls in a deletion or nothing but an insertion is
    /// left. A leading or trailing insertion is cut from the CIGAR, the reference flag is dropped
    /// when `ignore_ref_state` asks, and the score, k-mer size and the shifted alignment start are
    /// carried over.
    ///
    /// The reference's own note on the insertions (its two lines rewrapped here):
    ///
    /// ```java
    /// // note: trimCigarByReference does not remove leading or trailing indels, while
    /// // getBasesCoveringRefInterval does remove bases of leading and trailing insertions.  We
    /// // must remove leading and trailing insertions from the Cigar manually.
    /// ```
    ///
    /// The refusals come in the reference's order: no genome location, then `loc` not contained,
    /// then no cigar. The reference's `Utils.nonNull(loc, ...)` has no counterpart, because a Rust
    /// reference cannot be null. The trimmed haplotype has no event map, as the reference's has
    /// none.
    pub fn trim(
        &self,
        loc: &SimpleInterval,
        ignore_ref_state: bool,
    ) -> Result<Option<Haplotype>, HaplotypeTrimError> {
        let location = self
            .genome_location()
            .ok_or(HaplotypeTrimError::Missing("genomeLocation"))?;
        if !location.contains(loc) {
            return Err(HaplotypeTrimError::NotContained);
        }
        let cigar = self.cigar().ok_or(HaplotypeTrimError::Missing("cigar"))?;
        // Offsets into the haplotype's bases, not contig positions, as the reference's
        // `newStart` and `newStop` are.
        let new_start = loc.start - location.start;
        let new_stop = new_start + loc.end - loc.start;
        let bases = self.bases();
        // `getBasesCoveringRefInterval` returns `null` when an end falls in a deletion; that
        // `null` is the `None` the `else` branch takes, while a refusal is the error `?` passes on.
        let Some(new_bases) = bases_covering_ref_interval(new_start, new_stop, &bases, 0, cigar)
            .map_err(HaplotypeTrimError::Alignment)?
        else {
            return Ok(None);
        };
        if new_bases.is_empty() {
            return Ok(None);
        }
        let new_cigar = trim_cigar_by_reference(cigar, new_start, new_stop)
            .map_err(HaplotypeTrimError::Alignment)?
            .cigar;
        let leading_insertion = !new_cigar.elements[0].op.consumes_reference_bases();
        let trailing_insertion = !new_cigar.elements[new_cigar.elements.len() - 1]
            .op
            .consumes_reference_bases();
        // `usize::from(bool)` is 1 for `true` and 0 for `false`: the reference's
        // `leadingInsertion ? 1 : 0`, giving `firstIndexToKeepInclusive` and
        // `lastIndexToKeepExclusive`.
        let first = usize::from(leading_insertion);
        let last = new_cigar.elements.len() - usize::from(trailing_insertion);
        if last <= first {
            // The whole CIGAR is an insertion.
            return Ok(None);
        }
        // `elements[first..last]` is the half-open sub-range, `subList(first, last)`, and
        // `.to_vec()` copies it into a new list.
        let kept = if leading_insertion || trailing_insertion {
            Cigar::new(new_cigar.elements[first..last].to_vec())
        } else {
            new_cigar
        };
        let mut trimmed = Haplotype::new(&new_bases, !ignore_ref_state && self.is_reference())
            .map_err(HaplotypeTrimError::Allele)?;
        // `setCigar` runs the elements through a `CigarBuilder(false)`, as the reference's
        // `new CigarBuilder(false).addAll(...).make()` already has.
        trimmed
            .set_cigar(&kept)
            .map_err(HaplotypeTrimError::Cigar)?;
        trimmed.set_genome_location(loc.clone());
        trimmed.set_score(self.score());
        trimmed.set_kmer_size(self.kmer_size());
        trimmed.set_alignment_start_hap_wrt_ref(new_start + self.alignment_start_hap_wrt_ref());
        Ok(Some(trimmed))
    }

    /// `setKmerSize(int)`.
    pub fn set_kmer_size(&mut self, kmer_size: i32) {
        self.alignment.kmer_size = kmer_size;
    }
}
