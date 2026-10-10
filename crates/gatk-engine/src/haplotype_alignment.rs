//! Ported from `org.broadinstitute.hellbender.utils.haplotype.Haplotype` (GATK 4.6.2.0): the
//! fields the assembler sets on a haplotype after it is built, and their accessors.
//!
//! `Haplotype` itself keeps its identity (bases, reference flag, uniqueness value) in
//! `haplotype.rs`; this is the rest of its state, which `equals` never reads: the genome location,
//! the CIGAR against the reference haplotype, the offset of the haplotype in the padded reference
//! (`alignmentStartHapwrtRef`, an index into bases held in memory and not a contig position), the
//! score of the path it came from, and the k-mer size of that path's graph.

use htsjdk_bam::cigar::Cigar;

use crate::cigar_builder::{CigarBuilder, CigarError};
use crate::haplotype::Haplotype;
use crate::interval::SimpleInterval;

/// The assembler's state on a `Haplotype`.
#[derive(Debug, Clone)]
pub struct HaplotypeAlignment {
    genome_location: Option<SimpleInterval>,
    cigar: Option<Cigar>,
    alignment_start_hap_wrt_ref: i32,
    score: f64,
    kmer_size: i32,
}

impl Default for HaplotypeAlignment {
    /// The field initialisers: `score` starts at `Double.NaN`, not 0.
    fn default() -> Self {
        HaplotypeAlignment {
            genome_location: None,
            cigar: None,
            alignment_start_hap_wrt_ref: 0,
            score: f64::NAN,
            kmer_size: 0,
        }
    }
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
    pub fn set_cigar(&mut self, cigar: &Cigar) -> Result<(), HaplotypeCigarError> {
        let mut builder = CigarBuilder::new(false);
        for element in &cigar.elements {
            builder
                .add(*element)
                .map_err(HaplotypeCigarError::Builder)?;
        }
        let made = builder.make(false).map_err(HaplotypeCigarError::Builder)?;
        let length = self.as_allele().len();
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

    /// `setKmerSize(int)`.
    pub fn set_kmer_size(&mut self, kmer_size: i32) {
        self.alignment.kmer_size = kmer_size;
    }
}
