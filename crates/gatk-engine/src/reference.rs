//! Ported from `org.broadinstitute.hellbender.engine.ReferenceFileSource` and
//! `org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile` (GATK 4.6.2.0).
//!
//! # Why this is not "read the FASTA"
//!
//! The bases a GATK tool sees are **not** the bytes in the file. `CachingIndexedFastaSequenceFile`
//! defaults to `preserveCase = false` and `preserveIUPAC = false`, so every query comes back
//!
//!  * upper-cased, which erases soft-masking: `acgt` in the file is `ACGT` to every annotation;
//!  * with every IUPAC ambiguity code replaced by `N`, so `RYKMSWBDHV` all become `N`.
//!
//! Measured, not assumed: the golden shows `acgtNNNNacgt` in the file coming back as
//! `ACGTNNNNACGT`, and `ACGTRYKMSWBD` coming back as `ACGTNNNNNNNN`. A port that returned the
//! file's bytes, which is what any FASTA reader gives, would differ from the reference on every
//! soft-masked or ambiguous position in the genome, and those are not rare.
//!
//! # What is a dependency and what is a port
//!
//! The indexed-FASTA plumbing (`.fai` parsing, seeking to the right offset, skipping newlines)
//! is [`noodles_fasta`]: it is one well-tested implementation of a file format whose bytes are
//! unambiguous, and porting htsjdk's copy of it a second time would buy nothing measurable. The
//! *semantics* above are ported and measured, because they are GATK's and not the format's.
//!
//! Where the two disagree is also declared rather than smoothed over: the reference throws for a
//! query past the end of a contig, for `start > stop`, for a start below 1 and for an unknown
//! contig, and this returns [`ReferenceError`] in exactly those cases.

use std::path::Path;

use noodles_fasta as fasta;

/// What the reference throws rather than returning bases.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReferenceError {
    /// The dictionary does not declare this contig.
    UnknownContig(String),
    /// `start < 1`, `stop < start`, or a stop past the end of the contig.
    BadInterval,
    /// The FASTA or its index could not be read at all.
    Io(String),
}

/// `ReferenceFileSource`: an indexed FASTA queried by interval.
pub struct ReferenceFileSource {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<std::fs::File>>,
    lengths: Vec<(String, usize)>,
}

impl ReferenceFileSource {
    /// Open a FASTA and its `.fai`.
    pub fn open(path: &Path) -> Result<ReferenceFileSource, ReferenceError> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| ReferenceError::Io(e.to_string()))?;
        let lengths = reader
            .index()
            .as_ref()
            .iter()
            .map(|record| {
                (
                    String::from_utf8_lossy(record.name()).into_owned(),
                    record.length() as usize,
                )
            })
            .collect();
        Ok(ReferenceFileSource { reader, lengths })
    }

    fn length_of(&self, contig: &str) -> Option<usize> {
        self.lengths
            .iter()
            .find(|(name, _)| name == contig)
            .map(|(_, length)| *length)
    }

    /// `getSequenceDictionary().getSequence(contig).getSequenceLength()`, or `None` where the
    /// reference has no such contig and throws.
    pub fn sequence_length(&self, contig: &str) -> Option<usize> {
        self.length_of(contig)
    }

    /// `ReferenceDataSource.queryAndPrefetch(interval)`: the bases of `[start, stop]`, 1-based and
    /// inclusive, upper-cased and with IUPAC codes flattened to `N`.
    pub fn query(
        &mut self,
        contig: &str,
        start: i32,
        stop: i32,
    ) -> Result<Vec<u8>, ReferenceError> {
        let Some(length) = self.length_of(contig) else {
            return Err(ReferenceError::UnknownContig(contig.to_string()));
        };
        if start < 1 || stop < start || stop as usize > length {
            return Err(ReferenceError::BadInterval);
        }
        let region = format!("{contig}:{start}-{stop}")
            .parse()
            .map_err(|_| ReferenceError::BadInterval)?;
        let record = self
            .reader
            .query(&region)
            .map_err(|e| ReferenceError::Io(e.to_string()))?;
        let mut bases = record.sequence().as_ref().to_vec();
        upper_case_and_flatten_iupac(&mut bases);
        Ok(bases)
    }
}

/// `StringUtil.toUpperCase` then `BaseUtils.convertIUPACtoN`, which is what every query goes
/// through before a tool sees it.
///
/// The IUPAC set is the reference's own map: `N R Y M K W S B D H V`, in either case. Anything
/// that is not a base and not an ambiguity code is an error there (`errorOnBadReferenceBase`),
/// which no valid FASTA reaches.
pub fn upper_case_and_flatten_iupac(bases: &mut [u8]) {
    for base in bases.iter_mut() {
        *base = base.to_ascii_uppercase();
        if matches!(
            *base,
            b'N' | b'R' | b'Y' | b'M' | b'K' | b'W' | b'S' | b'B' | b'D' | b'H' | b'V'
        ) {
            *base = b'N';
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn soft_masking_and_ambiguity_do_not_survive_a_query() {
        let mut bases = b"acgtRYKMSWBDHVn".to_vec();
        upper_case_and_flatten_iupac(&mut bases);
        assert_eq!(&bases, b"ACGTNNNNNNNNNNN");
    }
}
