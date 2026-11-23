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
//! The indexed-FASTA plumbing (`.fai` parsing and the seek to the first base) is
//! [`htsjdk_bam::fasta_index`], the port's own reader. It was `noodles-fasta` until the
//! arithmetic was measured: `getSubsequenceAt` computes the byte offset from the `.fai`'s
//! bases-per-line and bytes-per-line columns and never scans for a newline, so two readers of the
//! format agree on a well-formed file and disagree on one whose index does not describe its own
//! lines. htsjdk-rs's `indexed-fasta` suite holds that measurement: sixteen answered queries and
//! three refusals over three files, one of them CRLF, which is read correctly only because the
//! index says the terminator is two bytes. The *semantics* above are ported and measured here,
//! because they are GATK's and not the format's.
//!
//! Where the two disagree is also declared rather than smoothed over: the reference throws for a
//! query past the end of a contig, for `start > stop`, for a start below 1 and for an unknown
//! contig, and this returns [`ReferenceError`] in exactly those cases.

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use htsjdk_bgzf::BgzfReader;

use htsjdk_bam::fasta_index::{FastaIndex, FastaIndexError, IndexedFasta};

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

/// An indexed FASTA whatever reader it reads through. `IndexedFasta::open` returns
/// `IndexedFasta<File>` up to htsjdk-rs#272 and `IndexedFasta<PositionalFile>` from it, so the
/// source names neither and compiles against either revision.
trait Fasta: Send + Sync {
    fn index(&self) -> &FastaIndex;
    fn query(&mut self, contig: &str, start: i64, stop: i64) -> Result<Vec<u8>, FastaIndexError>;
}

impl<R: Read + Seek + Send + Sync> Fasta for IndexedFasta<R> {
    fn index(&self) -> &FastaIndex {
        IndexedFasta::index(self)
    }

    fn query(&mut self, contig: &str, start: i64, stop: i64) -> Result<Vec<u8>, FastaIndexError> {
        IndexedFasta::query(self, contig, start, stop)
    }
}

/// `ReferenceFileSource`: an indexed FASTA queried by interval.
pub struct ReferenceFileSource {
    reader: Box<dyn Fasta>,
    lengths: Vec<(String, usize)>,
}

impl ReferenceFileSource {
    /// Open a FASTA and its `.fai`.
    pub fn open(path: &Path) -> Result<ReferenceFileSource, ReferenceError> {
        let reader: Box<dyn Fasta> = match BgzfFile::open(path)? {
            Some(bgzf) => Box::new(IndexedFasta::new(bgzf, index_beside(path)?)),
            None => {
                Box::new(IndexedFasta::open(path).map_err(|e| ReferenceError::Io(e.message()))?)
            }
        };
        let lengths = reader
            .index()
            .entries
            .iter()
            .map(|entry| (entry.name.clone(), entry.size as usize))
            .collect();
        Ok(ReferenceFileSource { reader, lengths })
    }

    fn length_of(&self, contig: &str) -> Option<usize> {
        self.lengths
            .iter()
            .find(|(name, _)| name == contig)
            .map(|(_, length)| *length)
    }

    /// `getSequenceDictionary().getSequences()`: every contig and its length, in the `.fai`'s
    /// order, which is the dictionary order every interval list is sorted into.
    pub fn sequences(&self) -> &[(String, usize)] {
        &self.lengths
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
        let mut bases = self
            .reader
            .query(contig, start as i64, stop as i64)
            .map_err(|e| ReferenceError::Io(e.message()))?;
        upper_case_and_flatten_iupac(&mut bases);
        Ok(bases)
    }
}

/// The `.fai` beside `path`, as `IndexedFasta::open` reads it.
fn index_beside(path: &Path) -> Result<FastaIndex, ReferenceError> {
    let fai = with_suffix(path, ".fai");
    let text = std::fs::read_to_string(&fai)
        .map_err(|e| ReferenceError::Io(format!("{}: {e}", fai.display())))?;
    FastaIndex::parse(&text).map_err(|e| ReferenceError::Io(e.message()))
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(suffix);
    PathBuf::from(name)
}

/// htsjdk's `BlockCompressedIndexedFastaSequenceFile`: a bgzipped FASTA read at the uncompressed
/// offsets its `.fai` gives, through the `.gzi` that maps each block's uncompressed start to its
/// compressed one. `ReferenceSequenceFileFactory` picks it whenever the FASTA is block-compressed.
struct BgzfFile {
    file: File,
    /// `(compressed, uncompressed)` start of every block, the implicit first one included.
    blocks: Vec<(u64, u64)>,
    position: u64,
    reader: Option<BgzfReader<File>>,
}

impl BgzfFile {
    /// `None` for a FASTA that is not block-compressed.
    fn open(path: &Path) -> Result<Option<BgzfFile>, ReferenceError> {
        let io = |e: io::Error| ReferenceError::Io(format!("{}: {e}", path.display()));
        let mut file = File::open(path).map_err(io)?;
        let mut magic = [0u8; 18];
        let n = file.read(&mut magic).map_err(io)?;
        if !htsjdk_bgzf::read::is_valid_block_header(&magic[..n]) {
            return Ok(None);
        }
        let gzi = with_suffix(path, ".gzi");
        let bytes = std::fs::read(&gzi)
            .map_err(|e| ReferenceError::Io(format!("{}: {e}", gzi.display())))?;
        let word = |i: usize| {
            bytes
                .get(8 * i..8 * i + 8)
                .map(|w| u64::from_le_bytes(w.try_into().expect("eight bytes")))
                .ok_or_else(|| ReferenceError::Io(format!("{}: truncated", gzi.display())))
        };
        let count = word(0)? as usize;
        let mut blocks = vec![(0, 0)];
        for i in 0..count {
            blocks.push((word(1 + 2 * i)?, word(2 + 2 * i)?));
        }
        Ok(Some(BgzfFile {
            file,
            blocks,
            position: 0,
            reader: None,
        }))
    }
}

impl Read for BgzfFile {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.reader.is_none() {
            let block = self
                .blocks
                .partition_point(|&(_, start)| start <= self.position)
                - 1;
            let (compressed, uncompressed) = self.blocks[block];
            let mut file = self.file.try_clone()?;
            file.seek(SeekFrom::Start(compressed))?;
            let mut reader = BgzfReader::new(file);
            let skip = self.position - uncompressed;
            if io::copy(&mut (&mut reader).take(skip), &mut io::sink())? != skip {
                return Err(io::ErrorKind::UnexpectedEof.into());
            }
            self.reader = Some(reader);
        }
        let n = self.reader.as_mut().expect("opened above").read(buf)?;
        self.position += n as u64;
        Ok(n)
    }
}

impl Seek for BgzfFile {
    /// Only `SeekFrom::Start`, the one `IndexedFasta::query` makes.
    fn seek(&mut self, to: SeekFrom) -> io::Result<u64> {
        let SeekFrom::Start(position) = to else {
            return Err(io::ErrorKind::Unsupported.into());
        };
        if position != self.position {
            self.reader = None;
            self.position = position;
        }
        Ok(position)
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

    /// A bgzipped FASTA spanning several blocks answers every query as its plain twin does.
    #[test]
    fn a_bgzipped_fasta_reads_as_its_plain_twin() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("gatk-rs-bgzf-fasta-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (mut text, mut fai) = (Vec::new(), String::new());
        let mut state = 7u64;
        for (name, length) in [("a", 150_000usize), ("b", 333)] {
            text.extend_from_slice(format!(">{name}\n").as_bytes());
            fai += &format!("{name}\t{length}\t{}\t60\t61\n", text.len());
            for i in 0..length {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                text.push(b"ACGTacgtN"[(state >> 33) as usize % 9]);
                if i % 60 == 59 || i + 1 == length {
                    text.push(b'\n');
                }
            }
        }
        let plain = dir.join("ref.fasta");
        let packed = dir.join("ref.fasta.gz");
        std::fs::write(&plain, &text).unwrap();
        std::fs::write(with_suffix(&plain, ".fai"), &fai).unwrap();
        std::fs::write(with_suffix(&packed, ".fai"), &fai).unwrap();
        let mut writer = htsjdk_bgzf::BgzfWriter::new(Vec::new());
        writer.write_all(&text).unwrap();
        let bytes = writer.into_inner().unwrap();
        std::fs::write(&packed, &bytes).unwrap();
        let (mut gzi, mut starts, mut uncompressed) = (Vec::new(), 0u64, 0u64);
        let mut blocks = htsjdk_bgzf::BgzfReader::new(&bytes[..]);
        while let Some(block) = blocks.next_block().unwrap() {
            if block.block_address > 0 && !block.data.is_empty() {
                gzi.extend_from_slice(&block.block_address.to_le_bytes());
                gzi.extend_from_slice(&uncompressed.to_le_bytes());
                starts += 1;
            }
            uncompressed += block.data.len() as u64;
        }
        assert!(starts >= 2, "the FASTA must span several blocks");
        std::fs::write(
            with_suffix(&packed, ".gzi"),
            [&starts.to_le_bytes()[..], &gzi].concat(),
        )
        .unwrap();

        let mut expected = ReferenceFileSource::open(&plain).unwrap();
        let mut actual = ReferenceFileSource::open(&packed).unwrap();
        for (contig, start, stop) in [
            ("b", 1, 333),
            ("a", 1, 150_000),
            ("a", 65_000, 66_000),
            ("a", 140_000, 140_001),
            ("a", 3, 3),
        ] {
            assert_eq!(
                actual.query(contig, start, stop),
                expected.query(contig, start, stop)
            );
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
