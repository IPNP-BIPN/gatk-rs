//! `PathSeqBuildKmers`, ported from the tool, `PSKmerUtils`, `PSKmerCollection`, `SVKmerShort` and
//! `SVKmerizer` (GATK 4.6.2.0).
//!
//! A host reference turned into the k-mer set PathSeq subtracts against. The containers the tool
//! writes, a hopscotch set or a Bloom filter, are not ported; what is here is what goes into them,
//! which is the whole of the tool's own arithmetic.
//!
//! # A k-mer is a long, and the first base is in the high bits
//!
//! Two bits per base, A=0 C=1 G=2 T=3, shifted up as the window advances:
//!
//! ```java
//! final long newV2 = ((valLow << 2) | (base.value & 3L)) & mask;
//! ```
//!
//! so the long is the whole k-mer and nothing else, and reading it back is a shift and a mask.
//!
//! # The set is canonical, which is why k must be odd
//!
//! ```java
//! if (((valLow >> kSize) & 1L) == 0) return this;
//! return reverseComplement(kSize);
//! ```
//!
//! A k-mer whose middle base is G or T is replaced by its reverse complement, so each strand of a
//! sequence produces the same set. The test is one bit of the middle base, which only means
//! anything when there IS a middle base, and an even size is refused by `canonical` itself rather
//! than by the argument parser: the tool accepts `--kmer-size 4` and dies inside the stream.
//!
//! # Masking happens after canonicalisation
//!
//! ```java
//! return val.canonical(kSize).mask(mask).getLong();
//! ```
//!
//! The mask is built once from the masked positions, counted from the START of the k-mer, and
//! ANDed in afterwards. So which bases a mask clears depends on which strand won the
//! canonicalisation, and two k-mers that differ only in a masked position collapse to one entry.
//!
//! # A bad base costs a whole window
//!
//! ```java
//! default: validBaseCount = -1;
//! ...
//! if ( ++validBaseCount == kSize ) return tmpKmer;
//! ```
//!
//! Anything that is not ACGT, in either case, sets the counter to -1, which the increment takes
//! back to zero: k good bases in a row are needed again. The k-mer's bits are NOT reset, but they
//! are all shifted out before the counter comes back, so it does not matter. `--kmer-spacing` works
//! on the same counter, the next window starting at `kSize - spacing` rather than at zero, so a
//! spacing of k gives non-overlapping windows.
//!
//! # And no empty set is ever written
//!
//! A k longer than every contig produces nothing, and the container refuses to be built from
//! nothing: the run dies with `Number of elements must be greater than 0`, which comes from the set
//! rather than from the tool.

use std::collections::BTreeSet;

/// `SVKmer.Base`: A=0, C=1, G=2, T=3.
fn base_value(base: u8) -> Option<u64> {
    match base {
        b'a' | b'A' => Some(0),
        b'c' | b'C' => Some(1),
        b'g' | b'G' => Some(2),
        b't' | b'T' => Some(3),
        _ => None,
    }
}

/// `SVKmer.reverseComplementByteValueAsLong`.
fn reverse_complement_byte(byte: u8) -> u64 {
    let value = byte as u64;
    (!(((value & 3) << 6)
        | (((value >> 2) & 3) << 4)
        | (((value >> 4) & 3) << 2)
        | ((value >> 6) & 3)))
        & 0xff
}

/// `SVKmer.reverseComplement(long)`, byte by byte.
pub fn reverse_complement_long(value: u64) -> u64 {
    let mut result = reverse_complement_byte((value & 0xff) as u8);
    let mut rest = value;
    for _ in 0..7 {
        rest >>= 8;
        result = (result << 8) | reverse_complement_byte((rest & 0xff) as u8);
    }
    result
}

/// What the run refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KmerError {
    /// `canonical`, which the tool reaches only once it is already streaming.
    EvenKmerSize,
    /// `PSUtils.parseMask`.
    InvalidMaskIndex { index: String },
    /// The container, when the reference produced no k-mer at all.
    EmptySet,
}

impl KmerError {
    pub fn java_class(&self) -> &str {
        "java.lang.IllegalArgumentException"
    }

    pub fn message(&self) -> String {
        match self {
            KmerError::EvenKmerSize => "Kmer length must be odd to canonicalize.".to_string(),
            KmerError::InvalidMaskIndex { index } => format!("Invalid kmer mask index: {index}"),
            KmerError::EmptySet => "Number of elements must be greater than 0".to_string(),
        }
    }
}

/// `SVKmerShort`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer(pub u64);

impl Kmer {
    /// `successor`, which drops the leading base and adds one at the end.
    pub fn successor(&self, base: u64, kmer_size: usize) -> Kmer {
        let window = (1u64 << (kmer_size * 2)) - 1;
        Kmer(((self.0 << 2) | (base & 3)) & window)
    }

    /// `reverseComplement`.
    pub fn reverse_complement(&self, kmer_size: usize) -> Kmer {
        let window = (1u64 << (kmer_size * 2)) - 1;
        let unused = 64 - kmer_size * 2;
        Kmer(reverse_complement_long(self.0 << unused) & window)
    }

    /// `canonical`, which is defined only for an odd size.
    pub fn canonical(&self, kmer_size: usize) -> Result<Kmer, KmerError> {
        if kmer_size.is_multiple_of(2) {
            return Err(KmerError::EvenKmerSize);
        }
        if (self.0 >> kmer_size) & 1 == 0 {
            Ok(*self)
        } else {
            Ok(self.reverse_complement(kmer_size))
        }
    }

    /// `mask`, which is a plain AND.
    pub fn mask(&self, mask: Kmer) -> Kmer {
        Kmer(self.0 & mask.0)
    }

    /// `toString(kSize)`, which spells the bases back out.
    pub fn bases(&self, kmer_size: usize) -> String {
        let mut out: Vec<u8> = Vec::with_capacity(kmer_size);
        let mut value = self.0;
        for _ in 0..kmer_size {
            out.push(b"ACGT"[(value & 3) as usize]);
            value >>= 2;
        }
        out.reverse();
        String::from_utf8(out).expect("ascii")
    }
}

/// `SVKmerShort.getMask`, whose positions are counted from the start of the k-mer.
pub fn get_mask(positions: &[u8], kmer_size: usize) -> Kmer {
    let mut mask: u64 = 0;
    for position in positions {
        mask |= 3u64 << (2 * (kmer_size - *position as usize - 1));
    }
    Kmer(!mask)
}

/// `PSUtils.parseMask`, whose empty string is no mask at all.
pub fn parse_mask(argument: &str, kmer_size: usize) -> Result<Vec<u8>, KmerError> {
    if argument.is_empty() {
        return Ok(Vec::new());
    }
    let mut positions = Vec::new();
    for field in argument.split(',') {
        let index = field
            .parse::<i32>()
            .map_err(|_| KmerError::InvalidMaskIndex {
                index: field.to_string(),
            })?;
        if index < 0 || index >= kmer_size as i32 {
            return Err(KmerError::InvalidMaskIndex {
                index: field.to_string(),
            });
        }
        positions.push(index as u8);
    }
    Ok(positions)
}

/// `SVKmerizer`, which emits a k-mer once it has seen k good bases since the last restart.
///
/// `spacing` is `--kmer-spacing`: the counter is rewound to `kmer_size - spacing` after each
/// k-mer rather than to zero.
pub fn kmerize(bases: &[u8], kmer_size: usize, spacing: usize) -> Vec<Kmer> {
    let advance = kmer_size as i32 - spacing as i32;
    let mut out = Vec::new();
    let mut kmer = Kmer(0);
    let mut valid = 0i32;
    let mut index = 0usize;
    loop {
        let mut emitted = false;
        while index < bases.len() {
            match base_value(bases[index]) {
                Some(value) => kmer = kmer.successor(value, kmer_size),
                // A bad base restarts the count, and the increment below takes -1 back to zero.
                None => valid = -1,
            }
            index += 1;
            valid += 1;
            if valid == kmer_size as i32 {
                out.push(kmer);
                emitted = true;
                break;
            }
        }
        if !emitted {
            return out;
        }
        valid = advance;
    }
}

/// `PSKmerCollection.canonicalizeAndMask` over every k-mer of one sequence.
pub fn masked_kmers(
    bases: &[u8],
    kmer_size: usize,
    spacing: usize,
    mask: Kmer,
) -> Result<Vec<u64>, KmerError> {
    kmerize(bases, kmer_size, spacing)
        .into_iter()
        .map(|kmer| Ok(kmer.canonical(kmer_size)?.mask(mask).0))
        .collect()
}

/// The whole tool: the contigs in, the distinct longs out.
pub fn build(
    contigs: &[Vec<u8>],
    kmer_size: usize,
    spacing: usize,
    mask_argument: &str,
) -> Result<BTreeSet<u64>, KmerError> {
    let positions = parse_mask(mask_argument, kmer_size)?;
    let mask = get_mask(&positions, kmer_size);
    let mut set = BTreeSet::new();
    for bases in contigs {
        for value in masked_kmers(bases, kmer_size, spacing, mask)? {
            set.insert(value);
        }
    }
    if set.is_empty() {
        // The container refuses to be built from nothing, so the run never gets a file.
        return Err(KmerError::EmptySet);
    }
    Ok(set)
}
