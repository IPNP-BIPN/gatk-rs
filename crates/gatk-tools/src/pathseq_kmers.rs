//! `PathSeqBuildKmers`, ported from the tool, `PSKmerUtils`, `PSKmerCollection`, `SVKmerShort` and
//! `SVKmerizer` (GATK 4.7.0.0; every one of those Java files is unchanged since 4.6.2.0).
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
//!
//! # Reading this file without reading Rust
//!
//! - `u64` is an unsigned 64-bit integer. Java's k-mer is a signed `long`; the port holds the same
//!   64 bits unsigned, which changes only what `>>` does: on a `u64` it shifts zeros in, where
//!   Java's `>>` on a `long` copies the sign bit. Every shift below is followed by a mask that
//!   discards the bits where the two differ, which the item comments point out where it matters.
//! - `Option<u64>` is "a value or nothing": `Some(v)` or `None`, and the compiler makes the caller
//!   handle both. `Result<T, KmerError>` is "a value or one of the refusals below"; `Err(..)` is
//!   the Java throw and `Ok(..)` the return.
//! - `?` after a call that returns a `Result` hands its `Err` straight to this function's caller,
//!   like an exception propagating out of a Java method, and otherwise unwraps the `Ok`.
//! - `x as i8`, `x as u64` and the like are casts. Between integer types they never fail: they
//!   truncate or reinterpret the bits, exactly as Java's `(byte)` and `(long)` casts do.
//! - `&[u8]` is a borrowed, read-only view of a `byte[]`; the caller keeps ownership.
//! - A function's last expression, written without `;`, is its return value.

use std::collections::BTreeSet;

/// `SVKmer.Base`: A=0, C=1, G=2, T=3.
///
/// **What**: the two-bit code of one base, or `None` for anything that is not ACGT in either
/// case. It is the `switch` of `SVKmerizer.nextKmer`, whose `default` arm is the `None` here:
///
/// ```java
/// case 'a': case 'A': tmpKmer = tmpKmer.successor(SVKmer.Base.A, kSize); break;
/// ```
fn base_value(base: u8) -> Option<u64> {
    match base {
        b'a' | b'A' => Some(0),
        b'c' | b'C' => Some(1),
        b'g' | b'G' => Some(2),
        b't' | b'T' => Some(3),
        _ => None,
    }
}

/// `SVKmer.reverseComplementByteValueAsLong`: the reverse complement of the four bases packed in
/// one byte.
///
/// ```java
/// return ~(((bIn & 3) << 6) | (((bIn >> 2) & 3) << 4) | (((bIn >> 4) & 3) << 2)
///         | ((bIn >> 6) & 3)) & 0xffL;
/// ```
///
/// The four two-bit pairs are reversed, then every bit is inverted, which complements each base
/// because A=0 and T=3, C=1 and G=2 are bitwise complements. Java computes the 256 answers once
/// into a table; the port computes each one on demand, which gives the same bytes. `!` is Rust's
/// bitwise NOT, Java's `~`.
fn reverse_complement_byte(byte: u8) -> u64 {
    let value = byte as u64;
    (!(((value & 3) << 6)
        | (((value >> 2) & 3) << 4)
        | (((value >> 4) & 3) << 2)
        | ((value >> 6) & 3)))
        & 0xff
}

/// `SVKmer.reverseComplement(long)`, byte by byte.
///
/// **How**: the low byte's reverse complement becomes the high byte of the result, and so on for
/// all eight bytes, which reverses the 32 bases of the whole long.
///
/// **Why the unsigned shift is safe**: Java writes `val >>= 8;` on a signed `long`, which fills
/// the top with copies of the sign bit, where this `u64` fills it with zeros. Only the low byte
/// is ever read after a shift, and seven shifts of eight bits never bring a filled bit down into
/// it, so the two agree on every input.
pub fn reverse_complement_long(value: u64) -> u64 {
    let mut result = reverse_complement_byte((value & 0xff) as u8);
    let mut rest = value;
    for _ in 0..7 {
        rest >>= 8;
        result = (result << 8) | reverse_complement_byte((rest & 0xff) as u8);
    }
    result
}

/// What the run refuses, each variant carrying what its Java message quotes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KmerError {
    /// `canonical`, which the tool reaches only once it is already streaming.
    EvenKmerSize,
    /// `PSUtils.parseMask`, over an index that parsed.
    InvalidMaskIndex { index: String },
    /// `Integer.parseInt`, over a field of the mask that is not an int.
    NotAnInt { field: String },
    /// The container, when the reference produced no k-mer at all.
    EmptySet,
}

impl KmerError {
    /// The exception class the reference throws. `Integer.parseInt` throws a
    /// `NumberFormatException`; `Utils.validateArg` throws an `IllegalArgumentException`, and every
    /// other refusal here is one.
    pub fn java_class(&self) -> &'static str {
        match self {
            KmerError::NotAnInt { .. } => "java.lang.NumberFormatException",
            _ => "java.lang.IllegalArgumentException",
        }
    }

    /// The exception's message, as the reference builds it.
    pub fn message(&self) -> String {
        match self {
            KmerError::EvenKmerSize => "Kmer length must be odd to canonicalize.".to_string(),
            KmerError::InvalidMaskIndex { index } => format!("Invalid kmer mask index: {index}"),
            KmerError::NotAnInt { field } => format!("For input string: \"{field}\""),
            KmerError::EmptySet => "Number of elements must be greater than 0".to_string(),
        }
    }
}

/// `SVKmerShort`: up to 31 bases in the low bits of one long, first base highest.
///
/// The struct is a single `u64` with no size attached, as Java's is a single `long`: like
/// `SVKmerShort`, every method that needs k is told it. The derived `Ord` compares the `u64`,
/// which agrees with Java's `Long.compare` in `compareTo` for a k-mer of at most 31 bases, whose
/// top bits are clear; a mask from `get_mask` has them set and would order differently.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer(pub u64);

impl Kmer {
    /// `successor`, which drops the leading base and adds one at the end.
    ///
    /// The window mask is k two-bit pairs of ones, Java's `(1L << kSize * 2) - 1L`.
    pub fn successor(&self, base: u64, kmer_size: usize) -> Kmer {
        let window = (1u64 << (kmer_size * 2)) - 1;
        Kmer(((self.0 << 2) | (base & 3)) & window)
    }

    /// `reverseComplement`: the k-mer is moved up to the top of the long, the whole long is
    /// reversed and complemented, and the result masked back to k bases.
    ///
    /// ```java
    /// final int compK = 64 - kSize * 2;
    /// // move the significant bits up to the top, reverse complement, and mask to kSize bits.
    /// return new SVKmerShort(reverseComplement(valLow << compK) & mask);
    /// ```
    pub fn reverse_complement(&self, kmer_size: usize) -> Kmer {
        let window = (1u64 << (kmer_size * 2)) - 1;
        let unused = 64 - kmer_size * 2;
        Kmer(reverse_complement_long(self.0 << unused) & window)
    }

    /// `canonical`, which is defined only for an odd size.
    ///
    /// The size is checked first, and then the one bit is tested (the Java comment line between
    /// the two is left out):
    ///
    /// ```java
    /// Utils.validateArg( (kSize & 1) != 0, "Kmer length must be odd to canonicalize.");
    /// if (((valLow >> kSize) & 1L) == 0) return this;
    /// ```
    ///
    /// The middle base of an odd k occupies bits `k - 1` and `k`, so bit `k` is its high bit: 0 for
    /// A and C, 1 for G and T. Because the size check comes first, an even size is refused even
    /// for a k-mer that would not need turning round.
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

    /// `mask`, which is a plain AND with the value `get_mask` built.
    pub fn mask(&self, mask: Kmer) -> Kmer {
        Kmer(self.0 & mask.0)
    }

    /// `toString(kSize)`, which spells the bases back out.
    ///
    /// Like the Java, it reads the last base first, two bits at a time, and reverses at the end.
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
///
/// ```java
/// mask |= 3L << 2*(kSize - pos - 1);
/// ```
///
/// The result is the complement of the masked pairs, so every bit above the k-mer is also set;
/// a k-mer has nothing there, so the AND in [`Kmer::mask`] leaves those bits zero either way.
/// `parse_mask` has already refused any position outside `0..k`, so the shift cannot underflow.
pub fn get_mask(positions: &[u8], kmer_size: usize) -> Kmer {
    let mut mask: u64 = 0;
    for position in positions {
        mask |= 3u64 << (2 * (kmer_size - *position as usize - 1));
    }
    Kmer(!mask)
}

/// `PSUtils.parseMask`, whose empty string is no mask at all.
///
/// Three things are Java's rather than the tool's: `String.split` drops TRAILING empty fields, so
/// `"5,"` is one index and `","` is none; a field `Integer.parseInt` refuses is a
/// `NumberFormatException` rather than the tool's own refusal; and the int is cast to a BYTE before
/// it is checked, so `257` is index 1 while the message still quotes the field as written.
///
/// ```java
/// kmerMask[i] = (byte) Integer.parseInt(kmerMaskString[i]);
/// Utils.validateArg(kmerMask[i] >= 0 && kmerMask[i] < kSize,
///         "Invalid kmer mask index: " + kmerMaskString[i]);
/// ```
///
/// The split is `PSUtils.parseCommaDelimitedArgList`, which returns no fields for an empty
/// argument and otherwise `return arg.split(",");`.
pub fn parse_mask(argument: &str, kmer_size: usize) -> Result<Vec<u8>, KmerError> {
    // Rust's `split` keeps every empty field, including the one an empty string produces, so the
    // trailing empties Java drops are popped here. An empty argument becomes no fields at all,
    // which is also what `parseCommaDelimitedArgList` returns for it.
    let mut fields: Vec<&str> = argument.split(',').collect();
    while fields.last() == Some(&"") {
        fields.pop();
    }
    let mut positions = Vec::new();
    for field in fields {
        // `parse::<i32>` is `Integer.parseInt`; `map_err` turns its failure into the
        // `NumberFormatException` above, `?` returns it, and `as i8` is the `(byte)` cast.
        let index = field.parse::<i32>().map_err(|_| KmerError::NotAnInt {
            field: field.to_string(),
        })? as i8;
        if index < 0 || index as i32 >= kmer_size as i32 {
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
///
/// **How**: the outer `loop` is one call of `nextKmer` per k-mer, the inner `while` its scan.
/// The scan resumes at `index`, where the previous k-mer ended, with the k-mer bits carried over,
/// as Java's iterator keeps `idx` and passes the last k-mer back in:
///
/// ```java
/// this.kAdvance = kSize - kSpace;
/// ```
///
/// A spacing larger than k makes the restart negative, which only means more bases must be read
/// before the next k-mer, as it does in Java.
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
///
/// ```java
/// return val.canonical(kSize).mask(mask).getLong();
/// ```
///
/// The first even-size `canonical` refusal stops the whole sequence: collecting an iterator of
/// `Result`s into one `Result` returns the first `Err` it meets, as the Java stream throws at the
/// first k-mer.
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
///
/// `PathSeqBuildKmers.doWork`: parse the mask, build it, k-merize each contig. A `BTreeSet` is a
/// sorted set, Java's `TreeSet`; the port keeps the distinct values sorted, which is a property of
/// the port's container and not of the hopscotch set or Bloom filter the tool writes.
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
        // The container refuses to be built from nothing, so the run never gets a file. Both
        // `LargeLongHopscotchSet` and `LongBloomFilter` carry the same check:
        // `Utils.validateArg(numElements > 0, "Number of elements must be greater than 0");`
        return Err(KmerError::EmptySet);
    }
    Ok(set)
}
