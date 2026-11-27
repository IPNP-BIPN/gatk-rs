//! `SVKmerLong`, `SVKmerizer` and `SVDUSTFilteredKmerizer` (GATK 4.7.0.0; these and
//! `SVKmerShort`, `SVKmer` and `SVUtils`, all in `tools/spark/sv/utils`, are unchanged since
//! 4.6.2.0).
//!
//! # Two longs, split at the middle of the k-mer
//!
//! A k-mer of K bases is 2K bits held as two K-bit values, `valHigh` over `valLow`, so an odd K has
//! its middle base straddling them and `valHigh`'s low bit is the middle base's high bit: that is
//! the whole of `canonical`. `SVKmerLong(int)` refuses a K outside 1 to 63, but the no-argument
//! constructor the tools use checks nothing, and a Java shift counts modulo 64, so `--k-size 70`
//! runs and means whatever those shifts mean. `wrapping_shl` and `wrapping_shr` are the same
//! arithmetic, and `>>` stays signed.
//!
//! # A base that is not ACGT still rolls in, as an A
//!
//! The DUST kmerizer keeps its trimer counts honest by treating a bad base as an A and restarting
//! the count of valid bases, so a window can come back good with a score that remembers it. It
//! starts from the poly-A k-mer's score, `(K-2)(K-3)/2`, and returns a k-mer only at a valid count
//! of at least K with the score at or under the limit.
//!
//! # Reading this file without reading Rust
//!
//! - `i64` is Java's `long` and `i32` its `int`. Rust's `<<` and `>>` refuse a shift of 64 or
//!   more in a debug build, so the shifts that Java takes modulo 64 are `wrapping_shl` and
//!   `wrapping_shr`, which do the same thing Java does. On an `i64`, `>>` and `wrapping_shr` are
//!   Java's signed `>>`; Java's unsigned `>>>` is written as a cast to `u64` and a shift.
//! - `wrapping_add`, `wrapping_sub` and `wrapping_mul` overflow as Java's `int` and `long`
//!   arithmetic does, where a plain `+` would panic in a debug build.
//! - `as i32`, `as u32`, `as usize` reinterpret or truncate as a Java cast does: `(int)` of a
//!   `long` keeps its low 32 bits.
//! - `!x` on an integer is Java's `~x`.
//! - `impl Iterator for Kmerizer` makes the type usable in a `for` loop, as Java's
//!   `Iterator<SVKmer>`; `next` returns `None` where Java's `hasNext()` would be false.
//! - `Kmerizer<'a>` borrows the sequence it walks (`&'a [u8]`) rather than copying it.

/// `SVUtils.FNV64_DEFAULT_SEED`, `new BigInteger("14695981039346656037").longValue()`: the
/// unsigned constant's 64 bits read as a signed `long`, which the `u64 as i64` cast reproduces.
const FNV64_DEFAULT_SEED: i64 = 14695981039346656037u64 as i64;

/// `SVUtils.fnvLong64(start, toHash)`.
///
/// The Java unrolls the eight byte steps (`start ^= (toHash >> 56) & 0xffL; start *= mult;`
/// and so on down to `>> 0`); the loop over the shifts is the same sequence.
pub fn fnv_long_64(mut start: i64, to_hash: i64) -> i64 {
    const MULT: i64 = 1099511628211;
    for shift in [56, 48, 40, 32, 24, 16, 8, 0] {
        start ^= (to_hash >> shift) & 0xff;
        start = start.wrapping_mul(MULT);
    }
    start
}

/// `SVKmerLong`: a k-mer, two bits a base, split over two `long`s.
///
/// `Default` is the no-argument constructor, `valHigh = valLow = 0`, the poly-A k-mer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SvKmerLong {
    /// `valHigh`, the first half of the bases.
    pub high: i64,
    /// `valLow`, the second half.
    pub low: i64,
}

/// `(1L << kSize) - 1L`.
fn mask(k: i32) -> i64 {
    1i64.wrapping_shl(k as u32).wrapping_sub(1)
}

/// `SVKmer.reverseComplementByteValueAsLong`.
///
/// ```java
/// // this turns the 8 bits [b1 b2 b3 b4 b5 b6 b7 b8] into [~b7 ~b8 ~b5 ~b6 ~b3 ~b4 ~b1 ~b2]
/// ```
///
/// The Java argument is an `int`, so the arithmetic is done at `i32` and widened after the `~`.
/// The Java computes it once per byte value into `BYTEWISE_REVERSE_COMPLEMENT`; this computes it
/// at each call, with the same result.
fn reverse_complement_byte(value: i64) -> i64 {
    let value = value as i32;
    (!(((value & 3) << 6)
        | (((value >> 2) & 3) << 4)
        | (((value >> 4) & 3) << 2)
        | ((value >> 6) & 3)) as i64)
        & 0xff
}

/// `SVKmer.reverseComplement(long)`, a byte at a time.
///
/// Low byte first, then seven more, each shifted down with the signed `val >>= 8` and the
/// result shifted up: the Java's `while ( --nBytes != 0 )` runs seven times.
fn reverse_complement_long(mut value: i64) -> i64 {
    let mut result = reverse_complement_byte(value & 0xff);
    for _ in 0..7 {
        value >>= 8;
        result = (result << 8) | reverse_complement_byte(value & 0xff);
    }
    result
}

impl SvKmerLong {
    /// `successor`, with the base as its two bits.
    ///
    /// ```java
    /// final long newV1 = ((valHigh << 2) | (valLow >> (kSize-2))) & mask;
    /// ```
    ///
    /// `valLow >> (kSize-2)` moves the top base of the low half into the high half; the shift is
    /// signed, as the Java's is.
    pub fn successor(self, base: i64, k: i32) -> SvKmerLong {
        let mask = mask(k);
        SvKmerLong {
            high: (self.high.wrapping_shl(2) | self.low.wrapping_shr(k.wrapping_sub(2) as u32))
                & mask,
            low: (self.low.wrapping_shl(2) | (base & 3)) & mask,
        }
    }

    /// `reverseComplement(kSize)`. An even K swaps the halves, each shifted to the top and
    /// reversed. An odd K also moves the middle base's two bits, which straddle the halves: the
    /// Java XORs `((valHigh & 1L) ^ 1L)` into the new high half and
    /// `(valLow & (1L << (kSize-1)))` into the new low half.
    pub fn reverse_complement(self, k: i32) -> SvKmerLong {
        let mask = mask(k);
        let comp = 64i32.wrapping_sub(k);
        if k & 1 == 0 {
            return SvKmerLong {
                high: reverse_complement_long(self.low.wrapping_shl(comp as u32)) & mask,
                low: reverse_complement_long(self.high.wrapping_shl(comp as u32)) & mask,
            };
        }
        SvKmerLong {
            high: (reverse_complement_long(self.low.wrapping_shl(comp.wrapping_add(1) as u32))
                .wrapping_shl(1)
                & mask)
                ^ ((self.high & 1) ^ 1),
            low: (reverse_complement_long(self.high.wrapping_shl(comp as u32)) & mask)
                ^ (self.low & 1i64.wrapping_shl(k.wrapping_sub(1) as u32)),
        }
    }

    /// `canonical`, where an even K was refused before this is reached.
    pub fn canonical(self, k: i32) -> SvKmerLong {
        if self.high & 1 == 0 {
            self
        } else {
            self.reverse_complement(k)
        }
    }

    /// `firstTrimer`, whose shift is unsigned.
    ///
    /// `(int)(valHigh >>> (kSize-6))`: the top three bases. The `as i32` is the Java's `(int)`;
    /// the `as usize` then indexes the 64-entry count table.
    fn first_trimer(self, k: i32) -> usize {
        ((self.high as u64).wrapping_shr(k.wrapping_sub(6) as u32)) as i32 as usize
    }

    /// `lastTrimer`, `(int)valLow & 0x3F`: the last three bases.
    fn last_trimer(self) -> usize {
        (self.low as i32 & 0x3f) as usize
    }

    /// `hashCode`, `(int)SVUtils.fnvLong64(SVUtils.fnvLong64(valHigh), valLow)`; the one-argument
    /// `fnvLong64` starts from `FNV64_DEFAULT_SEED`.
    pub fn hash_code(self) -> i32 {
        fnv_long_64(fnv_long_64(FNV64_DEFAULT_SEED, self.high), self.low) as i32
    }

    /// `toString(kSize)`.
    ///
    /// The bases are read off the bottom two bits at a time, last base first, and the text is
    /// reversed at the end. For an odd K the middle base's high bit is the low bit of `high`,
    /// so `high` is shifted up one and ORed in before the second half is read.
    pub fn to_text(self, k: i32) -> String {
        const BASES: &[u8; 4] = b"ACGT";
        let mut text = Vec::with_capacity(k.max(0) as usize);
        let mut value = self.low;
        for _ in 0..k / 2 {
            text.push(BASES[(value & 3) as usize]);
            value >>= 2;
        }
        if k & 1 == 0 {
            value = self.high;
        } else {
            value |= self.high.wrapping_shl(1);
        }
        for _ in 0..(k + 1) / 2 {
            text.push(BASES[(value & 3) as usize]);
            value >>= 2;
        }
        text.reverse();
        String::from_utf8(text).expect("bases are ASCII")
    }
}

/// `SVKmerizer` over `SVKmerLong`: a bad base rolls nothing in and starts the count again. The
/// count after each k-mer is `kSize - kSpace`, so a `kSpace` of 0 gets past the test for K only
/// with a bad base's help: the adapter sequence gives one k-mer per run of good bases.
pub struct Kmerizer<'a> {
    seq: &'a [u8],
    k: i32,
    advance: i32,
    index: usize,
    next: Option<SvKmerLong>,
}

impl<'a> Kmerizer<'a> {
    /// `new SVKmerizer(seq, kSize, kSpace, new SVKmerLong())`: `kAdvance = kSize - kSpace`, and
    /// the first k-mer is found at once, as the Java constructor's `nextKmer(kmer, 0)` does.
    pub fn new(seq: &'a [u8], k: i32, k_space: i32) -> Kmerizer<'a> {
        let mut kmerizer = Kmerizer {
            seq,
            k,
            advance: k.wrapping_sub(k_space),
            index: 0,
            next: None,
        };
        kmerizer.next = kmerizer.step(SvKmerLong::default(), 0);
        kmerizer
    }

    /// `nextKmer(tmpKmer, validBaseCount)`: roll bases in until `validBaseCount` reaches K, or
    /// return `None` (the Java's `null`) when the sequence ends first. A bad base sets the
    /// count to -1, so the `++` that follows makes it 0.
    fn step(&mut self, mut kmer: SvKmerLong, mut valid: i32) -> Option<SvKmerLong> {
        while self.index < self.seq.len() {
            match self.seq[self.index] {
                b'a' | b'A' => kmer = kmer.successor(0, self.k),
                b'c' | b'C' => kmer = kmer.successor(1, self.k),
                b'g' | b'G' => kmer = kmer.successor(2, self.k),
                b't' | b'T' => kmer = kmer.successor(3, self.k),
                _ => valid = -1,
            }
            self.index += 1;
            valid = valid.wrapping_add(1);
            if valid == self.k {
                return Some(kmer);
            }
        }
        None
    }
}

impl Iterator for Kmerizer<'_> {
    type Item = SvKmerLong;

    /// `next()`: hand out the k-mer already found and look for the one after it, starting the
    /// count at `kAdvance`.
    fn next(&mut self) -> Option<SvKmerLong> {
        let result = self.next?;
        self.next = self.step(result, self.advance);
        Some(result)
    }
}

/// `SVKmerShort`, for the K of 31 `FermiLiteAssemblyHandler` looks for shadowed contigs with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SvKmerShort(pub i64);

impl SvKmerShort {
    /// `(1L << kSize * 2) - 1L`. The plain `<<` is safe for the K of 31 this type is used with.
    const fn mask(k: i32) -> i64 {
        (1i64 << (k * 2)) - 1
    }

    /// `successor`: shift up one base, OR the new one in, mask to K bases.
    pub fn successor(self, base: i64, k: i32) -> SvKmerShort {
        SvKmerShort(((self.0 << 2) | (base & 3)) & Self::mask(k))
    }

    /// `reverseComplement`: the bases moved to the top (`valLow << compK`), reversed, masked.
    pub fn reverse_complement(self, k: i32) -> SvKmerShort {
        SvKmerShort(reverse_complement_long(self.0 << (64 - k * 2)) & Self::mask(k))
    }

    /// `canonical`: `((valLow >> kSize) & 1L) == 0` is a middle base of A or C, which is kept.
    pub fn canonical(self, k: i32) -> SvKmerShort {
        if (self.0 >> k) & 1 == 0 {
            self
        } else {
            self.reverse_complement(k)
        }
    }

    /// `hashCode`, `SVUtils.fnvLong64(2166136261L, valLow)` cut to an int.
    pub fn hash_code(self) -> i32 {
        fnv_long_64(2166136261, self.0) as i32
    }

    /// `new SVKmerizer(seq, k, new SVKmerShort())` run out: the same rule as [`Kmerizer`].
    pub fn kmers(seq: &[u8], k: i32) -> Vec<SvKmerShort> {
        let mut kmers = Vec::new();
        let mut kmer = SvKmerShort(0);
        let mut valid = 0;
        for base in seq {
            match base {
                b'a' | b'A' => kmer = kmer.successor(0, k),
                b'c' | b'C' => kmer = kmer.successor(1, k),
                b'g' | b'G' => kmer = kmer.successor(2, k),
                b't' | b'T' => kmer = kmer.successor(3, k),
                _ => valid = -1,
            }
            valid += 1;
            if valid == k {
                kmers.push(kmer);
                // The 3-argument `SVKmerizer` constructor has `kSpace` 1, so the count restarts
                // at `kAdvance = kSize - 1`.
                valid = k - 1;
            }
        }
        kmers
    }
}

/// `SVDUSTFilteredKmerizer.canonicalStream`'s source, before `canonical`.
pub struct DustFilteredKmerizer<'a> {
    seq: &'a [u8],
    k: i32,
    /// `maxDUSTScore`.
    max_dust: i32,
    /// `trimerCounts`: how many times each of the 64 trimers occurs in the current window.
    trimer_counts: Vec<i32>,
    /// `curDUSTScore`, the sum over trimers of `count*(count-1)/2`.
    score: i32,
    index: usize,
    next: Option<SvKmerLong>,
}

impl<'a> DustFilteredKmerizer<'a> {
    /// The caller has refused a K under 4, which this constructor throws on.
    pub fn new(seq: &'a [u8], k: i32, max_dust: i32) -> DustFilteredKmerizer<'a> {
        let poly_a = k.wrapping_sub(2);
        let mut trimer_counts = vec![0; 64];
        trimer_counts[0] = poly_a;
        let mut kmerizer = DustFilteredKmerizer {
            seq,
            k,
            max_dust,
            trimer_counts,
            score: poly_a.wrapping_mul(poly_a.wrapping_sub(1)) / 2,
            index: 0,
            next: None,
        };
        kmerizer.next = kmerizer.advance(SvKmerLong::default(), 0);
        kmerizer
    }

    /// `nextKmer(initialKmer, validBaseCount)`: per base, the leading trimer's count is taken
    /// out of the score, the base is rolled in, and the new trailing trimer's count is added.
    ///
    /// ```java
    /// curDUSTScore -= --trimerCounts[result.firstTrimer(kSize)];
    /// ```
    ///
    /// The Java pre-decrements inside the expression; the two lines below split it into the
    /// decrement and then the subtraction of the decremented count. The `+=` with a
    /// post-increment after the roll is split the other way round.
    fn advance(&mut self, initial: SvKmerLong, mut valid: i32) -> Option<SvKmerLong> {
        let mut result = initial;
        while self.index < self.seq.len() {
            let leaving = result.first_trimer(self.k);
            self.trimer_counts[leaving] = self.trimer_counts[leaving].wrapping_sub(1);
            self.score = self.score.wrapping_sub(self.trimer_counts[leaving]);
            match self.seq[self.index] {
                b'a' | b'A' => result = result.successor(0, self.k),
                b'c' | b'C' => result = result.successor(1, self.k),
                b'g' | b'G' => result = result.successor(2, self.k),
                b't' | b'T' => result = result.successor(3, self.k),
                _ => {
                    result = result.successor(0, self.k);
                    valid = -1;
                }
            }
            let entering = result.last_trimer();
            self.score = self.score.wrapping_add(self.trimer_counts[entering]);
            self.trimer_counts[entering] = self.trimer_counts[entering].wrapping_add(1);
            self.index += 1;
            valid = valid.wrapping_add(1);
            if valid >= self.k && self.score <= self.max_dust {
                return Some(result);
            }
        }
        None
    }
}

impl Iterator for DustFilteredKmerizer<'_> {
    type Item = SvKmerLong;

    fn next(&mut self) -> Option<SvKmerLong> {
        let result = self.next?;
        // `kAdvance` is K less the one base a step moves.
        self.next = self.advance(result, self.k.wrapping_sub(1));
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kmer(text: &str) -> SvKmerLong {
        let k = text.len() as i32;
        text.bytes().fold(SvKmerLong::default(), |kmer, base| {
            kmer.successor(b"ACGT".iter().position(|b| *b == base).unwrap() as i64, k)
        })
    }

    #[test]
    fn text_round_trips_and_the_middle_base_chooses_the_strand() {
        assert_eq!(kmer("ACGTA").to_text(5), "ACGTA");
        // Middle base A: kept. Middle base G: the reverse complement.
        assert_eq!(kmer("TTACC").canonical(5).to_text(5), "TTACC");
        assert_eq!(kmer("TTGCC").canonical(5).to_text(5), "GGCAA");
        assert_eq!(kmer("ACGT").reverse_complement(4).to_text(4), "ACGT");
    }

    #[test]
    fn a_low_complexity_window_is_withheld() {
        let polya = vec![b'A'; 20];
        assert_eq!(DustFilteredKmerizer::new(&polya, 7, 5).count(), 0);
        assert_eq!(DustFilteredKmerizer::new(&polya, 7, 10).count(), 14);
    }

    #[test]
    fn a_bad_base_costs_a_whole_window() {
        let seq = b"ACGTTGCANACGTTGCA";
        let kmers: Vec<String> = DustFilteredKmerizer::new(seq, 8, 100)
            .map(|kmer| kmer.to_text(8))
            .collect();
        assert_eq!(kmers, ["ACGTTGCA", "ACGTTGCA"]);
    }
}
