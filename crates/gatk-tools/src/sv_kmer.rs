//! `SVKmerLong`, `SVKmerizer` and `SVDUSTFilteredKmerizer` (GATK 4.6.2.0).
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

/// `SVUtils.FNV64_DEFAULT_SEED`.
const FNV64_DEFAULT_SEED: i64 = 14695981039346656037u64 as i64;

/// `SVUtils.fnvLong64(start, toHash)`.
pub fn fnv_long_64(mut start: i64, to_hash: i64) -> i64 {
    const MULT: i64 = 1099511628211;
    for shift in [56, 48, 40, 32, 24, 16, 8, 0] {
        start ^= (to_hash >> shift) & 0xff;
        start = start.wrapping_mul(MULT);
    }
    start
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SvKmerLong {
    pub high: i64,
    pub low: i64,
}

/// `(1L << kSize) - 1L`.
fn mask(k: i32) -> i64 {
    1i64.wrapping_shl(k as u32).wrapping_sub(1)
}

/// `SVKmer.reverseComplementByteValueAsLong`.
fn reverse_complement_byte(value: i64) -> i64 {
    let value = value as i32;
    (!(((value & 3) << 6)
        | (((value >> 2) & 3) << 4)
        | (((value >> 4) & 3) << 2)
        | ((value >> 6) & 3)) as i64)
        & 0xff
}

/// `SVKmer.reverseComplement(long)`, a byte at a time.
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
    pub fn successor(self, base: i64, k: i32) -> SvKmerLong {
        let mask = mask(k);
        SvKmerLong {
            high: (self.high.wrapping_shl(2) | self.low.wrapping_shr(k.wrapping_sub(2) as u32))
                & mask,
            low: (self.low.wrapping_shl(2) | (base & 3)) & mask,
        }
    }

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
    fn first_trimer(self, k: i32) -> usize {
        ((self.high as u64).wrapping_shr(k.wrapping_sub(6) as u32)) as i32 as usize
    }

    fn last_trimer(self) -> usize {
        (self.low as i32 & 0x3f) as usize
    }

    /// `hashCode`.
    pub fn hash_code(self) -> i32 {
        fnv_long_64(fnv_long_64(FNV64_DEFAULT_SEED, self.high), self.low) as i32
    }

    /// `toString(kSize)`.
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

/// `SVDUSTFilteredKmerizer.canonicalStream`'s source, before `canonical`.
pub struct DustFilteredKmerizer<'a> {
    seq: &'a [u8],
    k: i32,
    max_dust: i32,
    trimer_counts: Vec<i32>,
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
