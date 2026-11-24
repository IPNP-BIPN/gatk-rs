//! Ported from `org.apache.commons.math3.random.Well19937c`, `org.apache.commons.math3.random.AbstractWell`
//! and `org.apache.commons.math3.random.BitsStreamGenerator` (commons-math3 3.5, Apache 2.0,
//! verified from the licence header of the sources jar in the resolved dependency).
//!
//! GATK carries **two** static generators, not one:
//!
//! ```text
//! private static final Random randomGenerator = new Random(GATK_RANDOM_SEED);
//! private static final RandomDataGenerator randomDataGenerator =
//!         new RandomDataGenerator(new Well19937c(GATK_RANDOM_SEED));
//! ```
//!
//! They are different algorithms with different streams, and `Utils.resetRandomGenerator()` resets
//! both. [`crate::java_random`] is the first; this is the second. A port that routed both through
//! one generator would agree with the reference on neither.
//!
//! # Where it differs from `java.util.Random`, method by method
//!
//! The two share method *names* and almost nothing else. `BitsStreamGenerator` derives its methods
//! from `next(bits)` differently from `java.util.Random`, so even a generator whose `next` were
//! somehow identical would still disagree:
//!
//! | method | `java.util.Random` | `BitsStreamGenerator` |
//! |---|---|---|
//! | `nextDouble` | `((next(26) << 27) + next(27)) * 2^-53` | `((next(26) << 26) \| next(26)) * 2^-52` |
//! | `nextFloat` | `next(24) / 2^24` | `next(23) * 2^-23` |
//! | `nextLong` | `(next(32) << 32) + next(32)`, signed add | `(next(32) << 32) \| (next(32) & 0xffffffff)` |
//!
//! `nextInt(bound)` is the one that does agree: commons-math copied it from Apache Harmony's
//! `java.util.Random`, so the power-of-two path and the rejection loop are the same shape.
//!
//! # What is not here
//!
//! `nextGaussian` is deliberately absent. It is `r * FastMath.cos(alpha)` over commons-math's own
//! `FastMath`, whose `cos`, `sin`, `log` and `sqrt` are its implementations rather than the
//! platform's, so a port of the *shape* would be bit-wrong in the mantissa while looking right.
//! It belongs with jmath, and refusing is better than a plausible answer. Nothing GATK reaches
//! through `RandomDataGenerator.nextPermutation` needs it.
//!
//! # Reading this file without reading Rust
//!
//! - `i32` and `i64` are Java's `int` and `long`. `u32`, `u8` and `usize` are unsigned types Java
//!   lacks: a `u8` holds the same eight bits as a Java `byte`, and `usize` is the index type.
//! - `x as T` converts: to a narrower integer it keeps the low bits (Java's `(int)` or `(byte)`),
//!   to a wider signed one it sign-extends (Java's promotion of `int` to `long`), and between
//!   signed and unsigned of one width it reinterprets the same bits.
//! - Rust's `>>` on a signed integer is Java's arithmetic `>>`. Java's `>>>` has no operator here
//!   and is the private `ushr` at the end of the file.
//! - `wrapping_mul`, `wrapping_add` and `wrapping_sub` are Java's arithmetic, which wraps silently
//!   on overflow. Plain `*`, `+` and `-` would panic in a debug build instead.
//! - `Vec<i32>` is an `int[]` the struct owns. `&[i32]` and `&mut [u8]` are read-only and writable
//!   borrowed views of an array the caller keeps, and `v[a..b]` is the part from `a` up to but not
//!   including `b`.
//! - `&mut self` on a method means it advances this generator in place, not a copy of it.

/// Number of bits in the pool.
const K: usize = 19937;
/// `r`, the number of 32-bit blocks the pool is made of: `(K + 31) / 32`.
const R: usize = K.div_ceil(32);

// The recurrence's three tap offsets into the pool, read in `next` as `v[(index + M1) % R]` and
// so on.
const M1: usize = 70;
const M2: usize = 179;
const M3: usize = 449;

/// `Well19937c`.
///
/// The reference precomputes five index tables (`iRm1`, `iRm2`, `i1`, `i2`, `i3`) so the hot loop
/// avoids a modulo. That is an optimisation over `(j + r - 1) % r` and friends, and this computes
/// the same indices directly: same values, no state to keep in step.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Well19937c {
    /// The pool: `R` ints, used as a ring.
    v: Vec<i32>,
    /// Where the next draw's recurrence starts; it steps backwards by one per draw.
    index: usize,
}

impl Well19937c {
    /// `new Well19937c(long seed)`, which is the constructor GATK uses.
    ///
    /// The long is split into **two** ints, high word first, so it is not the same seeding as
    /// `new Well19937c((int) seed)` even when the long fits in an int: the array length changes
    /// how much of the pool the scrambler fills.
    pub fn from_long(seed: i64) -> Self {
        // `seed >> 32` is the high word (an arithmetic shift, Java's `>>`) and
        // `seed & 0xffff_ffff` the low one; `as i32` keeps each one's low 32 bits.
        Well19937c::from_int_array(&[(seed >> 32) as i32, (seed & 0xffff_ffff) as i32])
    }

    /// `new Well19937c(int seed)`.
    pub fn from_int(seed: i32) -> Self {
        Well19937c::from_int_array(&[seed])
    }

    /// `new Well19937c(int[] seed)`.
    pub fn from_int_array(seed: &[i32]) -> Self {
        let mut generator = Well19937c {
            v: vec![0; R],
            index: 0,
        };
        generator.set_seed_ints(seed);
        generator
    }

    /// `Utils.getRandomDataGenerator()`'s generator: `new Well19937c(47382911L)`.
    ///
    /// The static field is initialised with
    /// `new RandomDataGenerator(new Well19937c(GATK_RANDOM_SEED))`, and:
    ///
    /// ```java
    /// public static RandomDataGenerator getRandomDataGenerator() { return randomDataGenerator; }
    /// ```
    ///
    /// GATK reaches it through `MathUtils.sampleIndicesWithoutReplacement`
    /// (`return Utils.getRandomDataGenerator().nextPermutation(n, k);`, ported in
    /// [`crate::permutation`]). As with [`crate::java_random::JavaRandom::gatk`], the Java hands
    /// out one static instance and this returns a fresh one, so the port creates it once and
    /// passes the same `&mut Well19937c` to every consumer.
    pub fn gatk() -> Self {
        Well19937c::from_long(crate::java_random::GATK_RANDOM_SEED)
    }

    /// `setSeed(int[] seed)`.
    ///
    /// The seed is copied into the head of the pool and the tail is filled by a Knuth-style
    /// scrambler. The widening is the part worth being careful about: `v[i - seed.length]` is an
    /// `int` promoted to `long`, so a negative entry sign-extends and `l >> 30` is an *arithmetic*
    /// shift over 64 bits, not the logical shift over 32 that reading it as unsigned would give.
    pub fn set_seed_ints(&mut self, seed: &[i32]) {
        // Copy the first `taken` seed entries into the pool; a seed longer than the pool is cut.
        let taken = seed.len().min(R);
        self.v[..taken].copy_from_slice(&seed[..taken]);
        if seed.len() < R {
            // `a..b` is the half-open range a, a+1, ... b-1: `for (int i = a; i < b; i++)`.
            for i in seed.len()..R {
                let l = self.v[i - seed.len()] as i64;
                // `& 0xffff_ffff` then `as i32` keep the low 32 bits of the 64-bit result.
                self.v[i] = (1_812_433_253i64
                    .wrapping_mul(l ^ (l >> 30))
                    .wrapping_add(i as i64)
                    & 0xffff_ffff) as i32;
            }
        }
        self.index = 0;
    }

    /// `setSeed(long seed)`, which is what `RandomDataGenerator.reSeed(long)` calls and therefore
    /// what `Utils.resetRandomGenerator()` reaches:
    /// `randomDataGenerator.reSeed(GATK_RANDOM_SEED);`.
    pub fn set_seed_long(&mut self, seed: i64) {
        self.set_seed_ints(&[(seed >> 32) as i32, (seed & 0xffff_ffff) as i32]);
    }

    /// `protected int next(int bits)`: the WELL19937c recurrence, with the Matsumoto-Kurita
    /// tempering that makes it the "c" variant rather than `Well19937a`.
    pub fn next(&mut self, bits: u32) -> i32 {
        // `% R` wraps an index past the end of the ring; these are the values the reference's
        // precomputed tables hold (see the struct note).
        let index_rm1 = (self.index + R - 1) % R;
        let index_rm2 = (self.index + R - 2) % R;

        let v0 = self.v[self.index];
        let vm1 = self.v[(self.index + M1) % R];
        let vm2 = self.v[(self.index + M2) % R];
        let vm3 = self.v[(self.index + M3) % R];

        // `0x8000_0000u32 as i32` is Java's `int` literal `0x80000000`, which is negative. Rust
        // refuses an out-of-range `i32` literal, so it is written unsigned and reinterpreted:
        // the same bits.
        let z0 = (0x8000_0000u32 as i32 & self.v[index_rm1]) ^ (0x7fff_ffff & self.v[index_rm2]);
        let z1 = (v0 ^ (v0 << 25)) ^ (vm1 ^ ushr(vm1, 27));
        let z2 = ushr(vm2, 9) ^ (vm3 ^ ushr(vm3, 1));
        let z3 = z1 ^ z2;
        let mut z4 = z0 ^ (z1 ^ (z1 << 9)) ^ (z2 ^ (z2 << 21)) ^ (z3 ^ ushr(z3, 21));

        // Store the new state; the index then steps backwards, so the ring is walked downwards.
        self.v[self.index] = z3;
        self.v[index_rm1] = z4;
        self.v[index_rm2] &= 0x8000_0000u32 as i32;
        self.index = index_rm1;

        // The tempering. `z4` is a copy of what was just stored, so it alters the output only,
        // never the pool.
        z4 ^= (z4 << 7) & 0xe46e_1700u32 as i32;
        z4 ^= (z4 << 15) & 0x9b86_8000u32 as i32;

        // The top `bits` bits, by an unsigned shift.
        ushr(z4, 32 - bits)
    }

    /// `nextBoolean()`.
    pub fn next_boolean(&mut self) -> bool {
        self.next(1) != 0
    }

    /// `nextInt()`.
    pub fn next_int(&mut self) -> i32 {
        self.next(32)
    }

    /// `nextInt(int n)`, copied upstream from Apache Harmony's `java.util.Random`, so this one *is*
    /// the same shape as [`crate::java_random::JavaRandom::next_int_bound`]: a power of two takes
    /// the high bits of a single draw, anything else loops until the value is unbiased and so can
    /// consume more than one.
    ///
    /// The power-of-two product is a plain `*` where `JavaRandom` writes `wrapping_mul`: both
    /// factors are below 2^31, so it cannot overflow a long either way.
    pub fn next_int_bound(&mut self, n: i32) -> i32 {
        // `assert!` panics (stops the program) when false, standing in for the reference's
        // rejection of a non-positive `n`. `{n}` puts the value into the message.
        assert!(n > 0, "{n} is smaller than, or equal to, the minimum (0)");
        if (n & -n) == n {
            return (((n as i64) * (self.next(31) as i64)) >> 31) as i32;
        }
        loop {
            let bits = self.next(31);
            let val = bits % n;
            if bits.wrapping_sub(val).wrapping_add(n - 1) >= 0 {
                return val;
            }
        }
    }

    /// `nextLong()`: two draws or-ed, the low one masked. Not the signed add `java.util.Random`
    /// does.
    pub fn next_long(&mut self) -> i64 {
        let high = (self.next(32) as i64) << 32;
        let low = (self.next(32) as i64) & 0xffff_ffff;
        high | low
    }

    /// `nextLong(long n)`.
    ///
    /// The same rejection loop as `next_int_bound`, over 63 bits: the first draw is 31 bits, so
    /// `bits` is never negative.
    pub fn next_long_bound(&mut self, n: i64) -> i64 {
        assert!(n > 0, "{n} is smaller than, or equal to, the minimum (0)");
        loop {
            let mut bits = (self.next(31) as i64) << 32;
            bits |= (self.next(32) as i64) & 0xffff_ffff;
            let val = bits % n;
            if bits.wrapping_sub(val).wrapping_add(n - 1) >= 0 {
                return val;
            }
        }
    }

    /// `nextDouble()`: 26 bits shifted by 26 and or-ed with 26 more, scaled by 2^-52.
    pub fn next_double(&mut self) -> f64 {
        let high = (self.next(26) as i64) << 26;
        let low = self.next(26) as i64;
        ((high | low) as f64) * f64::from_bits(0x3cb0_0000_0000_0000) // 0x1.0p-52
    }

    /// `nextFloat()`: 23 bits scaled by 2^-23.
    pub fn next_float(&mut self) -> f32 {
        (self.next(23) as f32) * f32::from_bits(0x3400_0000) // 0x1.0p-23f
    }

    /// `nextBytes(byte[])`.
    ///
    /// The tail is the part that is not obvious: after the four-at-a-time loop the reference takes
    /// **one more full draw** even when nothing is left to fill, so the number of draws depends on
    /// the length in a way a naive port gets wrong for lengths that are a multiple of four.
    pub fn next_bytes(&mut self, bytes: &mut [u8]) {
        let mut i = 0usize;
        // `saturating_sub` is subtraction that stops at zero: the length is unsigned, so a plain
        // `- 3` on a slice shorter than three would panic. Either way the loop below runs only
        // while at least four bytes remain.
        let end = bytes.len().saturating_sub(3);
        while i < end {
            let random = self.next(32);
            bytes[i] = (random & 0xff) as u8;
            bytes[i + 1] = (ushr(random, 8) & 0xff) as u8;
            bytes[i + 2] = (ushr(random, 16) & 0xff) as u8;
            bytes[i + 3] = (ushr(random, 24) & 0xff) as u8;
            i += 4;
        }
        // Taken unconditionally, even when `i == bytes.len()`.
        let mut random = self.next(32);
        // At most three bytes remain, so `random` is shifted at most twice before its last use
        // and the sign bits a `>>` copies in never reach the low byte.
        while i < bytes.len() {
            bytes[i] = (random & 0xff) as u8;
            i += 1;
            random >>= 8;
        }
    }
}

/// Java's `>>>` over an `i32`.
///
/// Reinterpreting as unsigned first makes Rust's `>>` the logical shift, then the result is
/// reinterpreted back.
fn ushr(value: i32, bits: u32) -> i32 {
    ((value as u32) >> bits) as i32
}
