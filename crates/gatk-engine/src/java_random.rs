//! `java.util.Random`, implemented from its **specification** rather than from its implementation.
//!
//! This is the opposite case to [`crate::java_hash`], and the difference is the whole reason both
//! files exist. `HashMap`'s iteration order is documented as unspecified, so its only description
//! is GPL2 source and it had to be treated as a measured observable. `Random` is the other way
//! round: its Javadoc states the algorithm in full, constants included, and declares that all its
//! methods must produce exactly that sequence. Implementing a published contract is not
//! translation, and `docs/an-unspecified-order-that-reaches-the-output.md` records the distinction.
//!
//! The contract, as published:
//!
//! ```text
//! seed = (seed ^ 0x5DEECE66D) & ((1 << 48) - 1)                       // the constructor
//! seed = (seed * 0x5DEECE66D + 0xB) & ((1 << 48) - 1)                 // next(bits)
//! return (int)(seed >>> (48 - bits))
//! ```
//!
//! Every other method is defined in terms of `next`, and the derived ones are spelled out too, so
//! `nextInt(bound)`'s rejection loop and `nextDouble`'s two draws are part of the contract rather
//! than of an implementation.
//!
//! # Why GATK needs it to be exact
//!
//! `Utils.getRandomGenerator()` is a single static `new Random(47382911L)`, seeded once at class
//! initialisation. Every draw in a run therefore comes from one deterministic stream, and its
//! *position* in that stream depends on how many draws everything before it took. A downsampler
//! that consumes one value too few or too many does not produce a differently-shuffled result, it
//! produces a differently-shuffled result **for every later consumer as well**.
//!
//! That is why this is a separate slice with its own conformance suite: the sequence has to be
//! right before anything is allowed to draw from it.

/// `Utils.GATK_RANDOM_SEED`, the constant GATK seeds its static generator with.
pub const GATK_RANDOM_SEED: i64 = 47_382_911;

const MULTIPLIER: i64 = 0x5DEE_CE66D;
const ADDEND: i64 = 0xB;
const MASK: i64 = (1 << 48) - 1;

/// `java.util.Random`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JavaRandom {
    seed: i64,
    /// `nextNextGaussian` while `haveNextNextGaussian` holds, kept as bits so the type stays `Eq`.
    next_next_gaussian: Option<u64>,
}

impl JavaRandom {
    /// `new Random(long seed)`, whose scramble is part of the specification.
    pub fn new(seed: i64) -> Self {
        JavaRandom {
            seed: (seed ^ MULTIPLIER) & MASK,
            next_next_gaussian: None,
        }
    }

    /// `Utils.getRandomGenerator()`: the one GATK seeds at class initialisation.
    pub fn gatk() -> Self {
        JavaRandom::new(GATK_RANDOM_SEED)
    }

    /// `setSeed`, which is what `Utils.resetRandomGenerator` calls.
    pub fn set_seed(&mut self, seed: i64) {
        self.seed = (seed ^ MULTIPLIER) & MASK;
        // `haveNextNextGaussian = false`: a reseed discards the cached second deviate.
        self.next_next_gaussian = None;
    }

    /// `protected int next(int bits)`.
    pub fn next(&mut self, bits: u32) -> i32 {
        self.seed = self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND) & MASK;
        (((self.seed as u64) >> (48 - bits)) as u32) as i32
    }

    /// `nextInt()`.
    pub fn next_int(&mut self) -> i32 {
        self.next(32)
    }

    /// `nextInt(int bound)`.
    ///
    /// Two things here are contract rather than optimisation: a power-of-two bound takes a
    /// different path and therefore a *different* value from the same seed than a nearby non-power
    /// of two would, and the general path loops until the value is unbiased, so it can consume
    /// more than one draw. Both are in the published algorithm.
    pub fn next_int_bound(&mut self, bound: i32) -> i32 {
        assert!(bound > 0, "bound must be positive");
        if (bound & -bound) == bound {
            // A power of two: the high bits of one draw, exactly.
            return (((bound as i64).wrapping_mul(self.next(31) as i64)) >> 31) as i32;
        }
        loop {
            let bits = self.next(31);
            let value = bits % bound;
            // The rejection test, which is what keeps the distribution uniform and what makes the
            // number of draws depend on the values themselves.
            if bits.wrapping_sub(value).wrapping_add(bound - 1) >= 0 {
                return value;
            }
        }
    }

    /// `nextLong()`: two 32-bit draws, the first shifted, and *added* rather than or-ed, which the
    /// specification notes can overflow.
    pub fn next_long(&mut self) -> i64 {
        ((self.next(32) as i64) << 32).wrapping_add(self.next(32) as i64)
    }

    /// `nextBoolean()`.
    pub fn next_boolean(&mut self) -> bool {
        self.next(1) != 0
    }

    /// `nextDouble()`: 26 bits then 27, scaled by 2^-53.
    pub fn next_double(&mut self) -> f64 {
        let high = (self.next(26) as i64) << 27;
        let low = self.next(27) as i64;
        (high.wrapping_add(low) as f64) * (1.0f64 / (1u64 << 53) as f64)
    }

    /// `nextGaussian()`: Marsaglia's polar method, as the specification states it.
    ///
    /// ```text
    /// do { v1 = 2 * nextDouble() - 1; v2 = 2 * nextDouble() - 1; s = v1 * v1 + v2 * v2; }
    /// while (s >= 1 || s == 0);
    /// multiplier = StrictMath.sqrt(-2 * StrictMath.log(s) / s);
    /// nextNextGaussian = v2 * multiplier; haveNextNextGaussian = true; return v1 * multiplier;
    /// ```
    ///
    /// Each accepted pair yields two deviates: the first is returned and the second is cached for
    /// the next call, so a consumer drawing one Gaussian moves the stream by a whole pair only
    /// every other time. The logarithm is `StrictMath`'s, which is FDLIBM and not the correctly
    /// rounded one ([`jmath::strict_math::log`], htsjdk-rs decision 0044).
    pub fn next_gaussian(&mut self) -> f64 {
        if let Some(bits) = self.next_next_gaussian.take() {
            return f64::from_bits(bits);
        }
        loop {
            let v1 = 2.0 * self.next_double() - 1.0;
            let v2 = 2.0 * self.next_double() - 1.0;
            let s = v1 * v1 + v2 * v2;
            if s >= 1.0 || s == 0.0 {
                continue;
            }
            let multiplier = jmath::strict_math::sqrt(-2.0 * jmath::strict_math::log(s) / s);
            self.next_next_gaussian = Some((v2 * multiplier).to_bits());
            return v1 * multiplier;
        }
    }

    /// `nextFloat()`: 24 bits scaled by 2^-24.
    pub fn next_float(&mut self) -> f32 {
        self.next(24) as f32 / (1u32 << 24) as f32
    }
}
