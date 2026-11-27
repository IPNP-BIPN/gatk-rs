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
//! initialisation. Every draw taken through it comes from one deterministic stream, and its
//! *position* in that stream depends on how many draws everything before it took. (Some tools
//! build their own `new Random(seed)` instead, `KernelSegmenter` and `DownsampleByDuplicateSet`
//! among them; each of those is a separate stream with the same property.) A downsampler
//! that consumes one value too few or too many does not produce a differently-shuffled result, it
//! produces a differently-shuffled result **for every later consumer as well**.
//!
//! That is why this is a separate slice with its own conformance suite: the sequence has to be
//! right before anything is allowed to draw from it.
//!
//! # Reading this file without reading Rust
//!
//! - `i64` is Java's `long` and `i32` is `int`. `u32` and `u64` are unsigned types Java lacks;
//!   they appear only on the way to Java's `>>>` (see `next`).
//! - `wrapping_mul`, `wrapping_add` and `wrapping_sub` are Java's `long` and `int` arithmetic,
//!   which wraps silently on overflow. Plain `*`, `+` and `-` would panic in a debug build
//!   instead, and the multiply in `next` overflows a long for most states.
//! - `x as T` converts: to a narrower integer it keeps the low bits (Java's `(int)` cast), to a
//!   wider signed one it sign-extends (Java's promotion of `int` to `long`), and between signed
//!   and unsigned of one width it reinterprets the same bits.
//! - `&mut self` on a method means it advances this generator in place. A caller holding a
//!   `&mut JavaRandom` draws from that one stream, not from a copy of it.
//! - `Option<T>` is "a value or nothing", checked by the compiler.
//! - A function's last expression is its return value; `return` is written only for early exits.

/// `Utils.GATK_RANDOM_SEED`, the constant GATK seeds its static generator with.
///
/// ```java
/// private static final long GATK_RANDOM_SEED = 47382911L;
/// ```
///
/// Both of `Utils`'s static generators are seeded with it; [`crate::well19937c::Well19937c::gatk`]
/// reads it from here.
pub const GATK_RANDOM_SEED: i64 = 47_382_911;

const MULTIPLIER: i64 = 0x5DEE_CE66D;
const ADDEND: i64 = 0xB;
const MASK: i64 = (1 << 48) - 1;

/// `java.util.Random`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JavaRandom {
    /// The 48-bit state, stored after the constructor's scramble, so never the caller's seed.
    seed: i64,
    /// `nextNextGaussian` while `haveNextNextGaussian` holds, kept as bits so the type stays `Eq`.
    next_next_gaussian: Option<u64>,
}

impl JavaRandom {
    /// `new Random(long seed)`, whose scramble is part of the specification.
    ///
    /// Not every GATK `Random` is the `Utils` one: several tools build their own stream, which
    /// no other consumer shares, for example `DownsampleByDuplicateSet`:
    ///
    /// ```java
    /// private final Random rng = new Random(RANDOM_SEED);
    /// ```
    pub fn new(seed: i64) -> Self {
        JavaRandom {
            seed: (seed ^ MULTIPLIER) & MASK,
            next_next_gaussian: None,
        }
    }

    /// `Utils.getRandomGenerator()`: the one GATK seeds at class initialisation.
    ///
    /// ```java
    /// private static final Random randomGenerator = new Random(GATK_RANDOM_SEED);
    /// public static Random getRandomGenerator() { return randomGenerator; }
    /// ```
    ///
    /// The Java hands every caller the same static instance. This returns a *fresh* generator at
    /// the start of that stream on each call, so the port keeps the single stream by creating it
    /// once and passing the same `&mut JavaRandom` to every consumer. A second call mid-run
    /// restarts the stream, which the Java does only through `resetRandomGenerator`.
    pub fn gatk() -> Self {
        JavaRandom::new(GATK_RANDOM_SEED)
    }

    /// `setSeed`, which is what `Utils.resetRandomGenerator` calls.
    ///
    /// ```java
    /// public static void resetRandomGenerator() {
    ///     randomGenerator.setSeed(GATK_RANDOM_SEED);
    ///     randomDataGenerator.reSeed(GATK_RANDOM_SEED);
    /// }
    /// ```
    ///
    /// The second line is [`crate::well19937c::Well19937c::set_seed_long`]: a reset reaches both.
    pub fn set_seed(&mut self, seed: i64) {
        self.seed = (seed ^ MULTIPLIER) & MASK;
        // `haveNextNextGaussian = false`: a reseed discards the cached second deviate.
        self.next_next_gaussian = None;
    }

    /// `protected int next(int bits)`.
    pub fn next(&mut self, bits: u32) -> i32 {
        // The multiply wraps exactly as Java's `long` does; `& MASK` keeps the 48-bit state.
        self.seed = self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND) & MASK;
        // Java's `(int)(seed >>> (48 - bits))`. The state is non-negative after `& MASK`, so
        // `>>` and `>>>` agree here; `as u64` only makes the logical shift explicit. The `as u32`
        // then `as i32` keep the low 32 bits and reinterpret them as an `int`, Java's `(int)`.
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
        // `assert!` panics (stops the program) when its condition is false. It stands in for the
        // reference's rejection of a non-positive bound.
        assert!(bound > 0, "bound must be positive");
        if bound.isolate_lowest_one() == bound {
            // A power of two: the high bits of one draw, exactly.
            return (((bound as i64).wrapping_mul(self.next(31) as i64)) >> 31) as i32;
        }
        // `loop` is Java's `while (true)`: only the `return` leaves it, one draw per pass.
        loop {
            let bits = self.next(31);
            let value = bits % bound;
            // The rejection test, which is what keeps the distribution uniform and what makes the
            // number of draws depend on the values themselves.
            // The sum wraps negative exactly when `bits - value` starts the last, incomplete run
            // of `bound` values below 2^31, and those draws are the ones rejected.
            if bits.wrapping_sub(value).wrapping_add(bound - 1) >= 0 {
                return value;
            }
        }
    }

    /// `nextLong()`: two 32-bit draws, the first shifted, and *added* rather than or-ed, which the
    /// specification notes can overflow.
    pub fn next_long(&mut self) -> i64 {
        // `as i64` sign-extends the second draw like Java's `int` to `long` promotion, so a
        // negative low word borrows one from the high word: that is the add, not an or.
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
        // `take()` reads the cached deviate and clears it in one step, which is the
        // `haveNextNextGaussian` test and its reset to false. `if let Some(bits) = ...` runs the
        // block only when a value was cached. The cache holds the double's bit pattern
        // (`to_bits` / `from_bits`), an exact round trip.
        if let Some(bits) = self.next_next_gaussian.take() {
            return f64::from_bits(bits);
        }
        loop {
            let v1 = 2.0 * self.next_double() - 1.0;
            let v2 = 2.0 * self.next_double() - 1.0;
            let s = v1 * v1 + v2 * v2;
            // The specification's `do ... while`: `continue` draws a fresh pair.
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
