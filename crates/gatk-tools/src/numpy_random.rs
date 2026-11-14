//! `numpy.random.RandomState` as scikit-learn reaches it: MT19937 seeded by an integer, and the
//! legacy draws its bagging and its trees ask for.
//! Also `Generator(PCG64(seed))` as pymc's `reseed_rngs` builds it from a `SeedSequence` and spawns,
//! with the ziggurat normal pytensor's `NormalRV` draws (`random/bit_generator.pyx`,
//! `random/src/pcg64/pcg64.h`).
//!
//! Ported from numpy 1.26.2 (BSD 3-Clause), `random/src/mt19937/mt19937.c`,
//! `random/src/distributions/distributions.c` and `random/mtrand.pyx`, and checked by the stream
//! the isolation forest draws from it: the seeds in a scorer the reference pickled are this
//! generator's output.

const STATE: usize = 624;
const PERIOD: usize = 397;

/// `RandomState(seed)`, for a seed `int` in `[0, 2^32)`.
#[derive(Debug, Clone)]
pub struct RandomState {
    key: [u32; STATE],
    position: usize,
}

impl RandomState {
    /// `mt19937_seed`.
    pub fn new(seed: u32) -> RandomState {
        let mut key = [0u32; STATE];
        let mut value = seed;
        for (position, slot) in key.iter_mut().enumerate() {
            *slot = value;
            value = 1812433253u32
                .wrapping_mul(value ^ (value >> 30))
                .wrapping_add(position as u32 + 1);
        }
        RandomState {
            key,
            position: STATE,
        }
    }

    /// `mt19937_gen`, the whole state at a time.
    fn generate(&mut self) {
        const MATRIX: u32 = 0x9908_b0df;
        const UPPER: u32 = 0x8000_0000;
        const LOWER: u32 = 0x7fff_ffff;
        let key = &mut self.key;
        for i in 0..STATE - PERIOD {
            let y = (key[i] & UPPER) | (key[i + 1] & LOWER);
            key[i] = key[i + PERIOD] ^ (y >> 1) ^ (if y & 1 == 1 { MATRIX } else { 0 });
        }
        for i in STATE - PERIOD..STATE - 1 {
            let y = (key[i] & UPPER) | (key[i + 1] & LOWER);
            key[i] = key[i + PERIOD - STATE] ^ (y >> 1) ^ (if y & 1 == 1 { MATRIX } else { 0 });
        }
        let y = (key[STATE - 1] & UPPER) | (key[0] & LOWER);
        key[STATE - 1] = key[PERIOD - 1] ^ (y >> 1) ^ (if y & 1 == 1 { MATRIX } else { 0 });
        self.position = 0;
    }

    /// `mt19937_next32`.
    pub fn next_u32(&mut self) -> u32 {
        if self.position == STATE {
            self.generate();
        }
        let mut y = self.key[self.position];
        self.position += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// `mt19937_next_double`, which `RandomState.uniform(size=n)` returns as it stands over
    /// `[0, 1)`.
    pub fn next_double(&mut self) -> f64 {
        let a = i64::from(self.next_u32() >> 5);
        let b = i64::from(self.next_u32() >> 6);
        (a as f64 * 67108864.0 + b as f64) / 9007199254740992.0
    }

    /// The smallest mask of ones covering `value`.
    fn mask(value: u64) -> u64 {
        let mut mask = value;
        mask |= mask >> 1;
        mask |= mask >> 2;
        mask |= mask >> 4;
        mask |= mask >> 8;
        mask |= mask >> 16;
        mask | (mask >> 32)
    }

    /// `randint(low, high)` of the default dtype: `_rand_int64` with `_masked=True`, drawing 32
    /// bits at a time for a range under 2^32, which is every range scikit-learn asks for here.
    pub fn randint(&mut self, low: u64, high: u64) -> u64 {
        let range = high - low - 1;
        if range == 0 {
            return low;
        }
        debug_assert!(range < 0xffff_ffff);
        let mask = Self::mask(range);
        loop {
            let value = u64::from(self.next_u32()) & mask;
            if value <= range {
                return low + value;
            }
        }
    }

    /// `random_interval`, what `shuffle` draws its swaps with.
    fn interval(&mut self, maximum: u64) -> u64 {
        if maximum == 0 {
            return 0;
        }
        let mask = Self::mask(maximum);
        loop {
            let value = u64::from(self.next_u32()) & mask;
            if value <= maximum {
                return value;
            }
        }
    }

    /// `permutation(n)`: `arange(n)` shuffled from the last element down.
    pub fn permutation(&mut self, n: usize) -> Vec<usize> {
        let mut out: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = self.interval(i as u64) as usize;
            out.swap(i, j);
        }
        out
    }
}

/// `sklearn.utils.random.sample_without_replacement(population, samples, method="auto")`, of
/// scikit-learn 1.3.2 (BSD 3-Clause), `utils/_random.pyx`.
pub fn sample_without_replacement(
    population: usize,
    samples: usize,
    random: &mut RandomState,
) -> Vec<usize> {
    let ratio = if population == 0 {
        1.0
    } else {
        samples as f64 / population as f64
    };
    if ratio > 0.01 && ratio < 0.99 {
        let mut permuted = random.permutation(population);
        permuted.truncate(samples);
        return permuted;
    }
    if ratio < 0.2 {
        // `_sample_without_replacement_with_tracking_selection`.
        let mut selected = std::collections::BTreeSet::new();
        let mut out = Vec::with_capacity(samples);
        for _ in 0..samples {
            let mut j = random.randint(0, population as u64) as usize;
            while selected.contains(&j) {
                j = random.randint(0, population as u64) as usize;
            }
            selected.insert(j);
            out.push(j);
        }
        return out;
    }
    // `_sample_without_replacement_with_reservoir_sampling`.
    let mut out: Vec<usize> = (0..samples).collect();
    for i in samples..population {
        let j = random.randint(0, i as u64 + 1) as usize;
        if j < samples {
            out[j] = i;
        }
    }
    out
}

/// `SeedSequence(entropy, spawn_key=...)` with the default pool of four words, for an entropy and
/// spawn keys under 2^32 (one word each), which is what `reseed_rngs` passes.
#[derive(Debug, Clone)]
pub struct SeedSequence {
    entropy: u32,
    spawn_key: Vec<u32>,
    pool: [u32; 4],
    children: u32,
}

const INIT_A: u32 = 0x43b0_d7e5;
const MULT_A: u32 = 0x931e_8875;
const INIT_B: u32 = 0x8b51_f9dd;
const MULT_B: u32 = 0x58f3_8ded;
const MIX_MULT_L: u32 = 0xca01_f9dd;
const MIX_MULT_R: u32 = 0x4973_f715;

fn hashmix(value: u32, hash: &mut u32) -> u32 {
    let mut value = value ^ *hash;
    *hash = hash.wrapping_mul(MULT_A);
    value = value.wrapping_mul(*hash);
    value ^ (value >> 16)
}

fn mix(x: u32, y: u32) -> u32 {
    let result = MIX_MULT_L
        .wrapping_mul(x)
        .wrapping_sub(MIX_MULT_R.wrapping_mul(y));
    result ^ (result >> 16)
}

impl SeedSequence {
    pub fn new(entropy: u32) -> SeedSequence {
        SeedSequence::with_key(entropy, Vec::new())
    }

    fn with_key(entropy: u32, spawn_key: Vec<u32>) -> SeedSequence {
        // `get_assembled_entropy`: with a spawn key, the run entropy is padded to the pool size.
        let mut assembled = vec![entropy];
        if !spawn_key.is_empty() {
            assembled.resize(4, 0);
        }
        assembled.extend_from_slice(&spawn_key);
        let mut pool = [0u32; 4];
        let mut hash = INIT_A;
        for (i, slot) in pool.iter_mut().enumerate() {
            *slot = hashmix(assembled.get(i).copied().unwrap_or(0), &mut hash);
        }
        for source in 0..4 {
            for destination in 0..4 {
                if source != destination {
                    pool[destination] = mix(pool[destination], hashmix(pool[source], &mut hash));
                }
            }
        }
        for &word in assembled.iter().skip(4) {
            for slot in pool.iter_mut() {
                *slot = mix(*slot, hashmix(word, &mut hash));
            }
        }
        SeedSequence {
            entropy,
            spawn_key,
            pool,
            children: 0,
        }
    }

    pub fn spawn(&mut self, n: u32) -> Vec<SeedSequence> {
        let children = (self.children..self.children + n)
            .map(|i| {
                let mut key = self.spawn_key.clone();
                key.push(i);
                SeedSequence::with_key(self.entropy, key)
            })
            .collect();
        self.children += n;
        children
    }

    /// `generate_state(n, np.uint64)`.
    pub fn generate_state_u64(&self, n: usize) -> Vec<u64> {
        let mut hash = INIT_B;
        let words: Vec<u32> = (0..2 * n)
            .map(|i| {
                let mut value = self.pool[i % 4] ^ hash;
                hash = hash.wrapping_mul(MULT_B);
                value = value.wrapping_mul(hash);
                value ^ (value >> 16)
            })
            .collect();
        words
            .chunks(2)
            .map(|pair| u64::from(pair[0]) | (u64::from(pair[1]) << 32))
            .collect()
    }
}

const PCG_MULTIPLIER: u128 = (2549297995355413924u128 << 64) | 4865540595714422341u128;

/// `Generator(PCG64(seed_sequence))`: PCG XSL-RR 128/64.
#[derive(Debug, Clone)]
pub struct Pcg64 {
    pub state: u128,
    pub inc: u128,
}

impl Pcg64 {
    pub fn new(seed: &SeedSequence) -> Pcg64 {
        let words = seed.generate_state_u64(4);
        let initstate = (u128::from(words[0]) << 64) | u128::from(words[1]);
        let initseq = (u128::from(words[2]) << 64) | u128::from(words[3]);
        let mut rng = Pcg64 {
            state: 0,
            inc: (initseq << 1) | 1,
        };
        rng.step();
        rng.state = rng.state.wrapping_add(initstate);
        rng.step();
        rng
    }

    fn step(&mut self) {
        self.state = self
            .state
            .wrapping_mul(PCG_MULTIPLIER)
            .wrapping_add(self.inc);
    }

    pub fn next_u64(&mut self) -> u64 {
        self.step();
        let value = ((self.state >> 64) as u64) ^ (self.state as u64);
        value.rotate_right((self.state >> 122) as u32)
    }

    pub fn next_double(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9007199254740992.0)
    }

    /// `random_standard_normal`, the ziggurat; `exp` and `log1p` are glibc's, as numpy calls them.
    pub fn standard_normal(&mut self) -> f64 {
        use crate::glibc_libm::{exp, log1p};
        use crate::numpy_ziggurat::{
            FI_DOUBLE, KI_DOUBLE, WI_DOUBLE, ZIGGURAT_NOR_INV_R, ZIGGURAT_NOR_R,
        };
        loop {
            let mut r = self.next_u64();
            let idx = (r & 0xff) as usize;
            r >>= 8;
            let sign = r & 1;
            let rabs = (r >> 1) & 0x000f_ffff_ffff_ffff;
            let mut x = rabs as f64 * WI_DOUBLE[idx];
            if sign == 1 {
                x = -x;
            }
            if rabs < KI_DOUBLE[idx] {
                return x;
            }
            if idx == 0 {
                loop {
                    let xx = -ZIGGURAT_NOR_INV_R * log1p(-self.next_double());
                    let yy = -log1p(-self.next_double());
                    if yy + yy > xx * xx {
                        return if (rabs >> 8) & 1 == 1 {
                            -(ZIGGURAT_NOR_R + xx)
                        } else {
                            ZIGGURAT_NOR_R + xx
                        };
                    }
                }
            } else if (FI_DOUBLE[idx - 1] - FI_DOUBLE[idx]) * self.next_double() + FI_DOUBLE[idx]
                < exp(-0.5 * x * x)
            {
                return x;
            }
        }
    }
}

/// pymc's `reseed_rngs(rngs, seed)`: one `PCG64` per spawned child of `SeedSequence(seed)`.
pub fn reseed(seed: u32, n: u32) -> Vec<Pcg64> {
    SeedSequence::new(seed)
        .spawn(n)
        .iter()
        .map(Pcg64::new)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `numpy.random.RandomState(0)`'s first draws, from the published stream.
    #[test]
    fn the_first_doubles_of_seed_zero_are_numpys() {
        let mut random = RandomState::new(0);
        assert_eq!(random.next_double(), 0.5488135039273248);
        assert_eq!(random.next_double(), 0.7151893663724195);
        assert_eq!(random.next_double(), 0.6027633760716439);
    }

    #[test]
    fn the_sampler_keeps_a_whole_population_in_order_without_drawing() {
        let mut random = RandomState::new(1);
        let mut untouched = random.clone();
        assert_eq!(
            sample_without_replacement(4, 4, &mut random),
            vec![0, 1, 2, 3]
        );
        assert_eq!(random.next_u32(), untouched.next_u32());
    }

    /// `np.random.default_rng(42)`'s first draws, from numpy's documentation.
    #[test]
    fn default_rng_42_draws_numpys_normals() {
        let mut rng = Pcg64::new(&SeedSequence::new(42));
        let draws: Vec<f64> = (0..3).map(|_| rng.standard_normal()).collect();
        assert_eq!(format!("{:.8}", draws[0]), "0.30471708");
        assert_eq!(format!("{:.8}", draws[1]), "-1.03998411");
        assert_eq!(format!("{:.8}", draws[2]), "0.75045120");
        let mut rng = Pcg64::new(&SeedSequence::new(42));
        assert_eq!(rng.next_double(), 0.7739560485559633);
    }

    /// What the ploidy reference's trace shows: `approx.rng = RandomState(1984)` hands out
    /// `randint(2**30)` seeds, and `reseed_rngs(rngs, 515143279)` sets one PCG64 to this state.
    #[test]
    fn the_ploidy_references_seeds_and_states_are_reproduced() {
        let mut approx = RandomState::new(1984);
        assert_eq!(approx.randint(0, 1 << 30), 37605596);
        assert_eq!(approx.randint(0, 1 << 30), 515143279);
        let rng = &reseed(515143279, 1)[0];
        assert_eq!(rng.state, 0x898905ca812836b7abe908e5ed64eb1a);
        assert_eq!(rng.inc, 0x3d164b2ca944653a0f756cd48176a211);
    }
}
