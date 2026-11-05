//! `numpy.random.RandomState` as scikit-learn reaches it: MT19937 seeded by an integer, and the
//! legacy draws its bagging and its trees ask for.
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
}
