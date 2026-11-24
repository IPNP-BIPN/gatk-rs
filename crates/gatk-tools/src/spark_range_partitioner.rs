//! Spark 3.5.0's `RangePartitioner` bounds (Scala 2.13), sketch included: what `sortByKey` cuts a
//! sort into, and so where Disq's part files, and the BGZF blocks of a merged BAM, end.
//!
//! Each input partition is sketched by `SamplingUtils.reservoirSampleAndCount` with an
//! `XORShiftRandom` seeded from `byteswap32(index ^ (shift << 16))`, where `shift` is the id of the
//! RDD `rdd.map(_._1)` that `RangePartitioner` builds. That id is how many RDDs the tool made
//! before, so the caller passes it: the reference logs it, plus one, as the `MapPartitionsRDD`
//! of the `sortByKey` result stage.

use std::cmp::Ordering;

/// `scala.util.hashing.MurmurHash3.bytesHash` over a whole number of four-byte words.
fn bytes_hash(data: &[u8], seed: i32) -> i32 {
    let mut h = seed;
    for word in data.as_chunks::<4>().0 {
        let k = i32::from_le_bytes(*word);
        let k = k
            .wrapping_mul(0xcc9e_2d51_u32 as i32)
            .rotate_left(15)
            .wrapping_mul(0x1b87_3593);
        h = (h ^ k)
            .rotate_left(13)
            .wrapping_mul(5)
            .wrapping_add(0xe654_6b64_u32 as i32);
    }
    let mut h = (h ^ data.len() as i32) as u32;
    h ^= h >> 16;
    h = h.wrapping_mul(0x85eb_ca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2_ae35);
    h ^= h >> 16;
    h as i32
}

/// `scala.util.hashing.byteswap32`.
pub fn byteswap32(v: i32) -> i32 {
    let hc = v.wrapping_mul(0x9e37_75cd_u32 as i32).swap_bytes();
    hc.wrapping_mul(0x9e37_75cd_u32 as i32)
}

/// `org.apache.spark.util.random.XORShiftRandom`, with `java.util.Random.nextDouble` over it.
pub struct XorShiftRandom {
    seed: i64,
}

impl XorShiftRandom {
    /// `new XORShiftRandom(init)`: the seed is `hashSeed(init)`.
    pub fn new(init: i64) -> Self {
        let bytes = init.to_be_bytes();
        let low = bytes_hash(&bytes, 0x3c07_4a61);
        let high = bytes_hash(&bytes, low);
        XorShiftRandom {
            seed: (i64::from(high) << 32) | (i64::from(low) & 0xFFFF_FFFF),
        }
    }

    fn next(&mut self, bits: u32) -> i64 {
        let mut next = self.seed ^ (self.seed << 21);
        next ^= ((next as u64) >> 35) as i64;
        next ^= next << 4;
        self.seed = next;
        i64::from((next & ((1i64 << bits) - 1)) as i32)
    }

    pub fn next_double(&mut self) -> f64 {
        ((self.next(26) << 27) + self.next(27)) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

/// `SamplingUtils.reservoirSampleAndCount` for one input partition, fed in its order.
pub struct Reservoir<T> {
    k: usize,
    seen: u64,
    sample: Vec<T>,
    random: Option<XorShiftRandom>,
    seed: i64,
}

impl<T> Reservoir<T> {
    /// The sketch of partition `index`: `k` keys, seeded as `RangePartitioner.sketch` seeds it.
    pub fn new(k: usize, index: i32, shift: i32) -> Self {
        Reservoir {
            k,
            seen: 0,
            sample: Vec::new(),
            random: None,
            seed: i64::from(byteswap32(index ^ (shift << 16))),
        }
    }

    pub fn push(&mut self, item: T) {
        self.seen += 1;
        if self.sample.len() < self.k {
            self.sample.push(item);
            return;
        }
        let random = self
            .random
            .get_or_insert_with(|| XorShiftRandom::new(self.seed));
        let replacement = (random.next_double() * self.seen as f64) as i64;
        if replacement < self.k as i64 {
            self.sample[replacement as usize] = item;
        }
    }

    /// `(sample, n)`.
    pub fn finish(self) -> (Vec<T>, u64) {
        (self.sample, self.seen)
    }
}

/// `sampleSizePerPartition`, and `sampleSize` beside it, for `partitions` output partitions over
/// `inputs` input partitions.
pub fn sample_sizes(partitions: usize, inputs: usize) -> (usize, f64) {
    let sample_size = (20.0 * partitions as f64).min(1e6);
    let per_partition = (3.0 * sample_size / inputs.max(1) as f64).ceil() as usize;
    (per_partition, sample_size)
}

/// An input partition `fraction * n > sampleSizePerPartition`, which Spark samples again by
/// `BernoulliSampler`; the port does not.
#[derive(Debug)]
pub struct Imbalanced;

/// The bounds `RangePartitioner` takes for `partitions` output partitions, from each input
/// partition's sketch in index order: weight `n / sample.length` as a float, then
/// `determineBounds(candidates, min(partitions, candidates.size))`.
pub fn bounds<T: Clone>(
    sketches: Vec<(Vec<T>, u64)>,
    partitions: usize,
    compare: impl Fn(&T, &T) -> Ordering,
) -> Result<Vec<T>, Imbalanced> {
    if partitions <= 1 {
        return Ok(Vec::new());
    }
    let (per_partition, sample_size) = sample_sizes(partitions, sketches.len());
    let items: u64 = sketches.iter().map(|(_, n)| n).sum();
    if items == 0 {
        return Ok(Vec::new());
    }
    let fraction = (sample_size / items.max(1) as f64).min(1.0);
    let mut candidates: Vec<(T, f32)> = Vec::new();
    for (sample, n) in sketches {
        if fraction * n as f64 > per_partition as f64 {
            return Err(Imbalanced);
        }
        let weight = (n as f64 / sample.len() as f64) as f32;
        candidates.extend(sample.into_iter().map(|key| (key, weight)));
    }
    // `sortBy`, stable as Java's object sort is.
    candidates.sort_by(|a, b| compare(&a.0, &b.0));
    let partitions = partitions.min(candidates.len());
    let sum: f64 = candidates.iter().map(|(_, w)| f64::from(*w)).sum();
    let step = sum / partitions as f64;
    let mut cumulative = 0.0f64;
    let mut target = step;
    let mut found: Vec<T> = Vec::new();
    for (key, weight) in &candidates {
        if found.len() >= partitions - 1 {
            break;
        }
        cumulative += f64::from(*weight);
        if cumulative >= target && found.last().is_none_or(|last| compare(key, last).is_gt()) {
            found.push(key.clone());
            target += step;
        }
    }
    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A partition no larger than `k` is its own sample, unshuffled, with nothing drawn.
    #[test]
    fn a_small_partition_is_taken_whole() {
        let mut reservoir = Reservoir::new(5, 0, 15);
        (0..5).for_each(|i| reservoir.push(i));
        assert_eq!(reservoir.finish(), (vec![0, 1, 2, 3, 4], 5));
    }

    /// Every key sampled: the bounds of the unweighted walk the port took before the sketch.
    #[test]
    fn whole_sketches_cut_by_count() {
        let keys: Vec<i32> = (0..30).collect();
        let found = bounds(vec![(keys, 30)], 3, |a, b| a.cmp(b)).unwrap();
        assert_eq!(found, vec![9, 19]);
    }

    /// Past `k` every key is still counted, and the sample holds `k` of them, each seen.
    #[test]
    fn a_large_partition_is_sampled() {
        let mut reservoir = Reservoir::new(10, 3, 20);
        (0..1000).for_each(|i| reservoir.push(i));
        let (sample, n) = reservoir.finish();
        assert_eq!(n, 1000);
        assert_eq!(sample.len(), 10);
        assert!(sample.iter().any(|&i| i >= 10));
    }
}
