//! `SeriesStats`: the accumulator both ground-truth tools keep their numbers in.
//!
//! A sorted map from value to count rather than a list, and every statistic is computed off that
//! map. That is what makes several of them not the statistic their name suggests: the percentile is
//! an observed value rather than an interpolation, and `getUniq` counts bins rather than distinct
//! numbers, which is not the same thing once a negative zero is involved.

use std::collections::BTreeMap;

/// A `Double` key ordered the way `java.lang.Double.compareTo` orders it: -0.0 BEFORE 0.0, and NaN
/// after everything.
#[derive(Debug, Clone, Copy)]
pub struct Key(pub f64);

impl Key {
    /// `Double.compare`: the NUMERIC comparison first, and only when neither side is less or
    /// greater are the raw bits compared as SIGNED longs. That last step is what separates -0.0,
    /// whose bits are the most negative long, from 0.0, whose bits are zero, and what puts every
    /// NaN after everything else.
    fn compare(a: f64, b: f64) -> std::cmp::Ordering {
        if a < b {
            return std::cmp::Ordering::Less;
        }
        if a > b {
            return std::cmp::Ordering::Greater;
        }
        // `doubleToLongBits` collapses every NaN to one pattern, so two of them are equal.
        let bits = |value: f64| {
            if value.is_nan() {
                f64::NAN.to_bits() as i64
            } else {
                value.to_bits() as i64
            }
        };
        bits(a).cmp(&bits(b))
    }
}

impl PartialEq for Key {
    fn eq(&self, other: &Self) -> bool {
        Key::compare(self.0, other.0) == std::cmp::Ordering::Equal
    }
}
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        Key::compare(self.0, other.0)
    }
}

/// The accumulator.
#[derive(Debug, Clone, Default)]
pub struct SeriesStats {
    last: f64,
    count: i32,
    sum: f64,
    min: f64,
    max: f64,
    bins: BTreeMap<Key, i32>,
    /// How many values arrived through `add(int)`, which is what decides the CSV's format.
    int_count: i32,
}

impl SeriesStats {
    pub fn new() -> SeriesStats {
        SeriesStats {
            last: f64::NAN,
            count: 0,
            sum: 0.0,
            min: f64::NAN,
            max: f64::NAN,
            bins: BTreeMap::new(),
            int_count: 0,
        }
    }

    /// `add(int)`, which calls `add(double)` FIRST and only then counts the value as an integer.
    pub fn add_int(&mut self, value: i32) {
        self.add(value as f64);
        self.int_count += 1;
    }

    /// `add(double)`. `Math.min` and `Math.max` propagate a NaN, so one NaN poisons both bounds for
    /// every later value.
    pub fn add(&mut self, value: f64) {
        self.last = value;
        self.sum += value;
        if self.count > 0 {
            self.min = java_min(self.min, value);
            self.max = java_max(self.max, value);
        } else {
            self.min = value;
            self.max = value;
        }
        self.count += 1;
        *self.bins.entry(Key(value)).or_insert(0) += 1;
    }

    pub fn last(&self) -> f64 {
        self.last
    }

    pub fn count(&self) -> i32 {
        self.count
    }

    pub fn min(&self) -> f64 {
        if self.count != 0 {
            self.min
        } else {
            f64::NAN
        }
    }

    pub fn max(&self) -> f64 {
        if self.count != 0 {
            self.max
        } else {
            f64::NAN
        }
    }

    /// `getUniq`: the number of BINS, so a negative zero and a zero count as two.
    pub fn uniq(&self) -> usize {
        self.bins.len()
    }

    pub fn mean(&self) -> f64 {
        if self.count != 0 {
            self.sum / self.count as f64
        } else {
            f64::NAN
        }
    }

    pub fn median(&self) -> f64 {
        self.percentile(50.0)
    }

    /// `getPercentile`.
    ///
    /// The index is TRUNCATED and the walk returns a bin KEY, so the answer is always a value that
    /// was actually added. A single value short-circuits to the LAST added rather than to the only
    /// bin.
    pub fn percentile(&self, percentile: f64) -> f64 {
        if self.count == 0 {
            return f64::NAN;
        }
        if self.count == 1 {
            return self.last;
        }
        let target = (self.count as f64 * percentile / 100.0) as i32;
        let mut index = 0;
        for (key, size) in &self.bins {
            if target >= index && target < index + size {
                return key.0;
            }
            index += size;
        }
        // Past the last bin, which is where a percentile of 100 lands.
        self.bins
            .keys()
            .next_back()
            .map(|key| key.0)
            .unwrap_or(f64::NAN)
    }

    /// `getStd`, which divides by the COUNT rather than by the count less one: the population
    /// deviation, not the sample one.
    pub fn std(&self) -> f64 {
        if self.count == 0 {
            return f64::NAN;
        }
        let mean = self.mean();
        let mut sum = 0.0;
        for (key, size) in &self.bins {
            sum += (key.0 - mean).powi(2) * *size as f64;
        }
        (sum / self.count as f64).sqrt()
    }

    pub fn bins(&self) -> &BTreeMap<Key, i32> {
        &self.bins
    }

    /// `isIntKeys`: every value arrived through `add(int)`. A property of the ADD PATH rather than
    /// of the values, so two whole numbers added as doubles are NOT integer-keyed. An empty series
    /// is, because both counts are zero.
    pub fn is_int_keys(&self) -> bool {
        self.count == self.int_count
    }

    /// `csvWrite`'s body: `%d` under integer keys and `%f` otherwise, for the WHOLE file.
    pub fn csv(&self) -> String {
        let mut out = String::from("value,count\n");
        let integers = self.is_int_keys();
        for (key, size) in &self.bins {
            if integers {
                out.push_str(&format!("{},{size}\n", key.0 as i32));
            } else {
                out.push_str(&format!("{:.6},{size}\n", key.0));
            }
        }
        out
    }

    /// `toDigest`.
    ///
    /// Under integer keys the three statistics are cast with `(int)`, and `(int)Double.NaN` is
    /// ZERO in Java, so the digest of an EMPTY series claims a minimum of nought while every other
    /// reader of it is told NaN.
    pub fn to_digest(&self) -> String {
        if self.is_int_keys() {
            format!(
                "count={}, min={}, max={}, median={}, bin.count={}",
                self.count,
                java_double_to_int(self.min()),
                java_double_to_int(self.max()),
                java_double_to_int(self.median()),
                self.bins.len()
            )
        } else {
            format!(
                "count={}, min={:.6}, max={:.6}, median={:.6}, bin.count={}",
                self.count,
                self.min(),
                self.max(),
                self.median(),
                self.bins.len()
            )
        }
    }
}

/// `(int)` on a double: NaN becomes zero, and anything out of range saturates.
pub fn java_double_to_int(value: f64) -> i32 {
    if value.is_nan() {
        0
    } else if value >= i32::MAX as f64 {
        i32::MAX
    } else if value <= i32::MIN as f64 {
        i32::MIN
    } else {
        value as i32
    }
}

/// `Math.min`, which returns NaN if either argument is one.
fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.min(b)
    }
}

/// `Math.max`, the same.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}
