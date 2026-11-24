//! `CompressedDataList` and `Histogram`, ported from `org.broadinstitute.hellbender.utils`
//! (GATK 4.7.0.0; both Java files are unchanged since 4.6.2.0).
//!
//! The run-length encoding the allele-specific annotations carry through a gVCF. A list of doubles
//! is binned, the bins are counted, and the counts are written into an INFO field as
//! `value,count,value,count,...` so that combining two gVCFs is adding two count maps.
//!
//! This is **not** htsjdk's `Histogram`, which `htsjdk-rs` already has. That one keys on arbitrary
//! comparable values and reports a mean; this one keys on an integer bin index and reports a
//! median it computes by walking the bins.
//!
//! # The bin index is a rounded floor, with an epsilon in front of the division
//!
//! ```java
//! return Math.round(Math.floor((d+BIN_EPSILON*binSize)/binSize));
//! ```
//!
//! The comment says the epsilon is there "so values exactly on bin boundaries will stay in the same
//! bin", and it is one hundredth of a bin. So the bin edges sit slightly below the round numbers
//! they look like: with the default bin size of 0.1, a value of `-0.0009` bins to zero and not to
//! minus one. A bin index outside `int` range is a `GATKException` rather than a clamp.
//!
//! # The median walks the bins and can answer `None`
//!
//! ```java
//! int medianIndex = (numItems+1)/2;
//! ```
//!
//! Integer division, so for an even count the index is the **lower** of the two middle items. The
//! walk then returns the bin where the running count first exceeds that index, and for an even
//! count averages the two bins around it. An empty histogram answers null, which the reference
//! stores into a `Map<Allele, Double>` and later formats as the missing value.
//!
//! # An empty histogram prints as `NaN`
//!
//! ```java
//! if (keys.length == 0) { return Double.toString(Double.NaN); }
//! ```
//!
//! Not the empty string, and it does reach a record: the allele-specific rank sums print an empty
//! histogram for an allele with no reads, so a site with only reference reads carries
//! `AS_RAW_BaseQRankSum=|NaN` (`crates/gatk-annotation/tests/data/allele_specific_rank_sum.txt.gz`,
//! the `ref-only` rows).
//!
//! # Reading this file without reading Rust
//!
//! - `Result<(), HistogramError>` is a `void` method that may throw: `Ok(())` is the normal
//!   return and `Err(..)` is the `GATKException`. `Result<Option<i32>, _>` is the same for a
//!   method returning an `Integer` that may be `null`.
//! - `Option<T>` is "a value or nothing": `Some(x)` is a non-null `x`, `None` is `null`.
//! - `x as i32`, `x as f64`, `x as usize` are casts. Between integer types `as` keeps the low bits
//!   silently, as Java's `(int)` does on a `long`, so every narrowing cast below sits after a
//!   range check.
//! - `&self` is a method that only reads the object; `&mut self` is one that may change it.
//! - `impl std::fmt::Display for T` is Rust's `toString()`. `.to_string()` calls it, and inside it
//!   `write!(f, ...)` appends to the output.
//! - `*x` reads the value behind a reference, which is what iterating a map hands out.

use std::collections::BTreeMap;

/// `Histogram.BIN_EPSILON`.
const BIN_EPSILON: f64 = 0.01;

/// What the reference raises as `GATKException`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HistogramError {
    /// A bin index outside `int` range: "Histogram values are suspiciously extreme". Reachable
    /// from a value near the limits of a `double`, since the bin index is the value divided by the
    /// bin size.
    ValueTooExtreme,
    /// `add(value, count)` with a count below one.
    NonPositiveCount,
    /// `add(other)` where the two bin sizes differ.
    MismatchedBinSize,
}

/// `CompressedDataList<Integer>`: a map from value to count, iterated and printed in **sorted**
/// order however it was filled.
///
/// The reference's map is a `HashMap` and every read of it goes through a `TreeSet` of the keys or
/// an `Arrays.sort`, so the hash order is never observable. A `BTreeMap` is that, without the
/// intermediate copy.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CompressedDataList {
    value_counts: BTreeMap<i32, i32>,
}

impl CompressedDataList {
    /// An empty list. `Self` means "the type this block is for", so `Self::default()` is the
    /// derived empty value: an empty map.
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether anything has been added. Note that a value added with a count of zero still makes
    /// the list non-empty, because the key exists.
    pub fn is_empty(&self) -> bool {
        self.value_counts.is_empty()
    }

    /// `getValueCounts()`, in sorted order.
    pub fn value_counts(&self) -> &BTreeMap<i32, i32> {
        &self.value_counts
    }

    /// `add(val)`.
    pub fn add(&mut self, value: i32) {
        self.add_count(value, 1);
    }

    /// `add(val, count)`, which has **no** guard against a non-positive count. The guard is in
    /// `Histogram`, one level up, so a negative count reaching this class is stored.
    ///
    /// **How the one line works**: `.entry(value)` looks the key up once and hands back a handle to
    /// its slot, present or not; `.or_insert(0)` fills the slot with zero if it was absent and
    /// returns a writable reference either way; the leading `*` writes through that reference. Java
    /// would need a `get`, a null test and a `put`, and would look the key up twice.
    pub fn add_count(&mut self, value: i32, count: i32) {
        *self.value_counts.entry(value).or_insert(0) += count;
    }

    /// `add(CompressedDataList)`.
    pub fn add_all(&mut self, other: &CompressedDataList) {
        for (value, count) in &other.value_counts {
            self.add_count(*value, *count);
        }
    }

    /// The iteration order: each value repeated its count of times, values ascending.
    ///
    /// **What**: the run-length encoding expanded back out. A list holding `{2: 3, 5: 2}` yields
    /// `2, 2, 2, 5, 5`.
    ///
    /// **How**: `.flat_map(...)` turns each `(value, count)` pair into a small sequence and then
    /// flattens all of them into one. `repeat_n(v, n)` is that small sequence.
    ///
    /// **Why `.max(0)`**: a negative count could have been stored, since this class has no guard.
    /// Converting a negative number to an unsigned length would be a program error, so it is
    /// clamped and the entry simply yields nothing. That is the port's choice, not a reading of
    /// the Java iterator, which has no such clamp. No ported caller stores a count below one:
    /// `Histogram` refuses one and `add(val)` adds exactly one.
    ///
    /// `impl Iterator<Item = i32> + '_` is the return type: "some sequence of integers whose
    /// lifetime is tied to this list". The `'_` is what stops the sequence outliving the data it
    /// reads.
    pub fn iter(&self) -> impl Iterator<Item = i32> + '_ {
        self.value_counts
            .iter()
            .flat_map(|(value, count)| std::iter::repeat_n(*value, (*count).max(0) as usize))
    }
}

impl std::fmt::Display for CompressedDataList {
    /// `toString`: `value,count` pairs joined by commas, values ascending.
    ///
    /// ```java
    /// for (Object i: keys){
    ///     if(!str.isEmpty()) {
    ///         str += ",";
    ///     }
    ///     str+=(i+","+valueCounts.get(i));
    /// }
    /// ```
    ///
    /// The Java decides whether to write a comma by testing whether the string is still empty;
    /// the port's `first` flag is the same test, because every entry appends at least one
    /// character. Unlike `Histogram`, an empty list prints the empty string: golden
    /// `crates/gatk-annotation/tests/data/allele_specific_rank_sum.txt.gz`, row `cdl empty`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        for (value, count) in &self.value_counts {
            if !first {
                // `write!` appends to the formatter and can fail if the destination does; the `?`
                // hands that failure to the caller.
                write!(f, ",")?;
            }
            write!(f, "{value},{count}")?;
            first = false;
        }
        Ok(())
    }
}

/// `Histogram`: a `CompressedDataList` over bin indices, plus the bin size that gives them meaning.
#[derive(Debug, Clone, PartialEq)]
pub struct Histogram {
    /// `binSize`: a bin index times this is the value the bin stands for.
    bin_size: f64,
    /// `precisionFormat`: the number of decimals `toString` prints each bin's lower edge with.
    precision: usize,
    /// `dataList`: count per bin index.
    data: CompressedDataList,
}

/// Rust's conventional no-argument constructor, here the same as [`Histogram::new`].
impl Default for Histogram {
    fn default() -> Self {
        Self::new()
    }
}

impl Histogram {
    /// `new Histogram()`: bins of a tenth, printed with one decimal.
    pub fn new() -> Self {
        Histogram {
            bin_size: 0.1,
            precision: 1,
            data: CompressedDataList::new(),
        }
    }

    /// `new Histogram(binSize)`, whose precision is `Math.round(-Math.log10(binSize))` decimals.
    ///
    /// So a bin size that is not a power of a tenth still gets a whole number of decimals. That
    /// number goes negative only above the square root of ten, where minus its log rounds to -1
    /// or below; a bin size of one or two gets zero decimals. At -1 the Java format string is
    /// `%.-1f`, which `String.format` rejects when `toString` runs on a non-empty histogram (an
    /// empty one returns `NaN` before formatting anything).
    ///
    /// **Where the port departs**: it clamps a negative count to zero (`places.max(0)`), so it
    /// prints where the Java would throw. Nothing reaches that: production code builds only
    /// `Histogram::new()`, and the one caller of this constructor is the hundredth-wide histogram
    /// in `crates/gatk-annotation/tests/allele_specific_rank_sum.rs`.
    pub fn with_bin_size(bin_size: f64) -> Self {
        // `jmath::math::round` and `jmath::math::log10` are the port's `Math.round` and
        // `Math.log10`, not Rust's own.
        let places = jmath::math::round(-jmath::math::log10(bin_size));
        Histogram {
            bin_size,
            // Clamped to zero or more first, so the cast to an unsigned count cannot wrap.
            precision: places.max(0) as usize,
            data: CompressedDataList::new(),
        }
    }

    /// `isEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The bin size. The Java field is private with no getter; this accessor is the port's.
    pub fn bin_size(&self) -> f64 {
        self.bin_size
    }

    /// `getBinnedValue`: the epsilon goes in **before** the division, not after.
    ///
    /// **What**: which bin a value falls in, as a signed index that may be negative.
    ///
    /// **How**: shift by a hundredth of a bin, divide by the bin size, take the floor, then round.
    /// The rounding after a floor is redundant arithmetically and is kept because the reference has
    /// it: `Math.round` of an already-integral double is that double, but it also converts to a
    /// `long`, which is what the reference wants.
    ///
    /// **Why the floor and not a truncation**: the floor goes **down** for negative numbers, so a
    /// value of -1.23 with a tenth-wide bin lands in bin -13 and is reported as -1.3, away from
    /// zero. The binning is therefore not symmetric about zero, and the allele-specific rank sums,
    /// which store their Z scores through this, inherit that asymmetry.
    fn binned_value(&self, d: f64) -> i64 {
        jmath::math::round(((d + BIN_EPSILON * self.bin_size) / self.bin_size).floor())
    }

    /// `isValidBinKey`: whether a bin index fits an `int`. It is what makes each later
    /// `key as i32` exact, as it makes the Java's `(int) binKey` exact. `as i64` widens without
    /// loss.
    fn valid_bin_key(binned: i64) -> bool {
        binned <= i32::MAX as i64 && binned >= i32::MIN as i64
    }

    /// `add(Double d)`. A `NaN` is **dropped silently**, not binned and not an error.
    ///
    /// `AS_RankSumTest.outputSingletonValueAsHistogram` depends on that: a NaN Z score leaves its
    /// fresh histogram empty, and the empty histogram prints `NaN` (see the `Display` impl). In the
    /// golden `crates/gatk-annotation/tests/data/allele_specific_rank_sum.txt.gz` the `ref-only`
    /// case, where every read supports the reference, is written `AS_RAW_BaseQRankSum=|NaN`.
    pub fn add(&mut self, d: f64) -> Result<(), HistogramError> {
        if d.is_nan() {
            return Ok(());
        }
        let key = self.binned_value(d);
        if !Self::valid_bin_key(key) {
            return Err(HistogramError::ValueTooExtreme);
        }
        self.data.add(key as i32);
        Ok(())
    }

    /// `add(Double d, int count)`, which unlike the one-argument form does **not** drop a NaN: it
    /// bins it, and `Math.round(Math.floor(NaN))` is zero, so a NaN lands in bin zero.
    pub fn add_count(&mut self, d: f64, count: i32) -> Result<(), HistogramError> {
        if count < 1 {
            return Err(HistogramError::NonPositiveCount);
        }
        let key = self.binned_value(d);
        if !Self::valid_bin_key(key) {
            return Err(HistogramError::ValueTooExtreme);
        }
        self.data.add_count(key as i32, count);
        Ok(())
    }

    /// `add(Histogram h)`: the two count maps summed, refused unless the bin sizes are equal.
    ///
    /// ```java
    /// if (!this.binSize.equals(h.binSize)) {
    /// ```
    ///
    /// `binSize` is a boxed `Double`, so the Java compares with `Double.equals`, not `==`:
    /// `Double.equals` holds two NaNs equal and tells `0.0` from `-0.0`. The port's `!=` is the
    /// IEEE comparison, which does the opposite on both. Neither bin size occurs here (every
    /// histogram is a tenth or, in one test, a hundredth wide), so the plain comparison stands.
    pub fn add_histogram(&mut self, other: &Histogram) -> Result<(), HistogramError> {
        if self.bin_size != other.bin_size {
            return Err(HistogramError::MismatchedBinSize);
        }
        self.data.add_all(&other.data);
        Ok(())
    }

    /// `get(Double d)`: the count in the bin the value falls in, absent rather than zero.
    ///
    /// ```java
    /// return dataList.getValueCounts().get((int) binKey);
    /// ```
    ///
    /// A `Map.get`, so an empty bin is `null`, not `0`. The `bin` rows of the golden
    /// `crates/gatk-annotation/tests/data/allele_specific_rank_sum.txt.gz` print `get(0.0)` in
    /// their last column, and it reads `null` whenever the value landed in another bin.
    pub fn get(&self, d: f64) -> Result<Option<i32>, HistogramError> {
        let key = self.binned_value(d);
        if !Self::valid_bin_key(key) {
            return Err(HistogramError::ValueTooExtreme);
        }
        // The map's `.get` yields a reference to the count or nothing; `.copied()` turns that
        // into the count itself or nothing.
        Ok(self.data.value_counts().get(&(key as i32)).copied())
    }

    /// `median()`, walking the bins in ascending order. `None` for an empty histogram.
    ///
    /// **What**: the middle value, or the average of the two middle values when the count is even.
    ///
    /// **How**: accumulate the counts in key order until the running total reaches the middle
    /// position, then read off the bin.
    ///
    /// **Why the arithmetic looks off by one**: `medianIndex` is `(n + 1) / 2` in **integer**
    /// division, so for an even count it is the **lower** of the two middle positions. The walk
    /// then remembers that bin and averages it with the next one. For an odd count it returns
    /// immediately. The bins it picks are the ones a textbook median picks.
    ///
    /// **Why the average is taken on bin indices**: the reference averages the two indices and
    /// only then scales by the bin size,
    ///
    /// ```java
    /// return (firstMedian+key)/2.0*binSize;
    /// ```
    ///
    /// and that is not always the same double as averaging the two bin values: with tenth-wide
    /// bins, bins -30 and -26 give -2.8000000000000003 this way and -2.8 the other. The port keeps
    /// the reference's order.
    pub fn median(&self) -> Option<f64> {
        // Java's loop over `dataList.valueCounts.values()`, summing.
        let num_items: i32 = self.data.value_counts().values().sum();
        let odd = num_items % 2 != 0;
        // Integer division: for an even count this is the lower of the two middle positions.
        let median_index = (num_items + 1) / 2;

        // `counter` is the running total; `first_median` remembers the lower of the two middle
        // bins once it has been passed, and stays `None` for an odd count.
        let mut counter = 0i32;
        let mut first_median: Option<f64> = None;
        // Ascending key order, which the Java gets from
        // `new TreeSet<>(dataList.valueCounts.keySet())`.
        for (key, count) in self.data.value_counts() {
            counter += *count;
            if counter > median_index {
                // `match` on the `Option` is the Java's `if (firstMedian == null)`: `None` is the
                // null branch, and `Some(first)` names the stored value for the other.
                return Some(match first_median {
                    None => *key as f64 * self.bin_size,
                    Some(first) => (first + *key as f64) / 2.0 * self.bin_size,
                });
            }
            if counter == median_index {
                if odd {
                    return Some(*key as f64 * self.bin_size);
                }
                first_median = Some(*key as f64);
            }
        }
        None
    }
}

impl std::fmt::Display for Histogram {
    /// `toString`: `value,count` pairs, where a bin's value is its index times the bin size (its
    /// lower edge, since the index is a floor: golden row `bin 0.05 0.0,1 1`), printed with the
    /// bin size's precision, and `NaN` for an empty histogram.
    ///
    /// ```java
    /// if (keys.length == 0) {
    ///     return Double.toString(Double.NaN);
    /// }
    /// ```
    ///
    /// The empty case reaches output. `AS_RankSumTest.outputSingletonValueAsHistogram` prints a
    /// fresh histogram right after `h.add(rankSumValue);` with no `isEmpty()` test, and a NaN Z
    /// score leaves it empty. Golden
    /// `crates/gatk-annotation/tests/data/allele_specific_rank_sum.txt.gz`: the row
    /// `hist empty NaN null`, and the three `asraw` rows of the `ref-only` case, such as
    /// `AS_RAW_MQRankSum=|NaN`. `makeCombinedAnnotationString` does test
    /// `if (!alleleValue.isEmpty()) {`, so a combined raw string never carries it.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.data.is_empty() {
            // `Double.toString(Double.NaN)`, not the empty string.
            return write!(f, "NaN");
        }
        let mut first = true;
        for (key, count) in self.data.value_counts() {
            if !first {
                write!(f, ",")?;
            }
            // `String.format(precisionFormat,(double)(int)i*binSize)`: the index widened to a
            // double, times the bin size, in that order. See `format_fixed` for the rounding.
            let centre = *key as f64 * self.bin_size;
            write!(f, "{},{count}", format_fixed(centre, self.precision))?;
            first = false;
        }
        Ok(())
    }
}

/// `String.format("%.Nf", value)` for a bin's value, rounded half-up on the value's **exact**
/// decimal expansion (its first 30 fraction digits).
///
/// That is not the rule [`crate::java_format::format_decimals`] reproduces. That function, in this
/// same crate and the one `gatk-annotation` delegates to, rounds the shortest digits
/// `Double.toString` gives (its module header has the evidence), and the two part wherever the
/// first dropped digit is a 5 in one and not in the other. They agree on every value this
/// function is handed: production code uses only the default tenth-wide bins (one test uses a
/// hundredth), and `k * 0.1` or `k * 0.01` lies within a few ulps of `k / 10` or `k / 100`, so
/// the first dropped digit is a 0 or a 9 under both rules. A bin size of 0.05 printed with one
/// decimal would not agree: bin -763 is -38.15, which this prints as `-38.1` and
/// `format_decimals` as `-38.2`.
///
/// It handles at most 29 places: at 30, `fraction.as_bytes()[places]` is past the end of the
/// 30-digit fraction and the program stops. The tenth and the hundredth need one and two.
fn format_fixed(value: f64, places: usize) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    // `{:.*}` takes its precision from the argument list: the magnitude with 30 fraction digits.
    let text = format!("{:.*}", 30, value.abs());
    // `split_once` cuts at the first point into the two halves, and `expect` stops the program
    // with that message if there is none, which a non-zero precision rules out.
    let (whole, fraction) = text.split_once('.').expect("a decimal point");
    // The whole digits, then the kept fraction digits, as numbers 0 to 9: `b - b'0'` turns an
    // ASCII digit byte into its value.
    let mut digits: Vec<u8> = whole
        .bytes()
        .chain(fraction.bytes().take(places))
        .map(|b| b - b'0')
        .collect();
    // Half-up: the first dropped digit alone decides. Rounding up adds one at the last kept digit
    // and carries left through nines; a carry out of the first digit prepends a 1, so 9.96 at one
    // place becomes 10.0.
    if fraction.as_bytes()[places] >= b'5' {
        let mut index = digits.len();
        loop {
            if index == 0 {
                digits.insert(0, 1);
                break;
            }
            index -= 1;
            if digits[index] == 9 {
                digits[index] = 0;
            } else {
                digits[index] += 1;
                break;
            }
        }
    }
    let split = digits.len() - places;
    let whole: String = digits[..split].iter().map(|d| (d + b'0') as char).collect();
    let fraction: String = digits[split..].iter().map(|d| (d + b'0') as char).collect();
    let sign = if value.is_sign_negative() { "-" } else { "" };
    // `%.0f` prints no point at all.
    if places == 0 {
        return format!("{sign}{whole}");
    }
    format!("{sign}{whole}.{fraction}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_empty_histogram_has_no_median_and_prints_as_nan() {
        let histogram = Histogram::new();
        assert_eq!(histogram.median(), None);
        assert_eq!(histogram.to_string(), "NaN");
    }

    #[test]
    fn the_median_of_an_even_count_averages_the_two_middle_bins() {
        let mut histogram = Histogram::new();
        histogram.add(0.1).expect("a bin");
        histogram.add(0.2).expect("a bin");
        // Bins 1 and 2, averaged to 1.5, times a bin size of a tenth.
        let median = histogram.median().expect("a median");
        assert!((median - 0.15).abs() < 1e-12, "{median}");
    }

    #[test]
    fn a_value_just_below_a_bin_edge_stays_in_the_upper_bin() {
        let mut histogram = Histogram::new();
        // The epsilon is a hundredth of a bin, so -0.0009 still bins to zero.
        histogram.add(-0.0009).expect("a bin");
        assert_eq!(histogram.get(0.0), Ok(Some(1)));
    }

    #[test]
    fn the_rendering_is_the_bin_centre_and_its_count() {
        let mut histogram = Histogram::new();
        histogram.add_count(-1.2, 3).expect("a bin");
        histogram.add(0.5).expect("a bin");
        assert_eq!(histogram.to_string(), "-1.2,3,0.5,1");
    }
}
