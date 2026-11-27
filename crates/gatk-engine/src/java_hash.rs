//! The iteration order of the sample set, which a pileup's element order depends on.
//!
//! **This is not a port, and it must not become one.** The order comes from a `java.util.HashSet`,
//! and `java.util` is GPL2: the OpenJDK Assembly Exception grants permission to *link*, not to
//! translate. htsjdk-rs decision 0013 refused `FloatingDecimal` for that reason and
//! `docs/licence-compatibility-risk.md` records it as the programme's critical risk. The
//! provenance guard enforces it, and it caught the first version of this file, which claimed to be
//! ported from `java.util.HashMap`.
//!
//! # Why the order matters at all
//!
//! `AlignmentContextIteratorBuilder` collects the header's sample names with `Collectors.toSet()`.
//! `LocusIteratorByState` then creates one per-sample manager per element of that set, in iteration
//! order, and concatenates their elements in the same order to build every pileup. So the element
//! order of a multi-sample pileup is that set's iteration order: deterministic, and neither sorted
//! nor the header's. A port using a sorted map, an insertion-ordered vector or Rust's own hasher
//! agrees on single-sample data and diverges as soon as a second sample appears.
//!
//! # What this file stands on instead
//!
//! Two things, and neither is OpenJDK source:
//!
//!  * **`String.hashCode` is specified**, not implementation-defined. Its Javadoc states the value
//!    as `s[0]*31^(n-1) + s[1]*31^(n-2) + ... + s[n-1]`, using `int` arithmetic. Computing that is
//!    implementing a published contract;
//!  * **the bucket layout is not specified**, and the `HashMap` documentation says outright that
//!    iteration order is not guaranteed. So it is treated here as an *observable of the pinned
//!    oracle*: the conformance suite's golden records the order the reference produces for each
//!    probed name set, along with each name's `String.hashCode`, and that golden is this file's
//!    definition rather than a check on it. Where the two disagree, the measurement is right and
//!    this file is wrong.
//!
//! The consequence is a standing obligation rather than a one-off: any sample-name shape the suite
//! does not probe is unverified. The probe therefore includes a set large enough to cross the load
//! factor and a name whose hash is negative, and [`hash_set_order`] refuses outright rather than
//! guessing once a bucket grows past the point where the observed behaviour is known to change.
//!
//! # Reading this file without reading Rust
//!
//! - `i32` is Java's `int`, `u32` an unsigned 32-bit integer Java lacks, and `usize` an unsigned
//!   index type. `wrapping_mul` and `wrapping_add` are Java's `int` `*` and `+`, which wrap
//!   silently on overflow; plain `*` and `+` on an `i32` would panic in a debug build instead.
//! - `as` is a cast, and it does three different things here. `i32 as u32` keeps the 32 bits and
//!   reads them as unsigned, so `(h as u32) >> 16` is Java's `h >>> 16` (Rust's `>>` on an `i32`
//!   is Java's `>>`, and Rust has no `>>>`). `u8 as i8 as i32` reads a byte as signed and then
//!   sign-extends it, which is Java widening a `byte` to `int`. `f32 as usize` truncates toward
//!   zero, which is Java's `(int)` cast of a positive `float`.
//! - `Vec<T>` is an `ArrayList<T>`; a table here is a `Vec<Vec<..>>`, one inner list per bucket.
//!   `&[T]` is a borrowed, read-only view of an array: the caller keeps ownership.
//! - `Result<T, HashOrderError>` is "a value or the refusal below": `Ok(..)` is the return and
//!   `Err(..)` stands for a thrown exception. `Option<T>` is "a value or nothing", Java's `null`
//!   checked by the compiler.
//! - `|x| ...` is a lambda. `.iter().map(..).collect()` is a stream, map and collect;
//!   `.flatten()` is `flatMap` over the inner lists, in order.
//! - A body's last expression, written without `return` or `;`, is its value.

/// `String.hashCode`, over UTF-16 code units.
///
/// Rust strings are UTF-8, so a character outside the basic multilingual plane counts as its two
/// surrogates here, which is what Java sees. A sample name is ASCII in practice, and this is
/// written for the case where it is not.
///
/// The value is fixed by the `String.hashCode` Javadoc (the formula quoted in the module header),
/// so for a string key the hash is specified and only the layout is measured.
pub fn string_hash_code(text: &str) -> i32 {
    let mut hash: i32 = 0;
    // `encode_utf16` yields the string's UTF-16 code units, which are Java's `char`s; the loop is
    // the formula in Horner form, `h = 31 * h + c`, in wrapping `int` arithmetic.
    for unit in text.encode_utf16() {
        hash = hash.wrapping_mul(31).wrapping_add(unit as i32);
    }
    hash
}

/// The mixing step the measured orders are consistent with: the high bits folded into the low
/// ones, which is necessary because only the low bits select a bucket.
pub fn hash_map_hash(text: &str) -> i32 {
    let h = string_hash_code(text);
    // `h ^ (h >>> 16)` in Java's terms; see the module header for the cast.
    h ^ ((h as u32) >> 16) as i32
}

/// What this refuses rather than guessing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HashOrderError {
    /// A bucket reached eight entries, where the observed order stops following the simple
    /// layout this file reproduces. Nothing here has been measured beyond that point, so it
    /// refuses rather than answering confidently and wrongly.
    BucketTreeified { bucket: usize, length: usize },
}

/// The bucket length from which nothing here is measured, and the table size from which
/// `hash_set_order` and `hash_map_order` refuse such a bucket. Below 64 buckets they still answer
/// with the simple layout, a crowded bucket there being answered without being modelled;
/// `JavaHashMap` refuses a crowded bucket at any table size.
const TREEIFY_THRESHOLD: usize = 8;
const MIN_TREEIFY_CAPACITY: usize = 64;

/// The order a `HashSet<String>` built by `Collectors.toSet()` iterates in, as measured.
///
/// The input is the insertion order, which for a stream collector is the stream's order. Duplicates
/// are dropped, keeping the first.
///
/// The set the module header describes is this one, in
/// `org.broadinstitute.hellbender.utils.locusiterator.AlignmentContextIteratorBuilder`:
///
/// ```java
/// final Set<String> samples = header.getReadGroups().stream()
///         .map(SAMReadGroupRecord::getSample)
///         .collect(Collectors.toSet());
/// ```
///
/// Callers include `gatk-tools/src/locus_walker.rs` for that set and
/// `gatk-engine/tests/ac_builder.rs`, which checks it against the golden
/// `tests/data/ac_builder.txt.gz` (its `thirteen` probe crosses the growth point).
///
/// A crowded bucket is refused only once the table has at least `MIN_TREEIFY_CAPACITY` buckets;
/// below that the simple layout is answered.
pub fn hash_set_order(names: &[String]) -> Result<Vec<String>, HashOrderError> {
    // Sixteen buckets to start, which is what the measured orders are consistent with.
    let mut capacity: usize = 16;
    // `vec![x; n]` is a list of `n` copies of `x`: here sixteen empty buckets.
    let mut table: Vec<Vec<String>> = vec![Vec::new(); capacity];
    let mut size: usize = 0;

    for name in names {
        let hash = hash_map_hash(name);
        // The observed index rule: the low bits of the mixed hash select the bucket.
        // `capacity` is a power of two, so `capacity - 1` is a mask; casting the hash to `u32`
        // keeps a negative hash's bits, so the mask behaves like Java's `int` `&`.
        let index = ((capacity - 1) as u32 & hash as u32) as usize;
        // `Set.add` of a name already present changes nothing, so the first occurrence keeps its
        // place.
        if table[index].iter().any(|existing| existing == name) {
            continue;
        }
        table[index].push(name.clone());
        size += 1;
        if table[index].len() >= TREEIFY_THRESHOLD && capacity >= MIN_TREEIFY_CAPACITY {
            return Err(HashOrderError::BucketTreeified {
                bucket: index,
                length: table[index].len(),
            });
        }
        // The observed growth point: three quarters of the current bucket count.
        if size > capacity * 3 / 4 {
            capacity *= 2;
            let mut resized: Vec<Vec<String>> = vec![Vec::new(); capacity];
            // The growth preserves relative order within each bucket, which is what the
            // thirteen-name probe in the golden establishes.
            // `into_iter` consumes the old table, moving each name across rather than copying it.
            for bucket in table.into_iter() {
                for entry in bucket {
                    let h = hash_map_hash(&entry);
                    let to = ((capacity - 1) as u32 & h as u32) as usize;
                    resized[to].push(entry);
                }
            }
            table = resized;
        }
    }

    // Iteration is bucket 0 to the last, each bucket in its own order.
    Ok(table.into_iter().flatten().collect())
}

/// The order a `HashMap` iterates in, given each key's `hashCode` and its insertion order.
///
/// The same layout as [`hash_set_order`], which is the same structure underneath: a `HashSet` is a
/// `HashMap` with a constant value. This form takes the hashes because not every key is a string:
/// `Allele.hashCode` is `Arrays.hashCode(bases) * 31 + Boolean.hashCode(isRef)`, and the order it
/// produces decides which allele wins a tie in a marginalised likelihood matrix.
///
/// That formula is htsjdk 5.0.0's `SimpleAllele.hashCode`:
///
/// ```java
/// return Arrays.hashCode(bases) * 31 + Boolean.hashCode(isRef);
/// ```
///
/// and its caller is `marginalisation_order` in `gatk-annotation/src/depth_per_allele.rs`: the
/// map `DepthPerAlleleBySample` builds with `Collectors.toMap(a -> a, Arrays::asList)` reaches
/// `AlleleLikelihoods.marginalize`, which takes its keys as an array with
/// `newToOldAlleleMap.keySet().toArray(..)`. Others include
/// `gatk-engine/src/fragment.rs` and the copy-ratio segment means in
/// `gatk-tools/src/model_segments_tool.rs`.
///
/// **How**: the same steps as [`hash_set_order`], with `T: Clone + PartialEq` standing for a key
/// that can be copied into the table and compared with `equals`. Each entry carries its hash so
/// that a resize can rebucket it without knowing how to hash a `T`.
pub fn hash_map_order<T: Clone + PartialEq>(
    entries: &[(T, i32)],
) -> Result<Vec<T>, HashOrderError> {
    // The same fold as `hash_map_hash`, as a lambda over an already computed `hashCode`.
    let mix = |hash: i32| hash ^ ((hash as u32) >> 16) as i32;
    let mut capacity: usize = 16;
    let mut table: Vec<Vec<(T, i32)>> = vec![Vec::new(); capacity];
    let mut size: usize = 0;

    // `(key, hash)` takes each pair apart; both are borrowed from `entries`, hence the `*hash`
    // dereference and the `key.clone()` below.
    for (key, hash) in entries {
        let mixed = mix(*hash);
        let index = ((capacity - 1) as u32 & mixed as u32) as usize;
        if table[index].iter().any(|(existing, _)| existing == key) {
            continue;
        }
        table[index].push((key.clone(), *hash));
        size += 1;
        if table[index].len() >= TREEIFY_THRESHOLD && capacity >= MIN_TREEIFY_CAPACITY {
            return Err(HashOrderError::BucketTreeified {
                bucket: index,
                length: table[index].len(),
            });
        }
        if size > capacity * 3 / 4 {
            capacity *= 2;
            let mut resized: Vec<Vec<(T, i32)>> = vec![Vec::new(); capacity];
            for bucket in table.into_iter() {
                for (entry, hash) in bucket {
                    let to = ((capacity - 1) as u32 & mix(hash) as u32) as usize;
                    resized[to].push((entry, hash));
                }
            }
            table = resized;
        }
    }

    Ok(table.into_iter().flatten().map(|(key, _)| key).collect())
}

/// `java.util.Arrays.hashCode(byte[])`.
///
/// The seed is 1, not 0 as in [`string_hash_code`], so an empty array hashes to 1. Callers
/// include the `Allele.hashCode` reproduction `allele_hash_code` in
/// `gatk-annotation/src/depth_per_allele.rs`.
pub fn byte_array_hash_code(bytes: &[u8]) -> i32 {
    let mut hash: i32 = 1;
    for byte in bytes {
        // The element is widened to int, so a byte above 0x7f is negative.
        hash = hash.wrapping_mul(31).wrapping_add(*byte as i8 as i32);
    }
    hash
}

/// `String.compareTo`, which is **UTF-16 code-unit** order, not byte order and not Unicode
/// scalar order.
///
/// The two orders agree on ASCII and disagree above the BMP: a supplementary character encodes as
/// a surrogate pair whose first unit is in `0xD800..=0xDBFF`, so it sorts *before* every character
/// in `0xE000..=0xFFFF` even though its scalar value is larger. Read names are ASCII in practice,
/// but the comparator this feeds decides the order of an assembly region's reads, and a sort order
/// that is right in practice is exactly the kind of thing that is wrong once.
pub fn compare_strings(left: &str, right: &str) -> std::cmp::Ordering {
    let mut left_units = left.encode_utf16();
    let mut right_units = right.encode_utf16();
    loop {
        // `next()` gives `Some(unit)` or `None` at the end. Matching the pair covers all four
        // cases, and the compiler rejects a `match` that misses one. `if a != b` is a guard: the
        // arm applies only when it holds, and equal units fall through to `_`, the next unit.
        match (left_units.next(), right_units.next()) {
            (None, None) => return std::cmp::Ordering::Equal,
            // Java returns len1 - len2 when one is a prefix of the other.
            (None, Some(_)) => return std::cmp::Ordering::Less,
            (Some(_), None) => return std::cmp::Ordering::Greater,
            (Some(a), Some(b)) if a != b => return a.cmp(&b),
            _ => {}
        }
    }
}

/// A key whose Java `hashCode` is specified, which is what every map modelled here is keyed by.
///
/// A `trait` is a Java interface; `impl JavaHashCode for X` below is `X implements JavaHashCode`,
/// and because the trait is defined here it can be implemented for `i32` and `String`. Other crates
/// implement it for their own keys, for example `LibraryName` in `gatk-tools/src/read_metadata.rs`.
pub trait JavaHashCode {
    fn java_hash_code(&self) -> i32;
}

/// `Integer.hashCode`, which its Javadoc gives as the value itself.
impl JavaHashCode for i32 {
    fn java_hash_code(&self) -> i32 {
        *self
    }
}

/// `String.hashCode`, over UTF-16 code units.
impl JavaHashCode for String {
    fn java_hash_code(&self) -> i32 {
        string_hash_code(self)
    }
}

/// A `HashMap` whose iteration order follows the measured layout, including its capacity.
///
/// [`hash_map_order`] answers for a map built once with the default capacity and never shrunk.
/// `PathSeqBuildReferenceTaxonomy` does neither: it sizes maps at construction
/// (`new HashMap<>(n)`, `new HashSet<>(collection)`), removes from them (`retainAll`,
/// `removeChild`), and iterates them again afterwards. Each of those decides how many buckets the
/// table has when it is iterated, and the bucket count decides the order. So this keeps the table
/// itself, and every rule below is what the `pathseq-taxonomy-kryo` probes are consistent with
/// rather than a reading of any implementation:
///
///  * **the table exists from the first insertion**, sized by what the constructor was asked for:
///    sixteen by default, and the smallest power of two not below a requested capacity;
///  * **it doubles once the entry count exceeds three quarters of the bucket count**, and a
///    doubling keeps the relative order of entries that share a bucket;
///  * **a removal never shrinks it**, so a map that once held many entries iterates what is left
///    in the order the larger table gives;
///  * **replacing a key's value keeps its place**;
///  * **the bucket is the low bits of the hash with its high half folded in**, and within a bucket
///    the order is the order of insertion.
///
/// Iteration is therefore the insertion order, stably sorted by bucket at the final size. The one
/// state nothing here reproduces is a crowded bucket: from `TREEIFY_THRESHOLD` entries on the
/// observed layout changes shape, so the map records that it got there and [`Self::check`]
/// refuses, rather than answering an order that was never measured.
#[derive(Debug, Clone)]
pub struct JavaHashMap<K, V> {
    /// The bucket count the first insertion allocates.
    initial_capacity: usize,
    /// One list of entries per bucket, each in iteration order. Empty until the first insertion,
    /// which is how "the table exists from the first insertion" is kept.
    buckets: Vec<Vec<(K, V)>>,
    /// The entry count, Java's `size()`.
    len: usize,
    /// The first bucket that reached `TREEIFY_THRESHOLD` entries, kept so [`Self::check`] can
    /// refuse even after removals shrank it. Unlike [`hash_set_order`], this is recorded at any
    /// table size.
    crowded: Option<HashOrderError>,
}

/// `JavaHashMap::default()` is `new()`. `Default` is the trait generic code uses for "a fresh
/// empty value", a no-argument constructor in Java's terms.
impl<K: JavaHashCode + PartialEq + Clone, V> Default for JavaHashMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: JavaHashCode + PartialEq + Clone, V> JavaHashMap<K, V> {
    /// `new HashMap<>()`: sixteen buckets.
    pub fn new() -> Self {
        Self::allocating(16)
    }

    /// `new HashMap<>(capacity)`: the smallest power of two not below the request.
    pub fn with_capacity(capacity: usize) -> Self {
        Self::allocating(capacity.max(1).next_power_of_two())
    }

    /// The capacity `new HashSet<>(collection)` asks for: the collection's size over the load
    /// factor, in `float` arithmetic, plus one, and never below sixteen.
    ///
    /// `f32` is Java's `float`, so the division is done at `float` precision, not `double`. The
    /// `as usize` truncates.
    pub fn copy_capacity(size: usize) -> usize {
        ((size as f32 / 0.75f32) as usize + 1).max(16)
    }

    /// An empty map that will allocate `initial_capacity` buckets on its first insertion.
    fn allocating(initial_capacity: usize) -> Self {
        JavaHashMap {
            initial_capacity,
            buckets: Vec::new(),
            len: 0,
            crowded: None,
        }
    }

    /// The bucket a `hashCode` lands in at `capacity` buckets: the fold of [`hash_map_hash`], then
    /// the low bits.
    fn bucket_of(hash: i32, capacity: usize) -> usize {
        let mixed = hash ^ ((hash as u32) >> 16) as i32;
        // Through `u32` first, so a negative hash widens without sign extension.
        (mixed as u32 as usize) & (capacity - 1)
    }

    /// The bucket count the table has now, or would have at its first insertion.
    pub fn capacity(&self) -> usize {
        if self.buckets.is_empty() {
            self.initial_capacity
        } else {
            self.buckets.len()
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Where `key` is: its bucket and its index within the bucket, or `None` if absent.
    fn find(&self, key: &K) -> Option<(usize, usize)> {
        if self.buckets.is_empty() {
            return None;
        }
        let bucket = Self::bucket_of(key.java_hash_code(), self.buckets.len());
        // `position` is `indexOf` with a predicate; `.map` on the `Option` pairs the index with
        // its bucket and leaves `None` as it is.
        self.buckets[bucket]
            .iter()
            .position(|(existing, _)| existing == key)
            .map(|index| (bucket, index))
    }

    pub fn contains_key(&self, key: &K) -> bool {
        self.find(key).is_some()
    }

    pub fn get(&self, key: &K) -> Option<&V> {
        self.find(key)
            .map(|(bucket, index)| &self.buckets[bucket][index].1)
    }

    /// `get`, returning a writable borrow (`&mut V`) so the caller can update the value in place,
    /// as `gatk-tools/src/pathseq_score.rs` does to increment a count.
    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        self.find(key)
            .map(|(bucket, index)| &mut self.buckets[bucket][index].1)
    }

    /// `computeIfAbsent` on a key that is absent: unlike `put`, the node goes to the HEAD of its
    /// bucket, and the table grows on the way in rather than on the way out.
    ///
    /// The caller is `gatk-tools/src/pathseq_score.rs`, for this expression in GATK 4.7.0.0's
    /// `PSScorer`:
    ///
    /// ```java
    /// taxIds.collect(Collectors.groupingBy(e -> e, Collectors.counting()))
    /// ```
    pub fn insert_at_head(&mut self, key: K, value: V) {
        if self.buckets.is_empty() {
            // `(0..n).map(..).collect()` builds `n` empty buckets.
            self.buckets = (0..self.initial_capacity).map(|_| Vec::new()).collect();
        }
        // Tested against the count before this entry, so the growth comes first.
        if self.len > self.buckets.len() * 3 / 4 {
            self.grow();
        }
        let bucket = Self::bucket_of(key.java_hash_code(), self.buckets.len());
        // `Vec::insert(0, ..)` is `List.add(0, ..)`.
        self.buckets[bucket].insert(0, (key, value));
        self.len += 1;
        let length = self.buckets[bucket].len();
        if length >= TREEIFY_THRESHOLD && self.crowded.is_none() {
            self.crowded = Some(HashOrderError::BucketTreeified { bucket, length });
        }
    }

    /// `put`: a new key goes to the end of its bucket, an existing one keeps its place.
    ///
    /// Returns the previous value, as `put` does, with `None` for Java's `null`.
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        // `if let Some(x) = ..` runs the block only when there is a value, binding it to `x`.
        // `std::mem::replace` stores the new value and hands back the old one.
        if let Some((bucket, index)) = self.find(&key) {
            return Some(std::mem::replace(&mut self.buckets[bucket][index].1, value));
        }
        if self.buckets.is_empty() {
            self.buckets = (0..self.initial_capacity).map(|_| Vec::new()).collect();
        }
        let capacity = self.buckets.len();
        let bucket = Self::bucket_of(key.java_hash_code(), capacity);
        self.buckets[bucket].push((key, value));
        self.len += 1;
        let length = self.buckets[bucket].len();
        if length >= TREEIFY_THRESHOLD && self.crowded.is_none() {
            self.crowded = Some(HashOrderError::BucketTreeified { bucket, length });
        }
        // Tested after the entry is in, against the bucket count it went into: the opposite order
        // to `insert_at_head`.
        if self.len > capacity * 3 / 4 {
            self.grow();
        }
        None
    }

    /// The doubling, which keeps each bucket's relative order.
    fn grow(&mut self) {
        let capacity = self.buckets.len() * 2;
        let mut resized: Vec<Vec<(K, V)>> = (0..capacity).map(|_| Vec::new()).collect();
        // `std::mem::take` moves the old buckets out and leaves an empty list behind, so the
        // entries can be moved into `resized` while `self` is still borrowed. Walking the old
        // buckets in order and appending is what keeps each bucket's relative order.
        for bucket in std::mem::take(&mut self.buckets) {
            for (key, value) in bucket {
                let to = Self::bucket_of(key.java_hash_code(), capacity);
                resized[to].push((key, value));
            }
        }
        self.buckets = resized;
    }

    /// `remove`, which leaves the table its size.
    pub fn remove(&mut self, key: &K) -> Option<V> {
        // `?` on an `Option` returns `None` from this function at once when the key is absent.
        let (bucket, index) = self.find(key)?;
        self.len -= 1;
        Some(self.buckets[bucket].remove(index).1)
    }

    /// `keySet().retainAll(...)`, and anything else that removes by a predicate.
    ///
    /// `keep` is a lambda (`impl FnMut` accepts any closure, including one that updates state it
    /// captured); `Vec::retain` keeps the survivors in their order, and the table keeps its size.
    pub fn retain(&mut self, mut keep: impl FnMut(&K, &V) -> bool) {
        for bucket in &mut self.buckets {
            bucket.retain(|(key, value)| keep(key, value));
        }
        self.len = self.buckets.iter().map(Vec::len).sum();
    }

    /// The entries in iteration order.
    ///
    /// `impl Iterator` returns a lazy iterator without naming its type, like returning an
    /// `Iterator<Map.Entry<K, V>>`; it borrows the map, so the map cannot change while it is used.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        // The final `map` turns a borrowed pair `&(K, V)` into a pair of borrows `(&K, &V)`.
        self.buckets
            .iter()
            .flatten()
            .map(|(key, value)| (key, value))
    }

    /// The keys in iteration order.
    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.iter().map(|(key, _)| key)
    }

    /// Whether the order this answers was ever measured: a bucket that reached
    /// `TREEIFY_THRESHOLD` entries is refused, even if it later shrank.
    pub fn check(&self) -> Result<(), HashOrderError> {
        match &self.crowded {
            Some(error) => Err(error.clone()),
            None => Ok(()),
        }
    }
}

impl<K: JavaHashCode + PartialEq + Clone> JavaHashMap<K, ()> {
    /// `new HashSet<>(collection)`: sized by [`Self::copy_capacity`], filled in the source's
    /// iteration order.
    ///
    /// `JavaHashMap<K, ()>` is a set: `()` is an empty value, the constant a `HashSet` stores.
    /// `ExactSizeIterator` is an iterator that knows its length up front, standing in for
    /// `Collection.size()`; `'a` says the borrowed keys outlive the call, and each is cloned in.
    pub fn copy_of<'a>(keys: impl ExactSizeIterator<Item = &'a K>) -> Self
    where
        K: 'a,
    {
        let mut set = Self::with_capacity(Self::copy_capacity(keys.len()));
        for key in keys {
            set.insert(key.clone(), ());
        }
        set
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn order<V>(map: &JavaHashMap<i32, V>) -> Vec<i32> {
        map.keys().copied().collect()
    }

    #[test]
    fn a_default_map_orders_small_integers_by_value() {
        let mut map = JavaHashMap::new();
        for key in [9, 3, 15, 1] {
            map.insert(key, ());
        }
        assert_eq!(order(&map), vec![1, 3, 9, 15]);
    }

    #[test]
    fn two_keys_in_one_bucket_keep_their_insertion_order() {
        let mut map = JavaHashMap::new();
        for key in [33, 1, 17] {
            map.insert(key, ());
        }
        assert_eq!(order(&map), vec![33, 1, 17]);
    }

    #[test]
    fn the_thirteenth_entry_doubles_the_table() {
        let mut map = JavaHashMap::new();
        for key in 0..11 {
            map.insert(key, ());
        }
        map.insert(16, ());
        assert_eq!(map.capacity(), 16);
        map.insert(17, ());
        assert_eq!(map.capacity(), 32);
        // At thirty-two buckets sixteen is no longer in zero's bucket.
        assert_eq!(order(&map)[..2], [0, 1]);
    }

    #[test]
    fn a_removal_leaves_the_larger_table() {
        let mut map = JavaHashMap::new();
        for key in 0..13 {
            map.insert(key, ());
        }
        map.insert(32, ());
        map.retain(|key, _| *key == 32 || *key == 0 || *key == 1);
        // At thirty-two buckets 32 shares bucket 0 with 0, after it; at sixteen it would too, so
        // the capacity assertion is what shows the table did not shrink.
        assert_eq!(map.capacity(), 32);
        assert_eq!(order(&map), vec![0, 32, 1]);
    }

    #[test]
    fn a_requested_capacity_rounds_up_to_a_power_of_two() {
        assert_eq!(JavaHashMap::<i32, ()>::with_capacity(0).capacity(), 1);
        assert_eq!(JavaHashMap::<i32, ()>::with_capacity(2).capacity(), 2);
        assert_eq!(JavaHashMap::<i32, ()>::with_capacity(3).capacity(), 4);
        assert_eq!(JavaHashMap::<i32, ()>::with_capacity(17).capacity(), 32);
        assert_eq!(JavaHashMap::<i32, ()>::copy_capacity(12), 17);
        assert_eq!(JavaHashMap::<i32, ()>::copy_capacity(3), 16);
    }

    #[test]
    fn replacing_a_value_keeps_the_key_where_it_was() {
        let mut map = JavaHashMap::new();
        map.insert(33, 'a');
        map.insert(1, 'b');
        assert_eq!(map.insert(33, 'c'), Some('a'));
        assert_eq!(
            map.iter().collect::<Vec<_>>(),
            vec![(&33, &'c'), (&1, &'b')]
        );
    }

    #[test]
    fn a_crowded_bucket_is_refused() {
        let mut map = JavaHashMap::with_capacity(64);
        for index in 0..TREEIFY_THRESHOLD as i32 {
            map.insert(index * 64, ());
        }
        assert!(map.check().is_err());
    }
}
