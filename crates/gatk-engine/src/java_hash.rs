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

/// `String.hashCode`, over UTF-16 code units.
///
/// Rust strings are UTF-8, so a character outside the basic multilingual plane counts as its two
/// surrogates here, which is what Java sees. A sample name is ASCII in practice, and this is
/// written for the case where it is not.
pub fn string_hash_code(text: &str) -> i32 {
    let mut hash: i32 = 0;
    for unit in text.encode_utf16() {
        hash = hash.wrapping_mul(31).wrapping_add(unit as i32);
    }
    hash
}

/// The mixing step the measured orders are consistent with: the high bits folded into the low
/// ones, which is necessary because only the low bits select a bucket.
pub fn hash_map_hash(text: &str) -> i32 {
    let h = string_hash_code(text);
    h ^ ((h as u32) >> 16) as i32
}

/// What this refuses rather than guessing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HashOrderError {
    /// A bucket grew past eight entries, where the observed order stops following the simple
    /// layout this file reproduces. Nothing here has been measured beyond that point, so it
    /// refuses rather than answering confidently and wrongly.
    BucketTreeified { bucket: usize, length: usize },
}

/// The bucket length past which nothing here is measured, and the table size below which the
/// structure grows instead of changing shape.
const TREEIFY_THRESHOLD: usize = 8;
const MIN_TREEIFY_CAPACITY: usize = 64;

/// The order a `HashSet<String>` built by `Collectors.toSet()` iterates in, as measured.
///
/// The input is the insertion order, which for a stream collector is the stream's order. Duplicates
/// are dropped, keeping the first.
pub fn hash_set_order(names: &[String]) -> Result<Vec<String>, HashOrderError> {
    // Sixteen buckets to start, which is what the measured orders are consistent with.
    let mut capacity: usize = 16;
    let mut table: Vec<Vec<String>> = vec![Vec::new(); capacity];
    let mut size: usize = 0;

    for name in names {
        let hash = hash_map_hash(name);
        // The observed index rule: the low bits of the mixed hash select the bucket.
        let index = ((capacity - 1) as u32 & hash as u32) as usize;
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

    Ok(table.into_iter().flatten().collect())
}

/// The order a `HashMap` iterates in, given each key's `hashCode` and its insertion order.
///
/// The same layout as [`hash_set_order`], which is the same structure underneath: a `HashSet` is a
/// `HashMap` with a constant value. This form takes the hashes because not every key is a string:
/// `Allele.hashCode` is `Arrays.hashCode(bases) * 31 + Boolean.hashCode(isRef)`, and the order it
/// produces decides which allele wins a tie in a marginalised likelihood matrix.
pub fn hash_map_order<T: Clone + PartialEq>(
    entries: &[(T, i32)],
) -> Result<Vec<T>, HashOrderError> {
    let mix = |hash: i32| hash ^ ((hash as u32) >> 16) as i32;
    let mut capacity: usize = 16;
    let mut table: Vec<Vec<(T, i32)>> = vec![Vec::new(); capacity];
    let mut size: usize = 0;

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
/// state nothing here reproduces is a crowded bucket: past [`TREEIFY_THRESHOLD`] entries the
/// observed layout changes shape, so the map records that it got there and [`Self::check`]
/// refuses, rather than answering an order that was never measured.
#[derive(Debug, Clone)]
pub struct JavaHashMap<K, V> {
    /// The bucket count the first insertion allocates.
    initial_capacity: usize,
    buckets: Vec<Vec<(K, V)>>,
    len: usize,
    crowded: Option<HashOrderError>,
}

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
    pub fn copy_capacity(size: usize) -> usize {
        ((size as f32 / 0.75f32) as usize + 1).max(16)
    }

    fn allocating(initial_capacity: usize) -> Self {
        JavaHashMap {
            initial_capacity,
            buckets: Vec::new(),
            len: 0,
            crowded: None,
        }
    }

    fn bucket_of(hash: i32, capacity: usize) -> usize {
        let mixed = hash ^ ((hash as u32) >> 16) as i32;
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

    fn find(&self, key: &K) -> Option<(usize, usize)> {
        if self.buckets.is_empty() {
            return None;
        }
        let bucket = Self::bucket_of(key.java_hash_code(), self.buckets.len());
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

    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        self.find(key)
            .map(|(bucket, index)| &mut self.buckets[bucket][index].1)
    }

    /// `put`: a new key goes to the end of its bucket, an existing one keeps its place.
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
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
        if self.len > capacity * 3 / 4 {
            self.grow();
        }
        None
    }

    /// The doubling, which keeps each bucket's relative order.
    fn grow(&mut self) {
        let capacity = self.buckets.len() * 2;
        let mut resized: Vec<Vec<(K, V)>> = (0..capacity).map(|_| Vec::new()).collect();
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
        let (bucket, index) = self.find(key)?;
        self.len -= 1;
        Some(self.buckets[bucket].remove(index).1)
    }

    /// `keySet().retainAll(...)`, and anything else that removes by a predicate.
    pub fn retain(&mut self, mut keep: impl FnMut(&K, &V) -> bool) {
        for bucket in &mut self.buckets {
            bucket.retain(|(key, value)| keep(key, value));
        }
        self.len = self.buckets.iter().map(Vec::len).sum();
    }

    /// The entries in iteration order.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
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
    /// [`TREEIFY_THRESHOLD`] entries is refused, even if it later shrank.
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
