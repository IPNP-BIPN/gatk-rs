//! Ported from `java.util.HashMap` and `java.lang.String.hashCode` (JDK 17), which GATK's oracle
//! pins.
//!
//! This exists because an iteration order leaks into a pileup. `AlignmentContextIteratorBuilder`
//! collects the header's sample names with `Collectors.toSet()`, which is a `HashSet`, and hands
//! that set to `LocusIteratorByState`. The per-sample managers are then created in *that* order,
//! and `LocusIteratorByState` concatenates their elements in the same order to build each pileup.
//!
//! So the order of elements in every pileup is the bucket order of a Java `HashSet` over the sample
//! names. It is deterministic, and it is neither sorted nor the header's order. A port that used a
//! `BTreeMap`, an insertion-ordered vector, or Rust's `HashSet` would agree on single-sample data
//! and diverge as soon as a second sample appeared.
//!
//! What is reproduced here is exactly the part that decides that order:
//!
//!  * `String.hashCode`: `s[0]*31^(n-1) + s[1]*31^(n-2) + ... + s[n-1]`, over UTF-16 code units,
//!    wrapping on overflow;
//!  * `HashMap.hash`: `h ^ (h >>> 16)`, which mixes the high bits down because the index only uses
//!    the low ones;
//!  * the table: 16 buckets, doubling when `size` exceeds `capacity * 0.75`, index
//!    `(capacity - 1) & hash`;
//!  * the split on resize, which Java 8 and later made order-preserving: a bucket's entries are
//!    partitioned into a low list and a high list, each keeping its relative order.
//!
//! Only the parts an iteration order depends on are here. Treeification, which a bucket reaches at
//! eight collisions, is not: it would change the order within one bucket, and
//! [`hash_set_order`] refuses rather than guessing when a bucket gets that long.

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

/// `HashMap.hash(Object)`: spread the high bits into the low ones.
pub fn hash_map_hash(text: &str) -> i32 {
    let h = string_hash_code(text);
    h ^ ((h as u32) >> 16) as i32
}

/// What this port refuses rather than guessing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HashOrderError {
    /// A bucket reached `TREEIFY_THRESHOLD`, where the bucket becomes a red-black tree ordered by
    /// hash and then by comparable key. Reproducing that is a separate port, and answering without
    /// it would answer confidently and wrongly.
    BucketTreeified { bucket: usize, length: usize },
}

/// `TREEIFY_THRESHOLD`. A bucket is treeified when a *ninth* entry is added to it, and only when
/// the table has at least `MIN_TREEIFY_CAPACITY` (64) buckets; below that it resizes instead.
const TREEIFY_THRESHOLD: usize = 8;
const MIN_TREEIFY_CAPACITY: usize = 64;

/// The order a `HashSet<String>` built by `Collectors.toSet()` iterates in.
///
/// The input is the insertion order, which for a stream collector is the stream's order. Duplicates
/// are dropped, keeping the first, as `add` does.
pub fn hash_set_order(names: &[String]) -> Result<Vec<String>, HashOrderError> {
    // `Collectors.toSet()` uses `new HashSet<>()`, whose table is created lazily at 16.
    let mut capacity: usize = 16;
    let mut table: Vec<Vec<String>> = vec![Vec::new(); capacity];
    let mut size: usize = 0;

    for name in names {
        let hash = hash_map_hash(name);
        // `(n - 1) & hash` on the current table.
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
        // `++size > threshold` with `threshold = capacity * 0.75`.
        if size > capacity * 3 / 4 {
            capacity *= 2;
            let mut resized: Vec<Vec<String>> = vec![Vec::new(); capacity];
            // Java 8's order-preserving split: each old bucket's entries go to `index` or
            // `index + oldCapacity`, keeping their relative order in both.
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
