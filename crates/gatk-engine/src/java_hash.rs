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
