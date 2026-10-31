//! `HopscotchCollection` (GATK 4.6.2.0), whose iteration order is what a tool writes.
//!
//! The Spark SV and PathSeq tools collect k-mers through this table and print them in the order it
//! gives them, so its layout is part of the output: bucket by bucket from the first chain head,
//! each chain in link order, and where an entry lands depends on everything inserted before it.
//!
//! # The index survives the resize that invalidates it
//!
//! ```java
//! final int bucketIndex = hashToIndex(toKey().apply(entry));
//! try {
//!     return insert(entry, bucketIndex, collision);
//! } catch ( final IllegalStateException ise ) {
//!     resize();
//!     return insert(entry, bucketIndex, collision);
//! }
//! ```
//!
//! A failed hopscotch resizes and retries with the bucket computed for the OLD capacity, so the
//! entry heads a chain it does not hash to and `find` will not see it. A failed attempt also leaves
//! the moves it made behind. Both are kept here, since both change what is iterated.
//!
//! # A status byte is a Java byte
//!
//! The high bit marks a chain head and the low seven the offset to the next link, and `+=` on a
//! `byte` truncates. Hence `i8` and wrapping arithmetic throughout.

/// `SetSizeUtils.legalSizes`.
const LEGAL_SIZES: [i32; 47] = [
    251, 359, 509, 719, 1021, 1447, 2039, 2887, 4093, 5791, 8191, 11579, 16381, 23167, 32749,
    46337, 65521, 92681, 131071, 185363, 262139, 370723, 524287, 741431, 1048573, 1482907, 2097143,
    2965819, 4194301, 5931641, 8388593, 11863279, 16777213, 23726561, 33554393, 47453111, 67108859,
    94906249, 134217689, 189812507, 268435399, 379625047, 536870909, 759250111, 1073741789,
    1518500213, 2147483629,
];

const LOAD_FACTOR: f64 = 0.85;
const SPREADER: i32 = 241;

/// An entry's `toKey().apply(entry).hashCode()` and `entryCollides()`.
pub trait HopscotchEntry: Clone {
    fn key_hash(&self) -> i32;
    fn collides(&self, other: &Self) -> bool;
}

/// `IllegalStateException` out of a resize that could not place what it already held.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HopscotchFailed(pub String);

pub struct Hopscotch<T: HopscotchEntry> {
    capacity: i32,
    size: i32,
    buckets: Vec<Option<T>>,
    status: Vec<i8>,
}

impl<T: HopscotchEntry> Hopscotch<T> {
    /// `new HopscotchCollection(capacity)`, through `computeCapacity`.
    pub fn new(capacity: i32) -> Hopscotch<T> {
        let capacity = compute_capacity(capacity);
        Hopscotch {
            capacity,
            size: 0,
            buckets: vec![None; capacity as usize],
            status: vec![0; capacity as usize],
        }
    }

    pub fn len(&self) -> usize {
        self.size as usize
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    pub fn capacity(&self) -> i32 {
        self.capacity
    }

    /// `add`: false when an entry that collides was already there.
    pub fn add(&mut self, entry: T) -> Result<bool, HopscotchFailed> {
        if self.size == self.capacity {
            self.resize()?;
        }
        let bucket = self.hash_to_index(entry.key_hash());
        match self.insert(&entry, bucket, true) {
            Ok(added) => Ok(added),
            Err(_) => {
                self.resize()?;
                self.insert(&entry, bucket, true)
            }
        }
    }

    /// `find`: the chain the key hashes to, and no other.
    pub fn find_mut(&mut self, hash: i32, matches: impl Fn(&T) -> bool) -> Option<&mut T> {
        let mut bucket = self.hash_to_index(hash);
        if !self.is_chain_head(bucket) {
            return None;
        }
        loop {
            if matches(self.buckets[bucket as usize].as_ref()?) {
                return self.buckets[bucket as usize].as_mut();
            }
            let offset = self.offset(bucket);
            if offset == 0 {
                return None;
            }
            bucket = self.index(bucket, offset);
        }
    }

    /// `iterator`: every chain from its head, heads in bucket order.
    pub fn iter(&self) -> impl Iterator<Item = &T> + '_ {
        let mut order = Vec::with_capacity(self.size as usize);
        for head in 0..self.capacity {
            if !self.is_chain_head(head) {
                continue;
            }
            let mut bucket = head;
            loop {
                if let Some(entry) = &self.buckets[bucket as usize] {
                    order.push(entry);
                }
                let offset = self.offset(bucket);
                if offset == 0 {
                    break;
                }
                bucket = self.index(bucket, offset);
            }
        }
        order.into_iter()
    }

    fn hash_to_index(&self, hash: i32) -> i32 {
        let result = SPREADER.wrapping_mul(hash) % self.capacity;
        if result < 0 {
            result + self.capacity
        } else {
            result
        }
    }

    fn insert(&mut self, entry: &T, bucket: i32, test: bool) -> Result<bool, HopscotchFailed> {
        if self.buckets[bucket as usize].is_some() && !self.is_chain_head(bucket) {
            self.evict(bucket)?;
        }
        if self.buckets[bucket as usize].is_none() {
            self.buckets[bucket as usize] = Some(entry.clone());
            self.status[bucket as usize] = i8::MIN;
            self.size += 1;
            return Ok(true);
        }
        let mut end_of_chain = bucket;
        loop {
            if test
                && self.buckets[end_of_chain as usize]
                    .as_ref()
                    .is_some_and(|held| held.collides(entry))
            {
                return Ok(false);
            }
            let offset = self.offset(end_of_chain);
            if offset == 0 {
                break;
            }
            end_of_chain = self.index(end_of_chain, offset);
        }
        let empty = self.insert_into_chain(bucket, end_of_chain)?;
        self.buckets[empty as usize] = Some(entry.clone());
        self.size += 1;
        Ok(true)
    }

    fn insert_into_chain(
        &mut self,
        bucket: i32,
        end_of_chain: i32,
    ) -> Result<i32, HopscotchFailed> {
        let offset_to_end = self.index_diff(bucket, end_of_chain);
        let mut empty = self.find_empty_bucket(bucket);
        let max_offset = offset_to_end + i8::MAX as i32;
        let mut offset_to_empty;
        loop {
            offset_to_empty = self.index_diff(bucket, empty);
            if offset_to_empty <= max_offset {
                break;
            }
            empty = self.hopscotch(bucket, empty)?;
        }
        if offset_to_empty > offset_to_end {
            self.add_status(end_of_chain, offset_to_empty - offset_to_end);
        } else {
            self.link_into_chain(bucket, empty);
        }
        Ok(empty)
    }

    fn link_into_chain(&mut self, bucket: i32, empty: i32) {
        let mut offset_to_empty = self.index_diff(bucket, empty);
        let mut at = bucket;
        let mut offset;
        loop {
            offset = self.offset(at);
            if offset >= offset_to_empty {
                break;
            }
            at = self.index(at, offset);
            offset_to_empty -= offset;
        }
        offset -= offset_to_empty;
        self.add_status(at, -offset);
        self.status[empty as usize] = offset as i8;
    }

    fn evict(&mut self, to_evict: i32) -> Result<(), HopscotchFailed> {
        let hash = self.buckets[to_evict as usize]
            .as_ref()
            .expect("a squatter")
            .key_hash();
        let bucket = self.hash_to_index(hash);
        let offset_to_evictee = self.index_diff(bucket, to_evict);
        let mut empty = self.find_empty_bucket(bucket);
        let mut from = bucket;
        loop {
            while self.index_diff(bucket, empty) > offset_to_evictee {
                empty = self.hopscotch(from, empty)?;
            }
            if empty == to_evict {
                return Ok(());
            }
            from = empty;
            self.link_into_chain(bucket, empty);
            let mut previous = bucket;
            let mut offset_to_next = self.offset(previous);
            let mut next = self.index(previous, offset_to_next);
            loop {
                offset_to_next = self.offset(next);
                if offset_to_next == 0 {
                    break;
                }
                previous = next;
                next = self.index(next, offset_to_next);
            }
            self.buckets[empty as usize] = self.buckets[next as usize].take();
            self.status[next as usize] = 0;
            let offset = self.offset(previous);
            self.add_status(previous, -offset);
            empty = next;
        }
    }

    fn find_empty_bucket(&self, mut bucket: i32) -> i32 {
        loop {
            bucket = self.index(bucket, 1);
            if self.buckets[bucket as usize].is_none() {
                return bucket;
            }
        }
    }

    fn is_chain_head(&self, bucket: i32) -> bool {
        self.status[bucket as usize] < 0
    }

    fn offset(&self, bucket: i32) -> i32 {
        (self.status[bucket as usize] & i8::MAX) as i32
    }

    /// `status[bucket] += delta`, a compound assignment on a byte.
    fn add_status(&mut self, bucket: i32, delta: i32) {
        let status = &mut self.status[bucket as usize];
        *status = (*status as i32).wrapping_add(delta) as i8;
    }

    fn index(&self, bucket: i32, offset: i32) -> i32 {
        let result = bucket + offset;
        if result >= self.capacity {
            result - self.capacity
        } else if result < 0 {
            result + self.capacity
        } else {
            result
        }
    }

    fn index_diff(&self, first: i32, second: i32) -> i32 {
        let result = second - first;
        if result < 0 {
            result + self.capacity
        } else {
            result
        }
    }

    fn hopscotch(&mut self, from: i32, empty: i32) -> Result<i32, HopscotchFailed> {
        let from_to_empty = self.index_diff(from, empty);
        let mut offset_to_empty = i8::MAX as i32;
        while offset_to_empty > 1 {
            let bucket = self.index(empty, -offset_to_empty);
            let in_bucket = self.offset(bucket);
            if in_bucket != 0
                && in_bucket < offset_to_empty
                && offset_to_empty - in_bucket < from_to_empty
            {
                let to_move = self.index(bucket, in_bucket);
                self.relocate(bucket, to_move, empty);
                return Ok(to_move);
            }
            offset_to_empty -= 1;
        }
        Err(HopscotchFailed(format!(
            "Hopscotching failed at load factor {}",
            self.size as f64 / self.capacity as f64
        )))
    }

    /// `move`, a keyword here.
    fn relocate(&mut self, predecessor: i32, to_move: i32, empty: i32) {
        let mut predecessor = predecessor;
        let mut to_empty = self.index_diff(to_move, empty);
        let mut next_offset = self.offset(to_move);
        if next_offset == 0 || next_offset > to_empty {
            self.add_status(predecessor, to_empty);
        } else {
            self.add_status(predecessor, next_offset);
            to_empty -= next_offset;
            predecessor = self.index(to_move, next_offset);
            loop {
                next_offset = self.offset(predecessor);
                if next_offset == 0 || next_offset >= to_empty {
                    break;
                }
                to_empty -= next_offset;
                predecessor = self.index(predecessor, next_offset);
            }
            self.status[predecessor as usize] = to_empty as i8;
        }
        if next_offset != 0 {
            self.status[empty as usize] = (next_offset - to_empty) as i8;
        }
        self.buckets[empty as usize] = self.buckets[to_move as usize].take();
        self.status[to_move as usize] = 0;
    }

    fn resize(&mut self) -> Result<(), HopscotchFailed> {
        let old_capacity = self.capacity;
        let old_size = self.size;
        let new_capacity = legal_size_above(old_capacity);
        let old_buckets = std::mem::replace(&mut self.buckets, vec![None; new_capacity as usize]);
        let old_status = std::mem::replace(&mut self.status, vec![0; new_capacity as usize]);
        self.capacity = new_capacity;
        self.size = 0;
        // Not in bucket order: a stride of 127 around the old table.
        let mut index = 0;
        loop {
            if let Some(entry) = &old_buckets[index as usize] {
                let bucket = self.hash_to_index(entry.key_hash());
                if self.insert(entry, bucket, false).is_err() {
                    self.capacity = old_capacity;
                    self.size = old_size;
                    self.buckets = old_buckets;
                    self.status = old_status;
                    return Err(HopscotchFailed(format!(
                        "Hopscotching failed at load factor {}, and resizing didn't help.",
                        self.size as f64 / self.capacity as f64
                    )));
                }
            }
            index = ((index as i64 + 127) % old_capacity as i64) as i32;
            if index == 0 {
                break;
            }
        }
        Ok(())
    }
}

/// `computeCapacity`: the first legal size at or above `size / LOAD_FACTOR`.
fn compute_capacity(size: i32) -> i32 {
    if (size as f64) < LOAD_FACTOR * i32::MAX as f64 {
        let augmented = (size as f64 / LOAD_FACTOR) as i32;
        for legal in LEGAL_SIZES {
            if legal >= augmented {
                return legal;
            }
        }
    }
    LEGAL_SIZES[LEGAL_SIZES.len() - 1]
}

/// `SetSizeUtils.getLegalSizeAbove`, strictly above.
fn legal_size_above(elements: i32) -> i32 {
    LEGAL_SIZES
        .into_iter()
        .find(|&legal| legal > elements)
        .expect("no legal size large enough")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, PartialEq, Debug)]
    struct Key(i32);

    impl HopscotchEntry for Key {
        fn key_hash(&self) -> i32 {
            self.0
        }
        fn collides(&self, other: &Self) -> bool {
            self == other
        }
    }

    #[test]
    fn a_capacity_is_a_legal_size_above_the_load_factor() {
        // 208 / 0.85 is 244; the first legal size at or above it.
        assert_eq!(Hopscotch::<Key>::new(208).capacity(), 251);
        assert_eq!(Hopscotch::<Key>::new(0).capacity(), 251);
    }

    #[test]
    fn it_grows_and_holds_everything_once() {
        let mut table = Hopscotch::new(0);
        for value in 0..5000i32 {
            assert_eq!(table.add(Key(value.wrapping_mul(31))), Ok(true));
        }
        assert_eq!(table.add(Key(31)), Ok(false));
        assert_eq!(table.len(), 5000);
        let mut seen: Vec<i32> = table.iter().map(|key| key.0).collect();
        seen.sort_unstable();
        assert_eq!(seen, (0..5000).map(|v: i32| v * 31).collect::<Vec<_>>());
    }
}
