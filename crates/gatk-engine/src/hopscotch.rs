//! `LongHopscotchSet`, because the order it iterates in is the order PathSeq's file carries.
//!
//! # Why a set has to be ported at all
//!
//! `PathSeqBuildKmers` writes a `PSKmerSet` through Kryo, and that serializer writes the set's
//! CAPACITY, its size, and then every entry: the chain heads in bucket order and the squatters
//! after them. So the bytes depend on where each k-mer landed in the table, which depends on
//! FNV-1a, on the legal size above the capacity asked for, and on the hopscotch displacement that
//! ran while the set was filled. A port that holds the same k-mers in a `BTreeSet` answers a file
//! with the same numbers in a different order, which is a different file
//! (IPNP-BIPN/gatk-rs#1181).
//!
//! Nothing here is a Java collection: `LongHopscotchSet` is GATK's own, and this is a port of it
//! rather than of anything in the JDK.
//!
//! # The parts that decide the order
//!
//! - **the capacity is a legal size**, the first entry of a fixed table of primes strictly greater
//!   than the number asked for. Asking for eight and asking for sixty-four both give 251, which is
//!   what the `kryo-stream` golden measured;
//! - **the index is FNV-1a over the value's eight bytes**, truncated to an `int` and taken modulo
//!   the capacity, with a negative remainder brought back up;
//! - **a bucket holds its occupancy in the value's top bit**, so an entry must be non-negative, and
//!   a status byte holds "chain head" in its own top bit and the offset to the next entry of the
//!   chain in the low seven;
//! - **an entry that finds a squatter in its own bucket evicts it**, and a chain that cannot reach
//!   an empty bucket hopscotches one closer. Both move entries between buckets, so the order is a
//!   property of the whole insertion sequence and not of the values alone.

/// `SetSizeUtils.legalSizes`: the largest primes below each half power of two.
pub const LEGAL_SIZES: [i32; 47] = [
    251, 359, 509, 719, 1021, 1447, 2039, 2887, 4093, 5791, 8191, 11579, 16381, 23167, 32749,
    46337, 65521, 92681, 131071, 185363, 262139, 370723, 524287, 741431, 1048573, 1482907, 2097143,
    2965819, 4194301, 5931641, 8388593, 11863279, 16777213, 23726561, 33554393, 47453111, 67108859,
    94906249, 134217689, 189812507, 268435399, 379625047, 536870909, 759250111, 1073741789,
    1518500213, 2147483629,
];

/// `SetSizeUtils.getLegalSizeAbove(minElements)`, whose load factor defaults to one.
///
/// STRICTLY greater, which is why a request for 251 does not get 251.
pub fn legal_size_above(min_elements: i64) -> i32 {
    LEGAL_SIZES
        .iter()
        .copied()
        .find(|size| i64::from(*size) > min_elements)
        .expect("no legal size large enough")
}

/// `SVUtils.FNV64_DEFAULT_SEED`, the FNV-1a offset basis as a signed long.
const FNV64_DEFAULT_SEED: i64 = -3_750_763_034_362_895_579;

/// `SVUtils.fnvLong64(long)`: FNV-1a over the eight bytes, most significant first.
pub fn fnv_long64(to_hash: i64) -> i64 {
    const MULT: i64 = 1_099_511_628_211;
    let mut start = FNV64_DEFAULT_SEED;
    for shift in [56, 48, 40, 32, 24, 16, 8, 0] {
        start ^= (to_hash >> shift) & 0xff;
        start = start.wrapping_mul(MULT);
    }
    start
}

/// `LongHopscotchSet.longHash`: the low half of the FNV hash, as Java's narrowing cast gives it.
pub fn long_hash(value: i64) -> i32 {
    fnv_long64(value) as i32
}

/// A set of non-negative longs in GATK's own hopscotch table.
#[derive(Debug, Clone)]
pub struct LongHopscotchSet {
    capacity: i32,
    size: i32,
    /// The top bit says the bucket is occupied; the rest is the value.
    buckets: Vec<i64>,
    /// The top bit says the bucket is a chain head; the low seven are the offset to the next.
    status: Vec<i8>,
}

impl LongHopscotchSet {
    /// `new LongHopscotchSet(capacity)`, which rounds up to a legal size.
    pub fn with_capacity(capacity: i32) -> Self {
        let capacity = legal_size_above(i64::from(capacity));
        Self {
            capacity,
            size: 0,
            buckets: vec![0; capacity as usize],
            status: vec![0; capacity as usize],
        }
    }

    pub fn capacity(&self) -> i32 {
        self.capacity
    }

    pub fn size(&self) -> i32 {
        self.size
    }

    /// `add(entryValue)`: resize when full or when hopscotching gives up, then insert.
    ///
    /// # Panics
    ///
    /// On a negative value, which the reference refuses with "Tried to add negative entry to
    /// LongHopScotchSet": the top bit of a bucket is the occupancy flag, so a negative entry has
    /// nowhere to live.
    pub fn add(&mut self, value: i64) -> bool {
        assert!(
            value >= 0,
            "Tried to add negative entry to LongHopScotchSet"
        );
        if self.size == self.capacity {
            self.resize();
        }
        match self.insert(value) {
            Ok(added) => added,
            Err(Hopscotch::Failed) => {
                self.resize();
                self.insert(value).expect("a resized table has room")
            }
        }
    }

    /// The entries in the order `serialize` writes them: chain heads first, then squatters, each
    /// in bucket order.
    pub fn serialized_entries(&self) -> Vec<i64> {
        let mut entries = Vec::with_capacity(self.size as usize);
        for index in 0..self.capacity {
            if self.is_chain_head(index) {
                entries.push(value_of(self.buckets[index as usize]));
            }
        }
        for index in 0..self.capacity {
            let entry = self.buckets[index as usize];
            if !is_unused(entry) && !self.is_chain_head(index) {
                entries.push(value_of(entry));
            }
        }
        assert_eq!(
            entries.len() as i32,
            self.size,
            "Failed to serialize the expected number of objects"
        );
        entries
    }

    pub fn contains(&self, key: i64) -> bool {
        let mut index = self.hash_to_index(long_hash(key));
        if !self.is_chain_head(index) {
            return false;
        }
        if value_of(self.buckets[index as usize]) == key {
            return true;
        }
        loop {
            let offset = self.offset(index);
            if offset == 0 {
                return false;
            }
            index = self.index_at(index, offset);
            if value_of(self.buckets[index as usize]) == key {
                return true;
            }
        }
    }

    fn insert(&mut self, value: i64) -> Result<bool, Hopscotch> {
        let bucket_index = self.hash_to_index(long_hash(value));

        // A squatter in the bucket this entry belongs in is moved out of the way first.
        if !is_unused(self.buckets[bucket_index as usize]) && !self.is_chain_head(bucket_index) {
            self.evict(bucket_index)?;
        }

        if is_unused(self.buckets[bucket_index as usize]) {
            self.buckets[bucket_index as usize] = value;
            self.set_occupied(bucket_index);
            self.status[bucket_index as usize] = i8::MIN;
            self.size += 1;
            return Ok(true);
        }

        // Walk to the end of the chain, refusing a value the chain already holds.
        let mut end_of_chain = bucket_index;
        loop {
            if value_of(self.buckets[end_of_chain as usize]) == value {
                return Ok(false);
            }
            let offset = self.offset(end_of_chain);
            if offset == 0 {
                break;
            }
            end_of_chain = self.index_at(end_of_chain, offset);
        }

        let empty = self.insert_into_chain(bucket_index, end_of_chain)?;
        self.buckets[empty as usize] = value;
        self.set_occupied(empty);
        self.size += 1;
        Ok(true)
    }

    fn insert_into_chain(
        &mut self,
        bucket_index: i32,
        end_of_chain: i32,
    ) -> Result<i32, Hopscotch> {
        let offset_to_end = self.index_diff(bucket_index, end_of_chain);
        let mut empty = self.find_empty_bucket(bucket_index);
        let max_offset = offset_to_end + i32::from(i8::MAX);

        let mut offset_to_empty = self.index_diff(bucket_index, empty);
        while offset_to_empty > max_offset {
            empty = self.hopscotch(bucket_index, empty)?;
            offset_to_empty = self.index_diff(bucket_index, empty);
        }

        if offset_to_empty > offset_to_end {
            // Downstream of the chain's end, so it links straight on.
            let status = &mut self.status[end_of_chain as usize];
            *status = status.wrapping_add((offset_to_empty - offset_to_end) as i8);
        } else {
            self.link_into_chain(bucket_index, empty);
        }
        Ok(empty)
    }

    fn link_into_chain(&mut self, bucket_index: i32, empty: i32) {
        let mut offset_to_empty = self.index_diff(bucket_index, empty);
        let mut index = bucket_index;
        let mut offset = self.offset(index);
        while offset < offset_to_empty {
            index = self.index_at(index, offset);
            offset_to_empty -= offset;
            offset = self.offset(index);
        }
        offset -= offset_to_empty;
        self.status[index as usize] = self.status[index as usize].wrapping_sub(offset as i8);
        self.status[empty as usize] = offset as i8;
    }

    fn evict(&mut self, to_evict: i32) -> Result<(), Hopscotch> {
        let bucket_index = self.value_to_index(value_of(self.buckets[to_evict as usize]));
        let offset_to_evictee = self.index_diff(bucket_index, to_evict);
        let mut empty = self.find_empty_bucket(bucket_index);
        let mut from_index = bucket_index;
        loop {
            while self.index_diff(bucket_index, empty) > offset_to_evictee {
                empty = self.hopscotch(from_index, empty)?;
            }
            if empty == to_evict {
                return Ok(());
            }
            from_index = empty;
            self.link_into_chain(bucket_index, empty);

            // The chain's last entry moves into the hole, which is what frees the evictee's bucket
            // one step at a time.
            let mut previous = bucket_index;
            let mut offset_to_next = self.offset(previous);
            let mut next = self.index_at(previous, offset_to_next);
            loop {
                offset_to_next = self.offset(next);
                if offset_to_next == 0 {
                    break;
                }
                previous = next;
                next = self.index_at(next, offset_to_next);
            }
            self.buckets[empty as usize] = self.buckets[next as usize];
            self.buckets[next as usize] = 0;
            self.status[next as usize] = 0;
            self.status[previous as usize] =
                self.status[previous as usize].wrapping_sub(self.offset(previous) as i8);
            empty = next;
        }
    }

    fn find_empty_bucket(&self, mut index: i32) -> i32 {
        loop {
            index = self.index_at(index, 1);
            if is_unused(self.buckets[index as usize]) {
                return index;
            }
        }
    }

    /// `hopscotch`: pull an empty bucket closer by moving an entry into it.
    ///
    /// The reference throws `IllegalStateException` when it cannot, and `add` answers that by
    /// resizing and trying again, so the failure is a control-flow step rather than an error.
    fn hopscotch(&mut self, from_index: i32, empty: i32) -> Result<i32, Hopscotch> {
        let from_to_empty = self.index_diff(from_index, empty);
        let mut offset_to_empty = i32::from(i8::MAX);
        while offset_to_empty > 1 {
            let bucket_index = self.index_at(empty, -offset_to_empty);
            let offset_in_bucket = self.offset(bucket_index);
            if offset_in_bucket != 0
                && offset_in_bucket < offset_to_empty
                && offset_to_empty - offset_in_bucket < from_to_empty
            {
                let to_move = self.index_at(bucket_index, offset_in_bucket);
                self.move_entry(bucket_index, to_move, empty);
                return Ok(to_move);
            }
            offset_to_empty -= 1;
        }
        Err(Hopscotch::Failed)
    }

    fn move_entry(&mut self, predecessor: i32, to_move: i32, empty: i32) {
        let mut predecessor = predecessor;
        let mut to_empty = self.index_diff(to_move, empty);
        let mut next_offset = self.offset(to_move);
        if next_offset == 0 || next_offset > to_empty {
            self.status[predecessor as usize] =
                self.status[predecessor as usize].wrapping_add(to_empty as i8);
        } else {
            self.status[predecessor as usize] =
                self.status[predecessor as usize].wrapping_add(next_offset as i8);
            to_empty -= next_offset;
            predecessor = self.index_at(to_move, next_offset);
            loop {
                next_offset = self.offset(predecessor);
                if next_offset == 0 || next_offset >= to_empty {
                    break;
                }
                to_empty -= next_offset;
                predecessor = self.index_at(predecessor, next_offset);
            }
            self.status[predecessor as usize] = to_empty as i8;
        }
        if next_offset != 0 {
            self.status[empty as usize] = (next_offset - to_empty) as i8;
        }
        self.buckets[empty as usize] = self.buckets[to_move as usize];
        self.buckets[to_move as usize] = 0;
        self.status[to_move as usize] = 0;
    }

    /// `resize`: the next legal size, refilled by walking the old table in steps of 127.
    ///
    /// The step is not decoration. It spreads the old entries over the new table in an order that
    /// is neither the old bucket order nor the insertion order, and the entries' final places
    /// depend on it.
    fn resize(&mut self) {
        let old_capacity = self.capacity;
        let old_buckets = std::mem::take(&mut self.buckets);

        self.capacity = legal_size_above(i64::from(self.capacity));
        self.size = 0;
        self.buckets = vec![0; self.capacity as usize];
        self.status = vec![0; self.capacity as usize];

        let mut index: i32 = 0;
        loop {
            let entry = old_buckets[index as usize];
            if !is_unused(entry) {
                self.insert(value_of(entry))
                    .expect("a freshly resized table has room");
            }
            index = (index + 127) % old_capacity;
            if index == 0 {
                break;
            }
        }
    }

    fn value_to_index(&self, value: i64) -> i32 {
        self.hash_to_index(long_hash(value))
    }

    fn hash_to_index(&self, hash: i32) -> i32 {
        let mut result = hash % self.capacity;
        if result < 0 {
            result += self.capacity;
        }
        result
    }

    fn index_at(&self, index: i32, offset: i32) -> i32 {
        let mut result = index + offset;
        if result >= self.capacity {
            result -= self.capacity;
        } else if result < 0 {
            result += self.capacity;
        }
        result
    }

    /// The distance from the first bucket to the second, wrapping, and therefore never negative.
    fn index_diff(&self, first: i32, second: i32) -> i32 {
        let mut result = second - first;
        if result < 0 {
            result += self.capacity;
        }
        result
    }

    fn is_chain_head(&self, index: i32) -> bool {
        self.status[index as usize] & i8::MIN != 0
    }

    fn offset(&self, index: i32) -> i32 {
        i32::from(self.status[index as usize] & i8::MAX)
    }

    fn set_occupied(&mut self, index: i32) {
        self.buckets[index as usize] = (self.buckets[index as usize] as u64 | 1u64 << 63) as i64;
    }
}

/// `SetSizeUtils.getLegalSizeBelow(maxElements)`, which is what decides a partition count.
///
/// It refuses anything at or below the smallest legal size, and `LargeLongHopscotchSet` answers
/// that refusal with a single partition.
pub fn legal_size_below(max_elements: i64) -> Option<i32> {
    if max_elements <= i64::from(LEGAL_SIZES[0]) {
        return None;
    }
    for index in 1..LEGAL_SIZES.len() {
        if max_elements <= i64::from(LEGAL_SIZES[index]) {
            return Some(LEGAL_SIZES[index - 1]);
        }
    }
    Some(LEGAL_SIZES[LEGAL_SIZES.len() - 1])
}

/// `LargeLongHopscotchSet`: one hopscotch table per partition, chosen by the hash.
///
/// The partition count is the largest legal size below the SQUARE ROOT of the number of elements
/// the caller expects, which for anything small is one partition. Each partition is then asked for
/// `numElements / partitions + 1` and rounds that up to its own legal size.
#[derive(Debug, Clone)]
pub struct LargeLongHopscotchSet {
    sets: Vec<LongHopscotchSet>,
}

impl LargeLongHopscotchSet {
    pub fn new(num_elements: i64) -> Self {
        assert!(
            num_elements > 0,
            "Number of elements must be greater than 0"
        );
        let elements_per_partition = (num_elements as f64).sqrt() as i64;
        let partitions = legal_size_below(elements_per_partition).unwrap_or(1);
        let per_partition = (num_elements / i64::from(partitions)) + 1;
        Self {
            sets: (0..partitions)
                .map(|_| LongHopscotchSet::with_capacity(per_partition as i32))
                .collect(),
        }
    }

    /// `add`, which picks the partition by the hash read as UNSIGNED.
    pub fn add(&mut self, value: i64) -> bool {
        let hash = long_hash(value);
        let index = (hash as u32 % self.sets.len() as u32) as usize;
        self.sets[index].add(value)
    }

    pub fn contains(&self, value: i64) -> bool {
        let hash = long_hash(value);
        let index = (hash as u32 % self.sets.len() as u32) as usize;
        self.sets[index].contains(value)
    }

    pub fn sets(&self) -> &[LongHopscotchSet] {
        &self.sets
    }

    pub fn size(&self) -> i32 {
        self.sets.iter().map(LongHopscotchSet::size).sum()
    }
}

/// What `hopscotch` does when no entry can be moved: the reference throws, and `add` resizes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Hopscotch {
    Failed,
}

/// `getValue`: the entry with its occupancy bit cleared.
fn value_of(entry: i64) -> i64 {
    entry & i64::MAX
}

/// `isUnusedValue`: a bucket holding zero is empty, which is why zero is not a storable value.
fn is_unused(entry: i64) -> bool {
    entry == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_legal_size_is_strictly_above_the_request() {
        assert_eq!(legal_size_above(8), 251);
        assert_eq!(legal_size_above(64), 251);
        assert_eq!(legal_size_above(251), 359);
    }

    #[test]
    fn the_set_holds_what_it_was_given() {
        let mut set = LongHopscotchSet::with_capacity(8);
        for value in [1, 2, 3, 17, 1024, 65535, 1_048_577, 123_456_789] {
            assert!(set.add(value), "{value}");
        }
        assert_eq!(set.size(), 8);
        assert_eq!(set.capacity(), 251);
        for value in [1, 2, 3, 17, 1024, 65535, 1_048_577, 123_456_789] {
            assert!(set.contains(value), "{value}");
        }
        assert!(!set.contains(4));
        // The same value twice is one entry.
        assert!(!set.add(17));
        assert_eq!(set.size(), 8);
    }
}
