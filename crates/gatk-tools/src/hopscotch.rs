//! `HopscotchCollection` (GATK 4.7.0.0; `HopscotchCollection.java` is unchanged since 4.6.2.0),
//! whose iteration order is what a tool writes.
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
//!
//! # A table the size of the reference
//!
//! `KmerCleaner` asks for a table from twice the reference over the partition count, narrowed to an
//! `int`, whatever will go in it. The reference's cost is the arrays, four bytes of reference and a
//! status byte a bucket, which its heap holds or does not. Here the entries live elsewhere and a
//! bucket holds a number with 0 for empty, so the two arrays are zeroed pages that only what is
//! touched ever occupies. At the top of the legal sizes `bucket + offset` passes
//! `Integer.MAX_VALUE` and comes back negative, which `Hopscotch::index` does as Java does.
//!
//! # Reading this file without reading Rust
//!
//! - `Hopscotch<T: HopscotchEntry>` is a generic class whose element type must implement the
//!   [`HopscotchEntry`] interface (a Rust `trait`); `impl<T: ..> Hopscotch<T>` is its body;
//! - `Option<T>` is "a value or nothing" and `Result<T, E>` "a value or an error"; `?` after a
//!   call returns its error to the caller at once, which is how an `IllegalStateException`
//!   propagates here;
//! - `i8` is Java's signed `byte`; `x as i8` keeps the low eight bits, as Java's `(byte)` cast
//!   does, and `wrapping_add` / `wrapping_mul` are Java's `int` arithmetic (a plain `+` on Rust
//!   integers panics on overflow in a debug build);
//! - `&mut self` is a method that may modify the table; `&self` one that only reads it;
//! - `impl Fn(&T) -> bool` is a parameter of any lambda type taking an entry, Java's
//!   `Predicate<T>`.

/// `SetSizeUtils.legalSizes`: the sizes the table may take, all of them prime.
const LEGAL_SIZES: [i32; 47] = [
    251, 359, 509, 719, 1021, 1447, 2039, 2887, 4093, 5791, 8191, 11579, 16381, 23167, 32749,
    46337, 65521, 92681, 131071, 185363, 262139, 370723, 524287, 741431, 1048573, 1482907, 2097143,
    2965819, 4194301, 5931641, 8388593, 11863279, 16777213, 23726561, 33554393, 47453111, 67108859,
    94906249, 134217689, 189812507, 268435399, 379625047, 536870909, 759250111, 1073741789,
    1518500213, 2147483629,
];

/// `HopscotchCollection.LOAD_FACTOR`.
const LOAD_FACTOR: f64 = 0.85;
/// `HopscotchCollection.SPREADER`, the multiplier `hashToIndex` scatters a hash code with.
const SPREADER: i32 = 241;

/// An entry's `toKey().apply(entry).hashCode()` and `entryCollides()`.
///
/// `Clone` is required so an entry can be copied; `collides` is the `BiPredicate` `add` tests
/// against every entry of the chain, `(t1, t2) -> false` by default in the Java (a collection
/// that admits duplicates).
pub trait HopscotchEntry: Clone {
    fn key_hash(&self) -> i32;
    fn collides(&self, other: &Self) -> bool;
}

/// `IllegalStateException` out of a resize that could not place what it already held.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HopscotchFailed(pub String);

/// `HopscotchCollection<T>`.
///
/// The Java's `T[] buckets` holds the entries themselves; here `buckets` holds a number that
/// indexes `entries` (see the module note on the size of the table), with `free` the numbers a
/// removal gave back for reuse. `status` is the Java's `byte[] status`, one per bucket.
pub struct Hopscotch<T: HopscotchEntry> {
    capacity: i32,
    size: i32,
    /// One more than the entry's place in `entries`, and 0 for an empty bucket.
    buckets: Vec<u32>,
    status: Vec<i8>,
    entries: Vec<Option<T>>,
    free: Vec<u32>,
}

impl<T: HopscotchEntry> Hopscotch<T> {
    /// `new HopscotchCollection(capacity)`, through `computeCapacity`.
    pub fn new(capacity: i32) -> Hopscotch<T> {
        let capacity = compute_capacity(capacity);
        Hopscotch {
            capacity,
            size: 0,
            buckets: vec![0; capacity as usize],
            status: vec![0; capacity as usize],
            entries: Vec::new(),
            free: Vec::new(),
        }
    }

    /// `buckets[bucket]`, `None` for the Java's null.
    fn entry(&self, bucket: i32) -> Option<&T> {
        match self.buckets[bucket as usize] {
            0 => None,
            held => self.entries[held as usize - 1].as_ref(),
        }
    }

    fn is_empty_bucket(&self, bucket: i32) -> bool {
        self.buckets[bucket as usize] == 0
    }

    /// An entry, under the number a bucket will hold it by.
    fn store(&mut self, entry: T) -> u32 {
        match self.free.pop() {
            Some(held) => {
                self.entries[held as usize - 1] = Some(entry);
                held
            }
            None => {
                self.entries.push(Some(entry));
                self.entries.len() as u32
            }
        }
    }

    /// Drop the entry under `held` and give the number back.
    fn release(&mut self, held: u32) {
        self.entries[held as usize - 1] = None;
        self.free.push(held);
    }

    /// Move what a bucket holds to another, which leaves it empty.
    ///
    /// `buckets[to] = buckets[from]; buckets[from] = null;`. `std::mem::take` reads the value
    /// and leaves the default, 0, in its place.
    fn shift(&mut self, to: i32, from: i32) {
        self.buckets[to as usize] = std::mem::take(&mut self.buckets[from as usize]);
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
    ///
    /// The retry passes the `bucket` computed before the resize, as the Java does (see the module
    /// note). The entry is stored first so that `insert` can place it by number; if it was not
    /// added, the number is given back.
    pub fn add(&mut self, entry: T) -> Result<bool, HopscotchFailed> {
        if self.size == self.capacity {
            self.resize()?;
        }
        let bucket = self.hash_to_index(entry.key_hash());
        let held = self.store(entry);
        let added = match self.insert(held, bucket, true) {
            Ok(added) => Ok(added),
            Err(_) => match self.resize() {
                Ok(()) => self.insert(held, bucket, true),
                Err(failed) => Err(failed),
            },
        };
        if !matches!(added, Ok(true)) {
            self.release(held);
        }
        added
    }

    /// `find`: the chain the key hashes to, and no other.
    ///
    /// `hash` is the key's `hashCode()` and `matches` the Java's `equivalent(entry, key)`. The
    /// entry comes back as a writable borrow, so a caller can update it in place.
    pub fn find_mut(&mut self, hash: i32, matches: impl Fn(&T) -> bool) -> Option<&mut T> {
        let mut bucket = self.hash_to_index(hash);
        if !self.is_chain_head(bucket) {
            return None;
        }
        loop {
            if matches(self.entry(bucket)?) {
                let held = self.buckets[bucket as usize] as usize;
                return self.entries[held - 1].as_mut();
            }
            let offset = self.offset(bucket);
            if offset == 0 {
                return None;
            }
            bucket = self.index(bucket, offset);
        }
    }

    /// `iterator`: every chain from its head, heads in bucket order.
    ///
    /// `CompleteIterator`: `nextBucketHead` walks the buckets upward for the next chain head, and
    /// `next` follows the chain's offsets until one is zero. The order is collected into a list
    /// first and handed out as an iterator (`impl Iterator`, Java's `Iterator<T>`); the `+ '_`
    /// says it borrows the table.
    pub fn iter(&self) -> impl Iterator<Item = &T> + '_ {
        let mut order = Vec::with_capacity(self.size as usize);
        for head in 0..self.capacity {
            if !self.is_chain_head(head) {
                continue;
            }
            let mut bucket = head;
            loop {
                if let Some(entry) = self.entry(bucket) {
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

    /// `contains`, of whatever is in the bucket for this hash: `find(key) != null`.
    pub fn contains(&self, hash: i32, matches: impl Fn(&T) -> bool) -> bool {
        !self.find_each(hash, matches).is_empty()
    }

    /// `findEach`: what the chain holds of a key, in link order.
    ///
    /// The `ElementIterator` walk, collected into a list. `expect` aborts with its message if a
    /// chain link is empty, which the table's invariant rules out.
    pub fn find_each(&self, hash: i32, matches: impl Fn(&T) -> bool) -> Vec<&T> {
        let mut found = Vec::new();
        let mut bucket = self.hash_to_index(hash);
        if !self.is_chain_head(bucket) {
            return found;
        }
        loop {
            let entry = self.entry(bucket).expect("a chain holds no empty bucket");
            if matches(entry) {
                found.push(entry);
            }
            let offset = self.offset(bucket);
            if offset == 0 {
                return found;
            }
            bucket = self.index(bucket, offset);
        }
    }

    /// `removeEach`, an `ElementIterator` removing all it finds: a removal moves the chain's last
    /// entry into the hole, so the walk backs up to look at what arrived there.
    pub fn remove_each(&mut self, hash: i32, matches: impl Fn(&T) -> bool) -> bool {
        // `NO_ELEMENT_INDEX`.
        const NONE: i32 = -1;
        let bucket = self.hash_to_index(hash);
        if !self.is_chain_head(bucket) {
            return false;
        }
        let mut current = bucket;
        let mut previous = NONE;
        // `advanceUntilEquivalent`, as a closure handed the table and the two cursors, since a
        // closure that held them itself would block the removals below.
        let advance = |table: &Self, current: &mut i32, previous: &mut i32| loop {
            let offset = table.offset(*current);
            if offset == 0 {
                *current = NONE;
                return;
            }
            *previous = *current;
            *current = table.index(*current, offset);
            if matches(table.entry(*current).expect("a link")) {
                return;
            }
        };
        // The `ElementIterator` constructor's `ensureEquivalence`.
        if !matches(self.entry(current).expect("a head")) {
            advance(self, &mut current, &mut previous);
        }
        let mut removed = false;
        while current != NONE {
            // `next()` then `remove()`: note the element and advance, remove it, and if the
            // removal pulled the chain's last entry into the hole, back up onto it; then
            // `ensureEquivalence`.
            let (remove_at, remove_previous) = (current, previous);
            advance(self, &mut current, &mut previous);
            self.remove_at_index(remove_at, remove_previous);
            if !self.is_empty_bucket(remove_at) {
                current = remove_at;
                previous = remove_previous;
            }
            if current != NONE && !matches(self.entry(current).expect("a link")) {
                advance(self, &mut current, &mut previous);
            }
            removed = true;
        }
        removed
    }

    /// `Collection.removeIf`, through `CompleteIterator`.
    ///
    /// The same `next()` and `remove()` dance as [`Hopscotch::remove_each`], over every chain, with
    /// `next_head` the iterator's `nextBucketHead`.
    pub fn remove_if(&mut self, doomed: impl Fn(&T) -> bool) {
        const NONE: i32 = -1;
        let mut head = NONE;
        let mut current = NONE;
        let mut previous = NONE;
        let next_head = |table: &Self, head: &mut i32, current: &mut i32, previous: &mut i32| {
            *head += 1;
            while *head < table.capacity {
                if table.is_chain_head(*head) {
                    *current = *head;
                    *previous = NONE;
                    return;
                }
                *head += 1;
            }
            *head = NONE;
        };
        next_head(self, &mut head, &mut current, &mut previous);
        while head != NONE {
            let (remove_at, remove_previous) = (current, previous);
            let offset = self.offset(current);
            if offset == 0 {
                next_head(self, &mut head, &mut current, &mut previous);
            } else {
                previous = current;
                current = self.index(current, offset);
            }
            if doomed(self.entry(remove_at).expect("a link")) {
                self.remove_at_index(remove_at, remove_previous);
                if !self.is_empty_bucket(remove_at) {
                    current = remove_at;
                    previous = remove_previous;
                }
            }
        }
    }

    /// `removeAtIndex`.
    ///
    /// At the end of a chain the bucket is emptied and the predecessor's offset zeroed
    /// (`status[predecessorIndex] -= getOffset(predecessorIndex)`, which keeps the head bit). In
    /// the middle, the chain's LAST entry moves into the hole and the link to it is cut, so the
    /// chain's order changes: that is why the iterators back up after a removal.
    fn remove_at_index(&mut self, bucket: i32, predecessor: i32) {
        let offset = self.offset(bucket);
        self.release(self.buckets[bucket as usize]);
        if offset == 0 {
            self.buckets[bucket as usize] = 0;
            self.status[bucket as usize] = 0;
            if predecessor != -1 {
                let offset = self.offset(predecessor);
                self.add_status(predecessor, -offset);
            }
        } else {
            let mut previous = bucket;
            let mut next = self.index(previous, offset);
            loop {
                let offset_to_next = self.offset(next);
                if offset_to_next == 0 {
                    break;
                }
                previous = next;
                next = self.index(next, offset_to_next);
            }
            self.shift(bucket, next);
            let offset = self.offset(previous);
            self.add_status(previous, -offset);
        }
        self.size -= 1;
    }

    /// `hashToIndex`: `(SPREADER * entry.hashCode()) % capacity`, the product wrapping as an
    /// `int` does, and a negative remainder (Java's `%` keeps the dividend's sign, as Rust's does)
    /// moved up by the capacity.
    fn hash_to_index(&self, hash: i32) -> i32 {
        let result = SPREADER.wrapping_mul(hash) % self.capacity;
        if result < 0 {
            result + self.capacity
        } else {
            result
        }
    }

    /// `insert`, of an entry already stored under `held`.
    ///
    /// A squatter in the entry's own bucket is evicted first; an empty bucket takes the entry as a
    /// chain head (`status[bucketIndex] = Byte.MIN_VALUE`); otherwise the chain is walked to its
    /// end, testing each entry for a collision when `test` holds, and the entry goes into an empty
    /// bucket linked into the chain. `test` is false during a resize, the Java's
    /// `(t1, t2) -> false`.
    fn insert(&mut self, held: u32, bucket: i32, test: bool) -> Result<bool, HopscotchFailed> {
        if !self.is_empty_bucket(bucket) && !self.is_chain_head(bucket) {
            self.evict(bucket)?;
        }
        if self.is_empty_bucket(bucket) {
            self.buckets[bucket as usize] = held;
            self.status[bucket as usize] = i8::MIN;
            self.size += 1;
            return Ok(true);
        }
        let mut end_of_chain = bucket;
        loop {
            if test {
                let entry = self.entries[held as usize - 1].as_ref().expect("stored");
                if self
                    .entry(end_of_chain)
                    .is_some_and(|in_chain| in_chain.collides(entry))
                {
                    return Ok(false);
                }
            }
            let offset = self.offset(end_of_chain);
            if offset == 0 {
                break;
            }
            end_of_chain = self.index(end_of_chain, offset);
        }
        let empty = self.insert_into_chain(bucket, end_of_chain)?;
        self.buckets[empty as usize] = held;
        self.size += 1;
        Ok(true)
    }

    /// `insertIntoChain`: find an empty bucket after `bucket`, hopscotch it back while it is more
    /// than `Byte.MAX_VALUE` past the chain's end, then link it in after the end or, when it lies
    /// before the end, at its place in the chain. Returns the bucket.
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

    /// `linkIntoChain`: walk the chain to the link that jumps past `empty` and split that jump in
    /// two, `status[tmpIndex] -= offset; status[emptyBucketIndex] = (byte) offset;`.
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

    /// `evict`: free `to_evict` for its rightful head by moving its chain's entries forward,
    /// each time hopscotching an empty bucket into range and moving the chain's last entry into
    /// it, until the empty bucket is `to_evict` itself.
    fn evict(&mut self, to_evict: i32) -> Result<(), HopscotchFailed> {
        let hash = self.entry(to_evict).expect("a squatter").key_hash();
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
            self.shift(empty, next);
            self.status[next as usize] = 0;
            let offset = self.offset(previous);
            self.add_status(previous, -offset);
            empty = next;
        }
    }

    /// `findEmptyBucket`: the first empty bucket strictly after `bucket`, wrapping. It never
    /// fails to stop because `add` resizes a full table first.
    fn find_empty_bucket(&self, mut bucket: i32) -> i32 {
        loop {
            bucket = self.index(bucket, 1);
            if self.is_empty_bucket(bucket) {
                return bucket;
            }
        }
    }

    /// `(status[bucketIndex] & Byte.MIN_VALUE) != 0`: the high bit is the sign bit of an `i8`.
    fn is_chain_head(&self, bucket: i32) -> bool {
        self.status[bucket as usize] < 0
    }

    /// `status[bucketIndex] & Byte.MAX_VALUE`, the low seven bits.
    fn offset(&self, bucket: i32) -> i32 {
        (self.status[bucket as usize] & i8::MAX) as i32
    }

    /// `status[bucket] += delta`, a compound assignment on a byte.
    fn add_status(&mut self, bucket: i32, delta: i32) {
        let status = &mut self.status[bucket as usize];
        *status = (*status as i32).wrapping_add(delta) as i8;
    }

    /// `getIndex`, whose sum is an `int`.
    fn index(&self, bucket: i32, offset: i32) -> i32 {
        let result = bucket.wrapping_add(offset);
        if result >= self.capacity {
            result - self.capacity
        } else if result < 0 {
            result + self.capacity
        } else {
            result
        }
    }

    /// `getIndexDiff`: how far `second` lies after `first`, going round the end of the table.
    fn index_diff(&self, first: i32, second: i32) -> i32 {
        let result = second - first;
        if result < 0 {
            result + self.capacity
        } else {
            result
        }
    }

    /// `hopscotch`: bring the empty bucket closer by moving into it an entry from up to
    /// `Byte.MAX_VALUE - 1` buckets before it whose predecessor link can still reach it, scanning
    /// from the farthest. Returns the bucket emptied, or the `IllegalStateException` when no entry
    /// can move.
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
    ///
    /// Move the entry at `to_move` into `empty` and relink: the predecessor's offset now reaches
    /// past `to_move`, and `empty` is spliced into the chain where its distance falls.
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
        self.shift(empty, to_move);
        self.status[to_move as usize] = 0;
    }

    /// `resize`: the next legal size strictly above the capacity, every old entry reinserted
    /// without a collision test, and on failure the old table restored and the second
    /// `IllegalStateException` raised. `std::mem::replace` swaps the new arrays in and hands the
    /// old ones back. The Java's final `size != oldSize` check is not ported.
    ///
    /// The message's load factor is printed with Rust's `{}`, which writes a whole number without
    /// `.0`: a resize triggered by a full table restores `size == capacity`, so the port writes
    /// `load factor 1` where the Java's `"..."+1.*size/capacity` writes `load factor 1.0`.
    fn resize(&mut self) -> Result<(), HopscotchFailed> {
        let old_capacity = self.capacity;
        let old_size = self.size;
        let new_capacity = legal_size_above(old_capacity);
        let old_buckets = std::mem::replace(&mut self.buckets, vec![0; new_capacity as usize]);
        let old_status = std::mem::replace(&mut self.status, vec![0; new_capacity as usize]);
        self.capacity = new_capacity;
        self.size = 0;
        // Not in bucket order: a stride of 127 around the old table,
        // `while ( (idx = (idx+127)%oldCapacity) != 0 )`. The capacity is prime, so every bucket
        // is visited once. The sum is taken as `i64` so `idx + 127` cannot wrap.
        let mut index = 0;
        loop {
            let held = old_buckets[index as usize];
            if held != 0 {
                let entry = self.entries[held as usize - 1].as_ref().expect("stored");
                let bucket = self.hash_to_index(entry.key_hash());
                if self.insert(held, bucket, false).is_err() {
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
///
/// `as i32` truncates toward zero, as the Java's `(int)` cast does; the guard keeps the quotient
/// within range.
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
///
/// The load factor is 1.0 there, so the Java's `(long) (minElements / loadFactor)` is the count
/// itself. Past the largest legal size the Java throws `IllegalArgumentException`; this aborts.
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
