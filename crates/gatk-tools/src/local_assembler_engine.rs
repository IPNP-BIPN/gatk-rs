//! `LocalAssembler`'s traversal and assembly: the reads paired as `PairWalker` pairs them, the de
//! Bruijn graph they make, and the two files written from it.
//!
//! Ported from `org.broadinstitute.hellbender.tools.LocalAssembler`,
//! `org.broadinstitute.hellbender.tools.walkers.PairWalker` and
//! `org.broadinstitute.hellbender.utils.collections.HopscotchCollection`.
//!
//! What decides the output, and therefore what this reproduces rather than approximates:
//!
//!  * **the contig numbers are the kmer set's iteration order.** `buildContigs` numbers a contig
//!    as it meets its first kmer walking a `HopscotchSet`, so the set is ported bucket for bucket:
//!    the capacity `computeCapacity` picks from ten times the padded region, the `SPREADER`
//!    index, eviction, hopscotching and removal, each of which moves entries between buckets;
//!  * **the reads are kmerized in the order `PairWalker` hands them over**: a pair when its second
//!    read arrives, first-seen read first, and the reads whose mate never came at the end, in the
//!    iteration order of the pair buffer, which is another `HopscotchSet`, keyed by the read
//!    name's `String.hashCode`;
//!  * **gap fills are kmerized in `HashMap<String, Integer>` order**, which is a fixed function of
//!    the strings' hash codes and the table's size, so it is reproduced as well.
//!
//! What is NOT reproducible, and why: the traversal set is a `HashSet<Traversal>` whose hash is
//! the list hash of the contigs', and a contig's hash is its identity hash. The order of the GFA's
//! `O` lines, the order of the FASTA's records under `--no-scaffolding`, and which orientation a
//! traversal found from both of its ends is kept in, are therefore the JVM's and change from one
//! machine to the next. The transit map is an identity-keyed `HashMap` too; this port walks it in
//! insertion order. Scaffolds are sorted by contig id before they are written, so the default
//! FASTA is fully determined.

use std::collections::{HashMap, HashSet};

use gatk_engine::read;
use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::Tag;

/// `Kmer.KSIZE`.
pub const KSIZE: usize = 31;
const KMASK: i64 = (1i64 << (2 * KSIZE)) - 1;
/// `HopscotchCollection.SPREADER`.
const SPREADER: i32 = 241;
/// `HopscotchCollection.LOAD_FACTOR`.
const LOAD_FACTOR: f64 = 0.85;

/// A failure the tool raises, as the class GATK would print and its message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyError {
    pub class: &'static str,
    pub message: String,
}

impl AssemblyError {
    fn gatk(message: impl Into<String>) -> Self {
        AssemblyError {
            class: "org.broadinstitute.hellbender.exceptions.GATKException",
            message: message.into(),
        }
    }
}

// ---------------------------------------------------------------------------------------------
// HopscotchCollection
// ---------------------------------------------------------------------------------------------

/// `HopscotchCollection.computeCapacity`: the first legal size at or above `size / 0.85`.
pub fn compute_capacity(size: i32) -> i32 {
    if (size as f64) < LOAD_FACTOR * i32::MAX as f64 {
        let augmented = (size as f64 / LOAD_FACTOR) as i32;
        for legal in gatk_engine::hopscotch::LEGAL_SIZES {
            if legal >= augmented {
                return legal;
            }
        }
    }
    *gatk_engine::hopscotch::LEGAL_SIZES.last().unwrap()
}

/// A failed hopscotch, which `add` and `findOrAdd` answer with a resize.
struct HopscotchFailed;

/// `HopscotchCollection<T>` with set semantics (`HopscotchSet`), bucket for bucket.
///
/// `index_of` is the collection's `hashToIndex` applied to an entry's key, and `same_key` its
/// `equivalent`: the set's order is a property of both and of every insertion and removal.
pub struct Hopscotch<T: Clone> {
    capacity: i32,
    size: i32,
    buckets: Vec<Option<T>>,
    status: Vec<u8>,
    index_of: fn(&T, i32) -> i32,
}

impl<T: Clone> Hopscotch<T> {
    pub fn new(capacity: i32, index_of: fn(&T, i32) -> i32) -> Self {
        let capacity = compute_capacity(capacity);
        Hopscotch {
            capacity,
            size: 0,
            buckets: vec![None; capacity as usize],
            status: vec![0; capacity as usize],
            index_of,
        }
    }

    pub fn capacity(&self) -> i32 {
        self.capacity
    }

    fn is_chain_head(&self, index: i32) -> bool {
        self.status[index as usize] & 0x80 != 0
    }

    fn offset(&self, index: i32) -> i32 {
        (self.status[index as usize] & 0x7f) as i32
    }

    fn index(&self, index: i32, offset: i32) -> i32 {
        let result = index + offset;
        if result >= self.capacity {
            result - self.capacity
        } else if result < 0 {
            result + self.capacity
        } else {
            result
        }
    }

    fn diff(&self, first: i32, second: i32) -> i32 {
        let result = second - first;
        if result < 0 {
            result + self.capacity
        } else {
            result
        }
    }

    fn status_add(&mut self, index: i32, delta: i32) {
        let slot = &mut self.status[index as usize];
        *slot = (*slot as i8).wrapping_add(delta as i8) as u8;
    }

    /// `find`: the entry `key_index` names whose key `matches` accepts.
    pub fn find(&self, key_index: i32, matches: impl Fn(&T) -> bool) -> Option<&T> {
        let mut index = key_index;
        if !self.is_chain_head(index) {
            return None;
        }
        loop {
            let entry = self.buckets[index as usize].as_ref()?;
            if matches(entry) {
                return Some(entry);
            }
            let offset = self.offset(index);
            if offset == 0 {
                return None;
            }
            index = self.index(index, offset);
        }
    }

    /// `add`, with the set's collision rule: an entry whose key is already present is refused.
    pub fn add(&mut self, entry: T, same_key: impl Fn(&T, &T) -> bool) -> bool {
        if self.size == self.capacity {
            self.resize();
        }
        let bucket = (self.index_of)(&entry, self.capacity);
        match self.insert(entry.clone(), bucket, &same_key) {
            Ok(added) => added,
            Err(HopscotchFailed) => {
                self.resize();
                let bucket = (self.index_of)(&entry, self.capacity);
                self.insert(entry, bucket, &same_key)
                    .unwrap_or_else(|_| panic!("Hopscotching failed after a resize"))
            }
        }
    }

    /// `findOrAdd`: the entry with this key, or `produce()` inserted at the end of its chain.
    pub fn find_or_add(
        &mut self,
        key_index: impl Fn(i32) -> i32,
        matches: impl Fn(&T) -> bool,
        produce: impl Fn() -> T,
    ) -> T {
        match self.find_or_add_internal(key_index(self.capacity), &matches, &produce) {
            Ok(entry) => entry,
            Err(HopscotchFailed) => {
                self.resize();
                self.find_or_add_internal(key_index(self.capacity), &matches, &produce)
                    .unwrap_or_else(|_| panic!("Hopscotching failed after a resize"))
            }
        }
    }

    fn find_or_add_internal(
        &mut self,
        bucket: i32,
        matches: &impl Fn(&T) -> bool,
        produce: &impl Fn() -> T,
    ) -> Result<T, HopscotchFailed> {
        if !self.is_chain_head(bucket) {
            let entry = produce();
            self.insert(entry.clone(), bucket, &|_: &T, _: &T| false)?;
            return Ok(entry);
        }
        let mut end = bucket;
        loop {
            let entry = self.buckets[end as usize].as_ref().expect("a chain entry");
            if matches(entry) {
                return Ok(entry.clone());
            }
            let offset = self.offset(end);
            if offset == 0 {
                break;
            }
            end = self.index(end, offset);
        }
        let entry = produce();
        let empty = self.insert_into_chain(bucket, end)?;
        self.buckets[empty as usize] = Some(entry.clone());
        self.size += 1;
        Ok(entry)
    }

    /// `remove(key)`.
    pub fn remove(&mut self, key_index: i32, matches: impl Fn(&T) -> bool) -> bool {
        let mut index = key_index;
        if self.buckets[index as usize].is_none() || !self.is_chain_head(index) {
            return false;
        }
        let mut predecessor = -1;
        while !matches(
            self.buckets[index as usize]
                .as_ref()
                .expect("a chain entry"),
        ) {
            let offset = self.offset(index);
            if offset == 0 {
                return false;
            }
            predecessor = index;
            index = self.index(index, offset);
        }
        self.remove_at(index, predecessor);
        true
    }

    fn remove_at(&mut self, index: i32, predecessor: i32) {
        let offset = self.offset(index);
        if offset == 0 {
            self.buckets[index as usize] = None;
            self.status[index as usize] = 0;
            if predecessor != -1 {
                let off = self.offset(predecessor);
                self.status_add(predecessor, -off);
            }
        } else {
            let mut prev = index;
            let mut next = self.index(prev, offset);
            loop {
                let to_next = self.offset(next);
                if to_next == 0 {
                    break;
                }
                prev = next;
                next = self.index(next, to_next);
            }
            self.buckets[index as usize] = self.buckets[next as usize].take();
            let off = self.offset(prev);
            self.status_add(prev, -off);
        }
        self.size -= 1;
    }

    /// The `CompleteIterator` order: chain heads by bucket, each chain followed to its end.
    pub fn entries(&self) -> Vec<T> {
        let mut result = Vec::with_capacity(self.size as usize);
        for head in 0..self.capacity {
            if !self.is_chain_head(head) {
                continue;
            }
            let mut index = head;
            loop {
                result.push(self.buckets[index as usize].clone().expect("a chain entry"));
                let offset = self.offset(index);
                if offset == 0 {
                    break;
                }
                index = self.index(index, offset);
            }
        }
        result
    }

    fn insert(
        &mut self,
        entry: T,
        bucket: i32,
        collides: &impl Fn(&T, &T) -> bool,
    ) -> Result<bool, HopscotchFailed> {
        if self.buckets[bucket as usize].is_some() && !self.is_chain_head(bucket) {
            self.evict(bucket)?;
        }
        if self.buckets[bucket as usize].is_none() {
            self.buckets[bucket as usize] = Some(entry);
            self.status[bucket as usize] = 0x80;
            self.size += 1;
            return Ok(true);
        }
        let mut end = bucket;
        loop {
            if collides(
                self.buckets[end as usize].as_ref().expect("a chain entry"),
                &entry,
            ) {
                return Ok(false);
            }
            let offset = self.offset(end);
            if offset == 0 {
                break;
            }
            end = self.index(end, offset);
        }
        let empty = self.insert_into_chain(bucket, end)?;
        self.buckets[empty as usize] = Some(entry);
        self.size += 1;
        Ok(true)
    }

    fn insert_into_chain(&mut self, bucket: i32, end: i32) -> Result<i32, HopscotchFailed> {
        let to_end = self.diff(bucket, end);
        let mut empty = self.find_empty(bucket);
        let max_offset = to_end + 127;
        let mut to_empty;
        loop {
            to_empty = self.diff(bucket, empty);
            if to_empty <= max_offset {
                break;
            }
            empty = self.hopscotch(bucket, empty)?;
        }
        if to_empty > to_end {
            self.status_add(end, to_empty - to_end);
        } else {
            self.link_into_chain(bucket, empty);
        }
        Ok(empty)
    }

    fn link_into_chain(&mut self, bucket: i32, empty: i32) {
        let mut to_empty = self.diff(bucket, empty);
        let mut index = bucket;
        let mut offset;
        loop {
            offset = self.offset(index);
            if offset >= to_empty {
                break;
            }
            index = self.index(index, offset);
            to_empty -= offset;
        }
        offset -= to_empty;
        self.status_add(index, -offset);
        self.status[empty as usize] = offset as u8;
    }

    fn evict(&mut self, to_evict: i32) -> Result<(), HopscotchFailed> {
        let bucket = (self.index_of)(
            self.buckets[to_evict as usize].as_ref().expect("an entry"),
            self.capacity,
        );
        let to_evictee = self.diff(bucket, to_evict);
        let mut empty = self.find_empty(bucket);
        let mut from = bucket;
        loop {
            while self.diff(bucket, empty) > to_evictee {
                empty = self.hopscotch(from, empty)?;
            }
            if empty == to_evict {
                return Ok(());
            }
            from = empty;
            self.link_into_chain(bucket, empty);
            let mut prev = bucket;
            let mut to_next = self.offset(prev);
            let mut next = self.index(prev, to_next);
            loop {
                to_next = self.offset(next);
                if to_next == 0 {
                    break;
                }
                prev = next;
                next = self.index(next, to_next);
            }
            self.buckets[empty as usize] = self.buckets[next as usize].take();
            self.status[next as usize] = 0;
            let off = self.offset(prev);
            self.status_add(prev, -off);
            empty = next;
        }
    }

    fn find_empty(&self, mut index: i32) -> i32 {
        loop {
            index = self.index(index, 1);
            if self.buckets[index as usize].is_none() {
                return index;
            }
        }
    }

    fn hopscotch(&mut self, from: i32, empty: i32) -> Result<i32, HopscotchFailed> {
        let from_to_empty = self.diff(from, empty);
        let mut to_empty = 127;
        while to_empty > 1 {
            let bucket = self.index(empty, -to_empty);
            let in_bucket = self.offset(bucket);
            if in_bucket != 0 && in_bucket < to_empty && to_empty - in_bucket < from_to_empty {
                let to_move = self.index(bucket, in_bucket);
                self.move_entry(bucket, to_move, empty);
                return Ok(to_move);
            }
            to_empty -= 1;
        }
        Err(HopscotchFailed)
    }

    fn move_entry(&mut self, predecessor_arg: i32, to_move: i32, empty: i32) {
        let mut predecessor = predecessor_arg;
        let mut to_empty = self.diff(to_move, empty);
        let mut next_offset = self.offset(to_move);
        if next_offset == 0 || next_offset > to_empty {
            self.status_add(predecessor, to_empty);
        } else {
            self.status_add(predecessor, next_offset);
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
            self.status[predecessor as usize] = to_empty as u8;
        }
        if next_offset != 0 {
            self.status[empty as usize] = (next_offset - to_empty) as u8;
        }
        self.buckets[empty as usize] = self.buckets[to_move as usize].take();
        self.status[to_move as usize] = 0;
    }

    fn resize(&mut self) {
        let old_capacity = self.capacity;
        let old_size = self.size;
        let old_buckets = std::mem::take(&mut self.buckets);
        self.capacity = gatk_engine::hopscotch::legal_size_above(old_capacity as i64);
        self.size = 0;
        self.buckets = vec![None; self.capacity as usize];
        self.status = vec![0; self.capacity as usize];
        let mut index = 0usize;
        loop {
            if let Some(entry) = &old_buckets[index] {
                let bucket = (self.index_of)(entry, self.capacity);
                if self
                    .insert(entry.clone(), bucket, &|_: &T, _: &T| false)
                    .is_err()
                {
                    panic!("Hopscotching failed at load factor, and resizing didn't help.");
                }
            }
            index = (index + 127) % old_capacity as usize;
            if index == 0 {
                break;
            }
        }
        assert_eq!(self.size, old_size, "Lost some elements during resizing.");
    }
}

/// `String.hashCode`, over the UTF-16 units.
pub fn java_string_hash(text: &str) -> i32 {
    text.encode_utf16().fold(0i32, |hash, unit| {
        hash.wrapping_mul(31).wrapping_add(unit as i32)
    })
}

/// `HopscotchCollection.hashToIndex` for an entry whose `hashCode` is `hash`.
fn default_index(hash: i32, capacity: i32) -> i32 {
    let result = SPREADER.wrapping_mul(hash) % capacity;
    if result < 0 {
        result + capacity
    } else {
        result
    }
}

/// `KmerSet.hashToIndex`.
fn kmer_index(kval: i64, capacity: i32) -> i32 {
    let positive = (SPREADER as i64).wrapping_mul(kval) & i64::MAX;
    (positive % capacity as i64) as i32
}

/// The iteration order of a `java.util.HashMap<String, _>` filled in `keys` order.
///
/// A bucket is `(h ^ h >>> 16) & (n - 1)` for a table of `n`, which starts at sixteen and doubles
/// whenever the size passes three quarters of it; a resize keeps each bucket's relative order, so
/// the order is the buckets in turn, each in insertion order.
pub fn java_hash_map_order(keys: &[String]) -> Vec<usize> {
    let mut capacity = 16usize;
    while keys.len() > capacity * 3 / 4 {
        capacity *= 2;
    }
    let bucket = |key: &str| {
        let hash = java_string_hash(key);
        ((hash ^ ((hash as u32) >> 16) as i32) as u32 as usize) & (capacity - 1)
    };
    let mut order: Vec<usize> = (0..keys.len()).collect();
    order.sort_by_key(|&index| bucket(&keys[index]));
    order
}

// ---------------------------------------------------------------------------------------------
// Kmers
// ---------------------------------------------------------------------------------------------

fn is_canonical(kval: i64) -> bool {
    kval & (1i64 << KSIZE) == 0
}

fn initial_call(kval: i64) -> i32 {
    ((kval >> (KSIZE * 2 - 2)) & 3) as i32
}

fn final_call(kval: i64) -> i32 {
    (kval & 3) as i32
}

fn predecessor_val(kval: i64, call: i32) -> i64 {
    (kval >> 2) | ((call as i64) << (2 * (KSIZE - 1)))
}

fn successor_val(kval: i64, call: i32) -> i64 {
    ((kval << 2) & KMASK) | call as i64
}

/// `KmerAdjacency.reverseComplement`, a byte at a time through the same table.
pub fn reverse_complement_kval(mut val: i64) -> i64 {
    let table = |byte: i64| -> i64 {
        let b = byte & 0xff;
        !(((b & 3) << 6) | (((b >> 2) & 3) << 4) | (((b >> 4) & 3) << 2) | ((b >> 6) & 3)) & 0xff
    };
    let mut result = table(val);
    for _ in 0..7 {
        val >>= 8;
        result = (result << 8) | table(val);
    }
    ((result as u64) >> (64 - 2 * KSIZE)) as i64
}

fn kmer_string(kval: i64) -> Vec<u8> {
    let mut bases = Vec::with_capacity(KSIZE);
    let mut current = kval;
    for _ in 0..KSIZE {
        bases.push(b"ACGT"[(current & 3) as usize]);
        current >>= 2;
    }
    bases.reverse();
    bases
}

const COUNT_FOR_MASK: [i32; 16] = [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4];
const NIBREV: [u8; 16] = [
    0b0000, 0b1000, 0b0100, 0b1100, 0b0010, 0b1010, 0b0110, 0b1110, 0b0001, 0b1001, 0b0101, 0b1101,
    0b0011, 0b1011, 0b0111, 0b1111,
];

/// A `KmerAdjacency`: a canonical `KmerAdjacencyImpl`, or (`rc`) its `KmerAdjacencyRC`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct Adj {
    index: u32,
    rc: bool,
}

impl Adj {
    fn rc(self) -> Adj {
        Adj {
            index: self.index,
            rc: !self.rc,
        }
    }
}

/// A `Contig`: a `ContigImpl`, or (`rc`) its `ContigRCImpl`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct Ctg {
    index: u32,
    rc: bool,
}

impl Ctg {
    fn rc(self) -> Ctg {
        Ctg {
            index: self.index,
            rc: !self.rc,
        }
    }
    fn canonical(self) -> Ctg {
        Ctg {
            index: self.index,
            rc: false,
        }
    }
}

struct AdjData {
    kval: i64,
    rc_kval: i64,
    sole_predecessor: Option<Adj>,
    sole_successor: Option<Adj>,
    predecessor_mask: u8,
    successor_mask: u8,
    observations: i32,
    contig: Option<Ctg>,
    contig_offset: i32,
}

struct ContigData {
    id: i32,
    sequence: Vec<u8>,
    max_observations: i32,
    first: Adj,
    last: Adj,
    predecessors: Vec<Ctg>,
    successors: Vec<Ctg>,
    cyclic: bool,
    marked: bool,
}

/// An entry of the kmer set: the canonical kmer's value and its adjacency.
#[derive(Clone)]
struct KmerEntry {
    kval: i64,
    adj: u32,
}

fn kmer_entry_index(entry: &KmerEntry, capacity: i32) -> i32 {
    kmer_index(entry.kval, capacity)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Orientation {
    Fwd,
    Rev,
    Both,
}

#[derive(Clone)]
struct ContigEnd {
    kval: i64,
    contig: Ctg,
    orientation: Orientation,
}

fn contig_end_index(entry: &ContigEnd, capacity: i32) -> i32 {
    kmer_index(entry.kval, capacity)
}

/// Which of a contig's two lists a view reads, and whether it reflects it (`ContigListRC`).
#[derive(Clone, Copy)]
struct ListView {
    contig: u32,
    successors: bool,
    reflected: bool,
}

/// `PathPart`.
#[derive(Clone, Debug)]
enum PathPart {
    Gap(Vec<u8>),
    Contig { contig: Ctg, start: i32, stop: i32 },
}

impl PathPart {
    fn contig(&self) -> Option<Ctg> {
        match self {
            PathPart::Gap(_) => None,
            PathPart::Contig { contig, .. } => Some(*contig),
        }
    }
    fn is_gap(&self) -> bool {
        matches!(self, PathPart::Gap(_))
    }
}

/// `TransitPairCount`, one of a pair sharing its count with its reverse complement.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Tpc {
    pair: u32,
    rc: bool,
}

struct TpcData {
    prev: Ctg,
    next: Ctg,
    count: i32,
}

type Traversal = Vec<Ctg>;

/// `TraversalSet`: refuses a traversal whose reverse complement it holds, and explodes past its
/// limit. Its order is this port's insertion order (see the module note).
struct TraversalSet {
    items: Vec<Traversal>,
    members: HashSet<Traversal>,
    too_many: i32,
}

struct TooComplex;

impl TraversalSet {
    fn new(too_many: i32) -> Self {
        TraversalSet {
            items: Vec::new(),
            members: HashSet::new(),
            too_many,
        }
    }
}

/// The settings `LocalAssembler` declares.
#[derive(Clone, Debug)]
pub struct Settings {
    pub assembly_name: String,
    pub q_min: i8,
    pub min_thin_observations: i32,
    pub min_gapfill_count: i32,
    pub too_many_traversals: i32,
    pub too_many_scaffolds: i32,
    pub min_sv_size: i32,
    pub no_scaffolding: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            assembly_name: String::new(),
            q_min: 25,
            min_thin_observations: 4,
            min_gapfill_count: 3,
            too_many_traversals: 100_000,
            too_many_scaffolds: 50_000,
            min_sv_size: 50,
            no_scaffolding: false,
        }
    }
}

/// What the tool writes: the GFA (absent when the traversal was too complex) and the FASTA.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Assembly {
    pub gfa: Option<String>,
    pub fasta: String,
}

/// A read as the assembler keeps it: its calls and their qualities.
#[derive(Clone)]
pub struct AssemblyRead {
    pub bases: Vec<u8>,
    pub qualities: Vec<u8>,
}

struct Graph {
    adjs: Vec<AdjData>,
    contigs: Vec<ContigData>,
    kmers: Hopscotch<KmerEntry>,
    tpcs: Vec<TpcData>,
}

impl Graph {
    // ----- adjacency views -----

    fn kval(&self, a: Adj) -> i64 {
        let data = &self.adjs[a.index as usize];
        if a.rc {
            data.rc_kval
        } else {
            data.kval
        }
    }

    fn sole_predecessor(&self, a: Adj) -> Option<Adj> {
        let data = &self.adjs[a.index as usize];
        if a.rc {
            data.sole_successor.map(Adj::rc)
        } else {
            data.sole_predecessor
        }
    }

    fn sole_successor(&self, a: Adj) -> Option<Adj> {
        let data = &self.adjs[a.index as usize];
        if a.rc {
            data.sole_predecessor.map(Adj::rc)
        } else {
            data.sole_successor
        }
    }

    fn predecessor_mask(&self, a: Adj) -> u8 {
        let data = &self.adjs[a.index as usize];
        if a.rc {
            NIBREV[data.successor_mask as usize]
        } else {
            data.predecessor_mask
        }
    }

    fn successor_mask(&self, a: Adj) -> u8 {
        let data = &self.adjs[a.index as usize];
        if a.rc {
            NIBREV[data.predecessor_mask as usize]
        } else {
            data.successor_mask
        }
    }

    fn predecessor_count(&self, a: Adj) -> i32 {
        COUNT_FOR_MASK[self.predecessor_mask(a) as usize]
    }

    fn successor_count(&self, a: Adj) -> i32 {
        COUNT_FOR_MASK[self.successor_mask(a) as usize]
    }

    fn observations(&self, a: Adj) -> i32 {
        self.adjs[a.index as usize].observations
    }

    fn adj_contig(&self, a: Adj) -> Option<Ctg> {
        let contig = self.adjs[a.index as usize].contig;
        if a.rc {
            contig.map(Ctg::rc)
        } else {
            contig
        }
    }

    fn adj_contig_offset(&self, a: Adj) -> i32 {
        let data = &self.adjs[a.index as usize];
        if !a.rc {
            return data.contig_offset;
        }
        match data.contig {
            None => 0,
            Some(contig) => self.size(contig) - data.contig_offset - KSIZE as i32,
        }
    }

    fn set_contig_offset(&mut self, a: Adj, contig: Ctg, offset: i32) -> Result<(), AssemblyError> {
        let (contig, offset) = if a.rc {
            (contig.rc(), self.size(contig) - offset - KSIZE as i32)
        } else {
            (contig, offset)
        };
        let data = &mut self.adjs[a.index as usize];
        if data.contig.is_some() {
            return Err(AssemblyError::gatk(
                "Internal error: overwriting kmer contig and offset.",
            ));
        }
        data.contig = Some(contig);
        data.contig_offset = offset;
        Ok(())
    }

    fn clear_adj_contig(&mut self, a: Adj) {
        let data = &mut self.adjs[a.index as usize];
        data.contig = None;
        data.contig_offset = 0;
    }

    fn find_canonical(&self, kval: i64) -> Option<u32> {
        let key = kval;
        self.kmers
            .find(kmer_index(key, self.kmers.capacity()), |entry| {
                entry.kval == key
            })
            .map(|entry| entry.adj)
    }

    /// `KmerAdjacency.find`.
    fn find(&self, kval: i64) -> Option<Adj> {
        if is_canonical(kval) {
            self.find_canonical(kval & KMASK)
                .map(|index| Adj { index, rc: false })
        } else {
            self.find_canonical(reverse_complement_kval(kval))
                .map(|index| Adj { index, rc: true })
        }
    }

    /// `KmerAdjacency.findOrAdd`.
    fn find_or_add(&mut self, kval: i64) -> Adj {
        let (key, rc) = if is_canonical(kval) {
            (kval & KMASK, false)
        } else {
            (reverse_complement_kval(kval), true)
        };
        let next_index = self.adjs.len() as u32;
        let entry = self.kmers.find_or_add(
            |capacity| kmer_index(key, capacity),
            |entry| entry.kval == key,
            || KmerEntry {
                kval: key,
                adj: next_index,
            },
        );
        if entry.adj == next_index && self.adjs.len() as u32 == next_index {
            self.adjs.push(AdjData {
                kval: key,
                rc_kval: reverse_complement_kval(key),
                sole_predecessor: None,
                sole_successor: None,
                predecessor_mask: 0,
                successor_mask: 0,
                observations: 0,
                contig: None,
                contig_offset: 0,
            });
        }
        Adj {
            index: entry.adj,
            rc,
        }
    }

    fn observe(&mut self, a: Adj, predecessor: Option<Adj>, successor: Option<Adj>, count: i32) {
        if a.rc {
            return self.observe(
                a.rc(),
                successor.map(Adj::rc),
                predecessor.map(Adj::rc),
                count,
            );
        }
        if let Some(predecessor) = predecessor {
            let initial = initial_call(self.kval(predecessor));
            let bit = 1u8 << initial;
            let data = &mut self.adjs[a.index as usize];
            if bit & data.predecessor_mask == 0 {
                if data.predecessor_mask == 0 {
                    data.sole_predecessor = Some(predecessor);
                    data.predecessor_mask = bit;
                } else {
                    data.sole_predecessor = None;
                    data.predecessor_mask |= bit;
                }
            }
        }
        if let Some(successor) = successor {
            let final_ = final_call(self.kval(successor));
            let bit = 1u8 << final_;
            let data = &mut self.adjs[a.index as usize];
            if bit & data.successor_mask == 0 {
                if data.successor_mask == 0 {
                    data.sole_successor = Some(successor);
                    data.successor_mask = bit;
                } else {
                    data.sole_successor = None;
                    data.successor_mask |= bit;
                }
            }
        }
        self.adjs[a.index as usize].observations += count;
    }

    fn remove_predecessor(&mut self, a: Adj, call: i32) {
        if a.rc {
            return self.remove_successor(a.rc(), 3 - call);
        }
        let data = &mut self.adjs[a.index as usize];
        data.predecessor_mask &= !(1u8 << call);
        data.sole_predecessor = None;
        if COUNT_FOR_MASK[data.predecessor_mask as usize] == 1 {
            let mask = data.predecessor_mask;
            let kval = data.kval;
            for c in 0..4 {
                if (1u8 << c) & mask != 0 {
                    let found = self.find(predecessor_val(kval, c));
                    self.adjs[a.index as usize].sole_predecessor = found;
                    break;
                }
            }
        }
    }

    fn remove_successor(&mut self, a: Adj, call: i32) {
        if a.rc {
            return self.remove_predecessor(a.rc(), 3 - call);
        }
        let data = &mut self.adjs[a.index as usize];
        data.successor_mask &= !(1u8 << call);
        data.sole_successor = None;
        if COUNT_FOR_MASK[data.successor_mask as usize] == 1 {
            let mask = data.successor_mask;
            let kval = data.kval;
            for c in 0..4 {
                if (1u8 << c) & mask != 0 {
                    let found = self.find(successor_val(kval, c));
                    self.adjs[a.index as usize].sole_successor = found;
                    break;
                }
            }
        }
    }

    /// `KmerAdjacency.kmerize(calls, quals, qMin, set)`.
    fn kmerize_read(&mut self, calls: &[u8], quals: &[u8], q_min: i8) -> Result<(), AssemblyError> {
        let mut count = 0usize;
        let mut kval: i64 = 0;
        let mut prev: Option<Adj> = None;
        let mut current: Option<Adj> = None;
        for (index, &call) in calls.iter().enumerate() {
            let quality = *quals.get(index).ok_or_else(|| AssemblyError {
                class: "java.lang.ArrayIndexOutOfBoundsException",
                message: format!("Index {index} out of bounds for length {}", quals.len()),
            })? as i8;
            if quality < q_min {
                if let Some(current) = current {
                    self.observe(current, prev, None, 1);
                }
                count = 0;
                current = None;
                prev = None;
                continue;
            }
            kval <<= 2;
            match call {
                b'A' | b'a' => {}
                b'C' | b'c' => kval += 1,
                b'G' | b'g' => kval += 2,
                b'T' | b't' => kval += 3,
                _ => {
                    if let Some(current) = current {
                        self.observe(current, prev, None, 1);
                    }
                    count = 0;
                    current = None;
                    prev = None;
                    continue;
                }
            }
            count += 1;
            if count >= KSIZE {
                let next = self.find_or_add(kval);
                if let Some(current) = current {
                    self.observe(current, prev, Some(next), 1);
                }
                prev = current;
                current = Some(next);
            }
        }
        if let Some(current) = current {
            self.observe(current, prev, None, 1);
        }
        Ok(())
    }

    /// `KmerAdjacency.kmerize(sequence, nObservations, set)`, for gap fills.
    fn kmerize_fill(&mut self, sequence: &[u8], observations: i32) -> Result<(), AssemblyError> {
        let mut count = 0usize;
        let mut kval: i64 = 0;
        let mut obs = 0;
        let mut prev: Option<Adj> = None;
        let mut current: Option<Adj> = None;
        for &call in sequence {
            kval <<= 2;
            match call {
                b'A' | b'a' => {}
                b'C' | b'c' => kval += 1,
                b'G' | b'g' => kval += 2,
                b'T' | b't' => kval += 3,
                _ => {
                    return Err(AssemblyError::gatk(
                        "unexpected base call in string to kmerize.",
                    ));
                }
            }
            count += 1;
            if count >= KSIZE {
                let next = self.find_or_add(kval);
                if let Some(current) = current {
                    self.observe(current, prev, Some(next), obs);
                    obs = observations;
                }
                prev = current;
                current = Some(next);
            }
        }
        if let Some(current) = current {
            self.observe(current, prev, None, 0);
        }
        Ok(())
    }

    // ----- contig views -----

    fn contig_data(&self, c: Ctg) -> &ContigData {
        &self.contigs[c.index as usize]
    }

    fn id(&self, c: Ctg) -> i32 {
        let id = self.contig_data(c).id;
        if c.rc {
            !id
        } else {
            id
        }
    }

    fn size(&self, c: Ctg) -> i32 {
        self.contig_data(c).sequence.len() as i32
    }

    fn n_kmers(&self, c: Ctg) -> i32 {
        self.size(c) - KSIZE as i32 + 1
    }

    fn sequence(&self, c: Ctg) -> Vec<u8> {
        let sequence = &self.contig_data(c).sequence;
        if c.rc {
            reverse_complement_bases(sequence)
        } else {
            sequence.clone()
        }
    }

    fn base_at(&self, c: Ctg, index: i32) -> u8 {
        let sequence = &self.contig_data(c).sequence;
        if c.rc {
            complement(sequence[sequence.len() - 1 - index as usize])
        } else {
            sequence[index as usize]
        }
    }

    fn max_observations(&self, c: Ctg) -> i32 {
        self.contig_data(c).max_observations
    }

    fn first_kmer(&self, c: Ctg) -> Adj {
        let data = self.contig_data(c);
        if c.rc {
            data.last.rc()
        } else {
            data.first
        }
    }

    fn last_kmer(&self, c: Ctg) -> Adj {
        let data = self.contig_data(c);
        if c.rc {
            data.first.rc()
        } else {
            data.last
        }
    }

    fn is_cycle_member(&self, c: Ctg) -> bool {
        self.contig_data(c).cyclic
    }

    fn set_cycle_member(&mut self, c: Ctg, value: bool) {
        self.contigs[c.index as usize].cyclic = value;
    }

    fn is_marked(&self, c: Ctg) -> bool {
        self.contig_data(c).marked
    }

    fn set_marked(&mut self, c: Ctg, value: bool) {
        self.contigs[c.index as usize].marked = value;
    }

    fn name(&self, c: Ctg) -> String {
        let id = self.contig_data(c).id;
        if c.rc {
            format!("c{id}RC")
        } else {
            format!("c{id}")
        }
    }

    fn reference(&self, c: Ctg) -> String {
        let id = self.contig_data(c).id;
        format!("c{id}{}", if c.rc { "-" } else { "+" })
    }

    fn predecessors(c: Ctg) -> ListView {
        ListView {
            contig: c.index,
            successors: c.rc,
            reflected: c.rc,
        }
    }

    fn successors(c: Ctg) -> ListView {
        ListView {
            contig: c.index,
            successors: !c.rc,
            reflected: c.rc,
        }
    }

    fn raw(&self, view: ListView) -> &Vec<Ctg> {
        let data = &self.contigs[view.contig as usize];
        if view.successors {
            &data.successors
        } else {
            &data.predecessors
        }
    }

    fn raw_mut(&mut self, view: ListView) -> &mut Vec<Ctg> {
        let data = &mut self.contigs[view.contig as usize];
        if view.successors {
            &mut data.successors
        } else {
            &mut data.predecessors
        }
    }

    fn list(&self, view: ListView) -> Vec<Ctg> {
        let raw = self.raw(view);
        if view.reflected {
            raw.iter().rev().map(|c| c.rc()).collect()
        } else {
            raw.clone()
        }
    }

    fn list_len(&self, view: ListView) -> usize {
        self.raw(view).len()
    }

    fn list_contains(&self, view: ListView, contig: Ctg) -> bool {
        self.list(view).contains(&contig)
    }

    fn list_index_of(&self, view: ListView, contig: Ctg) -> Option<usize> {
        self.list(view).iter().position(|c| *c == contig)
    }

    fn list_set(&mut self, view: ListView, index: usize, contig: Ctg) {
        let reflected = view.reflected;
        let raw = self.raw_mut(view);
        let len = raw.len();
        if reflected {
            raw[len - 1 - index] = contig.rc();
        } else {
            raw[index] = contig;
        }
    }

    fn list_remove(&mut self, view: ListView, contig: Ctg) -> bool {
        match self.list_index_of(view, contig) {
            None => false,
            Some(index) => {
                let reflected = view.reflected;
                let raw = self.raw_mut(view);
                let len = raw.len();
                raw.remove(if reflected { len - 1 - index } else { index });
                true
            }
        }
    }

    fn list_push(&mut self, view: ListView, contig: Ctg) {
        // Only ever called on a ContigImpl's own lists.
        debug_assert!(!view.reflected);
        self.raw_mut(view).push(contig);
    }

    // ----- contigs -----

    /// `new ContigImpl(id, firstKmerAdjacency)`.
    fn new_contig(&mut self, id: i32, first: Adj) -> Ctg {
        let mut sequence = kmer_string(self.kval(first));
        let mut max_observations = self.observations(first);
        let mut last = first;
        let mut next = self.sole_successor(first);
        while let Some(kmer) = next {
            if kmer == first || self.predecessor_count(kmer) != 1 || kmer == last.rc() {
                break;
            }
            sequence.push(b"ACGT"[final_call(self.kval(kmer)) as usize]);
            max_observations = max_observations.max(self.observations(kmer));
            last = kmer;
            next = self.sole_successor(kmer);
        }
        let index = self.contigs.len() as u32;
        self.contigs.push(ContigData {
            id,
            sequence,
            max_observations,
            first,
            last,
            predecessors: Vec::new(),
            successors: Vec::new(),
            cyclic: false,
            marked: false,
        });
        Ctg { index, rc: false }
    }

    /// `new ContigImpl(id, predecessor, successor)`.
    fn join(&mut self, id: i32, predecessor: Ctg, successor: Ctg) -> Result<Ctg, AssemblyError> {
        if predecessor == successor || predecessor == successor.rc() {
            return Err(AssemblyError::gatk("can't self-join"));
        }
        let first_sequence = self.sequence(predecessor);
        let second_sequence = self.sequence(successor);
        if first_sequence[first_sequence.len() - KSIZE + 1..] != second_sequence[..KSIZE - 1] {
            return Err(AssemblyError::gatk("sequences can't be joined"));
        }
        let mut sequence = first_sequence;
        sequence.extend_from_slice(&second_sequence[KSIZE - 1..]);
        let index = self.contigs.len() as u32;
        let joined = Ctg { index, rc: false };
        self.contigs.push(ContigData {
            id,
            sequence,
            max_observations: self
                .max_observations(predecessor)
                .max(self.max_observations(successor)),
            first: self.first_kmer(predecessor),
            last: self.last_kmer(successor),
            predecessors: Vec::new(),
            successors: Vec::new(),
            cyclic: false,
            marked: false,
        });
        for pred_predecessor in self.list(Graph::predecessors(predecessor)) {
            if pred_predecessor == successor {
                self.list_push(Graph::predecessors(joined), joined);
            } else if pred_predecessor == predecessor.rc() {
                self.list_push(Graph::predecessors(joined), joined.rc());
            } else {
                self.list_push(Graph::predecessors(joined), pred_predecessor);
                let view = Graph::successors(pred_predecessor);
                let at = self
                    .list_index_of(view, predecessor)
                    .ok_or_else(|| AssemblyError {
                        class: "java.lang.ArrayIndexOutOfBoundsException",
                        message: "Index -1 out of bounds".to_string(),
                    })?;
                self.list_set(view, at, joined);
            }
        }
        for succ_successor in self.list(Graph::successors(successor)) {
            if succ_successor == predecessor {
                self.list_push(Graph::successors(joined), joined);
            } else if succ_successor == successor.rc() {
                self.list_push(Graph::successors(joined), joined.rc());
            } else {
                self.list_push(Graph::successors(joined), succ_successor);
                let view = Graph::predecessors(succ_successor);
                let at = self
                    .list_index_of(view, successor)
                    .ok_or_else(|| AssemblyError {
                        class: "java.lang.ArrayIndexOutOfBoundsException",
                        message: "Index -1 out of bounds".to_string(),
                    })?;
                self.list_set(view, at, joined);
            }
        }
        self.clear_kmer_contig(joined)?;
        self.set_kmer_contig(joined)?;
        Ok(joined)
    }

    fn clear_kmer_contig(&mut self, contig: Ctg) -> Result<(), AssemblyError> {
        let mut count = 0;
        let first = self.first_kmer(contig);
        let last = self.last_kmer(contig);
        let mut kmer = Some(first);
        while kmer != Some(last) {
            let Some(current) = kmer else {
                return Err(AssemblyError::gatk(
                    "contig does not have a flat pipeline of kmers",
                ));
            };
            if self.adj_contig(current).is_none() {
                return Err(AssemblyError::gatk(
                    "we've returned to a kmer we've already cleared",
                ));
            }
            self.clear_adj_contig(current);
            count += 1;
            kmer = self.sole_successor(current);
        }
        self.clear_adj_contig(last);
        if count + KSIZE as i32 != self.size(contig) {
            return Err(AssemblyError::gatk(
                "kmer chain length does not equal contig size",
            ));
        }
        Ok(())
    }

    fn set_kmer_contig(&mut self, contig: Ctg) -> Result<(), AssemblyError> {
        let mut offset = 0;
        let first = self.first_kmer(contig);
        let last = self.last_kmer(contig);
        let mut kmer = Some(first);
        while kmer != Some(last) {
            let Some(current) = kmer else {
                return Err(AssemblyError::gatk(
                    "contig does not have a flat pipeline of kmers",
                ));
            };
            if self.adj_contig(current).is_some() {
                return Err(AssemblyError::gatk(
                    "we've returned to a kmer we've already updated",
                ));
            }
            self.set_contig_offset(current, contig, offset)?;
            offset += 1;
            kmer = self.sole_successor(current);
        }
        self.set_contig_offset(last, contig, offset)?;
        if offset + KSIZE as i32 != self.size(contig) {
            return Err(AssemblyError::gatk(
                "kmer chain length does not equal contig size",
            ));
        }
        Ok(())
    }

    /// `buildContigs`.
    fn build_contigs(&mut self) -> Result<Vec<Ctg>, AssemblyError> {
        let mut contigs = Vec::new();
        let mut n_contigs = 0;
        let entries = self.kmers.entries();
        for entry in &entries {
            let kmer = Adj {
                index: entry.adj,
                rc: false,
            };
            if self.adj_contig(kmer).is_some() {
                continue;
            }
            let mut contig = None;
            let predecessor = self.sole_predecessor(kmer);
            let starts = match predecessor {
                None => true,
                Some(p) => self.successor_count(p) > 1 || p == kmer.rc(),
            };
            if starts {
                n_contigs += 1;
                contig = Some(self.new_contig(n_contigs, kmer));
            } else {
                let successor = self.sole_successor(kmer);
                let ends = match successor {
                    None => true,
                    Some(s) => self.predecessor_count(s) > 1 || s == kmer.rc(),
                };
                if ends {
                    n_contigs += 1;
                    contig = Some(self.new_contig(n_contigs, kmer.rc()));
                }
            }
            if let Some(contig) = contig {
                self.set_kmer_contig(contig)?;
                contigs.push(contig);
            }
        }
        for entry in &entries {
            let kmer = Adj {
                index: entry.adj,
                rc: false,
            };
            if self.adj_contig(kmer).is_none() {
                n_contigs += 1;
                let contig = self.new_contig(n_contigs, kmer);
                self.set_kmer_contig(contig)?;
                contigs.push(contig);
            }
        }
        Ok(contigs)
    }

    /// `connectContigs`.
    fn connect_contigs(&mut self, contigs: &[Ctg]) -> Result<(), AssemblyError> {
        let mut ends: Hopscotch<ContigEnd> =
            Hopscotch::new(2 * contigs.len() as i32, contig_end_index);
        let same = |a: &ContigEnd, b: &ContigEnd| a.kval == b.kval;
        for &contig in contigs {
            let fwd = self.first_kmer(contig);
            let rev = self.last_kmer(contig).rc();
            if fwd == rev {
                ends.add(
                    ContigEnd {
                        kval: self.kval(fwd),
                        contig,
                        orientation: Orientation::Both,
                    },
                    same,
                );
            } else {
                ends.add(
                    ContigEnd {
                        kval: self.kval(fwd),
                        contig,
                        orientation: Orientation::Fwd,
                    },
                    same,
                );
                ends.add(
                    ContigEnd {
                        kval: self.kval(rev),
                        contig,
                        orientation: Orientation::Rev,
                    },
                    same,
                );
            }
        }
        let find = |ends: &Hopscotch<ContigEnd>, kval: i64| {
            ends.find(kmer_index(kval, ends.capacity()), |entry| {
                entry.kval == kval
            })
            .cloned()
            .ok_or_else(|| AssemblyError::gatk("missing contig end kmer"))
        };
        for &contig in contigs {
            let start = self.first_kmer(contig);
            if self.predecessor_count(start) > 0 {
                let mask = self.predecessor_mask(start);
                for call in 0..4 {
                    if mask & (1 << call) != 0 {
                        let kval = reverse_complement_kval(predecessor_val(self.kval(start), call));
                        let end = find(&ends, kval)?;
                        let view = Graph::predecessors(contig);
                        match end.orientation {
                            Orientation::Fwd => self.list_push(view, end.contig.rc()),
                            Orientation::Rev => self.list_push(view, end.contig),
                            Orientation::Both => {
                                self.list_push(view, end.contig);
                                self.list_push(view, end.contig.rc());
                            }
                        }
                    }
                }
            }
            let last = self.last_kmer(contig);
            if self.successor_count(last) > 0 {
                let mask = self.successor_mask(last);
                for call in 0..4 {
                    if mask & (1 << call) != 0 {
                        let kval = successor_val(self.kval(last), call);
                        let end = find(&ends, kval)?;
                        let view = Graph::successors(contig);
                        match end.orientation {
                            Orientation::Fwd => self.list_push(view, end.contig),
                            Orientation::Rev => self.list_push(view, end.contig.rc()),
                            Orientation::Both => {
                                self.list_push(view, end.contig);
                                self.list_push(view, end.contig.rc());
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// `removeThinContigs`.
    fn remove_thin_contigs(
        &mut self,
        contigs: &mut Vec<Ctg>,
        min_thin_observations: i32,
    ) -> Result<(), AssemblyError> {
        contigs.sort_by_key(|c| self.max_observations(*c));
        let mut next_visit = 0;
        loop {
            let mut cut_data: HashMap<Ctg, (i32, i32)> = HashMap::new();
            for &contig in contigs.iter() {
                if cut_data.contains_key(&contig) {
                    continue;
                }
                next_visit += 1;
                cut_data.insert(contig, (next_visit, next_visit));
                let mut children = 0;
                for next in self.list(Graph::successors(contig)) {
                    if !cut_data.contains_key(&next) {
                        self.find_cuts(next, contig, &mut next_visit, &mut cut_data);
                        children += 1;
                    }
                }
                for next in self.list(Graph::predecessors(contig)) {
                    if !cut_data.contains_key(&next) {
                        self.find_cuts(next, contig, &mut next_visit, &mut cut_data);
                        children += 1;
                    }
                }
                if children >= 2 {
                    self.set_marked(contig, true);
                }
            }
            let mut removed = false;
            for position in 0..contigs.len() {
                let contig = contigs[position];
                if self.max_observations(contig) < min_thin_observations && !self.is_marked(contig)
                {
                    self.unlink_contig(contig)?;
                    contigs.remove(position);
                    removed = true;
                    break;
                }
            }
            if !removed {
                break;
            }
        }
        contigs.sort_by_key(|c| self.id(*c));
        Ok(())
    }

    fn find_cuts(
        &mut self,
        contig: Ctg,
        parent: Ctg,
        next_visit: &mut i32,
        cut_data: &mut HashMap<Ctg, (i32, i32)>,
    ) -> (i32, i32) {
        *next_visit += 1;
        let visit = *next_visit;
        let mut min_visit = visit;
        cut_data.insert(contig, (visit, visit));
        let neighbours: Vec<Ctg> = self
            .list(Graph::successors(contig))
            .into_iter()
            .chain(self.list(Graph::predecessors(contig)))
            .collect();
        for next in neighbours {
            if next == parent {
                continue;
            }
            match cut_data.get(&next).copied() {
                Some((next_visit_num, _)) => {
                    min_visit = min_visit.min(next_visit_num);
                    cut_data.insert(contig, (visit, min_visit));
                }
                None => {
                    let (_, next_min) = self.find_cuts(next, contig, next_visit, cut_data);
                    min_visit = min_visit.min(next_min);
                    cut_data.insert(contig, (visit, min_visit));
                    if next_min >= visit {
                        self.set_marked(contig, true);
                    }
                }
            }
        }
        (visit, min_visit)
    }

    /// `unlinkContig`.
    fn unlink_contig(&mut self, contig: Ctg) -> Result<(), AssemblyError> {
        let first = self.first_kmer(contig);
        let first_final = final_call(self.kval(first));
        for predecessor in self.list(Graph::predecessors(contig)) {
            if predecessor != contig && predecessor != contig.rc() {
                let last = self.last_kmer(predecessor);
                self.remove_successor(last, first_final);
                if !self.list_remove(Graph::successors(predecessor), contig) {
                    return Err(AssemblyError::gatk("failed to find predecessor link"));
                }
            }
        }
        let last = self.last_kmer(contig);
        let last_initial = initial_call(self.kval(last));
        for successor in self.list(Graph::successors(contig)) {
            if successor != contig && successor != contig.rc() {
                let first = self.first_kmer(successor);
                self.remove_predecessor(first, last_initial);
                if !self.list_remove(Graph::predecessors(successor), contig) {
                    return Err(AssemblyError::gatk("failed to find successor link"));
                }
            }
        }
        let mut next = Some(first);
        loop {
            let Some(kmer) = next else {
                return Err(AssemblyError {
                    class: "java.lang.NullPointerException",
                    message: "null".to_string(),
                });
            };
            next = self.sole_successor(kmer);
            let key = self.adjs[kmer.index as usize].kval;
            self.kmers
                .remove(kmer_index(key, self.kmers.capacity()), |entry| {
                    entry.kval == key
                });
            if kmer == last {
                break;
            }
        }
        Ok(())
    }

    /// `weldPipes`.
    fn weld_pipes(&mut self, contigs: &mut Vec<Ctg>) -> Result<(), AssemblyError> {
        let mut index: isize = 0;
        while (index as usize) < contigs.len() {
            let contig = contigs[index as usize];
            let successors = self.list(Graph::successors(contig));
            if successors.len() == 1 {
                let successor = successors[0];
                if successor != contig
                    && successor != contig.rc()
                    && self.list_len(Graph::predecessors(successor)) == 1
                {
                    let joined = self.join(self.id(contig), contig, successor)?;
                    contigs[index as usize] = joined;
                    let target = successor.canonical();
                    match contigs.iter().position(|c| *c == target) {
                        Some(at) => {
                            contigs.remove(at);
                        }
                        None => return Err(AssemblyError::gatk("successor linkage is messed up")),
                    }
                    // `contigIdx -= 1; continue;`: the loop's increment brings it back, so the
                    // joined contig is considered again.
                    continue;
                }
            }
            let predecessors = self.list(Graph::predecessors(contig));
            if predecessors.len() == 1 {
                let predecessor = predecessors[0];
                if predecessor != contig
                    && predecessor != contig.rc()
                    && self.list_len(Graph::successors(predecessor)) == 1
                {
                    let joined = self.join(self.id(contig), predecessor, contig)?;
                    contigs[index as usize] = joined;
                    let target = predecessor.canonical();
                    match contigs.iter().position(|c| *c == target) {
                        Some(at) => {
                            contigs.remove(at);
                        }
                        None => {
                            return Err(AssemblyError::gatk("predecessor linkage is messed up"));
                        }
                    }
                    index -= 1;
                }
            }
            index += 1;
        }
        Ok(())
    }

    /// `markCycles`.
    fn mark_cycles(&mut self, contigs: &[Ctg]) {
        for &contig in contigs {
            self.set_cycle_member(contig, false);
        }
        let mut deque: Vec<Ctg> = Vec::new();
        let mut walk: HashMap<Ctg, (i32, i32)> = HashMap::new();
        let mut next_visit = 0;
        for &contig in contigs {
            if !walk.contains_key(&contig) {
                self.mark_cycles_recursion(contig, &mut deque, &mut next_visit, &mut walk);
            }
        }
    }

    fn mark_cycles_recursion(
        &mut self,
        contig: Ctg,
        deque: &mut Vec<Ctg>,
        next_visit: &mut i32,
        walk: &mut HashMap<Ctg, (i32, i32)>,
    ) -> i32 {
        *next_visit += 1;
        let visit = *next_visit;
        let mut min_visit = visit;
        walk.insert(contig, (visit, visit));
        deque.push(contig);
        for successor in self.list(Graph::successors(contig)) {
            match walk.get(&successor).copied() {
                None => {
                    let recursion = self.mark_cycles_recursion(successor, deque, next_visit, walk);
                    min_visit = min_visit.min(recursion);
                }
                Some((successor_visit, _)) => {
                    min_visit = min_visit.min(successor_visit);
                }
            }
            let (own_visit, _) = walk[&contig];
            walk.insert(contig, (own_visit, min_visit));
        }
        let (own_visit, _) = walk[&contig];
        if own_visit == min_visit {
            let mut tig = deque.pop().expect("a deque entry");
            if tig == contig {
                let entry = walk.get_mut(&tig).expect("walk data");
                entry.0 = i32::MAX;
                if self.list_contains(Graph::successors(tig), tig) {
                    self.set_cycle_member(tig, true);
                }
            } else {
                loop {
                    walk.get_mut(&tig).expect("walk data").0 = i32::MAX;
                    self.set_cycle_member(tig, true);
                    if tig == contig {
                        break;
                    }
                    tig = deque.pop().expect("a deque entry");
                }
            }
        }
        min_visit
    }

    // ----- paths -----

    fn gap_last_call(sequence: &[u8]) -> u8 {
        sequence[sequence.len() - KSIZE + 1]
    }

    fn part_first_call(&self, part: &PathPart) -> u8 {
        match part {
            PathPart::Gap(sequence) => sequence[KSIZE - 1],
            PathPart::Contig { contig, start, .. } => {
                self.base_at(*contig, start + KSIZE as i32 - 1)
            }
        }
    }

    fn part_last_call(&self, part: &PathPart) -> u8 {
        match part {
            PathPart::Gap(sequence) => Graph::gap_last_call(sequence),
            PathPart::Contig { contig, stop, .. } => self.base_at(*contig, stop - 1),
        }
    }

    fn part_rc(&self, part: &PathPart) -> PathPart {
        match part {
            PathPart::Gap(sequence) => PathPart::Gap(reverse_complement_bases(sequence)),
            PathPart::Contig {
                contig,
                start,
                stop,
            } => {
                let rev_base = self.size(*contig) - KSIZE as i32 + 1;
                PathPart::Contig {
                    contig: contig.rc(),
                    start: rev_base - stop,
                    stop: rev_base - start,
                }
            }
        }
    }

    fn path_rc(&self, path: &[PathPart]) -> Vec<PathPart> {
        path.iter().rev().map(|part| self.part_rc(part)).collect()
    }

    fn zero_length_gap(&self, current: &PathPart) -> PathPart {
        let PathPart::Contig { contig, stop, .. } = current else {
            unreachable!("a zero-length gap follows a contig part");
        };
        let sequence = self.sequence(*contig);
        PathPart::Gap(sequence[*stop as usize..*stop as usize + KSIZE - 1].to_vec())
    }

    /// `PathBuilder.processCalls`.
    fn path(&self, calls: &[u8]) -> Vec<PathPart> {
        let mut parts: Vec<PathPart> = Vec::new();
        let mut kval: i64 = 0;
        let mut count = 0usize;
        let mut current: Option<usize> = None;
        for &call in calls {
            kval <<= 2;
            match call {
                b'C' | b'c' => kval += 1,
                b'G' | b'g' => kval += 2,
                b'T' | b't' => kval += 3,
                _ => {}
            }
            count += 1;
            if count < KSIZE {
                continue;
            }
            match self.find(kval) {
                None => {
                    let extend = current.is_some_and(|index| parts[index].is_gap());
                    if extend {
                        if let PathPart::Gap(sequence) = &mut parts[current.unwrap()] {
                            sequence.push(call);
                        }
                    } else {
                        parts.push(PathPart::Gap(kmer_string(kval)));
                        current = Some(parts.len() - 1);
                    }
                }
                Some(kmer) => {
                    let contig = self.adj_contig(kmer).expect("a kmer on a contig");
                    let offset = self.adj_contig_offset(kmer);
                    match current {
                        None => {
                            parts.push(PathPart::Contig {
                                contig,
                                start: offset,
                                stop: offset + 1,
                            });
                            current = Some(parts.len() - 1);
                        }
                        Some(index) if parts[index].contig() == Some(contig) => {
                            let PathPart::Contig { stop, .. } = parts[index] else {
                                unreachable!()
                            };
                            if offset == stop {
                                if let PathPart::Contig { stop, .. } = &mut parts[index] {
                                    *stop += 1;
                                }
                            } else if offset == 0 && self.n_kmers(contig) == stop {
                                parts.push(PathPart::Contig {
                                    contig,
                                    start: 0,
                                    stop: 1,
                                });
                                current = Some(parts.len() - 1);
                            } else {
                                let gap = self.zero_length_gap(&parts[index]);
                                parts.push(gap);
                                parts.push(PathPart::Contig {
                                    contig,
                                    start: offset,
                                    stop: offset + 1,
                                });
                                current = Some(parts.len() - 1);
                            }
                        }
                        Some(index) => {
                            if let PathPart::Gap(sequence) = &parts[index] {
                                let gap_length = sequence.len() - KSIZE + 1;
                                if gap_length == KSIZE && parts.len() >= 2 {
                                    let previous = parts.len() - 2;
                                    if let Some(squashed) =
                                        self.squash(&mut parts, previous, contig, offset)
                                    {
                                        current = Some(squashed);
                                        continue;
                                    }
                                }
                            } else {
                                let PathPart::Contig {
                                    contig: part_contig,
                                    stop,
                                    ..
                                } = parts[index]
                                else {
                                    unreachable!()
                                };
                                let stops_at_end =
                                    stop + KSIZE as i32 - 1 == self.size(part_contig);
                                if !stops_at_end
                                    || offset != 0
                                    || !self.list_contains(Graph::successors(part_contig), contig)
                                {
                                    let gap = self.zero_length_gap(&parts[index]);
                                    parts.push(gap);
                                }
                            }
                            parts.push(PathPart::Contig {
                                contig,
                                start: offset,
                                stop: offset + 1,
                            });
                            current = Some(parts.len() - 1);
                        }
                    }
                }
            }
        }
        parts
    }

    /// `gapIsSquashed`: the index of the new current part when the gap was squashed.
    fn squash(
        &self,
        parts: &mut Vec<PathPart>,
        previous: usize,
        contig: Ctg,
        offset: i32,
    ) -> Option<usize> {
        let PathPart::Contig {
            contig: previous_contig,
            start: previous_start,
            stop: previous_stop,
        } = parts[previous].clone()
        else {
            // `prevPart.getContig()` of a gap is null and its size() would throw; a gap never
            // follows a gap, so this is not reached.
            return None;
        };
        let previous_max_stop = self.size(previous_contig) - KSIZE as i32 + 1;
        let new_stop = offset + 1;
        if previous_contig == contig {
            if offset - previous_stop == KSIZE as i32 {
                parts[previous] = PathPart::Contig {
                    contig: previous_contig,
                    start: previous_start,
                    stop: new_stop,
                };
                parts.remove(previous + 1);
                return Some(previous);
            }
        } else if previous_max_stop - previous_stop + offset == KSIZE as i32 {
            parts[previous] = PathPart::Contig {
                contig: previous_contig,
                start: previous_start,
                stop: previous_max_stop,
            };
            parts[previous + 1] = PathPart::Contig {
                contig,
                start: 0,
                stop: new_stop,
            };
            return Some(previous + 1);
        }
        None
    }

    /// `fillGaps`.
    fn fill_gaps(
        &mut self,
        min_gapfill_count: i32,
        reads: &[AssemblyRead],
    ) -> Result<bool, AssemblyError> {
        let mut keys: Vec<String> = Vec::new();
        let mut counts: HashMap<String, i32> = HashMap::new();
        for read in reads {
            let parts = self.path(&read.bases);
            if parts.len() < 3 {
                continue;
            }
            let last = parts.len() - 1;
            for index in 1..last {
                if let PathPart::Gap(sequence) = &parts[index] {
                    let prev_call = self.part_last_call(&parts[index - 1]);
                    let next_call = self.part_first_call(&parts[index + 1]);
                    let mut fill = vec![prev_call];
                    fill.extend_from_slice(sequence);
                    fill.push(next_call);
                    let fill_rc = reverse_complement_bases(&fill);
                    let upper: Vec<u8> = fill.iter().map(|b| b.to_ascii_uppercase()).collect();
                    let fill = if fill_rc < upper { fill_rc } else { fill };
                    let key = String::from_utf8_lossy(&fill).into_owned();
                    match counts.get_mut(&key) {
                        Some(count) => *count += 1,
                        None => {
                            counts.insert(key.clone(), 1);
                            keys.push(key);
                        }
                    }
                }
            }
        }
        let mut new_kmers = false;
        for index in java_hash_map_order(&keys) {
            let key = &keys[index];
            let observations = counts[key];
            if observations >= min_gapfill_count {
                self.kmerize_fill(key.as_bytes(), observations)?;
                new_kmers = true;
            }
        }
        if new_kmers {
            for entry in self.kmers.entries() {
                self.clear_adj_contig(Adj {
                    index: entry.adj,
                    rc: false,
                });
            }
        }
        Ok(new_kmers)
    }

    fn create_assembly(&mut self, min_thin_observations: i32) -> Result<Vec<Ctg>, AssemblyError> {
        let mut contigs = self.build_contigs()?;
        self.connect_contigs(&contigs)?;
        self.remove_thin_contigs(&mut contigs, min_thin_observations)?;
        self.weld_pipes(&mut contigs)?;
        Ok(contigs)
    }

    // ----- transits and traversals -----

    fn tpc_prev(&self, t: Tpc) -> Ctg {
        let data = &self.tpcs[t.pair as usize];
        if t.rc {
            data.next.rc()
        } else {
            data.prev
        }
    }

    fn tpc_next(&self, t: Tpc) -> Ctg {
        let data = &self.tpcs[t.pair as usize];
        if t.rc {
            data.prev.rc()
        } else {
            data.next
        }
    }

    fn tpc_count(&self, t: Tpc) -> i32 {
        self.tpcs[t.pair as usize].count
    }

    fn tpc_reset(&mut self, t: Tpc) {
        self.tpcs[t.pair as usize].count = 0;
    }

    /// `collectTransitPairCounts`, as an insertion-ordered map.
    fn collect_transits(&mut self, paths: &[Vec<PathPart>]) -> TransitMap {
        let mut map = TransitMap::default();
        for path in paths {
            if path.len() < 3 {
                continue;
            }
            let last = path.len() - 1;
            let mut index = 1;
            while index < last {
                let Some(prev) = path[index - 1].contig() else {
                    index += 1;
                    continue;
                };
                let Some(current) = path[index].contig() else {
                    index += 2;
                    continue;
                };
                let Some(next) = path[index + 1].contig() else {
                    index += 3;
                    continue;
                };
                let found = map.get(current).and_then(|list| {
                    list.iter()
                        .copied()
                        .find(|t| self.tpc_prev(*t) == prev && self.tpc_next(*t) == next)
                });
                match found {
                    Some(t) => self.tpcs[t.pair as usize].count += 1,
                    None => {
                        let pair = self.tpcs.len() as u32;
                        self.tpcs.push(TpcData {
                            prev,
                            next,
                            count: 1,
                        });
                        map.entry(current).push(Tpc { pair, rc: false });
                        map.entry(current.rc()).push(Tpc { pair, rc: true });
                    }
                }
                index += 1;
            }
        }
        map
    }

    fn add_traversal(
        &self,
        set: &mut TraversalSet,
        traversal: Traversal,
    ) -> Result<(), TooComplex> {
        let rc = traversal_rc(&traversal);
        if set.members.contains(&rc) {
            return Ok(());
        }
        if set.items.len() as i32 >= set.too_many {
            return Err(TooComplex);
        }
        if set.members.insert(traversal.clone()) {
            set.items.push(traversal);
        }
        Ok(())
    }

    /// `traverseAllPaths`.
    fn traverse_all_paths(
        &mut self,
        contigs: &[Ctg],
        paths: &[Vec<PathPart>],
        too_many: i32,
        map: &TransitMap,
    ) -> Result<Vec<Traversal>, TooComplex> {
        let mut set = TraversalSet::new(too_many);
        let mut list: Vec<Ctg> = Vec::new();
        for &contig in contigs {
            if map.get(contig).is_none() {
                let successors = self.list(Graph::successors(contig));
                let predecessors = self.list(Graph::predecessors(contig));
                if successors.is_empty() && predecessors.is_empty() {
                    self.add_traversal(&mut set, vec![contig])?;
                } else {
                    for successor in successors {
                        self.traverse(successor, contig, &mut list, paths, map, &mut set)?;
                    }
                    for predecessor in predecessors {
                        self.traverse(
                            predecessor.rc(),
                            contig.rc(),
                            &mut list,
                            paths,
                            map,
                            &mut set,
                        )?;
                    }
                }
            }
        }
        for (contig, transits) in map.iter() {
            if self.is_cycle_member(contig) {
                continue;
            }
            for &t in transits {
                if self.tpc_count(t) > 0 {
                    self.tpc_reset(t);
                    let mut fwd = TraversalSet::new(too_many);
                    self.traverse(self.tpc_next(t), contig, &mut list, paths, map, &mut fwd)?;
                    let mut rev = TraversalSet::new(too_many);
                    self.traverse(
                        self.tpc_prev(t).rc(),
                        contig.rc(),
                        &mut list,
                        paths,
                        map,
                        &mut rev,
                    )?;
                    for rev_traversal in &rev.items {
                        let rev_rc = traversal_rc(rev_traversal);
                        for fwd_traversal in &fwd.items {
                            let combined =
                                combine(&rev_rc, fwd_traversal, 1).map_err(|_| TooComplex)?;
                            self.add_traversal(&mut set, combined)?;
                        }
                    }
                }
            }
        }
        for (contig, transits) in map.iter() {
            for &t in transits {
                if self.tpc_count(t) > 0 {
                    self.tpc_reset(t);
                    let last = list.len();
                    list.push(self.tpc_prev(t));
                    self.traverse(self.tpc_next(t), contig, &mut list, paths, map, &mut set)?;
                    list[last] = self.tpc_next(t).rc();
                    self.traverse(
                        self.tpc_prev(t).rc(),
                        contig.rc(),
                        &mut list,
                        paths,
                        map,
                        &mut set,
                    )?;
                    list.remove(last);
                }
            }
        }
        Ok(set.items)
    }

    fn traverse(
        &mut self,
        contig: Ctg,
        predecessor: Ctg,
        list: &mut Vec<Ctg>,
        paths: &[Vec<PathPart>],
        map: &TransitMap,
        set: &mut TraversalSet,
    ) -> Result<(), TooComplex> {
        list.push(predecessor);
        if self.is_cycle_member(contig) {
            self.traverse_cycle(contig, list, paths, map, set)?;
            list.pop();
            return Ok(());
        }
        let mut found = false;
        if let Some(transits) = map.get(contig) {
            for &t in transits {
                if self.tpc_prev(t) == predecessor {
                    let successor = self.tpc_next(t);
                    if predecessor == contig.rc() {
                        let n = list.len();
                        if n > 1 && successor.rc() == list[n - 2] {
                            continue;
                        }
                    }
                    self.tpc_reset(t);
                    self.traverse(successor, contig, list, paths, map, set)?;
                    found = true;
                }
            }
        }
        if !found {
            list.push(contig);
            let traversal = list.clone();
            list.pop();
            self.add_traversal(set, traversal)?;
        }
        list.pop();
        Ok(())
    }

    fn traverse_cycle(
        &mut self,
        contig: Ctg,
        list: &mut Vec<Ctg>,
        paths: &[Vec<PathPart>],
        map: &TransitMap,
        set: &mut TraversalSet,
    ) -> Result<(), TooComplex> {
        list.push(contig);
        let n = list.len();
        let to_match: Vec<Ctg> = if n <= 2 {
            list.clone()
        } else {
            list[n - 2..].to_vec()
        };
        let longest = self.find_longest_paths(&to_match, paths);
        if longest.is_empty() {
            self.add_traversal(set, list.clone())?;
        } else {
            for path in longest {
                if path.is_empty() {
                    self.add_traversal(set, list.clone())?;
                    continue;
                }
                let mut extended = list.clone();
                if self.is_cycle_member(*path.last().unwrap()) {
                    extended.extend_from_slice(&path);
                    self.add_traversal(set, extended.clone())?;
                } else {
                    for &current in &path {
                        if self.is_cycle_member(current) {
                            extended.push(current);
                        } else {
                            let prev = extended.pop().expect("a previous contig");
                            self.traverse(current, prev, &mut extended, paths, map, set)?;
                            extended.push(prev);
                            break;
                        }
                    }
                }
                self.clear_transit_pairs(map, &extended);
            }
        }
        list.pop();
        Ok(())
    }

    fn clear_transit_pairs(&mut self, map: &TransitMap, list: &[Ctg]) {
        if list.len() < 3 {
            return;
        }
        let last = list.len() - 1;
        for index in 1..last {
            if let Some(transits) = map.get(list[index]) {
                let predecessor = list[index - 1];
                let successor = list[index + 1];
                for &t in transits {
                    if self.tpc_prev(t) == predecessor && self.tpc_next(t) == successor {
                        self.tpc_reset(t);
                        break;
                    }
                }
            }
        }
    }

    fn find_longest_paths(&self, to_match: &[Ctg], paths: &[Vec<PathPart>]) -> Vec<Vec<Ctg>> {
        let mut longest: Vec<Vec<Ctg>> = Vec::new();
        for path in paths {
            self.test_path(path, to_match, &mut longest);
            let rc = self.path_rc(path);
            self.test_path(&rc, to_match, &mut longest);
        }
        longest
    }

    fn test_path(&self, path: &[PathPart], to_match: &[Ctg], longest: &mut Vec<Vec<Ctg>>) {
        let contigs: Vec<Option<Ctg>> = path.iter().map(PathPart::contig).collect();
        let want: Vec<Option<Ctg>> = to_match.iter().copied().map(Some).collect();
        let found = if want.len() > contigs.len() {
            None
        } else {
            (0..=contigs.len() - want.len())
                .find(|&start| contigs[start..start + want.len()] == want[..])
        };
        if let Some(match_index) = found {
            let suffix = match_index + want.len();
            if suffix < contigs.len() {
                let sub = self.extract_sub_path(&contigs, suffix);
                add_sub_path(sub, longest);
            }
        }
    }

    fn extract_sub_path(&self, contigs: &[Option<Ctg>], suffix: usize) -> Vec<Ctg> {
        let mut prev = contigs[suffix - 1].expect("a matched contig");
        let mut sub = Vec::new();
        for tig in &contigs[suffix..] {
            match tig {
                Some(tig) if self.list_contains(Graph::successors(prev), *tig) => {
                    sub.push(*tig);
                    prev = *tig;
                }
                _ => break,
            }
        }
        sub
    }

    // ----- scaffolds -----

    fn sequence_length(&self, traversal: &[Ctg]) -> i32 {
        traversal.iter().map(|c| self.n_kmers(*c)).sum::<i32>() + KSIZE as i32 - 1
    }

    fn min_max_observations(&self, traversal: &[Ctg]) -> i32 {
        traversal
            .iter()
            .map(|c| self.max_observations(*c))
            .min()
            .unwrap_or(i32::MAX)
    }

    fn traversal_sequence(&self, traversal: &[Ctg]) -> Vec<u8> {
        let Some(first) = traversal.first() else {
            return Vec::new();
        };
        let mut sequence = self.sequence(*first)[..KSIZE - 1].to_vec();
        for contig in traversal {
            sequence.extend_from_slice(&self.sequence(*contig)[KSIZE - 1..]);
        }
        sequence
    }

    fn traversal_name(&self, traversal: &[Ctg]) -> String {
        traversal
            .iter()
            .map(|c| self.name(*c))
            .collect::<Vec<_>>()
            .join("+")
    }

    /// `TraversalEndpointComparator`.
    fn compare_endpoints(&self, first: &[Ctg], second: &[Ctg]) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        let cmp = self.id(first[0]).cmp(&self.id(second[0]));
        if cmp != Ordering::Equal {
            return cmp;
        }
        let last1 = first.len() - 1;
        let last2 = second.len() - 1;
        let cmp = self.id(first[last1]).cmp(&self.id(second[last2]));
        if cmp != Ordering::Equal {
            return cmp;
        }
        let cmp = self
            .min_max_observations(second)
            .cmp(&self.min_max_observations(first));
        if cmp != Ordering::Equal {
            return cmp;
        }
        let end = last1.min(last2);
        for index in 1..end {
            let cmp = self.id(first[index]).cmp(&self.id(second[index]));
            if cmp != Ordering::Equal {
                return cmp;
            }
        }
        last1.cmp(&last2)
    }

    fn sorted_insert(&self, sorted: &mut Vec<Traversal>, traversal: Traversal) {
        match sorted.binary_search_by(|probe| self.compare_endpoints(probe, &traversal)) {
            Ok(_) => {}
            Err(at) => sorted.insert(at, traversal),
        }
    }

    fn sorted_contains(&self, sorted: &[Traversal], traversal: &[Ctg]) -> bool {
        sorted
            .binary_search_by(|probe| self.compare_endpoints(probe, traversal))
            .is_ok()
    }

    fn remove_trivially_different(
        &self,
        traversals: &mut Vec<Traversal>,
        min_sv_size: i32,
    ) -> Result<(), AssemblyError> {
        if traversals.is_empty() {
            return Ok(());
        }
        let mut sorted: Vec<Traversal> = Vec::new();
        for traversal in traversals.iter() {
            self.sorted_insert(&mut sorted, traversal.clone());
            self.sorted_insert(&mut sorted, traversal_rc(traversal));
        }
        let mut kept: Vec<Traversal> = Vec::with_capacity(sorted.len());
        let mut iter = sorted.into_iter();
        let mut prev = iter.next().expect("one traversal");
        kept.push(prev.clone());
        for current in iter {
            if self.is_trivially_different(&prev, &current, min_sv_size)? {
                continue;
            }
            kept.push(current.clone());
            prev = current;
        }
        let mut index = 0;
        while index < kept.len() {
            let rc = traversal_rc(&kept[index]);
            if self.sorted_contains(&kept, &rc) {
                kept.remove(index);
            } else {
                index += 1;
            }
        }
        *traversals = kept;
        Ok(())
    }

    fn is_trivially_different(
        &self,
        first: &[Ctg],
        second: &[Ctg],
        min_sv_size: i32,
    ) -> Result<bool, AssemblyError> {
        let (first1, last1) = (first[0], first[first.len() - 1]);
        let (first2, last2) = (second[0], second[second.len() - 1]);
        if first1 != first2 || last1 != last2 {
            return Ok(false);
        }
        let interior1 = self.sequence_length(first) - self.size(first1) - self.size(last1);
        let interior2 = self.sequence_length(second) - self.size(first2) - self.size(last2);
        if (interior1 - interior2).abs() >= min_sv_size {
            return Ok(false);
        }
        let max_interior = interior1.max(interior2);
        if max_interior < min_sv_size {
            return Ok(true);
        }
        let out_of_bounds = || AssemblyError {
            class: "java.lang.ArrayIndexOutOfBoundsException",
            message: "Index -1 out of bounds for length 0".to_string(),
        };
        let row_len = first.len() - 1;
        if row_len == 0 {
            return Err(out_of_bounds());
        }
        let mut rows = [vec![0i32; row_len], vec![0i32; row_len]];
        let mut pair = 0usize;
        let n_rows = second.len() - 1;
        let mut idx2 = 1;
        while idx2 != n_rows {
            if idx2 >= second.len() {
                return Err(out_of_bounds());
            }
            let (current, previous) = if pair == 0 {
                let (a, b) = rows.split_at_mut(1);
                (&mut a[0], &b[0])
            } else {
                let (a, b) = rows.split_at_mut(1);
                (&mut b[0], &a[0])
            };
            pair ^= 1;
            let id2 = self.id(second[idx2]);
            for idx1 in 1..row_len {
                let tig1 = first[idx1];
                if self.id(tig1) == id2 {
                    let extend = self.id(first[idx1 - 1]) == self.id(second[idx2 - 1]);
                    current[idx1] = previous[idx1 - 1]
                        + if extend {
                            self.n_kmers(tig1)
                        } else {
                            self.size(tig1)
                        };
                } else {
                    current[idx1] = current[idx1 - 1].max(previous[idx1]);
                }
            }
            idx2 += 1;
        }
        let common = rows[pair ^ 1][row_len - 1];
        Ok(max_interior - common < min_sv_size)
    }

    /// `createScaffolds`; `Err(None)` is `AssemblyTooComplexException`.
    fn create_scaffolds(
        &self,
        traversals: &mut Vec<Traversal>,
        too_many: i32,
        min_sv_size: i32,
    ) -> Result<Vec<Traversal>, Option<AssemblyError>> {
        self.remove_trivially_different(traversals, min_sv_size)
            .map_err(Some)?;
        let mut by_first: HashMap<Ctg, Vec<i64>> = HashMap::new();
        for (index, traversal) in traversals.iter().enumerate() {
            by_first.entry(traversal[0]).or_default().push(index as i64);
            let rc = traversal_rc(traversal);
            by_first.entry(rc[0]).or_default().push(!(index as i64));
        }
        let mut scaffolds = Vec::new();
        let mut touched = vec![false; traversals.len()];
        for index in 0..traversals.len() {
            if !touched[index] {
                let traversal = traversals[index].clone();
                touched[index] = true;
                let mut starting = HashSet::new();
                let mut down = Vec::new();
                self.walk_traversals(
                    traversal.clone(),
                    &mut touched,
                    &mut starting,
                    &by_first,
                    traversals,
                    &mut down,
                )
                .map_err(Some)?;
                let mut up = Vec::new();
                self.walk_traversals(
                    traversal_rc(&traversal),
                    &mut touched,
                    &mut starting,
                    &by_first,
                    traversals,
                    &mut up,
                )
                .map_err(Some)?;
                for d in &down {
                    for u in &up {
                        if scaffolds.len() as i32 >= too_many {
                            return Err(None);
                        }
                        scaffolds
                            .push(combine(&traversal_rc(u), d, traversal.len()).map_err(Some)?);
                    }
                }
            }
        }
        Ok(scaffolds)
    }

    fn walk_traversals(
        &self,
        traversal: Traversal,
        touched: &mut [bool],
        starting: &mut HashSet<Ctg>,
        by_first: &HashMap<Ctg, Vec<i64>>,
        traversals: &[Traversal],
        extensions: &mut Vec<Traversal>,
    ) -> Result<(), AssemblyError> {
        let first = traversal[0];
        let last = *traversal.last().unwrap();
        let list = if starting.contains(&first) || self.is_cycle_member(last) {
            None
        } else {
            by_first.get(&last)
        };
        let Some(list) = list else {
            extensions.push(traversal);
            return Ok(());
        };
        starting.insert(first);
        for &index in list {
            let extension = if index >= 0 {
                touched[index as usize] = true;
                traversals[index as usize].clone()
            } else {
                let rc_index = (!index) as usize;
                touched[rc_index] = true;
                traversal_rc(&traversals[rc_index])
            };
            let combined = combine(&traversal, &extension, 1)?;
            self.walk_traversals(
                combined, touched, starting, by_first, traversals, extensions,
            )?;
        }
        starting.remove(&first);
        Ok(())
    }

    // ----- writers -----

    fn write_gfa(&mut self, contigs: &[Ctg], traversals: &[Traversal]) -> String {
        for &contig in contigs {
            self.set_marked(contig, false);
        }
        let mut out = String::from("H\tVN:Z:2.0\n");
        for &contig in contigs {
            if !self.is_marked(contig) {
                self.write_contig(contig, &mut out);
            }
        }
        for traversal in traversals {
            out.push_str("O\t*\t");
            out.push_str(
                &traversal
                    .iter()
                    .map(|c| self.reference(*c))
                    .collect::<Vec<_>>()
                    .join(" "),
            );
            out.push('\n');
        }
        out
    }

    fn write_contig(&mut self, contig: Ctg, out: &mut String) {
        let canonical = contig.canonical();
        self.set_marked(canonical, true);
        let sequence = self.sequence(canonical);
        out.push_str(&format!(
            "S\t{}\t{}\t{}\tMO:i:{}\tFO:i:{}\tLO:i:{}\n",
            self.name(canonical),
            sequence.len(),
            String::from_utf8_lossy(&sequence),
            self.max_observations(canonical),
            self.observations(self.first_kmer(canonical)),
            self.observations(self.last_kmer(canonical)),
        ));
        for successor in self.list(Graph::successors(contig)) {
            if !self.is_marked(successor) {
                self.write_contig(successor, out);
            }
            let length = self.size(contig);
            out.push_str(&format!(
                "E\t*\t{}\t{}\t{}\t{}$\t0\t{}\t{}M\n",
                self.reference(contig),
                self.reference(successor),
                length - KSIZE as i32 + 1,
                length,
                KSIZE - 1,
                KSIZE - 1
            ));
        }
        for predecessor in self.list(Graph::predecessors(contig)) {
            if !self.is_marked(predecessor) {
                self.write_contig(predecessor, out);
            }
        }
    }

    fn write_traversals(&self, assembly_name: &str, traversals: &[Traversal]) -> String {
        let mut out = String::new();
        for (index, traversal) in traversals.iter().enumerate() {
            out.push_str(&format!(
                ">{assembly_name}_t{} {}\n",
                index + 1,
                self.traversal_name(traversal)
            ));
            out.push_str(&String::from_utf8_lossy(
                &self.traversal_sequence(traversal),
            ));
            out.push('\n');
        }
        out
    }
}

/// An insertion-ordered `Map<Contig, List<TransitPairCount>>`.
#[derive(Default)]
struct TransitMap {
    keys: Vec<Ctg>,
    values: HashMap<Ctg, Vec<Tpc>>,
}

impl TransitMap {
    fn get(&self, key: Ctg) -> Option<&Vec<Tpc>> {
        self.values.get(&key)
    }
    fn entry(&mut self, key: Ctg) -> &mut Vec<Tpc> {
        if !self.values.contains_key(&key) {
            self.keys.push(key);
        }
        self.values.entry(key).or_default()
    }
    fn iter(&self) -> impl Iterator<Item = (Ctg, &Vec<Tpc>)> {
        self.keys.iter().map(|k| (*k, &self.values[k]))
    }
}

fn add_sub_path(sub: Vec<Ctg>, longest: &mut Vec<Vec<Ctg>>) {
    for test in longest.iter_mut() {
        if is_prefix(&sub, test) {
            return;
        }
        if is_prefix(test, &sub) {
            *test = sub;
            return;
        }
    }
    longest.push(sub);
}

fn is_prefix(first: &[Ctg], second: &[Ctg]) -> bool {
    first.len() <= second.len() && first == &second[..first.len()]
}

fn traversal_rc(traversal: &[Ctg]) -> Traversal {
    traversal.iter().rev().map(|c| c.rc()).collect()
}

fn combine(first: &[Ctg], second: &[Ctg], overlap: usize) -> Result<Traversal, AssemblyError> {
    let len1 = first.len();
    if len1 < overlap || second.len() < overlap || first[len1 - overlap..] != second[..overlap] {
        return Err(AssemblyError::gatk("combining non-overlapping traversals"));
    }
    let mut combined = first.to_vec();
    combined.extend_from_slice(&second[overlap..]);
    Ok(combined)
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// `SequenceRC`, materialized: upper-cased, complemented, reversed, anything else an `N`.
fn reverse_complement_bases(sequence: &[u8]) -> Vec<u8> {
    sequence.iter().rev().map(|b| complement(*b)).collect()
}

/// `LocalAssembler.onTraversalSuccess`, after `PairWalker` has handed over every read.
///
/// `region_size` is the summed size of the traversal intervals, which `PairWalker` has padded
/// and merged, and which sizes the kmer set.
pub fn assemble(
    reads: &[AssemblyRead],
    region_size: i32,
    settings: &Settings,
) -> Result<Assembly, AssemblyError> {
    let mut graph = Graph {
        adjs: Vec::new(),
        contigs: Vec::new(),
        kmers: Hopscotch::new(region_size.wrapping_mul(10), kmer_entry_index),
        tpcs: Vec::new(),
    };
    for read in reads {
        graph.kmerize_read(&read.bases, &read.qualities, settings.q_min)?;
    }
    let mut contigs = graph.create_assembly(settings.min_thin_observations)?;
    if graph.fill_gaps(settings.min_gapfill_count, reads)? {
        contigs = graph.create_assembly(settings.min_thin_observations)?;
    }
    graph.mark_cycles(&contigs);
    let paths: Vec<Vec<PathPart>> = reads.iter().map(|read| graph.path(&read.bases)).collect();
    let map = graph.collect_transits(&paths);
    match graph.traverse_all_paths(&contigs, &paths, settings.too_many_traversals, &map) {
        Ok(mut traversals) => {
            contigs.sort_by_key(|c| graph.id(*c));
            let gfa = graph.write_gfa(&contigs, &traversals);
            if settings.no_scaffolding {
                let fasta = graph.write_traversals(&settings.assembly_name, &traversals);
                return Ok(Assembly {
                    gfa: Some(gfa),
                    fasta,
                });
            }
            let fasta = match graph.create_scaffolds(
                &mut traversals,
                settings.too_many_scaffolds,
                settings.min_sv_size,
            ) {
                Ok(scaffolds) => graph.write_traversals(&settings.assembly_name, &scaffolds),
                Err(None) => graph.write_traversals(&settings.assembly_name, &traversals),
                Err(Some(error)) => return Err(error),
            };
            Ok(Assembly {
                gfa: Some(gfa),
                fasta,
            })
        }
        Err(TooComplex) => {
            let singles: Vec<Traversal> = contigs.iter().map(|c| vec![*c]).collect();
            Ok(Assembly {
                gfa: None,
                fasta: graph.write_traversals(&settings.assembly_name, &singles),
            })
        }
    }
}

// ---------------------------------------------------------------------------------------------
// PairWalker
// ---------------------------------------------------------------------------------------------

/// An interval as `PairWalker` holds it: contig, one-based start and end.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Span {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

/// `PairWalker.transformTraversalIntervals`: each interval padded within its contig, then
/// merged into the one before it whenever both are on the same contig.
///
/// The merge test compares the previous START with the current END twice, so it holds for every
/// pair of intervals on one contig, and `mergeWithContiguous` then refuses the two that do not
/// touch: `-L` intervals further apart than twice the padding are a `GATKException`.
pub fn pad_intervals(
    intervals: &[Span],
    padding: i32,
    contig_length: impl Fn(&str) -> Option<i32>,
) -> Result<Vec<Span>, AssemblyError> {
    let mut padded = Vec::new();
    let mut previous: Option<Span> = None;
    for interval in intervals {
        let length = contig_length(&interval.contig).ok_or_else(|| AssemblyError {
            class: "java.lang.IllegalArgumentException",
            message: format!(
                "Contig {} not found in provided dictionary",
                interval.contig
            ),
        })?;
        if padding < 0 {
            return Err(AssemblyError {
                class: "java.lang.IllegalArgumentException",
                message: "padding must be >= 0".to_string(),
            });
        }
        let current = Span {
            contig: interval.contig.clone(),
            start: (interval.start - padding).max(1),
            end: (interval.end + padding).min(length),
        };
        match previous.take() {
            None => previous = Some(current),
            Some(prev) if prev.contig == current.contig && prev.start <= current.end + 1 => {
                let contiguous = prev.start <= current.end + 1 && current.start <= prev.end + 1;
                if !contiguous {
                    return Err(AssemblyError::gatk(format!(
                        "The two intervals need to be contiguous: {}:{}-{} {}:{}-{}",
                        prev.contig,
                        prev.start,
                        prev.end,
                        current.contig,
                        current.start,
                        current.end
                    )));
                }
                previous = Some(Span {
                    contig: prev.contig,
                    start: prev.start.min(current.start),
                    end: prev.end.max(current.end),
                });
            }
            Some(prev) => {
                padded.push(prev);
                previous = Some(current);
            }
        }
    }
    if let Some(prev) = previous {
        padded.push(prev);
    }
    Ok(padded)
}

/// `PairWalker.RegionChecker`.
struct RegionChecker<'a> {
    intervals: &'a [Span],
    next: usize,
    current: Option<&'a Span>,
}

impl<'a> RegionChecker<'a> {
    fn new(intervals: &'a [Span]) -> Self {
        RegionChecker {
            intervals,
            next: 1,
            current: intervals.first(),
        }
    }

    fn advance(&mut self) -> Option<&'a Span> {
        let result = self.intervals.get(self.next);
        self.next += 1;
        result
    }

    fn is_in_interval(
        &mut self,
        read: &BamRecord,
        contig_of: &dyn Fn(i32) -> Option<String>,
        index_of: &dyn Fn(&str) -> i32,
    ) -> bool {
        let Some(mut current) = self.current else {
            return false;
        };
        if read::is_unmapped(read) {
            return false;
        }
        let read_contig = contig_of(read.reference_index).unwrap_or_default();
        if current.contig != read_contig {
            let read_id = index_of(&read_contig);
            let mut current_id;
            loop {
                current_id = index_of(&current.contig);
                if current_id >= read_id {
                    break;
                }
                match self.advance() {
                    None => {
                        self.current = None;
                        return false;
                    }
                    Some(next) => {
                        current = next;
                        self.current = Some(next);
                    }
                }
            }
            if current_id > read_id {
                return false;
            }
        }
        let start = read.alignment_start;
        while current.end < start {
            match self.advance() {
                None => {
                    self.current = None;
                    return false;
                }
                Some(next) => {
                    current = next;
                    self.current = Some(next);
                    if current.contig != read_contig {
                        return false;
                    }
                }
            }
        }
        current.contig == read_contig
            && current.start <= read.alignment_end()
            && start <= current.end
    }
}

struct PairBufferEntry {
    read: BamRecord,
    in_interval: bool,
    distant: bool,
}

/// A slot of the pair buffer: the read name's hash code and where its entry is kept. The set is
/// sized for a million reads, as `PairWalker`'s is, so it holds a handle and not the read.
#[derive(Clone)]
struct PairSlot {
    hash: i32,
    entry: u32,
}

fn pair_slot_index(slot: &PairSlot, capacity: i32) -> i32 {
    default_index(slot.hash, capacity)
}

fn assembly_read(read: &BamRecord) -> AssemblyRead {
    AssemblyRead {
        bases: read.read_bases.clone(),
        qualities: read.base_qualities.clone(),
    }
}

/// `PairWalker`'s traversal over the reads that reached `apply`, and `LocalAssembler.apply`'s
/// `trimOverruns`: the reads in the order `LocalAssembler` stores them.
pub fn pair_reads(
    reads: Vec<BamRecord>,
    intervals: &[Span],
    contig_of: &dyn Fn(i32) -> Option<String>,
    index_of: &dyn Fn(&str) -> i32,
) -> Vec<AssemblyRead> {
    let mut checker = RegionChecker::new(intervals);
    let mut buffer: Hopscotch<PairSlot> = Hopscotch::new(1_000_000, pair_slot_index);
    let mut entries: Vec<Option<PairBufferEntry>> = Vec::new();
    let mut distant_processed: HashSet<String> = HashSet::new();
    let mut out = Vec::new();
    for read in reads {
        if !read::is_paired(&read)
            || read::is_secondary_alignment(&read)
            || read::is_supplementary_alignment(&read)
        {
            out.push(assembly_read(&read));
            continue;
        }
        let in_interval = checker.is_in_interval(&read, contig_of, index_of);
        let distant = read.tags.get(Tag::new(b"DM")).is_some();
        let read = if distant {
            undo_distant_mate(read, index_of)
        } else {
            read
        };
        let name = read.read_name.clone();
        let hash = java_string_hash(&name);
        let entry = PairBufferEntry {
            read,
            in_interval,
            distant,
        };
        let same_name = |slot: &PairSlot| {
            entries[slot.entry as usize]
                .as_ref()
                .is_some_and(|held| held.read.read_name == name)
        };
        let mate = buffer
            .find(default_index(hash, buffer.capacity()), same_name)
            .map(|slot| slot.entry);
        match mate {
            None => {
                let slot = PairSlot {
                    hash,
                    entry: entries.len() as u32,
                };
                entries.push(Some(entry));
                let names = |a: &PairSlot, b: &PairSlot| {
                    entries[a.entry as usize]
                        .as_ref()
                        .map(|e| &e.read.read_name)
                        == entries[b.entry as usize]
                            .as_ref()
                            .map(|e| &e.read.read_name)
                };
                buffer.add(slot, names);
            }
            Some(mate_slot) => {
                let mate = entries[mate_slot as usize].take().expect("a buffered mate");
                if entry.in_interval || mate.in_interval {
                    if !entry.distant && !mate.distant {
                        push_pair(mate.read, entry.read, &mut out);
                    } else if !distant_processed.remove(&name) {
                        distant_processed.insert(name.clone());
                        push_pair(mate.read, entry.read, &mut out);
                    }
                }
                let capacity = buffer.capacity();
                buffer.remove(default_index(hash, capacity), |slot| {
                    slot.entry == mate_slot
                });
            }
        }
    }
    for slot in buffer.entries() {
        if let Some(entry) = &entries[slot.entry as usize] {
            if entry.in_interval {
                out.push(assembly_read(&entry.read));
            }
        }
    }
    out
}

fn push_pair(mut read: BamRecord, mut mate: BamRecord, out: &mut Vec<AssemblyRead>) {
    trim_overruns(&mut read, &mut mate);
    out.push(assembly_read(&read));
    out.push(assembly_read(&mate));
}

/// `PrintDistantMates.undoDistantMateAlterations`: the alignment the `OA` tag remembers.
fn undo_distant_mate(mut read: BamRecord, index_of: &dyn Fn(&str) -> i32) -> BamRecord {
    let Some(htsjdk_bam::tag::TagValue::Str(oa)) = read.tags.get(Tag::new(b"OA")).cloned() else {
        return read;
    };
    read.tags.remove(Tag::new(b"DM"));
    read.tags.remove(Tag::new(b"OA"));
    let tokens: Vec<&str> = oa.split(',').collect();
    if tokens.len() >= 5 {
        read.reference_index = index_of(tokens[0]);
        read.alignment_start = tokens[1].parse().unwrap_or(read.alignment_start);
        if tokens[2] == "-" {
            read.flags |= read::flags::READ_REVERSE_STRAND;
        } else {
            read.flags &= !read::flags::READ_REVERSE_STRAND;
        }
        if let Some(cigar) = parse_cigar(tokens[3]) {
            read.cigar = cigar;
        }
        read.mapping_quality = tokens[4].parse().unwrap_or(read.mapping_quality);
        read.flags &= !read::flags::READ_UNMAPPED;
    }
    read
}

/// A CIGAR string as `TextCigarCodec.decode` reads it, `None` where it would refuse.
fn parse_cigar(text: &str) -> Option<htsjdk_bam::cigar::Cigar> {
    let mut elements = Vec::new();
    let mut length = 0u32;
    for character in text.bytes() {
        if character.is_ascii_digit() {
            length = length * 10 + (character - b'0') as u32;
            continue;
        }
        let op = match character {
            b'M' => Op::M,
            b'I' => Op::I,
            b'D' => Op::D,
            b'N' => Op::N,
            b'S' => Op::S,
            b'H' => Op::H,
            b'P' => Op::P,
            b'=' => Op::Eq,
            b'X' => Op::X,
            _ => return None,
        };
        elements.push(CigarElement { length, op });
        length = 0;
    }
    Some(htsjdk_bam::cigar::Cigar { elements })
}

/// `LocalAssembler.trimOverruns`.
fn trim_overruns(read: &mut BamRecord, mate: &mut BamRecord) {
    if read::is_unmapped(read)
        || read::is_unmapped(mate)
        || read::is_reverse_strand(read) == read::is_reverse_strand(mate)
    {
        return;
    }
    if (read.alignment_start - read.mate_alignment_start).abs() > 1 {
        return;
    }
    let read_length = read.cigar.reference_length() as i64;
    let mate_length = mate.cigar.reference_length() as i64;
    if (read_length - mate_length).abs() > 1 {
        return;
    }
    if read::is_reverse_strand(mate) {
        trim_clips(read, mate);
    } else {
        trim_clips(mate, read);
    }
}

fn trim_clips(fwd: &mut BamRecord, rev: &mut BamRecord) {
    let (Some(fwd_first), Some(fwd_last)) = (fwd.cigar.elements.first(), fwd.cigar.elements.last())
    else {
        return;
    };
    let (Some(rev_first), Some(rev_last)) = (rev.cigar.elements.first(), rev.cigar.elements.last())
    else {
        return;
    };
    if !(fwd_first.op == Op::M
        && fwd_last.op == Op::S
        && rev_first.op == Op::S
        && rev_last.op == Op::M)
    {
        return;
    }
    let last_length = fwd_last.length as usize;
    let first_length = rev_first.length as usize;
    let keep = fwd.read_bases.len() - last_length;
    fwd.read_bases.truncate(keep);
    if !fwd.base_qualities.is_empty() {
        let keep = fwd.base_qualities.len() - last_length;
        fwd.base_qualities.truncate(keep);
    }
    let last_index = fwd.cigar.elements.len() - 1;
    fwd.cigar.elements[last_index] = CigarElement {
        length: last_length as u32,
        op: Op::H,
    };
    rev.read_bases.drain(..first_length);
    if !rev.base_qualities.is_empty() {
        rev.base_qualities.drain(..first_length);
    }
    rev.cigar.elements[0] = CigarElement {
        length: first_length as u32,
        op: Op::H,
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_matches_the_table() {
        // AAAA...A (31) is all zero; its reverse complement is all T, 3 in every position.
        assert_eq!(reverse_complement_kval(0), KMASK);
        let kval = 0b01_10_11; // ...ACGT tail: C G T at the low end
        let rc = reverse_complement_kval(kval);
        assert_eq!(reverse_complement_kval(rc), kval);
    }

    #[test]
    fn capacity_is_the_first_legal_size_above_the_load() {
        assert_eq!(compute_capacity(1_000_000), 1_482_907);
        assert_eq!(compute_capacity(0), 251);
    }

    #[test]
    fn java_string_hash_is_javas() {
        assert_eq!(java_string_hash("abc"), 96354);
        assert_eq!(java_string_hash(""), 0);
    }
}
