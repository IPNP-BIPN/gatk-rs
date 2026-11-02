//! `FindBreakpointEvidenceSpark.gatherEvidenceAndWriteContigSamFile` (GATK 4.6.2.0): the evidence,
//! its clusters, the template names and k-mers of each interval, and each interval's assembly.
//!
//! # One partition, and what Spark does to it
//!
//! Everything runs as the reference runs it under `--spark-master local[1]` on an input smaller than
//! a split. `PartitionCrossingChecker` is then never on a boundary, since `PartitionBounds.UNMAPPED`
//! is `Integer.MAX_VALUE`. The `reduceByKey`s still shuffle, and the order that leaves a shuffle is
//! what `KmerCleaner` fills a hopscotch table with, and so what `--kmer-intervals` prints: see
//! [`reduce_by_key`].
//!
//! # Hopscotch all the way down
//!
//! The template names and the k-mers live in hopscotch tables through every stage, with entries
//! removed in between, and each table's layout decides what the next stage finds first: the order of
//! a name's intervals decides the order a read reaches them in, which is the order of the FASTQ. So
//! the tables here are the reference's, and nothing is sorted that it does not sort.
//!
//! # A broadcast is the same object
//!
//! In a local master `TorrentBroadcast` hands back the object it was given, so a table is not
//! rebuilt by Kryo on its way to a task, and the cached `evidenceRDD` holds the objects the density
//! filter sets `validated` on. A `collect` does go through Kryo, and only the target links have
//! anything in them that shows it: see [`NameSet`].
//!
//! # What this port refuses
//!
//! Fermi-lite aborts the JVM when no read of an interval survives its filters, and corrupts its heap
//! on a read holding a Phred 0, whose quality string ends at the NUL. Neither gives an answer to be
//! the same as, so each stops at a refusal of the port's own.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::Instant;

use gatk_engine::java_hash::{compare_strings, string_hash_code, JavaHashMap};
use gatk_engine::read_utils;
use gatk_fermilite::{FmlError, FmlOpt, MAG_F_AGGRESSIVE, MAG_F_POPOPEN};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op as CigarOp};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Value, VariantContext};

use crate::bwa_aligner::{ContigAligner, ContigAlignment};
use crate::hopscotch::{Hopscotch, HopscotchEntry, HopscotchFailed};
use crate::read_metadata::ReadMetadata;
use crate::sv_discovery::{append_only_group, Reference};
use crate::sv_evidence::{
    self, contained_in_region_to_ignore, fragment_ordinal, is_mapped, is_primary_line, Evidence,
    EvidenceArguments, EvidenceKind, ReadKind, SvInterval,
};
use crate::sv_kmer::{DustFilteredKmerizer, Kmerizer, SvKmerLong, SvKmerShort};

const DUPLICATE: u16 = 0x400;
const QC_FAIL: u16 = 0x200;
const REVERSE: u16 = 0x10;

/// The largest table this port will allocate up front. `KmerCleaner` sizes its table from twice the
/// reference, which is the reference's choice and a dozen gigabytes on a human genome.
const MAX_TABLE: i64 = 1 << 27;

/// How a stage stops.
#[derive(Debug)]
pub enum Failure {
    /// Thrown on the driver and not a `UserException`.
    Thrown { class: &'static str, message: String },
    /// Thrown in a task, as `class: message`, or `class` alone where there is no message.
    Task(String),
    /// `UserException` and its subclasses.
    User(String),
    /// `BucketUtils.createFile`.
    Create { path: String, error: std::io::Error },
    /// `BucketUtils.openFile`.
    Read { path: String, error: std::io::Error },
    /// The port's own, which the runner reports as one.
    Limitation(String),
    /// What the port itself could not do with no reference behaviour behind it.
    Port(String),
}

const GATK_EXCEPTION: &str = "org.broadinstitute.hellbender.exceptions.GATKException";
const ILLEGAL_STATE: &str = "java.lang.IllegalStateException";
const ILLEGAL_ARGUMENT: &str = "java.lang.IllegalArgumentException";

impl From<HopscotchFailed> for Failure {
    fn from(failed: HopscotchFailed) -> Failure {
        Failure::Task(format!("{ILLEGAL_STATE}: {}", failed.0))
    }
}

fn hash_order_limit(error: gatk_engine::java_hash::HashOrderError) -> Failure {
    Failure::Limitation(format!(
        "a Java HashSet grew past what this port has measured ({error:?})."
    ))
}

/// `FindBreakpointEvidenceSparkArgumentCollection`.
#[derive(Debug, Clone)]
pub struct Arguments {
    pub evidence: EvidenceArguments,
    pub k_size: i32,
    pub max_dust_score: i32,
    pub min_evidence_weight_per_coverage: f64,
    pub min_coherent_evidence_weight_per_coverage: f64,
    pub min_kmers_per_interval: i32,
    pub cleaner_max_intervals: i32,
    pub cleaner_min_kmer_count: i32,
    pub cleaner_max_copy_number: i32,
    pub assembly_to_mapped_size_ratio_guess: i32,
    pub max_fastq_size: i32,
    pub exclusion_interval_padding: i32,
    pub include_mapping_location: bool,
    pub interval_only_assembly: bool,
    pub adapter_sequence: Option<String>,
    pub write_gfas: bool,
    pub pop_variant_bubbles: bool,
    pub remove_shadowed_contigs: bool,
    pub expand_assembly_graph: bool,
    pub z_dropoff: i32,
    pub aligner_index_image: Option<String>,
    pub evidence_dir: Option<String>,
    pub unfiltered_evidence_dir: Option<String>,
    pub interval_file: Option<String>,
    pub qnames_mapped_file: Option<String>,
    pub kmer_file: Option<String>,
    pub qnames_assembly_file: Option<String>,
    pub fastq_dir: Option<String>,
    pub target_link_file: Option<String>,
    pub exclusion_intervals_file: Option<String>,
    pub kmers_to_ignore_file: Option<String>,
}

fn write_file(path: &str, bytes: &[u8]) -> Result<(), Failure> {
    std::fs::write(path, bytes).map_err(|error| Failure::Create {
        path: path.to_string(),
        error,
    })
}

fn read_file(path: &str) -> Result<Vec<u8>, Failure> {
    std::fs::read(path).map_err(|error| Failure::Read {
        path: path.to_string(),
        error,
    })
}

/// `JavaRDD.saveAsTextFile` of one partition, which an existing directory refuses before anything
/// runs.
fn save_as_text_file(directory: &str, lines: impl Iterator<Item = String>) -> Result<(), Failure> {
    let path = std::path::Path::new(directory);
    if path.exists() {
        let absolute = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
        return Err(Failure::Thrown {
            class: "org.apache.hadoop.mapred.FileAlreadyExistsException",
            message: format!("Output directory file:{} already exists", absolute.display()),
        });
    }
    let mut text = String::new();
    for line in lines {
        text.push_str(&line);
        text.push('\n');
    }
    std::fs::create_dir_all(path)
        .and_then(|()| std::fs::write(path.join("part-00000"), &text))
        .and_then(|()| std::fs::write(path.join("_SUCCESS"), b""))
        .map_err(|error| Failure::Port(format!("{directory}: {error}")))
}

/// `MurmurHash3.productHash` of Scala 2.13's `Tuple2(0, key)`, the partition being 0: the seed, the
/// name of the product, the two elements, the arity.
fn partitioned_key_hash(key_hash: i32) -> i32 {
    fn mix(hash: u32, data: u32) -> u32 {
        let data = data
            .wrapping_mul(0xcc9e_2d51)
            .rotate_left(15)
            .wrapping_mul(0x1b87_3593);
        (hash ^ data)
            .rotate_left(13)
            .wrapping_mul(5)
            .wrapping_add(0xe654_6b64)
    }
    // `MurmurHash3.productSeed` and `"Tuple2".hashCode`.
    let mut hash = (-889275714i32) as u32;
    hash = mix(hash, (-1778623094i32) as u32);
    hash = mix(hash, 0);
    hash = mix(hash, key_hash as u32);
    let mut hash = hash ^ 2;
    hash ^= hash >> 16;
    hash = hash.wrapping_mul(0x85eb_ca6b);
    hash ^= hash >> 13;
    hash = hash.wrapping_mul(0xc2b2_ae35);
    (hash ^ (hash >> 16)) as i32
}

/// `reduceByKey` of one map partition into one reduce partition. On the map side `ExternalSorter`
/// combines in a `PartitionedAppendOnlyMap` keyed by `(partition, key)`, walks it in slot order, and
/// sorts stably by `key.hashCode` as a signed int; the reducer combines in an `AppendOnlyMap` keyed
/// by the key alone and walks that in slot order.
fn reduce_by_key<K: PartialEq, V>(
    pairs: Vec<(K, V)>,
    hash_of: impl Fn(&K) -> i32,
    reduce: impl Fn(V, V) -> V,
) -> Vec<(K, V)> {
    let combine = |values: Vec<V>| {
        values
            .into_iter()
            .reduce(&reduce)
            .expect("a group holds a value")
    };
    let mut mapped: Vec<(K, V)> = append_only_group(pairs, |key| partitioned_key_hash(hash_of(key)))
        .into_iter()
        .map(|(key, values)| (key, combine(values)))
        .collect();
    mapped.sort_by_key(|(key, _)| hash_of(key));
    append_only_group(mapped, &hash_of)
        .into_iter()
        .map(|(key, values)| (key, combine(values)))
        .collect()
}

/// `SVIntervalTree`: keys in `SVInterval` order, and `overlappers` found by a range wide enough for
/// the longest interval ever put.
struct IntervalMap<V> {
    map: BTreeMap<SvInterval, V>,
    longest: i32,
}

impl<V> IntervalMap<V> {
    fn new() -> IntervalMap<V> {
        IntervalMap {
            map: BTreeMap::new(),
            longest: 0,
        }
    }

    fn put(&mut self, interval: SvInterval, value: V) {
        self.longest = self.longest.max(interval.end.wrapping_sub(interval.start));
        self.map.insert(interval, value);
    }

    fn overlappers(&self, query: &SvInterval) -> Vec<SvInterval> {
        let from = SvInterval {
            contig: query.contig,
            start: query.start.saturating_sub(self.longest),
            end: i32::MIN,
        };
        self.map
            .range(from..)
            .take_while(|(key, _)| key.contig == query.contig && key.start < query.end)
            .filter(|(key, _)| key.overlaps(query))
            .map(|(key, _)| *key)
            .collect()
    }
}

type Stranded = (SvInterval, bool);

/// `PairedStrandedIntervalTree`. Its overlap iterator moves to the next left end on every `hasNext`,
/// so a query finds at most one right end per left end, and a left end keeps the strand it first came
/// with whatever is added under it.
struct PairedTree<V> {
    left: IntervalMap<(bool, IntervalMap<(bool, V)>)>,
}

impl<V> PairedTree<V> {
    fn new() -> PairedTree<V> {
        PairedTree {
            left: IntervalMap::new(),
        }
    }

    fn contains(&self, left: &Stranded, right: &Stranded) -> bool {
        match self.left.map.get(&left.0) {
            Some((strand, rights)) if *strand == left.1 => rights
                .map
                .get(&right.0)
                .is_some_and(|(strand, _)| *strand == right.1),
            _ => false,
        }
    }

    fn put(&mut self, left: Stranded, right: Stranded, value: V) -> bool {
        if self.contains(&left, &right) {
            return false;
        }
        if let Some((_, rights)) = self.left.map.get_mut(&left.0) {
            rights.put(right.0, (right.1, value));
        } else {
            let mut rights = IntervalMap::new();
            rights.put(right.0, (right.1, value));
            self.left.put(left.0, (left.1, rights));
        }
        true
    }

    /// The pairs `overlappers` walks to, as (left, right) keys.
    fn overlappers(&self, left: &Stranded, right: &Stranded) -> Vec<(SvInterval, SvInterval)> {
        let mut found = Vec::new();
        for key in self.left.overlappers(&left.0) {
            let (strand, rights) = &self.left.map[&key];
            if *strand != left.1 {
                continue;
            }
            for right_key in rights.overlappers(&right.0) {
                if rights.map[&right_key].0 == right.1 {
                    found.push((key, right_key));
                    break;
                }
            }
        }
        found
    }

    /// The iterator's `remove`: the left end goes with its last right end.
    fn remove(&mut self, key: &(SvInterval, SvInterval)) -> Option<(bool, V)> {
        let (_, rights) = self.left.map.get_mut(&key.0)?;
        let removed = rights.map.remove(&key.1);
        if rights.map.is_empty() {
            self.left.map.remove(&key.0);
        }
        removed
    }

    /// `iterator`: left ends in order, each one's right ends in order, with the strands as held.
    fn entries(&self) -> Vec<(Stranded, Stranded, &V)> {
        let mut entries = Vec::new();
        for (left, (left_strand, rights)) in &self.left.map {
            for (right, (right_strand, value)) in &rights.map {
                entries.push(((*left, *left_strand), (*right, *right_strand), value));
            }
        }
        entries
    }
}

/// An `EvidenceTargetLink`'s template names: `Collections.singleton` until two links merge, after
/// which a `HashSet` copied with each merge. `Utils.join` walks them, so the order is the file's.
#[derive(Clone)]
enum NameSet {
    Single(String),
    Hashed(JavaHashMap<String, ()>),
}

impl NameSet {
    fn empty() -> NameSet {
        NameSet::Hashed(JavaHashMap::new())
    }

    fn contains(&self, name: &str) -> bool {
        match self {
            NameSet::Single(held) => held == name,
            NameSet::Hashed(set) => set.contains_key(&name.to_string()),
        }
    }

    fn names(&self) -> Vec<&String> {
        match self {
            NameSet::Single(held) => vec![held],
            NameSet::Hashed(set) => set.keys().collect(),
        }
    }

    /// `new HashSet<>(this)`.
    fn copied(&self) -> JavaHashMap<String, ()> {
        let names = self.names();
        JavaHashMap::copy_of(names.into_iter())
    }

    /// What `collect` leaves after Kryo's `CollectionSerializer`: a `new HashSet()` filled in this
    /// one's order, and a singleton still a singleton.
    fn collected(&self) -> Result<NameSet, Failure> {
        Ok(match self {
            NameSet::Single(held) => NameSet::Single(held.clone()),
            NameSet::Hashed(set) => {
                set.check().map_err(hash_order_limit)?;
                let mut fresh = JavaHashMap::new();
                for name in set.keys() {
                    fresh.insert(name.clone(), ());
                }
                fresh.check().map_err(hash_order_limit)?;
                NameSet::Hashed(fresh)
            }
        })
    }
}

/// `EvidenceTargetLink`.
#[derive(Clone)]
pub struct TargetLink {
    pub source: (SvInterval, bool),
    pub target: (SvInterval, bool),
    pub split_reads: i32,
    pub read_pairs: i32,
    read_pair_names: NameSet,
    split_read_names: NameSet,
}

impl TargetLink {
    fn new(
        source: Stranded,
        target: Stranded,
        split_reads: i32,
        read_pairs: i32,
        read_pair_names: NameSet,
        split_read_names: NameSet,
    ) -> TargetLink {
        let (source, target) = if source.0.is_upstream_of(&target.0) {
            (source, target)
        } else {
            (target, source)
        };
        TargetLink {
            source,
            target,
            split_reads,
            read_pairs,
            read_pair_names,
            split_read_names,
        }
    }

    /// `toBedpeString`.
    fn bedpe(&self, metadata: &ReadMetadata) -> String {
        let (source, target) = (self.source.0, self.target.0);
        let name = |contig: i32| &metadata.contigs[contig as usize];
        let join = |names: &NameSet| {
            names
                .names()
                .into_iter()
                .map(String::as_str)
                .collect::<Vec<_>>()
                .join(",")
        };
        let id = format!(
            "{}_{}_{}_{}_{}_{}_{}{}_{}_{}",
            name(source.contig),
            source.start - 1,
            source.end,
            name(target.contig),
            target.start - 1,
            target.end,
            if self.source.1 { "P" } else { "M" },
            if self.target.1 { "P" } else { "M" },
            self.split_reads,
            self.read_pairs
        );
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{id}\t{}\t{}\t{}\tSR:{}\tRP:{}",
            name(source.contig),
            source.start - 1,
            source.end,
            name(target.contig),
            target.start - 1,
            target.end,
            self.read_pairs.wrapping_add(self.split_reads),
            if self.source.1 { "+" } else { "-" },
            if self.target.1 { "+" } else { "-" },
            join(&self.split_read_names),
            join(&self.read_pair_names)
        )
    }
}

/// `SVInterval.intersect`, of what an overlap search found.
fn intersect(one: &SvInterval, two: &SvInterval) -> SvInterval {
    SvInterval {
        contig: one.contig,
        start: one.start.max(two.start),
        end: one.end.min(two.end),
    }
}

/// `EvidenceTargetLinkClusterer.cluster`.
fn cluster_links(evidence: &[Evidence]) -> Vec<TargetLink> {
    let mut current: PairedTree<TargetLink> = PairedTree::new();
    for one in evidence {
        let (Some(distal), EvidenceKind::Read(read)) = (&one.distal, &one.kind) else {
            continue;
        };
        let name = &read.template_name;
        let is_split = matches!(read.kind, ReadKind::SplitRead { .. });
        let next = TargetLink::new(
            (one.location, read.forward),
            *distal,
            i32::from(is_split),
            i32::from(!is_split),
            if is_split {
                NameSet::empty()
            } else {
                NameSet::Single(name.clone())
            },
            if is_split {
                NameSet::Single(name.clone())
            } else {
                NameSet::empty()
            },
        );
        let mut updated = None;
        if let Some(key) = current.overlappers(&next.source, &next.target).first() {
            let (old_rights_strand, old) = {
                let rights = &current.left.map[&key.0].1;
                let (strand, link) = &rights.map[&key.1];
                (*strand, link.clone())
            };
            let _ = old_rights_strand;
            let seen = if is_split {
                old.split_read_names.contains(name)
            } else {
                old.read_pair_names.contains(name)
            };
            if seen {
                continue;
            }
            current.remove(key);
            let mut read_pair_names = old.read_pair_names.copied();
            let mut split_read_names = old.split_read_names.copied();
            if is_split {
                split_read_names.insert(name.clone(), ());
            } else {
                read_pair_names.insert(name.clone(), ());
            }
            updated = Some(TargetLink::new(
                (intersect(&old.source.0, &next.source.0), old.source.1),
                (intersect(&old.target.0, &next.target.0), old.target.1),
                old.split_reads + i32::from(is_split),
                old.read_pairs + i32::from(!is_split),
                NameSet::Hashed(read_pair_names),
                NameSet::Hashed(split_read_names),
            ));
        }
        let link = updated.unwrap_or(next);
        current.put(link.source, link.target, link.clone());
    }
    current
        .entries()
        .into_iter()
        .map(|(_, _, link)| link.clone())
        .collect()
}

/// `BreakpointDensityFilter` run out, on no partition boundary: what passes, in tree order, with
/// `validated` set on it.
fn density_filter(evidence: Vec<Evidence>, min_weight: f64, min_coherent: f64) -> Vec<Evidence> {
    let mut tree: IntervalMap<Vec<Evidence>> = IntervalMap::new();
    for one in evidence {
        match tree.map.get_mut(&one.location) {
            Some(list) => list.push(one),
            None => tree.put(one.location, vec![one]),
        }
    }
    let mut passed = Vec::new();
    let keys: Vec<SvInterval> = tree.map.keys().copied().collect();
    let has_enough = |tree: &IntervalMap<Vec<Evidence>>, interval: &SvInterval| -> bool {
        let mut targets: PairedTree<()> = PairedTree::new();
        let mut weight = 0i32;
        for overlapper in tree.overlappers(interval) {
            let list = &tree.map[&overlapper];
            weight = weight.wrapping_add(list.iter().fold(0, |sum, one| sum.wrapping_add(one.weight)));
            if f64::from(weight) >= min_weight {
                return true;
            }
            for one in list {
                if let (Some(distal), EvidenceKind::Read(read)) = (&one.distal, &one.kind) {
                    targets.put((one.location, read.forward), *distal, ());
                }
            }
        }
        targets.entries().iter().any(|(left, right, _)| {
            targets.overlappers(left, right).len() as f64 >= min_coherent
        })
    };
    for key in keys {
        let pass = tree.map[&key].iter().any(|one| one.validated) || has_enough(&tree, &key);
        if pass {
            let list = tree.map.get_mut(&key).expect("a key just listed");
            for one in list.iter_mut() {
                one.validated = true;
            }
            passed.extend(list.iter().cloned());
        }
    }
    passed
}

/// `BreakpointEvidenceClusterer` over a stream and the sentinel past the last contig.
fn cluster_evidence(evidence: &[Evidence], gap_size: i32, contigs: i32) -> Vec<Evidence> {
    let sentinel = (
        SvInterval {
            contig: contigs,
            start: 1,
            end: 1,
        },
        0,
    );
    let mut clusters = Vec::new();
    let mut current: Option<(SvInterval, i32)> = None;
    for (interval, weight) in evidence
        .iter()
        .map(|one| (one.location, one.weight))
        .chain(std::iter::once(sentinel))
    {
        current = Some(match current {
            None => (interval, weight),
            Some((held, sum)) if held.gap_len(&interval) < gap_size => {
                (held.join(&interval), sum.wrapping_add(weight))
            }
            Some((held, sum)) => {
                clusters.push(Evidence {
                    location: held,
                    weight: sum,
                    validated: true,
                    kind: EvidenceKind::Cluster,
                    distal: None,
                });
                (interval, weight)
            }
        });
    }
    clusters
}

/// `QNameAndInterval`, keyed by the name, colliding on both.
#[derive(Clone)]
struct QNameAndInterval {
    name: String,
    interval: i32,
}

impl HopscotchEntry for QNameAndInterval {
    fn key_hash(&self) -> i32 {
        string_hash_code(&self.name)
    }
    fn collides(&self, other: &Self) -> bool {
        self.interval == other.interval && self.name == other.name
    }
}

/// `KmerAndInterval`, keyed by the k-mer, colliding on both.
#[derive(Clone, PartialEq)]
struct KmerAndInterval {
    kmer: SvKmerLong,
    interval: i32,
}

impl HopscotchEntry for KmerAndInterval {
    fn key_hash(&self) -> i32 {
        self.kmer.hash_code()
    }
    fn collides(&self, other: &Self) -> bool {
        self == other
    }
}

/// `HopscotchSet<SVKmer>`.
#[derive(Clone)]
struct Kmer(SvKmerLong);

impl HopscotchEntry for Kmer {
    fn key_hash(&self) -> i32 {
        self.0.hash_code()
    }
    fn collides(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

/// `KmerAndCount` in a `HopscotchMap`.
#[derive(Clone)]
struct KmerAndCount {
    kmer: SvKmerLong,
    count: i32,
}

impl HopscotchEntry for KmerAndCount {
    fn key_hash(&self) -> i32 {
        self.kmer.hash_code()
    }
    fn collides(&self, other: &Self) -> bool {
        self.kmer == other.kmer
    }
}

fn not_junk(record: &BamRecord) -> bool {
    record.flags & (DUPLICATE | QC_FAIL) == 0
}

/// `Math.round(float)`, which is `floor(x + 1/2)` without the double rounding of the sum, and which
/// `f64` cannot get wrong since a float and a half both fit it exactly.
fn java_round(value: f32) -> i32 {
    (f64::from(value) + 0.5).floor() as i32
}

/// `ReadMetadata.getAccurateKmerCoverage`. `Math.exp` and `Math.pow` are intrinsics that can stand
/// an ulp off `StrictMath`, which the narrowing to float hides: measured, the two forms agree at
/// every one of 34310 points across K and mean qualities.
fn accurate_kmer_coverage(metadata: &ReadMetadata, k: i32) -> f32 {
    let quality = f64::from(metadata.mean_base_quality);
    let error = jmath::strict_math::pow(10.0, quality / -10.0);
    let probability = jmath::strict_math::exp(f64::from(k) * (-error).ln_1p()) as f32;
    metadata.coverage * probability
}

fn canonical(kmer: SvKmerLong, k: i32) -> Result<SvKmerLong, Failure> {
    if k & 1 == 0 {
        return Err(Failure::Task(format!(
            "{ILLEGAL_ARGUMENT}: Kmer length must be odd to canonicalize."
        )));
    }
    Ok(kmer.canonical(k))
}

/// `SVFileUtils.readKmersFile`, then the adapter's k-mers.
fn kmer_kill_set(arguments: &Arguments) -> Result<Hopscotch<Kmer>, Failure> {
    let k = arguments.k_size;
    let Some(path) = &arguments.kmers_to_ignore_file else {
        return Err(Failure::Thrown {
            class: "java.lang.IllegalArgumentException",
            message: "provided path for file containing kmers is null".to_string(),
        });
    };
    if k <= 0 {
        return Err(Failure::Thrown {
            class: ILLEGAL_ARGUMENT,
            message: format!("provided k-size is non positive: {k}"),
        });
    }
    let bytes = read_file(path)?;
    let mut kills = Hopscotch::new((bytes.len() as i64 / (i64::from(k) + 1)) as i32);
    let gatk = |message: String| Failure::Thrown {
        class: GATK_EXCEPTION,
        message,
    };
    let text = String::from_utf8_lossy(&bytes);
    for line in java_lines(&text) {
        if line.chars().count() != k as usize {
            return Err(gatk(format!(
                "SVKmer kill set contains a line of length {} but we were expecting K = {k}",
                line.chars().count()
            )));
        }
        if !(1..64).contains(&k) {
            return Err(Failure::Thrown {
                class: ILLEGAL_ARGUMENT,
                message: "Kmer length must be between 1 and 63.".to_string(),
            });
        }
        let Some(kmer) = Kmerizer::new(line.as_bytes(), k, 1).next() else {
            return Err(gatk(format!(
                "Unable to kmerize the kmer kill set string '{line}'."
            )));
        };
        kills.add(Kmer(kmer))?;
    }
    if let Some(adapter) = &arguments.adapter_sequence {
        for kmer in Kmerizer::new(adapter.as_bytes(), k, 0) {
            if k & 1 == 0 {
                return Err(Failure::Thrown {
                    class: ILLEGAL_ARGUMENT,
                    message: "Kmer length must be odd to canonicalize.".to_string(),
                });
            }
            kills.add(Kmer(kmer.canonical(k)))?;
        }
    }
    Ok(kills)
}

/// `BufferedReader.readLine` over a whole file: `\n`, `\r` and `\r\n`, and no line after a last
/// terminator.
fn java_lines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let bytes = text.as_bytes();
    let (mut from, mut at) = (0, 0);
    while at < bytes.len() {
        if bytes[at] == b'\n' || bytes[at] == b'\r' {
            lines.push(&text[from..at]);
            if bytes[at] == b'\r' && bytes.get(at + 1) == Some(&b'\n') {
                at += 1;
            }
            at += 1;
            from = at;
        } else {
            at += 1;
        }
    }
    if from < bytes.len() {
        lines.push(&text[from..]);
    }
    lines
}

/// `String.split("\t")`, which drops what is empty at the end.
fn java_split_tabs(line: &str) -> Vec<&str> {
    let mut tokens: Vec<&str> = line.split('\t').collect();
    while tokens.len() > 1 && tokens.last().is_some_and(|token| token.is_empty()) {
        tokens.pop();
    }
    tokens
}

/// `SVFileUtils.readIntervalsFile`.
fn read_intervals_file(path: &str, metadata: &ReadMetadata) -> Result<Vec<SvInterval>, Failure> {
    let bytes = read_file(path)?;
    let text = String::from_utf8_lossy(&bytes);
    let mut intervals = Vec::new();
    for (index, line) in java_lines(&text).into_iter().enumerate() {
        let number = index + 1;
        if line.starts_with('#') {
            continue;
        }
        let tokens = java_split_tabs(line);
        if tokens.len() != 3 {
            return Err(Failure::Thrown {
                class: GATK_EXCEPTION,
                message: format!(
                    "Interval file {path} line {number} did not contain 3 columns: {line}"
                ),
            });
        }
        let parsed = (|| {
            let contig = *metadata.contig_ids.get(&tokens[0].to_string())?;
            let start: i32 = tokens[1].parse().ok()?;
            let end: i32 = tokens[2].parse().ok()?;
            (start >= 0 && end >= 0).then_some(SvInterval { contig, start, end })
        })();
        match parsed {
            Some(interval) => intervals.push(interval),
            None => {
                return Err(Failure::Thrown {
                    class: GATK_EXCEPTION,
                    message: format!(
                        "Unable to parse interval file {path} line {number}: {line}"
                    ),
                })
            }
        }
    }
    Ok(intervals)
}

/// `removeIntervalsNearGaps`, through a `TreeSet`'s neighbours on either side.
fn remove_intervals_near_gaps(
    intervals: Vec<SvInterval>,
    padding: i32,
    path: &str,
    metadata: &ReadMetadata,
) -> Result<Vec<SvInterval>, Failure> {
    let gaps: std::collections::BTreeSet<SvInterval> =
        read_intervals_file(path, metadata)?.into_iter().collect();
    Ok(intervals
        .into_iter()
        .filter(|interval| {
            if let Some(gap) = gaps.range(..*interval).next_back() {
                if gap.gap_len(interval) < padding {
                    return false;
                }
            }
            if let Some(gap) = gaps.range(*interval..).next() {
                if interval.gap_len(gap) < padding {
                    return false;
                }
            }
            true
        })
        .collect())
}

/// `QNameFinder` over the partition: `intervalsIndex` only ever moves forward.
fn mapped_qnames(
    records: &[BamRecord],
    intervals: &[SvInterval],
    ignored: &[SvInterval],
) -> Vec<QNameAndInterval> {
    let mut found = Vec::new();
    let mut index = 0;
    for record in records {
        if !is_mapped(record) {
            continue;
        }
        let contig = record.reference_index;
        while index < intervals.len() {
            let interval = &intervals[index];
            if interval.contig > contig
                || (interval.contig == contig && interval.end > record.alignment_start)
            {
                break;
            }
            index += 1;
        }
        if index >= intervals.len() {
            continue;
        }
        let unclipped = SvInterval {
            contig,
            start: read_utils::unclipped_start(record).max(0),
            end: read_utils::unclipped_end(record),
        };
        if !intervals[index].overlaps(&unclipped) {
            continue;
        }
        let clipped = SvInterval {
            contig,
            start: record.alignment_start,
            end: read_utils::end(record),
        };
        if contained_in_region_to_ignore(&clipped, ignored) {
            continue;
        }
        found.push(QNameAndInterval {
            name: record.read_name.clone(),
            interval: index as i32,
        });
    }
    found
}

fn write_qnames(path: &str, table: &Hopscotch<QNameAndInterval>) -> Result<(), Failure> {
    let mut text = String::new();
    for entry in table.iter() {
        text.push_str(&format!("{} {}\n", entry.name, entry.interval));
    }
    write_file(path, text.as_bytes())
}

/// `getKmerIntervals`: the k-mers of each interval's templates, cleaned, and the intervals with too
/// few of them excused.
fn kmer_intervals(
    records: &[BamRecord],
    arguments: &Arguments,
    metadata: &ReadMetadata,
    qnames: &mut Hopscotch<QNameAndInterval>,
    n_intervals: usize,
    kills: &Hopscotch<Kmer>,
    dispositions: &mut Vec<(i32, Disposition)>,
) -> Result<Vec<KmerAndInterval>, Failure> {
    let k = arguments.k_size;
    let kmers_per_partition = metadata.ref_bases.wrapping_mul(2) as i32;
    let diploid = accurate_kmer_coverage(metadata, k);
    let haploid = diploid / 2.0;
    let deviation = f64::from(haploid).sqrt() as f32;
    let min_coverage = haploid - 3.0 * deviation;
    let min_kmers = arguments.cleaner_min_kmer_count.max(java_round(min_coverage));
    let max_kmers = java_round(arguments.cleaner_max_copy_number as f32 * diploid);

    // `QNameKmerizer`, read by read.
    let mut pairs: Vec<(KmerAndInterval, i32)> = Vec::new();
    for record in records {
        if !not_junk(record) || !is_primary_line(record) {
            continue;
        }
        let name = &record.read_name;
        let hash = string_hash_code(name);
        let found: Vec<i32> = qnames
            .find_each(hash, |entry| entry.name == *name)
            .iter()
            .map(|entry| entry.interval)
            .collect();
        for interval in found {
            if k < 4 {
                return Err(Failure::Task(format!(
                    "{GATK_EXCEPTION}: kmer size must be at least 4 for this filter to work \
                     properly."
                )));
            }
            for kmer in DustFilteredKmerizer::new(&record.read_bases, k, arguments.max_dust_score) {
                let kmer = canonical(kmer, k)?;
                if kills.contains(kmer.hash_code(), |held| held.0 == kmer) {
                    continue;
                }
                pairs.push((KmerAndInterval { kmer, interval }, 1));
            }
        }
    }
    let reduced = reduce_by_key(
        pairs,
        |key: &KmerAndInterval| key.kmer.hash_code(),
        |one: i32, two: i32| one.wrapping_add(two),
    );

    // `KmerCleaner`.
    if i64::from(kmers_per_partition) > MAX_TABLE {
        return Err(Failure::Limitation(format!(
            "KmerCleaner sizes its table from the reference, here for {kmers_per_partition} \
             entries, and this port stops short of allocating it."
        )));
    }
    let mut table: Hopscotch<KmerAndInterval> = Hopscotch::new(kmers_per_partition);
    for (key, count) in reduced {
        if count >= min_kmers && count <= max_kmers {
            table.add(key)?;
        }
    }
    let mut unique: Hopscotch<Kmer> = Hopscotch::new(table.len() as i32);
    for entry in table.iter() {
        unique.add(Kmer(entry.kmer))?;
    }
    let ubiquitous: Vec<SvKmerLong> = unique
        .iter()
        .filter(|kmer| {
            table
                .find_each(kmer.0.hash_code(), |entry| entry.kmer == kmer.0)
                .len() as i64
                > i64::from(arguments.cleaner_max_intervals)
        })
        .map(|kmer| kmer.0)
        .collect();
    for kmer in ubiquitous {
        table.remove_each(kmer.hash_code(), |entry| entry.kmer == kmer);
    }
    let kmer_intervals: Vec<KmerAndInterval> = table.iter().cloned().collect();

    let mut counts = vec![0i32; n_intervals];
    for entry in &kmer_intervals {
        counts[entry.interval as usize] += 1;
    }
    let mut kill = HashSet::new();
    for (index, count) in counts.iter().enumerate() {
        if *count < arguments.min_kmers_per_interval {
            kill.insert(index as i32);
            dispositions.push((
                index as i32,
                Disposition::Excuse("FASTQ not written -- too few kmers".to_string()),
            ));
        }
    }
    qnames.remove_if(|entry| kill.contains(&entry.interval));
    let filtered: Vec<KmerAndInterval> = kmer_intervals
        .into_iter()
        .filter(|entry| !kill.contains(&entry.interval))
        .collect();
    if let Some(path) = &arguments.kmer_file {
        let mut text = String::new();
        for entry in &filtered {
            text.push_str(&format!("{} {}\n", entry.kmer.to_text(k), entry.interval));
        }
        write_file(path, text.as_bytes())?;
    }
    Ok(filtered)
}

/// `removeUbiquitousKmers`.
fn remove_ubiquitous_kmers(
    records: &[BamRecord],
    arguments: &Arguments,
    metadata: &ReadMetadata,
    table: &mut Hopscotch<KmerAndInterval>,
) -> Result<(), Failure> {
    let k = arguments.k_size;
    let span = metadata.partition.span;
    if span == 0 {
        return Err(Failure::Thrown {
            class: "java.lang.ArithmeticException",
            message: "/ by zero".to_string(),
        });
    }
    let per_locus = (span
        .wrapping_add(metadata.avg_read_len)
        .wrapping_sub(k)
        .wrapping_div(span)) as f32;
    let max_partition_count = (arguments.cleaner_max_copy_number as f32 * per_locus) as i32;
    let max_count =
        (arguments.cleaner_max_copy_number as f32 * accurate_kmer_coverage(metadata, k)) as i32;

    // `KmerCounter`.
    let mut counts: Hopscotch<KmerAndCount> = Hopscotch::new(table.len() as i32);
    for record in records {
        if !not_junk(record) || !is_primary_line(record) {
            continue;
        }
        for kmer in Kmerizer::new(&record.read_bases, k, 1) {
            let kmer = canonical(kmer, k)?;
            let hash = kmer.hash_code();
            if !table.contains(hash, |entry| entry.kmer == kmer) {
                continue;
            }
            match counts.find_mut(hash, |entry| entry.kmer == kmer) {
                Some(entry) => entry.count = entry.count.wrapping_add(1),
                None => {
                    counts.add(KmerAndCount { kmer, count: 1 })?;
                }
            }
        }
    }
    let pairs: Vec<(SvKmerLong, (i32, i32))> = counts
        .iter()
        .map(|entry| (entry.kmer, (1, entry.count)))
        .collect();
    let reduced = reduce_by_key(
        pairs,
        |kmer: &SvKmerLong| kmer.hash_code(),
        |one: (i32, i32), two: (i32, i32)| {
            (one.0.wrapping_add(two.0), one.1.wrapping_add(two.1))
        },
    );
    for (kmer, (partitions, count)) in reduced {
        if partitions > max_partition_count || count > max_count {
            table.remove_each(kmer.hash_code(), |entry| entry.kmer == kmer);
        }
    }
    Ok(())
}

/// `QNameIntervalFinder`, whose intervals come in the order the read's k-mers reach them.
fn assembly_qnames(
    records: &[BamRecord],
    k: i32,
    table: &Hopscotch<KmerAndInterval>,
) -> Result<Vec<QNameAndInterval>, Failure> {
    let mut found = Vec::new();
    for record in records {
        if !not_junk(record) || !is_primary_line(record) {
            continue;
        }
        let mut intervals: Vec<i32> = Vec::new();
        for kmer in Kmerizer::new(&record.read_bases, k, 1) {
            let kmer = canonical(kmer, k)?;
            for entry in table.find_each(kmer.hash_code(), |entry| entry.kmer == kmer) {
                if !intervals.contains(&entry.interval) {
                    intervals.push(entry.interval);
                }
            }
        }
        for interval in intervals {
            found.push(QNameAndInterval {
                name: record.read_name.clone(),
                interval,
            });
        }
    }
    Ok(found)
}

/// `SVFastqUtils.FastqRead`.
struct FastqRead {
    header: String,
    bases: Vec<u8>,
    quals: Vec<u8>,
}

/// `SequenceUtil.complement`.
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        other => other,
    }
}

fn reverse_complement(bases: &mut [u8]) {
    bases.reverse();
    for base in bases {
        *base = complement(*base);
    }
}

/// `SVFastqUtils.Mapping.toString(GATKRead)`.
fn mapping_text(record: &BamRecord, contigs: &[String]) -> String {
    let padded: u32 = record
        .cigar
        .elements
        .iter()
        .filter(|element| element.op.consumes_reference_bases() || element.op == CigarOp::P)
        .map(|element| element.length)
        .sum();
    if unmapped(record) || padded == 0 {
        return "*".to_string();
    }
    let int_tag = |tag: &[u8; 2]| match record.tags.get(Tag::new(tag)) {
        Some(TagValue::Int(value)) => Some(*value),
        _ => None,
    };
    let mut text = format!(
        "{},{},{},{},{}",
        contigs[record.reference_index as usize],
        record.alignment_start,
        if record.flags & REVERSE != 0 { "-" } else { "+" },
        record.cigar.to_text(),
        record.mapping_quality
    );
    match (int_tag(b"NM"), int_tag(b"AS")) {
        (None, None) => {}
        (None, Some(score)) => text.push_str(&format!(",.,{score}")),
        (Some(mismatches), None) => text.push_str(&format!(",{mismatches}")),
        (Some(mismatches), Some(score)) => text.push_str(&format!(",{mismatches},{score}")),
    }
    if let Some(TagValue::Str(sa)) = record.tags.get(Tag::new(b"SA")) {
        text.push(';');
        text.push_str(sa);
    }
    if text.ends_with(';') {
        text.pop();
    }
    text
}

fn unmapped(record: &BamRecord) -> bool {
    record.flags & 0x4 != 0 || record.reference_index < 0 || record.alignment_start == 0
}

impl FastqRead {
    fn new(record: &BamRecord, metadata: &ReadMetadata, include_mapping: bool) -> FastqRead {
        let mut header = format!("@{}{}", record.read_name, fragment_ordinal(record));
        if include_mapping {
            header.push_str("\tmapping=");
            header.push_str(&mapping_text(record, &metadata.contigs));
        }
        let mut bases = record.read_bases.clone();
        let mut quals = record.base_qualities.clone();
        if !unmapped(record) && record.flags & REVERSE != 0 {
            reverse_complement(&mut bases);
            quals.reverse();
        }
        FastqRead {
            header,
            bases,
            quals,
        }
    }
}

/// `SVFastqUtils.writeFastqFile`.
fn write_fastq(path: &str, reads: &[&FastqRead]) -> Result<(), Failure> {
    let mut bytes = Vec::new();
    for read in reads {
        if read.header.contains(' ') {
            return Err(Failure::Task(format!(
                "{ILLEGAL_STATE}: Blank found: {}",
                read.header
            )));
        }
        bytes.extend_from_slice(read.header.as_bytes());
        bytes.push(b'\n');
        bytes.extend_from_slice(&read.bases);
        bytes.extend_from_slice(b"\n+\n");
        for quality in &read.quals {
            // `SAMUtils.phredToFastq`.
            if *quality as i8 > 93 || (*quality as i8) < 0 {
                return Err(Failure::Task(format!(
                    "{ILLEGAL_ARGUMENT}: Cannot encode phred score: {}",
                    *quality as i8
                )));
            }
            bytes.push(quality + 33);
        }
        bytes.push(b'\n');
    }
    write_file(path, &bytes)
}

/// `FermiLiteAssembly.Connection`. Java compares these by identity, so each has an id of its own.
#[derive(Clone, Copy)]
struct Connection {
    target: usize,
    overlap: i32,
    is_rc: bool,
    is_target_rc: bool,
    id: u64,
}

/// `FermiLiteAssembly.Contig`.
struct Contig {
    sequence: Vec<u8>,
    supporting_reads: i32,
    connections: Vec<Connection>,
}

/// Every contig an assembly's revision makes, by index: contigs are compared by identity too.
struct Graph {
    contigs: Vec<Contig>,
    next_connection: u64,
}

impl Graph {
    fn connection(&mut self, target: usize, overlap: i32, is_rc: bool, is_target_rc: bool) -> Connection {
        self.next_connection += 1;
        Connection {
            target,
            overlap,
            is_rc,
            is_target_rc,
            id: self.next_connection,
        }
    }

    /// `rcConnection`, called with the contig among whose connections this one appears.
    fn reversed(&mut self, connection: &Connection, contig: usize) -> Connection {
        self.connection(
            contig,
            connection.overlap,
            !connection.is_target_rc,
            !connection.is_rc,
        )
    }

    /// `getSingletonConnection`.
    fn singleton(&self, contig: usize, is_rc: bool) -> Option<Connection> {
        let mut singleton = None;
        for connection in &self.contigs[contig].connections {
            if connection.is_rc == is_rc {
                if singleton.is_some() {
                    return None;
                }
                singleton = Some(*connection);
            }
        }
        singleton
    }

    /// `removeShadowedContigs`: a contig another contig covers with few mismatches goes. The order
    /// a k-mer's locations come in cannot show, since every one of them is tried.
    fn remove_shadowed(&mut self, assembly: Vec<usize>) -> Result<Vec<usize>, Failure> {
        const K: i32 = 31;

        #[derive(Clone)]
        struct KmerLocation {
            kmer: SvKmerShort,
            contig: usize,
            offset: i32,
            canonical: bool,
        }
        impl HopscotchEntry for KmerLocation {
            fn key_hash(&self) -> i32 {
                self.kmer.hash_code()
            }
            fn collides(&self, _: &Self) -> bool {
                false
            }
        }

        let capacity = assembly.iter().fold(0i32, |sum, contig| {
            sum.wrapping_add(self.contigs[*contig].sequence.len() as i32 - K + 1)
        });
        let mut map: Hopscotch<KmerLocation> = Hopscotch::new(capacity);
        for contig in &assembly {
            for (offset, kmer) in SvKmerShort::kmers(&self.contigs[*contig].sequence, K)
                .into_iter()
                .enumerate()
            {
                let canonical = kmer.canonical(K);
                map.add(KmerLocation {
                    kmer: canonical,
                    contig: *contig,
                    offset: offset as i32,
                    canonical: kmer == canonical,
                })?;
            }
        }

        let mut to_remove: HashSet<usize> = HashSet::new();
        for contig in &assembly {
            let contig = *contig;
            let bases = &self.contigs[contig].sequence;
            let max_mismatches = (bases.len() as f64 * 0.01) as i32;
            let mut tested: HashSet<(usize, i32, bool)> = HashSet::new();
            for (offset, kmer) in SvKmerShort::kmers(bases, K).into_iter().enumerate() {
                let offset = offset as i32;
                let canonical = kmer.canonical(K);
                let is_canonical = kmer == canonical;
                for location in map.find_each(canonical.hash_code(), |entry| entry.kmer == canonical)
                {
                    let other = location.contig;
                    if contig == other || to_remove.contains(&other) {
                        continue;
                    }
                    let other_bases = &self.contigs[other].sequence;
                    let is_rc = is_canonical != location.canonical;
                    let other_offset = if is_rc {
                        other_bases.len() as i32 - location.offset - K
                    } else {
                        location.offset
                    };
                    if offset > other_offset {
                        continue;
                    }
                    if bases.len() as i32 - offset > other_bases.len() as i32 - other_offset {
                        continue;
                    }
                    let other_start = other_offset - offset;
                    if !tested.insert((other, other_start, is_rc)) {
                        continue;
                    }
                    let mut mismatches = 0;
                    for (at, base) in bases.iter().enumerate() {
                        let differs = if is_rc {
                            let from = other_bases.len() - other_start as usize - 1 - at;
                            *base != simple_complement(other_bases[from])
                        } else {
                            *base != other_bases[other_start as usize + at]
                        };
                        if differs {
                            mismatches += 1;
                            if mismatches > max_mismatches {
                                break;
                            }
                        }
                    }
                    if mismatches <= max_mismatches {
                        to_remove.insert(contig);
                        break;
                    }
                }
            }
        }

        let surviving: Vec<usize> = assembly
            .iter()
            .copied()
            .filter(|contig| !to_remove.contains(contig))
            .collect();
        let mut stale: HashSet<usize> = HashSet::new();
        for contig in &to_remove {
            for connection in &self.contigs[*contig].connections {
                stale.insert(connection.target);
            }
        }
        for contig in stale {
            let connections = &mut self.contigs[contig].connections;
            if connections.is_empty() {
                return Err(Failure::Task(format!(
                    "{ILLEGAL_ARGUMENT}: Illegal Capacity: -1"
                )));
            }
            connections.retain(|connection| !to_remove.contains(&connection.target));
        }
        Ok(surviving)
    }

    /// `removeUnbranchedConnections`.
    fn remove_unbranched(&mut self, assembly: Vec<usize>) -> Result<Vec<usize>, Failure> {
        let mut list = Vec::with_capacity(assembly.len());
        let mut examined: HashSet<usize> = HashSet::new();
        for contig in assembly {
            if !examined.insert(contig) {
                continue;
            }
            let mut contig = contig;
            for is_rc in [true, false] {
                while let Some(connection) = self.singleton(contig, is_rc) {
                    if examined.contains(&connection.target) {
                        break;
                    }
                    let Some(back) = self.singleton(connection.target, !connection.is_target_rc)
                    else {
                        break;
                    };
                    examined.insert(connection.target);
                    contig = self.join_with_connections(contig, &connection, &back)?;
                }
            }
            list.push(contig);
        }
        Ok(list)
    }

    /// `joinContigsWithConnections`.
    fn join_with_connections(
        &mut self,
        first: usize,
        connection: &Connection,
        back: &Connection,
    ) -> Result<usize, Failure> {
        let joined = self.join(first, std::slice::from_ref(connection))?;
        let last = connection.target;
        let mut connections = Vec::new();
        for (from, skip, is_rc) in [(first, connection.id, true), (last, back.id, false)] {
            for conn in self.contigs[from].connections.clone() {
                if conn.id == skip {
                    continue;
                }
                let fresh = self.connection(conn.target, conn.overlap, is_rc, conn.is_target_rc);
                let old = self.reversed(&conn, from);
                let new = self.reversed(&fresh, joined);
                self.replace_connection(conn.target, &old, new);
                connections.push(fresh);
            }
        }
        self.contigs[joined].connections = connections;
        Ok(joined)
    }

    /// `replaceConnection`.
    fn replace_connection(&mut self, contig: usize, old: &Connection, new: Connection) {
        for conn in &mut self.contigs[contig].connections {
            if conn.target == old.target
                && conn.is_rc == old.is_rc
                && conn.is_target_rc == old.is_target_rc
            {
                *conn = new;
            }
        }
    }

    /// `joinContigs`.
    fn join(&mut self, first: usize, path: &[Connection]) -> Result<usize, Failure> {
        if path.is_empty() {
            return Ok(first);
        }
        let mut sequence = self.contigs[first].sequence.clone();
        let mut reads = self.contigs[first].supporting_reads;
        if path[0].is_rc {
            reverse_complement(&mut sequence);
        }
        for connection in path {
            let target = &self.contigs[connection.target];
            reads = reads.wrapping_add(target.supporting_reads);
            let overlap = connection.overlap;
            let length = target.sequence.len() as i32 - overlap;
            if overlap < 0 && connection.is_target_rc || length < 0 {
                return Err(Failure::Limitation(
                    "an assembly joined across a connection its contig cannot hold.".to_string(),
                ));
            }
            if connection.is_target_rc {
                let mut piece = target.sequence[..length as usize].to_vec();
                reverse_complement(&mut piece);
                sequence.extend_from_slice(&piece);
            } else {
                if overlap < 0 {
                    return Err(Failure::Limitation(
                        "an assembly joined across a connection its contig cannot hold."
                            .to_string(),
                    ));
                }
                sequence.extend_from_slice(&target.sequence[overlap as usize..]);
            }
        }
        self.contigs.push(Contig {
            sequence,
            supporting_reads: reads,
            connections: Vec::new(),
        });
        Ok(self.contigs.len() - 1)
    }

    fn predecessors(&self, contig: usize) -> usize {
        self.contigs[contig]
            .connections
            .iter()
            .filter(|connection| connection.is_rc)
            .count()
    }

    /// `expandAssemblyGraph`: every path from every source and sink.
    fn expand(&mut self, assembly: Vec<usize>) -> Result<Vec<usize>, Failure> {
        let mut walk = Walk {
            list: Vec::with_capacity(assembly.len()),
            examined: HashSet::new(),
            visited: HashSet::new(),
        };
        for contig in &assembly {
            let contig = *contig;
            if walk.examined.contains(&contig) {
                continue;
            }
            if self.contigs[contig].connections.is_empty() {
                walk.list.push(contig);
                walk.examined.insert(contig);
            } else {
                let predecessors = self.predecessors(contig);
                let successors = self.contigs[contig].connections.len() - predecessors;
                if predecessors == 0 {
                    self.trace_paths(contig, false, &mut walk)?;
                } else if successors == 0 {
                    self.trace_paths(contig, true, &mut walk)?;
                }
            }
        }
        // What is left is a smooth cycle.
        for contig in &assembly {
            if !walk.examined.contains(contig) {
                self.trace_paths(*contig, false, &mut walk)?;
            }
        }
        Ok(walk.list)
    }

    /// `tracePaths`.
    fn trace_paths(&mut self, contig: usize, is_rc: bool, walk: &mut Walk) -> Result<(), Failure> {
        walk.examined.insert(contig);
        let mut path = Vec::new();
        let is_cycle = !walk.visited.insert((contig, is_rc));
        for connection in self.contigs[contig].connections.clone() {
            if connection.is_rc == is_rc && connection.overlap >= 0 {
                self.extend_path(contig, connection, &mut path, walk)?;
            }
        }
        if !is_cycle {
            walk.visited.remove(&(contig, is_rc));
        }
        Ok(())
    }

    /// `extendPath`.
    fn extend_path(
        &mut self,
        first: usize,
        connection: Connection,
        path: &mut Vec<Connection>,
        walk: &mut Walk,
    ) -> Result<(), Failure> {
        path.push(connection);
        let strand = (connection.target, connection.is_target_rc);
        let is_cycle = !walk.visited.insert(strand);
        let mut at_end = true;
        if !is_cycle {
            let target = connection.target;
            let predecessors = self.predecessors(target);
            let successors = self.contigs[target].connections.len() - predecessors;
            if predecessors > 1 && successors > 1 {
                if walk.examined.insert(target) {
                    self.trace_paths(target, connection.is_target_rc, walk)?;
                }
            } else {
                walk.examined.insert(target);
                for next in self.contigs[target].connections.clone() {
                    if next.is_rc == connection.is_target_rc && next.overlap >= 0 {
                        self.extend_path(first, next, path, walk)?;
                        at_end = false;
                    }
                }
            }
        }
        if at_end {
            let joined = self.join(first, path)?;
            walk.list.push(joined);
        }
        if !is_cycle {
            walk.visited.remove(&strand);
        }
        path.pop();
        Ok(())
    }

    /// `FermiLiteAssembly.writeGFA`.
    fn gfa(&self, assembly: &[usize]) -> Result<String, Failure> {
        let ids: HashMap<usize, usize> = assembly
            .iter()
            .enumerate()
            .map(|(id, contig)| (*contig, id))
            .collect();
        let mut text = String::from("H\tVN:Z:1.0\n");
        for contig in assembly {
            let id = ids[contig];
            let tig = &self.contigs[*contig];
            text.push_str(&format!(
                "S\ttig{id}\t{}\tLN:i:{}\tRC:i:{}\n",
                String::from_utf8_lossy(&tig.sequence),
                tig.sequence.len(),
                tig.supporting_reads
            ));
            for connection in &tig.connections {
                let Some(target) = ids.get(&connection.target) else {
                    return Err(Failure::Task(
                        "java.lang.NullPointerException: Cannot invoke \
                         \"java.lang.Integer.intValue()\" because the return value of \
                         \"java.util.HashMap.get(Object)\" is null"
                            .to_string(),
                    ));
                };
                if id <= *target {
                    let overlap = connection.overlap;
                    text.push_str(&format!(
                        "L\ttig{id}\t{}\ttig{target}\t{}\t{}\n",
                        if connection.is_rc { "-" } else { "+" },
                        if connection.is_target_rc { "-" } else { "+" },
                        if overlap < 0 {
                            format!("{}H", -overlap)
                        } else {
                            format!("{overlap}M")
                        }
                    ));
                }
            }
        }
        Ok(text)
    }
}

/// `BaseUtils.simpleComplement`.
fn simple_complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'.',
    }
}

struct Walk {
    list: Vec<usize>,
    examined: HashSet<usize>,
    visited: HashSet<(usize, bool)>,
}

/// `AlignedAssemblyOrExcuse`.
#[derive(Debug, Clone)]
pub enum Disposition {
    Excuse(String),
    Aligned {
        contigs: Vec<Vec<u8>>,
        seconds: i32,
        alignments: Vec<Vec<ContigAlignment>>,
    },
}

/// `AssembledEvidenceResults`.
pub struct Results {
    /// Whether `gather` got past its two early returns, which write neither the interval file nor
    /// the contigs.
    pub finished: bool,
    pub intervals: Vec<SvInterval>,
    pub assemblies: Vec<(i32, Disposition)>,
    pub links: Vec<TargetLink>,
}

/// `FermiLiteAssemblyHandler.apply`.
fn assemble_interval(
    id: i32,
    reads: &[std::rc::Rc<FastqRead>],
    arguments: &Arguments,
    aligner: &mut Option<ContigAligner>,
) -> Result<Disposition, Failure> {
    let name = format!("asm{id:06}");
    let size = reads
        .iter()
        .fold(0i32, |sum, read| sum.wrapping_add(read.bases.len() as i32));
    if size > arguments.max_fastq_size {
        return Ok(Disposition::Excuse(format!(
            "no assembly -- too big ({size} bytes)."
        )));
    }
    if let Some(directory) = &arguments.fastq_dir {
        let mut sorted: Vec<&FastqRead> = reads.iter().map(|read| read.as_ref()).collect();
        sorted.sort_by(|one, two| compare_strings(&one.header, &two.header));
        write_fastq(&format!("{directory}/{name}.fastq"), &sorted)?;
    }

    let mut options = FmlOpt::default();
    if arguments.pop_variant_bubbles {
        options = options.with_cleaning_flag(MAG_F_AGGRESSIVE | MAG_F_POPOPEN);
    }
    let pairs: Vec<(&[u8], &[u8])> = reads
        .iter()
        .map(|read| (read.bases.as_slice(), read.quals.as_slice()))
        .collect();
    let started = Instant::now();
    let initial = gatk_fermilite::assemble(&options, &pairs).map_err(|error| match error {
        FmlError::Abort(assertion) => Failure::Limitation(format!(
            "fermi-lite aborts on interval {id} ({assertion}), which ends the JVM in the \
             reference, so there is no output to be the same as."
        )),
        FmlError::Undefined(what) => Failure::Limitation(format!(
            "the assembly of interval {id} reads memory fermi-lite does not own ({what}), so the \
             reference's answer depends on its heap."
        )),
        FmlError::JavaBufferOverflow => {
            Failure::Task("java.nio.BufferOverflowException".to_string())
        }
    })?;
    let seconds = ((started.elapsed().as_millis() as i64 + 500) / 1000) as i32;
    if initial.is_empty() {
        return Ok(Disposition::Excuse(
            "no assembly -- no contigs produced by assembler.".to_string(),
        ));
    }

    let mut graph = Graph {
        contigs: Vec::with_capacity(initial.len()),
        next_connection: 0,
    };
    for contig in &initial {
        graph.contigs.push(Contig {
            sequence: contig.seq.clone(),
            supporting_reads: contig.nsr,
            connections: Vec::new(),
        });
    }
    for (index, contig) in initial.iter().enumerate() {
        let connections = contig
            .ovlp
            .iter()
            .map(|overlap| {
                graph.connection(
                    overlap.target(),
                    overlap.overlap_len(),
                    overlap.is_rc(),
                    overlap.is_target_rc(),
                )
            })
            .collect();
        graph.contigs[index].connections = connections;
    }
    let mut assembly: Vec<usize> = (0..initial.len()).collect();
    if arguments.remove_shadowed_contigs {
        assembly = graph.remove_shadowed(assembly)?;
    }
    if arguments.expand_assembly_graph {
        assembly = graph.remove_unbranched(assembly)?;
        assembly = graph.expand(assembly)?;
    }
    if arguments.fastq_dir.is_some() && arguments.write_gfas {
        let directory = arguments.fastq_dir.as_deref().expect("just tested");
        write_file(
            &format!("{directory}/{name}.gfa"),
            graph.gfa(&assembly)?.as_bytes(),
        )?;
    }

    if aligner.is_none() {
        let path = arguments.aligner_index_image.as_deref().unwrap_or("");
        let image = std::fs::read(path)
            .ok()
            .filter(|bytes| !bytes.is_empty())
            .ok_or_else(|| {
                Failure::Task(format!(
                    "org.broadinstitute.hellbender.utils.bwa.BwaMemException: cannot read the \
                     image file '{path}': is empty or is not readable"
                ))
            })?;
        *aligner = Some(
            ContigAligner::from_image(&image, arguments.z_dropoff).map_err(Failure::Port)?,
        );
    }
    let contigs: Vec<Vec<u8>> = assembly
        .iter()
        .map(|contig| graph.contigs[*contig].sequence.clone())
        .collect();
    let alignments = aligner.as_ref().expect("just opened").align(&contigs);
    Ok(Disposition::Aligned {
        contigs,
        seconds,
        alignments,
    })
}

/// `getIntervalsAndEvidenceTargetLinks`.
fn intervals_and_links(
    records: &[BamRecord],
    arguments: &Arguments,
    metadata: &ReadMetadata,
    external: Vec<Evidence>,
    ignored: &[SvInterval],
) -> Result<(Vec<SvInterval>, Vec<TargetLink>), Failure> {
    let evidence = sv_evidence::evidence(records, metadata, ignored, &arguments.evidence.metadata)
        .map_err(|message| Failure::Task(format!("{ILLEGAL_STATE}: {message}")))?;
    if let Some(directory) = &arguments.unfiltered_evidence_dir {
        save_as_text_file(directory, evidence.iter().map(Evidence::string_rep))?;
    }

    let mut links = Vec::new();
    for link in cluster_links(&evidence) {
        if link.read_pairs >= 2 || link.split_reads >= 1 {
            links.push(TargetLink {
                read_pair_names: link.read_pair_names.collected()?,
                split_read_names: link.split_read_names.collected()?,
                ..link
            });
        }
    }
    if let Some(path) = &arguments.target_link_file {
        let mut text = String::new();
        for link in &links {
            text.push_str(&link.bedpe(metadata));
            text.push('\n');
        }
        write_file(path, text.as_bytes())?;
    }

    let coverage = f64::from(metadata.coverage);
    let min_weight = arguments.min_evidence_weight_per_coverage * coverage;
    let min_coherent = arguments.min_coherent_evidence_weight_per_coverage * coverage;
    let mut stream = evidence;
    stream.extend(external);
    let filtered = density_filter(stream, min_weight, min_coherent);
    if let Some(directory) = &arguments.evidence_dir {
        save_as_text_file(
            directory,
            filtered.iter().filter(|one| one.validated).map(Evidence::to_text),
        )?;
    }

    let max_fragment = metadata
        .libraries
        .iter()
        .map(|(_, statistics)| statistics.median)
        .max()
        .unwrap_or(0);
    let contigs = metadata.contigs.len() as i32;
    let collected = cluster_evidence(&filtered, max_fragment, contigs);
    let again = density_filter(collected, min_weight, min_coherent);
    if let Some(directory) = &arguments.evidence_dir {
        // The cross-partition file, which one partition has nothing to put in.
        write_file(&format!("{directory}/part-xxxxx"), b"")?;
    }
    let intervals = cluster_evidence(&again, max_fragment, contigs)
        .into_iter()
        .map(|one| one.location)
        .collect();
    Ok((intervals, links))
}

/// `AlignedAssemblyOrExcuse.writeIntervalFile`.
fn write_interval_file(
    path: &str,
    intervals: &[SvInterval],
    assemblies: &[(i32, Disposition)],
    header: &SamHeader,
) -> Result<(), Failure> {
    let mut text = String::new();
    for (index, interval) in intervals.iter().enumerate() {
        let disposition = match assemblies.iter().find(|(id, _)| *id == index as i32) {
            None => "unknown".to_string(),
            Some((_, Disposition::Excuse(excuse))) => excuse.clone(),
            Some((
                _,
                Disposition::Aligned {
                    contigs, seconds, ..
                },
            )) => format!("produced {} contigs in {seconds} secs.", contigs.len()),
        };
        text.push_str(&format!(
            "{index}\t{}:{}-{}\t{disposition}\n",
            header.sequences[interval.contig as usize].name, interval.start, interval.end
        ));
    }
    write_file(path, text.as_bytes())
}

/// `gatherEvidenceAndWriteContigSamFile`, past the metadata, the external evidence and the
/// high-coverage intervals, up to the file of aligned contigs the caller writes.
pub fn gather(
    records: &[BamRecord],
    arguments: &Arguments,
    metadata: &ReadMetadata,
    header: &SamHeader,
    external: Vec<Evidence>,
    ignored: &[SvInterval],
) -> Result<Results, Failure> {
    let (intervals, links) = intervals_and_links(records, arguments, metadata, external, ignored)?;
    if intervals.is_empty() {
        return Ok(Results {
            finished: false,
            intervals,
            assemblies: Vec::new(),
            links,
        });
    }

    let mut intervals = intervals;
    if let Some(path) = &arguments.exclusion_intervals_file {
        intervals = remove_intervals_near_gaps(
            intervals,
            arguments.exclusion_interval_padding,
            path,
            metadata,
        )?;
    }
    intervals.retain(|interval| !contained_in_region_to_ignore(interval, ignored));

    let mapped = mapped_qnames(records, &intervals, ignored);
    let mut qnames: Hopscotch<QNameAndInterval> = Hopscotch::new(
        arguments
            .assembly_to_mapped_size_ratio_guess
            .wrapping_mul(mapped.len() as i32),
    );
    for entry in mapped {
        qnames.add(entry)?;
    }
    if let Some(path) = &arguments.qnames_mapped_file {
        write_qnames(path, &qnames)?;
    }
    if intervals.is_empty() {
        return Ok(Results {
            finished: false,
            intervals,
            assemblies: Vec::new(),
            links,
        });
    }

    let n_intervals = intervals.len();
    let mut assemblies: Vec<(i32, Disposition)> = Vec::new();
    if !arguments.interval_only_assembly {
        let kills = kmer_kill_set(arguments)?;
        let kmers = kmer_intervals(
            records,
            arguments,
            metadata,
            &mut qnames,
            n_intervals,
            &kills,
            &mut assemblies,
        )?;
        let mut table: Hopscotch<KmerAndInterval> = Hopscotch::new(kmers.len() as i32);
        for entry in kmers {
            table.add(entry)?;
        }
        remove_ubiquitous_kmers(records, arguments, metadata, &mut table)?;
        for entry in assembly_qnames(records, arguments.k_size, &table)? {
            qnames.add(entry)?;
        }
        if let Some(path) = &arguments.qnames_assembly_file {
            write_qnames(path, &qnames)?;
        }
    }

    // `ReadsForQNamesFinder`, every interval's reads sharing one `FastqRead`.
    let mut interval_reads: Vec<Vec<std::rc::Rc<FastqRead>>> = Vec::new();
    interval_reads.resize_with(n_intervals, Vec::new);
    for record in records {
        if !is_primary_line(record) {
            continue;
        }
        let name = &record.read_name;
        let found = qnames.find_each(string_hash_code(name), |entry| entry.name == *name);
        if found.is_empty() {
            continue;
        }
        let read = std::rc::Rc::new(FastqRead::new(
            record,
            metadata,
            arguments.include_mapping_location,
        ));
        for entry in found {
            interval_reads[entry.interval as usize].push(read.clone());
        }
    }

    // The handler is a map over the populated intervals, whichever partition each landed in, and
    // the driver sorts what comes back.
    let mut aligner = None;
    for (id, reads) in interval_reads.iter().enumerate() {
        if reads.is_empty() {
            continue;
        }
        let disposition = assemble_interval(id as i32, reads, arguments, &mut aligner)?;
        assemblies.push((id as i32, disposition));
    }
    assemblies.sort_by_key(|(id, _)| *id);
    if let Some(path) = &arguments.interval_file {
        write_interval_file(path, &intervals, &assemblies, header)?;
    }
    Ok(Results {
        finished: true,
        intervals,
        assemblies,
        links,
    })
}

/// `AlignedAssemblyOrExcuse.toSAMStreamForAlignmentsOfThisAssembly` for every assembly, in assembly
/// order, which runs on the driver once everything else is written. `ref_names` is the BAM's
/// dictionary, and the JNI's reference id is an index into it whichever image it came from.
pub fn contig_records(
    assemblies: &[(i32, Disposition)],
    ref_names: &[String],
) -> Result<Vec<BamRecord>, Failure> {
    let mut records = Vec::new();
    for (id, disposition) in assemblies {
        let Disposition::Aligned {
            contigs,
            alignments,
            ..
        } = disposition
        else {
            continue;
        };
        for (index, lines) in alignments.iter().enumerate() {
            let name = format!("asm{id:06}:tig{index:05}");
            let tags = sa_tags(lines, ref_names)?;
            for (line, sa) in lines.iter().zip(tags) {
                records.push(apply_alignment(&name, &contigs[index], line, sa, ref_names)?);
            }
        }
    }
    Ok(records)
}

fn ref_name(ref_names: &[String], index: i32) -> Result<&String, Failure> {
    ref_names.get(index as usize).ok_or_else(|| Failure::Thrown {
        class: "java.lang.IndexOutOfBoundsException",
        message: format!("Index {index} out of bounds for length {}", ref_names.len()),
    })
}

/// `BwaMemAlignmentUtils.createSATags`.
fn sa_tags(
    lines: &[ContigAlignment],
    ref_names: &[String],
) -> Result<Vec<Option<String>>, Failure> {
    if lines.len() < 2 {
        return Ok(vec![None; lines.len()]);
    }
    let mut selves = Vec::with_capacity(lines.len());
    for line in lines {
        selves.push(format!(
            "{},{},{},{},{},{};",
            ref_name(ref_names, line.ref_id)?,
            line.start + 1,
            if line.reverse { "-" } else { "+" },
            line.cigar.replace('H', "S"),
            line.mapq,
            line.mismatches
        ));
    }
    Ok((0..lines.len())
        .map(|index| {
            let others: Vec<&str> = selves
                .iter()
                .enumerate()
                .filter(|(other, _)| *other != index)
                .map(|(_, tag)| tag.as_str())
                .collect();
            Some(others.concat())
        })
        .collect())
}

/// `BwaMemAlignmentUtils.applyAlignment` with no qualities, not soft-clipping, and the tags of
/// `toSAMStreamForRead` after it.
fn apply_alignment(
    name: &str,
    sequence: &[u8],
    line: &ContigAlignment,
    sa: Option<String>,
    ref_names: &[String],
) -> Result<BamRecord, Failure> {
    let mut record = BamRecord {
        read_name: name.to_string(),
        ..BamRecord::default()
    };
    let mut bases = sequence.to_vec();
    if line.ref_id < 0 {
        record.flags = 0x4;
    } else {
        ref_name(ref_names, line.ref_id)?;
        record.reference_index = line.ref_id;
        record.alignment_start = line.start as i32 + 1;
        record.mapping_quality = line.mapq as u8;
        if line.reverse {
            record.flags |= REVERSE;
            for base in bases.iter_mut() {
                *base = simple_complement(*base);
            }
            bases.reverse();
        }
        let mut elements = crate::pathseq_score::text_cigar(&line.cigar)
            .map_err(|error| Failure::Port(format!("{}: {error:?}", line.cigar)))?;
        if line.supplementary {
            record.flags |= 0x800;
            let clip = |element: Option<&CigarElement>| {
                element.filter(|one| one.op == CigarOp::S).map_or(0, |one| one.length as usize)
            };
            let (from, to) = (clip(elements.first()), bases.len() - clip(elements.last()));
            bases = bases[from..to].to_vec();
            for element in &mut elements {
                if element.op == CigarOp::S {
                    element.op = CigarOp::H;
                }
            }
        }
        record.cigar = Cigar::new(elements);
        let int = |value: i32| TagValue::Int(i64::from(value));
        record.tags.insert(Tag::new(b"NM"), int(line.mismatches));
        record.tags.insert(Tag::new(b"AS"), int(line.score));
        record.tags.insert(Tag::new(b"XS"), int(line.suboptimal));
        record.tags.insert(Tag::new(b"MD"), TagValue::Str(line.md.clone()));
        if let Some(xa) = &line.xa {
            record.tags.insert(Tag::new(b"XA"), TagValue::Str(xa.clone()));
        }
    }
    if let Some(sa) = sa {
        record.tags.insert(Tag::new(b"SA"), TagValue::Str(sa));
    }
    record
        .tags
        .insert(Tag::new(b"RG"), TagValue::Str(CONTIG_READ_GROUP.to_string()));
    record.read_bases = bases;
    Ok(record)
}

/// `SVInterval.midpoint`.
fn midpoint(interval: &SvInterval) -> i32 {
    interval.start.wrapping_add(interval.end) / 2
}

/// `SVContext.composePaddedInterval`, as the half-open `SVInterval` `getPairedStrandedIntervals`
/// makes of it.
fn padded(contig: i32, length: i32, start: i32, end: i32, padding: i32) -> SvInterval {
    SvInterval {
        contig,
        start: (if padding > 0 { start - padding + 1 } else { start }).max(1),
        end: end.saturating_add(padding).min(length) + 1,
    }
}

/// `processEvidenceTargetLinks`: the tree the links make, each deletion taking the links that
/// overlap it, then an imprecise deletion from each link left over with enough to it.
pub fn annotate_with_links(
    variants: Vec<VariantContext>,
    links: &[TargetLink],
    metadata: &ReadMetadata,
    header: &SamHeader,
    reference: &mut dyn Reference,
    threshold: i32,
    max_deletion: i32,
) -> Result<Vec<VariantContext>, Failure> {
    let mut tree: PairedTree<TargetLink> = PairedTree::new();
    for link in links {
        tree.put(link.source, link.target, link.clone());
    }
    let int = |value: i32| Value::Int(i64::from(value));
    let padding = metadata
        .libraries
        .iter()
        .map(|(_, statistics)| statistics.median)
        .max()
        .unwrap_or(0)
        / 2;
    let mut annotated = Vec::with_capacity(variants.len());
    for mut variant in variants {
        let is_deletion = variant
            .attributes
            .iter()
            .any(|(name, value)| name == "SVTYPE" && *value == Value::Str("DEL".to_string()));
        if is_deletion {
            let contig = header
                .sequences
                .iter()
                .position(|sequence| sequence.name == variant.contig)
                .expect("a variant on the dictionary");
            let length = header.sequences[contig].length;
            let homology = variant
                .attributes
                .iter()
                .find_map(|(name, value)| match value {
                    Value::Int(homology) if name == "HOMLEN" => Some(*homology as i32),
                    _ => None,
                })
                .unwrap_or(0);
            let (start, end) = (variant.start as i32, variant.stop as i32);
            let left = padded(contig as i32, length, start + 1, start + 1 + homology, padding);
            let right = padded(contig as i32, length, end, end + homology, padding);
            let (mut read_pairs, mut split_reads) = (0i32, 0i32);
            for key in tree.overlappers(&(left, true), &(right, false)) {
                if let Some((_, link)) = tree.remove(&key) {
                    read_pairs = read_pairs.wrapping_add(link.read_pairs);
                    split_reads = split_reads.wrapping_add(link.split_reads);
                }
            }
            if read_pairs > 0 {
                variant.attributes.push(("READ_PAIR_SUPPORT".to_string(), int(read_pairs)));
            }
            if split_reads > 0 {
                variant.attributes.push(("SPLIT_READ_SUPPORT".to_string(), int(split_reads)));
            }
        }
        annotated.push(variant);
    }

    // `callImpreciseDeletionFromEvidenceLinks`, over what the deletions left.
    for (_, _, link) in tree.entries() {
        let (left, right) = (link.source, link.target);
        let same_contig = left.0.contig == right.0.contig;
        let distance = if same_contig {
            right.0.start.wrapping_sub(left.0.end)
        } else {
            -1
        };
        if !(same_contig && left.1 && !right.1)
            || distance >= max_deletion
            || link.read_pairs.wrapping_add(link.split_reads) <= threshold
        {
            continue;
        }
        let contig = &metadata.contigs[left.0.contig as usize];
        let (start, end) = (midpoint(&left.0), midpoint(&right.0));
        let bases = reference
            .bases(contig, start, start)
            .map_err(|message| Failure::Task(message))?;
        let reference_allele = Allele::create(&bases, true)
            .map_err(|error| Failure::Port(format!("{error:?}")))?;
        let alternate =
            Allele::from_str("<DEL>", false).map_err(|error| Failure::Port(format!("{error:?}")))?;
        let mut variant = VariantContext::new(
            contig,
            i64::from(start),
            vec![reference_allele, alternate],
        );
        variant.stop = i64::from(end);
        variant.id = format!(
            "DEL_IMPRECISE_{contig}_{}_{}_{}_{}",
            left.0.start, left.0.end, right.0.start, right.0.end
        );
        for (name, value) in [
            ("END", int(end)),
            ("SVTYPE", Value::Str("DEL".to_string())),
            ("SVLEN", int(start.wrapping_sub(end))),
            ("IMPRECISE", Value::Bool(true)),
            (
                "CIPOS",
                Value::Str(format!("{},{}", left.0.start - start, left.0.end - start)),
            ),
            (
                "CIEND",
                Value::Str(format!("{},{}", right.0.start - end, right.0.end - end)),
            ),
            ("READ_PAIR_SUPPORT", int(link.read_pairs)),
            ("SPLIT_READ_SUPPORT", int(link.split_reads)),
        ] {
            variant.attributes.push((name.to_string(), value));
        }
        annotated.push(variant);
    }
    Ok(annotated)
}

/// `SVUtils.GATKSV_CONTIG_ALIGNMENTS_READ_GROUP_ID`.
pub const CONTIG_READ_GROUP: &str = "GATKSVContigAlignments";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_partitioned_key_hashes_as_scala_2_13_does() {
        // `new scala.Tuple2<>(0, key).hashCode()` in the oracle's Scala 2.13.14.
        for (key, expected) in [
            (0, -737015216),
            (1, 1084890308),
            (-1, 1475012417),
            (5, 179566373),
            (123456789, 1207535907),
            (i32::MIN, -1466823222),
            (i32::MAX, -249071175),
        ] {
            assert_eq!(partitioned_key_hash(key), expected);
        }
    }

    #[test]
    fn a_shuffle_leaves_the_hash_order_of_its_slots() {
        let pairs: Vec<(i32, i32)> = (0..200).map(|key| (key % 40, 1)).collect();
        let reduced = reduce_by_key(pairs, |key| *key, |one, two| one + two);
        assert_eq!(reduced.len(), 40);
        assert!(reduced.iter().all(|(_, count)| *count == 5));
    }

    #[test]
    fn an_overlap_search_takes_one_right_end_per_left_end() {
        let at = |start, end| SvInterval {
            contig: 0,
            start,
            end,
        };
        let mut tree = PairedTree::new();
        assert!(tree.put((at(10, 20), true), (at(100, 110), false), 1));
        assert!(!tree.put((at(10, 20), true), (at(100, 110), false), 2));
        assert!(tree.put((at(10, 20), true), (at(105, 115), false), 3));
        assert!(tree.put((at(15, 25), true), (at(100, 110), false), 4));
        let found = tree.overlappers(&(at(12, 18), true), &(at(102, 108), false));
        assert_eq!(
            found,
            vec![(at(10, 20), at(100, 110)), (at(15, 25), at(100, 110))]
        );
        tree.remove(&found[0]);
        assert_eq!(tree.entries().len(), 2);
    }

    fn graph_of(contigs: &[&str], connections: &[(usize, usize, i32, bool, bool)]) -> Graph {
        let mut graph = Graph {
            contigs: contigs
                .iter()
                .map(|text| Contig {
                    sequence: text.as_bytes().to_vec(),
                    supporting_reads: 1,
                    connections: Vec::new(),
                })
                .collect(),
            next_connection: 0,
        };
        for (from, to, overlap, is_rc, is_target_rc) in connections {
            let connection = graph.connection(*to, *overlap, *is_rc, *is_target_rc);
            graph.contigs[*from].connections.push(connection);
        }
        graph
    }

    #[test]
    fn two_contigs_with_nothing_else_to_join_become_one() {
        let mut graph = graph_of(
            &["AACCGG", "GGTTAA"],
            &[(0, 1, 2, false, false), (1, 0, 2, true, false)],
        );
        let list = graph.remove_unbranched(vec![0, 1]).expect("joined");
        assert_eq!(list.len(), 1);
        assert_eq!(graph.contigs[list[0]].sequence, b"AACCGGTTAA");
        assert_eq!(graph.contigs[list[0]].supporting_reads, 2);
    }

    #[test]
    fn a_bubble_expands_into_both_of_its_paths() {
        let mut graph = graph_of(
            &["AAAC", "ACGG", "ACTT", "GGAA"],
            &[
                (0, 1, 2, false, false),
                (0, 2, 2, false, false),
                (1, 0, 2, true, false),
                (2, 0, 2, true, false),
            ],
        );
        let list = graph.expand(vec![0, 1, 2, 3]).expect("expanded");
        let sequences: Vec<&[u8]> = list
            .iter()
            .map(|contig| graph.contigs[*contig].sequence.as_slice())
            .collect();
        assert_eq!(sequences, vec![&b"AAACGG"[..], &b"AAACTT"[..], &b"GGAA"[..]]);
    }

    #[test]
    fn lines_and_tabs_read_as_java_reads_them() {
        assert_eq!(java_lines("a\nb\r\nc"), vec!["a", "b", "c"]);
        assert_eq!(java_lines("a\n"), vec!["a"]);
        assert_eq!(java_split_tabs("a\tb\t\t"), vec!["a", "b"]);
    }
}
