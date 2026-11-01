//! `PSFilter` and its transformers (GATK 4.6.2.0): the first step of PathSeq, without its BWA step.
//!
//! # One partition's order is three tables
//!
//! The file order is not a sort. `filterDuplicateSequences` is a `groupByKey`, whose reducer is
//! Spark's `AppendOnlyMap`, walked in slot order, the slot being the murmur3 of the key's
//! `hashCode`, probed quadratically and doubled at seven tenths. `getPairedAndUnpairedLists` then
//! walks the reads: a read whose name is waiting joins the paired list ahead of its mate, and those
//! still waiting come out in `HashMap` order, a table sized by `--filter-reads-per-partition`.
//!
//! # A read is rebuilt twice, and only the second time is what is written
//!
//! `clearReadAlignment` keeps the name, bases, qualities, pairing and read group, turning a
//! reverse-strand read around, and drops every other tag. It runs again after the pairing flags are
//! set, so a read that lost its mate loses the mate flags too.

use std::collections::{HashSet, VecDeque};

use gatk_engine::java_hash::JavaHashMap;
use gatk_engine::java_random::JavaRandom;
use gatk_engine::read;
use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sequence::complement;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};

use crate::pathseq_kmers::{kmerize, Kmer};
use crate::sv_kmer::fnv_long_64;

/// `ADAPTER_SEQUENCES`, each `IlluminaAdapterPair`'s five and three prime in the reference's order.
const ADAPTERS: [&str; 6] = [
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
];

/// `PSFilterArgumentCollection`, without the BWA arguments.
#[derive(Debug, Clone)]
pub struct FilterArguments {
    pub aligned_input: bool,
    pub contigs_to_ignore: Vec<String>,
    pub min_identity: i32,
    pub skip_filters: bool,
    pub max_adapter_mismatches: i32,
    pub min_adapter_length: i32,
    pub dust_mask: i32,
    pub dust_window: i32,
    pub dust_t: f64,
    pub read_trim_threshold: i32,
    pub min_read_length: i32,
    pub quality_threshold: i32,
    pub max_ambiguous_bases: i32,
    pub host_kmer_threshold: i32,
    pub filter_duplicates: bool,
    pub reads_per_partition: i32,
}

/// What the tool hands the writer and the logger.
pub struct Filtered {
    pub paired: Vec<BamRecord>,
    pub unpaired: Vec<BamRecord>,
    pub metrics: FilterMetrics,
}

/// `PSFilterMetrics`.
#[derive(Debug, Clone, Default)]
pub struct FilterMetrics {
    pub primary: i64,
    pub after_prealigned_host: i64,
    pub after_quality: i64,
    pub after_host: i64,
    pub after_deduplication: i64,
    pub final_paired: i64,
}

impl MetricBean for FilterMetrics {
    fn class_name(&self) -> &str {
        "org.broadinstitute.hellbender.tools.spark.pathseq.loggers.PSFilterMetrics"
    }

    fn columns(&self) -> &[&'static str] {
        &[
            "PRIMARY_READS",
            "READS_AFTER_PREALIGNED_HOST_FILTER",
            "READS_AFTER_QUALITY_AND_COMPLEXITY_FILTER",
            "READS_AFTER_HOST_FILTER",
            "READS_AFTER_DEDUPLICATION",
            "FINAL_PAIRED_READS",
            "FINAL_UNPAIRED_READS",
            "FINAL_TOTAL_READS",
            "LOW_QUALITY_OR_LOW_COMPLEXITY_READS_FILTERED",
            "HOST_READS_FILTERED",
            "DUPLICATE_READS_FILTERED",
        ]
    }

    fn values(&self) -> Vec<Value> {
        // `computeDerivedMetrics`.
        [
            self.primary,
            self.after_prealigned_host,
            self.after_quality,
            self.after_host,
            self.after_deduplication,
            self.final_paired,
            self.after_deduplication - self.final_paired,
            self.after_deduplication,
            self.after_prealigned_host - self.after_quality,
            self.primary - self.after_prealigned_host + self.after_quality - self.after_host,
            self.after_host - self.after_deduplication,
        ]
        .into_iter()
        .map(Value::Long)
        .collect()
    }
}

/// `PSFilterFileLogger.close`'s file.
pub fn metrics_file(metrics: &FilterMetrics, headers: &[String]) -> String {
    let mut file = MetricsFile::default();
    for line in headers {
        file.add_header(line);
    }
    file.add_metric(metrics);
    file.write()
}

/// `PSKmerCollection`, read: membership is all the filter asks of it.
pub enum KmerLibrary {
    Set {
        kmer_size: usize,
        mask: Kmer,
        values: HashSet<i64>,
    },
    Bloom {
        kmer_size: usize,
        mask: Kmer,
        filter: LongBloomFilter,
    },
}

/// `LongBloomFilter`, read.
pub struct LongBloomFilter {
    total_bits: i64,
    num_hashes: i32,
    bucket_array_size: i32,
    buckets: Vec<Vec<u8>>,
}

impl LongBloomFilter {
    fn contains(&self, key: i64) -> bool {
        const HASH_SEED_2: i64 = 0x6cebe6dca7f118a6;
        let hash1 = fnv_long_64(FNV64_DEFAULT_SEED, key);
        let hash2 = fnv_long_64(HASH_SEED_2, key);
        for i in 0..self.num_hashes {
            let result = hash1.wrapping_add((i as i64).wrapping_mul(hash2)) % self.total_bits;
            let bit = if result < 0 {
                result + self.total_bits
            } else {
                result
            };
            let byte = (bit as u64 >> 3) as i64;
            let array = (byte / self.bucket_array_size as i64) as usize;
            let index = (byte % self.bucket_array_size as i64) as usize;
            if self.buckets[array][index] & (1u8 << (bit & 7)) == 0 {
                return false;
            }
        }
        true
    }
}

/// `SVUtils.FNV64_DEFAULT_SEED`.
const FNV64_DEFAULT_SEED: i64 = 14695981039346656037u64 as i64;

/// Kryo's fixed-width reads and its reference marker, over the bytes of a k-mer file.
struct Input<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Input<'_> {
    fn take(&mut self, count: usize) -> Result<&[u8], String> {
        let end = self
            .at
            .checked_add(count)
            .filter(|end| *end <= self.bytes.len());
        let end = end.ok_or_else(|| "Buffer underflow.".to_string())?;
        let taken = &self.bytes[self.at..end];
        self.at = end;
        Ok(taken)
    }

    fn int(&mut self) -> Result<i32, String> {
        let bytes = self.take(4)?;
        Ok(i32::from_be_bytes(bytes.try_into().expect("four bytes")))
    }

    fn long(&mut self) -> Result<i64, String> {
        let bytes = self.take(8)?;
        Ok(i64::from_be_bytes(bytes.try_into().expect("eight bytes")))
    }

    fn marker(&mut self) -> Result<(), String> {
        self.take(1).map(|_| ())
    }
}

/// `PSKmerUtils.readKmerFilter`: the extension chooses the class.
pub fn read_kmer_library(path: &str, bytes: &[u8]) -> Result<KmerLibrary, String> {
    let mut input = Input { bytes, at: 0 };
    if path.ends_with(".hss") {
        input.marker()?;
        let kmer_size = input.int()? as usize;
        let mask = Kmer(input.long()? as u64);
        input.marker()?;
        let mut values = HashSet::new();
        for _ in 0..input.int()? {
            let _capacity = input.int()?;
            for _ in 0..input.int()? {
                values.insert(input.long()?);
            }
        }
        return Ok(KmerLibrary::Set {
            kmer_size,
            mask,
            values,
        });
    }
    input.marker()?;
    let kmer_size = input.int()? as usize;
    let mask = Kmer(input.long()? as u64);
    input.marker()?;
    let total_bits = input.long()?;
    let _total_buckets = input.long()?;
    let num_bucket_arrays = input.int()?;
    let bucket_array_size = input.int()?;
    let final_bucket_array_size = input.int()?;
    let num_hashes = input.int()?;
    let mut buckets = Vec::new();
    for i in 0..num_bucket_arrays {
        let size = if i + 1 == num_bucket_arrays {
            final_bucket_array_size
        } else {
            bucket_array_size
        };
        buckets.push(input.take(size as usize)?.to_vec());
    }
    Ok(KmerLibrary::Bloom {
        kmer_size,
        mask,
        filter: LongBloomFilter {
            total_bits,
            num_hashes,
            bucket_array_size,
            buckets,
        },
    })
}

impl KmerLibrary {
    fn kmer_size(&self) -> usize {
        match self {
            KmerLibrary::Set { kmer_size, .. } | KmerLibrary::Bloom { kmer_size, .. } => *kmer_size,
        }
    }

    fn contains(&self, kmer: Kmer) -> bool {
        let size = self.kmer_size();
        // `canonicalizeAndMask`. A file holds an odd size, which is the only one `canonical` takes.
        let canonical = |mask: Kmer| kmer.canonical(size).map(|kmer| kmer.mask(mask).0 as i64);
        match self {
            KmerLibrary::Set { mask, values, .. } => {
                canonical(*mask).is_ok_and(|key| values.contains(&key))
            }
            KmerLibrary::Bloom { mask, filter, .. } => {
                canonical(*mask).is_ok_and(|key| filter.contains(key))
            }
        }
    }
}

/// `PSUtils.getMatchesLessDeletions`.
pub fn matches_less_deletions_of(elements: &[CigarElement], nm: i32) -> i32 {
    let (mut aligned, mut deletions, mut insertions) = (0i32, 0i32, 0i32);
    for element in elements {
        match element.op {
            Op::M | Op::Eq | Op::X => aligned += element.length as i32,
            Op::D => deletions += element.length as i32,
            Op::I => insertions += element.length as i32,
            _ => {}
        }
    }
    let substitutions = nm - deletions - insertions;
    let matches = aligned - substitutions;
    if substitutions < 0 || matches < 0 {
        return 0;
    }
    matches - deletions
}

/// `HostAlignmentReadFilter.test`.
fn survives_host_alignment(
    record: &BamRecord,
    contigs: &[String],
    ignored: &[String],
    min_identity: i32,
) -> bool {
    if read::is_unmapped(record) {
        return true;
    }
    let contig = usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| contigs.get(index));
    if contig.is_some_and(|contig| ignored.contains(contig)) {
        return true;
    }
    match record.tags.get(Tag::new(b"NM")) {
        Some(TagValue::Int(nm)) => {
            matches_less_deletions_of(&record.cigar.elements, *nm as i32) < min_identity
        }
        _ => true,
    }
}

/// `SequenceUtil.reverseComplement`, in place over a copy.
fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|base| complement(*base)).collect()
}

/// `clearReadAlignment`.
fn cleared(record: &BamRecord) -> BamRecord {
    // A reverse read stays as it is: the reference reverses `getBases()`, which is a copy.
    let mut out = BamRecord {
        read_name: record.read_name.clone(),
        read_bases: record.read_bases.clone(),
        base_qualities: record.base_qualities.clone(),
        ..BamRecord::default()
    };
    let mut flags = 0x4u16;
    if read::is_paired(record) {
        flags |= 0x1 | 0x8;
        if read::is_first_of_pair(record) {
            flags |= 0x40;
        } else if read::is_second_of_pair(record) {
            flags |= 0x80;
        }
    }
    out.flags = flags;
    if let Some(group) = record.tags.get(Tag::new(b"RG")) {
        out.tags.insert(Tag::new(b"RG"), group.clone());
    }
    out
}

/// `StripMateNumberTransformer`: `String.trim`, a mate suffix, and `trim` again.
fn strip_mate_number(record: &mut BamRecord) {
    let is_trimmed = |c: char| c <= ' ';
    let name = record.read_name.trim_matches(is_trimmed);
    if name.ends_with("/1") || name.ends_with("/2") {
        let stripped = name[..name.len() - 2].trim_matches(is_trimmed).to_string();
        record.read_name = stripped;
    }
}

/// `ReadClipper` with `HARDCLIP_BASES`, over a read with no alignment: `stop` is inclusive, and past
/// the end it means the end. An `emptyRead` is the same flags and the read group.
fn hard_clip(record: &mut BamRecord, start: usize, stop: usize) {
    let length = record.read_bases.len();
    if start >= length {
        return;
    }
    let stop = stop.min(length - 1);
    let keep = |values: &mut Vec<u8>| {
        if values.len() == length {
            if start == 0 {
                values.drain(..=stop);
            } else {
                values.truncate(start);
            }
        }
    };
    keep(&mut record.read_bases);
    keep(&mut record.base_qualities);
}

/// `AdapterTrimTransformer`.
fn trim_adapters(record: &mut BamRecord, max_mismatches: i32, min_clip_length: i32) {
    let bases = &record.read_bases;
    if (bases.len() as i32) < min_clip_length {
        return;
    }
    let mut best_mismatches = max_mismatches;
    let mut best_start = bases.len() as i32;
    let mut found = false;
    for adapter in ADAPTERS {
        let adapter = adapter.as_bytes();
        let mut alignment_length = min_clip_length;
        let mut start = bases.len() as i32 - min_clip_length;
        while start >= 0 {
            let mut mismatches = 0;
            for j in 0..alignment_length {
                // `charAt` past the adapter, for a --min-adapter-length longer than it, throws.
                let wanted = adapter.get(j as usize).copied().unwrap_or(b'N');
                let is_no_call = matches!(wanted, b'N' | b'n' | b'.');
                if !is_no_call && bases[(start + j) as usize] != wanted {
                    mismatches += 1;
                    if mismatches > max_mismatches {
                        break;
                    }
                }
            }
            if mismatches < best_mismatches || (mismatches == best_mismatches && start < best_start)
            {
                best_mismatches = mismatches;
                best_start = start;
                found = true;
            }
            if (alignment_length as usize) < adapter.len() {
                alignment_length += 1;
            }
            start -= 1;
        }
    }
    if found {
        let length = record.read_bases.len();
        hard_clip(record, best_start as usize, length);
    }
}

/// `SimpleRepeatMaskTransformer`.
fn mask_simple_repeats(record: &mut BamRecord, thresh_at: i32, thresh_gc: i32, window: usize) {
    let original = record.read_bases.clone();
    let mut counts = [0i8; 256];
    let is_repeat = |counts: &[i8; 256]| {
        let count = |base: u8| counts[base as usize] as i32;
        let n = count(b'N') + count(b'n');
        let at = count(b'A') + count(b'a') + count(b'T') + count(b't') + n;
        let gc = count(b'G') + count(b'g') + count(b'C') + count(b'c') + n;
        at >= thresh_at || gc >= thresh_gc
    };
    let mut window_start = 0;
    let mut window_end = window.min(original.len());
    for base in &original[..window_end] {
        counts[*base as usize] = counts[*base as usize].wrapping_add(1);
    }
    let (mut mask_start, mut mask_end) = (0, 0);
    if is_repeat(&counts) {
        mask_end = window_end;
    }
    while window_end < original.len() {
        let leaving = original[window_start] as usize;
        counts[leaving] = counts[leaving].wrapping_sub(1);
        window_start += 1;
        let entering = original[window_end] as usize;
        counts[entering] = counts[entering].wrapping_add(1);
        window_end += 1;
        if is_repeat(&counts) {
            if window_end == mask_end + 1 {
                mask_end += 1;
            } else {
                mask(record, mask_start, mask_end);
                mask_start = window_start;
                mask_end = window_end;
            }
        }
    }
    mask(record, mask_start, mask_end);
}

/// `WRITE_NS_Q0S` over `start .. end`.
fn mask(record: &mut BamRecord, start: usize, end: usize) {
    for at in start..end.min(record.read_bases.len()) {
        record.read_bases[at] = b'N';
    }
    for at in start..end.min(record.base_qualities.len()) {
        record.base_qualities[at] = 0;
    }
}

#[derive(Clone, Copy)]
struct Perfect {
    start: i32,
    finish: i32,
    score: f64,
}

fn save_masked_regions(perfect: &mut Vec<Perfect>, intervals: &mut Vec<(i32, i32)>, wstart: i32) {
    let Some(last) = perfect.last().copied() else {
        return;
    };
    if last.start >= wstart {
        return;
    }
    match intervals.last_mut() {
        Some(previous) if last.start < previous.1 + 1 => {
            previous.1 = previous.1.max(last.finish);
        }
        _ => intervals.push((last.start, last.finish)),
    }
    while perfect.last().is_some_and(|p| p.start < wstart) {
        perfect.pop();
    }
}

/// `DUSTReadTransformer.runDUST`, the sDUST pseudocode as the reference wrote it. `None` is the
/// `NoSuchElementException` a window of two or less pops from an empty list.
fn run_dust(bases: &[u8], window: i32, tscore: f64) -> Option<Vec<(i32, i32)>> {
    let mut rng: Option<JavaRandom> = None;
    let q: Vec<i32> = bases
        .iter()
        .map(|base| match base {
            b'A' | b'a' => 0,
            b'T' | b't' => 1,
            b'C' | b'c' => 2,
            b'G' | b'g' => 3,
            _ => rng
                .get_or_insert_with(|| JavaRandom::new(0))
                .next_int_bound(4),
        })
        .collect();
    let mut intervals = Vec::new();
    let mut perfect: Vec<Perfect> = Vec::new();
    let mut count_suffix = [0i16; 64];
    let mut count_window = [0i16; 64];
    let (mut score_suffix, mut score_window) = (0i32, 0i32);
    let mut triplets: VecDeque<usize> = VecDeque::new();
    let mut l = 0i32;
    let length = q.len() as i32;

    let mut wfinish = 2;
    while wfinish < length {
        let wstart = (wfinish - window + 1).max(0);
        save_masked_regions(&mut perfect, &mut intervals, wstart);
        let at = wfinish as usize;
        let mut t = (q[at - 2] | (q[at - 1] << 2) | (q[at] << 4)) as usize;
        if triplets.len() as i32 >= window - 2 {
            let s = triplets.pop_front()?;
            count_window[s] -= 1;
            score_window -= count_window[s] as i32;
            if l > triplets.len() as i32 {
                l -= 1;
                count_suffix[s] -= 1;
                score_suffix -= count_suffix[s] as i32;
            }
        }
        triplets.push_back(t);
        l += 1;
        score_window += count_window[t] as i32;
        count_window[t] += 1;
        score_suffix += count_suffix[t] as i32;
        count_suffix[t] += 1;
        if count_suffix[t] as f64 * 10.0 > 2.0 * tscore {
            loop {
                let s = triplets[triplets.len() - l as usize];
                count_suffix[s] -= 1;
                score_suffix -= count_suffix[s] as i32;
                l -= 1;
                if s == t {
                    break;
                }
            }
        }
        if (score_window * 10) as f64 > l as f64 * tscore {
            let mut c = count_suffix;
            let mut r = score_suffix;
            // A `ListIterator`: the cursor is the index the next `next()` answers.
            let mut cursor = 0usize;
            let advance = |cursor: &mut usize, perfect: &Vec<Perfect>| {
                let next = perfect.get(*cursor).copied();
                if next.is_some() {
                    *cursor += 1;
                }
                next
            };
            let mut perf = advance(&mut cursor, &perfect);
            let mut max_score = 0.0f64;
            let mut i = triplets.len() as i32 - l - 1;
            while i >= 0 {
                t = triplets[i as usize];
                r += c[t] as i32;
                c[t] += 1;
                let new_score = r as f64 / (triplets.len() as i32 - i - 1) as f64;
                if new_score * 10.0 > tscore {
                    while let Some(current) = perf {
                        if current.start < i + wstart {
                            break;
                        }
                        max_score = max_score.max(current.score);
                        perf = advance(&mut cursor, &perfect);
                    }
                    if new_score >= max_score {
                        max_score = new_score;
                        if cursor > 0 {
                            cursor -= 1;
                        }
                        perfect.insert(
                            cursor,
                            Perfect {
                                start: i + wstart,
                                finish: triplets.len() as i32 + 1 + wstart,
                                score: new_score,
                            },
                        );
                        cursor += 1;
                        if cursor < perfect.len() {
                            cursor += 1;
                        }
                    }
                }
                i -= 1;
            }
        }
        wfinish += 1;
    }
    let mut wstart = (length - window + 1).max(0);
    while !perfect.is_empty() {
        save_masked_regions(&mut perfect, &mut intervals, wstart);
        wstart += 1;
    }
    Some(intervals)
}

/// `DUSTReadTransformer.apply`.
fn dust(record: &mut BamRecord, mask_quality: i32, window: i32, tscore: f64) -> Result<(), String> {
    let intervals = run_dust(&record.read_bases, window, tscore)
        .ok_or_else(|| "java.util.NoSuchElementException".to_string())?;
    let with_qualities = record.base_qualities.len() == record.read_bases.len();
    for (start, end) in intervals {
        for at in start..=end {
            record.read_bases[at as usize] = b'N';
            if with_qualities {
                record.base_qualities[at as usize] = mask_quality as u8;
            }
        }
    }
    Ok(())
}

/// `BaseQualityClipReadTransformer`.
fn clip_low_quality(record: &mut BamRecord, threshold: i32) {
    let point = |qualities: &mut dyn Iterator<Item = (usize, u8)>| {
        let (mut sum, mut last_max, mut clip) = (0i32, -1i32, None);
        for (at, quality) in qualities {
            sum += threshold - quality as i8 as i32;
            if sum >= 0 && sum >= last_max {
                last_max = sum;
                clip = Some(at);
            }
        }
        clip
    };
    let qualities = record.base_qualities.clone();
    let right = point(&mut qualities.iter().copied().enumerate().rev());
    if let Some(right) = right {
        let length = record.read_bases.len();
        hard_clip(record, right, length);
    }
    let qualities = record.base_qualities.clone();
    if let Some(left) = point(&mut qualities.iter().copied().enumerate()) {
        hard_clip(record, 0, left);
    }
}

/// `SVUtils.fnvByteArray64(byte[])`, whose seed the reference notes is the wrong constant.
fn fnv_byte_array_64(bytes: &[u8]) -> i64 {
    let mut start = 1099511628211i64;
    for chunk in bytes.chunks(8) {
        let mut value = 0i64;
        for byte in chunk {
            value = (value << 8) | *byte as i8 as i64;
        }
        start = fnv_long_64(start, value);
    }
    start
}

/// `Hashing.murmur3_32().hashInt`, which is what `AppendOnlyMap.rehash` calls.
fn murmur3_32_int(input: i32) -> i32 {
    let mut k = (input as u32).wrapping_mul(0xcc9e2d51);
    k = k.rotate_left(15).wrapping_mul(0x1b873593);
    let mut h = k;
    h = h.rotate_left(13).wrapping_mul(5).wrapping_add(0xe6546b64);
    h ^= 4;
    h ^= h >> 16;
    h = h.wrapping_mul(0x85ebca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2ae35);
    h ^= h >> 16;
    h as i32
}

/// `AppendOnlyMap` keyed on a `Long`, holding the first value and whatever `first` keeps, walked in
/// slot order as `ExternalAppendOnlyMap.iterator` walks it when nothing has spilled.
fn group_by_key<T>(pairs: Vec<(i64, T)>) -> Vec<Vec<T>> {
    let slot_hash =
        |key: i64| murmur3_32_int((key ^ ((key as u64) >> 32) as i64) as i32) as u32 as usize;
    let mut capacity = 64usize;
    let mut slots: Vec<Option<usize>> = vec![None; capacity];
    let mut keys: Vec<i64> = Vec::new();
    let mut groups: Vec<Vec<T>> = Vec::new();
    for (key, value) in pairs {
        let mask = capacity - 1;
        let mut position = slot_hash(key) & mask;
        let mut delta = 1;
        loop {
            match slots[position] {
                None => {
                    slots[position] = Some(groups.len());
                    keys.push(key);
                    groups.push(vec![value]);
                    if groups.len() > (0.7 * capacity as f64) as usize {
                        let doubled = capacity * 2;
                        let mut grown = vec![None; doubled];
                        for slot in slots.iter().flatten() {
                            let mut at = slot_hash(keys[*slot]) & (doubled - 1);
                            let mut step = 1;
                            while grown[at].is_some() {
                                at = (at + step) & (doubled - 1);
                                step += 1;
                            }
                            grown[at] = Some(*slot);
                        }
                        slots = grown;
                        capacity = doubled;
                    }
                    break;
                }
                Some(group) if keys[group] == key => {
                    groups[group].push(value);
                    break;
                }
                Some(_) => {
                    position = (position + delta) & mask;
                    delta += 1;
                }
            }
        }
    }
    let order: Vec<usize> = slots.into_iter().flatten().collect();
    let mut taken: Vec<Option<Vec<T>>> = groups.into_iter().map(Some).collect();
    order
        .into_iter()
        .map(|group| taken[group].take().expect("each group once"))
        .collect()
}

/// `getPairedAndUnpairedLists`: a name seen twice joins the paired list, ahead of the read waiting.
fn paired_and_unpaired(reads: Vec<BamRecord>, guess: i32) -> (Vec<BamRecord>, Vec<BamRecord>) {
    let mut paired = Vec::new();
    let mut waiting: JavaHashMap<String, BamRecord> = JavaHashMap::with_capacity(guess as usize);
    for record in reads {
        match waiting.remove(&record.read_name) {
            Some(mate) => {
                paired.push(record);
                paired.push(mate);
            }
            None => {
                waiting.insert(record.read_name.clone(), record);
            }
        }
    }
    let unpaired = waiting.iter().map(|(_, read)| read.clone()).collect();
    (paired, unpaired)
}

/// `setPairFlags`.
fn set_pair_flags(reads: Vec<BamRecord>, guess: i32) -> Vec<BamRecord> {
    let (mut all, unpaired) = paired_and_unpaired(reads, guess);
    for mut record in unpaired {
        // `setIsPaired(false)`: the proper-pair flag goes with it, the mate flags stay.
        record.flags &= !(0x1 | 0x2);
        all.push(record);
    }
    all
}

/// `filterDuplicateSequences`: the first unpaired read of a sequence, else its first read.
fn filter_duplicate_sequences(reads: Vec<BamRecord>) -> Vec<BamRecord> {
    let pairs = reads
        .into_iter()
        .map(|record| {
            let forward = fnv_byte_array_64(&record.read_bases);
            let reverse = fnv_byte_array_64(&reverse_complement(&record.read_bases));
            (forward.min(reverse), record)
        })
        .collect();
    group_by_key(pairs)
        .into_iter()
        .map(|mut group| {
            let at = group
                .iter()
                .position(|record| !read::is_paired(record))
                .unwrap_or(0);
            group.swap_remove(at)
        })
        .collect()
}

/// `PSFilter.doFilter`, with `load` opening the k-mer file on the first read that needs it, as the
/// task's filter does. `ambiguous` and the rest stay arguments for the caller's read filters.
pub fn do_filter(
    reads: Vec<BamRecord>,
    contigs: &[String],
    arguments: &FilterArguments,
    load_kmers: Option<&dyn Fn() -> Result<KmerLibrary, String>>,
) -> Result<Filtered, String> {
    let mut metrics = FilterMetrics::default();
    let mut reads: Vec<BamRecord> = reads
        .into_iter()
        .filter(|record| !read::is_secondary_alignment(record))
        .filter(|record| !read::is_supplementary_alignment(record))
        .collect();
    metrics.primary = reads.len() as i64;
    if arguments.aligned_input {
        reads.retain(|record| {
            survives_host_alignment(
                record,
                contigs,
                &arguments.contigs_to_ignore,
                arguments.min_identity,
            )
        });
    }
    metrics.after_prealigned_host = reads.len() as i64;
    let mut reads: Vec<BamRecord> = reads.iter().map(cleared).collect();
    reads.iter_mut().for_each(strip_mate_number);

    if !arguments.skip_filters {
        let mut kept = Vec::with_capacity(reads.len());
        for mut record in reads {
            trim_adapters(
                &mut record,
                arguments.max_adapter_mismatches,
                arguments.min_adapter_length,
            );
            mask_simple_repeats(&mut record, 29, 29, 30);
            mask_simple_repeats(&mut record, 87, 89, 100);
            dust(
                &mut record,
                arguments.dust_mask,
                arguments.dust_window,
                arguments.dust_t,
            )?;
            clip_low_quality(&mut record, arguments.read_trim_threshold);
            if (record.read_bases.len() as i32) < arguments.min_read_length {
                continue;
            }
            // `BaseQualityReadTransformer`.
            if record.base_qualities.len() == record.read_bases.len() {
                for at in 0..record.read_bases.len() {
                    if (record.base_qualities[at] as i8 as i32) < arguments.quality_threshold {
                        record.read_bases[at] = b'N';
                    }
                }
            }
            let ambiguous = gatk_readfilter::Parameterized::AmbiguousBase {
                max_fraction: 0.05,
                max_bases: Some(arguments.max_ambiguous_bases),
            };
            if ambiguous.decide(&record).unwrap_or(false) {
                kept.push(record);
            }
        }
        reads = kept;
    }
    metrics.after_quality = reads.len() as i64;

    if let Some(load) = load_kmers {
        if !reads.is_empty() {
            let library = load()?;
            let size = library.kmer_size();
            reads.retain(|record| {
                let mut found = 0;
                for kmer in kmerize(&record.read_bases, size, 1) {
                    if library.contains(kmer) {
                        found += 1;
                        if found >= arguments.host_kmer_threshold {
                            return false;
                        }
                    }
                }
                true
            });
        }
    }
    metrics.after_host = reads.len() as i64;

    let guess = arguments.reads_per_partition;
    if arguments.filter_duplicates {
        reads = filter_duplicate_sequences(set_pair_flags(reads, guess));
    }
    metrics.after_deduplication = reads.len() as i64;

    let reads: Vec<BamRecord> = set_pair_flags(reads, guess).iter().map(cleared).collect();
    let (paired, unpaired) = paired_and_unpaired(reads, guess);
    metrics.final_paired = paired.len() as i64;
    Ok(Filtered {
        paired,
        unpaired,
        metrics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn murmur_matches_guava_on_zero() {
        // `Hashing.murmur3_32().hashInt(0)`.
        assert_eq!(murmur3_32_int(0) as u32, 0x2362f9de);
    }

    #[test]
    fn a_mate_suffix_goes_and_only_a_mate_suffix() {
        let mut record = BamRecord {
            read_name: "a/1".to_string(),
            ..BamRecord::default()
        };
        strip_mate_number(&mut record);
        assert_eq!(record.read_name, "a");
        record.read_name = "a/3".to_string();
        strip_mate_number(&mut record);
        assert_eq!(record.read_name, "a/3");
    }

    #[test]
    fn dust_masks_a_poly_a() {
        let bases = vec![b'A'; 60];
        let intervals = run_dust(&bases, 64, 20.0).unwrap();
        assert_eq!(intervals, [(0, 59)]);
        let plain = b"ACGTTGCAAGGCTTAACCGATCGATGCATCGGATTACGCTAGCTAGGATCCTAGT";
        assert_eq!(run_dust(plain, 64, 20.0).unwrap(), []);
    }
}
