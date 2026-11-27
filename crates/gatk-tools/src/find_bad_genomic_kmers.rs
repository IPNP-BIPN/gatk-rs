//! `FindBadGenomicKmersSpark` (GATK 4.7.0.0; `FindBadGenomicKmersSpark.java` is unchanged since
//! 4.6.2.0): the k-mers a reference holds more than three times.
//!
//! # The order is the partitions' and the tables'
//!
//! The file is written from `collect()`, so nothing sorts it. `getReferenceBasesRDD` cuts each
//! contig into records of 10000 bases overlapping by K-1 and `parallelize`s them into
//! `records / 104 + 1` slices. Each slice counts its k-mers in its own hopscotch map, which it
//! iterates into a `HashPartitioner` over the same number of partitions, and each reducer sums what
//! arrives, map by map, in another map, keeps the counts over three, and iterates that. The
//! partitions come back in order. So what decides a k-mer's place is the order of the slices and of the two
//! tables, all of which are reproduced here.
//!
//! # The second pass is another table
//!
//! `--high-copy-fasta` adds every k-mer of that file, repeats and all, through `SVUtils.uniquify`:
//! a hopscotch set sized for both lists, holding the first list then the second, and iterated.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, KmerFailure>` is "a value or an error"; each [`KmerFailure`] variant is one of the
//!   exceptions the Java would throw. `?` after a call hands its error straight to the caller,
//!   Java's `throws` at that one call site.
//! - `impl From<HopscotchFailed> for KmerFailure` lets `?` turn a table failure into a
//!   `KmerFailure` on the way out, with no explicit conversion at the call.
//! - `&[Vec<u8>]` is a borrowed, read-only list of byte arrays; `&mut Vec<..>` is a borrowed list
//!   the function appends to, which the caller keeps.
//! - `as` between integer types converts without a check: `as i32` from a `usize` drops the high
//!   bits, `as usize` from a non-negative `i64` keeps the value.

use crate::hopscotch::{Hopscotch, HopscotchEntry, HopscotchFailed};
use crate::sv_kmer::{DustFilteredKmerizer, SvKmerLong};

/// What stops the tool once a k-mer exists.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KmerFailure {
    /// `canonical`'s `validateArg`: an even K runs until its first k-mer.
    EvenSize,
    /// `charAt(0)` on an empty line of the high-copy FASTA.
    EmptyLine,
    /// The hopscotch table's own `IllegalStateException`, out of a resize.
    Table(HopscotchFailed),
}

impl From<HopscotchFailed> for KmerFailure {
    fn from(failed: HopscotchFailed) -> KmerFailure {
        KmerFailure::Table(failed)
    }
}

/// `SVKmerLong.canonical(kSize)`, with its argument check made a refusal:
///
/// ```java
///         Utils.validateArg( (kSize & 1) != 0, "Kmer length must be odd to canonicalize.");
/// ```
fn canonical(kmer: SvKmerLong, k: i32) -> Result<SvKmerLong, KmerFailure> {
    if k & 1 == 0 {
        return Err(KmerFailure::EvenSize);
    }
    Ok(kmer.canonical(k))
}

/// `REF_RECORD_LEN`: the length of the records `getReferenceBasesRDD` cuts a contig into.
pub const REF_RECORD_LEN: usize = 10000;
/// `REF_RECORDS_PER_PARTITION`: `1024*1024 / REF_RECORD_LEN`, integer division, so 104.
pub const REF_RECORDS_PER_PARTITION: usize = 1024 * 1024 / REF_RECORD_LEN;
/// `MAX_KMER_FREQ`: a k-mer is kept when its count is strictly greater than this.
pub const MAX_KMER_FREQ: i32 = 3;

/// `KmerAndCount`, keyed on the k-mer alone: the hash and the equality the hopscotch map sees are
/// the k-mer's, so a count never splits one k-mer into two entries.
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

/// An `SVKmer` as `HopscotchSetSpark` holds it in `SVUtils.uniquify`: keyed on itself.
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

/// `SVReferenceUtils.getReferenceBasesRDD`'s chunks of one contig, appended to `chunks` so that
/// the contigs follow one another in dictionary order:
///
/// ```java
///         final int effectiveRecLen = refRecordLen - kSize + 1;
/// ```
///
/// ```java
///                 for ( int start = 0; start < seqLen; start += effectiveRecLen ) {
/// ```
///
/// and each record is `bases[start, min(start + refRecordLen, seqLen))`, so records overlap by
/// K-1 bases and the last one is short.
///
/// A step that does not advance is the reference's unbounded loop; the caller has refused it.
/// The arithmetic is in `i64` so that a K past the record length gives a step of zero or less
/// rather than an unsigned wrap.
pub fn chunk_contig(bases: &[u8], k: i32, chunks: &mut Vec<Vec<u8>>) {
    let step = REF_RECORD_LEN as i64 - k as i64 + 1;
    let mut start = 0i64;
    while start < bases.len() as i64 {
        let from = start as usize;
        let to = (start as usize + REF_RECORD_LEN).min(bases.len());
        chunks.push(bases[from..to].to_vec());
        start += step;
    }
}

/// `kmerCounts.find(kmer)`, then `add(new KmerAndCount(kmer, by))` or `entry.bumpCount(by)`.
///
/// `bumpCount` is `count += extra` on a Java `int`, which wraps on overflow; `wrapping_add` is
/// that wrap written out, where Rust's plain `+` would stop a debug build.
fn bump(
    table: &mut Hopscotch<KmerAndCount>,
    kmer: SvKmerLong,
    by: i32,
) -> Result<(), HopscotchFailed> {
    match table.find_mut(kmer.hash_code(), |entry| entry.kmer == kmer) {
        Some(entry) => entry.count = entry.count.wrapping_add(by),
        None => {
            table.add(KmerAndCount { kmer, count: by })?;
        }
    }
    Ok(())
}

/// `collectUbiquitousKmersInReference`, with Spark's two stages run one after the other.
///
/// **What**: the canonical k-mers whose count over all chunks is above `max_frequency`.
///
/// **How**: the first loop is the first `mapPartitions`: one hopscotch map per slice, iterated
/// into the target partition of `HashPartitioner(nPartitions)`. The second loop is the second
/// `mapPartitions`, each partition summing what arrived into a fresh map and iterating it.
///
/// **Why it is not one global count**: the output has no sort, so its order is `collect()`'s:
/// partition by partition and, inside each, the iteration order of that partition's map. One
/// global map would hold the same k-mers in a different order.
pub fn collect_ubiquitous_kmers(
    k: i32,
    max_dust: i32,
    max_frequency: i32,
    chunks: &[Vec<u8>],
) -> Result<Vec<SvKmerLong>, KmerFailure> {
    // `ctx.parallelize(sequenceChunks, sequenceChunks.size()/refRecordsPerPartition+1)`.
    let partitions = chunks.len() / REF_RECORDS_PER_PARTITION + 1;
    let hash_size = (2 * REF_RECORDS_PER_PARTITION) as i32;

    let mut shuffled: Vec<Vec<KmerAndCount>> = vec![Vec::new(); partitions];
    for slice in 0..partitions {
        // `ParallelCollectionRDD.slice`'s positions.
        let start = slice * chunks.len() / partitions;
        let end = (slice + 1) * chunks.len() / partitions;
        let mut counts = Hopscotch::new(hash_size);
        for chunk in &chunks[start..end] {
            // `SVDUSTFilteredKmerizer.canonicalStream`: the kmerizer yields the k-mers that pass
            // the DUST filter, and `canonical` is applied to each.
            for kmer in DustFilteredKmerizer::new(chunk, k, max_dust) {
                bump(&mut counts, canonical(kmer, k)?, 1)?;
            }
        }
        for entry in counts.iter() {
            // `Partitioner.nonNegativeMod`.
            let target = entry.kmer.hash_code().rem_euclid(partitions as i32) as usize;
            shuffled[target].push(entry.clone());
        }
    }

    // Each partition in turn, as `collect()` concatenates them.
    let mut kept = Vec::new();
    for arrived in &shuffled {
        let mut counts = Hopscotch::new(hash_size);
        for entry in arrived {
            bump(&mut counts, entry.kmer, entry.count)?;
        }
        kept.extend(
            counts
                .iter()
                .filter(|entry| entry.count > max_frequency)
                .map(|entry| entry.kmer),
        );
    }
    Ok(kept)
}

/// `processFasta` over what `readLine` handed back, whose first character decides each line:
///
/// ```java
///                 if ( line.charAt(0) != '>' ) sb.append(line);
///                 else if ( sb.length() > 0 ) {
/// ```
///
/// So the sequence lines of one record are joined before they are cut into k-mers, a header ends
/// the record, and an empty line throws on `charAt(0)` ([`KmerFailure::EmptyLine`]).
pub fn process_fasta(
    lines: &[Vec<u8>],
    k: i32,
    max_dust: i32,
) -> Result<Vec<SvKmerLong>, KmerFailure> {
    let mut kmers = Vec::new();
    let mut sequence = Vec::new();
    // A closure, Java's lambda: the `canonicalStream(...).forEach(kmers::add)` both the header
    // branch and the end of the file run.
    let flush = |sequence: &Vec<u8>, kmers: &mut Vec<SvKmerLong>| -> Result<(), KmerFailure> {
        for kmer in DustFilteredKmerizer::new(sequence, k, max_dust) {
            kmers.push(canonical(kmer, k)?);
        }
        Ok(())
    };
    // `line.first()` is `None` for an empty line, and `ok_or(..)?` turns that into the refusal.
    for line in lines {
        match line.first().ok_or(KmerFailure::EmptyLine)? {
            b'>' => {
                if !sequence.is_empty() {
                    flush(&sequence, &mut kmers)?;
                    sequence.clear();
                }
            }
            _ => sequence.extend_from_slice(line),
        }
    }
    if !sequence.is_empty() {
        flush(&sequence, &mut kmers)?;
    }
    Ok(kmers)
}

/// `SVUtils.uniquify`: a `HopscotchSetSpark` sized `coll1.size() + coll2.size()`, given the first
/// list then the second, and returned as the set itself, so its iteration order is the output
/// order. A k-mer already held is not added again.
pub fn uniquify(
    first: &[SvKmerLong],
    second: &[SvKmerLong],
) -> Result<Vec<SvKmerLong>, HopscotchFailed> {
    let mut set = Hopscotch::new((first.len() + second.len()) as i32);
    for kmer in first.iter().chain(second) {
        set.add(Kmer(*kmer))?;
    }
    Ok(set.iter().map(|kmer| kmer.0).collect())
}

/// `SVFileUtils.writeKmersFile`: one k-mer per line, as K bases of text, in the order given.
pub fn kmers_file(kmers: &[SvKmerLong], k: i32) -> String {
    let mut text = String::new();
    for kmer in kmers {
        text.push_str(&kmer.to_text(k));
        text.push('\n');
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn records_overlap_by_k_less_one() {
        let bases = vec![b'A'; 25000];
        let mut chunks = Vec::new();
        chunk_contig(&bases, 51, &mut chunks);
        let lengths: Vec<usize> = chunks.iter().map(Vec::len).collect();
        // Starts 0, 9950, 19900.
        assert_eq!(lengths, [10000, 10000, 5100]);
    }

    #[test]
    fn a_kmer_needs_more_than_the_frequency() {
        let unit = b"ACGTTGCAAGGCTTAACCGATCGAT".to_vec();
        let mut bases = Vec::new();
        for _ in 0..3 {
            bases.extend_from_slice(&unit);
            bases.extend_from_slice(b"NNNN");
        }
        let few = collect_ubiquitous_kmers(21, 100, 3, std::slice::from_ref(&bases)).unwrap();
        assert!(few.is_empty());
        let many = collect_ubiquitous_kmers(21, 100, 2, std::slice::from_ref(&bases)).unwrap();
        assert_eq!(many.len(), 5);
    }

    #[test]
    fn an_empty_line_is_the_reference_throwing() {
        assert_eq!(
            process_fasta(&[b">a".to_vec(), Vec::new()], 5, 3),
            Err(KmerFailure::EmptyLine)
        );
    }
}
