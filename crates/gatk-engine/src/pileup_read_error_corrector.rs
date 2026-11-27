//! `--error-correction-log-odds`, ported from
//! `tools.walkers.haplotypecaller.PileupReadErrorCorrector` (GATK 4.7.0.0;
//! `PileupReadErrorCorrector.java` is unchanged since 4.6.2.0).
//!
//! Before assembly, a copy of the region's reads is piled up locus by locus. At each locus the
//! plurality of A, C, G and T (the first of a tie, in that order) stands for the reference, and
//! when Mutect2's log likelihood ratio of the other bases' qualities against it is below the
//! threshold, every other base that is not a deletion, not beside an indel and not just after a
//! soft clip is a candidate correction to it (the reference tests `isAfterSoftClip` and not
//! `isBeforeSoftClip`). A read takes its candidates, set to quality 30, except at either end where
//! three fall within fifteen bases, which may be an indel the aligner missed. Only the assembly
//! sees the corrections.
//!
//! The pileup concatenates samples in `getSamplesFromHeader`'s sorted order, and the quality sum
//! inside the ratio is taken in that order.
//!
//! # Reading this file without reading Rust
//!
//! - `&[BamRecord]` is a borrowed, read-only array of reads; the caller keeps them, and
//!   `reads.to_vec()` is the copy that gets corrected (the reference's `GATKRead::deepCopy`).
//! - `*const BamRecord` is a raw pointer, used here only as an identity: two pointers are equal
//!   when they point at the same read, which is how a pileup element is traced back to its index
//!   in `reads`. It is never dereferenced.
//! - `Result<_, String>` is "a value or an error message", and `.map_err(..)?` turns the
//!   callee's error into that message and returns it at once.
//! - `i64` loop counters below are where the reference's `int` index can go negative; an index
//!   (`usize`) cannot, so the backward scan counts in `i64` and casts back only to index.

use std::collections::{BTreeSet, HashMap};

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use crate::locus_iterator::{self, LocusIteratorOptions};
use crate::read_states::{DownsamplingInfo, ReadStateManager};

/// `GOOD_QUAL`: the quality a corrected base is given.
const GOOD_QUAL: u8 = 30;
/// `INDEL_SPAN` and `INDEL_MISMATCHES`: three candidate corrections within fifteen bases of each
/// other, at either end of a read's list, are taken for a missed indel and left alone.
const INDEL_SPAN: i32 = 15;
const INDEL_MISMATCHES: usize = 3;

/// `correctReads(reads)`: the reads, copied and corrected, in their order.
///
/// The reference iterates its `LocusIteratorByState` lazily and gathers each read's candidate
/// edits in a map keyed by the read; this collects every locus first and keeps the edits in a
/// list per read index. The edits of one read are appended in locus order either way.
///
/// `log_odds_threshold` is `--error-correction-log-odds` (the argument collection's
/// `pileupErrorCorrectionLogOdds`). An `Err` is the locus iterator's refusal, rendered as text.
pub fn correct_reads(
    reads: &[BamRecord],
    header: &SamHeader,
    log_odds_threshold: f64,
) -> Result<Vec<BamRecord>, String> {
    // `ReadUtils.getSamplesFromHeader`: a `TreeSet` of the read groups' samples. `filter_map`
    // drops a read group with no `SM`, as the reference's `if ( sample != null )` does, and the
    // `BTreeSet` is the sorted, de-duplicated set; `.map(Some)` wraps each name in an `Option`
    // because the state manager also takes a `None` sample.
    let samples: Vec<Option<String>> = header
        .read_groups
        .iter()
        .filter_map(|group| group.attributes.get("SM").map(str::to_string))
        .collect::<BTreeSet<String>>()
        .into_iter()
        .map(Some)
        .collect();
    // `DownsamplingMethod.NONE`.
    let states = ReadStateManager::new(samples.clone(), DownsamplingInfo::NONE)
        .map_err(|e| format!("{e:?}"))?;
    // The reference's iterator:
    // `new LocusIteratorByState(reads.iterator(), DownsamplingMethod.NONE,
    // ReadUtils.getSamplesFromHeader(header), header, false)`. Its last argument is
    // `includeReadsWithDeletionAtLoci`, and the reference passes false. The default options
    // passed here have `include_deletions: true`, so a read inside a deletion is in this pileup
    // and its Q16 enters `alt_quals` below, where the reference leaves the read out.
    let contexts = locus_iterator::contexts(
        reads,
        samples,
        header,
        LocusIteratorOptions::default(),
        states,
    )
    .map_err(|e| format!("{e:?}"))?;
    // Each read's index, keyed by its address: a pileup element borrows its read from `reads`,
    // so the element's pointer finds the read it came from.
    let index_of: HashMap<*const BamRecord, usize> = reads
        .iter()
        .enumerate()
        .map(|(i, read)| (read as *const BamRecord, i))
        .collect();
    // One list of `(offset, base)` per read: the reference's `List<Pair<Integer, Byte>>`.
    let mut edits: Vec<Vec<(i32, u8)>> = vec![Vec::new(); reads.len()];
    for context in &contexts {
        let pileup = &context.pileup.elements;
        // `Nucleotide.Counter` over `pe.getBase()`, for the four standard bases only. The test
        // here is upper case only, where `Nucleotide.decode` also maps `a`, `c`, `g` and `t`.
        let mut counts = [0i64; 4];
        for element in pileup {
            if let Some(i) = b"ACGT".iter().position(|b| *b == element.base()) {
                counts[i] += 1;
            }
        }
        // `max(comparingLong(counter::get))` keeps the first of equals. Strictly greater, so a
        // tie goes to the earlier of A, C, G, T, and an all-deletion pileup answers A.
        let mut best = 0;
        for i in 1..4 {
            if counts[i] > counts[best] {
                best = i;
            }
        }
        let reference = b"ACGT"[best];
        let mut ref_count = 0;
        let mut alt_quals: Vec<u8> = Vec::new();
        for element in pileup {
            if element.base() == reference {
                ref_count += 1;
            } else {
                alt_quals.push(element.qual());
            }
        }
        // The two-argument `Mutect2Engine.logLikelihoodRatio(refCount, altQualBuffer)`, which
        // calls the three-argument one with a repeat factor of 1.
        let log_odds = crate::mutect_engine::log_likelihood_ratio(ref_count, &alt_quals, 1);
        if log_odds < log_odds_threshold {
            for element in pileup {
                if element.base() != reference
                    && !(element.is_deletion()
                        || element.is_before_insertion()
                        || element.is_after_deletion_end()
                        || element.is_before_deletion_start()
                        || element.is_after_insertion()
                        || element.is_after_soft_clip())
                {
                    // `potentialCorrections.get(pe.getRead())
                    // .add(ImmutablePair.of(pe.getOffset(), ref))`.
                    let read = index_of[&(element.read as *const BamRecord)];
                    edits[read].push((element.offset, reference));
                }
            }
        }
    }
    let mut corrected = reads.to_vec();
    // `.iter_mut().zip(&edits)` walks each copied read beside its own edit list; `iter_mut`
    // hands out a writable reference, so the assignments below change the copy in place.
    for (read, edits) in corrected.iter_mut().zip(&edits) {
        let size = edits.len();
        // The reference's `for` loops, written as `while`: skip leading edits for as long as each
        // window of three starts within `INDEL_SPAN` of where it ends. `.0` is the pair's first
        // member, the offset.
        let mut first_edit = 0usize;
        let mut n = 0usize;
        while n + INDEL_MISMATCHES < size
            && edits[n + INDEL_MISMATCHES - 1].0 - edits[n].0 < INDEL_SPAN
        {
            first_edit = n + INDEL_MISMATCHES;
            n += 1;
        }
        // The same from the far end. `lastEdit` starts at `size - 1`, which is -1 for a read with
        // no edits, and `n - INDEL_MISMATCHES` can go negative too, hence `i64`.
        let mut last_edit = size as i64 - 1;
        let mut n = size as i64 - 1;
        while n >= INDEL_MISMATCHES as i64 - 1
            && edits[n as usize].0 - edits[(n - INDEL_MISMATCHES as i64 + 1) as usize].0
                < INDEL_SPAN
        {
            last_edit = n - INDEL_MISMATCHES as i64;
            n -= 1;
        }
        // Apply what is left between the two cut points: the base and a quality of 30.
        let mut n = first_edit as i64;
        while n <= last_edit {
            let (offset, base) = edits[n as usize];
            read.read_bases[offset as usize] = base;
            read.base_qualities[offset as usize] = GOOD_QUAL;
            n += 1;
        }
    }
    Ok(corrected)
}
