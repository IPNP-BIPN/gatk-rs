//! `--error-correction-log-odds`, ported from `tools.walkers.haplotypecaller.PileupReadErrorCorrector`
//! (GATK 4.6.2.0).
//!
//! Before assembly, a copy of the region's reads is piled up locus by locus. At each locus the
//! plurality of A, C, G and T (the first of a tie, in that order) stands for the reference, and
//! when Mutect2's log likelihood ratio of the other bases' qualities against it is below the
//! threshold, every other base not next to an indel or a soft clip is a candidate correction to
//! it. A read takes its candidates, set to quality 30, except at either end where three fall within
//! fifteen bases, which may be an indel the aligner missed. Only the assembly sees the corrections.
//!
//! The pileup concatenates samples in `getSamplesFromHeader`'s sorted order, and the quality sum
//! inside the ratio is taken in that order.

use std::collections::{BTreeSet, HashMap};

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use crate::locus_iterator::{self, LocusIteratorOptions};
use crate::read_states::{DownsamplingInfo, ReadStateManager};

const GOOD_QUAL: u8 = 30;
const INDEL_SPAN: i32 = 15;
const INDEL_MISMATCHES: usize = 3;

/// `correctReads(reads)`: the reads, copied and corrected, in their order.
pub fn correct_reads(
    reads: &[BamRecord],
    header: &SamHeader,
    log_odds_threshold: f64,
) -> Result<Vec<BamRecord>, String> {
    // `ReadUtils.getSamplesFromHeader`: a `TreeSet` of the read groups' samples.
    let samples: Vec<Option<String>> = header
        .read_groups
        .iter()
        .filter_map(|group| group.attributes.get("SM").map(str::to_string))
        .collect::<BTreeSet<String>>()
        .into_iter()
        .map(Some)
        .collect();
    let states = ReadStateManager::new(samples.clone(), DownsamplingInfo::NONE)
        .map_err(|e| format!("{e:?}"))?;
    let contexts = locus_iterator::contexts(
        reads,
        samples,
        header,
        LocusIteratorOptions::default(),
        states,
    )
    .map_err(|e| format!("{e:?}"))?;
    let index_of: HashMap<*const BamRecord, usize> = reads
        .iter()
        .enumerate()
        .map(|(i, read)| (read as *const BamRecord, i))
        .collect();
    let mut edits: Vec<Vec<(i32, u8)>> = vec![Vec::new(); reads.len()];
    for context in &contexts {
        let pileup = &context.pileup.elements;
        let mut counts = [0i64; 4];
        for element in pileup {
            if let Some(i) = b"ACGT".iter().position(|b| *b == element.base()) {
                counts[i] += 1;
            }
        }
        // `max(comparingLong(counter::get))` keeps the first of equals.
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
                    let read = index_of[&(element.read as *const BamRecord)];
                    edits[read].push((element.offset, reference));
                }
            }
        }
    }
    let mut corrected = reads.to_vec();
    for (read, edits) in corrected.iter_mut().zip(&edits) {
        let size = edits.len();
        let mut first_edit = 0usize;
        let mut n = 0usize;
        while n + INDEL_MISMATCHES < size
            && edits[n + INDEL_MISMATCHES - 1].0 - edits[n].0 < INDEL_SPAN
        {
            first_edit = n + INDEL_MISMATCHES;
            n += 1;
        }
        let mut last_edit = size as i64 - 1;
        let mut n = size as i64 - 1;
        while n >= INDEL_MISMATCHES as i64 - 1
            && edits[n as usize].0 - edits[(n - INDEL_MISMATCHES as i64 + 1) as usize].0
                < INDEL_SPAN
        {
            last_edit = n - INDEL_MISMATCHES as i64;
            n -= 1;
        }
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
