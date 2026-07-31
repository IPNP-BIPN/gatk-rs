//! `Fragment` and `AlleleLikelihoods.groupEvidence`, ported from GATK 4.6.2.0.
//!
//! A fragment is "all available evidence coming from a single biological fragment", which is one
//! read or a read pair. Two annotations count fragments rather than reads, and the difference is
//! not cosmetic: a read pair that straddles a variant votes **once**.
//!
//! # Grouping **sums** the log likelihoods, it does not average them
//!
//! ```java
//! newLikelihoodValues[s][a][newEvidenceIndex] += oldSampleValues[a][oldEvidenceIndex];
//! ```
//!
//! "corresponding to an independent evidence assumption. Since this container's likelihoods
//! generally pertain to sequencing only (and not sample prep etc) this is usually a good
//! assumption." So a pair of reads that each support an allele at -1 give the fragment a -2, and the
//! informativeness threshold, which is an absolute difference, is therefore easier to clear after
//! grouping than before.
//!
//! # The fragment order is a `HashMap`'s
//!
//! ```java
//! new ArrayList<>(sampleEvidence(s).stream().collect(Collectors.groupingBy(groupingFunction)).values())
//! ```
//!
//! `Collectors.groupingBy` builds a `HashMap`, so the new evidence order is hash order over the
//! grouping key, which for these annotations is the read name. The order **within** each group is
//! the stream's, so the reads of a pair keep the sample's order and "the first read of each
//! fragment" is well defined. [`gatk_engine::java_hash::hash_map_order`] reproduces the outer one.
//!
//! # More than two reads with one name is a warning, not an error
//!
//! `Fragment.createAndAvoidFailure` drops duplicates, secondary and supplementary alignments; if
//! more than two survive it logs "Using two reads randomly to combine as a fragment" and takes the
//! **first two**, which is not random but is the sublist. If none survive it takes the first read of
//! the original list, supplementary or not.

use crate::java_hash::{hash_map_order, string_hash_code};
use htsjdk_bam::record::BamRecord;

/// `Fragment`: one read or a pair, with the interval spanning them.
#[derive(Debug, Clone, PartialEq)]
pub struct Fragment {
    pub contig_index: i32,
    pub start: i32,
    pub end: i32,
    /// One or two reads, in the order they appeared in the sample.
    pub reads: Vec<BamRecord>,
}

/// What fragment construction refuses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FragmentError {
    /// `Utils.validateArg(!reads.isEmpty(), "Need one or two reads to construct a fragment")`.
    NoReads,
    /// `Utils.validateArg(reads.size() <= 2, "Cannot construct fragment from more than two reads")`,
    /// which only `create` raises: `createAndAvoidFailure` trims instead.
    TooManyReads { found: usize },
}

/// `GATKRead.getEnd()` for an aligned read: the last reference position the cigar covers.
fn read_end(read: &BamRecord) -> i32 {
    read.alignment_end()
}

impl Fragment {
    /// `Fragment.create(List<GATKRead>)`.
    ///
    /// The interval takes `min` and `max` of the two reads' starts and ends and then `min`/`max`
    /// again, because a read whose end precedes its start (an unmapped one) would otherwise build an
    /// interval the constructor rejects.
    pub fn create(reads: &[BamRecord]) -> Result<Fragment, FragmentError> {
        match reads.len() {
            0 => Err(FragmentError::NoReads),
            1 => {
                let read = &reads[0];
                let end = read_end(read);
                Ok(Fragment {
                    contig_index: read.reference_index,
                    start: read.alignment_start.min(end),
                    end: read.alignment_start.max(end),
                    reads: reads.to_vec(),
                })
            }
            2 => {
                let (left, right) = (&reads[0], &reads[1]);
                let start = left.alignment_start.min(right.alignment_start);
                let end = read_end(left).max(read_end(right));
                Ok(Fragment {
                    contig_index: left.reference_index,
                    start: start.min(end),
                    end: start.max(end),
                    reads: reads.to_vec(),
                })
            }
            found => Err(FragmentError::TooManyReads { found }),
        }
    }

    /// `Fragment.createAndAvoidFailure`.
    pub fn create_and_avoid_failure(reads: &[BamRecord]) -> Result<Fragment, FragmentError> {
        if reads.len() <= 2 {
            return Fragment::create(reads);
        }
        let primary: Vec<BamRecord> = reads
            .iter()
            .filter(|read| {
                let flags = read.flags;
                // Duplicate, secondary and supplementary respectively.
                flags & 0x400 == 0 && flags & 0x100 == 0 && flags & 0x800 == 0
            })
            .cloned()
            .collect();
        if primary.len() > 2 {
            // "Using two reads randomly to combine as a fragment", which takes the first two.
            return Fragment::create(&primary[..2]);
        }
        if primary.is_empty() {
            return Fragment::create(&reads[..1]);
        }
        Fragment::create(&primary)
    }

    /// `ReadUtils.isF2R1`: reverse-strand and first-of-pair agree.
    ///
    /// A read that is not paired at all has `isFirstOfPair()` false, so a **forward** unpaired read
    /// is `F2R1` (false equals false) and a reverse one is `F1R2`. The orientation is defined
    /// whether or not there is a mate, and for an unpaired read it is the opposite of what the
    /// names suggest.
    pub fn is_f2r1(read: &BamRecord) -> bool {
        let reverse = read.flags & 0x10 != 0;
        let first_of_pair = read.flags & 0x40 != 0;
        reverse == first_of_pair
    }
}

/// The grouping key and the resulting order of `AlleleLikelihoods.groupEvidence` by read name.
///
/// Returns, per sample, the groups in the order the reference's `HashMap` iteration produces, each
/// group holding the indices of its reads within that sample in the sample's own order.
pub fn group_by_read_name(reads: &[BamRecord]) -> Vec<Vec<usize>> {
    // Insertion order first, as `groupingBy` fills the map.
    let mut names: Vec<String> = Vec::new();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (index, read) in reads.iter().enumerate() {
        match names.iter().position(|name| *name == read.read_name) {
            Some(position) => groups[position].push(index),
            None => {
                names.push(read.read_name.clone());
                groups.push(vec![index]);
            }
        }
    }
    let entries: Vec<(usize, i32)> = (0..names.len())
        .map(|index| (index, string_hash_code(&names[index])))
        .collect();
    match hash_map_order(&entries) {
        Ok(order) => order
            .into_iter()
            .map(|index| groups[index].clone())
            .collect(),
        // A bucket past the treeify threshold, which needs eight names colliding. Unreachable for
        // read names, and the insertion order is the honest fallback.
        Err(_) => groups,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(name: &str, flags: u16, start: i32) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            flags,
            reference_index: 0,
            alignment_start: start,
            mapping_quality: 60,
            cigar: htsjdk_bam::text_parse::parse_cigar("10M").expect("a cigar"),
            read_bases: vec![b'A'; 10],
            base_qualities: vec![30; 10],
            ..Default::default()
        }
    }

    #[test]
    fn a_pair_spans_both_reads() {
        let fragment =
            Fragment::create(&[read("r", 0x41, 100), read("r", 0x81, 200)]).expect("a fragment");
        assert_eq!((fragment.start, fragment.end), (100, 209));
    }

    #[test]
    fn an_unpaired_forward_read_is_f2r1() {
        // false == false, so the unpaired forward read lands in the F2R1 bucket.
        assert!(Fragment::is_f2r1(&read("r", 0, 100)));
        assert!(!Fragment::is_f2r1(&read("r", 0x10, 100)));
        // A paired first-of-pair forward read is F1R2, which is the ordinary case.
        assert!(!Fragment::is_f2r1(&read("r", 0x41, 100)));
    }

    #[test]
    fn a_pair_is_one_group() {
        let groups = group_by_read_name(&[
            read("a", 0x41, 100),
            read("b", 0x41, 100),
            read("a", 0x81, 200),
        ]);
        assert_eq!(groups.len(), 2);
        assert!(groups.contains(&vec![0, 2]));
        assert!(groups.contains(&vec![1]));
    }
}
