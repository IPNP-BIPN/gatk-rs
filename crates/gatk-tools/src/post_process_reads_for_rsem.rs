//! Ported from `org.broadinstitute.hellbender.tools.walkers.qc.PostProcessReadsForRSEM` and its
//! `ReadPair` (GATK 4.6.2.0).
//!
//! The twelfth whole tool of the record-transform archetype, the third that is not a walker, and the
//! first that groups the traversal into query-name runs of its own accord.
//!
//! # A fourth `getDefaultReadFilters` pattern
//!
//! ```java
//! public List<ReadFilter> getDefaultReadFilters() {
//!     return Collections.singletonList(ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT);
//! }
//! ```
//!
//! [`crate::print_reads`] takes `GATKTool`'s default, [`crate::unmark_duplicates`] replaces it with
//! `ALLOW_ALL_READS`, [`crate::print_distant_mates`] extends it with four more. This one replaces
//! the whole list with a **single** filter that is not `Wellformed`, so a supplementary alignment
//! never reaches the tool and the `supplementaryAlignments` list `ReadPair` maintains for it is
//! dead code on this path. Reproduced anyway, because the reference maintains it.
//!
//! # Two of its own guards dereference null
//!
//! ```java
//! if (read1 == null || read2 == null){
//!     logger.warn("read1 or read2 is null. This read will not be output. " + read1.getName());
//!     return false;
//! }
//! ```
//!
//! The branch exists because `read1` may be null and dereferences it. A query-name group holding
//! only a second-of-pair therefore kills the run:
//! `Cannot invoke "GATKRead.getName()" because "read1" is null`.
//!
//! `groupSecondaryReads` has the same shape. It collects into
//! `Collectors.groupingBy(GATKRead::isFirstOfPair)`, which produces no `false` key at all when every
//! secondary is a first-of-pair, and then calls `read2Reads.size()` on the line before the guard
//! that would have caught it: `Cannot invoke "List.size()" because "read2Reads" is null`.
//!
//! Both are [`RsemError::NullDereference`] here rather than a panic. A port that returned an empty
//! result, or skipped the pair, would write a file the reference never wrote: the reference writes
//! **nothing at all**, because the exception escapes the traversal.
//!
//! # The output order is not the input order
//!
//! Primary pair first, then each secondary pair, and within each pair first-of-pair before
//! second-of-pair. A secondary finds its mate by contig plus `getStart() == mate.getMateStart()` and
//! `getMateStart() == mate.getStart()`, not by anything in the record that says the two belong
//! together.
//!
//! # Three reasons drop both reads, and one of them drops more
//!
//! Either read unmapped, the two on different contigs, or a cigar that is not exactly one `M`.
//! `100=` is a single element and is still refused, because the test is on the operator. A primary
//! pair that fails takes its secondary alignments with it, since the secondary loop sits inside the
//! primary's `if`; a secondary pair that fails on its own drops only itself.

use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;

use gatk_engine::read;
use gatk_engine::reads::{ReadsDataSource, ReadsError};

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK PostProcessReadsForRSEM";

/// `getDefaultReadFilters()`: one filter, and it is not `Wellformed`.
pub const DEFAULT_READ_FILTERS: [&str; 1] = ["NotSupplementaryAlignmentReadFilter"];

/// `ReadFilterLibrary.NOT_SUPPLEMENTARY_ALIGNMENT`, which is this tool's whole filter chain.
pub fn default_read_filter(record: &BamRecord) -> bool {
    !read::is_supplementary_alignment(record)
}

/// What ends a run, rather than dropping a pair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RsemError {
    /// `UserException("Input must be query-name sorted.")`, from `onTraversalStart`.
    NotQueryNameSorted,
    /// `UserException("Read names do not match: ...")`, from `ReadPair.add`.
    NameMismatch { expected: String, found: String },
    /// `Utils.validate`: a second primary first-of-pair or second-of-pair in one group.
    PrimaryAlreadySet { name: String, first: bool },
    /// `UserException("Unknown read type: " + read.getContig())`, which `ReadPair.add` raises for a
    /// primary read that is neither first nor second of pair.
    UnknownReadType { contig: String },
    /// One of the two `NullPointerException`s the reference throws out of its own null guards. The
    /// message is the JVM's, and it is what the golden carries.
    NullDereference(&'static str),
}

impl RsemError {
    /// The message the reference raises, for the ones worth comparing against a golden.
    pub fn message(&self) -> String {
        match self {
            RsemError::NotQueryNameSorted => "Input must be query-name sorted.".to_string(),
            RsemError::NameMismatch { expected, found } => {
                format!("Read names do not match: {expected} vs {found}")
            }
            RsemError::PrimaryAlreadySet { name, first } => format!(
                "The primary {} is already set. Read = {name}",
                if *first {
                    "firstOfPair"
                } else {
                    "secondOfPair"
                }
            ),
            RsemError::UnknownReadType { contig } => format!("Unknown read type: {contig}"),
            RsemError::NullDereference(what) => (*what).to_string(),
        }
    }
}

/// The JVM's message for `read1.getName()` inside the guard that exists because `read1` may be null.
pub const READ1_IS_NULL: &str =
    "Cannot invoke \"org.broadinstitute.hellbender.utils.read.GATKRead.getName()\" \
     because \"read1\" is null";

/// The JVM's message for `read2Reads.size()` when `groupingBy` produced no `false` key.
pub const READ2_READS_IS_NULL: &str =
    "Cannot invoke \"java.util.List.size()\" because \"read2Reads\" is null";

/// And for `read1Reads.size()` when it produced no `true` key.
///
/// Two messages rather than one, because `read1Reads.size() != read2Reads.size()` evaluates the
/// left operand first: which list is missing decides which name the JVM prints. Measured, not
/// guessed; the golden carries both.
pub const READ1_READS_IS_NULL: &str =
    "Cannot invoke \"java.util.List.size()\" because \"read1Reads\" is null";

/// `ReadPair`: one query-name group, split the way the reference splits it.
#[derive(Debug, Clone, Default)]
pub struct ReadPair {
    pub query_name: String,
    pub first_of_pair: Option<BamRecord>,
    pub second_of_pair: Option<BamRecord>,
    pub secondary_alignments: Vec<BamRecord>,
    /// Always empty on this tool's path, because the read filter removes every supplementary
    /// alignment before the traversal reaches `add`. Kept because the reference keeps it.
    pub supplementary_alignments: Vec<BamRecord>,
}

impl ReadPair {
    pub fn new(read: &BamRecord) -> Result<ReadPair, RsemError> {
        let mut pair = ReadPair {
            query_name: read.read_name.clone(),
            ..ReadPair::default()
        };
        pair.add(read)?;
        Ok(pair)
    }

    /// `ReadPair.add`, including the two `Utils.validate` calls and the final `else`.
    pub fn add(&mut self, read: &BamRecord) -> Result<(), RsemError> {
        if self.query_name != read.read_name {
            return Err(RsemError::NameMismatch {
                expected: self.query_name.clone(),
                found: read.read_name.clone(),
            });
        }
        let primary =
            !read::is_secondary_alignment(read) && !read::is_supplementary_alignment(read);
        if primary && read::is_first_of_pair(read) {
            if self.first_of_pair.is_some() {
                return Err(RsemError::PrimaryAlreadySet {
                    name: read.read_name.clone(),
                    first: true,
                });
            }
            self.first_of_pair = Some(read.clone());
        } else if primary && read::is_second_of_pair(read) {
            if self.second_of_pair.is_some() {
                return Err(RsemError::PrimaryAlreadySet {
                    name: read.read_name.clone(),
                    first: false,
                });
            }
            self.second_of_pair = Some(read.clone());
        } else if read::is_secondary_alignment(read) {
            self.secondary_alignments.push(read.clone());
        } else if read::is_supplementary_alignment(read) {
            self.supplementary_alignments.push(read.clone());
        } else {
            // A primary read that is neither first nor second of pair. The message reports the
            // contig rather than the name, which is the reference's choice and not a useful one.
            return Err(RsemError::UnknownReadType {
                contig: read.reference_index.to_string(),
            });
        }
        Ok(())
    }
}

/// `passesRSEMFilter`, including the null guard that dereferences null.
pub fn passes_rsem_filter(
    read1: Option<&BamRecord>,
    read2: Option<&BamRecord>,
) -> Result<bool, RsemError> {
    match (read1, read2) {
        // `read1.getName()` runs whichever of the two is null, so a null `read1` throws and a null
        // `read2` beside a present `read1` merely warns.
        (None, _) => return Err(RsemError::NullDereference(READ1_IS_NULL)),
        (Some(_), None) => return Ok(false),
        _ => {}
    }
    let (read1, read2) = (read1.expect("checked"), read2.expect("checked"));

    if read::is_unmapped(read1) || read::is_unmapped(read2) {
        return Ok(false);
    }
    // `contigsMatch`: the reference indexes, and an unmapped read has none. Both are mapped here.
    if read1.reference_index != read2.reference_index {
        return Ok(false);
    }
    if read1.cigar.num_elements() != 1 || read2.cigar.num_elements() != 1 {
        return Ok(false);
    }
    let single_m = |record: &BamRecord| record.cigar.elements[0].op == Op::M;
    Ok(single_m(read1) && single_m(read2))
}

/// `groupSecondaryReads`: the secondary alignments, paired up by their mate positions.
///
/// Returns `Err` for the shape the reference dereferences null on: every secondary a first-of-pair,
/// so `groupingBy` produced no `false` key.
pub fn group_secondary_reads(
    secondary: &[BamRecord],
) -> Result<Vec<(BamRecord, BamRecord)>, RsemError> {
    if secondary.is_empty() {
        return Ok(Vec::new());
    }
    let read1s: Vec<&BamRecord> = secondary
        .iter()
        .filter(|r| read::is_first_of_pair(r))
        .collect();
    let read2s: Vec<&BamRecord> = secondary
        .iter()
        .filter(|r| !read::is_first_of_pair(r))
        .collect();
    // `groupedByRead1.get(...)` is null when nothing grouped under that key, and the next line asks
    // both for their size. `read1Reads.size() != read2Reads.size()` evaluates the left operand
    // first, so an absent `true` key is reported before an absent `false` one.
    if read1s.is_empty() {
        return Err(RsemError::NullDereference(READ1_READS_IS_NULL));
    }
    if read2s.is_empty() {
        return Err(RsemError::NullDereference(READ2_READS_IS_NULL));
    }
    if read1s.len() != read2s.len() {
        // `logger.warn` then an empty list, so the whole secondary set is dropped rather than the
        // odd one out.
        return Ok(Vec::new());
    }

    let mut result = Vec::with_capacity(read1s.len());
    for read1 in &read1s {
        let mates: Vec<&&BamRecord> = read2s
            .iter()
            .filter(|r| {
                r.reference_index == read1.reference_index
                    && r.alignment_start == read1.mate_alignment_start
                    && r.mate_alignment_start == read1.alignment_start
            })
            .collect();
        // Exactly one mate is kept; none and more than one are both a warning and a skip.
        if mates.len() == 1 {
            result.push(((*read1).clone(), (**mates[0]).clone()));
        }
    }
    Ok(result)
}

/// `writeReads`: the records one query-name group contributes, in the order it contributes them.
pub fn write_reads(pair: &ReadPair) -> Result<Vec<BamRecord>, RsemError> {
    let mut out = Vec::new();
    if !passes_rsem_filter(pair.first_of_pair.as_ref(), pair.second_of_pair.as_ref())? {
        // The secondary loop is inside this `if`, so a failing primary takes them with it.
        return Ok(out);
    }
    out.push(pair.first_of_pair.clone().expect("checked by the filter"));
    out.push(pair.second_of_pair.clone().expect("checked by the filter"));

    for (read1, read2) in group_secondary_reads(&pair.secondary_alignments)? {
        if passes_rsem_filter(Some(&read1), Some(&read2))? {
            out.push(read1);
            out.push(read2);
        }
    }
    Ok(out)
}

/// What a run produces: the output BAM, and no index.
pub type RunResult = Result<Result<(Vec<u8>, Option<Vec<u8>>), RsemError>, ReadsError>;

/// `PostProcessReadsForRSEM`: the query-name groups that survive, reordered for RSEM.
pub fn post_process_reads_for_rsem(source: &ReadsDataSource, options: &Options) -> RunResult {
    if source.header().attributes.get("SO") != Some("queryname") {
        return Ok(Err(RsemError::NotQueryNameSorted));
    }

    let reads: Vec<BamRecord> = source
        .iter_all()?
        .into_iter()
        .filter(default_read_filter)
        .collect();

    let mut records = Vec::new();
    let mut current: Option<ReadPair> = None;
    for read in &reads {
        match &mut current {
            None => {
                current = match ReadPair::new(read) {
                    Ok(pair) => Some(pair),
                    Err(error) => return Ok(Err(error)),
                }
            }
            Some(pair) if pair.query_name != read.read_name => {
                match write_reads(pair) {
                    Ok(written) => records.extend(written),
                    Err(error) => return Ok(Err(error)),
                }
                current = match ReadPair::new(read) {
                    Ok(pair) => Some(pair),
                    Err(error) => return Ok(Err(error)),
                };
            }
            Some(pair) => {
                if let Err(error) = pair.add(read) {
                    return Ok(Err(error));
                }
            }
        }
    }
    if let Some(pair) = &current {
        match write_reads(pair) {
            Ok(written) => records.extend(written),
            Err(error) => return Ok(Err(error)),
        }
    }

    let header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    // `createSAMWriter(outSam, true)`: presorted, so nothing is re-ordered on the way out, and a
    // queryname header has no index.
    Ok(Ok(write_records(&header, &records, false)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::text_parse::parse_cigar;

    fn read(name: &str, flags: u16, contig: i32, start: i32, cigar: &str, mate: i32) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            flags,
            reference_index: contig,
            alignment_start: start,
            mate_reference_index: contig,
            mate_alignment_start: mate,
            cigar: parse_cigar(cigar).unwrap(),
            mapping_quality: 60,
            ..BamRecord::default()
        }
    }

    const PAIRED: u16 = 0x1;
    const UNMAPPED: u16 = 0x4;
    const FIRST: u16 = 0x40;
    const SECOND: u16 = 0x80;
    const SECONDARY: u16 = 0x100;
    const SUPPLEMENTARY: u16 = 0x800;

    #[test]
    fn a_missing_first_of_pair_throws_out_of_the_guard_that_checks_for_it() {
        let second = read("q1", PAIRED | SECOND, 0, 100, "100M", 300);
        assert_eq!(
            passes_rsem_filter(None, Some(&second)),
            Err(RsemError::NullDereference(READ1_IS_NULL))
        );
        // The other way round only warns, because `read1.getName()` succeeds.
        assert_eq!(passes_rsem_filter(Some(&second), None), Ok(false));
    }

    #[test]
    fn one_sided_secondary_alignments_throw_too() {
        let only_first = vec![read("r1", PAIRED | FIRST | SECONDARY, 0, 500, "100M", 700)];
        assert_eq!(
            group_secondary_reads(&only_first),
            Err(RsemError::NullDereference(READ2_READS_IS_NULL))
        );
        // The mirror image names the other list, because the comparison evaluates left first.
        let only_second = vec![read("t1", PAIRED | SECOND | SECONDARY, 0, 700, "100M", 500)];
        assert_eq!(
            group_secondary_reads(&only_second),
            Err(RsemError::NullDereference(READ1_READS_IS_NULL))
        );
        assert_eq!(group_secondary_reads(&[]), Ok(Vec::new()), "empty is fine");
    }

    #[test]
    fn a_single_element_cigar_that_is_not_m_is_refused() {
        let first = read("p5", PAIRED | FIRST, 0, 1700, "100=", 1900);
        let second = read("p5", PAIRED | SECOND, 0, 1900, "100M", 1700);
        assert_eq!(passes_rsem_filter(Some(&first), Some(&second)), Ok(false));

        let ok_first = read("p1", PAIRED | FIRST, 0, 100, "100M", 300);
        let ok_second = read("p1", PAIRED | SECOND, 0, 300, "100M", 100);
        assert_eq!(
            passes_rsem_filter(Some(&ok_first), Some(&ok_second)),
            Ok(true)
        );
    }

    #[test]
    fn unmapped_and_chimeric_pairs_are_refused() {
        let mapped = read("p2", PAIRED | FIRST, 0, 900, "100M", 900);
        let mut unmapped = read("p2", PAIRED | SECOND | UNMAPPED, 0, 900, "100M", 900);
        unmapped.cigar = parse_cigar("").unwrap();
        assert_eq!(
            passes_rsem_filter(Some(&mapped), Some(&unmapped)),
            Ok(false)
        );

        let here = read("p3", PAIRED | FIRST, 0, 1100, "100M", 100);
        let there = read("p3", PAIRED | SECOND, 1, 100, "100M", 1100);
        assert_eq!(passes_rsem_filter(Some(&here), Some(&there)), Ok(false));
    }

    #[test]
    fn a_secondary_pair_is_matched_by_its_mate_positions() {
        let secondary = vec![
            read("p1", PAIRED | FIRST | SECONDARY, 0, 500, "100M", 700),
            read("p1", PAIRED | SECOND | SECONDARY, 0, 700, "100M", 500),
        ];
        let paired = group_secondary_reads(&secondary).expect("both sides present");
        assert_eq!(paired.len(), 1);
        assert_eq!(paired[0].0.alignment_start, 500);
        assert_eq!(paired[0].1.alignment_start, 700);
    }

    #[test]
    fn the_output_order_is_primary_then_each_secondary() {
        let mut pair = ReadPair::new(&read("p1", PAIRED | FIRST, 0, 100, "100M", 300)).unwrap();
        pair.add(&read("p1", PAIRED | SECOND, 0, 300, "100M", 100))
            .unwrap();
        pair.add(&read("p1", PAIRED | FIRST | SECONDARY, 0, 500, "100M", 700))
            .unwrap();
        pair.add(&read(
            "p1",
            PAIRED | SECOND | SECONDARY,
            0,
            700,
            "100M",
            500,
        ))
        .unwrap();
        let written = write_reads(&pair).expect("no null dereference");
        assert_eq!(
            written
                .iter()
                .map(|r| r.alignment_start)
                .collect::<Vec<_>>(),
            vec![100, 300, 500, 700]
        );
    }

    #[test]
    fn a_failing_primary_takes_its_secondary_alignments_with_it() {
        // The primary is chimeric, so nothing at all is written.
        let mut pair = ReadPair::new(&read("p3", PAIRED | FIRST, 0, 1100, "100M", 100)).unwrap();
        pair.add(&read("p3", PAIRED | SECOND, 1, 100, "100M", 1100))
            .unwrap();
        pair.add(&read("p3", PAIRED | FIRST | SECONDARY, 0, 500, "100M", 700))
            .unwrap();
        pair.add(&read(
            "p3",
            PAIRED | SECOND | SECONDARY,
            0,
            700,
            "100M",
            500,
        ))
        .unwrap();
        assert!(write_reads(&pair).unwrap().is_empty());
    }

    #[test]
    fn the_filter_is_one_filter_and_it_is_not_wellformed() {
        assert_eq!(DEFAULT_READ_FILTERS.len(), 1);
        let supplementary = read("p6", PAIRED | FIRST | SUPPLEMENTARY, 0, 2500, "100M", 2300);
        assert!(!default_read_filter(&supplementary));
        // No read group and no mate, which Wellformed would reject and this filter does not.
        assert!(default_read_filter(&read(
            "p6",
            PAIRED | FIRST,
            0,
            2100,
            "100M",
            2300
        )));
    }
}
