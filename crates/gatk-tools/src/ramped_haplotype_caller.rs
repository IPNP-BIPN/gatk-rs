//! `RampedHaplotypeCaller`'s ramps: the state file an off ramp writes and an on ramp restarts from.
//!
//! The tool is `HaplotypeCaller` with a different engine, and the engine is milestone G3's. What is
//! the tool's own is the ramp file's shape, the two orderings its contents are compared under, and
//! what the comparison refuses. Those are ported; the caller is not.

/// `RampBase.Type`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RampType {
    Off,
    On,
}

/// A genomic interval, as the ramp names entries from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

/// `getLocFilenameSuffix`: a directory name made out of coordinates.
pub fn loc_filename_suffix(loc: &Interval) -> String {
    format!("{}-{}-{}", loc.contig, loc.start, loc.end)
}

/// `addEntry`'s name: the region's suffix as a directory, or the bare name at the root.
pub fn entry_name(loc: Option<&Interval>, name: &str) -> String {
    match loc {
        Some(loc) => format!("{}/{}", loc_filename_suffix(loc), name),
        None => name.to_string(),
    }
}

/// `getReadSuppName`: the read's name with a supplementary flag appended.
pub fn read_supp_name(name: &str, supplementary: bool) -> String {
    format!("{name},{}", if supplementary { "1" } else { "0" })
}

/// `getBamIndexPath`, which is `String.replace` rather than a suffix change: EVERY `.bam` in the
/// path is replaced, including one inside a directory name, and a path holding none comes back
/// unchanged.
pub fn bam_index_path(path: &str) -> String {
    path.replace(".bam", ".bai")
}

/// `readInfo`: the supplementary flag, then the record's own text.
pub fn read_info(supplementary: bool, record_text: &str) -> String {
    format!("{},{record_text}", if supplementary { 1 } else { 0 })
}

/// One haplotype, reduced to what the table and the comparator read.
#[derive(Debug, Clone, PartialEq)]
pub struct Haplotype {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub is_reference: bool,
    pub cigar: String,
    pub bases: String,
    pub score: f64,
    pub alignment_start_hap_wrt_ref: i32,
}

/// The header `addHaplotypes` writes above the rows.
pub const HAPLOTYPE_HEADER: &str = "contig,start,end,ref,cigar,bases,score,alignmentStartHapwrtRef";

/// `addHaplotypes`' table. The reference column is 1 or 0 and the score is `Double.toString`.
pub fn haplotype_table(haplotypes: &[Haplotype]) -> String {
    let mut out = String::from(HAPLOTYPE_HEADER);
    out.push('\n');
    for haplotype in haplotypes {
        out.push_str(&format!(
            "{},{},{},{},{},{},{},{}\n",
            haplotype.contig,
            haplotype.start,
            haplotype.end,
            if haplotype.is_reference { 1 } else { 0 },
            haplotype.cigar,
            haplotype.bases,
            gatk_engine::tsv_table::java_double_to_string(haplotype.score),
            haplotype.alignment_start_hap_wrt_ref
        ));
    }
    out
}

/// `RampUtils.HaplotypeComparator`.
///
/// A REFERENCE haplotype sorts LAST, because the first key is the difference of the reference
/// flags and a true is one. The score is compared by SIGN rather than by difference, so two scores
/// a hair apart are still ordered.
pub fn compare_haplotypes(a: &Haplotype, b: &Haplotype) -> i32 {
    let delta = i32::from(a.is_reference) - i32::from(b.is_reference);
    if delta != 0 {
        return delta;
    }
    let delta = match a.contig.cmp(&b.contig) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
    };
    if delta != 0 {
        return delta;
    }
    let delta = a.start - b.start;
    if delta != 0 {
        return delta;
    }
    let delta = a.end - b.end;
    if delta != 0 {
        return delta;
    }
    let difference = a.score - b.score;
    let delta = if difference < 0.0 {
        -1
    } else if difference > 0.0 {
        1
    } else {
        0
    };
    if delta != 0 {
        return delta;
    }
    match a.bases.cmp(&b.bases) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
    }
}

/// One read, reduced to what the comparator reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub reverse_strand: bool,
    /// `commonToString`, which the comparator uses before the bases.
    pub common_text: String,
    pub bases: String,
    /// The raw quality bytes, compared as a string.
    pub base_qualities: Vec<u8>,
    pub soft_start: i32,
    pub soft_end: i32,
    pub start: i32,
    pub end: i32,
    pub unclipped_start: i32,
    pub unclipped_end: i32,
}

/// `RampUtils.GATKReadComparator`, which starts from the STRAND: a reverse read sorts after a
/// forward one whatever their positions.
pub fn compare_reads(a: &Read, b: &Read) -> i32 {
    let delta = i32::from(a.reverse_strand) - i32::from(b.reverse_strand);
    if delta != 0 {
        return delta;
    }
    let text = |ordering: std::cmp::Ordering| match ordering {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
    };
    for delta in [
        text(a.common_text.cmp(&b.common_text)),
        text(a.bases.cmp(&b.bases)),
        text(a.base_qualities.cmp(&b.base_qualities)),
        a.soft_start - b.soft_start,
        a.soft_end - b.soft_end,
        a.start - b.start,
        a.end - b.end,
        a.unclipped_start - b.unclipped_start,
        a.unclipped_end - b.unclipped_end,
    ] {
        if delta != 0 {
            return delta;
        }
    }
    0
}

/// What the verification refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VerificationError {
    HaplotypeSize { left: usize, right: usize },
    HaplotypeIndex { index: usize },
    ReadSize { left: usize, right: usize },
    ReadIndex { index: usize },
}

impl VerificationError {
    pub fn message(&self) -> String {
        match self {
            VerificationError::HaplotypeSize { left, right } => {
                format!("haplotype size verification failed: {left} vs {right}")
            }
            VerificationError::HaplotypeIndex { index } => {
                format!("haplotype failed verification on index {index}")
            }
            VerificationError::ReadSize { left, right } => {
                format!("reads size verification failed: {left} vs {right}")
            }
            VerificationError::ReadIndex { index } => {
                format!("reads failed verification on index {index}")
            }
        }
    }
}

/// `compareHaplotypes`: the SIZE is checked before anything is compared, and each side is sorted
/// before the pairwise walk, so the order the two collections arrived in does not matter.
pub fn verify_haplotypes(left: &[Haplotype], right: &[Haplotype]) -> Result<(), VerificationError> {
    if left.len() != right.len() {
        return Err(VerificationError::HaplotypeSize {
            left: left.len(),
            right: right.len(),
        });
    }
    let mut sorted_left = left.to_vec();
    let mut sorted_right = right.to_vec();
    sorted_left.sort_by(|a, b| compare_haplotypes(a, b).cmp(&0));
    sorted_right.sort_by(|a, b| compare_haplotypes(a, b).cmp(&0));
    for index in 0..sorted_left.len() {
        if compare_haplotypes(&sorted_left[index], &sorted_right[index]) != 0 {
            return Err(VerificationError::HaplotypeIndex { index });
        }
    }
    Ok(())
}

/// `compareReads`, the same shape.
pub fn verify_reads(left: &[Read], right: &[Read]) -> Result<(), VerificationError> {
    if left.len() != right.len() {
        return Err(VerificationError::ReadSize {
            left: left.len(),
            right: right.len(),
        });
    }
    let mut sorted_left = left.to_vec();
    let mut sorted_right = right.to_vec();
    sorted_left.sort_by(|a, b| compare_reads(a, b).cmp(&0));
    sorted_right.sort_by(|a, b| compare_reads(a, b).cmp(&0));
    for index in 0..sorted_left.len() {
        if compare_reads(&sorted_left[index], &sorted_right[index]) != 0 {
            return Err(VerificationError::ReadIndex { index });
        }
    }
    Ok(())
}
