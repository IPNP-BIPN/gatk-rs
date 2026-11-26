//! `RampedHaplotypeCaller`'s ramps: the state file an off ramp writes and an on ramp restarts from.
//!
//! The tool is `HaplotypeCaller` with a different engine, and the engine is milestone G3's. What is
//! the tool's own is the ramp file's shape, the two orderings its contents are compared under, and
//! what the comparison refuses. Those are ported; the caller is not.
//!
//! The Java is `RampBase`, `OffRampBase` and `RampUtils` under
//! `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps` (GATK 4.7.0.0).
//!
//! # The comparators return a sign, the Java returns a difference
//!
//! `RampUtils`' comparators chain `String.compareTo`, which returns the difference of the first
//! two unequal characters, not -1 or 1. The port maps every string comparison to -1, 0 or 1. Both
//! callers here only sort with the result and test it against zero, so only its sign is
//! observable; a caller that read the magnitude would see a different number.
//!
//! # Reading this file without reading Rust
//!
//! - `&str` is a borrowed, read-only view of a string, `&[T]` the same for an array; the caller
//!   keeps ownership. `String` and `Vec<T>` are the owned, growable versions;
//! - `Option<&Interval>` is "an interval or nothing", Java's nullable reference checked by the
//!   compiler, and `match` on it must handle both cases;
//! - `Result<(), E>` is "nothing, or an error `E`": how this codebase models the Java's
//!   `throw new RuntimeException(...)`;
//! - `a.cmp(&b)` returns an `Ordering` (`Less`, `Equal`, `Greater`) rather than an `int`; on `str`
//!   it compares bytes, which for the ASCII these fields hold is the order `String.compareTo`
//!   gives;
//! - `i32::from(flag)` turns a `bool` into 0 or 1, the Java's `flag ? 1 : 0`;
//! - a function's last expression is its return value; `return` is only written for early exits.

/// `RampBase.Type`, whose constants are `OffRamp` and `OnRamp`.
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

/// `RampBase.getLocFilenameSuffix`: a directory name made out of coordinates.
///
/// ```java
/// return String.format("%s-%d-%d", loc.getContig(), loc.getStart(), loc.getEnd());
/// ```
pub fn loc_filename_suffix(loc: &Interval) -> String {
    format!("{}-{}-{}", loc.contig, loc.start, loc.end)
}

/// `OffRampBase.addEntry`'s name: the region's suffix as a directory, or the bare name at the
/// root. The on ramp's `getEntry` builds the name it looks up the same way.
///
/// ```java
/// final String      prefix = loc != null ? getLocFilenameSuffix(loc) + "/" : "";
/// ```
pub fn entry_name(loc: Option<&Interval>, name: &str) -> String {
    match loc {
        Some(loc) => format!("{}/{}", loc_filename_suffix(loc), name),
        None => name.to_string(),
    }
}

/// `RampBase.getReadSuppName`: the read's name with a supplementary flag appended, so a primary
/// and a supplementary record of one read are different keys.
///
/// ```java
/// return name + (supp ? ",1" : ",0");
/// ```
pub fn read_supp_name(name: &str, supplementary: bool) -> String {
    format!("{name},{}", if supplementary { "1" } else { "0" })
}

/// `getBamIndexPath`, which is `String.replace` rather than a suffix change: EVERY `.bam` in the
/// path is replaced, including one inside a directory name, and a path holding none comes back
/// unchanged.
///
/// ```java
/// return (new File(bamPath.toFile().getAbsolutePath().replace(".bam", ".bai"))).toPath();
/// ```
///
/// The Java makes the path absolute first, so a `.bam` in the working directory's own path is
/// replaced too; this function replaces only in the text it is given.
pub fn bam_index_path(path: &str) -> String {
    path.replace(".bam", ".bai")
}

/// `RampUtils.readInfo`: the supplementary flag, then the record's own text. The Java logs it for
/// each side of a read pair that fails verification.
///
/// ```java
/// final String      prefix = read.isSupplementaryAlignment() ? "1," : "0,";
/// ```
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

/// The header `OffRampBase.addHaplotypes` writes above the rows.
pub const HAPLOTYPE_HEADER: &str = "contig,start,end,ref,cigar,bases,score,alignmentStartHapwrtRef";

/// `OffRampBase.addHaplotypes`' table. The reference column is 1 or 0 and the score is
/// `Double.toString`, which is why it goes through `java_double_to_string` rather than Rust's own
/// float formatting (the two differ, for instance `1.0` against `1`).
///
/// ```java
/// pw.println(String.format("%s,%d,%d,%d,%s,%s,%s,%d",
/// ```
///
/// `println` ends each row with the line separator, `\n` on the Linux the reference runs on.
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
///
/// ```java
/// double      ddelta = o1.getScore() - o2.getScore();
/// if ( ddelta < 0 )
///     delta = -1;
/// else if ( ddelta > 0 )
///     delta = 1;
/// ```
///
/// Order of keys: reference flag, contig, start, end, score, bases. The CIGAR is not a key, so two
/// haplotypes differing only in CIGAR compare equal.
pub fn compare_haplotypes(a: &Haplotype, b: &Haplotype) -> i32 {
    let delta = i32::from(a.is_reference) - i32::from(b.is_reference);
    if delta != 0 {
        return delta;
    }
    // `String.compareTo`, reduced to its sign (see the module header).
    let delta = match a.contig.cmp(&b.contig) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
    };
    if delta != 0 {
        return delta;
    }
    // Integer subtraction, as the Java's `o1.getStart() - o2.getStart()`. Positions are far from
    // `i32`'s range, so the overflow Java would wrap through cannot arise.
    let delta = a.start - b.start;
    if delta != 0 {
        return delta;
    }
    let delta = a.end - b.end;
    if delta != 0 {
        return delta;
    }
    // A NaN difference is neither below nor above zero, so it falls to the bases, as in the Java.
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
    /// The raw quality bytes, compared as a string. The Java decodes them with
    /// `new String(o1.getBaseQualitiesNoCopy())`; for qualities below 128, which is every Phred
    /// score a SAM can carry, that decoding is one character per byte and byte order is the same.
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
///
/// ```java
/// int     delta = (o1.isReverseStrand() ? 1 : 0) - (o2.isReverseStrand() ? 1 : 0);
/// ```
///
/// Then, in order: `commonToString`, bases, qualities, soft start, soft end, start, end,
/// unclipped start, unclipped end. The Java nests one `if ( delta == 0 )` per key; the port lists
/// the keys in an array and returns the first that is not zero.
pub fn compare_reads(a: &Read, b: &Read) -> i32 {
    let delta = i32::from(a.reverse_strand) - i32::from(b.reverse_strand);
    if delta != 0 {
        return delta;
    }
    // A closure (an inline function bound to a name, Java's lambda) reducing `Ordering` to a sign.
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

/// What the verification refuses: the four `RuntimeException`s of `RampUtils.compareHaplotypes`
/// and `compareReads`, `message` giving each one's text.
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
///
/// `LinkedList.sort` is stable and so is Rust's `sort_by`, so equal elements keep their order on
/// both sides.
pub fn verify_haplotypes(left: &[Haplotype], right: &[Haplotype]) -> Result<(), VerificationError> {
    if left.len() != right.len() {
        return Err(VerificationError::HaplotypeSize {
            left: left.len(),
            right: right.len(),
        });
    }
    // `to_vec` copies the borrowed slice, the Java's `new LinkedList<>(c1)`. `sort_by` wants an
    // `Ordering`, so the comparator's sign is turned into one by comparing it with zero.
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

/// `RampUtils.compareReads`, the same shape. The Java also logs `readInfo` of both reads before
/// throwing; the port returns only the error.
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
