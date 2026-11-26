//! `CollectSVEvidence`: what one BAM contributes to structural-variant calling.
//!
//! One traversal writes four evidence files, and each has its own rule for which reads it will
//! look at: discordant pairs, split-read positions, per-site allele depths and per-interval read
//! counts.
//!
//! Reading the BAM, the VCF and the interval file is not ported, nor are the streaming buffers
//! that keep each writer's memory bounded. Which read contributes what, and what each contribution
//! is, are.
//!
//! The reference is `org.broadinstitute.hellbender.tools.walkers.sv.CollectSVEvidence` (GATK
//! 4.7.0.0; `CollectSVEvidence.java` is unchanged since 4.6.2.0), with the encoders of
//! `DiscordantPairEvidenceCodec` and `SplitReadEvidenceCodec` for the zero-based positions.
//!
//! # Two output orders that are not the reference's
//!
//! - Split reads: the Java flushes its `PriorityQueue<SplitPos>` whenever the contig changes, so
//!   the file follows the traversal's contig order. [`split_reads`] sorts the whole run by contig
//!   NAME, so a BAM with split reads on `chr2` and then `chr10` writes `chr10` first here.
//! - Discordant pairs: at every new start the Java sorts the buffered reads with
//!   `DiscordantReadComparator` (contig index, start, read strand, mate contig index, mate start,
//!   mate strand) before writing them. [`discordant_pairs`] keeps traversal order, so two pairs
//!   at one start whose mates are at 500 and then 300 come out in the other order.
//!
//! # Reading this file without reading Rust
//!
//! - `&[Read]` is a borrowed, read-only view of a list; the caller keeps ownership.
//! - `Option<T>` is "a value or nothing", Java's nullable reference with the check enforced;
//!   `let Some(x) = y else { .. };` unwraps it or takes the early exit.
//! - `match` must cover every case; `Some(('M', _))` matches a first cigar element whose operator
//!   is `M`, whatever its length.
//! - `(i32, Side)` is a tuple, two values returned together without a holder class.
//! - `.iter().filter(..).map(..).sum()` is a Java stream ending in `sum()`; `.count() as i32`
//!   narrows the count to an `int`.
//! - `x.then_some(v)` is `x ? Optional.of(v) : Optional.empty()`.

/// One read, reduced to what the four writers read off it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    /// `getName()`.
    pub name: String,
    /// The contig's index in the sequence dictionary, `getSequenceIndex(read.getContig())`.
    pub contig_index: usize,
    /// `getContig()`.
    pub contig: String,
    /// `getStart()`, one-based.
    pub start: i32,
    /// `getMappingQuality()`.
    pub mapping_quality: i32,
    /// Cigar as (operator, length), the operator being one of `M`, `S`, `D`, `I`, `N`, `H`.
    pub cigar: Vec<(char, i32)>,
    // The flags below are the `GATKRead` accessors of the same names (`isPaired()`,
    // `isProperlyPaired()`, `mateIsUnmapped()` and so on). The mate's fields are `None` where the
    // read has no mate position to report.
    pub paired: bool,
    pub properly_paired: bool,
    pub mate_unmapped: bool,
    pub mate_contig_index: Option<usize>,
    pub mate_contig: Option<String>,
    pub mate_start: Option<i32>,
    pub reverse_strand: bool,
    pub mate_reverse_strand: bool,
    pub supplementary: bool,
    pub secondary: bool,
    pub duplicate: bool,
    pub unmapped: bool,
    /// The bases and their qualities, for the site-depth counter.
    pub bases: Vec<u8>,
    /// One quality per base.
    pub base_qualities: Vec<i32>,
}

/// The two default read filters this tool adds on top of the walker's own:
/// `ReadFilterLibrary.MAPPED` and `ReadFilterLibrary.NOT_DUPLICATE`.
pub fn passes_default_filters(read: &Read) -> bool {
    !read.unmapped && !read.duplicate
}

/// The outer test in `apply`: a read whose MATE is unmapped, a supplementary alignment and a
/// secondary alignment are all skipped by both the split and the discordant writer, and by neither
/// of the two depth counters.
///
/// ```java
/// if ( !(read.isPaired() && read.mateIsUnmapped()) &&
///         !read.isSupplementaryAlignment() &&
///         !read.isSecondaryAlignment() ) {
/// ```
pub fn contributes_to_pair_evidence(read: &Read) -> bool {
    !(read.paired && read.mate_unmapped) && !read.supplementary && !read.secondary
}

/// `isSoftClipped`: EXACTLY one end. A read clipped at both ends is not a split read at all.
///
/// ```java
/// return (firstOperator == CigarOperator.SOFT_CLIP && lastOperator != CigarOperator.SOFT_CLIP) ||
///         (firstOperator != CigarOperator.SOFT_CLIP && lastOperator == CigarOperator.SOFT_CLIP);
/// ```
///
/// which is "exactly one of the two" and is written here as `!=` on the two tests.
pub fn is_soft_clipped(read: &Read) -> bool {
    let Some((first, _)) = read.cigar.first() else {
        return false;
    };
    let Some((last, _)) = read.cigar.last() else {
        return false;
    };
    (*first == 'S') != (*last == 'S')
}

/// Which side of the split a position marks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Side {
    /// `POSITION.LEFT`: the clip is on the left, the position is the read's start.
    Left,
    /// `POSITION.RIGHT`: the clip is on the right, the position is past the alignment's end.
    Right,
    /// Neither end is a match, so the read gives no position.
    Middle,
}

impl Side {
    /// The text the codec writes.
    pub fn name(self) -> &'static str {
        match self {
            Side::Left => "left",
            Side::Right => "right",
            Side::Middle => "middle",
        }
    }
}

/// `CigarOperator.consumesReferenceBases()`.
fn consumes_reference(operator: char) -> bool {
    matches!(operator, 'M' | 'D' | 'N' | '=' | 'X')
}

/// `getSplitPosition`.
///
/// The FIRST element being a match wins, whatever the last one is, and the position it gives is the
/// start plus every reference-consuming length: a deletion inside the alignment moves it out by its
/// own length. Only when the first is not a match is the last one asked.
///
/// ```java
/// if (read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.M) {
/// ```
///
/// Only `M`: a cigar written with `=` or `X` is neither side.
pub fn split_position(read: &Read) -> (i32, Side) {
    match read.cigar.first() {
        Some(('M', _)) => {
            let match_length: i32 = read
                .cigar
                .iter()
                .filter(|(operator, _)| consumes_reference(*operator))
                .map(|(_, length)| *length)
                .sum();
            (read.start + match_length, Side::Right)
        }
        _ => match read.cigar.last() {
            Some(('M', _)) => (read.start, Side::Left),
            _ => (-1, Side::Middle),
        },
    }
}

/// One split-read record: a position, a side and how many reads were clipped there.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplitRead {
    /// The read's contig.
    pub contig: String,
    /// The codec writes it zero-based.
    pub position: i32,
    /// Which side of the clip the position marks.
    pub side: Side,
    /// How many reads gave this position and side.
    pub count: i32,
}

/// Every split-read record one BAM produces, in output order.
///
/// Identical (position, side) pairs are counted together; the SAME position reached from the two
/// sides stays two records, because the side is part of what is compared.
pub fn split_reads(reads: &[Read]) -> Vec<SplitRead> {
    let mut positions: Vec<(String, i32, Side)> = Vec::new();
    for read in reads {
        if !passes_default_filters(read) || !contributes_to_pair_evidence(read) {
            continue;
        }
        if !is_soft_clipped(read) {
            continue;
        }
        let (position, side) = split_position(read);
        if side == Side::Middle {
            continue;
        }
        positions.push((read.contig.clone(), position, side));
    }
    // The buffer is a priority queue over position then side, flushed in that order. The contig
    // in front of them is compared by name, which is not the Java's order: see the module note.
    positions.sort_by(|a, b| (&a.0, a.1, a.2).cmp(&(&b.0, b.1, b.2)));
    let mut out: Vec<SplitRead> = Vec::new();
    for (contig, position, side) in positions {
        // `last_mut` borrows the last record writably so its count can be bumped in place; the
        // `if` after the pattern is a guard, so a different position falls to `_` and is pushed.
        match out.last_mut() {
            Some(last)
                if last.contig == contig && last.position == position - 1 && last.side == side =>
            {
                last.count += 1;
            }
            _ => out.push(SplitRead {
                contig,
                // The codec writes zero-based positions.
                position: position - 1,
                side,
                count: 1,
            }),
        }
    }
    out
}

/// One discordant-pair record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscordantPair {
    /// The read's contig.
    pub contig: String,
    /// Zero-based, as the codec writes it.
    pub position: i32,
    /// `!isReverseStrand`. The codec prints `+` for true, so the two negations cancel and the
    /// output shows the read's own strand.
    pub strand: bool,
    /// The mate's contig.
    pub mate_contig: String,
    /// The mate's start, zero-based.
    pub mate_position: i32,
    /// `!isMateReverseStrand`, as for `strand`.
    pub mate_strand: bool,
}

/// `getReportableDiscordantReadPair`: a pair is written ONCE, by whichever of its two reads has the
/// smaller contig index, or the smaller start, and at an EQUAL start by the first of the two SEEN.
///
/// The last case is tracked by NAME within the locus, so two different pairs at one position are
/// still two records: the name set only ever suppresses a read's own mate.
///
/// ```java
/// final boolean seenBefore = observedDiscordantNamesAtThisLocus.remove(read.getName());
/// ```
///
/// Records are returned in traversal order, without the Java's per-start sort; see the module
/// note.
pub fn discordant_pairs(reads: &[Read]) -> Vec<DiscordantPair> {
    let mut out = Vec::new();
    let mut seen_at_locus: Vec<String> = Vec::new();
    let mut current_position: i32 = -1;
    for read in reads {
        if !passes_default_filters(read) || !contributes_to_pair_evidence(read) {
            continue;
        }
        // An UNPAIRED read reports properly_paired as false, and the writer then asks it for a
        // mate it does not have: see `crashes_on_an_unpaired_read`.
        if read.properly_paired {
            continue;
        }
        if read.start != current_position {
            current_position = read.start;
            seen_at_locus.clear();
        }
        let Some(mate_index) = read.mate_contig_index else {
            continue;
        };
        let mate_start = read.mate_start.unwrap_or(0);
        let reportable = if read.contig_index < mate_index {
            true
        } else if read.contig_index == mate_index {
            if read.start < mate_start {
                true
            } else if read.start == mate_start {
                match seen_at_locus.iter().position(|name| *name == read.name) {
                    // Seen before at this locus: removed and dropped.
                    Some(at) => {
                        seen_at_locus.remove(at);
                        false
                    }
                    None => {
                        seen_at_locus.push(read.name.clone());
                        true
                    }
                }
            } else {
                false
            }
        } else {
            false
        };
        if reportable {
            out.push(DiscordantPair {
                contig: read.contig.clone(),
                position: read.start - 1,
                strand: !read.reverse_strand,
                mate_contig: read.mate_contig.clone().unwrap_or_default(),
                mate_position: mate_start - 1,
                mate_strand: !read.mate_reverse_strand,
            });
        }
    }
    out
}

/// Asking an unpaired read for its mate is what the discordant writer does, and it throws.
pub const UNPAIRED_MATE_MESSAGE: &str = "Cannot get mate information for an unpaired read";

/// True when the discordant writer would ask this read for a mate it does not have.
pub fn crashes_on_an_unpaired_read(read: &Read) -> bool {
    passes_default_filters(read) && contributes_to_pair_evidence(read) && !read.paired
}

/// One candidate site-depth locus.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    /// The VCF record's contig.
    pub contig: String,
    /// Its one-based position.
    pub position: i32,
    /// The reference allele's bases.
    pub reference: String,
    /// The alternate alleles' bases, in record order.
    pub alternates: Vec<String>,
}

impl Site {
    /// `isSNP`, which here also requires the site to be biallelic.
    pub fn is_biallelic_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternates.len() == 1 && self.alternates[0].len() == 1
    }
}

/// `BAFSiteIterator`: biallelic SNPs at NEW loci only, so a repeated position, a triallelic site
/// and an indel are all walked past.
pub fn baf_sites(sites: &[Site]) -> Vec<Site> {
    let mut out: Vec<Site> = Vec::new();
    for site in sites {
        if !site.is_biallelic_snp() {
            continue;
        }
        let new_locus = match out.last() {
            None => true,
            Some(last) => last.contig != site.contig || last.position < site.position,
        };
        if new_locus {
            out.push(site.clone());
        }
    }
    out
}

/// `MIN_SITE_DEPTH_MAPQ`, `MIN_SITE_DEPTH_BASEQ` and `MIN_DEPTH_EVIDENCE_MAPQ`, whose three
/// different defaults are what let one read be counted for depth and not for site depth.
///
/// `public int minMapQ = 30;`
pub const DEFAULT_SITE_DEPTH_MIN_MAPQ: i32 = 30;
/// `public int minQ = 20;`
pub const DEFAULT_SITE_DEPTH_MIN_BASEQ: i32 = 20;
/// `public int minDepthEvidenceMapQ = 0;`
pub const DEFAULT_DEPTH_EVIDENCE_MIN_MAPQ: i32 = 0;

/// One site-depth record: the four base counts, in `A`, `C`, `G`, `T` order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteDepth {
    /// The site's contig.
    pub contig: String,
    /// Zero-based, as the codec writes it.
    pub position: i32,
    /// Reads carrying `A`, `C`, `G`, `T` at the site.
    pub counts: [i32; 4],
}

/// The column a base is counted in, either case; any other base counts nowhere.
fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// The site depths one BAM produces over the given loci.
///
/// The mapping-quality floor removes a whole READ; the base-quality floor removes a single BASE.
/// A locus no read covers is still written, with four zeros.
pub fn site_depths(
    reads: &[Read],
    sites: &[Site],
    min_mapping_quality: i32,
    min_base_quality: i32,
) -> Vec<SiteDepth> {
    baf_sites(sites)
        .into_iter()
        .map(|site| {
            let mut counts = [0; 4];
            for read in reads {
                if !passes_default_filters(read) || read.contig != site.contig {
                    continue;
                }
                if read.mapping_quality < min_mapping_quality {
                    continue;
                }
                // Only the aligned matches are walked, so the offset is into the leading match.
                let Some(offset) = reference_offset(read, site.position) else {
                    continue;
                };
                if read.base_qualities.get(offset).copied().unwrap_or(0) < min_base_quality {
                    continue;
                }
                if let Some(index) = read.bases.get(offset).copied().and_then(base_index) {
                    counts[index] += 1;
                }
            }
            SiteDepth {
                contig: site.contig,
                position: site.position - 1,
                counts,
            }
        })
        .collect()
}

/// Where a reference position falls in a read's bases, walking its cigar. `None` when the position
/// is outside the read or falls in a deletion.
pub fn reference_offset(read: &Read, position: i32) -> Option<usize> {
    let mut reference = read.start;
    let mut query = 0usize;
    for (operator, length) in &read.cigar {
        match operator {
            'M' | '=' | 'X' => {
                if position >= reference && position < reference + length {
                    return Some(query + (position - reference) as usize);
                }
                reference += length;
                query += *length as usize;
            }
            'D' | 'N' => {
                if position >= reference && position < reference + length {
                    return None;
                }
                reference += length;
            }
            'I' | 'S' => query += *length as usize,
            _ => {}
        }
    }
    None
}

/// One depth-evidence interval and the reads that started inside it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DepthEvidence {
    /// The interval's contig.
    pub contig: String,
    /// The interval file's own coordinates, written back unchanged.
    pub start: i32,
    /// The interval's end, as written in the interval file.
    pub end: i32,
    /// Reads whose start falls in the interval.
    pub count: i32,
}

/// The read counts one BAM produces over the given intervals.
///
/// A read is counted for the interval its START falls in, and an interval no read reaches is still
/// written with a count of zero.
pub fn depth_evidence(
    reads: &[Read],
    intervals: &[(String, i32, i32)],
    min_mapping_quality: i32,
) -> Vec<DepthEvidence> {
    intervals
        .iter()
        .map(|(contig, start, end)| DepthEvidence {
            contig: contig.clone(),
            start: *start,
            end: *end,
            count: reads
                .iter()
                .filter(|read| {
                    passes_default_filters(read)
                        && read.contig == *contig
                        && read.mapping_quality >= min_mapping_quality
                        // The BED start is exclusive of the first base, so the interval covers
                        // start + 1 through end.
                        && read.start > *start
                        && read.start <= *end
                })
                .count() as i32,
        })
        .collect()
}

/// The message when no output file was asked for at all.
pub const NO_OUTPUT_MESSAGE: &str = "You must supply at least one output file: PE, SR, SD, or RD";

/// Each writer refuses a file name it could not read back, with its own wording and its own three
/// extensions.
pub fn bad_name_message(kind: &str, filename: &str) -> String {
    let (what, extensions) = match kind {
        "pe" => (
            "discordant pair evidence",
            ".pe.txt\", \".pe.txt.gz\", or \".pe.bci",
        ),
        "sr" => (
            "split read evidence",
            ".sr.txt\", \".sr.txt.gz\", or \".sr.bci",
        ),
        "sd" => (
            "site depth evidence",
            ".sd.txt\", \".sd.txt.gz\", or \".sd.bci",
        ),
        "rd" => ("depth evidence", ".rd.txt\", \".rd.txt.gz\", or \".rd.bci"),
        other => panic!("an unknown writer {other}"),
    };
    format!(
        "Attempting to write {what} to a file that can't be read as {what}: {filename}.  The \
         file name should end with \"{extensions}\"."
    )
}

/// How a writer lays out the file its name selects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Encoding {
    /// Tab-separated text, block compressed when the name says so.
    Text { block_compressed: bool },
    /// The binary container, `.bci`.
    Bci,
}

/// `FileExtensions.BLOCK_COMPRESSED`, which `IOUtil.hasBlockCompressedExtension` tests.
const BLOCK_COMPRESSED: [&str; 4] = [".gz", ".gzip", ".bgz", ".bgzf"];

/// Which of the `kind` writer's two codecs claims `filename`, or `None` where neither does and the
/// writer refuses it with [`bad_name_message`].
///
/// The writer asks the binary codec first: it lower-cases the name and tests `.<kind>.bci`. The
/// text codec lower-cases it, strips ONE block-compressed extension, and tests `.<kind>.txt`, so
/// `x.PE.TXT.BGZ` is block-compressed discordant pair evidence and `x.pe.bci.gz` is neither
/// codec's.
pub fn encoding(kind: &str, filename: &str) -> Option<Encoding> {
    let lower = filename.to_lowercase();
    if lower.ends_with(&format!(".{kind}.bci")) {
        return Some(Encoding::Bci);
    }
    let block_compressed = BLOCK_COMPRESSED
        .iter()
        .any(|extension| lower.ends_with(extension));
    let stripped = if block_compressed {
        &lower[..lower.rfind('.').unwrap_or(lower.len())]
    } else {
        &lower[..]
    };
    stripped
        .ends_with(&format!(".{kind}.txt"))
        .then_some(Encoding::Text { block_compressed })
}

/// The message when the interval file holds nothing.
pub fn empty_intervals_message(filename: &str) -> String {
    format!("{filename} contains no intervals.")
}
