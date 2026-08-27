//! `CollectSVEvidence`: what one BAM contributes to structural-variant calling.
//!
//! One traversal writes four evidence files, and each has its own rule for which reads it will
//! look at: discordant pairs, split-read positions, per-site allele depths and per-interval read
//! counts.
//!
//! Reading the BAM, the VCF and the interval file is not ported, nor are the streaming buffers
//! that keep each writer's memory bounded. Which read contributes what, and what each contribution
//! is, are.

/// One read, reduced to what the four writers read off it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub name: String,
    pub contig_index: usize,
    pub contig: String,
    pub start: i32,
    pub mapping_quality: i32,
    /// Cigar as (operator, length), the operator being one of `M`, `S`, `D`, `I`, `N`, `H`.
    pub cigar: Vec<(char, i32)>,
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
    pub base_qualities: Vec<i32>,
}

/// The two default read filters this tool adds on top of the walker's own.
pub fn passes_default_filters(read: &Read) -> bool {
    !read.unmapped && !read.duplicate
}

/// The outer test in `apply`: a read whose MATE is unmapped, a supplementary alignment and a
/// secondary alignment are all skipped by both the split and the discordant writer, and by neither
/// of the two depth counters.
pub fn contributes_to_pair_evidence(read: &Read) -> bool {
    !(read.paired && read.mate_unmapped) && !read.supplementary && !read.secondary
}

/// `isSoftClipped`: EXACTLY one end. A read clipped at both ends is not a split read at all.
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
    Left,
    Right,
    /// Neither end is a match, so the read gives no position.
    Middle,
}

impl Side {
    pub fn name(self) -> &'static str {
        match self {
            Side::Left => "left",
            Side::Right => "right",
            Side::Middle => "middle",
        }
    }
}

fn consumes_reference(operator: char) -> bool {
    matches!(operator, 'M' | 'D' | 'N' | '=' | 'X')
}

/// `getSplitPosition`.
///
/// The FIRST element being a match wins, whatever the last one is, and the position it gives is the
/// start plus every reference-consuming length: a deletion inside the alignment moves it out by its
/// own length. Only when the first is not a match is the last one asked.
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
    pub contig: String,
    /// The codec writes it zero-based.
    pub position: i32,
    pub side: Side,
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
    // The buffer is a priority queue over position then side, flushed in that order.
    positions.sort_by(|a, b| (&a.0, a.1, a.2).cmp(&(&b.0, b.1, b.2)));
    let mut out: Vec<SplitRead> = Vec::new();
    for (contig, position, side) in positions {
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
    pub contig: String,
    /// Zero-based, as the codec writes it.
    pub position: i32,
    /// `!isReverseStrand`. The codec prints `+` for true, so the two negations cancel and the
    /// output shows the read's own strand.
    pub strand: bool,
    pub mate_contig: String,
    pub mate_position: i32,
    pub mate_strand: bool,
}

/// `getReportableDiscordantReadPair`: a pair is written ONCE, by whichever of its two reads has the
/// smaller contig index, or the smaller start, and at an EQUAL start by the first of the two SEEN.
///
/// The last case is tracked by NAME within the locus, so two different pairs at one position are
/// still two records: the name set only ever suppresses a read's own mate.
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
    pub contig: String,
    pub position: i32,
    pub reference: String,
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
pub const DEFAULT_SITE_DEPTH_MIN_MAPQ: i32 = 30;
pub const DEFAULT_SITE_DEPTH_MIN_BASEQ: i32 = 20;
pub const DEFAULT_DEPTH_EVIDENCE_MIN_MAPQ: i32 = 0;

/// One site-depth record: the four base counts, in `A`, `C`, `G`, `T` order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteDepth {
    pub contig: String,
    /// Zero-based, as the codec writes it.
    pub position: i32,
    pub counts: [i32; 4],
}

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
    pub contig: String,
    /// The interval file's own coordinates, written back unchanged.
    pub start: i32,
    pub end: i32,
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

/// The message when the interval file holds nothing.
pub fn empty_intervals_message(filename: &str) -> String {
    format!("{filename} contains no intervals.")
}
