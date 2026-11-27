//! Ported from `org.broadinstitute.hellbender.tools.walkers.rnaseq.ASEReadCounter`
//! (GATK 4.7.0.0; `ASEReadCounter.java` is unchanged since 4.6.2.0).
//!
//! Reference and alternate counts at heterozygous sites. The counting is trivial; the order of the
//! tests that decide which reads are counted is not, and neither is what "one read" means when two
//! mates overlap the same locus.
//!
//! # The cascade is ordered, and every test but the last one `continue`s
//!
//! Improper pair, then low mapping quality, then low base quality, then other-base. A read that
//! would fail two of them is charged to the FIRST and to nothing else, so the five depth columns
//! partition the pileup rather than overlapping. `rawDepth` is the only one that sees every
//! element, and it sees the pileup **after** the overlap filter has run.
//!
//! # Three answers to one overlapping pair
//!
//! `COUNT_FRAGMENTS_REQUIRE_SAME_BASE`, the default, keeps one element per read name and DISCARDS
//! THE PAIR ENTIRELY when the two disagree. `COUNT_FRAGMENTS` keeps the better-quality element
//! either way. `COUNT_READS` keeps both, so a disagreeing mate is counted as an other-base. The
//! golden has one site where all three differ.
//!
//! # A discarded pair stays discarded
//!
//! The filter keeps a set of the names it has deleted, so a third element with that read's name
//! later in the pileup does not resurrect it. That cannot happen with two mates, and is written out
//! because the reference wrote it out.
//!
//! ```java
//! if (discardDiscordant && existing.getBase() != p.getBase()) {
//!     filteredPileup.remove(readName);
//!     readNamesDeleted.add(readName);
//! } else if (tieBreaker.compare(existing, p) < 0) {
//!     filteredPileup.put(readName, p);
//! }
//! ```
//!
//! (`ReadPileup.filterSingleSampleForOverlaps`, `ReadPileup.java`, GATK 4.7.0.0.)
//!
//! # Reading this file without reading Rust
//!
//! - `&ReadPileup<'a>` is a borrowed, read-only pileup; the `'a` says its elements borrow their
//!   reads from something that lives at least as long, so a filtered pileup can hand the same
//!   reads on without copying them.
//! - `Option<String>` is "a line or nothing", Java's `null` return checked by the compiler.
//! - `continue` and `match` are as in Java, except that a `match` must cover every enum variant.
//! - `i32::from(x)` widens a byte to an `int` without sign extension: a Rust `u8` is unsigned,
//!   where Java's `byte` is signed.
//! - `flags & 0x1 != 0` tests one bit; Rust's `&` binds tighter than `!=`.

use gatk_engine::pileup::PileupElement;
use gatk_engine::read_pileup::ReadPileup;
use gatk_readfilter as filters;
use htsjdk_bam::record::BamRecord;

/// `CountPileupType`. `#[default]` marks the value an unset argument takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CountType {
    /// Every element, mates included.
    CountReads,
    /// One element per read name, the better quality winning a disagreement.
    CountFragments,
    /// One element per read name, a disagreement discarding both.
    #[default]
    CountFragmentsRequireSameBase,
}

/// `OUTPUT_FORMAT`, which decides the separator and nothing else.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutputFormat {
    /// Tab separated.
    Table,
    /// Tab separated, and identical to [`OutputFormat::Table`] in every byte.
    #[default]
    RTable,
    /// Comma separated.
    Csv,
}

impl OutputFormat {
    /// The separator this format writes.
    pub fn separator(&self) -> &'static str {
        match self {
            OutputFormat::Csv => ",",
            OutputFormat::Table | OutputFormat::RTable => "\t",
        }
    }
}

/// The eight filters `getDefaultReadFilters` returns, in its order
/// (`ASEReadCounter.java:153-160`).
pub const DEFAULT_READ_FILTERS: [&str; 8] = [
    "ValidAlignmentStartReadFilter",
    "ValidAlignmentEndReadFilter",
    "HasReadGroupReadFilter",
    "MatchingBasesAndQualsReadFilter",
    "SeqIsStoredReadFilter",
    "NotDuplicateReadFilter",
    "NotSecondaryAlignmentReadFilter",
    "MappedReadFilter",
];

/// The conjunction those eight names make.
///
/// The walker's own defaults are NOT included: this tool builds the list from an empty one, so
/// `WellformedReadFilter` never runs and a read this tool accepts may be one another tool refuses.
pub fn default_read_filter(read: &BamRecord) -> bool {
    filters::valid_alignment_start(read)
        && filters::valid_alignment_end(read)
        && filters::has_read_group(read)
        && filters::matching_bases_and_quals(read)
        && filters::seq_is_stored(read)
        && filters::not_duplicate(read)
        && filters::not_secondary_alignment(read)
        && filters::mapped(read)
}

/// The counts one site produces, in the order the table prints them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SiteCounts {
    /// Elements carrying the reference base and passing every test.
    pub reference: i32,
    /// Elements carrying the alternate base and passing every test.
    pub alternate: i32,
    /// `refCount + altCount`, which is what the depth threshold is compared against.
    pub total: i32,
    /// Elements stopped by the mapping quality threshold.
    pub low_mapping_quality: i32,
    /// Elements stopped by the base quality threshold.
    pub low_base_quality: i32,
    /// Every element in the filtered pileup, whatever became of it.
    pub raw: i32,
    /// Elements carrying neither the reference nor the alternate base.
    pub other_bases: i32,
    /// Elements whose read is paired and not properly paired, or whose mate is unmapped.
    pub improper_pairs: i32,
}

/// `calculateLineForSite`'s counting, without the formatting.
///
/// ```java
/// if (base.getMappingQual() < minMappingQuality) {
///     lowMAPQDepth++;
///     continue;
/// }
///
/// if (base.getQual() < minBaseQuality) {
/// ```
///
/// (`ASEReadCounter.java:261-266`.) `getQual()` is a Java `byte`, compared after sign extension;
/// here `qual()` is a `u8` compared unsigned. The two agree for every quality below 128.
pub fn count_site(
    pileup: &ReadPileup,
    reference_allele: u8,
    alternate_allele: u8,
    minimum_mapping_quality: i32,
    minimum_base_quality: u8,
) -> SiteCounts {
    let mut counts = SiteCounts::default();
    // `&pileup.elements` iterates by reference: the loop reads the elements and leaves them in the
    // pileup.
    for element in &pileup.elements {
        counts.raw += 1;
        let read = element.read;
        // `isPaired() && (mateIsUnmapped() || !isProperlyPaired())`, which is the first bucket and
        // therefore the one an improperly paired read is charged to whatever else is wrong with it.
        if is_paired(read) && (mate_is_unmapped(read) || !is_properly_paired(read)) {
            counts.improper_pairs += 1;
            continue;
        }
        if i32::from(element.mapping_qual()) < minimum_mapping_quality {
            counts.low_mapping_quality += 1;
            continue;
        }
        if element.qual() < minimum_base_quality {
            counts.low_base_quality += 1;
            continue;
        }
        // A base that is neither the reference nor the alternate is charged to `otherBases` and
        // skips `totalCount`, which only the two alleles reach.
        if element.base() == reference_allele {
            counts.reference += 1;
        } else if element.base() == alternate_allele {
            counts.alternate += 1;
        } else {
            counts.other_bases += 1;
            continue;
        }
        counts.total += 1;
    }
    counts
}

/// `SAMFlag.READ_PAIRED`.
fn is_paired(read: &BamRecord) -> bool {
    read.flags & 0x1 != 0
}

/// `SAMFlag.PROPER_PAIR`, read only on a paired read: `SAMRecordToGATKReadAdapter.isProperlyPaired`
/// is `isPaired() && samRecord.getProperPairFlag()`.
fn is_properly_paired(read: &BamRecord) -> bool {
    is_paired(read) && read.flags & 0x2 != 0
}

/// `SAMFlag.MATE_UNMAPPED`.
fn mate_is_unmapped(read: &BamRecord) -> bool {
    read.flags & 0x8 != 0
}

/// `getOverlappingFragmentFilteredPileup`, restricted to one sample.
///
/// The reference splits the pileup by sample with `splitBySample` and filters each sample's
/// elements in a map of its own; here the whole pileup is one map. The two agree whenever the
/// pileup holds one sample, or no read name is shared between samples.
///
/// One element per read name. A second element with a name already kept either replaces it (when
/// its base quality is higher) or, when `discard_discordant` and the bases differ, DELETES BOTH and
/// records the name so a third element cannot bring it back.
///
/// The reference collects the survivors out of a `HashMap`, so their order is the map's rather than
/// the pileup's. Nothing this tool does with the result depends on the order, which is why the port
/// keeps the pileup's instead of reproducing a Java hash table.
pub fn filter_overlapping_fragments<'a>(
    pileup: &ReadPileup<'a>,
    discard_discordant: bool,
) -> ReadPileup<'a> {
    // `kept` is the `HashMap<String, PileupElement>` as a list of (name, element) pairs in
    // first-seen order, and `deleted` is `readNamesDeleted`. Linear lookups, since a pileup's
    // names are few.
    let mut kept: Vec<(String, PileupElement<'a>)> = Vec::new();
    let mut deleted: Vec<String> = Vec::new();

    for element in &pileup.elements {
        let name = element.read.read_name.clone();
        match kept.iter().position(|(existing, _)| *existing == name) {
            None => {
                if !deleted.contains(&name) {
                    kept.push((name, element.clone()));
                }
            }
            Some(index) => {
                let existing = kept[index].1.clone();
                if discard_discordant && existing.base() != element.base() {
                    kept.remove(index);
                    deleted.push(name);
                // Strictly less: on equal qualities the element kept first stays.
                } else if base_qual_tie_breaker(&existing, element) < 0 {
                    kept[index].1 = element.clone();
                }
            }
        }
    }

    ReadPileup::new(
        &pileup.contig,
        pileup.start,
        kept.into_iter().map(|(_, element)| element).collect(),
    )
}

/// `ReadPileup.baseQualTieBreaker`: compare the two base qualities.
///
/// ```java
/// return Byte.compare(o1.getQual(), o2.getQual());
/// ```
///
/// `Byte.compare` is signed and this subtraction is of unsigned bytes, so the two order a quality
/// of 128 or more differently; only the sign of the answer is read.
fn base_qual_tie_breaker(left: &PileupElement, right: &PileupElement) -> i32 {
    i32::from(left.qual()) - i32::from(right.qual())
}

/// `filterPileup`: the overlap handling, then the deletions.
///
/// ```java
/// return pileupWithDeletions.makeFilteredPileup(p -> !p.isDeletion());
/// ```
///
/// The deletions go after the overlap filter, so a mate whose element is a deletion still takes
/// part in the overlap test: under the default it can discard a pair whose other mate carries a
/// base, because a deletion's base is not that base.
pub fn filter_pileup<'a>(pileup: &ReadPileup<'a>, count_type: CountType) -> ReadPileup<'a> {
    let with_deletions = match count_type {
        // `COUNT_READS` hands the pileup through untouched, which is one `ReadPileup` rebuilt
        // rather than cloned: the elements borrow their reads and the struct is not `Clone`.
        CountType::CountReads => {
            ReadPileup::new(&pileup.contig, pileup.start, pileup.elements.clone())
        }
        CountType::CountFragments => filter_overlapping_fragments(pileup, false),
        CountType::CountFragmentsRequireSameBase => filter_overlapping_fragments(pileup, true),
    };
    with_deletions.filtered(|element| !element.is_deletion())
}

/// The header line, which is printed before the traversal and therefore even for an empty run.
pub fn header(format: OutputFormat) -> String {
    [
        "contig",
        "position",
        "variantID",
        "refAllele",
        "altAllele",
        "refCount",
        "altCount",
        "totalCount",
        "lowMAPQDepth",
        "lowBaseQDepth",
        "rawDepth",
        "otherBases",
        "improperPairs",
    ]
    .join(format.separator())
}

/// One site's line, or nothing when it is below the depth threshold: `calculateLineForSite`'s
/// `StringBuilder`, the thirteen columns joined by the format's separator. The alleles are
/// written as characters, as the reference's `char refAllele` appends.
#[allow(clippy::too_many_arguments)]
pub fn line(
    contig: &str,
    position: i32,
    site_id: &str,
    reference_allele: u8,
    alternate_allele: u8,
    counts: SiteCounts,
    minimum_depth: i32,
    format: OutputFormat,
) -> Option<String> {
    // The comparison is against `totalCount`, the reference plus alternate, and not against the
    // raw depth: a site can be deep and still produce no line.
    if counts.total < minimum_depth {
        return None;
    }
    Some(
        [
            contig.to_string(),
            position.to_string(),
            site_id.to_string(),
            (reference_allele as char).to_string(),
            (alternate_allele as char).to_string(),
            counts.reference.to_string(),
            counts.alternate.to_string(),
            counts.total.to_string(),
            counts.low_mapping_quality.to_string(),
            counts.low_base_quality.to_string(),
            counts.raw.to_string(),
            counts.other_bases.to_string(),
            counts.improper_pairs.to_string(),
        ]
        .join(format.separator()),
    )
}
