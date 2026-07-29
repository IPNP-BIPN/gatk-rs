//! Ported from `org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary` (GATK 4.6.2.0).
//!
//! The read filters are the first thing ported in this repository, on purpose. They are stateless,
//! they touch no floating point, and every tool that reads reads runs a chain of them, so the 55
//! of them are both the cheapest thing to get right and the widest-reaching. A wrong filter does
//! not produce a wrong number, it produces a different set of reads, and every number downstream
//! inherits that.
//!
//! # The part that is easy to get wrong
//!
//! A filter reads its predicate off `GATKRead`, not off the SAM flags, and the two are not the
//! same thing. `SAMRecordToGATKReadAdapter.isUnmapped` is:
//!
//! ```java
//! samRecord.getReadUnmappedFlag()
//!     || samRecord.getReferenceName() == null
//!     || samRecord.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)
//!     || samRecord.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START
//! ```
//!
//! Three criteria, of which the 0x4 flag is one. A record with the flag clear, a reference index of
//! -1 and a start of 0 is mapped by the flag and unmapped by GATK, and a port that tested the flag
//! alone would keep it. The same shape applies to `mateIsUnmapped`, and `isFirstOfPair` is
//! `isPaired() && firstOfPairFlag`, not the 0x40 flag on its own.
//!
//! Each filter below therefore records which GATKRead accessor it goes through, and the accessors
//! are implemented once, in [`read`], rather than inlined per filter.

use htsjdk_bam::record::BamRecord;

/// The SAM flag bits, as `SAMFlag` defines them.
pub mod flags {
    pub const READ_PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const READ_UNMAPPED: u16 = 0x4;
    pub const MATE_UNMAPPED: u16 = 0x8;
    pub const READ_REVERSE_STRAND: u16 = 0x10;
    pub const MATE_REVERSE_STRAND: u16 = 0x20;
    pub const FIRST_OF_PAIR: u16 = 0x40;
    pub const SECOND_OF_PAIR: u16 = 0x80;
    pub const NOT_PRIMARY_ALIGNMENT: u16 = 0x100;
    pub const READ_FAILS_VENDOR_QUALITY_CHECK: u16 = 0x200;
    pub const DUPLICATE_READ: u16 = 0x400;
    pub const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;
}

/// `SAMRecord.NO_ALIGNMENT_START`.
pub const NO_ALIGNMENT_START: i32 = 0;
/// The index htsjdk uses for `SAMRecord.NO_ALIGNMENT_REFERENCE_NAME` ("*").
pub const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;
/// `QualityUtils.MAPPING_QUALITY_UNAVAILABLE`.
pub const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;

/// The `GATKRead` accessors the filters are written against.
///
/// They live here rather than inside each filter because their definitions are where the
/// divergences hide: three of them are not the flag test they look like.
pub mod read {
    use super::*;

    pub fn is_paired(read: &BamRecord) -> bool {
        read.flags & flags::READ_PAIRED != 0
    }

    /// `SAMRecordToGATKReadAdapter.isUnmapped`: the flag, an absent reference, or a zero start.
    pub fn is_unmapped(read: &BamRecord) -> bool {
        read.flags & flags::READ_UNMAPPED != 0
            || read.reference_index == NO_ALIGNMENT_REFERENCE_INDEX
            || read.alignment_start == NO_ALIGNMENT_START
    }

    /// `SAMRecordToGATKReadAdapter.mateIsUnmapped`, same three criteria applied to the mate.
    ///
    /// The Java asserts `isPaired()` first and throws otherwise; here that is the caller's job, and
    /// every caller in this crate is a filter that has already tested pairing.
    pub fn mate_is_unmapped(read: &BamRecord) -> bool {
        read.flags & flags::MATE_UNMAPPED != 0
            || read.mate_reference_index == NO_ALIGNMENT_REFERENCE_INDEX
            || read.mate_alignment_start == NO_ALIGNMENT_START
    }

    /// `isPaired() && firstOfPairFlag`, not the 0x40 flag alone.
    pub fn is_first_of_pair(read: &BamRecord) -> bool {
        is_paired(read) && read.flags & flags::FIRST_OF_PAIR != 0
    }

    /// `isPaired() && secondOfPairFlag`.
    pub fn is_second_of_pair(read: &BamRecord) -> bool {
        is_paired(read) && read.flags & flags::SECOND_OF_PAIR != 0
    }

    pub fn is_proper_pair(read: &BamRecord) -> bool {
        is_paired(read) && read.flags & flags::PROPER_PAIR != 0
    }

    pub fn is_duplicate(read: &BamRecord) -> bool {
        read.flags & flags::DUPLICATE_READ != 0
    }

    pub fn is_secondary_alignment(read: &BamRecord) -> bool {
        read.flags & flags::NOT_PRIMARY_ALIGNMENT != 0
    }

    pub fn is_supplementary_alignment(read: &BamRecord) -> bool {
        read.flags & flags::SUPPLEMENTARY_ALIGNMENT != 0
    }

    pub fn fails_vendor_quality_check(read: &BamRecord) -> bool {
        read.flags & flags::READ_FAILS_VENDOR_QUALITY_CHECK != 0
    }

    pub fn is_reverse_strand(read: &BamRecord) -> bool {
        read.flags & flags::READ_REVERSE_STRAND != 0
    }

    pub fn mate_is_reverse_strand(read: &BamRecord) -> bool {
        read.flags & flags::MATE_REVERSE_STRAND != 0
    }

    pub fn fragment_length(read: &BamRecord) -> i32 {
        read.inferred_insert_size
    }

    /// `GATKRead.getLength()`: the number of bases, which is what the record carries.
    pub fn length(read: &BamRecord) -> usize {
        read.read_bases.len()
    }

    /// `GATKRead.getBaseQualityCount()`.
    ///
    /// htsjdk encodes "qualities absent" as a run of 0xFF of the read's length, and htsjdk-rs
    /// represents that as an empty vector, so an absent quality array counts as zero here exactly
    /// as `SAMRecord.getBaseQualities()` returns an empty array for it.
    pub fn base_quality_count(read: &BamRecord) -> usize {
        read.base_qualities.len()
    }

    pub fn mapping_quality(read: &BamRecord) -> u8 {
        read.mapping_quality
    }

    /// Whether the record carries an `RG` tag at all, which is what `getReadGroup() != null` tests.
    ///
    /// Only presence, deliberately. `SAMRecord.getReadGroup()` resolves the tag against the
    /// header's `@RG` lines and returns null when the header has no such group, so a record whose
    /// `RG` names a group the header does not declare is filtered out by GATK and kept here. That
    /// gap closes when the ported filters take a header; it is recorded rather than left implicit,
    /// because it is the kind of difference that only shows on a malformed file.
    pub fn has_read_group(read: &BamRecord) -> bool {
        read.tags.iter().any(|(tag, _)| tag.name() == *b"RG")
    }
}

/// `htsjdk.samtools.Cigar` validity and the two GATK predicates built on it.
///
/// `GoodCigarReadFilter` calls `CigarUtils.isGood`, which is `Cigar.isValid(null, -1) == null`
/// plus two GATK rules. The htsjdk half is the long one and the easy one to approximate: it is ported
/// here in full rather than reduced to "looks sensible", because every rule it applies is a way a
/// read can be kept by the reference and dropped by the port.
pub mod cigar {
    use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

    fn is_clipping(op: Op) -> bool {
        matches!(op, Op::S | Op::H)
    }

    /// `Cigar.isRealOperator`: M, =, X, I, D, N.
    fn is_real(op: Op) -> bool {
        matches!(op, Op::M | Op::Eq | Op::X | Op::I | Op::D | Op::N)
    }

    fn is_indel(op: Op) -> bool {
        matches!(op, Op::I | Op::D)
    }

    fn is_padding(op: Op) -> bool {
        matches!(op, Op::P)
    }

    /// `Cigar.isValid(null, -1) == null`: no validation error of any kind.
    ///
    /// An empty cigar returns `null` in the Java, meaning *valid*, which is not the reading the
    /// name suggests: an unmapped read with `*` for its cigar passes.
    pub fn is_valid(c: &Cigar) -> bool {
        let elements: &[CigarElement] = &c.elements;
        if elements.is_empty() {
            return true;
        }
        let mut seen_real = false;
        for (i, element) in elements.iter().enumerate() {
            if element.length == 0 {
                return false;
            }
            let op = element.op;
            if is_clipping(op) {
                if op == Op::H {
                    // Hard clips only at either end.
                    if i != 0 && i != elements.len() - 1 {
                        return false;
                    }
                } else if i == 0 || i == elements.len() - 1 {
                    // A soft clip at either end is fine.
                } else if i == 1 {
                    // The special case the Java calls funky: S is both one from the beginning and
                    // one from the end.
                    let funky = elements.len() == 3 && elements[2].op == Op::H;
                    if !funky && elements[0].op != Op::H {
                        return false;
                    }
                } else if i == elements.len() - 2 {
                    if elements[elements.len() - 1].op != Op::H {
                        return false;
                    }
                } else {
                    return false;
                }
            } else if is_real(op) {
                seen_real = true;
                if is_indel(op) {
                    // There must be an M, N or P between any pair of the *same* indel operator.
                    for next in &elements[i + 1..] {
                        let next_op = next.op;
                        if (is_real(next_op) && !is_indel(next_op)) || is_padding(next_op) {
                            break;
                        }
                        if is_indel(next_op) && next_op == op {
                            return false;
                        }
                    }
                }
            } else if is_padding(op) && i != 0 {
                // Position 0 is allowed: a read starting inside a pad needs leading padding, and
                // the Java carries a comment saying the restriction was removed deliberately.
                // Everything else must sit between two real operators, which the last position
                // cannot do. The two rejections are one condition here because they are one
                // condition in effect; the Java writes them apart only to word two messages.
                let at_end = i == elements.len() - 1;
                if at_end || !is_real(elements[i - 1].op) || !is_real(elements[i + 1].op) {
                    return false;
                }
            }
        }
        seen_real
    }

    /// `CigarUtils.hasConsecutiveIndels`.
    fn has_consecutive_indels(elements: &[CigarElement]) -> bool {
        let mut previous_indel = false;
        for element in elements {
            let indel = is_indel(element.op);
            if previous_indel && indel {
                return true;
            }
            previous_indel = indel;
        }
        false
    }

    /// `CigarUtils.startsOrEndsWithDeletionIgnoringClips`.
    fn starts_or_ends_with_deletion_ignoring_clips(elements: &[CigarElement]) -> bool {
        for from_left in [true, false] {
            let iter: Box<dyn Iterator<Item = &CigarElement>> = if from_left {
                Box::new(elements.iter())
            } else {
                Box::new(elements.iter().rev())
            };
            for element in iter {
                if element.op == Op::D {
                    return true;
                } else if !is_clipping(element.op) {
                    break;
                }
            }
        }
        false
    }

    /// `CigarUtils.isGood`.
    pub fn is_good(c: &Cigar) -> bool {
        if !is_valid(c) {
            return false;
        }
        let elements = &c.elements;
        !(has_consecutive_indels(elements) || starts_or_ends_with_deletion_ignoring_clips(elements))
    }

    /// `CigarUtils.containsNOperator`.
    pub fn contains_n_operator(c: &Cigar) -> bool {
        c.elements.iter().any(|e| e.op == Op::N)
    }
}

/// A read filter: `true` keeps the read, exactly as `ReadFilter.test` does.
pub type ReadFilter = fn(&BamRecord) -> bool;

/// `ReadFilterLibrary.AllowAllReadsReadFilter`: "Do not filter out any read."
pub fn allow_all_reads(_read: &BamRecord) -> bool {
    true
}

/// `ReadFilterLibrary.MappedReadFilter`: filter out unmapped reads.
pub fn mapped(read: &BamRecord) -> bool {
    !read::is_unmapped(read)
}

/// `ReadFilterLibrary.MappingQualityAvailableReadFilter`: filter out MAPQ 255.
pub fn mapping_quality_available(read: &BamRecord) -> bool {
    read::mapping_quality(read) != MAPPING_QUALITY_UNAVAILABLE
}

/// `ReadFilterLibrary.MappingQualityNotZeroReadFilter`.
pub fn mapping_quality_not_zero(read: &BamRecord) -> bool {
    read::mapping_quality(read) != 0
}

/// `ReadFilterLibrary.NotDuplicateReadFilter`.
pub fn not_duplicate(read: &BamRecord) -> bool {
    !read::is_duplicate(read)
}

/// `ReadFilterLibrary.NotSecondaryAlignmentReadFilter`.
pub fn not_secondary_alignment(read: &BamRecord) -> bool {
    !read::is_secondary_alignment(read)
}

/// `ReadFilterLibrary.NotSupplementaryAlignmentReadFilter`.
pub fn not_supplementary_alignment(read: &BamRecord) -> bool {
    !read::is_supplementary_alignment(read)
}

/// `ReadFilterLibrary.PassesVendorQualityCheckReadFilter`.
pub fn passes_vendor_quality_check(read: &BamRecord) -> bool {
    !read::fails_vendor_quality_check(read)
}

/// `ReadFilterLibrary.PairedReadFilter`.
pub fn paired(read: &BamRecord) -> bool {
    read::is_paired(read)
}

/// `ReadFilterLibrary.NotProperlyPairedReadFilter`.
///
/// `read.isPaired() && !read.isProperlyPaired()`, which is **not** the negation of
/// `ProperlyPairedReadFilter`: an unpaired read is filtered out by both. The first version here
/// was the negation, and it kept every unpaired read; the oracle's decision matrix caught it on
/// the first run, on five records of a nineteen-record corpus.
pub fn not_properly_paired(read: &BamRecord) -> bool {
    read::is_paired(read) && !read::is_proper_pair(read)
}

/// `ReadFilterLibrary.ProperlyPairedReadFilter`, for the contrast with the filter above.
pub fn properly_paired(read: &BamRecord) -> bool {
    read::is_proper_pair(read)
}

/// `ReadFilterLibrary.FirstOfPairReadFilter`.
pub fn first_of_pair(read: &BamRecord) -> bool {
    read::is_first_of_pair(read)
}

/// `ReadFilterLibrary.SecondOfPairReadFilter`.
pub fn second_of_pair(read: &BamRecord) -> bool {
    read::is_second_of_pair(read)
}

/// `ReadFilterLibrary.NonZeroFragmentLengthReadFilter`.
pub fn non_zero_fragment_length(read: &BamRecord) -> bool {
    read::fragment_length(read) != 0
}

/// `ReadFilterLibrary.MatchingBasesAndQualsReadFilter`.
pub fn matching_bases_and_quals(read: &BamRecord) -> bool {
    read::length(read) == read::base_quality_count(read)
}

/// `ReadFilterLibrary.HasReadGroupReadFilter`.
pub fn has_read_group(read: &BamRecord) -> bool {
    read::has_read_group(read)
}

/// `ReadFilterLibrary.MateDifferentStrandReadFilter`.
///
/// Keep only paired reads whose mate maps to the opposite strand, both ends mapped.
pub fn mate_different_strand(read: &BamRecord) -> bool {
    read::is_paired(read)
        && !read::is_unmapped(read)
        && !read::mate_is_unmapped(read)
        && read::mate_is_reverse_strand(read) != read::is_reverse_strand(read)
}

/// `ReadFilterLibrary.MateOnSameContigOrNoMappedMateReadFilter`.
///
/// Keep a read that is unpaired, whose mate is unmapped, or whose mate is on this read's contig.
pub fn mate_on_same_contig_or_no_mapped_mate(read: &BamRecord) -> bool {
    !read::is_paired(read)
        || read::mate_is_unmapped(read)
        || read.mate_reference_index == read.reference_index
}

/// `ReadFilterLibrary.CigarContainsNoNOperator`.
pub fn cigar_contains_no_n_operator(read: &BamRecord) -> bool {
    !cigar::contains_n_operator(&read.cigar)
}

/// `ReadFilterLibrary.GoodCigarReadFilter`.
pub fn good_cigar(read: &BamRecord) -> bool {
    cigar::is_good(&read.cigar)
}

/// `ReadFilterLibrary.NonZeroReferenceLengthAlignmentReadFilter`.
pub fn non_zero_reference_length_alignment(read: &BamRecord) -> bool {
    read.cigar
        .elements
        .iter()
        .any(|e| e.op.consumes_reference_bases() && e.length > 0)
}

/// `ReadFilterLibrary.PrimaryLineReadFilter`: neither secondary nor supplementary.
pub fn primary_line(read: &BamRecord) -> bool {
    !read::is_secondary_alignment(read) && !read::is_supplementary_alignment(read)
}

/// `ReadFilterLibrary.ReadLengthEqualsCigarLengthReadFilter`.
///
/// Unmapped reads pass unconditionally, which matters: their cigar is `*`, so the comparison would
/// otherwise reject every one of them.
pub fn read_length_equals_cigar_length(read: &BamRecord) -> bool {
    read::is_unmapped(read) || read::length(read) as u32 == read.cigar.read_length()
}

/// `ReadFilterLibrary.SeqIsStoredReadFilter`.
pub fn seq_is_stored(read: &BamRecord) -> bool {
    read::length(read) > 0
}

/// `ReadFilterLibrary.ValidAlignmentStartReadFilter`.
///
/// `GATKRead.getStart()` returns `ReadConstants.UNSET_POSITION` for an unmapped read, so the
/// filter's first clause is what keeps this from testing a sentinel.
pub fn valid_alignment_start(read: &BamRecord) -> bool {
    read::is_unmapped(read) || read.alignment_start > 0
}

/// `ReadFilterLibrary.ValidAlignmentEndReadFilter`: `end - start + 1 >= 0`.
///
/// The condition looks vacuous and is not: `getEnd()` is derived from the cigar's reference
/// length, and a cigar that consumes no reference gives an end before the start.
pub fn valid_alignment_end(read: &BamRecord) -> bool {
    read::is_unmapped(read) || (read.alignment_end() - read.alignment_start + 1) >= 0
}

/// `ReadFilterLibrary.MateUnmappedAndUnmappedReadFilter`.
///
/// Keeps a read only when neither it nor (if paired) its mate is unmapped.
pub fn mate_unmapped_and_unmapped(read: &BamRecord) -> bool {
    !(read::is_unmapped(read) || (read::is_paired(read) && read::mate_is_unmapped(read)))
}

/// `ReadFilterLibrary.NonChimericOriginalAlignmentReadFilter`.
///
/// Compares the contig in the `OA` tag with the mate contig `AddOriginalAlignmentTags` writes.
///
/// That tag is **`XM`**, not `MC`. The first version of this port assumed `MC`, the SAM standard
/// mate-cigar-adjacent tag, and so did the conformance corpus, so the filter returned "no tags,
/// pass" on every record and the two sides agreed perfectly while testing nothing. A corpus
/// written from the same assumption as the port confirms the assumption, not the behaviour. The
/// constant is `AddOriginalAlignmentTags.MATE_CONTIG_TAG_NAME`.
///
/// A read missing either tag passes.
pub fn non_chimeric_original_alignment(read: &BamRecord) -> bool {
    let oa = read.tags.iter().find(|(t, _)| t.name() == *b"OA");
    let mate_contig = read.tags.iter().find(|(t, _)| t.name() == *b"XM");
    match (oa, mate_contig) {
        (Some((_, oa_value)), Some((_, mate_value))) => {
            // AddOriginalAlignmentTags.getOAContig takes the first comma-separated field of OA.
            let oa_text = tag_text(oa_value);
            let contig = oa_text.split(',').next().unwrap_or("");
            contig == tag_text(mate_value)
        }
        _ => true,
    }
}

fn tag_text(value: &htsjdk_bam::tag::TagValue) -> &str {
    match value {
        htsjdk_bam::tag::TagValue::Str(text) => text,
        _ => "",
    }
}

/// The filters that take arguments.
///
/// They are modelled as data rather than as functions because their behaviour is a function of
/// their parameters, and the parameters are what the command line supplies. Porting the *logic*
/// here keeps it independent of the Barclay argument layer, which is a separate slice: a filter
/// can be proved byte-identical against the reference long before anything can parse
/// `--read-filter MappingQualityReadFilter --minimum-mapping-quality 20`.
///
/// The conformance golden names each instance with its parameters, and [`Parameterized::parse`]
/// rebuilds it from that name, so the reference's own instantiation is what the port is tested
/// against rather than a second list that could drift.
#[derive(Debug, Clone, PartialEq)]
pub enum Parameterized {
    /// `MappingQualityReadFilter`. The maximum is optional and absent by default.
    MappingQuality { min: i32, max: Option<i32> },
    /// `ReadLengthReadFilter`.
    ReadLength { min: i32, max: i32 },
    /// `FragmentLengthReadFilter`.
    FragmentLength { min: i32, max: i32 },
    /// `MateDistantReadFilter`.
    MateDistant { threshold: i32 },
    /// `ReadNameReadFilter`.
    ReadName { names: Vec<String> },
    /// `ReadStrandFilter`.
    ReadStrand { keep_reverse: bool },
    /// `AmbiguousBaseReadFilter`: an absolute count when given, otherwise a fraction of the read.
    AmbiguousBase {
        max_bases: Option<i32>,
        max_fraction: f64,
    },
}

impl Parameterized {
    pub fn test(&self, read: &BamRecord) -> bool {
        match self {
            Parameterized::MappingQuality { min, max } => {
                let mq = read::mapping_quality(read) as i32;
                mq >= *min && max.is_none_or(|limit| mq <= limit)
            }
            Parameterized::ReadLength { min, max } => {
                let length = read::length(read) as i32;
                length >= *min && length <= *max
            }
            Parameterized::FragmentLength { min, max } => {
                // An unpaired read passes: the fragment length means nothing for it.
                if !read::is_paired(read) {
                    return true;
                }
                // Negative when the mate maps before the read, hence the absolute value.
                let length = read::fragment_length(read).abs();
                length <= *max && length >= *min
            }
            Parameterized::MateDistant { threshold } => {
                read::is_paired(read)
                    && !read::is_unmapped(read)
                    && !read::mate_is_unmapped(read)
                    && ((read.alignment_start - read.mate_alignment_start).abs() >= *threshold
                        || read.reference_index != read.mate_reference_index)
            }
            Parameterized::ReadName { names } => names.contains(&read.read_name),
            Parameterized::ReadStrand { keep_reverse } => {
                read::is_reverse_strand(read) == *keep_reverse
            }
            Parameterized::AmbiguousBase {
                max_bases,
                max_fraction,
            } => {
                // `(int)(length * fraction)` truncates towards zero in Java, which `as i32` also
                // does for the non-negative values a read length can take.
                let max_n =
                    max_bases.unwrap_or_else(|| (read::length(read) as f64 * max_fraction) as i32);
                let mut ambiguous = 0;
                for base in &read.read_bases {
                    if !is_regular_base(*base) {
                        ambiguous += 1;
                        if ambiguous > max_n {
                            return false;
                        }
                    }
                }
                true
            }
        }
    }

    /// Rebuild an instance from the golden's `Name(key=value,...)` label.
    pub fn parse(spec: &str) -> Option<Parameterized> {
        let (name, rest) = spec.split_once('(')?;
        let args = rest.strip_suffix(')')?;
        let mut values: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
        for pair in args.split(',').filter(|p| !p.is_empty()) {
            let (key, value) = pair.split_once('=')?;
            values.insert(key, value);
        }
        let int = |key: &str| -> Option<i32> { values.get(key)?.parse().ok() };
        let maybe_int = |key: &str| -> Option<Option<i32>> {
            match *values.get(key)? {
                "null" => Some(None),
                other => other.parse().ok().map(Some),
            }
        };
        Some(match name {
            "MappingQualityReadFilter" => Parameterized::MappingQuality {
                min: int("min")?,
                max: maybe_int("max")?,
            },
            "ReadLengthReadFilter" => Parameterized::ReadLength {
                min: int("min")?,
                max: int("max")?,
            },
            "FragmentLengthReadFilter" => Parameterized::FragmentLength {
                min: int("min")?,
                max: int("max")?,
            },
            "MateDistantReadFilter" => Parameterized::MateDistant {
                threshold: int("threshold")?,
            },
            "ReadNameReadFilter" => Parameterized::ReadName {
                names: values
                    .get("names")?
                    .split('+')
                    .map(str::to_string)
                    .collect(),
            },
            "ReadStrandFilter" => Parameterized::ReadStrand {
                keep_reverse: *values.get("keepReverse")? == "true",
            },
            "AmbiguousBaseReadFilter" => Parameterized::AmbiguousBase {
                max_bases: maybe_int("maxBases")?,
                max_fraction: values.get("maxFraction")?.parse().ok()?,
            },
            _ => return None,
        })
    }
}

/// `BaseUtils.isRegularBase`.
///
/// The trap: `*` is a **regular** base. `BaseUtils.baseIndexMap['*']` is set to A's index with the
/// comment "the wildcard character counts as an A", so a read full of `*` has no ambiguous bases
/// at all. A port that tested `matches!(base, b'A' | b'C' | b'G' | b'T')` would count them.
pub fn is_regular_base(base: u8) -> bool {
    matches!(
        base,
        b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'*'
    )
}

/// The filters by the name GATK exposes on the command line (`--read-filter <Name>`).
pub fn by_name(name: &str) -> Option<ReadFilter> {
    Some(match name {
        "AllowAllReadsReadFilter" => allow_all_reads as ReadFilter,
        "MappedReadFilter" => mapped,
        "MappingQualityAvailableReadFilter" => mapping_quality_available,
        "MappingQualityNotZeroReadFilter" => mapping_quality_not_zero,
        "NotDuplicateReadFilter" => not_duplicate,
        "NotSecondaryAlignmentReadFilter" => not_secondary_alignment,
        "NotSupplementaryAlignmentReadFilter" => not_supplementary_alignment,
        "PassesVendorQualityCheckReadFilter" => passes_vendor_quality_check,
        "PairedReadFilter" => paired,
        "NotProperlyPairedReadFilter" => not_properly_paired,
        "FirstOfPairReadFilter" => first_of_pair,
        "SecondOfPairReadFilter" => second_of_pair,
        "NonZeroFragmentLengthReadFilter" => non_zero_fragment_length,
        "MatchingBasesAndQualsReadFilter" => matching_bases_and_quals,
        "HasReadGroupReadFilter" => has_read_group,
        "MateDifferentStrandReadFilter" => mate_different_strand,
        "MateOnSameContigOrNoMappedMateReadFilter" => mate_on_same_contig_or_no_mapped_mate,
        "ProperlyPairedReadFilter" => properly_paired,
        "CigarContainsNoNOperator" => cigar_contains_no_n_operator,
        "GoodCigarReadFilter" => good_cigar,
        "NonZeroReferenceLengthAlignmentReadFilter" => non_zero_reference_length_alignment,
        "PrimaryLineReadFilter" => primary_line,
        "ReadLengthEqualsCigarLengthReadFilter" => read_length_equals_cigar_length,
        "SeqIsStoredReadFilter" => seq_is_stored,
        "ValidAlignmentStartReadFilter" => valid_alignment_start,
        "ValidAlignmentEndReadFilter" => valid_alignment_end,
        "MateUnmappedAndUnmappedReadFilter" => mate_unmapped_and_unmapped,
        "NonChimericOriginalAlignmentReadFilter" => non_chimeric_original_alignment,
        _ => return None,
    })
}

/// Every filter ported so far, by name. The conformance harness iterates this, so a filter that is
/// added here is exercised against the oracle without touching the harness.
pub const PORTED: &[&str] = &[
    "AllowAllReadsReadFilter",
    "CigarContainsNoNOperator",
    "FirstOfPairReadFilter",
    "GoodCigarReadFilter",
    "HasReadGroupReadFilter",
    "MappedReadFilter",
    "MappingQualityAvailableReadFilter",
    "MappingQualityNotZeroReadFilter",
    "MatchingBasesAndQualsReadFilter",
    "MateDifferentStrandReadFilter",
    "MateOnSameContigOrNoMappedMateReadFilter",
    "MateUnmappedAndUnmappedReadFilter",
    "NonChimericOriginalAlignmentReadFilter",
    "NonZeroFragmentLengthReadFilter",
    "NonZeroReferenceLengthAlignmentReadFilter",
    "NotDuplicateReadFilter",
    "NotProperlyPairedReadFilter",
    "NotSecondaryAlignmentReadFilter",
    "NotSupplementaryAlignmentReadFilter",
    "PairedReadFilter",
    "PassesVendorQualityCheckReadFilter",
    "PrimaryLineReadFilter",
    "ProperlyPairedReadFilter",
    "ReadLengthEqualsCigarLengthReadFilter",
    "SecondOfPairReadFilter",
    "SeqIsStoredReadFilter",
    "ValidAlignmentEndReadFilter",
    "ValidAlignmentStartReadFilter",
];

#[cfg(test)]
mod tests {
    use super::*;

    fn mapped_read() -> BamRecord {
        BamRecord {
            reference_index: 0,
            alignment_start: 100,
            mapping_quality: 60,
            read_bases: vec![b'A'; 10],
            base_qualities: vec![30; 10],
            ..BamRecord::default()
        }
    }

    #[test]
    fn unmapped_is_three_criteria_not_one() {
        // The flag is clear, so a port that tested the flag alone would call this mapped.
        let mut read = mapped_read();
        read.reference_index = NO_ALIGNMENT_REFERENCE_INDEX;
        assert!(read::is_unmapped(&read));
        assert!(!mapped(&read));

        let mut read = mapped_read();
        read.alignment_start = NO_ALIGNMENT_START;
        assert!(read::is_unmapped(&read));

        let mut read = mapped_read();
        read.flags |= flags::READ_UNMAPPED;
        assert!(read::is_unmapped(&read));

        assert!(mapped(&mapped_read()));
    }

    #[test]
    fn first_of_pair_requires_pairing() {
        // 0x40 without 0x1 is not first-of-pair to GATK, whatever the flag says.
        let mut read = mapped_read();
        read.flags |= flags::FIRST_OF_PAIR;
        assert!(!first_of_pair(&read));

        read.flags |= flags::READ_PAIRED;
        assert!(first_of_pair(&read));
    }

    #[test]
    fn matching_bases_and_quals_counts_absent_qualities_as_zero() {
        let mut read = mapped_read();
        assert!(matching_bases_and_quals(&read));
        read.base_qualities.clear();
        assert!(!matching_bases_and_quals(&read));
    }

    #[test]
    fn mate_on_same_contig_or_no_mapped_mate() {
        let mut read = mapped_read();
        // Unpaired: kept.
        assert!(super::mate_on_same_contig_or_no_mapped_mate(&read));

        read.flags |= flags::READ_PAIRED;
        read.mate_reference_index = 1;
        read.mate_alignment_start = 50;
        assert!(!super::mate_on_same_contig_or_no_mapped_mate(&read));

        read.mate_reference_index = 0;
        assert!(super::mate_on_same_contig_or_no_mapped_mate(&read));
    }

    #[test]
    fn every_ported_name_resolves() {
        for name in PORTED {
            assert!(by_name(name).is_some(), "{name} is listed but not wired");
        }
    }
}
