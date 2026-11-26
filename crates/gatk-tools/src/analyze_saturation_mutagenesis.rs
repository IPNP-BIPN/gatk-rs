//! `AnalyzeSaturationMutagenesis`: reads over an ORF, counted by what each one turned out to be.
//!
//! What is ported is the census the tool writes and the decisions behind it: the quality trim, the
//! report type each read lands in, the flank test a variant has to pass, and the codon and
//! amino-acid names a variant row carries. The alignment itself is not ported.
//!
//! Ported from `org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis` in GATK 4.7.0.0;
//! `AnalyzeSaturationMutagenesis.java` is unchanged since 4.6.2.0.
//!
//! # This is the measured model, not the tool's own path
//!
//! The command line runs [`crate::analyze_saturation_mutagenesis_engine`], which carries its
//! own copy of each rule written against the reads. The functions here take the rules one at a
//! time, on plain values, so the `readfilter-conformance` dump
//! (`tests/data/analyze_saturation_mutagenesis.txt.gz`) can measure each of them alone. Where
//! the two copies differ (each difference is noted on its item below), the engine's is what the
//! tool does.
//!
//! # Reading this file without reading Rust
//!
//! - `&[u8]` is a borrowed, read-only view of a byte array; `u8` is an UNSIGNED byte, where
//!   Java's `byte` is signed, which matters wherever a quality is compared.
//! - `Option<T>` is "a value or nothing"; `?` inside a function returning `Option` answers
//!   nothing as soon as one step does. `Result<(), String>` is "nothing, or an error message".
//! - `usize` is an unsigned index type, so `saturating_sub` (subtract, but stop at zero) stands
//!   in for a Java subtraction that could go negative.
//! - a function's last expression, without a `;`, is its return value.

/// The classes a read or a pair is counted under.
///
/// The variants are NOT in the Java's declaration order, which is the order
/// `writeReportTypeCounts` walks (`for ( final ReportType reportType : ReportType.values() )`):
///
/// ```java
/// WILD_TYPE(null, "Wild type"),
/// CALLED_VARIANT(null, "Called variants"),
/// INCONSISTENT("inconsistent", "Inconsistent pair"),
/// IGNORED_MATE("ignoredMate", "Mate ignored"),
/// LOW_Q_VAR("lowQVar", "Low quality variation"),
/// NO_FLANK("noFlank", "Insufficient flank");
/// ```
///
/// [`category_block`] writes the rows in the order a [`Counts`] holds them instead. The
/// engine's `ReportType` follows the declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ReportType {
    Unmapped,
    LowQuality,
    Evaluable,
    WildType,
    LowQualityVariant,
    NoFlank,
    Inconsistent,
    IgnoredMate,
    CalledVariant,
}

impl ReportType {
    /// The label the census writes for each, which is not the enum's own name.
    ///
    /// `match` is a `switch` the compiler checks covers every variant.
    pub fn label(self) -> &'static str {
        match self {
            ReportType::Unmapped => "Unmapped Reads",
            ReportType::LowQuality => "LowQ Reads",
            ReportType::Evaluable => "Evaluable Reads",
            ReportType::WildType => "Wild type",
            ReportType::LowQualityVariant => "LowQ variant",
            ReportType::NoFlank => "Insufficient flank",
            ReportType::Inconsistent => "Inconsistent pair",
            ReportType::IgnoredMate => "Mate ignored",
            ReportType::CalledVariant => "Called variants",
        }
    }
}

/// The arguments that decide what a read is.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Arguments {
    /// `--min-q`: the quality a base needs to count as high quality.
    pub min_quality: u8,
    /// `--min-length`: how many consecutive high-quality bases a read needs.
    pub min_length: usize,
    /// `--min-flanking-length`: wild-type calls needed on each side of a variant.
    pub min_flanking_length: usize,
    /// `--min-mapq`: below this a read is counted UNMAPPED, not badly mapped.
    pub min_mapping_quality: i32,
    /// `--min-variant-obs`: how many observations a variant needs to be reported.
    pub min_variant_observations: u64,
    /// `--paired-mode`.
    pub paired_mode: bool,
    /// The negation of `--dont-ignore-disjoint-pairs`.
    pub ignore_disjoint_pairs: bool,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            min_quality: 30,
            min_length: 15,
            min_flanking_length: 2,
            min_mapping_quality: 4,
            min_variant_observations: 3,
            paired_mode: true,
            ignore_disjoint_pairs: true,
        }
    }
}

/// A half-open interval of read offsets, which is what the trim is: the Java's `Interval`,
/// without its constructor's refusal of a reversed interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
}

impl Interval {
    pub const NULL: Interval = Interval { start: 0, end: 0 };

    pub fn size(&self) -> usize {
        self.end.saturating_sub(self.start)
    }
}

/// `calculateQualityTrim`: the first and last runs of `min_length` bases at or above `min_q`.
///
/// The walk forwards stops at the END of the first such run and then steps back over it, so the
/// interval starts at the run's first base. A read with no such run at all yields the null
/// interval, whose size is zero and so under any minimum.
///
/// The backward walk here counts `end` from the length down, one past the Java's `readEnd`, so
/// `end - 1 + min_length` is the Java's `readEnd += minLength`. The Java compares a signed
/// `byte` (`if ( quals[readStart] < minQ )`), and this compares the unsigned `u8`, so the two
/// part on a quality above 127; the engine casts to `i8` first.
pub fn quality_trim(qualities: &[u8], arguments: &Arguments) -> Interval {
    let mut start = 0usize;
    let mut high = 0usize;
    while start < qualities.len() {
        if qualities[start] < arguments.min_quality {
            high = 0;
        } else {
            high += 1;
            if high == arguments.min_length {
                break;
            }
        }
        start += 1;
    }
    if start == qualities.len() {
        return Interval::NULL;
    }
    let start = start + 1 - arguments.min_length;

    let mut end = qualities.len();
    let mut high = 0usize;
    while end > 0 {
        if qualities[end - 1] < arguments.min_quality {
            high = 0;
        } else {
            high += 1;
            if high == arguments.min_length {
                break;
            }
        }
        end -= 1;
    }
    let end = end - 1 + arguments.min_length;
    Interval { start, end }
}

/// `calculateShortFragmentTrim`: the trim cut back to the fragment's own ends.
///
/// The fragment length is read from the record's TLEN, so a properly-paired read whose TLEN is
/// ZERO has its whole trim cut away and is counted LOW_QUALITY however good its bases are. That
/// is not a guard against a missing TLEN: it is the guard against a fragment shorter than the
/// read, applied to a fragment length of nothing.
///
/// Where the Java returns `Interval.NULL_INTERVAL` because too few calls are left, this returns
/// the cut-back interval itself; both are shorter than `min_length`, so [`classify_read`]
/// answers LOW_QUALITY for either. `unsigned_abs` is not Java's `Math.abs`, which leaves `int`'s
/// minimum negative.
pub fn fragment_trim(
    trim: Interval,
    read_length: usize,
    is_properly_paired: bool,
    is_reverse: bool,
    fragment_length: i32,
    arguments: &Arguments,
) -> Interval {
    if trim.size() < arguments.min_length || !is_properly_paired {
        return trim;
    }
    let fragment = fragment_length.unsigned_abs() as usize;
    if is_reverse {
        // The read has been reverse-complemented, so its beginning is what runs past the end.
        let minimum_start = read_length.saturating_sub(fragment);
        Interval {
            start: trim.start.max(minimum_start),
            end: trim.end.max(trim.start.max(minimum_start)),
        }
    } else {
        Interval {
            start: trim.start.min(fragment),
            end: trim.end.min(fragment),
        }
    }
}

/// Which report type a read lands in before any variant is looked at.
///
/// The order is the reference's: an unmapped, duplicate, vendor-failed or badly-mapped read is
/// UNMAPPED, and only then is the trim computed and measured against the minimum.
pub fn classify_read(
    is_unmapped: bool,
    is_duplicate: bool,
    fails_vendor_check: bool,
    mapping_quality: i32,
    trim: Interval,
    arguments: &Arguments,
) -> ReportType {
    if is_unmapped
        || is_duplicate
        || fails_vendor_check
        || mapping_quality < arguments.min_mapping_quality
    {
        return ReportType::Unmapped;
    }
    if trim.size() < arguments.min_length {
        return ReportType::LowQuality;
    }
    ReportType::Evaluable
}

/// Whether a variant has enough wild-type calls on each side of it to be called.
///
/// Both sides are measured inside the read's own evaluable span, so a variant at the very edge of
/// a trimmed read fails however much reference lies beyond it.
///
/// The reference's `hasCleanLeftFlank` clamps the flank at the amplicon's start, and the right
/// one clamps at its end; this test does not clamp:
///
/// ```java
/// Math.max(0, snvList.get(0).getRefIndex() - minFlankingLength) >= getFirstRefIndex();
/// ```
pub fn has_sufficient_flank(
    variant_offset: usize,
    coverage: Interval,
    arguments: &Arguments,
) -> bool {
    variant_offset >= coverage.start + arguments.min_flanking_length
        && variant_offset + arguments.min_flanking_length < coverage.end
}

// ================================================================================================
// The codons.
// ================================================================================================

/// The default `--codon-translation`, one amino-acid code per codon.
pub const DEFAULT_CODON_TRANSLATION: &str =
    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF";

/// The base order the translation string is indexed in: `A`, `C`, `G`, `T`.
pub const BASE_ORDER: [u8; 4] = *b"ACGT";

/// A codon's index into the translation string, which is its three bases in base-four.
///
/// A base that is not one of the four yields nothing, so an `N` in the reference has no codon.
/// `position` is Java's `indexOf`, answering `None` for -1, and the `?` after it returns that
/// `None` from the whole function.
pub fn codon_value(bases: &[u8]) -> Option<usize> {
    if bases.len() != 3 {
        return None;
    }
    let mut value = 0usize;
    for base in bases {
        let index = BASE_ORDER.iter().position(|known| known == base)?;
        value = value * 4 + index;
    }
    Some(value)
}

/// The amino acid a codon translates to: `codonTranslation.charAt(codonValue)`.
pub fn translate(bases: &[u8], translation: &str) -> Option<char> {
    codon_value(bases).and_then(|value| translation.chars().nth(value))
}

/// The refusal a translation string of the wrong length produces.
pub const TRANSLATION_LENGTH_MESSAGE: &str =
    "codon-translation string must contain exactly 64 characters";

/// The refusal an ORF whose length does not divide by three produces.
pub const ORF_LENGTH_MESSAGE: &str = "ORF length must be divisible by 3.";

/// The `onTraversalStart` refusal of a translation that is not 64 characters long.
pub fn check_translation(translation: &str) -> Result<(), String> {
    if translation.chars().count() != 64 {
        return Err(TRANSLATION_LENGTH_MESSAGE.to_string());
    }
    Ok(())
}

/// One interval of the ORF, one-based and inclusive as the argument writes them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OrfInterval {
    pub start: i32,
    pub end: i32,
}

/// `--orf`, parsed: `134-180,214-238`, with no spaces.
///
/// `collect()` into `Option<Vec<_>>` answers `None` if any one pair failed to parse.
pub fn parse_orf(text: &str) -> Option<Vec<OrfInterval>> {
    text.split(',')
        .map(|part| {
            let (start, end) = part.split_once('-')?;
            Some(OrfInterval {
                start: start.parse().ok()?,
                end: end.parse().ok()?,
            })
        })
        .collect()
}

/// The ORF's total length, which is what has to divide by three.
///
/// The intervals are spliced before translation, so a codon may straddle two of them.
pub fn orf_length(intervals: &[OrfInterval]) -> i32 {
    intervals
        .iter()
        .map(|interval| interval.end - interval.start + 1)
        .sum()
}

/// Two of `getExons`'s refusals: an end past the reference, and a length not divisible by three.
///
/// The Java's message ends with the whole `--orf` text, where this one ends with the offending
/// end coordinate alone:
///
/// ```java
/// throw new UserException("Found ORF end coordinate larger than reference length: " + orfCoords);
/// ```
pub fn check_orf(intervals: &[OrfInterval], reference_length: i32) -> Result<(), String> {
    if let Some(past) = intervals.iter().find(|i| i.end > reference_length) {
        return Err(format!(
            "Found ORF end coordinate larger than reference length: {}",
            past.end
        ));
    }
    if orf_length(intervals) % 3 != 0 {
        return Err(ORF_LENGTH_MESSAGE.to_string());
    }
    Ok(())
}

// ================================================================================================
// The census.
// ================================================================================================

/// The counts one category holds, as `(type, count)` pairs in the order they were first bumped.
///
/// The Java's `ReportTypeCounts` is an array indexed by ordinal, so its rows always come out in
/// declaration order; this list comes out in first-bumped order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Counts {
    pub entries: Vec<(ReportType, u64)>,
}

impl Counts {
    pub fn get(&self, kind: ReportType) -> u64 {
        self.entries
            .iter()
            .find(|(name, _)| *name == kind)
            .map_or(0, |(_, count)| *count)
    }

    pub fn total(&self) -> u64 {
        self.entries.iter().map(|(_, count)| count).sum()
    }

    pub fn bump(&mut self, kind: ReportType) {
        match self.entries.iter_mut().find(|(name, _)| *name == kind) {
            Some((_, count)) => *count += 1,
            None => self.entries.push((kind, 1)),
        }
    }
}

/// `new DecimalFormat("0.000")`, which is what every percentage in the census is written with.
///
/// It rounds half to EVEN, which is the Java default and not the half-up a naive port would give.
/// The engine formats through `gatk_annotation::decimal_format::DecimalFormat` instead.
pub fn percentage(part: u64, whole: u64) -> String {
    decimal_format(100.0 * part as f64 / whole as f64)
}

/// `DecimalFormat("0.000")` over one value, including its NaN and infinity spellings.
pub fn decimal_format(value: f64) -> String {
    if value.is_nan() {
        return "\u{fffd}".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "\u{221e}" } else { "-\u{221e}" }.to_string();
    }
    let scaled = value * 1000.0;
    let rounded = round_half_even(scaled) / 1000.0;
    format!("{rounded:.3}")
}

/// Half-to-even rounding, which is what `DecimalFormat` uses unless told otherwise.
///
/// The tie is detected on `value * 1000` after the multiplication has rounded, within an
/// epsilon, not on the double's exact decimal expansion. `100.0 * 1 / 8000` is
/// `0.01250000000000000069...`, a little above the tie, and scaled it rounds to exactly
/// `12.5`, so this answers `0.012`.
fn round_half_even(value: f64) -> f64 {
    let floor = value.floor();
    let difference = value - floor;
    if (difference - 0.5).abs() < f64::EPSILON * value.abs().max(1.0) {
        if (floor as i64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else {
        value.round()
    }
}

/// One line of the census: a depth of `>` markers, a label, a count and a percentage.
pub fn census_line(depth: usize, label: &str, count: u64, whole: u64) -> String {
    format!(
        "{}{label}:\t{count}\t{}%",
        ">".repeat(depth),
        percentage(count, whole)
    )
}

/// The three top lines of the census, whose denominator is every read.
pub fn top_lines(reads: &Counts) -> Vec<String> {
    let total = reads.total();
    let mut lines = vec![format!("Total Reads:\t{total}\t100.000%")];
    for kind in [
        ReportType::Unmapped,
        ReportType::LowQuality,
        ReportType::Evaluable,
    ] {
        lines.push(census_line(1, kind.label(), reads.get(kind), total));
    }
    lines
}

/// One category's block: its own line against the evaluable reads, then its rows against itself.
///
/// The overlapping category's line counts READS and so is twice its pair count, while the rows
/// beneath it count pairs: the two levels are not in the same unit.
pub fn category_block(
    label: &str,
    counts: &Counts,
    reads_in_category: u64,
    evaluable: u64,
) -> Vec<String> {
    let mut lines = vec![format!(
        ">>{label}:\t{reads_in_category}\t{}%",
        percentage(reads_in_category, evaluable)
    )];
    let total = counts.total();
    for (kind, count) in &counts.entries {
        if *count != 0 {
            lines.push(census_line(3, kind.label(), *count, total));
        }
    }
    lines
}

/// Whether a variant reaches the report at all: `entry.getCount() >= minVariantObservations`.
pub fn is_reported(observations: u64, arguments: &Arguments) -> bool {
    observations >= arguments.min_variant_observations
}
