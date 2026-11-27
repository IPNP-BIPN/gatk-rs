//! `AnalyzeSaturationMutagenesis`, end to end: the reads walked, the variants counted and the eight
//! reports (and the rejected-reads BAM) written.
//!
//! [`crate::analyze_saturation_mutagenesis`] carries the census rules the `readfilter-conformance`
//! dump measures one at a time. This module is the tool those rules live in, written as the
//! reference writes it: one pass over the reads in file order, primary lines only, with the pairs
//! found by adjacency; every read reported as a list of single-base differences from the amplicon
//! and the reference span it covered; and the reports built from the counts once the pass ends.
//!
//! What a command line observes and this reproduces:
//!
//!  * the traversal is the tool's own `getTransformedReadStream(PRIMARY_LINE)`, so the read
//!    filters a command line resolves are validated and never applied, and `-L` bounds nothing;
//!  * a read's variants are counted only if they are clean: every call at or above `--min-q`, one
//!    of `-ACGT`, and `--min-flanking-length` wild-type calls on each side;
//!  * the reports are written in `onTraversalSuccess`'s order, and the rejected reads, when asked
//!    for, carry an `XX` tag naming why and are written as the pass meets them, after any change
//!    `--find-large-deletions` made to their alignment;
//!  * every failure inside the pass reaches the handler wrapped as `Caught unexpected exception on
//!    read N`, where N is the number of reads counted so far and the name is the read held as the
//!    first of a pair at that moment.
//!
//! Ported from `org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis` (GATK 4.7.0.0;
//! `AnalyzeSaturationMutagenesis.java` is unchanged since 4.6.2.0).
//!
//! # The tool's state is static, and this is one `Run`
//!
//! The arguments, the amplicon and every counter the reference keeps are static fields of the
//! tool class. `reference` and `codonTracker` are assigned in `onTraversalStart`, but the counts
//! are initialised once per JVM and never reset:
//!
//! ```java
//! private static final ReportTypeCounts readCounts = new ReportTypeCounts();
//! ```
//!
//! so two runs in one JVM add up. A command line runs the tool once per JVM, so a [`Run`] that
//! starts from zero is what it observes. The static helpers that read `minQ`, `minLength` and
//! the other static arguments are methods here that read [`Settings`] instead.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, AsmError>` is "a `T` or the exception the reference would throw"; [`AsmError`]
//!   names the Java class. `?` after a call hands that exception straight back to the caller,
//!   which is Java's unchecked exception propagating out of the call site.
//! - `Option<T>` is "a `T` or nothing", the compiler-checked form of a reference that may be
//!   `null`. `Some(x)` holds a value and `None` is the `null`.
//! - `let Some(x) = y else { ... };` takes the value out of `y`, or runs the block, which must
//!   leave the function or loop. Java would write `if (y == null) { ... }` and then use `y`.
//! - `match` is a `switch` the compiler checks for completeness; an arm such as
//!   `Some(present) if present.ref_index == ref_index =>` matches only when the guard after `if`
//!   is true, and `_ =>` is `default:`.
//! - `x as i8` reinterprets a byte as signed, which is what Java's `byte` already is: the read's
//!   bases and qualities are stored here as unsigned `u8`, so every comparison that Java makes on
//!   a signed `byte` casts first. `as usize` turns a Java `int` into an index.
//! - `&mut self` is a method that may change the run; `&self` one that only reads it.
//! - `|x| ...` is a lambda. A function's last expression, without a `;`, is its return value.
//! - `wrapping_add` is Java's `int` addition, which wraps on overflow; Rust's `+` is allowed to
//!   stop the program instead.

use std::collections::{BTreeMap, HashMap};

use gatk_engine::java_format::format_decimals;
use gatk_engine::read;
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// `getToolName()`, which is what the rejected BAM's `@PG` record carries.
pub const TOOL_NAME: &str = "GATK AnalyzeSaturationMutagenesis";

/// The default `--codon-translation`.
pub const DEFAULT_CODON_TRANSLATION: &str =
    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF";

/// `GATKException`, which is what every failure inside the pass is wrapped in.
pub const GATK_EXCEPTION: &str = "org.broadinstitute.hellbender.exceptions.GATKException";

/// `NO_CALL`: the call an `Snv` carries on the side that has no base, the read's for a
/// deletion and the reference's for an insertion.
const NO_CALL: u8 = b'-';
/// `UPPERCASE_MASK`, `0xDF`: clearing bit 5 turns an ASCII lower-case letter into its upper
/// case, so a soft-masked `a` in the reference or the read compares equal to `A`.
const UPPERCASE_MASK: u8 = 0xDF;
/// `CodonTracker`'s row layout: one count per codon value (`0..64`), then the in-frame indels
/// (`NFS`) and the frame shifts (`FS`), the last two columns of `.codonCounts`.
const N_REGULAR_CODONS: usize = 64;
const FRAME_PRESERVING_INDEL_INDEX: usize = 64;
const FRAME_SHIFTING_INDEL_INDEX: usize = 65;
const CODON_COUNT_ROW_SIZE: usize = 66;
/// `CodonTracker.NO_FRAME_SHIFT_CODON`: what `findFrameShift` answers when the read stays in
/// frame.
const NO_FRAME_SHIFT_CODON: i32 = -1;

/// `LABEL_FOR_CODON_VALUE`: a codon value is its three bases in base four over `ACGT`, first
/// base most significant, so value 14 (`0b00_11_10`) is `ATG`.
const LABEL_FOR_CODON_VALUE: [&str; 64] = [
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA",
    "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC",
    "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG",
    "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT",
    "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT",
];

/// The tool's own arguments, at the values the command line gave them.
///
/// The static `@Argument` fields, one per command-line argument, with the reference's
/// defaults. The integer types are the Java fields' own (`int`, `long`), so the comparisons
/// and the overflow behave as the reference's do.
#[derive(Debug, Clone)]
pub struct Settings {
    pub min_q: i32,
    pub min_length: i32,
    pub min_flanking_length: i32,
    pub min_mapq: i32,
    /// `--orf`, kept as text: [`Run::start`] parses it, and most of its refusals quote it whole.
    pub orf: String,
    pub min_variant_observations: i64,
    pub find_large_deletions: bool,
    /// `--min-alt-length`: the clip a primary alignment needs at one end before
    /// `--find-large-deletions` looks at its `SA` tag.
    pub min_alt_length: i32,
    /// `--codon-translation`, as the UTF-16 units `charAt` indexes.
    pub codon_translation: Vec<u16>,
    pub paired_mode: bool,
    /// `--dont-ignore-disjoint-pairs`: combine the mates of a pair that do not overlap, rather
    /// than count the first of pair alone and reject the other as `ignoredMate`.
    pub dont_ignore_disjoint_pairs: bool,
    /// `--write-rejected-reads`: keep the rejected reads for the `.rejected.bam`.
    pub write_rejected_reads: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            min_q: 30,
            min_length: 15,
            min_flanking_length: 2,
            min_mapq: 4,
            orf: String::new(),
            min_variant_observations: 3,
            find_large_deletions: false,
            min_alt_length: 15,
            codon_translation: DEFAULT_CODON_TRANSLATION.encode_utf16().collect(),
            paired_mode: true,
            dont_ignore_disjoint_pairs: false,
            write_rejected_reads: false,
        }
    }
}

/// Why a run stopped: the reference exception's class and message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AsmError {
    /// A `UserException`, whose message the handler prints in its banner.
    User(String),
    /// Anything else, named by the class the reference throws.
    Internal {
        class: &'static str,
        message: String,
    },
}

impl AsmError {
    fn gatk(message: impl Into<String>) -> Self {
        AsmError::Internal {
            class: GATK_EXCEPTION,
            message: message.into(),
        }
    }

    fn user(message: impl Into<String>) -> Self {
        AsmError::User(message.into())
    }
}

/// An array index past the end, which Java raises as `ArrayIndexOutOfBoundsException`.
///
/// The message is the JVM's own wording for that exception. The tool lets that exception
/// reach the traversal's `catch`, which wraps it (see [`Run::traverse`]).
fn out_of_bounds(index: i64, length: usize) -> AsmError {
    AsmError::Internal {
        class: "java.lang.ArrayIndexOutOfBoundsException",
        message: format!("Index {index} out of bounds for length {length}"),
    }
}

/// `array[index]` with Java's bounds check: Rust's own indexing would stop the program instead
/// of raising an exception the traversal can wrap. `&[T]` is a borrowed, read-only view of an
/// array of any element type `T`; `T: Copy` limits it to plain values such as `u8` and `i32`.
fn at<T: Copy>(array: &[T], index: i64) -> Result<T, AsmError> {
    if index < 0 || index as usize >= array.len() {
        return Err(out_of_bounds(index, array.len()));
    }
    Ok(array[index as usize])
}

// ================================================================================================
// Report types.
// ================================================================================================

/// `ReportType`, in declaration order, which is the order the census lists them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportType {
    Unmapped,
    LowQuality,
    Evaluable,
    WildType,
    CalledVariant,
    Inconsistent,
    IgnoredMate,
    LowQVar,
    NoFlank,
}

/// `ReportType.values()`, which is the order `writeReportTypeCounts` walks. The census module's
/// enum has its own order; this one is the Java declaration's.
const REPORT_TYPES: [ReportType; 9] = [
    ReportType::Unmapped,
    ReportType::LowQuality,
    ReportType::Evaluable,
    ReportType::WildType,
    ReportType::CalledVariant,
    ReportType::Inconsistent,
    ReportType::IgnoredMate,
    ReportType::LowQVar,
    ReportType::NoFlank,
];

impl ReportType {
    /// `attributeValue`: the `XX` value a rejected read is tagged with, or none for the types
    /// that are not rejections.
    pub fn attribute_value(self) -> Option<&'static str> {
        match self {
            ReportType::Unmapped => Some("unmapped"),
            ReportType::LowQuality => Some("lowQ"),
            ReportType::Inconsistent => Some("inconsistent"),
            ReportType::IgnoredMate => Some("ignoredMate"),
            ReportType::LowQVar => Some("lowQVar"),
            ReportType::NoFlank => Some("noFlank"),
            ReportType::Evaluable | ReportType::WildType | ReportType::CalledVariant => None,
        }
    }

    /// `label`: the name the `.readCounts` census writes, which is not the constant's.
    /// `&'static str` is a string fixed into the program, which is Java's string literal.
    pub fn label(self) -> &'static str {
        match self {
            ReportType::Unmapped => "Unmapped Reads",
            ReportType::LowQuality => "LowQ Reads",
            ReportType::Evaluable => "Evaluable Reads",
            ReportType::WildType => "Wild type",
            ReportType::CalledVariant => "Called variants",
            ReportType::Inconsistent => "Inconsistent pair",
            ReportType::IgnoredMate => "Mate ignored",
            ReportType::LowQVar => "Low quality variation",
            ReportType::NoFlank => "Insufficient flank",
        }
    }

    /// `Enum.ordinal()`: the position in `REPORT_TYPES`. Rust has no built-in ordinal for an
    /// enum used this way, so the table is searched.
    fn ordinal(self) -> usize {
        REPORT_TYPES
            .iter()
            .position(|kind| *kind == self)
            .expect("every type is listed")
    }
}

/// `ReportTypeCounts`.
///
/// `long[] counts = new long[ReportType.values().length]`, indexed by `ordinal()`.
#[derive(Debug, Clone, Default)]
struct ReportTypeCounts {
    counts: [i64; 9],
}

impl ReportTypeCounts {
    fn bump(&mut self, kind: ReportType) {
        self.counts[kind.ordinal()] += 1;
    }
    fn get(&self, kind: ReportType) -> i64 {
        self.counts[kind.ordinal()]
    }
    fn total(&self) -> i64 {
        self.counts.iter().sum()
    }
}

// ================================================================================================
// Intervals, SNVs and the reference.
// ================================================================================================

/// `Interval`: 0-based, half-open.
///
/// The constructor's refusal of a negative start or a reversed interval is kept, because it is
/// a `GATKException` the traversal wraps; the fields are `int`, as in the Java.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Interval {
    start: i32,
    end: i32,
}

impl Interval {
    const NULL: Interval = Interval { start: 0, end: 0 };

    fn new(start: i32, end: i32) -> Result<Interval, AsmError> {
        if start < 0 || end < start {
            return Err(AsmError::gatk(format!("Illegal interval: [{start},{end})")));
        }
        Ok(Interval { start, end })
    }

    fn size(&self) -> i32 {
        self.end - self.start
    }

    fn overlap_length(&self, that: &Interval) -> i32 {
        self.end.min(that.end) - self.start.max(that.start)
    }
}

/// `SNV`: one base call that differs from the amplicon. The quality is not part of its identity.
///
/// `quality` is a raw `byte`, compared as signed wherever the Java compares it.
#[derive(Debug, Clone, Copy)]
struct Snv {
    ref_index: i32,
    ref_call: u8,
    variant_call: u8,
    quality: u8,
}

impl Snv {
    fn key(&self) -> (i32, i8, i8) {
        (self.ref_index, self.ref_call as i8, self.variant_call as i8)
    }

    fn same_as(&self, that: &Snv) -> bool {
        self.key() == that.key()
    }

    fn text(&self) -> String {
        format!(
            "{}:{}>{}",
            self.ref_index + 1,
            self.ref_call as char,
            self.variant_call as char
        )
    }
}

/// `IntervalCounter`: molecules counted by where they start and how far they reach.
///
/// A triangular matrix indexed first by start and then by length: row `start` has room for
/// every length from `0` to `reference_length - start`. `Vec<Vec<i64>>` is Java's `long[][]`,
/// a list of rows of different lengths.
struct IntervalCounter {
    counts: Vec<Vec<i64>>,
}

impl IntervalCounter {
    fn new(reference_length: usize) -> Self {
        IntervalCounter {
            counts: (0..reference_length)
                .map(|row| vec![0; reference_length - row + 1])
                .collect(),
        }
    }

    fn add_count(&mut self, span: Interval) -> Result<(), AsmError> {
        if span.end as usize > self.counts.len() {
            return Err(AsmError::gatk(format!(
                "illegal span: [{},{})",
                span.start, span.end
            )));
        }
        let row = self
            .counts
            .get_mut(span.start as usize)
            .ok_or_else(|| out_of_bounds(span.start as i64, 0))?;
        let length = row.len();
        let cell = row
            .get_mut((span.end - span.start) as usize)
            .ok_or_else(|| out_of_bounds((span.end - span.start) as i64, length))?;
        *cell += 1;
        Ok(())
    }

    /// `countSpanners`: the molecules whose span contains `[start, end)`, which is every row that
    /// starts at or before `start` and every length in that row that reaches `end`.
    fn count_spanners(&self, start: i32, end: i32) -> Result<i64, AsmError> {
        if start < 0 || end < start || end as usize > self.counts.len() {
            return Err(AsmError::gatk(format!("illegal span: [{start},{end})")));
        }
        let mut total = 0;
        for row_index in 0..=start {
            let row = self
                .counts
                .get(row_index as usize)
                .ok_or_else(|| out_of_bounds(row_index as i64, self.counts.len()))?;
            let mut span = (end - row_index) as usize;
            while span < row.len() {
                total += row[span];
                span += 1;
            }
        }
        Ok(total)
    }
}

/// `Reference`: the amplicon, upper-cased, and what the molecules covered of it.
///
/// `coverage` counts molecules per reference base; `coverage_size_histogram` counts molecules
/// by how many reference bases they covered, which is `.coverageLengthCounts`.
struct Reference {
    sequence: Vec<u8>,
    coverage: Vec<i64>,
    coverage_size_histogram: Vec<i64>,
    interval_counter: IntervalCounter,
}

impl Reference {
    fn new(sequence: Vec<u8>) -> Self {
        let length = sequence.len();
        Reference {
            sequence,
            coverage: vec![0; length],
            coverage_size_histogram: vec![0; length + 1],
            interval_counter: IntervalCounter::new(length),
        }
    }

    /// `updateCoverage`: one more molecule over each base of each interval, and one more in the
    /// histogram bucket for the sum of the intervals' sizes, which it also returns.
    fn update_coverage(&mut self, intervals: &[Interval]) -> Result<i32, AsmError> {
        let mut coverage_length = 0;
        for interval in intervals {
            coverage_length += interval.end - interval.start;
            for index in interval.start..interval.end {
                let length = self.coverage.len();
                *self
                    .coverage
                    .get_mut(index as usize)
                    .ok_or_else(|| out_of_bounds(index as i64, length))? += 1;
            }
        }
        let length = self.coverage_size_histogram.len();
        *self
            .coverage_size_histogram
            .get_mut(coverage_length as usize)
            .ok_or_else(|| out_of_bounds(coverage_length as i64, length))? += 1;
        Ok(coverage_length)
    }
}

/// The amplicon a reference's single contig holds, refused the way the reference refuses it.
///
/// `contigs` is how many sequences the reference declares; `bases` are the first one's, already
/// upper-cased and with IUPAC codes flattened to `N` the way `queryAndPrefetch` returns them.
///
/// The reference upper-cases in place (`refSeq[idx] &= UPPERCASE_MASK`) inside the `switch`, and
/// refuses at the first base that is not then one of `ACGT`. `iter_mut()` walks the copy with a
/// writable reference to each byte, and `*base` is the byte it points at.
pub fn amplicon(contigs: usize, bases: &[u8]) -> Result<Vec<u8>, AsmError> {
    if contigs != 1 {
        return Err(AsmError::user(format!(
            "Expecting a reference with a single contig. The supplied reference has {contigs} contigs."
        )));
    }
    let mut sequence = bases.to_vec();
    for base in sequence.iter_mut() {
        *base &= UPPERCASE_MASK;
        if !matches!(*base, b'A' | b'C' | b'G' | b'T') {
            return Err(AsmError::user(
                "Reference sequence contains something other than A, C, G, and T.",
            ));
        }
    }
    Ok(sequence)
}

// ================================================================================================
// Codons.
// ================================================================================================

/// `CodonVariationType`. Only equality is ever asked of it, so its declaration order does not
/// matter here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CodonVariationType {
    Frameshift,
    Insertion,
    Deletion,
    Modification,
}

/// `CodonVariation`: one altered codon. `codon_value` is `-1` for a frame shift or a deletion,
/// as the Java's factory methods set it.
#[derive(Debug, Clone, Copy)]
struct CodonVariation {
    codon_id: i32,
    codon_value: i32,
    kind: CodonVariationType,
}

impl CodonVariation {
    fn frameshift(codon_id: i32) -> Self {
        CodonVariation {
            codon_id,
            codon_value: -1,
            kind: CodonVariationType::Frameshift,
        }
    }
    fn insertion(codon_id: i32, codon_value: i32) -> Self {
        CodonVariation {
            codon_id,
            codon_value,
            kind: CodonVariationType::Insertion,
        }
    }
    fn deletion(codon_id: i32) -> Self {
        CodonVariation {
            codon_id,
            codon_value: -1,
            kind: CodonVariationType::Deletion,
        }
    }
    fn modification(codon_id: i32, codon_value: i32) -> Self {
        CodonVariation {
            codon_id,
            codon_value,
            kind: CodonVariationType::Modification,
        }
    }
}

/// `CodonTracker.isStop`: `TAA`, `TAG` and `TGA` as codon values.
///
/// ```java
/// final int STOP_OCH = 0x30;
/// final int STOP_AMB = 0x32;
/// final int STOP_OPA = 0x38;
/// ```
fn is_stop(codon_value: i32) -> bool {
    codon_value == 0x30 || codon_value == 0x32 || codon_value == 0x38
}

/// `"ACGT".indexOf(base)`.
///
/// The `-1` for anything else is kept: OR-ed into a codon value it sets every bit, as the Java's
/// `| "ACGT".indexOf(...)` does. The callers only reach it with `ACGT`.
fn base_value(base: u8) -> i32 {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => -1,
    }
}

/// `String.split(separator)` with Java's rules: trailing empty strings are dropped, and a string
/// with no separator in it is returned whole.
fn java_split(text: &str, separator: char) -> Vec<&str> {
    if !text.contains(separator) {
        return vec![text];
    }
    let mut parts: Vec<&str> = text.split(separator).collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }
    parts
}

/// `Integer.parseInt`.
///
/// An optional sign and ASCII digits, nothing else; out-of-range values fail as they do in Java.
/// `strip_prefix` drops one leading `+` or `-` when there is one, and `.ok()` turns Rust's parse
/// error into `None`, the `NumberFormatException` the callers catch.
fn parse_int(text: &str) -> Option<i32> {
    if text.is_empty() {
        return None;
    }
    let digits = text.strip_prefix(['+', '-']).unwrap_or(text);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    text.parse().ok()
}

/// `CodonTracker`: the ORF, its wild-type codons and what was observed at each.
///
/// `codon_counts` is `long[nCodons][CODON_COUNT_ROW_SIZE]`, written out as `.codonCounts`;
/// `ref_codon_values` is the wild-type codon value of each codon of the ORF.
struct CodonTracker {
    reference: Vec<u8>,
    exons: Vec<Interval>,
    codon_counts: Vec<[i64; CODON_COUNT_ROW_SIZE]>,
    ref_codon_values: Vec<i32>,
}

impl CodonTracker {
    /// The `CodonTracker` constructor: the ORF parsed and checked, then the reference read into
    /// codons, in the Java's order, so the first refusal raised is the Java's.
    fn new(orf: &str, reference: &[u8]) -> Result<Self, AsmError> {
        let exons = Self::exons(orf, reference.len() as i32)?;
        let codons = (exons.iter().map(Interval::size).sum::<i32>() / 3) as usize;
        let ref_codon_values = Self::parse_reference_into_codons(reference, &exons)?;
        Ok(CodonTracker {
            reference: reference.to_vec(),
            exons,
            codon_counts: vec![[0; CODON_COUNT_ROW_SIZE]; codons],
            ref_codon_values,
        })
    }

    /// `getExons`: `--orf`'s one-based inclusive pairs as zero-based half-open intervals.
    ///
    /// The checks run per pair in the Java's order (shape, start an integer and at least one, end
    /// an integer, end not before start, end inside the reference), and the sortedness check runs
    /// inside the loop after each pair is added, so an unsorted list is refused before a later
    /// malformed pair is read. Every per-pair message but the 1-based one carries the whole
    /// `--orf` text.
    fn exons(orf: &str, reference_length: i32) -> Result<Vec<Interval>, AsmError> {
        let mut exons: Vec<Interval> = Vec::new();
        for pair in java_split(orf, ',') {
            let coordinates = java_split(pair, '-');
            if coordinates.len() != 2 {
                return Err(AsmError::user(format!(
                    "Can't interpret ORF as list of pairs of coords: {orf}"
                )));
            }
            let not_integers =
                || AsmError::user(format!("Can't interpret ORF coords as integers: {orf}"));
            let start = parse_int(coordinates[0]).ok_or_else(not_integers)?;
            if start < 1 {
                return Err(AsmError::user("Coordinates of ORF are 1-based."));
            }
            let end = parse_int(coordinates[1]).ok_or_else(not_integers)?;
            if end < start {
                return Err(AsmError::user(format!(
                    "Found ORF end coordinate less than start: {orf}"
                )));
            }
            if end > reference_length {
                return Err(AsmError::user(format!(
                    "Found ORF end coordinate larger than reference length: {orf}"
                )));
            }
            exons.push(Interval::new(start - 1, end)?);
            for index in 1..exons.len() {
                if exons[index - 1].end >= exons[index].start {
                    return Err(AsmError::user(format!(
                        "ORF coordinates are not sorted: {orf}"
                    )));
                }
            }
        }
        let length: i32 = exons.iter().map(Interval::size).sum();
        if length % 3 != 0 {
            return Err(AsmError::user("ORF length must be divisible by 3."));
        }
        Ok(exons)
    }

    /// `parseReferenceIntoCodons`: the wild-type codon values, two bits per base, with a stop
    /// codon anywhere but last refused.
    fn parse_reference_into_codons(
        reference: &[u8],
        exons: &[Interval],
    ) -> Result<Vec<i32>, AsmError> {
        let codons = (exons.iter().map(Interval::size).sum::<i32>() / 3) as usize;
        let mut values = vec![0; codons];
        let mut codon_id = 0usize;
        let mut phase = 0;
        let mut value = 0;
        for exon in exons {
            for index in exon.start..exon.end {
                value = (value << 2) | base_value(reference[index as usize]);
                phase += 1;
                if phase == 3 {
                    if is_stop(value) && codon_id != codons - 1 {
                        return Err(AsmError::user(format!(
                            "There is an upstream stop codon at reference index {}.",
                            index + 1
                        )));
                    }
                    values[codon_id] = value;
                    value = 0;
                    phase = 0;
                    codon_id += 1;
                }
            }
        }
        // The start and stop codons are only warned about, on the log.
        Ok(values)
    }

    /// `isExonic`: whether a reference index lies inside one of the (sorted) exons.
    fn is_exonic(&self, index: i32) -> bool {
        for exon in &self.exons {
            if exon.start > index {
                return false;
            }
            if exon.end > index {
                return true;
            }
        }
        false
    }

    /// `exonicBaseIndex`: how many ORF bases lie before a reference index. For an intronic index
    /// that is the first base of the next exon, so dividing by three gives the next codon.
    fn exonic_base_index(&self, index: i32) -> i32 {
        let mut count = 0;
        for exon in &self.exons {
            if index >= exon.end {
                count += exon.size();
            } else {
                if index > exon.start {
                    count += index - exon.start;
                }
                break;
            }
        }
        count
    }

    /// `findFrameShift`: the codon at which the read's exonic indels leave it out of frame, or
    /// `NO_FRAME_SHIFT_CODON` when the deletions and insertions come back into frame. `lead_lag`
    /// counts inserted bases up and deleted bases down, and resets at plus or minus three.
    fn find_frame_shift(&self, snvs: &[Snv]) -> i32 {
        let mut codon_id = NO_FRAME_SHIFT_CODON;
        let mut lead_lag = 0;
        for snv in snvs {
            if self.is_exonic(snv.ref_index) {
                if snv.variant_call == NO_CALL {
                    if lead_lag == 0 {
                        codon_id = self.exonic_base_index(snv.ref_index) / 3;
                    }
                    lead_lag -= 1;
                    if lead_lag == -3 {
                        lead_lag = 0;
                    }
                } else if snv.ref_call == NO_CALL {
                    if lead_lag == 0 {
                        codon_id = self.exonic_base_index(snv.ref_index) / 3;
                    }
                    lead_lag += 1;
                    if lead_lag == 3 {
                        lead_lag = 0;
                    }
                }
                if lead_lag == 0 {
                    codon_id = NO_FRAME_SHIFT_CODON;
                }
            }
        }
        codon_id
    }

    /// `encodeSNVsAsCodons`: a read's base differences as codon differences.
    ///
    /// From each exonic variant, the codon it falls in is rebuilt base by base from the read
    /// (variant bases) and the amplicon (the rest), two bits at a time, until the read is back in
    /// frame at a codon boundary; each finished codon that differs from the wild type is one
    /// variation. Written as the Java is, step for step, because the early returns (a stop codon,
    /// the last codon, the end of the amplicon) and the `codon_id -= 1` after an insertion decide
    /// which codons the `.variantCounts` row names.
    fn encode_snvs_as_codons(&self, snvs: &[Snv]) -> Result<Vec<CodonVariation>, AsmError> {
        let mut variations = Vec::new();
        let mut next = 0usize;
        // `snvs` is walked through one index, `next`, shared by this search and the advance inside
        // the loop below, as the Java shares one `snvIterator`. The first variant taken must be
        // exonic and not an insertion before the ORF's first base.
        let mut snv: Option<Snv> = None;
        let orf_start = self.exons[0].start;
        while next < snvs.len() {
            let test = snvs[next];
            next += 1;
            if (test.ref_index != orf_start || test.ref_call != NO_CALL)
                && self.is_exonic(test.ref_index)
            {
                snv = Some(test);
                break;
            }
        }

        let mut frame_shift_codon_id = self.find_frame_shift(snvs);
        let last_exon_end = self.exons[self.exons.len() - 1].end;
        let codon_total = self.ref_codon_values.len() as i32;
        // `while let Some(current) = snv` loops while `snv` holds a variant: `while ( snv != null
        // )`.
        while let Some(current) = snv {
            let mut ref_index = current.ref_index;
            if ref_index >= last_exon_end {
                break;
            }
            let mut exon_cursor = 0usize;
            let mut current_exon = None;
            while exon_cursor < self.exons.len() {
                let test = self.exons[exon_cursor];
                exon_cursor += 1;
                if test.start <= ref_index && test.end > ref_index {
                    current_exon = Some(test);
                    break;
                }
            }
            // `ok_or_else` turns `None` (no exon found) into the Java's `GATKException`, and `?`
            // throws it.
            let mut current_exon = current_exon.ok_or_else(|| {
                AsmError::gatk("can't find current exon, even though refIndex should be exonic.")
            })?;

            let mut codon_id = self.exonic_base_index(ref_index);
            let mut codon_phase = codon_id % 3;
            codon_id /= 3;

            let mut codon_value = at(&self.ref_codon_values, codon_id as i64)?;
            if codon_phase == 0 {
                codon_value = 0;
            } else if codon_phase == 1 {
                codon_value >>= 4;
            } else {
                codon_value >>= 2;
            }

            let mut lead_lag = 0;
            loop {
                let mut codon_value_altered = false;
                let mut bump_ref_index = false;
                // The Java's `if` is the other way round; the arm below with the guard is its
                // `else`, the variant sitting at this reference base.
                match snv {
                    Some(present) if present.ref_index == ref_index => {
                        if present.variant_call == NO_CALL {
                            if codon_id == frame_shift_codon_id {
                                variations.push(CodonVariation::frameshift(codon_id));
                                frame_shift_codon_id = NO_FRAME_SHIFT_CODON;
                            }
                            lead_lag -= 1;
                            if lead_lag == -3 {
                                variations.push(CodonVariation::deletion(codon_id));
                                codon_id += 1;
                                if codon_id == codon_total {
                                    return Ok(variations);
                                }
                                lead_lag = 0;
                            }
                            bump_ref_index = true;
                        } else if present.ref_call == NO_CALL {
                            lead_lag += 1;
                            codon_value = (codon_value << 2) | base_value(present.variant_call);
                            codon_value_altered = true;
                        } else {
                            codon_value = (codon_value << 2) | base_value(present.variant_call);
                            codon_value_altered = true;
                            bump_ref_index = true;
                        }
                        snv = None;
                        while next < snvs.len() {
                            let test = snvs[next];
                            next += 1;
                            if test.ref_index >= last_exon_end || self.is_exonic(test.ref_index) {
                                snv = Some(test);
                                break;
                            }
                        }
                    }
                    // `snv == null || snv.getRefIndex() != refIndex`: no variant at this base, so
                    // the wild-type base goes into the codon.
                    _ => {
                        let base = at(&self.reference, ref_index as i64)?;
                        codon_value = (codon_value << 2) | base_value(base);
                        codon_value_altered = true;
                        bump_ref_index = true;
                    }
                }
                if bump_ref_index {
                    ref_index += 1;
                    // `exon_cursor` is where the exon search above stopped, so this is the Java's
                    // `exonIterator.next()` continuing the same iterator.
                    if ref_index == current_exon.end && exon_cursor < self.exons.len() {
                        current_exon = self.exons[exon_cursor];
                        exon_cursor += 1;
                        ref_index = current_exon.start;
                    }
                    if ref_index as usize == self.reference.len() {
                        return Ok(variations);
                    }
                }
                if codon_value_altered {
                    codon_phase += 1;
                    if codon_phase == 3 {
                        if codon_id == frame_shift_codon_id {
                            variations.push(CodonVariation::frameshift(codon_id));
                            frame_shift_codon_id = NO_FRAME_SHIFT_CODON;
                        }
                        if lead_lag >= 3 {
                            variations.push(CodonVariation::insertion(codon_id, codon_value));
                            lead_lag -= 3;
                            codon_id -= 1;
                        } else if codon_value != at(&self.ref_codon_values, codon_id as i64)? {
                            variations.push(CodonVariation::modification(codon_id, codon_value));
                        }
                        if is_stop(codon_value) {
                            return Ok(variations);
                        }
                        codon_id += 1;
                        if codon_id == codon_total {
                            return Ok(variations);
                        }
                        codon_phase = 0;
                        codon_value = 0;
                    }
                }
                // The Java's `do { ... } while ( leadLag != 0 || codonPhase != 0 );`: Rust has no
                // do-while, so the loop runs at least once and the condition is tested here, at its
                // end.
                if lead_lag == 0 && codon_phase == 0 {
                    break;
                }
            }
        }
        Ok(variations)
    }

    /// `reportVariantCodonCounts`: one count per codon the molecule fully covers, the observed
    /// value where a variation lies and the wild type elsewhere.
    ///
    /// `starting` rounds the first covered ORF base up to a codon boundary and `ending` rounds the
    /// end down, so a partly covered codon is not counted. A codon with any in-frame indel counts
    /// once in `NFS`, however many indels it holds.
    fn report_variant_codon_counts(
        &mut self,
        coverage: Interval,
        variations: &[CodonVariation],
    ) -> Result<(), AsmError> {
        let starting = (self.exonic_base_index(coverage.start) + 2) / 3;
        let ending = self.exonic_base_index(coverage.end) / 3;
        let mut cursor = 0usize;
        // The Java's `hasNext() ? next() : null` over the variations, as a lambda that advances
        // `cursor`. It takes the cursor as `&mut usize`, a writable borrow of the local, so the
        // loop below can still read `cursor` between calls.
        let next = |cursor: &mut usize| {
            let found = variations.get(*cursor).copied();
            if found.is_some() {
                *cursor += 1;
            }
            found
        };
        let mut variation = next(&mut cursor);
        for codon_id in starting..ending {
            while let Some(present) = variation {
                if present.codon_id >= codon_id {
                    break;
                }
                variation = next(&mut cursor);
            }
            let row_length = self.codon_counts.len();
            let row = self
                .codon_counts
                .get_mut(codon_id as usize)
                .ok_or_else(|| out_of_bounds(codon_id as i64, row_length))?;
            match variation {
                Some(present) if present.codon_id == codon_id => {
                    let mut frame_preserving_indel = false;
                    let mut current = Some(present);
                    while let Some(here) = current {
                        match here.kind {
                            CodonVariationType::Frameshift => {
                                row[FRAME_SHIFTING_INDEL_INDEX] += 1;
                            }
                            CodonVariationType::Deletion | CodonVariationType::Insertion => {
                                frame_preserving_indel = true;
                            }
                            CodonVariationType::Modification => {
                                row[here.codon_value as usize] += 1;
                            }
                        }
                        current = next(&mut cursor);
                        // `is_none_or` is true for `None` or for a value the test accepts: the
                        // Java's `while ( codonVariation != null && codonVariation.getCodonId() ==
                        // codonId )`, negated.
                        if current.is_none_or(|c| c.codon_id != codon_id) {
                            break;
                        }
                    }
                    variation = current;
                    if frame_preserving_indel {
                        row[FRAME_PRESERVING_INDEL_INDEX] += 1;
                    }
                }
                _ => {
                    let value = self.ref_codon_values[codon_id as usize];
                    row[value as usize] += 1;
                }
            }
        }
        Ok(())
    }

    /// `reportWildCodonCounts`: the wild-type value of every codon the molecule fully covers.
    fn report_wild_codon_counts(&mut self, coverage: Interval) -> Result<(), AsmError> {
        let starting = (self.exonic_base_index(coverage.start) + 2) / 3;
        let ending = self.exonic_base_index(coverage.end) / 3;
        for codon_id in starting..ending {
            let value = at(&self.ref_codon_values, codon_id as i64)?;
            self.codon_counts[codon_id as usize][value as usize] += 1;
        }
        Ok(())
    }
}

/// `CodonVariationGroup`: codon variations that one HGVS term describes.
///
/// `'a` is a lifetime: the group borrows the wild-type codon values and the translation from
/// the tracker and the settings rather than copying them, and the compiler checks that it does
/// not outlive them. `alt_calls` is the Java's `StringBuilder` of amino-acid `char`s.
struct CodonVariationGroup<'a> {
    ref_codon_values: &'a [i32],
    translation: &'a [u16],
    alt_calls: Vec<u16>,
    starting_codon: i32,
    ending_codon: i32,
    is_frame_shift: bool,
    ins_count: i32,
    del_count: i32,
    sub_count: i32,
}

impl<'a> CodonVariationGroup<'a> {
    fn new(ref_codon_values: &'a [i32], translation: &'a [u16], first: CodonVariation) -> Self {
        let mut group = CodonVariationGroup {
            ref_codon_values,
            translation,
            alt_calls: Vec::new(),
            starting_codon: first.codon_id,
            ending_codon: first.codon_id,
            is_frame_shift: false,
            ins_count: 0,
            del_count: 0,
            sub_count: 0,
        };
        match first.kind {
            CodonVariationType::Frameshift => group.is_frame_shift = true,
            CodonVariationType::Insertion => {
                group.ins_count = 1;
                group
                    .alt_calls
                    .push(translation[first.codon_value as usize]);
            }
            CodonVariationType::Deletion => group.del_count = 1,
            CodonVariationType::Modification => {
                group.sub_count = 1;
                group
                    .alt_calls
                    .push(translation[first.codon_value as usize]);
            }
        }
        group
    }

    /// `isEmpty`: a group begun by a frame shift that has gathered nothing but synonymous codons.
    fn is_empty(&self) -> bool {
        self.sub_count + self.ins_count + self.del_count == 0
    }

    /// `addVariation`: whether the variation joins this group. A skipped codon ends a group that is
    /// not a frame shift; a frame shift always starts a new one; a synonymous modification ends a
    /// non-frame-shift group and is otherwise only appended once the group holds something.
    fn add(&mut self, variation: CodonVariation) -> bool {
        let codon_id = variation.codon_id;
        if codon_id > self.ending_codon + 1 && !self.is_frame_shift {
            return false;
        }
        match variation.kind {
            CodonVariationType::Frameshift => return false,
            CodonVariationType::Insertion => {
                self.ins_count += 1;
                self.alt_calls
                    .push(self.translation[variation.codon_value as usize]);
                if self.is_frame_shift && self.is_empty() {
                    self.starting_codon = codon_id;
                }
            }
            CodonVariationType::Deletion => self.del_count += 1,
            CodonVariationType::Modification => {
                let aa = self.translation[variation.codon_value as usize];
                if aa == self.translation[self.ref_codon_values[codon_id as usize] as usize] {
                    if self.is_frame_shift {
                        if !self.is_empty() {
                            self.alt_calls.push(aa);
                        }
                    } else {
                        return false;
                    }
                } else {
                    if self.is_frame_shift && self.is_empty() {
                        self.starting_codon = codon_id;
                    }
                    self.sub_count += 1;
                    self.alt_calls.push(aa);
                }
            }
        }
        self.ending_codon = codon_id;
        true
    }

    fn amino_acid(&self, codon: i32) -> String {
        utf16(&[self.translation[self.ref_codon_values[codon as usize] as usize]])
    }

    /// `asHGVSString`: the group in HGVS protein nomenclature.
    ///
    /// It takes `&mut self` because, like the Java, it moves `starting_codon` back one for a pure
    /// insertion, and the moved value is what the rest of the string reads.
    fn hgvs(&mut self) -> String {
        let alts = utf16(&self.alt_calls);
        let mut text = String::new();
        if self.is_frame_shift && !self.alt_calls.is_empty() {
            text.push_str(&self.amino_acid(self.starting_codon));
            text.push_str(&(self.starting_codon + 1).to_string());
            text.push_str(&utf16(&self.alt_calls[..1]));
            let length = self.ending_codon - self.starting_codon + 1;
            if length > 1 {
                text.push_str("fs*");
                if self.alt_calls[self.alt_calls.len() - 1] == u16::from(b'X') {
                    text.push_str(&length.to_string());
                } else {
                    text.push('?');
                }
            }
        } else {
            if self.ins_count != 0 && self.del_count == 0 && self.sub_count == 0 {
                self.starting_codon -= 1;
            }
            text.push_str(&self.amino_acid(self.starting_codon));
            text.push_str(&(self.starting_codon + 1).to_string());
            if self.starting_codon != self.ending_codon {
                text.push('_');
                text.push_str(&self.amino_acid(self.ending_codon));
                text.push_str(&(self.ending_codon + 1).to_string());
            }
            if self.sub_count == 0 && self.ins_count == 0 {
                text.push_str("del");
            } else if self.sub_count == 0 && self.del_count == 0 {
                text.push_str("ins");
            } else if self.sub_count + self.del_count + self.ins_count > 1 {
                text.push_str("insdel");
            }
            text.push_str(&alts);
        }
        text
    }
}

/// A Java `char` sequence back to text. The translation string is kept as UTF-16 units so that
/// `charAt` indexes the same unit the Java does; `from_utf16_lossy` replaces a lone surrogate
/// with U+FFFD rather than failing.
fn utf16(units: &[u16]) -> String {
    String::from_utf16_lossy(units)
}

// ================================================================================================
// Read reports.
// ================================================================================================

/// `ReadReport`: the reference a read (or a pair) covered, and how its calls differed.
///
/// `snvs` is `None` for a pair whose mates disagree over their overlap, which is the reference's
/// null list.
///
/// `#[derive(Default)]` gives an empty coverage and `snvs: None`; every report the run builds
/// sets both fields itself.
#[derive(Debug, Clone, Default)]
struct ReadReport {
    coverage: Vec<Interval>,
    snvs: Option<Vec<Snv>>,
}

impl ReadReport {
    /// `ReadReport.NULL_REPORT`: no coverage and an empty, not null, variant list.
    fn null() -> Self {
        ReadReport {
            coverage: Vec::new(),
            snvs: Some(Vec::new()),
        }
    }

    fn first_ref_index(&self) -> i32 {
        self.coverage[0].start
    }

    fn last_ref_index(&self) -> i32 {
        self.coverage[self.coverage.len() - 1].end
    }

    /// `getVariations()`, with the null list read as empty: `combineVariations` asks only for
    /// emptiness and iterates, and both see an empty list the same way.
    fn variations(&self) -> &[Snv] {
        self.snvs.as_deref().unwrap_or(&[])
    }
}

/// `CigarUtils.countClippedBases(cigar, tail)`: soft and hard clips at one end.
///
/// The Java counts soft clips and hard clips in two passes from the same end, each stopping at
/// the first element that is not a clip; one pass summing both is the same total and raises the
/// same `IllegalArgumentException`s.
fn count_clipped_bases(cigar: &Cigar, left: bool) -> Result<i32, AsmError> {
    let is_clip = |op: Op| matches!(op, Op::S | Op::H);
    let size = cigar.elements.len();
    if size < 2 {
        if size == 1 && !is_clip(cigar.elements[0].op) {
            return Ok(0);
        }
        return Err(AsmError::Internal {
            class: "java.lang.IllegalArgumentException",
            message: "cigar is empty or completely clipped.".to_string(),
        });
    }
    let mut total = 0;
    for step in 0..size {
        let element = cigar.elements[if left { step } else { size - step - 1 }];
        if !is_clip(element.op) {
            return Ok(total);
        }
        total += element.length as i32;
    }
    Err(AsmError::Internal {
        class: "java.lang.IllegalArgumentException",
        message: format!("Input cigar {} is completely clipped.", cigar.to_text()),
    })
}

/// `GATKRead.getAttributeAsString`: the tag's value as text, or `None` (the Java's `null`) when
/// it is absent. Only the `SA` tag is asked for, and it is a string.
fn string_tag(record: &BamRecord, name: &[u8; 2]) -> Option<String> {
    match record.tags.get(Tag::new(name))? {
        TagValue::Str(text) => Some(text.clone()),
        TagValue::Char(c) => Some((*c as char).to_string()),
        TagValue::Int(value) => Some(value.to_string()),
        TagValue::Float(value) => Some(value.to_string()),
        _ => None,
    }
}

/// `ReadReport.replaceCigar`: `first` minus its last element (the clip), a `D` of
/// `deletion_length`, then `second` minus its first element, with the match after it shortened
/// by `overlap_length` when the two alignments overlap on the read. `&first[..n]` is a borrowed
/// slice of the first `n` elements, Java's `subList(0, n)`.
fn replace_cigar(
    first: &[CigarElement],
    overlap_length: i32,
    deletion_length: i32,
    second: &[CigarElement],
) -> Cigar {
    let mut elements: Vec<CigarElement> = first[..first.len() - 1].to_vec();
    elements.push(CigarElement {
        length: deletion_length as u32,
        op: Op::D,
    });
    if overlap_length == 0 {
        elements.extend_from_slice(&second[1..]);
    } else {
        let first_match = second[1];
        elements.push(CigarElement {
            length: (first_match.length as i32 - overlap_length) as u32,
            op: Op::M,
        });
        elements.extend_from_slice(&second[2..]);
    }
    Cigar::new(elements)
}

/// The state of one run: the arguments, the amplicon, the counts, and the rejected reads.
///
/// `variation_counts` is the reference's `HopscotchMap<SNVCollectionCount, ...>`: one entry per
/// distinct list of variants, with its count and the sum of its coverage (`int`, so it wraps
/// as the Java's does). The hash map's own order is never observable, because
/// `getVariationEntries` sorts before writing, so a list in first-seen order with an index by
/// key is enough. `rejected` is `None` unless `--write-rejected-reads` was given.
pub struct Run {
    settings: Settings,
    reference: Reference,
    codon_tracker: CodonTracker,
    variation_counts: Vec<(Vec<Snv>, i64, i32)>,
    variation_index: HashMap<Vec<(i32, i8, i8)>, usize>,
    total_base_calls: i64,
    read_counts: ReportTypeCounts,
    unpaired_counts: ReportTypeCounts,
    disjoint_pair_counts: ReportTypeCounts,
    overlapping_pair_counts: ReportTypeCounts,
    rejected: Option<Vec<BamRecord>>,
    /// The copies written while the current read (or pair) is processed, each with the address of
    /// the record it was taken from; see [`Run::settle_step`].
    step_writes: Vec<(usize, usize)>,
}

/// Which census a processed report is counted in.
/// Which of `unpairedCounts` and `disjointPairCounts` `processReport` is handed; the
/// overlapping census is only ever bumped by `updateCountsForPair` directly.
#[derive(Clone, Copy)]
enum Census {
    Unpaired,
    Disjoint,
}

impl Run {
    /// `onTraversalStart`, from the refusal of a coordinate-sorted input on: `coordinate_sorted` is
    /// the reads header's sort order, and `reference` the amplicon [`amplicon`] accepted.
    pub fn start(
        settings: Settings,
        coordinate_sorted: bool,
        contigs: usize,
        reference_bases: &[u8],
    ) -> Result<Run, AsmError> {
        if settings.paired_mode && coordinate_sorted {
            return Err(AsmError::user(
                "In paired mode the BAM cannot be coordinate sorted.  Mates must be adjacent.",
            ));
        }
        if settings.codon_translation.len() != N_REGULAR_CODONS {
            return Err(AsmError::user(
                "codon-translation string must contain exactly 64 characters",
            ));
        }
        let sequence = amplicon(contigs, reference_bases)?;
        let codon_tracker = CodonTracker::new(&settings.orf, &sequence)?;
        let rejected = settings.write_rejected_reads.then(Vec::new);
        Ok(Run {
            settings,
            reference: Reference::new(sequence),
            codon_tracker,
            variation_counts: Vec::new(),
            variation_index: HashMap::new(),
            total_base_calls: 0,
            read_counts: ReportTypeCounts::default(),
            unpaired_counts: ReportTypeCounts::default(),
            disjoint_pair_counts: ReportTypeCounts::default(),
            overlapping_pair_counts: ReportTypeCounts::default(),
            rejected,
            step_writes: Vec::new(),
        })
    }

    /// `read.setAttribute("XX", reportType.attributeValue)` and `rejectedReadsBAMWriter.addRead`,
    /// when there is a writer. A `null` value removes the tag, as `setAttribute` does.
    ///
    /// The address of the record (`read as *const BamRecord as usize`, the pointer as a number) is
    /// kept with the copy so that `Run::settle_step` can tell two writes of the same Java object
    /// apart from two different reads. It is only compared, never followed.
    fn write_rejected(&mut self, read: &mut BamRecord, kind: ReportType) {
        if let Some(rejected) = &mut self.rejected {
            match kind.attribute_value() {
                Some(value) => read
                    .tags
                    .insert(Tag::new(b"XX"), TagValue::Str(value.to_string())),
                None => read.tags.remove(Tag::new(b"XX")),
            }
            self.step_writes
                .push((rejected.len(), read as *const BamRecord as usize));
            rejected.push(read.clone());
        }
    }

    /// The copies of one read written while one read (or pair) was processed, given the tag the
    /// read ended that step with.
    ///
    /// The writer is htsjdk's asynchronous one (`GATKConfig` sets `use_async_io_write_samtools`),
    /// whose queue holds the record OBJECT and encodes it on its own thread. A read rejected by
    /// `getReadReport` and then written again by `processReport` is the same object queued twice,
    /// retagged in between, and the writer thread reaches the first copy only after the tool's
    /// thread has moved on: both records carry the second tag (`lowQ` twice, where a synchronous
    /// writer would write `unmapped` then `lowQ`).
    fn settle_step(&mut self) {
        let writes = std::mem::take(&mut self.step_writes);
        let Some(rejected) = &mut self.rejected else {
            return;
        };
        for (index, (position, address)) in writes.iter().enumerate() {
            if let Some((last, _)) = writes[index + 1..]
                .iter()
                .rev()
                .find(|(_, later)| later == address)
            {
                let tags = rejected[*last].tags.clone();
                rejected[*position].tags = tags;
            }
        }
    }

    /// `rejectRead`: counted in the read census, written when asked, and an empty report back.
    fn reject_read(&mut self, read: &mut BamRecord, kind: ReportType) -> ReadReport {
        self.read_counts.bump(kind);
        self.write_rejected(read, kind);
        ReadReport::null()
    }

    /// `calculateQualityTrim`.
    fn quality_trim(&self, qualities: &[u8]) -> Result<Interval, AsmError> {
        let min_q = self.settings.min_q;
        let min_length = self.settings.min_length;
        let mut start = 0usize;
        let mut high = 0;
        while start < qualities.len() {
            // `quals[readStart] < minQ` compares a signed `byte`: `as i8` reinterprets the stored
            // byte, then `as i32` widens it to `int`, so a quality above 127 is low.
            if (qualities[start] as i8 as i32) < min_q {
                high = 0;
            } else {
                high += 1;
                if high == min_length {
                    break;
                }
            }
            start += 1;
        }
        if start == qualities.len() {
            return Ok(Interval::NULL);
        }
        let start = start as i32 - (min_length - 1);
        let mut end = qualities.len() as i32 - 1;
        high = 0;
        while end >= 0 {
            if (qualities[end as usize] as i8 as i32) < min_q {
                high = 0;
            } else {
                high += 1;
                if high == min_length {
                    break;
                }
            }
            end -= 1;
        }
        Interval::new(start, end + min_length)
    }

    /// `calculateShortFragmentTrim`.
    fn fragment_trim(&self, read: &BamRecord, trim: Interval) -> Result<Interval, AsmError> {
        let min_length = self.settings.min_length;
        if trim.size() < min_length {
            return Ok(trim);
        }
        if read::is_proper_pair(read) {
            // `Math.abs(read.getFragmentLength())`: `wrapping_abs` keeps Java's answer for `int`'s
            // minimum, which is the minimum itself, still negative.
            let fragment = read::fragment_length(read).wrapping_abs();
            if read::is_reverse_strand(read) {
                let minimum_start = read::length(read) as i32 - fragment;
                if trim.start < minimum_start {
                    if trim.end - minimum_start < min_length {
                        return Ok(Interval::NULL);
                    }
                    return Interval::new(minimum_start, trim.end);
                }
            } else if fragment < trim.end {
                if fragment - trim.start < min_length {
                    return Ok(Interval::NULL);
                }
                return Interval::new(trim.start, fragment);
            }
        }
        Ok(trim)
    }

    /// `findLargeDeletions`: a supplementary alignment that continues the primary one across a gap
    /// on the reference, turned into one alignment with a deletion in the middle.
    ///
    /// The READ is changed, not a copy of its cigar, so a rejected read is written with the
    /// alignment this made.
    ///
    /// Each `continue` is one of the Java's: a malformed `SA` entry, a mapq that is low or not an
    /// integer, the other strand, a cigar that does not parse or has fewer than two elements, read
    /// intervals that gap or overlap by more than two bases, or a start that is not an integer.
    fn find_large_deletions(&self, read: &mut BamRecord) -> Result<(), AsmError> {
        let cigar = read.cigar.clone();
        if cigar.elements.len() < 2 {
            return Ok(());
        }
        let initial_clip = count_clipped_bases(&cigar, true)?;
        let final_clip = count_clipped_bases(&cigar, false)?;
        let min_alt_length = self.settings.min_alt_length;
        if initial_clip < min_alt_length && final_clip < min_alt_length {
            return Ok(());
        }
        let Some(sa) = string_tag(read, b"SA") else {
            return Ok(());
        };
        let read_length = read::length(read) as i32;
        let primary = Interval::new(initial_clip, read_length - final_clip)?;
        for alt in java_split(&sa, ';') {
            let fields = java_split(alt, ',');
            if fields.len() != 6 {
                continue;
            }
            match parse_int(fields[4]) {
                Some(mapq) if mapq < self.settings.min_mapq => continue,
                Some(_) => {}
                None => continue,
            }
            let strand = if read::is_reverse_strand(read) {
                "-"
            } else {
                "+"
            };
            if strand != fields[2] {
                continue;
            }
            let Ok(alt_cigar) = htsjdk_bam::text_parse::parse_cigar(fields[3]) else {
                continue;
            };
            if alt_cigar.elements.len() < 2 {
                continue;
            }
            let initial_alt_clip = count_clipped_bases(&alt_cigar, true)?;
            let final_alt_clip = count_clipped_bases(&alt_cigar, false)?;
            let alt_interval = Interval::new(initial_alt_clip, read_length - final_alt_clip)?;
            let overlap = primary.overlap_length(&alt_interval);
            if overlap.abs() > 2 {
                continue;
            }
            let Some(alt_start) = parse_int(fields[1]).map(|start| start - 1) else {
                continue;
            };
            if initial_clip < initial_alt_clip {
                let deletion = alt_start - read.alignment_end() + overlap;
                if deletion > 2 {
                    read.cigar =
                        replace_cigar(&cigar.elements, overlap, deletion, &alt_cigar.elements);
                    break;
                }
            } else {
                let deletion = read.alignment_start
                    - (alt_start + alt_cigar.reference_length() as i32 + 1)
                    + overlap;
                if deletion > 2 {
                    read.cigar =
                        replace_cigar(&alt_cigar.elements, overlap, deletion, &cigar.elements);
                    read.alignment_start = alt_start + 1;
                    break;
                }
            }
        }
        Ok(())
    }

    /// The `ReadReport(GATKRead, Interval, byte[])` constructor.
    ///
    /// One walk along the cigar from the alignment start, stepping the read for every operator but
    /// `D` and the reference for every one but `I`, recording a difference at each base inside the
    /// trim and on the reference. A leading soft clip is treated as a match placed before the
    /// alignment start (`refIndex -= cigarElementCount; // this can go negative and that's ok.`),
    /// so a read starting just before the amplicon still reports its first bases. Any other
    /// operator (`N`, `H`, `P`, `=`, `X`) reached inside the trim is a `GATKException`.
    fn read_report(&self, read: &mut BamRecord, trim: Interval) -> Result<ReadReport, AsmError> {
        if self.settings.find_large_deletions {
            self.find_large_deletions(read)?;
        }
        let reference = &self.reference.sequence;
        let reference_length = reference.len() as i32;
        let mut snvs = Vec::new();
        let mut coverage = Vec::new();
        let elements = &read.cigar.elements;
        let mut element_index = 0usize;
        let first = *elements.first().ok_or_else(|| AsmError::Internal {
            class: "java.util.NoSuchElementException",
            message: String::new(),
        })?;
        let mut op = first.op;
        let mut remaining = first.length as i32;
        let bases = &read.read_bases;
        let qualities = &read.base_qualities;
        let mut ref_index = read.alignment_start - 1;
        let mut read_index = 0;
        if op == Op::S {
            ref_index -= remaining;
        }
        let mut coverage_begin = -1;
        let mut coverage_end = -1;
        loop {
            if read_index >= trim.start && ref_index >= 0 {
                if coverage_begin == -1 {
                    coverage_begin = ref_index;
                    coverage_end = ref_index;
                }
                match op {
                    Op::D => snvs.push(Snv {
                        ref_index,
                        ref_call: at(reference, ref_index as i64)?,
                        variant_call: NO_CALL,
                        quality: at(qualities, read_index as i64)?,
                    }),
                    Op::I => snvs.push(Snv {
                        ref_index,
                        ref_call: NO_CALL,
                        variant_call: at(bases, read_index as i64)? & UPPERCASE_MASK,
                        quality: at(qualities, read_index as i64)?,
                    }),
                    Op::M | Op::S => {
                        let call = at(bases, read_index as i64)? & UPPERCASE_MASK;
                        let reference_call = at(reference, ref_index as i64)?;
                        if call != reference_call {
                            snvs.push(Snv {
                                ref_index,
                                ref_call: reference_call,
                                variant_call: call,
                                quality: at(qualities, read_index as i64)?,
                            });
                        }
                        if ref_index == coverage_end {
                            coverage_end += 1;
                        } else {
                            coverage.push(Interval::new(coverage_begin, coverage_end)?);
                            coverage_begin = ref_index;
                            coverage_end = ref_index + 1;
                        }
                    }
                    other => {
                        return Err(AsmError::gatk(format!(
                            "unanticipated cigar operator: {}",
                            other.to_char() as char
                        )));
                    }
                }
            }
            if op != Op::D {
                read_index += 1;
                if read_index == trim.end {
                    break;
                }
            }
            if op != Op::I {
                ref_index += 1;
                if ref_index == reference_length {
                    break;
                }
            }
            remaining -= 1;
            if remaining == 0 {
                element_index += 1;
                let element = *elements
                    .get(element_index)
                    .ok_or_else(|| AsmError::gatk("unexpectedly exhausted cigar iterator"))?;
                op = element.op;
                remaining = element.length as i32;
            }
        }
        if coverage_begin < coverage_end {
            coverage.push(Interval::new(coverage_begin, coverage_end)?);
        }
        Ok(ReadReport {
            coverage,
            snvs: Some(snvs),
        })
    }

    /// `getReadReport`.
    ///
    /// Every read's length counts toward `Total base calls`, rejected or not, and a read with
    /// a long enough trim that still covers no reference base is LOW_QUALITY too.
    fn get_read_report(&mut self, read: &mut BamRecord) -> Result<ReadReport, AsmError> {
        self.total_base_calls += read::length(read) as i64;
        if read::is_unmapped(read)
            || read::is_duplicate(read)
            || read::fails_vendor_quality_check(read)
            || (read::mapping_quality(read) as i32) < self.settings.min_mapq
        {
            return Ok(self.reject_read(read, ReportType::Unmapped));
        }
        let quality = self.quality_trim(&read.base_qualities)?;
        let trim = self.fragment_trim(read, quality)?;
        if trim.size() < self.settings.min_length {
            return Ok(self.reject_read(read, ReportType::LowQuality));
        }
        let report = self.read_report(read, trim)?;
        if report.coverage.is_empty() {
            return Ok(self.reject_read(read, ReportType::LowQuality));
        }
        self.read_counts.bump(ReportType::Evaluable);
        Ok(report)
    }

    /// `hasCleanFlanks`: `minFlankingLength` covered bases before the first variant and after the
    /// last, clamped to the amplicon's ends, so a variant near an end needs only what exists.
    fn has_clean_flanks(&self, report: &ReadReport, snvs: &[Snv]) -> bool {
        let flank = self.settings.min_flanking_length;
        let left = snvs.is_empty() || 0.max(snvs[0].ref_index - flank) >= report.first_ref_index();
        let right = snvs.is_empty()
            || (self.reference.sequence.len() as i32 - 1)
                .min(snvs[snvs.len() - 1].ref_index + flank)
                < report.last_ref_index();
        left && right
    }

    /// `ReadReport.updateCounts`.
    ///
    /// The tests run in the Java's order, and only a report that passes all of them changes the
    /// coverage and the codon counts. `snv.quality as i8` is the Java's signed `byte` against
    /// `minQ`, so a quality above 127 fails as negative.
    fn update_counts(&mut self, report: &ReadReport) -> Result<ReportType, AsmError> {
        if report.coverage.is_empty() {
            return Ok(ReportType::LowQuality);
        }
        let Some(snvs) = &report.snvs else {
            return Ok(ReportType::Inconsistent);
        };
        if snvs.iter().any(|snv| {
            (snv.quality as i8 as i32) < self.settings.min_q
                || !matches!(snv.variant_call, b'-' | b'A' | b'C' | b'G' | b'T')
        }) {
            return Ok(ReportType::LowQVar);
        }
        if !self.has_clean_flanks(report, snvs) {
            return Ok(ReportType::NoFlank);
        }
        let coverage = self.reference.update_coverage(&report.coverage)?;
        let total = Interval::new(report.first_ref_index(), report.last_ref_index())?;
        self.reference.interval_counter.add_count(total)?;
        if snvs.is_empty() {
            self.codon_tracker.report_wild_codon_counts(total)?;
            return Ok(ReportType::WildType);
        }
        let variations = self.codon_tracker.encode_snvs_as_codons(snvs)?;
        self.codon_tracker
            .report_variant_codon_counts(total, &variations)?;
        let key: Vec<(i32, i8, i8)> = snvs.iter().map(Snv::key).collect();
        match self.variation_index.get(&key) {
            Some(&index) => {
                let entry = &mut self.variation_counts[index];
                entry.1 += 1;
                // `bumpCount`: one more observation, and `totalRefCoverage += refCoverage` on an
                // `int`.
                entry.2 = entry.2.wrapping_add(coverage);
            }
            None => {
                self.variation_index
                    .insert(key, self.variation_counts.len());
                self.variation_counts.push((snvs.clone(), 1, coverage));
            }
        }
        Ok(ReportType::CalledVariant)
    }

    /// The counts a `Census` names, borrowed writably so the caller can bump them.
    fn census(&mut self, census: Census) -> &mut ReportTypeCounts {
        match census {
            Census::Unpaired => &mut self.unpaired_counts,
            Census::Disjoint => &mut self.disjoint_pair_counts,
        }
    }

    /// `processReport`: counted in the given census, and written to the rejected reads when
    /// the type is a rejection.
    fn process_report(
        &mut self,
        read: &mut BamRecord,
        report: &ReadReport,
        census: Census,
    ) -> Result<(), AsmError> {
        let kind = self.update_counts(report)?;
        self.census(census).bump(kind);
        if kind.attribute_value().is_some() {
            self.write_rejected(read, kind);
        }
        Ok(())
    }

    /// `updateCountsForPair`.
    ///
    /// A mate with no coverage leaves the other to be counted alone, in the disjoint census.
    /// Mates whose spans overlap or are adjacent (`overlapStart <= overlapEnd`, the end being
    /// exclusive) are combined and counted once in the overlapping census; disjoint mates are
    /// combined only under `--dont-ignore-disjoint-pairs`, and then nothing is written to the
    /// rejected reads whatever the outcome, because that branch of the Java writes nothing.
    fn update_counts_for_pair(
        &mut self,
        read1: &mut BamRecord,
        report1: &ReadReport,
        read2: &mut BamRecord,
        report2: &ReadReport,
    ) -> Result<(), AsmError> {
        if report1.coverage.is_empty() {
            if !report2.coverage.is_empty() {
                self.process_report(read2, report2, Census::Disjoint)?;
            }
            return Ok(());
        }
        if report2.coverage.is_empty() {
            return self.process_report(read1, report1, Census::Disjoint);
        }
        let overlap_start = report1.first_ref_index().max(report2.first_ref_index());
        let overlap_end = report1.last_ref_index().min(report2.last_ref_index());
        if overlap_start <= overlap_end {
            let combined = combine(report1, report2)?;
            let kind = self.update_counts(&combined)?;
            self.overlapping_pair_counts.bump(kind);
            if kind.attribute_value().is_some() {
                self.write_rejected(read1, kind);
                self.write_rejected(read2, kind);
            }
        } else if self.settings.dont_ignore_disjoint_pairs {
            let combined = combine(report1, report2)?;
            let kind = self.update_counts(&combined)?;
            self.disjoint_pair_counts.bump(kind);
        } else {
            if read::is_first_of_pair(read1) {
                self.process_report(read1, report1, Census::Disjoint)?;
                self.write_rejected(read2, ReportType::IgnoredMate);
            } else {
                self.process_report(read2, report2, Census::Disjoint)?;
                self.write_rejected(read1, ReportType::IgnoredMate);
            }
            self.disjoint_pair_counts.bump(ReportType::IgnoredMate);
        }
        Ok(())
    }

    /// `processReport(read, getReadReport(read), unpairedCounts)`, which is also `processRead1`
    /// apart from its `Read ... has no mate.` warning on the log.
    fn process_alone(&mut self, read: &mut BamRecord) -> Result<(), AsmError> {
        let report = self.get_read_report(read)?;
        self.process_report(read, &report, Census::Unpaired)
    }

    /// `traverse`: the primary lines, in file order, paired by adjacency in paired mode.
    ///
    /// A failure anywhere in the pass is wrapped with the count of reads the census holds at that
    /// moment and the name of the read the tool was holding as `read1`.
    pub fn traverse(&mut self, reads: Vec<BamRecord>) -> Result<(), AsmError> {
        let mut held_name = String::new();
        let result = self.traverse_inner(reads, &mut held_name);
        // The Java passes the cause to the `GATKException`; what reaches the user is the message,
        // so the cause is dropped here. `map_err` rewrites the error and passes a success through.
        result.map_err(|_cause| {
            AsmError::gatk(format!(
                "Caught unexpected exception on read {}: {}",
                self.read_counts.total(),
                held_name
            ))
        })
    }

    /// The body of `traverse`'s `try`. `held_name` is written before every step so that the
    /// wrapper sees the name `read1` holds when the failure happens; `filter` is the
    /// `PRIMARY_LINE` read filter.
    fn traverse_inner(
        &mut self,
        reads: Vec<BamRecord>,
        held_name: &mut String,
    ) -> Result<(), AsmError> {
        let primary = reads.into_iter().filter(|read| {
            !read::is_secondary_alignment(read) && !read::is_supplementary_alignment(read)
        });
        if !self.settings.paired_mode {
            for mut read in primary {
                *held_name = read.read_name.clone();
                self.process_alone(&mut read)?;
                self.settle_step();
            }
            return Ok(());
        }
        let mut read1: Option<BamRecord> = None;
        for mut read in primary {
            if !read::is_paired(&read) {
                if let Some(mut held) = read1.take() {
                    *held_name = held.read_name.clone();
                    self.process_alone(&mut held)?;
                }
                *held_name = read.read_name.clone();
                self.process_alone(&mut read)?;
            } else {
                match read1.take() {
                    None => {
                        *held_name = read.read_name.clone();
                        read1 = Some(read);
                    }
                    Some(mut held) if held.read_name != read.read_name => {
                        *held_name = held.read_name.clone();
                        self.process_alone(&mut held)?;
                        *held_name = read.read_name.clone();
                        read1 = Some(read);
                    }
                    Some(mut held) => {
                        *held_name = held.read_name.clone();
                        let report1 = self.get_read_report(&mut held)?;
                        let report2 = self.get_read_report(&mut read)?;
                        self.update_counts_for_pair(&mut held, &report1, &mut read, &report2)?;
                    }
                }
            }
            self.settle_step();
        }
        if let Some(mut held) = read1 {
            *held_name = held.read_name.clone();
            self.process_alone(&mut held)?;
            self.settle_step();
        }
        Ok(())
    }

    /// The rejected reads, in the order they were written, when `--write-rejected-reads` asked.
    pub fn rejected_reads(&self) -> Option<&[BamRecord]> {
        self.rejected.as_deref()
    }

    /// `onTraversalSuccess`'s eight reports, as `(suffix, text)` in the order they are written.
    pub fn reports(&self) -> Result<Vec<(&'static str, String)>, AsmError> {
        Ok(vec![
            (".variantCounts", self.variation_counts_report()?),
            (".refCoverage", self.ref_coverage_report()),
            (".codonCounts", self.codon_counts_report()),
            (".codonFractions", self.codon_fractions_report()),
            (".aaCounts", self.aa_counts_report()),
            (".aaFractions", self.aa_fractions_report()),
            (".readCounts", self.read_counts_report()),
            (".coverageLengthCounts", self.coverage_length_report()),
        ])
    }

    /// `codonTranslation.charAt(value)`, as a one-character string.
    fn translation(&self, value: i32) -> String {
        utf16(&[self.settings.codon_translation[value as usize]])
    }

    /// `getVariationEntries` and `writeVariationCounts`: the variant lists seen at least
    /// `--min-variant-obs` times, sorted, one row each.
    ///
    /// The sort is `SNVCollectionCount.compareTo`: variant by variant on reference index, then
    /// reference call, then variant call (the calls as signed `byte`s), with a list that is a
    /// prefix of another first. A `Vec` of tuples compares exactly that way. The mean coverage is
    /// `float` arithmetic (`1.f * totalRefCoverage / count`), widened to `double` only for
    /// `DecimalFormat("0.0")`, which is why it is computed in `f32`.
    fn variation_counts_report(&self) -> Result<String, AsmError> {
        let mut entries: Vec<&(Vec<Snv>, i64, i32)> = self
            .variation_counts
            .iter()
            .filter(|entry| entry.1 >= self.settings.min_variant_observations)
            .collect();
        entries.sort_by(|a, b| {
            let left: Vec<(i32, i8, i8)> = a.0.iter().map(Snv::key).collect();
            let right: Vec<(i32, i8, i8)> = b.0.iter().map(Snv::key).collect();
            left.cmp(&right)
        });
        let reference_length = self.reference.sequence.len() as i32;
        let flank = self.settings.min_flanking_length;
        let mut out = String::new();
        for (snvs, count, total_coverage) in entries {
            out.push_str(&count.to_string());
            out.push('\t');
            let start = 0.max(snvs[0].ref_index - flank);
            let end = reference_length.min(snvs[snvs.len() - 1].ref_index + flank);
            out.push_str(
                &self
                    .reference
                    .interval_counter
                    .count_spanners(start, end)?
                    .to_string(),
            );
            out.push('\t');
            let mean = (*total_coverage as f32) / (*count as f32);
            out.push_str(&decimal_format(f64::from(mean), 1));
            out.push('\t');
            out.push_str(&snvs.len().to_string());
            let mut separator = "\t";
            for snv in snvs {
                out.push_str(separator);
                separator = ", ";
                out.push_str(&snv.text());
            }
            self.describe_variants_as_codons(&mut out, snvs)?;
            out.push('\n');
        }
        Ok(out)
    }

    /// `isSynonymous`: a modification whose amino acid is the wild type's.
    fn is_synonymous(&self, variation: &CodonVariation) -> bool {
        let translation = &self.settings.codon_translation;
        variation.kind == CodonVariationType::Modification
            && translation[variation.codon_value as usize]
                == translation
                    [self.codon_tracker.ref_codon_values[variation.codon_id as usize] as usize]
    }

    /// `describeVariantsAsCodons`: the row's codon columns. The count, then each variation as DNA,
    /// then as an amino acid (`S` synonymous, `N` nonsense, `M` missense), then the HGVS groups
    /// separated by `;`. A group still open at the end is written only if it holds something.
    fn describe_variants_as_codons(&self, out: &mut String, snvs: &[Snv]) -> Result<(), AsmError> {
        let variations = self.codon_tracker.encode_snvs_as_codons(snvs)?;
        if variations.is_empty() {
            out.push_str("\t0");
            return Ok(());
        }
        out.push('\t');
        out.push_str(&variations.len().to_string());
        let ref_values = &self.codon_tracker.ref_codon_values;

        let mut separator = "\t";
        for variation in &variations {
            out.push_str(separator);
            separator = ", ";
            out.push_str(&(variation.codon_id + 1).to_string());
            out.push(':');
            if variation.kind == CodonVariationType::Frameshift {
                out.push_str("FS");
            } else {
                out.push_str(if variation.kind == CodonVariationType::Insertion {
                    "---"
                } else {
                    LABEL_FOR_CODON_VALUE[ref_values[variation.codon_id as usize] as usize]
                });
                out.push('>');
                out.push_str(if variation.kind == CodonVariationType::Deletion {
                    "---"
                } else {
                    LABEL_FOR_CODON_VALUE[variation.codon_value as usize]
                });
            }
        }

        separator = "\t";
        for variation in &variations {
            out.push_str(separator);
            separator = ", ";
            match variation.kind {
                CodonVariationType::Frameshift => out.push_str("FS"),
                CodonVariationType::Insertion => {
                    out.push_str("I:->");
                    out.push_str(&self.translation(variation.codon_value));
                }
                CodonVariationType::Deletion => {
                    out.push_str("D:");
                    out.push_str(&self.translation(ref_values[variation.codon_id as usize]));
                    out.push_str(">-");
                }
                CodonVariationType::Modification => {
                    let from = self.settings.codon_translation
                        [ref_values[variation.codon_id as usize] as usize];
                    let to = self.settings.codon_translation[variation.codon_value as usize];
                    let label = if from == to {
                        'S'
                    } else if is_stop(variation.codon_value) {
                        'N'
                    } else {
                        'M'
                    };
                    out.push(label);
                    out.push(':');
                    out.push_str(&utf16(&[from]));
                    out.push('>');
                    out.push_str(&utf16(&[to]));
                }
            }
        }

        separator = "\t";
        let translation = &self.settings.codon_translation;
        let mut group: Option<CodonVariationGroup> = None;
        for variation in &variations {
            match &mut group {
                None => {
                    if !self.is_synonymous(variation) {
                        group = Some(CodonVariationGroup::new(
                            ref_values,
                            translation,
                            *variation,
                        ));
                    }
                }
                Some(current) => {
                    if !current.add(*variation) {
                        out.push_str(separator);
                        separator = ";";
                        out.push_str(&current.hgvs());
                        group = None;
                        if !self.is_synonymous(variation) {
                            group = Some(CodonVariationGroup::new(
                                ref_values,
                                translation,
                                *variation,
                            ));
                        }
                    }
                }
            }
        }
        if let Some(mut current) = group {
            if !current.is_empty() {
                out.push_str(separator);
                out.push_str(&current.hgvs());
            }
        }
        Ok(())
    }

    /// `writeRefCoverage`.
    fn ref_coverage_report(&self) -> String {
        let mut out = String::from("RefPos\tCoverage\n");
        for (index, coverage) in self.reference.coverage.iter().enumerate() {
            out.push_str(&format!("{}\t{coverage}\n", index + 1));
        }
        out
    }

    /// `writeCodonCounts`: each data row ends with a tab before its total, as the header does.
    fn codon_counts_report(&self) -> String {
        let mut out = String::new();
        for label in LABEL_FOR_CODON_VALUE {
            out.push_str(label);
            out.push('\t');
        }
        out.push_str("NFS\tFS\tTotal\n");
        for row in &self.codon_tracker.codon_counts {
            let mut total = 0;
            for count in row {
                out.push_str(&count.to_string());
                out.push('\t');
                total += count;
            }
            out.push_str(&total.to_string());
            out.push('\n');
        }
        out
    }

    /// `writeCodonFractions`: `String.format("%6.2f", 100. * count / total)` per cell, so a codon
    /// no molecule covered divides zero by zero.
    fn codon_fractions_report(&self) -> String {
        let mut out = String::from("Codon");
        for label in LABEL_FOR_CODON_VALUE {
            out.push_str("   ");
            out.push_str(label);
        }
        out.push_str("   NFS    FS    Total\n");
        for (codon_id, row) in self.codon_tracker.codon_counts.iter().enumerate() {
            out.push_str(&format!("{:>5}", codon_id + 1));
            let total: i64 = row.iter().sum();
            for count in row {
                out.push_str(&percent_field(100.0 * *count as f64 / total as f64));
            }
            out.push_str(&format!("{total:>9}\n"));
        }
        out
    }

    /// The per-codon amino-acid tallies, keyed and ordered the way a `TreeMap<Character, Long>`
    /// orders them.
    fn aa_counts(&self, row: &[i64; CODON_COUNT_ROW_SIZE]) -> BTreeMap<u16, i64> {
        let mut counts = BTreeMap::new();
        for (value, count) in row.iter().take(N_REGULAR_CODONS).enumerate() {
            *counts
                .entry(self.settings.codon_translation[value])
                .or_insert(0) += count;
        }
        counts
    }

    /// `writeAACounts`: the header comes from the first codon's keys; every row has the same keys,
    /// since they are the translation string's distinct characters.
    fn aa_counts_report(&self) -> String {
        let mut out = String::new();
        for (codon_id, row) in self.codon_tracker.codon_counts.iter().enumerate() {
            let counts = self.aa_counts(row);
            if codon_id == 0 {
                let header: Vec<String> = counts.keys().map(|key| utf16(&[*key])).collect();
                out.push_str(&header.join("\t"));
                out.push('\n');
            }
            let values: Vec<String> = counts.values().map(|count| count.to_string()).collect();
            out.push_str(&values.join("\t"));
            out.push('\n');
        }
        out
    }

    /// `writeAAFractions`: the fractions are over the whole row's total, `NFS` and `FS` included,
    /// so a row's amino-acid fractions need not add to 100.
    fn aa_fractions_report(&self) -> String {
        let mut out = String::new();
        for (codon_id, row) in self.codon_tracker.codon_counts.iter().enumerate() {
            let counts = self.aa_counts(row);
            if codon_id == 0 {
                out.push_str("Codon");
                for key in counts.keys() {
                    out.push_str("     ");
                    out.push_str(&utf16(&[*key]));
                }
                out.push_str("    Total\n");
            }
            out.push_str(&format!("{:>5}", codon_id + 1));
            let total: i64 = row.iter().sum();
            for count in counts.values() {
                out.push_str(&percent_field(100.0 * *count as f64 / total as f64));
            }
            out.push_str(&format!("{total:>9}\n"));
        }
        out
    }

    /// `writeReadCounts`: the census, with each level's percentage over a different denominator.
    /// The overlapping line counts reads, `2 * overlappingPairCounts.totalCounts()`, while the rows
    /// under it count pairs.
    fn read_counts_report(&self) -> String {
        let percent = |part: i64, whole: i64| decimal_format(100.0 * part as f64 / whole as f64, 3);
        let mut out = String::new();
        let total_reads = self.read_counts.total();
        out.push_str(&format!("Total Reads:\t{total_reads}\t100.000%\n"));
        for kind in [
            ReportType::Unmapped,
            ReportType::LowQuality,
            ReportType::Evaluable,
        ] {
            let count = self.read_counts.get(kind);
            out.push_str(&format!(
                ">{}:\t{count}\t{}%\n",
                kind.label(),
                percent(count, total_reads)
            ));
        }
        let evaluable = self.read_counts.get(ReportType::Evaluable);

        let unpaired = self.unpaired_counts.total();
        if unpaired > 0 {
            out.push_str(&format!(
                ">>Unpaired reads:\t{unpaired}\t{}%\n",
                percent(unpaired, evaluable)
            ));
            type_counts(&mut out, &self.unpaired_counts);
        }
        let disjoint = self.disjoint_pair_counts.total();
        out.push_str(&format!(
            ">>Reads in disjoint pairs evaluated separately:\t{disjoint}\t{}%\n",
            percent(disjoint, evaluable)
        ));
        type_counts(&mut out, &self.disjoint_pair_counts);
        let overlapping = 2 * self.overlapping_pair_counts.total();
        out.push_str(&format!(
            ">>Reads in overlapping pairs evaluated together:\t{overlapping}\t{}%\n",
            percent(overlapping, evaluable)
        ));
        type_counts(&mut out, &self.overlapping_pair_counts);

        let total_bases = self.total_base_calls;
        out.push_str(&format!("Total base calls:\t{total_bases}\t100.000%\n"));
        let total_coverage: i64 = self.reference.coverage.iter().sum();
        out.push_str(&format!(
            ">Base calls evaluated for variants:\t{total_coverage}\t{}%\n",
            percent(total_coverage, total_bases)
        ));
        let unevaluated = total_bases - total_coverage;
        out.push_str(&format!(
            ">Base calls unevaluated:\t{unevaluated}\t{}%\n",
            percent(unevaluated, total_bases)
        ));
        out
    }

    /// `writeCoverageSizeHistogram`: the empty tail trimmed (`peel off the empty tail -- but always
    /// print the 1st 100 entries`), and bucket zero never written.
    fn coverage_length_report(&self) -> String {
        let histogram = &self.reference.coverage_size_histogram;
        let mut length = histogram.len();
        while length > 101 {
            if histogram[length - 1] > 0 {
                break;
            }
            length -= 1;
        }
        let mut out = String::new();
        for (index, count) in histogram.iter().enumerate().take(length).skip(1) {
            out.push_str(&format!("{index}\t{count}\n"));
        }
        out
    }
}

/// `writeReportTypeCounts`: the non-zero types, in declaration order.
fn type_counts(out: &mut String, counts: &ReportTypeCounts) {
    let total = counts.total();
    for kind in REPORT_TYPES {
        let count = counts.get(kind);
        if count != 0 {
            out.push_str(&format!(
                ">>>{}:\t{count}\t{}%\n",
                kind.label(),
                decimal_format(100.0 * count as f64 / total as f64, 3)
            ));
        }
    }
}

/// `combineCoverage` and `combineVariations`: a pair's two reports as one.
fn combine(report1: &ReadReport, report2: &ReadReport) -> Result<ReadReport, AsmError> {
    Ok(ReadReport {
        coverage: combine_coverage(report1, report2)?,
        snvs: combine_variations(report1, report2),
    })
}

/// `combineCoverage`: the two coverage lists merged by start, with intervals that touch or
/// overlap joined (`curInterval.getEnd() < testInterval.getStart()` is what keeps two apart).
/// On equal starts the second report's interval is taken first, as the Java's `<` does.
fn combine_coverage(report1: &ReadReport, report2: &ReadReport) -> Result<Vec<Interval>, AsmError> {
    let first = &report1.coverage;
    let second = &report2.coverage;
    if first.is_empty() {
        return Ok(second.clone());
    }
    if second.is_empty() {
        return Ok(first.clone());
    }
    let mut combined = Vec::with_capacity(first.len() + second.len());
    let (mut i, mut j) = (0usize, 0usize);
    let mut current;
    if first[0].start < second[0].start {
        current = first[0];
        i = 1;
    } else {
        current = second[0];
        j = 1;
    }
    while i < first.len() || j < second.len() {
        let test = if i >= first.len() {
            j += 1;
            second[j - 1]
        } else if j >= second.len() || first[i].start < second[j].start {
            i += 1;
            first[i - 1]
        } else {
            j += 1;
            second[j - 1]
        };
        if current.end < test.start {
            combined.push(current);
            current = test;
        } else {
            current = Interval::new(current.start, current.end.max(test.end))?;
        }
    }
    combined.push(current);
    Ok(combined)
}

/// `combineVariations`: the two variant lists merged by reference index, or `None` (the Java's
/// `null`, read later as INCONSISTENT) when they disagree inside the overlap.
///
/// A variant inside `[overlapStart, overlapEnd)` that the other mate lacks is a disagreement,
/// and so is a different variant at the same index. Where both agree the higher-quality copy
/// is kept, the second on a tie. When either list is empty the other is returned unchecked, so
/// a mate with no variants never makes a pair inconsistent.
fn combine_variations(report1: &ReadReport, report2: &ReadReport) -> Option<Vec<Snv>> {
    if report1.variations().is_empty() {
        return report2.snvs.clone();
    }
    if report2.variations().is_empty() {
        return report1.snvs.clone();
    }
    let overlap_start = report1.first_ref_index().max(report2.first_ref_index());
    let overlap_end = report1.last_ref_index().min(report2.last_ref_index());
    let inside = |index: i32| index >= overlap_start && index < overlap_end;
    let list1 = report1.variations();
    let list2 = report2.variations();
    let (mut i, mut j) = (0usize, 0usize);
    let mut combined = Vec::new();
    while i < list1.len() || j < list2.len() {
        let next = if i >= list1.len() {
            let snv = list2[j];
            j += 1;
            if inside(snv.ref_index) {
                return None;
            }
            snv
        } else if j >= list2.len() {
            let snv = list1[i];
            i += 1;
            if inside(snv.ref_index) {
                return None;
            }
            snv
        } else {
            let (snv1, snv2) = (list1[i], list2[j]);
            if snv1.ref_index < snv2.ref_index {
                i += 1;
                if inside(snv1.ref_index) {
                    return None;
                }
                snv1
            } else if snv2.ref_index < snv1.ref_index {
                j += 1;
                if inside(snv2.ref_index) {
                    return None;
                }
                snv2
            } else if !snv1.same_as(&snv2) {
                return None;
            } else {
                i += 1;
                j += 1;
                if (snv1.quality as i8) > (snv2.quality as i8) {
                    snv1
                } else {
                    snv2
                }
            }
        };
        combined.push(next);
    }
    Some(combined)
}

/// `new DecimalFormat(pattern)` for `0.0` and `0.000`: HALF_EVEN, and exactly `places` fraction
/// digits, which is what a `0` in the pattern means where `#` would drop trailing zeros.
pub fn decimal_format(value: f64, places: usize) -> String {
    let body = gatk_annotation::decimal_format::DecimalFormat::new(places).format(value);
    if !value.is_finite() {
        return body;
    }
    let (whole, fraction) = body.split_once('.').unwrap_or((&body, ""));
    format!("{whole}.{fraction:0<places$}")
}

/// `String.format("%6.2f", value)`.
///
/// The scale-2 rounding is [`format_decimals`]'s, and `{:>6}` right-aligns it to width 6 as the
/// `6` in `%6.2f` does.
fn percent_field(value: f64) -> String {
    format!("{:>6}", format_decimals(value, 2))
}
