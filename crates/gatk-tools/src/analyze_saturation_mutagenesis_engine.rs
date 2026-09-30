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
//! Ported from `org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis` (GATK 4.6.2.0).

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

const NO_CALL: u8 = b'-';
const UPPERCASE_MASK: u8 = 0xDF;
const N_REGULAR_CODONS: usize = 64;
const FRAME_PRESERVING_INDEL_INDEX: usize = 64;
const FRAME_SHIFTING_INDEL_INDEX: usize = 65;
const CODON_COUNT_ROW_SIZE: usize = 66;
const NO_FRAME_SHIFT_CODON: i32 = -1;

const LABEL_FOR_CODON_VALUE: [&str; 64] = [
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA",
    "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC",
    "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG",
    "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT",
    "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT",
];

/// The tool's own arguments, at the values the command line gave them.
#[derive(Debug, Clone)]
pub struct Settings {
    pub min_q: i32,
    pub min_length: i32,
    pub min_flanking_length: i32,
    pub min_mapq: i32,
    pub orf: String,
    pub min_variant_observations: i64,
    pub find_large_deletions: bool,
    pub min_alt_length: i32,
    /// `--codon-translation`, as the UTF-16 units `charAt` indexes.
    pub codon_translation: Vec<u16>,
    pub paired_mode: bool,
    pub dont_ignore_disjoint_pairs: bool,
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

/// Why a run stopped.
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
fn out_of_bounds(index: i64, length: usize) -> AsmError {
    AsmError::Internal {
        class: "java.lang.ArrayIndexOutOfBoundsException",
        message: format!("Index {index} out of bounds for length {length}"),
    }
}

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
    /// The `XX` value a rejected read is tagged with, or none for the types that are not rejections.
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

    fn ordinal(self) -> usize {
        REPORT_TYPES
            .iter()
            .position(|kind| *kind == self)
            .expect("every type is listed")
    }
}

/// `ReportTypeCounts`.
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CodonVariationType {
    Frameshift,
    Insertion,
    Deletion,
    Modification,
}

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

fn is_stop(codon_value: i32) -> bool {
    codon_value == 0x30 || codon_value == 0x32 || codon_value == 0x38
}

/// `"ACGT".indexOf(base)`.
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
struct CodonTracker {
    reference: Vec<u8>,
    exons: Vec<Interval>,
    codon_counts: Vec<[i64; CODON_COUNT_ROW_SIZE]>,
    ref_codon_values: Vec<i32>,
}

impl CodonTracker {
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

    /// `encodeSNVsAsCodons`.
    fn encode_snvs_as_codons(&self, snvs: &[Snv]) -> Result<Vec<CodonVariation>, AsmError> {
        let mut variations = Vec::new();
        let mut next = 0usize;
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
                    _ => {
                        let base = at(&self.reference, ref_index as i64)?;
                        codon_value = (codon_value << 2) | base_value(base);
                        codon_value_altered = true;
                        bump_ref_index = true;
                    }
                }
                if bump_ref_index {
                    ref_index += 1;
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
                if lead_lag == 0 && codon_phase == 0 {
                    break;
                }
            }
        }
        Ok(variations)
    }

    fn report_variant_codon_counts(
        &mut self,
        coverage: Interval,
        variations: &[CodonVariation],
    ) -> Result<(), AsmError> {
        let starting = (self.exonic_base_index(coverage.start) + 2) / 3;
        let ending = self.exonic_base_index(coverage.end) / 3;
        let mut cursor = 0usize;
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

    fn is_empty(&self) -> bool {
        self.sub_count + self.ins_count + self.del_count == 0
    }

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
#[derive(Debug, Clone, Default)]
struct ReadReport {
    coverage: Vec<Interval>,
    snvs: Option<Vec<Snv>>,
}

impl ReadReport {
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

    fn variations(&self) -> &[Snv] {
        self.snvs.as_deref().unwrap_or(&[])
    }
}

/// `CigarUtils.countClippedBases(cigar, tail)`: soft and hard clips at one end.
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

fn string_tag(record: &BamRecord, name: &[u8; 2]) -> Option<String> {
    match record.tags.get(Tag::new(name))? {
        TagValue::Str(text) => Some(text.clone()),
        TagValue::Char(c) => Some((*c as char).to_string()),
        TagValue::Int(value) => Some(value.to_string()),
        TagValue::Float(value) => Some(value.to_string()),
        _ => None,
    }
}

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

    fn census(&mut self, census: Census) -> &mut ReportTypeCounts {
        match census {
            Census::Unpaired => &mut self.unpaired_counts,
            Census::Disjoint => &mut self.disjoint_pair_counts,
        }
    }

    /// `processReport`.
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
        result.map_err(|_cause| {
            AsmError::gatk(format!(
                "Caught unexpected exception on read {}: {}",
                self.read_counts.total(),
                held_name
            ))
        })
    }

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

    fn translation(&self, value: i32) -> String {
        utf16(&[self.settings.codon_translation[value as usize]])
    }

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

    fn is_synonymous(&self, variation: &CodonVariation) -> bool {
        let translation = &self.settings.codon_translation;
        variation.kind == CodonVariationType::Modification
            && translation[variation.codon_value as usize]
                == translation
                    [self.codon_tracker.ref_codon_values[variation.codon_id as usize] as usize]
    }

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

    fn ref_coverage_report(&self) -> String {
        let mut out = String::from("RefPos\tCoverage\n");
        for (index, coverage) in self.reference.coverage.iter().enumerate() {
            out.push_str(&format!("{}\t{coverage}\n", index + 1));
        }
        out
    }

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
fn percent_field(value: f64) -> String {
    format!("{:>6}", format_decimals(value, 2))
}
