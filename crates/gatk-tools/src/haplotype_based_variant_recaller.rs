//! `HaplotypeBasedVariantRecaller`: every allele a haplotype carries, scored against every read
//! that haplotype spans.
//!
//! The PairHMM that produces the likelihoods is not ported. What is ported is everything around
//! it: which haplotype group a variant is scored against, how a matrix line is built and sorted,
//! and the ways a line is dropped or comes out wrong.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.variantrecalling.HaplotypeRegionWalker` and
//! `org.broadinstitute.hellbender.tools.walkers.variantrecalling.VariantRecallerResultWriter`
//! in GATK 4.6.2.0.

/// A half-open-free interval, closed at both ends, which is what every span here is.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Span {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Span {
    pub fn new(contig: &str, start: i32, end: i32) -> Span {
        Span {
            contig: contig.to_string(),
            start,
            end,
        }
    }

    /// `SimpleInterval.contains`, which needs the whole of the other interval.
    pub fn contains(&self, other: &Span) -> bool {
        self.contig == other.contig && self.start <= other.start && self.end >= other.end
    }

    pub fn overlaps(&self, other: &Span) -> bool {
        self.contig == other.contig && self.start <= other.end && self.end >= other.start
    }
}

impl std::fmt::Display for Span {
    /// `SimpleInterval.toString`, which is what the matrix's header line carries.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start, self.end)
    }
}

// ================================================================================================
// The haplotype groups.
// ================================================================================================

/// The prefix that makes a record in the haplotype BAM a haplotype.
///
/// Any other record is passed over however well it fits, so the file may hold anything else
/// alongside them.
pub const HAPLOTYPE_NAME_PREFIX: &str = "HC_";

pub fn is_haplotype_record(name: &str) -> bool {
    name.starts_with(HAPLOTYPE_NAME_PREFIX)
}

/// One record of the haplotype BAM, reduced to what the walk reads off it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HaplotypeRecord {
    pub name: String,
    pub span: Span,
}

/// The groups a query yields, in the order the reader hands them over.
///
/// A record whose span differs from the group's first CLOSES the group and opens a new one, so
/// two runs of the same span separated by a third are two groups rather than one.
pub fn groups(records: &[HaplotypeRecord]) -> Vec<Vec<HaplotypeRecord>> {
    let mut groups: Vec<Vec<HaplotypeRecord>> = Vec::new();
    let mut current: Vec<HaplotypeRecord> = Vec::new();
    for record in records
        .iter()
        .filter(|record| is_haplotype_record(&record.name))
    {
        if let Some(first) = current.first() {
            if first.span != record.span {
                groups.push(std::mem::take(&mut current));
            }
        }
        current.push(record.clone());
    }
    if !current.is_empty() {
        groups.push(current);
    }
    groups
}

/// `fitnessScore`: one less twice the distance of the variant from the group's halfway point.
///
/// Both gaps are floored at one before the ratio is taken, so a variant flush against either end
/// does not divide by zero and does not score zero either. An empty group scores zero.
pub fn fitness_score(location: &Span, group: &[HaplotypeRecord]) -> f64 {
    let Some(first) = group.first() else {
        return 0.0;
    };
    let before = std::cmp::max(1, location.start - first.span.start) as f64;
    let after = std::cmp::max(1, first.span.end - location.end) as f64;
    1.0 - 2.0 * (0.5 - before / (before + after)).abs()
}

/// `forBest`: the group with the highest fitness, ties going to the FIRST.
///
/// The comparison is strict, so a later group has to beat the one held rather than match it.
pub fn best_group<'a>(
    location: &Span,
    groups: &'a [Vec<HaplotypeRecord>],
) -> Option<&'a Vec<HaplotypeRecord>> {
    let mut best: Option<&Vec<HaplotypeRecord>> = None;
    for group in groups {
        match best {
            None => best = Some(group),
            Some(held) if fitness_score(location, group) > fitness_score(location, held) => {
                best = Some(group)
            }
            Some(_) => {}
        }
    }
    best.filter(|group| !group.is_empty())
}

// ================================================================================================
// The cigar walk.
// ================================================================================================

/// One cigar element: an operator and a length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarElement {
    pub operator: char,
    pub length: i32,
}

impl CigarElement {
    pub fn consumes_read_bases(self) -> bool {
        matches!(self.operator, 'M' | 'I' | 'S' | '=' | 'X')
    }

    pub fn consumes_reference_bases(self) -> bool {
        matches!(self.operator, 'M' | 'D' | 'N' | '=' | 'X')
    }
}

/// A cigar string, as its elements.
pub fn parse_cigar(text: &str) -> Vec<CigarElement> {
    let mut elements = Vec::new();
    let mut length = 0i32;
    for character in text.chars() {
        if let Some(digit) = character.to_digit(10) {
            length = length * 10 + digit as i32;
        } else {
            elements.push(CigarElement {
                operator: character,
                length,
            });
            length = 0;
        }
    }
    elements
}

/// `getOffsetOnRead`: the offset in the read of the base at `offset` positions into the alignment.
///
/// The walk has a hole in it. A read-consuming element returns as soon as the remaining offset
/// fits inside it, and the reference-consuming subtraction happens AFTERWARDS, so a deletion
/// drives the offset negative and the very next match element then returns
/// `read_offset + offset` with a negative offset: an index that many bases too early rather than
/// the refusal a variant inside a deletion should be. Only an offset that runs off the end of the
/// read returns nothing.
pub fn offset_on_read(cigar: &[CigarElement], initial: i32) -> Option<i32> {
    let mut read_offset = 0i32;
    let mut offset = initial;
    for element in cigar {
        if element.consumes_read_bases() {
            if offset < element.length {
                return Some(read_offset + offset);
            }
            read_offset += element.length;
        }
        if element.consumes_reference_bases() {
            offset -= element.length;
        }
    }
    None
}

// ================================================================================================
// The matrix lines.
// ================================================================================================

/// One read, reduced to what a matrix line carries.
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    pub name: String,
    pub span: Span,
    pub cigar: Vec<CigarElement>,
    pub bases: Vec<u8>,
    pub is_duplicate: bool,
    pub is_reverse: bool,
    pub mapping_quality: i32,
    /// The flow key length, which is zero for a read that is not flow-based.
    pub key_length: i32,
    pub sample: String,
    pub unclipped_start: i32,
    pub unclipped_end: i32,
}

/// `Double.toString` for the likelihood columns, which is Rust's `{:?}` for these values.
fn java_double(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    format!("{value:?}")
}

/// The variant's type, which decides only whether its end is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantKind {
    Mixed,
    Other,
}

/// The header line one variant gets: its position, the haplotype span and the alleles.
///
/// The end is omitted for a one-base variant AND for a MIXED one however long it is, which is the
/// one place the variant's type is consulted at all.
pub fn header_line(
    contig: &str,
    start: i32,
    end: i32,
    kind: VariantKind,
    haplotype_span: &Span,
    alleles: &[String],
) -> String {
    let mut line = format!("#{contig}:{start}");
    if kind != VariantKind::Mixed && end != start {
        line.push_str(&format!("-{end}"));
    }
    line.push_str(&format!(" {haplotype_span}"));
    for allele in alleles {
        line.push(' ');
        line.push_str(allele);
    }
    line
}

/// One matrix line, and the key it is sorted by.
#[derive(Debug, Clone, PartialEq)]
pub struct MatrixLine {
    pub sort_key: f64,
    pub text: String,
}

/// The bases a read carries over the variant, and the unclipped offset of the first of them.
///
/// Nothing comes back when the read does not span the whole variant, and nothing comes back when
/// the cigar walk runs off the end of the read: both are lines that are never added.
pub fn variant_bases(read: &Read, variant: &Span) -> Option<(String, i32)> {
    if !read.span.contains(variant) {
        return None;
    }
    let offset = variant.start - read.span.start;
    let length = variant.end - variant.start + 1;
    let mut bases = String::new();
    let mut first_unclipped = 0;
    for i in 0..length {
        let read_offset = offset_on_read(&read.cigar, offset + i)?;
        bases.push(*read.bases.get(read_offset as usize)? as char);
        first_unclipped = if read.is_reverse {
            (read.bases.len() as i32 - read_offset - 1) + (read.unclipped_end - read.span.end)
        } else {
            read_offset + (read.span.start - read.unclipped_start)
        };
    }
    Some((bases, first_unclipped))
}

/// One line of the matrix, or nothing.
///
/// Two things drop a line. Every likelihood being negative infinity is an alignment failure, and
/// the read not yielding bases over the variant leaves the line unfinished. The SORT KEY is the
/// LAST allele's likelihood rather than the best of them: the loop that collects the columns
/// overwrites it each time round.
pub fn matrix_line(read: &Read, variant: &Span, likelihoods: &[f64]) -> Option<MatrixLine> {
    if likelihoods.iter().all(|value| *value == f64::NEG_INFINITY) {
        return None;
    }
    let sort_key = *likelihoods.last().unwrap_or(&f64::NEG_INFINITY);
    let (bases, first_unclipped) = variant_bases(read, variant)?;
    if bases.is_empty() {
        return None;
    }
    let columns: Vec<String> = likelihoods
        .iter()
        .map(|value| java_double(*value))
        .collect();
    Some(MatrixLine {
        sort_key,
        text: format!(
            "{} {} {} {} {} {} {} {} {}",
            read.name,
            read.key_length,
            if read.is_duplicate { 1 } else { 0 },
            if read.is_reverse { 1 } else { 0 },
            read.mapping_quality,
            columns.join(" "),
            bases,
            first_unclipped,
            read.sample
        ),
    })
}

/// The lines of one variant, sorted by their key, descending.
///
/// The sort is `-Double.compare(a, b)`, which is stable, so two reads with the same last column
/// keep the order they were evidenced in.
pub fn sorted_lines(lines: &[MatrixLine]) -> Vec<String> {
    let mut lines = lines.to_vec();
    lines.sort_by(|a, b| compare_double(b.sort_key, a.sort_key));
    lines.into_iter().map(|line| line.text).collect()
}

/// `Double.compare`: numeric first, then the signed bit pattern.
fn compare_double(a: f64, b: f64) -> std::cmp::Ordering {
    if a < b {
        return std::cmp::Ordering::Less;
    }
    if a > b {
        return std::cmp::Ordering::Greater;
    }
    let bits = |value: f64| {
        if value.is_nan() {
            f64::NAN.to_bits() as i64
        } else {
            value.to_bits() as i64
        }
    };
    bits(a).cmp(&bits(b))
}

/// One variant's whole block: its header line, then its sorted matrix lines.
pub fn variant_block(header: &str, lines: &[MatrixLine]) -> String {
    let mut text = String::from(header);
    text.push('\n');
    for line in sorted_lines(lines) {
        text.push_str(&line);
        text.push('\n');
    }
    text
}
