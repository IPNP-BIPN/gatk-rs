//! `GroundTruthReadsBuilder`: every read scored against the haplotype its two ancestral
//! references give it.
//!
//! The flow-based scoring engine is not ported. What is ported is the translation from the
//! aligned contig to the two ancestral ones, the filters that decide which reads survive it, and
//! the shape of the row each survivor becomes.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.groundtruth.SingleFileLocationTranslator`,
//! `org.broadinstitute.hellbender.tools.walkers.groundtruth.AncestralContigLocationTranslator`
//! and `org.broadinstitute.hellbender.tools.walkers.groundtruth.GroundTruthReadsBuilder`
//! in GATK 4.6.2.0.

/// The two ancestor names, which are what the translated contig and the CSV file are named for.
pub const MATERNAL: &str = "maternal";
pub const PATERNAL: &str = "paternal";

/// The fill values the tool writes into a flow key it could not read.
pub const DEFAULT_FILL_VALUE: i32 = -65;
pub const NONREF_FILL_VALUE: i32 = -80;
pub const UNKNOWN_FILL_VALUE: i32 = -85;
pub const SOFTCLIP_FILL_VALUE: i32 = -83;

/// One translation table: a position and the offset that applies from it on.
///
/// The first line of the file is IGNORED whatever it holds, so a table without a header loses its
/// first row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Translator {
    pub positions: Vec<i32>,
    pub offsets: Vec<i32>,
}

impl Translator {
    /// The file's rows, its first line dropped.
    pub fn parse(text: &str) -> Translator {
        let mut positions = Vec::new();
        let mut offsets = Vec::new();
        for line in text.lines().skip(1) {
            if line.is_empty() {
                continue;
            }
            let mut parts = line.split(',');
            let position = parts.next().and_then(|v| v.parse().ok());
            let offset = parts.next().and_then(|v| v.parse().ok());
            if let (Some(position), Some(offset)) = (position, offset) {
                positions.push(position);
                offsets.push(offset);
            }
        }
        Translator { positions, offsets }
    }

    /// The position, translated.
    ///
    /// A position between two rows takes the EARLIER row's offset, the search falling back on the
    /// insertion point less two. A position BEFORE the first row therefore indexes at minus one,
    /// which is why the file is documented as starting at position one: nothing checks it.
    pub fn translate(&self, from: i32) -> Option<i32> {
        match self.positions.binary_search(&from) {
            Ok(index) => Some(from + self.offsets[index]),
            Err(insertion) => {
                // `-index - 2` in the reference, where `index` is `-insertion - 1`.
                let earlier = insertion as i64 - 1;
                if earlier < 0 {
                    return None;
                }
                Some(from + self.offsets[earlier as usize])
            }
        }
    }
}

/// A closed interval on one of the ancestral contigs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

/// The refusal a read whose translated span collapses produces.
///
/// It is caught by the traversal and counted rather than propagated, so a read that hits it is
/// skipped and the run carries on.
pub fn translation_failure(
    contig: &str,
    start: i32,
    end: i32,
    ancestor: &str,
    translated_start: i32,
    translated_end: i32,
) -> String {
    format!(
        "location {contig}:{start}-{end} failed to translate for {ancestor}, \
         start:{translated_start} ,end:{translated_end}"
    )
}

/// One read's span on one ancestral contig.
///
/// The contig is the read's own name with the ancestor appended, so the reference file has to
/// carry `<contig>_maternal` and `<contig>_paternal` rather than the aligned name. The end must be
/// STRICTLY past the start, so a translation that collapses a read is a failure and a read of one
/// base never translates at all.
pub fn translate_span(
    translator: &Translator,
    contig: &str,
    start: i32,
    end: i32,
    ancestor: &str,
) -> Result<Interval, String> {
    let translated_start = translator
        .translate(start)
        .ok_or_else(|| translation_failure(contig, start, end, ancestor, 0, 0))?;
    let translated_end = translator
        .translate(end)
        .ok_or_else(|| translation_failure(contig, start, end, ancestor, translated_start, 0))?;
    if translated_end > translated_start {
        Ok(Interval {
            contig: format!("{contig}_{ancestor}"),
            start: translated_start,
            end: translated_end,
        })
    } else {
        Err(translation_failure(
            contig,
            start,
            end,
            ancestor,
            translated_start,
            translated_end,
        ))
    }
}

/// One cigar element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarElement {
    pub operator: char,
    pub length: i32,
}

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

/// `isEndSoftclipped`: whether the read's LAST cigar element is a soft clip.
///
/// It is the last element and not either one, so a read clipped only at its front is not
/// soft-clipped as far as this filter is concerned.
pub fn is_end_softclipped(cigar: &[CigarElement]) -> bool {
    cigar.last().is_some_and(|element| element.operator == 'S')
}

/// Whether a soft clip is poly-T, which is what spares it from the discard.
pub fn is_polyt(bases: &[u8]) -> bool {
    !bases.is_empty() && bases.iter().all(|base| *base == b'T')
}

/// The arguments that decide which reads survive.
#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub min_mapping_quality: i32,
    pub max_read_quality: Option<i32>,
    pub discard_non_polyt_softclipped_reads: bool,
    pub include_supplementary_alignments: bool,
    /// Zero means the filter is off, not a threshold of zero.
    pub min_haplotype_score: f64,
    pub min_haplotype_score_delta: f64,
    pub max_output_reads: Option<usize>,
    pub prepend_sequence: String,
    pub append_sequence: String,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            min_mapping_quality: 0,
            max_read_quality: None,
            discard_non_polyt_softclipped_reads: true,
            include_supplementary_alignments: false,
            min_haplotype_score: 0.0,
            min_haplotype_score_delta: 0.0,
            max_output_reads: None,
            prepend_sequence: String::new(),
            append_sequence: String::new(),
        }
    }
}

/// Whether the two score filters keep a read.
///
/// Both are off when zero rather than being a threshold of zero, and both compare with a strict
/// GREATER-THAN against a value the reference's own comment doubts: the scores are negative, so
/// `--min-haplotype-score` keeps the reads whose worse haplotype scores at or BELOW it.
pub fn keeps_scores(maternal: f64, paternal: f64, arguments: &Arguments) -> bool {
    if arguments.min_haplotype_score != 0.0
        && maternal.min(paternal) > arguments.min_haplotype_score
    {
        return false;
    }
    if arguments.min_haplotype_score_delta != 0.0
        && (maternal - paternal).abs() > arguments.min_haplotype_score_delta
    {
        return false;
    }
    true
}

/// The columns the CSV carries, in the order the tool holds them.
pub const CSV_FIELD_ORDER: [&str; 22] = [
    "ReadName",
    "ReadChrom",
    "ReadStart",
    "ReadEnd",
    "PaternalHaplotypeScore",
    "MaternalHaplotypeScore",
    "RefHaplotypeScore",
    "ReadKey",
    "BestHaplotypeKey",
    "ConsensusHaplotypeKey",
    "tm",
    "mapq",
    "flags",
    "ReadCigar",
    "ReadSequence",
    "PaternalHaplotypeSequence",
    "MaternalHaplotypeSequence",
    "BestHaplotypeSequence",
    "ReadUnclippedStart",
    "ReadUnclippedEnd",
    "PaternalHaplotypeInterval",
    "MaternalHaplotypeInterval",
];

/// The header line, which is the column order joined by commas.
pub fn header() -> String {
    CSV_FIELD_ORDER.join(",")
}

/// One CSV field, quoted when it holds a comma.
///
/// The flow keys hold commas of their own, so a reader that splits on the comma alone reads the
/// columns out of step.
pub fn quote(value: &str) -> String {
    if value.contains(',') || value.contains('"') {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

/// One CSV row, split on commas that are not inside quotes.
pub fn split_row(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut quoted = false;
    for character in line.chars() {
        match character {
            '"' => quoted = !quoted,
            ',' if !quoted => fields.push(std::mem::take(&mut current)),
            _ => current.push(character),
        }
    }
    fields.push(current);
    fields
}
