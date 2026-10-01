//! `GroundTruthReadsBuilder`: every read scored against the haplotype its two ancestral
//! references give it.
//!
//! The translation from the aligned contig to the two ancestral ones, the filters that decide
//! which reads survive it, and the shape of the row each survivor becomes: the haplotype keys, the
//! false SNP compensation and the read's key as the CSV writes it. The score is
//! `FlowFeatureMapper.computeLikelihoodLocal`, and the walk is the `gatk-cli` runner's.
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

/// `SequenceUtil.reverseComplement` over bytes: the four bases in either case complemented, every
/// other byte kept, and the order reversed.
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|base| match *base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        })
        .collect()
}

/// `reverseComplement(bases, isReversed)`: the bases as they are for a forward read.
pub fn oriented(bases: &[u8], reversed: bool) -> Vec<u8> {
    if reversed {
        reverse_complement(bases)
    } else {
        bases.to_vec()
    }
}

/// What the reference raises instead of answering, by class and message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Thrown {
    pub class: &'static str,
    pub message: String,
}

fn out_of_bounds(index: i64, length: usize) -> Thrown {
    Thrown {
        class: "java.lang.ArrayIndexOutOfBoundsException",
        message: format!("Index {index} out of bounds for length {length}"),
    }
}

/// `new Haplotype(bases, ...)`, which refuses an empty allele.
pub fn check_haplotype(bases: &[u8]) -> Result<(), Thrown> {
    if bases.is_empty() {
        Err(Thrown {
            class: "java.lang.IllegalArgumentException",
            message: "Null alleles are not supported".to_string(),
        })
    } else {
        Ok(())
    }
}

/// `detectFalseSNP`: how many haplotype bases to skip so that the read's leading homopolymer and
/// the five bases after the skip line up, at most five, zero when no skip does.
pub fn detect_false_snp(haplotype: &[u8], read: &[u8]) -> Result<usize, Thrown> {
    const MAX_SKIP: usize = 5;
    const REMAINING: usize = 5;
    let first_h = *haplotype.first().ok_or_else(|| out_of_bounds(0, 0))?;
    let first_r = *read.first().ok_or_else(|| out_of_bounds(0, 0))?;
    if first_h == first_r {
        return Ok(0);
    }
    let homo = read.iter().take_while(|base| **base == first_r).count();
    for skip in 1..=MAX_SKIP {
        if skip + REMAINING > haplotype.len() {
            break;
        }
        if (skip + 1).max(REMAINING) > read.len() {
            break;
        }
        if skip + 1 > homo {
            break;
        }
        let mut equal = true;
        for i in 0..REMAINING {
            let r = *read
                .get(skip + i)
                .ok_or_else(|| out_of_bounds((skip + i) as i64, read.len()))?;
            if r != haplotype[skip + i] {
                equal = false;
                break;
            }
        }
        if equal {
            return Ok(skip);
        }
    }
    Ok(0)
}

/// `buildHaplotypeKey`: the sequence's key in the flow order the read was synthesised in, its
/// leading zero flows dropped, and zeros prepended for the flows between the order's first `T` and
/// the sequence's first base, unless that base is a `T` or an `N`.
pub fn haplotype_key(
    sequence: &[u8],
    flow_order: &str,
    reversed: bool,
) -> Result<Vec<i32>, Thrown> {
    let seq = oriented(sequence, reversed);
    check_haplotype(&seq)?;
    let order = if reversed {
        String::from_utf8_lossy(&reverse_complement(flow_order.as_bytes())).into_owned()
    } else {
        flow_order.to_string()
    };
    let haplotype = crate::flow_pairhmm_align_reads_to_haplotypes::FlowHaplotype::new(&seq, &order)
        .ok_or_else(|| Thrown {
            class: "org.broadinstitute.hellbender.exceptions.GATKException",
            message: format!(
                "baseArrayToKey periodGuard tripped, on {}, flowOrder: {order} This probably indicates the presence of a base (value) in the sequence that is not included in the provided flow order",
                String::from_utf8_lossy(&seq)
            ),
        })?;
    let mut key: &[i32] = &haplotype.key;
    while *key.first().ok_or_else(|| out_of_bounds(0, 0))? == 0 {
        key = &key[1..];
    }
    let mut zeros = 0usize;
    if seq[0] != b'T' && seq[0] != b'N' {
        let flows = &haplotype.flow_order;
        let mut offset = 0usize;
        while *flows
            .get(offset)
            .ok_or_else(|| out_of_bounds(offset as i64, flows.len()))?
            != b'T'
        {
            offset += 1;
        }
        while flows[offset] != seq[0] {
            zeros += 1;
            offset = (offset + 1) % flows.len();
        }
    }
    let mut out = vec![0; zeros];
    out.extend_from_slice(key);
    Ok(out)
}

/// `keyBases`: the bases a key spells, its fill values not counted.
pub fn key_bases(key: &[i32]) -> usize {
    key.iter().filter(|v| **v > 0).map(|v| *v as usize).sum()
}

/// `buildConsensusKey`: flow by flow over the shorter key, the value where the two agree and -72
/// where they do not.
pub fn consensus_key(a: &[i32], b: &[i32]) -> Vec<i32> {
    a.iter()
        .zip(b)
        .map(|(x, y)| if x == y { *x } else { -72 })
        .collect()
}

/// `flowKeyAsCsvString(key)`: the values comma-joined inside quotes.
pub fn key_csv(key: &[i32]) -> String {
    let joined: Vec<String> = key.iter().map(|v| v.to_string()).collect();
    format!("\"{}\"", joined.join(","))
}

/// `flowKeyAsCsvString(key, seq, flowOrder)`: the read's key with its leading zero flows dropped
/// and zeros written for the flows from the order's first `T` to the read's first base.
pub fn read_key_csv(key: &[i32], sequence: &[u8], flow_order: &[u8]) -> Result<String, Thrown> {
    let mut key: &[i32] = key;
    while *key.first().ok_or_else(|| out_of_bounds(0, 0))? == 0 {
        key = &key[1..];
    }
    let mut out = String::from("\"");
    let first = *sequence.first().ok_or_else(|| Thrown {
        class: "java.lang.StringIndexOutOfBoundsException",
        message: "index 0, length 0".to_string(),
    })?;
    if first != b'T' && first != b'N' {
        let mut offset = 0usize;
        while *flow_order.get(offset).ok_or_else(|| Thrown {
            class: "java.lang.StringIndexOutOfBoundsException",
            message: format!("index {offset}, length {}", flow_order.len()),
        })? != b'T'
        {
            offset += 1;
        }
        while flow_order[offset] != first {
            out.push_str("0,");
            offset = (offset + 1) % flow_order.len();
        }
    }
    let joined: Vec<String> = key.iter().map(|v| v.to_string()).collect();
    out.push_str(&joined.join(","));
    out.push('"');
    Ok(out)
}
