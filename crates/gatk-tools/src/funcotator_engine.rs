//! The GENCODE annotation `Funcotator` and `FuncotateSegments` share: a GTF read into genes, a
//! variant or a segment placed on their transcripts, and the funcotations that come out.
//!
//! The strings the funcotations are built out of (codon, protein and cDNA changes) are
//! [`crate::funcotator`]'s, measured by the `funcotator` suite. What is here is the walk around
//! them: which transcripts a variant is placed on, which region of each it falls in, the reference
//! windows each field is cut from, the order the transcripts are sorted into and the other
//! transcripts each one names.
//!
//! Substitutions (a SNP or an equal-length ONP) and symbolic alleles are ported; an insertion or a
//! deletion is refused as the port's limitation by [`EngineError::Limitation`], because the
//! reference's indel branches are not.
//!
//! Ported from `org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory`,
//! `GencodeFuncotation`, `FuncotatorUtils`, `ProteinChangeInfo`, `TranscriptSelectionMode`,
//! `SegmentExonUtils` and `org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfCodec` in
//! GATK 4.6.2.0.

use crate::funcotator as fu;
use std::collections::BTreeMap;

// ================================================================================================
// The GTF.
// ================================================================================================

/// A genomic strand, as a GTF writes it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Positive,
    Negative,
}

impl Strand {
    /// `Strand.toString()`, which is the GTF's own symbol.
    pub fn symbol(self) -> &'static str {
        match self {
            Strand::Positive => "+",
            Strand::Negative => "-",
        }
    }
}

/// A closed genomic interval on the one contig everything here shares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: i32,
    pub end: i32,
}

impl Span {
    pub fn new(start: i32, end: i32) -> Self {
        Span { start, end }
    }
    /// `Locatable.overlaps`.
    pub fn overlaps(&self, other: &Span) -> bool {
        self.start <= other.end && other.start <= self.end
    }
    /// `Locatable.contains`.
    pub fn contains(&self, other: &Span) -> bool {
        self.start <= other.start && other.end <= self.end
    }
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }
}

/// An exon, with the coding pieces the codec attached to it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Exon {
    pub span: Span,
    pub number: i32,
    pub gene_type: String,
    pub cds: Option<Span>,
    pub start_codon: Option<Span>,
    pub stop_codon: Option<Span>,
}

/// A transcript, its exons in file order and its UTRs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Transcript {
    pub contig: String,
    pub span: Span,
    pub strand: Strand,
    pub transcript_id: String,
    pub transcript_type: String,
    pub gene_name: String,
    pub gene_type: String,
    pub level: Option<i32>,
    pub tags: Vec<String>,
    pub exons: Vec<Exon>,
    pub utrs: Vec<Span>,
}

/// A gene, which is the feature the GTF source is queried for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gene {
    pub contig: String,
    pub span: Span,
    pub transcripts: Vec<Transcript>,
}

/// One GTF line, as the codec's base data reads it.
#[derive(Debug, Clone)]
struct GtfLine {
    contig: String,
    feature: String,
    span: Span,
    strand: Strand,
    attributes: Vec<(String, String)>,
}

impl GtfLine {
    fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
    }
}

/// What reading a GENCODE GTF can go wrong with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GtfError(pub String);

fn parse_gtf_line(line: &str) -> Result<GtfLine, GtfError> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() != 9 {
        return Err(GtfError(format!(
            "Found an invalid number of columns in the given GTF file: {line}"
        )));
    }
    let number = |text: &str| {
        text.parse::<i32>()
            .map_err(|_| GtfError(format!("Could not parse a position: {line}")))
    };
    let strand = match fields[6] {
        "+" => Strand::Positive,
        "-" => Strand::Negative,
        _ => return Err(GtfError(format!("Unsupported strand: {line}"))),
    };
    let mut attributes = Vec::new();
    for part in fields[8].split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (key, value) = part.split_once(' ').unwrap_or((part, ""));
        attributes.push((key.to_string(), value.trim().trim_matches('"').to_string()));
    }
    Ok(GtfLine {
        contig: fields[0].to_string(),
        feature: fields[2].to_string(),
        span: Span::new(number(fields[3])?, number(fields[4])?),
        strand,
        attributes,
    })
}

/// A transcript still being read: itself, its exons, and the leaf features not yet attached.
type OpenTranscript = (Transcript, Vec<Exon>, Vec<GtfLine>);

/// `GencodeGtfCodec.decode` over a whole file: the five header lines, then each gene with its
/// transcripts, each transcript with its exons, and each CDS, codon and UTR attached to the exon
/// that CONTAINS it (a UTR goes to the transcript).
pub fn parse_gencode_gtf(text: &str) -> Result<Vec<Gene>, GtfError> {
    let mut genes: Vec<Gene> = Vec::new();
    let mut current: Option<(Gene, Option<OpenTranscript>)> = None;
    let flush_transcript = |gene: &mut Gene, transcript: Option<OpenTranscript>| {
        if let Some((mut transcript, exons, mut leaves)) = transcript {
            for mut exon in exons {
                leaves.retain(|leaf| {
                    if !exon.span.contains(&leaf.span) {
                        return true;
                    }
                    match leaf.feature.as_str() {
                        "CDS" => exon.cds = Some(leaf.span),
                        "start_codon" => exon.start_codon = Some(leaf.span),
                        "stop_codon" => exon.stop_codon = Some(leaf.span),
                        "UTR" => transcript.utrs.push(leaf.span),
                        _ => {}
                    }
                    false
                });
                transcript.exons.push(exon);
            }
            gene.transcripts.push(transcript);
        }
    };
    for line in text.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let parsed = parse_gtf_line(line)?;
        match parsed.feature.as_str() {
            "gene" => {
                if let Some((mut gene, transcript)) = current.take() {
                    flush_transcript(&mut gene, transcript);
                    genes.push(gene);
                }
                current = Some((
                    Gene {
                        contig: parsed.contig.clone(),
                        span: parsed.span,
                        transcripts: Vec::new(),
                    },
                    None,
                ));
            }
            "transcript" => {
                let Some((gene, transcript)) = current.as_mut() else {
                    return Err(GtfError(format!("A transcript before any gene: {line}")));
                };
                flush_transcript(gene, transcript.take());
                let tags = parsed
                    .attributes
                    .iter()
                    .filter(|(key, _)| key == "tag")
                    .map(|(_, value)| value.clone())
                    .collect();
                *transcript = Some((
                    Transcript {
                        contig: parsed.contig.clone(),
                        span: parsed.span,
                        strand: parsed.strand,
                        transcript_id: parsed.attribute("transcript_id").unwrap_or("").to_string(),
                        transcript_type: parsed
                            .attribute("transcript_type")
                            .unwrap_or("")
                            .to_string(),
                        gene_name: parsed.attribute("gene_name").unwrap_or("").to_string(),
                        gene_type: parsed.attribute("gene_type").unwrap_or("").to_string(),
                        level: parsed.attribute("level").and_then(|v| v.parse().ok()),
                        tags,
                        exons: Vec::new(),
                        utrs: Vec::new(),
                    },
                    Vec::new(),
                    Vec::new(),
                ));
            }
            "exon" => {
                let Some((_, Some((_, exons, _)))) = current.as_mut() else {
                    return Err(GtfError(format!("An exon before any transcript: {line}")));
                };
                exons.push(Exon {
                    span: parsed.span,
                    number: parsed
                        .attribute("exon_number")
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(-1),
                    gene_type: parsed.attribute("gene_type").unwrap_or("").to_string(),
                    cds: None,
                    start_codon: None,
                    stop_codon: None,
                });
            }
            _ => {
                let Some((_, Some((_, _, leaves)))) = current.as_mut() else {
                    return Err(GtfError(format!("A feature before any transcript: {line}")));
                };
                leaves.push(parsed);
            }
        }
    }
    if let Some((mut gene, transcript)) = current.take() {
        flush_transcript(&mut gene, transcript);
        genes.push(gene);
    }
    Ok(genes)
}

// ================================================================================================
// The transcript FASTA.
// ================================================================================================

/// `MappedTranscriptIdInfo`: where a transcript's CDS and UTRs sit in its FASTA record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MappedTranscript {
    pub map_key: String,
    pub coding_sequence_start: i32,
    pub coding_sequence_end: i32,
    pub has_five_prime_utr: bool,
    pub five_prime_utr_start: i32,
    pub five_prime_utr_end: i32,
}

/// The transcript FASTA, by record name and by every `|`-separated field of every name.
#[derive(Debug, Clone, Default)]
pub struct TranscriptFasta {
    pub sequences: BTreeMap<String, Vec<u8>>,
    pub ids: BTreeMap<String, MappedTranscript>,
}

impl TranscriptFasta {
    /// `createTranscriptIdMap` over a FASTA's records: each record's name split at `|`, every
    /// piece mapped to the record's coordinates. A later record claiming a piece wins it.
    pub fn parse(text: &str) -> TranscriptFasta {
        let mut out = TranscriptFasta::default();
        let mut name: Option<String> = None;
        let mut bases: Vec<u8> = Vec::new();
        let finish = |name: Option<String>, bases: &mut Vec<u8>, out: &mut TranscriptFasta| {
            if let Some(name) = name {
                let mut mapped = MappedTranscript {
                    map_key: name.clone(),
                    coding_sequence_start: 0,
                    coding_sequence_end: 0,
                    has_five_prime_utr: false,
                    five_prime_utr_start: 0,
                    five_prime_utr_end: 0,
                };
                let range = |field: &str| -> Option<(i32, i32)> {
                    let (_, numbers) = field.split_once(':')?;
                    let (a, b) = numbers.split_once('-')?;
                    Some((a.parse().ok()?, b.parse().ok()?))
                };
                for field in name.split('|') {
                    if field.len() > 4 && field.starts_with("UTR5:") {
                        if let Some((a, b)) = range(field) {
                            mapped.five_prime_utr_start = a;
                            mapped.five_prime_utr_end = b;
                            mapped.has_five_prime_utr = true;
                        }
                    } else if field.len() > 3 && field.starts_with("CDS:") {
                        if let Some((a, b)) = range(field) {
                            mapped.coding_sequence_start = a;
                            mapped.coding_sequence_end = b;
                        }
                    }
                }
                if mapped.coding_sequence_start == 0 {
                    mapped.coding_sequence_start = 1;
                }
                if mapped.coding_sequence_end == 0 {
                    mapped.coding_sequence_end = bases.len() as i32;
                }
                for field in name.split('|') {
                    out.ids.insert(field.to_string(), mapped.clone());
                }
                out.sequences.insert(name, std::mem::take(bases));
            }
        };
        for line in text.lines() {
            if let Some(header) = line.strip_prefix('>') {
                finish(name.take(), &mut bases, &mut out);
                name = Some(header.split_whitespace().next().unwrap_or("").to_string());
            } else {
                bases.extend(line.trim_end().bytes());
            }
        }
        finish(name.take(), &mut bases, &mut out);
        out
    }

    /// `queryAndPrefetch` over a record, one-based and inclusive.
    fn query(&self, key: &str, start: i32, end: i32) -> String {
        let sequence = &self.sequences[key];
        let from = (start - 1).max(0) as usize;
        let to = (end.max(0) as usize).min(sequence.len());
        String::from_utf8_lossy(&sequence[from.min(to)..to]).into_owned()
    }
}

// ================================================================================================
// The reference.
// ================================================================================================

/// The reference bases a factory reads: one contig, whole, and its length.
pub trait Reference {
    /// The bases over `[start, end]`, one-based, clipped to the contig.
    fn bases(&self, contig: &str, start: i64, end: i64) -> Vec<u8>;
    fn length(&self, contig: &str) -> i64;
}

/// `ReferenceContext.getBases(leading, trailing)` around `[start, end]`: the window expanded and
/// then clipped to the contig.
fn window_bases(
    reference: &dyn Reference,
    contig: &str,
    start: i32,
    end: i32,
    leading: i32,
    trailing: i32,
) -> Vec<u8> {
    let length = reference.length(contig);
    let from = (i64::from(start) - i64::from(leading)).max(1);
    let to = (i64::from(end) + i64::from(trailing)).min(length);
    reference.bases(contig, from, to)
}

fn reverse_complement(bases: &str) -> String {
    bases
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'a' => 't',
            b'c' => 'g',
            b'g' => 'c',
            b't' => 'a',
            other => other as char,
        })
        .collect()
}

// ================================================================================================
// The variant.
// ================================================================================================

/// An allele as the funcotations key on it: bases, or a symbolic name such as `<INS>`.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Allele(pub String);

impl Allele {
    pub fn is_symbolic(&self) -> bool {
        self.0.starts_with('<') || self.0.contains('[') || self.0.contains(']') || self.0 == "*"
    }
    pub fn len(&self) -> usize {
        if self.is_symbolic() {
            0
        } else {
            self.0.len()
        }
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// `getBaseString()`, which is empty for a symbolic allele.
    pub fn bases(&self) -> &str {
        if self.is_symbolic() {
            ""
        } else {
            &self.0
        }
    }
}

/// A variant context as the factory reads it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub reference: Allele,
    pub alternates: Vec<Allele>,
}

impl Variant {
    fn span(&self) -> Span {
        Span::new(self.start, self.end)
    }
}

/// What the factory refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EngineError {
    /// A shape of variant whose branches are not ported.
    Limitation(String),
    /// A `UserException` the reference throws.
    User {
        class: &'static str,
        message: String,
    },
}

/// `FuncotatorUtils.TranscriptCodingSequenceException`, which the factory catches and turns into
/// a `COULD_NOT_DETERMINE` funcotation.
#[derive(Debug)]
enum Failure {
    Coding,
    Engine(EngineError),
}

impl From<EngineError> for Failure {
    fn from(error: EngineError) -> Self {
        Failure::Engine(error)
    }
}

/// A real insertion or deletion: two base alleles of different lengths, which is the shape whose
/// branches are not ported.
fn is_base_indel(reference: &Allele, alternate: &Allele) -> bool {
    !reference.is_symbolic() && !alternate.is_symbolic() && reference.len() != alternate.len()
}

/// `GATKVariantContextUtils.isIndel(Allele, Allele)`: the lengths differ, a symbolic allele
/// counting as none, so a `<DEL>` against one base IS an indel here.
fn is_indel(reference: &Allele, alternate: &Allele) -> bool {
    reference.len() != alternate.len()
}

/// `isInsertion` / `isDeletion` over `Allele.length()`, with the same symbolic rule.
fn is_insertion(reference: &Allele, alternate: &Allele) -> bool {
    reference.len() < alternate.len()
}
fn is_deletion(reference: &Allele, alternate: &Allele) -> bool {
    reference.len() > alternate.len()
}

// ================================================================================================
// The funcotations.
// ================================================================================================

/// `GencodeFuncotation.VariantClassification`, with its default severity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum VariantClassification {
    CouldNotDetermine,
    Intron,
    FivePrimeUtr,
    ThreePrimeUtr,
    Igr,
    FivePrimeFlank,
    ThreePrimeFlank,
    Missense,
    Nonsense,
    Nonstop,
    Silent,
    SpliceSite,
    InFrameDel,
    InFrameIns,
    FrameShiftIns,
    FrameShiftDel,
    StartCodonSnp,
    StartCodonIns,
    StartCodonDel,
    DeNovoStartInFrame,
    DeNovoStartOutFrame,
    Rna,
    Lincrna,
}

impl VariantClassification {
    pub const ALL: [VariantClassification; 23] = [
        VariantClassification::CouldNotDetermine,
        VariantClassification::Intron,
        VariantClassification::FivePrimeUtr,
        VariantClassification::ThreePrimeUtr,
        VariantClassification::Igr,
        VariantClassification::FivePrimeFlank,
        VariantClassification::ThreePrimeFlank,
        VariantClassification::Missense,
        VariantClassification::Nonsense,
        VariantClassification::Nonstop,
        VariantClassification::Silent,
        VariantClassification::SpliceSite,
        VariantClassification::InFrameDel,
        VariantClassification::InFrameIns,
        VariantClassification::FrameShiftIns,
        VariantClassification::FrameShiftDel,
        VariantClassification::StartCodonSnp,
        VariantClassification::StartCodonIns,
        VariantClassification::StartCodonDel,
        VariantClassification::DeNovoStartInFrame,
        VariantClassification::DeNovoStartOutFrame,
        VariantClassification::Rna,
        VariantClassification::Lincrna,
    ];

    pub fn name(self) -> &'static str {
        match self {
            VariantClassification::CouldNotDetermine => "COULD_NOT_DETERMINE",
            VariantClassification::Intron => "INTRON",
            VariantClassification::FivePrimeUtr => "FIVE_PRIME_UTR",
            VariantClassification::ThreePrimeUtr => "THREE_PRIME_UTR",
            VariantClassification::Igr => "IGR",
            VariantClassification::FivePrimeFlank => "FIVE_PRIME_FLANK",
            VariantClassification::ThreePrimeFlank => "THREE_PRIME_FLANK",
            VariantClassification::Missense => "MISSENSE",
            VariantClassification::Nonsense => "NONSENSE",
            VariantClassification::Nonstop => "NONSTOP",
            VariantClassification::Silent => "SILENT",
            VariantClassification::SpliceSite => "SPLICE_SITE",
            VariantClassification::InFrameDel => "IN_FRAME_DEL",
            VariantClassification::InFrameIns => "IN_FRAME_INS",
            VariantClassification::FrameShiftIns => "FRAME_SHIFT_INS",
            VariantClassification::FrameShiftDel => "FRAME_SHIFT_DEL",
            VariantClassification::StartCodonSnp => "START_CODON_SNP",
            VariantClassification::StartCodonIns => "START_CODON_INS",
            VariantClassification::StartCodonDel => "START_CODON_DEL",
            VariantClassification::DeNovoStartInFrame => "DE_NOVO_START_IN_FRAME",
            VariantClassification::DeNovoStartOutFrame => "DE_NOVO_START_OUT_FRAME",
            VariantClassification::Rna => "RNA",
            VariantClassification::Lincrna => "LINCRNA",
        }
    }

    pub fn from_name(name: &str) -> Option<VariantClassification> {
        Self::ALL.into_iter().find(|vc| vc.name() == name)
    }

    pub fn default_severity(self) -> i32 {
        match self {
            VariantClassification::CouldNotDetermine => 99,
            VariantClassification::Intron => 10,
            VariantClassification::FivePrimeUtr | VariantClassification::ThreePrimeUtr => 6,
            VariantClassification::Igr => 20,
            VariantClassification::FivePrimeFlank => 15,
            VariantClassification::ThreePrimeFlank => 16,
            VariantClassification::Missense => 1,
            VariantClassification::Nonsense | VariantClassification::Nonstop => 0,
            VariantClassification::Silent => 5,
            VariantClassification::SpliceSite => 4,
            VariantClassification::InFrameDel | VariantClassification::InFrameIns => 1,
            VariantClassification::FrameShiftIns | VariantClassification::FrameShiftDel => 2,
            VariantClassification::StartCodonSnp
            | VariantClassification::StartCodonIns
            | VariantClassification::StartCodonDel => 3,
            VariantClassification::DeNovoStartInFrame => 1,
            VariantClassification::DeNovoStartOutFrame => 0,
            VariantClassification::Rna | VariantClassification::Lincrna => 4,
        }
    }
}

/// `GencodeFuncotation.VariantType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    Ins,
    Del,
    Snp,
    Dnp,
    Tnp,
    Onp,
    Na,
}

impl VariantType {
    pub fn name(self) -> &'static str {
        match self {
            VariantType::Ins => "INS",
            VariantType::Del => "DEL",
            VariantType::Snp => "SNP",
            VariantType::Dnp => "DNP",
            VariantType::Tnp => "TNP",
            VariantType::Onp => "ONP",
            VariantType::Na => "NA",
        }
    }
}

/// `getVariantType`.
pub fn variant_type(reference: &Allele, alternate: &Allele) -> VariantType {
    if alternate.len() > reference.len() {
        VariantType::Ins
    } else if alternate.0 == "*" || alternate.0 == "." || alternate.is_symbolic() {
        VariantType::Na
    } else if alternate.len() < reference.len() {
        VariantType::Del
    } else {
        match reference.len() {
            1 => VariantType::Snp,
            2 => VariantType::Dnp,
            3 => VariantType::Tnp,
            _ => VariantType::Onp,
        }
    }
}

/// The APPRIS tags, in the enum's declaration order, which is the order the ranks sort in.
const APPRIS_RANKS: [&str; 14] = [
    "appris_principal",
    "appris_principal_1",
    "appris_principal_2",
    "appris_principal_3",
    "appris_principal_4",
    "appris_principal_5",
    "appris_alternative_1",
    "appris_alternative_2",
    "appris_candidate_highest_score",
    "appris_candidate_longest_ccds",
    "appris_candidate_ccds",
    "appris_candidate_longest_seq",
    "appris_candidate_longest",
    "appris_candidate",
];

/// The names of the twenty-two GENCODE fields, after the `<name>_<version>_` prefix.
pub const GENCODE_FIELDS: [&str; 22] = [
    "hugoSymbol",
    "ncbiBuild",
    "chromosome",
    "start",
    "end",
    "variantClassification",
    "secondaryVariantClassification",
    "variantType",
    "refAllele",
    "tumorSeqAllele1",
    "tumorSeqAllele2",
    "genomeChange",
    "annotationTranscript",
    "transcriptStrand",
    "transcriptExon",
    "transcriptPos",
    "cDnaChange",
    "codonChange",
    "proteinChange",
    "gcContent",
    "referenceContext",
    "otherTranscripts",
];

/// `GencodeFuncotation`: the fields, and the serialization overrides a user set.
#[derive(Debug, Clone, PartialEq)]
pub struct GencodeFuncotation {
    pub hugo_symbol: Option<String>,
    pub ncbi_build: Option<String>,
    pub chromosome: Option<String>,
    pub start: i32,
    pub end: i32,
    pub variant_classification: Option<VariantClassification>,
    pub secondary_variant_classification: Option<VariantClassification>,
    pub variant_type: Option<VariantType>,
    pub ref_allele: Option<String>,
    pub tumor_seq_allele2: Option<String>,
    pub genome_change: Option<String>,
    pub annotation_transcript: Option<String>,
    pub transcript_strand: Option<String>,
    pub transcript_exon: Option<i32>,
    pub transcript_start_pos: Option<i32>,
    pub transcript_end_pos: Option<i32>,
    pub cdna_change: Option<String>,
    pub codon_change: Option<String>,
    pub protein_change: Option<String>,
    pub gc_content: Option<f64>,
    pub reference_context: Option<String>,
    pub other_transcripts: Option<Vec<String>>,
    pub data_source_name: String,
    pub version: String,
    pub locus_level: Option<i32>,
    pub appris_rank: Option<usize>,
    pub transcript_length: Option<i32>,
    pub gene_transcript_type: Option<String>,
    /// Overrides, by short field name.
    pub overrides: BTreeMap<String, String>,
}

impl GencodeFuncotation {
    fn empty() -> Self {
        GencodeFuncotation {
            hugo_symbol: None,
            ncbi_build: None,
            chromosome: None,
            start: 0,
            end: 0,
            variant_classification: None,
            secondary_variant_classification: None,
            variant_type: None,
            ref_allele: None,
            tumor_seq_allele2: None,
            genome_change: None,
            annotation_transcript: None,
            transcript_strand: None,
            transcript_exon: None,
            transcript_start_pos: None,
            transcript_end_pos: None,
            cdna_change: None,
            codon_change: None,
            protein_change: None,
            gc_content: None,
            reference_context: None,
            other_transcripts: None,
            data_source_name: String::new(),
            version: String::new(),
            locus_level: None,
            appris_rank: None,
            transcript_length: None,
            gene_transcript_type: None,
            overrides: BTreeMap::new(),
        }
    }

    fn prefix(&self) -> String {
        format!("{}_{}_", self.data_source_name, self.version)
    }

    /// `getFieldNames`.
    pub fn field_names(&self) -> Vec<String> {
        let prefix = self.prefix();
        GENCODE_FIELDS
            .iter()
            .map(|field| format!("{prefix}{field}"))
            .collect()
    }

    /// `getAltAllele`: the tumor allele as bases, which is how the map compares it.
    pub fn alt_allele(&self) -> Allele {
        Allele(self.tumor_seq_allele2.clone().unwrap_or_default())
    }

    fn transcript_pos_string(&self) -> String {
        match (self.transcript_start_pos, self.transcript_end_pos) {
            (None, _) => String::new(),
            (Some(start), Some(end)) if start == end => start.to_string(),
            (Some(start), Some(end)) => format!("{start}_{end}"),
            (Some(start), None) => format!("{start}_null"),
        }
    }

    /// `getField`, the override winning over the value.
    pub fn field(&self, name: &str) -> Option<String> {
        let prefix = self.prefix();
        let short = name.strip_prefix(&prefix).unwrap_or(name);
        if !GENCODE_FIELDS.contains(&short) {
            return None;
        }
        if let Some(value) = self.overrides.get(short) {
            return Some(value.clone());
        }
        let text = |value: &Option<String>| value.clone().unwrap_or_default();
        Some(match short {
            "hugoSymbol" => match &self.hugo_symbol {
                Some(symbol) if !symbol.is_empty() => symbol.clone(),
                _ => fu::UNKNOWN_GENE.to_string(),
            },
            "ncbiBuild" => text(&self.ncbi_build),
            "chromosome" => text(&self.chromosome),
            "start" => self.start.to_string(),
            "end" => self.end.to_string(),
            "variantClassification" => self
                .variant_classification
                .map(|vc| vc.name().to_string())
                .unwrap_or_default(),
            "secondaryVariantClassification" => self
                .secondary_variant_classification
                .map(|vc| vc.name().to_string())
                .unwrap_or_default(),
            "variantType" => self
                .variant_type
                .map(|vt| vt.name().to_string())
                .unwrap_or_default(),
            "refAllele" | "tumorSeqAllele1" => text(&self.ref_allele),
            "tumorSeqAllele2" => text(&self.tumor_seq_allele2),
            "genomeChange" => text(&self.genome_change),
            "annotationTranscript" => text(&self.annotation_transcript),
            "transcriptStrand" => text(&self.transcript_strand),
            "transcriptExon" => self
                .transcript_exon
                .map(|exon| exon.to_string())
                .unwrap_or_default(),
            "transcriptPos" => self.transcript_pos_string(),
            "cDnaChange" => text(&self.cdna_change),
            "codonChange" => text(&self.codon_change),
            "proteinChange" => text(&self.protein_change),
            "gcContent" => self
                .gc_content
                .map(gatk_engine::tsv_table::java_double_to_string)
                .unwrap_or_default(),
            "referenceContext" => text(&self.reference_context),
            "otherTranscripts" => self
                .other_transcripts
                .as_ref()
                .map(|others| others.join("/"))
                .unwrap_or_default(),
            _ => unreachable!("the field list is closed"),
        })
    }

    /// `setFieldSerializationOverrideValue`, which strips the prefix from the given name.
    pub fn set_override(&mut self, name: &str, value: &str) -> Result<(), EngineError> {
        let prefix = self.prefix();
        let short = name.strip_prefix(&prefix).unwrap_or(name);
        if !GENCODE_FIELDS.contains(&short) {
            return Err(EngineError::User {
                class: "org.broadinstitute.hellbender.exceptions.UserException",
                message: format!(
                    "Attempted to override invalid field in this GencodeFuncotation: {name} (value was: {value})"
                ),
            });
        }
        self.overrides.insert(short.to_string(), value.to_string());
        Ok(())
    }

    fn severity(&self, severities: &Severities) -> i32 {
        severities.of(self
            .variant_classification
            .expect("a sorted funcotation carries a classification"))
    }
}

/// `VariantClassification.getSeverity()`, which a custom severity file can reassign.
#[derive(Debug, Clone, Default)]
pub struct Severities {
    pub custom: BTreeMap<VariantClassification, i32>,
}

impl Severities {
    pub fn of(&self, vc: VariantClassification) -> i32 {
        self.custom
            .get(&vc)
            .copied()
            .unwrap_or_else(|| vc.default_severity())
    }
}

/// `TranscriptSelectionMode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TranscriptSelectionMode {
    BestEffect,
    Canonical,
    All,
}

/// The chained comparators the two modes build, over the RAW user transcript set.
fn compare_funcotations(
    mode: TranscriptSelectionMode,
    user: &[String],
    severities: &Severities,
    a: &GencodeFuncotation,
    b: &GencodeFuncotation,
) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let in_list = |f: &GencodeFuncotation| match &f.annotation_transcript {
        Some(id) => {
            let bare = fu::transcript_id_without_version_number(id);
            user.iter()
                .any(|u| fu::transcript_id_without_version_number(u) == bare)
        }
        None => false,
    };
    let by_user = || match (in_list(a), in_list(b)) {
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        _ => Ordering::Equal,
    };
    let is_igr =
        |f: &GencodeFuncotation| f.variant_classification == Some(VariantClassification::Igr);
    let by_igr = || match (is_igr(a), is_igr(b)) {
        (false, true) => Ordering::Less,
        (true, false) => Ordering::Greater,
        _ => Ordering::Equal,
    };
    let by_classification = || a.severity(severities).cmp(&b.severity(severities));
    let coding =
        |f: &GencodeFuncotation| f.gene_transcript_type.as_deref() == Some("protein_coding");
    let by_coding = || match (coding(a), coding(b)) {
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        _ => Ordering::Equal,
    };
    let by_level = || match (a.locus_level, b.locus_level) {
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => Ordering::Equal,
        (Some(x), Some(y)) => x.cmp(&y),
    };
    let by_appris = || match (a.appris_rank, b.appris_rank) {
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (Some(x), Some(y)) if x != y => x.cmp(&y),
        _ => Ordering::Equal,
    };
    let by_length = || match (a.transcript_length, b.transcript_length) {
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (Some(x), Some(y)) if x != y => y.cmp(&x),
        _ => Ordering::Equal,
    };
    let by_name = || match (&a.annotation_transcript, &b.annotation_transcript) {
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => Ordering::Less,
        (Some(x), Some(y)) => java_string_compare(x, y),
    };
    match mode {
        TranscriptSelectionMode::BestEffect => by_user()
            .then_with(by_igr)
            .then_with(by_classification)
            .then_with(by_coding)
            .then_with(by_level)
            .then_with(by_appris)
            .then_with(by_length)
            .then_with(by_name),
        _ => by_user()
            .then_with(by_coding)
            .then_with(by_level)
            .then_with(by_appris)
            .then_with(by_igr)
            .then_with(by_classification)
            .then_with(by_length)
            .then_with(by_name),
    }
}

/// `String.compareTo`, by UTF-16 code unit.
fn java_string_compare(a: &str, b: &str) -> std::cmp::Ordering {
    a.encode_utf16().cmp(b.encode_utf16())
}

// ================================================================================================
// The factory.
// ================================================================================================

/// The factory's settings, from the data source's config and the command line.
#[derive(Debug, Clone)]
pub struct Factory {
    pub name: String,
    pub version: String,
    pub ncbi_build: String,
    pub mode: TranscriptSelectionMode,
    /// The user's transcripts as given, which is what the comparator reads.
    pub user_transcripts: Vec<String>,
    pub five_prime_flank: i32,
    pub three_prime_flank: i32,
    pub splice_site_window: i32,
    pub prefer_mane: bool,
    pub segment_funcotation: bool,
    pub min_bases_for_segment: i32,
    /// The overrides this factory supports, by full field name.
    pub overrides: Vec<(String, String)>,
    pub severities: Severities,
    pub transcripts: TranscriptFasta,
}

const REFERENCE_WINDOW: i32 = 10;
const GC_WINDOW: i32 = 200;
const LEADING_UTR_BASES: i32 = 2;
const TRAILING_UTR_BASES: i32 = 3;

/// The segment fields' suffixes, in the order the metadata declares them.
pub const SEGMENT_SUFFIXES: [&str; 7] = [
    "_genes",
    "_start_gene",
    "_end_gene",
    "_start_exon",
    "_end_exon",
    "_alt_allele",
    "_ref_allele",
];

/// A funcotation of the TABLE kind: named fields, one alternate, a data source.
#[derive(Debug, Clone, PartialEq)]
pub struct TableFuncotation {
    pub fields: Vec<(String, String)>,
    pub alt_allele: Allele,
    pub data_source_name: String,
}

/// Either kind of funcotation.
#[derive(Debug, Clone, PartialEq)]
pub enum Funcotation {
    Gencode(Box<GencodeFuncotation>),
    Table(TableFuncotation),
}

impl Funcotation {
    pub fn field_names(&self) -> Vec<String> {
        match self {
            Funcotation::Gencode(g) => g.field_names(),
            Funcotation::Table(t) => t.fields.iter().map(|(name, _)| name.clone()).collect(),
        }
    }
    pub fn field(&self, name: &str) -> Option<String> {
        match self {
            Funcotation::Gencode(g) => g.field(name),
            Funcotation::Table(t) => t
                .fields
                .iter()
                .find(|(field, _)| field == name)
                .map(|(_, value)| value.clone()),
        }
    }
    pub fn alt_allele(&self) -> Allele {
        match self {
            Funcotation::Gencode(g) => g.alt_allele(),
            Funcotation::Table(t) => t.alt_allele.clone(),
        }
    }
    pub fn data_source_name(&self) -> &str {
        match self {
            Funcotation::Gencode(g) => &g.data_source_name,
            Funcotation::Table(t) => &t.data_source_name,
        }
    }
}

impl Factory {
    /// `getSupportedFuncotationFields`.
    pub fn supported_fields(&self) -> Vec<String> {
        GENCODE_FIELDS
            .iter()
            .map(|field| format!("{}_{}_{field}", self.name, self.version))
            .collect()
    }

    /// `getSupportedFuncotationFieldsForSegments`.
    pub fn segment_fields(&self) -> Vec<String> {
        if !self.segment_funcotation {
            return Vec::new();
        }
        SEGMENT_SUFFIXES
            .iter()
            .map(|suffix| format!("{}_{}{suffix}", self.name, self.version))
            .collect()
    }

    /// `getInfoString`.
    pub fn info_string(&self) -> String {
        let mode = match self.mode {
            TranscriptSelectionMode::BestEffect => "BEST_EFFECT",
            TranscriptSelectionMode::Canonical => "CANONICAL",
            TranscriptSelectionMode::All => "ALL",
        };
        format!("{} {} {mode}", self.name, self.version)
    }

    /// `transformFeatureQueryInterval`: the query padded by the larger flank.
    pub fn query_span(&self, span: Span) -> Span {
        let padding = self.five_prime_flank.max(self.three_prime_flank);
        let start = span.start - padding;
        Span::new(if start < 1 { 1 } else { start }, span.end + padding)
    }

    /// `createFuncotations` for one variant, over the genes its query returned.
    pub fn create_funcotations(
        &self,
        variant: &Variant,
        reference: &dyn Reference,
        genes: &[&Gene],
    ) -> Result<Vec<Funcotation>, EngineError> {
        for alternate in &variant.alternates {
            if is_base_indel(&variant.reference, alternate) {
                return Err(EngineError::Limitation(format!(
                    "the variant at {}:{} is an insertion or a deletion, whose GENCODE branches this port does not carry.",
                    variant.contig, variant.start
                )));
            }
        }
        if genes.is_empty() {
            return self.default_funcotations(variant, reference);
        }
        let mut output: Vec<Funcotation> = if is_segment_variant(
            variant,
            self.min_bases_for_segment,
        ) && self.segment_funcotation
        {
            self.segment_funcotations(variant, reference, genes)?
        } else {
            self.funcotations_on_variant(variant, reference, genes)?
                .into_iter()
                .map(|g| Funcotation::Gencode(Box::new(g)))
                .collect()
        };
        for funcotation in &mut output {
            for (name, value) in &self.overrides {
                match funcotation {
                    Funcotation::Gencode(g) => g.set_override(name, value)?,
                    Funcotation::Table(t) => {
                        if let Some(slot) = t.fields.iter_mut().find(|(field, _)| field == name) {
                            slot.1 = value.clone();
                        } else {
                            return Err(EngineError::User {
                                class: "org.broadinstitute.hellbender.exceptions.GATKException",
                                message: format!(
                                    "Attempted to override a field that is not contained in this TableFuncotation: {name} is not one of [{}]",
                                    t.fields.iter().map(|(f, _)| f.as_str()).collect::<Vec<_>>().join(",")
                                ),
                            });
                        }
                    }
                }
            }
        }
        if output.is_empty() {
            return self.default_funcotations(variant, reference);
        }
        Ok(output)
    }

    /// `createDefaultFuncotationsOnVariant`: a segment with no gene, or an IGR per alternate.
    fn default_funcotations(
        &self,
        variant: &Variant,
        reference: &dyn Reference,
    ) -> Result<Vec<Funcotation>, EngineError> {
        if is_segment_variant(variant, self.min_bases_for_segment) {
            return Ok(self.segment_table(variant, &[], None, None, None, None));
        }
        Ok(variant
            .alternates
            .iter()
            .map(|alternate| {
                Funcotation::Gencode(Box::new(self.igr(variant, alternate, reference)))
            })
            .collect())
    }

    /// `createFuncotationsOnVariant`: per alternate, every transcript's funcotation, or an IGR.
    fn funcotations_on_variant(
        &self,
        variant: &Variant,
        reference: &dyn Reference,
        genes: &[&Gene],
    ) -> Result<Vec<GencodeFuncotation>, EngineError> {
        let mut out = Vec::new();
        for alternate in &variant.alternates {
            let list = self.by_all_transcripts(variant, alternate, reference, genes)?;
            if list.is_empty() {
                out.push(self.igr(variant, alternate, reference));
            } else {
                out.extend(list);
            }
        }
        Ok(out)
    }

    fn sort(&self, list: &mut [GencodeFuncotation]) {
        list.sort_by(|a, b| {
            compare_funcotations(self.mode, &self.user_transcripts, &self.severities, a, b)
        });
    }

    /// `createGencodeFuncotationsByAllTranscripts`.
    fn by_all_transcripts(
        &self,
        variant: &Variant,
        alternate: &Allele,
        reference: &dyn Reference,
        genes: &[&Gene],
    ) -> Result<Vec<GencodeFuncotation>, EngineError> {
        let mut list = Vec::new();
        for gene in genes {
            let transcripts: Vec<&Transcript> = if self.prefer_mane {
                mane_transcripts(&gene.transcripts)
            } else {
                gene.transcripts.iter().filter(|t| is_basic(t)).collect()
            };
            list.extend(self.helper(variant, alternate, reference, &transcripts)?);
        }
        if !list.is_empty() {
            self.sort(&mut list);
            // `populateOtherTranscriptsMapForFuncotation`: a LinkedHashMap keyed by transcript, so a
            // repeated transcript keeps its first place and its last value.
            let mut condensed: Vec<(String, String)> = Vec::new();
            for f in &list {
                let key = f.annotation_transcript.clone().unwrap_or_default();
                let value = condense(f);
                match condensed.iter_mut().find(|(k, _)| *k == key) {
                    Some(slot) => slot.1 = value,
                    None => condensed.push((key, value)),
                }
            }
            for f in &mut list {
                let own = f.annotation_transcript.clone().unwrap_or_default();
                f.other_transcripts = Some(
                    condensed
                        .iter()
                        .filter(|(k, _)| *k != own)
                        .map(|(_, v)| v.clone())
                        .collect(),
                );
            }
        }
        if self.mode != TranscriptSelectionMode::All && !list.is_empty() {
            list.truncate(1);
        }
        Ok(list)
    }

    /// `createFuncotationsHelper` over a transcript list: each one's funcotation, a coding
    /// failure turned into a `COULD_NOT_DETERMINE` one.
    fn helper(
        &self,
        variant: &Variant,
        alternate: &Allele,
        reference: &dyn Reference,
        transcripts: &[&Transcript],
    ) -> Result<Vec<GencodeFuncotation>, EngineError> {
        let mut out = Vec::new();
        for transcript in transcripts {
            match self.single_transcript(variant, alternate, reference, transcript) {
                Ok(Some(f)) => out.push(f),
                Ok(None) => {}
                Err(Failure::Coding) => {
                    out.push(self.problem(variant, alternate, reference, transcript))
                }
                Err(Failure::Engine(error)) => return Err(error),
            }
        }
        Ok(out)
    }

    /// `createDefaultFuncotationsOnProblemVariant`.
    fn problem(
        &self,
        variant: &Variant,
        alternate: &Allele,
        reference: &dyn Reference,
        transcript: &Transcript,
    ) -> GencodeFuncotation {
        let mut f = self.trivial(variant, alternate, transcript);
        f.version = self.version.clone();
        let snippet = reference_snippet(variant, alternate, reference, transcript.strand);
        f.reference_context = Some(snippet.positive());
        f.gc_content = Some(gc_content(variant, alternate, reference));
        f.variant_classification = Some(VariantClassification::CouldNotDetermine);
        f.data_source_name = self.name.clone();
        f
    }

    /// `createGencodeFuncotationBuilderWithTrivialFieldsPopulated`.
    fn trivial(
        &self,
        variant: &Variant,
        alternate: &Allele,
        transcript: &Transcript,
    ) -> GencodeFuncotation {
        let mut f = GencodeFuncotation::empty();
        f.ref_allele = Some(variant.reference.0.clone());
        f.transcript_strand = Some(transcript.strand.symbol().to_string());
        f.hugo_symbol = Some(transcript.gene_name.clone());
        f.ncbi_build = Some(self.ncbi_build.clone());
        f.chromosome = Some(transcript.contig.clone());
        f.start = variant.start;
        f.gene_transcript_type = Some(transcript.transcript_type.clone());
        f.end = variant.end;
        f.variant_type = Some(variant_type(&variant.reference, alternate));
        f.tumor_seq_allele2 = Some(alternate.bases().to_string());
        f.genome_change = Some(genome_change(variant, alternate));
        f.annotation_transcript = Some(transcript.transcript_id.clone());
        f.locus_level = Some(transcript.level.unwrap_or(0));
        f.appris_rank = transcript
            .tags
            .iter()
            .filter_map(|tag| APPRIS_RANKS.iter().position(|rank| rank == tag))
            .min();
        f.transcript_length = Some(transcript.exons.iter().map(|e| e.span.length()).sum());
        // The builder's data source and version stay unset until the caller sets them.
        f
    }

    /// `createGencodeFuncotationOnSingleTranscript`.
    fn single_transcript(
        &self,
        variant: &Variant,
        alternate: &Allele,
        reference: &dyn Reference,
        transcript: &Transcript,
    ) -> Result<Option<GencodeFuncotation>, Failure> {
        let span = variant.span();
        let Some(region) = containing_region(&span, transcript) else {
            if transcript.span.overlaps(&span) {
                return Err(Failure::Coding);
            }
            if is_five_prime_flank(variant, transcript, self.five_prime_flank) {
                return Ok(Some(self.flank(
                    variant,
                    alternate,
                    transcript,
                    reference,
                    VariantClassification::FivePrimeFlank,
                )));
            }
            if is_three_prime_flank(variant, transcript, self.three_prime_flank) {
                return Ok(Some(self.flank(
                    variant,
                    alternate,
                    transcript,
                    reference,
                    VariantClassification::ThreePrimeFlank,
                )));
            }
            return Ok(None);
        };
        if alternate.is_symbolic() || alternate.0 == "*" {
            return Ok(Some(self.symbolic(
                variant,
                alternate,
                transcript,
                reference,
                VariantClassification::CouldNotDetermine,
            )));
        }
        if variant.reference.bases().contains('N') || alternate.bases().contains('N') {
            return Ok(Some(self.masked(variant, alternate, transcript, reference)));
        }
        let start_in_transcript = start_position_in_transcript(
            &span,
            &sorted_cds_and_codons(transcript),
            transcript.strand,
        );
        match region {
            Region::Exon(index) => {
                if start_in_transcript == -1 {
                    return Err(Failure::Coding);
                }
                let exon = &transcript.exons[index];
                if transcript.gene_type == "protein_coding" {
                    self.coding(variant, alternate, reference, transcript, exon)
                        .map(Some)
                } else {
                    Ok(Some(self.non_coding(
                        variant, alternate, reference, transcript, exon,
                    )))
                }
            }
            Region::Utr(index) => Ok(Some(self.utr(
                variant,
                alternate,
                reference,
                transcript,
                transcript.utrs[index],
            )?)),
            Region::Transcript => Ok(Some(self.intron(variant, alternate, reference, transcript))),
        }
    }

    /// `createCodingRegionFuncotationForNonProteinCodingFeature`.
    fn non_coding(
        &self,
        variant: &Variant,
        alternate: &Allele,
        reference: &dyn Reference,
        transcript: &Transcript,
        exon: &Exon,
    ) -> GencodeFuncotation {
        let positions = sorted_cds_and_codons(transcript);
        let mut f = self.trivial(variant, alternate, transcript);
        f.transcript_exon = Some(exon.number);
        f.version = self.version.clone();
        let start = transcript_allele_start(variant, &transcript.exons, transcript.strand);
        set_transcript_position(variant, alternate, start, &mut f);
        let snippet = reference_snippet(variant, alternate, reference, transcript.strand);
        f.reference_context = Some(snippet.positive());
        let overlap = overlapping_exon_positions(variant, &positions);
        f.gc_content = Some(gc_content(variant, alternate, reference));
        let (ref_bases, alt_bases) = match transcript.strand {
            Strand::Positive => (variant.reference.0.clone(), alternate.0.clone()),
            Strand::Negative => (
                reverse_complement(&variant.reference.0),
                reverse_complement(&alternate.0),
            ),
        };
        let _ = overlap;
        f.cdna_change = Some(fu::coding_sequence_change_string_for_xnp(
            start_position_in_transcript(&variant.span(), &positions, transcript.strand),
            &ref_bases,
            &alt_bases,
        ));
        f.variant_classification = Some(lincrna_or_rna(&exon.gene_type));
        f.data_source_name = self.name.clone();
        f
    }

    /// `createCodingRegionFuncotationForProteinCodingFeature`.
    fn coding(
        &self,
        variant: &Variant,
        alternate: &Allele,
        reference: &dyn Reference,
        transcript: &Transcript,
        exon: &Exon,
    ) -> Result<GencodeFuncotation, Failure> {
        let positions = sorted_cds_and_codons(transcript);
        let kind = variant_type(&variant.reference, alternate);
        let mut f = self.trivial(variant, alternate, transcript);
        f.transcript_exon = Some(exon.number);
        f.version = self.version.clone();
        let comparison = sequence_comparison(
            variant,
            alternate,
            reference,
            transcript,
            &positions,
            &self.transcripts,
        )?;
        set_transcript_position(
            variant,
            alternate,
            comparison.transcript_allele_start,
            &mut f,
        );
        f.reference_context = Some(match comparison.strand {
            Strand::Positive => comparison.reference_bases.clone(),
            Strand::Negative => reverse_complement(&comparison.reference_bases),
        });
        f.gc_content = Some(comparison.gc_content);
        f.cdna_change = Some(fu::coding_sequence_change_string_for_xnp(
            comparison.coding_sequence_allele_start,
            &comparison.reference_allele,
            &comparison.alternate_allele,
        ));
        match &comparison.sequence {
            Some(sequence) => {
                f.codon_change = Some(fu::codon_change_string_for_onp(
                    &sequence.aligned_coding_reference,
                    &sequence.aligned_coding_alternate,
                    comparison.aligned_coding_sequence_allele_start,
                    comparison.aligned_reference_allele_stop,
                ));
                f.protein_change = Some(fu::protein_change_string_for_onp(
                    &sequence.protein.ref_aa,
                    &sequence.protein.alt_aa,
                    sequence.protein.start,
                    sequence.protein.end,
                ));
                let class = classification(
                    variant,
                    alternate,
                    kind,
                    exon,
                    transcript.exons.len(),
                    &comparison,
                    self.splice_site_window,
                );
                f.variant_classification = Some(class);
                if class == VariantClassification::SpliceSite {
                    f.secondary_variant_classification =
                        Some(coding_region_classification(&comparison));
                }
            }
            None => f.variant_classification = Some(lincrna_or_rna(&exon.gene_type)),
        }
        f.data_source_name = self.name.clone();
        Ok(f)
    }

    /// `createUtrFuncotation`.
    fn utr(
        &self,
        variant: &Variant,
        alternate: &Allele,
        reference: &dyn Reference,
        transcript: &Transcript,
        utr: Span,
    ) -> Result<GencodeFuncotation, Failure> {
        let mut f = self.trivial(variant, alternate, transcript);
        if let Some(exon) = transcript.exons.iter().find(|e| e.span.contains(&utr)) {
            f.transcript_exon = Some(exon.number);
        }
        f.gc_content = Some(gc_content(variant, alternate, reference));
        let strand = transcript.strand;
        let corrected_alt = match strand {
            Strand::Positive => alternate.0.clone(),
            Strand::Negative => reverse_complement(&alternate.0),
        };
        let snippet = reference_snippet(variant, alternate, reference, strand);
        f.reference_context = Some(snippet.positive());
        if is_five_prime_utr(utr, transcript) {
            f.variant_classification = Some(VariantClassification::FivePrimeUtr);
            if let Some(mapped) = self.transcripts.ids.get(&transcript.transcript_id) {
                let reference_length = variant.reference.len() as i32;
                let trailing = if reference_length < TRAILING_UTR_BASES {
                    TRAILING_UTR_BASES
                } else {
                    reference_length + 1
                };
                let utr_sequence = if mapped.has_five_prime_utr {
                    self.transcripts.query(
                        &mapped.map_key,
                        mapped.five_prime_utr_start,
                        mapped.five_prime_utr_end + trailing,
                    )
                } else {
                    String::new()
                };
                let coding_start = start_position_in_transcript(
                    &variant.span(),
                    &transcript.exons.iter().map(|e| e.span).collect::<Vec<_>>(),
                    strand,
                );
                let indel_offset = if variant.reference.len() != alternate.len() {
                    1
                } else {
                    0
                };
                let front = if strand == Strand::Positive {
                    indel_offset
                } else {
                    0
                };
                let back = if strand == Strand::Negative {
                    indel_offset
                } else {
                    0
                };
                let bases = snippet.bases.as_str();
                let cut = |from: i32, to: i32| -> Result<&str, Failure> {
                    bases.get(from as usize..to as usize).ok_or(Failure::Engine(
                        EngineError::Limitation(
                            "a 5' UTR window past the reference snippet".to_string(),
                        ),
                    ))
                };
                let raw = format!(
                    "{}{}{}",
                    cut(
                        REFERENCE_WINDOW - LEADING_UTR_BASES + front,
                        REFERENCE_WINDOW
                    )?,
                    corrected_alt,
                    cut(
                        REFERENCE_WINDOW + reference_length,
                        REFERENCE_WINDOW + trailing + back
                    )?
                );
                let mut found = false;
                let mut offset = front - LEADING_UTR_BASES;
                let mut i = 0usize;
                while i + 3 < raw.len() {
                    found = fu::eukaryotic_amino_acid_by_codon(&raw[i..i + 3]) == Some("M");
                    if found {
                        offset += i as i32;
                        break;
                    }
                    i += 1;
                }
                if found {
                    f.variant_classification = Some(
                        if in_frame_with_end_of_region(
                            coding_start + offset,
                            utr_sequence.len() as i32,
                        ) {
                            VariantClassification::DeNovoStartInFrame
                        } else {
                            VariantClassification::DeNovoStartOutFrame
                        },
                    );
                }
            }
        } else {
            f.variant_classification = Some(VariantClassification::ThreePrimeUtr);
        }
        f.version = self.version.clone();
        f.data_source_name = self.name.clone();
        Ok(f)
    }

    /// `createIntronFuncotation`.
    fn intron(
        &self,
        variant: &Variant,
        alternate: &Allele,
        reference: &dyn Reference,
        transcript: &Transcript,
    ) -> GencodeFuncotation {
        let (corrected_ref, corrected_alt) = match transcript.strand {
            Strand::Positive => (variant.reference.0.clone(), alternate.0.clone()),
            Strand::Negative => (
                reverse_complement(&variant.reference.0),
                reverse_complement(&alternate.0),
            ),
        };
        let mut f = self.trivial(variant, alternate, transcript);
        let snippet = reference_snippet(variant, alternate, reference, transcript.strand);
        f.reference_context = Some(snippet.positive());
        f.variant_classification = Some(if transcript.gene_type == "protein_coding" {
            VariantClassification::Intron
        } else {
            lincrna_or_rna(&transcript.gene_type)
        });
        f.gc_content = Some(gc_content(variant, alternate, reference));
        if let Some(exon) = splice_site_exon(variant, transcript, self.splice_site_window) {
            f.variant_classification = Some(VariantClassification::SpliceSite);
            f.secondary_variant_classification = Some(VariantClassification::Intron);
            f.codon_change = Some(splice_site_codon_change(
                variant.start,
                exon.number,
                exon.span.start,
                exon.span.end,
                transcript.strand,
                0,
            ));
        }
        let exons: Vec<(i32, i32)> = transcript
            .exons
            .iter()
            .map(|e| (e.span.start, e.span.end))
            .collect();
        f.cdna_change = Some(fu::intronic_cdna_string(
            variant.start,
            &exons,
            &corrected_ref,
            &corrected_alt,
        ));
        f.version = self.version.clone();
        f.data_source_name = self.name.clone();
        f
    }

    /// `createIgrFuncotation`.
    pub fn igr(
        &self,
        variant: &Variant,
        alternate: &Allele,
        reference: &dyn Reference,
    ) -> GencodeFuncotation {
        let mut f = GencodeFuncotation::empty();
        f.gc_content = Some(gc_content(variant, alternate, reference));
        f.variant_classification = Some(VariantClassification::Igr);
        f.ref_allele = Some(variant.reference.0.clone());
        f.tumor_seq_allele2 = Some(if alternate.bases().is_empty() {
            alternate.0.clone()
        } else {
            alternate.bases().to_string()
        });
        f.start = variant.start;
        f.end = variant.end;
        f.variant_type = Some(variant_type(&variant.reference, alternate));
        f.chromosome = Some(variant.contig.clone());
        f.annotation_transcript = Some(fu::NO_TRANSCRIPT.to_string());
        if !alternate.is_symbolic() && alternate.0 != "*" {
            f.genome_change = Some(genome_change(variant, alternate));
        }
        f.ncbi_build = Some(self.ncbi_build.clone());
        f.reference_context =
            Some(reference_snippet(variant, alternate, reference, Strand::Positive).bases);
        f.version = self.version.clone();
        f.data_source_name = self.name.clone();
        f
    }

    /// `createFlankFuncotation`.
    fn flank(
        &self,
        variant: &Variant,
        alternate: &Allele,
        transcript: &Transcript,
        reference: &dyn Reference,
        kind: VariantClassification,
    ) -> GencodeFuncotation {
        if alternate.is_symbolic() || alternate.0 == "*" {
            return self.symbolic(variant, alternate, transcript, reference, kind);
        }
        let mut f = GencodeFuncotation::empty();
        let snippet = reference_snippet(variant, alternate, reference, Strand::Positive);
        f.hugo_symbol = Some(transcript.gene_name.clone());
        f.chromosome = Some(variant.contig.clone());
        f.start = variant.start;
        f.end = variant.end;
        f.transcript_strand = Some(transcript.strand.symbol().to_string());
        f.variant_classification = Some(kind);
        f.variant_type = Some(variant_type(&variant.reference, alternate));
        f.ref_allele = Some(variant.reference.0.clone());
        f.tumor_seq_allele2 = Some(alternate.bases().to_string());
        f.genome_change = Some(genome_change(variant, alternate));
        f.annotation_transcript = Some(transcript.transcript_id.clone());
        f.reference_context = Some(snippet.positive());
        f.gc_content = Some(gc_content(variant, alternate, reference));
        f.ncbi_build = Some(self.ncbi_build.clone());
        f.gene_transcript_type = Some(transcript.transcript_type.clone());
        f.version = self.version.clone();
        f.data_source_name = self.name.clone();
        f
    }

    /// `createFuncotationForSymbolicAltAllele`.
    fn symbolic(
        &self,
        variant: &Variant,
        alternate: &Allele,
        transcript: &Transcript,
        reference: &dyn Reference,
        kind: VariantClassification,
    ) -> GencodeFuncotation {
        let mut f = GencodeFuncotation::empty();
        f.gc_content = Some(gc_content(variant, alternate, reference));
        f.variant_classification = Some(kind);
        f.ref_allele = Some(variant.reference.0.clone());
        f.transcript_strand = Some(Strand::Positive.symbol().to_string());
        f.tumor_seq_allele2 = Some(if alternate.bases().is_empty() {
            alternate.0.clone()
        } else {
            alternate.bases().to_string()
        });
        f.start = variant.start;
        f.end = variant.end;
        f.variant_type = Some(variant_type(&variant.reference, alternate));
        f.chromosome = Some(variant.contig.clone());
        f.other_transcripts = Some(Vec::new());
        f.annotation_transcript = Some(transcript.transcript_id.clone());
        f.hugo_symbol = Some(transcript.gene_name.clone());
        f.ncbi_build = Some(self.ncbi_build.clone());
        f.reference_context =
            Some(reference_snippet(variant, &variant.reference, reference, Strand::Positive).bases);
        f.version = self.version.clone();
        f.data_source_name = self.name.clone();
        f
    }

    /// `createFuncotationForMaskedBases`.
    fn masked(
        &self,
        variant: &Variant,
        alternate: &Allele,
        transcript: &Transcript,
        reference: &dyn Reference,
    ) -> GencodeFuncotation {
        let mut f = self.symbolic(
            variant,
            alternate,
            transcript,
            reference,
            VariantClassification::CouldNotDetermine,
        );
        f.reference_context = Some(
            reference_snippet(variant, &variant.reference, reference, Strand::Positive).positive(),
        );
        f
    }

    // --------------------------------------------------------------------------------------------
    // Segments.
    // --------------------------------------------------------------------------------------------

    /// `createSegmentFuncotations` over the genes the segment's query returned.
    fn segment_funcotations(
        &self,
        segment: &Variant,
        reference: &dyn Reference,
        genes: &[&Gene],
    ) -> Result<Vec<Funcotation>, EngineError> {
        if self.mode == TranscriptSelectionMode::All {
            return Err(EngineError::User {
                class: "java.lang.IllegalArgumentException",
                message: "Cannot create funcotations on segments if the selection mode is ALL"
                    .to_string(),
            });
        }
        let span = segment.span();
        let overlapping: Vec<&Transcript> = genes
            .iter()
            .flat_map(|gene| gene.transcripts.iter())
            .filter(|t| is_basic(t))
            .filter(|t| t.span.overlaps(&span))
            .collect();
        let mut names: Vec<String> = Vec::new();
        for t in &overlapping {
            if !names.contains(&t.gene_name) {
                names.push(t.gene_name.clone());
            }
        }
        names.sort_by(|a, b| java_string_compare(a, b));

        let point = |position: i32| Variant {
            contig: segment.contig.clone(),
            start: position,
            end: position,
            reference: segment.reference.clone(),
            alternates: segment.alternates.clone(),
        };
        let pick = |position: i32| -> Result<(Option<GencodeFuncotation>, Option<Transcript>), EngineError> {
            let sub = point(position);
            let at: Vec<&Transcript> = overlapping.iter().copied().filter(|t| t.span.overlaps(&sub.span())).collect();
            let mut list = self.helper(&sub, &sub.reference, reference, &at)?;
            self.sort(&mut list);
            let first = list.into_iter().next();
            let chosen = first.as_ref().and_then(|f| {
                at.iter()
                    .find(|t| Some(&t.transcript_id) == f.annotation_transcript.as_ref())
                    .map(|t| (*t).clone())
            });
            Ok((first, chosen))
        };
        let (start_f, start_t) = pick(segment.start)?;
        let (end_f, end_t) = pick(segment.end)?;
        Ok(self.segment_table(
            segment,
            &names,
            start_f.as_ref(),
            end_f.as_ref(),
            start_t.as_ref(),
            end_t.as_ref(),
        ))
    }

    /// The seven segment fields, one table funcotation per alternate.
    fn segment_table(
        &self,
        segment: &Variant,
        genes: &[String],
        start: Option<&GencodeFuncotation>,
        end: Option<&GencodeFuncotation>,
        start_transcript: Option<&Transcript>,
        end_transcript: Option<&Transcript>,
    ) -> Vec<Funcotation> {
        let hugo = |f: Option<&GencodeFuncotation>| {
            f.map(|f| {
                f.field(&format!("{}hugoSymbol", f.prefix()))
                    .unwrap_or_default()
            })
            .unwrap_or_default()
        };
        let start_exon = start_transcript
            .map(|t| segment_exon_position(t, segment.span()).0)
            .unwrap_or_default();
        let end_exon = end_transcript
            .map(|t| segment_exon_position(t, segment.span()).1)
            .unwrap_or_default();
        let values = [
            genes.join(","),
            hugo(start),
            hugo(end),
            start_exon,
            end_exon,
            String::new(),
            String::new(),
        ];
        let names = self.segment_fields_always();
        segment
            .alternates
            .iter()
            .map(|alternate| {
                Funcotation::Table(TableFuncotation {
                    fields: names.iter().cloned().zip(values.iter().cloned()).collect(),
                    alt_allele: alternate.clone(),
                    data_source_name: self.name.clone(),
                })
            })
            .collect()
    }

    fn segment_fields_always(&self) -> Vec<String> {
        SEGMENT_SUFFIXES
            .iter()
            .map(|suffix| format!("{}_{}{suffix}", self.name, self.version))
            .collect()
    }
}

/// `isBasic`.
fn is_basic(transcript: &Transcript) -> bool {
    transcript.tags.iter().any(|tag| tag == "basic")
}

/// `retreiveMANESelectModeTranscriptsCriteria`.
fn mane_transcripts(transcripts: &[Transcript]) -> Vec<&Transcript> {
    let plus: Vec<&Transcript> = transcripts
        .iter()
        .filter(|t| t.tags.iter().any(|tag| tag == "MANE_Plus_Clinical"))
        .collect();
    if !plus.is_empty() {
        return plus;
    }
    let select: Vec<&Transcript> = transcripts
        .iter()
        .filter(|t| t.tags.iter().any(|tag| tag == "MANE_Select"))
        .collect();
    if !select.is_empty() {
        return select;
    }
    transcripts.iter().filter(|t| is_basic(t)).collect()
}

/// `condenseGencodeFuncotation`.
fn condense(f: &GencodeFuncotation) -> String {
    if f.variant_classification != Some(VariantClassification::Igr) {
        let mut out = format!(
            "{}_{}_{}",
            f.hugo_symbol.clone().unwrap_or_else(|| "null".to_string()),
            f.annotation_transcript
                .clone()
                .unwrap_or_else(|| "null".to_string()),
            f.variant_classification
                .map(|vc| vc.name())
                .unwrap_or("null")
        );
        let intronic = f.variant_classification == Some(VariantClassification::Intron)
            || f.secondary_variant_classification == Some(VariantClassification::Intron);
        if let Some(protein) = &f.protein_change {
            if !intronic {
                out.push('_');
                out.push_str(protein);
            }
        }
        out
    } else {
        "IGR_ANNOTATON".to_string()
    }
}

/// `convertGeneTranscriptTypeToVariantClassification`.
fn lincrna_or_rna(gene_type: &str) -> VariantClassification {
    if gene_type == "lincRNA" || gene_type == "macro_lncRNA" {
        VariantClassification::Lincrna
    } else {
        VariantClassification::Rna
    }
}

/// `isSegmentVariantContext`.
pub fn is_segment_variant(variant: &Variant, minimum: i32) -> bool {
    let acceptable = variant
        .alternates
        .iter()
        .any(|a| crate::funcotate_segments::ACCEPTABLE_ALTERNATES.contains(&a.0.as_str()));
    acceptable && (variant.end - variant.start + 1) > minimum
}

/// Where in a transcript a variant sits.
enum Region {
    Exon(usize),
    Utr(usize),
    Transcript,
}

/// `getContainingGtfSubfeature`: a UTR it overlaps, then an exon whose CDS contains it or whose
/// codons it overlaps (the LAST such wins), else the transcript itself.
fn containing_region(span: &Span, transcript: &Transcript) -> Option<Region> {
    if !transcript.span.contains(span) {
        return None;
    }
    let mut found: Option<Region> = None;
    for (index, utr) in transcript.utrs.iter().enumerate() {
        if utr.overlaps(span) {
            found = Some(Region::Utr(index));
        }
    }
    for (index, exon) in transcript.exons.iter().enumerate() {
        let hit = exon.cds.is_some_and(|cds| cds.contains(span))
            || exon.start_codon.is_some_and(|c| c.overlaps(span))
            || exon.stop_codon.is_some_and(|c| c.overlaps(span));
        if hit {
            found = Some(Region::Exon(index));
        }
    }
    Some(found.unwrap_or(Region::Transcript))
}

/// `getSortedCdsAndStartStopPositions`: exons by number, each one's start codon (outside its CDS),
/// CDS and stop codon (outside its CDS), reversed by start on the negative strand.
fn sorted_cds_and_codons(transcript: &Transcript) -> Vec<Span> {
    let mut exons: Vec<&Exon> = transcript.exons.iter().collect();
    exons.sort_by_key(|e| e.number);
    let mut out = Vec::new();
    for exon in exons {
        if let Some(cds) = exon.cds {
            if let Some(start) = exon.start_codon {
                if !cds.contains(&start) {
                    out.push(start);
                }
            }
            out.push(cds);
            if let Some(stop) = exon.stop_codon {
                if !cds.contains(&stop) {
                    out.push(stop);
                }
            }
        } else if let Some(start) = exon.start_codon {
            out.push(start);
        } else if let Some(stop) = exon.stop_codon {
            out.push(stop);
        }
    }
    if transcript.strand == Strand::Negative {
        out.sort_by_key(|span| std::cmp::Reverse(span.start));
    }
    out
}

/// `getStartPositionInTranscript`.
fn start_position_in_transcript(variant: &Span, pieces: &[Span], strand: Strand) -> i32 {
    let locus = match strand {
        Strand::Positive => variant.start,
        Strand::Negative => variant.end,
    };
    let mut position = 1;
    for piece in pieces {
        if piece.start <= locus && locus <= piece.end {
            position += match strand {
                Strand::Positive => locus - piece.start,
                Strand::Negative => piece.end - locus,
            };
            return position;
        }
        position += piece.end - piece.start + 1;
    }
    -1
}

/// `getTranscriptAlleleStartPosition`.
fn transcript_allele_start(variant: &Variant, exons: &[Exon], strand: Strand) -> i32 {
    let spans: Vec<Span> = exons.iter().map(|e| e.span).collect();
    transcript_allele_start_spans(variant, &spans, strand)
}

fn transcript_allele_start_spans(variant: &Variant, exons: &[Span], strand: Strand) -> i32 {
    let mut filtered: Vec<Span>;
    let mut position;
    match strand {
        Strand::Positive => {
            filtered = exons
                .iter()
                .copied()
                .filter(|e| e.start <= variant.start)
                .collect();
            filtered.sort_by_key(|e| e.start);
            position = variant.start - filtered.last().map_or(0, |e| e.start) + 1;
        }
        Strand::Negative => {
            filtered = exons
                .iter()
                .copied()
                .filter(|e| e.end >= variant.start)
                .collect();
            filtered.sort_by_key(|span| std::cmp::Reverse(span.start));
            position = filtered.last().map_or(0, |e| e.end) - variant.start + 1;
        }
    }
    for exon in filtered.iter().take(filtered.len().saturating_sub(1)) {
        position += exon.end - exon.start + 1;
    }
    position
}

/// `setTranscriptPosition`.
fn set_transcript_position(
    variant: &Variant,
    alternate: &Allele,
    start: i32,
    f: &mut GencodeFuncotation,
) {
    f.transcript_start_pos = Some(start);
    if variant.reference.len() == 1 {
        f.transcript_end_pos = Some(start);
    } else {
        let indel = if is_indel(&variant.reference, alternate) {
            1
        } else {
            0
        };
        f.transcript_end_pos = Some(start + variant.reference.len() as i32 - indel - 1);
    }
}

/// `getOverlappingExonPositions`: the LAST piece the variant overlaps.
fn overlapping_exon_positions(variant: &Variant, pieces: &[Span]) -> Option<Span> {
    let span = variant.span();
    pieces.iter().rev().find(|p| span.overlaps(p)).copied()
}

/// `getGenomeChangeString`'s substitution branches.
fn genome_change(variant: &Variant, alternate: &Allele) -> String {
    fu::genome_change_string_for_xnp(
        &variant.contig,
        variant.start,
        variant.reference.bases(),
        alternate.bases(),
    )
}

/// `StrandCorrectedReferenceBases`: the bases as stored, and the strand they were stored for.
struct Snippet {
    bases: String,
    strand: Strand,
}

impl Snippet {
    /// `getBaseString(Strand.POSITIVE)`.
    fn positive(&self) -> String {
        match self.strand {
            Strand::Positive => self.bases.clone(),
            Strand::Negative => reverse_complement(&self.bases),
        }
    }
}

/// `createReferenceSnippet` with a window of ten.
fn reference_snippet(
    variant: &Variant,
    alternate: &Allele,
    reference: &dyn Reference,
    strand: Strand,
) -> Snippet {
    let adjustment = if is_indel(&variant.reference, alternate) {
        1
    } else {
        0
    };
    let start = (i64::from(variant.start) - i64::from(REFERENCE_WINDOW) + adjustment).max(1);
    let end = i64::from(variant.end) + i64::from(REFERENCE_WINDOW);
    let bases = String::from_utf8_lossy(&reference.bases(&variant.contig, start, end)).into_owned();
    match strand {
        Strand::Positive => Snippet { bases, strand },
        Strand::Negative => Snippet {
            bases: reverse_complement(&bases),
            strand,
        },
    }
}

/// `calculateGcContent` with a window of two hundred.
fn gc_content(variant: &Variant, alternate: &Allele, reference: &dyn Reference) -> f64 {
    let leading = if is_insertion(&variant.reference, alternate)
        || is_deletion(&variant.reference, alternate)
    {
        GC_WINDOW - 1
    } else {
        GC_WINDOW
    };
    let bases = window_bases(
        reference,
        &variant.contig,
        variant.start,
        variant.end,
        leading,
        GC_WINDOW,
    );
    let gc = bases
        .iter()
        .filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    gc as f64 / bases.len() as f64
}

/// `is5PrimeUtr`.
fn is_five_prime_utr(utr: Span, transcript: &Transcript) -> bool {
    let start_codon = transcript.exons.iter().find_map(|e| e.start_codon);
    match (transcript.strand, start_codon) {
        (Strand::Positive, Some(codon)) => utr.start < codon.start && utr.end < codon.start,
        (Strand::Negative, Some(codon)) => utr.start > codon.end && utr.end > codon.end,
        _ => false,
    }
}

/// `isInFrameWithEndOfRegion`.
fn in_frame_with_end_of_region(start: i32, length: i32) -> bool {
    if start > 0 {
        (length - start + 1) % 3 == 0
    } else {
        let pre = 1 - start;
        ((length + pre) - (start + pre) + 1) % 3 == 0
    }
}

fn is_five_prime_flank(variant: &Variant, transcript: &Transcript, size: i32) -> bool {
    match transcript.strand {
        Strand::Positive => left_flank(variant, transcript, size),
        Strand::Negative => right_flank(variant, transcript, size),
    }
}

fn is_three_prime_flank(variant: &Variant, transcript: &Transcript, size: i32) -> bool {
    match transcript.strand {
        Strand::Positive => right_flank(variant, transcript, size),
        Strand::Negative => left_flank(variant, transcript, size),
    }
}

fn left_flank(variant: &Variant, transcript: &Transcript, size: i32) -> bool {
    variant.contig == transcript.contig
        && variant.end < transcript.span.start
        && transcript.span.start - variant.end <= size
}

fn right_flank(variant: &Variant, transcript: &Transcript, size: i32) -> bool {
    variant.contig == transcript.contig
        && variant.start > transcript.span.end
        && variant.start - transcript.span.end <= size
}

/// `getExonWithinSpliceSiteWindow`, for a substitution, whose changed bases are the variant's own.
fn splice_site_exon<'a>(
    variant: &Variant,
    transcript: &'a Transcript,
    window: i32,
) -> Option<&'a Exon> {
    let changed = variant.span();
    transcript.exons.iter().find(|exon| {
        let near =
            |position: i32| changed.start - window <= position && position <= changed.end + window;
        near(exon.span.start) || near(exon.span.end)
    })
}

/// `createSpliceSiteCodonChange`.
fn splice_site_codon_change(
    variant_start: i32,
    exon_number: i32,
    exon_start: i32,
    exon_end: i32,
    strand: Strand,
    adjustment: i32,
) -> String {
    let mut sign = '-';
    let mut offset = exon_start - variant_start;
    if offset.abs() > (variant_start - exon_end).abs() {
        offset = variant_start - exon_end;
        sign = '+';
    }
    offset = offset.abs();
    if strand == Strand::Negative {
        sign = if sign == '+' { '-' } else { '+' };
    }
    if sign == '+' {
        offset += adjustment;
    } else {
        offset -= adjustment;
    }
    if offset < 0 {
        offset = -offset;
        sign = if sign == '+' { '-' } else { '+' };
    }
    format!("c.e{exon_number}{sign}{offset}")
}

/// `SegmentExonUtils.determineSegmentExonPosition`: the start's and the end's exon strings.
fn segment_exon_position(transcript: &Transcript, segment: Span) -> (String, String) {
    let mut exons: Vec<&Exon> = transcript.exons.iter().collect();
    if transcript.strand == Strand::Negative {
        exons.reverse();
    }
    if exons.is_empty() {
        return (String::new(), String::new());
    }
    let find = |point: i32, is_start: bool| -> i32 {
        // `IntervalTree.minOverlapper`: the overlapping exon with the smallest start.
        let exon = exons
            .iter()
            .filter(|e| e.span.start <= point && point <= e.span.end)
            .min_by_key(|e| e.span.start);
        let extent = Span::new(exons[0].span.start, exons[exons.len() - 1].span.end);
        let is_intron = exon.is_none() && extent.start <= point && point <= extent.end;
        if exon.is_none() && !is_intron {
            -1
        } else if is_intron {
            let increment = match (transcript.strand, is_start) {
                (Strand::Positive, true) => 0,
                (Strand::Positive, false) => -1,
                (Strand::Negative, true) => -1,
                (Strand::Negative, false) => 0,
            };
            let mut result = -1;
            for i in 1..exons.len() {
                let current = exons[i];
                let previous = exons[i - 1];
                let intron = Span::new(previous.span.end + 1, current.span.start - 1);
                if intron.start <= point && point <= intron.end {
                    result = match transcript.strand {
                        Strand::Negative => exons.len() as i32 - i as i32 + increment,
                        Strand::Positive => i as i32 + increment,
                    };
                }
            }
            result
        } else {
            exon.map(|e| e.number - 1).unwrap_or(-1)
        }
    };
    let direction = |is_start: bool| {
        if is_start ^ (transcript.strand == Strand::Positive) {
            "-"
        } else {
            "+"
        }
    };
    let start = find(segment.start, true);
    let end = find(segment.end, false);
    (
        if start != -1 {
            format!("{start}{}", direction(true))
        } else {
            String::new()
        },
        if end != -1 {
            format!("{end}{}", direction(false))
        } else {
            String::new()
        },
    )
}

// ================================================================================================
// The sequence comparison.
// ================================================================================================

struct ProteinChange {
    start: i32,
    end: i32,
    ref_aa: String,
    alt_aa: String,
}

struct SequenceInfo {
    aligned_coding_reference: String,
    aligned_coding_alternate: String,
    protein: ProteinChange,
}

struct SequenceComparison {
    strand: Strand,
    reference_bases: String,
    gc_content: f64,
    reference_allele: String,
    alternate_allele: String,
    transcript_allele_start: i32,
    coding_sequence_allele_start: i32,
    aligned_coding_sequence_allele_start: i32,
    aligned_reference_allele_stop: i32,
    sequence: Option<SequenceInfo>,
}

/// `createSequenceComparison` for a substitution.
fn sequence_comparison(
    variant: &Variant,
    alternate: &Allele,
    reference: &dyn Reference,
    transcript: &Transcript,
    positions: &[Span],
    transcripts: &TranscriptFasta,
) -> Result<SequenceComparison, Failure> {
    let strand = transcript.strand;
    let (ref_allele, alt_allele) = match strand {
        Strand::Positive => (variant.reference.0.clone(), alternate.0.clone()),
        Strand::Negative => (
            reverse_complement(&variant.reference.0),
            reverse_complement(&alternate.0),
        ),
    };
    let snippet = reference_snippet(variant, alternate, reference, strand);
    let gc = gc_content(variant, &Allele(alt_allele.clone()), reference);
    let transcript_allele_start = transcript_allele_start(variant, &transcript.exons, strand);
    let coding_start = start_position_in_transcript(&variant.span(), positions, strand);
    let aligned_start = fu::aligned_position(coding_start);
    let aligned_ref_stop = fu::aligned_end_position(coding_start + ref_allele.len() as i32 - 1);
    let aligned_ref = aligned_ref_allele(
        &snippet.bases,
        REFERENCE_WINDOW,
        &ref_allele,
        coding_start,
        aligned_start,
    )?;
    let aligned_ref_start = coding_start - aligned_start + 1;
    // `getAlternateSequence` over the aligned reference, which is computed for its refusal.
    alternate_sequence(&aligned_ref, aligned_ref_start, &ref_allele, &alt_allele)?;

    let mut sequence = None;
    if let Some(mapped) = transcripts.ids.get(&transcript.transcript_id) {
        // No tail padding: that is for an indel near the transcript's end.
        let raw = transcripts.query(
            &mapped.map_key,
            mapped.coding_sequence_start,
            mapped.coding_sequence_end,
        );
        if (coding_start - 1 + ref_allele.len() as i32) > raw.len() as i32 {
            return Err(Failure::Coding);
        }
        let from = (coding_start - 1) as usize;
        let corrected = format!(
            "{}{}{}",
            &raw[..from],
            ref_allele,
            &raw[from + ref_allele.len()..]
        );
        let aligned_coding_ref = aligned_coding_sequence_allele(
            &corrected,
            aligned_start,
            aligned_ref_stop,
            &ref_allele,
            coding_start,
        )?;
        let aligned_coding_alt = alternate_sequence(
            &aligned_coding_ref,
            aligned_ref_start,
            &ref_allele,
            &alt_allele,
        )?;
        let protein = protein_change(
            &ref_allele,
            &alt_allele,
            coding_start,
            aligned_start,
            &corrected,
        );
        sequence = Some(SequenceInfo {
            aligned_coding_reference: aligned_coding_ref,
            aligned_coding_alternate: aligned_coding_alt,
            protein,
        });
    }
    Ok(SequenceComparison {
        strand,
        reference_bases: snippet.bases,
        gc_content: gc,
        reference_allele: ref_allele,
        alternate_allele: alt_allele,
        transcript_allele_start,
        coding_sequence_allele_start: coding_start,
        aligned_coding_sequence_allele_start: aligned_start,
        aligned_reference_allele_stop: aligned_ref_stop,
        sequence,
    })
}

/// `getAlignedRefAllele` for a substitution: the codon-aligned reference bases out of the snippet,
/// with the given reference allele substituted in when the snippet disagrees with it.
fn aligned_ref_allele(
    snippet: &str,
    padding: i32,
    ref_allele: &str,
    coding_start: i32,
    aligned_start: i32,
) -> Result<String, Failure> {
    let extra = coding_start - aligned_start;
    let start = (padding - extra).max(0);
    let end = start + ((f64::from(extra + ref_allele.len() as i32) / 3.0).ceil() as i32) * 3;
    let cut = |text: &str, from: i32, to: i32| -> Result<String, Failure> {
        text.get(from as usize..to as usize)
            .map(str::to_string)
            .ok_or(Failure::Engine(EngineError::Limitation(
                "a codon window past the reference snippet".to_string(),
            )))
    };
    let aligned = cut(snippet, start, end)?;
    let computed = cut(&aligned, extra, extra + ref_allele.len() as i32)?;
    if computed != ref_allele {
        let substituted = alternate_sequence(snippet, padding + 1, ref_allele, ref_allele)?;
        return cut(&substituted, start, end);
    }
    Ok(aligned)
}

/// `getAlternateSequence`.
fn alternate_sequence(
    sequence: &str,
    allele_start: i32,
    ref_allele: &str,
    alt_allele: &str,
) -> Result<String, Failure> {
    let index = (allele_start - 1).unsigned_abs() as usize;
    let end = index + ref_allele.len();
    if end > sequence.len() {
        return Err(Failure::Engine(EngineError::Limitation(
            "an allele past the end of its sequence".to_string(),
        )));
    }
    Ok(format!(
        "{}{}{}",
        &sequence[..index],
        alt_allele,
        &sequence[end..]
    ))
}

/// `getAlignedCodingSequenceAllele` over the POSITIVE-strand coding sequence.
fn aligned_coding_sequence_allele(
    coding: &str,
    aligned_start: i32,
    aligned_stop: i32,
    ref_allele: &str,
    ref_start: i32,
) -> Result<String, Failure> {
    let slice = |text: &str| -> Result<String, Failure> {
        let start = (aligned_start - 1) as usize;
        let end = aligned_stop as usize;
        if end > text.len() {
            return Err(Failure::Coding);
        }
        Ok(text[start..end].to_string())
    };
    let aligned = slice(coding)?;
    let from = (ref_start - 1) as usize;
    let expected = coding.get(from..from + ref_allele.len()).unwrap_or("");
    if expected != ref_allele {
        let substituted = alternate_sequence(coding, ref_start, ref_allele, ref_allele)?;
        return slice(&substituted);
    }
    Ok(aligned)
}

/// `createAminoAcidSequence`.
fn amino_acids(coding: &str) -> String {
    let max = coding.len() - coding.len() % 3;
    let mut out = String::new();
    let mut i = 0;
    while i < max {
        out.push_str(fu::eukaryotic_amino_acid_by_codon(&coding[i..i + 3]).unwrap_or("?"));
        i += 3;
    }
    out
}

/// `ProteinChangeInfo` for a substitution.
fn protein_change(
    ref_allele: &str,
    alt_allele: &str,
    coding_start: i32,
    aligned_start: i32,
    coding: &str,
) -> ProteinChange {
    let from = (coding_start - 1) as usize;
    let alt_coding = format!(
        "{}{}{}",
        &coding[..from],
        alt_allele,
        &coding[from + ref_allele.len()..]
    );
    let reference_protein = amino_acids(coding);
    let alternate_protein = amino_acids(&alt_coding);
    let r: Vec<char> = reference_protein.chars().collect();
    let a: Vec<char> = alternate_protein.chars().collect();
    let mut start = ((aligned_start - 1) / 3) as usize;
    for i in 0..r.len().max(a.len()) {
        if i >= r.len() || i >= a.len() || r[i] != a[i] {
            start = i;
            break;
        }
    }
    let mut i = start;
    while i < r.len() && i < a.len() && r[i] != a[i] {
        i += 1;
    }
    if start == i {
        let aa: String = r.get(start).map(|c| c.to_string()).unwrap_or_default();
        ProteinChange {
            start: start as i32 + 1,
            end: start as i32 + 1,
            ref_aa: aa.clone(),
            alt_aa: aa,
        }
    } else {
        ProteinChange {
            start: start as i32 + 1,
            end: i as i32,
            ref_aa: r[start..i].iter().collect(),
            alt_aa: a[start..i].iter().collect(),
        }
    }
}

/// `createVariantClassification` for a substitution.
fn classification(
    variant: &Variant,
    alternate: &Allele,
    kind: VariantType,
    exon: &Exon,
    exon_count: usize,
    comparison: &SequenceComparison,
    window: i32,
) -> VariantClassification {
    let _ = alternate;
    let changed = variant.span();
    if let Some(stop) = exon.stop_codon {
        if stop.overlaps(&changed) {
            let aligned = &comparison
                .sequence
                .as_ref()
                .expect("classified with sequence information")
                .aligned_coding_alternate;
            let mut found = false;
            let mut i = 0;
            while i + 2 < aligned.len() {
                if fu::eukaryotic_amino_acid_by_codon(&aligned[i..i + 3]) == Some("*") {
                    found = true;
                    break;
                }
                i += 3;
            }
            if !found {
                return VariantClassification::Nonstop;
            }
        }
    }
    let internal = exon.number != 1 && exon.number != exon_count as i32;
    let strand = comparison.strand;
    let left_check = internal
        || (strand == Strand::Negative && exon.number == 1)
        || (strand == Strand::Positive && exon.number == exon_count as i32);
    let right_check = internal
        || (strand == Strand::Positive && exon.number == 1)
        || (strand == Strand::Negative && exon.number == exon_count as i32);
    let adjuster = if window > 0 { 1 } else { 0 };
    let left = left_check
        && Span::new(
            exon.span.start - window,
            exon.span.start + (window - adjuster),
        )
        .overlaps(&changed);
    let right = right_check
        && Span::new(
            exon.span.end - window + 1,
            exon.span.end + (window - adjuster) + 1,
        )
        .overlaps(&changed);
    if left || right {
        VariantClassification::SpliceSite
    } else if exon.start_codon.is_some_and(|c| c.overlaps(&changed)) {
        match kind {
            VariantType::Ins => VariantClassification::StartCodonIns,
            VariantType::Del => VariantClassification::StartCodonDel,
            _ => VariantClassification::StartCodonSnp,
        }
    } else {
        coding_region_classification(comparison)
    }
}

/// `getVarClassFromEqualLengthCodingRegions`.
fn coding_region_classification(comparison: &SequenceComparison) -> VariantClassification {
    let protein = &comparison
        .sequence
        .as_ref()
        .expect("classified with sequence information")
        .protein;
    let r: Vec<char> = protein.ref_aa.chars().collect();
    let mut class = VariantClassification::Silent;
    for (i, alt) in protein.alt_aa.chars().enumerate() {
        if Some(&alt) != r.get(i) {
            if alt == '*' {
                return VariantClassification::Nonsense;
            }
            class = VariantClassification::Missense;
        }
    }
    class
}

// ================================================================================================
// The funcotation map.
// ================================================================================================

/// `FuncotationMap.NO_TRANSCRIPT_AVAILABLE_KEY`.
pub const NO_TRANSCRIPT_KEY: &str = "no_transcript";

/// `FuncotationMap`: funcotations by transcript, each list a `LinkedHashSet`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FuncotationMap {
    pub entries: Vec<(String, Vec<Funcotation>)>,
}

impl FuncotationMap {
    /// `createFromGencodeFuncotations`: each under its own transcript.
    pub fn from_gencode(list: &[Funcotation]) -> FuncotationMap {
        let mut map = FuncotationMap::default();
        for f in list {
            let key = match f {
                Funcotation::Gencode(g) => g.annotation_transcript.clone().unwrap_or_default(),
                Funcotation::Table(_) => NO_TRANSCRIPT_KEY.to_string(),
            };
            map.add(&key, std::slice::from_ref(f));
        }
        map
    }

    /// `createNoTranscriptInfo`.
    pub fn no_transcript(list: &[Funcotation]) -> FuncotationMap {
        let mut map = FuncotationMap::default();
        map.add(NO_TRANSCRIPT_KEY, list);
        map
    }

    pub fn add(&mut self, transcript: &str, list: &[Funcotation]) {
        let index = match self.entries.iter().position(|(key, _)| key == transcript) {
            Some(index) => index,
            None => {
                self.entries.push((transcript.to_string(), Vec::new()));
                self.entries.len() - 1
            }
        };
        for f in list {
            if !self.entries[index].1.contains(f) {
                self.entries[index].1.push(f.clone());
            }
        }
    }

    pub fn transcripts(&self) -> Vec<String> {
        self.entries.iter().map(|(key, _)| key.clone()).collect()
    }

    pub fn get(&self, transcript: &str) -> &[Funcotation] {
        self.entries
            .iter()
            .find(|(key, _)| key == transcript)
            .map(|(_, list)| list.as_slice())
            .unwrap_or(&[])
    }

    /// `getFieldNames(txId)`.
    pub fn field_names(&self, transcript: &str) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for f in self.get(transcript) {
            for name in f.field_names() {
                if !out.contains(&name) {
                    out.push(name);
                }
            }
        }
        out
    }

    /// `getAlleles(txId)`, in first-seen order.
    pub fn alleles(&self, transcript: &str) -> Vec<Allele> {
        let mut out: Vec<Allele> = Vec::new();
        for f in self.get(transcript) {
            let allele = f.alt_allele();
            if !out.contains(&allele) {
                out.push(allele);
            }
        }
        out
    }

    /// `getFieldValue`: the one value this field has for this allele, or a refusal when two
    /// funcotations disagree.
    pub fn field_value(
        &self,
        transcript: &str,
        field: &str,
        allele: &Allele,
    ) -> Result<Option<String>, EngineError> {
        let mut values: Vec<String> = Vec::new();
        for f in self.get(transcript) {
            if f.alt_allele() != *allele {
                continue;
            }
            if let Some(value) = f.field(field) {
                if !values.contains(&value) {
                    values.push(value);
                }
            }
        }
        if values.len() > 1 {
            return Err(EngineError::User {
                class: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
                message: format!(
                    "Bad input: Found more than one unique value for the tuple {{{transcript}, {}, {field}}}: {}",
                    allele.0,
                    values.join(", ")
                ),
            });
        }
        Ok(values.pop())
    }

    fn all_field_names(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for (key, _) in &self.entries {
            for name in self.field_names(key) {
                if !out.contains(&name) {
                    out.push(name);
                }
            }
        }
        out
    }

    /// `doAllTxAlleleCombinationsHaveTheSameFields`.
    pub fn same_fields_everywhere(&self) -> bool {
        let mut all = self.all_field_names();
        all.sort();
        self.entries.iter().all(|(key, list)| {
            self.alleles(key).iter().all(|allele| {
                let mut names: Vec<String> = Vec::new();
                for f in list.iter().filter(|f| f.alt_allele() == *allele) {
                    for name in f.field_names() {
                        if !names.contains(&name) {
                            names.push(name);
                        }
                    }
                }
                names.sort();
                names == all
            })
        })
    }
}

/// A table funcotation over the given fields, for one allele.
pub fn table(fields: &[(String, String)], allele: &Allele, source: &str) -> Funcotation {
    Funcotation::Table(TableFuncotation {
        fields: fields.to_vec(),
        alt_allele: allele.clone(),
        data_source_name: source.to_string(),
    })
}

// ================================================================================================
// The VCF rendering.
// ================================================================================================

/// `VcfOutputRenderer.write`'s `FUNCOTATION` value: per alternate, per transcript, the bracketed
/// funcotations joined by `|`, transcripts joined by `#` and alternates by `,`.
pub fn vcf_funcotation(
    alternates: &[Allele],
    existing: Option<&str>,
    map: &FuncotationMap,
    manual: &[(String, String)],
    included: &[String],
) -> String {
    let existing: Vec<&str> = existing
        .map(|text| text.split(',').collect())
        .unwrap_or_default();
    let mut out = String::new();
    for (index, alternate) in alternates.iter().enumerate() {
        if let Some(previous) = existing.get(index) {
            out.push_str(previous);
            out.push('|');
        }
        for transcript in map.transcripts() {
            out.push('[');
            let manual_f = table(manual, alternate, "UnaccountedManualAnnotations");
            let rendered: Vec<String> = map
                .get(&transcript)
                .iter()
                .chain(std::iter::once(&manual_f))
                .filter(|f| f.alt_allele() == *alternate)
                .filter(|f| !f.field_names().is_empty())
                .filter(|f| f.data_source_name() != "INPUT_VCF")
                .map(|f| {
                    let fields: Vec<(String, String)> = f
                        .field_names()
                        .into_iter()
                        .map(|name| {
                            let value = f.field(&name).unwrap_or_default();
                            (name, value)
                        })
                        .collect();
                    fu::render_sanitized_funcotation_for_vcf(&fields, included)
                })
                .collect();
            out.push_str(&rendered.join("|"));
            out.push_str("]#");
        }
        out.pop();
        out.push(',');
    }
    out.pop();
    out
}

// ================================================================================================
// The MAF rendering.
// ================================================================================================

/// `outputFieldNameMap`: each MAF column and the funcotation fields it takes, in order.
pub const MAF_ALIASES: &[(&str, &[&str])] = &[
    (
        "Hugo_Symbol",
        &[
            "Hugo_Symbol",
            "Gencode_19_hugoSymbol",
            "Gencode_27_hugoSymbol",
            "Gencode_28_hugoSymbol",
            "Gencode_34_hugoSymbol",
            "Gencode_43_hugoSymbol",
            "gene",
            "Gene",
        ],
    ),
    (
        "Entrez_Gene_Id",
        &[
            "Entrez_Gene_Id",
            "HGNC_Entrez_Gene_ID",
            "HGNC_Entrez Gene ID",
            "HGNC_Entrez_Gene_ID(supplied_by_NCBI)",
            "HGNC_Entrez Gene ID(supplied by NCBI)",
            "entrez_id",
            "gene_id",
        ],
    ),
    ("Center", &["Center", "center"]),
    (
        "NCBI_Build",
        &[
            "NCBI_Build",
            "Gencode_19_ncbiBuild",
            "Gencode_27_ncbiBuild",
            "Gencode_28_ncbiBuild",
            "Gencode_34_ncbiBuild",
            "Gencode_43_ncbiBuild",
            "ncbi_build",
        ],
    ),
    (
        "Chromosome",
        &[
            "Chromosome",
            "Gencode_19_chromosome",
            "Gencode_27_chromosome",
            "Gencode_28_chromosome",
            "Gencode_34_chromosome",
            "Gencode_43_chromosome",
            "chr",
            "contig",
            "chromosome",
            "chrom",
            "Chrom",
        ],
    ),
    (
        "Start_Position",
        &[
            "Start_Position",
            "Start_position",
            "Gencode_19_start",
            "Gencode_27_start",
            "Gencode_28_start",
            "Gencode_34_start",
            "Gencode_43_start",
            "start",
            "Start",
            "start_pos",
            "pos",
        ],
    ),
    (
        "End_Position",
        &[
            "End_Position",
            "End_position",
            "Gencode_19_end",
            "Gencode_27_end",
            "Gencode_28_end",
            "Gencode_34_end",
            "Gencode_43_end",
            "end",
            "End",
            "end_pos",
        ],
    ),
    ("Strand", &["Strand"]),
    (
        "Variant_Classification",
        &[
            "Variant_Classification",
            "Gencode_19_variantClassification",
            "Gencode_27_variantClassification",
            "Gencode_28_variantClassification",
            "Gencode_34_variantClassification",
            "Gencode_43_variantClassification",
            "variant_classification",
        ],
    ),
    (
        "Variant_Type",
        &[
            "Variant_Type",
            "Gencode_19_variantType",
            "Gencode_27_variantType",
            "Gencode_28_variantType",
            "Gencode_34_variantType",
            "Gencode_43_variantType",
            "variant_type",
        ],
    ),
    (
        "Reference_Allele",
        &[
            "Reference_Allele",
            "Gencode_19_refAllele",
            "Gencode_27_refAllele",
            "Gencode_28_refAllele",
            "Gencode_34_refAllele",
            "Gencode_43_refAllele",
            "ref",
            "ref_allele",
            "reference_allele",
        ],
    ),
    (
        "Tumor_Seq_Allele1",
        &[
            "Tumor_Seq_Allele1",
            "Gencode_19_tumorSeqAllele1",
            "Gencode_27_tumorSeqAllele1",
            "Gencode_28_tumorSeqAllele1",
            "Gencode_34_tumorSeqAllele1",
            "Gencode_43_tumorSeqAllele1",
            "ref",
            "ref_allele",
            "reference_allele",
        ],
    ),
    (
        "Tumor_Seq_Allele2",
        &[
            "Tumor_Seq_Allele2",
            "Gencode_19_tumorSeqAllele2",
            "Gencode_27_tumorSeqAllele2",
            "Gencode_28_tumorSeqAllele2",
            "Gencode_34_tumorSeqAllele2",
            "Gencode_43_tumorSeqAllele2",
            "alt",
            "alt_allele",
            "alt2",
            "alt_allele2",
            "alternate_allele2",
            "observed_allele2",
            "alternate_allele",
            "observed_allele",
            "alt1",
            "alt_allele1",
            "alternate_allele1",
            "observed_allele1",
        ],
    ),
    ("dbSNP_RS", &["dbSNP_RS", "dbsnp_rs", "dbSNP_RSPOS"]),
    (
        "dbSNP_Val_Status",
        &[
            "dbSNP_Val_Status",
            "custom_dbsnp_val_status",
            "dbsnp_val_status",
            "dbSNP_VLD",
        ],
    ),
    (
        "Tumor_Sample_Barcode",
        &[
            "Tumor_Sample_Barcode",
            "tumor_barcode",
            "tumor_id",
            "case_barcode",
            "case_id",
            "tumor_name",
        ],
    ),
    (
        "Matched_Norm_Sample_Barcode",
        &[
            "Matched_Norm_Sample_Barcode",
            "normal_barcode",
            "normal_id",
            "control_barcode",
            "control_id",
            "normal_name",
            "sample_name",
        ],
    ),
    (
        "Match_Norm_Seq_Allele1",
        &["Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele1"],
    ),
    (
        "Match_Norm_Seq_Allele2",
        &["Match_Norm_Seq_Allele2", "Match_Norm_Seq_Allele2"],
    ),
    (
        "Tumor_Validation_Allele1",
        &["Tumor_Validation_Allele1", "Tumor_Validation_Allele1"],
    ),
    (
        "Tumor_Validation_Allele2",
        &["Tumor_Validation_Allele2", "Tumor_Validation_Allele2"],
    ),
    (
        "Match_Norm_Validation_Allele1",
        &[
            "Match_Norm_Validation_Allele1",
            "Match_Norm_Validation_Allele1",
        ],
    ),
    (
        "Match_Norm_Validation_Allele2",
        &[
            "Match_Norm_Validation_Allele2",
            "Match_Norm_Validation_Allele2",
        ],
    ),
    (
        "Verification_Status",
        &["Verification_Status", "Verification_Status"],
    ),
    (
        "Validation_Status",
        &["Validation_Status", "validation_status"],
    ),
    ("Mutation_Status", &["Mutation_Status", "status"]),
    ("Sequencing_Phase", &["Sequencing_Phase", "phase"]),
    ("Sequence_Source", &["Sequence_Source", "source"]),
    (
        "Validation_Method",
        &["Validation_Method", "Validation_Method"],
    ),
    ("Score", &["Score", "Score"]),
    ("BAM_File", &["BAM_File", "BAM_file", "bam", "bam_file"]),
    ("Sequencer", &["Sequencer", "sequencer", "platform"]),
    (
        "Tumor_Sample_UUID",
        &[
            "Tumor_Sample_UUID",
            "tumor_uuid",
            "case_uuid",
            "tumor_barcode",
            "tumor_id",
            "case_barcode",
            "case_id",
            "tumor_name",
            "Tumor_Sample_Barcode",
        ],
    ),
    (
        "Matched_Norm_Sample_UUID",
        &[
            "Matched_Norm_Sample_UUID",
            "normal_uuid",
            "control_uuid",
            "normal_barcode",
            "normal_id",
            "control_barcode",
            "control_id",
            "normal_name",
            "sample_name",
            "Matched_Norm_Sample_Barcode",
        ],
    ),
    (
        "Genome_Change",
        &[
            "Genome_Change",
            "Gencode_19_genomeChange",
            "Gencode_27_genomeChange",
            "Gencode_28_genomeChange",
            "Gencode_34_genomeChange",
            "Gencode_43_genomeChange",
            "genome_change",
        ],
    ),
    (
        "Annotation_Transcript",
        &[
            "Annotation_Transcript",
            "Gencode_19_annotationTranscript",
            "Gencode_27_annotationTranscript",
            "Gencode_28_annotationTranscript",
            "Gencode_34_annotationTranscript",
            "Gencode_43_annotationTranscript",
            "annotation_transcript",
            "transcript_id",
        ],
    ),
    (
        "Transcript_Strand",
        &[
            "Transcript_Strand",
            "Gencode_19_transcriptStrand",
            "Gencode_27_transcriptStrand",
            "Gencode_28_transcriptStrand",
            "Gencode_34_transcriptStrand",
            "Gencode_43_transcriptStrand",
            "transcript_strand",
        ],
    ),
    (
        "Transcript_Exon",
        &[
            "Transcript_Exon",
            "Gencode_19_transcriptExon",
            "Gencode_27_transcriptExon",
            "Gencode_28_transcriptExon",
            "Gencode_34_transcriptExon",
            "Gencode_43_transcriptExon",
            "transcript_exon",
        ],
    ),
    (
        "Transcript_Position",
        &[
            "Transcript_Position",
            "Gencode_19_transcriptPos",
            "Gencode_27_transcriptPos",
            "Gencode_28_transcriptPos",
            "Gencode_34_transcriptPos",
            "Gencode_43_transcriptPos",
            "transcript_position",
        ],
    ),
    (
        "cDNA_Change",
        &[
            "cDNA_Change",
            "Gencode_19_cDnaChange",
            "Gencode_27_cDnaChange",
            "Gencode_28_cDnaChange",
            "Gencode_34_cDnaChange",
            "Gencode_43_cDnaChange",
            "transcript_change",
        ],
    ),
    (
        "Codon_Change",
        &[
            "Codon_Change",
            "Gencode_19_codonChange",
            "Gencode_27_codonChange",
            "Gencode_28_codonChange",
            "Gencode_34_codonChange",
            "Gencode_43_codonChange",
            "codon_change",
        ],
    ),
    (
        "Protein_Change",
        &[
            "Protein_Change",
            "Gencode_19_proteinChange",
            "Gencode_27_proteinChange",
            "Gencode_28_proteinChange",
            "Gencode_34_proteinChange",
            "Gencode_43_proteinChange",
            "protein_change",
        ],
    ),
    (
        "Other_Transcripts",
        &[
            "Other_Transcripts",
            "Gencode_19_otherTranscripts",
            "Gencode_27_otherTranscripts",
            "Gencode_28_otherTranscripts",
            "Gencode_34_otherTranscripts",
            "Gencode_43_otherTranscripts",
            "other_transcripts",
        ],
    ),
    (
        "Refseq_mRNA_Id",
        &[
            "Refseq_mRNA_Id",
            "Gencode_XRefSeq_mRNA_id",
            "gencode_xref_refseq_mRNA_id",
            "ENSEMBL_RefSeq_mRNA_accession",
            "RefSeq_mRNA_Id",
            "HGNC_RefSeq IDs",
        ],
    ),
    (
        "Refseq_prot_Id",
        &[
            "Refseq_prot_Id",
            "Gencode_XRefSeq_prot_acc",
            "gencode_xref_refseq_prot_acc",
            "ENSEMBL_RefSeq_protein_accession",
            "RefSeq_prot_Id",
        ],
    ),
    (
        "SwissProt_acc_Id",
        &[
            "SwissProt_acc_Id",
            "Simple_Uniprot_uniprot_accession",
            "uniprot_accession",
            "UniProt_uniprot_accession",
        ],
    ),
    (
        "SwissProt_entry_Id",
        &[
            "SwissProt_entry_Id",
            "Simple_Uniprot_uniprot_entry_name",
            "uniprot_entry_name",
            "UniProt_uniprot_entry_name",
        ],
    ),
    (
        "Description",
        &[
            "Description",
            "RefSeq_Description",
            "HGNC_Approved_Name",
            "HGNC_Approved Name",
        ],
    ),
    (
        "UniProt_AApos",
        &["UniProt_AApos", "UniProt_AAxform_aapos", "uniprot_AA_pos"],
    ),
    ("UniProt_Region", &["UniProt_Region", "UniProt_AA_region"]),
    ("UniProt_Site", &["UniProt_Site", "UniProt_AA_site"]),
    (
        "UniProt_Natural_Variations",
        &["UniProt_Natural_Variations", "UniProt_AA_natural_variation"],
    ),
    (
        "UniProt_Experimental_Info",
        &["UniProt_Experimental_Info", "UniProt_AA_experimental_info"],
    ),
    (
        "GO_Biological_Process",
        &[
            "GO_Biological_Process",
            "Simple_Uniprot_GO_Biological_Process",
            "UniProt_GO_Biological_Process",
        ],
    ),
    (
        "GO_Cellular_Component",
        &[
            "GO_Cellular_Component",
            "Simple_Uniprot_GO_Cellular_Component",
            "UniProt_GO_Cellular_Component",
        ],
    ),
    (
        "GO_Molecular_Function",
        &[
            "GO_Molecular_Function",
            "Simple_Uniprot_GO_Molecular_Function",
            "UniProt_GO_Molecular_Function",
        ],
    ),
    (
        "COSMIC_overlapping_mutations",
        &[
            "COSMIC_overlapping_mutations",
            "Cosmic_overlapping_mutations",
            "COSMIC_overlapping_mutations",
            "COSMIC_overlapping_mutation_AAs",
        ],
    ),
    (
        "COSMIC_fusion_genes",
        &[
            "COSMIC_fusion_genes",
            "CosmicFusion_fusion_genes",
            "COSMIC_FusionGenes_fusion_genes",
        ],
    ),
    (
        "COSMIC_tissue_types_affected",
        &[
            "COSMIC_tissue_types_affected",
            "CosmicTissue_tissue_types_affected",
            "COSMIC_tissue_types_affected",
            "COSMIC_Tissue_tissue_types_affected",
        ],
    ),
    (
        "COSMIC_total_alterations_in_gene",
        &[
            "COSMIC_total_alterations_in_gene",
            "CosmicTissue_total_alterations_in_gene",
            "COSMIC_total_alterations_in_gene",
            "COSMIC_Tissue_total_alterations_in_gene",
        ],
    ),
    (
        "Tumorscape_Amplification_Peaks",
        &[
            "Tumorscape_Amplification_Peaks",
            "TUMORScape_Amplification_Peaks",
        ],
    ),
    (
        "Tumorscape_Deletion_Peaks",
        &["Tumorscape_Deletion_Peaks", "TUMORScape_Deletion_Peaks"],
    ),
    (
        "TCGAscape_Amplification_Peaks",
        &[
            "TCGAscape_Amplification_Peaks",
            "TCGAScape_Amplification_Peaks",
        ],
    ),
    (
        "TCGAscape_Deletion_Peaks",
        &["TCGAscape_Deletion_Peaks", "TCGAScape_Deletion_Peaks"],
    ),
    (
        "DrugBank",
        &["DrugBank", "Simple_Uniprot_DrugBank", "UniProt_DrugBank"],
    ),
    (
        "ref_context",
        &[
            "ref_context",
            "Gencode_19_referenceContext",
            "Gencode_27_referenceContext",
            "Gencode_28_referenceContext",
            "Gencode_34_referenceContext",
            "Gencode_43_referenceContext",
            "ref_context",
        ],
    ),
    (
        "gc_content",
        &[
            "gc_content",
            "Gencode_19_gcContent",
            "Gencode_27_gcContent",
            "Gencode_28_gcContent",
            "Gencode_34_gcContent",
            "Gencode_43_gcContent",
            "gc_content",
        ],
    ),
    (
        "CCLE_ONCOMAP_overlapping_mutations",
        &[
            "CCLE_ONCOMAP_overlapping_mutations",
            "CCLE_By_GP_overlapping_mutations",
        ],
    ),
    (
        "CCLE_ONCOMAP_total_mutations_in_gene",
        &[
            "CCLE_ONCOMAP_total_mutations_in_gene",
            "CCLE_By_Gene_total_mutations_in_gene",
        ],
    ),
    (
        "CGC_Mutation_Type",
        &["CGC_Mutation_Type", "CGC_Mutation Type"],
    ),
    (
        "CGC_Translocation_Partner",
        &["CGC_Translocation_Partner", "CGC_Translocation Partner"],
    ),
    (
        "CGC_Tumor_Types_Somatic",
        &[
            "CGC_Tumor_Types_Somatic",
            "CGC_Tumour Types  (Somatic Mutations)",
            "CGC_Tumour_Types__(Somatic_Mutations)",
        ],
    ),
    (
        "CGC_Tumor_Types_Germline",
        &[
            "CGC_Tumor_Types_Germline",
            "CGC_Tumour Types (Germline Mutations)",
            "CGC_Tumour_Types_(Germline_Mutations)",
        ],
    ),
    (
        "CGC_Other_Diseases",
        &[
            "CGC_Other_Diseases",
            "CGC_Other Syndrome/Disease",
            "CGC_Other_Syndrome/Disease",
        ],
    ),
    (
        "DNARepairGenes_Activity_linked_to_OMIM",
        &["DNARepairGenes_Activity_linked_to_OMIM"],
    ),
    (
        "FamilialCancerDatabase_Syndromes",
        &[
            "FamilialCancerDatabase_Syndromes",
            "Familial_Cancer_Genes_Syndrome",
        ],
    ),
    (
        "MUTSIG_Published_Results",
        &[
            "MUTSIG_Published_Results",
            "MutSig Published Results_Published_Results",
        ],
    ),
    (
        "OREGANNO_ID",
        &["OREGANNO_ID", "Oreganno_ID", "ORegAnno_ID"],
    ),
    (
        "OREGANNO_Values",
        &["OREGANNO_Values", "Oreganno_Values", "ORegAnno_Values"],
    ),
    ("tumor_f", &["tumor_f", "sample_allelic_fraction"]),
    ("t_alt_count", &["t_alt_count"]),
    ("t_ref_count", &["t_ref_count"]),
    ("n_alt_count", &["n_alt_count"]),
    ("n_ref_count", &["n_ref_count"]),
];

/// `MafOutputRendererConstants.VariantClassificationMap`.
fn maf_variant_classification(value: &str) -> Option<&'static str> {
    Some(match value {
        "IN_FRAME_DEL" => "In_Frame_Del",
        "IN_FRAME_INS" => "In_Frame_Ins",
        "FRAME_SHIFT_INS" => "Frame_Shift_Ins",
        "FRAME_SHIFT_DEL" => "Frame_Shift_Del",
        "MISSENSE" => "Missense_Mutation",
        "NONSENSE" => "Nonsense_Mutation",
        "SILENT" => "Silent",
        "SPLICE_SITE" => "Splice_Site",
        "START_CODON_DEL" => "Translation_Start_Site",
        "NONSTOP" => "Nonstop_Mutation",
        "FIVE_PRIME_UTR" => "5'UTR",
        "THREE_PRIME_UTR" => "3'UTR",
        "FIVE_PRIME_FLANK" => "5'Flank",
        "INTRON" => "Intron",
        "LINCRNA" => "RNA",
        _ => return None,
    })
}

/// The classifications `mafTransform` rewrites inside the other transcripts, in its order.
const ORDERED_GENCODE_CLASSIFICATIONS: [&str; 19] = [
    "IN_FRAME_DEL",
    "IN_FRAME_INS",
    "FRAME_SHIFT_INS",
    "FRAME_SHIFT_DEL",
    "MISSENSE",
    "NONSENSE",
    "SILENT",
    "SPLICE_SITE",
    "DE_NOVO_START_IN_FRAME",
    "DE_NOVO_START_OUT_FRAME",
    "START_CODON_SNP",
    "START_CODON_INS",
    "START_CODON_DEL",
    "NONSTOP",
    "FIVE_PRIME_UTR",
    "THREE_PRIME_UTR",
    "FIVE_PRIME_FLANK",
    "INTRON",
    "LINCRNA",
];

/// `StringUtils.replaceEachRepeatedly` with the classification lists.
fn replace_classifications(value: &str) -> String {
    // A pair whose replacement is null (a classification the map does not carry) is skipped.
    let pairs: Vec<(&str, &str)> = ORDERED_GENCODE_CLASSIFICATIONS
        .iter()
        .filter_map(|c| maf_variant_classification(c).map(|m| (*c, m)))
        .collect();
    let searches: Vec<&str> = pairs.iter().map(|(s, _)| *s).collect();
    let replacements: Vec<&str> = pairs.iter().map(|(_, r)| *r).collect();
    let mut text = value.to_string();
    // `replaceEach` finds, at each position, the earliest match over all search strings (the
    // lowest index among those with the smallest position), replaces it and moves past it; the
    // repeat runs that pass again until nothing changes, with a depth guard the lists never reach.
    for _ in 0..searches.len() + 1 {
        let mut out = String::new();
        let mut rest = text.as_str();
        let mut changed = false;
        loop {
            let mut best: Option<(usize, usize)> = None;
            for (index, search) in searches.iter().enumerate() {
                if let Some(position) = rest.find(search) {
                    if best.is_none_or(|(p, _)| position < p) {
                        best = Some((position, index));
                    }
                }
            }
            let Some((position, index)) = best else {
                out.push_str(rest);
                break;
            };
            out.push_str(&rest[..position]);
            out.push_str(replacements[index]);
            rest = &rest[position + searches[index].len()..];
            changed = true;
        }
        if !changed {
            break;
        }
        text = out;
    }
    text
}

/// `mafTransform`.
fn maf_transform(key: &str, value: &str, reference_version: &str) -> String {
    match key {
        "Variant_Classification" => {
            if let Some(mapped) = maf_variant_classification(value) {
                return mapped.to_string();
            }
        }
        "Chromosome" => {
            if value == "chrM" {
                return "MT".to_string();
            }
            let version = reference_version.to_lowercase();
            if value.to_lowercase().starts_with("chr") && (version == "hg19" || version == "b37") {
                let trimmed = &value[3..];
                let numbered = matches!(trimmed, "X" | "Y")
                    || trimmed
                        .parse::<u32>()
                        .is_ok_and(|n| (1..=22).contains(&n) && trimmed == n.to_string());
                if numbered {
                    return trimmed.to_string();
                }
            }
        }
        "Other_Transcripts" => return replace_classifications(value),
        _ => {}
    }
    value.to_string()
}

/// The MAF renderer's state: the default map, the aliases, and whether the header is out.
#[derive(Debug, Clone)]
pub struct MafRenderer {
    pub default_map: Vec<(String, String)>,
    /// The default annotations no column absorbed, which the header lists and no row carries.
    pub manual: Vec<(String, String)>,
    pub overrides: Vec<(String, String)>,
    pub excluded: Vec<String>,
    pub reference_version: String,
    pub header_written: bool,
}

fn put(map: &mut Vec<(String, String)>, key: &str, value: String) {
    match map.iter_mut().find(|(k, _)| k == key) {
        Some(slot) => slot.1 = value,
        None => map.push((key.to_string(), value)),
    }
}

impl MafRenderer {
    /// The constructor: the columns, the supported fields no alias takes, the default
    /// annotations folded in, and `Score` and `BAM_File` set to `NA`.
    pub fn new(
        supported_fields: &[String],
        defaults: &[(String, String)],
        overrides: &[(String, String)],
        excluded: &[String],
        reference_version: &str,
    ) -> MafRenderer {
        let mut default_map: Vec<(String, String)> = MAF_ALIASES
            .iter()
            .map(|(column, _)| {
                let value = if *column == "Strand" {
                    "+"
                } else {
                    "__UNKNOWN__"
                };
                (column.to_string(), value.to_string())
            })
            .collect();
        for field in supported_fields {
            if default_map.iter().any(|(k, _)| k == field) {
                continue;
            }
            let aliased = MAF_ALIASES
                .iter()
                .any(|(_, aliases)| aliases.contains(&field.as_str()));
            if !aliased {
                default_map.push((field.clone(), "__UNKNOWN__".to_string()));
            }
        }
        let mut manual: Vec<(String, String)> = Vec::new();
        for (key, value) in defaults {
            if default_map.iter().any(|(k, _)| k == key) {
                put(&mut default_map, key, value.clone());
                continue;
            }
            let mut aliased = false;
            for (column, aliases) in MAF_ALIASES {
                if aliases.contains(&key.as_str()) {
                    aliased = true;
                    put(&mut default_map, column, value.clone());
                }
            }
            if !aliased {
                manual.push((key.clone(), value.clone()));
            }
        }
        put(&mut default_map, "Score", "NA".to_string());
        put(&mut default_map, "BAM_File", "NA".to_string());
        MafRenderer {
            default_map,
            manual,
            overrides: overrides.to_vec(),
            excluded: excluded.to_vec(),
            reference_version: reference_version.to_string(),
            header_written: false,
        }
    }

    /// `createMafCompliantOutputMap`.
    pub fn row(&self, alternate: &Allele, funcotations: &[Funcotation]) -> Vec<(String, String)> {
        let mut output = self.default_map.clone();
        let mut extra: Vec<(String, String)> = Vec::new();
        for f in funcotations {
            if f.alt_allele() == *alternate {
                for name in f.field_names() {
                    let value = f.field(&name).unwrap_or_else(|| "__UNKNOWN__".to_string());
                    put(&mut extra, &name, value);
                }
            }
        }
        for (key, value) in &self.overrides {
            put(&mut extra, key, value.clone());
        }
        for (column, aliases) in MAF_ALIASES {
            for alias in *aliases {
                if let Some(index) = extra.iter().position(|(k, _)| k == alias) {
                    let (_, value) = extra.remove(index);
                    put(&mut output, column, value);
                    break;
                }
            }
        }
        for (key, value) in extra {
            put(&mut output, &key, value);
        }
        let mut transformed: Vec<(String, String)> = output
            .into_iter()
            .map(|(key, value)| {
                let value = maf_transform(&key, &value, &self.reference_version);
                (key, value)
            })
            .collect();
        if let Some(slot) = transformed
            .iter_mut()
            .find(|(k, _)| k == "Other_Transcripts")
        {
            slot.1 = slot.1.replace('/', "|");
        }
        transformed
            .into_iter()
            .filter(|(key, _)| !self.excluded.contains(key))
            .map(|(key, value)| (key, fu::sanitize_funcotation_field_for_maf(&value)))
            .collect()
    }

    /// `writeHeader`.
    pub fn header(
        &self,
        columns: &[String],
        input_lines: &[String],
        tool_lines: &[String],
        tool_version: &str,
        date: &str,
        data_sources: &str,
    ) -> String {
        let mut out = String::from("#version 2.4\n##\n");
        for line in input_lines.iter().chain(tool_lines) {
            out.push_str("## ");
            out.push_str(line);
            out.push('\n');
        }
        out.push_str(&format!(
            "##  Funcotator {tool_version} | Date {date} | {data_sources}\n"
        ));
        out.push_str(&columns.join("\t"));
        if self.manual.is_empty() {
            out.push('\n');
        } else {
            out.push('\t');
            out.push_str(
                &self
                    .manual
                    .iter()
                    .map(|(k, _)| k.as_str())
                    .collect::<Vec<_>>()
                    .join("\t"),
            );
            out.push('\n');
        }
        out
    }
}

// ================================================================================================
// The simple TSV and the gene list.
// ================================================================================================

/// The SEG output's columns and aliases, `simple_funcotator_seg_file.config`.
pub const SEG_ALIASES: &[(&str, &[&str])] = &[
    (
        "alt_allele",
        &[
            "Gencode_19_alt_allele",
            "Gencode_27_alt_allele",
            "Gencode_28_alt_allele",
        ],
    ),
    (
        "end_gene",
        &[
            "Gencode_19_end_gene",
            "Gencode_27_end_gene",
            "Gencode_28_end_gene",
        ],
    ),
    (
        "end",
        &[
            "END",
            "End",
            "End_Position",
            "end_position",
            "chromEnd",
            "segment_end",
            "End_position",
            "target_end",
            "stop",
            "Stop",
            "Position",
            "position",
            "pos",
            "POS",
            "segment_end",
        ],
    ),
    (
        "start_gene",
        &[
            "Gencode_19_start_gene",
            "Gencode_27_start_gene",
            "Gencode_28_start_gene",
        ],
    ),
    ("Segment_Mean", &["MEAN_LOG2_COPY_RATIO"]),
    (
        "genes",
        &["Gencode_19_genes", "Gencode_27_genes", "Gencode_28_genes"],
    ),
    ("Sample", &["sample", "sample_id"]),
    (
        "start",
        &[
            "START",
            "Start",
            "Start_Position",
            "start_position",
            "chromStart",
            "segment_start",
            "Start_position",
            "target_start",
            "Position",
            "position",
            "pos",
            "POS",
            "segment_start",
        ],
    ),
    (
        "chr",
        &[
            "CONTIG",
            "contig",
            "Chromosome",
            "chrom",
            "chromosome",
            "Chrom",
            "seqname",
            "seqnames",
            "CHROM",
            "target_contig",
            "segment_contig",
        ],
    ),
    ("build", &[]),
    ("Num_Probes", &["NUM_POINTS_COPY_RATIO"]),
    (
        "start_exon",
        &[
            "Gencode_19_start_exon",
            "Gencode_27_start_exon",
            "Gencode_28_start_exon",
        ],
    ),
    (
        "end_exon",
        &[
            "Gencode_19_end_exon",
            "Gencode_27_end_exon",
            "Gencode_28_end_exon",
        ],
    ),
    (
        "ref_allele",
        &[
            "Gencode_19_ref_allele",
            "Gencode_27_ref_allele",
            "Gencode_28_ref_allele",
        ],
    ),
    ("Segment_Call", &[]),
];

/// The gene list's columns and aliases, `gene_list_output.config`.
pub const GENE_LIST_ALIASES: &[(&str, &[&str])] = &[
    ("gene", &[]),
    ("exon", &[]),
    (
        "segment_contig",
        &[
            "CONTIG",
            "chr",
            "contig",
            "Chromosome",
            "chrom",
            "chromosome",
            "Chrom",
            "seqname",
            "seqnames",
            "CHROM",
            "target_contig",
        ],
    ),
    (
        "segment_start",
        &[
            "START",
            "start",
            "Start",
            "Start_Position",
            "start_position",
            "chromStart",
            "segment_start",
            "Start_position",
            "target_start",
            "Position",
            "position",
            "pos",
            "POS",
            "segment_start",
        ],
    ),
    (
        "segment_start_gene",
        &[
            "start_gene",
            "Gencode_19_start_gene",
            "Gencode_27_start_gene",
            "Gencode_28_start_gene",
        ],
    ),
    (
        "segment_start_exon",
        &[
            "start_exon",
            "Gencode_19_start_exon",
            "Gencode_27_start_exon",
            "Gencode_28_start_exon",
        ],
    ),
    (
        "segment_end",
        &[
            "END",
            "end",
            "End",
            "End_Position",
            "end_position",
            "chromEnd",
            "segment_end",
            "End_position",
            "target_end",
            "stop",
            "Stop",
            "Position",
            "position",
            "pos",
            "POS",
            "segment_end",
        ],
    ),
    (
        "segment_end_gene",
        &[
            "end_gene",
            "Gencode_19_end_gene",
            "Gencode_27_end_gene",
            "Gencode_28_end_gene",
        ],
    ),
    (
        "segment_end_exon",
        &[
            "end_exon",
            "Gencode_19_end_exon",
            "Gencode_27_end_exon",
            "Gencode_28_end_exon",
        ],
    ),
    (
        "segment_num_probes",
        &["Num_Probes", "NUM_POINTS_COPY_RATIO"],
    ),
    ("segment_mean", &["Segment_Mean", "MEAN_LOG2_COPY_RATIO"]),
    ("segment_call", &["Segment_Call"]),
    ("build", &[]),
    ("sample", &["Sample", "sample_id"]),
];

/// `SimpleNaturalComparator`: digit runs compared by value, everything else by character.
pub fn natural_compare(a: &str, b: &str) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i].is_ascii_digit() && b[j].is_ascii_digit() {
            let si = i;
            while i < a.len() && a[i].is_ascii_digit() {
                i += 1;
            }
            let sj = j;
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }
            let x: String = a[si..i]
                .iter()
                .collect::<String>()
                .trim_start_matches('0')
                .to_string();
            let y: String = b[sj..j]
                .iter()
                .collect::<String>()
                .trim_start_matches('0')
                .to_string();
            let order = x.len().cmp(&y.len()).then_with(|| x.cmp(&y));
            if order != Ordering::Equal {
                return order;
            }
        } else {
            if a[i] != b[j] {
                return a[i].cmp(&b[j]);
            }
            i += 1;
            j += 1;
        }
    }
    (a.len() - i).cmp(&(b.len() - j))
}

/// `SimpleTsvOutputRenderer`: the columns fixed at the first row, then a row per allele.
#[derive(Debug, Clone)]
pub struct SimpleTsv {
    pub aliases: &'static [(&'static str, &'static [&'static str])],
    pub defaults: Vec<(String, String)>,
    pub overrides: Vec<(String, String)>,
    pub excluded: Vec<String>,
    pub write_fields_not_in_config: bool,
    pub columns: Option<Vec<(String, String)>>,
    pub text: String,
}

impl SimpleTsv {
    pub fn new(
        aliases: &'static [(&'static str, &'static [&'static str])],
        defaults: &[(String, String)],
        overrides: &[(String, String)],
        excluded: &[String],
        write_fields_not_in_config: bool,
    ) -> SimpleTsv {
        SimpleTsv {
            aliases,
            defaults: defaults.to_vec(),
            overrides: overrides.to_vec(),
            excluded: excluded.to_vec(),
            write_fields_not_in_config,
            columns: None,
            text: String::new(),
        }
    }

    /// `createColumnNameToFieldNameMap`.
    fn column_map(
        &self,
        map: &FuncotationMap,
        transcript: &str,
        write_leftovers: bool,
    ) -> Vec<(String, String)> {
        let names = map.field_names(transcript);
        let mut result: Vec<(String, String)> = self
            .aliases
            .iter()
            .map(|(column, aliases)| {
                let found = std::iter::once(*column)
                    .chain(aliases.iter().copied())
                    .find(|candidate| names.iter().any(|n| n == candidate))
                    .unwrap_or("");
                (column.to_string(), found.to_string())
            })
            .collect();
        let used: Vec<String> = result.iter().map(|(_, v)| v.clone()).collect();
        let mut leftovers: Vec<String> = Vec::new();
        if write_leftovers {
            let mut extra: Vec<String> = names
                .iter()
                .filter(|n| !used.contains(n))
                .cloned()
                .collect();
            extra.sort();
            leftovers.extend(extra);
        }
        for source in [&self.defaults, &self.overrides] {
            let mut extra: Vec<String> = source
                .iter()
                .map(|(k, _)| k.clone())
                .filter(|k| !leftovers.contains(k))
                .collect();
            extra.sort();
            extra.dedup();
            leftovers.extend(extra);
        }
        leftovers.sort_by(|a, b| natural_compare(a, b));
        for column in leftovers {
            if !result.iter().any(|(k, _)| *k == column) {
                result.push((column.clone(), column));
            }
        }
        result.retain(|(k, _)| !self.excluded.contains(k));
        result
    }

    fn initialize(&mut self, columns: Vec<(String, String)>) -> Result<(), EngineError> {
        if columns.is_empty() {
            return Err(EngineError::User {
                class: "java.lang.IllegalArgumentException",
                message: "TSV output renderer has been configured to produce a blank file.  This is usually a user error.  Please check excluded columns.".to_string(),
            });
        }
        self.text.push_str(
            &columns
                .iter()
                .map(|(k, _)| k.as_str())
                .collect::<Vec<_>>()
                .join("\t"),
        );
        self.text.push('\n');
        self.columns = Some(columns);
        Ok(())
    }

    /// `write`: the locatable fields added to every transcript and allele, then the rows.
    pub fn write(
        &mut self,
        contig: &str,
        start: i32,
        end: i32,
        map: &mut FuncotationMap,
    ) -> Result<(), EngineError> {
        if !map.same_fields_everywhere() {
            return Err(EngineError::User {
                class: "org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException",
                message: "The funcotation map cannot be written by this simple output renderer.  The fields in the funcotation map do not match across transcript-allele combinations.  This is almost certainly an issue for the GATK development team.".to_string(),
            });
        }
        for transcript in map.transcripts() {
            for allele in map.alleles(&transcript) {
                let locatable = table(
                    &[
                        ("CONTIG".to_string(), contig.to_string()),
                        ("START".to_string(), start.to_string()),
                        ("END".to_string(), end.to_string()),
                    ],
                    &allele,
                    "SIMPLE_TSV_OUTPUT_RENDERER",
                );
                map.add(&transcript, &[locatable]);
            }
        }
        if self.columns.is_none() {
            let first = map.transcripts().first().cloned().unwrap_or_default();
            let columns = self.column_map(map, &first, self.write_fields_not_in_config);
            self.initialize(columns)?;
        }
        let columns = self.columns.clone().unwrap_or_default();
        for transcript in map.transcripts() {
            for allele in map.alleles(&transcript) {
                let mut values = Vec::new();
                for (column, field) in &columns {
                    if self.excluded.contains(column) {
                        continue;
                    }
                    let value = map.field_value(&transcript, field, &allele)?;
                    let value = match self.overrides.iter().find(|(k, _)| k == column) {
                        Some((_, v)) => v.clone(),
                        None => match value {
                            Some(v) => v,
                            None => self
                                .defaults
                                .iter()
                                .find(|(k, _)| k == column)
                                .map(|(_, v)| v.clone())
                                .unwrap_or_else(|| "__UNKNOWN__".to_string()),
                        },
                    };
                    values.push(value);
                }
                self.text.push_str(&values.join("\t"));
                self.text.push('\n');
            }
        }
        Ok(())
    }

    /// `close`: a file with no row still gets the columns it would have had.
    pub fn close(&mut self) -> Result<String, EngineError> {
        if self.columns.is_none() {
            let columns = self.column_map(&FuncotationMap::default(), "DUMMY", false);
            self.initialize(columns)?;
        }
        Ok(self.text.clone())
    }
}

/// A gene and one of its exons, or the empty exon for the whole gene.
pub type GeneExon = (String, String);

/// A segment's locus and the funcotation map it was written with.
pub type SegmentRecord = (String, i32, i32, FuncotationMap);

/// `GeneListOutputRenderer`: the genes and exons each segment covers, flushed at the end.
#[derive(Debug, Clone)]
pub struct GeneList {
    pub tsv: SimpleTsv,
    pub min_bases: i32,
    /// `(gene, exon)` to the segment and its map, natural-ordered on the pair.
    pub entries: Vec<(GeneExon, SegmentRecord)>,
}

impl GeneList {
    pub fn new(
        defaults: &[(String, String)],
        overrides: &[(String, String)],
        excluded: &[String],
        min_bases: i32,
    ) -> GeneList {
        GeneList {
            tsv: SimpleTsv::new(GENE_LIST_ALIASES, defaults, overrides, excluded, false),
            min_bases,
            entries: Vec::new(),
        }
    }

    fn put(&mut self, key: GeneExon, value: SegmentRecord) {
        match self.entries.binary_search_by(|(k, _)| {
            natural_compare(&k.0, &key.0).then_with(|| natural_compare(&k.1, &key.1))
        }) {
            Ok(index) => self.entries[index].1 = value,
            Err(index) => self.entries.insert(index, (key, value)),
        }
    }

    /// `write`, which validates the segment and records its genes.
    pub fn write(
        &mut self,
        variant: &Variant,
        vc_string: &str,
        map: &FuncotationMap,
    ) -> Result<(), EngineError> {
        let bad = |message: String| EngineError::User {
            class: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
            message: format!("Bad input: {message}"),
        };
        const IS_SEG: &str = "Is the input a file of segment variant contexts?";
        if !is_segment_variant(variant, self.min_bases) {
            return Err(bad(format!(
                "{IS_SEG}  Variant context does not represent a copy number segment: {vc_string}"
            )));
        }
        let transcripts = map.transcripts();
        if transcripts.len() != 1 {
            return Err(bad(format!(
                "{IS_SEG}  Need exactly one transcript ID: {}",
                transcripts.join(",")
            )));
        }
        if transcripts[0] != NO_TRANSCRIPT_KEY {
            return Err(bad(format!(
                "{IS_SEG}  Invalid transcript ID seen  (must be no transcript available dummy ID): {}",
                transcripts[0]
            )));
        }
        let alleles = map.alleles(NO_TRANSCRIPT_KEY);
        if alleles.len() != 1 {
            return Err(bad(format!(
                "{IS_SEG}  Only one alternate allele per variant context is accepted."
            )));
        }
        let allele = &alleles[0];
        let names = map.field_names(NO_TRANSCRIPT_KEY);
        let field = |suffix: &str| -> Result<String, EngineError> {
            let matches: Vec<&String> = names
                .iter()
                .filter(|n| {
                    n.strip_prefix("Gencode_")
                        .and_then(|rest| rest.split_once('_'))
                        .is_some_and(|(number, tail)| {
                            !number.is_empty()
                                && number.bytes().all(|b| b.is_ascii_digit())
                                && format!("_{tail}") == suffix
                        })
                })
                .collect();
            if matches.len() != 1 {
                return Err(bad(format!(
                    "Could not find exactly one gencode field match in the funcotation map: {}",
                    matches
                        .iter()
                        .map(|s| s.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                )));
            }
            Ok(map
                .field_value(NO_TRANSCRIPT_KEY, matches[0], allele)?
                .unwrap_or_default())
        };
        let genes_value = field("_genes")?;
        let start_gene = field("_start_gene")?;
        let start_exon = field("_start_exon")?;
        let end_gene = field("_end_gene")?;
        let end_exon = field("_end_exon")?;
        let record = (
            variant.contig.clone(),
            variant.start,
            variant.end,
            map.clone(),
        );
        for (gene, exon) in [
            (start_gene.clone(), start_exon),
            (end_gene.clone(), end_exon),
        ] {
            if !gene.is_empty() {
                self.put((gene, exon), record.clone());
            }
        }
        for gene in genes_value.split(',').filter(|g| !g.is_empty()) {
            if gene != start_gene && gene != end_gene {
                self.put((gene.to_string(), String::new()), record.clone());
            }
        }
        Ok(())
    }

    /// `close`: every gene and exon written, then the file.
    pub fn close(&mut self) -> Result<String, EngineError> {
        let entries = std::mem::take(&mut self.entries);
        for ((gene, exon), (contig, start, end, mut map)) in entries {
            let allele = map.get(NO_TRANSCRIPT_KEY)[0].alt_allele();
            let f = table(
                &[("gene".to_string(), gene), ("exon".to_string(), exon)],
                &allele,
                "GENE_LIST_OUTPUT_RENDERER",
            );
            map.add(NO_TRANSCRIPT_KEY, &[f]);
            self.tsv.write(&contig, start, end, &mut map)?;
        }
        self.tsv.close()
    }
}
