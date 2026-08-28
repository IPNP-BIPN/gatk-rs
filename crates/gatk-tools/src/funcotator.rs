//! `Funcotator`: the strings a GENCODE funcotation is built out of.
//!
//! Reading the data sources, walking the transcripts and rendering the VCF are not ported. What is
//! ported is the arithmetic and the formatting the annotation is made of: the genetic code, the
//! codon alignment, and the four change strings a coding variant carries.
//!
//! Ported from `org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils`,
//! `org.broadinstitute.hellbender.tools.funcotator.AminoAcid` and
//! `org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory`.

/// `AminoAcid.CODON_LENGTH`.
pub const CODON_LENGTH: i32 = 3;
/// `VcfOutputRenderer.FIELD_DELIMITER`, which is also the header-listed one.
pub const FIELD_DELIMITER: char = '|';
/// `VcfOutputRenderer.ALL_TRANSCRIPT_DELIMITER`.
pub const ALL_TRANSCRIPT_DELIMITER: char = '#';
/// `VcfOutputRenderer.DESCRIPTION_PREAMBLE_DELIMITER`.
pub const DESCRIPTION_PREAMBLE_DELIMITER: &str = ": ";
/// `GencodeFuncotation`'s gene name where no transcript covers the variant.
pub const UNKNOWN_GENE: &str = "Unknown";
/// And its transcript name in the same place.
pub const NO_TRANSCRIPT: &str = "no_transcript";

/// `AminoAcid`, whose codon lists carry the IUPAC ambiguity codes beside the plain ones.
///
/// The order is the enum's own, which is alphabetical by name with the undecodable one last.
pub const AMINO_ACIDS: &[(&str, &str, &str, &[&str])] = &[
    ("Alanine", "Ala", "A", &["GCA", "GCC", "GCG", "GCT", "GCN"]),
    (
        "Arganine",
        "Arg",
        "R",
        &[
            "AGA", "AGG", "CGA", "CGC", "CGG", "CGT", "CGN", "AGR", "CGY", "MGR",
        ],
    ),
    ("Asparagine", "Asn", "N", &["AAC", "AAT", "AAY"]),
    ("Aspartic acid", "Asp", "D", &["GAT", "GAC", "GAY"]),
    ("Cysteine", "Cys", "C", &["TGC", "TGT", "TGY"]),
    ("Glutamic acid", "Glu", "E", &["GAA", "GAG", "GAR"]),
    ("Glutamine", "Gln", "Q", &["CAA", "CAG", "CAR"]),
    ("Glycine", "Gly", "G", &["GGA", "GGC", "GGG", "GGT", "GGN"]),
    ("Histidine", "His", "H", &["CAC", "CAT", "CAY"]),
    ("Isoleucine", "Ile", "I", &["ATA", "ATC", "ATT", "ATH"]),
    (
        "Leucine",
        "Leu",
        "L",
        &[
            "CTA", "CTC", "CTG", "CTT", "TTA", "TTG", "CTN", "CTY", "TTR", "YTR",
        ],
    ),
    ("Lysine", "Lys", "K", &["AAA", "AAG", "AAR"]),
    ("Methionine", "Met", "M", &["ATG"]),
    ("Phenylalanine", "Phe", "F", &["TTC", "TTT", "TTY"]),
    ("Proline", "Pro", "P", &["CCA", "CCC", "CCG", "CCT", "CCN"]),
    (
        "Serine",
        "Ser",
        "S",
        &["AGC", "AGT", "TCA", "TCC", "TCG", "TCT", "TCN", "AGY"],
    ),
    (
        "Stop codon",
        "Stop",
        "*",
        &["TAA", "TAG", "TGA", "TRA", "TAR"],
    ),
    (
        "Threonine",
        "Thr",
        "T",
        &["ACA", "ACC", "ACG", "ACT", "ACN"],
    ),
    ("Tryptophan", "Trp", "W", &["TGG"]),
    ("Tyrosine", "Tyr", "Y", &["TAC", "TAT", "TAY"]),
    ("Valine", "Val", "V", &["GTA", "GTC", "GTG", "GTT", "GTN"]),
    ("Undecodable Amino Acid", "UNDECODABLE", "?", &[]),
];

/// `getEukaryoticAminoAcidByCodon`, which UPPER-CASES the codon before it looks it up and answers
/// nothing at all for a codon no amino acid claims.
pub fn eukaryotic_amino_acid_by_codon(codon: &str) -> Option<&'static str> {
    let codon = codon.to_uppercase();
    AMINO_ACIDS
        .iter()
        .find(|(_, _, _, codons)| codons.contains(&codon.as_str()))
        .map(|(_, _, letter, _)| *letter)
}

/// `isPositionInFrame`, which is a position 1-based on the coding sequence.
pub fn is_position_in_frame(position: i32) -> bool {
    (position - 1) % CODON_LENGTH == 0
}

/// `getAlignedPosition`: the next LOWEST position a codon holding this one would start at.
///
/// A position at or before the start of the sequence is folded back the other way, which is what
/// an upstream UTR or flank needs.
pub fn aligned_position(position: i32) -> i32 {
    if position > 0 {
        position - ((position - 1) % CODON_LENGTH)
    } else {
        let adjusted = 1 - position;
        -(adjusted - ((adjusted - 1) % CODON_LENGTH) + 1)
    }
}

/// `getAlignedEndPosition`: the end of the codon the given end position falls in.
pub fn aligned_end_position(allele_end_position: i32) -> i32 {
    ((f64::from(allele_end_position) / f64::from(CODON_LENGTH)).ceil() as i32) * CODON_LENGTH
}

/// `getProteinChangePosition`, the amino acid a coding position sits in.
pub fn protein_change_position(aligned_coding_sequence_allele_start: i32) -> i32 {
    ((aligned_coding_sequence_allele_start - 1) / CODON_LENGTH) + 1
}

/// `getProteinChangeEndPosition`, which adds the ALLELE's own length in amino acids, less one.
pub fn protein_change_end_position(
    protein_change_start_position: i32,
    aligned_alternate_allele_length: i32,
) -> i32 {
    protein_change_start_position + protein_change_position(aligned_alternate_allele_length) - 1
}

/// `getCodonChangeStringForOnp`: the aligned codons with the CHANGED bases in upper case and the
/// unchanged ones in lower.
///
/// A one-base codon window prints one position and a wider one prints the range, which is why a
/// SNP in the middle of a codon still reads `c.(61-63)Gcg>Acg`.
pub fn codon_change_string_for_onp(
    aligned_reference_allele: &str,
    aligned_alternate_allele: &str,
    aligned_coding_sequence_allele_start: i32,
    aligned_reference_allele_stop: i32,
) -> String {
    let mut reference = String::new();
    let mut alternate = String::new();
    for (r, a) in aligned_reference_allele
        .chars()
        .zip(aligned_alternate_allele.chars())
    {
        if r != a {
            reference.extend(r.to_uppercase());
            alternate.extend(a.to_uppercase());
        } else {
            let lower: String = r.to_lowercase().collect();
            reference.push_str(&lower);
            alternate.push_str(&lower);
        }
    }
    if aligned_coding_sequence_allele_start == aligned_reference_allele_stop {
        format!("c.({aligned_coding_sequence_allele_start}){reference}>{alternate}")
    } else {
        format!(
            "c.({aligned_coding_sequence_allele_start}-{aligned_reference_allele_stop}){reference}>{alternate}"
        )
    }
}

/// `renderProteinChangeString`'s last branch, which is the one an ONP takes.
///
/// A change confined to one amino acid is written `p.<ref><position><alt>`, with the position in
/// the MIDDLE. A wider one is written `p.<start>_<end><ref>><alt>`, with the positions first.
pub fn protein_change_string_for_onp(
    reference_amino_acids: &str,
    alternate_amino_acids: &str,
    start_position: i32,
    end_position: i32,
) -> String {
    if start_position != end_position {
        format!("p.{start_position}_{end_position}{reference_amino_acids}>{alternate_amino_acids}")
    } else {
        format!("p.{reference_amino_acids}{start_position}{alternate_amino_acids}")
    }
}

/// `getCodingSequenceChangeString`'s XNP branch.
///
/// A single base prints one position; a longer ONP prints the range. The bases are the
/// strand-corrected ones, so a negative-strand transcript has already reverse-complemented them.
pub fn coding_sequence_change_string_for_xnp(
    coding_sequence_allele_start: i32,
    reference_allele: &str,
    alternate_allele: &str,
) -> String {
    if alternate_allele.len() > 1 {
        let stop = coding_sequence_allele_start + reference_allele.len() as i32 - 1;
        format!("c.{coding_sequence_allele_start}_{stop}{reference_allele}>{alternate_allele}")
    } else {
        format!("c.{coding_sequence_allele_start}{reference_allele}>{alternate_allele}")
    }
}

/// `getClosestExonIndex`: the exon the fewest bases from the variant's start, by the SMALLER of
/// the two distances to its ends.
pub fn closest_exon_index(variant_start: i32, exons: &[(i32, i32)]) -> Option<usize> {
    let mut best: Option<(usize, i32)> = None;
    for (index, (start, end)) in exons.iter().enumerate() {
        let distance = (variant_start - start)
            .abs()
            .min((variant_start - end).abs());
        if best.is_none_or(|(_, closest)| distance < closest) {
            best = Some((index, distance));
        }
    }
    best.map(|(index, _)| index)
}

/// `createIntronicCDnaString`, whose offset is signed from the NEARER end of the closest exon.
///
/// The exon is named by its one-based index, so `c.e1+300T>A` is three hundred bases past the end
/// of the first exon. A transcript with no exon at all answers the literal `NA`.
pub fn intronic_cdna_string(
    variant_start: i32,
    exons: &[(i32, i32)],
    reference_allele: &str,
    alternate_allele: &str,
) -> String {
    let Some(index) = closest_exon_index(variant_start, exons) else {
        return "NA".to_string();
    };
    let (start, end) = exons[index];
    let start_difference = variant_start - start;
    let end_difference = variant_start - end;
    let offset = if start_difference.abs() <= end_difference.abs() {
        start_difference
    } else {
        end_difference
    };
    let sign = if offset < 0 { '-' } else { '+' };
    let magnitude = offset.abs();
    format!(
        "c.e{}{sign}{magnitude}{reference_allele}>{alternate_allele}",
        index + 1
    )
}

/// `getGenomeChangeString`'s SNP and ONP branches.
pub fn genome_change_string_for_xnp(
    contig: &str,
    start: i32,
    reference_allele: &str,
    alternate_allele: &str,
) -> String {
    if reference_allele.len() == 1 {
        format!("g.{contig}:{start}{reference_allele}>{alternate_allele}")
    } else {
        let stop = start + reference_allele.len() as i32 - 1;
        format!("g.{contig}:{start}_{stop}{reference_allele}>{alternate_allele}")
    }
}

/// `getTranscriptIdWithoutVersionNumber`, which strips a trailing dot and digits and nothing else.
pub fn transcript_id_without_version_number(transcript_id: &str) -> String {
    match transcript_id.rfind('.') {
        Some(dot)
            if dot + 1 < transcript_id.len()
                && transcript_id[dot + 1..].bytes().all(|b| b.is_ascii_digit()) =>
        {
            transcript_id[..dot].to_string()
        }
        _ => transcript_id.to_string(),
    }
}

/// The bad letters `sanitizeFuncotationFieldForVcf` encodes, in the order it lists them.
const BAD_LETTERS: &[&str] = &[",", ";", "=", "\t", "|", " ", "\n", "#"];

/// `sanitizeFuncotationFieldForVcf`: each of eight characters becomes `_%<HEX>_`, upper case, with
/// a leading zero under sixteen.
pub fn sanitize_funcotation_field_for_vcf(field: &str) -> String {
    let mut sanitized = field.to_string();
    for letter in BAD_LETTERS {
        let byte = letter.as_bytes()[0];
        let hex = if byte < 16 {
            format!("_%0{:X}_", byte)
        } else {
            format!("_%{:X}_", byte)
        };
        sanitized = sanitized.replace(letter, &hex);
    }
    sanitized
}

/// `sanitizeFuncotationFieldForMaf`, which encodes only the tab and the newline.
pub fn sanitize_funcotation_field_for_maf(field: &str) -> String {
    field.replace('\t', "_%09_").replace('\n', "_%0A_")
}

/// `renderSanitizedFuncotationForVcf`: the included fields' values, sanitized, joined by `|`.
///
/// The order is the FUNCOTATION's own field order and not the included list's, and an empty
/// included list renders nothing at all.
pub fn render_sanitized_funcotation_for_vcf(
    fields: &[(String, String)],
    included: &[String],
) -> String {
    if included.is_empty() {
        return String::new();
    }
    fields
        .iter()
        .filter(|(name, _)| included.contains(name))
        .map(|(_, value)| sanitize_funcotation_field_for_vcf(value))
        .collect::<Vec<_>>()
        .join(&FIELD_DELIMITER.to_string())
}

/// `extractFuncotatorKeysFromHeaderDescription`: everything past the first `": "`, split on `|`.
///
/// The split is on the WHOLE separator and preserves empty tokens, so a description ending on the
/// delimiter yields a trailing empty key, which is what the reference's own header does.
pub fn extract_funcotator_keys_from_header_description(description: &str) -> Vec<String> {
    let Some((_, rest)) = description.split_once(DESCRIPTION_PREAMBLE_DELIMITER) else {
        return Vec::new();
    };
    rest.split(FIELD_DELIMITER).map(str::to_string).collect()
}
