//! `Funcotator`: the strings a GENCODE funcotation is built out of.
//!
//! Reading the data sources, walking the transcripts and rendering the VCF are not ported. What is
//! ported is the arithmetic and the formatting the annotation is made of: the genetic code, the
//! codon alignment, and the four change strings a coding variant carries.
//!
//! Ported from `org.broadinstitute.hellbender.tools.funcotator.FuncotatorUtils`,
//! `org.broadinstitute.hellbender.tools.funcotator.AminoAcid` and
//! `org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory`.
//!
//! # Which GATK
//!
//! GATK 4.7.0.0 for everything here except [`transcript_id_without_version_number`], which keeps
//! the 4.6.2.0 body (see its note): `FuncotatorUtils.java` changed between the two releases in
//! that one method. `AminoAcid.java` is unchanged since 4.6.2.0.
//!
//! # The reference's argument checks are not ported
//!
//! Most of these Java methods open with `ParamUtils.isPositive(...)` or `Utils.nonNull(...)`,
//! which throw `IllegalArgumentException` on a non-positive position. The functions below take
//! any value and compute with it; their callers pass positions that are already one-based.
//!
//! # Reading this file without reading Rust
//!
//! - `&str` is a borrowed, read-only string; `String` is an owned one. `format!` is Java's
//!   `String.format`, with `{name}` substituting the variable `name`.
//! - `Option<T>` is "a value or nothing": `Some(x)` or `None`, Java's nullable return with the
//!   null check enforced by the compiler.
//! - `.iter().map(..).collect()` and friends are Java's stream, map, collect; a closure
//!   `|x| ...` is a Java lambda `x -> ...`.
//! - `i32` is Java's `int`; `usize` is an unsigned index type Java lacks. `x as i32` converts
//!   between them, truncating (wrapping) a value that does not fit.
//! - A function's last expression is its return value; `return` is written only for early exits.

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
/// Each row is `(name, three-letter code, one-letter code, codons)`, copied from the enum's
/// constructor arguments, the misspelt `"Arganine"` included. The order is the enum's own, which
/// is alphabetical by name with the undecodable one last.
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
///
/// The Java answers from a `HashMap` built once from the enum (`tableByCodon`); a linear search of
/// [`AMINO_ACIDS`] finds the same answer, because no codon is listed under two amino acids.
pub fn eukaryotic_amino_acid_by_codon(codon: &str) -> Option<&'static str> {
    let codon = codon.to_uppercase();
    AMINO_ACIDS
        .iter()
        .find(|(_, _, _, codons)| codons.contains(&codon.as_str()))
        .map(|(_, _, letter, _)| *letter)
}

/// `isPositionInFrame`, which is a position 1-based on the coding sequence.
///
/// Rust's `%` truncates toward zero like Java's, so a negative position gives the same remainder
/// in both; the reference refuses one before it gets here.
pub fn is_position_in_frame(position: i32) -> bool {
    (position - 1) % CODON_LENGTH == 0
}

/// `getAlignedPosition`: the next LOWEST position a codon holding this one would start at.
///
/// A position at or before the start of the sequence is folded back the other way, which is what
/// an upstream UTR or flank needs. Both branches are the Java's expressions unchanged:
///
/// ```java
/// final int adjustedPos = 1 - position;
/// return -(adjustedPos - ((adjustedPos - 1) % AminoAcid.CODON_LENGTH) + 1);
/// ```
pub fn aligned_position(position: i32) -> i32 {
    if position > 0 {
        position - ((position - 1) % CODON_LENGTH)
    } else {
        let adjusted = 1 - position;
        -(adjusted - ((adjusted - 1) % CODON_LENGTH) + 1)
    }
}

/// `getAlignedEndPosition`: the end of the codon the given end position falls in.
///
/// ```java
/// (int)(Math.ceil(alleleEndPosition / ((double)AminoAcid.CODON_LENGTH)) * AminoAcid.CODON_LENGTH)
/// ```
///
/// The division is in `double`, as in the Java, and `as i32` is the `(int)` cast. One difference
/// in order: the Java multiplies the `double` and then casts, this casts and then multiplies in
/// `i32`. They agree until the product passes `Integer.MAX_VALUE` (only at an end position of
/// 2147483647), where the Java's `(int)` saturates and this multiplication overflows.
pub fn aligned_end_position(allele_end_position: i32) -> i32 {
    ((f64::from(allele_end_position) / f64::from(CODON_LENGTH)).ceil() as i32) * CODON_LENGTH
}

/// `getProteinChangePosition`, the amino acid a coding position sits in.
///
/// Integer division, as the Java's `int` division: it truncates toward zero.
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
///
/// The Java walks `alignedRefAllele.length()` characters and reads the alternate at the same
/// index, so an alternate shorter than the reference throws there. `zip` stops at the shorter of
/// the two instead; callers pass two alleles of one length.
pub fn codon_change_string_for_onp(
    aligned_reference_allele: &str,
    aligned_alternate_allele: &str,
    aligned_coding_sequence_allele_start: i32,
    aligned_reference_allele_stop: i32,
) -> String {
    let mut reference = String::new();
    let mut alternate = String::new();
    // `zip` walks the two strings in step, pairing their characters, Java's loop over one index.
    // `to_uppercase()` yields an iterator of characters (one, for every base letter), which
    // `extend` appends; `Character.toUpperCase` is the Java.
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
///
/// Both are `String.format(PROTEIN_CHANGE_FORMAT_STRING, ...)` with `"p.%s%s%s%s"`: the wide one
/// passes `protChangeStartPos + "_" + protChangeEndPos, refAaSeq , ">", altAaSeq`, the narrow one
/// `"", refAaSeq , protChangeStartPos, altAaSeq`.
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
///
/// The test is on the ALTERNATE's length while the range is computed from the REFERENCE's, as in
/// the Java (`if (altAllele.length() > 1)`, then `refAllele.length() - 1`); for an XNP the two
/// are equal. `.len()` is a byte count, which for base letters is the character count.
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
///
/// The Java loops over `Arrays.asList(start, end)` per exon and replaces the best on a strictly
/// smaller distance; taking the per-exon minimum first and comparing once gives the same winner,
/// and the strict `<` keeps the EARLIER exon on a tie, as the Java does. `None` here is the Java's
/// `-1`.
pub fn closest_exon_index(variant_start: i32, exons: &[(i32, i32)]) -> Option<usize> {
    let mut best: Option<(usize, i32)> = None;
    for (index, (start, end)) in exons.iter().enumerate() {
        let distance = (variant_start - start)
            .abs()
            .min((variant_start - end).abs());
        // `is_none_or` is true when nothing is held yet, or when the closure is true of what is
        // held: Java's `best == null || distance < best.distance`.
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
///
/// On a tie between the two distances the START wins (`Math.abs(startDiff) <= Math.abs(endDiff)`
/// in the Java), and an offset of zero prints as `+0`.
pub fn intronic_cdna_string(
    variant_start: i32,
    exons: &[(i32, i32)],
    reference_allele: &str,
    alternate_allele: &str,
) -> String {
    // `let Some(x) = .. else { .. }` takes the value out of an `Option` or runs the `else` block,
    // which must leave the function. Here it is the Java's `exonIndex != -1` test.
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
///
/// The insertion and deletion branches before them (decided on the two alleles' lengths) are not
/// ported. The SNP branch is chosen by `variant.getReference().length() == 1`, which is the test
/// here.
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

/// `getTranscriptIdWithoutVersionNumber` AS GATK 4.6.2.0 WROTE IT: strip one trailing dot and
/// digits and nothing else, which is `transcriptId.replaceAll("\\.\\d+$", "")`.
///
/// GATK 4.7.0.0 cuts at the FIRST dot instead, whatever follows it:
///
/// ```java
/// final int end = transcriptId.indexOf('.');
/// if(end != -1) {
///     return transcriptId.substring(0, end);
/// }
/// ```
///
/// So the two disagree on an id with a second dot or a non-digit suffix: GENCODE's
/// `ENST00000381192.10_PAR_Y` is `ENST00000381192` to 4.7.0.0 and is left whole here, and the
/// suite's own `ENST1.2.3` is `ENST1` to 4.7.0.0 and `ENST1.2` here.
pub fn transcript_id_without_version_number(transcript_id: &str) -> String {
    // `rfind` is `lastIndexOf`. The guard (`if ...` after the pattern) requires at least one
    // character after the dot and every one of them an ASCII digit, which is `\.\d+$`.
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
///
/// The Java makes one `StringUtils.replaceEach` call over all eight; this replaces them one after
/// another. The result is the same because no replacement text (`_`, `%`, hex digits) contains
/// any of the eight letters, so a later pass never rewrites an earlier pass's output.
pub fn sanitize_funcotation_field_for_vcf(field: &str) -> String {
    let mut sanitized = field.to_string();
    for letter in BAD_LETTERS {
        // The Java's `s.getBytes(StandardCharsets.US_ASCII)[0]`; every letter here is one ASCII
        // byte. `{:X}` is upper-case hex, Java's `Integer.toHexString(..).toUpperCase()`.
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
///
/// The Java is one `StringUtils.replaceEach` with `new String[]{"\t", "\n"}` and
/// `new String[]{"_%09_", "_%0A_"}`; neither replacement contains a tab or a newline, so two
/// `replace` calls give the same text.
pub fn sanitize_funcotation_field_for_maf(field: &str) -> String {
    field.replace('\t', "_%09_").replace('\n', "_%0A_")
}

/// `renderSanitizedFuncotationForVcf`: the included fields' values, sanitized, joined by `|`.
///
/// The order is the FUNCOTATION's own field order and not the included list's, and an empty
/// included list renders nothing at all.
///
/// ```java
/// return funcotation.getFieldNames().stream()
///         .filter(includedFields::contains)
/// ```
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

/// `extractFuncotatorKeysFromHeaderDescription`: the text after the first `": "`, split on `|`.
///
/// The split on `|` preserves empty tokens, so a description ending on the delimiter yields a
/// trailing empty key, which is what the reference's own header does.
///
/// The Java splits the WHOLE description on `VcfOutputRenderer.DESCRIPTION_PREAMBLE_DELIMITER`
/// with `StringUtils.splitByWholeSeparatorPreserveAllTokens` and then splits only
/// `descriptionSplit[1]`, the second piece, on `VcfOutputRenderer.FIELD_DELIMITER`.
///
/// So the two differ on a description holding `": "` twice (`a: B|C: D` gives `B`, `C` there and
/// `B`, `C: D` here), and on one holding no `": "` at all, which the Java's `descriptionSplit[1]`
/// throws on and this answers with no keys.
pub fn extract_funcotator_keys_from_header_description(description: &str) -> Vec<String> {
    // `split_once` splits at the first occurrence only and answers `None` when there is none.
    let Some((_, rest)) = description.split_once(DESCRIPTION_PREAMBLE_DELIMITER) else {
        return Vec::new();
    };
    rest.split(FIELD_DELIMITER).map(str::to_string).collect()
}
