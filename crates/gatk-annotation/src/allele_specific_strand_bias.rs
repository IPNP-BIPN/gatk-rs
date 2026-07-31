//! `AS_FisherStrand`, `AS_StrandOddsRatio` and the `AS_StrandBiasTest` / `StrandBiasUtils`
//! machinery under them, ported from GATK 4.6.2.0.
//!
//! `AS_FS`, `AS_SOR` and the `AS_SB_TABLE` they are both derived from: one forward/reverse pair per
//! allele, summed across samples, carried through a gVCF as a pipe-separated list.
//!
//! # The raw string has an entry for the **reference**, unlike the rank sums'
//!
//! ```java
//! for (final Allele a : vcAlleles) {
//!     if (!annotationString.isEmpty()) { annotationString += ALLELE_SPECIFIC_RAW_DELIM; }
//! ```
//!
//! `AS_SB_TABLE` writes every allele including the reference and puts the delimiter **between**
//! entries, so it does not start with one. The rank sums skip the reference's value but keep its
//! slot, so theirs does. Two allele-specific families, two incompatible conventions for the same
//! delimiter, and the parser here refuses a count mismatch outright:
//!
//! ```java
//! throw new IllegalStateException("Number of alleles and number of allele-specific entries do not match.")
//! ```
//!
//! # A sample contributes nothing at all unless it has **more than two** informative reads
//!
//! ```java
//! return readCount > minCount;   // minCount == 2
//! ```
//!
//! Strictly greater, and the threshold is on the sample's whole table rather than on any one
//! allele. A sample with two informative reads is dropped entirely, so its strand counts do not
//! appear in `AS_SB_TABLE` even for the allele that had both of them.
//!
//! # `AS_SOR` computes a value for the reference allele and then never prints it
//!
//! ```java
//! for (final Allele a : perAlleleData.keySet()) {          // AS_StrandOddsRatio: no filter
//! for (final Allele a : perAlleleData.keySet()) {
//!     if(!a.equals(combinedData.getRefAllele(),true)) {    // AS_FisherStrand: filtered
//! ```
//!
//! `AS_StrandOddsRatio.calculateReducedData` iterates every allele, so it takes a symmetric odds
//! ratio of the reference against itself. `makeReducedAnnotationString` then walks only the
//! alternates, so that value is computed and discarded. `AS_FisherStrand` filters the reference out
//! first. The two siblings differ in exactly that line.
//!
//! # `AS_FS` floors its p-value at `1e-320`, which is a subnormal
//!
//! ```java
//! QualityUtils.phredScaleErrorRate(Math.max(FisherStrand.pValueForContingencyTable(refAltTable), MIN_PVALUE))
//! ```
//!
//! `1.0E-320` is below `Double.MIN_NORMAL`, so the floor is a denormalised double and the phred
//! scale of it is about 3199. The clamp inside `phredScaleLog10ErrorRate` is at
//! `log10(Double.MIN_VALUE)`, about -323.3, which is lower still, so the `MIN_PVALUE` floor is the
//! one that binds and the value never reaches 3233.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

use crate::rank_sum::format_decimals;
use crate::strand_bias::{calculate_sor, p_value_for_contingency_table};

/// `GATKVCFConstants.AS_SB_TABLE_KEY`.
pub const AS_SB_TABLE_KEY: &str = "AS_SB_TABLE";
/// `GATKVCFConstants.AS_FISHER_STRAND_KEY`.
pub const AS_FISHER_STRAND_KEY: &str = "AS_FS";
/// `GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY`.
pub const AS_STRAND_ODDS_RATIO_KEY: &str = "AS_SOR";

/// `StrandBiasUtils.FORWARD` and `REVERSE`, which are indices into a two-element list.
pub const FORWARD: usize = 0;
pub const REVERSE: usize = 1;
/// `StrandBiasUtils.MIN_COUNT`.
pub const MIN_COUNT: i32 = 2;
/// `AS_StrandBiasTest.MIN_PVALUE`.
pub const MIN_PVALUE: f64 = 1.0E-320;
/// `AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM`.
const RAW_DELIM: char = '|';
/// `AS_StrandBiasTest.REDUCED_DELIM`.
const REDUCED_DELIM: char = ',';
/// `VCFConstants.MISSING_VALUE_v4`.
const MISSING_VALUE: &str = ".";

/// What this family refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum AsStrandBiasError {
    /// `IllegalStateException`: "Number of alleles and number of allele-specific entries do not
    /// match. Allele-specific annotations should have an entry for each allele including the
    /// reference."
    AlleleCountMismatch { entries: usize, alleles: usize },
    /// `NumberFormatException` from `Integer.parseInt` on a non-integer count.
    MalformedCount { raw: String },
    /// No reference allele, or more than one.
    ReferenceAlleleCount { found: usize },
}

/// The two members, which differ in their key and in whether they filter the reference out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsStrandBias {
    Fisher,
    OddsRatio,
}

impl AsStrandBias {
    pub fn vcf_key(&self) -> &'static str {
        match self {
            AsStrandBias::Fisher => AS_FISHER_STRAND_KEY,
            AsStrandBias::OddsRatio => AS_STRAND_ODDS_RATIO_KEY,
        }
    }

    /// `getPrimaryRawKey()`, which is the **same** key for both: they share one raw table.
    pub fn raw_key(&self) -> &'static str {
        AS_SB_TABLE_KEY
    }

    /// `getEmptyRawValue()`, which is a pair of zeros rather than the empty string the rank sums
    /// use.
    pub fn empty_raw_value(&self) -> &'static str {
        "0,0"
    }
}

/// `AlleleSpecificAnnotationData`'s reference-allele check.
fn reference_allele(alleles: &[Allele]) -> Result<&Allele, AsStrandBiasError> {
    let references: Vec<&Allele> = alleles.iter().filter(|a| a.is_reference()).collect();
    match references.len() {
        1 => Ok(references[0]),
        found => Err(AsStrandBiasError::ReferenceAlleleCount { found }),
    }
}

/// `StrandBiasUtils.getStrandCountsFromLikelihoodMap`: one forward/reverse pair per allele.
///
/// Absent rather than zero: an allele no read was assigned to has **no** entry, which is what makes
/// the difference between `0,0` and an empty entry in the raw string.
pub fn strand_counts(
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    min_count: i32,
) -> Vec<(Allele, Option<[i32; 2]>)> {
    let mut combined: Vec<(Allele, Option<[i32; 2]>)> = vc
        .alleles
        .iter()
        .map(|allele| (allele.clone(), None))
        .collect();
    let Some(likelihoods) = likelihoods else {
        return combined;
    };
    let reference = vc.alleles.iter().find(|a| a.is_reference());

    for sample_index in 0..likelihoods.number_of_samples() {
        let mut sample_table: Vec<(Allele, Option<[i32; 2]>)> = vc
            .alleles
            .iter()
            .map(|allele| (allele.clone(), None))
            .collect();
        for best in likelihoods.best_alleles_breaking_ties_for_sample(sample_index, None) {
            if !best.is_informative() {
                continue;
            }
            let Some(allele) = best.allele.as_ref() else {
                continue;
            };
            // "can happen if a read's most likely allele has been removed when
            // --max_alternate_alleles is exceeded": neither the reference nor a declared
            // alternate, so the read is dropped rather than counted somewhere.
            let matches_reference = reference == Some(allele);
            let matches_alternate = vc
                .alternate_alleles()
                .iter()
                .any(|alternate| alternate == allele);
            if !(matches_reference || matches_alternate) {
                continue;
            }
            let Some(read) = likelihoods
                .sample_evidence(sample_index)
                .and_then(|reads| reads.get(best.evidence_index))
            else {
                continue;
            };
            let strand = if is_reverse_strand(read) {
                REVERSE
            } else {
                FORWARD
            };
            if let Some(slot) = sample_table.iter_mut().find(|(a, _)| a == allele) {
                let counts = slot.1.get_or_insert([0, 0]);
                counts[strand] += 1;
            }
        }
        // The threshold is on the sample's whole table, and it is strictly greater than.
        let read_count: i32 = sample_table
            .iter()
            .filter_map(|(_, counts)| counts.as_ref())
            .map(|counts| counts[FORWARD] + counts[REVERSE])
            .sum();
        if read_count > min_count {
            for (allele, counts) in &sample_table {
                let Some(counts) = counts else { continue };
                if let Some(slot) = combined.iter_mut().find(|(a, _)| a == allele) {
                    match &mut slot.1 {
                        Some(existing) => {
                            existing[FORWARD] += counts[FORWARD];
                            existing[REVERSE] += counts[REVERSE];
                        }
                        None => slot.1 = Some(*counts),
                    }
                }
            }
        }
    }
    combined
}

/// `GATKRead.isReverseStrand()`: the `0x10` flag.
fn is_reverse_strand(read: &BamRecord) -> bool {
    read.flags & 0x10 != 0
}

/// `StrandBiasUtils.makeRawAnnotationString`: the delimiter goes **between** entries, and a missing
/// allele is written as the shared `ZERO_LIST`.
pub fn make_raw_annotation_string(per_allele: &[(Allele, Option<[i32; 2]>)]) -> String {
    let mut out = String::new();
    for (_, counts) in per_allele {
        if !out.is_empty() {
            out.push(RAW_DELIM);
        }
        match counts {
            // `ZERO_LIST` is a `static` mutable `ArrayList` shared by every call site; nothing
            // writes through it, so the sharing is safe and remains one edit away from not being.
            None => out.push_str("0,0"),
            Some(counts) => out.push_str(&format!("{},{}", counts[FORWARD], counts[REVERSE])),
        }
    }
    out
}

/// `AS_StrandBiasTest.annotateRawData`.
pub fn annotate_raw_data(
    annotation: AsStrandBias,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Option<(String, String)> {
    // "for allele-specific annotations we only call from HC and we only use likelihoods".
    likelihoods?;
    let counts = strand_counts(vc, likelihoods, MIN_COUNT);
    Some((
        annotation.raw_key().to_string(),
        make_raw_annotation_string(&counts),
    ))
}

/// `AS_StrandBiasTest.parseRawDataString`.
///
/// An entry may hold **any** number of comma-separated integers, or none; only the count of
/// pipe-separated entries is checked, and it must equal the allele count exactly.
pub fn parse_raw_data_string(
    alleles: &[Allele],
    raw: &str,
) -> Result<Vec<(Allele, Vec<i32>)>, AsStrandBiasError> {
    // `getAlleleLengthListOfStringFromRawData`: a **leading** bracket triggers stripping both ends
    // and removing every whitespace character, and split keeps trailing empties.
    let cleaned = if raw.starts_with('[') {
        raw[1..raw.len() - 1]
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<String>()
    } else {
        raw.to_string()
    };
    let entries: Vec<&str> = cleaned.split(RAW_DELIM).collect();
    if entries.len() != alleles.len() {
        return Err(AsStrandBiasError::AlleleCountMismatch {
            entries: entries.len(),
            alleles: alleles.len(),
        });
    }
    let mut out = Vec::with_capacity(alleles.len());
    for (allele, entry) in alleles.iter().zip(entries) {
        let mut counts = Vec::new();
        for piece in entry.split(',') {
            if piece.is_empty() {
                continue;
            }
            let value =
                piece
                    .trim()
                    .parse::<i32>()
                    .map_err(|_| AsStrandBiasError::MalformedCount {
                        raw: piece.to_string(),
                    })?;
            counts.push(value);
        }
        out.push((allele.clone(), counts));
    }
    Ok(out)
}

/// `AS_StrandBiasTest.combineRawData`: the tables added componentwise.
pub fn combine_raw_data(
    alleles: &[Allele],
    raw_strings: &[String],
) -> Result<String, AsStrandBiasError> {
    let mut combined: Vec<(Allele, Option<[i32; 2]>)> = alleles
        .iter()
        .map(|allele| (allele.clone(), None))
        .collect();
    for raw in raw_strings {
        let parsed = parse_raw_data_string(alleles, raw)?;
        for (allele, counts) in &parsed {
            let Some(slot) = combined.iter_mut().find(|(a, _)| a == allele) else {
                continue;
            };
            match &mut slot.1 {
                Some(existing) => {
                    // `combineAttributeMap` reads index 0 and 1 unconditionally once the entry
                    // exists, so a one-element entry here is an IndexOutOfBounds upstream. Only an
                    // empty entry is representable, and it is handled below.
                    if counts.len() >= 2 {
                        existing[FORWARD] += counts[FORWARD];
                        existing[REVERSE] += counts[REVERSE];
                    }
                }
                None => {
                    slot.1 = if counts.len() >= 2 {
                        Some([counts[FORWARD], counts[REVERSE]])
                    } else {
                        // An empty entry puts an **empty list** in the map, which is not the same
                        // as an absent one: it renders as `0,0` here but as the missing value once
                        // reduced.
                        Some([0, 0])
                    };
                }
            }
        }
    }
    Ok(make_raw_annotation_string(&combined))
}

/// `QualityUtils.phredScaleErrorRate`.
///
/// The `Math.abs` is there because an error rate of exactly one gives a log of `-0.0`, and `-10 *
/// -0.0` is `0.0` while `-10 * 0.0` is `-0.0`, which would print with a sign.
pub fn phred_scale_error_rate(error_rate: f64) -> f64 {
    let log10 = jmath::math::log10(error_rate);
    // `MIN_LOG10_SCALED_QUAL = Math.log10(Double.MIN_VALUE)`, about -323.3.
    let min_log10 = jmath::math::log10(f64::from_bits(1));
    (-10.0 * gatk_engine::math_utils::java_max(log10, min_log10)).abs()
}

/// `AS_FisherStrand.calculateReducedData` and `AS_StrandOddsRatio.calculateReducedData`.
///
/// `None` in the value position is the reference's null, which becomes the missing value. The two
/// members differ only in whether the reference allele is skipped; see the module note.
pub fn calculate_reduced_data(
    annotation: AsStrandBias,
    per_allele: &[(Allele, Vec<i32>)],
    reference: &Allele,
) -> Vec<(Allele, Option<f64>)> {
    let reference_counts: Vec<i32> = per_allele
        .iter()
        .find(|(a, _)| a == reference)
        .map(|(_, counts)| counts.clone())
        .unwrap_or_default();
    per_allele
        .iter()
        .filter(|(allele, _)| annotation == AsStrandBias::OddsRatio || allele != reference)
        .map(|(allele, counts)| {
            if counts.is_empty() {
                return (allele.clone(), None);
            }
            let table = [
                [reference_counts[FORWARD], reference_counts[REVERSE]],
                [counts[FORWARD], counts[REVERSE]],
            ];
            let value = match annotation {
                AsStrandBias::Fisher => phred_scale_error_rate(gatk_engine::math_utils::java_max(
                    p_value_for_contingency_table(table),
                    MIN_PVALUE,
                )),
                AsStrandBias::OddsRatio => calculate_sor(table),
            };
            (allele.clone(), Some(value))
        })
        .collect()
}

/// `AS_StrandBiasTest.makeReducedAnnotationString`.
///
/// An alternate the annotation data does not carry is **skipped entirely** rather than written as
/// the missing value, so the field can come out with fewer entries than the variant has alternates.
/// That is the opposite of what the rank sums do with the same situation, and it is what the
/// reference's own log line calls out as an error without acting on it.
pub fn make_reduced_annotation_string(
    alternates: &[Allele],
    per_allele: &[(Allele, Option<f64>)],
) -> String {
    let mut out = String::new();
    for allele in alternates {
        match per_allele.iter().find(|(a, _)| a == allele) {
            // "ERROR: VC allele not found in annotation alleles -- maybe there was trimming?",
            // logged and then nothing appended, not even a delimiter.
            None => continue,
            Some((_, value)) => {
                if !out.is_empty() {
                    out.push(REDUCED_DELIM);
                }
                match value {
                    None => out.push_str(MISSING_VALUE),
                    Some(value) => out.push_str(&format_decimals(*value, 3)),
                }
            }
        }
    }
    out
}

/// `AS_StrandBiasTest.finalizeRawData`: the reduced value and the raw table, both.
///
/// The raw string it writes back is built over the **current** variant's alleles while the data was
/// parsed over the original's, so an allele the merge dropped disappears from the table and one it
/// added is written as `0,0`.
pub fn finalize_raw_data(
    annotation: AsStrandBias,
    vc_alleles: &[Allele],
    original_alleles: &[Allele],
    raw: Option<&str>,
) -> Result<Option<(String, String, String, String)>, AsStrandBiasError> {
    let Some(raw) = raw else { return Ok(None) };
    let reference = reference_allele(original_alleles)?.clone();
    let parsed = parse_raw_data_string(original_alleles, raw)?;
    let reduced_values = calculate_reduced_data(annotation, &parsed, &reference);

    let alternates: Vec<Allele> = vc_alleles
        .iter()
        .filter(|a| !a.is_reference())
        .cloned()
        .collect();
    let reduced = make_reduced_annotation_string(&alternates, &reduced_values);

    let raw_out = make_raw_annotation_string_preserving_empties(&parsed, vc_alleles);
    Ok(Some((
        annotation.vcf_key().to_string(),
        reduced,
        annotation.raw_key().to_string(),
        raw_out,
    )))
}

/// `makeRawAnnotationString` over a parsed map, where an entry that exists but is empty renders as
/// the empty string and only an **absent** entry falls back to `0,0`.
fn make_raw_annotation_string_preserving_empties(
    parsed: &[(Allele, Vec<i32>)],
    vc_alleles: &[Allele],
) -> String {
    let mut out = String::new();
    for allele in vc_alleles {
        if !out.is_empty() {
            out.push(RAW_DELIM);
        }
        match parsed.iter().find(|(a, _)| a == allele) {
            None => out.push_str("0,0"),
            Some((_, counts)) => {
                let text: Vec<String> = counts.iter().map(|c| c.to_string()).collect();
                out.push_str(&text.join(","));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alleles() -> Vec<Allele> {
        vec![
            Allele::from_str("A", true).expect("an allele"),
            Allele::from_str("C", false).expect("an allele"),
        ]
    }

    #[test]
    fn the_raw_string_has_an_entry_for_the_reference_and_no_leading_delimiter() {
        let alleles = alleles();
        let counts = vec![
            (alleles[0].clone(), Some([10, 8])),
            (alleles[1].clone(), Some([3, 4])),
        ];
        assert_eq!(make_raw_annotation_string(&counts), "10,8|3,4");
    }

    #[test]
    fn a_count_mismatch_is_refused_rather_than_padded() {
        assert_eq!(
            parse_raw_data_string(&alleles(), "10,8"),
            Err(AsStrandBiasError::AlleleCountMismatch {
                entries: 1,
                alleles: 2
            })
        );
    }

    #[test]
    fn the_odds_ratio_keeps_the_reference_and_fisher_drops_it() {
        let alleles = alleles();
        let parsed = vec![
            (alleles[0].clone(), vec![10, 8]),
            (alleles[1].clone(), vec![3, 4]),
        ];
        let fisher = calculate_reduced_data(AsStrandBias::Fisher, &parsed, &alleles[0]);
        let sor = calculate_reduced_data(AsStrandBias::OddsRatio, &parsed, &alleles[0]);
        assert_eq!(fisher.len(), 1);
        // The reference's own row survives, and is a symmetric table against itself.
        assert_eq!(sor.len(), 2);
    }
}
