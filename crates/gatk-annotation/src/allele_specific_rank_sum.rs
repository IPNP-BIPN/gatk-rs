//! `AS_BaseQualityRankSumTest`, `AS_MappingQualityRankSumTest` and `AS_ReadPosRankSumTest`, and the
//! `AS_RankSumTest` machinery under them, ported from GATK 4.6.2.0.
//!
//! `AS_BaseQRankSum`, `AS_MQRankSum` and `AS_ReadPosRankSum`, plus the three `AS_RAW_` keys the
//! gVCF path writes so that combining two samples is adding two histograms.
//!
//! # The direct path of an allele-specific annotation is **not** allele-specific
//!
//! ```java
//! final MannWhitneyU.Result result = mannWhitneyU.test(
//!         Doubles.toArray(altQuals), Doubles.toArray(refQuals), MannWhitneyU.TestType.FIRST_DOMINATES);
//! ```
//!
//! `AS_RankSumTest.annotate` pools every alternate allele's reads into one series and reports one Z
//! score, exactly as its non-allele-specific parent does. The two differ only in the key they write
//! it under, so `AS_MQRankSum` and `MQRankSum` on the same site carry the same number. The
//! allele-specific part lives entirely in [`annotate_raw_data`] and [`finalize_raw_data`], which the
//! HaplotypeCaller reaches in gVCF mode and nothing else does.
//!
//! # The raw string starts with its delimiter, and the parser depends on that
//!
//! ```java
//! if (!vcAlleles.get(i).isReference()) {
//!     if (i != 0) { //strings will always start with a printDelim ...
//!         annotationString += AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM;
//!     }
//! ```
//!
//! The reference allele is skipped but its **slot is not**: the writer emits a leading `|`, so
//! splitting the string yields an empty first token that the parser assigns back to the reference.
//! Position, not content, is what carries the allele identity through a gVCF, and a raw string
//! whose leading delimiter were trimmed would silently shift every allele by one.
//!
//! # A site with no unambiguous reference read produces **no** raw annotation
//!
//! ```java
//! if (perAlleleValues.get(ref).isEmpty()) { return perAltRankSumResults; }
//! ```
//!
//! Not a missing value per allele: an empty map, which makes the raw string a bare run of
//! delimiters. The finalising step then has nothing to take a median of.
//!
//! # The opposite case writes the four characters `NaN` into the record
//!
//! A site with reference reads and **no** alternate read takes a rank sum of an empty series
//! against a non-empty one, which is `NaN`. `Histogram.add` drops a `NaN` silently, so the
//! histogram is empty, and an empty histogram renders as `Double.toString(Double.NaN)`. The raw
//! field is then `|NaN`, and the golden carries that row. The guard in
//! [`make_combined_annotation_string`] keeps those four characters out of the *combined* string,
//! but there is no such guard on the way in.
//!
//! # Each Z score is stored as a one-entry histogram, so the raw form is lossy
//!
//! ```java
//! public String outputSingletonValueAsHistogram(final Double rankSumValue) {
//!     Histogram h = new Histogram(); h.add(rankSumValue); return h.toString();
//! }
//! ```
//!
//! The Z score is binned to a tenth before it is written, so `-1.23` and `-1.27` both come back as
//! `-1.3`: the bin index is a **floor**, so a negative score bins away from zero and the binning is
//! not symmetric about it. Combining gVCFs then takes the **median of the bins**, which is why
//! `AS_MQRankSum` from a joint call is not the `MQRankSum` of the pooled reads and cannot be.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::histogram::{CompressedDataList, Histogram, HistogramError};
use gatk_engine::mann_whitney::{self, TestType};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

use crate::rank_sum::{
    format_decimals, BaseQualityRankSumTest, MappingQualityRankSumTest, RankSumTest,
    ReadPosRankSumTest, INVALID_ELEMENT_FROM_READ,
};

/// `GATKVCFConstants.AS_BASE_QUAL_RANK_SUM_KEY` and its raw form.
pub const AS_BASE_QUAL_RANK_SUM_KEY: &str = "AS_BaseQRankSum";
pub const AS_RAW_BASE_QUAL_RANK_SUM_KEY: &str = "AS_RAW_BaseQRankSum";
/// `GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY` and its raw form.
pub const AS_MAP_QUAL_RANK_SUM_KEY: &str = "AS_MQRankSum";
pub const AS_RAW_MAP_QUAL_RANK_SUM_KEY: &str = "AS_RAW_MQRankSum";
/// `GATKVCFConstants.AS_READ_POS_RANK_SUM_KEY` and its raw form.
pub const AS_READ_POS_RANK_SUM_KEY: &str = "AS_ReadPosRankSum";
pub const AS_RAW_READ_POS_RANK_SUM_KEY: &str = "AS_RAW_ReadPosRankSum";

/// `AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM`.
pub const ALLELE_SPECIFIC_RAW_DELIM: char = '|';
/// `AnnotationUtils.ALLELE_SPECIFIC_REDUCED_DELIM`.
pub const ALLELE_SPECIFIC_REDUCED_DELIM: char = ',';
/// `VCFConstants.MISSING_VALUE_v4`.
pub const MISSING_VALUE: &str = ".";
/// `AS_RankSumTest.RAW_DELIM`.
const RAW_DELIM: char = ',';

/// What this family refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum AsRankSumError {
    /// `IllegalStateException`: the raw path needs a variant context with exactly one sample, "as
    /// in a gVCF".
    NotExactlyOneSample { found: usize },
    /// `IllegalArgumentException` from `AlleleSpecificAnnotationData`: no reference allele, or more
    /// than one.
    ReferenceAlleleCount { found: usize },
    /// `ArrayIndexOutOfBoundsException`: the raw string carries more allele slots than the variant
    /// has alleles.
    TooManyAlleleSlots { slots: usize, alleles: usize },
    /// A bin index the histogram refuses.
    Histogram(HistogramError),
}

/// The three members, which differ only in their keys and in the element they take from a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsRankSum {
    BaseQuality,
    MappingQuality,
    ReadPosition,
}

impl AsRankSum {
    /// `getPrimaryRawKey()`.
    pub fn raw_key(&self) -> &'static str {
        match self {
            AsRankSum::BaseQuality => AS_RAW_BASE_QUAL_RANK_SUM_KEY,
            AsRankSum::MappingQuality => AS_RAW_MAP_QUAL_RANK_SUM_KEY,
            AsRankSum::ReadPosition => AS_RAW_READ_POS_RANK_SUM_KEY,
        }
    }

    /// `getEmptyRawValue()`, which is the empty string and not the missing value.
    pub fn empty_raw_value(&self) -> &'static str {
        ""
    }
}

impl RankSumTest for AsRankSum {
    fn vcf_key(&self) -> &'static str {
        match self {
            AsRankSum::BaseQuality => AS_BASE_QUAL_RANK_SUM_KEY,
            AsRankSum::MappingQuality => AS_MAP_QUAL_RANK_SUM_KEY,
            AsRankSum::ReadPosition => AS_READ_POS_RANK_SUM_KEY,
        }
    }

    fn element_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<f64> {
        match self {
            AsRankSum::BaseQuality => BaseQualityRankSumTest.element_for_read(read, vc),
            AsRankSum::MappingQuality => MappingQualityRankSumTest.element_for_read(read, vc),
            AsRankSum::ReadPosition => ReadPosRankSumTest.element_for_read(read, vc),
        }
    }

    fn is_usable_read(&self, read: &BamRecord, vc: &VariantContext) -> bool {
        match self {
            // Only `AS_ReadPosRankSumTest` overrides the usability test, and it overrides it with
            // the same one its non-allele-specific sibling uses.
            AsRankSum::ReadPosition => ReadPosRankSumTest.is_usable_read(read, vc),
            AsRankSum::BaseQuality => BaseQualityRankSumTest.is_usable_read(read, vc),
            AsRankSum::MappingQuality => MappingQualityRankSumTest.is_usable_read(read, vc),
        }
    }
}

/// `AlleleSpecificAnnotationData`'s reference-allele check, which runs before any data is stored.
fn reference_allele(alleles: &[Allele]) -> Result<&Allele, AsRankSumError> {
    let references: Vec<&Allele> = alleles.iter().filter(|a| a.is_reference()).collect();
    match references.len() {
        1 => Ok(references[0]),
        found => Err(AsRankSumError::ReferenceAlleleCount { found }),
    }
}

/// `AS_RankSumTest.calculateRawData`: one run-length list of **truncated** values per allele.
///
/// The cast to `int` is a truncation towards zero, so a read position of `-0.5` becomes zero and a
/// base quality is unaffected. Every value this family produces is already integral except the
/// clipping one, which is not a member here.
pub fn calculate_raw_data(
    annotation: AsRankSum,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Result<Vec<(Allele, CompressedDataList)>, AsRankSumError> {
    let samples: std::collections::BTreeSet<&str> = vc
        .genotypes
        .iter()
        .map(|g| g.sample_name.as_str())
        .collect();
    if samples.len() != 1 {
        return Err(AsRankSumError::NotExactlyOneSample {
            found: samples.len(),
        });
    }
    let mut per_allele: Vec<(Allele, CompressedDataList)> = vc
        .alleles
        .iter()
        .map(|allele| (allele.clone(), CompressedDataList::new()))
        .collect();
    let Some(likelihoods) = likelihoods else {
        return Ok(per_allele);
    };

    for best in likelihoods.best_alleles_breaking_ties(None) {
        if !best.is_informative() {
            continue;
        }
        let Some(read) = likelihoods
            .sample_evidence(likelihoods.index_of_sample(&best.sample).unwrap_or(0))
            .and_then(|reads| reads.get(best.evidence_index))
        else {
            continue;
        };
        if !annotation.is_usable_read(read, vc) {
            continue;
        }
        let Some(value) = annotation.element_for_read(read, vc) else {
            continue;
        };
        if value == INVALID_ELEMENT_FROM_READ {
            continue;
        }
        let Some(allele) = best.allele else { continue };
        // `perAlleleValues.containsKey(bestAllele.allele)`: an allele the matrix holds and the
        // variant does not is dropped rather than counted into a new bucket.
        if let Some(slot) = per_allele.iter_mut().find(|(a, _)| *a == allele) {
            slot.1.add(value as i32);
        }
    }
    Ok(per_allele)
}

/// `AS_RankSumTest.calculateRankSum`: one Z score per alternate, or nothing at all.
pub fn calculate_rank_sum(
    per_allele: &[(Allele, CompressedDataList)],
    reference: &Allele,
) -> Vec<(Allele, f64)> {
    let ref_values: Vec<f64> = per_allele
        .iter()
        .find(|(a, _)| a == reference)
        .map(|(_, list)| list.iter().map(|v| v as f64).collect())
        .unwrap_or_default();
    // "shortcut to not try to calculate rank sum if there are no reads that unambiguously support
    // the ref": every alternate is dropped, not just the ones with no reads.
    if ref_values.is_empty() {
        return Vec::new();
    }
    per_allele
        .iter()
        .filter(|(allele, _)| allele != reference)
        .map(|(allele, list)| {
            let alts: Vec<f64> = list.iter().map(|v| v as f64).collect();
            let result = mann_whitney::test(&alts, &ref_values, TestType::FirstDominates);
            (allele.clone(), result.z)
        })
        .collect()
}

/// `AS_RankSumTest.outputSingletonValueAsHistogram`: a Z score binned to a tenth.
pub fn singleton_histogram(z: f64) -> Result<String, AsRankSumError> {
    let mut histogram = Histogram::new();
    histogram.add(z).map_err(AsRankSumError::Histogram)?;
    Ok(histogram.to_string())
}

/// `AS_RankSumTest.makeRawAnnotationString`: the leading delimiter is deliberate. See the module
/// note.
pub fn make_raw_annotation_string(
    alleles: &[Allele],
    per_allele: &[(Allele, f64)],
) -> Result<String, AsRankSumError> {
    let mut out = String::new();
    for (i, allele) in alleles.iter().enumerate() {
        if allele.is_reference() {
            continue;
        }
        if i != 0 {
            out.push(ALLELE_SPECIFIC_RAW_DELIM);
        }
        // Null when `calculateRankSum` returned nothing, in which case only the delimiters are
        // written and the field is a run of pipes.
        if let Some((_, z)) = per_allele.iter().find(|(a, _)| a == allele) {
            out.push_str(&singleton_histogram(*z)?);
        }
    }
    Ok(out)
}

/// `AS_RankSumTest.annotateRawData`.
pub fn annotate_raw_data(
    annotation: AsRankSum,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Result<Option<(String, String)>, AsRankSumError> {
    // `if (likelihoods == null) return emptyMap()` happens **before** the one-sample check, so a
    // null matrix is an empty answer and not an exception.
    if likelihoods.is_none() {
        return Ok(None);
    }
    let reference = reference_allele(&vc.alleles)?.clone();
    let per_allele = calculate_raw_data(annotation, vc, likelihoods)?;
    let rank_sums = calculate_rank_sum(&per_allele, &reference);
    let text = make_raw_annotation_string(&vc.alleles, &rank_sums)?;
    Ok(Some((annotation.raw_key().to_string(), text)))
}

/// `AS_RankSumTest.parseRawDataString`: positional, over the alleles the caller supplies.
///
/// The raw string's tokens line up with the allele list index for index, the reference's token
/// being the empty one the leading delimiter creates. More tokens than alleles is the reference's
/// `ArrayIndexOutOfBoundsException`.
pub fn parse_raw_data_string(
    alleles: &[Allele],
    raw: &str,
) -> Result<Vec<(Allele, Histogram)>, AsRankSumError> {
    reference_allele(alleles)?;
    let mut per_allele: Vec<(Allele, Histogram)> = alleles
        .iter()
        .map(|allele| (allele.clone(), Histogram::new()))
        .collect();

    // "Map gives back list with []": only a **leading** bracket is tested, and then both ends are
    // trimmed, so a string with a trailing bracket and no leading one keeps it.
    let without_brackets = if raw.starts_with('[') {
        &raw[1..raw.len() - 1]
    } else {
        raw
    };

    let tokens = split_dropping_trailing_empties(without_brackets, ALLELE_SPECIFIC_RAW_DELIM);
    if tokens.len() > alleles.len() {
        return Err(AsRankSumError::TooManyAlleleSlots {
            slots: tokens.len(),
            alleles: alleles.len(),
        });
    }
    for (index, token) in tokens.iter().enumerate() {
        let entries: Vec<&str> = token.split(RAW_DELIM).collect();
        let mut j = 0usize;
        while j < entries.len() {
            if !entries[j].is_empty() {
                let Ok(value) = entries[j].trim().parse::<f64>() else {
                    // `Double.parseDouble` throws a NumberFormatException, which the reference does
                    // not catch. No writer of this field produces one.
                    j += 2;
                    continue;
                };
                // A NaN bin is skipped here, which is the only reason `Histogram.add(d, count)`
                // never sees one despite having no guard of its own.
                if !value.is_nan() && j + 1 < entries.len() && !entries[j + 1].is_empty() {
                    let count: i32 = entries[j + 1].trim().parse().unwrap_or(0);
                    per_allele[index]
                        .1
                        .add_count(value, count)
                        .map_err(AsRankSumError::Histogram)?;
                }
            }
            j += 2;
        }
    }
    Ok(per_allele)
}

/// `String.split(regex)`, which drops trailing empty tokens but keeps leading and interior ones.
fn split_dropping_trailing_empties(text: &str, delimiter: char) -> Vec<String> {
    if text.is_empty() {
        // Java's split on an empty input yields one empty token, not none.
        return vec![String::new()];
    }
    let mut parts: Vec<String> = text.split(delimiter).map(|s| s.to_string()).collect();
    while parts.last().is_some_and(|part| part.is_empty()) {
        parts.pop();
    }
    parts
}

/// `AS_RankSumTest.combineAttributeMap` over a list of raw strings.
pub fn combine_raw_data(
    alleles: &[Allele],
    raw_strings: &[String],
) -> Result<String, AsRankSumError> {
    let mut combined: Vec<(Allele, Histogram)> = alleles
        .iter()
        .map(|allele| (allele.clone(), Histogram::new()))
        .collect();
    for raw in raw_strings {
        let parsed = parse_raw_data_string(alleles, raw)?;
        for (allele, histogram) in &parsed {
            if let Some(slot) = combined.iter_mut().find(|(a, _)| a == allele) {
                slot.1
                    .add_histogram(histogram)
                    .map_err(AsRankSumError::Histogram)?;
            }
        }
    }
    make_combined_annotation_string(alleles, &combined)
}

/// `AS_RankSumTest.makeCombinedAnnotationString`.
pub fn make_combined_annotation_string(
    alleles: &[Allele],
    per_allele: &[(Allele, Histogram)],
) -> Result<String, AsRankSumError> {
    let mut out = String::new();
    for (i, allele) in alleles.iter().enumerate() {
        if allele.is_reference() {
            continue;
        }
        if i != 0 {
            out.push(ALLELE_SPECIFIC_RAW_DELIM);
        }
        if let Some((_, histogram)) = per_allele.iter().find(|(a, _)| a == allele) {
            // An empty histogram would render as `NaN`, so the guard is what keeps those four
            // characters out of the field.
            if !histogram.is_empty() {
                out.push_str(&histogram.to_string());
            }
        }
    }
    Ok(out)
}

/// `AS_RankSumTest.finalizeRawData`: the per-allele medians, and the combined raw string.
///
/// `None` when the raw key is absent. The reduced string is written over the **current** variant's
/// alternates while the histograms were parsed over the **original** variant's alleles, so an
/// alternate that survived trimming under a different representation is written as the missing
/// value rather than dropped, keeping the field's arity equal to the alternate count.
pub fn finalize_raw_data(
    annotation: AsRankSum,
    vc_alternates: &[Allele],
    original_alleles: &[Allele],
    raw: Option<&str>,
) -> Result<Option<(String, String, String, String)>, AsRankSumError> {
    let Some(raw) = raw else { return Ok(None) };
    let reference = reference_allele(original_alleles)?.clone();
    let per_allele = parse_raw_data_string(original_alleles, raw)?;

    let medians: Vec<(Allele, Option<f64>)> = per_allele
        .iter()
        .filter(|(allele, _)| *allele != reference)
        .map(|(allele, histogram)| (allele.clone(), histogram.median()))
        .collect();
    // "shortcut for no ref values": an empty map here is an empty answer, both keys absent.
    if medians.is_empty() {
        return Ok(None);
    }

    let mut reduced = String::new();
    for allele in vc_alternates {
        if !reduced.is_empty() {
            reduced.push(ALLELE_SPECIFIC_REDUCED_DELIM);
        }
        match medians.iter().find(|(a, _)| a == allele) {
            // "VC allele not found in annotation alleles -- maybe there was trimming?"
            None => reduced.push_str(MISSING_VALUE),
            Some((_, None)) => reduced.push_str(MISSING_VALUE),
            Some((_, Some(median))) => reduced.push_str(&format_decimals(*median, 3)),
        }
    }
    let combined = make_combined_annotation_string(original_alleles, &per_allele)?;
    Ok(Some((
        annotation.vcf_key().to_string(),
        reduced,
        annotation.raw_key().to_string(),
        combined,
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alleles() -> Vec<Allele> {
        vec![
            Allele::from_str("A", true).expect("an allele"),
            Allele::from_str("C", false).expect("an allele"),
            Allele::from_str("G", false).expect("an allele"),
        ]
    }

    #[test]
    fn the_raw_string_starts_with_its_delimiter_and_round_trips() {
        let alleles = alleles();
        let raw = make_raw_annotation_string(
            &alleles,
            &[(alleles[1].clone(), -1.23), (alleles[2].clone(), 0.4)],
        )
        .expect("a raw string");
        // Binned to a tenth by a **floor**, so a negative Z score bins away from zero: -1.23
        // becomes -1.3 and not -1.2. The binning is not symmetric about zero.
        assert_eq!(raw, "|-1.3,1|0.4,1");
        let parsed = parse_raw_data_string(&alleles, &raw).expect("a parse");
        assert!(parsed[0].1.is_empty(), "the reference keeps an empty slot");
        let median = parsed[1].1.median().expect("a median");
        assert!((median + 1.3).abs() < 1e-12, "{median}");
    }

    #[test]
    fn a_site_with_no_reference_read_produces_only_delimiters() {
        let alleles = alleles();
        let raw = make_raw_annotation_string(&alleles, &[]).expect("a raw string");
        assert_eq!(raw, "||");
    }

    #[test]
    fn the_reference_allele_count_is_checked_before_anything_is_stored() {
        let two_references = vec![
            Allele::from_str("A", true).expect("an allele"),
            Allele::from_str("C", true).expect("an allele"),
        ];
        assert_eq!(
            parse_raw_data_string(&two_references, "|"),
            Err(AsRankSumError::ReferenceAlleleCount { found: 2 })
        );
    }
}
