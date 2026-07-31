//! `StrandBiasTest` and its three members, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator` (GATK 4.6.2.0).
//!
//! `FS`, `SOR` and `SB`: the same 2x2 contingency table of forward and reverse reads on the
//! reference and the alternate, read three ways.
//!
//! # The table is built per sample and only then added up
//!
//! ```java
//! for (final String sample : samples) {
//!     final int[] sampleTable = new int[ARRAY_SIZE];
//!     ...
//!     if (passesMinimumThreshold(sampleTable, minCount)) { copyToMainTable(sampleTable, table); }
//! }
//! ```
//!
//! A sample under the threshold contributes **nothing**, rather than contributing its counts to a
//! pooled total that would pass. So the same reads split across two samples can produce a
//! different table from the same reads in one, and `FS` and `SOR` disagree about the threshold:
//! `FisherStrand` uses `MIN_COUNT = 2` and `StrandOddsRatio` uses `0`.
//!
//! # The genotype field wins over the likelihoods
//!
//! If any genotype carries `SB`, the annotation is computed from those per-sample arrays and the
//! likelihood matrix is not consulted at all. That is what makes `SB` (`StrandBiasBySample`) a
//! load-bearing annotation rather than a diagnostic: writing it changes what `FS` and `SOR`
//! compute later.
//!
//! # `FS` normalises the table, and the normalisation truncates
//!
//! ```java
//! final double normFactor = sum / TARGET_TABLE_SIZE;
//! return new int[][]{{(int) (table[0][0] / normFactor), ...}};
//! ```
//!
//! Above 400 total reads the counts are scaled to about 200 and **truncated** to int, so `FS` on a
//! deep site is computed from a table that no longer sums to the coverage. `SOR` does not
//! normalise at all.
//!
//! # Both values are strings, and `FS` is phred-scaled off a floor
//!
//! ```java
//! String.format("%.3f", QualityUtils.phredScaleErrorRate(Math.max(pValue, MIN_PVALUE)))
//! ```
//!
//! `MIN_PVALUE = 1e-320` exists "to prevent INFINITYs", so the largest `FS` a site can report is
//! about 3200 whatever the evidence.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::fisher_exact;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Value, VariantContext};

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};

/// `GATKVCFConstants.FISHER_STRAND_KEY`.
pub const FISHER_STRAND_KEY: &str = "FS";
/// `GATKVCFConstants.STRAND_ODDS_RATIO_KEY`.
pub const STRAND_ODDS_RATIO_KEY: &str = "SOR";
/// `GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY`.
pub const STRAND_BIAS_BY_SAMPLE_KEY: &str = "SB";

/// `FisherStrand.MIN_PVALUE`, which caps `FS` at about 3200.
const MIN_PVALUE: f64 = 1e-320;
/// `FisherStrand.TARGET_TABLE_SIZE`.
const TARGET_TABLE_SIZE: f64 = 200.0;
/// `StrandOddsRatio.PSEUDOCOUNT`.
const PSEUDOCOUNT: f64 = 1.0;

/// `StrandBiasTest.getContingencyTable`, per sample and then summed.
///
/// The rows are reference and alternate; the columns are forward and reverse.
pub fn contingency_table(
    likelihoods: &AlleleLikelihoods<BamRecord>,
    vc: &VariantContext,
    min_count: i32,
) -> [[i32; 2]; 2] {
    let mut table = [[0i32; 2]; 2];
    let Some(reference) = vc.alleles.iter().find(|a| a.is_reference()) else {
        return table;
    };
    let alternates: Vec<&Allele> = vc.alleles.iter().filter(|a| !a.is_reference()).collect();

    for sample_index in 0..likelihoods.number_of_samples() {
        let mut sample_table = [0i32; 4];
        for best in likelihoods.best_alleles_breaking_ties_for_sample(sample_index, None) {
            if !best.is_informative() {
                continue;
            }
            let Some(read) = likelihoods
                .sample_evidence(sample_index)
                .and_then(|reads| reads.get(best.evidence_index))
            else {
                continue;
            };
            let Some(allele) = best.allele else { continue };
            // `allele.equals(ref, true)` ignores the reference flag, so an allele with the same
            // bases counts as the reference whichever way it was constructed.
            let matches_ref = allele.display_string() == reference.display_string();
            let matches_alt = alternates.iter().any(|alt| **alt == allele);
            if !matches_ref && !matches_alt {
                continue;
            }
            let offset = if matches_ref { 0 } else { 2 };
            let is_forward = !is_reverse_strand(read);
            sample_table[offset + usize::from(!is_forward)] += 1;
        }
        // A sample under the threshold contributes nothing at all.
        if passes_minimum_threshold(&sample_table, min_count) {
            table[0][0] += sample_table[0];
            table[0][1] += sample_table[1];
            table[1][0] += sample_table[2];
            table[1][1] += sample_table[3];
        }
    }
    table
}

/// `SAMRecord.getReadNegativeStrandFlag`, which is bit 0x10.
fn is_reverse_strand(read: &BamRecord) -> bool {
    read.flags & 0x10 != 0
}

/// `passesMinimumThreshold`: the **total** of all four cells, strictly greater than the count.
fn passes_minimum_threshold(data: &[i32; 4], min_count: i32) -> bool {
    data[0] + data[1] + data[2] + data[3] > min_count
}

/// `StrandBiasTest.getTableFromSamples`, which reads the per-sample `SB` fields.
///
/// `None` is the reference's null: no genotype carried `SB`, so the caller falls through to the
/// likelihoods.
pub fn table_from_samples(vc: &VariantContext, min_count: i32) -> Option<[[i32; 2]; 2]> {
    let mut totals = [0i32; 4];
    let mut found = false;
    for genotype in &vc.genotypes {
        let Some(counts) = strand_counts(genotype) else {
            continue;
        };
        found = true;
        if passes_minimum_threshold(&counts, min_count) {
            for (index, value) in counts.iter().enumerate() {
                totals[index] += value;
            }
        }
    }
    if !found {
        return None;
    }
    // `decodeSBBS`: the flat array becomes the 2x2 table in order.
    Some([[totals[0], totals[1]], [totals[2], totals[3]]])
}

/// `getStrandCounts`: the `SB` field, which the reference accepts as a String, a list of Integers
/// or a list of Strings, and refuses as anything else.
fn strand_counts(genotype: &htsjdk_vcf::variant::Genotype) -> Option<[i32; 4]> {
    let value = genotype
        .extended
        .iter()
        .find(|(key, _)| key == STRAND_BIAS_BY_SAMPLE_KEY)
        .map(|(_, value)| value)?;
    let parts: Vec<String> = match value {
        Value::Str(text) => text.split(',').map(|part| part.to_string()).collect(),
        Value::Int(number) => vec![number.to_string()],
        Value::List(values) => values.iter().map(render_value).collect(),
        other => vec![render_value(other)],
    };
    let mut counts = [0i32; 4];
    for (index, slot) in counts.iter_mut().enumerate() {
        // `Integer.parseInt(tokenizer.nextToken())`: fewer than four tokens is a
        // NoSuchElementException, which nothing here can produce from a written record.
        *slot = parts.get(index)?.trim().parse().ok()?;
    }
    Some(counts)
}

fn render_value(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        Value::Int(number) => number.to_string(),
        Value::Double(number) => number.to_string(),
        Value::Bool(flag) => flag.to_string(),
        Value::Missing => ".".to_string(),
        Value::List(values) => values
            .iter()
            .map(render_value)
            .collect::<Vec<_>>()
            .join(","),
    }
}

/// `FisherStrand.normalizeContingencyTable`, which truncates rather than rounds.
fn normalize_contingency_table(table: [[i32; 2]; 2]) -> [[i32; 2]; 2] {
    let sum = table[0][0] as i64 + table[0][1] as i64 + table[1][0] as i64 + table[1][1] as i64;
    if sum as f64 <= TARGET_TABLE_SIZE * 2.0 {
        return table;
    }
    let norm_factor = sum as f64 / TARGET_TABLE_SIZE;
    [
        [
            (table[0][0] as f64 / norm_factor) as i32,
            (table[0][1] as f64 / norm_factor) as i32,
        ],
        [
            (table[1][0] as f64 / norm_factor) as i32,
            (table[1][1] as f64 / norm_factor) as i32,
        ],
    ]
}

/// `FisherStrand.pValueForContingencyTable`.
pub fn p_value_for_contingency_table(table: [[i32; 2]; 2]) -> f64 {
    fisher_exact::two_sided_p_value(normalize_contingency_table(table))
}

/// `QualityUtils.phredScaleErrorRate`, with the floor `MIN_LOG10_SCALED_QUAL` applies.
fn phred_scale_error_rate(error_rate: f64) -> f64 {
    let log10 = jmath::math::log10(error_rate);
    let floor = jmath::math::log10(f64::MIN_POSITIVE * f64::EPSILON / 2.0); // log10(Double.MIN_VALUE)
                                                                            // `abs` is in the reference, "for edge base with errorRateLog10 = 0 producing -0.0 doubles".
    (-10.0 * log10.max(floor)).abs()
}

/// `StrandOddsRatio.calculateSOR`, over the table with one pseudocount added to every cell.
pub fn calculate_sor(table: [[i32; 2]; 2]) -> f64 {
    let t00 = table[0][0] as f64 + PSEUDOCOUNT;
    let t01 = table[0][1] as f64 + PSEUDOCOUNT;
    let t11 = table[1][1] as f64 + PSEUDOCOUNT;
    let t10 = table[1][0] as f64 + PSEUDOCOUNT;
    let ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);
    let ref_ratio = t00.min(t01) / t00.max(t01);
    let alt_ratio = t10.min(t11) / t10.max(t11);
    // `Math.log`, which is the JDK's and is correctly rounded, not commons-math3's.
    jmath::math::log(ratio) + jmath::math::log(ref_ratio) - jmath::math::log(alt_ratio)
}

/// `String.format("%.3f", value)`, half-up on the decimal expansion.
fn format_three_decimals(value: f64) -> String {
    crate::rank_sum::format_three_decimals(value)
}

/// The shared template: not a variant site is nothing, then the `SB` fields, then the likelihoods.
fn annotate_table<F>(
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    min_count: i32,
    key: &str,
    value_of: F,
) -> Vec<(String, AnnotationValue)>
where
    F: Fn([[i32; 2]; 2]) -> String,
{
    if !vc.is_variant() {
        return Vec::new();
    }
    // The genotype field wins: if any genotype carries SB, the matrix is never consulted.
    if vc.genotypes.iter().any(|g| {
        g.extended
            .iter()
            .any(|(k, _)| k == STRAND_BIAS_BY_SAMPLE_KEY)
    }) {
        return match table_from_samples(vc, min_count) {
            Some(table) => vec![(key.to_string(), AnnotationValue::Str(value_of(table)))],
            // `null` from the template method, which the caller turns into an absent key.
            None => Vec::new(),
        };
    }
    let Some(likelihoods) = likelihoods else {
        return Vec::new();
    };
    let table = contingency_table(likelihoods, vc, min_count);
    vec![(key.to_string(), AnnotationValue::Str(value_of(table)))]
}

/// `FisherStrand`: `FS`, the phred-scaled two-sided Fisher p-value.
pub struct FisherStrand;

impl InfoFieldAnnotation for FisherStrand {
    fn key_names(&self) -> Vec<&'static str> {
        vec![FISHER_STRAND_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        // `MIN_COUNT = ARRAY_DIM`, which is 2.
        annotate_table(vc, likelihoods, 2, FISHER_STRAND_KEY, |table| {
            let p_value = p_value_for_contingency_table(table);
            format_three_decimals(phred_scale_error_rate(p_value.max(MIN_PVALUE)))
        })
    }
}

/// `StrandOddsRatio`: `SOR`.
pub struct StrandOddsRatio;

impl InfoFieldAnnotation for StrandOddsRatio {
    fn key_names(&self) -> Vec<&'static str> {
        vec![STRAND_ODDS_RATIO_KEY]
    }

    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Vec<(String, AnnotationValue)> {
        // `MIN_COUNT = 0`, so a sample with a single read still contributes.
        annotate_table(vc, likelihoods, 0, STRAND_ODDS_RATIO_KEY, |table| {
            format_three_decimals(calculate_sor(table))
        })
    }
}

/// `StrandBiasBySample`: `SB`, the four counts themselves, on the genotype rather than the site.
pub struct StrandBiasBySample;

impl StrandBiasBySample {
    /// `annotate(...)`, which writes nothing when the genotype is not called or there are no
    /// likelihoods, and computes the table over **that sample alone**.
    pub fn counts(
        &self,
        vc: &VariantContext,
        sample: &str,
        likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    ) -> Option<Vec<i32>> {
        let likelihoods = likelihoods?;
        let sample_index = likelihoods.index_of_sample(sample)?;
        let mut table = [[0i32; 2]; 2];
        let reference = vc.alleles.iter().find(|a| a.is_reference())?;
        let alternates: Vec<&Allele> = vc.alleles.iter().filter(|a| !a.is_reference()).collect();
        let mut sample_table = [0i32; 4];
        for best in likelihoods.best_alleles_breaking_ties_for_sample(sample_index, None) {
            if !best.is_informative() {
                continue;
            }
            let Some(read) = likelihoods
                .sample_evidence(sample_index)
                .and_then(|reads| reads.get(best.evidence_index))
            else {
                continue;
            };
            let Some(allele) = best.allele else { continue };
            let matches_ref = allele.display_string() == reference.display_string();
            let matches_alt = alternates.iter().any(|alt| **alt == allele);
            if !matches_ref && !matches_alt {
                continue;
            }
            let offset = if matches_ref { 0 } else { 2 };
            sample_table[offset + usize::from(is_reverse_strand(read))] += 1;
        }
        // The threshold is zero here, so any read at all is counted.
        if passes_minimum_threshold(&sample_table, 0) {
            table[0][0] += sample_table[0];
            table[0][1] += sample_table[1];
            table[1][0] += sample_table[2];
            table[1][1] += sample_table[3];
        }
        Some(vec![table[0][0], table[0][1], table[1][0], table[1][1]])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_normalisation_truncates_and_only_above_four_hundred() {
        assert_eq!(
            normalize_contingency_table([[100, 100], [100, 100]]),
            [[100, 100], [100, 100]]
        );
        let big = normalize_contingency_table([[1000, 1000], [1000, 1000]]);
        let sum: i32 = big.iter().flatten().sum();
        assert!(sum <= 200, "normalised to {sum}");
    }

    #[test]
    fn the_phred_floor_caps_fs_at_about_three_thousand() {
        let capped = phred_scale_error_rate(MIN_PVALUE);
        assert!((3190.0..3210.0).contains(&capped), "capped = {capped}");
    }
}
