//! `StrandBiasTest`, its two members `FisherStrand` and `StrandOddsRatio`, and
//! `StrandBiasBySample`, ported from `org.broadinstitute.hellbender.tools.walkers.annotator`
//! (GATK 4.7.0.0; the four Java files are unchanged since 4.6.2.0). `StrandBiasBySample` is a
//! `GenotypeAnnotation` rather than a `StrandBiasTest`, and is here because it writes the table
//! the other two read back.
//!
//! `FS`, `SOR` and `SB`: the same 2x2 contingency table of forward and reverse reads on the
//! reference and the alternate, read three ways.
//!
//! # The table is built per sample and only then added up
//!
//! `StrandBiasTest.java:194-202`, with the stream that fills `sampleTable` elided:
//!
//! ```java
//! for (final String sample : samples) {
//!     final int[] sampleTable = new int[ARRAY_SIZE];
//!     ...
//!     if (passesMinimumThreshold(sampleTable, minCount)) {
//!         copyToMainTable(sampleTable, table);
//!     }
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
//! `FisherStrand.java:106-110`:
//!
//! ```java
//! final double normFactor = sum / TARGET_TABLE_SIZE;
//!
//! return new int[][]{
//!         {(int) (table[0][0] / normFactor), (int) (table[0][1] / normFactor)},
//!         {(int) (table[1][0] / normFactor), (int) (table[1][1] / normFactor)},
//! ```
//!
//! Above 400 total reads the counts are scaled to about 200 and **truncated** to int, so `FS` on a
//! deep site is computed from a table that no longer sums to the coverage. `SOR` does not
//! normalise at all.
//!
//! # Both values are strings, and `FS` is phred-scaled off a floor
//!
//! `FisherStrand.java:85` (one line in the source, wrapped here):
//!
//! ```java
//! return String.format("%.3f", QualityUtils.phredScaleErrorRate(Math.max(pValue, MIN_PVALUE)));
//!         // prevent INFINITYs
//! ```
//!
//! `MIN_PVALUE = 1e-320` is there, in the reference's trailing comment's words, to "prevent
//! INFINITYs", so the largest `FS` a site can report is about 3200 whatever the evidence.
//!
//! # Reading this file without reading Rust
//!
//! - `[[i32; 2]; 2]` is a fixed-size two-by-two array of `int`, Java's `int[2][2]`, and
//!   `[i32; 4]` is `int[4]`. They are values: passing one copies it.
//! - `Option<T>` is a value or nothing; `None` stands for the reference's `null` here.
//! - `let Some(x) = expr else { continue; };` unwraps a present value or runs the `else` branch.
//! - `?` after an `Option` returns `None` from the whole function when the value is absent.
//! - `&VariantContext` is a borrowed, read-only reference; the caller keeps the object.
//! - `x as f64` and `x as i32` are Java's casts; `as i32` from a double truncates toward zero,
//!   as Java's `(int)` does.
//! - `|table| ...` is a lambda; `F: Fn(..) -> String` says a parameter is one.
//! - `match` is a `switch` the compiler checks covers every case.

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

/// `StrandBiasTest.getContingencyTable(likelihoods, vc, minCount)`, per sample and then summed,
/// with `updateTable` and `copyToMainTable` inlined.
///
/// The rows are reference and alternate; the columns are forward and reverse. A variant with no
/// reference allele gives the all-zero table here.
///
/// `updateTable` matches the reference allele as `StrandBiasTest.java:222` (one line in the
/// source, wrapped here):
///
/// ```java
/// final boolean matchesRef = (ref instanceof  InverseAllele) ? (ref).equals(allele)
///         : allele.equals(ref, true);
/// ```
///
/// This port has only the second branch, compared through the display strings.
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

    // `likelihoods.samples()` in matrix order; `0..n` is the half-open range 0, 1, ... n - 1.
    for sample_index in 0..likelihoods.number_of_samples() {
        let mut sample_table = [0i32; 4];
        // `bestAllelesBreakingTies(sample).stream().filter(ba -> ba.isInformative())`.
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
            // `table[offset + (isForward ? 0 : 1)]++`: `usize::from(bool)` is 1 for true and
            // 0 for false.
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
/// `None` is the reference's null: no genotype carried `SB`. The reference sets `foundData`
/// before parsing; here `found` is set only once `strand_counts` succeeds, so a genotype whose
/// `SB` this port cannot parse is skipped where the reference would have thrown or counted it.
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
///
/// Every shape is flattened to tokens and the first four parsed. On shapes no writer here
/// produces this answers `None` where the reference does otherwise: a list of fewer than four
/// Integers is zero-padded by `encodeSBBS(List)`, and a bare Integer is a `GATKException`
/// ("Unexpected SB type").
fn strand_counts(genotype: &htsjdk_vcf::variant::Genotype) -> Option<[i32; 4]> {
    // `g.getAnyAttribute(SB)`; the trailing `?` returns `None` when the genotype has none.
    let value = genotype
        .extended
        .iter()
        .find(|(key, _)| key == STRAND_BIAS_BY_SAMPLE_KEY)
        .map(|(_, value)| value)?;
    // A field read back from a file can also arrive as a list holding the one comma-separated
    // string the codec kept, which is the reference's String under another shape: every piece is
    // split on the comma, so the three shapes give the same four tokens.
    let parts: Vec<String> = match value {
        Value::Str(text) => text.split(',').map(|part| part.to_string()).collect(),
        Value::Int(number) => vec![number.to_string()],
        Value::List(values) => values
            .iter()
            .flat_map(|value| {
                render_value(value)
                    .split(',')
                    .map(str::to_string)
                    .collect::<Vec<_>>()
            })
            .collect(),
        other => vec![render_value(other)],
    };
    let mut counts = [0i32; 4];
    for (index, slot) in counts.iter_mut().enumerate() {
        // `Integer.parseInt(tokenizer.nextToken())`: for the String shape, fewer than four
        // tokens is a NoSuchElementException, which nothing here can produce from a written
        // record. `*slot = ..` writes through the mutable reference `iter_mut` handed out.
        *slot = parts.get(index)?.trim().parse().ok()?;
    }
    Some(counts)
}

/// A parsed field value back as text, so that every shape can be split on the comma; the
/// reference's `o.toString()` on a list element.
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

/// `FisherStrand.normalizeContingencyTable`, which truncates rather than rounds (see the module
/// note for the Java).
///
/// The reference's sum is `int` through `addExact`, which throws `ArithmeticException` past
/// 2^31 - 1 (`FisherStrand.java:101`); here it is `i64` and does not, so the two differ only on a
/// table of more than two billion reads.
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

/// `QualityUtils.phredScaleErrorRate`, with the floor `MIN_LOG10_SCALED_QUAL` it applies.
///
/// `QualityUtils.java:322`, after `phredScaleErrorRate` takes the base-ten logarithm of
/// `errorRate`:
///
/// ```java
/// return Math.abs(-10.0 * Math.max(errorRateLog10, MIN_LOG10_SCALED_QUAL));
/// ```
///
/// `MIN_LOG10_SCALED_QUAL` is the base-ten logarithm of `Double.MIN_VALUE`
/// (`QualityUtils.java:25`), and `Double.MIN_VALUE` is the smallest subnormal, 2^-1074. The
/// `floor` below is meant to be that (its trailing comment says so), but
/// `MIN_POSITIVE * EPSILON` is already 2^-1074 and the further `/ 2.0` rounds to zero, so `floor`
/// is negative infinity rather than about -323.3. Both callers pass `p_value.max(MIN_PVALUE)`,
/// whose logarithm is about -320 or more, so neither
/// floor is ever the larger and the output does not change. The reference's
/// `Utils.validateArg(MathUtils.isValidLog10Probability(errorRateLog10), ..)` at line 319 is not
/// ported: an error rate above one would throw there and is answered here.
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
    // `Math.log` (`StrandOddsRatio.java:161`), the JDK's and not commons-math3's `FastMath`.
    // `jmath::math::log` rounds the true logarithm once, because `Math.log` was measured to be
    // correctly rounded (see that function).
    jmath::math::log(ratio) + jmath::math::log(ref_ratio) - jmath::math::log(alt_ratio)
}

/// `String.format("%.3f", value)`, half-up on the digits `Double.toString` would produce (see
/// `rank_sum::format_three_decimals`).
fn format_three_decimals(value: f64) -> String {
    crate::rank_sum::format_three_decimals(value)
}

/// The shared template: not a variant site is nothing, then the `SB` fields, then the likelihoods.
///
/// `StrandBiasTest.annotate`, with the two abstract `calculateAnnotationFrom...` methods passed in
/// as `min_count` and `value_of`: `F` is any lambda from a table to the formatted value, which is
/// how each member supplies its override.
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
            // `null` from the template method, which the caller turns into an absent key. Some
            // genotype carries `SB` here, so the reference never returns null at this point;
            // this port does when no `SB` could be parsed (see `table_from_samples`).
            None => Vec::new(),
        };
    }
    let Some(likelihoods) = likelihoods else {
        return Vec::new();
    };
    let table = contingency_table(likelihoods, vc, min_count);
    vec![(key.to_string(), AnnotationValue::Str(value_of(table)))]
}

/// `FisherStrand.makeValueObjectForAnnotation(pValueForContingencyTable(table))`: `FS` as written.
pub fn fisher_strand_value(table: [[i32; 2]; 2]) -> String {
    let p_value = p_value_for_contingency_table(table);
    format_three_decimals(phred_scale_error_rate(p_value.max(MIN_PVALUE)))
}

/// `StrandOddsRatio.formattedValue(calculateSOR(table))`: `SOR` as written.
pub fn strand_odds_ratio_value(table: [[i32; 2]; 2]) -> String {
    format_three_decimals(calculate_sor(table))
}

/// `StrandBiasTest.decodeSBBS`: four counts as the two-by-two table, forward then reverse.
///
/// The reference returns null for an array whose length is not four (`StrandBiasTest.java:281`);
/// this reads a missing count as zero and ignores any past the fourth. `.get(index)` is a bounds
/// checked read that gives `None` past the end, and `.copied()` turns the borrowed count into a
/// value.
pub fn decode_sbbs(counts: &[i32]) -> [[i32; 2]; 2] {
    let at = |index: usize| counts.get(index).copied().unwrap_or(0);
    [[at(0), at(1)], [at(2), at(3)]]
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
    ///
    /// Only the likelihoods half of that guard is here (`likelihoods?`, and `?` again for a
    /// sample or reference allele the matrix and variant lack). The not-called half,
    /// `StrandBiasBySample.java:77`, is the caller's: `variant_annotator_engine.rs` in
    /// `gatk-tools` skips this annotation for a genotype that is not called. The body is
    /// `contingency_table` for one sample with `minCount` 0, as the reference's
    /// `FisherStrand.getContingencyTable(likelihoods, vc, 0, Arrays.asList(g.getSampleName()))`.
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
