//! `Concordance`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.Concordance` (GATK 4.6.2.0), with
//! `ConcordanceSummaryRecord` beside it.
//!
//! The tool [`gatk_engine::concordance_walker`] was written for. Every step of the traversal is
//! tallied by its state into one of two buckets, and the two buckets are the two rows of a
//! six-column table. This module is the table alone: `--summary`.
//!
//! # An empty callset reports zero rather than NaN
//!
//! ```java
//! public double getSensitivity() { return (double) truePositives / (truePositives + falseNegatives); }
//! ...
//! .set(SENSITIVITY_COLUMN_NAME, record.getSensitivity(), 3)
//! ```
//!
//! Both rates are `0/0` for a bucket nothing landed in, which is NaN, and the three-decimal setter
//! is `MathUtils.roundToNDecimalPlaces`:
//!
//! ```java
//! return Math.round( (in+Math.ulp(in))*mult )/mult;
//! ```
//!
//! `Math.ulp(NaN)` is NaN and `Math.round(NaN)` is **0**, so the column is written `0.0`. The same
//! `0/0` that [`crate::evaluate_info_field_concordance`] prints as `NaN` prints as a rate of zero
//! here, and the two tools differ only in whether the double reaches the rounding.
//!
//! # A filtered eval record alone leaves no trace
//!
//! ```java
//! snpCounts.get(ConcordanceState.FALSE_NEGATIVE).longValue() + snpCounts.get(ConcordanceState.FILTERED_FALSE_NEGATIVE).longValue()
//! ```
//!
//! Five states are counted and four are read: the FN column folds `FILTERED_FALSE_NEGATIVE` into
//! `FALSE_NEGATIVE`, the FP column is `FALSE_POSITIVE` alone, and `FILTERED_TRUE_NEGATIVE` is
//! incremented and then never looked at. A filtered eval record at a truth locus is a miss without
//! being a false positive; the same record with no truth beside it moves neither rate, and the
//! golden's filtered run reports a precision of `1.0` with an unmatched record in the eval file.
//!
//! # Everything that is not a SNP is an indel
//!
//! ```java
//! if (truthVersusEval.getTruthIfPresentElseEval().isSNP()) { snpCounts... } else { indelCounts... }
//! ```
//!
//! One test, and the row it does not choose is labelled `INDEL` whatever the record was: an MNP, a
//! mixed record and a symbolic one all land there. Which record is asked is truth's when truth is
//! present, so only the two eval-only states are stratified by the eval record.
//!
//! # The truth side drops symbolic records and the eval side keeps them
//!
//! `makeTruthVariantFilter` is `!vc.isFiltered() && !vc.isSymbolicOrSV()` against the base class's
//! `vc -> true`, so a symbolic record written into both files survives on one side only and comes
//! out a false positive rather than a true positive.
//!
//! # Agreement needs the same number of alternates but only truth's first
//!
//! ```java
//! (truth.getAlternateAlleles().size() == eval.getAlternateAlleles().size()) &&
//! ((truth.getAlternateAlleles().size() > 0) &&
//!         eval.getAlternateAlleles().contains(truth.getAlternateAllele(0)))
//! ```
//!
//! The size test asks for nothing but the count and the membership test asks for one allele, so
//! truth `A/C,G` against eval `A/G,C` agrees while truth `A/C` against eval `A/C,G` does not.

use gatk_engine::base_recalibration_engine::round_to_n_decimal_places;
use gatk_engine::concordance_walker::ConcordanceState;
use gatk_engine::tsv_table::{java_double_to_string, write_table};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK Concordance";

/// `SUMMARY_TABLE_COLUMN_HEADER`.
pub const COLUMNS: [&str; 6] = ["type", "TP", "FP", "FN", "RECALL", "PRECISION"];

/// The `EnumMap<ConcordanceState, MutableLong>` of one variant type, all five states of it.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Counts {
    pub true_positive: i64,
    pub false_positive: i64,
    pub false_negative: i64,
    /// Counted here and read nowhere.
    pub filtered_true_negative: i64,
    pub filtered_false_negative: i64,
}

impl Counts {
    /// `apply`: one state, one counter.
    pub fn increment(&mut self, state: ConcordanceState) {
        match state {
            ConcordanceState::TruePositive => self.true_positive += 1,
            ConcordanceState::FalsePositive => self.false_positive += 1,
            ConcordanceState::FalseNegative => self.false_negative += 1,
            ConcordanceState::FilteredTrueNegative => self.filtered_true_negative += 1,
            ConcordanceState::FilteredFalseNegative => self.filtered_false_negative += 1,
        }
    }

    /// The FN column: the two negative states that had a truth record, summed.
    pub fn false_negatives(&self) -> i64 {
        self.false_negative + self.filtered_false_negative
    }

    /// `getSensitivity()`, which is `0/0` for an empty bucket.
    pub fn sensitivity(&self) -> f64 {
        self.true_positive as f64 / (self.true_positive + self.false_negatives()) as f64
    }

    /// `getPrecision()`, whose denominator never sees a filtered record.
    pub fn precision(&self) -> f64 {
        self.true_positive as f64 / (self.true_positive + self.false_positive) as f64
    }
}

/// The two buckets, in the order the table writes them.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Summary {
    pub snp: Counts,
    pub indel: Counts,
}

impl Summary {
    /// `apply`, whose one test is `getTruthIfPresentElseEval().isSNP()`.
    pub fn add(&mut self, state: ConcordanceState, is_snp: bool) {
        if is_snp {
            self.snp.increment(state);
        } else {
            self.indel.increment(state);
        }
    }

    /// The whole summary table: the column line and one row per type, whatever the counts are.
    pub fn table(&self) -> String {
        let row = |name: &str, counts: &Counts| {
            vec![
                name.to_string(),
                counts.true_positive.to_string(),
                counts.false_positive.to_string(),
                counts.false_negatives().to_string(),
                rate(counts.sensitivity()),
                rate(counts.precision()),
            ]
        };
        write_table(
            &COLUMNS,
            &[row("SNP", &self.snp), row("INDEL", &self.indel)],
            &[],
        )
    }
}

/// `DataLine.set(column, value, 3)`: rounded to three places, then `Double.toString`.
///
/// The rounding is what turns NaN into zero, so it cannot be skipped for a rate that happens to be
/// exact: `Math.round` saturates NaN to `0L` and the division by 1000 gives `0.0`.
pub fn rate(value: f64) -> String {
    let rounded = round_to_n_decimal_places(value, 3).expect("three places is more than zero");
    java_double_to_string(rounded)
}

/// `makeTruthVariantFilter()`: what the truth side keeps.
///
/// The eval side has no filter of its own, the base class's being `vc -> true`.
pub fn truth_variant_filter(is_filtered: bool, is_symbolic_or_sv: bool) -> bool {
    !is_filtered && !is_symbolic_or_sv
}

/// `areVariantsAtSameLocusConcordant`, as the two allele lists.
pub fn variants_at_same_locus_are_concordant(
    truth_reference: &str,
    truth_alternates: &[String],
    eval_reference: &str,
    eval_alternates: &[String],
) -> bool {
    let same_reference_allele = truth_reference == eval_reference;
    let contains_alternate_allele = truth_alternates.len() == eval_alternates.len()
        && !truth_alternates.is_empty()
        && eval_alternates.contains(&truth_alternates[0]);
    same_reference_allele && contains_alternate_allele
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alleles(list: &[&str]) -> Vec<String> {
        list.iter().map(|allele| allele.to_string()).collect()
    }

    #[test]
    fn an_empty_callset_reports_zero_rather_than_nan() {
        let summary = Summary::default();
        assert!(summary.snp.sensitivity().is_nan());
        assert!(summary.snp.precision().is_nan());
        assert_eq!(
            summary.table(),
            "type\tTP\tFP\tFN\tRECALL\tPRECISION\n\
             SNP\t0\t0\t0\t0.0\t0.0\n\
             INDEL\t0\t0\t0\t0.0\t0.0\n"
        );
    }

    #[test]
    fn a_filtered_eval_record_alone_leaves_no_trace() {
        // The golden's filtered run: one true positive, one filtered false negative at a truth
        // locus, and one filtered true negative of its own.
        let mut summary = Summary::default();
        summary.add(ConcordanceState::FilteredFalseNegative, true);
        summary.add(ConcordanceState::TruePositive, true);
        summary.add(ConcordanceState::FilteredTrueNegative, true);
        assert_eq!(summary.snp.filtered_true_negative, 1);
        assert_eq!(summary.snp.false_negatives(), 1);
        // The unmatched filtered record is in neither column, so the precision is still exactly one.
        assert_eq!(
            summary.table().lines().nth(1).expect("the SNP row"),
            "SNP\t1\t0\t1\t0.5\t1.0"
        );
    }

    #[test]
    fn the_rounding_is_half_up_after_an_ulp_is_added() {
        // Two true positives against one miss and one false alarm: 2/3 on both rates.
        let mut summary = Summary::default();
        summary.add(ConcordanceState::TruePositive, true);
        summary.add(ConcordanceState::TruePositive, true);
        summary.add(ConcordanceState::FalseNegative, true);
        summary.add(ConcordanceState::FalsePositive, true);
        assert_eq!(rate(summary.snp.sensitivity()), "0.667");
        assert_eq!(rate(summary.snp.precision()), "0.667");
        assert_eq!(rate(0.5), "0.5");
        assert_eq!(rate(1.0), "1.0");
        assert_eq!(rate(0.0), "0.0");
    }

    #[test]
    fn everything_that_is_not_a_snp_shares_one_row() {
        let mut summary = Summary::default();
        // An MNP and a symbolic record, neither of which is a SNP.
        summary.add(ConcordanceState::TruePositive, false);
        summary.add(ConcordanceState::FalsePositive, false);
        assert_eq!(summary.indel.true_positive, 1);
        assert_eq!(summary.indel.false_positive, 1);
        assert_eq!(
            summary.table().lines().nth(2).expect("the INDEL row"),
            "INDEL\t1\t1\t0\t1.0\t0.5"
        );
    }

    #[test]
    fn agreement_needs_the_same_count_but_only_truths_first_alternate() {
        // A reordered alternate list still agrees.
        assert!(variants_at_same_locus_are_concordant(
            "A",
            &alleles(&["C", "G"]),
            "A",
            &alleles(&["G", "C"])
        ));
        // One more alternate on the eval side does not.
        assert!(!variants_at_same_locus_are_concordant(
            "A",
            &alleles(&["C"]),
            "A",
            &alleles(&["C", "G"])
        ));
        // Nor does a different reference allele.
        assert!(!variants_at_same_locus_are_concordant(
            "AT",
            &alleles(&["A"]),
            "ACC",
            &alleles(&["A"])
        ));
        // A truth record with no alternate at all agrees with nothing.
        assert!(!variants_at_same_locus_are_concordant("A", &[], "A", &[]));
    }

    #[test]
    fn the_truth_side_drops_symbolic_records_and_the_eval_side_keeps_them() {
        assert!(truth_variant_filter(false, false));
        assert!(!truth_variant_filter(false, true));
        assert!(!truth_variant_filter(true, false));
    }
}
