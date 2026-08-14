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
//!
//! # The filter analysis counts without the flag that asks for it
//!
//! ```java
//! if (filterAnalysis != null && concordanceState == ConcordanceState.FILTERED_TRUE_NEGATIVE || concordanceState == ConcordanceState.FILTERED_FALSE_NEGATIVE) {
//! ```
//!
//! `&&` binds tighter than `||`, so this is `(flag && FTN) || FFN`: a filtered false negative walks
//! into the map whether or not `--filter-analysis` was given. Nothing is written when it was not, so
//! the only way the difference reaches a user is as a crash. The map is keyed by the **eval
//! header's** FILTER lines and the lookup is unguarded, so a filter the header does not declare is a
//! null, and the increment on it is a `NullPointerException` on a run that asked for no such table.
//! [`FilterAnalysis::apply`] keeps the condition in that shape rather than the one the layout
//! suggests, and returns the null as an error rather than reaching for a record that is not there.
//!
//! # One step writes two different records
//!
//! ```java
//! tryToWrite(truePositivesAndFalseNegativesVcfWriter, annotateWithConcordanceState(truthVersusEval.getTruth(), state));
//! tryToWrite(truePositivesAndFalsePositivesVcfWriter, annotateWithConcordanceState(truthVersusEval.getEval(), state));
//! ```
//!
//! Three of the five states reach two of the three optional VCFs, and never with the same record:
//! [`writes`] is that routing. A true positive is truth's record in `-tpfn` and eval's in `-tpfp`,
//! and a filtered false negative is truth's in `-tpfn` and eval's in `-ftnfn`, labelled `FFN` in
//! both rather than `FN` in the file documented as false negatives. Only `-tpfn` is written against
//! the truth header, so its sample column is the truth file's.

use gatk_engine::base_recalibration_engine::round_to_n_decimal_places;
use gatk_engine::concordance_walker::ConcordanceState;
use gatk_engine::java_hash::{hash_map_order, string_hash_code, HashOrderError};
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

/// `FilterAnalysisTableColumn.COLUMNS`.
pub const FILTER_ANALYSIS_COLUMNS: [&str; 5] = ["filter", "tn", "fn", "uniq_tn", "uniq_fn"];

/// What the filter analysis can fail with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FilterAnalysisError {
    /// `filterAnalysisRecords.get(filter)` answered null and was incremented anyway.
    UndeclaredFilter {
        /// The increment the reference reached for, which is the state's own.
        method: &'static str,
    },
    /// A bucket nothing here has measured, from [`hash_map_order`].
    HashOrder(HashOrderError),
}

impl FilterAnalysisError {
    /// The exception class the reference throws.
    pub fn class(&self) -> &'static str {
        match self {
            FilterAnalysisError::UndeclaredFilter { .. } => "java.lang.NullPointerException",
            FilterAnalysisError::HashOrder(_) => "unmeasured",
        }
    }

    /// The message, which is the JDK's own helpful wording rather than anything GATK writes.
    pub fn message(&self) -> String {
        match self {
            FilterAnalysisError::UndeclaredFilter { method } => format!(
                "Cannot invoke \"org.broadinstitute.hellbender.tools.walkers.validation.FilterAnalysisRecord.{method}()\" because \"record\" is null"
            ),
            FilterAnalysisError::HashOrder(HashOrderError::BucketTreeified { bucket, length }) => {
                format!("bucket {bucket} holds {length} filters, which is past what is measured")
            }
        }
    }
}

/// `FilterAnalysisRecord`: one filter, and the four counters it accumulates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilterAnalysisRecord {
    pub filter: String,
    pub true_negative: i32,
    pub false_negative: i32,
    pub unique_true_negative: i32,
    pub unique_false_negative: i32,
}

/// The `HashMap<String, FilterAnalysisRecord>` of `onTraversalStart`, in declaration order.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FilterAnalysis {
    records: Vec<FilterAnalysisRecord>,
}

impl FilterAnalysis {
    /// One record per FILTER line of the **eval** header, whatever any record carries.
    ///
    /// The truth header's filters are never looked at, and a declared filter nothing carries still
    /// gets a row of zeroes.
    pub fn new(declared_by_the_eval_header: &[String]) -> Self {
        FilterAnalysis {
            records: declared_by_the_eval_header
                .iter()
                .map(|filter| FilterAnalysisRecord {
                    filter: filter.clone(),
                    true_negative: 0,
                    false_negative: 0,
                    unique_true_negative: 0,
                    unique_false_negative: 0,
                })
                .collect(),
        }
    }

    /// The tail of `apply`, condition and all.
    ///
    /// `requested` is `filterAnalysis != null`, and it is deliberately consulted for one of the two
    /// states only: that is what the reference's precedence says.
    pub fn apply(
        &mut self,
        state: ConcordanceState,
        filters_on_the_eval_record: &[String],
        requested: bool,
    ) -> Result<(), FilterAnalysisError> {
        let filtered_true_negative = state == ConcordanceState::FilteredTrueNegative;
        let filtered_false_negative = state == ConcordanceState::FilteredFalseNegative;
        if !((requested && filtered_true_negative) || filtered_false_negative) {
            return Ok(());
        }
        // `filters.size() == 1`: a property of the record, decided once and handed to each of its
        // filters, so a record carrying two of them makes neither of the two unique.
        let unique = filters_on_the_eval_record.len() == 1;
        let method = if filtered_true_negative {
            "incrementTrueNegative"
        } else {
            "incrementFalseNegative"
        };
        for filter in filters_on_the_eval_record {
            let Some(record) = self
                .records
                .iter_mut()
                .find(|record| &record.filter == filter)
            else {
                // The unguarded `filterAnalysisRecords.get(filter)`, on a filter the eval header
                // never declared.
                return Err(FilterAnalysisError::UndeclaredFilter { method });
            };
            if filtered_true_negative {
                record.true_negative += 1;
                if unique {
                    record.unique_true_negative += 1;
                }
            } else {
                record.false_negative += 1;
                if unique {
                    record.unique_false_negative += 1;
                }
            }
        }
        Ok(())
    }

    /// The records in the order `filterAnalysisRecords.values()` hands them over.
    pub fn ordered(&self) -> Result<Vec<&FilterAnalysisRecord>, FilterAnalysisError> {
        let entries: Vec<(String, i32)> = self
            .records
            .iter()
            .map(|record| (record.filter.clone(), string_hash_code(&record.filter)))
            .collect();
        let order = hash_map_order(&entries).map_err(FilterAnalysisError::HashOrder)?;
        Ok(order
            .into_iter()
            .map(|filter| {
                self.records
                    .iter()
                    .find(|record| record.filter == filter)
                    .expect("the order is of these very keys")
            })
            .collect())
    }

    /// The whole filter-analysis table, written only when the flag asked for it.
    pub fn table(&self) -> Result<String, FilterAnalysisError> {
        let rows: Vec<Vec<String>> = self
            .ordered()?
            .into_iter()
            .map(|record| {
                vec![
                    record.filter.clone(),
                    record.true_negative.to_string(),
                    record.false_negative.to_string(),
                    record.unique_true_negative.to_string(),
                    record.unique_false_negative.to_string(),
                ]
            })
            .collect();
        Ok(write_table(&FILTER_ANALYSIS_COLUMNS, &rows, &[]))
    }
}

/// `TRUTH_STATUS_VCF_ATTRIBUTE`.
pub const TRUTH_STATUS_VCF_ATTRIBUTE: &str = "STATUS";

/// `TRUTH_STATUS_HEADER_LINE`, whose description names three of the five states it can hold.
pub const TRUTH_STATUS_HEADER_LINE: &str = "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Truth status: TP/FP/FN for true positive/false positive/false negative.\">";

/// The three optional record outputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotatedVcf {
    /// `-tpfn`, everything truth had.
    TruePositivesAndFalseNegatives,
    /// `-tpfp`, everything eval called.
    TruePositivesAndFalsePositives,
    /// `-ftnfn`, everything eval filtered.
    FilteredTrueNegativesAndFalseNegatives,
}

/// Which of a step's two records is the one written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    Truth,
    Eval,
}

impl AnnotatedVcf {
    /// Which header the file is written against.
    ///
    /// `-tpfn` is the only one built on the **truth** header, so its sample column and its INFO
    /// declarations are the truth file's while the other two are the eval file's.
    pub fn header(&self) -> Side {
        match self {
            AnnotatedVcf::TruePositivesAndFalseNegatives => Side::Truth,
            AnnotatedVcf::TruePositivesAndFalsePositives
            | AnnotatedVcf::FilteredTrueNegativesAndFalseNegatives => Side::Eval,
        }
    }
}

/// `writeTruePositive` and its four siblings: which files a step reaches, and with which record.
///
/// Three of the five states write twice, and never the same record: a true positive is truth's
/// record in `-tpfn` and eval's in `-tpfp`, and a filtered false negative is truth's in `-tpfn` and
/// eval's in `-ftnfn`. Both copies carry the same STATUS, so the two files agree on the label and
/// disagree on everything the two records disagree on.
pub fn writes(state: ConcordanceState) -> Vec<(AnnotatedVcf, Side)> {
    match state {
        ConcordanceState::TruePositive => vec![
            (AnnotatedVcf::TruePositivesAndFalseNegatives, Side::Truth),
            (AnnotatedVcf::TruePositivesAndFalsePositives, Side::Eval),
        ],
        ConcordanceState::FalsePositive => {
            vec![(AnnotatedVcf::TruePositivesAndFalsePositives, Side::Eval)]
        }
        ConcordanceState::FalseNegative => {
            vec![(AnnotatedVcf::TruePositivesAndFalseNegatives, Side::Truth)]
        }
        // Labelled FFN in both, rather than FN in the file documented as false negatives.
        ConcordanceState::FilteredFalseNegative => vec![
            (AnnotatedVcf::TruePositivesAndFalseNegatives, Side::Truth),
            (
                AnnotatedVcf::FilteredTrueNegativesAndFalseNegatives,
                Side::Eval,
            ),
        ],
        ConcordanceState::FilteredTrueNegative => vec![(
            AnnotatedVcf::FilteredTrueNegativesAndFalseNegatives,
            Side::Eval,
        )],
    }
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

    fn declared() -> Vec<String> {
        ["weak", "shallow", "noisy", "unused"]
            .iter()
            .map(|filter| filter.to_string())
            .collect()
    }

    fn filters(list: &[&str]) -> Vec<String> {
        list.iter().map(|filter| filter.to_string()).collect()
    }

    #[test]
    fn a_filtered_false_negative_counts_without_the_flag_that_asks_for_it() {
        let mut analysis = FilterAnalysis::new(&declared());
        // No flag, and a filtered false negative still walks into the map.
        analysis
            .apply(
                ConcordanceState::FilteredFalseNegative,
                &filters(&["weak"]),
                false,
            )
            .expect("weak is declared");
        // No flag, and a filtered true negative does not.
        analysis
            .apply(
                ConcordanceState::FilteredTrueNegative,
                &filters(&["weak"]),
                false,
            )
            .expect("nothing happens at all");
        let weak = analysis
            .records
            .iter()
            .find(|record| record.filter == "weak")
            .expect("declared");
        assert_eq!((weak.false_negative, weak.true_negative), (1, 0));
    }

    #[test]
    fn an_undeclared_filter_is_the_references_null_pointer() {
        let mut analysis = FilterAnalysis::new(&declared());
        // At a truth locus, with no flag: the crash happens on a run that asked for no table.
        let error = analysis
            .apply(
                ConcordanceState::FilteredFalseNegative,
                &filters(&["ghost"]),
                false,
            )
            .expect_err("ghost is not declared");
        assert_eq!(error.class(), "java.lang.NullPointerException");
        assert_eq!(
            error.message(),
            "Cannot invoke \"org.broadinstitute.hellbender.tools.walkers.validation.FilterAnalysisRecord.incrementFalseNegative()\" because \"record\" is null"
        );
        // Standing alone, with no flag: the guard keeps it away from the map entirely.
        analysis
            .apply(
                ConcordanceState::FilteredTrueNegative,
                &filters(&["ghost"]),
                false,
            )
            .expect("the state never reaches the lookup");
        // Standing alone, with the flag: the same null, reached through the other increment.
        let error = analysis
            .apply(
                ConcordanceState::FilteredTrueNegative,
                &filters(&["ghost"]),
                true,
            )
            .expect_err("ghost is not declared");
        assert!(error.message().contains("incrementTrueNegative()"));
    }

    #[test]
    fn a_record_with_two_filters_increments_neither_unique_column() {
        let mut analysis = FilterAnalysis::new(&declared());
        analysis
            .apply(
                ConcordanceState::FilteredFalseNegative,
                &filters(&["weak", "shallow"]),
                true,
            )
            .expect("both declared");
        for name in ["weak", "shallow"] {
            let record = analysis
                .records
                .iter()
                .find(|record| record.filter == name)
                .expect("declared");
            assert_eq!(
                (record.false_negative, record.unique_false_negative),
                (1, 0)
            );
        }
    }

    #[test]
    fn the_row_order_is_a_hash_maps_and_every_declared_filter_has_one() {
        // The golden's baseline: two filtered false negatives, one of them carrying two filters,
        // and two filtered true negatives, one of them carrying two.
        let mut analysis = FilterAnalysis::new(&declared());
        for (state, on_the_record) in [
            (ConcordanceState::FilteredFalseNegative, filters(&["weak"])),
            (
                ConcordanceState::FilteredFalseNegative,
                filters(&["weak", "shallow"]),
            ),
            (ConcordanceState::FilteredTrueNegative, filters(&["weak"])),
            (
                ConcordanceState::FilteredTrueNegative,
                filters(&["shallow", "noisy"]),
            ),
        ] {
            analysis
                .apply(state, &on_the_record, true)
                .expect("all declared");
        }
        assert_eq!(
            analysis.table().expect("four filters in sixteen buckets"),
            "filter\ttn\tfn\tuniq_tn\tuniq_fn\n\
             shallow\t1\t1\t0\t0\n\
             unused\t0\t0\t0\t0\n\
             noisy\t1\t0\t0\t0\n\
             weak\t1\t2\t1\t1\n"
        );
    }

    #[test]
    fn three_of_the_five_states_write_two_different_records() {
        assert_eq!(
            writes(ConcordanceState::TruePositive),
            vec![
                (AnnotatedVcf::TruePositivesAndFalseNegatives, Side::Truth),
                (AnnotatedVcf::TruePositivesAndFalsePositives, Side::Eval),
            ]
        );
        // A filtered false negative reaches the file documented as false negatives, as truth's
        // record, and the filtered file as eval's.
        assert_eq!(
            writes(ConcordanceState::FilteredFalseNegative),
            vec![
                (AnnotatedVcf::TruePositivesAndFalseNegatives, Side::Truth),
                (
                    AnnotatedVcf::FilteredTrueNegativesAndFalseNegatives,
                    Side::Eval
                ),
            ]
        );
        // The two that write once.
        assert_eq!(writes(ConcordanceState::FalseNegative).len(), 1);
        assert_eq!(writes(ConcordanceState::FilteredTrueNegative).len(), 1);
    }

    #[test]
    fn only_the_first_file_is_written_against_the_truth_header() {
        assert_eq!(
            AnnotatedVcf::TruePositivesAndFalseNegatives.header(),
            Side::Truth
        );
        assert_eq!(
            AnnotatedVcf::TruePositivesAndFalsePositives.header(),
            Side::Eval
        );
        assert_eq!(
            AnnotatedVcf::FilteredTrueNegativesAndFalseNegatives.header(),
            Side::Eval
        );
    }

    #[test]
    fn the_status_of_a_filtered_state_is_its_own_abbreviation() {
        // What the STATUS attribute holds, which the header line's description does not list.
        assert_eq!(
            ConcordanceState::FilteredFalseNegative.abbreviation(),
            "FFN"
        );
        assert_eq!(ConcordanceState::FilteredTrueNegative.abbreviation(), "FTN");
        assert!(TRUTH_STATUS_HEADER_LINE.contains("ID=STATUS,Number=1,Type=String"));
    }

    #[test]
    fn the_truth_side_drops_symbolic_records_and_the_eval_side_keeps_them() {
        assert!(truth_variant_filter(false, false));
        assert!(!truth_variant_filter(false, true));
        assert!(!truth_variant_filter(true, false));
    }
}
