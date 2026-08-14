//! `EvaluateInfoFieldConcordance`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.EvaluateInfoFieldConcordance`
//! (GATK 4.6.2.0).
//!
//! The first tool written on [`gatk_engine::concordance_walker`]: for every true positive it takes
//! the difference between one INFO field of the eval record and one of the truth record, and reports
//! a mean and a standard deviation for SNPs and for indels.
//!
//! # A record whose key is absent is counted but contributes nothing
//!
//! ```java
//! if (truthVersusEval.getEval().isSNP()) { snpCount++; } else if (...isIndel()) { indelCount++; }
//! this.infoDifference(truthVersusEval.getEval(), truthVersusEval.getTruth());
//! ```
//!
//! `infoDifference` returns immediately unless **both** records carry their key, while the counter
//! was already incremented. The mean is therefore a sum of deltas divided by a count that can be
//! larger than the number of deltas: the golden's baseline has four deltas of 1.0 over five counted
//! true positives and reports `0.8`.
//!
//! # The standard deviation can be NaN by arithmetic alone
//!
//! ```java
//! final double snpVariance = (sumDeltaSquared - sumDelta * sumDelta / count) / count;
//! ```
//!
//! The cancelling form, with no guard. Equal deltas make the two terms equal, and the subtraction
//! of two rounded doubles can land just below zero, whose square root is NaN. The port keeps the
//! expression exactly as written, including the order of the operations, because any algebraically
//! equivalent rearrangement is a different double.
//!
//! # An empty bucket is two NaN columns
//!
//! `count` of zero makes the mean `0/0` and the variance `(0 - 0/0)/0`, so a run with no indel among
//! its true positives still writes an INDEL row and that row is `NaN NaN`. A run where nothing
//! agrees writes two of them.
//!
//! # The mean is of absolute differences, computed the long way
//!
//! `sumDelta += Math.sqrt(deltaSquared)` rather than `Math.abs(delta)`. The two agree on the
//! golden's inputs; they are not the same function of a double in general, since `delta * delta`
//! rounds before `sqrt` does.
//!
//! # Only true positives are looked at
//!
//! The other four states fall through the switch untouched, and this walker drops filtered records
//! on both sides, so the two filtered states cannot occur at all.

use gatk_engine::tsv_table::{java_double_to_string, write_table};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK EvaluateInfoFieldConcordance";

/// `INFO_CONCORDANCE_COLUMN_HEADER`.
pub const COLUMNS: [&str; 5] = [
    "type",
    "eval_info_key",
    "true_info_key",
    "mean_difference",
    "std_difference",
];

/// What this tool refuses, before any record is read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvaluateInfoFieldConcordanceError {
    /// A key the eval header does not declare.
    MissingEvalKey { key: String, file: String },
    /// A key the truth header does not declare.
    MissingTruthKey { key: String, file: String },
}

impl EvaluateInfoFieldConcordanceError {
    /// A plain `UserException`.
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    /// The message, whose colons have no spaces around them.
    pub fn message(&self) -> String {
        match self {
            EvaluateInfoFieldConcordanceError::MissingEvalKey { key, file } => {
                format!("Missing key:{key} in Eval VCF:{file}")
            }
            EvaluateInfoFieldConcordanceError::MissingTruthKey { key, file } => {
                format!("Missing key:{key} in Truth VCF:{file}")
            }
        }
    }
}

/// `onTraversalStart`: both headers are checked, eval first.
pub fn check_keys(
    eval_declares: bool,
    eval_key: &str,
    eval_file: &str,
    truth_declares: bool,
    truth_key: &str,
    truth_file: &str,
) -> Result<(), EvaluateInfoFieldConcordanceError> {
    if !eval_declares {
        return Err(EvaluateInfoFieldConcordanceError::MissingEvalKey {
            key: eval_key.to_string(),
            file: eval_file.to_string(),
        });
    }
    if !truth_declares {
        return Err(EvaluateInfoFieldConcordanceError::MissingTruthKey {
            key: truth_key.to_string(),
            file: truth_file.to_string(),
        });
    }
    Ok(())
}

/// The four counters the traversal keeps, per variant type.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Bucket {
    /// Every true positive of this type, whether or not it carried the key.
    pub count: i64,
    /// The sum of `sqrt(delta * delta)`.
    pub sum_delta: f64,
    /// The sum of `delta * delta`.
    pub sum_delta_squared: f64,
}

impl Bucket {
    /// `snpSumDelta / snpCount`, which is `0/0` for an untouched bucket.
    pub fn mean(&self) -> f64 {
        self.sum_delta / self.count as f64
    }

    /// `(sumSq - sum * sum / n) / n`, in that order, and its square root.
    pub fn standard_deviation(&self) -> f64 {
        let count = self.count as f64;
        let variance = (self.sum_delta_squared - self.sum_delta * self.sum_delta / count) / count;
        variance.sqrt()
    }
}

/// The two buckets, in the order the table writes them.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Concordance {
    pub snp: Bucket,
    pub indel: Bucket,
}

/// Which bucket a true positive falls in, decided by the **eval** record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvalType {
    Snp,
    Indel,
    /// Anything else: counted nowhere, though its delta is still not taken.
    Other,
}

impl Concordance {
    /// `apply` for one true positive: the counter first, the delta only if both keys are there.
    pub fn add(&mut self, eval_type: EvalType, eval_value: Option<f64>, truth_value: Option<f64>) {
        match eval_type {
            EvalType::Snp => self.snp.count += 1,
            EvalType::Indel => self.indel.count += 1,
            EvalType::Other => {}
        }
        let (Some(eval_value), Some(truth_value)) = (eval_value, truth_value) else {
            return;
        };
        let delta = eval_value - truth_value;
        let delta_squared = delta * delta;
        // `Math.sqrt(deltaSquared)` where `Math.abs(delta)` was meant.
        let absolute = delta_squared.sqrt();
        match eval_type {
            EvalType::Snp => {
                self.snp.sum_delta += absolute;
                self.snp.sum_delta_squared += delta_squared;
            }
            EvalType::Indel => {
                self.indel.sum_delta += absolute;
                self.indel.sum_delta_squared += delta_squared;
            }
            // The type test is made twice, and a record that is neither adds to no bucket even
            // though its difference was computed.
            EvalType::Other => {}
        }
    }

    /// The whole summary table: the header line and one row per type.
    pub fn table(&self, eval_key: &str, truth_key: &str) -> String {
        let row = |name: &str, bucket: &Bucket| {
            vec![
                name.to_string(),
                eval_key.to_string(),
                truth_key.to_string(),
                java_double_to_string(bucket.mean()),
                java_double_to_string(bucket.standard_deviation()),
            ]
        };
        write_table(
            &COLUMNS,
            &[row("SNP", &self.snp), row("INDEL", &self.indel)],
            &[],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_record_whose_key_is_absent_is_counted_but_contributes_nothing() {
        // The golden's baseline: four deltas of 1.0 and a fifth true positive with no truth key.
        let mut concordance = Concordance::default();
        for _ in 0..4 {
            concordance.add(EvalType::Snp, Some(2.0), Some(1.0));
        }
        concordance.add(EvalType::Snp, Some(1.0), None);
        assert_eq!(concordance.snp.count, 5);
        assert_eq!(concordance.snp.sum_delta, 4.0);
        assert_eq!(java_double_to_string(concordance.snp.mean()), "0.8");
        assert_eq!(
            java_double_to_string(concordance.snp.standard_deviation()),
            "0.39999999999999997"
        );
    }

    #[test]
    fn an_empty_bucket_is_two_nan_columns() {
        let concordance = Concordance::default();
        assert!(concordance.indel.mean().is_nan());
        assert!(concordance.indel.standard_deviation().is_nan());
        assert_eq!(
            concordance.table("SCORE", "SCORE"),
            "type\teval_info_key\ttrue_info_key\tmean_difference\tstd_difference\n\
             SNP\tSCORE\tSCORE\tNaN\tNaN\n\
             INDEL\tSCORE\tSCORE\tNaN\tNaN\n"
        );
    }

    #[test]
    fn the_goldens_spread_run_is_reproduced() {
        // Deltas of 0.5, -0.75 and 100.0.
        let mut concordance = Concordance::default();
        concordance.add(EvalType::Snp, Some(1.5), Some(1.0));
        concordance.add(EvalType::Snp, Some(0.25), Some(1.0));
        concordance.add(EvalType::Snp, Some(101.0), Some(1.0));
        assert_eq!(java_double_to_string(concordance.snp.mean()), "33.75");
        assert_eq!(
            java_double_to_string(concordance.snp.standard_deviation()),
            "46.84593543378835"
        );
    }

    #[test]
    fn one_delta_has_a_standard_deviation_of_zero() {
        let mut concordance = Concordance::default();
        concordance.add(EvalType::Indel, Some(2.5), Some(1.0));
        assert_eq!(java_double_to_string(concordance.indel.mean()), "1.5");
        assert_eq!(
            java_double_to_string(concordance.indel.standard_deviation()),
            "0.0"
        );
    }

    #[test]
    fn both_refusals_carry_the_references_words() {
        let missing_eval = check_keys(false, "SCORE", "no-key.vcf", true, "SCORE", "truth.vcf")
            .expect_err("the eval header lacks it");
        assert_eq!(
            missing_eval.message(),
            "Missing key:SCORE in Eval VCF:no-key.vcf"
        );
        let missing_truth = check_keys(true, "SCORE", "eval.vcf", false, "SCORE", "no-key.vcf")
            .expect_err("the truth header lacks it");
        assert_eq!(
            missing_truth.message(),
            "Missing key:SCORE in Truth VCF:no-key.vcf"
        );
    }
}
