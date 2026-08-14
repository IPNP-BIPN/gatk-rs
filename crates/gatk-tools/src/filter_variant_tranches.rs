//! `FilterVariantTranches`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.FilterVariantTranches` (GATK 4.6.2.0).
//!
//! Two passes over the same VCF: the first collects the **input** record's own score wherever it
//! overlaps a resource, the second turns those scores into cutoffs and filters against them.
//!
//! # The cutoff is a truncated index into the resource scores
//!
//! ```java
//! Collections.sort(resourceSNPScores, Collections.reverseOrder());
//! ...
//! int snpIndex = (int) ((t / 100.0) * (double) (resourceSNPScores.size() - 1));
//! snpCutoffs.add(resourceSNPScores.get(snpIndex));
//! ```
//!
//! Descending order and a truncating cast, so this is not a quantile by interpolation and not a
//! median by any usual definition: five scores and a tranche of 50 give index 2. Because the order
//! is descending, a **higher** tranche means a **lower** cutoff, which is what makes a higher
//! tranche more sensitive.
//!
//! # Membership is decided by the first cutoff alone
//!
//! ```java
//! private boolean isTrancheFiltered(double score, List<Double> cutoffs) {
//!     return score <= cutoffs.get(0);
//! }
//! ```
//!
//! The tranches were sorted ascending, so `cutoffs.get(0)` belongs to the smallest of them, and it
//! alone decides filtered-or-not. The comparison is `<=`, so a score sitting exactly on the cutoff
//! is filtered. Only after that does [`filter_string_from_score`] walk the cutoffs to name a band.
//!
//! # The score is the input's, and it is taken once
//!
//! The first pass matches an input record against the resources and then reads
//! `variant.getAttribute(infoKey)` — the **input's** value, not the resource's — and `return`s, so a
//! record overlapping two resources contributes one score rather than two.

use gatk_engine::tsv_table::java_double_to_string;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK FilterVariantTranches";

/// `SNPString` and `INDELString`.
pub const SNP_STRING: &str = "SNP";
pub const INDEL_STRING: &str = "INDEL";

/// The default tranches, which maximise the F1 score on whole-genome human data.
pub const DEFAULT_SNP_TRANCHE: f64 = 99.95;
pub const DEFAULT_INDEL_TRANCHE: f64 = 99.4;

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FilterVariantTranchesError {
    /// `validateTranches`, which is a command-line error rather than a user one.
    BadTranches,
    /// An info key the input header does not declare.
    InfoKeyNotInHeader(String),
    /// Nothing in the input carried the key.
    NothingScored(String),
    /// No resource overlapped anything.
    NoOverlap,
    /// SNPs with no SNP resource, or indels with no indel resource.
    NoResourceFor(&'static str),
}

impl FilterVariantTranchesError {
    /// The exception class, which is not the same for all of them.
    pub fn class(&self) -> &'static str {
        match self {
            FilterVariantTranchesError::BadTranches => {
                "org.broadinstitute.barclay.argparser.CommandLineException"
            }
            FilterVariantTranchesError::InfoKeyNotInHeader(_) => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            _ => "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        }
    }

    /// The message. `BadInput` prefixes its own with `Bad input: `.
    pub fn message(&self) -> String {
        match self {
            FilterVariantTranchesError::BadTranches => "At least 1 tranche value must be given and all tranches must be greater than 0 and less than 100.".to_string(),
            FilterVariantTranchesError::InfoKeyNotInHeader(key) => {
                format!("Input VCF does not contain a header line for specified info key:{key}")
            }
            FilterVariantTranchesError::NothingScored(key) => format!(
                "Bad input: VCF contains no variants or no variants with INFO score key \"{key}\""
            ),
            FilterVariantTranchesError::NoOverlap => "Bad input: Neither SNP nor indel resource contains variants overlapping input.  Filtering cannot be performed.".to_string(),
            FilterVariantTranchesError::NoResourceFor(what) => format!(
                "Bad input: {what}s are present in input VCF, but cannot be filtered because no overlapping {what}s were found in the resources."
            ),
        }
    }
}

/// `validateTranches`: at least one, each in `[0, 100)`, deduplicated and sorted.
///
/// The message says "greater than 0" and the test is `d < 0`, so zero itself is allowed.
pub fn validate_tranches(tranches: &[f64]) -> Result<Vec<f64>, FilterVariantTranchesError> {
    if tranches.is_empty() || tranches.iter().any(|value| *value < 0.0 || *value >= 100.0) {
        return Err(FilterVariantTranchesError::BadTranches);
    }
    let mut distinct: Vec<f64> = Vec::new();
    for value in tranches {
        // `stream().distinct()` keeps the first of each, which the sort then reorders anyway.
        if !distinct.contains(value) {
            distinct.push(*value);
        }
    }
    distinct.sort_by(f64::total_cmp);
    Ok(distinct)
}

/// `afterFirstPass`: the cutoffs, or the refusal the counts imply.
///
/// `scores` are the input's own values at resource sites, in the order the traversal found them.
pub fn cutoffs(
    snp_scores: &[f64],
    indel_scores: &[f64],
    scored_snps: usize,
    scored_indels: usize,
    snp_tranches: &[f64],
    indel_tranches: &[f64],
    info_key: &str,
) -> Result<(Vec<f64>, Vec<f64>), FilterVariantTranchesError> {
    if scored_snps == 0 && scored_indels == 0 {
        return Err(FilterVariantTranchesError::NothingScored(
            info_key.to_string(),
        ));
    }
    if snp_scores.is_empty() && indel_scores.is_empty() {
        return Err(FilterVariantTranchesError::NoOverlap);
    }
    if scored_snps > 0 && snp_scores.is_empty() {
        return Err(FilterVariantTranchesError::NoResourceFor(SNP_STRING));
    }
    if scored_indels > 0 && indel_scores.is_empty() {
        return Err(FilterVariantTranchesError::NoResourceFor("indel"));
    }
    Ok((
        cutoffs_of(snp_scores, snp_tranches),
        cutoffs_of(indel_scores, indel_tranches),
    ))
}

/// One class's cutoffs: the scores sorted descending, indexed by a truncated fraction.
pub fn cutoffs_of(scores: &[f64], tranches: &[f64]) -> Vec<f64> {
    if scores.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<f64> = scores.to_vec();
    // `Collections.sort(scores, reverseOrder())`.
    sorted.sort_by(|left, right| right.total_cmp(left));
    tranches
        .iter()
        .map(|tranche| {
            // `(int)` truncates towards zero, and the multiplication is by size - 1.
            let index = ((tranche / 100.0) * (sorted.len() - 1) as f64) as usize;
            sorted[index]
        })
        .collect()
}

/// `isTrancheFiltered`: the smallest tranche's cutoff, inclusive.
pub fn is_tranche_filtered(score: f64, cutoffs: &[f64]) -> bool {
    !cutoffs.is_empty() && score <= cutoffs[0]
}

/// `filterKeyFromTranches`, whose two numbers are `%.2f`.
pub fn filter_key(info_key: &str, class: &str, from: f64, to: f64) -> String {
    format!("{info_key}_{class}_Tranche_{from:.2}_{to:.2}")
}

/// `filterDescriptionFromTranches`.
pub fn filter_description(info_key: &str, class: &str, from: f64, to: f64) -> String {
    format!(
        "{class} truth resource sensitivity between {from:.2} and {to:.2} for info key {info_key}"
    )
}

/// `filterStringFromScore`: which band a filtered score falls in.
///
/// The walk stops at the first cutoff the score is above, and names the band **ending** there. A
/// score below every cutoff falls out of the loop and takes the last tranche's band, which always
/// runs to 100. The `i == 0` case is the reference's own `GATKException`, unreachable because
/// [`is_tranche_filtered`] has already excluded it.
pub fn filter_string_from_score(
    info_key: &str,
    class: &str,
    score: f64,
    tranches: &[f64],
    cutoffs: &[f64],
) -> String {
    for (index, cutoff) in cutoffs.iter().enumerate() {
        if score > *cutoff {
            if index == 0 {
                // "Trying to add a filter to a passing variant."
                return String::new();
            }
            return filter_key(info_key, class, tranches[index - 1], tranches[index]);
        }
    }
    filter_key(info_key, class, tranches[tranches.len() - 1], 100.0)
}

/// `addTrancheHeaderFields`: one line per gap between tranches, plus one to 100.
pub fn tranche_header_lines(
    info_key: &str,
    class: &str,
    tranches: &[f64],
) -> Vec<(String, String)> {
    let mut lines = Vec::new();
    if tranches.len() >= 2 {
        for pair in tranches.windows(2) {
            lines.push((
                filter_key(info_key, class, pair[0], pair[1]),
                filter_description(info_key, class, pair[0], pair[1]),
            ));
        }
    }
    let last = tranches[tranches.len() - 1];
    lines.push((
        filter_key(info_key, class, last, 100.0),
        filter_description(info_key, class, last, 100.0),
    ));
    lines
}

/// The FILTER column of one record of the second pass.
///
/// `previous` are the record's own filters, which `--invalidate-previous-filters` drops. Anything
/// left with no filter at all is written `PASS`, including a record the tool never scored.
pub fn filters_of(
    previous: &[String],
    invalidate_previous: bool,
    added: Option<String>,
) -> Vec<String> {
    let mut filters: Vec<String> = if invalidate_previous {
        Vec::new()
    } else {
        previous.to_vec()
    };
    if let Some(added) = added {
        filters.push(added);
    }
    if filters.is_empty() {
        return vec!["PASS".to_string()];
    }
    filters
}

/// A score as the tool reads it: `Double.parseDouble` of the attribute's own text.
pub fn score(text: &str) -> Option<f64> {
    text.trim().parse::<f64>().ok()
}

/// A score rendered the way this module's messages render one.
pub fn rendered(value: f64) -> String {
    java_double_to_string(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_cutoff_is_a_truncated_index_into_a_descending_list() {
        // The golden's five SNP scores, whatever order they arrived in.
        let scores = [3.0, 5.0, 1.0, 4.0, 2.0];
        // 50% of four is two: the third score down.
        assert_eq!(cutoffs_of(&scores, &[50.0]), vec![3.0]);
        // 99% of four truncates to three, so a higher tranche has a lower cutoff.
        assert_eq!(cutoffs_of(&scores, &[50.0, 99.0]), vec![3.0, 2.0]);
        // The extremes: 0 is the highest score and anything under 100 stops one short of the last.
        assert_eq!(cutoffs_of(&scores, &[0.0]), vec![5.0]);
        assert_eq!(cutoffs_of(&scores, &[99.9]), vec![2.0]);
        // Three indels, half of two is one.
        assert_eq!(cutoffs_of(&[9.0, 8.0, 7.0], &[50.0]), vec![8.0]);
    }

    #[test]
    fn a_score_on_the_cutoff_is_filtered() {
        let cutoffs = [8.0];
        assert!(is_tranche_filtered(8.0, &cutoffs));
        assert!(is_tranche_filtered(7.0, &cutoffs));
        assert!(!is_tranche_filtered(9.0, &cutoffs));
        // With no cutoffs at all nothing is filtered, which is the `size != 0` guard.
        assert!(!is_tranche_filtered(0.0, &[]));
    }

    #[test]
    fn the_band_is_the_first_cutoff_the_score_is_above() {
        let tranches = [50.0, 99.0];
        let cutoffs = [3.0, 2.0];
        // Above the second cutoff: the band between the two tranches.
        assert_eq!(
            filter_string_from_score("SCORE", SNP_STRING, 3.0, &tranches, &cutoffs),
            "SCORE_SNP_Tranche_50.00_99.00"
        );
        // Below every cutoff: the last band, which runs to 100.
        assert_eq!(
            filter_string_from_score("SCORE", SNP_STRING, 2.0, &tranches, &cutoffs),
            "SCORE_SNP_Tranche_99.00_100.00"
        );
        assert_eq!(
            filter_string_from_score("SCORE", SNP_STRING, 1.0, &tranches, &cutoffs),
            "SCORE_SNP_Tranche_99.00_100.00"
        );
        // One tranche: every filtered score lands in the same band.
        assert_eq!(
            filter_string_from_score("SCORE", INDEL_STRING, 7.0, &[50.0], &[8.0]),
            "SCORE_INDEL_Tranche_50.00_100.00"
        );
    }

    #[test]
    fn the_tranches_are_deduplicated_and_sorted_and_bounded() {
        assert_eq!(
            validate_tranches(&[99.0, 50.0, 99.0]).expect("three values, two of them the same"),
            vec![50.0, 99.0]
        );
        // Zero is allowed by the test even though the message says otherwise.
        assert_eq!(validate_tranches(&[0.0]).expect("zero"), vec![0.0]);
        assert_eq!(
            validate_tranches(&[100.0]).expect_err("a hundred").message(),
            "At least 1 tranche value must be given and all tranches must be greater than 0 and less than 100."
        );
        assert_eq!(
            validate_tranches(&[]).expect_err("nothing at all"),
            FilterVariantTranchesError::BadTranches
        );
        assert_eq!(
            validate_tranches(&[100.0]).expect_err("a hundred").class(),
            "org.broadinstitute.barclay.argparser.CommandLineException"
        );
    }

    #[test]
    fn the_four_refusals_of_the_first_pass_each_say_something_different() {
        let none: [f64; 0] = [];
        assert_eq!(
            cutoffs(&none, &none, 0, 0, &[50.0], &[50.0], "SCORE")
                .expect_err("nothing scored")
                .message(),
            "Bad input: VCF contains no variants or no variants with INFO score key \"SCORE\""
        );
        assert_eq!(
            cutoffs(&none, &none, 5, 3, &[50.0], &[50.0], "SCORE")
                .expect_err("no overlap")
                .message(),
            "Bad input: Neither SNP nor indel resource contains variants overlapping input.  Filtering cannot be performed."
        );
        assert_eq!(
            cutoffs(&none, &[9.0], 5, 3, &[50.0], &[50.0], "SCORE")
                .expect_err("indels only")
                .message(),
            "Bad input: SNPs are present in input VCF, but cannot be filtered because no overlapping SNPs were found in the resources."
        );
        assert_eq!(
            cutoffs(&[5.0], &none, 5, 3, &[50.0], &[50.0], "SCORE")
                .expect_err("snps only")
                .message(),
            "Bad input: indels are present in input VCF, but cannot be filtered because no overlapping indels were found in the resources."
        );
    }

    #[test]
    fn a_record_the_tool_did_not_filter_passes() {
        let weak = vec!["weak".to_string()];
        // Nothing added and nothing there: PASS.
        assert_eq!(filters_of(&[], false, None), vec!["PASS"]);
        // A record the tool never scored keeps its own filters rather than passing.
        assert_eq!(filters_of(&weak, false, None), vec!["weak"]);
        // And the new filter joins the old one unless the old ones were invalidated.
        assert_eq!(
            filters_of(
                &weak,
                false,
                Some("SCORE_SNP_Tranche_99.00_100.00".to_string())
            ),
            vec!["weak", "SCORE_SNP_Tranche_99.00_100.00"]
        );
        assert_eq!(
            filters_of(
                &weak,
                true,
                Some("SCORE_SNP_Tranche_99.00_100.00".to_string())
            ),
            vec!["SCORE_SNP_Tranche_99.00_100.00"]
        );
        assert_eq!(filters_of(&weak, true, None), vec!["PASS"]);
    }

    #[test]
    fn the_header_lines_are_one_per_gap_plus_one_to_a_hundred() {
        let lines = tranche_header_lines("SCORE", SNP_STRING, &[50.0, 99.0]);
        assert_eq!(
            lines.iter().map(|(id, _)| id.clone()).collect::<Vec<_>>(),
            vec![
                "SCORE_SNP_Tranche_50.00_99.00",
                "SCORE_SNP_Tranche_99.00_100.00"
            ]
        );
        assert_eq!(
            lines[0].1,
            "SNP truth resource sensitivity between 50.00 and 99.00 for info key SCORE"
        );
        // One tranche is one line.
        assert_eq!(
            tranche_header_lines("SCORE", INDEL_STRING, &[50.0]).len(),
            1
        );
        assert_eq!(rendered(50.0), "50.0");
        assert_eq!(score("3.0"), Some(3.0));
    }
}
