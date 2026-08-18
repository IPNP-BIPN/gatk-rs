//! `StrandArtifactFilter`'s E step, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.StrandArtifactFilter`
//! (GATK 4.6.2.0).
//!
//! "Do the reads supporting this allele all come from one strand?" Three hypotheses are weighed
//! against each other for every alternate: a forward-strand artifact, a reverse-strand artifact, and
//! neither. The filter reports the first two added together.
//!
//! The M step, which re-estimates the prior and the beta shape between passes by brute-force
//! optimisation, needs a Brent optimiser and is not here.
//!
//! # The strand counts come out of a string
//!
//! `AS_SB_TABLE` is read with `getAttributeAsString` and split on `|` and then on `,`, with brackets
//! stripped and each field trimmed before `Integer.parseInt`. So:
//!
//!  * a non-integer field, an empty field included, is a `NumberFormatException` out of a filter;
//!  * an absent or empty annotation is an empty list;
//!  * and **a table with one entry drops the whole filter**, `sbs.size() <= 1` being a guard on the
//!    filter rather than on one allele.
//!
//! # Two branches answer a hard zero
//!
//! ```java
//! if (altSB.stream().mapToInt(Integer::intValue).sum() == 0 || altIndelSize > LONGEST_STRAND_ARTIFACT_INDEL_SIZE) {
//!     return new EStep(0, 0, totalFwd, totalRev, altSB.get(0), altSB.get(1));
//! }
//! ```
//!
//! An alternate with no reads on either strand, and an indel longer than four, are both an `EStep`
//! with both responsibilities zero. The filter still reports for that allele: zero is an answer, not
//! an abstention.
//!
//! # The normalisation is in log10 and the likelihoods are natural
//!
//! The three log-likelihoods are multiplied by `MathUtils.LOG10_E` on the way into
//! `normalizeLog10`, a conversion by multiplication rather than a change of base inside the
//! logarithm. The port keeps the multiplication where the reference has it.
//!
//! # A prior of one is not a probability, and it does not show
//!
//! At `strandArtifactPrior = 1` the none hypothesis has a log prior of negative infinity, the
//! forward responsibility normalises to exactly `1.0`, and the reverse one to `1.3E-25`. The three
//! sum to more than one in the reals; the two the filter reports sum back to exactly `1.0` in
//! doubles, `1.3E-25` being far below the ulp of one. The golden pins both numbers, so the
//! discrepancy is visible in the responsibilities even though it cannot be in their sum.
//!
//! # What is not modelled, and why
//!
//! `indelSizes` is computed over **every** alternate allele while the strand table has had its
//! symbolic entries removed, and the two are then indexed by the same `i`. A symbolic allele that is
//! not the last alternate therefore pairs each remaining table entry with the wrong allele's indel
//! size. Nothing measured reaches it -- the golden's symbolic allele is last, where the shift cannot
//! show -- so [`calculate_artifact_probabilities`] takes the sizes already paired with the table and
//! refuses to guess at the shifted case.

use crate::beta_binomial::{BetaBinomialDistribution, BetaBinomialError};
use crate::math_utils::{log10_sum_log10, pow10};
use crate::somatic_clustering_model::AlternateAllele;
use jmath::combinatorics::{binomial_coefficient_log, CombinatoricsError};

/// `StrandArtifactFilter`'s identity.
pub const FILTER_NAME: &str = "strand_bias";

/// `phredScaledPosteriorAnnotationName`, `STRAND_QUAL_KEY`.
pub const ANNOTATION: &str = "STRANDQ";

/// `INITIAL_ALPHA_STRAND`, the beta prior on the artifact allele fraction.
pub const INITIAL_ALPHA_STRAND: f64 = 1.0;

/// `INITIAL_BETA_STRAND`.
pub const INITIAL_BETA_STRAND: f64 = 20.0;

/// `INITIAL_STRAND_ARTIFACT_PRIOR`.
pub const INITIAL_STRAND_ARTIFACT_PRIOR: f64 = 0.001;

/// `ALPHA_SEQ`, the beta prior on the sequencing-error allele fraction.
pub const ALPHA_SEQ: f64 = 1.0;

/// `BETA_SEQ_SNV`.
pub const BETA_SEQ_SNV: f64 = 1000.0;

/// `BETA_SEQ_SHORT_INDEL`.
pub const BETA_SEQ_SHORT_INDEL: f64 = 5000.0;

/// `BETA_SEQ_LONG_INDEL`.
pub const BETA_SEQ_LONG_INDEL: f64 = 50000.0;

/// `LONG_INDEL_SIZE`, the size from which the long-indel prior applies.
pub const LONG_INDEL_SIZE: i32 = 3;

/// `LONGEST_STRAND_ARTIFACT_INDEL_SIZE`, past which the filter answers zero.
pub const LONGEST_STRAND_ARTIFACT_INDEL_SIZE: i32 = 4;

/// What this filter refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum StrandArtifactError {
    /// `Integer.parseInt` on a field of `AS_SB_TABLE` that is not one.
    NumberFormat { input: String },
    /// A table entry without both a forward and a reverse count, which the reference reaches by
    /// `get(0)`/`get(1)` and nothing measured produces.
    TableEntryTooShort { index: usize, length: usize },
    /// The beta binomial refused.
    BetaBinomial(BetaBinomialError),
    /// A binomial coefficient the reference could compute and `jmath` has not measured.
    Combinatorics(CombinatoricsError),
}

impl StrandArtifactError {
    pub fn class(&self) -> Option<&'static str> {
        match self {
            StrandArtifactError::NumberFormat { .. } => Some("java.lang.NumberFormatException"),
            _ => None,
        }
    }

    /// `NumberFormatException`'s message, which quotes the input.
    pub fn message(&self) -> Option<String> {
        match self {
            StrandArtifactError::NumberFormat { input } => {
                Some(format!("For input string: \"{input}\""))
            }
            _ => None,
        }
    }
}

impl From<BetaBinomialError> for StrandArtifactError {
    fn from(error: BetaBinomialError) -> Self {
        StrandArtifactError::BetaBinomial(error)
    }
}

impl From<CombinatoricsError> for StrandArtifactError {
    fn from(error: CombinatoricsError) -> Self {
        StrandArtifactError::Combinatorics(error)
    }
}

/// `StrandArtifactFilter.EStep`: two responsibilities and the counts they were computed from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EStep {
    pub forward_artifact_responsibility: f64,
    pub reverse_artifact_responsibility: f64,
    pub forward_count: i32,
    pub reverse_count: i32,
    pub forward_alt_count: i32,
    pub reverse_alt_count: i32,
}

impl EStep {
    /// `getArtifactProbability`, which is what the filter reports.
    pub fn artifact_probability(&self) -> f64 {
        self.forward_artifact_responsibility + self.reverse_artifact_responsibility
    }
}

/// `StrandBiasUtils.getSBsForAlleles`: `AS_SB_TABLE` split on `|` and then on `,`.
///
/// `None` for an absent or empty annotation, which is the reference's empty list.
pub fn parse_strand_bias_table(table: &str) -> Result<Vec<Vec<i32>>, StrandArtifactError> {
    if table.is_empty() {
        return Ok(Vec::new());
    }
    // `replaceAll(BRACKET_REGEX, "")` strips every bracket anywhere, not only at the ends.
    let stripped: String = table.chars().filter(|c| *c != '[' && *c != ']').collect();
    // `splitByWholeSeparatorPreserveAllTokens`, so an empty field is kept and then refused.
    let mut rows = Vec::new();
    for entry in stripped.split('|') {
        let mut counts = Vec::new();
        for field in entry.split(',') {
            let field = field.trim();
            counts.push(
                field
                    .parse::<i32>()
                    .map_err(|_| StrandArtifactError::NumberFormat {
                        input: field.to_string(),
                    })?,
            );
        }
        rows.push(counts);
    }
    Ok(rows)
}

/// `Math.abs(vc.getReference().length() - alt.length())`.
pub fn indel_size(reference_length: i32, alternate: AlternateAllele) -> i32 {
    (reference_length - alternate.length).abs()
}

/// `calculateArtifactProbabilities`, one `EStep` per alternate allele.
///
/// `strand_bias_table` is what [`parse_strand_bias_table`] produced, **with the symbolic entries
/// already removed**, and `indel_sizes` is one size per remaining entry: see the module note on why
/// the pairing is the caller's to make.
pub fn calculate_artifact_probabilities(
    strand_bias_table: &[Vec<i32>],
    indel_sizes: &[i32],
    strand_artifact_prior: f64,
    alpha_strand: f64,
    beta_strand: f64,
) -> Result<Vec<EStep>, StrandArtifactError> {
    // `sbs == null || sbs.isEmpty() || sbs.size() <= 1`: one entry drops the whole filter.
    if strand_bias_table.len() <= 1 {
        return Ok(Vec::new());
    }
    let mut total_forward = 0;
    let mut total_reverse = 0;
    for (index, entry) in strand_bias_table.iter().enumerate() {
        total_forward += *at(entry, 0, index)?;
        total_reverse += *at(entry, 1, index)?;
    }

    let mut steps = Vec::with_capacity(strand_bias_table.len() - 1);
    for (index, entry) in strand_bias_table[1..].iter().enumerate() {
        let forward_alt = *at(entry, 0, index + 1)?;
        let reverse_alt = *at(entry, 1, index + 1)?;
        let size = indel_sizes[index];
        if forward_alt + reverse_alt == 0 || size > LONGEST_STRAND_ARTIFACT_INDEL_SIZE {
            steps.push(EStep {
                forward_artifact_responsibility: 0.0,
                reverse_artifact_responsibility: 0.0,
                forward_count: total_forward,
                reverse_count: total_reverse,
                forward_alt_count: forward_alt,
                reverse_alt_count: reverse_alt,
            });
        } else {
            steps.push(strand_artifact_probability(
                strand_artifact_prior,
                total_forward,
                total_reverse,
                forward_alt,
                reverse_alt,
                size,
                alpha_strand,
                beta_strand,
            )?);
        }
    }
    Ok(steps)
}

/// `calculateErrorProbabilityForAlleles`: the two responsibilities added, or an empty list.
pub fn error_probabilities(steps: &[EStep]) -> Vec<f64> {
    steps.iter().map(EStep::artifact_probability).collect()
}

/// `strandArtifactProbability`, the three-hypothesis comparison.
#[allow(clippy::too_many_arguments)]
pub fn strand_artifact_probability(
    strand_artifact_prior: f64,
    forward_count: i32,
    reverse_count: i32,
    forward_alt_count: i32,
    reverse_alt_count: i32,
    indel_size: i32,
    alpha_strand: f64,
    beta_strand: f64,
) -> Result<EStep, StrandArtifactError> {
    let forward_log_likelihood =
        artifact_strand_log_likelihood(
            forward_count,
            forward_alt_count,
            alpha_strand,
            beta_strand,
        )? + non_artifact_strand_log_likelihood(reverse_count, reverse_alt_count, indel_size)?;
    let reverse_log_likelihood =
        artifact_strand_log_likelihood(
            reverse_count,
            reverse_alt_count,
            alpha_strand,
            beta_strand,
        )? + non_artifact_strand_log_likelihood(forward_count, forward_alt_count, indel_size)?;
    // Three binomial coefficients and a beta binomial, in one expression.
    let none_log_likelihood =
        binomial_coefficient_log(i64::from(forward_count), i64::from(forward_alt_count))?
            + binomial_coefficient_log(i64::from(reverse_count), i64::from(reverse_alt_count))?
            - binomial_coefficient_log(
                i64::from(forward_count) + i64::from(reverse_count),
                i64::from(forward_alt_count) + i64::from(reverse_alt_count),
            )?
            + BetaBinomialDistribution::new(1.0, 1.0, forward_count + reverse_count)?
                .log_probability(forward_alt_count + reverse_alt_count)?;

    // `Math.log`, and the artifact prior is halved between the two strands.
    let forward_log_prior = jmath::math::log(strand_artifact_prior / 2.0);
    let reverse_log_prior = jmath::math::log(strand_artifact_prior / 2.0);
    let none_log_prior = jmath::math::log(1.0 - strand_artifact_prior);

    let log10_e = log10_e();
    let unnormalized = [
        (forward_log_likelihood + forward_log_prior) * log10_e,
        (reverse_log_likelihood + reverse_log_prior) * log10_e,
        (none_log_likelihood + none_log_prior) * log10_e,
    ];
    let probabilities = normalize_log10(&unnormalized);

    Ok(EStep {
        forward_artifact_responsibility: probabilities[0],
        reverse_artifact_responsibility: probabilities[1],
        forward_count,
        reverse_count,
        forward_alt_count,
        reverse_alt_count,
    })
}

/// `MathUtils.LOG10_E`, which is `Math.log10(Math.E)`.
fn log10_e() -> f64 {
    jmath::math::log10(std::f64::consts::E)
}

/// `MathUtils.normalizeLog10(array, false, true)`: subtract the log10 sum, then raise ten to it.
fn normalize_log10(values: &[f64; 3]) -> [f64; 3] {
    let log10_sum = log10_sum_log10(values);
    [
        pow10(values[0] - log10_sum),
        pow10(values[1] - log10_sum),
        pow10(values[2] - log10_sum),
    ]
}

/// `artifactStrandLogLikelihood(strandCount, strandAltCount, alpha, beta)`.
fn artifact_strand_log_likelihood(
    strand_count: i32,
    strand_alt_count: i32,
    alpha: f64,
    beta: f64,
) -> Result<f64, StrandArtifactError> {
    Ok(BetaBinomialDistribution::new(alpha, beta, strand_count)?
        .log_probability(strand_alt_count)?)
}

/// `nonArtifactStrandLogLikelihood`, whose beta takes three steps in the indel size.
fn non_artifact_strand_log_likelihood(
    strand_count: i32,
    strand_alt_count: i32,
    indel_size: i32,
) -> Result<f64, StrandArtifactError> {
    let beta_seq = if indel_size == 0 {
        BETA_SEQ_SNV
    } else if indel_size < LONG_INDEL_SIZE {
        BETA_SEQ_SHORT_INDEL
    } else {
        BETA_SEQ_LONG_INDEL
    };
    Ok(
        BetaBinomialDistribution::new(ALPHA_SEQ, beta_seq, strand_count)?
            .log_probability(strand_alt_count)?,
    )
}

/// One count of a table entry, refusing where the reference would be out of bounds.
fn at(entry: &[i32], position: usize, index: usize) -> Result<&i32, StrandArtifactError> {
    entry
        .get(position)
        .ok_or(StrandArtifactError::TableEntryTooShort {
            index,
            length: entry.len(),
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The annotation is a string, and every shape of it that is not a table is a refusal or an
    /// empty list rather than a zero.
    #[test]
    fn the_table_is_parsed_and_not_read() {
        assert_eq!(
            parse_strand_bias_table("50,50|20,0").expect("parsed"),
            vec![vec![50, 50], vec![20, 0]]
        );
        // Brackets anywhere, and whitespace around every field.
        assert_eq!(
            parse_strand_bias_table("[50, 50 | 20, 0]").expect("parsed"),
            vec![vec![50, 50], vec![20, 0]]
        );
        // An empty annotation is an empty list, not an entry of nothing.
        assert_eq!(
            parse_strand_bias_table("").expect("parsed"),
            Vec::<Vec<i32>>::new()
        );
        // An empty FIELD is a refusal, because `Integer.parseInt("")` is one.
        assert_eq!(
            parse_strand_bias_table("50,50|,0"),
            Err(StrandArtifactError::NumberFormat {
                input: String::new()
            })
        );
        assert_eq!(
            parse_strand_bias_table("50,50|twenty,0"),
            Err(StrandArtifactError::NumberFormat {
                input: "twenty".to_string()
            })
        );
    }

    /// A table with one entry drops the whole filter, and the empty list is what
    /// `ErrorProbabilities` discards.
    #[test]
    fn a_one_entry_table_drops_the_filter() {
        let table = parse_strand_bias_table("50,50").expect("parsed");
        let steps = calculate_artifact_probabilities(
            &table,
            &[0],
            INITIAL_STRAND_ARTIFACT_PRIOR,
            INITIAL_ALPHA_STRAND,
            INITIAL_BETA_STRAND,
        )
        .expect("answered");
        assert!(steps.is_empty());
        assert!(error_probabilities(&steps).is_empty());
    }

    /// The two branches that answer zero rather than abstaining.
    #[test]
    fn no_reads_and_a_long_indel_both_answer_a_hard_zero() {
        let table = parse_strand_bias_table("50,50|0,0").expect("parsed");
        let steps = calculate_artifact_probabilities(
            &table,
            &[0],
            INITIAL_STRAND_ARTIFACT_PRIOR,
            INITIAL_ALPHA_STRAND,
            INITIAL_BETA_STRAND,
        )
        .expect("answered");
        assert_eq!(
            error_probabilities(&steps),
            vec![0.0],
            "no reads on either strand"
        );

        // Twenty forward reads, and an indel one longer than the filter considers.
        let table = parse_strand_bias_table("50,50|20,0").expect("parsed");
        let long = calculate_artifact_probabilities(
            &table,
            &[LONGEST_STRAND_ARTIFACT_INDEL_SIZE + 1],
            INITIAL_STRAND_ARTIFACT_PRIOR,
            INITIAL_ALPHA_STRAND,
            INITIAL_BETA_STRAND,
        )
        .expect("answered");
        assert_eq!(
            error_probabilities(&long),
            vec![0.0],
            "one past the longest"
        );
        // And at the longest itself it is a probability, not a zero.
        let at_the_edge = calculate_artifact_probabilities(
            &table,
            &[LONGEST_STRAND_ARTIFACT_INDEL_SIZE],
            INITIAL_STRAND_ARTIFACT_PRIOR,
            INITIAL_ALPHA_STRAND,
            INITIAL_BETA_STRAND,
        )
        .expect("answered");
        assert!(error_probabilities(&at_the_edge)[0] > 0.9);
    }

    /// A prior of one leaves the two responsibilities summing to more than one, and a prior of zero
    /// makes both artifact hypotheses impossible.
    #[test]
    fn the_priors_ends_are_answers() {
        let one = strand_artifact_probability(
            1.0,
            50,
            50,
            20,
            0,
            0,
            INITIAL_ALPHA_STRAND,
            INITIAL_BETA_STRAND,
        )
        .expect("answered");
        assert_eq!(one.forward_artifact_responsibility, 1.0);
        assert!(one.reverse_artifact_responsibility > 0.0);
        // The three normalised responsibilities sum to more than one in the reals, and the two the
        // filter reports sum back to exactly 1.0 in doubles: 1.3e-25 is far below the ulp of one.
        assert_eq!(one.artifact_probability(), 1.0);

        let zero = strand_artifact_probability(
            0.0,
            50,
            50,
            20,
            0,
            0,
            INITIAL_ALPHA_STRAND,
            INITIAL_BETA_STRAND,
        )
        .expect("answered");
        assert_eq!(zero.artifact_probability(), 0.0);
    }
}
