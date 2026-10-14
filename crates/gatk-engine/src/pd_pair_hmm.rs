//! `LoglessPDPairHMM`: the PairHMM over a partially determined haplotype.
//!
//! Ported from `org.broadinstitute.hellbender.utils.pairhmm.LoglessPDPairHMM`, `N2MemoryPDPairHMM`
//! and `PDPairHMM` (GATK 4.6.2.0).
//!
//! The recursion is [`crate::pair_hmm`]'s with two additions, both read off the haplotype's
//! alternate bytes. A base flagged `SNP` also matches the read wherever its flags name the read's
//! base. And a deletion that may or may not be there, flagged `DEL_START` before its first base
//! and `DEL_END` on its last, is skipped as well as walked: a second set of matrices carries the
//! state from before the deletion across it, and after its end each cell takes the larger of the
//! two paths.
//!
//! # What the reference's state machine does that a reader would not expect
//!
//! The state (normal, inside a deletion, just after one) is advanced at the end of each column and
//! is **not reset at the start of a row**: a row that ends inside a deletion starts the next one in
//! that state. The port keeps it that way.
//!
//! # Why every call starts fresh
//!
//! The reference keeps its matrices between calls and only writes the cells a call reads: the
//! deletion row seeded with the initial condition over the haplotype's length (re-seeded only when
//! the length changes, to the same values), every cell of the read's rows over the haplotype's
//! columns, and column zero never, so it stays at zero. Allocating zeroed matrices per call is the
//! same state.

use crate::haplotype::pd_bytes;
use crate::math_utils::{qual_to_error_prob, qual_to_prob};
use crate::pair_hmm::{
    initial_condition, initial_condition_log10, qual_to_trans_probs, DELETION_TO_DELETION,
    INDEL_TO_MATCH, INSERTION_TO_INSERTION, MATCH_TO_DELETION, MATCH_TO_INSERTION, MATCH_TO_MATCH,
    TRANS_PROB_ARRAY_LENGTH, TRISTATE_CORRECTION,
};

/// `HMMState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HmmState {
    Normal,
    InsideDel,
    AfterDel,
}

/// `isBasePDMatching`: whether a base flagged `SNP` may be the read's base. A read base other than
/// the four is the reference's `RuntimeException`.
fn is_base_pd_matching(read_base: u8, pd: u8) -> Result<bool, String> {
    if pd & pd_bytes::SNP == 0 {
        return Ok(false);
    }
    let flag = match read_base {
        b'A' | b'a' => pd_bytes::A,
        b'C' | b'c' => pd_bytes::C,
        b'T' | b't' => pd_bytes::T,
        b'G' | b'g' => pd_bytes::G,
        _ => {
            return Err(
                "java.lang.RuntimeException: Found unexpected base in alt alleles".to_string(),
            )
        }
    };
    Ok(pd & flag != 0)
}

/// `computeReadLikelihoodGivenHaplotypeLog10` for a partially determined haplotype, with the
/// reference's two checks on the answer.
#[allow(clippy::too_many_arguments)]
pub fn pd_read_likelihood_given_haplotype_log10(
    haplotype: &[u8],
    haplotype_pd: &[u8],
    read: &[u8],
    read_quals: &[u8],
    insertion_quals: &[u8],
    deletion_quals: &[u8],
    gap_continuation: &[u8],
) -> Result<f64, String> {
    let padded_read = read.len() + 1;
    let padded_haplotype = haplotype.len() + 1;
    let zeros = || vec![vec![0.0f64; padded_haplotype]; padded_read];
    let (mut matches, mut insertions, mut deletions) = (zeros(), zeros(), zeros());
    let (mut branch_matches, mut branch_insertions, mut branch_deletions) =
        (zeros(), zeros(), zeros());

    let initial = initial_condition() / haplotype.len() as f64;
    for cell in deletions[0].iter_mut() {
        *cell = initial;
    }

    let mut transitions = vec![[0.0f64; TRANS_PROB_ARRAY_LENGTH]; padded_read];
    for index in 0..read.len() {
        transitions[index + 1] = qual_to_trans_probs(
            insertion_quals[index],
            deletion_quals[index],
            gap_continuation[index],
        );
    }

    // `initializePriors`.
    let mut priors = zeros();
    for i in 0..read.len() {
        let x = read[i];
        let qual = read_quals[i];
        for j in 0..haplotype.len() {
            let y = haplotype[j];
            let matching =
                x == y || x == b'N' || y == b'N' || is_base_pd_matching(x, haplotype_pd[j])?;
            priors[i + 1][j + 1] = if matching {
                qual_to_prob(f64::from(qual))
            } else {
                qual_to_error_prob(f64::from(qual)) / TRISTATE_CORRECTION
            };
        }
    }

    let del_end = |j: usize| haplotype_pd[j - 1] & pd_bytes::DEL_END == pd_bytes::DEL_END;
    let mut state = HmmState::Normal;
    for i in 1..padded_read {
        let t = transitions[i];
        for j in 1..padded_haplotype {
            match state {
                HmmState::Normal | HmmState::InsideDel => {
                    if state == HmmState::Normal {
                        branch_matches[i][j] = matches[i][j - 1];
                        branch_deletions[i][j] = deletions[i][j - 1];
                        branch_insertions[i][j] = insertions[i][j - 1];
                    } else {
                        branch_matches[i][j] = branch_matches[i][j - 1];
                        branch_deletions[i][j] = branch_deletions[i][j - 1];
                        branch_insertions[i][j] = branch_insertions[i][j - 1];
                    }
                    matches[i][j] = priors[i][j]
                        * (matches[i - 1][j - 1] * t[MATCH_TO_MATCH]
                            + insertions[i - 1][j - 1] * t[INDEL_TO_MATCH]
                            + deletions[i - 1][j - 1] * t[INDEL_TO_MATCH]);
                    deletions[i][j] = matches[i][j - 1] * t[MATCH_TO_DELETION]
                        + deletions[i][j - 1] * t[DELETION_TO_DELETION];
                }
                HmmState::AfterDel => {
                    branch_matches[i][j] = branch_matches[i][j - 1].max(matches[i][j - 1]);
                    branch_deletions[i][j] = branch_deletions[i][j - 1].max(deletions[i][j - 1]);
                    branch_insertions[i][j] = branch_insertions[i][j - 1].max(insertions[i][j - 1]);
                    matches[i][j] = priors[i][j]
                        * (branch_matches[i - 1][j - 1].max(matches[i - 1][j - 1])
                            * t[MATCH_TO_MATCH]
                            + branch_insertions[i - 1][j - 1].max(insertions[i - 1][j - 1])
                                * t[INDEL_TO_MATCH]
                            + branch_deletions[i - 1][j - 1].max(deletions[i - 1][j - 1])
                                * t[INDEL_TO_MATCH]);
                    deletions[i][j] = branch_matches[i][j - 1].max(matches[i][j - 1])
                        * t[MATCH_TO_DELETION]
                        + branch_deletions[i][j - 1].max(deletions[i][j - 1])
                            * t[DELETION_TO_DELETION];
                }
            }
            insertions[i][j] = if del_end(j) {
                branch_matches[i - 1][j].max(matches[i - 1][j]) * t[MATCH_TO_INSERTION]
                    + branch_insertions[i - 1][j].max(insertions[i - 1][j])
                        * t[INSERTION_TO_INSERTION]
            } else {
                matches[i - 1][j] * t[MATCH_TO_INSERTION]
                    + insertions[i - 1][j] * t[INSERTION_TO_INSERTION]
            };
            if state == HmmState::AfterDel {
                state = HmmState::Normal;
            }
            if haplotype_pd[j - 1] & pd_bytes::DEL_START == pd_bytes::DEL_START {
                state = HmmState::InsideDel;
            }
            if del_end(j) {
                state = HmmState::AfterDel;
            }
        }
    }

    let end = padded_read - 1;
    let mut total = 0.0;
    for j in 1..padded_haplotype {
        total += matches[end][j] + insertions[end][j];
    }
    let result = jmath::math::log10(total) - initial_condition_log10();
    // `Utils.validate(result <= 0.0, ...)`, then `MathUtils.isValidLog10Probability`, which a NaN
    // fails.
    if result > 0.0 {
        return Err(format!(
            "java.lang.IllegalStateException: PairHMM Log Probability cannot be greater than 0: \
             haplotype: {}, read: {}, result: {result:.6}, PairHMM: LoglessPDPairHMM",
            String::from_utf8_lossy(haplotype),
            String::from_utf8_lossy(read),
        ));
    }
    if result.is_nan() {
        return Err("java.lang.IllegalStateException: Invalid Log Probability: NaN".to_string());
    }
    Ok(result)
}
