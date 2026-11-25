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
//!
//! The other cache the reference has, starting a call at the first column where this haplotype
//! differs from the last one, is switched off in the reference itself (`PDPairHMM.java`):
//!
//! ```java
//! int nextHapStartIndex = 0; //disable the optimization for now until its confirmed to be correct
//! ```
//!
//! so `hapStartIndex` is always zero and every column from one is computed, as here.
//!
//! # Which kernel, and why the arithmetic is written out
//!
//! This is the Java `LoglessPDPairHMM`, not a native one. `PDPairHMM.java` changed in 4.7.0.0:
//! its `FASTEST_AVAILABLE` now tries the native `VectorLoglessPairPDHMM` (`AVX_LOGLESS_CACHING`)
//! and falls back to `LoglessPDPairHMM` on `HardwareFeatureException`.
//! `docs/the-move-to-4.7.0.0-measured.md` counts the PDHMM's work in that move as "arguments (not
//! its kernel: no suite reaches a DRAGEN PDHMM path)"; the `hc-pdhmm` golden (`hc_pdhmm.txt.gz`) is
//! what pins these bits.
//!
//! Every cell is the reference's expression in the reference's order: products summed left to
//! right, `Math.max` taken where it takes it (`f64::max` differs from `Math.max` only when
//! an operand is NaN). Floating-point addition is not associative, and
//! Rust does not fuse a multiply and an add unless asked, so writing the same expression gives
//! the same rounding. The final logarithm is `jmath::math::log10`, Java's `Math.log10`, not the
//! platform's, for the reason `docs/fgkl-and-lorikeet.md` measured on the PairHMM: the host
//! `log10` differs from Java's in the last bits.
//!
//! # Reading this file without reading Rust
//!
//! - `&[u8]` is a borrowed, read-only `byte[]`; `Vec<Vec<f64>>` is `double[][]`.
//! - `Result<f64, String>` is "a likelihood or the reference's exception as text", and `?` after
//!   a call returns that error from this function at once.
//! - `let zeros = || vec![...]` is a lambda that builds a fresh zeroed matrix each time it is
//!   called.
//! - `match state { A | B => ..., C => ... }` is a `switch` whose arms must cover every case; the
//!   reference's `NORMAL` and `INSIDE_DEL` arms differ only in the three copies, so they share one
//!   arm here with an `if` inside.

use crate::haplotype::pd_bytes;
use crate::math_utils::{qual_to_error_prob, qual_to_prob};
use crate::pair_hmm::{
    initial_condition, initial_condition_log10, qual_to_trans_probs, DELETION_TO_DELETION,
    INDEL_TO_MATCH, INSERTION_TO_INSERTION, MATCH_TO_DELETION, MATCH_TO_INSERTION, MATCH_TO_MATCH,
    TRANS_PROB_ARRAY_LENGTH, TRISTATE_CORRECTION,
};

/// `HMMState`: what the column just finished said about the next one.
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

    // `INITIAL_CONDITION / haplotypeBases.length`: free deletions before the read starts, spread
    // over the haplotype's own length, so each haplotype is seeded with its own value.
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
            // `||` short-circuits as Java's does, so a read base outside ACGT is refused only where
            // the haplotype base does not already match it (or is `N`).
            let matching =
                x == y || x == b'N' || y == b'N' || is_base_pd_matching(x, haplotype_pd[j])?;
            priors[i + 1][j + 1] = if matching {
                qual_to_prob(f64::from(qual))
            } else {
                qual_to_error_prob(f64::from(qual)) / TRISTATE_CORRECTION
            };
        }
    }

    // `(haplotypePDBases[j-OFF] & PartiallyDeterminedHaplotype.DEL_END) == ...DEL_END`, with
    // `OFF` one. Rust's `&` binds tighter than `==`, so no parentheses are needed.
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
            // The same in all three of the reference's arms: at a deletion's end the insertion
            // takes the larger of the two paths from the row above.
            insertions[i][j] = if del_end(j) {
                branch_matches[i - 1][j].max(matches[i - 1][j]) * t[MATCH_TO_INSERTION]
                    + branch_insertions[i - 1][j].max(insertions[i - 1][j])
                        * t[INSERTION_TO_INSERTION]
            } else {
                matches[i - 1][j] * t[MATCH_TO_INSERTION]
                    + insertions[i - 1][j] * t[INSERTION_TO_INSERTION]
            };
            // The reference sets `NORMAL` at the end of its `AFTER_DEL` arm, then lets a
            // `DEL_START` and then a `DEL_END` on this column override it, in that order: "Being
            // after a deltion overrides being inside a deletion".
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

    // Only the match and insertion states of the last row are summed, in column order: paths
    // ending in a deletion are left out, as in the reference.
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
