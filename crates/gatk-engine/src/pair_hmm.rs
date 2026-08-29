//! `LoglessPairHMM`: the likelihood of a read given a haplotype, in linear space.
//!
//! The model is a three-state HMM over the read and the haplotype: match, insertion and deletion.
//! What makes it "logless" is that the recursion runs in LINEAR probability space and is kept from
//! underflowing by an initial condition of `2^1020`, subtracted back off at the end. A log-space
//! version exists in the reference and gives slightly different bits, which is why the golden
//! records both and this ports the one the callers use.
//!
//! Ported from `org.broadinstitute.hellbender.utils.pairhmm.LoglessPairHMM`,
//! `org.broadinstitute.hellbender.utils.pairhmm.PairHMMModel` and the two helpers they reach:
//! `QualityUtils.qualToProb`, `QualityUtils.qualToErrorProb` and
//! `MathUtils.approximateLog10SumLog10`.

use crate::math_utils::{qual_to_error_prob, qual_to_prob};

/// `LoglessPairHMM.INITIAL_CONDITION`, which is `2^1020`: large enough that the recursion cannot
/// underflow a double, and a power of two so that dividing it out costs no precision.
pub fn initial_condition() -> f64 {
    2f64.powi(1020)
}
/// Its base-ten logarithm, subtracted from the answer. Computed rather than written down: a
/// transcribed decimal is a different double from `Math.log10(Math.pow(2, 1020))`, and the
/// difference lands in the seventh digit of every likelihood.
pub fn initial_condition_log10() -> f64 {
    initial_condition().log10()
}
/// `TRISTATE_CORRECTION`: a mismatch may be any of the three other bases, so the error probability
/// is divided by three.
pub const TRISTATE_CORRECTION: f64 = 3.0;
/// `QualityUtils.MAX_QUAL`.
pub const MAX_QUAL: usize = 254;

/// The six transitions the model carries, in the reference's own order.
pub const MATCH_TO_MATCH: usize = 0;
pub const INDEL_TO_MATCH: usize = 1;
pub const MATCH_TO_INSERTION: usize = 2;
pub const INSERTION_TO_INSERTION: usize = 3;
pub const MATCH_TO_DELETION: usize = 4;
pub const DELETION_TO_DELETION: usize = 5;
pub const TRANS_PROB_ARRAY_LENGTH: usize = 6;

/// `MathUtils.JacobianLogTable`, which is a quantised correction rather than a computed one.
///
/// The table holds `log10(1 + 10^-k*step)` for `k` from zero to eighty thousand, and a difference
/// beyond `MAX_TOLERANCE` contributes nothing at all. Computing the correction exactly instead
/// would be a different function, and the two disagree in the last bits.
pub const JACOBIAN_MAX_TOLERANCE: f64 = 8.0;
pub const JACOBIAN_TABLE_STEP: f64 = 0.0001;

/// `MathUtils.fastRound`, which rounds away from zero rather than to even.
pub fn fast_round(value: f64) -> i32 {
    if value > 0.0 {
        (value + 0.5) as i32
    } else {
        (value - 0.5) as i32
    }
}

/// `JacobianLogTable.get`, computed at the table's own quantisation.
///
/// The reference builds the table once and reads it by index; reading it by index is what this
/// reproduces, so the value is the table's entry and not `log10(1 + 10^-difference)`.
pub fn jacobian_log(difference: f64) -> f64 {
    let index = fast_round(difference / JACOBIAN_TABLE_STEP);
    let quantised = f64::from(index) * JACOBIAN_TABLE_STEP;
    (1.0 + 10f64.powf(-quantised)).log10()
}

/// `MathUtils.approximateLog10SumLog10`.
pub fn approximate_log10_sum_log10(a: f64, b: f64) -> f64 {
    if a > b {
        return approximate_log10_sum_log10(b, a);
    }
    if a == f64::NEG_INFINITY {
        return b;
    }
    let difference = b - a;
    b + if difference < JACOBIAN_MAX_TOLERANCE {
        jacobian_log(difference)
    } else {
        0.0
    }
}

/// `PairHMMModel.matchToMatchProb`, which is one minus the probability that EITHER an insertion or
/// a deletion opened.
///
/// The reference reads a triangular cache built at class-init with `log1p`, and computes the same
/// expression only for a quality past `MAX_QUAL`. Both paths are written out here because the
/// cached one is not `1 - 10^log10sum`: it is `exp(log1p(-min(1, 10^log10sum)))`, which is a
/// different double for a small sum.
pub fn match_to_match_prob(insertion_qual: i32, deletion_qual: i32) -> f64 {
    let (min_qual, max_qual) = if insertion_qual <= deletion_qual {
        (insertion_qual, deletion_qual)
    } else {
        (deletion_qual, insertion_qual)
    };
    let log10_sum =
        approximate_log10_sum_log10(-0.1 * f64::from(min_qual), -0.1 * f64::from(max_qual));
    if max_qual as usize > MAX_QUAL {
        return 1.0 - 10f64.powf(log10_sum);
    }
    // `Math.log1p(-Math.min(1, Math.pow(10, log10Sum))) * INV_LN10`, then `Math.pow(10, that)`.
    let log10 = (-(10f64.powf(log10_sum)).min(1.0)).ln_1p() * std::f64::consts::LOG10_E;
    10f64.powf(log10)
}

/// `PairHMMModel.qualToTransProbs` for one base.
pub fn qual_to_trans_probs(
    insertion_qual: u8,
    deletion_qual: u8,
    gap_continuation: u8,
) -> [f64; TRANS_PROB_ARRAY_LENGTH] {
    let mut dest = [0.0; TRANS_PROB_ARRAY_LENGTH];
    dest[MATCH_TO_MATCH] = match_to_match_prob(i32::from(insertion_qual), i32::from(deletion_qual));
    dest[MATCH_TO_INSERTION] = qual_to_error_prob(f64::from(insertion_qual));
    dest[MATCH_TO_DELETION] = qual_to_error_prob(f64::from(deletion_qual));
    dest[INDEL_TO_MATCH] = qual_to_prob(f64::from(gap_continuation));
    dest[INSERTION_TO_INSERTION] = qual_to_error_prob(f64::from(gap_continuation));
    dest[DELETION_TO_DELETION] = dest[INSERTION_TO_INSERTION];
    dest
}

/// `initializePriors`: the probability of the read's base given the haplotype's.
///
/// An `N` on either side matches, which is what keeps an uncalled base from costing a mismatch.
pub fn prior(read_base: u8, haplotype_base: u8, quality: u8) -> f64 {
    if read_base == haplotype_base || read_base == b'N' || haplotype_base == b'N' {
        qual_to_prob(f64::from(quality))
    } else {
        qual_to_error_prob(f64::from(quality)) / TRISTATE_CORRECTION
    }
}

/// `computeReadLikelihoodGivenHaplotypeLog10`, the whole recursion.
///
/// The three matrices are padded by one in each direction, the deletion row is seeded with the
/// initial condition spread over the haplotype's length, and the answer is the log of the last
/// row's match and insertion terms with that condition taken back off.
pub fn read_likelihood_given_haplotype_log10(
    haplotype: &[u8],
    read: &[u8],
    read_quals: &[u8],
    insertion_quals: &[u8],
    deletion_quals: &[u8],
    gap_continuation: &[u8],
) -> f64 {
    let padded_read = read.len() + 1;
    let padded_haplotype = haplotype.len() + 1;
    let mut matches = vec![vec![0.0f64; padded_haplotype]; padded_read];
    let mut insertions = vec![vec![0.0f64; padded_haplotype]; padded_read];
    let mut deletions = vec![vec![0.0f64; padded_haplotype]; padded_read];

    // The seed is spread over the haplotype, so a longer haplotype starts from a smaller value and
    // the answer does not depend on its length.
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

    for i in 1..padded_read {
        for j in 1..padded_haplotype {
            let p = prior(read[i - 1], haplotype[j - 1], read_quals[i - 1]);
            matches[i][j] = p
                * (matches[i - 1][j - 1] * transitions[i][MATCH_TO_MATCH]
                    + insertions[i - 1][j - 1] * transitions[i][INDEL_TO_MATCH]
                    + deletions[i - 1][j - 1] * transitions[i][INDEL_TO_MATCH]);
            insertions[i][j] = matches[i - 1][j] * transitions[i][MATCH_TO_INSERTION]
                + insertions[i - 1][j] * transitions[i][INSERTION_TO_INSERTION];
            deletions[i][j] = matches[i][j - 1] * transitions[i][MATCH_TO_DELETION]
                + deletions[i][j - 1] * transitions[i][DELETION_TO_DELETION];
        }
    }

    let end = padded_read - 1;
    let mut total = 0.0;
    for j in 1..padded_haplotype {
        total += matches[end][j] + insertions[end][j];
    }
    total.log10() - initial_condition_log10()
}
