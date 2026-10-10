//! `SmithWatermanJavaAligner`: the pairwise aligner `HaplotypeCaller` and `Mutect2` realign
//! haplotypes and reads with.
//!
//! It is not textbook Smith-Waterman. Only the rightmost column and the bottom row are searched for
//! the best cell, which is what lets a read overhang the reference, and four overhang strategies
//! decide what the overhang becomes. The affine gap is evaluated in linear time per cell by keeping
//! the best open gap per column and per row, which is valid only because the gap penalty is linear
//! in its length.
//!
//! What a port must keep, because every one of them moves a CIGAR:
//!
//!  * the traceback priority is diagonal, then a gap in the row direction (an insertion), then a
//!    gap in the column direction (a deletion), with `>=` so a tie stays on the diagonal;
//!  * the rightmost-column scan takes the LAST best row (`>=`), the bottom-row scan takes the first
//!    one that is strictly better, or equal and closer to the diagonal;
//!  * an alternate that occurs exactly in the reference is answered by a substring search that
//!    takes the LAST occurrence, for `SOFTCLIP` and `IGNORE` only, before the matrix is built;
//!  * the matrix is floored at `-1e8`, and an empty sequence is refused.
//!
//! The Intel aligner (GKL, JNI) is a different implementation behind the same interface; which one
//! `FASTEST_AVAILABLE` resolves to depends on the host, and the `smith-waterman` conformance suite
//! records which one loaded on the runner. This is the Java one.
//!
//! Ported from `org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner`,
//! `SmithWatermanAlignmentConstants`, the `SWParameters` and `SWOverhangStrategy` of the GATK
//! native bindings, and `Utils.lastIndexOf`.

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

use crate::alignment_utils::AlignmentError;

/// `SWParameters`: the four weights of an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SwParameters {
    pub match_value: i32,
    pub mismatch_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
}

impl SwParameters {
    pub const fn new(
        match_value: i32,
        mismatch_penalty: i32,
        gap_open_penalty: i32,
        gap_extend_penalty: i32,
    ) -> Self {
        SwParameters {
            match_value,
            mismatch_penalty,
            gap_open_penalty,
            gap_extend_penalty,
        }
    }
}

/// `SmithWatermanAlignmentConstants.ORIGINAL_DEFAULT`, which only test code uses.
pub const ORIGINAL_DEFAULT: SwParameters = SwParameters::new(3, -1, -4, -3);
/// `SmithWatermanAlignmentConstants.STANDARD_NGS`: dangling head and tail recovery.
pub const STANDARD_NGS: SwParameters = SwParameters::new(25, -50, -110, -6);
/// `SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS`: haplotype to reference, and heavily
/// favouring indels over substitutions.
pub const NEW_SW_PARAMETERS: SwParameters = SwParameters::new(200, -150, -260, -11);
/// `SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS`: read to haplotype.
pub const ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS: SwParameters =
    SwParameters::new(10, -15, -30, -5);

/// `SWOverhangStrategy`: what to do with the part of the alternate that hangs off the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SwOverhangStrategy {
    /// Treat overhangs as insertions and deletions, starting from the corner of the matrix.
    Indel,
    /// Soft-clip them.
    SoftClip,
    /// Allow an indel at the start but not at the end.
    LeadingIndel,
    /// Count them in the first element, and report a negative offset.
    Ignore,
}

/// `SmithWatermanAlignment`: a CIGAR and where on the reference it starts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SwAlignment {
    pub cigar: Cigar,
    pub alignment_offset: i32,
}

/// `SmithWatermanAligner`.
pub trait SmithWatermanAligner {
    fn align(
        &self,
        reference: &[u8],
        alternate: &[u8],
        parameters: &SwParameters,
        strategy: SwOverhangStrategy,
    ) -> Result<SwAlignment, AlignmentError>;
}

/// `SmithWatermanJavaAligner`, stateless.
#[derive(Debug, Clone, Copy, Default)]
pub struct SmithWatermanJavaAligner;

/// `MATRIX_MIN_CUTOFF`: `(int) -1.0e8`, the floor no cell drops below.
const MATRIX_MIN_CUTOFF: i32 = -100_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Match,
    Insertion,
    Deletion,
    Clip,
}

/// `Utils.lastIndexOf`: the last position of `query` in `reference`, searched from the right.
pub fn last_index_of(reference: &[u8], query: &[u8]) -> Option<usize> {
    if query.len() > reference.len() {
        return None;
    }
    (0..=reference.len() - query.len())
        .rev()
        .find(|&r| reference[r..r + query.len()] == *query)
}

impl SmithWatermanAligner for SmithWatermanJavaAligner {
    fn align(
        &self,
        reference: &[u8],
        alternate: &[u8],
        parameters: &SwParameters,
        strategy: SwOverhangStrategy,
    ) -> Result<SwAlignment, AlignmentError> {
        if reference.is_empty() || alternate.is_empty() {
            return Err(AlignmentError::EmptySequence);
        }

        // Skip the matrix when the alternate occurs exactly in the reference. Only valid for the
        // two strategies that do not account for overhangs as indels.
        if matches!(
            strategy,
            SwOverhangStrategy::SoftClip | SwOverhangStrategy::Ignore
        ) {
            if let Some(index) = last_index_of(reference, alternate) {
                return Ok(SwAlignment {
                    cigar: Cigar::new(vec![CigarElement {
                        length: alternate.len() as u32,
                        op: Op::M,
                    }]),
                    alignment_offset: index as i32,
                });
            }
        }

        let n = reference.len() + 1;
        let m = alternate.len() + 1;
        let mut sw = vec![vec![0i32; m]; n];
        let mut btrack = vec![vec![0i32; m]; n];
        calculate_matrix(
            reference,
            alternate,
            &mut sw,
            &mut btrack,
            strategy,
            parameters,
        );
        Ok(calculate_cigar(&sw, &btrack, strategy))
    }
}

/// `calculateMatrix`: fill the score matrix and the back track matrix.
fn calculate_matrix(
    reference: &[u8],
    alternate: &[u8],
    sw: &mut [Vec<i32>],
    btrack: &mut [Vec<i32>],
    strategy: SwOverhangStrategy,
    parameters: &SwParameters,
) {
    let ncol = sw[0].len();
    let nrow = sw.len();

    let low_init_value = i32::MIN / 2;
    let mut best_gap_v = vec![low_init_value; ncol + 1];
    let mut gap_size_v = vec![0i32; ncol + 1];
    let mut best_gap_h = vec![low_init_value; nrow + 1];
    let mut gap_size_h = vec![0i32; nrow + 1];

    let w_open = parameters.gap_open_penalty;
    let w_extend = parameters.gap_extend_penalty;
    let w_match = parameters.match_value;
    let w_mismatch = parameters.mismatch_penalty;

    // Gap penalties on the first row and column, to keep track of indels at the edges.
    if matches!(
        strategy,
        SwOverhangStrategy::Indel | SwOverhangStrategy::LeadingIndel
    ) {
        sw[0][1] = w_open;
        let mut current = w_open;
        for cell in sw[0].iter_mut().skip(2) {
            current = current.wrapping_add(w_extend);
            *cell = current;
        }
        sw[1][0] = w_open;
        current = w_open;
        for row in sw.iter_mut().skip(2) {
            current = current.wrapping_add(w_extend);
            row[0] = current;
        }
    }

    for i in 1..nrow {
        let a_base = reference[i - 1];
        let (last_rows, cur_rows) = sw.split_at_mut(i);
        let last_row = &last_rows[i - 1];
        let cur_row = &mut cur_rows[0];
        let cur_back_track_row = &mut btrack[i];

        for j in 1..ncol {
            let b_base = alternate[j - 1];
            let step_diag = last_row[j - 1].wrapping_add(if a_base == b_base {
                w_match
            } else {
                w_mismatch
            });

            // Best gap arriving from above: opened just now, or an earlier one extended by one.
            let mut prev_gap = last_row[j].wrapping_add(w_open);
            best_gap_v[j] = best_gap_v[j].wrapping_add(w_extend);
            if prev_gap > best_gap_v[j] {
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            } else {
                gap_size_v[j] += 1;
            }
            let step_down = best_gap_v[j];
            let kd = gap_size_v[j];

            // Best gap arriving from the left.
            prev_gap = cur_row[j - 1].wrapping_add(w_open);
            best_gap_h[i] = best_gap_h[i].wrapping_add(w_extend);
            if prev_gap > best_gap_h[i] {
                best_gap_h[i] = prev_gap;
                gap_size_h[i] = 1;
            } else {
                gap_size_h[i] += 1;
            }
            let step_right = best_gap_h[i];
            let ki = gap_size_h[i];

            // Priority: diagonal, then right, then down.
            let diag_highest_or_equal = step_diag >= step_down && step_diag >= step_right;
            if diag_highest_or_equal {
                cur_row[j] = MATRIX_MIN_CUTOFF.max(step_diag);
                cur_back_track_row[j] = 0;
            } else if step_right >= step_down {
                cur_row[j] = MATRIX_MIN_CUTOFF.max(step_right);
                cur_back_track_row[j] = -ki;
            } else {
                cur_row[j] = MATRIX_MIN_CUTOFF.max(step_down);
                cur_back_track_row[j] = kd;
            }
        }
    }
}

fn make_element(state: State, length: i32) -> CigarElement {
    let op = match state {
        State::Match => Op::M,
        State::Insertion => Op::I,
        State::Deletion => Op::D,
        State::Clip => Op::S,
    };
    CigarElement {
        length: length as u32,
        op,
    }
}

/// `calculateCigar`: walk the back track matrix from the best cell to the origin.
fn calculate_cigar(
    sw: &[Vec<i32>],
    btrack: &[Vec<i32>],
    strategy: SwOverhangStrategy,
) -> SwAlignment {
    let mut p1: i32 = 0;
    let mut p2: i32;

    let ref_length = (sw.len() - 1) as i32;
    let alt_length = (sw[0].len() - 1) as i32;

    let mut max_score = i32::MIN;
    let mut segment_length: i32 = 0;

    if strategy == SwOverhangStrategy::Indel {
        p1 = ref_length;
        p2 = alt_length;
    } else {
        // The largest score on the rightmost column; `>=` plus the traversal direction picks the
        // one closer to the diagonal on a tie.
        p2 = alt_length;
        for (i, row) in sw.iter().enumerate().skip(1) {
            let cur_score = row[alt_length as usize];
            if cur_score >= max_score {
                p1 = i as i32;
                max_score = cur_score;
            }
        }
        // A larger score on the bottom-most row.
        if strategy != SwOverhangStrategy::LeadingIndel {
            let bottom_row = &sw[ref_length as usize];
            for (j, &cur_score) in bottom_row.iter().enumerate().skip(1) {
                let j = j as i32;
                if cur_score > max_score
                    || (cur_score == max_score && (ref_length - j).abs() < (p1 - p2).abs())
                {
                    p1 = ref_length;
                    p2 = j;
                    max_score = cur_score;
                    segment_length = alt_length - j;
                }
            }
        }
    }

    let mut lce: Vec<CigarElement> = Vec::with_capacity(5);
    if segment_length > 0 && strategy == SwOverhangStrategy::SoftClip {
        lce.push(make_element(State::Clip, segment_length));
        segment_length = 0;
    }

    // Insertions and deletions are placed in the alternate, so the states are named for it.
    let mut state = State::Match;
    loop {
        let btr = btrack[p1 as usize][p2 as usize];
        let mut step_length = 1;
        let new_state = if btr > 0 {
            step_length = btr;
            State::Deletion
        } else if btr < 0 {
            step_length = -btr;
            State::Insertion
        } else {
            State::Match
        };

        match new_state {
            State::Match => {
                p1 -= 1;
                p2 -= 1;
            }
            State::Insertion => p2 -= step_length,
            State::Deletion => p1 -= step_length,
            State::Clip => {}
        }

        if new_state == state {
            segment_length += step_length;
        } else {
            if segment_length > 0 {
                lce.push(make_element(state, segment_length));
            }
            segment_length = step_length;
            state = new_state;
        }
        if !(p1 > 0 && p2 > 0) {
            break;
        }
    }

    let alignment_offset = match strategy {
        SwOverhangStrategy::SoftClip => {
            lce.push(make_element(state, segment_length));
            if p2 > 0 {
                lce.push(make_element(State::Clip, p2));
            }
            p1
        }
        SwOverhangStrategy::Ignore => {
            lce.push(make_element(state, segment_length + p2));
            p1 - p2
        }
        SwOverhangStrategy::Indel | SwOverhangStrategy::LeadingIndel => {
            lce.push(make_element(state, segment_length));
            if p1 > 0 {
                lce.push(make_element(State::Deletion, p1));
            } else if p2 > 0 {
                lce.push(make_element(State::Insertion, p2));
            }
            0
        }
    };

    lce.reverse();
    SwAlignment {
        cigar: Cigar::new(lce),
        alignment_offset,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn align(reference: &str, alternate: &str, strategy: SwOverhangStrategy) -> (i32, String) {
        let result = SmithWatermanJavaAligner
            .align(
                reference.as_bytes(),
                alternate.as_bytes(),
                &NEW_SW_PARAMETERS,
                strategy,
            )
            .expect("an alignment");
        (result.alignment_offset, result.cigar.to_text())
    }

    #[test]
    fn exact_substring_takes_the_last_occurrence_for_softclip() {
        let (offset, cigar) = align(
            "TTGACCAGTCATTGACCAGTCATTTGAC",
            "TTGACCAGTCA",
            SwOverhangStrategy::SoftClip,
        );
        assert_eq!((offset, cigar.as_str()), (11, "11M"));
    }

    #[test]
    fn indel_strategy_skips_the_substring_search() {
        let (offset, cigar) = align(
            "TTGACCAGTCATTGACCAGTCATTTGAC",
            "TTGACCAGTCA",
            SwOverhangStrategy::Indel,
        );
        assert_eq!((offset, cigar.as_str()), (0, "11M17D"));
    }

    #[test]
    fn overhang_is_clipped_or_counted() {
        assert_eq!(
            align(
                "ACGTCCAGTTGACCAT",
                "GGGGGACGTCCAGTTGACCAT",
                SwOverhangStrategy::SoftClip
            ),
            (0, "5S16M".to_string())
        );
        assert_eq!(
            align(
                "ACGTCCAGTTGACCAT",
                "GGGGGACGTCCAGTTGACCAT",
                SwOverhangStrategy::Ignore
            ),
            (-5, "21M".to_string())
        );
        assert_eq!(
            align(
                "ACGTCCAGTTGACCAT",
                "GGGGGACGTCCAGTTGACCAT",
                SwOverhangStrategy::Indel
            ),
            (0, "5I16M".to_string())
        );
    }

    #[test]
    fn empty_input_is_refused() {
        assert_eq!(
            SmithWatermanJavaAligner.align(
                b"",
                b"A",
                &NEW_SW_PARAMETERS,
                SwOverhangStrategy::Indel
            ),
            Err(AlignmentError::EmptySequence)
        );
    }
}
