//! `BaseUtils`, ported from `org.broadinstitute.hellbender.utils.BaseUtils` (GATK 4.6.2.0).
//!
//! The four functions that turn a base into a number and back. They live in one module because a
//! second copy of them is a second definition of what a base is, and three of them have a case a
//! reader would not guess:
//!
//!  * [`simple_base_to_base_index`] maps `*` to `A`, so a wildcard is a base as far as every caller
//!    is concerned, while `N` is not;
//!  * [`base_index_to_simple_base`] answers `.` for anything outside `0..=3` rather than failing,
//!    so a key with a bad length nibble decodes to a string of dots rather than to an error;
//!  * [`simple_complement`] **uppercases**: the complement of `a` is `T`, not `t`. htsjdk's
//!    `SequenceUtil.complement` preserves case, and the two are not interchangeable.

/// `BaseUtils.simpleBaseToBaseIndex`.
///
/// `A`, `C`, `G` and `T` in either case, plus one entry that is not a base at all: the wildcard `*`
/// maps to `A`. Everything else, `N` included, is `-1`.
pub fn simple_base_to_base_index(base: u8) -> i32 {
    match base {
        b'A' | b'a' | b'*' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => -1,
    }
}

/// `BaseUtils.baseIndexToSimpleBase`: the inverse, with `.` for every index it does not know.
///
/// It does not fail. `ContextCovariate::context_from_key` walks a length taken from the key's low
/// four bits, so a key claiming more bases than it holds decodes to dots rather than to an error,
/// and the golden carries that: key 4095 answers `TTTTAAAAAAAAAAA`.
pub fn base_index_to_simple_base(base_index: i32) -> u8 {
    match base_index {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'.',
    }
}

/// `BaseUtils.simpleComplement`: **uppercasing**, and the identity on anything else.
pub fn simple_complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        other => other,
    }
}

/// `BaseUtils.simpleReverseComplement`.
pub fn simple_reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|b| simple_complement(*b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_wildcard_is_an_a_and_an_n_is_nothing() {
        assert_eq!(simple_base_to_base_index(b'*'), 0);
        assert_eq!(simple_base_to_base_index(b'A'), 0);
        assert_eq!(simple_base_to_base_index(b'N'), -1);
        assert_eq!(simple_base_to_base_index(b'n'), -1);
    }

    #[test]
    fn an_unknown_index_is_a_dot_and_not_a_failure() {
        assert_eq!(base_index_to_simple_base(3), b'T');
        assert_eq!(base_index_to_simple_base(4), b'.');
        assert_eq!(base_index_to_simple_base(-1), b'.');
    }

    #[test]
    fn the_complement_uppercases() {
        assert_eq!(simple_complement(b'a'), b'T');
        assert_eq!(simple_complement(b'N'), b'N');
        assert_eq!(simple_reverse_complement(b"acgt"), b"ACGT".to_vec());
    }
}
