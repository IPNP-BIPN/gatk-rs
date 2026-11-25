//! `BaseUtils`, ported from `org.broadinstitute.hellbender.utils.BaseUtils` (GATK 4.7.0.0;
//! `BaseUtils.java` is unchanged since 4.6.2.0).
//!
//! The four functions that turn a base into a number and back. They live in one module because a
//! second copy of them is a second definition of what a base is, and three of them have a case a
//! reader would not guess:
//!
//!  * [`simple_base_to_base_index`] maps `*` to `A`, so a wildcard is a base as far as every caller
//!    is concerned, while `N` is not;
//!  * [`base_index_to_simple_base`] answers `.` for anything outside `0..=3` rather than failing;
//!  * [`simple_complement`] **uppercases**: the complement of `a` is `T`, not `t`. htsjdk's
//!    `SequenceUtil.complement` preserves case, and the two are not interchangeable.
//!
//! # Reading this file without reading Rust
//!
//! - `u8` is an unsigned byte. Java's `byte` is signed, which is why the reference masks with
//!   `& 0xff` before it indexes its table; a `u8` needs no mask.
//! - `b'A'` is the byte value of the character `A`, Java's `(byte) 'A'`.
//! - `match` is Java's `switch` without fall-through. `b'A' | b'a'` is two case labels sharing one
//!   arm, and `_` (or a bare name such as `other`) is `default`. The compiler refuses a `match`
//!   that leaves a value uncovered, so every function below answers for all 256 bytes.
//! - `i32` is Java's `int`; a function's last expression is its return value.

/// `BaseUtils.simpleBaseToBaseIndex`.
///
/// `A`, `C`, `G` and `T` in either case, plus one entry that is not a base at all: the wildcard `*`
/// maps to `A`. Everything else, `N` included, is `-1`.
///
/// The reference is a 256-entry table filled with `-1` and then set for these nine bytes, one of
/// them with its own note:
///
/// ```java
/// baseIndexMap['*'] = Base.A.ordinal();    // the wildcard character counts as an A
/// ```
///
/// A `match` over the same nine bytes is that table written out. The class's other table,
/// `baseIndexWithIupacMap`, sets `*` back to `-1` and maps the IUPAC codes to `N`; only
/// `convertIUPACtoN` reads it, and nothing here does.
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
/// It does not fail, and the `.` is only reached by an index outside `0..=3`.
/// `ContextCovariate::context_from_key` never produces one: it reads each index through a two-bit
/// mask, so a key whose low four bits claim more bases than were written decodes the extra ones
/// as `A`, not as `.`. Its test pins that: key 4095 answers `TTTTAAAAAAAAAAA`, with no dot in it.
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
///
/// The reference's `switch` returns `'T'` for both `'A'` and `'a'`, and so on for the other three;
/// `default: return base;` is the `other => other` arm, so `N`, `n` and `*` come back unchanged and
/// in their own case.
pub fn simple_complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        other => other,
    }
}

/// `BaseUtils.simpleReverseComplement`: a new array, the input read backwards and complemented.
///
/// The reference fills `rcbases[i] = simpleComplement(bases[bases.length - 1 - i])`. Here
/// `.iter().rev()` walks the input from its last byte, `.map(..)` complements each, and
/// `.collect()` gathers them into a new `Vec<u8>` (Java's `byte[]`, growable). `|b| ...` is a
/// closure, Java's lambda, and `*b` reads the byte the iterator lends rather than the reference to
/// it. `&[u8]` is a borrowed, read-only view of the caller's array, which is left untouched.
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
