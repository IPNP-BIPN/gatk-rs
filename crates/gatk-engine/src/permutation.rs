//! Ported from `org.apache.commons.math3.random.RandomDataGenerator.nextPermutation`,
//! `org.apache.commons.math3.util.MathArrays.shuffle` and
//! `org.apache.commons.math3.distribution.UniformIntegerDistribution.sample` (commons-math3 3.5,
//! Apache 2.0), which together are what
//! `org.broadinstitute.hellbender.utils.MathUtils.sampleIndicesWithoutReplacement` is:
//!
//! ```text
//! public static int[] sampleIndicesWithoutReplacement(final int n, final int k) {
//!     //No error checking : RandomDataGenetator.nextPermutation does it
//!     return Utils.getRandomDataGenerator().nextPermutation(n, k);
//! }
//! ```
//!
//! Three things about it decide the answer and none is inherent to "sample k of n":
//!
//!  * **the whole array is shuffled, then truncated.** `nextPermutation(n, k)` shuffles all `n`
//!    entries and returns the first `k`, so it consumes `n - 1` draws whatever `k` is. Sampling 2
//!    of 1000 costs 999 draws from the shared stream, and everything drawn later in the run moves
//!    with it;
//!  * **the shuffle runs downwards.** Fisher-Yates from `list.length - 1` to `start`, with the
//!    target drawn from `[start, i]`, so the draw bounds *shrink*. Running it upwards produces a
//!    valid permutation from the same draws and a different one;
//!  * **the last step takes no draw.** At `i == start` the target is `start` without a call, so
//!    the count is `n - 1` and not `n`. One extra draw here shifts every later consumer.
//!
//! The returned indices are the shuffled array's head, so they are **not sorted**, and the caller
//! (`LevelingDownsampler`) only ever uses them as set membership. That is why the golden compares
//! the raw index list too: a port that sorted them would still keep the right reads and would
//! still be wrong here.
//!
//! # Reading this file without reading Rust
//!
//! - `&mut [i32]` is a **writable borrowed view** of an `int[]`: the function may change the
//!   elements in place but does not own the array, and while it holds the view nothing else may
//!   touch it. `&mut Well19937c` is the same for the generator, which is how a draw advances the
//!   caller's stream rather than a copy of it.
//! - `Result<Vec<i32>, PermutationError>` is "an `int[]` or one of the two refusals below"; it is
//!   how this codebase models a method that throws. `Err(..)` is the throw, `Ok(..)` the return.
//! - `wrapping_sub` and `wrapping_add` are Java's `int` arithmetic, which wraps silently on
//!   overflow. Plain `-` and `+` on an `i32` would panic in a debug build instead.

use crate::well19937c::Well19937c;

/// The two ways `nextPermutation` refuses, kept apart because they are different exception types
/// upstream and a caller can tell them apart.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PermutationError {
    /// `NumberIsTooLargeException`: `k > n`.
    SizeExceedsN { n: i32, k: i32 },
    /// `NotStrictlyPositiveException`: `k <= 0`. `LevelingDownsampler` can reach this with
    /// `minElementsPerStack` of zero, which is why it is a refusal and not a clamp.
    NotStrictlyPositive { k: i32 },
}

impl std::fmt::Display for PermutationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PermutationError::SizeExceedsN { n, k } => {
                write!(f, "permutation k ({k}) exceeds n ({n})")
            }
            PermutationError::NotStrictlyPositive { k } => {
                write!(f, "permutation size ({k}) must be positive")
            }
        }
    }
}

impl std::error::Error for PermutationError {}

/// `MathArrays.natural(n)`: the array `[0, 1, ..., n - 1]`, the identity before the shuffle.
///
/// `(0..n as i32)` is the half-open range `0` to `n - 1` and `.collect()` gathers it into a new
/// `Vec`, which is Java's `int[]` with a length that may grow.
pub fn natural(n: usize) -> Vec<i32> {
    (0..n as i32).collect()
}

/// `UniformIntegerDistribution(rng, lower, upper).sample()`: one draw, uniform over the closed
/// range `[lower, upper]`.
///
/// **How**: when the width `upper - lower + 1` fits in an `int`, one bounded draw
/// (`nextInt(max)`) offset by `lower`; otherwise unbounded draws until one lands in range.
///
/// The `max <= 0` branch is the overflow case, when the range is wider than `i32` can count. It
/// cannot be reached from `shuffle`, whose lower bound is `start` and whose upper is an index, but
/// it is kept because it is what makes the method total and because leaving it out would look like
/// a decision rather than an omission.
fn uniform_sample(rng: &mut Well19937c, lower: i32, upper: i32) -> i32 {
    let max = upper.wrapping_sub(lower).wrapping_add(1);
    if max <= 0 {
        loop {
            let r = rng.next_int();
            if r >= lower && r <= upper {
                return r;
            }
        }
    } else {
        lower.wrapping_add(rng.next_int_bound(max))
    }
}

/// `MathArrays.shuffle(int[] list, RandomGenerator rng)`, which is `shuffle(list, 0, TAIL, rng)`:
/// the whole array, in place.
pub fn shuffle(list: &mut [i32], rng: &mut Well19937c) {
    shuffle_tail(list, 0, rng);
}

/// `MathArrays.shuffle(int[] list, int start, Position.TAIL, RandomGenerator rng)`: Fisher-Yates
/// over `list[start..]`, in place.
///
/// **How**: for `i` from the last index down to `start`, swap `list[i]` with `list[target]`, where
/// `target` is drawn uniformly from `[start, i]`.
///
/// **Why a `loop` with a `break` rather than a `for` over a reversed range**: the index is
/// unsigned (`usize`), so a countdown that stops *after* `start` cannot be written as Java's
/// `i >= start; i--` when `start` is zero, because `i` would have to reach -1. Breaking at
/// `i == start` before the decrement visits the same indices in the same order.
pub fn shuffle_tail(list: &mut [i32], start: usize, rng: &mut Well19937c) {
    // An empty array has no last index to start from. The reference's loop would simply not run.
    if list.is_empty() {
        return;
    }
    let mut i = list.len() - 1;
    loop {
        // At `i == start` the target is `start` and no draw is taken. That is the reason the
        // shuffle costs `n - 1` draws rather than `n`.
        let target = if i == start {
            start
        } else {
            uniform_sample(rng, start as i32, i as i32) as usize
        };
        list.swap(target, i);
        if i == start {
            break;
        }
        i -= 1;
    }
}

/// `RandomDataGenerator.nextPermutation(int n, int k)`: `k` distinct indices from `[0, n)`, in the
/// order the shuffle left them.
///
/// **How**: the two refusals, `k > n` checked first, so `k = 5, n = 3` reports the size and not
/// the sign; then shuffle the whole identity array and keep its first `k`.
///
/// **Why the whole array is shuffled**: see the module note. It is the draw count, `n - 1` for any
/// `k`, that the shared stream sees.
pub fn next_permutation(
    n: i32,
    k: i32,
    rng: &mut Well19937c,
) -> Result<Vec<i32>, PermutationError> {
    if k > n {
        return Err(PermutationError::SizeExceedsN { n, k });
    }
    if k <= 0 {
        return Err(PermutationError::NotStrictlyPositive { k });
    }
    // `n.max(0)` only keeps the cast to an unsigned length total. After the two checks above
    // `n >= k >= 1`, so it never changes `n`.
    let mut index = natural(n.max(0) as usize);
    shuffle(&mut index, rng);
    // `truncate` keeps the first `k` elements and drops the rest, in place.
    // The last line, with no `return` and no semicolon, is the function's value.
    index.truncate(k as usize);
    Ok(index)
}

/// `MathUtils.sampleIndicesWithoutReplacement`, which is `nextPermutation` under another name and
/// with its error checking left to it.
pub fn sample_indices_without_replacement(
    n: i32,
    k: i32,
    rng: &mut Well19937c,
) -> Result<Vec<i32>, PermutationError> {
    next_permutation(n, k, rng)
}
