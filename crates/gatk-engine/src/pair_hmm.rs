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
    jmath::math::log10(initial_condition())
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

/// `JacobianLogTable.get`: the table's entry at the difference's quantisation.
///
/// The reference builds the table once and reads it by index; so does this, and every entry is
/// the reference's expression, so the value is the table's entry and not `log10(1 + 10^-difference)`.
/// An index past the table, which `approximate_log10_sum_log10` never asks for, is computed.
pub fn jacobian_log(difference: f64) -> f64 {
    // `fastRound(difference * INV_STEP)`, a product by `1.0 / TABLE_STEP` rather than a quotient
    // by the step.
    let index = fast_round(difference * (1.0 / JACOBIAN_TABLE_STEP));
    static TABLE: std::sync::LazyLock<Vec<f64>> = std::sync::LazyLock::new(|| {
        // `(int) (MAX_TOLERANCE / TABLE_STEP) + 1` entries.
        let size = (JACOBIAN_MAX_TOLERANCE / JACOBIAN_TABLE_STEP) as i32 + 1;
        (0..size).map(jacobian_entry).collect()
    });
    match usize::try_from(index).ok().and_then(|i| TABLE.get(i)) {
        Some(&entry) => entry,
        None => jacobian_entry(index),
    }
}

/// One entry: `Math.log10(1.0 + Math.pow(10.0, -k * TABLE_STEP))`, whose `Math.log10` is correctly
/// rounded.
fn jacobian_entry(index: i32) -> f64 {
    let quantised = f64::from(-index) * JACOBIAN_TABLE_STEP;
    jmath::math::log10(1.0 + std::hint::black_box(10f64).powf(quantised))
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
    // `INV_LN10` is `1.0 / Math.log(10)`, which is one ulp below `LOG10_E`.
    let inv_ln10 = 1.0 / std::f64::consts::LN_10;
    let log10 = (-(10f64.powf(log10_sum)).min(1.0)).ln_1p() * inv_ln10;
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

/// Every quality's three probabilities, computed once by the functions above.
struct Tables {
    /// `qual_to_prob(q)`.
    prob: [f64; 256],
    /// `qual_to_error_prob(q)`.
    error: [f64; 256],
    /// `qual_to_error_prob(q) / TRISTATE_CORRECTION`, what [`prior`] gives a mismatch.
    mismatch: [f64; 256],
    /// `match_to_match_prob(insertion, deletion)`, indexed `insertion * 256 + deletion`.
    match_to_match: Vec<f64>,
}

static TABLES: std::sync::LazyLock<Tables> = std::sync::LazyLock::new(|| {
    let prob = std::array::from_fn(|q| qual_to_prob(q as f64));
    let error = std::array::from_fn(|q| qual_to_error_prob(q as f64));
    let mismatch = std::array::from_fn(|q| qual_to_error_prob(q as f64) / TRISTATE_CORRECTION);
    let match_to_match = (0..256 * 256)
        .map(|k: i32| match_to_match_prob(k / 256, k % 256))
        .collect();
    Tables {
        prob,
        error,
        mismatch,
        match_to_match,
    }
});

/// `computeReadLikelihoodGivenHaplotypeLog10`, the whole recursion.
///
/// The three matrices are padded by one in each direction, the deletion row is seeded with the
/// initial condition spread over the haplotype's length, and the answer is the log of the last
/// row's match and insertion terms with that condition taken back off.
///
/// Each cell does the reference's operations in the reference's order, so the bits are its own;
/// what is saved is the rest. Only two rows of each matrix are kept, and the transitions and
/// priors come from [`TABLES`], built by [`qual_to_trans_probs`]' and [`prior`]'s own functions,
/// instead of a `pow` per cell.
pub fn read_likelihood_given_haplotype_log10(
    haplotype: &[u8],
    read: &[u8],
    read_quals: &[u8],
    insertion_quals: &[u8],
    deletion_quals: &[u8],
    gap_continuation: &[u8],
) -> f64 {
    let tables = &*TABLES;
    let width = haplotype.len() + 1;
    // Row i - 1 and row i of the match, insertion and deletion matrices.
    let mut rows = vec![0.0f64; 6 * width];
    let (previous, current) = rows.split_at_mut(3 * width);
    let (mut previous_m, rest) = previous.split_at_mut(width);
    let (mut previous_x, mut previous_y) = rest.split_at_mut(width);
    let (mut current_m, rest) = current.split_at_mut(width);
    let (mut current_x, mut current_y) = rest.split_at_mut(width);

    // The seed is spread over the haplotype, so a longer haplotype starts from a smaller value and
    // the answer does not depend on its length. Row zero's column zero carries it too.
    let initial = initial_condition() / haplotype.len() as f64;
    previous_y.fill(initial);

    for i in 0..read.len() {
        let (insertion, deletion, gap) = (
            usize::from(insertion_quals[i]),
            usize::from(deletion_quals[i]),
            usize::from(gap_continuation[i]),
        );
        let match_to_match = tables.match_to_match[insertion * 256 + deletion];
        let indel_to_match = tables.prob[gap];
        let match_to_insertion = tables.error[insertion];
        let insertion_to_insertion = tables.error[gap];
        let match_to_deletion = tables.error[deletion];
        let deletion_to_deletion = insertion_to_insertion;
        let quality = usize::from(read_quals[i]);
        let (matched, mismatched) = (tables.prob[quality], tables.mismatch[quality]);
        let read_base = read[i];

        // Column zero of every row past the first is zero.
        current_m[0] = 0.0;
        current_x[0] = 0.0;
        current_y[0] = 0.0;
        let (mut left_m, mut left_y) = (0.0, 0.0);
        for j in 1..width {
            let haplotype_base = haplotype[j - 1];
            let p = if read_base == haplotype_base || read_base == b'N' || haplotype_base == b'N' {
                matched
            } else {
                mismatched
            };
            let m = p
                * (previous_m[j - 1] * match_to_match
                    + previous_x[j - 1] * indel_to_match
                    + previous_y[j - 1] * indel_to_match);
            current_x[j] =
                previous_m[j] * match_to_insertion + previous_x[j] * insertion_to_insertion;
            let y = left_m * match_to_deletion + left_y * deletion_to_deletion;
            current_m[j] = m;
            current_y[j] = y;
            left_m = m;
            left_y = y;
        }
        std::mem::swap(&mut previous_m, &mut current_m);
        std::mem::swap(&mut previous_x, &mut current_x);
        std::mem::swap(&mut previous_y, &mut current_y);
    }

    // After the last swap the read's last row is in `previous`. A read with no bases reads row
    // zero, as the reference does.
    let mut total = 0.0;
    for j in 1..width {
        total += previous_m[j] + previous_x[j];
    }
    jmath::math::log10(total) - initial_condition_log10()
}

/// One read's bases and its four quality arrays, all of the same length, as
/// [`read_likelihood_given_haplotype_log10`] takes them.
#[derive(Debug, Clone, Copy)]
pub struct HmmRead<'a> {
    pub bases: &'a [u8],
    pub quals: &'a [u8],
    pub insertion_quals: &'a [u8],
    pub deletion_quals: &'a [u8],
    pub gap_continuation: &'a [u8],
}

/// One column of a row: the match, insertion and deletion values of every lane, side by side.
/// Each lane is one read and does that read's own multiplies and adds, never fused, so computing
/// the lanes in one vector instruction cannot move a bit.
#[derive(Clone, Copy)]
struct Cell<const N: usize> {
    m: [f64; N],
    x: [f64; N],
    y: [f64; N],
}

impl<const N: usize> Cell<N> {
    const ZERO: Self = Cell {
        m: [0.0; N],
        x: [0.0; N],
        y: [0.0; N],
    };
}

/// The two rows of each width, kept per thread between calls: allocating them per haplotype was a
/// tenth of a call.
type Rows<const N: usize> = std::cell::RefCell<(Vec<Cell<N>>, Vec<Cell<N>>)>;
thread_local! {
    static ROWS_4: Rows<4> = const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
    static ROWS_8: Rows<8> = const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

/// How [`reads_likelihoods_given_haplotype_log10`] runs its lanes. Every kernel does each lane's
/// scalar operations in the scalar order (`tools/audit/simd_kernels.json` declares this file), so
/// they differ in speed only; the test forces each one the host has against the single pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kernel {
    /// Four lanes, in whatever the target's baseline vectors are: SSE2 on x86-64, NEON on arm64.
    Portable,
    /// Four lanes in one 256-bit register.
    Avx2,
    /// Eight lanes in one 512-bit register.
    Avx512,
}

impl Kernel {
    /// Every kernel this host can run, the portable one first.
    pub fn available() -> Vec<Kernel> {
        #[cfg(target_arch = "x86_64")]
        {
            let mut kernels = vec![Kernel::Portable];
            if std::arch::is_x86_feature_detected!("avx2") {
                kernels.push(Kernel::Avx2);
            }
            if std::arch::is_x86_feature_detected!("avx512f") {
                kernels.push(Kernel::Avx512);
            }
            kernels
        }
        #[cfg(not(target_arch = "x86_64"))]
        vec![Kernel::Portable]
    }

    /// The widest kernel this host can run, detected once.
    pub fn best() -> Kernel {
        static BEST: std::sync::OnceLock<Kernel> = std::sync::OnceLock::new();
        *BEST.get_or_init(|| *Kernel::available().last().expect("the portable kernel"))
    }
}

/// [`read_likelihood_given_haplotype_log10`] for every read of `reads` against one haplotype,
/// into `out` in read order, with the same bits, on the widest kernel the host has.
pub fn reads_likelihoods_given_haplotype_log10(
    haplotype: &[u8],
    reads: &[HmmRead<'_>],
    out: &mut [f64],
) {
    reads_likelihoods_given_haplotype_log10_with(Kernel::best(), haplotype, reads, out);
}

/// [`reads_likelihoods_given_haplotype_log10`] on a named kernel, which is how the test reaches
/// each one.
///
/// # Panics
///
/// If the host cannot run `kernel` (it is not in [`Kernel::available`]).
pub fn reads_likelihoods_given_haplotype_log10_with(
    kernel: Kernel,
    haplotype: &[u8],
    reads: &[HmmRead<'_>],
    out: &mut [f64],
) {
    assert!(
        Kernel::available().contains(&kernel),
        "{kernel:?} is not available here"
    );
    match kernel {
        Kernel::Portable => lanes::<4>(&ROWS_4, haplotype, reads, out),
        #[cfg(target_arch = "x86_64")]
        // SAFETY: `available` detected the feature on this host.
        Kernel::Avx2 => unsafe { avx2(haplotype, reads, out) },
        #[cfg(target_arch = "x86_64")]
        // SAFETY: as above.
        Kernel::Avx512 => unsafe { avx512(haplotype, reads, out) },
        #[cfg(not(target_arch = "x86_64"))]
        Kernel::Avx2 | Kernel::Avx512 => unreachable!("refused by the assertion"),
    }
}

/// [`lanes`] compiled with AVX2 enabled and nothing else: the same source, so the same operations
/// per lane, in 256-bit registers. Fused multiply-add is deliberately not enabled (Rust would not contract
/// `a * b + c` anyway, and the guard refuses it here).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn avx2(haplotype: &[u8], reads: &[HmmRead<'_>], out: &mut [f64]) {
    lanes::<4>(&ROWS_4, haplotype, reads, out);
}

/// [`lanes`] with eight lanes and AVX-512F enabled.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn avx512(haplotype: &[u8], reads: &[HmmRead<'_>], out: &mut [f64]) {
    lanes::<8>(&ROWS_8, haplotype, reads, out);
}

/// The body every kernel shares, `N` reads at a time.
///
/// Reads are taken in order of length, `N` at a time, so a lane rarely idles. A lane whose read
/// has ended runs on zero transitions and priors, and its answer was taken at its last row. The
/// prior is looked up, not chosen: each distinct haplotype byte gets a code, and every row holds
/// the prior of each code for each lane, which keeps the inner loop free of branches.
/// `inline(always)` is what puts this body, and [`sweep`], inside each `target_feature` wrapper.
#[inline(always)]
fn lanes<const N: usize>(
    rows: &'static std::thread::LocalKey<Rows<N>>,
    haplotype: &[u8],
    reads: &[HmmRead<'_>],
    out: &mut [f64],
) {
    assert_eq!(reads.len(), out.len(), "one answer per read");
    let width = haplotype.len() + 1;
    let initial = initial_condition() / haplotype.len() as f64;
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|&r| reads[r].bases.len());

    let mut distinct: Vec<u8> = Vec::new();
    let codes: Vec<u8> = haplotype
        .iter()
        .map(|&base| match distinct.iter().position(|&d| d == base) {
            Some(code) => code as u8,
            None => {
                distinct.push(base);
                (distinct.len() - 1) as u8
            }
        })
        .collect();
    let mut priors = vec![[0.0f64; N]; distinct.len()];

    // Taken out and put back rather than borrowed through a closure: a closure is compiled once
    // for every caller of `lanes::<N>`, without the caller's `target_feature`, so the AVX2 path
    // would run the portable body.
    let (mut previous, mut current) = rows.take();
    previous.resize(width, Cell::ZERO);
    current.resize(width, Cell::ZERO);
    sweep(
        reads,
        out,
        &order,
        &codes,
        &distinct,
        &mut priors,
        &mut previous[..width],
        &mut current[..width],
        initial,
    );
    rows.set((previous, current));
}

/// The rows of [`lanes`], `N` reads at a time.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn sweep<'a, const N: usize>(
    reads: &[HmmRead<'_>],
    out: &mut [f64],
    order: &[usize],
    codes: &[u8],
    distinct: &[u8],
    priors: &mut [[f64; N]],
    mut previous: &'a mut [Cell<N>],
    mut current: &'a mut [Cell<N>],
    initial: f64,
) {
    let tables = &*TABLES;

    // The last row's match and insertion terms, summed in the reference's order.
    let answer = |row: &[Cell<N>], lane: usize| -> f64 {
        let mut total = 0.0;
        for cell in &row[1..] {
            total += cell.m[lane] + cell.x[lane];
        }
        jmath::math::log10(total) - initial_condition_log10()
    };

    for chunk in order.chunks(N) {
        previous.fill(Cell {
            y: [initial; N],
            ..Cell::ZERO
        });
        for (lane, &r) in chunk.iter().enumerate() {
            if reads[r].bases.is_empty() {
                out[r] = answer(previous, lane);
            }
        }

        let longest = chunk
            .iter()
            .map(|&r| reads[r].bases.len())
            .max()
            .unwrap_or(0);
        for i in 0..longest {
            let mut match_to_match = [0.0; N];
            let mut indel_to_match = [0.0; N];
            let mut match_to_insertion = [0.0; N];
            let mut insertion_to_insertion = [0.0; N];
            let mut match_to_deletion = [0.0; N];
            for prior in priors.iter_mut() {
                *prior = [0.0; N];
            }
            for (lane, &r) in chunk.iter().enumerate() {
                let read = &reads[r];
                if i >= read.bases.len() {
                    continue;
                }
                let insertion = usize::from(read.insertion_quals[i]);
                let deletion = usize::from(read.deletion_quals[i]);
                let gap = usize::from(read.gap_continuation[i]);
                let quality = usize::from(read.quals[i]);
                match_to_match[lane] = tables.match_to_match[insertion * 256 + deletion];
                indel_to_match[lane] = tables.prob[gap];
                match_to_insertion[lane] = tables.error[insertion];
                insertion_to_insertion[lane] = tables.error[gap];
                match_to_deletion[lane] = tables.error[deletion];
                // `prior`, for each haplotype byte.
                let read_base = read.bases[i];
                for (prior, &base) in priors.iter_mut().zip(distinct) {
                    prior[lane] = if read_base == base || read_base == b'N' || base == b'N' {
                        tables.prob[quality]
                    } else {
                        tables.mismatch[quality]
                    };
                }
            }
            let deletion_to_deletion = insertion_to_insertion;

            current[0] = Cell::ZERO;
            let mut left = Cell::ZERO;
            for ((cell, window), &code) in
                current[1..].iter_mut().zip(previous.windows(2)).zip(codes)
            {
                let (diagonal, up) = (&window[0], &window[1]);
                let p = &priors[usize::from(code)];
                let mut next = Cell::ZERO;
                for lane in 0..N {
                    next.m[lane] = p[lane]
                        * (diagonal.m[lane] * match_to_match[lane]
                            + diagonal.x[lane] * indel_to_match[lane]
                            + diagonal.y[lane] * indel_to_match[lane]);
                    next.x[lane] = up.m[lane] * match_to_insertion[lane]
                        + up.x[lane] * insertion_to_insertion[lane];
                    next.y[lane] = left.m[lane] * match_to_deletion[lane]
                        + left.y[lane] * deletion_to_deletion[lane];
                }
                *cell = next;
                left = next;
            }
            std::mem::swap(&mut previous, &mut current);
            for (lane, &r) in chunk.iter().enumerate() {
                if reads[r].bases.len() == i + 1 {
                    out[r] = answer(previous, lane);
                }
            }
        }
    }
}
