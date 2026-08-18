//! Ported from `org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter`.
//!
//! `CalculateContamination` reaches this class through `ContaminationSegmenter`, and it is the only
//! thing in the port that needs a singular value decomposition. This module carries the part the
//! `kernel-segmentation` golden pins directly: the subsampled kernel matrix that
//! [`crate::singular_value_decomposition`] is handed.
//!
//! # The subsample is seeded, and the seed is the whole comparison
//!
//! `new Random(1216)` through `RandomGeneratorFactory.createRandomGenerator`, then
//! `rng.nextInt(data.size())` once per subsampled point. The factory's `nextInt(n)` delegates
//! straight to the wrapped `java.util.Random`, so [`crate::java_random`] is the right generator and
//! not [`crate::well19937c`]. A port that subsampled differently would decompose a different matrix
//! and every comparison downstream would be meaningless rather than wrong.
//!
//! Note the reference draws **with replacement** and does not sort: the same data point can appear
//! twice in the subsample, and the matrix then has repeated rows, which is where the rank
//! deficiency the decomposition reports comes from.
//!
//! # Which `exp` the kernel uses depends on the caller
//!
//! `ContaminationSegmenter.SEGMENTATION_KERNEL` calls `FastMath.exp`, which is commons-math's own
//! and is [`jmath::fast_math::exp`] here. The `kernel-segmentation` dump calls `Math.exp`, so the
//! golden's matrices are the platform's. The kernel is a parameter for exactly that reason: this
//! module does not choose one, the caller does, and the two are not the same function.

use crate::java_random::JavaRandom;

/// `ContaminationSegmenter.SEGMENTATION_KERNEL_VARIANCE`.
pub const SEGMENTATION_KERNEL_VARIANCE: f64 = 0.025;

/// `ContaminationSegmenter.SEGMENTATION_KERNEL`, on alt fractions rather than on `PileupSummary`.
///
/// The reference folds each fraction to its minor allele fraction first, `min(af, 1 - af)`, so a
/// site at 0.8 and a site at 0.2 are the same point to the kernel. `FastMath.min` propagates NaN,
/// which `f64::min` does not, so the fold is written out.
pub fn segmentation_kernel(first: f64, second: f64) -> f64 {
    let maf1 = java_min(first, 1.0 - first);
    let maf2 = java_min(second, 1.0 - second);
    jmath::fast_math::exp(-square(maf1 - maf2) / (2.0 * SEGMENTATION_KERNEL_VARIANCE))
}

/// `MathUtils.square`.
fn square(value: f64) -> f64 {
    value * value
}

/// `FastMath.min`, which returns NaN when either argument is NaN.
fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        return f64::NAN;
    }
    if a <= b {
        a
    } else {
        b
    }
}

/// The subsampled kernel matrix, as `calculateReducedObservationMatrix` builds it.
///
/// `dimension` is `kernelApproximationDimension`. When it is at least the data size the reference
/// uses the data itself and draws nothing, which is why the generator is created before the branch
/// but not always used: a port that created it lazily would still agree here, but not if the
/// reference ever drew before the branch.
///
/// The matrix is filled the reference's way, lower triangle first with each value mirrored, then
/// the diagonal. That matters: the diagonal is `kernel(x, x)` computed separately rather than
/// assumed to be one, so a kernel that is not exactly one on the diagonal keeps its value.
pub fn sub_kernel_matrix<K>(data: &[f64], dimension: usize, kernel: K) -> Vec<Vec<f64>>
where
    K: Fn(f64, f64) -> f64,
{
    let mut rng = JavaRandom::new(1216);
    let num_subsample = dimension.min(data.len());
    let subsample: Vec<f64> = if num_subsample == data.len() {
        data.to_vec()
    } else {
        (0..num_subsample)
            .map(|_| data[rng.next_int_bound(data.len() as i32) as usize])
            .collect()
    };

    let mut matrix = vec![vec![0.0; num_subsample]; num_subsample];
    for i in 0..num_subsample {
        for j in 0..i {
            let value = kernel(subsample[i], subsample[j]);
            matrix[i][j] = value;
            matrix[j][i] = value;
        }
        matrix[i][i] = kernel(subsample[i], subsample[i]);
    }
    matrix
}

/// `KernelSegmenter.EPSILON`, which is not the decomposition's epsilon.
pub const EPSILON: f64 = 1e-10;

/// `ChangepointSortOrder`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChangepointSortOrder {
    /// Keep the order backward selection produced.
    BackwardSelection,
    /// Sort by increasing index.
    Index,
}

/// The N by p reduced observation matrix, `Z` in the paper.
///
/// This is where the decomposition is used, and it is used through
/// `1 / (sqrt(s) + EPSILON)`. A singular value of exactly zero becomes `1e10`, which is why
/// `docs/what-the-kernel-segmenter-needs-from-the-decomposition.md` exists: the columns that get
/// amplified the most are the ones the decomposition is least entitled to.
///
/// The two matrix walks are the reference's, and so is the multiplication order: commons-math
/// accumulates `sum += row[i] * column[i]` in index order, so the sum is written the same way here.
pub fn reduced_observation_matrix<K>(data: &[f64], dimension: usize, kernel: K) -> Vec<Vec<f64>>
where
    K: Fn(f64, f64) -> f64 + Copy,
{
    let num_subsample = dimension.min(data.len());
    // The reference draws the subsample inside this function, from a generator `findChangepoints`
    // made a moment earlier and uses nowhere else, so rebuilding it here draws the same points.
    let mut rng = JavaRandom::new(RANDOM_SEED);
    let subsample: Vec<f64> = if num_subsample == data.len() {
        data.to_vec()
    } else {
        (0..num_subsample)
            .map(|_| data[rng.next_int_bound(data.len() as i32) as usize])
            .collect()
    };

    let sub_kernel = sub_kernel_matrix(data, dimension, kernel);
    let svd = crate::singular_value_decomposition::SingularValueDecomposition::new(&sub_kernel);
    let inv_sqrt_singular_values: Vec<f64> = svd
        .singular_values
        .iter()
        .map(|value| 1.0 / (value.sqrt() + EPSILON))
        .collect();

    // `subKernelUMatrix`: U with each column scaled by its own inverse root singular value.
    let mut scaled_u = vec![vec![0.0; num_subsample]; num_subsample];
    for (i, row) in scaled_u.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = svd.u[i][j] * inv_sqrt_singular_values[j];
        }
    }

    // `reducedKernelMatrix`: every data point against every subsampled point.
    let mut reduced = vec![vec![0.0; num_subsample]; data.len()];
    for (i, row) in reduced.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = kernel(data[i], subsample[j]);
        }
    }

    let mut product = vec![vec![0.0; num_subsample]; data.len()];
    for row in 0..data.len() {
        for column in 0..num_subsample {
            let mut sum = 0.0;
            for i in 0..num_subsample {
                sum += reduced[row][i] * scaled_u[i][column];
            }
            product[row][column] = sum;
        }
    }
    product
}

/// `KernelSegmenter.RANDOM_SEED`.
pub const RANDOM_SEED: i64 = 1216;

/// `calculateKernelApproximationDiagonal`: the squared norm of each row.
///
/// The reference takes `MathUtils.square(getRowVector(i).getNorm())`, so it is the square of a
/// square root rather than a plain sum of squares, and the two are not the same double.
pub fn kernel_approximation_diagonal(reduced: &[Vec<f64>]) -> Vec<f64> {
    reduced
        .iter()
        .map(|row| {
            // `ArrayRealVector.getNorm`, which accumulates the squares in index order.
            let mut sum = 0.0;
            for value in row {
                sum += value * value;
            }
            let norm = sum.sqrt();
            norm * norm
        })
        .collect()
}

/// The quantities `calculateSegmentCost` carries between calls.
#[derive(Debug, Clone)]
struct Cost {
    /// The diagonal term.
    d: f64,
    /// The running row sum used for the off-diagonal terms.
    w: Vec<f64>,
    /// The off-diagonal accumulation.
    v: f64,
    /// The total cost.
    c: f64,
}

/// `calculateSegmentCost`, both indices inclusive, wrapping when `start > end`.
fn segment_cost(start: usize, end: usize, reduced: &[Vec<f64>], diagonal: &[f64]) -> Cost {
    let n = reduced.len();
    let p = if n == 0 { 0 } else { reduced[0].len() };

    let mut d = diagonal[start];
    let mut w = reduced[start].clone();
    let mut v: f64 = w.iter().map(|value| value * value).sum();

    // The reference wraps around the beginning of the data when the segment does.
    let indices: Vec<usize> = if start <= end {
        ((start + 1)..(end + 1)).collect()
    } else {
        ((start + 1)..n).chain(0..(end + 1)).collect()
    };

    for tau_prime in &indices {
        d += diagonal[*tau_prime];
        let mut z_dot_w = 0.0;
        for j in 0..p {
            z_dot_w += reduced[*tau_prime][j] * w[j];
            w[j] += reduced[*tau_prime][j];
        }
        v += 2.0 * z_dot_w + diagonal[*tau_prime];
    }
    let c = d - v / (indices.len() + 1) as f64;
    Cost { d, w, v, c }
}

/// `calculateWindowCosts`: the cost of a changepoint at each index, for one window size.
///
/// The segments wrap, so the first and last points are given costs against the other end of the
/// data rather than being skipped. `findChangepointCandidates` removes those two afterwards.
fn window_costs(reduced: &[Vec<f64>], diagonal: &[f64], window_size: usize) -> Vec<f64> {
    let n = reduced.len();
    let p = if n == 0 { 0 } else { reduced[0].len() };

    let mut center = 0usize;
    let mut start = (center + n - window_size + 1) % n;
    let mut end = (center + window_size) % n;

    let left = segment_cost(start, center, reduced, diagonal);
    let right = segment_cost(center + 1, end, reduced, diagonal);
    let total = segment_cost(start, end, reduced, diagonal);

    let (mut left_d, mut left_w, mut left_v, mut left_c) = (left.d, left.w, left.v, left.c);
    let (mut right_d, mut right_w, mut right_v, mut right_c) = (right.d, right.w, right.v, right.c);
    let (mut total_d, mut total_w, mut total_v, mut total_c) = (total.d, total.w, total.v, total.c);

    let mut costs = vec![0.0; n];
    costs[center] = left_c + right_c - total_c;

    let window_size_reciprocal = 1.0 / window_size as f64;

    // Slide the three segments one point at a time, updating by recurrence rather than recomputing.
    for c in 0..n {
        center = c;
        let center_next = (center + 1) % n;
        let end_next = (end + 1) % n;

        // The left segment loses its first point and gains the next centre.
        left_d -= diagonal[start];
        let mut z_dot_w = 0.0;
        for j in 0..p {
            z_dot_w += reduced[start][j] * left_w[j];
            left_w[j] -= reduced[start][j];
        }
        left_v += -2.0 * z_dot_w + diagonal[start];

        left_d += diagonal[center_next];
        z_dot_w = 0.0;
        for j in 0..p {
            z_dot_w += reduced[center_next][j] * left_w[j];
            left_w[j] += reduced[center_next][j];
        }
        left_v += 2.0 * z_dot_w + diagonal[center_next];

        left_c = left_d - left_v * window_size_reciprocal;

        // The right segment loses the next centre and gains the point past its end.
        right_d -= diagonal[center_next];
        z_dot_w = 0.0;
        for j in 0..p {
            z_dot_w += reduced[center_next][j] * right_w[j];
            right_w[j] -= reduced[center_next][j];
        }
        right_v += -2.0 * z_dot_w + diagonal[center_next];

        right_d += diagonal[end_next];
        z_dot_w = 0.0;
        for j in 0..p {
            z_dot_w += reduced[end_next][j] * right_w[j];
            right_w[j] += reduced[end_next][j];
        }
        right_v += 2.0 * z_dot_w + diagonal[end_next];

        right_c = right_d - right_v * window_size_reciprocal;

        // The total segment loses the first point and gains the point past the end. Note the half:
        // it spans two windows, so its reciprocal is halved rather than doubled.
        total_d -= diagonal[start];
        z_dot_w = 0.0;
        for j in 0..p {
            z_dot_w += reduced[start][j] * total_w[j];
            total_w[j] -= reduced[start][j];
        }
        total_v += -2.0 * z_dot_w + diagonal[start];

        total_d += diagonal[end_next];
        z_dot_w = 0.0;
        for j in 0..p {
            z_dot_w += reduced[end_next][j] * total_w[j];
            total_w[j] += reduced[end_next][j];
        }
        total_v += 2.0 * z_dot_w + diagonal[end_next];

        total_c = total_d - 0.5 * total_v * window_size_reciprocal;

        costs[center_next] = left_c + right_c - total_c;

        start = (start + 1) % n;
        end = end_next;
    }
    costs
}

/// `findChangepointCandidates`: the local minima of the window costs, over every window size.
///
/// A window size wider than half the data is skipped with a warning rather than being an error, so
/// a short contig simply contributes fewer candidates.
fn changepoint_candidates(
    data_size: usize,
    reduced: &[Vec<f64>],
    diagonal: &[f64],
    max_num_changepoints: usize,
    window_sizes: &[usize],
) -> Vec<usize> {
    let mut candidates = Vec::with_capacity(window_sizes.len() * max_num_changepoints);
    for window_size in window_sizes {
        if 2 * window_size > data_size {
            continue;
        }
        let costs = window_costs(reduced, diagonal, *window_size);
        let minima = crate::persistence_optimizer::persistence_optimizer(&costs)
            .expect("the window costs are never empty")
            .minima_indices;
        let mut minima: Vec<usize> = minima;
        // `List.remove(Object)` removes the first occurrence, and there is only ever one.
        if let Some(position) = minima.iter().position(|index| *index == 0) {
            minima.remove(position);
        }
        if let Some(position) = minima.iter().position(|index| *index == data_size - 1) {
            minima.remove(position);
        }
        candidates.extend(minima.into_iter().take(max_num_changepoints));
    }
    candidates
}

/// `calculateChangepointPenalty`: `A * C + B * C * log(N / (C + EPSILON))`.
///
/// The epsilon keeps the logarithm finite at zero changepoints, where the whole term is multiplied
/// by zero anyway; `Math.log` of a large number is the platform's.
fn changepoint_penalty(
    num_changepoints: usize,
    linear_factor: f64,
    log_linear_factor: f64,
    num_data: usize,
) -> f64 {
    linear_factor * num_changepoints as f64
        + log_linear_factor
            * num_changepoints as f64
            * (num_data as f64 / (num_changepoints as f64 + EPSILON)).ln()
}

/// One segment during backward selection.
#[derive(Debug, Clone)]
struct Segment {
    /// Inclusive start index.
    start: usize,
    /// Inclusive end index.
    end: usize,
    /// The segment's cost.
    cost: f64,
}

/// `selectChangepoints`: merge the cheapest adjacent pair until one segment is left.
///
/// The changepoints come out in the order they were merged away, which is increasing cost, and the
/// penalty then decides how many of them to keep.
fn select_changepoints(
    candidates: &[usize],
    max_num_changepoints: usize,
    linear_factor: f64,
    log_linear_factor: f64,
    reduced: &[Vec<f64>],
    diagonal: &[f64],
) -> Vec<usize> {
    let mut changepoints: Vec<usize> = Vec::with_capacity(candidates.len());
    let num_data = reduced.len();

    let penalties: Vec<f64> = (0..=max_num_changepoints)
        .map(|count| changepoint_penalty(count, linear_factor, log_linear_factor, num_data))
        .collect();

    // `sorted().distinct()`, then the starts shifted by one and clamped to the last index.
    let mut sorted: Vec<usize> = candidates.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let mut starts: Vec<usize> = sorted.iter().map(|i| (i + 1).min(num_data - 1)).collect();
    starts.insert(0, 0);
    let mut ends: Vec<usize> = sorted.clone();
    ends.push(num_data - 1);

    let num_segments = starts.len();
    let mut segments: Vec<Segment> = (0..num_segments)
        .map(|i| Segment {
            start: starts[i],
            end: ends[i],
            cost: segment_cost(starts[i], ends[i], reduced, diagonal).c,
        })
        .collect();

    let mut total_costs: Vec<f64> = vec![segments.iter().map(|s| s.cost).sum()];
    let mut costs_for_pairs: Vec<f64> = (0..num_segments - 1)
        .map(|i| segments[i].cost + segments[i + 1].cost)
        .collect();
    let mut costs_for_merged: Vec<f64> = (0..num_segments - 1)
        .map(|i| segment_cost(starts[i], ends[i + 1], reduced, diagonal).c)
        .collect();
    let mut costs_for_merging: Vec<f64> = (0..num_segments - 1)
        .map(|i| costs_for_pairs[i] - costs_for_merged[i])
        .collect();

    for _ in 0..(num_segments - 1) {
        // `Collections.max` keeps the first of equal values, and `indexOf` then finds that first
        // one, so ties go to the leftmost pair.
        let maximum = costs_for_merging
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, |best, value| {
                if crate::persistence_optimizer::java_compare(value, best)
                    == std::cmp::Ordering::Greater
                {
                    value
                } else {
                    best
                }
            });
        // `List.indexOf` compares with `Double.equals`, which is bit equality on the canonical
        // NaN rather than `==`: it holds for two NaNs and fails between the two zeros.
        let index = costs_for_merging
            .iter()
            .position(|value| value.to_bits() == maximum.to_bits())
            .expect("the maximum is in the list");

        let new_cost = costs_for_merged[index];
        let new_start = segments[index].start;
        let mergepoint = segments[index].end;
        let new_end = segments[index + 1].end;

        segments.remove(index);
        segments.remove(index);
        segments.insert(
            index,
            Segment {
                start: new_start,
                end: new_end,
                cost: new_cost,
            },
        );

        costs_for_pairs.remove(index);
        costs_for_merged.remove(index);
        costs_for_merging.remove(index);

        // The pair to the left of the merge, if there is one.
        if index > 0 {
            costs_for_pairs[index - 1] = segments[index - 1].cost + segments[index].cost;
            costs_for_merged[index - 1] =
                segment_cost(segments[index - 1].start, new_end, reduced, diagonal).c;
            costs_for_merging[index - 1] = costs_for_pairs[index - 1] - costs_for_merged[index - 1];
        }
        // And the pair to the right.
        if index < segments.len() - 1 {
            costs_for_pairs[index] = segments[index].cost + segments[index + 1].cost;
            costs_for_merged[index] =
                segment_cost(new_start, segments[index + 1].end, reduced, diagonal).c;
            costs_for_merging[index] = costs_for_pairs[index] - costs_for_merged[index];
        }

        total_costs.insert(0, segments.iter().map(|s| s.cost).sum());
        changepoints.insert(0, mergepoint);
    }

    // The penalty decides how many to keep.
    let effective_max = max_num_changepoints.min(changepoints.len());
    let with_penalties: Vec<f64> = (0..=effective_max)
        .map(|i| total_costs[i] + penalties[i])
        .collect();
    let minimum = with_penalties
        .iter()
        .copied()
        .fold(f64::INFINITY, |best, value| {
            if crate::persistence_optimizer::java_compare(value, best) == std::cmp::Ordering::Less {
                value
            } else {
                best
            }
        });
    let optimal = with_penalties
        .iter()
        .position(|value| value.to_bits() == minimum.to_bits())
        .expect("the minimum is in the list");

    changepoints.truncate(optimal);
    changepoints
}

/// `KernelSegmenter.findChangepoints`, on a series of doubles.
///
/// The eight arguments are the reference's eight, in its order, which is worth more to a reader
/// checking the transcription than a struct would be.
///
/// The argument checks the reference makes are the caller's business here: this port is reached
/// from `ContaminationSegmenter`, which passes constants. The two early returns are kept, because
/// they are the reference's answers rather than its refusals: no changepoints requested, or no
/// data, and the list is empty.
#[allow(clippy::too_many_arguments)]
pub fn find_changepoints<K>(
    data: &[f64],
    max_num_changepoints: usize,
    kernel: K,
    kernel_approximation_dimension: usize,
    window_sizes: &[usize],
    linear_factor: f64,
    log_linear_factor: f64,
    sort_order: ChangepointSortOrder,
) -> Vec<usize>
where
    K: Fn(f64, f64) -> f64 + Copy,
{
    if max_num_changepoints == 0 || data.is_empty() {
        return Vec::new();
    }

    let reduced = reduced_observation_matrix(data, kernel_approximation_dimension, kernel);
    let diagonal = kernel_approximation_diagonal(&reduced);
    let candidates = changepoint_candidates(
        data.len(),
        &reduced,
        &diagonal,
        max_num_changepoints,
        window_sizes,
    );
    // No early return when the candidate list is empty: the reference logs a warning and carries
    // on into backward selection, which then has one segment, merges nothing and keeps nothing.
    let mut selected = select_changepoints(
        &candidates,
        max_num_changepoints,
        linear_factor,
        log_linear_factor,
        &reduced,
        &diagonal,
    );
    // `BACKWARD_SELECTION` keeps the order backward selection produced, so only `INDEX` sorts.
    if sort_order == ChangepointSortOrder::Index {
        selected.sort_unstable();
    }
    selected
}
