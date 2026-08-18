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
