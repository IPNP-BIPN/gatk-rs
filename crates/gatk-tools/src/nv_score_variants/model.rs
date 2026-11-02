//! `scorevariants/models/gatk_1d.py` and `gatk_2d.py`, as the pickled modules in the jar build them.
//!
//! The dropouts are in evaluation mode and do nothing, so they are not here. Every step stores what
//! it produced under the name its module has in the reference, for the comparison that
//! `tools/nvscorevariants/dump_layers.py` is the other half of.

use super::kernels::{
    batch_norm, conv2d, linear, max_pool_height, relu, BatchNorm, ConvWeights, Shape, Terms,
};
use super::weights::Weights;
use std::collections::HashMap;

/// What the layers produced, by module name. Empty unless the caller asks for them.
#[derive(Default)]
pub struct Trace {
    pub keep: bool,
    pub layers: HashMap<String, Vec<f32>>,
}

impl Trace {
    fn store(&mut self, name: &str, values: &[f32]) {
        if self.keep {
            self.layers.insert(name.to_string(), values.to_vec());
        }
    }
}

fn convolution(
    weights: &Weights,
    prefix: &str,
    src: &[f32],
    shape: Shape,
    threads: usize,
) -> (Vec<f32>, Shape) {
    convolve(weights, prefix, src, shape, threads, false)
}

/// The 2D model's convolutions, which oneDNN gives to a different implementation.
fn brg_convolution(
    weights: &Weights,
    prefix: &str,
    src: &[f32],
    shape: Shape,
    threads: usize,
) -> (Vec<f32>, Shape) {
    convolve(weights, prefix, src, shape, threads, true)
}

fn convolve(
    weights: &Weights,
    prefix: &str,
    src: &[f32],
    shape: Shape,
    threads: usize,
    bias_last: bool,
) -> (Vec<f32>, Shape) {
    let weight_name = format!("{prefix}.weight");
    let dimensions = weights.shape(&weight_name);
    // A `Conv1d`'s weights have no height, and it runs as a convolution of height one.
    let (height, width) = match dimensions[2..] {
        [height, width] => (height, width),
        _ => (1, dimensions[2]),
    };
    let parameters = ConvWeights {
        outputs: dimensions[0],
        inputs: dimensions[1],
        height,
        width,
        weight: weights.get(&weight_name),
        bias: weights.get(&format!("{prefix}.bias")),
        bias_last,
    };
    conv2d(src, shape, &parameters, threads)
}

fn normalize(weights: &Weights, prefix: &str, src: &[f32], eps: f64) -> Vec<f32> {
    batch_norm(
        src,
        &BatchNorm {
            weight: weights.get(&format!("{prefix}.weight")),
            bias: weights.get(&format!("{prefix}.bias")),
            mean: weights.get(&format!("{prefix}.running_mean")),
            var: weights.get(&format!("{prefix}.running_var")),
            eps,
        },
    )
}

fn dense(
    weights: &Weights,
    prefix: &str,
    src: &[f32],
    threads: usize,
    terms: &(dyn Fn(usize, usize) -> Terms + Sync),
) -> Vec<f32> {
    linear(
        src,
        weights.get(&format!("{prefix}.weight")),
        weights.get(&format!("{prefix}.bias")),
        terms,
        threads,
    )
}

/// `(B, H, W, C)` as `(B, W, C)` flattened, which is what the two models' `permute` and `flatten`
/// come to.
fn channels_last(src: &[f32], shape: Shape) -> Vec<f32> {
    let mut out = vec![0f32; src.len()];
    let mut at = 0;
    for b in 0..shape.batch {
        let origin = b * shape.channels * shape.height * shape.width;
        for h in 0..shape.height {
            for w in 0..shape.width {
                for c in 0..shape.channels {
                    let plane = origin + c * shape.height * shape.width;
                    out[at] = src[plane + h * shape.width + w];
                    at += 1;
                }
            }
        }
    }
    out
}

/// Each row of `a` followed by the same row of `b`.
fn concatenate(a: &[f32], a_row: usize, b: &[f32], b_row: usize) -> Vec<f32> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    for (x, y) in a.chunks(a_row).zip(b.chunks(b_row)) {
        out.extend_from_slice(x);
        out.extend_from_slice(y);
    }
    out
}

const BEST_PRACTICES: usize = 7;

/// `GATK1D.forward`: the logits of `reference`, `(B, 128, 4)`, and `best_practices`, `(B, 7)`.
pub fn forward_1d(
    weights: &Weights,
    reference: &[f32],
    best_practices: &[f32],
    threads: usize,
    trace: &mut Trace,
) -> Vec<f32> {
    let batch = best_practices.len() / BEST_PRACTICES;
    let width = reference.len() / batch / 4;
    let shape = Shape {
        batch,
        channels: 4,
        height: 1,
        width,
    };
    let mut planar = vec![0f32; reference.len()];
    for b in 0..batch {
        for w in 0..width {
            for c in 0..4 {
                planar[(b * 4 + c) * width + w] = reference[(b * width + w) * 4 + c];
            }
        }
    }
    let (mut x, shape) = convolution(weights, "conv.0", &planar, shape, threads);
    trace.store("conv.0", &x);
    relu(&mut x);
    let (mut x, shape) = convolution(weights, "conv.2", &x, shape, threads);
    trace.store("conv.2", &x);
    relu(&mut x);
    let (mut x, shape) = convolution(weights, "conv.5", &x, shape, threads);
    trace.store("conv.5", &x);
    relu(&mut x);
    let flattened = channels_last(&x, shape);
    let row = flattened.len() / batch;

    let normalized = normalize(weights, "normalize", best_practices, 0.001);
    trace.store("normalize", &normalized);
    let classes = |outputs: usize, k: usize| {
        move |output: usize, column: usize| gemm_terms(outputs, batch, k, output, column)
    };
    let mut first = dense(
        weights,
        "dense_1.0",
        &normalized,
        threads,
        &classes(40, BEST_PRACTICES),
    );
    trace.store("dense_1.0", &first);
    relu(&mut first);

    let joined = concatenate(&flattened, row, &first, 40);
    let mut second = dense(
        weights,
        "dense_2.0",
        &joined,
        threads,
        &classes(40, row + 40),
    );
    trace.store("dense_2.0", &second);
    relu(&mut second);

    let last = concatenate(&second, 40, &normalized, BEST_PRACTICES);
    let logits = dense(
        weights,
        "dense_3",
        &last,
        threads,
        &classes(4, 40 + BEST_PRACTICES),
    );
    trace.store("dense_3", &logits);
    logits
}

/// `GATK2D.forward`: the logits of `read_tensor`, `(B, 128, 128, 15)`, and `best_practices`.
pub fn forward_2d(
    weights: &Weights,
    read_tensor: &[f32],
    best_practices: &[f32],
    threads: usize,
    trace: &mut Trace,
) -> Vec<f32> {
    let batch = best_practices.len() / BEST_PRACTICES;
    let (height, width, channels) = (128, 128, 15);
    let mut planar = vec![0f32; read_tensor.len()];
    for b in 0..batch {
        for h in 0..height {
            for w in 0..width {
                for c in 0..channels {
                    let to = (((b * channels + c) * height) + h) * width + w;
                    planar[to] = read_tensor[((b * height + h) * width + w) * channels + c];
                }
            }
        }
    }
    let shape = Shape {
        batch,
        channels,
        height,
        width,
    };
    let (mut x, shape) = brg_convolution(weights, "conv.0", &planar, shape, threads);
    trace.store("conv.0", &x);
    relu(&mut x);
    let (mut x, shape) = brg_convolution(weights, "conv.3", &x, shape, threads);
    trace.store("conv.3", &x);
    relu(&mut x);
    let (mut x, shape) = brg_convolution(weights, "conv.6", &x, shape, threads);
    trace.store("conv.6", &x);
    relu(&mut x);
    let (x, shape) = max_pool_height(&x, shape, 3);
    let (mut x, shape) = brg_convolution(weights, "conv.10", &x, shape, threads);
    trace.store("conv.10", &x);
    relu(&mut x);
    let (x, shape) = max_pool_height(&x, shape, 3);
    let flattened = channels_last(&x, shape);
    let row = flattened.len() / batch;

    let classes = |outputs: usize, k: usize| {
        move |output: usize, column: usize| gemm_terms(outputs, batch, k, output, column)
    };
    let normalized = normalize(weights, "dense_1.0", best_practices, 0.001);
    trace.store("dense_1.0", &normalized);
    let mut first = dense(
        weights,
        "dense_1.1",
        &normalized,
        threads,
        &classes(64, BEST_PRACTICES),
    );
    trace.store("dense_1.1", &first);
    relu(&mut first);

    let joined = concatenate(&flattened, row, &first, 64);
    let mut second = dense(
        weights,
        "dense_2.0",
        &joined,
        threads,
        &classes(24, row + 64),
    );
    trace.store("dense_2.0", &second);
    relu(&mut second);
    let logits = dense(weights, "dense_2.3", &second, threads, &classes(4, 24));
    trace.store("dense_2.3", &logits);
    logits
}

/// Which sum `SGEMM` gives a position, as `tools/nvscorevariants/probe_gemm.py` maps it under AVX2
/// for the sizes these two models use: `outputs` of 4, 24, 40 and 64, with `batch` from 1 to 69 and
/// 127 to 129. Outputs in sixteens get the one-accumulator kernel over batch rows in sixes, and a
/// tail of four or five. Everything else gets the even and odd one.
///
/// A batch of one goes another way through MKL that has not been decoded; `Chain` is the nearest and
/// `docs/nvscorevariants.md` gives its distance.
pub fn gemm_terms(outputs: usize, batch: usize, _k: usize, output: usize, column: usize) -> Terms {
    if batch == 1 {
        return Terms::Chain;
    }
    let in_sixteens = output < outputs - outputs % 16;
    let in_sixes = column < batch - batch % 6 || batch % 6 >= 4;
    // Measured, and not understood: from 56 rows of a batch, 64 outputs over 7 terms are all chain.
    if in_sixteens && (in_sixes || (outputs == 64 && batch >= 56)) {
        Terms::Chain
    } else {
        Terms::EvenOdd
    }
}
