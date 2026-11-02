//! The arithmetic of the reference's forward pass, in the order its CPU backend does it.
//!
//! Each kernel here answers to a layer of the dump `tools/nvscorevariants/dump_layers.py` writes,
//! compared bit for bit. `docs/nvscorevariants.md` has the numbers and what they do not cover. The
//! order is the one a single-threaded run of the image takes on an AVX2 dispatch, and that is not
//! one order among equals by accident: `dense_2` moves its last bits with the thread count in the
//! reference itself.
//!
//! * **Convolutions** in the 1D model are oneDNN 3.1.1's `jit:avx2`: one accumulator for each output,
//!   starting from the bias, taking `vfmadd231ps` over input channel blocks of 8, then the kernel's
//!   height and width, then the channels of the block. In the 2D model they are `brgconv:avx2`: from
//!   zero over the kernel's height and width, then all the input channels, with the bias added last.
//! * **Dense layers** are MKL's `SGEMM`, whose source is not public, so its order was read off from
//!   outside by `tools/nvscorevariants/probe_gemm.py`: which terms have met before which others is
//!   observable by adding a huge number and two ones. Terms are taken 256 at a time, the last
//!   block and a bit split in two, and each block's sum, started from zero, is added to the output.
//! * **Batch norm** is `batch_norm_cpu_contiguous_impl`, where GCC's contraction fused both
//!   `bias - mean * alpha` and `x * alpha + beta`.
//! * **Softmax** is `_vec_softmax_lastdim` over SLEEF's `expf8_u10`.

/// `vfmadd231ps` on one lane: `a * b + c` rounded once.
///
/// The product of two floats is exact in a double. The sum rounds again there, which only matters
/// when that rounding lands exactly on a midpoint between two floats, and then the float it goes to
/// depends on which side the exact sum was. Everywhere else the double's rounding cannot have
/// crossed a midpoint, since every midpoint is itself a double.
#[inline(always)]
pub fn fused(a: f32, b: f32, c: f32) -> f32 {
    let product = f64::from(a) * f64::from(b);
    let addend = f64::from(c);
    let sum = product + addend;
    let bits = sum.to_bits();
    // The midpoint between two normal floats is 29 low bits of 1 followed by zeros; a result too
    // small to be a normal float has its midpoints elsewhere and goes the long way too.
    if bits & 0x1FFF_FFFF != 0x1000_0000 && sum.abs() >= 1.1754943508222875e-38 {
        return sum as f32;
    }
    rounded_once(product, addend, sum)
}

/// Round to odd in the double, which a later rounding to fewer bits cannot be fooled by.
#[cold]
fn rounded_once(product: f64, addend: f64, sum: f64) -> f32 {
    let from_product = sum - addend;
    let error = (product - from_product) + (addend - (sum - from_product));
    let bits = sum.to_bits();
    if error == 0.0 || bits & 1 == 1 {
        return sum as f32;
    }
    let outward = (error > 0.0) == (sum > 0.0);
    f64::from_bits(if outward { bits + 1 } else { bits - 1 }) as f32
}

/// The shape of an NCHW tensor.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Shape {
    pub batch: usize,
    pub channels: usize,
    pub height: usize,
    pub width: usize,
}

impl Shape {
    pub fn len(&self) -> usize {
        self.batch * self.channels * self.height * self.width
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A convolution's weights, `(output channels, input channels, height, width)`, and its bias.
pub struct ConvWeights<'a> {
    pub outputs: usize,
    pub inputs: usize,
    pub height: usize,
    pub width: usize,
    pub weight: &'a [f32],
    pub bias: &'a [f32],
    /// `brgconv:avx2`, which the 2D model gets, and not `jit:avx2`: from zero with every input
    /// channel at each kernel position, and the bias added at the end. Measured on 4000 positions
    /// of each of its layers, with no difference, by the layer test's first batch.
    pub bias_last: bool,
}

/// `jit_avx2_convolution_fwd_t` on a NCHW tensor, stride 1 and no padding.
///
/// The blocks are of 8 input channels unless there are fewer than 8 in all, which is oneDNN's flat
/// first convolution. A partial last block is padded with zeros in both the input and the weights,
/// which leaves every sum where it was except a negative zero. The threads share out outputs and
/// not terms, so no accumulator is ever split: this is the same number of threads as the reference
/// has or not, and the bits do not care.
pub fn conv2d(src: &[f32], input: Shape, conv: &ConvWeights, threads: usize) -> (Vec<f32>, Shape) {
    let out = Shape {
        batch: input.batch,
        channels: conv.outputs,
        height: input.height + 1 - conv.height,
        width: input.width + 1 - conv.width,
    };
    let mut dst = vec![0f32; out.len()];
    let block = if conv.inputs < 8 || conv.bias_last {
        conv.inputs
    } else {
        8
    };
    let padded = !conv.inputs.is_multiple_of(block);
    let plane = out.height * out.width;
    let each_channel = |batch_channel: usize, slab: &mut [f32]| {
        let batch = batch_channel / conv.outputs;
        let channel = batch_channel % conv.outputs;
        let image = &src[batch * input.channels * input.height * input.width..];
        let weights = &conv.weight[channel * conv.inputs * conv.height * conv.width..];
        for y in 0..out.height {
            for x in 0..out.width {
                let mut sum = if conv.bias_last {
                    0.0
                } else {
                    conv.bias[channel]
                };
                let mut first = 0;
                while first < conv.inputs {
                    let last = (first + block).min(conv.inputs);
                    for ky in 0..conv.height {
                        for kx in 0..conv.width {
                            let at_k = ky * conv.width + kx;
                            let at_i = (y + ky) * input.width + x + kx;
                            let span = input.height * input.width;
                            let span_w = conv.height * conv.width;
                            for c in first..last {
                                sum =
                                    fused(image[c * span + at_i], weights[c * span_w + at_k], sum);
                            }
                        }
                    }
                    first = last;
                }
                if padded {
                    sum += 0.0;
                }
                if conv.bias_last {
                    sum += conv.bias[channel];
                }
                slab[y * out.width + x] = sum;
            }
        }
    };
    let chunks: Vec<&mut [f32]> = dst.chunks_mut(plane).collect();
    run_over(chunks, threads, &each_channel);
    (dst, out)
}

/// Hand each piece to a thread, in order within a thread and without any sharing between them.
fn run_over<T: Send>(items: Vec<T>, threads: usize, work: &(dyn Fn(usize, T) + Sync)) {
    let threads = threads.max(1);
    if threads == 1 || items.len() < 2 {
        for (at, item) in items.into_iter().enumerate() {
            work(at, item);
        }
        return;
    }
    let mut buckets: Vec<Vec<(usize, T)>> = (0..threads).map(|_| Vec::new()).collect();
    for (at, item) in items.into_iter().enumerate() {
        buckets[at % threads].push((at, item));
    }
    std::thread::scope(|scope| {
        for bucket in buckets {
            scope.spawn(move || {
                for (at, item) in bucket {
                    work(at, item);
                }
            });
        }
    });
}

/// `clamp_min(0)`: what `ReLU` does, and a NaN stays one.
pub fn relu(values: &mut [f32]) {
    for value in values {
        if *value <= 0.0 {
            *value = 0.0;
        }
    }
}

/// `MaxPool2d` with a stride equal to its height and a width of one.
pub fn max_pool_height(src: &[f32], input: Shape, height: usize) -> (Vec<f32>, Shape) {
    let out = Shape {
        height: input.height / height,
        ..input
    };
    let mut dst = vec![0f32; out.len()];
    for b in 0..input.batch {
        for c in 0..input.channels {
            let image = (b * input.channels + c) * input.height * input.width;
            let target = (b * out.channels + c) * out.height * out.width;
            for y in 0..out.height {
                for x in 0..out.width {
                    let mut best = src[image + (y * height) * input.width + x];
                    for k in 1..height {
                        let value = src[image + (y * height + k) * input.width + x];
                        // The kernel's `(val > maxval) || isnan(val)`.
                        if value > best || value.is_nan() {
                            best = value;
                        }
                    }
                    dst[target + y * out.width + x] = best;
                }
            }
        }
    }
    (dst, out)
}

/// A batch norm's parameters, in evaluation mode.
pub struct BatchNorm<'a> {
    pub weight: &'a [f32],
    pub bias: &'a [f32],
    pub mean: &'a [f32],
    pub var: &'a [f32],
    pub eps: f64,
}

/// `batch_norm_cpu_contiguous_impl`.
pub fn batch_norm(src: &[f32], parameters: &BatchNorm) -> Vec<f32> {
    let channels = parameters.weight.len();
    let eps = parameters.eps as f32;
    let terms: Vec<(f32, f32)> = (0..channels)
        .map(|c| {
            let invstd = 1.0 / (parameters.var[c] + eps).sqrt();
            let alpha = invstd * parameters.weight[c];
            (alpha, fused(-parameters.mean[c], alpha, parameters.bias[c]))
        })
        .collect();
    src.chunks(channels)
        .flat_map(|row| {
            row.iter()
                .zip(&terms)
                .map(|(x, (alpha, beta))| fused(*x, *alpha, *beta))
                .collect::<Vec<_>>()
        })
        .collect()
}

/// How MKL's `SGEMM` sums the terms of one output within one block of `K`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Terms {
    /// One accumulator, in order.
    Chain,
    /// One accumulator for the even terms and one for the odd, added together, and a term left
    /// over when there is an odd number of them goes on afterwards.
    EvenOdd,
}

/// `SGEMM`'s block of `K`.
pub const GEMM_BLOCK: usize = 256;

/// How much of the `remaining` of `K` the next block takes: a whole block while two would still
/// fit, half of what is left when one and a bit are, so the tail is never one long and one short.
fn gemm_block(remaining: usize) -> usize {
    if remaining >= 2 * GEMM_BLOCK {
        GEMM_BLOCK
    } else if remaining > GEMM_BLOCK {
        remaining.div_ceil(2)
    } else {
        remaining
    }
}

fn block_sum(row: &[f32], column: &[f32], terms: Terms) -> f32 {
    match terms {
        Terms::Chain => row
            .iter()
            .zip(column)
            .fold(0f32, |sum, (a, b)| fused(*a, *b, sum)),
        Terms::EvenOdd => {
            let (mut even, mut odd) = (0f32, 0f32);
            let pairs = row.len() / 2;
            for t in 0..pairs {
                even = fused(row[2 * t], column[2 * t], even);
                odd = fused(row[2 * t + 1], column[2 * t + 1], odd);
            }
            let mut sum = even + odd;
            if row.len() % 2 == 1 {
                sum = fused(row[row.len() - 1], column[row.len() - 1], sum);
            }
            sum
        }
    }
}

/// `Linear` on `(batch, K)`: `at::addmm(bias, input, weight.t())`, so `SGEMM` with `M` outputs, `N`
/// rows of the batch and `beta` 1 over a copy of the bias. `terms(output, row)` is which sum that
/// position gets, which is a property of `M`, `N` and `K` only.
pub fn linear(
    src: &[f32],
    weight: &[f32],
    bias: &[f32],
    terms: &(dyn Fn(usize, usize) -> Terms + Sync),
    threads: usize,
) -> Vec<f32> {
    let outputs = bias.len();
    let k = weight.len() / outputs;
    let batch = src.len() / k;
    let mut dst = vec![0f32; batch * outputs];
    let each_row = |row: usize, out: &mut [f32]| {
        let column = &src[row * k..(row + 1) * k];
        for (o, value) in out.iter_mut().enumerate() {
            let line = &weight[o * k..(o + 1) * k];
            let kind = terms(o, row);
            let mut sum = bias[o];
            let mut start = 0;
            while start < k {
                let end = start + gemm_block(k - start);
                sum += block_sum(&line[start..end], &column[start..end], kind);
                start = end;
            }
            *value = sum;
        }
    };
    let chunks: Vec<&mut [f32]> = dst.chunks_mut(outputs).collect();
    run_over(chunks, threads, &each_row);
    dst
}

/// SLEEF's `xexpf`, on one lane of the AVX2 build, where `vmla` is `vfmadd`.
#[allow(clippy::excessive_precision, clippy::approx_constant)]
pub fn sleef_expf(d: f32) -> f32 {
    const R_LN2: f32 = 1.442_695_040_888_963_4;
    const L2U: f32 = 0.693_145_751_953_125;
    const L2L: f32 = 1.428_606_765_330_187_0e-6;
    // `vrint_vi2_vf`, `_mm256_cvtps_epi32`: to nearest, ties to even.
    let q = (d * R_LN2).round_ties_even();
    let integer = q as i32;
    let mut s = fused(q, -L2U, d);
    s = fused(q, -L2L, s);
    let mut u = 0.000_198_527_62_f32;
    u = fused(u, s, 0.001_393_043_6);
    u = fused(u, s, 0.008_333_361);
    u = fused(u, s, 0.041_666_485);
    u = fused(u, s, 0.166_666_67);
    u = fused(u, s, 0.5);
    u = 1.0 + fused(s * s, u, s);
    // `vldexp2`: two multiplications by powers of two, so a subnormal rounds in the second.
    let half = integer >> 1;
    u = u * power_of_two(half) * power_of_two(integer - half);
    if d < -104.0 {
        return 0.0;
    }
    if d > 100.0 {
        return f32::INFINITY;
    }
    u
}

/// `vpow2i`.
fn power_of_two(exponent: i32) -> f32 {
    f32::from_bits(((exponent.wrapping_add(0x7f)) as u32) << 23)
}

/// `_vec_softmax_lastdim` over four outputs.
pub fn softmax(logits: &[f32]) -> Vec<f32> {
    let mut out = Vec::with_capacity(logits.len());
    for row in logits.chunks(4) {
        let mut most = row[0];
        for value in &row[1..] {
            // `vec::maximum`: a NaN wins.
            if most.is_nan() || value.is_nan() {
                most = f32::NAN;
            } else if *value > most {
                most = *value;
            }
        }
        let exps: Vec<f32> = row.iter().map(|x| sleef_expf(*x - most)).collect();
        let sum = exps[0] + exps[1] + exps[2] + exps[3];
        let inverse = 1.0 / sum;
        out.extend(exps.iter().map(|x| *x * inverse));
    }
    out
}

/// `VariantLabel`.
const NOT_SNP: usize = 0;
const NOT_INDEL: usize = 1;
const SNP: usize = 2;
const INDEL: usize = 3;

/// `torch.log`, which is MKL's `vsLn` at `VML_HA`. Its source is not public and its table is not
/// reproduced: this is the correctly rounded log of the double, rounded to a float, and the distance
/// from the reference is what `docs/nvscorevariants.md` measures.
fn log(value: f32) -> f32 {
    if value.is_nan() {
        return value;
    }
    jmath::math::log(f64::from(value)) as f32
}

/// `predictions_to_score`, on the softmax and the encoded `VariantType`.
pub fn scores(probabilities: &[f32], types: &[i64]) -> Vec<f32> {
    const EPS: f32 = 1e-7;
    probabilities
        .chunks(4)
        .zip(types)
        .map(|(p, kind)| {
            let snp = log(EPS + p[SNP] / (p[NOT_SNP] + EPS));
            let indel = log(EPS + p[INDEL] / (p[NOT_INDEL] + EPS));
            match kind {
                0 => snp,
                1 => indel,
                // `torch.maximum`, where a NaN wins.
                _ if snp.is_nan() || indel.is_nan() => f32::NAN,
                _ => {
                    if indel > snp {
                        indel
                    } else {
                        snp
                    }
                }
            }
        })
        .collect()
}
