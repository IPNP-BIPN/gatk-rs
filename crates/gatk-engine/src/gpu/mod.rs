//! Milestone GPU: a second implementation of a path that is already byte-identical, which ships
//! only while it produces the CPU path's bytes. docs/gpu/eligibility.md says what may run here.
//!
//! The device path is taken when the crate is built with the `gpu` feature, the host has a Metal
//! device (macOS) or a CUDA one (Linux), and `GATK_RS_GPU` is set to anything but `0`. Otherwise every entry point answers
//! `None` and the caller runs its CPU code, which is the reference either way.

#[cfg(all(feature = "gpu", target_os = "linux"))]
mod cuda;
#[cfg(all(feature = "gpu", target_os = "macos"))]
mod metal;
#[cfg(all(feature = "gpu", target_os = "linux"))]
use cuda as backend;
#[cfg(all(feature = "gpu", target_os = "macos"))]
use metal as backend;

/// What `read_likelihood_given_haplotype_log10` reads off one read.
pub struct ReadInputs<'a> {
    pub bases: &'a [u8],
    pub qualities: &'a [u8],
    pub insertion: &'a [u8],
    pub deletion: &'a [u8],
    pub gap_continuation: &'a [u8],
}

/// Whether the environment asks for the device path.
pub fn requested() -> bool {
    std::env::var_os("GATK_RS_GPU").is_some_and(|v| !v.is_empty() && v != "0")
}

/// The device the path would run on, when it is requested and present.
pub fn device_name() -> Option<String> {
    #[cfg(all(feature = "gpu", any(target_os = "macos", target_os = "linux")))]
    if requested() {
        return backend::device_name();
    }
    None
}

/// `read_likelihood_given_haplotype_log10` for every read against every haplotype, read-major
/// (`reads.len() * haplotypes.len()` values), or `None` when the device path is not taken.
pub fn pair_hmm_log10(haplotypes: &[Vec<u8>], reads: &[ReadInputs]) -> Option<Vec<f64>> {
    #[cfg(all(feature = "gpu", any(target_os = "macos", target_os = "linux")))]
    if requested() {
        return backend::pair_hmm_log10(haplotypes, reads);
    }
    let _ = (haplotypes, reads);
    None
}

/// The kernel's own binary64 addition (`operation` 0) or multiplication (1) over IEEE bit
/// patterns, so a test can hold it to the host's. `None` when the device path is not taken.
#[doc(hidden)]
pub fn f64_ops(operation: u32, a: &[u64], b: &[u64]) -> Option<Vec<u64>> {
    #[cfg(all(feature = "gpu", any(target_os = "macos", target_os = "linux")))]
    if requested() {
        return backend::f64_ops(operation, a, b);
    }
    let _ = (operation, a, b);
    None
}

/// A convolution for [`conv2d_fused`]: an NCHW input, `(outputs, inputs, kernel_height,
/// kernel_width)` weights, stride 1, no padding, in `kernels::conv2d`'s order of NVScoreVariants.
pub struct ConvShape {
    pub batch: usize,
    pub inputs: usize,
    pub height: usize,
    pub width: usize,
    pub outputs: usize,
    pub kernel_height: usize,
    pub kernel_width: usize,
    /// Input channels per block of the sum.
    pub block: usize,
    /// The sum starts from zero and the bias is added last, rather than starting from the bias.
    pub bias_last: bool,
    /// The last block is partial, so the sum gets the `+ 0.0` its zero padding gives it.
    pub padded: bool,
}

/// NVScoreVariants' convolution on the device, or `None` when the device path is not taken.
pub fn conv2d_fused(
    src: &[f32],
    weight: &[f32],
    bias: &[f32],
    shape: &ConvShape,
) -> Option<Vec<f32>> {
    #[cfg(all(feature = "gpu", any(target_os = "macos", target_os = "linux")))]
    if requested() {
        return backend::conv2d_fused(src, weight, bias, shape);
    }
    let _ = (src, weight, bias, shape);
    None
}

/// The kernel's own `fused(a, b, c)` (`operation` 0) or binary32 addition of `a` and `b` (1) over
/// IEEE bit patterns, so a test can hold it to the host's. `None` when the device path is not taken.
#[doc(hidden)]
pub fn f32_ops(operation: u32, a: &[u32], b: &[u32], c: &[u32]) -> Option<Vec<u32>> {
    #[cfg(all(feature = "gpu", any(target_os = "macos", target_os = "linux")))]
    if requested() {
        return backend::f32_ops(operation, a, b, c);
    }
    let _ = (operation, a, b, c);
    None
}

/// `pair_hmm`'s inputs, packed once for either backend. Everything that reaches a transcendental
/// stays on the host, in the CPU path's own code: the transitions, the priors and the initial
/// condition go down as bits, the sum comes back as bits, and the logarithm is taken here.
#[cfg(all(feature = "gpu", any(target_os = "macos", target_os = "linux")))]
struct PairHmmBatch {
    read_bases: Vec<u8>,
    read_spans: Vec<[u32; 2]>,
    read_rows: Vec<u64>,
    haplotype_bases: Vec<u8>,
    haplotype_spans: Vec<[u32; 2]>,
    haplotype_initial: Vec<u64>,
    jobs: Vec<[u32; 2]>,
    /// The words of scratch one pair needs: three rows of the longest haplotype plus one.
    stride: usize,
}

#[cfg(all(feature = "gpu", any(target_os = "macos", target_os = "linux")))]
impl PairHmmBatch {
    fn pack(haplotypes: &[Vec<u8>], reads: &[ReadInputs]) -> Self {
        use crate::pair_hmm::{initial_condition, prior, qual_to_trans_probs};
        let mut read_bases = Vec::new();
        let mut read_spans = Vec::with_capacity(reads.len());
        let mut read_rows = Vec::new();
        for read in reads {
            read_spans.push([read_bases.len() as u32, read.bases.len() as u32]);
            read_bases.extend_from_slice(read.bases);
            for k in 0..read.bases.len() {
                let t = qual_to_trans_probs(
                    read.insertion[k],
                    read.deletion[k],
                    read.gap_continuation[k],
                );
                read_rows.extend(t.iter().map(|v| v.to_bits()));
                // The two values `prior` can take for this base: the kernel picks one per cell.
                read_rows.push(prior(0, 0, read.qualities[k]).to_bits());
                read_rows.push(prior(0, 1, read.qualities[k]).to_bits());
            }
        }
        let mut haplotype_bases = Vec::new();
        let mut haplotype_spans = Vec::with_capacity(haplotypes.len());
        let mut haplotype_initial = Vec::with_capacity(haplotypes.len());
        for haplotype in haplotypes {
            haplotype_spans.push([haplotype_bases.len() as u32, haplotype.len() as u32]);
            haplotype_bases.extend_from_slice(haplotype);
            haplotype_initial.push((initial_condition() / haplotype.len() as f64).to_bits());
        }
        let jobs = (0..reads.len() as u32)
            .flat_map(|r| (0..haplotypes.len() as u32).map(move |h| [r, h]))
            .collect();
        let longest = haplotypes.iter().map(Vec::len).max().unwrap_or(0);
        PairHmmBatch {
            read_bases,
            read_spans,
            read_rows,
            haplotype_bases,
            haplotype_spans,
            haplotype_initial,
            jobs,
            stride: 3 * (longest + 1),
        }
    }

    fn log10(totals: Vec<u64>) -> Vec<f64> {
        let offset = crate::pair_hmm::initial_condition_log10();
        totals
            .into_iter()
            .map(|bits| jmath::math::log10(f64::from_bits(bits)) - offset)
            .collect()
    }
}

/// How many outputs a [`conv2d_fused`] writes.
#[cfg(all(feature = "gpu", any(target_os = "macos", target_os = "linux")))]
fn conv_outputs(shape: &ConvShape) -> usize {
    shape.batch
        * shape.outputs
        * (shape.height + 1 - shape.kernel_height)
        * (shape.width + 1 - shape.kernel_width)
}

/// The kernel's `ConvShape` words for the `count` outputs from `first`, and `count` as a word.
#[cfg(all(feature = "gpu", any(target_os = "macos", target_os = "linux")))]
fn conv_words(shape: &ConvShape, first: usize, count: usize) -> ([u32; 11], u32) {
    let word = |v: usize| u32::try_from(v).expect("a convolution past 2^32 elements");
    let words = [
        shape.batch,
        shape.inputs,
        shape.height,
        shape.width,
        shape.outputs,
        shape.kernel_height,
        shape.kernel_width,
        shape.block,
        usize::from(shape.bias_last),
        usize::from(shape.padded),
        first,
    ]
    .map(word);
    (words, word(count))
}
