//! Milestone GPU: a second implementation of a path that is already byte-identical, which ships
//! only while it produces the CPU path's bytes. docs/gpu/eligibility.md says what may run here.
//!
//! The device path is taken when the crate is built with the `gpu` feature, the host has a Metal
//! device, and `GATK_RS_GPU` is set to anything but `0`. Otherwise every entry point answers
//! `None` and the caller runs its CPU code, which is the reference either way.

#[cfg(all(feature = "gpu", target_os = "macos"))]
mod metal;

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
    #[cfg(all(feature = "gpu", target_os = "macos"))]
    if requested() {
        return metal::device_name();
    }
    None
}

/// `read_likelihood_given_haplotype_log10` for every read against every haplotype, read-major
/// (`reads.len() * haplotypes.len()` values), or `None` when the device path is not taken.
pub fn pair_hmm_log10(haplotypes: &[Vec<u8>], reads: &[ReadInputs]) -> Option<Vec<f64>> {
    #[cfg(all(feature = "gpu", target_os = "macos"))]
    if requested() {
        return metal::pair_hmm_log10(haplotypes, reads);
    }
    let _ = (haplotypes, reads);
    None
}

/// The kernel's own binary64 addition (`operation` 0) or multiplication (1) over IEEE bit
/// patterns, so a test can hold it to the host's. `None` when the device path is not taken.
#[doc(hidden)]
pub fn f64_ops(operation: u32, a: &[u64], b: &[u64]) -> Option<Vec<u64>> {
    #[cfg(all(feature = "gpu", target_os = "macos"))]
    if requested() {
        return metal::f64_ops(operation, a, b);
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
    #[cfg(all(feature = "gpu", target_os = "macos"))]
    if requested() {
        return metal::conv2d_fused(src, weight, bias, shape);
    }
    let _ = (src, weight, bias, shape);
    None
}

/// The kernel's own `fused(a, b, c)` (`operation` 0) or binary32 addition of `a` and `b` (1) over
/// IEEE bit patterns, so a test can hold it to the host's. `None` when the device path is not taken.
#[doc(hidden)]
pub fn f32_ops(operation: u32, a: &[u32], b: &[u32], c: &[u32]) -> Option<Vec<u32>> {
    #[cfg(all(feature = "gpu", target_os = "macos"))]
    if requested() {
        return metal::f32_ops(operation, a, b, c);
    }
    let _ = (operation, a, b, c);
    None
}
