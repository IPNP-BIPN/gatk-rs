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
