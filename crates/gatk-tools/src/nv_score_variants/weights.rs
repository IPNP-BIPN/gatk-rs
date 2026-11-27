//! The two models' weights, pinned by hash.
//!
//! `tools/nvscorevariants/extract_weights.py` takes them out of the pickled modules in the GATK
//! jar, once. They are about 8 MB and live outside the repository; a file whose hash is not the one
//! here was not extracted from the release this port reproduces, and is refused rather than scored
//! with.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing", and `?` after an `Option` returns `None` from the
//!   enclosing function at once. So every read below that would run past the end of the file makes
//!   `parse` answer `None` instead of throwing, as `ArrayIndexOutOfBoundsException` would in Java;
//! - `Result<T, String>` is "a value or an error message", and `?` after one returns the error;
//! - `&mut usize` is a cursor the callee may move: `take_u32` reads at `*at` and advances it;
//! - `bytes.get(a..b)` is a range lookup that answers `None` rather than failing when the range
//!   does not fit.

use sha2::{Digest, Sha256};
use std::collections::HashMap;

/// The environment variable naming the directory that holds the two files.
pub const DIRECTORY_VARIABLE: &str = "GATK_RS_NVSCORE_WEIGHTS";

/// The file name and the SHA-256 of the extracted file, for the 1D model taken out of
/// `large/nvscorevariants/1d_cnn_mix_train_full_bn.pt`, whose hash is f1d3d913…98fb0.
pub const ONE_D: (&str, &str) = (
    "nvscore-1d.weights",
    "89dcb0a22a103a1773af23d70c06386ef988fdb6520134dad5edf528acd8bd12",
);
/// The file name and the SHA-256 of the extracted file, for the 2D model taken out of
/// `large/nvscorevariants/small_2d.pt`, whose hash is 2682dc9e…d12.
pub const TWO_D: (&str, &str) = (
    "nvscore-2d.weights",
    "2ee26bb9db2b58eabd942448f977d45dae7e600ce67534095b31694269faa441",
);

/// A file's tensors by name, flat.
///
/// Each value is the tensor's shape and its values in row-major order, the order
/// `extract_weights.py` writes them in, so a convolution's weight at
/// `(output, input, height, width)` is one flat index computed from the four.
pub struct Weights {
    tensors: HashMap<String, (Vec<usize>, Vec<f32>)>,
}

/// A little-endian `u32` at `*at`, moving the cursor past it, or `None` when the file is short.
fn take_u32(bytes: &[u8], at: &mut usize) -> Option<u32> {
    let value = u32::from_le_bytes(bytes.get(*at..*at + 4)?.try_into().ok()?);
    *at += 4;
    Some(value)
}

/// A little-endian `u64` at `*at`, as `take_u32` does for four bytes.
fn take_u64(bytes: &[u8], at: &mut usize) -> Option<u64> {
    let value = u64::from_le_bytes(bytes.get(*at..*at + 8)?.try_into().ok()?);
    *at += 8;
    Some(value)
}

impl Weights {
    /// Read and check one of [`ONE_D`] or [`TWO_D`] from `directory`.
    ///
    /// The hash is checked before anything is parsed, so a file that is not the pinned one is
    /// refused even when it would parse.
    pub fn load(directory: &str, (name, expected): (&str, &str)) -> Result<Weights, String> {
        let path = format!("{directory}/{name}");
        let bytes = std::fs::read(&path).map_err(|e| format!("{path}: {e}"))?;
        // `.iter().map(..).collect()` is Java's stream, map, collect: the digest's bytes as lower
        // case hexadecimal, joined into one string.
        let digest: String = Sha256::digest(&bytes)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        if digest != expected {
            return Err(format!(
                "{path} has SHA-256 {digest}, and the port is pinned to {expected}"
            ));
        }
        Weights::parse(&bytes).ok_or_else(|| format!("{path} is not a weights file"))
    }

    /// The layout `extract_weights.py` writes: the six bytes `NVSW1` and a newline, a `u32` count
    /// of tensors, then for each a `u32` name length, the UTF-8 name, a `u32` number of dimensions,
    /// one `u64` per dimension, and the product of the dimensions in little-endian `f32`s.
    fn parse(bytes: &[u8]) -> Option<Weights> {
        if bytes.get(..6)? != b"NVSW1\n" {
            return None;
        }
        let mut at = 6;
        let mut tensors = HashMap::new();
        for _ in 0..take_u32(bytes, &mut at)? {
            let length = take_u32(bytes, &mut at)? as usize;
            let name = String::from_utf8(bytes.get(at..at + length)?.to_vec()).ok()?;
            at += length;
            let dimensions = take_u32(bytes, &mut at)?;
            let mut shape = Vec::new();
            for _ in 0..dimensions {
                shape.push(take_u64(bytes, &mut at)? as usize);
            }
            let count: usize = shape.iter().product();
            // `as_chunks::<4>()` views the bytes as four-byte groups; `.0` is the whole groups and
            // drops the remainder, which the range above has already made empty.
            let data = bytes
                .get(at..at + 4 * count)?
                .as_chunks::<4>()
                .0
                .iter()
                .map(|b| f32::from_le_bytes(*b))
                .collect();
            at += 4 * count;
            tensors.insert(name, (shape, data));
        }
        Some(Weights { tensors })
    }

    /// A tensor's values, which the hash has already vouched for.
    ///
    /// `self.tensors[name]` panics on a name that is not there, rather than answering `None`: the
    /// names are the pickled modules' own, fixed by the hash, so a missing one is a bug here.
    pub fn get(&self, name: &str) -> &[f32] {
        &self.tensors[name].1
    }

    /// A tensor's shape.
    pub fn shape(&self, name: &str) -> &[usize] {
        &self.tensors[name].0
    }
}
