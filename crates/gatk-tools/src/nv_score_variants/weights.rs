//! The two models' weights, pinned by hash.
//!
//! `tools/nvscorevariants/extract_weights.py` takes them out of the pickled modules in the GATK jar,
//! once. They are about 8 MB and live outside the repository; a file whose hash is not the one here
//! was not extracted from the release this port reproduces, and is refused rather than scored with.

use sha2::{Digest, Sha256};
use std::collections::HashMap;

/// The environment variable naming the directory that holds the two files.
pub const DIRECTORY_VARIABLE: &str = "GATK_RS_NVSCORE_WEIGHTS";

/// `large/nvscorevariants/1d_cnn_mix_train_full_bn.pt`, whose hash is f1d3d913…98fb0.
pub const ONE_D: (&str, &str) = (
    "nvscore-1d.weights",
    "89dcb0a22a103a1773af23d70c06386ef988fdb6520134dad5edf528acd8bd12",
);
/// `large/nvscorevariants/small_2d.pt`, whose hash is 2682dc9e…d12.
pub const TWO_D: (&str, &str) = (
    "nvscore-2d.weights",
    "2ee26bb9db2b58eabd942448f977d45dae7e600ce67534095b31694269faa441",
);

/// A file's tensors by name, flat.
pub struct Weights {
    tensors: HashMap<String, (Vec<usize>, Vec<f32>)>,
}

fn take_u32(bytes: &[u8], at: &mut usize) -> Option<u32> {
    let value = u32::from_le_bytes(bytes.get(*at..*at + 4)?.try_into().ok()?);
    *at += 4;
    Some(value)
}

fn take_u64(bytes: &[u8], at: &mut usize) -> Option<u64> {
    let value = u64::from_le_bytes(bytes.get(*at..*at + 8)?.try_into().ok()?);
    *at += 8;
    Some(value)
}

impl Weights {
    /// Read and check one of [`ONE_D`] or [`TWO_D`] from `directory`.
    pub fn load(directory: &str, (name, expected): (&str, &str)) -> Result<Weights, String> {
        let path = format!("{directory}/{name}");
        let bytes = std::fs::read(&path).map_err(|e| format!("{path}: {e}"))?;
        let digest: String = Sha256::digest(&bytes).iter().map(|b| format!("{b:02x}")).collect();
        if digest != expected {
            return Err(format!("{path} has SHA-256 {digest}, and the port is pinned to {expected}"));
        }
        Weights::parse(&bytes).ok_or_else(|| format!("{path} is not a weights file"))
    }

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
            let data = bytes
                .get(at..at + 4 * count)?
                .chunks_exact(4)
                .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
                .collect();
            at += 4 * count;
            tensors.insert(name, (shape, data));
        }
        Some(Weights { tensors })
    }

    /// A tensor's values, which the hash has already vouched for.
    pub fn get(&self, name: &str) -> &[f32] {
        &self.tensors[name].1
    }

    /// A tensor's shape.
    pub fn shape(&self, name: &str) -> &[usize] {
        &self.tensors[name].0
    }
}
