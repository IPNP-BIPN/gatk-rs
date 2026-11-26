//! The groups and the double and string datasets `HDF5File` reads back, read through hdf5-pure.
//!
//! Two kinds of file reach it. The reference's are written by the HDF5 C library in its default
//! layout: a version 0 superblock, version 1 object headers, groups held in a B-tree of symbol
//! table nodes over a local heap. The port's are [`crate::hdf5_writer`]'s, in the "latest" layout.
//! hdf5-pure reads both and checks each object header's checksum where the layout has one
//! (docs/hdf5-second-reader.md).
//!
//! The reference reads every dataset here with `HDF5File.readDoubleMatrix`, `readDoubleArray` or
//! `readStringArray`, which ask the library for doubles or strings. A dataset of another class
//! reads as absent, as in the hand-written reader this replaced, rather than as a converted value.

/// The 8 bytes every HDF5 file begins with, which `IOUtils.isHDF5File` compares.
pub const SIGNATURE: &[u8] = b"\x89HDF\r\n\x1a\n";

/// `IOUtils.isHDF5File`: whether the file's first 8 bytes are the signature. The reference decides
/// by the bytes rather than the name, so a counts file named `.tsv` that holds HDF5 is read as
/// HDF5, and one named `.hdf5` that holds text is read as text.
pub fn is_hdf5_file(path: &str) -> bool {
    use std::io::Read;
    let mut head = [0u8; 8];
    std::fs::File::open(path)
        .and_then(|mut file| file.read_exact(&mut head))
        .is_ok_and(|()| head == SIGNATURE)
}

/// A file's bytes, opened.
pub struct Hdf5Reader {
    file: hdf5_pure::File,
}

impl std::fmt::Debug for Hdf5Reader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Hdf5Reader")
    }
}

impl Hdf5Reader {
    /// The file, when the bytes are one hdf5-pure can open.
    pub fn open(bytes: Vec<u8>) -> Option<Self> {
        if bytes.get(..8)? != SIGNATURE {
            return None;
        }
        hdf5_pure::File::from_bytes(bytes)
            .ok()
            .map(|file| Self { file })
    }

    /// The doubles at `path` in row-major order, with their shape.
    pub fn doubles(&self, path: &str) -> Option<(Vec<u64>, Vec<f64>)> {
        let dataset = self.file.dataset(path).ok()?;
        if dataset.dtype().ok()? != hdf5_pure::DType::F64 {
            return None;
        }
        Some((dataset.shape().ok()?, dataset.read_f64().ok()?))
    }

    /// The strings at `path`, fixed-length or variable-length, whichever the library chose.
    pub fn strings(&self, path: &str) -> Option<Vec<String>> {
        let dataset = self.file.dataset(path).ok()?;
        match dataset.dtype().ok()? {
            hdf5_pure::DType::String | hdf5_pure::DType::VariableLengthString => {
                dataset.read_string().ok()
            }
            _ => None,
        }
    }

    /// Whether a dataset is at `path`: `HDF5File.isPresent`.
    pub fn is_present(&self, path: &str) -> bool {
        self.file.dataset(path).is_ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hdf5_writer::Hdf5File;

    #[test]
    fn what_the_writer_wrote_reads_back() {
        let mut file = Hdf5File::new();
        file.double_array("/labels/snp", &[1.0, 0.0]);
        file.double_matrix(
            "/annotations/chunk_0",
            &[vec![1.5, f64::NAN], vec![3.0, 4.0]],
        );
        file.string_array("/annotations/names", &["QD".to_string(), "MQ".to_string()]);
        let reader = Hdf5Reader::open(file.to_bytes()).expect("an HDF5 file");
        assert_eq!(
            reader.doubles("/labels/snp"),
            Some((vec![2], vec![1.0, 0.0]))
        );
        let (shape, values) = reader.doubles("/annotations/chunk_0").expect("a matrix");
        assert_eq!(shape, vec![2, 2]);
        assert!(values[1].is_nan());
        assert_eq!(values[2], 3.0);
        assert_eq!(
            reader.strings("/annotations/names"),
            Some(vec!["QD".to_string(), "MQ".to_string()])
        );
        assert_eq!(reader.doubles("/labels/training"), None);
        // A string dataset is not doubles, and doubles are not strings.
        assert_eq!(reader.doubles("/annotations/names"), None);
        assert_eq!(reader.strings("/labels/snp"), None);
        assert!(reader.is_present("/labels/snp"));
        assert!(!reader.is_present("/labels/training"));
    }

    #[test]
    fn other_bytes_are_not_a_file() {
        assert!(Hdf5Reader::open(b"This is not an HDF5 file.\n".to_vec()).is_none());
    }
}
