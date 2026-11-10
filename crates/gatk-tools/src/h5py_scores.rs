//! The bytes of the scores file `isolation-forest.py` writes with h5py 3.10.0 and HDF5 1.14:
//! `f.create_dataset('data/scores', (n,), dtype='d')`.
//!
//! Unlike the file `HDF5File` writes through the Java binding, this one holds no time, and two
//! runs of the reference write the same bytes, measured in `docs/python-tools.md`. Its metadata is
//! 2432 bytes long whatever `n` is, and the only fields `n` moves are the end of the file, the
//! dataspace's two dimensions and the layout's size: the rest is the template below, as the
//! reference leaves it.

/// The runs of the template that are not zero, as `(offset, bytes)`.
const TEMPLATE: [(usize, &[u8]); 7] = [
    (0x0000, b"\x89HDF\r\n\x1a\n\0\0\0\0\0\x08\x08\0\x04\0\x10\0\0\0\0\0\0\0\0\0\0\0\0\0\xff\xff\xff\xff\xff\xff\xff\xff\0\0\0\0\0\0\0\0\xff\xff\xff\xff\xff\xff\xff\xff\0\0\0\0\0\0\0\0\x60\0\0\0\0\0\0\0\x01\0\0\0\0\0\0\0\x88\0\0\0\0\0\0\0\xa8\x02\0\0\0\0\0\0\x01\0\x01\0\x01\0\0\0\x18\0\0\0\0\0\0\0\x11\0\x10\0\0\0\0\0\x88\0\0\0\0\0\0\0\xa8\x02\0\0\0\0\0\0TREE\0\0\x01\0\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\0\0\0\0\0\0\0\0\xe0\x05\0\0\0\0\0\0\x08"),
    (0x02a8, b"HEAP\0\0\0\0\x58\0\0\0\0\0\0\0\x10\0\0\0\0\0\0\0\xc8\x02\0\0\0\0\0\0\0\0\0\0\0\0\0\0data\0\0\0\0\x01\0\0\0\0\0\0\0\x48"),
    (0x0320, b"\x01\0\x01\0\x01\0\0\0\x18\0\0\0\0\0\0\0\x11\0\x10\0\0\0\0\0\x48\x03\0\0\0\0\0\0\x68\x05\0\0\0\0\0\0TREE\0\0\x01\0\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\0\0\0\0\0\0\0\0\x38\x08\0\0\0\0\0\0\x08"),
    (0x0568, b"HEAP\0\0\0\0\x58\0\0\0\0\0\0\0\x10\0\0\0\0\0\0\0\x88\x05\0\0\0\0\0\0\0\0\0\0\0\0\0\0scores\0\0\x01\0\0\0\0\0\0\0\x48"),
    (0x05e0, b"SNOD\x01\0\x01\0\x08\0\0\0\0\0\0\0\x20\x03\0\0\0\0\0\0\x01\0\0\0\0\0\0\0\x48\x03\0\0\0\0\0\0\x68\x05"),
    (0x0728, b"\x01\0\x05\0\x01\0\0\0\0\x01\0\0\0\0\0\0\x01\0\x18\0\0\0\0\0\x01\x01\x01\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\x03\0\x18\0\x01\0\0\0\x11\x20\x3f\0\x08\0\0\0\0\0\x40\0\x34\x0b\0\x34\xff\x03\0\0\0\0\0\0\x05\0\x08\0\x01\0\0\0\x02\x02\x02\x01\0\0\0\0\x08\0\x18\0\0\0\0\0\x03\x01\x80\x09\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\x88\0\0"),
    (0x0838, b"SNOD\x01\0\x01\0\x08\0\0\0\0\0\0\0\x28\x07"),
];

const METADATA: usize = 0x980;
const END_OF_FILE: usize = 0x28;
const DIMENSION: usize = 0x748;
const MAXIMUM_DIMENSION: usize = 0x750;
const LAYOUT_SIZE: usize = 0x79a;

/// The file.
pub fn write(scores: &[f64]) -> Vec<u8> {
    let mut out = vec![0u8; METADATA];
    for (at, bytes) in TEMPLATE {
        out[at..at + bytes.len()].copy_from_slice(bytes);
    }
    let n = scores.len() as u64;
    let data = n * 8;
    let end = METADATA as u64 + data;
    out[END_OF_FILE..END_OF_FILE + 8].copy_from_slice(&end.to_le_bytes());
    out[DIMENSION..DIMENSION + 8].copy_from_slice(&n.to_le_bytes());
    out[MAXIMUM_DIMENSION..MAXIMUM_DIMENSION + 8].copy_from_slice(&n.to_le_bytes());
    out[LAYOUT_SIZE..LAYOUT_SIZE + 8].copy_from_slice(&data.to_le_bytes());
    for score in scores {
        out.extend_from_slice(&score.to_le_bytes());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hdf5_reader::Hdf5Reader;

    #[test]
    fn what_is_written_reads_back_through_the_reader_of_the_reference_files() {
        let written = write(&[-0.25, -0.5, -0.75]);
        assert_eq!(written.len(), 2432 + 24);
        let reader = Hdf5Reader::open(written).expect("an HDF5 file");
        let (shape, values) = reader.doubles("/data/scores").expect("the scores");
        assert_eq!(shape, vec![3]);
        assert_eq!(values, vec![-0.25, -0.5, -0.75]);
    }
}
