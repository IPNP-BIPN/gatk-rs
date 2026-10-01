//! A minimal HDF5 writer: groups, and contiguous datasets of doubles or fixed-length strings.
//!
//! GATK writes its HDF5 files through the HDF5 C library, and two runs of the reference on one
//! input already write different bytes (an object header carries the time it was made), so no
//! claim about the BYTES of such a file is available to anyone. What a tool decides is the tree:
//! which groups and datasets exist, each dataset's type, shape and values. That is what the
//! covering-array harness reads back from both sides, and that is all this writer has to get
//! right, in a file any HDF5 reader opens.
//!
//! The file is the "latest format" layout the library itself writes with `H5F_LIBVER_LATEST`:
//!
//! * a version 2 superblock, 48 bytes, pointing at the root group's object header;
//! * version 2 object headers (`OHDR`), each checksummed with Jenkins' `lookup3` as the format
//!   specification requires;
//! * groups in compact storage: a Link Info message, a Group Info message and one hard Link
//!   message per child;
//! * datasets with a version 2 dataspace, an IEEE little-endian double or a null-padded
//!   fixed-length ASCII string type, a version 3 fill value message, and a contiguous layout
//!   pointing at the raw data.
//!
//! Written from "HDF5 File Format Specification Version 3.0" (The HDF Group), sections II.C
//! (superblock), IV.A.1.b (version 2 object headers) and IV.A.2 (messages); no HDF5 source was
//! consulted.

use std::collections::BTreeMap;

/// `UNDEF_ADDRESS`, all ones.
const UNDEFINED: u64 = u64::MAX;

/// What a dataset holds.
#[derive(Debug, Clone, PartialEq)]
pub enum Data {
    /// Doubles in row-major order, with their shape.
    Doubles { shape: Vec<u64>, values: Vec<f64> },
    /// A one-dimensional array of strings, stored at the longest one's length.
    Strings(Vec<String>),
}

/// A file being built: datasets by absolute path. Every group a path names is created.
#[derive(Debug, Clone, Default)]
pub struct Hdf5File {
    datasets: BTreeMap<String, Data>,
}

/// A group's children, by name.
#[derive(Default)]
struct Group {
    groups: BTreeMap<String, Group>,
    datasets: BTreeMap<String, Data>,
}

impl Hdf5File {
    pub fn new() -> Self {
        Self::default()
    }

    /// A one-dimensional array of doubles.
    pub fn double_array(&mut self, path: &str, values: &[f64]) {
        self.datasets.insert(
            path.to_string(),
            Data::Doubles {
                shape: vec![values.len() as u64],
                values: values.to_vec(),
            },
        );
    }

    /// A matrix of doubles, one row after another.
    pub fn double_matrix(&mut self, path: &str, rows: &[Vec<f64>]) {
        let columns = rows.first().map_or(0, Vec::len) as u64;
        self.datasets.insert(
            path.to_string(),
            Data::Doubles {
                shape: vec![rows.len() as u64, columns],
                values: rows.iter().flatten().copied().collect(),
            },
        );
    }

    pub fn string_array(&mut self, path: &str, values: &[String]) {
        self.datasets
            .insert(path.to_string(), Data::Strings(values.to_vec()));
    }

    /// The file's bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut root = Group::default();
        for (path, data) in &self.datasets {
            let parts: Vec<&str> = path.split('/').filter(|p| !p.is_empty()).collect();
            let Some((name, parents)) = parts.split_last() else {
                continue;
            };
            let mut group = &mut root;
            for parent in parents {
                group = group.groups.entry((*parent).to_string()).or_default();
            }
            group.datasets.insert((*name).to_string(), data.clone());
        }
        let mut out = vec![0u8; 48];
        let root_address = write_group(&root, &mut out);
        let end = out.len() as u64;
        let mut superblock = Vec::with_capacity(48);
        superblock.extend_from_slice(b"\x89HDF\r\n\x1a\n");
        superblock.extend_from_slice(&[2, 8, 8, 0]);
        superblock.extend_from_slice(&0u64.to_le_bytes());
        superblock.extend_from_slice(&UNDEFINED.to_le_bytes());
        superblock.extend_from_slice(&end.to_le_bytes());
        superblock.extend_from_slice(&root_address.to_le_bytes());
        let checksum = lookup3(&superblock);
        superblock.extend_from_slice(&checksum.to_le_bytes());
        out[..48].copy_from_slice(&superblock);
        out
    }
}

/// One message: its type and its body.
struct Message {
    kind: u8,
    body: Vec<u8>,
}

/// A version 2 object header holding `messages`, appended to `out`; its address is returned.
fn write_object_header(messages: &[Message], out: &mut Vec<u8>) -> u64 {
    let address = out.len() as u64;
    let size: usize = messages.iter().map(|m| 4 + m.body.len()).sum();
    let mut header = Vec::with_capacity(12 + size);
    header.extend_from_slice(b"OHDR");
    // Version 2; flags 2: the size of chunk 0 is a four-byte field, and nothing optional follows.
    header.extend_from_slice(&[2, 2]);
    header.extend_from_slice(&(size as u32).to_le_bytes());
    for message in messages {
        header.push(message.kind);
        header.extend_from_slice(&(message.body.len() as u16).to_le_bytes());
        header.push(0);
        header.extend_from_slice(&message.body);
    }
    let checksum = lookup3(&header);
    header.extend_from_slice(&checksum.to_le_bytes());
    out.extend_from_slice(&header);
    address
}

/// A group, its children first so that their addresses are known when its links are written.
fn write_group(group: &Group, out: &mut Vec<u8>) -> u64 {
    let mut children: BTreeMap<&str, u64> = BTreeMap::new();
    for (name, child) in &group.groups {
        children.insert(name, write_group(child, out));
    }
    for (name, data) in &group.datasets {
        children.insert(name, write_dataset(data, out));
    }
    let mut messages = Vec::new();
    // Link Info: version 0, no flags, and neither a fractal heap nor a name index.
    let mut link_info = vec![0u8, 0];
    link_info.extend_from_slice(&UNDEFINED.to_le_bytes());
    link_info.extend_from_slice(&UNDEFINED.to_le_bytes());
    messages.push(Message {
        kind: 0x02,
        body: link_info,
    });
    // Group Info: version 0, no flags.
    messages.push(Message {
        kind: 0x0A,
        body: vec![0, 0],
    });
    for (name, address) in children {
        // Link: version 1; flags 0 means a one-byte name length, a hard link and ASCII.
        let mut link = vec![1u8, 0, name.len() as u8];
        link.extend_from_slice(name.as_bytes());
        link.extend_from_slice(&address.to_le_bytes());
        messages.push(Message {
            kind: 0x06,
            body: link,
        });
    }
    write_object_header(&messages, out)
}

/// A dataset: the raw data first, then the header that points at it.
fn write_dataset(data: &Data, out: &mut Vec<u8>) -> u64 {
    let (shape, datatype, raw): (Vec<u64>, Vec<u8>, Vec<u8>) = match data {
        Data::Doubles { shape, values } => {
            let mut datatype = vec![0x11, 0x20, 0x3f, 0x00];
            datatype.extend_from_slice(&8u32.to_le_bytes());
            datatype.extend_from_slice(&0u16.to_le_bytes());
            datatype.extend_from_slice(&64u16.to_le_bytes());
            datatype.extend_from_slice(&[52, 11, 0, 52]);
            datatype.extend_from_slice(&1023u32.to_le_bytes());
            let raw = values.iter().flat_map(|v| v.to_le_bytes()).collect();
            (shape.clone(), datatype, raw)
        }
        Data::Strings(values) => {
            let width = values.iter().map(String::len).max().unwrap_or(0).max(1);
            // Class 3, version 1; null padded, ASCII.
            let mut datatype = vec![0x13, 0x01, 0x00, 0x00];
            datatype.extend_from_slice(&(width as u32).to_le_bytes());
            let mut raw = Vec::with_capacity(width * values.len());
            for value in values {
                raw.extend_from_slice(value.as_bytes());
                raw.resize(raw.len() + width - value.len(), 0);
            }
            (vec![values.len() as u64], datatype, raw)
        }
    };
    let address = if raw.is_empty() {
        UNDEFINED
    } else {
        let at = out.len() as u64;
        out.extend_from_slice(&raw);
        at
    };
    let mut dataspace = vec![2u8, shape.len() as u8, 0, 1];
    for dimension in &shape {
        dataspace.extend_from_slice(&dimension.to_le_bytes());
    }
    let mut layout = vec![3u8, 1];
    layout.extend_from_slice(&address.to_le_bytes());
    layout.extend_from_slice(&(raw.len() as u64).to_le_bytes());
    let messages = [
        Message {
            kind: 0x01,
            body: dataspace,
        },
        Message {
            kind: 0x03,
            body: datatype,
        },
        // Fill value, version 3: allocated early, never written, none defined.
        Message {
            kind: 0x05,
            body: vec![3, 0x05],
        },
        Message {
            kind: 0x08,
            body: layout,
        },
    ];
    write_object_header(&messages, out)
}

/// Bob Jenkins' `hashlittle` from `lookup3.c` with an initial value of zero, which is the
/// checksum the HDF5 format puts on every version 2 structure.
pub fn lookup3(key: &[u8]) -> u32 {
    fn rot(x: u32, k: u32) -> u32 {
        x.rotate_left(k)
    }
    let length = key.len();
    let mut a: u32 = 0xdeadbeef_u32.wrapping_add(length as u32);
    let mut b = a;
    let mut c = a;
    let word = |chunk: &[u8]| -> u32 {
        let mut w = 0u32;
        for (i, byte) in chunk.iter().enumerate() {
            w |= (*byte as u32) << (8 * i);
        }
        w
    };
    let mut rest = key;
    while rest.len() > 12 {
        a = a.wrapping_add(word(&rest[0..4]));
        b = b.wrapping_add(word(&rest[4..8]));
        c = c.wrapping_add(word(&rest[8..12]));
        // mix(a, b, c)
        a = a.wrapping_sub(c);
        a ^= rot(c, 4);
        c = c.wrapping_add(b);
        b = b.wrapping_sub(a);
        b ^= rot(a, 6);
        a = a.wrapping_add(c);
        c = c.wrapping_sub(b);
        c ^= rot(b, 8);
        b = b.wrapping_add(a);
        a = a.wrapping_sub(c);
        a ^= rot(c, 16);
        c = c.wrapping_add(b);
        b = b.wrapping_sub(a);
        b ^= rot(a, 19);
        a = a.wrapping_add(c);
        c = c.wrapping_sub(b);
        c ^= rot(b, 4);
        b = b.wrapping_add(a);
        rest = &rest[12..];
    }
    if rest.is_empty() {
        return c;
    }
    let mut tail = [0u8; 12];
    tail[..rest.len()].copy_from_slice(rest);
    a = a.wrapping_add(word(&tail[0..4]));
    b = b.wrapping_add(word(&tail[4..8]));
    c = c.wrapping_add(word(&tail[8..12]));
    // final(a, b, c)
    c ^= b;
    c = c.wrapping_sub(rot(b, 14));
    a ^= c;
    a = a.wrapping_sub(rot(c, 11));
    b ^= a;
    b = b.wrapping_sub(rot(a, 25));
    c ^= b;
    c = c.wrapping_sub(rot(b, 16));
    a ^= c;
    a = a.wrapping_sub(rot(c, 4));
    b ^= a;
    b = b.wrapping_sub(rot(a, 14));
    c ^= b;
    c = c.wrapping_sub(rot(b, 24));
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The test vectors `lookup3.c` prints from its own driver.
    #[test]
    fn lookup3_matches_the_published_vectors() {
        assert_eq!(lookup3(b""), 0xdeadbeef);
        assert_eq!(lookup3(b"Four score and seven years ago"), 0x17770551);
    }

    #[test]
    fn a_file_starts_with_the_signature() {
        let mut file = Hdf5File::new();
        file.double_array("/labels/snp", &[1.0, 0.0]);
        let bytes = file.to_bytes();
        assert_eq!(&bytes[..8], b"\x89HDF\r\n\x1a\n");
        assert_eq!(bytes[8], 2);
    }
}
