//! A minimal HDF5 reader: the groups and the double and string datasets `HDF5File` reads back.
//!
//! Two kinds of file reach it. The reference's are written by the HDF5 C library in its default
//! layout: a version 0 superblock, version 1 object headers, groups held in a B-tree of symbol
//! table nodes over a local heap. The port's are [`crate::hdf5_writer`]'s, in the "latest" layout:
//! a version 2 superblock, version 2 object headers and links held in the header itself.
//!
//! Written from "HDF5 File Format Specification Version 3.0" (The HDF Group); no HDF5 source was
//! consulted. Anything outside what those two writers make (a filter on a chunk, big-endian
//! values, a dense group) reads as absent rather than as a guess.

/// `UNDEF_ADDRESS`, all ones.
const UNDEFINED: u64 = u64::MAX;

const SIGNATURE: &[u8] = b"\x89HDF\r\n\x1a\n";

const DATASPACE: u16 = 0x01;
const DATATYPE: u16 = 0x03;
const LINK: u16 = 0x06;
const LAYOUT: u16 = 0x08;
const CONTINUATION: u16 = 0x10;
const SYMBOL_TABLE: u16 = 0x11;

/// A file's bytes, opened.
#[derive(Debug, Clone)]
pub struct Hdf5Reader {
    bytes: Vec<u8>,
    root: u64,
}

/// One object header's messages, by type.
type Messages = Vec<(u16, Vec<u8>)>;

impl Hdf5Reader {
    /// Whether the bytes begin as an HDF5 file does, and the root group's address.
    pub fn open(bytes: Vec<u8>) -> Option<Self> {
        if bytes.get(..8)? != SIGNATURE {
            return None;
        }
        let root = match *bytes.get(8)? {
            // The root group's symbol table entry is at 56, its header's address eight bytes in.
            0 => u64_at(&bytes, 64)?,
            2 | 3 => u64_at(&bytes, 12 + 24)?,
            _ => return None,
        };
        Some(Self { bytes, root })
    }

    /// The doubles at `path` in row-major order, with their shape.
    pub fn doubles(&self, path: &str) -> Option<(Vec<u64>, Vec<f64>)> {
        let messages = self.object(path)?;
        let datatype = message(&messages, DATATYPE)?;
        // Class 1, floating point, eight bytes, and little-endian: bit 0 of the first class field.
        if datatype.first()? & 0x0f != 1 || u32_at(datatype, 4)? != 8 || datatype.get(1)? & 1 != 0 {
            return None;
        }
        let shape = shape_of(message(&messages, DATASPACE)?)?;
        let count: u64 = shape.iter().product();
        let raw = self.raw(&messages, count.checked_mul(8)?)?;
        let values = raw
            .as_chunks::<8>()
            .0
            .iter()
            .map(|chunk| f64::from_le_bytes(*chunk))
            .collect();
        Some((shape, values))
    }

    /// The strings at `path`, fixed-length or variable-length, whichever the library chose.
    pub fn strings(&self, path: &str) -> Option<Vec<String>> {
        let messages = self.object(path)?;
        let datatype = message(&messages, DATATYPE)?;
        let shape = shape_of(message(&messages, DATASPACE)?)?;
        let count: u64 = shape.iter().product();
        let size = u64::from(u32_at(datatype, 4)?);
        let raw = self.raw(&messages, count.checked_mul(size)?)?;
        match datatype.first()? & 0x0f {
            3 => Some(
                raw.chunks_exact(usize::try_from(size).ok()?.max(1))
                    .map(|chunk| {
                        let end = chunk.iter().position(|b| *b == 0).unwrap_or(chunk.len());
                        String::from_utf8_lossy(&chunk[..end]).into_owned()
                    })
                    .collect(),
            ),
            9 => raw
                .chunks_exact(usize::try_from(size).ok()?.max(1))
                .map(|chunk| self.global_heap_string(chunk))
                .collect(),
            _ => None,
        }
    }

    /// One variable-length element: its length, its heap's address and its index in that heap.
    fn global_heap_string(&self, element: &[u8]) -> Option<String> {
        let length = usize::try_from(u32_at(element, 0)?).ok()?;
        let heap = usize::try_from(u64_at(element, 4)?).ok()?;
        let index = u32_at(element, 12)?;
        if self.bytes.get(heap..heap + 4)? != b"GCOL" {
            return None;
        }
        let end = heap + usize::try_from(u64_at(&self.bytes, heap + 8)?).ok()?;
        let mut at = heap + 16;
        while at + 16 <= end {
            let object = u16_at(&self.bytes, at)?;
            let size = usize::try_from(u64_at(&self.bytes, at + 8)?).ok()?;
            if u32::from(object) == index {
                let data = self.bytes.get(at + 16..at + 16 + length)?;
                return Some(String::from_utf8_lossy(data).into_owned());
            }
            if object == 0 {
                return None;
            }
            at += 16 + size.div_ceil(8) * 8;
        }
        None
    }

    /// The object header at the end of `path`, walked from the root group.
    fn object(&self, path: &str) -> Option<Messages> {
        let mut address = self.root;
        for name in path.split('/').filter(|part| !part.is_empty()) {
            address = self.child(address, name)?;
        }
        self.messages(address)
    }

    /// A group's child by name, through its links or its symbol table.
    fn child(&self, group: u64, name: &str) -> Option<u64> {
        let messages = self.messages(group)?;
        for (kind, body) in &messages {
            if *kind == LINK {
                if let Some((child, address)) = link(body) {
                    if child == name {
                        return Some(address);
                    }
                }
            }
        }
        let table = message(&messages, SYMBOL_TABLE)?;
        let tree = u64_at(table, 0)?;
        let heap = u64_at(table, 8)?;
        let heap = usize::try_from(heap).ok()?;
        if self.bytes.get(heap..heap + 4)? != b"HEAP" {
            return None;
        }
        let names = usize::try_from(u64_at(&self.bytes, heap + 24)?).ok()?;
        self.search_group_tree(tree, names, name)
    }

    /// A version 1 B-tree of symbol table nodes, every node visited.
    fn search_group_tree(&self, tree: u64, names: usize, name: &str) -> Option<u64> {
        let at = usize::try_from(tree).ok()?;
        if self.bytes.get(at..at + 4)? != b"TREE" {
            return None;
        }
        let level = *self.bytes.get(at + 5)?;
        let used = usize::from(u16_at(&self.bytes, at + 6)?);
        for entry in 0..used {
            // After the 24-byte node header, keys and children alternate, eight bytes each.
            let child = u64_at(&self.bytes, at + 24 + 8 + entry * 16)?;
            let found = if level > 0 {
                self.search_group_tree(child, names, name)
            } else {
                self.search_symbol_node(child, names, name)
            };
            if found.is_some() {
                return found;
            }
        }
        None
    }

    fn search_symbol_node(&self, node: u64, names: usize, name: &str) -> Option<u64> {
        let at = usize::try_from(node).ok()?;
        if self.bytes.get(at..at + 4)? != b"SNOD" {
            return None;
        }
        let count = usize::from(u16_at(&self.bytes, at + 6)?);
        for entry in 0..count {
            let base = at + 8 + entry * 40;
            let offset = usize::try_from(u64_at(&self.bytes, base)?).ok()?;
            let text = self.bytes.get(names + offset..)?;
            let end = text.iter().position(|b| *b == 0)?;
            if &text[..end] == name.as_bytes() {
                return u64_at(&self.bytes, base + 8);
            }
        }
        None
    }

    /// An object header's messages, followed through its continuation blocks.
    fn messages(&self, address: u64) -> Option<Messages> {
        let at = usize::try_from(address).ok()?;
        let mut out = Vec::new();
        if self.bytes.get(at..at + 4)? == b"OHDR" {
            let flags = *self.bytes.get(at + 5)?;
            let mut cursor = at + 6;
            if flags & 0x20 != 0 {
                cursor += 16;
            }
            if flags & 0x10 != 0 {
                cursor += 4;
            }
            let width = 1usize << (flags & 3);
            let mut size = [0u8; 8];
            size[..width].copy_from_slice(self.bytes.get(cursor..cursor + width)?);
            let size = usize::try_from(u64::from_le_bytes(size)).ok()?;
            cursor += width;
            self.messages_v2(cursor, cursor + size, flags, &mut out)?;
        } else {
            if *self.bytes.get(at)? != 1 {
                return None;
            }
            let size = usize::try_from(u32_at(&self.bytes, at + 8)?).ok()?;
            self.messages_v1(at + 16, at + 16 + size, &mut out)?;
        }
        Some(out)
    }

    fn messages_v1(&self, mut at: usize, end: usize, out: &mut Messages) -> Option<()> {
        while at + 8 <= end {
            let kind = u16_at(&self.bytes, at)?;
            let size = usize::from(u16_at(&self.bytes, at + 2)?);
            let body = self.bytes.get(at + 8..at + 8 + size)?;
            if kind == CONTINUATION {
                let block = usize::try_from(u64_at(body, 0)?).ok()?;
                let length = usize::try_from(u64_at(body, 8)?).ok()?;
                self.messages_v1(block, block + length, out)?;
            } else {
                out.push((kind, body.to_vec()));
            }
            at += 8 + size;
        }
        Some(())
    }

    fn messages_v2(&self, mut at: usize, end: usize, flags: u8, out: &mut Messages) -> Option<()> {
        let prefix = if flags & 0x04 != 0 { 6 } else { 4 };
        while at + prefix <= end {
            let kind = u16::from(*self.bytes.get(at)?);
            let size = usize::from(u16_at(&self.bytes, at + 1)?);
            let body = self.bytes.get(at + prefix..at + prefix + size)?;
            if kind == CONTINUATION {
                let block = usize::try_from(u64_at(body, 0)?).ok()?;
                let length = usize::try_from(u64_at(body, 8)?).ok()?;
                if self.bytes.get(block..block + 4)? != b"OCHK" {
                    return None;
                }
                self.messages_v2(block + 4, block + length - 4, flags, out)?;
            } else if kind != 0 {
                out.push((kind, body.to_vec()));
            }
            at += prefix + size;
        }
        Some(())
    }

    /// A dataset's raw bytes, `expected` of them.
    fn raw(&self, messages: &Messages, expected: u64) -> Option<Vec<u8>> {
        let layout = message(messages, LAYOUT)?;
        let expected = usize::try_from(expected).ok()?;
        if expected == 0 {
            return Some(Vec::new());
        }
        if *layout.first()? != 3 {
            return None;
        }
        match *layout.get(1)? {
            0 => {
                let size = usize::from(u16_at(layout, 2)?);
                Some(layout.get(4..4 + size)?.to_vec())
            }
            1 => {
                let address = u64_at(layout, 2)?;
                if address == UNDEFINED {
                    return None;
                }
                let at = usize::try_from(address).ok()?;
                Some(self.bytes.get(at..at + expected)?.to_vec())
            }
            2 => {
                let rank = usize::from(*layout.get(2)?);
                let tree = u64_at(layout, 3)?;
                let chunk: Vec<usize> = (0..rank)
                    .map(|i| u32_at(layout, 11 + i * 4).and_then(|v| usize::try_from(v).ok()))
                    .collect::<Option<_>>()?;
                let shape = shape_of(message(messages, DATASPACE)?)?;
                let mut out = vec![0u8; expected];
                self.copy_chunks(tree, &chunk, &shape, &mut out)?;
                Some(out)
            }
            _ => None,
        }
    }

    /// Every chunk of a version 1 chunk B-tree into its place, unfiltered chunks only.
    fn copy_chunks(&self, tree: u64, chunk: &[usize], shape: &[u64], out: &mut [u8]) -> Option<()> {
        let at = usize::try_from(tree).ok()?;
        if self.bytes.get(at..at + 4)? != b"TREE" {
            return None;
        }
        let level = *self.bytes.get(at + 5)?;
        let used = usize::from(u16_at(&self.bytes, at + 6)?);
        let rank = chunk.len();
        let key = 8 + rank * 8;
        for entry in 0..used {
            let base = at + 24 + entry * (key + 8);
            let child = u64_at(&self.bytes, base + key)?;
            if level > 0 {
                self.copy_chunks(child, chunk, shape, out)?;
                continue;
            }
            if u32_at(&self.bytes, base + 4)? != 0 {
                return None;
            }
            let size = usize::try_from(u32_at(&self.bytes, base)?).ok()?;
            let data = self
                .bytes
                .get(usize::try_from(child).ok()?..)?
                .get(..size)?;
            let offsets: Vec<usize> = (0..rank - 1)
                .map(|i| {
                    u64_at(&self.bytes, base + 8 + i * 8).and_then(|v| usize::try_from(v).ok())
                })
                .collect::<Option<_>>()?;
            copy_chunk(data, &offsets, chunk, shape, out)?;
        }
        Some(())
    }
}

/// One chunk's elements into the row-major whole. A chunk may overhang the dataset's edge.
fn copy_chunk(
    data: &[u8],
    offsets: &[usize],
    chunk: &[usize],
    shape: &[u64],
    out: &mut [u8],
) -> Option<()> {
    let rank = offsets.len();
    let element = *chunk.get(rank)?;
    let chunk_elements: usize = chunk[..rank].iter().product();
    for index in 0..chunk_elements {
        let mut remaining = index;
        let mut position = vec![0usize; rank];
        for dimension in (0..rank).rev() {
            position[dimension] = remaining % chunk[dimension];
            remaining /= chunk[dimension];
        }
        let mut target = 0usize;
        let mut inside = true;
        for dimension in 0..rank {
            let value = offsets[dimension] + position[dimension];
            let extent = usize::try_from(shape[dimension]).ok()?;
            if value >= extent {
                inside = false;
                break;
            }
            target = target * extent + value;
        }
        if inside {
            let from = data.get(index * element..(index + 1) * element)?;
            out.get_mut(target * element..(target + 1) * element)?
                .copy_from_slice(from);
        }
    }
    Some(())
}

fn message(messages: &Messages, kind: u16) -> Option<&[u8]> {
    messages
        .iter()
        .find(|(found, _)| *found == kind)
        .map(|(_, body)| body.as_slice())
}

/// A hard link's target by name, in a version 1 Link message.
fn link(body: &[u8]) -> Option<(String, u64)> {
    if *body.first()? != 1 {
        return None;
    }
    let flags = *body.get(1)?;
    let mut at = 2;
    let kind = if flags & 0x08 != 0 {
        at += 1;
        *body.get(at - 1)?
    } else {
        0
    };
    if flags & 0x04 != 0 {
        at += 8;
    }
    if flags & 0x10 != 0 {
        at += 1;
    }
    let width = 1usize << (flags & 3);
    let mut length = [0u8; 8];
    length[..width].copy_from_slice(body.get(at..at + width)?);
    let length = usize::try_from(u64::from_le_bytes(length)).ok()?;
    at += width;
    let name = String::from_utf8_lossy(body.get(at..at + length)?).into_owned();
    at += length;
    (kind == 0).then_some(())?;
    Some((name, u64_at(body, at)?))
}

/// A Dataspace message's dimensions, in either version.
fn shape_of(body: &[u8]) -> Option<Vec<u64>> {
    let rank = usize::from(*body.get(1)?);
    let start = match *body.first()? {
        1 => 8,
        2 => 4,
        _ => return None,
    };
    (0..rank).map(|i| u64_at(body, start + i * 8)).collect()
}

fn u16_at(bytes: &[u8], at: usize) -> Option<u16> {
    Some(u16::from_le_bytes(bytes.get(at..at + 2)?.try_into().ok()?))
}

fn u32_at(bytes: &[u8], at: usize) -> Option<u32> {
    Some(u32::from_le_bytes(bytes.get(at..at + 4)?.try_into().ok()?))
}

fn u64_at(bytes: &[u8], at: usize) -> Option<u64> {
    Some(u64::from_le_bytes(bytes.get(at..at + 8)?.try_into().ok()?))
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
    }

    #[test]
    fn other_bytes_are_not_a_file() {
        assert!(Hdf5Reader::open(b"This is not an HDF5 file.\n".to_vec()).is_none());
    }
}
