//! A stable sort whose memory is a budget, not the input: sorted runs spill to a directory once
//! the budget is spent, and are merged back in one pass.
//!
//! Stable as Java's `TimSort` is (Spark's `sortByKey`, htsjdk's `SortingCollection`): items that
//! compare equal come out in the order they went in. Each run is a stable in-memory sort of a
//! contiguous stretch of the input, and the merge breaks ties by run index and then by position
//! in the run, so the merged order is exactly that of one stable sort of the whole input.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::io::{BufRead, BufReader, BufWriter, Seek, Write};
use std::path::PathBuf;

/// How an item goes to a run on disk and comes back.
pub trait Spill: Sized {
    /// About what the item holds in memory, which is what the budget counts.
    fn footprint(&self) -> usize;
    /// Appends the item's bytes; an `InvalidData` error carries the item's own refusal.
    fn spill(&self, out: &mut Vec<u8>) -> std::io::Result<()>;
    /// The next item, or `None` at the end of the run.
    fn unspill(input: &mut impl BufRead) -> std::io::Result<Option<Self>>;
}

/// Items pushed in input order, sorted with `compare` when [`Self::finish`] is called.
pub struct ExternalSorter<T, F> {
    compare: F,
    budget: usize,
    dir: PathBuf,
    buffer: Vec<T>,
    used: usize,
    runs: Vec<std::fs::File>,
    len: u64,
}

impl<T: Spill, F: Fn(&T, &T) -> Ordering + Copy> ExternalSorter<T, F> {
    /// Runs go to `dir` as unnamed files, gone when the sorter or its output is dropped.
    pub fn new(compare: F, budget: usize, dir: impl Into<PathBuf>) -> Self {
        ExternalSorter {
            compare,
            budget,
            dir: dir.into(),
            buffer: Vec::new(),
            used: 0,
            runs: Vec::new(),
            len: 0,
        }
    }

    pub fn push(&mut self, item: T) -> std::io::Result<()> {
        self.used += item.footprint() + std::mem::size_of::<T>();
        self.buffer.push(item);
        self.len += 1;
        if self.used > self.budget {
            self.spill()?;
        }
        Ok(())
    }

    /// How many items were pushed.
    pub fn len(&self) -> u64 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn spill(&mut self) -> std::io::Result<()> {
        self.buffer.sort_by(self.compare);
        let run = write_run(&self.dir, self.buffer.drain(..).map(Ok))?;
        self.runs.push(run);
        self.used = 0;
        if self.runs.len() >= MAX_RUNS {
            // shortcut: every run so far merges into one, rewriting what earlier merges wrote;
            // a merge by levels would bound that if inputs of thousands of runs become common.
            let runs = std::mem::take(&mut self.runs);
            let merged = merge(runs, Vec::new(), self.compare)?;
            let run = write_run(&self.dir, merged)?;
            self.runs.push(run);
        }
        Ok(())
    }

    /// The items in order. Nothing spilled: the buffer, sorted in place. Otherwise a k-way merge of
    /// the runs and of what is left in memory, which is the last run.
    pub fn finish(mut self) -> std::io::Result<Sorted<T, F>> {
        self.buffer.sort_by(self.compare);
        merge(self.runs, self.buffer, self.compare)
    }
}

/// A sorted run to disk, rewound for the merge.
fn write_run<T: Spill>(
    dir: &std::path::Path,
    items: impl Iterator<Item = std::io::Result<T>>,
) -> std::io::Result<std::fs::File> {
    let mut file = BufWriter::with_capacity(1 << 20, tempfile::tempfile_in(dir)?);
    let mut bytes = Vec::new();
    for item in items {
        bytes.clear();
        item?.spill(&mut bytes)?;
        file.write_all(&bytes)?;
    }
    let mut file = file.into_inner().map_err(|error| error.into_error())?;
    file.rewind()?;
    Ok(file)
}

/// Runs open at once before they are merged into one, which keeps a sort of any size within the
/// descriptors a process is given.
const MAX_RUNS: usize = 128;

/// The runs, in input order, and a sorted last run still in memory, merged.
fn merge<T: Spill, F: Fn(&T, &T) -> Ordering + Copy>(
    runs: Vec<std::fs::File>,
    memory: Vec<T>,
    compare: F,
) -> std::io::Result<Sorted<T, F>> {
    let mut sources: Vec<Source<T>> = Vec::with_capacity(runs.len() + 1);
    for run in runs {
        sources.push(Source::Run(BufReader::with_capacity(1 << 16, run)));
    }
    sources.push(Source::Memory(memory.into_iter()));
    let mut heap = BinaryHeap::with_capacity(sources.len());
    for (run, source) in sources.iter_mut().enumerate() {
        if let Some(item) = source.next()? {
            heap.push(Head { item, run, compare });
        }
    }
    Ok(Sorted { sources, heap })
}

enum Source<T> {
    Run(BufReader<std::fs::File>),
    Memory(std::vec::IntoIter<T>),
}

impl<T: Spill> Source<T> {
    fn next(&mut self) -> std::io::Result<Option<T>> {
        match self {
            Source::Run(reader) => T::unspill(reader),
            Source::Memory(items) => Ok(items.next()),
        }
    }
}

/// A run's next item. `BinaryHeap` pops the greatest, so the order is reversed: the least item,
/// and among equal ones the earliest run.
struct Head<T, F> {
    item: T,
    run: usize,
    compare: F,
}

impl<T, F: Fn(&T, &T) -> Ordering> Ord for Head<T, F> {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.compare)(&other.item, &self.item).then(other.run.cmp(&self.run))
    }
}

impl<T, F: Fn(&T, &T) -> Ordering> PartialOrd for Head<T, F> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T, F: Fn(&T, &T) -> Ordering> PartialEq for Head<T, F> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl<T, F: Fn(&T, &T) -> Ordering> Eq for Head<T, F> {}

/// The merged items, read from the runs as they are taken.
pub struct Sorted<T, F> {
    sources: Vec<Source<T>>,
    heap: BinaryHeap<Head<T, F>>,
}

impl<T: Spill, F: Fn(&T, &T) -> Ordering + Copy> Iterator for Sorted<T, F> {
    type Item = std::io::Result<T>;

    fn next(&mut self) -> Option<Self::Item> {
        let Head { item, run, compare } = self.heap.pop()?;
        match self.sources[run].next() {
            Ok(Some(next)) => self.heap.push(Head {
                item: next,
                run,
                compare,
            }),
            Ok(None) => {}
            Err(error) => {
                self.heap.clear();
                return Some(Err(error));
            }
        }
        Some(Ok(item))
    }
}

/// A read goes to a run as its BAM record, length prefix included.
impl Spill for htsjdk_bam::record::BamRecord {
    fn footprint(&self) -> usize {
        // The heap behind the fields, and an allocator's header for each of them.
        self.read_name.len()
            + self.read_bases.len()
            + self.base_qualities.len()
            + 8 * self.cigar.elements.len()
            + 64 * self.tags.len()
            + 64
    }

    fn spill(&self, out: &mut Vec<u8>) -> std::io::Result<()> {
        self.encode_into(out).map_err(|error| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, format!("{error:?}"))
        })
    }

    fn unspill(input: &mut impl BufRead) -> std::io::Result<Option<Self>> {
        let mut length = [0u8; 4];
        if input.fill_buf()?.is_empty() {
            return Ok(None);
        }
        input.read_exact(&mut length)?;
        let mut bytes = length.to_vec();
        bytes.resize(4 + u32::from_le_bytes(length) as usize, 0);
        input.read_exact(&mut bytes[4..])?;
        match htsjdk_bam::record::BamRecord::decode(&bytes) {
            Ok(Some((record, _))) => Ok(Some(record)),
            Ok(None) => Err(std::io::ErrorKind::UnexpectedEof.into()),
            Err(error) => Err(std::io::Error::other(format!("{error:?}"))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    /// A key with few values and the item's input position, which the comparison ignores.
    #[derive(Debug, Clone, PartialEq)]
    struct Item {
        key: u8,
        position: u32,
    }

    impl Spill for Item {
        fn footprint(&self) -> usize {
            0
        }

        fn spill(&self, out: &mut Vec<u8>) -> std::io::Result<()> {
            out.push(self.key);
            out.extend_from_slice(&self.position.to_le_bytes());
            Ok(())
        }

        fn unspill(input: &mut impl BufRead) -> std::io::Result<Option<Self>> {
            let mut bytes = [0u8; 5];
            if input.fill_buf()?.is_empty() {
                return Ok(None);
            }
            input.read_exact(&mut bytes)?;
            Ok(Some(Item {
                key: bytes[0],
                position: u32::from_le_bytes(bytes[1..].try_into().unwrap()),
            }))
        }
    }

    fn by_key(a: &Item, b: &Item) -> Ordering {
        a.key.cmp(&b.key)
    }

    fn external(items: &[Item], budget: usize) -> Vec<Item> {
        let dir = tempfile::tempdir().unwrap();
        let mut sorter = ExternalSorter::new(by_key, budget, dir.path());
        for item in items {
            sorter.push(item.clone()).unwrap();
        }
        assert_eq!(sorter.len(), items.len() as u64);
        sorter.finish().unwrap().map(Result::unwrap).collect()
    }

    #[test]
    fn equals_a_stable_in_memory_sort_on_many_equal_keys() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let item = std::mem::size_of::<Item>();
        for (count, keys) in [
            (0, 1),
            (1, 1),
            (7, 1),
            (1000, 1),
            (1000, 3),
            (4099, 17),
            (5000, 256),
        ] {
            let items: Vec<Item> = (0..count)
                .map(|position| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    Item {
                        key: (state % keys) as u8,
                        position,
                    }
                })
                .collect();
            let mut stable = items.clone();
            stable.sort_by(by_key);
            // A run of one item, of a few, of many, and none spilled at all.
            for budget in [0, item, 3 * item, 64 * item, usize::MAX] {
                assert_eq!(
                    external(&items, budget),
                    stable,
                    "{count} items, {keys} keys, {budget}"
                );
            }
        }
    }

    #[test]
    fn descending_and_constant_inputs_keep_their_order_among_equals() {
        let descending: Vec<Item> = (0..300)
            .map(|position| Item {
                key: (299 - position) as u8 / 10,
                position,
            })
            .collect();
        let mut stable = descending.clone();
        stable.sort_by(by_key);
        assert_eq!(
            external(&descending, 7 * std::mem::size_of::<Item>()),
            stable
        );
        let constant: Vec<Item> = (0..300).map(|position| Item { key: 4, position }).collect();
        assert_eq!(
            external(&constant, 7 * std::mem::size_of::<Item>()),
            constant
        );
    }

    #[test]
    fn a_read_comes_back_from_a_run_as_it_went() {
        let mut record = htsjdk_bam::record::BamRecord {
            read_name: "q1".to_string(),
            read_bases: b"ACGT".to_vec(),
            base_qualities: vec![30, 31, 32, 33],
            ..Default::default()
        };
        record.tags.insert(
            htsjdk_bam::tag::Tag::new(b"RG"),
            htsjdk_bam::tag::TagValue::Str("rg".to_string()),
        );
        let mut bytes = Vec::new();
        record.spill(&mut bytes).unwrap();
        record.spill(&mut bytes).unwrap();
        let mut input = std::io::BufReader::new(bytes.as_slice());
        for _ in 0..2 {
            let back = htsjdk_bam::record::BamRecord::unspill(&mut input).unwrap();
            assert_eq!(back.as_ref(), Some(&record));
        }
        assert!(htsjdk_bam::record::BamRecord::unspill(&mut input)
            .unwrap()
            .is_none());
        assert_eq!(input.bytes().count(), 0);
    }
}
