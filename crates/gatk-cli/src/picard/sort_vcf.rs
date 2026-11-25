//! `SortVcf.doWork`, over [`gatk_tools::sort_vcf`].
//!
//! The library sorts texts; what is here is what the tool does with paths first: every input and
//! the dictionary asserted readable, the dictionary read through `SamReaderFactory`, and handed to
//! each input that has none of its own. An input that declares one keeps it, and the given one is
//! only checked against it.

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::picard_gather_vcfs::{assert_same_dictionary, dictionary, Sequence};
use gatk_tools::sort_vcf::{contig_indexes, sort, Checks, SortError};
use std::io::{BufRead, Write};

use super::{
    absolute, assert_readable, thrown, vcf_lines, vcf_text, with_contigs, write_vcf, VcfStream,
};
use crate::runners::{argument, arguments, dictionary_from_any, flag};
use crate::vcf_stream::{Body, Records};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let inputs = arguments(parser, "INPUT");
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    for input in &inputs {
        assert_readable(input)?;
    }
    let dictionary = match argument(parser, "SEQUENCE_DICTIONARY") {
        None => None,
        Some(path) => {
            assert_readable(&path)?;
            Some(dictionary_from_any(&path)?)
        }
    };

    let create_index = flag(parser, "CREATE_INDEX");
    if streamed(
        parser,
        &inputs,
        &output,
        dictionary.as_deref(),
        create_index,
    )
    .is_some()
    {
        return Ok(0);
    }

    let mut texts = Vec::new();
    for input in &inputs {
        let text = vcf_text(input)?;
        texts.push(match (&dictionary, has_contigs(&text)) {
            // `header.setSequenceDictionary(samSequenceDictionary)`.
            (Some(sequences), false) => with_contigs(&text, sequences),
            // `samSequenceDictionary.assertSameDictionary(dict)`, rethrown as the cause of an
            // `IllegalArgumentException`.
            (Some(sequences), true) => {
                if let Err(message) = assert_same_dictionary(&given(sequences), &own(&text)) {
                    return Err(thrown(
                        "java.lang.IllegalArgumentException",
                        format!("java.lang.AssertionError: {message}"),
                    ));
                }
                text
            }
            _ => text,
        });
    }

    let sorted = sort(&texts).map_err(|error| {
        // The library masks the path its message names; the tool prints the input's absolute one.
        let culprit = match error {
            SortError::MissingDictionary => {
                texts.iter().position(|text| !text.contains("\n##contig="))
            }
            SortError::DifferentSamples => {
                let first = samples(&texts[0]);
                texts.iter().position(|text| samples(text) != first)
            }
            _ => None,
        };
        let message = match culprit {
            Some(index) => error
                .message()
                .replace("<masked>", &absolute(&inputs[index])),
            None => error.message(),
        };
        thrown(error.java_class(), message)
    })?;
    write_vcf(parser, &output, &sorted, create_index)?;
    Ok(0)
}

fn has_contigs(text: &str) -> bool {
    text.lines()
        .take_while(|line| line.starts_with("##"))
        .any(|line| line.starts_with("##contig="))
}

/// The run above with the records read a file at a time, each encoded as it is read, and sorted
/// outside memory: in runs of [`RUN`] bytes, each sorted and written beside the output, then
/// merged. `None` on any refusal, having written nothing, so the run above refuses in its own
/// words and order.
///
/// The output is the run above's: its header, and its records as `write_vcf` encodes them under
/// it, in the order of a stable sort on the comparator's key. Each run is sorted stably and holds
/// records read after every earlier run's, and the merge takes the lowest run on a tie, so a tie
/// keeps its input order across runs too.
fn streamed(
    parser: &Parser,
    inputs: &[String],
    output: &str,
    dictionary: Option<&[htsjdk_bam::header::SequenceRecord]>,
    create_index: bool,
) -> Option<()> {
    let mut checks = Checks::new();
    let (mut readers, mut files) = (Vec::new(), Vec::new());
    for input in inputs {
        let mut lines = vcf_lines(input)?;
        let (text, header_lines) = lines.read_header().ok()?;
        let text = match (dictionary, has_contigs(&text)) {
            (Some(sequences), false) => with_contigs(&text, sequences),
            (Some(sequences), true) => {
                assert_same_dictionary(&given(sequences), &own(&text)).ok()?;
                text
            }
            _ => text,
        };
        let file = htsjdk_vcf::reader::read_vcf(&text).ok()?;
        checks.add(file.header.clone()).ok()?;
        readers.push((lines, header_lines));
        files.push(file);
    }
    let (header, _comparator) = checks.finish().ok()?;
    let indexes = contig_indexes(&header);

    let mut stream = VcfStream::create(parser, output, &header, create_index, None)?;
    let mut sorter = Sorter::new(output, RUN);
    for ((lines, header_lines), file) in readers.into_iter().zip(&files) {
        for record in Records::new(Body::Stream(lines), file, header_lines) {
            let record = record.ok()?;
            let contig = *indexes.get(&record.contig)?;
            let from = sorter.text.len();
            stream
                .encoder()
                .encode_into(&record, &mut sorter.text)
                .ok()?;
            sorter.push(contig, record.start, from)?;
        }
    }
    sorter.drain(&mut |line| stream.write_line(line))?;
    stream.finish(parser)
}

/// How much encoded text a run holds before it is sorted and written out.
const RUN: usize = 128 << 20;

/// The records' lines and their keys, and the runs already written beside the output, which are
/// removed when it is dropped.
struct Sorter {
    text: String,
    /// Contig index, start, and where the line (with its newline) is in `text`.
    keys: Vec<(i32, i64, usize, usize)>,
    base: String,
    runs: Vec<String>,
    run: usize,
}

impl Sorter {
    fn new(output: &str, run: usize) -> Self {
        Sorter {
            text: String::new(),
            keys: Vec::new(),
            base: output.to_string(),
            runs: Vec::new(),
            run,
        }
    }

    /// The line encoded into `text` from `from` on.
    fn push(&mut self, contig: i32, start: i64, from: usize) -> Option<()> {
        self.text.push('\n');
        self.keys.push((contig, start, from, self.text.len()));
        if self.text.len() >= self.run {
            self.spill()?;
        }
        Some(())
    }

    /// shortcut: the key is compared as a tuple, which is the comparator's subtraction for any
    /// start below 2^31; a longer contig would need its `as i32` wrap reproduced.
    fn sort(&mut self) {
        self.keys
            .sort_by_key(|&(contig, start, _, _)| (contig, start));
    }

    fn spill(&mut self) -> Option<()> {
        self.sort();
        let path = format!("{}.sort{}.tmp", self.base, self.runs.len());
        self.runs.push(path.clone());
        let mut out = std::io::BufWriter::new(std::fs::File::create(&path).ok()?);
        for &(contig, start, from, to) in &self.keys {
            out.write_all(&contig.to_le_bytes()).ok()?;
            out.write_all(&start.to_le_bytes()).ok()?;
            out.write_all(&((to - from) as u64).to_le_bytes()).ok()?;
            out.write_all(&self.text.as_bytes()[from..to]).ok()?;
        }
        out.into_inner().ok()?;
        self.text.clear();
        self.keys.clear();
        Some(())
    }

    /// Every line, in order, to `sink`.
    fn drain(&mut self, sink: &mut dyn FnMut(&str) -> Option<()>) -> Option<()> {
        if self.runs.is_empty() {
            self.sort();
            for &(_, _, from, to) in &self.keys {
                sink(&self.text[from..to])?;
            }
            return Some(());
        }
        if !self.keys.is_empty() {
            self.spill()?;
        }
        let mut readers = Vec::new();
        for path in &self.runs {
            readers.push(std::io::BufReader::new(std::fs::File::open(path).ok()?));
        }
        let mut heads: Vec<Option<(i32, i64, String)>> = Vec::new();
        for reader in &mut readers {
            heads.push(next_line(reader)?);
        }
        let mut heap = std::collections::BinaryHeap::new();
        for (run, head) in heads.iter().enumerate() {
            if let Some((contig, start, _)) = head {
                heap.push(std::cmp::Reverse((*contig, *start, run)));
            }
        }
        while let Some(std::cmp::Reverse((_, _, run))) = heap.pop() {
            let (_, _, line) = heads[run].take()?;
            sink(&line)?;
            heads[run] = next_line(&mut readers[run])?;
            if let Some((contig, start, _)) = &heads[run] {
                heap.push(std::cmp::Reverse((*contig, *start, run)));
            }
        }
        Some(())
    }
}

/// A run's next line and its key: `Some(None)` at its end, `None` when it cannot be read.
fn next_line(reader: &mut impl BufRead) -> Option<Option<(i32, i64, String)>> {
    let mut contig = [0u8; 4];
    if reader.fill_buf().ok()?.is_empty() {
        return Some(None);
    }
    reader.read_exact(&mut contig).ok()?;
    let (mut start, mut length) = ([0u8; 8], [0u8; 8]);
    reader.read_exact(&mut start).ok()?;
    reader.read_exact(&mut length).ok()?;
    let mut line = vec![0u8; u64::from_le_bytes(length) as usize];
    reader.read_exact(&mut line).ok()?;
    Some(Some((
        i32::from_le_bytes(contig),
        i64::from_le_bytes(start),
        String::from_utf8(line).ok()?,
    )))
}

impl Drop for Sorter {
    fn drop(&mut self) {
        for path in &self.runs {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// `getSampleNamesInOrder()`, which is sorted.
fn samples(text: &str) -> Vec<String> {
    let mut names: Vec<String> = text
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .map(|line| line.split('\t').skip(9).map(str::to_string).collect())
        .unwrap_or_default();
    names.sort();
    names
}

/// The input's own dictionary, from its contig lines.
fn own(text: &str) -> Vec<Sequence> {
    htsjdk_vcf::reader::read_vcf(text)
        .ok()
        .and_then(|file| dictionary(&file.header))
        .unwrap_or_default()
}

/// The given dictionary, as `SAMSequenceRecord`s print.
fn given(sequences: &[htsjdk_bam::header::SequenceRecord]) -> Vec<Sequence> {
    sequences
        .iter()
        .enumerate()
        .map(|(index, record)| Sequence {
            name: record.name.clone(),
            length: i64::from(record.length),
            index,
            assembly: record.attributes.get("AS").map(str::to_string),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Whatever the run size, the lines come back in a stable sort's order, ties across runs
    /// included, and the runs are removed.
    #[test]
    fn runs_merge_as_a_stable_sort() {
        let keys = [(1, 5), (0, 9), (1, 5), (0, 1), (1, 2), (0, 9), (1, 5)];
        let mut expected: Vec<(usize, (i32, i64))> = keys.iter().copied().enumerate().collect();
        expected.sort_by_key(|&(_, key)| key);
        let expected: Vec<String> = expected.iter().map(|(i, _)| format!("{i}\n")).collect();
        let base = std::env::temp_dir().join(format!("sortvcf-{}", std::process::id()));
        let base = base.to_str().unwrap().to_string();
        for run in [1, 5, 1 << 20] {
            let mut sorter = Sorter::new(&base, run);
            for (i, (contig, start)) in keys.iter().enumerate() {
                let from = sorter.text.len();
                sorter.text.push_str(&i.to_string());
                sorter.push(*contig, *start, from).unwrap();
            }
            let mut got = Vec::new();
            sorter
                .drain(&mut |line| {
                    got.push(line.to_string());
                    Some(())
                })
                .unwrap();
            assert_eq!(got, expected, "run size {run}");
            let runs = sorter.runs.clone();
            drop(sorter);
            assert!(runs.iter().all(|path| !std::path::Path::new(path).exists()));
        }
    }
}
