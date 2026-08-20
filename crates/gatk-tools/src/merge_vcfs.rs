//! `MergeVcfs`, ported from `picard.vcf.MergeVcfs` (Picard 3.4.0).
//!
//! Several already-sorted VCFs merged by a heap, under a header smart-merged from all of them. The
//! sibling of [`crate::sort_vcf`], and three of its differences are worth naming.
//!
//! # An input that is not sorted is a refusal, not a repair
//!
//! `MergingIterator` checks each iterator's next element against the last one it emitted and throws
//! when the order breaks. The message names the comparator class and no file at all, so an input of
//! twenty files gives no clue which one was wrong.
//!
//! # A tie is not decided by the order the inputs were given
//!
//! `SortVcf` adds every record to one collection and its sort is stable, so its ties follow the
//! input order. Here the heap is over ITERATORS, and the tie falls out of the heap's own discipline
//! rather than out of the input order: both orders of the same pair write the same file, which the
//! golden's `two-files` and `reversed` pin.
//!
//! Reproducing that needs the heap itself, not a stable merge. This port carries a small
//! `java.util.PriorityQueue`: an array-backed binary heap whose `offer` sifts up while the new
//! element is STRICTLY smaller than its parent and whose `poll` moves the last element to the root
//! and sifts down. On a tie the element already nearer the root stays there, which is why the
//! stream that reached the root first keeps winning whatever order the files were given in.
//!
//! # Two comments collapse into one
//!
//! Every `CO=` is added as a header line keyed `MergeVcfs.comment`, and the smart merge keys
//! unstructured lines by their key, so a second note silently replaces nothing and disappears.
//!
//! # The contig check is about indices
//!
//! `isCompatible` requires every contig a later file declares to sit at the SAME index as in the
//! first. A file declaring a subset therefore fails as surely as one declaring a reordering: the
//! subset shifts the index of everything after what it dropped.

use htsjdk_vcf::comparator::VariantContextComparator;
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::merge::{smart_merge_headers, Source};
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::VariantContext;
use htsjdk_vcf::vcf_file::write_vcf;

/// The key every `CO=` becomes.
pub const COMMENT_KEY: &str = "MergeVcfs.comment";

/// What the tool refuses. The two that name a path carry the golden's mask in its place.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MergeError {
    /// An input with no contig lines and no dictionary supplied.
    MissingDictionary,
    /// A later input whose contigs are not at the same indices as the first's.
    IncompatibleContigs,
    /// A later input whose sorted sample names are not the first's.
    DifferentSamples,
    /// An input whose own records are out of order.
    UnsortedInput,
    /// The reader, the merge or the writer refused.
    Vcf(String, String),
}

impl MergeError {
    pub fn java_class(&self) -> &str {
        match self {
            MergeError::MissingDictionary
            | MergeError::IncompatibleContigs
            | MergeError::DifferentSamples => "java.lang.IllegalArgumentException",
            MergeError::UnsortedInput => "java.lang.IllegalStateException",
            MergeError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            MergeError::MissingDictionary => "A sequence dictionary must be available (either \
                 through the input file or by setting it explicitly)."
                .to_string(),
            MergeError::IncompatibleContigs => {
                "The contig entries in input path <masked> are not compatible with the others."
                    .to_string()
            }
            MergeError::DifferentSamples => {
                "Input path <masked> has sample entries that don't match the other files."
                    .to_string()
            }
            MergeError::UnsortedInput => {
                "The elements of the input Iterators are not sorted \
                 according to the comparator htsjdk.variant.variantcontext.VariantContextComparator"
                    .to_string()
            }
            MergeError::Vcf(_, message) => message.clone(),
        }
    }
}

/// The contig lines of a header, in order.
fn contig_lines(header: &VcfHeader) -> Vec<HeaderLine> {
    header
        .lines
        .iter()
        .filter(|line| matches!(line, HeaderLine::Contig { .. }))
        .cloned()
        .collect()
}

/// The `(id, index)` pairs a contig list declares.
fn contig_indices(lines: &[HeaderLine]) -> Vec<(String, i32)> {
    lines
        .iter()
        .filter_map(|line| match line {
            HeaderLine::Contig { index, fields } => {
                let id = fields
                    .iter()
                    .find(|(key, _)| key == "ID")
                    .map(|(_, value)| value.clone())
                    .unwrap_or_default();
                Some((id, *index))
            }
            _ => None,
        })
        .collect()
}

/// `isCompatible`: every contig this file declares must sit at the same index as in the first.
fn is_compatible(first: &[(String, i32)], later: &[(String, i32)]) -> bool {
    later.iter().all(|(id, index)| {
        first
            .iter()
            .any(|(known, known_index)| known == id && known_index == index)
    })
}

/// `doWork()`: every input's text in, one text out.
pub fn merge(inputs: &[String], comments: &[String]) -> Result<String, MergeError> {
    let mut headers: Vec<VcfHeader> = Vec::new();
    let mut streams: Vec<Vec<VariantContext>> = Vec::new();
    let mut samples: Option<Vec<String>> = None;
    let mut first_contigs: Option<Vec<(String, i32)>> = None;

    for input in inputs {
        let file = read_vcf(input).map_err(|failure| {
            MergeError::Vcf(failure.error.class().to_string(), failure.error.message())
        })?;
        let contigs = contig_lines(&file.header);
        if contigs.is_empty() {
            return Err(MergeError::MissingDictionary);
        }
        let indices = contig_indices(&contigs);
        match &first_contigs {
            None => first_contigs = Some(indices),
            Some(first) => {
                if !is_compatible(first, &indices) {
                    return Err(MergeError::IncompatibleContigs);
                }
            }
        }
        let mut names = file.header.samples.clone();
        names.sort();
        match &samples {
            None => samples = Some(names),
            Some(first) if *first != names => return Err(MergeError::DifferentSamples),
            Some(_) => {}
        }
        streams.push(file.records.clone());
        // `new LinkedHashSet<>(...)`: an input whose header is one the run already has adds
        // nothing to the merge.
        if !headers.contains(&file.header) {
            headers.push(file.header);
        }
    }

    // `COMMENT.forEach(... addMetaDataLine)` on the FIRST header only. Two comments share a key,
    // and the merge keeps one line per key, so only the first reaches the output.
    if let Some(header) = headers.first_mut() {
        for comment in comments {
            header.lines.push(HeaderLine::Unstructured {
                key: COMMENT_KEY.to_string(),
                value: comment.clone(),
            });
        }
    }

    let sources: Vec<Source> = headers
        .iter()
        .map(|header| Source {
            header,
            version: None,
        })
        .collect();
    let (merged, _warnings) = smart_merge_headers(&sources, false)
        .map_err(|error| MergeError::Vcf(error.class().to_string(), error.message().to_string()))?;
    let mut header = VcfHeader {
        lines: merged,
        samples: samples.unwrap_or_default(),
    };
    header.samples.sort();

    let comparator = VariantContextComparator::from_header_lines(&contig_lines(&header))
        .map_err(|error| MergeError::Vcf(error.class().to_string(), error.message().to_string()))?;

    // Each input must already be sorted, which the merging iterator checks as it goes.
    for stream in &streams {
        for pair in stream.windows(2) {
            let order = comparator
                .compare(&pair[0], &pair[1])
                .map_err(|_| MergeError::UnsortedInput)?;
            if order > 0 {
                return Err(MergeError::UnsortedInput);
            }
        }
    }

    // `new MergingIterator<>(comparator, iterators)`: a PriorityQueue of iterators, ordered by
    // their peeked head. The queue is Java's, and so is the tie-breaking, which is why it is
    // written out rather than replaced by a stable merge.
    let mut positions = vec![0usize; streams.len()];
    let mut heap: Vec<usize> = Vec::new();
    let head = |stream: usize, positions: &[usize]| -> Option<VariantContext> {
        streams[stream].get(positions[stream]).cloned()
    };
    let order = |left: usize, right: usize, positions: &[usize]| -> Result<i32, MergeError> {
        let first = head(left, positions).expect("a live stream");
        let second = head(right, positions).expect("a live stream");
        comparator
            .compare(&first, &second)
            .map_err(|_| MergeError::UnsortedInput)
    };

    // `offer`, which sifts up while the new element is strictly less than its parent.
    for stream in 0..streams.len() {
        if head(stream, &positions).is_none() {
            continue;
        }
        heap.push(stream);
        let mut child = heap.len() - 1;
        while child > 0 {
            let parent = (child - 1) / 2;
            if order(heap[child], heap[parent], &positions)? >= 0 {
                break;
            }
            heap.swap(child, parent);
            child = parent;
        }
    }

    let mut records = Vec::new();
    while !heap.is_empty() {
        let stream = heap[0];
        records.push(streams[stream][positions[stream]].clone());
        positions[stream] += 1;

        // `poll` then `add`: the root leaves, the last element takes its place and sifts down, and
        // the advanced iterator is offered again if it still has a record.
        let last = heap.pop().expect("a non-empty heap");
        if !heap.is_empty() {
            heap[0] = last;
            let mut parent = 0;
            loop {
                let left = 2 * parent + 1;
                if left >= heap.len() {
                    break;
                }
                let right = left + 1;
                let mut smallest = left;
                if right < heap.len() && order(heap[right], heap[left], &positions)? < 0 {
                    smallest = right;
                }
                if order(heap[smallest], heap[parent], &positions)? >= 0 {
                    break;
                }
                heap.swap(parent, smallest);
                parent = smallest;
            }
        }
        if head(stream, &positions).is_some() {
            heap.push(stream);
            let mut child = heap.len() - 1;
            while child > 0 {
                let parent = (child - 1) / 2;
                if order(heap[child], heap[parent], &positions)? >= 0 {
                    break;
                }
                heap.swap(child, parent);
                child = parent;
            }
        }
    }

    write_vcf(&header, &records).map_err(|error| {
        MergeError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!("{error:?}"),
        )
    })
}
