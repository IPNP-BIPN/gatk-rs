//! `MultiVariantWalker` and `MultiVariantDataSource`, ported from
//! `org.broadinstitute.hellbender.engine` (GATK 4.7.0.0, unchanged since 4.6.2.0), with the two
//! grouping walkers on top of them, `MultiVariantWalkerGroupedOnStart` and
//! `MultiVariantWalkerGroupedByOverlap`.
//!
//! Eleven tools stand on this: CombineGVCFs, ReblockGVCF, VariantEval, ApplyVQSR,
//! VariantRecalibrator, FilterAlignmentArtifacts, VCFComparator, SVStratify, SVCluster,
//! GroupedSVCluster and JointGermlineCNVSegmentation. What they share is everything that happens
//! before their `apply` sees a record, and that is what is here. The `multi-variant-walker` suite
//! measures it through three probe walkers.
//!
//! # Whose dictionary
//!
//! Each input's dictionary is its header's `##contig` lines, or, when it has none, the contig
//! names of its Tribble index at length 0, or nothing. The data source cross-checks them pairwise
//! (`validateAllSequenceDictionaries`): the first input to name a contig owns it, and every later
//! input naming it must agree on the length unless either side is 0. **The md5 branch is
//! unreachable from a VCF**: a `##contig` line's `md5` is not carried into the record the check
//! reads, so two inputs whose `md5`s differ are accepted, which the suite shows.
//!
//! The merged dictionary is the merged header's contig lines in `smartMergeHeaders`' first-seen
//! order, so two inputs naming {chr1, chr2} and {chr1, chr3} merge to chr1, chr2, chr3 and the
//! same two the other way round to chr1, chr3, chr2. An index-derived dictionary is written into
//! its header before the merge, so its length-0 lines can win.
//!
//! # The merged header loses contig lines when it is written
//!
//! The contig lines keep the index they had in their own file, and `getMetaDataInSortedOrder` is a
//! `TreeSet` whose contig comparison is by index alone. So chr2 from one file and chr3 from another,
//! both at index 1, are the same element: the first one inserted stays and the other vanishes from
//! every header written from the merge, while the dictionary still holds both. [`sorted_metadata`]
//! reproduces that.
//!
//! # The merge order
//!
//! htsjdk's `MergingIterator` over `VariantContextComparator`: contig index in the merged
//! dictionary, then start, through a `java.util.PriorityQueue` of iterators, which is a binary heap
//! and not stable. Three inputs with a record at the same position come out 1, 3, 2. A contig the
//! merged dictionary does not name is a `NullPointerException` at the first comparison that meets
//! it, and an input that goes backwards is an `IllegalStateException` when the record that went
//! backwards is handed out. Both arrive after whatever was handed out before them.

use std::collections::BTreeSet;

use htsjdk_vcf::header::{HeaderLine, VcfHeader};

use crate::sequence_dictionary;

/// One contig of a dictionary: what the cross-check reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    /// `SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH` is 0.
    pub length: i64,
}

/// What the data source and the traversal refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// `UserException.BadInput`.
    BadInput(String),
    /// `UserException.IncompatibleSequenceDictionaries`, the whole message.
    Incompatible(String),
    /// More than one input and no dictionary anywhere.
    NoDictionary,
    /// Anything else, by the class the reference throws.
    Other { class: String, message: String },
}

impl Refusal {
    pub fn java_class(&self) -> &str {
        match self {
            Refusal::BadInput(_) => "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
            Refusal::Incompatible(_) => {
                "org.broadinstitute.hellbender.exceptions.UserException$IncompatibleSequenceDictionaries"
            }
            Refusal::NoDictionary => "org.broadinstitute.hellbender.exceptions.UserException",
            Refusal::Other { class, .. } => class,
        }
    }

    /// Whether the reference reports it as a user error.
    pub fn is_user(&self) -> bool {
        !matches!(self, Refusal::Other { .. })
    }

    pub fn message(&self) -> String {
        match self {
            Refusal::BadInput(message) => format!("Bad input: {message}"),
            Refusal::Incompatible(message) => message.clone(),
            Refusal::NoDictionary => "No sequence dictionary was found for any input. When using \
                multiple inputs, at least one input must have a sequence dictionary, or an index \
                from which a sequence dictionary can be derived."
                .to_string(),
            Refusal::Other { message, .. } => message.clone(),
        }
    }
}

/// One `--variant`, as `FeatureInput` holds it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureInput {
    /// The `-V:tag` name, when the user gave one.
    pub tag: Option<String>,
    /// The argument as written, which equality compares.
    pub raw: String,
    /// The path made absolute, which is the name when there is no tag.
    pub absolute: String,
}

impl FeatureInput {
    /// `getName`: the tag, or the absolute path.
    pub fn name(&self) -> &str {
        self.tag.as_deref().unwrap_or(&self.absolute)
    }

    /// `toString`: `tag:path`, or the path.
    pub fn display(&self) -> String {
        match &self.tag {
            Some(tag) => format!("{tag}:{}", self.absolute),
            None => self.absolute.clone(),
        }
    }
}

/// `initializeDrivingVariants`: a repeated input is refused before any is opened. The same file
/// under two tags is two inputs, and is read twice.
pub fn require_unique(inputs: &[FeatureInput]) -> Result<(), Refusal> {
    for (index, input) in inputs.iter().enumerate() {
        if inputs[..index]
            .iter()
            .any(|other| other.tag == input.tag && other.raw == input.raw)
        {
            return Err(Refusal::BadInput(format!(
                "Feature inputs must be unique: {}",
                input.display()
            )));
        }
    }
    Ok(())
}

/// The `##contig` lines of a header as a dictionary, in the order they appear.
pub fn dictionary_of(lines: &[HeaderLine]) -> Vec<Contig> {
    lines
        .iter()
        .filter_map(|line| match line {
            HeaderLine::Contig { fields, .. } => {
                let field = |key: &str| {
                    fields
                        .iter()
                        .find(|(k, _)| k == key)
                        .map(|(_, value)| value.as_str())
                };
                Some(Contig {
                    name: field("ID")?.to_string(),
                    length: field("length")
                        .and_then(|value| value.parse().ok())
                        .unwrap_or(0),
                })
            }
            _ => None,
        })
        .collect()
}

fn as_records(dictionary: &[Contig]) -> Vec<htsjdk_bam::header::SequenceRecord> {
    dictionary
        .iter()
        .map(|contig| htsjdk_bam::header::SequenceRecord::new(&contig.name, contig.length as i32))
        .collect()
}

/// `validateAllSequenceDictionaries`: every contig checked against the first input that named it.
///
/// `sources` is each input's name and the dictionary its data source reports, `None` where it has
/// none (which the reference only warns about).
pub fn cross_validate(sources: &[(String, Option<Vec<Contig>>)]) -> Result<(), Refusal> {
    // contig name -> the source that first named it
    let mut owner: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for (index, (name, dictionary)) in sources.iter().enumerate() {
        let Some(dictionary) = dictionary else {
            continue;
        };
        for contig in dictionary {
            let Some(&previous) = owner.get(contig.name.as_str()) else {
                owner.insert(&contig.name, index);
                continue;
            };
            let (previous_name, previous_dictionary) = &sources[previous];
            let previous_dictionary = previous_dictionary.as_deref().unwrap_or(&[]);
            let target = previous_dictionary
                .iter()
                .find(|other| other.name == contig.name)
                .expect("the owner names the contig");
            // `sequenceRecordsAreEquivalent`: the same name, and the same length unless either is
            // unknown. The md5 arm after it never fires for a VCF, whose records carry none.
            if contig.length != target.length && contig.length != 0 && target.length != 0 {
                return Err(Refusal::Incompatible(sequence_dictionary::incompatible_message(
                    &format!(
                        "Incompatible sequences found ({}: {}) and ({}: {})",
                        contig.name, contig.length, target.name, target.length
                    ),
                    name,
                    &as_records(dictionary),
                    previous_name,
                    &as_records(previous_dictionary),
                )));
            }
        }
    }
    Ok(())
}

/// One input's header as the data source holds it.
pub struct Source<'a> {
    pub header: &'a VcfHeader,
    /// `getVCFHeaderVersion()`.
    pub version: Option<&'a str>,
    /// The contig names of its Tribble index, used only when the header has no `##contig` line.
    pub index_contigs: Option<Vec<String>>,
}

impl Source<'_> {
    /// `getHeaderWithUpdatedSequenceDictionary`: the header, with the index's dictionary written
    /// into it when it has none of its own.
    fn updated(&self) -> VcfHeader {
        let mut header = self.header.clone();
        let has_contigs = header
            .lines
            .iter()
            .any(|line| matches!(line, HeaderLine::Contig { .. }));
        if let (false, Some(contigs)) = (has_contigs, &self.index_contigs) {
            for (index, contig) in contigs.iter().enumerate() {
                header
                    .lines
                    .push(HeaderLine::contig(contig, 0, index as i32));
            }
        }
        header
    }

    /// The dictionary this input's data source reports: its header's, else its index's.
    pub fn dictionary(&self) -> Option<Vec<Contig>> {
        let own = dictionary_of(&self.header.lines);
        if !own.is_empty() {
            return Some(own);
        }
        self.index_contigs.as_ref().map(|contigs| {
            contigs
                .iter()
                .map(|name| Contig {
                    name: name.clone(),
                    length: 0,
                })
                .collect()
        })
    }
}

/// What the data source reports before the first record.
#[derive(Debug, Clone)]
pub struct Merged {
    /// The merged header's lines, in the merge's order, every contig line kept.
    pub lines: Vec<HeaderLine>,
    /// `getSequenceDictionary`, `None` when no input has one.
    pub dictionary: Option<Vec<Contig>>,
    /// `getSamples`: every input's samples, sorted, each once.
    pub samples: Vec<String>,
}

impl Merged {
    /// The contig names in dictionary order, which the merge compares by.
    pub fn contigs(&self) -> Vec<String> {
        self.dictionary
            .iter()
            .flatten()
            .map(|contig| contig.name.clone())
            .collect()
    }
}

/// The data source's constructor after the cross-check: the merged header, the samples, and the
/// refusal when several inputs leave no dictionary at all.
pub fn merge(sources: &[Source]) -> Result<Merged, Refusal> {
    let updated: Vec<VcfHeader> = sources.iter().map(Source::updated).collect();
    let lines = if updated.len() > 1 {
        let merge_sources: Vec<htsjdk_vcf::merge::Source> = updated
            .iter()
            .zip(sources)
            .map(|(header, source)| htsjdk_vcf::merge::Source {
                header,
                version: source.version,
            })
            .collect();
        htsjdk_vcf::merge::smart_merge_headers(&merge_sources, true)
            .map_err(|error| Refusal::Other {
                class: error.class().to_string(),
                message: error.message().to_string(),
            })?
            .0
    } else {
        updated
            .into_iter()
            .next()
            .map(|header| header.lines)
            .unwrap_or_default()
    };
    let dictionary = dictionary_of(&lines);
    let dictionary = (!dictionary.is_empty()).then_some(dictionary);
    if dictionary.is_none() && sources.len() > 1 {
        return Err(Refusal::NoDictionary);
    }
    let samples: BTreeSet<String> = sources
        .iter()
        .flat_map(|source| source.header.samples.iter().cloned())
        .collect();
    Ok(Merged {
        lines,
        dictionary,
        samples: samples.into_iter().collect(),
    })
}

/// `getMetaDataInSortedOrder` without its prepended `fileformat` line: the writer's order, where a
/// contig line sharing its index with an earlier one is dropped, as the `TreeSet` drops it.
pub fn sorted_metadata(lines: &[HeaderLine]) -> Vec<HeaderLine> {
    let header = VcfHeader {
        lines: lines.to_vec(),
        samples: Vec::new(),
    };
    let mut seen: BTreeSet<i32> = BTreeSet::new();
    header
        .sorted_lines()
        .into_iter()
        .filter(|line| match line {
            HeaderLine::Contig { index, .. } => seen.insert(*index),
            _ => true,
        })
        .cloned()
        .collect()
}

/// The merged lines with the contig lines a writer would drop already dropped, so that a tool
/// writing them through `VcfHeader::write` writes what the reference writes.
pub fn written_lines(lines: &[HeaderLine]) -> Vec<HeaderLine> {
    let kept = sorted_metadata(lines);
    lines
        .iter()
        .filter(|line| match line {
            HeaderLine::Contig { .. } => kept.contains(line),
            _ => true,
        })
        .cloned()
        .collect()
}

/// A record's place in the merge: its contig's index in the merged dictionary, `None` when the
/// dictionary does not name it, and its start.
pub type Key = (Option<i64>, i64);

/// `VariantContextComparator.compare`, which unboxes a null for a contig it does not know.
fn compare(a: Key, b: Key) -> Result<std::cmp::Ordering, Refusal> {
    match (a.0, b.0) {
        (Some(x), Some(y)) => Ok(x.cmp(&y).then(a.1.cmp(&b.1))),
        _ => Err(Refusal::Other {
            class: "java.lang.NullPointerException".to_string(),
            message: "Cannot invoke \"java.lang.Integer.intValue()\" because the return value of \
                      \"java.util.Map.get(Object)\" is null"
                .to_string(),
        }),
    }
}

/// The merge run to exhaustion or to its first failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Merge {
    /// `(input, index)` of every record handed out, in order.
    pub order: Vec<(usize, usize)>,
    /// What stopped it, after `order` was handed out.
    pub failure: Option<Refusal>,
}

/// htsjdk's `MergingIterator` over one key list per input.
///
/// A single input is not merged at all: `MultiVariantDataSource` hands its iterator over as it is,
/// so nothing is compared and neither refusal can happen.
pub fn merging_order(keys: &[Vec<Key>]) -> Merge {
    let mut order = Vec::new();
    if keys.len() == 1 {
        order.extend((0..keys[0].len()).map(|index| (0, index)));
        return Merge {
            order,
            failure: None,
        };
    }
    let failure = merge_into(keys, &mut order).err();
    Merge { order, failure }
}

fn merge_into(keys: &[Vec<Key>], order: &mut Vec<(usize, usize)>) -> Result<(), Refusal> {
    let mut cursors = vec![0usize; keys.len()];
    let mut queue = JavaQueue::default();
    for (input, list) in keys.iter().enumerate() {
        if let Some(key) = list.first() {
            queue.offer((*key, input))?;
        }
    }
    let mut last: Option<Key> = None;
    while let Some((key, input)) = queue.poll()? {
        let index = cursors[input];
        cursors[input] += 1;
        if let Some(previous) = last {
            if compare(previous, key)? == std::cmp::Ordering::Greater {
                return Err(Refusal::Other {
                    class: "java.lang.IllegalStateException".to_string(),
                    message: "The elements of the input Iterators are not sorted according to the \
                              comparator htsjdk.variant.variantcontext.VariantContextComparator"
                        .to_string(),
                });
            }
        }
        if let Some(next) = keys[input].get(cursors[input]) {
            queue.offer((*next, input))?;
        }
        last = Some(key);
        order.push((input, index));
    }
    Ok(())
}

/// `java.util.PriorityQueue` of iterators, ordered by the key alone: `siftUp` and `siftDown`
/// transcribed, since the order of equal keys is the heap's shape.
#[derive(Default)]
struct JavaQueue {
    queue: Vec<(Key, usize)>,
}

impl JavaQueue {
    /// `offer`: append, then `siftUp`, which stops at an equal parent.
    fn offer(&mut self, element: (Key, usize)) -> Result<(), Refusal> {
        let mut k = self.queue.len();
        self.queue.push(element);
        while k > 0 {
            let parent = (k - 1) >> 1;
            if compare(element.0, self.queue[parent].0)? != std::cmp::Ordering::Less {
                break;
            }
            self.queue[k] = self.queue[parent];
            k = parent;
        }
        self.queue[k] = element;
        Ok(())
    }

    /// `poll`: take the head, move the last element to the root, then `siftDown`, which prefers
    /// the left child on a tie and stops at an equal child.
    fn poll(&mut self) -> Result<Option<(Key, usize)>, Refusal> {
        if self.queue.is_empty() {
            return Ok(None);
        }
        let result = self.queue[0];
        let last = self.queue.pop().expect("a non-empty queue");
        let size = self.queue.len();
        if size > 0 {
            let mut k = 0;
            let half = size >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let right = child + 1;
                if right < size
                    && compare(self.queue[child].0, self.queue[right].0)?
                        == std::cmp::Ordering::Greater
                {
                    child = right;
                }
                if compare(last.0, self.queue[child].0)? != std::cmp::Ordering::Greater {
                    break;
                }
                self.queue[k] = self.queue[child];
                k = child;
            }
            self.queue[k] = last;
        }
        Ok(Some(result))
    }
}

/// The key of a record under a merged dictionary.
pub fn key_of(contigs: &[String], contig: &str, start: i64) -> Key {
    (
        contigs
            .iter()
            .position(|name| name == contig)
            .map(|index| index as i64),
        start,
    )
}

/// `MultiVariantWalker.traverse`'s wrapper around whatever `apply` throws.
pub fn apply_failure_message(contig: &str, start: i64, record: &str) -> String {
    format!("Exception thrown at {contig}:{start} {record}")
}

/// A record as the grouping walkers see it.
pub trait Located {
    fn contig(&self) -> &str;
    fn start(&self) -> i64;
    fn end(&self) -> i64;
}

/// `isWithinInterval` over a one-base interval at the start: true with no intervals at all.
pub fn starts_within(intervals: Option<&[(String, i64, i64)]>, contig: &str, start: i64) -> bool {
    intervals.is_none_or(|list| {
        list.iter()
            .any(|(c, s, e)| c == contig && *s <= start && start <= *e)
    })
}

/// `MultiVariantWalkerGroupedOnStart.apply(variant, ...)`: records whose starts chain within
/// `distance` of the last one and lie within `max_distance` of the first.
pub struct GroupedOnStart<T> {
    pub distance: i64,
    pub max_distance: i64,
    current: Vec<T>,
    first_start: i64,
    last_start: i64,
}

impl<T: Located> GroupedOnStart<T> {
    pub fn new(distance: i64, max_distance: i64) -> Self {
        GroupedOnStart {
            distance,
            max_distance,
            current: Vec::new(),
            first_start: 0,
            last_start: 0,
        }
    }

    /// Takes one record; returns the group it closed, if it closed one.
    pub fn accept(&mut self, variant: T) -> Option<Vec<T>> {
        let start = variant.start();
        let mut closed = None;
        if self.current.is_empty() {
            self.first_start = start;
        } else if self.current[0].contig() != variant.contig()
            || self.last_start < start - self.distance
            || self.first_start < start - self.max_distance
        {
            closed = Some(std::mem::take(&mut self.current));
            self.first_start = start;
        }
        self.current.push(variant);
        self.last_start = start;
        closed
    }

    /// `afterTraverse`: the last group, when there is one.
    pub fn finish(&mut self) -> Option<Vec<T>> {
        (!self.current.is_empty()).then(|| std::mem::take(&mut self.current))
    }
}

/// `MultiVariantWalkerGroupedByOverlap.apply(variant, ...)`: a record joins the group while it
/// starts at or before the furthest end seen on the contig. That end is never reset between
/// groups, only when the contig changes.
pub struct GroupedByOverlap<T> {
    current: Vec<T>,
    last_contig: Option<String>,
    last_end: i64,
}

impl<T: Located> Default for GroupedByOverlap<T> {
    fn default() -> Self {
        GroupedByOverlap {
            current: Vec::new(),
            last_contig: None,
            last_end: 0,
        }
    }
}

impl<T: Located> GroupedByOverlap<T> {
    /// Takes one record; returns the group it closed, if it closed one.
    pub fn accept(&mut self, variant: T) -> Option<Vec<T>> {
        let mut closed = None;
        if self.current.is_empty() {
            self.last_contig = Some(variant.contig().to_string());
        } else if self.current[0].contig() != variant.contig() || variant.start() > self.last_end {
            closed = Some(std::mem::take(&mut self.current));
        }
        if variant.end() > self.last_end || self.last_contig.as_deref() != Some(variant.contig()) {
            self.last_end = variant.end();
            self.last_contig = Some(variant.contig().to_string());
        }
        self.current.push(variant);
        closed
    }

    /// `afterTraverse`: the last group, when there is one.
    pub fn finish(&mut self) -> Option<Vec<T>> {
        (!self.current.is_empty()).then(|| std::mem::take(&mut self.current))
    }
}

/// `makeSpanningReferenceContext`: the group's span, padded by `padding` on each side and clipped
/// to the contig, as `ReferenceContext.setWindow` clips it.
pub fn spanning_window<T: Located>(group: &[T], padding: i64, contig_length: i64) -> (i64, i64) {
    let start = group.iter().map(Located::start).min().expect("a group");
    let end = group.iter().map(Located::end).max().expect("a group");
    ((start - padding).max(1), (end + padding).min(contig_length))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(contig: i64, start: i64) -> Key {
        (Some(contig), start)
    }

    #[test]
    fn three_equal_heads_come_out_first_third_second() {
        let keys = vec![vec![key(0, 100)], vec![key(0, 100)], vec![key(0, 100)]];
        assert_eq!(merging_order(&keys).order, vec![(0, 0), (2, 0), (1, 0)]);
    }

    #[test]
    fn an_unknown_contig_fails_at_the_first_comparison() {
        let keys = vec![vec![key(0, 100)], vec![(None, 5)]];
        let merge = merging_order(&keys);
        assert!(merge.order.is_empty());
        assert_eq!(
            merge.failure.map(|f| f.java_class().to_string()),
            Some("java.lang.NullPointerException".to_string())
        );
    }

    #[test]
    fn a_single_input_is_not_compared() {
        let keys = vec![vec![key(0, 300), (None, 5), key(0, 100)]];
        assert_eq!(merging_order(&keys).order, vec![(0, 0), (0, 1), (0, 2)]);
    }

    struct At(&'static str, i64, i64);
    impl Located for At {
        fn contig(&self) -> &str {
            self.0
        }
        fn start(&self) -> i64 {
            self.1
        }
        fn end(&self) -> i64 {
            self.2
        }
    }

    #[test]
    fn overlap_keeps_the_furthest_end_across_groups() {
        let mut grouping = GroupedByOverlap::default();
        let mut groups: Vec<Vec<i64>> = Vec::new();
        for at in [At("c", 100, 109), At("c", 105, 105), At("c", 110, 110), At("c", 111, 111)] {
            if let Some(group) = grouping.accept(at) {
                groups.push(group.iter().map(|a| a.1).collect());
            }
        }
        groups.extend(grouping.finish().map(|g| g.iter().map(|a| a.1).collect()));
        assert_eq!(groups, vec![vec![100, 105], vec![110], vec![111]]);
    }
}
