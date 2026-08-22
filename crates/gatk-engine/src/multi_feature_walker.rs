//! `MultiFeatureWalker`, ported from `org.broadinstitute.hellbender.engine.MultiFeatureWalker`
//! (GATK 4.6.2.0).
//!
//! Several locus-sorted feature files presented as one sorted stream. It is the class under
//! `SiteDepthtoBAF`, `CondenseDepthEvidence`, `PrintSVEvidence` and the rest of the SV family, so
//! its merge order is their file order.
//!
//! # The merge is a heap, and a heap is not stable
//!
//! `MergingIterator` holds a `java.util.PriorityQueue`. Two features that compare equal come out
//! in the order the heap's shape decides, which is neither the order their files were named in nor
//! anything the comparator says: three files each holding the same interval, offered in the order
//! they were named, come out **first, third, second**. [`JavaHeap`] is `siftUp` and `siftDown`
//! transcribed, for the same reason `OverhangFixingManager` needed them.
//!
//! # A contig the dictionary does not name is not refused, it is misdiagnosed
//!
//! ```java
//! this.contigIndex = context.getDictionary().getSequenceIndex(feature.getContig());
//! ```
//!
//! `getSequenceIndex` answers `-1` for a contig it does not hold, and `-1` sorts before every named
//! contig. So a file holding chr1 then chr2, under a dictionary naming chr1 and chr3, is refused
//! with `inputs are not sorted at chr2:101`. The file is sorted. The dictionary is incomplete. The
//! message says neither.
//!
//! That path needs the input to hold a named contig too: with no overlap at all the run never
//! reaches here, the engine's own dictionary comparison refusing it first as
//! `IncompatibleSequenceDictionaries`, which is not this class's doing and is not ported here.
//!
//! # The sort check fires late, and against the same input
//!
//! `next()` polls, then draws the replacement from the **same** context and compares it against
//! the entry it just returned. So a file that goes backwards is caught only when its next record
//! is pulled, and the message names the locus of the new feature rather than of the one it should
//! have followed. It is also unreachable through an indexed file, since an unsorted file cannot be
//! indexed at all.
//!
//! # On equal sizes the new dictionary is the smaller one
//!
//! ```java
//! if ( newDict.getDictionary().size() <= curDict.getDictionary().size() ) { smallDict = newDict; ... }
//! ```
//!
//! The `<=` decides which source each half of the two refusal messages names when the two
//! dictionaries are the same size, and the reference is offered as `newDict`, so it is the one
//! called small.

use std::cmp::Ordering;

/// One feature, reduced to what the merge reads and what the tool prints.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Located {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    /// `Feature.toString()`, which is what `ExampleMultiFeatureWalker.apply` prints.
    pub text: String,
}

/// A dictionary and the argument it was read from, which the refusals name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DictSource {
    pub contigs: Vec<String>,
    pub source: String,
}

impl DictSource {
    /// `SAMSequenceDictionary.getSequenceIndex`, which answers -1 rather than refusing.
    pub fn sequence_index(&self, contig: &str) -> i32 {
        self.contigs
            .iter()
            .position(|name| name == contig)
            .map_or(-1, |index| index as i32)
    }
}

/// What the walk refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WalkerError {
    /// No dictionary anywhere: a `.rd.txt` header carries sample names and a null dictionary.
    NoDictionary,
    /// A contig of the smaller dictionary is absent from the larger.
    ContigAbsent {
        contig: String,
        small_source: String,
        large_source: String,
    },
    /// The two dictionaries hold the same contigs in a different order.
    ContigsOutOfOrder {
        contig: String,
        previous: String,
        large_source: String,
        small_source: String,
    },
    /// An input went backwards, or a contig it names is missing from the dictionary.
    NotSorted { contig: String, start: i32 },
}

impl WalkerError {
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    pub fn message(&self) -> String {
        match self {
            WalkerError::NoDictionary => {
                // Two spaces after the full stop, as the reference writes it.
                "No dictionary found.  Provide one as --sequence-dictionary or --reference."
                    .to_string()
            }
            WalkerError::ContigAbsent {
                contig,
                small_source,
                large_source,
            } => format!(
                "Contig {contig} in the dictionary read from {small_source} does not appear in \
                 the larger dictionary read from {large_source}"
            ),
            WalkerError::ContigsOutOfOrder {
                contig,
                previous,
                large_source,
                small_source,
            } => format!(
                "Contigs out of order: Contig {contig} comes before contig {previous} in the \
                 dictionary read from {large_source}, but follows it in the dictionary read from \
                 {small_source}"
            ),
            WalkerError::NotSorted { contig, start } => {
                format!("inputs are not sorted at {contig}:{start}")
            }
        }
    }
}

/// `betterDictionary`: the larger of the two, once the smaller is known to be a subset of it in
/// the same relative order.
pub fn better_dictionary(
    new_dict: Option<DictSource>,
    current: Option<DictSource>,
) -> Result<Option<DictSource>, WalkerError> {
    let (new_dict, current) = match (new_dict, current) {
        (new_dict, None) => return Ok(new_dict),
        (None, current) => return Ok(current),
        (Some(new_dict), Some(current)) => (new_dict, current),
    };
    // The `<=` is what decides which source each refusal names when the two are the same size.
    let (small, large) = if new_dict.contigs.len() <= current.contigs.len() {
        (new_dict, current)
    } else {
        (current, new_dict)
    };
    let mut last_index: i32 = -1;
    for contig in &small.contigs {
        let index = large.sequence_index(contig);
        if index == -1 {
            return Err(WalkerError::ContigAbsent {
                contig: contig.clone(),
                small_source: small.source.clone(),
                large_source: large.source.clone(),
            });
        }
        if index <= last_index {
            return Err(WalkerError::ContigsOutOfOrder {
                contig: contig.clone(),
                previous: large.contigs[last_index as usize].clone(),
                large_source: large.source.clone(),
                small_source: small.source.clone(),
            });
        }
        last_index = index;
    }
    Ok(Some(large))
}

/// `setDictionaryAndSamples`, for the sources this walk can be given: the master dictionary, then
/// the reference, in that order.
pub fn choose_dictionary(
    master: Option<DictSource>,
    reference: Option<DictSource>,
) -> Result<DictSource, WalkerError> {
    let chosen = better_dictionary(reference, master)?;
    chosen.ok_or(WalkerError::NoDictionary)
}

/// One entry of the queue: a feature, its contig's index, and the input it came from.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Entry {
    input: usize,
    contig_index: i32,
    feature: Located,
}

impl Entry {
    /// `PQEntry.compareTo`: contig index, then start, then end. The index rather than the name, so
    /// the dictionary's order is the file's.
    fn compare(&self, other: &Entry) -> Ordering {
        self.contig_index
            .cmp(&other.contig_index)
            .then(self.feature.start.cmp(&other.feature.start))
            .then(self.feature.end.cmp(&other.feature.end))
    }
}

/// `java.util.PriorityQueue`, transcribed for this element type.
///
/// The order two equal elements come out in is a property of the heap's shape, so a stable merge
/// would produce a different, equally sorted stream and no golden could hold.
#[derive(Debug, Default)]
struct JavaHeap {
    queue: Vec<Entry>,
}

impl JavaHeap {
    fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }

    /// `offer(e)`: append, then sift it up. Equal stops, which is why the heap is not stable.
    fn push(&mut self, element: Entry) {
        let mut k = self.queue.len();
        self.queue.push(element);
        while k > 0 {
            let parent = (k - 1) >> 1;
            if self.queue[k].compare(&self.queue[parent]) != Ordering::Less {
                break;
            }
            self.queue.swap(k, parent);
            k = parent;
        }
    }

    /// `poll()`: take the head, move the last element to the root and sift it down.
    fn poll(&mut self) -> Option<Entry> {
        if self.queue.is_empty() {
            return None;
        }
        let head = self.queue.swap_remove(0);
        self.sift_down(0);
        Some(head)
    }

    fn sift_down(&mut self, mut k: usize) {
        let size = self.queue.len();
        let half = size >> 1;
        while k < half {
            let mut child = 2 * k + 1;
            let right = child + 1;
            if right < size && self.queue[child].compare(&self.queue[right]) == Ordering::Greater {
                child = right;
            }
            if self.queue[k].compare(&self.queue[child]) != Ordering::Greater {
                break;
            }
            self.queue.swap(k, child);
            k = child;
        }
    }
}

/// `MergingIterator` run to exhaustion: every feature in the order the walk hands it over.
///
/// The inputs are offered in the order given, which is `FeatureManager.getAllInputs()`'s order.
pub fn merge(
    inputs: &[Vec<Located>],
    dictionary: &DictSource,
) -> Result<Vec<Located>, WalkerError> {
    let mut cursors: Vec<usize> = vec![0; inputs.len()];
    let mut heap = JavaHeap::default();
    for (index, input) in inputs.iter().enumerate() {
        if let Some(feature) = input.first() {
            cursors[index] = 1;
            heap.push(entry_for(index, feature, dictionary));
        }
    }
    let mut written = Vec::new();
    while !heap.is_empty() {
        let entry = heap.poll().expect("a non-empty heap");
        // The replacement is drawn from the SAME input, and compared against what was just
        // returned rather than against the head of the queue.
        let input = entry.input;
        if let Some(feature) = inputs[input].get(cursors[input]) {
            cursors[input] += 1;
            let replacement = entry_for(input, feature, dictionary);
            let compared = replacement.compare(&entry);
            heap.push(replacement);
            if compared == Ordering::Less {
                return Err(WalkerError::NotSorted {
                    contig: feature.contig.clone(),
                    start: feature.start,
                });
            }
        }
        written.push(entry.feature);
    }
    Ok(written)
}

fn entry_for(input: usize, feature: &Located, dictionary: &DictSource) -> Entry {
    Entry {
        input,
        contig_index: dictionary.sequence_index(&feature.contig),
        feature: feature.clone(),
    }
}
