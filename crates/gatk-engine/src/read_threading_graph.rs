//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading`
//! (`AbstractReadThreadingGraph`, `ReadThreadingGraph`, `MultiDeBruijnVertex`) and
//! `haplotypecaller.graphs` (`BaseGraph`, `BaseEdge`, `MultiSampleEdge`), GATK 4.6.2.0: the read
//! threading graph as `buildGraphIfNecessary` leaves it. Pruning, dangling-end recovery and the
//! conversion to a sequence graph are later stages and are not here.
//!
//! # What is built
//!
//! The reference haplotype is threaded first, then each read's runs of usable bases, sample by
//! sample in the order the samples were first seen. A k-mer that occurs twice within any ONE
//! sequence is non-unique: it is never a merge point, so each occurrence gets its own vertex. Every
//! other k-mer is one vertex, shared by all the sequences that contain it, and an edge's
//! multiplicity counts the sequences that walked it.
//!
//! # Order is part of the answer
//!
//! The reference graph is a JGraphT `DefaultDirectedGraph` (jgrapht-core 1.1.0). Vertices and edges
//! compare by identity, the vertex and edge sets are `LinkedHashMap`s, and each vertex's incoming
//! and outgoing edges are an array-backed set in insertion order. `extendChainByOne` follows the
//! FIRST outgoing edge whose target ends in the next base, and the backwards count increase walks
//! the incoming edges in their order, so insertion order decides the graph. Here every collection
//! is a `Vec` in insertion order and a vertex or an edge is its index.
//!
//! The one hash-ordered collection, the set of non-unique k-mers, is only ever asked `contains`.
//!
//! # Two quirks kept on purpose
//!
//! `determineNonUniqueKmers` counts k-mers from the start of the WHOLE sequence, not from the
//! sequence's own `start`, so the bases of a read before its first usable run still decide which
//! k-mers are non-unique. And `findStart` searches `start .. stop - k` exclusive, so the last
//! k-mer of a sequence is never its threading start. Both are the reference's behaviour.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};

/// `ANONYMOUS_SAMPLE`: the sample every sequence added without one belongs to, the reference's
/// included.
pub const ANONYMOUS_SAMPLE: &str = "XXX_UNNAMED_XXX";

/// A refusal the reference throws while building.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphError {
    /// `IllegalStateException("Found two refSources! ...")`.
    TwoRefSources,
    /// `Utils.validate(!(isRef && uniqueMergeVertex != null), "Found a unique vertex to merge into
    /// the reference graph ...")`.
    UniqueVertexOnReference,
    /// `IllegalStateException("Attempting to add sequence to a graph that has already been built")`.
    AlreadyBuilt,
    /// `SequenceForKmers`'s argument checks.
    InvalidSequence(String),
}

/// `MultiDeBruijnVertex`: a k-mer's bases and the debugging text `buildGraphIfNecessary` appends to.
#[derive(Debug, Clone)]
pub struct Vertex {
    pub sequence: Vec<u8>,
    pub additional_info: String,
}

impl Vertex {
    /// `getSuffix()`: the last base.
    pub fn suffix(&self) -> u8 {
        self.sequence[self.sequence.len() - 1]
    }
}

/// `MultiSampleEdge`: a `BaseEdge` (multiplicity, reference flag) that also keeps the largest
/// per-sample multiplicities, up to `numPruningSamples` of them, for the pruner.
#[derive(Debug, Clone)]
pub struct Edge {
    pub source: usize,
    pub target: usize,
    pub is_ref: bool,
    pub multiplicity: i32,
    current_single_sample: i32,
    /// `PriorityQueue<Integer>`, a min-heap: the smallest is what `poll` drops and `peek` reads.
    single_sample: BinaryHeap<Reverse<i32>>,
    capacity: usize,
}

impl Edge {
    fn new(source: usize, target: usize, is_ref: bool, multiplicity: i32, capacity: usize) -> Edge {
        let mut single_sample = BinaryHeap::with_capacity(capacity + 1);
        single_sample.push(Reverse(multiplicity));
        Edge {
            source,
            target,
            is_ref,
            multiplicity,
            current_single_sample: multiplicity,
            single_sample,
            capacity,
        }
    }

    /// `incMultiplicity`: the total and the current sample's count.
    fn inc_multiplicity(&mut self, increment: i32) {
        self.multiplicity += increment;
        self.current_single_sample += increment;
    }

    /// `flushSingleSampleMultiplicity`: the current sample's count joins the queue, which keeps
    /// only the largest `capacity` of them.
    fn flush_single_sample_multiplicity(&mut self) {
        self.single_sample.push(Reverse(self.current_single_sample));
        if self.single_sample.len() == self.capacity + 1 {
            self.single_sample.pop();
        }
        self.current_single_sample = 0;
    }

    /// `getPruningMultiplicity()`: `singleSampleMultiplicities.peek()`, the smallest kept.
    pub fn pruning_multiplicity(&self) -> i32 {
        self.single_sample.peek().map(|value| value.0).unwrap_or(0)
    }
}

/// `SequenceForKmers`.
#[derive(Debug, Clone)]
struct SequenceForKmers {
    sequence: Vec<u8>,
    start: usize,
    stop: usize,
    count: i32,
    is_ref: bool,
}

/// `ReadThreadingGraph`.
#[derive(Debug, Clone)]
pub struct ReadThreadingGraph {
    kmer_size: usize,
    min_base_quality: u8,
    pruning_samples: usize,
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    outgoing: Vec<Vec<usize>>,
    incoming: Vec<Vec<usize>>,
    /// `pending`, a `LinkedHashMap` from sample to its sequences.
    pending: Vec<(String, Vec<SequenceForKmers>)>,
    /// `kmerToVertexMap`, a `LinkedHashMap`: the entries in insertion order, and an index into them.
    kmer_entries: Vec<(Vec<u8>, usize)>,
    kmer_index: HashMap<Vec<u8>, usize>,
    non_unique: HashSet<Vec<u8>>,
    reference_path: Option<Vec<usize>>,
    ref_source: Option<Vec<u8>>,
    start_only_at_existing_vertex: bool,
    increase_counts_through_branches: bool,
    already_built: bool,
}

impl ReadThreadingGraph {
    /// `ReadThreadingGraph(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly,
    /// numPruningSamples, numDanglingMatchingPrefixBases)`, without the debugging switch and the
    /// dangling-end argument, which this stage does not read.
    pub fn new(
        kmer_size: usize,
        min_base_quality: u8,
        pruning_samples: usize,
    ) -> ReadThreadingGraph {
        ReadThreadingGraph {
            kmer_size,
            min_base_quality,
            pruning_samples,
            vertices: Vec::new(),
            edges: Vec::new(),
            outgoing: Vec::new(),
            incoming: Vec::new(),
            pending: Vec::new(),
            kmer_entries: Vec::new(),
            kmer_index: HashMap::new(),
            non_unique: HashSet::new(),
            reference_path: None,
            ref_source: None,
            start_only_at_existing_vertex: false,
            increase_counts_through_branches: false,
            already_built: false,
        }
    }

    /// `setThreadingStartOnlyAtExistingVertex`.
    pub fn set_threading_start_only_at_existing_vertex(&mut self, value: bool) {
        self.start_only_at_existing_vertex = value;
    }

    /// `setIncreaseCountsThroughBranches`.
    pub fn set_increase_counts_through_branches(&mut self, value: bool) {
        self.increase_counts_through_branches = value;
    }

    pub fn kmer_size(&self) -> usize {
        self.kmer_size
    }

    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }

    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    pub fn outgoing_edges(&self, vertex: usize) -> &[usize] {
        &self.outgoing[vertex]
    }

    pub fn incoming_edges(&self, vertex: usize) -> &[usize] {
        &self.incoming[vertex]
    }

    /// `referencePath`, once the reference has been threaded.
    pub fn reference_path(&self) -> Option<&[usize]> {
        self.reference_path.as_deref()
    }

    /// `kmerToVertexMap`, in its order.
    pub fn kmer_to_vertex(&self) -> &[(Vec<u8>, usize)] {
        &self.kmer_entries
    }

    /// `getNonUniqueKmers()`.
    pub fn non_unique_kmers(&self) -> &HashSet<Vec<u8>> {
        &self.non_unique
    }

    /// `addSequence(seqName, sequence, isRef)`: the whole sequence, count one, the anonymous sample.
    pub fn add_sequence(&mut self, sequence: &[u8], is_ref: bool) -> Result<(), GraphError> {
        self.add_sequence_for_sample(ANONYMOUS_SAMPLE, sequence, 0, sequence.len(), 1, is_ref)
    }

    /// `addSequence(seqName, sampleName, sequence, start, stop, count, isRef)`.
    pub fn add_sequence_for_sample(
        &mut self,
        sample: &str,
        sequence: &[u8],
        start: usize,
        stop: usize,
        count: i32,
        is_ref: bool,
    ) -> Result<(), GraphError> {
        if self.already_built {
            return Err(GraphError::AlreadyBuilt);
        }
        if stop < start {
            return Err(GraphError::InvalidSequence(format!("Invalid stop {stop}")));
        }
        if count <= 0 {
            return Err(GraphError::InvalidSequence(format!(
                "Invalid count {count}"
            )));
        }
        let entry = SequenceForKmers {
            sequence: sequence.to_vec(),
            start,
            stop,
            count,
            is_ref,
        };
        match self.pending.iter_mut().find(|(name, _)| name == sample) {
            Some((_, list)) => list.push(entry),
            None => self.pending.push((sample.to_string(), vec![entry])),
        }
        Ok(())
    }

    /// `addRead(read, header)`: each maximal run of usable bases at least `k` long, as its own
    /// sequence of the read's sample, the whole read's bases carried along.
    pub fn add_read(
        &mut self,
        sample: &str,
        bases: &[u8],
        qualities: &[u8],
    ) -> Result<(), GraphError> {
        let mut last_good: Option<usize> = None;
        for end in 0..=bases.len() {
            if end == bases.len() || !self.base_is_usable_for_assembly(bases[end], qualities[end]) {
                if let Some(start) = last_good {
                    if end - start >= self.kmer_size {
                        self.add_sequence_for_sample(sample, bases, start, end, 1, false)?;
                    }
                }
                last_good = None;
            } else if last_good.is_none() {
                last_good = Some(end);
            }
        }
        Ok(())
    }

    /// `baseIsUsableForAssembly`: not an `N` (upper case only, as `BaseUtils.Base.N.base` is), and
    /// at the minimum quality.
    fn base_is_usable_for_assembly(&self, base: u8, quality: u8) -> bool {
        base != b'N' && quality >= self.min_base_quality
    }

    /// `buildGraphIfNecessary`.
    pub fn build_graph_if_necessary(&mut self) -> Result<(), GraphError> {
        if self.already_built {
            return Ok(());
        }
        self.non_unique = self.determine_non_uniques();

        let pending = std::mem::take(&mut self.pending);
        for (_, sequences) in &pending {
            for sequence in sequences {
                self.thread_sequence(sequence)?;
            }
            for edge in &mut self.edges {
                edge.flush_single_sample_multiplicity();
            }
        }
        // `shouldRemoveReadsAfterGraphConstruction()` is true: the pending pile is cleared.
        self.already_built = true;
        for (_, vertex) in &self.kmer_entries {
            self.vertices[*vertex].additional_info.push('+');
        }
        Ok(())
    }

    /// `determineNonUniques`: every k-mer seen twice within any one pending sequence.
    fn determine_non_uniques(&self) -> HashSet<Vec<u8>> {
        let k = self.kmer_size;
        let mut non_unique = HashSet::new();
        for (_, sequences) in &self.pending {
            for sequence in sequences {
                // `determineNonUniqueKmers` starts at 0, not at `start`.
                let mut seen = HashSet::new();
                if sequence.stop >= k {
                    for i in 0..=sequence.stop - k {
                        let kmer = &sequence.sequence[i..i + k];
                        if !seen.insert(kmer) {
                            non_unique.insert(kmer.to_vec());
                        }
                    }
                }
            }
        }
        non_unique
    }

    /// `findStart`: 0 for the reference, else the first position in `start .. stop - k` whose
    /// k-mer may start a thread.
    fn find_start(&self, sequence: &SequenceForKmers) -> Option<usize> {
        if sequence.is_ref {
            return Some(0);
        }
        let k = self.kmer_size;
        let end = sequence.stop.checked_sub(k)?;
        (sequence.start..end).find(|&i| self.is_threading_start(&sequence.sequence[i..i + k]))
    }

    /// `isThreadingStart`.
    fn is_threading_start(&self, kmer: &[u8]) -> bool {
        if self.start_only_at_existing_vertex {
            self.kmer_index.contains_key(kmer)
        } else {
            !self.non_unique.contains(kmer)
        }
    }

    /// `threadSequence`.
    fn thread_sequence(&mut self, sequence: &SequenceForKmers) -> Result<(), GraphError> {
        let Some(start) = self.find_start(sequence) else {
            return Ok(());
        };
        let k = self.kmer_size;
        let starting = self.get_or_create_kmer_vertex(&sequence.sequence[start..start + k]);

        // `INCREASE_COUNTS_BACKWARDS`.
        let original = self.vertices[starting].sequence.clone();
        if k >= 2 {
            self.increase_counts_in_matched_kmers(
                sequence.count,
                starting,
                &original,
                k as isize - 2,
            );
        }

        if sequence.is_ref {
            if self.ref_source.is_some() {
                return Err(GraphError::TwoRefSources);
            }
            self.reference_path = Some(vec![starting]);
            self.ref_source = Some(sequence.sequence[sequence.start..sequence.start + k].to_vec());
        }

        let mut vertex = starting;
        if sequence.stop >= k {
            for i in start + 1..=sequence.stop - k {
                vertex = self.extend_chain_by_one(
                    vertex,
                    &sequence.sequence,
                    i,
                    sequence.count,
                    sequence.is_ref,
                )?;
                if sequence.is_ref {
                    if let Some(path) = self.reference_path.as_mut() {
                        path.push(vertex);
                    }
                }
            }
        }
        Ok(())
    }

    /// `increaseCountsInMatchedKmers`: walk back along the incoming edges whose source ends in the
    /// base before, adding the sequence's count, while the vertex has one way in (or always, if
    /// counts may go through branches).
    fn increase_counts_in_matched_kmers(
        &mut self,
        count: i32,
        vertex: usize,
        original: &[u8],
        offset: isize,
    ) {
        if offset == -1 {
            return;
        }
        let incoming = self.incoming[vertex].clone();
        for edge in incoming {
            let previous = self.edges[edge].source;
            let suffix = self.vertices[previous].suffix();
            if suffix == original[offset as usize]
                && (self.increase_counts_through_branches || self.incoming[vertex].len() == 1)
            {
                self.edges[edge].inc_multiplicity(count);
                self.increase_counts_in_matched_kmers(count, previous, original, offset - 1);
            }
        }
    }

    /// `getOrCreateKmerVertex`.
    fn get_or_create_kmer_vertex(&mut self, kmer: &[u8]) -> usize {
        match self.kmer_vertex(kmer, true) {
            Some(vertex) => vertex,
            None => self.create_vertex(kmer),
        }
    }

    /// `getKmerVertex(kmer, allowRefSource)`.
    fn kmer_vertex(&self, kmer: &[u8], allow_ref_source: bool) -> Option<usize> {
        if !allow_ref_source && self.ref_source.as_deref() == Some(kmer) {
            return None;
        }
        self.kmer_index
            .get(kmer)
            .map(|&entry| self.kmer_entries[entry].1)
    }

    /// `createVertex`: a new vertex, tracked unless its k-mer is non-unique or already tracked.
    fn create_vertex(&mut self, kmer: &[u8]) -> usize {
        let vertex = self.vertices.len();
        self.vertices.push(Vertex {
            sequence: kmer.to_vec(),
            additional_info: String::new(),
        });
        self.outgoing.push(Vec::new());
        self.incoming.push(Vec::new());
        if !self.non_unique.contains(kmer) && !self.kmer_index.contains_key(kmer) {
            self.kmer_index
                .insert(kmer.to_vec(), self.kmer_entries.len());
            self.kmer_entries.push((kmer.to_vec(), vertex));
        }
        vertex
    }

    /// `extendChainByOne`.
    fn extend_chain_by_one(
        &mut self,
        previous: usize,
        sequence: &[u8],
        kmer_start: usize,
        count: i32,
        is_ref: bool,
    ) -> Result<usize, GraphError> {
        let next_base = sequence[kmer_start + self.kmer_size - 1];
        for position in 0..self.outgoing[previous].len() {
            let edge = self.outgoing[previous][position];
            let target = self.edges[edge].target;
            if self.vertices[target].suffix() == next_base {
                self.edges[edge].inc_multiplicity(count);
                return Ok(target);
            }
        }

        let kmer = &sequence[kmer_start..kmer_start + self.kmer_size];
        // `getNextKmerVertexForChainExtension`.
        let merge = self.kmer_vertex(kmer, false);
        if is_ref && merge.is_some() {
            return Err(GraphError::UniqueVertexOnReference);
        }
        let next = match merge {
            Some(vertex) => vertex,
            None => self.create_vertex(kmer),
        };
        self.add_edge(previous, next, is_ref, count);
        Ok(next)
    }

    fn add_edge(&mut self, source: usize, target: usize, is_ref: bool, multiplicity: i32) {
        let edge = self.edges.len();
        self.edges.push(Edge::new(
            source,
            target,
            is_ref,
            multiplicity,
            self.pruning_samples,
        ));
        self.outgoing[source].push(edge);
        self.incoming[target].push(edge);
    }

    /// `isLowQualityGraph()`: more than a quarter as many non-unique k-mers as tracked ones.
    pub fn is_low_quality_graph(&self) -> bool {
        self.non_unique.len() * 4 > self.kmer_entries.len()
    }

    /// `isRefSource(v)`: no incoming reference edge and an outgoing one, or the only vertex.
    pub fn is_ref_source(&self, vertex: usize) -> bool {
        if self.incoming[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return false;
        }
        if self.outgoing[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return true;
        }
        self.vertices.len() == 1
    }

    /// `isRefSink(v)`.
    pub fn is_ref_sink(&self, vertex: usize) -> bool {
        if self.outgoing[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return false;
        }
        if self.incoming[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return true;
        }
        self.vertices.len() == 1
    }

    /// `getReferenceSourceVertex()`: the first vertex, in vertex order, that is a reference source.
    pub fn reference_source_vertex(&self) -> Option<usize> {
        (0..self.vertices.len()).find(|&v| self.is_ref_source(v))
    }

    /// `getReferenceSinkVertex()`.
    pub fn reference_sink_vertex(&self) -> Option<usize> {
        (0..self.vertices.len()).find(|&v| self.is_ref_sink(v))
    }

    /// `hasCycles()`: JGraphT's `CycleDetector.detectCycles()`, a self-loop included. Only the
    /// answer is observable, so any traversal will do.
    pub fn has_cycles(&self) -> bool {
        // 0 unvisited, 1 on the current path, 2 finished.
        let mut state = vec![0u8; self.vertices.len()];
        for root in 0..self.vertices.len() {
            if state[root] != 0 {
                continue;
            }
            let mut stack: Vec<(usize, usize)> = vec![(root, 0)];
            state[root] = 1;
            while let Some(&mut (vertex, ref mut next)) = stack.last_mut() {
                if *next < self.outgoing[vertex].len() {
                    let target = self.edges[self.outgoing[vertex][*next]].target;
                    *next += 1;
                    match state[target] {
                        1 => return true,
                        0 => {
                            state[target] = 1;
                            stack.push((target, 0));
                        }
                        _ => {}
                    }
                } else {
                    state[vertex] = 2;
                    stack.pop();
                }
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_reference_alone_is_one_chain_of_reference_edges() {
        let mut graph = ReadThreadingGraph::new(4, 10, 1);
        graph.add_sequence(b"ACGTTGCA", true).unwrap();
        graph.build_graph_if_necessary().unwrap();
        assert_eq!(graph.vertices().len(), 5);
        assert_eq!(graph.edges().len(), 4);
        assert!(graph
            .edges()
            .iter()
            .all(|e| e.is_ref && e.multiplicity == 1));
        assert_eq!(graph.reference_path().unwrap(), &[0, 1, 2, 3, 4]);
        assert_eq!(graph.reference_source_vertex(), Some(0));
        assert_eq!(graph.reference_sink_vertex(), Some(4));
        assert!(!graph.has_cycles());
    }

    #[test]
    fn a_repeated_kmer_is_never_merged_into() {
        // ACGT occurs twice in the reference: two vertices, neither in the k-mer map.
        let mut graph = ReadThreadingGraph::new(4, 10, 1);
        graph.add_sequence(b"ACGTACGTT", true).unwrap();
        graph.build_graph_if_necessary().unwrap();
        let acgt = graph
            .vertices()
            .iter()
            .filter(|v| v.sequence == b"ACGT")
            .count();
        assert_eq!(acgt, 2);
        assert!(graph.kmer_to_vertex().iter().all(|(k, _)| k != b"ACGT"));
    }

    #[test]
    fn the_pruning_multiplicity_keeps_the_largest_per_sample_counts() {
        let mut edge = Edge::new(0, 1, false, 3, 2);
        edge.flush_single_sample_multiplicity();
        edge.inc_multiplicity(1);
        edge.flush_single_sample_multiplicity();
        // Queue was [3, 3], then 1 joined and the smallest left: [3, 3].
        assert_eq!(edge.pruning_multiplicity(), 3);
        assert_eq!(edge.multiplicity, 4);
    }
}
