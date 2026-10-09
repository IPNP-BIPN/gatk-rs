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
    /// `IllegalStateException("Graph must have ref source and sink vertices")`.
    NoReferenceSourceOrSink,
    /// `IllegalStateException("Should have eliminated all but the reference sink, ...")`.
    MoreThanOneSink,
    /// `IllegalStateException("Should have eliminated all but the reference source, ...")`.
    MoreThanOneSource,
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
    /// Whether each vertex and edge is still in the graph. A removal keeps the others' order, as
    /// removing from a `LinkedHashMap` or an array-backed set does.
    vertex_alive: Vec<bool>,
    edge_alive: Vec<bool>,
    outgoing: Vec<Vec<usize>>,
    incoming: Vec<Vec<usize>>,
    /// `pending`, a `LinkedHashMap` from sample to its sequences.
    pending: Vec<(String, Vec<SequenceForKmers>)>,
    /// `kmerToVertexMap`, a `LinkedHashMap`: the entries in insertion order (`None` once removed),
    /// and an index into them.
    kmer_entries: Vec<Option<(Vec<u8>, usize)>>,
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
            vertex_alive: Vec::new(),
            edge_alive: Vec::new(),
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

    /// `vertexSet()`, in its order: the ids of the vertices still in the graph.
    pub fn vertex_ids(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.vertices.len()).filter(|&v| self.vertex_alive[v])
    }

    /// `edgeSet()`, in its order.
    pub fn edge_ids(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.edges.len()).filter(|&e| self.edge_alive[e])
    }

    pub fn vertex(&self, id: usize) -> &Vertex {
        &self.vertices[id]
    }

    pub fn edge(&self, id: usize) -> &Edge {
        &self.edges[id]
    }

    pub fn vertex_count(&self) -> usize {
        self.vertex_ids().count()
    }

    pub fn edge_count(&self) -> usize {
        self.edge_ids().count()
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
    pub fn kmer_to_vertex(&self) -> impl Iterator<Item = &(Vec<u8>, usize)> + '_ {
        self.kmer_entries.iter().flatten()
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
        for (_, vertex) in self.kmer_entries.iter().flatten() {
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
            .and_then(|&entry| self.kmer_entries[entry].as_ref())
            .map(|(_, vertex)| *vertex)
    }

    /// `createVertex`: a new vertex, tracked unless its k-mer is non-unique or already tracked.
    fn create_vertex(&mut self, kmer: &[u8]) -> usize {
        let vertex = self.vertices.len();
        self.vertices.push(Vertex {
            sequence: kmer.to_vec(),
            additional_info: String::new(),
        });
        self.vertex_alive.push(true);
        self.outgoing.push(Vec::new());
        self.incoming.push(Vec::new());
        if !self.non_unique.contains(kmer) && !self.kmer_index.contains_key(kmer) {
            self.kmer_index
                .insert(kmer.to_vec(), self.kmer_entries.len());
            self.kmer_entries.push(Some((kmer.to_vec(), vertex)));
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
        self.edge_alive.push(true);
        self.outgoing[source].push(edge);
        self.incoming[target].push(edge);
    }

    /// `isLowQualityGraph()`: more than a quarter as many non-unique k-mers as tracked ones.
    pub fn is_low_quality_graph(&self) -> bool {
        self.non_unique.len() * 4 > self.kmer_to_vertex().count()
    }

    /// `isRefSource(v)`: no incoming reference edge and an outgoing one, or the only vertex.
    pub fn is_ref_source(&self, vertex: usize) -> bool {
        if self.incoming[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return false;
        }
        if self.outgoing[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return true;
        }
        self.vertex_count() == 1
    }

    /// `isRefSink(v)`.
    pub fn is_ref_sink(&self, vertex: usize) -> bool {
        if self.outgoing[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return false;
        }
        if self.incoming[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return true;
        }
        self.vertex_count() == 1
    }

    /// `getReferenceSourceVertex()`: the first vertex, in vertex order, that is a reference source.
    pub fn reference_source_vertex(&self) -> Option<usize> {
        self.vertex_ids().find(|&v| self.is_ref_source(v))
    }

    /// `getReferenceSinkVertex()`.
    pub fn reference_sink_vertex(&self) -> Option<usize> {
        self.vertex_ids().find(|&v| self.is_ref_sink(v))
    }

    /// `hasCycles()`: JGraphT's `CycleDetector.detectCycles()`, a self-loop included. Only the
    /// answer is observable, so any traversal will do.
    pub fn has_cycles(&self) -> bool {
        // 0 unvisited, 1 on the current path, 2 finished.
        let mut state = vec![0u8; self.vertices.len()];
        for root in self.vertex_ids() {
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

    /// `inDegreeOf(v)`.
    pub fn in_degree(&self, vertex: usize) -> usize {
        self.incoming[vertex].len()
    }

    /// `outDegreeOf(v)`.
    pub fn out_degree(&self, vertex: usize) -> usize {
        self.outgoing[vertex].len()
    }

    /// `getSources()`: the vertices with no incoming edge, in vertex order.
    pub fn sources(&self) -> Vec<usize> {
        self.vertex_ids()
            .filter(|&v| self.in_degree(v) == 0)
            .collect()
    }

    /// `getSinks()`.
    pub fn sinks(&self) -> Vec<usize> {
        self.vertex_ids()
            .filter(|&v| self.out_degree(v) == 0)
            .collect()
    }

    /// `removeEdge(e)`: out of the edge set and both endpoints' lists, the rest in order.
    pub fn remove_edge(&mut self, edge: usize) {
        if !self.edge_alive[edge] {
            return;
        }
        self.edge_alive[edge] = false;
        let Edge { source, target, .. } = self.edges[edge];
        self.outgoing[source].retain(|&e| e != edge);
        self.incoming[target].retain(|&e| e != edge);
    }

    /// `AbstractReadThreadingGraph.removeVertex`: JGraphT's (its edges first, then the vertex),
    /// and the vertex's sequence out of the k-mer map, whichever vertex the entry pointed to.
    pub fn remove_vertex(&mut self, vertex: usize) {
        if !self.vertex_alive[vertex] {
            return;
        }
        let touching: Vec<usize> = self.outgoing[vertex]
            .iter()
            .chain(self.incoming[vertex].iter())
            .copied()
            .collect();
        for edge in touching {
            self.remove_edge(edge);
        }
        self.vertex_alive[vertex] = false;
        if let Some(entry) = self.kmer_index.remove(&self.vertices[vertex].sequence) {
            self.kmer_entries[entry] = None;
        }
    }

    /// `AbstractReadThreadingGraph.removeSingletonOrphanVertices`: every vertex with no edge at
    /// all. Unlike `BaseGraph`'s, it does not spare an isolated reference source.
    pub fn remove_singleton_orphan_vertices(&mut self) {
        let orphans: Vec<usize> = self
            .vertex_ids()
            .filter(|&v| self.in_degree(v) == 0 && self.out_degree(v) == 0)
            .collect();
        for vertex in orphans {
            self.remove_vertex(vertex);
        }
    }

    /// `removePathsNotConnectedToRef`: keep only the vertices reachable forward from the reference
    /// source and backward from the reference sink. The reference removes them from a `HashSet`,
    /// in hash order, but the result of removing a set does not depend on the order.
    pub fn remove_paths_not_connected_to_ref(&mut self) -> Result<(), GraphError> {
        let (Some(source), Some(sink)) =
            (self.reference_source_vertex(), self.reference_sink_vertex())
        else {
            return Err(GraphError::NoReferenceSourceOrSink);
        };
        let forward = self.reachable(source, true);
        let backward = self.reachable(sink, false);
        let doomed: Vec<usize> = self
            .vertex_ids()
            .filter(|&v| !(forward[v] && backward[v]))
            .collect();
        for vertex in doomed {
            self.remove_vertex(vertex);
        }
        if self.sinks().len() > 1 {
            return Err(GraphError::MoreThanOneSink);
        }
        if self.sources().len() > 1 {
            return Err(GraphError::MoreThanOneSource);
        }
        Ok(())
    }

    /// The vertices reachable from `start` along outgoing edges (or incoming ones), `start`
    /// included: `BaseGraphIterator`, of which only the membership is read.
    fn reachable(&self, start: usize, forward: bool) -> Vec<bool> {
        let mut seen = vec![false; self.vertices.len()];
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(vertex) = stack.pop() {
            let edges = if forward {
                &self.outgoing[vertex]
            } else {
                &self.incoming[vertex]
            };
            for &edge in edges {
                let next = if forward {
                    self.edges[edge].target
                } else {
                    self.edges[edge].source
                };
                if !seen[next] {
                    seen[next] = true;
                    stack.push(next);
                }
            }
        }
        seen
    }
}

/// A `Path` as `ChainPruner.findChain` builds it: its edges in order, and the vertex it ends on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chain {
    pub edges: Vec<usize>,
    pub last_vertex: usize,
}

impl Chain {
    /// `getVertices()`: the first edge's source, then every edge's target.
    pub fn vertices(&self, graph: &ReadThreadingGraph) -> Vec<usize> {
        let mut out = vec![graph.edge(self.edges[0]).source];
        out.extend(self.edges.iter().map(|&e| graph.edge(e).target));
        out
    }
}

/// `ChainPruner.findAllChains`: from every source, and then from the end of every chain found,
/// the maximal linear chain behind each outgoing edge. A chain stops at a vertex with other than
/// one way out, more than one way in, or that is the chain's own start.
pub fn find_all_chains(graph: &ReadThreadingGraph) -> Vec<Chain> {
    let mut starts: std::collections::VecDeque<usize> = graph.sources().into();
    let mut seen: HashSet<usize> = starts.iter().copied().collect();
    let mut chains = Vec::new();
    while let Some(start) = starts.pop_front() {
        for &edge in graph.outgoing_edges(start) {
            let chain = find_chain(graph, edge);
            if seen.insert(chain.last_vertex) {
                starts.push_back(chain.last_vertex);
            }
            chains.push(chain);
        }
    }
    chains
}

fn find_chain(graph: &ReadThreadingGraph, start: usize) -> Chain {
    let mut edges = vec![start];
    let first = graph.edge(start).source;
    let mut last = graph.edge(start).target;
    loop {
        if graph.out_degree(last) != 1 || graph.in_degree(last) > 1 || last == first {
            break;
        }
        let next = graph.outgoing_edges(last)[0];
        edges.push(next);
        last = graph.edge(next).target;
    }
    Chain {
        edges,
        last_vertex: last,
    }
}

/// `LowWeightChainPruner.pruneLowWeightChains`: remove every chain whose edges are all
/// non-reference with a pruning multiplicity under `prune_factor`, then the orphans.
pub fn prune_low_weight_chains(graph: &mut ReadThreadingGraph, prune_factor: i32) {
    let chains = find_all_chains(graph);
    let doomed: Vec<Chain> = chains
        .into_iter()
        .filter(|chain| {
            chain.edges.iter().all(|&e| {
                let edge = graph.edge(e);
                edge.pruning_multiplicity() < prune_factor && !edge.is_ref
            })
        })
        .collect();
    for chain in doomed {
        for edge in chain.edges {
            graph.remove_edge(edge);
        }
    }
    graph.remove_singleton_orphan_vertices();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_reference_alone_is_one_chain_of_reference_edges() {
        let mut graph = ReadThreadingGraph::new(4, 10, 1);
        graph.add_sequence(b"ACGTTGCA", true).unwrap();
        graph.build_graph_if_necessary().unwrap();
        assert_eq!(graph.vertex_count(), 5);
        assert_eq!(graph.edge_count(), 4);
        assert!(graph
            .edge_ids()
            .all(|e| graph.edge(e).is_ref && graph.edge(e).multiplicity == 1));
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
            .vertex_ids()
            .filter(|&v| graph.vertex(v).sequence == b"ACGT")
            .count();
        assert_eq!(acgt, 2);
        assert!(graph.kmer_to_vertex().all(|(k, _)| k != b"ACGT"));
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

    #[test]
    fn a_weight_one_bubble_is_pruned_and_the_reference_kept() {
        let reference = b"ACGTTGCATGTCGCATGATGCATGAGAG";
        let mut read = reference.to_vec();
        read[14] = b'T';
        let mut graph = ReadThreadingGraph::new(10, 10, 1);
        graph.add_sequence(reference, true).unwrap();
        graph.add_read("s1", &read, &vec![30; read.len()]).unwrap();
        graph.build_graph_if_necessary().unwrap();
        assert!(graph.vertex_count() > reference.len() - 9);
        prune_low_weight_chains(&mut graph, 2);
        assert_eq!(graph.vertex_count(), reference.len() - 9);
        assert!(graph.edge_ids().all(|e| graph.edge(e).is_ref));
    }
}
