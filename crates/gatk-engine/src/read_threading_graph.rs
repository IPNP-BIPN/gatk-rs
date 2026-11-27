//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading`
//! (`AbstractReadThreadingGraph`, `ReadThreadingGraph`, `MultiDeBruijnVertex`) and
//! `haplotypecaller.graphs` (`BaseGraph`, `BaseEdge`, `MultiSampleEdge`, `ChainPruner`,
//! `LowWeightChainPruner`, `AdaptiveChainPruner`), GATK 4.7.0.0; every one of those `.java` files
//! is unchanged since 4.6.2.0. It holds the read threading graph as `buildGraphIfNecessary` leaves
//! it, the chain pruners, and (in `dangling_ends`) the recovery of dangling tails and heads. The
//! conversion to a sequence graph is `seq_graph`.
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
//!
//! # Where an order reaches the output
//!
//! Every golden that reads this graph compares it line for line, vertices, edges and k-mer map
//! entries in their order (`tests/read_threading_graph.rs` against `read_threading_graph.txt.gz`,
//! `tests/chain_pruner.rs` against `chain_pruner.txt.gz`, `tests/adaptive_chain_pruner.rs`
//! against `adaptive_chain_pruner.txt.gz`). The orders that matter:
//!
//! * vertex order is creation order: `vertexSet()` of the reference's `LinkedHashMap`. It decides
//!   which vertex `getReferenceSourceVertex` finds first and the order `toSequenceGraph` copies;
//! * each vertex's edge lists are in creation order and keep it across removals, which is what
//!   `extendChainByOne`'s "first matching edge" and `findAllChains`' order read;
//! * `kmerToVertexMap` is a `LinkedHashMap`, and the goldens' `kmer` rows list it in that order.
//!   Which vertices it holds also shows in their `vertex` rows and in the dot files, because
//!   `buildGraphIfNecessary` appends `+` to exactly those vertices' debugging text.
//!
//! # Reading this file without reading Rust
//!
//! - A vertex or an edge is a `usize` (an unsigned index type Java lacks) into `vertices` or
//!   `edges`. Nothing is ever freed: a removal clears a flag in `vertex_alive` or `edge_alive`, so
//!   an index stays a stable identity, the way a Java object reference does.
//! - `Option<T>` is "a value or nothing", Java's nullable reference checked by the compiler;
//!   `Some(x)` is the value and `None` the `null`. `Result<T, GraphError>` is "a value or one of
//!   the refusals below": `Err(..)` is the throw, `Ok(..)` the return, and `?` after a call passes
//!   an `Err` (or, in a function answering `Option`, a `None`) straight up to the caller, as an
//!   uncaught exception would.
//! - `let Some(x) = y else { ... };` unwraps `y` or takes the early exit in the braces.
//! - `&T` is a read-only borrowed reference, `&mut T` a writable one, and while a `&mut` exists
//!   nothing else may touch the value. That is why some loops first `.clone()` a list they read
//!   while they change the graph: the copy is for the compiler, not for the algorithm.
//! - `impl Iterator<Item = usize> + '_` is a lazy stream (Java's `Stream<Integer>`) that borrows
//!   the graph; `.filter`, `.map`, `.find`, `.any`, `.all` and `.collect()` are the stream
//!   operations of the same names.
//! - `|a, b| ...` is a lambda. A function's last expression, written without `return`, is its
//!   value.
//! - `x as i8` reinterprets or truncates a number the way a Java cast does; `isize` and `i32` are
//!   signed, `usize` is not, so a Java index that can be -1 becomes `isize` or `Option<usize>`.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};

mod dangling_ends;
pub mod junction_trees;

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
    /// `IllegalStateException("Attempting to add sequence to a graph that has already been
    /// built")`.
    AlreadyBuilt,
    /// `SequenceForKmers`'s argument checks.
    InvalidSequence(String),
    /// `IllegalStateException("Graph must have ref source and sink vertices")`.
    NoReferenceSourceOrSink,
    /// `IllegalStateException("Should have eliminated all but the reference sink, ...")`.
    MoreThanOneSink,
    /// `IllegalStateException("Should have eliminated all but the reference source, ...")`.
    MoreThanOneSource,
    /// `Utils.validateArg` in dangling-end recovery: a negative prune factor or branch length, or a
    /// CIGAR whose first or last element is not the M the merge expects.
    IllegalArgument(&'static str),
    /// `IllegalStateException`: recovery before the graph is built, or a dangling end that is no
    /// longer dangling when its turn comes.
    IllegalState(&'static str),
    /// `IndexOutOfBoundsException` or `ArrayIndexOutOfBoundsException`: an index the reference
    /// computes and then reads without a check.
    IndexOutOfBounds,
    /// `NoSuchElementException`: the heaviest edge of a vertex that has none.
    NoSuchElement,
    /// The Smith-Waterman aligner's own refusal.
    Alignment(crate::alignment_utils::AlignmentError),
}

/// `MultiDeBruijnVertex`: a k-mer's bases and the debugging text `buildGraphIfNecessary` appends
/// to.
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
    /// `currentSingleSampleMultiplicity`: what the sample being threaded has added so far.
    current_single_sample: i32,
    /// `PriorityQueue<Integer>`, a min-heap: the smallest is what `poll` drops and `peek` reads.
    /// Rust's `BinaryHeap` is a max-heap; wrapping each value in `Reverse` flips its ordering, so
    /// the largest `Reverse` is the smallest count, as in Java's queue.
    single_sample: BinaryHeap<Reverse<i32>>,
    /// `singleSampleCapacity`, the graph's `numPruningSamples`.
    capacity: usize,
}

impl Edge {
    /// `new MultiSampleEdge(isRef, multiplicity, singleSampleCapacity)`. The queue starts with the
    /// initial multiplicity already in it, so that entry competes with the per-sample counts:
    ///
    /// ```java
    /// singleSampleMultiplicities = new PriorityQueue<>(singleSampleCapacity);
    /// singleSampleMultiplicities.add(multiplicity);
    /// ```
    ///
    /// (`MultiSampleEdge.java:48-49`.) The unit test
    /// `the_pruning_multiplicity_keeps_the_largest_per_sample_counts` pins that.
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
    ///
    /// The queue always holds at least the constructor's entry, so the `unwrap_or(0)` (the value
    /// taken when the `Option` from `peek` is `None`) is never reached.
    pub fn pruning_multiplicity(&self) -> i32 {
        self.single_sample.peek().map(|value| value.0).unwrap_or(0)
    }
}

/// `SequenceForKmers`: the whole of a read's (or the reference's) bases, and the run of them to
/// thread. The name the reference also keeps is only for debugging output, so it is not kept.
#[derive(Debug, Clone)]
struct SequenceForKmers {
    sequence: Vec<u8>,
    /// First base to thread, inclusive.
    start: usize,
    /// Last base to thread, exclusive.
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
    /// Per vertex, `outgoingEdgesOf(v)` in the reference's order: creation order, kept across
    /// removals.
    outgoing: Vec<Vec<usize>>,
    /// Per vertex, `incomingEdgesOf(v)`, likewise.
    incoming: Vec<Vec<usize>>,
    /// `pending`, a `LinkedHashMap` from sample to its sequences. A `Vec` of pairs searched
    /// linearly keeps the samples in the order first seen, which is the order they are threaded.
    pending: Vec<(String, Vec<SequenceForKmers>)>,
    /// `kmerToVertexMap`, a `LinkedHashMap`: the entries in insertion order (`None` once removed),
    /// and an index into them.
    kmer_entries: Vec<Option<(Vec<u8>, usize)>>,
    /// The hash half of `kmerToVertexMap`: k-mer to its position in `kmer_entries`. Only looked
    /// up, never iterated, so its own order is not observable.
    kmer_index: HashMap<Vec<u8>, usize>,
    /// `nonUniqueKmers`. Only asked `contains`, and stays empty for the linked graph.
    non_unique: HashSet<Vec<u8>>,
    /// `referencePath`: the reference's vertices in threading order, `None` until it is threaded.
    reference_path: Option<Vec<usize>>,
    /// `refSource`: the reference's first k-mer, which `getKmerVertex(kmer, false)` refuses.
    ref_source: Option<Vec<u8>>,
    start_only_at_existing_vertex: bool,
    increase_counts_through_branches: bool,
    already_built: bool,
    /// `minMatchingBasesToDanglingEndRecovery`: -1 takes the legacy dangling-head merge.
    min_matching_bases: i32,
    /// `JunctionTreeLinkedDeBruijnGraph` rather than `ReadThreadingGraph`: no non-unique k-mers, a
    /// k-mer merges wherever it is already tracked, the pending sequences are kept for threading
    /// junction trees, and the reference ends are the first and last vertex of the reference path.
    linked: bool,
}

impl ReadThreadingGraph {
    /// `ReadThreadingGraph(kmerSize, debugGraphTransformations, minBaseQualityToUseInAssembly,
    /// numPruningSamples, numDanglingMatchingPrefixBases)`, without the debugging switch. The
    /// dangling-end argument starts at -1 and is set with
    /// `set_min_matching_bases_to_dangling_end_recovery`.
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
            min_matching_bases: -1,
            linked: false,
        }
    }

    /// `JunctionTreeLinkedDeBruijnGraph(kmerSize, debugGraphTransformations,
    /// minBaseQualityToUseInAssembly, numPruningSamples, numDanglingMatchingPrefixBases)`, without
    /// the debugging switch.
    pub fn new_linked(
        kmer_size: usize,
        min_base_quality: u8,
        pruning_samples: usize,
    ) -> ReadThreadingGraph {
        ReadThreadingGraph {
            linked: true,
            // `..x` takes every field not named above from `x`: a plain graph with the flag set.
            ..ReadThreadingGraph::new(kmer_size, min_base_quality, pruning_samples)
        }
    }

    /// Whether this is the linked de Bruijn graph.
    pub fn is_linked(&self) -> bool {
        self.linked
    }

    /// `setThreadingStartOnlyAtExistingVertex`.
    pub fn set_threading_start_only_at_existing_vertex(&mut self, value: bool) {
        self.start_only_at_existing_vertex = value;
    }

    /// The constructor's `numDanglingMatchingPrefixBases`, or
    /// `setMinMatchingBasesToDanglingEndRecovery`.
    pub fn set_min_matching_bases_to_dangling_end_recovery(&mut self, value: i32) {
        self.min_matching_bases = value;
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

    /// `addSequence(seqName, sequence, isRef)`: the whole sequence, count one, the anonymous
    /// sample.
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
        // `pending.computeIfAbsent(sampleName, ...)`: a new sample goes to the end.
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
        // The reference's `int lastGood = -1`: `None` is the -1.
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
        // `preprocessReads()`: the linked graph tracks no non-unique k-mers.
        if !self.linked {
            self.non_unique = self.determine_non_uniques();
        }

        // `shouldRemoveReadsAfterGraphConstruction()`: the linked graph keeps them, to thread
        // junction trees.
        let pending = if self.linked {
            self.pending.clone()
        } else {
            std::mem::take(&mut self.pending)
        };
        // Sample by sample, in the order first seen. Every edge closes the sample's count once
        // the sample is done; the reference walks `edgeSet()`, the live edges, and nothing is
        // removed while building, so walking every edge here is the same set.
        for (_, sequences) in &pending {
            for sequence in sequences {
                self.thread_sequence(sequence)?;
            }
            for edge in &mut self.edges {
                edge.flush_single_sample_multiplicity();
            }
        }
        self.already_built = true;
        // `.flatten()` over the `Option` entries skips the removed (`None`) ones.
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
                // `determineNonUniqueKmers` starts at 0, not at `start` (ReadThreadingGraph.java
                // lines 161-162):
                // `final int stopPosition = seqForKmers.stop - kmerSize;`
                // `for (int i = 0; i <= stopPosition; i++) {`
                // `seen.insert` answers false for a k-mer already in the set, like `Set.add`.
                let mut seen = HashSet::new();
                // A negative `stopPosition` runs no iteration in Java; `usize` cannot go below
                // zero (the subtraction would panic), hence the guard.
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
    /// k-mer may start a thread. `None` is the reference's -1.
    ///
    /// ```java
    /// for (int i = seqForKmers.start; i < seqForKmers.stop - kmerSize; i++) {
    /// ```
    ///
    /// (`AbstractReadThreadingGraph.java:205`.) The bound is exclusive, so the last k-mer never
    /// starts a thread; and the reference starts at 0 even when its `start` is not 0.
    fn find_start(&self, sequence: &SequenceForKmers) -> Option<usize> {
        if sequence.is_ref {
            return Some(0);
        }
        let k = self.kmer_size;
        // `checked_sub` is `None` where `stop - k` would be negative, where Java's loop runs no
        // iteration and answers -1; the `?` returns that `None`.
        let end = sequence.stop.checked_sub(k)?;
        (sequence.start..end).find(|&i| self.is_threading_start(&sequence.sequence[i..i + k]))
    }

    /// `isThreadingStart`, in `ReadThreadingGraph` and in `JunctionTreeLinkedDeBruijnGraph`.
    fn is_threading_start(&self, kmer: &[u8]) -> bool {
        if self.linked {
            // `return !startThreadingOnlyAtExistingVertex || kmers.contains(kmer);`
            // (JunctionTreeLinkedDeBruijnGraph.java:101). Its set `kmers` is never filled, so a
            // thread that may start only at an existing vertex never starts.
            return !self.start_only_at_existing_vertex;
        }
        if self.start_only_at_existing_vertex {
            self.kmer_index.contains_key(kmer)
        } else {
            !self.non_unique.contains(kmer)
        }
    }

    /// `threadSequence`: find the start, then extend the chain one base at a time to the last
    /// k-mer, adding the sequence's count to every edge walked or made.
    fn thread_sequence(&mut self, sequence: &SequenceForKmers) -> Result<(), GraphError> {
        let Some(start) = self.find_start(sequence) else {
            return Ok(());
        };
        let k = self.kmer_size;
        // `&v[a..b]` is a borrowed view of elements `a` to `b - 1`: Java's `new Kmer(seq, a, k)`.
        let starting = self.get_or_create_kmer_vertex(&sequence.sequence[start..start + k]);

        // `INCREASE_COUNTS_BACKWARDS`, always true:
        // `increaseCountsInMatchedKmers(seqForKmers, startingVertex,
        // startingVertex.getSequence(), kmerSize - 2);` (AbstractReadThreadingGraph.java:280).
        // The k-mer is cloned only so the graph may be changed while it is read. For k = 1 the
        // reference's offset is -1 and the walk returns at once; the guard skips the same case.
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
            // The k-mer at `start`, not at the threading start:
            // `refSource = new Kmer(seqForKmers.sequence, seqForKmers.start, kmerSize);`
            // (AbstractReadThreadingGraph.java:294).
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
        // In the reference's order. Copied because the loop changes edges while reading the
        // list; only multiplicities change, so the copy and the live list stay the same.
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

    /// `getKmerVertex(kmer, allowRefSource)`: the tracked vertex of the k-mer, unless it is the
    /// reference's first k-mer and that is not allowed.
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
    ///
    /// That is `ReadThreadingGraph.trackKmer`; the linked graph's is
    /// `kmerToVertexMap.putIfAbsent(kmer, newVertex);` (JunctionTreeLinkedDeBruijnGraph.java:126),
    /// the same test because its non-unique set is empty. A k-mer tracked again after its entry
    /// was removed goes to the end of the map, as a `LinkedHashMap` re-insertion does.
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

    /// `extendChainByOne`: follow the FIRST outgoing edge, in creation order, whose target ends
    /// in the next base, adding the count; else merge into a tracked vertex or make a new one,
    /// and join it with a new edge carrying the count.
    fn extend_chain_by_one(
        &mut self,
        previous: usize,
        sequence: &[u8],
        kmer_start: usize,
        count: i32,
        is_ref: bool,
    ) -> Result<usize, GraphError> {
        let next_base = sequence[kmer_start + self.kmer_size - 1];
        // By position rather than by iterator: the body changes an edge, which the compiler
        // would not allow while an iterator borrows the list.
        for position in 0..self.outgoing[previous].len() {
            let edge = self.outgoing[previous][position];
            let target = self.edges[edge].target;
            if self.vertices[target].suffix() == next_base {
                self.edges[edge].inc_multiplicity(count);
                return Ok(target);
            }
        }

        let kmer = &sequence[kmer_start..kmer_start + self.kmer_size];
        // `getNextKmerVertexForChainExtension`: the linked graph merges into any tracked k-mer,
        // the reference source's included, and never refuses.
        if self.linked {
            let next = match self.kmer_vertex(kmer, true) {
                Some(vertex) => vertex,
                None => self.create_vertex(kmer),
            };
            self.add_edge(previous, next, is_ref, count);
            return Ok(next);
        }
        // `ReadThreadingGraph.getNextKmerVertexForChainExtension`: never the reference source,
        // and the reference itself may never merge into a tracked vertex.
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

    /// `addEdge(source, target, createEdge(isRef, multiplicity))`, appended to the edge set and to
    /// both endpoints' lists.
    ///
    /// Unlike JGraphT it does not check for an existing `source -> target` edge, because no
    /// caller can have one: `extend_chain_by_one` only gets here when no outgoing edge's target
    /// ends in the next base (an existing edge to the target would), and the other callers join a
    /// vertex just made. `dangling_ends` goes through `add_edge_unless_present` instead.
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

    /// `isLowQualityGraph()`: more than a quarter as many non-unique k-mers as tracked ones,
    /// `return nonUniqueKmers.size() * 4 > kmerToVertexMap.size();` (ReadThreadingGraph.java:84).
    /// The linked graph always answers false.
    pub fn is_low_quality_graph(&self) -> bool {
        !self.linked && self.non_unique.len() * 4 > self.kmer_to_vertex().count()
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
    /// Vertex order decides between candidates. The linked graph overrides it with the first
    /// vertex of the stored reference path.
    pub fn reference_source_vertex(&self) -> Option<usize> {
        if self.linked {
            return self
                .reference_path
                .as_ref()
                .and_then(|path| path.first().copied());
        }
        self.vertex_ids().find(|&v| self.is_ref_source(v))
    }

    /// `getReferenceSinkVertex()`.
    pub fn reference_sink_vertex(&self) -> Option<usize> {
        if self.linked {
            return self
                .reference_path
                .as_ref()
                .and_then(|path| path.last().copied());
        }
        self.vertex_ids().find(|&v| self.is_ref_sink(v))
    }

    /// `hasCycles()`: JGraphT's `CycleDetector.detectCycles()`, a self-loop included. Only the
    /// answer is observable, so any traversal will do.
    pub fn has_cycles(&self) -> bool {
        // A depth-first search with an explicit stack: an edge back to a vertex still on the
        // current path is a cycle. 0 unvisited, 1 on the current path, 2 finished.
        let mut state = vec![0u8; self.vertices.len()];
        for root in self.vertex_ids() {
            if state[root] != 0 {
                continue;
            }
            // Each frame is a vertex and the position of the next outgoing edge to try.
            let mut stack: Vec<(usize, usize)> = vec![(root, 0)];
            state[root] = 1;
            // `while let` loops while the top frame exists; the pattern copies its vertex and
            // binds `next` as a writable reference into the frame, so `*next += 1` advances it.
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
        // Destructuring: copy out the two named fields, ignoring the rest (`..`).
        let Edge { source, target, .. } = self.edges[edge];
        // `retain` is `removeIf` with the test inverted; the survivors keep their order.
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
        // By k-mer, not by vertex: `kmerToVertexMap.remove(kmer);`
        // (AbstractReadThreadingGraph.java:370), so removing an untracked copy of a k-mer drops
        // the tracked vertex's entry too.
        if let Some(entry) = self.kmer_index.remove(&self.vertices[vertex].sequence) {
            self.kmer_entries[entry] = None;
        }
    }

    /// `AbstractReadThreadingGraph.removeSingletonOrphanVertices`: every vertex with no edge at
    /// all. Unlike `BaseGraph`'s, it does not spare an isolated reference source; its test is
    /// `if (inDegreeOf(v) == 0 && outDegreeOf(v) == 0) {` (AbstractReadThreadingGraph.java:380).
    /// The orphans are listed first, in vertex order, then removed.
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
    /// in hash order, but the result of removing a set does not depend on the order: the
    /// survivors keep their relative order, and the k-mer entries dropped are those of every
    /// removed vertex's sequence whichever goes first.
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

impl ReadThreadingGraph {
    /// `getAdditionalSequence(v)`: the whole k-mer at a source, else its last base, which is
    /// `return source ? super.getAdditionalSequence(source) : getSuffixAsArray();`
    /// (MultiDeBruijnVertex.java:136) with `source` being `isSource(v)`, in-degree zero.
    pub fn additional_sequence(&self, vertex: usize) -> Vec<u8> {
        if self.in_degree(vertex) == 0 {
            self.vertices[vertex].sequence.clone()
        } else {
            vec![self.vertices[vertex].suffix()]
        }
    }

    /// `JunctionTreeLinkedDeBruijnGraph.getReferenceBytes(from, to, true, true)`: from the first
    /// occurrence of `from` on the reference path to the last of `to`:
    ///
    /// ```java
    /// int fromIndex = referencePath.indexOf(fromVertex);
    /// int toIndex = referencePath.lastIndexOf(toVertex);
    /// ```
    ///
    /// (`JunctionTreeLinkedDeBruijnGraph.java:370-371`.) The middle loop runs
    /// `for (int i = fromIndex + 1; i < toIndex; i++) {` (line 376), so a `from` that is not on
    /// the path (-1) starts it at the path's first vertex; the arithmetic is kept signed for that.
    pub fn linked_reference_bytes(&self, from: usize, to: usize) -> Vec<u8> {
        let path = self.reference_path.as_deref().unwrap_or(&[]);
        // `position` is `indexOf` and `rposition` is `lastIndexOf`, `None` for the -1.
        let from_index = path.iter().position(|&v| v == from);
        let to_index = path.iter().rposition(|&v| v == to);
        let mut bytes = self.additional_sequence_as_source(from);
        let start = from_index.map(|i| i as isize).unwrap_or(-1) + 1;
        let end = to_index.map(|i| i as isize).unwrap_or(-1);
        let mut i = start;
        while i < end {
            bytes.extend(self.additional_sequence(path[i as usize]));
            i += 1;
        }
        bytes.extend(self.additional_sequence(to));
        bytes
    }

    /// `getAdditionalSequence(v, true)`: the whole k-mer.
    fn additional_sequence_as_source(&self, vertex: usize) -> Vec<u8> {
        self.vertices[vertex].sequence.clone()
    }

    /// `addVertex` and `addEdge(source, target)` for the symbolic end the junction trees need: the
    /// edge factory's default edge, multiplicity one and not on the reference. The `_` vertex is
    /// added directly, so it never enters the k-mer map. Answers the new edge, which the junction
    /// trees compare against to recognise a read that ended on the reference sink.
    fn add_symbolic_end(&mut self, sink: usize) -> usize {
        let vertex = self.vertices.len();
        self.vertices.push(Vertex {
            sequence: vec![b'_'],
            additional_info: String::new(),
        });
        self.vertex_alive.push(true);
        self.outgoing.push(Vec::new());
        self.incoming.push(Vec::new());
        let edge = self.edges.len();
        self.add_edge(sink, vertex, false, 1);
        edge
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
///
/// The chain order is the answer's order, and everything downstream reads it: the adaptive
/// pruner's tie-breaks and its "chains left out, in chain order", and the `chain` rows of
/// `chain_pruner.txt.gz`. It comes from the sources in vertex order, a FIFO of chain ends
/// (`chainStarts.pop()` takes from the front of a `LinkedList`, `add` appends), and each vertex's
/// outgoing edges in their order.
pub fn find_all_chains(graph: &ReadThreadingGraph) -> Vec<Chain> {
    let mut starts: std::collections::VecDeque<usize> = graph.sources().into();
    // `alreadySeen`, a `HashSet` only asked membership.
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

/// `ChainPruner.findChain`: from the start edge, follow the only way out until (ChainPruner.java
/// line 52)
///
/// ```java
/// if (outEdges.size() != 1 || graph.inDegreeOf(lastVertex) > 1 || lastVertex.equals(firstVertex))
/// ```
///
/// is true. The test is made before the step, so a chain always ends ON the
/// vertex that stopped it.
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
///
/// Every chain is judged before any is removed, as `chainsToRemove(chains)` is computed whole
/// before `removeAllEdges` runs, so one removal cannot change another chain's verdict.
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

/// `Path.getBases()`: the first vertex's whole sequence, then each following vertex's last base
/// (`MultiDeBruijnVertex.getAdditionalSequence(false)` is the suffix).
pub fn chain_bases(graph: &ReadThreadingGraph, chain: &Chain) -> Vec<u8> {
    let mut bases = graph
        .vertex(graph.edge(chain.edges[0]).source)
        .sequence
        .clone();
    for &edge in &chain.edges {
        bases.push(graph.vertex(graph.edge(edge).target).suffix());
    }
    bases
}

/// `BaseUtils.BASES_COMPARATOR`: signed bytes in order, then length. The bytes are compared as
/// Java's signed `byte`, `final int cmp = Byte.compare(o1[i],o2[i]);` (BaseUtils.java:317),
/// hence `as i8`; for bases, which are all below 128, signed and unsigned agree.
///
/// `pub(crate)` makes it visible to this crate only, Java's package-private.
pub(crate) fn compare_bases(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    // `zip` pairs the two arrays element by element and stops at the shorter: `minLength`.
    for (x, y) in a.iter().zip(b) {
        let order = (*x as i8).cmp(&(*y as i8));
        if order != std::cmp::Ordering::Equal {
            return order;
        }
    }
    a.len().cmp(&b.len())
}

/// `Double.compare`: `-0.0` below `0.0`, NaN above everything, NaN equal to NaN.
pub(crate) fn java_double_compare(a: f64, b: f64) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    if a < b {
        return Ordering::Less;
    }
    if a > b {
        return Ordering::Greater;
    }
    // Equal, or a NaN involved: the bit patterns as signed 64-bit integers decide, which puts
    // `-0.0` (sign bit set, so negative) below `0.0` and the canonical NaN above every number.
    let (a, b) = (java_double_bits(a), java_double_bits(b));
    a.cmp(&b)
}

/// `Double.doubleToLongBits`: every NaN collapsed to the canonical one.
fn java_double_bits(value: f64) -> i64 {
    if value.is_nan() {
        0x7ff8_0000_0000_0000
    } else {
        value.to_bits() as i64
    }
}

/// `java.util.PriorityQueue` with a comparator: the same binary heap, `siftUp` on `offer` and
/// `siftDown` on `poll`, so that elements the comparator calls equal come out in the reference's
/// order and not merely in some order.
///
/// Its comparators can tie: two JT paths of exactly the same score, adaptive-pruner entries of
/// equal log odds and equal bases, and `createMapOfPivotalEdgesInTopologicalOrder`, whose every
/// entry compares equal. Rust's own `BinaryHeap` breaks ties its own way.
///
/// `F: Fn(&T, &T) -> Ordering` is the `Comparator<T>`: any function or lambda of that shape.
pub(crate) struct JavaPriorityQueue<T, F: Fn(&T, &T) -> std::cmp::Ordering> {
    heap: Vec<T>,
    compare: F,
}

impl<T: Clone, F: Fn(&T, &T) -> std::cmp::Ordering> JavaPriorityQueue<T, F> {
    pub(crate) fn new(compare: F) -> Self {
        JavaPriorityQueue {
            heap: Vec::new(),
            compare,
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    /// `size()`.
    pub(crate) fn len(&self) -> usize {
        self.heap.len()
    }

    /// `offer` / `siftUpUsingComparator`.
    pub(crate) fn add(&mut self, item: T) {
        let mut k = self.heap.len();
        self.heap.push(item.clone());
        while k > 0 {
            let parent = (k - 1) >> 1;
            if (self.compare)(&item, &self.heap[parent]) != std::cmp::Ordering::Less {
                break;
            }
            self.heap[k] = self.heap[parent].clone();
            k = parent;
        }
        self.heap[k] = item;
    }

    /// `poll` / `siftDownUsingComparator`: answer the root, move the last element into its place
    /// and sift it down, stopping as soon as it is not greater than the smaller child. With every
    /// element equal nothing moves, so the order is root, then the last element, and so on.
    pub(crate) fn poll(&mut self) -> Option<T> {
        // `?` on an empty heap answers `None`, Java's `null` from `poll()`.
        let last = self.heap.pop()?;
        if self.heap.is_empty() {
            return Some(last);
        }
        let result = self.heap[0].clone();
        let n = self.heap.len();
        let mut k = 0;
        let half = n >> 1;
        while k < half {
            let mut child = 2 * k + 1;
            let right = child + 1;
            if right < n
                && (self.compare)(&self.heap[child], &self.heap[right])
                    == std::cmp::Ordering::Greater
            {
                child = right;
            }
            if (self.compare)(&last, &self.heap[child]) != std::cmp::Ordering::Greater {
                break;
            }
            self.heap[k] = self.heap[child].clone();
            k = child;
        }
        self.heap[k] = last;
        Some(result)
    }
}

/// The parameters of `AdaptiveChainPruner`.
#[derive(Debug, Clone, Copy)]
pub struct AdaptivePruning {
    /// `initialErrorProbability`: the error rate of the first pass.
    pub initial_error_probability: f64,
    /// `logOddsThreshold`: the bar for a chain to extend the good subgraph.
    pub log_odds_threshold: f64,
    /// `seedingLogOddsThreshold`: the stricter bar, on both sides, for a chain to seed it.
    pub seeding_log_odds_threshold: f64,
    /// `maxUnprunedVariants`.
    pub max_unpruned_variants: i32,
}

/// `AdaptiveChainPruner.chainLogOdds`: the log odds that the chain is real against the edges
/// leaving its first vertex (left) and entering its last (right), zero at a source or a sink.
///
/// `None` where `QualityUtils.errorProbToQual` throws on the error rate.
pub fn chain_log_odds(
    graph: &ReadThreadingGraph,
    chain: &Chain,
    error_rate: f64,
) -> Option<(f64, f64)> {
    let first = graph.edge(chain.edges[0]).source;
    let last = chain.last_vertex;
    // `MathUtils.sumIntFunction(edges, E::getMultiplicity)`.
    let multiplicity =
        |edges: &[usize]| -> i32 { edges.iter().map(|&e| graph.edge(e).multiplicity).sum() };
    let left_total = multiplicity(graph.outgoing_edges(first));
    let right_total = multiplicity(graph.incoming_edges(last));
    let left = graph.edge(chain.edges[0]).multiplicity;
    let right = graph
        .edge(*chain.edges.last().expect("a chain has an edge"))
        .multiplicity;
    // `graph.isSource(first)` and, below, `graph.isSink(last)`: no evidence to weigh, zero.
    let left_odds = if graph.in_degree(first) == 0 {
        0.0
    } else {
        crate::mutect_engine::log_likelihood_ratio_from_error(left_total - left, left, error_rate)?
    };
    let right_odds = if graph.out_degree(last) == 0 {
        0.0
    } else {
        crate::mutect_engine::log_likelihood_ratio_from_error(
            right_total - right,
            right,
            error_rate,
        )?
    };
    Some((left_odds, right_odds))
}

/// `AdaptiveChainPruner.likelyErrorChains`: grow a subgraph of good chains from seeds, along chains
/// whose log odds pass the threshold, and answer every chain left out, in chain order:
/// `return chains.stream().filter(c -> !goodChains.contains(c)).collect(Collectors.toList());`
/// (AdaptiveChainPruner.java:132). A chain is its index in `chains`.
///
/// The queue orders by log odds, highest first, then by the first vertex's bases, then by the
/// chain's bases, both with `BaseUtils.BASES_COMPARATOR` (lines 79-82). `None` where the error
/// rate is refused.
pub fn likely_error_chains(
    graph: &ReadThreadingGraph,
    chains: &[Chain],
    error_rate: f64,
    params: &AdaptivePruning,
) -> Option<Vec<usize>> {
    // `chainLogOdds`. Collecting into `Option<Vec<_>>` gives `None` as soon as any chain's odds
    // are `None`, and the `?` returns it: one refusal refuses the whole pass.
    let odds: Vec<(f64, f64)> = chains
        .iter()
        .map(|chain| chain_log_odds(graph, chain, error_rate))
        .collect::<Option<_>>()?;
    let first_vertex = |c: usize| graph.edge(chains[c].edges[0]).source;
    let starts_on_ref = |c: usize| graph.edge(chains[c].edges[0]).is_ref;

    // `ArrayListMultimap`s: per vertex, its chains in chain order.
    let n = graph.vertices.len();
    let mut seedable: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut good_incoming: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut good_outgoing: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (c, chain) in chains.iter().enumerate() {
        let (left, right) = odds[c];
        if right >= params.log_odds_threshold || starts_on_ref(c) {
            good_incoming[chain.last_vertex].push(c);
        }
        if left >= params.log_odds_threshold || starts_on_ref(c) {
            good_outgoing[first_vertex(c)].push(c);
        }
        if right >= params.seeding_log_odds_threshold && left >= params.seeding_log_odds_threshold {
            seedable[first_vertex(c)].push(c);
            seedable[chain.last_vertex].push(c);
        }
    }

    let bases: Vec<Vec<u8>> = chains
        .iter()
        .map(|chain| chain_bases(graph, chain))
        .collect();
    let mut queue = JavaPriorityQueue::new(|a: &(usize, f64), b: &(usize, f64)| {
        java_double_compare(-a.1, -b.1)
            .then_with(|| {
                compare_bases(
                    &graph.vertex(first_vertex(a.0)).sequence,
                    &graph.vertex(first_vertex(b.0)).sequence,
                )
            })
            .then_with(|| compare_bases(&bases[a.0], &bases[b.0]))
    });

    queue.add((max_weight_chain(graph, chains), f64::INFINITY));
    // The seeding loop walks the multimap's keys in hash order; it only offers chains to the
    // queue and marks vertices processed, and neither depends on the order of the offers.
    let mut processed = vec![false; n];
    for vertex in graph.vertex_ids() {
        if seedable[vertex].len() > 2 {
            for &c in &good_outgoing[vertex] {
                queue.add((c, odds[c].0));
            }
            for &c in &good_incoming[vertex] {
                queue.add((c, odds[c].1));
            }
            processed[vertex] = true;
        }
    }

    // `goodChains` and `verticesThatAlreadyHaveOutgoingGoodChains`, `LinkedHashSet`s whose order
    // is never read: flags per chain and per vertex.
    let mut good = vec![false; chains.len()];
    let mut has_good_outgoing = vec![false; n];
    let mut variant_count = 0;
    while !queue.is_empty() && variant_count <= params.max_unpruned_variants {
        let (c, _) = queue.poll().expect("not empty");
        // `if (!goodChains.add(chain)) continue;`: a chain queued from both sides counts once.
        if good[c] {
            continue;
        }
        good[c] = true;
        let first = first_vertex(c);
        // `!verticesThatAlreadyHaveOutgoingGoodChains.add(first)`: a second good chain out of
        // the same vertex is a variant.
        let new_variant = has_good_outgoing[first];
        has_good_outgoing[first] = true;
        if new_variant {
            variant_count += 1;
        }
        if new_variant && variant_count > params.max_unpruned_variants {
            continue;
        }
        for vertex in [first, chains[c].last_vertex] {
            if !processed[vertex] {
                for &other in &good_outgoing[vertex] {
                    queue.add((other, odds[other].0));
                }
                for &other in &good_incoming[vertex] {
                    queue.add((other, odds[other].1));
                }
                processed[vertex] = true;
            }
        }
    }
    Some((0..chains.len()).filter(|&c| !good[c]).collect())
}

/// `getMaxWeightChain`: the chain holding the heaviest edge, then the longest, then the first
/// vertex's bases; `Stream.max` keeps the FIRST of equals.
fn max_weight_chain(graph: &ReadThreadingGraph, chains: &[Chain]) -> usize {
    let key = |c: usize| {
        let heaviest = chains[c]
            .edges
            .iter()
            .map(|&e| graph.edge(e).multiplicity)
            .max()
            .unwrap_or(0);
        (heaviest, chains[c].edges.len())
    };
    let mut best = 0;
    for c in 1..chains.len() {
        let order = key(c).cmp(&key(best)).then_with(|| {
            compare_bases(
                &graph.vertex(graph.edge(chains[c].edges[0]).source).sequence,
                &graph
                    .vertex(graph.edge(chains[best].edges[0]).source)
                    .sequence,
            )
        });
        // `Stream.max` keeps the running best unless a later chain compares strictly greater,
        // so the first of equal chains wins; hence `Greater` and not `!= Less`.
        if order == std::cmp::Ordering::Greater {
            best = c;
        }
    }
    best
}

/// `AdaptiveChainPruner.chainsToRemove`: a first pass at the initial error rate estimates the
/// graph's own error rate (the multiplicity of the probable error chains' last edges over every
/// chain edge's), and a second pass at that rate decides; a chain touching the reference is kept.
pub fn adaptive_chains_to_remove(
    graph: &ReadThreadingGraph,
    chains: &[Chain],
    params: &AdaptivePruning,
) -> Result<Vec<usize>, GraphError> {
    if chains.is_empty() {
        return Ok(Vec::new());
    }
    let refused = || GraphError::InvalidSequence("errorRate must be good probability".to_string());
    let probable = likely_error_chains(graph, chains, params.initial_error_probability, params)
        .ok_or_else(refused)?;
    let error_count: i32 = probable
        .iter()
        .map(|&c| {
            graph
                .edge(*chains[c].edges.last().expect("an edge"))
                .multiplicity
        })
        .sum();
    let total_bases: i32 = chains
        .iter()
        .map(|chain| {
            chain
                .edges
                .iter()
                .map(|&e| graph.edge(e).multiplicity)
                .sum::<i32>()
        })
        .sum();
    // `(double) errorCount / totalBases`: both widened to `double` before dividing, so a zero
    // total gives NaN, which `errorProbToQual` refuses wherever the second pass uses the rate.
    let error_rate = f64::from(error_count) / f64::from(total_bases);
    let likely = likely_error_chains(graph, chains, error_rate, params).ok_or_else(refused)?;
    Ok(likely
        .into_iter()
        .filter(|&c| !chains[c].edges.iter().any(|&e| graph.edge(e).is_ref))
        .collect())
}

/// `AdaptiveChainPruner.pruneLowWeightChains`: `ChainPruner`'s removal with the adaptive choice.
pub fn prune_adaptive(
    graph: &mut ReadThreadingGraph,
    params: &AdaptivePruning,
) -> Result<(), GraphError> {
    let chains = find_all_chains(graph);
    let doomed = adaptive_chains_to_remove(graph, &chains, params)?;
    for c in doomed {
        for &edge in &chains[c].edges {
            graph.remove_edge(edge);
        }
    }
    graph.remove_singleton_orphan_vertices();
    Ok(())
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
