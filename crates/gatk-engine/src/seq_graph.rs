//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs` (`SeqGraph`,
//! `SeqVertex`, `BaseEdge`, `BaseGraph`'s clean-ups, `MergeDiamonds`, `MergeTails`,
//! `SplitCommonSuffices`, `MergeCommonSuffices`, `SharedVertexSequenceSplitter`,
//! `CommonSuffixSplitter`, `SharedSequenceMerger`, `GraphUtils`) and the sequence-graph half of
//! `readthreading.ReadThreadingAssembler` (`cleanupSeqGraph`), GATK 4.6.2.0.
//!
//! # The graph underneath
//!
//! The reference is a JGraphT 1.1.0 `DefaultDirectedGraph`, and three of its rules decide what the
//! transforms produce:
//!
//! * vertices (`SeqVertex` compares by identity) and edges keep insertion order, in the graph and in
//!   each vertex's incoming and outgoing lists, and a removal keeps the others' order;
//! * `addEdge(u, v, e)` does nothing when an edge `u -> v` already exists: no multiple edges. So
//!   when two predecessors that share a source are merged, the second edge to the merged vertex is
//!   refused and its multiplicity is lost, rather than added;
//! * `addEdge(u, v, e)` also does nothing when the edge OBJECT `e` is already in the graph, which
//!   is how `SharedVertexSequenceSplitter`, whose split graph lends its edges to the outer one, never
//!   adds the direct prefix-to-suffix edge twice.
//!
//! Here a vertex or an edge is its index, `add_edge` answers `None` where JGraphT answers `false`,
//! and the splitter keeps the split graph's edges as identities in its own map.

use crate::read_threading_graph::ReadThreadingGraph;

/// A refusal the reference throws.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeqGraphError {
    /// `IllegalStateException("Graph must have ref source and sink vertices")`.
    NoReferenceSourceOrSink,
    /// `removePathsNotConnectedToRef`'s sanity checks.
    MoreThanOneSink,
    MoreThanOneSource,
    /// `IllegalStateException("Infinite loop detected in simplification routines ...")`.
    InfiniteSimplification,
    /// `assertVertexExist`: an edge to a vertex the graph does not hold.
    NoSuchVertex,
    /// `KBestHaplotypeFinder`: "could not find any path from the source vertex to the sink vertex
    /// after removing cycles".
    NoPathAfterRemovingCycles,
    /// `KBestHaplotypeFinder`: "cannot find a way to remove the cycles".
    CannotRemoveCycles,
}

/// `BaseEdge`: a multiplicity and a reference flag.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeqEdge {
    pub source: usize,
    pub target: usize,
    pub is_ref: bool,
    pub multiplicity: i32,
}

/// `SeqGraph`.
#[derive(Debug, Clone, Default)]
pub struct SeqGraph {
    sequences: Vec<Vec<u8>>,
    vertex_alive: Vec<bool>,
    edges: Vec<SeqEdge>,
    edge_alive: Vec<bool>,
    outgoing: Vec<Vec<usize>>,
    incoming: Vec<Vec<usize>>,
}

impl SeqGraph {
    pub fn new() -> SeqGraph {
        SeqGraph::default()
    }

    /// `addVertex(new SeqVertex(sequence))`: always a new vertex, identity being the vertex.
    pub fn add_vertex(&mut self, sequence: &[u8]) -> usize {
        self.sequences.push(sequence.to_vec());
        self.vertex_alive.push(true);
        self.outgoing.push(Vec::new());
        self.incoming.push(Vec::new());
        self.sequences.len() - 1
    }

    /// `addEdge(source, target, new BaseEdge(isRef, multiplicity))`: `None` if an edge between the
    /// two already exists.
    pub fn add_edge(
        &mut self,
        source: usize,
        target: usize,
        is_ref: bool,
        multiplicity: i32,
    ) -> Option<usize> {
        if self.edge_between(source, target).is_some() {
            return None;
        }
        let edge = self.edges.len();
        self.edges.push(SeqEdge {
            source,
            target,
            is_ref,
            multiplicity,
        });
        self.edge_alive.push(true);
        self.outgoing[source].push(edge);
        self.incoming[target].push(edge);
        Some(edge)
    }

    pub fn vertex_ids(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.sequences.len()).filter(|&v| self.vertex_alive[v])
    }

    pub fn edge_ids(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.edges.len()).filter(|&e| self.edge_alive[e])
    }

    pub fn sequence(&self, vertex: usize) -> &[u8] {
        &self.sequences[vertex]
    }

    pub fn edge(&self, edge: usize) -> &SeqEdge {
        &self.edges[edge]
    }

    pub fn contains_vertex(&self, vertex: usize) -> bool {
        vertex < self.vertex_alive.len() && self.vertex_alive[vertex]
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

    pub fn in_degree(&self, vertex: usize) -> usize {
        self.incoming[vertex].len()
    }

    pub fn out_degree(&self, vertex: usize) -> usize {
        self.outgoing[vertex].len()
    }

    /// `getEdge(source, target)`.
    pub fn edge_between(&self, source: usize, target: usize) -> Option<usize> {
        self.outgoing[source]
            .iter()
            .copied()
            .find(|&e| self.edges[e].target == target)
    }

    /// `outgoingVerticesOf`: a `LinkedHashSet` of targets, in edge order.
    pub fn outgoing_vertices(&self, vertex: usize) -> Vec<usize> {
        let mut out: Vec<usize> = Vec::new();
        for &e in &self.outgoing[vertex] {
            let target = self.edges[e].target;
            if !out.contains(&target) {
                out.push(target);
            }
        }
        out
    }

    /// `incomingVerticesOf`.
    pub fn incoming_vertices(&self, vertex: usize) -> Vec<usize> {
        let mut out: Vec<usize> = Vec::new();
        for &e in &self.incoming[vertex] {
            let source = self.edges[e].source;
            if !out.contains(&source) {
                out.push(source);
            }
        }
        out
    }

    pub fn remove_edge(&mut self, edge: usize) {
        if !self.edge_alive[edge] {
            return;
        }
        self.edge_alive[edge] = false;
        let SeqEdge { source, target, .. } = self.edges[edge];
        self.outgoing[source].retain(|&e| e != edge);
        self.incoming[target].retain(|&e| e != edge);
    }

    /// `removeVertex`: its edges, then the vertex.
    pub fn remove_vertex(&mut self, vertex: usize) {
        if !self.contains_vertex(vertex) {
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
    }

    /// `isReferenceNode(v)`: a reference edge either way, or the only vertex.
    pub fn is_reference_node(&self, vertex: usize) -> bool {
        self.outgoing[vertex]
            .iter()
            .chain(self.incoming[vertex].iter())
            .any(|&e| self.edges[e].is_ref)
            || self.vertex_count() == 1
    }

    pub fn is_ref_source(&self, vertex: usize) -> bool {
        if self.incoming[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return false;
        }
        if self.outgoing[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return true;
        }
        self.vertex_count() == 1
    }

    pub fn is_ref_sink(&self, vertex: usize) -> bool {
        if self.outgoing[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return false;
        }
        if self.incoming[vertex].iter().any(|&e| self.edges[e].is_ref) {
            return true;
        }
        self.vertex_count() == 1
    }

    pub fn reference_source_vertex(&self) -> Option<usize> {
        self.vertex_ids().find(|&v| self.is_ref_source(v))
    }

    pub fn reference_sink_vertex(&self) -> Option<usize> {
        self.vertex_ids().find(|&v| self.is_ref_sink(v))
    }

    pub fn sources(&self) -> Vec<usize> {
        self.vertex_ids()
            .filter(|&v| self.in_degree(v) == 0)
            .collect()
    }

    pub fn sinks(&self) -> Vec<usize> {
        self.vertex_ids()
            .filter(|&v| self.out_degree(v) == 0)
            .collect()
    }

    /// `BaseGraph.toSequenceGraph()`: a sequence vertex per k-mer vertex, a source keeping its
    /// whole k-mer and every other vertex its last base, and a copy of every edge, both in the read
    /// threading graph's order.
    pub fn from_read_threading_graph(graph: &ReadThreadingGraph) -> SeqGraph {
        let mut seq = SeqGraph::new();
        let mut map = vec![usize::MAX; graph.vertex_ids().max().map_or(0, |m| m + 1)];
        for v in graph.vertex_ids() {
            let sequence = &graph.vertex(v).sequence;
            map[v] = if graph.in_degree(v) == 0 {
                seq.add_vertex(sequence)
            } else {
                seq.add_vertex(&sequence[sequence.len() - 1..])
            };
        }
        for e in graph.edge_ids() {
            let edge = graph.edge(e);
            seq.add_edge(
                map[edge.source],
                map[edge.target],
                edge.is_ref,
                edge.multiplicity,
            );
        }
        seq
    }

    /// `cleanNonRefPaths`: remove the non-reference edges reachable backward from the reference
    /// source, and forward from the reference sink, through non-reference edges, then the
    /// orphans. The reference walks a `HashSet` in identity-hash order; what it removes is that
    /// closure whatever the order.
    pub fn clean_non_ref_paths(&mut self) {
        if self.reference_source_vertex().is_none() || self.reference_sink_vertex().is_none() {
            return;
        }
        let source = self.reference_source_vertex().expect("checked");
        let mut pending: Vec<usize> = self.incoming[source].clone();
        while let Some(e) = pending.pop() {
            if self.edge_alive[e] && !self.edges[e].is_ref {
                pending.extend(self.incoming[self.edges[e].source].iter().copied());
                self.remove_edge(e);
            }
        }
        if let Some(sink) = self.reference_sink_vertex() {
            let mut pending: Vec<usize> = self.outgoing[sink].clone();
            while let Some(e) = pending.pop() {
                if self.edge_alive[e] && !self.edges[e].is_ref {
                    pending.extend(self.outgoing[self.edges[e].target].iter().copied());
                    self.remove_edge(e);
                }
            }
        }
        self.remove_singleton_orphan_vertices();
    }

    /// `BaseGraph.removeSingletonOrphanVertices`: no edge at all, and not the reference source.
    pub fn remove_singleton_orphan_vertices(&mut self) {
        let orphans: Vec<usize> = self
            .vertex_ids()
            .filter(|&v| {
                self.in_degree(v) == 0 && self.out_degree(v) == 0 && !self.is_ref_source(v)
            })
            .collect();
        for v in orphans {
            self.remove_vertex(v);
        }
    }

    fn reachable(&self, start: usize, backward: bool, forward: bool) -> Vec<bool> {
        let mut seen = vec![false; self.sequences.len()];
        let mut stack = vec![start];
        while let Some(v) = stack.pop() {
            if seen[v] {
                continue;
            }
            seen[v] = true;
            if backward {
                stack.extend(self.incoming[v].iter().map(|&e| self.edges[e].source));
            }
            if forward {
                stack.extend(self.outgoing[v].iter().map(|&e| self.edges[e].target));
            }
        }
        seen
    }

    /// `removeVerticesNotConnectedToRefRegardlessOfEdgeDirection`.
    pub fn remove_vertices_not_connected_to_ref_regardless_of_edge_direction(&mut self) {
        let keep = match self.reference_source_vertex() {
            Some(source) => self.reachable(source, true, true),
            None => vec![false; self.sequences.len()],
        };
        let doomed: Vec<usize> = self.vertex_ids().filter(|&v| !keep[v]).collect();
        for v in doomed {
            self.remove_vertex(v);
        }
    }

    /// `removePathsNotConnectedToRef`.
    pub fn remove_paths_not_connected_to_ref(&mut self) -> Result<(), SeqGraphError> {
        let (Some(source), Some(sink)) =
            (self.reference_source_vertex(), self.reference_sink_vertex())
        else {
            return Err(SeqGraphError::NoReferenceSourceOrSink);
        };
        let forward = self.reachable(source, false, true);
        let backward = self.reachable(sink, true, false);
        let doomed: Vec<usize> = self
            .vertex_ids()
            .filter(|&v| !(forward[v] && backward[v]))
            .collect();
        for v in doomed {
            self.remove_vertex(v);
        }
        if self.sinks().len() > 1 {
            return Err(SeqGraphError::MoreThanOneSink);
        }
        if self.sources().len() > 1 {
            return Err(SeqGraphError::MoreThanOneSource);
        }
        Ok(())
    }

    // ----------------------------------------------------------------------------------------
    // Zipping.

    /// `isLinearChainStart`.
    fn is_linear_chain_start(&self, v: usize) -> bool {
        self.out_degree(v) == 1
            && (self.in_degree(v) != 1 || self.out_degree(self.incoming_vertices(v)[0]) > 1)
    }

    /// `traceLinearChain`.
    fn trace_linear_chain(&self, start: usize) -> Vec<usize> {
        let mut chain = vec![start];
        let mut last_is_ref = self.is_reference_node(start);
        let mut last = start;
        loop {
            if self.out_degree(last) != 1 {
                break;
            }
            let target = self.edges[self.outgoing[last][0]].target;
            if self.in_degree(target) != 1 || last == target {
                break;
            }
            let target_is_ref = self.is_reference_node(target);
            if last_is_ref != target_is_ref {
                break;
            }
            chain.push(target);
            last = target;
            last_is_ref = target_is_ref;
        }
        chain
    }

    /// `zipLinearChains`.
    pub fn zip_linear_chains(&mut self) -> bool {
        let starts: Vec<usize> = self
            .vertex_ids()
            .filter(|&v| self.is_linear_chain_start(v))
            .collect();
        if starts.is_empty() {
            return false;
        }
        let mut merged = false;
        for start in starts {
            let chain = self.trace_linear_chain(start);
            merged |= self.merge_linear_chain(&chain);
        }
        merged
    }

    /// `mergeLinearChainVertex`.
    fn merge_linear_chain(&mut self, chain: &[usize]) -> bool {
        let first = chain[0];
        let last = *chain.last().expect("a chain has a vertex");
        if first == last {
            return false;
        }
        let sequence: Vec<u8> = chain
            .iter()
            .flat_map(|&v| self.sequences[v].clone())
            .collect();
        let added = self.add_vertex(&sequence);
        for e in self.outgoing[last].clone() {
            let SeqEdge {
                target,
                is_ref,
                multiplicity,
                ..
            } = self.edges[e];
            self.add_edge(added, target, is_ref, multiplicity);
        }
        for e in self.incoming[first].clone() {
            let SeqEdge {
                source,
                is_ref,
                multiplicity,
                ..
            } = self.edges[e];
            self.add_edge(source, added, is_ref, multiplicity);
        }
        for &v in chain {
            self.remove_vertex(v);
        }
        true
    }

    // ----------------------------------------------------------------------------------------
    // Simplification.

    /// `simplifyGraph()`.
    pub fn simplify_graph(&mut self) -> Result<(), SeqGraphError> {
        self.zip_linear_chains();
        let mut previous: Option<SeqGraph> = None;
        let mut i = 0;
        loop {
            if i > 100 {
                return Err(SeqGraphError::InfiniteSimplification);
            }
            if !self.simplify_graph_once()? {
                break;
            }
            if i > 5 {
                if let Some(prev) = &previous {
                    if graph_equals(prev, self) {
                        break;
                    }
                }
                previous = Some(self.clone());
            }
            i += 1;
        }
        Ok(())
    }

    fn simplify_graph_once(&mut self) -> Result<bool, SeqGraphError> {
        let mut did = false;
        did |= self.transform_until_complete(Transform::MergeDiamonds)?;
        did |= self.transform_until_complete(Transform::MergeTails)?;
        did |= self.transform_until_complete(Transform::SplitCommonSuffices)?;
        did |= self.transform_until_complete(Transform::MergeCommonSuffices)?;
        did |= self.zip_linear_chains();
        Ok(did)
    }

    /// `VertexBasedTransformer.transformUntilComplete`: try every vertex in order, and after the
    /// first that transforms, start again from the first vertex.
    fn transform_until_complete(&mut self, transform: Transform) -> Result<bool, SeqGraphError> {
        let mut did = false;
        // `SplitCommonSuffices.alreadySplit`, which lives as long as the transformer.
        let mut already_split = vec![false; 0];
        loop {
            let mut found = false;
            let vertices: Vec<usize> = self.vertex_ids().collect();
            for v in vertices {
                found = match transform {
                    Transform::MergeDiamonds => self.merge_diamonds(v)?,
                    Transform::MergeTails => self.merge_tails(v)?,
                    Transform::SplitCommonSuffices => {
                        if already_split.len() <= v {
                            already_split.resize(self.sequences.len().max(v + 1), false);
                        }
                        if already_split[v] {
                            false
                        } else {
                            already_split[v] = true;
                            self.split_common_suffix(v)
                        }
                    }
                    Transform::MergeCommonSuffices => self.merge_shared_sequence(v),
                };
                if found {
                    did = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }
        Ok(did)
    }

    /// `MergeDiamonds.tryToTransform`.
    fn merge_diamonds(&mut self, top: usize) -> Result<bool, SeqGraphError> {
        let middles = self.outgoing_vertices(top);
        if middles.len() <= 1 {
            return Ok(false);
        }
        let mut bottom: Option<usize> = None;
        for &m in &middles {
            if self.out_degree(m) < 1 || self.in_degree(m) != 1 {
                return Ok(false);
            }
            for t in self.outgoing_vertices(m) {
                match bottom {
                    None => bottom = Some(t),
                    Some(b) if b != t => return Ok(false),
                    _ => {}
                }
            }
        }
        let bottom = bottom.expect("every middle has an outgoing vertex");
        if self.in_degree(bottom) != middles.len() {
            return Ok(false);
        }
        let splitter = Splitter::new(self, &middles);
        // `meetsMinMergableSequenceForEitherPrefixOrSuffix(1)`.
        if splitter.prefix.is_empty() && splitter.suffix.is_empty() {
            return Ok(false);
        }
        splitter.split_and_update(self, Some(top), Some(bottom))?;
        Ok(true)
    }

    /// `MergeTails.tryToTransform`.
    fn merge_tails(&mut self, top: usize) -> Result<bool, SeqGraphError> {
        let tails = self.outgoing_vertices(top);
        if tails.len() <= 1 {
            return Ok(false);
        }
        for &t in &tails {
            if self.out_degree(t) != 0 || self.in_degree(t) > 1 {
                return Ok(false);
            }
        }
        let splitter = Splitter::new(self, &tails);
        // `MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES`.
        if splitter.suffix.len() < 10 {
            return Ok(false);
        }
        splitter.split_and_update(self, Some(top), None)?;
        Ok(true)
    }

    /// `CommonSuffixSplitter.split`.
    fn split_common_suffix(&mut self, v: usize) -> bool {
        let to_split = self.incoming_vertices(v);
        let Some(suffix) = self.common_suffix_template(v, &to_split) else {
            return false;
        };
        for &mid in &to_split {
            let suffix_v = self.add_vertex(&suffix);
            let prefix_size = self.sequences[mid].len() - suffix.len();
            let out = self.outgoing[mid][0];
            let out_edge = self.edges[out];
            let incoming_target = if prefix_size > 0 {
                let prefix = self.sequences[mid][..prefix_size].to_vec();
                let prefix_v = self.add_vertex(&prefix);
                self.add_edge(prefix_v, suffix_v, out_edge.is_ref, 1);
                prefix_v
            } else {
                suffix_v
            };
            self.add_edge(
                suffix_v,
                out_edge.target,
                out_edge.is_ref,
                out_edge.multiplicity,
            );
            for e in self.incoming[mid].clone() {
                let SeqEdge {
                    source,
                    is_ref,
                    multiplicity,
                    ..
                } = self.edges[e];
                self.add_edge(source, incoming_target, is_ref, multiplicity);
            }
        }
        for &mid in &to_split {
            self.remove_vertex(mid);
        }
        true
    }

    /// `CommonSuffixSplitter.commonSuffix(graph, v, toSplit)`.
    fn common_suffix_template(&self, v: usize, to_split: &[usize]) -> Option<Vec<u8>> {
        if to_split.len() < 2 || !self.safe_to_split(v, to_split) {
            return None;
        }
        let kmers: Vec<&[u8]> = to_split
            .iter()
            .map(|&m| self.sequences[m].as_slice())
            .collect();
        let min = kmers.iter().map(|k| k.len()).min().unwrap_or(0);
        let suffix_len = common_maximum_suffix_length(&kmers, min);
        let first = kmers[0];
        let suffix = first[first.len() - suffix_len..].to_vec();
        if suffix.is_empty() {
            return None;
        }
        // `wouldEliminateRefSource`: the first reference source among them decides.
        if let Some(&source) = to_split.iter().find(|&&m| self.is_ref_source(m)) {
            if self.sequences[source].len() == suffix.len() {
                return None;
            }
        }
        if to_split
            .iter()
            .all(|&m| self.sequences[m].len() == suffix.len())
        {
            return None;
        }
        Some(suffix)
    }

    /// `CommonSuffixSplitter.safeToSplit`.
    fn safe_to_split(&self, bottom: usize, to_merge: &[usize]) -> bool {
        let outgoing_of_bottom = self.outgoing_vertices(bottom);
        for &m in to_merge {
            if m == bottom
                || self.out_degree(m) != 1
                || !self.outgoing_vertices(m).contains(&bottom)
            {
                return false;
            }
            if outgoing_of_bottom.contains(&m) {
                return false;
            }
        }
        true
    }

    /// `SharedSequenceMerger.merge`.
    fn merge_shared_sequence(&mut self, v: usize) -> bool {
        let prevs = self.incoming_vertices(v);
        if !self.can_merge(v, &prevs) {
            return false;
        }
        let mut sequence = self.sequences[prevs[0]].clone();
        sequence.extend_from_slice(&self.sequences[v]);
        let new_v = self.add_vertex(&sequence);
        for &prev in &prevs {
            for e in self.incoming[prev].clone() {
                let SeqEdge {
                    source,
                    is_ref,
                    multiplicity,
                    ..
                } = self.edges[e];
                // Refused, and its multiplicity lost, when `source -> new_v` already exists.
                self.add_edge(source, new_v, is_ref, multiplicity);
            }
        }
        for e in self.outgoing[v].clone() {
            let SeqEdge {
                target,
                is_ref,
                multiplicity,
                ..
            } = self.edges[e];
            self.add_edge(new_v, target, is_ref, multiplicity);
        }
        for &prev in &prevs {
            self.remove_vertex(prev);
        }
        self.remove_vertex(v);
        true
    }

    /// `SharedSequenceMerger.canMerge`.
    fn can_merge(&self, v: usize, incoming: &[usize]) -> bool {
        let Some(&first) = incoming.first() else {
            return false;
        };
        for &prev in incoming {
            if self.sequences[prev] != self.sequences[first] {
                return false;
            }
            let outs = self.outgoing_vertices(prev);
            if outs.len() != 1 || outs[0] != v || self.in_degree(prev) == 0 {
                return false;
            }
        }
        true
    }
}

#[derive(Debug, Clone, Copy)]
enum Transform {
    MergeDiamonds,
    MergeTails,
    SplitCommonSuffices,
    MergeCommonSuffices,
}

/// `GraphUtils.commonMaximumPrefixLength`.
fn common_maximum_prefix_length(kmers: &[&[u8]]) -> usize {
    let min = kmers.iter().map(|k| k.len()).min().unwrap_or(0);
    for i in 0..min {
        let b = kmers[0][i];
        if kmers[1..].iter().any(|k| k[i] != b) {
            return i;
        }
    }
    min
}

/// `GraphUtils.commonMaximumSuffixLength`.
fn common_maximum_suffix_length(kmers: &[&[u8]], min: usize) -> usize {
    for suffix_len in 0..min {
        let first = kmers[0];
        let b = first[first.len() - suffix_len - 1];
        if kmers[1..].iter().any(|k| k[k.len() - suffix_len - 1] != b) {
            return suffix_len;
        }
    }
    min
}

/// `BaseGraph.graphEquals`: as many vertices and edges, every vertex's sequence found among the
/// other's, and every edge matched by one between vertices of the same sequences, both ways.
pub fn graph_equals(g1: &SeqGraph, g2: &SeqGraph) -> bool {
    if g1.vertex_count() != g2.vertex_count() || g1.edge_count() != g2.edge_count() {
        return false;
    }
    if !g1
        .vertex_ids()
        .all(|a| g2.vertex_ids().any(|b| g1.sequence(a) == g2.sequence(b)))
    {
        return false;
    }
    let same = |x: &SeqGraph, ex: usize, y: &SeqGraph, ey: usize| {
        x.sequence(x.edge(ex).source) == y.sequence(y.edge(ey).source)
            && x.sequence(x.edge(ex).target) == y.sequence(y.edge(ey).target)
    };
    g1.edge_ids()
        .all(|a| g2.edge_ids().any(|b| same(g1, a, g2, b)))
        && g2
            .edge_ids()
            .all(|b| g1.edge_ids().any(|a| same(g2, b, g1, a)))
}

/// `SharedVertexSequenceSplitter`: the shared prefix and suffix of the vertices to split, factored
/// out into their own vertices.
struct Splitter {
    to_splits: Vec<usize>,
    prefix: Vec<u8>,
    suffix: Vec<u8>,
}

/// A vertex of the split graph: the prefix, the suffix, or a new middle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SplitVertex {
    Prefix,
    Suffix,
    Middle(usize),
}

/// An edge of the split graph. Its index in the split graph's list is its identity, which the
/// outer graph remembers once it has taken the edge, so taking it again is a no-op.
#[derive(Debug, Clone, Copy)]
struct SplitEdge {
    source: SplitVertex,
    target: SplitVertex,
    is_ref: bool,
    multiplicity: i32,
}

impl Splitter {
    /// `new SharedVertexSequenceSplitter(graph, toSplits)`: `commonPrefixAndSuffixOfVertices`.
    fn new(graph: &SeqGraph, to_splits: &[usize]) -> Splitter {
        let kmers: Vec<&[u8]> = to_splits.iter().map(|&v| graph.sequence(v)).collect();
        let min = kmers.iter().map(|k| k.len()).min().unwrap_or(0);
        let prefix_len = common_maximum_prefix_length(&kmers);
        let suffix_len = common_maximum_suffix_length(&kmers, min - prefix_len);
        let first = kmers[0];
        Splitter {
            to_splits: to_splits.to_vec(),
            prefix: first[..prefix_len].to_vec(),
            suffix: first[first.len() - suffix_len..].to_vec(),
        }
    }

    /// `splitAndUpdate(top, bottom)`: `split()` then `updateGraph(top, bottom)`.
    fn split_and_update(
        &self,
        outer: &mut SeqGraph,
        top: Option<usize>,
        bottom: Option<usize>,
    ) -> Result<(), SeqGraphError> {
        // `split()`. The split graph: prefix, suffix, then each remaining middle, and its edges
        // in the order they were added; `addOrUpdateEdge` adds to an existing prefix -> suffix edge.
        let mut middles: Vec<Vec<u8>> = Vec::new();
        let mut edges: Vec<SplitEdge> = Vec::new();
        for &mid in &self.to_splits {
            let to_mid = Self::process_edge(outer, mid, true);
            let from_mid = Self::process_edge(outer, mid, false);
            let sequence = outer.sequence(mid);
            let start = self.prefix.len();
            let length =
                sequence.len() as isize - self.suffix.len() as isize - self.prefix.len() as isize;
            if length > 0 {
                let remaining = sequence[start..start + length as usize].to_vec();
                middles.push(remaining);
                let m = SplitVertex::Middle(middles.len() - 1);
                edges.push(SplitEdge {
                    source: SplitVertex::Prefix,
                    target: m,
                    is_ref: to_mid.0,
                    multiplicity: to_mid.1,
                });
                edges.push(SplitEdge {
                    source: m,
                    target: SplitVertex::Suffix,
                    is_ref: from_mid.0,
                    multiplicity: from_mid.1,
                });
            } else {
                // `toMid.copy().add(fromMid)`: the sum, and the OR of the reference flags.
                let combined = (to_mid.0 || from_mid.0, to_mid.1 + from_mid.1);
                match edges
                    .iter_mut()
                    .find(|e| e.source == SplitVertex::Prefix && e.target == SplitVertex::Suffix)
                {
                    Some(existing) => {
                        existing.multiplicity += combined.1;
                        existing.is_ref = existing.is_ref || combined.0;
                    }
                    None => edges.push(SplitEdge {
                        source: SplitVertex::Prefix,
                        target: SplitVertex::Suffix,
                        is_ref: combined.0,
                        multiplicity: combined.1,
                    }),
                }
            }
        }

        // `updateGraph(top, bot)`.
        for &v in &self.to_splits {
            outer.remove_vertex(v);
        }
        let middle_ids: Vec<usize> = middles.iter().map(|m| outer.add_vertex(m)).collect();

        let prefix_out = |edges: &[SplitEdge]| -> Vec<usize> {
            (0..edges.len())
                .filter(|&i| edges[i].source == SplitVertex::Prefix)
                .collect()
        };
        let suffix_in = |edges: &[SplitEdge]| -> Vec<usize> {
            (0..edges.len())
                .filter(|&i| edges[i].target == SplitVertex::Suffix)
                .collect()
        };
        let has_prefix_suffix_edge = edges
            .iter()
            .any(|e| e.source == SplitVertex::Prefix && e.target == SplitVertex::Suffix);
        let has_only_prefix_suffix_edges = has_prefix_suffix_edge && prefix_out(&edges).len() == 1;
        let need_prefix =
            !self.prefix.is_empty() || (top.is_none() && !has_only_prefix_suffix_edges);
        let need_suffix =
            !self.suffix.is_empty() || (bottom.is_none() && !has_only_prefix_suffix_edges);

        let mut prefix_id: Option<usize> = None;
        let mut suffix_id: Option<usize> = None;
        if need_prefix {
            let id = outer.add_vertex(&self.prefix);
            prefix_id = Some(id);
            if let Some(top) = top {
                let any_ref = prefix_out(&edges).iter().any(|&i| edges[i].is_ref);
                outer.add_edge(top, id, any_ref, 1);
            }
        }
        if need_suffix {
            let id = outer.add_vertex(&self.suffix);
            suffix_id = Some(id);
            if let Some(bottom) = bottom {
                let any_ref = suffix_in(&edges).iter().any(|&i| edges[i].is_ref);
                outer.add_edge(id, bottom, any_ref, 1);
            }
        }
        let top_for_connect = if need_prefix { prefix_id } else { top };
        let bottom_for_connect = if need_suffix { suffix_id } else { bottom };

        // Which split edges the outer graph already holds: JGraphT's `containsEdge(e)`.
        let mut taken = vec![false; edges.len()];
        let resolve = |v: SplitVertex| -> Option<usize> {
            match v {
                SplitVertex::Prefix => prefix_id,
                SplitVertex::Suffix => suffix_id,
                SplitVertex::Middle(i) => Some(middle_ids[i]),
            }
        };
        let add = |outer: &mut SeqGraph,
                   taken: &mut Vec<bool>,
                   i: usize,
                   source: Option<usize>,
                   target: Option<usize>|
         -> Result<(), SeqGraphError> {
            if taken[i] {
                return Ok(());
            }
            let (Some(source), Some(target)) = (source, target) else {
                return Err(SeqGraphError::NoSuchVertex);
            };
            if outer
                .add_edge(source, target, edges[i].is_ref, edges[i].multiplicity)
                .is_some()
            {
                taken[i] = true;
            }
            Ok(())
        };
        if let Some(top) = top_for_connect {
            for i in prefix_out(&edges) {
                if edges[i].target == SplitVertex::Suffix {
                    if let Some(bottom) = bottom_for_connect {
                        add(outer, &mut taken, i, Some(top), Some(bottom))?;
                    }
                } else {
                    add(outer, &mut taken, i, Some(top), resolve(edges[i].target))?;
                }
            }
        }
        if let Some(bottom) = bottom_for_connect {
            for i in suffix_in(&edges) {
                add(outer, &mut taken, i, resolve(edges[i].source), Some(bottom))?;
            }
        }
        Ok(())
    }

    /// `processEdgeToRemove(v, incomingEdgeOf(v))` (or the outgoing one): a copy of the single
    /// edge, or a new zero-multiplicity edge, reference if the vertex is a reference node.
    fn process_edge(outer: &SeqGraph, v: usize, incoming: bool) -> (bool, i32) {
        let edges = if incoming {
            outer.incoming_edges(v)
        } else {
            outer.outgoing_edges(v)
        };
        match edges.first() {
            Some(&e) => (outer.edge(e).is_ref, outer.edge(e).multiplicity),
            None => (outer.is_reference_node(v), 0),
        }
    }
}
