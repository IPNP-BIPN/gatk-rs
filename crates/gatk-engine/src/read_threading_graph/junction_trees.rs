//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.
//! JunctionTreeLinkedDeBruijnGraph` (its junction trees) and `haplotypecaller.graphs`
//! (`JTBestHaplotype`, `JunctionTreeKBestHaplotypeFinder`), GATK 4.6.2.0: the linked de Bruijn
//! graph's memory of which way each read went at each fork, and the search that follows it.
//!
//! # The trees
//!
//! A vertex any of whose outgoing edges leads into a vertex with more than one way in holds a
//! tree. Each read is threaded again through the finished graph; at every such vertex it passes it
//! opens (or re-enters) that vertex's tree, and at every fork it then takes, each tree it has open
//! records the edge as a child of its current node and moves to it. A read that ends on the
//! reference sink records the symbolic end edge the trees are given for that purpose.
//!
//! # One order the reference does not fix
//!
//! A node's children are a `HashMap` keyed by edge, and an edge hashes by identity, so the
//! reference walks them in an order that depends on the JVM's identity hashes. Here they are kept
//! in the order they were first recorded. The order reaches the search only through the queue,
//! where it decides which of two paths of exactly the same score is taken first.

use std::collections::{HashMap, HashSet};

use super::{java_double_compare, JavaPriorityQueue, ReadThreadingGraph};

/// `DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE`.
pub const DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD: i32 = 3;
/// `DEFAULT_MINIMUM_WEIGHT_FOR_JT_BRANCH_TO_NOT_BE_PRUNED`.
const DEFAULT_MINIMUM_WEIGHT_FOR_JT_BRANCH: i32 = 1;
/// `DEFAULT_MAX_ACCEPTABLE_DECISION_EDGES_WITHOUT_JT_GUIDANCE`.
const MAX_DECISION_EDGES_WITHOUT_JT_GUIDANCE: i32 = 5;
/// `DEFAULT_MAX_ACCEPTABLE_REPETITIONS_OF_A_KMER_IN_A_PATH`.
const MAX_REPETITIONS_OF_A_KMER_IN_A_PATH: usize = 1;
/// `DEFAULT_MAX_PATHS_TO_CONSIDER_WITHOUT_RESULT`.
const MAX_PATHS_WITHOUT_RESULT: usize = 1000;
/// `DEFAULT_MAX_PATHS_TO_EVER_CONSIDER`.
const MAX_PATHS_EVER: usize = 10000;

/// `ThreadingNode`: the edge that led here (none at a root), how many threads passed, and the
/// children by edge, in the order first recorded.
#[derive(Debug, Clone)]
struct ThreadingNode {
    prev_edge: Option<usize>,
    count: i32,
    children: Vec<(usize, usize)>,
}

/// The junction trees of a linked graph, once `postProcessForHaplotypeFinding` has run.
#[derive(Debug, Clone)]
pub struct JunctionTrees {
    nodes: Vec<ThreadingNode>,
    /// `readThreadingJunctionTrees`: each tree's vertex and root node.
    trees: Vec<(usize, usize)>,
    by_vertex: HashMap<usize, usize>,
    symbolic_end_edge: usize,
}

/// A refusal of the search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum JunctionTreeError {
    /// `GATKException("While constructing graph, there was an incongruity ...")`.
    Incongruity,
    /// `GATKException.ShouldNeverReachHereException`: a finished path that does not start on the
    /// reference source.
    NotFromSource,
}

impl JunctionTrees {
    fn child(&self, node: usize, edge: usize) -> Option<usize> {
        self.nodes[node]
            .children
            .iter()
            .find(|(e, _)| *e == edge)
            .map(|(_, n)| *n)
    }

    /// `ThreadingNode.addEdge`: the child for the edge, made if new, its count incremented.
    fn add_edge(&mut self, node: usize, edge: usize) -> usize {
        let next = match self.child(node, edge) {
            Some(next) => next,
            None => {
                let next = self.nodes.len();
                self.nodes.push(ThreadingNode {
                    prev_edge: Some(edge),
                    count: 0,
                    children: Vec::new(),
                });
                self.nodes[node].children.push((edge, next));
                next
            }
        };
        self.nodes[next].count += 1;
        next
    }

    /// `getJunctionTreeForNode`, after pruning: a tree whose root has a child.
    pub fn tree_for_vertex(&self, vertex: usize) -> Option<usize> {
        self.by_vertex
            .get(&vertex)
            .copied()
            .filter(|&tree| !self.nodes[self.trees[tree].1].children.is_empty())
    }

    /// `getExtraGraphFileLines()`: each tree (only those with a child once `pruned`) as its vertex's
    /// edge to a point, then its nodes depth first. The reference walks hash maps keyed by
    /// identity here, so its order is the JVM's; a node is named by its index in hex.
    pub fn dot_lines(&self, graph: &ReadThreadingGraph, pruned: bool) -> Vec<String> {
        let mut out = Vec::new();
        for &(vertex, root) in &self.trees {
            if pruned && self.nodes[root].children.is_empty() {
                continue;
            }
            out.push(format!(
                "\t{} -> TreadingNode_{root:x} [color=blue];",
                crate::graph_dot::read_threading_vertex_name(graph, vertex)
            ));
            out.push(format!("\tTreadingNode_{root:x} [shape=point];"));
            self.node_dot_lines(graph, root, &mut out);
        }
        out
    }

    /// `edgesForNodeRecursive`.
    fn node_dot_lines(&self, graph: &ReadThreadingGraph, node: usize, out: &mut Vec<String>) {
        for &(edge, child) in &self.nodes[node].children {
            out.push(format!(
                "\tTreadingNode_{node:x} -> TreadingNode_{child:x}[color=blue,label=\"{}\"]; ",
                self.nodes[child].count
            ));
            let target = graph.edge(edge).target;
            out.push(format!(
                "\tTreadingNode_{child:x} [label=\"{}\",shape=plaintext]",
                char::from(graph.vertex(target).suffix())
            ));
            self.node_dot_lines(graph, child, out);
        }
    }

    fn root(&self, tree: usize) -> usize {
        self.trees[tree].1
    }

    fn total_out(&self, node: usize) -> i32 {
        self.nodes[node]
            .children
            .iter()
            .map(|(_, n)| self.nodes[*n].count)
            .sum()
    }

    fn is_symbolic_end(&self, node: usize) -> bool {
        self.nodes[node].prev_edge == Some(self.symbolic_end_edge)
    }
}

/// `JunctionTreeThreadingHelper`.
struct ThreadingHelper {
    tracked: Vec<usize>,
}

impl ReadThreadingGraph {
    /// `postProcessForHaplotypeFinding`: the edges' reference indices (which nothing reads once
    /// `annotatePathBasedOnGraph` is disabled), the junction trees, and their pruning.
    pub fn post_process_for_haplotype_finding(&mut self) -> JunctionTrees {
        let sink = self
            .reference_sink_vertex()
            .expect("a built linked graph has a reference path");
        // `generateJunctionTrees`: the symbolic end, then every pending sequence threaded again.
        let symbolic_end_edge = self.add_symbolic_end(sink);
        let mut trees = JunctionTrees {
            nodes: Vec::new(),
            trees: Vec::new(),
            by_vertex: HashMap::new(),
            symbolic_end_edge,
        };
        let pending: Vec<_> = self
            .pending
            .iter()
            .flat_map(|(_, sequences)| sequences.iter().cloned())
            .collect();
        for sequence in &pending {
            self.thread_for_junction_trees(&mut trees, sequence, sink);
        }
        // `pruneJunctionTrees(1)` keeps every child (each was walked at least once); a tree whose
        // root has none is dropped, which `tree_for_vertex` applies.
        debug_assert_eq!(DEFAULT_MINIMUM_WEIGHT_FOR_JT_BRANCH, 1);
        trees
    }

    /// `threadSequenceForJunctionTree`.
    fn thread_for_junction_trees(
        &self,
        trees: &mut JunctionTrees,
        sequence: &super::SequenceForKmers,
        sink: usize,
    ) {
        if sequence.is_ref {
            return;
        }
        let k = self.kmer_size;
        let mut helper = ThreadingHelper {
            tracked: Vec::new(),
        };
        // `findStartForJunctionThreading`: `start .. stop - k`, exclusive.
        let Some(end) = sequence.stop.checked_sub(k) else {
            return;
        };
        let Some(start) = (sequence.start..end)
            .find(|&i| self.kmer_index.contains_key(&sequence.sequence[i..i + k]))
        else {
            return;
        };
        let mut last = self
            .kmer_vertex(&sequence.sequence[start..start + k], true)
            .expect("a tracked k-mer");
        let mut rediscover = false;
        for i in start + 1..=end {
            let mut vertex = if !rediscover {
                self.extend_junction_threading_by_one(
                    last,
                    &sequence.sequence,
                    i,
                    Some((&mut helper, &mut *trees)),
                )
            } else {
                self.kmer_vertex(&sequence.sequence[i..i + k], true)
            };
            if vertex.is_none() && !rediscover {
                // Recover the path if the way forward is unambiguous for up to k bases.
                let outgoing = &self.outgoing[last];
                let mut tentative = if outgoing.len() == 1 {
                    Some(self.edges[outgoing[0]].target)
                } else {
                    None
                };
                let mut j = i + 1;
                while j <= end && j <= i + k {
                    let Some(t) = tentative else { break };
                    tentative =
                        self.extend_junction_threading_by_one(t, &sequence.sequence, j, None);
                    j += 1;
                }
                if tentative.is_some() {
                    vertex = Some(self.edges[outgoing[0]].target);
                }
            }
            match vertex {
                Some(v) => {
                    last = v;
                    rediscover = false;
                }
                None => {
                    helper.tracked.clear();
                    rediscover = true;
                }
            }
        }
        if last == sink {
            let symbolic = trees.symbolic_end_edge;
            add_edge_to_tracked(trees, &mut helper, symbolic);
        }
    }

    /// `extendJunctionThreadingByOne`: the target of the first outgoing edge that ends in the next
    /// base; when trees are being altered, the vertex left opens its tree if it warrants one, and
    /// at a fork every open tree records the edge.
    fn extend_junction_threading_by_one(
        &self,
        previous: usize,
        sequence: &[u8],
        kmer_start: usize,
        alter: Option<(&mut ThreadingHelper, &mut JunctionTrees)>,
    ) -> Option<usize> {
        let outgoing = &self.outgoing[previous];
        let next_base = sequence[kmer_start + self.kmer_size - 1];
        for &edge in outgoing {
            let target = self.edges[edge].target;
            if self.vertices[target].suffix() == next_base {
                if let Some((helper, trees)) = alter {
                    // `addTreeIfNecessary`.
                    let warrants = outgoing
                        .iter()
                        .any(|&e| self.in_degree(self.edges[e].target) > 1);
                    if warrants {
                        let tree = match trees.by_vertex.get(&previous) {
                            Some(&tree) => tree,
                            None => {
                                let root = trees.nodes.len();
                                trees.nodes.push(ThreadingNode {
                                    prev_edge: None,
                                    count: 0,
                                    children: Vec::new(),
                                });
                                trees.trees.push((previous, root));
                                trees.by_vertex.insert(previous, trees.trees.len() - 1);
                                trees.trees.len() - 1
                            }
                        };
                        let root = trees.trees[tree].1;
                        trees.nodes[root].count += 1;
                        helper.tracked.push(root);
                    }
                    if outgoing.len() != 1 {
                        add_edge_to_tracked(trees, helper, edge);
                    }
                }
                return Some(target);
            }
        }
        None
    }
}

/// `addEdgeToJunctionTreeNodes`: every tracked node moves to its child for the edge.
fn add_edge_to_tracked(trees: &mut JunctionTrees, helper: &mut ThreadingHelper, edge: usize) {
    let tracked = std::mem::take(&mut helper.tracked);
    helper.tracked = tracked
        .into_iter()
        .map(|node| trees.add_edge(node, edge))
        .collect();
}

/// `JTBestHaplotype`: a path (its first vertex and edges), its score, the junction tree nodes it
/// follows and the trees it has visited, and the forks taken since a tree last guided it.
#[derive(Debug, Clone)]
pub struct JtPath {
    first: usize,
    pub edges: Vec<usize>,
    last: usize,
    pub score: f64,
    active: Vec<usize>,
    visited: HashSet<usize>,
    decisions_without_jt: i32,
}

impl JtPath {
    fn new(vertex: usize) -> JtPath {
        JtPath {
            first: vertex,
            edges: Vec::new(),
            last: vertex,
            score: 0.0,
            active: Vec::new(),
            visited: HashSet::new(),
            decisions_without_jt: 0,
        }
    }

    /// `getVertices()`.
    pub fn vertices(&self, graph: &ReadThreadingGraph) -> Vec<usize> {
        let mut out = vec![self.first];
        out.extend(self.edges.iter().map(|&e| graph.edge(e).target));
        out
    }

    /// `getBases()`: the first vertex's k-mer, then each target's last base.
    pub fn bases(&self, graph: &ReadThreadingGraph) -> Vec<u8> {
        let mut out = graph.vertex(self.first).sequence.clone();
        out.extend(
            self.edges
                .iter()
                .map(|&e| graph.vertex(graph.edge(e).target).suffix()),
        );
        out
    }

    fn contains_vertex(&self, graph: &ReadThreadingGraph, vertex: usize) -> bool {
        self.first == vertex || self.edges.iter().any(|&e| graph.edge(e).target == vertex)
    }

    fn occurrences(&self, graph: &ReadThreadingGraph, vertex: usize) -> usize {
        self.vertices(graph)
            .iter()
            .filter(|&&v| v == vertex)
            .count()
    }

    /// `JTBestHaplotype(previousPath, edgesToExtend, edgePenalty)`: the trees are copied as they
    /// are, not walked along the new edges.
    fn extended(&self, graph: &ReadThreadingGraph, edges: &[usize], penalty: f64) -> JtPath {
        let mut out = self.clone();
        out.edges.extend_from_slice(edges);
        out.last = graph.edge(*edges.last().expect("edges to extend")).target;
        out.score = self.score + penalty;
        out.decisions_without_jt = if out.active.is_empty() {
            self.decisions_without_jt
        } else {
            0
        };
        out
    }

    /// The private constructor: scored by the edge's share of the outgoing weight, every open
    /// tree walked along the chain's last edge.
    fn extended_at_fork(
        &self,
        graph: &ReadThreadingGraph,
        trees: &JunctionTrees,
        chain: &[usize],
        multiplicity: i32,
        total: i32,
        based_on_jt: bool,
    ) -> JtPath {
        let mut out = self.clone();
        out.edges.extend_from_slice(chain);
        let last_edge = *chain.last().expect("a chain");
        out.last = graph.edge(last_edge).target;
        // `computeLogPenaltyScore`, added as one term.
        out.score = self.score
            + (jmath::math::log10(f64::from(multiplicity)) - jmath::math::log10(f64::from(total)));
        // `traverseEdgeForAllTrees`.
        out.active = self
            .active
            .iter()
            .filter_map(|&node| trees.child(node, last_edge))
            .filter(|&node| !trees.nodes[node].children.is_empty())
            .collect();
        out.decisions_without_jt = if based_on_jt {
            0
        } else {
            self.decisions_without_jt + 1
        };
        out
    }

    /// `addJunctionTree`: a tree not yet visited and with evidence is opened.
    fn add_junction_tree(&mut self, trees: &JunctionTrees, tree: usize) {
        let root = trees.root(tree);
        if !self.visited.contains(&tree) && !trees.nodes[root].children.is_empty() {
            self.visited.insert(tree);
            self.active.push(root);
            self.decisions_without_jt = 0;
        }
    }

    /// `removeEmptyNodesAndReturnIterator`.
    fn remove_empty_nodes(&mut self, trees: &JunctionTrees) -> Vec<usize> {
        self.active.retain(|&node| trees.total_out(node) > 0);
        self.active.clone()
    }

    /// `hasStoppingEvidence`.
    fn has_stopping_evidence(&mut self, trees: &JunctionTrees, threshold: i32) -> bool {
        for node in self.remove_empty_nodes(trees) {
            let total = trees.total_out(node);
            if trees.nodes[node]
                .children
                .iter()
                .any(|&(_, child)| trees.is_symbolic_end(child))
            {
                return true;
            }
            if total >= threshold {
                return false;
            }
        }
        true
    }

    /// `getApplicableNextEdgesBasedOnJunctionTrees`.
    fn next_paths(
        &mut self,
        graph: &ReadThreadingGraph,
        trees: &JunctionTrees,
        chain: &[usize],
        outgoing: &[usize],
        threshold: i32,
    ) -> Result<Vec<JtPath>, JunctionTreeError> {
        let mut accounted: Vec<usize> = Vec::new();
        let mut output = Vec::new();
        for node in self.remove_empty_nodes(trees) {
            let total = trees.total_out(node);
            for &(edge, child) in &trees.nodes[node].children {
                if !outgoing.contains(&edge) {
                    return Err(JunctionTreeError::Incongruity);
                }
                if !trees.is_symbolic_end(child) && !accounted.contains(&edge) {
                    accounted.push(edge);
                    let mut next_chain = chain.to_vec();
                    next_chain.push(edge);
                    output.push(self.extended_at_fork(
                        graph,
                        trees,
                        &next_chain,
                        trees.nodes[child].count,
                        total,
                        true,
                    ));
                }
            }
            if total >= threshold {
                return Ok(output);
            }
        }
        let total: i32 = outgoing.iter().map(|&e| graph.edge(e).multiplicity).sum();
        for &edge in outgoing {
            let multiplicity = graph.edge(edge).multiplicity;
            if !accounted.contains(&edge) && total != 0 && multiplicity != 0 {
                if graph.edge(edge).is_ref
                    && multiplicity == 1
                    && !(accounted.is_empty() && outgoing.len() < 2)
                {
                    continue;
                }
                let mut next_chain = chain.to_vec();
                next_chain.push(edge);
                output.push(self.extended_at_fork(
                    graph,
                    trees,
                    &next_chain,
                    multiplicity,
                    total,
                    false,
                ));
            }
        }
        Ok(output)
    }
}

/// `JunctionTreeKBestHaplotypeFinder.findBestHaplotypes`, from the reference source to the
/// reference sink, `recover` being `experimentalEndRecoveryMode`.
pub fn find_best_haplotypes(
    graph: &ReadThreadingGraph,
    trees: &JunctionTrees,
    source: usize,
    sink: usize,
    max_haplotypes: usize,
    recover: bool,
) -> Result<Vec<JtPath>, JunctionTreeError> {
    let threshold = DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD;
    let mut unvisited: Vec<usize> = if recover {
        pivotal_edges_in_topological_order(graph, source)
    } else {
        Vec::new()
    };
    let mut result: Vec<JtPath> = Vec::new();
    let mut queue =
        JavaPriorityQueue::new(|a: &JtPath, b: &JtPath| java_double_compare(b.score, a.score));
    queue.add(JtPath::new(source));
    let mut chain_cache: HashMap<usize, Vec<usize>> = HashMap::new();
    while result.len() < max_haplotypes && (!queue.is_empty() || !unvisited.is_empty()) {
        let limit = if result.is_empty() {
            MAX_PATHS_WITHOUT_RESULT
        } else {
            MAX_PATHS_EVER
        };
        if queue.len() > limit {
            break;
        }
        if queue.is_empty() {
            enqueue_next_pivotal_edge(graph, trees, &mut unvisited, &result, &mut queue);
            continue;
        }
        let mut path = queue.poll().expect("a path");
        if path.decisions_without_jt > MAX_DECISION_EDGES_WITHOUT_JT_GUIDANCE {
            continue;
        }

        let mut vertex = path.last;
        let mut outgoing = graph.outgoing_edges(vertex).to_vec();
        let chain = match chain_cache.get(&path.last).filter(|c| !c.is_empty()) {
            Some(cached) => {
                vertex = graph.edge(*cached.last().expect("non-empty")).target;
                outgoing = graph.outgoing_edges(vertex).to_vec();
                cached.clone()
            }
            None => {
                let mut chain = Vec::new();
                while outgoing.len() == 1
                    && trees.tree_for_vertex(vertex).is_none()
                    && vertex != sink
                {
                    let edge = outgoing[0];
                    if chain.contains(&edge) {
                        break;
                    }
                    chain.push(edge);
                    vertex = graph.edge(edge).target;
                    outgoing = graph.outgoing_edges(vertex).to_vec();
                }
                chain_cache.insert(path.last, chain.clone());
                chain
            }
        };

        if let Some(tree) = trees.tree_for_vertex(vertex) {
            path.add_junction_tree(trees, tree);
        }
        if vertex == sink && path.has_stopping_evidence(trees, threshold) {
            let finished = if chain.is_empty() {
                path.clone()
            } else {
                path.extended(graph, &chain, 0.0)
            };
            // `reconcilePathMissingReferenceStartPositions`.
            if finished.first != source {
                return Err(JunctionTreeError::NotFromSource);
            }
            result.push(finished);
            for edge in &path.edges {
                if let Some(at) = unvisited.iter().position(|e| e == edge) {
                    unvisited.remove(at);
                }
            }
        }

        if outgoing.len() > 1 {
            let paths = path.next_paths(graph, trees, &chain, &outgoing, threshold)?;
            for candidate in paths {
                let keep = !candidate.active.is_empty()
                    || candidate.decisions_without_jt == 0
                    || candidate.occurrences(graph, candidate.last)
                        <= MAX_REPETITIONS_OF_A_KMER_IN_A_PATH;
                if keep {
                    queue.add(candidate);
                }
            }
        } else if !outgoing.is_empty() {
            if path.active.is_empty()
                && path.occurrences(graph, vertex) > MAX_REPETITIONS_OF_A_KMER_IN_A_PATH
            {
                // A loop with no tree to vouch for it.
            } else {
                let mut next_chain = chain.clone();
                next_chain.push(outgoing[0]);
                queue.add(path.extended(graph, &next_chain, 0.0));
            }
        }
    }
    Ok(result)
}

/// `enqueueNextPivotalEdge`: the first unvisited pivotal edge, appended to the best finished path
/// through its source, cut after that path's last arrival there.
fn enqueue_next_pivotal_edge<F: Fn(&JtPath, &JtPath) -> std::cmp::Ordering>(
    graph: &ReadThreadingGraph,
    trees: &JunctionTrees,
    unvisited: &mut Vec<usize>,
    result: &[JtPath],
    queue: &mut JavaPriorityQueue<JtPath, F>,
) {
    let first_edge = unvisited.remove(0);
    let pivotal = graph.edge(first_edge).source;
    // `Stream.max` keeps the first of equals.
    let mut best: Option<&JtPath> = None;
    for path in result.iter().filter(|p| p.contains_vertex(graph, pivotal)) {
        best = match best {
            Some(b) if java_double_compare(b.score, path.score) != std::cmp::Ordering::Less => {
                Some(b)
            }
            _ => Some(path),
        };
    }
    let Some(best) = best else {
        return;
    };
    let incoming: Vec<usize> = best
        .edges
        .iter()
        .copied()
        .filter(|&e| graph.edge(e).target == pivotal)
        .collect();
    let Some(&last_incoming) = incoming.last() else {
        return;
    };
    let cut = best
        .edges
        .iter()
        .rposition(|&e| e == last_incoming)
        .expect("present");
    let mut edges: Vec<usize> = best.edges[..=cut].to_vec();
    edges.push(first_edge);
    let mut path = JtPath::new(best.first).extended(graph, &edges, best.score);
    for vertex in path.vertices(graph) {
        if let Some(tree) = trees.tree_for_vertex(vertex) {
            path.visited.insert(tree);
        }
    }
    queue.add(path);
}

/// `createMapOfPivotalEdgesInTopologicalOrder`: a breadth-first walk from the source whose queue
/// scores every edge 0, so that it pops in the heap's order; at every fork, the edges out not yet
/// visited, bar reference edges walked only by the reference.
fn pivotal_edges_in_topological_order(graph: &ReadThreadingGraph, source: usize) -> Vec<usize> {
    let mut visited: HashSet<usize> = HashSet::new();
    let mut queue = JavaPriorityQueue::new(|_: &usize, _: &usize| std::cmp::Ordering::Equal);
    let mut output: Vec<usize> = Vec::new();
    for &edge in graph.outgoing_edges(source) {
        queue.add(edge);
    }
    while let Some(edge) = queue.poll() {
        let outgoing = graph.outgoing_edges(graph.edge(edge).target).to_vec();
        if outgoing.len() > 1 {
            for &e in &outgoing {
                let e_edge = graph.edge(e);
                if !(e_edge.is_ref && e_edge.multiplicity == 1)
                    && !visited.contains(&e)
                    && !output.contains(&e)
                {
                    output.push(e);
                }
            }
        }
        for &e in &outgoing {
            if !visited.contains(&e) {
                queue.add(e);
                visited.insert(e);
            }
        }
    }
    output
}
