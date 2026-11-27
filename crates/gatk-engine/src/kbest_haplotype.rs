//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs`
//! (`KBestHaplotypeFinder`, `GraphBasedKBestHaplotypeFinder`, `KBestHaplotype`, `Path`), GATK
//! 4.7.0.0 (all four files are unchanged since 4.6.2.0): the k best haplotypes of a sequence graph.
//!
//! Dijkstra's k-shortest-paths: a path leaves each source, each step adds `log10(multiplicity) -
//! log10(total outgoing multiplicity)` to its score, a `java.util.PriorityQueue` hands out the
//! highest score first and, between equal scores, the larger bases first (`BASES_COMPARATOR`
//! reversed), a vertex is extended at most `k` times, and a path reaching a sink is a haplotype.
//! The queue is the reference's own heap, so paths the comparator calls equal leave it in the
//! reference's order.
//!
//! Three things a reader would not guess:
//!
//! * `KBestHaplotype.isReference` starts `false` and is only ever `&=`-ed, so it is `false` for
//!   every haplotype; the assembler marks the reference haplotype elsewhere;
//! * a zero multiplicity makes a step `-Infinity`, an all-zero vertex makes it `NaN`, and
//!   `Double.compare` ranks NaN above everything, so reversed, a NaN path is handed out FIRST;
//! * a graph with a cycle is first cut by `removeCyclesAndVerticesThatDontLeadToSinks`, a
//!   depth-first walk whose set of "parent" vertices is never emptied on the way back up, so it
//!   cuts every edge into an already-visited vertex: cross edges as well as back edges.
//!
//! # Reading this file without reading Rust
//!
//! - vertices and edges are `usize` ids into the [`SeqGraph`], not objects; a Java
//!   `Set<V>` of visited vertices becomes a `Vec<bool>` indexed by id;
//! - `&[usize]` is a borrowed, read-only view of a list; `&mut Vec<bool>` a borrowed, writable
//!   one, the only way a recursive call can add to the caller's set;
//! - `Result<T, SeqGraphError>` is "a value or an error", how the reference's exceptions are
//!   modelled, and a trailing `?` returns the error to the caller at once;
//! - `graph.vertex_ids().max().map_or(0, |m| m + 1)` is "one more than the largest id, or zero
//!   for an empty graph": the length a by-id table needs;
//! - a function's last expression is its return value; `return` is only written for early exits.

use crate::read_threading_graph::{compare_bases, java_double_compare, JavaPriorityQueue};
use crate::seq_graph::{SeqGraph, SeqGraphError};

/// `KBestHaplotype`: a path, its score, and the reference flag the reference never sets.
///
/// ```java
/// public KBestHaplotype(final KBestHaplotype<V, E> p, final E edge,
///         final int totalOutgoingMultiplicity) {
///     super(p, edge);
///     score = p.score
///             + computeLogPenaltyScore( edge.getMultiplicity(), totalOutgoingMultiplicity);
///     isReference &= edge.isRef();
/// }
/// ```
///
/// (Quotes in this file are broken at token boundaries to fit its width.)
///
/// `isReference &=` reads the new object's own field, which starts `false`, so it stays `false`.
#[derive(Debug, Clone)]
pub struct KBestHaplotype {
    /// Every vertex of the path in order, the first being its source. `Path` keeps only the
    /// edges and the last vertex; the list is the same information, kept so it need not be
    /// rebuilt.
    pub vertices: Vec<usize>,
    /// `Path.edgesInOrder`.
    pub edges: Vec<usize>,
    /// `KBestHaplotype.score`: the sum of the step penalties, 0 at a source.
    pub score: f64,
    pub is_reference: bool,
    /// `Path.getBases()`, extended one vertex at a time as the path grows rather than recomputed
    /// on each call as the reference does; the queue's comparator reads it on every comparison.
    bases: Vec<u8>,
}

impl KBestHaplotype {
    /// `Path.getBases()`: a `SeqVertex` adds its whole sequence wherever it is on the path.
    ///
    /// `getBases` asks each vertex for `getAdditionalSequence(isSource)`, and `BaseVertex`
    /// answers `getSequence()` whatever the flag, so the bases are the concatenated sequences.
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    /// `Path.getLastVertex()`. A path is never empty, so `expect` (stop the program with this
    /// message on `None`) is never reached.
    pub fn last_vertex(&self) -> usize {
        *self.vertices.last().expect("a path has a vertex")
    }
}

/// `CycleDetector.detectCycles()` over the whole graph, a self-loop included.
///
/// The reference is jgrapht's `org.jgrapht.alg.CycleDetector`, whose source is not in the
/// pinned clones; only its yes-or-no answer is used, so this is a plain depth-first search with
/// three states per vertex: 0 not yet reached, 1 on the current path, 2 finished. An edge into
/// a vertex in state 1 closes a cycle.
///
/// The search keeps an explicit stack instead of recursing: each entry is a vertex and the index
/// of the next outgoing edge to try from it.
pub fn has_cycles(graph: &SeqGraph) -> bool {
    let n = graph.vertex_ids().max().map_or(0, |m| m + 1);
    let mut state = vec![0u8; n];
    for root in graph.vertex_ids() {
        if state[root] != 0 {
            continue;
        }
        let mut stack: Vec<(usize, usize)> = vec![(root, 0)];
        state[root] = 1;
        // `while let Some(..) = stack.last_mut()` loops while the stack has a top entry, binding
        // `v` to a copy of its vertex and `next` to a writable reference to its edge index.
        while let Some(&mut (v, ref mut next)) = stack.last_mut() {
            if *next < graph.outgoing_edges(v).len() {
                let target = graph.edge(graph.outgoing_edges(v)[*next]).target;
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
                state[v] = 2;
                stack.pop();
            }
        }
    }
    false
}

/// `findGuiltyVerticesAndEdgesToRemoveCycles`: depth first along outgoing edges in order; an edge
/// into a vertex already in `parents` is cut, and a vertex from which no sink is reached is
/// removed. `parents` is never emptied on the way back up.
///
/// ```java
/// final Set<E> outgoingEdges = graph.outgoingEdgesOf(currentVertex);
/// parentVertices.add(currentVertex);
/// ```
///
/// Nothing ever removes from `parentVertices`, so despite its name it is every vertex this
/// source's walk has entered; that is why cross edges are cut too (see the module note).
fn find_guilty(
    graph: &SeqGraph,
    current: usize,
    sinks: &[usize],
    edges_to_remove: &mut Vec<bool>,
    vertices_to_remove: &mut Vec<bool>,
    parents: &mut Vec<bool>,
) -> bool {
    if sinks.contains(&current) {
        return true;
    }
    parents[current] = true;
    let mut reaches_sink = false;
    for &edge in graph.outgoing_edges(current) {
        let child = graph.edge(edge).target;
        if parents[child] {
            edges_to_remove[edge] = true;
        } else {
            let child_reaches = find_guilty(
                graph,
                child,
                sinks,
                edges_to_remove,
                vertices_to_remove,
                parents,
            );
            reaches_sink = reaches_sink || child_reaches;
        }
    }
    if !reaches_sink {
        vertices_to_remove[current] = true;
    }
    reaches_sink
}

/// `removeCyclesAndVerticesThatDontLeadToSinks`: a copy of the graph without the guilty edges and
/// vertices. `IllegalStateException` where no source reaches a sink, or nothing would be cut.
///
/// ```java
/// foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(original, source, sinks,
///         edgesToRemove, vertexToRemove, parentVertices) || foundSomePath;
/// ```
///
/// The call is on the left of `||`, so it runs for every source even once a path is found; the
/// port keeps that order. Each source gets a fresh `parentVertices`, while the two removal sets
/// are shared across sources. They are `HashSet`s upstream, and only membership is read, so a
/// flag per id carries the same information.
fn remove_cycles(
    graph: &SeqGraph,
    sources: &[usize],
    sinks: &[usize],
) -> Result<SeqGraph, SeqGraphError> {
    let edges = graph.edge_ids().max().map_or(0, |m| m + 1);
    let vertices = graph.vertex_ids().max().map_or(0, |m| m + 1);
    let mut edges_to_remove = vec![false; edges];
    let mut vertices_to_remove = vec![false; vertices];
    let mut found = false;
    for &source in sources {
        let mut parents = vec![false; vertices];
        found = find_guilty(
            graph,
            source,
            sinks,
            &mut edges_to_remove,
            &mut vertices_to_remove,
            &mut parents,
        ) || found;
    }
    if !found {
        return Err(SeqGraphError::NoPathAfterRemovingCycles);
    }
    if !edges_to_remove.iter().any(|&x| x) && !vertices_to_remove.iter().any(|&x| x) {
        return Err(SeqGraphError::CannotRemoveCycles);
    }
    let mut result = graph.clone();
    for (edge, &remove) in edges_to_remove.iter().enumerate() {
        if remove {
            result.remove_edge(edge);
        }
    }
    for (vertex, &remove) in vertices_to_remove.iter().enumerate() {
        if remove {
            result.remove_vertex(vertex);
        }
    }
    Ok(result)
}

/// `new GraphBasedKBestHaplotypeFinder(graph, sources, sinks).findBestHaplotypes(k)`.
///
/// The constructor removes cycles first when the graph has any (`keepCycles()` is `false` for
/// this finder); the search then runs on that copy. The search itself:
///
/// ```java
/// if (vertexCounts.get(vertexToExtend).getAndIncrement() < maxNumberOfHaplotypes) {
/// ```
///
/// `getAndIncrement` answers the count from BEFORE the increment, hence `seen` below is read
/// before `counts[last]` is bumped. A path that reaches a sink is a result and is never
/// extended, even if the sink has outgoing edges.
pub fn find_best_haplotypes(
    graph: &SeqGraph,
    sources: &[usize],
    sinks: &[usize],
    k: usize,
) -> Result<Vec<KBestHaplotype>, SeqGraphError> {
    // `cut` is declared outside the `if` so the copy outlives it; `graph` is then rebound to
    // borrow either the copy or the original, which is how a Java reassignment of a parameter is
    // written when the new value has to be owned somewhere.
    let cut;
    let graph = if has_cycles(graph) {
        cut = remove_cycles(graph, sources, sinks)?;
        &cut
    } else {
        graph
    };

    let mut queue = JavaPriorityQueue::new(|a: &KBestHaplotype, b: &KBestHaplotype| {
        // `comparingDouble(score).reversed().thenComparing(getBases, BASES_COMPARATOR.reversed())`.
        java_double_compare(b.score, a.score).then_with(|| compare_bases(&b.bases, &a.bases))
    });
    for &source in sources {
        // `new KBestHaplotype<>(source, graph)` reaches `Path`'s
        // `Utils.validateArg(graph.containsVertex(initialVertex), ...)`, asked of the graph after
        // the cut, which may have removed a source that reaches no sink.
        if !graph.contains_vertex(source) {
            return Err(SeqGraphError::NoSuchVertex);
        }
        queue.add(KBestHaplotype {
            vertices: vec![source],
            edges: Vec::new(),
            score: 0.0,
            is_reference: false,
            bases: graph.sequence(source).to_vec(),
        });
    }

    let mut counts = vec![0usize; graph.vertex_ids().max().map_or(0, |m| m + 1)];
    let mut result = Vec::new();
    while !queue.is_empty() && result.len() < k {
        let path = queue.poll().expect("not empty");
        let last = path.last_vertex();
        if sinks.contains(&last) {
            result.push(path);
            continue;
        }
        let seen = counts[last];
        counts[last] += 1;
        if seen < k {
            let outgoing = graph.outgoing_edges(last);
            let total: i32 = outgoing.iter().map(|&e| graph.edge(e).multiplicity).sum();
            for &edge in outgoing {
                let target = graph.edge(edge).target;
                // `computeLogPenaltyScore`: the base-ten logarithm of the multiplicity minus that
                // of the total, both through jmath, as the reference takes both through `Math`'s
                // `log10` (`KBestHaplotype.computeLogPenaltyScore`).
                let step = jmath::math::log10(f64::from(graph.edge(edge).multiplicity))
                    - jmath::math::log10(f64::from(total));
                let mut vertices = path.vertices.clone();
                vertices.push(target);
                let mut edges = path.edges.clone();
                edges.push(edge);
                let mut bases = path.bases.clone();
                bases.extend_from_slice(graph.sequence(target));
                queue.add(KBestHaplotype {
                    vertices,
                    edges,
                    score: path.score + step,
                    is_reference: false,
                    bases,
                });
            }
        }
    }
    Ok(result)
}
