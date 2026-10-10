//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs`
//! (`KBestHaplotypeFinder`, `GraphBasedKBestHaplotypeFinder`, `KBestHaplotype`, `Path`), GATK
//! 4.6.2.0: the k best haplotypes of a sequence graph.
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

use crate::read_threading_graph::{compare_bases, java_double_compare, JavaPriorityQueue};
use crate::seq_graph::{SeqGraph, SeqGraphError};

/// `KBestHaplotype`: a path, its score, and the reference flag the reference never sets.
#[derive(Debug, Clone)]
pub struct KBestHaplotype {
    pub vertices: Vec<usize>,
    pub edges: Vec<usize>,
    pub score: f64,
    pub is_reference: bool,
    bases: Vec<u8>,
}

impl KBestHaplotype {
    /// `Path.getBases()`: a `SeqVertex` adds its whole sequence wherever it is on the path.
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    pub fn last_vertex(&self) -> usize {
        *self.vertices.last().expect("a path has a vertex")
    }
}

/// `CycleDetector.detectCycles()` over the whole graph, a self-loop included.
pub fn has_cycles(graph: &SeqGraph) -> bool {
    let n = graph.vertex_ids().max().map_or(0, |m| m + 1);
    let mut state = vec![0u8; n];
    for root in graph.vertex_ids() {
        if state[root] != 0 {
            continue;
        }
        let mut stack: Vec<(usize, usize)> = vec![(root, 0)];
        state[root] = 1;
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
pub fn find_best_haplotypes(
    graph: &SeqGraph,
    sources: &[usize],
    sinks: &[usize],
    k: usize,
) -> Result<Vec<KBestHaplotype>, SeqGraphError> {
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
                // `computeLogPenaltyScore`: the base-ten logarithm of the multiplicity minus that of
                // the total, both through jmath.
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
