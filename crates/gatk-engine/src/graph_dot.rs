//! `BaseGraph.printGraph`, ported for `--graph-output` and `--debug-graph-transformations`
//! (GATK 4.6.2.0).
//!
//! The reference names each vertex by `System.identityHashCode`, so no two JVM runs agree on the
//! numbers; here a vertex's number is its index, and a comparison masks them. Everything else is
//! the reference's text: the edges in the graph's order, each a `print` with no newline unless it
//! is a reference edge, which is printed again in red with one; then each vertex with its added
//! sequence, its number and its debugging text; then a linked graph's junction trees.

use crate::read_threading_graph::ReadThreadingGraph;
use crate::seq_graph::SeqGraph;

/// `MultiDeBruijnVertex.toString()`.
pub fn read_threading_vertex_name(graph: &ReadThreadingGraph, vertex: usize) -> String {
    format!(
        "MultiDeBruijnVertex_id_{vertex}_seq_{}",
        String::from_utf8_lossy(&graph.vertex(vertex).sequence)
    )
}

/// `SeqVertex.toString()`.
fn seq_vertex_name(graph: &SeqGraph, vertex: usize) -> String {
    format!(
        "SeqVertex_id_{vertex}_seq_{}",
        String::from_utf8_lossy(graph.sequence(vertex))
    )
}

/// One edge: `print(edgeString)`, `print(label)`, and for a reference edge
/// `println(edgeString + " [color=red];")`.
fn edge_text(
    out: &mut String,
    source: &str,
    target: &str,
    multiplicity: i32,
    label: &str,
    prune_factor: i32,
    is_ref: bool,
) {
    let edge = format!("\t{source} -> {target} ");
    out.push_str(&edge);
    if multiplicity > 0 && multiplicity < prune_factor {
        out.push_str(&format!("[style=dotted,color=grey,label=\"{label}\"];"));
    } else {
        out.push_str(&format!("[label=\"{label}\"];"));
    }
    if is_ref {
        out.push_str(&edge);
        out.push_str(" [color=red];\n");
    }
}

/// `printGraph(writer, writeHeader, pruneFactor)` for a read threading graph, with its extra lines.
pub fn read_threading_graph(
    graph: &ReadThreadingGraph,
    write_header: bool,
    prune_factor: i32,
    extra: &[String],
) -> String {
    let mut out = String::new();
    if write_header {
        out.push_str("digraph assemblyGraphs {\n");
    }
    for e in graph.edge_ids() {
        let edge = graph.edge(e);
        edge_text(
            &mut out,
            &read_threading_vertex_name(graph, edge.source),
            &read_threading_vertex_name(graph, edge.target),
            edge.multiplicity,
            &format!("{}/{}", edge.multiplicity, edge.pruning_multiplicity()),
            prune_factor,
            edge.is_ref,
        );
    }
    for v in graph.vertex_ids() {
        let vertex = graph.vertex(v);
        // `getAdditionalSequence(isSource(v))`: a source's whole k-mer, any other vertex's last base.
        let added = if graph.in_degree(v) == 0 {
            String::from_utf8_lossy(&vertex.sequence).into_owned()
        } else {
            char::from(vertex.suffix()).to_string()
        };
        out.push_str(&format!(
            "\t{} [label=\"{added} ({v}) {}\",shape=box]\n",
            read_threading_vertex_name(graph, v),
            vertex.additional_info
        ));
    }
    for line in extra {
        out.push_str(line);
        out.push('\n');
    }
    if write_header {
        out.push_str("}\n");
    }
    out
}

/// `printGraph(writer, writeHeader, pruneFactor)` for a sequence graph.
pub fn seq_graph(graph: &SeqGraph, write_header: bool, prune_factor: i32) -> String {
    let mut out = String::new();
    if write_header {
        out.push_str("digraph assemblyGraphs {\n");
    }
    for e in graph.edge_ids() {
        let edge = graph.edge(e);
        let label = match edge.pruning {
            Some(pruning) => format!("{}/{pruning}", edge.multiplicity),
            None => edge.multiplicity.to_string(),
        };
        edge_text(
            &mut out,
            &seq_vertex_name(graph, edge.source),
            &seq_vertex_name(graph, edge.target),
            edge.multiplicity,
            &label,
            prune_factor,
            edge.is_ref,
        );
    }
    for v in graph.vertex_ids() {
        out.push_str(&format!(
            "\t{} [label=\"{} ({v}) {}\",shape=box]\n",
            seq_vertex_name(graph, v),
            String::from_utf8_lossy(graph.sequence(v)),
            graph.info(v)
        ));
    }
    if write_header {
        out.push_str("}\n");
    }
    out
}
