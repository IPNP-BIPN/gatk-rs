//! `LowWeightChainPruner` and `removePathsNotConnectedToRef` on the read threading graph, against
//! the golden `tools/readfilter-conformance/ChainPrunerDump.java` wrote in the pinned container.
//!
//! Each case is built as the dump builds it, then printed after the build, after pruning and after
//! the removal of paths not connected to the reference, with the chains `findAllChains` found in
//! between. A vertex's index is its position in the current vertex set, as in the dump.

use gatk_corpus as corpus;
use gatk_engine::read_threading_graph::{
    find_all_chains, prune_low_weight_chains, GraphError, ReadThreadingGraph,
};

const REF: &str = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";

struct Case {
    label: &'static str,
    k: usize,
    pruning_samples: usize,
    prune_factor: i32,
    reference: String,
    reads: Vec<(&'static str, String, usize)>,
}

fn case(
    label: &'static str,
    k: usize,
    pruning_samples: usize,
    prune_factor: i32,
    reference: &str,
    reads: Vec<(&'static str, String, usize)>,
) -> Case {
    Case {
        label,
        k,
        pruning_samples,
        prune_factor,
        reference: reference.to_string(),
        reads,
    }
}

fn cases() -> Vec<Case> {
    let snp = format!("{}T{}", &REF[..20], &REF[21..]);
    let snp2 = format!("{}A{}", &REF[..12], &REF[13..]);
    let tail = format!("{}TTTTTTTT", &REF[..32]);
    let head = format!("GGGGGGGG{}", &REF[8..]);
    let before = format!("TTCCAAGGTT{}", &REF[..20]);
    let rearranged = format!("{}{}", &REF[25..], &REF[..15]);
    vec![
        case(
            "snp-once-factor-2",
            10,
            1,
            2,
            REF,
            vec![("s1", snp.clone(), 1)],
        ),
        case(
            "snp-twice-factor-2",
            10,
            1,
            2,
            REF,
            vec![("s1", snp.clone(), 2)],
        ),
        case(
            "snp-once-factor-0",
            10,
            1,
            0,
            REF,
            vec![("s1", snp.clone(), 1)],
        ),
        case(
            "snp-once-factor-1",
            10,
            1,
            1,
            REF,
            vec![("s1", snp.clone(), 1)],
        ),
        case(
            "two-bubbles-factor-3",
            10,
            1,
            3,
            REF,
            vec![("s1", snp.clone(), 3), ("s1", snp2, 2)],
        ),
        case(
            "dangling-tail",
            10,
            1,
            2,
            REF,
            vec![("s1", tail.clone(), 1)],
        ),
        case("dangling-head", 10, 1, 2, REF, vec![("s1", head, 1)]),
        case("heavy-dangling-tail", 10, 1, 2, REF, vec![("s1", tail, 4)]),
        case(
            "heavy-prefix-before-reference",
            10,
            1,
            2,
            REF,
            vec![("s1", before, 3)],
        ),
        case(
            "two-samples-two-pruning-samples",
            10,
            2,
            2,
            REF,
            vec![("s1", snp.clone(), 1), ("s2", snp.clone(), 1)],
        ),
        case(
            "two-samples-one-pruning-sample",
            10,
            1,
            2,
            REF,
            vec![("s1", snp.clone(), 1), ("s2", snp.clone(), 1)],
        ),
        case(
            "mixed-weights",
            10,
            1,
            2,
            REF,
            vec![("s1", snp.clone(), 1), ("s1", snp[..28].to_string(), 1)],
        ),
        case("ref-only", 10, 1, 2, REF, vec![]),
        case(
            "rearranged-read",
            5,
            1,
            2,
            REF,
            vec![("s1", rearranged.clone(), 1)],
        ),
        case(
            "rearranged-read-heavy",
            5,
            1,
            2,
            REF,
            vec![("s1", rearranged, 3)],
        ),
        case("k5-snp-once", 5, 1, 2, REF, vec![("s1", snp, 1)]),
    ]
}

fn positions(graph: &ReadThreadingGraph) -> Vec<Option<usize>> {
    let mut index = vec![None; graph.vertex_ids().max().map_or(0, |m| m + 1)];
    for (position, id) in graph.vertex_ids().enumerate() {
        index[id] = Some(position);
    }
    index
}

fn print(out: &mut Vec<String>, label: &str, stage: &str, graph: &ReadThreadingGraph) {
    let index = positions(graph);
    let at = |id: usize| index.get(id).copied().flatten();
    for id in graph.vertex_ids() {
        let v = graph.vertex(id);
        out.push(format!(
            "vertex\t{label}\t{stage}\t{}\t{}\t{}",
            at(id).unwrap(),
            String::from_utf8_lossy(&v.sequence),
            v.additional_info
        ));
    }
    for id in graph.edge_ids() {
        let e = graph.edge(id);
        out.push(format!(
            "edge\t{label}\t{stage}\t{}\t{}\t{}\t{}\t{}",
            at(e.source).unwrap(),
            at(e.target).unwrap(),
            e.multiplicity,
            e.pruning_multiplicity(),
            e.is_ref
        ));
    }
    for (kmer, vertex) in graph.kmer_to_vertex() {
        let shown = at(*vertex).map_or("removed".to_string(), |p| p.to_string());
        out.push(format!(
            "kmer\t{label}\t{stage}\t{}\t{shown}",
            String::from_utf8_lossy(kmer)
        ));
    }
    let show = |v: Option<usize>| v.and_then(at).map_or("-".to_string(), |p| p.to_string());
    out.push(format!(
        "summary\t{label}\t{stage}\tvertices={}\tedges={}\tcycles={}\trefsource={}\trefsink={}",
        graph.vertex_count(),
        graph.edge_count(),
        graph.has_cycles(),
        show(graph.reference_source_vertex()),
        show(graph.reference_sink_vertex())
    ));
}

fn render(c: &Case) -> Vec<String> {
    let mut out = vec![format!(
        "case\t{}\tk={}\tminbq=10\tpruning={}\tprunefactor={}",
        c.label, c.k, c.pruning_samples, c.prune_factor
    )];
    let mut graph = ReadThreadingGraph::new(c.k, 10, c.pruning_samples);
    graph.add_sequence(c.reference.as_bytes(), true).unwrap();
    for (sample, bases, copies) in &c.reads {
        for _ in 0..*copies {
            graph
                .add_read(sample, bases.as_bytes(), &vec![30; bases.len()])
                .unwrap();
        }
    }
    graph.build_graph_if_necessary().unwrap();
    print(&mut out, c.label, "built", &graph);

    let index = positions(&graph);
    for chain in find_all_chains(&graph) {
        let vertices: Vec<String> = chain
            .vertices(&graph)
            .iter()
            .map(|&v| index[v].unwrap().to_string())
            .collect();
        let multiplicities: Vec<String> = chain
            .edges
            .iter()
            .map(|&e| graph.edge(e).multiplicity.to_string())
            .collect();
        let pruning: Vec<String> = chain
            .edges
            .iter()
            .map(|&e| graph.edge(e).pruning_multiplicity().to_string())
            .collect();
        let removed = chain.edges.iter().all(|&e| {
            graph.edge(e).pruning_multiplicity() < c.prune_factor && !graph.edge(e).is_ref
        });
        out.push(format!(
            "chain\t{}\t{}\t{}\t{}\t{removed}",
            c.label,
            vertices.join(","),
            multiplicities.join(","),
            pruning.join(",")
        ));
    }
    prune_low_weight_chains(&mut graph, c.prune_factor);
    print(&mut out, c.label, "pruned", &graph);

    match graph.remove_paths_not_connected_to_ref() {
        Ok(()) => print(&mut out, c.label, "connected", &graph),
        Err(error) => out.push(format!(
            "error\t{}\tconnected\t{}",
            c.label,
            describe(&error)
        )),
    }
    out
}

fn describe(error: &GraphError) -> String {
    format!("IllegalStateException: {error:?}")
}

#[test]
fn every_case_prunes_and_disconnects_as_the_reference_does() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/chain_pruner.txt.gz"),
    );
    let expected: Vec<&str> = golden.lines().filter(|l| !l.starts_with('#')).collect();
    let actual: Vec<String> = cases().iter().flat_map(render).collect();
    for (i, (want, got)) in expected.iter().zip(&actual).enumerate() {
        assert_eq!(got, want, "line {i}");
    }
    assert_eq!(actual.len(), expected.len(), "line count");
}
