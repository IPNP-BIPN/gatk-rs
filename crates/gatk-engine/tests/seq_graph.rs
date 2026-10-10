//! `toSequenceGraph`, `cleanupSeqGraph` and `SeqGraph.simplifyGraph`, against the golden
//! `tools/readfilter-conformance/SeqGraphDump.java` wrote in the pinned container.
//!
//! Part one threads reads, prunes, disconnects and converts, then prints the sequence graph after
//! every stage of `cleanupSeqGraph`; part two builds graphs by hand to reach each transform.

use gatk_corpus as corpus;
use gatk_engine::read_threading_graph::{prune_low_weight_chains, ReadThreadingGraph};
use gatk_engine::seq_graph::{SeqGraph, SeqGraphError};

const REF: &str = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";

fn print(out: &mut Vec<String>, label: &str, stage: &str, g: &SeqGraph) {
    let mut index = vec![usize::MAX; g.vertex_ids().max().map_or(0, |m| m + 1)];
    for (position, id) in g.vertex_ids().enumerate() {
        index[id] = position;
        out.push(format!(
            "vertex\t{label}\t{stage}\t{position}\t{}",
            String::from_utf8_lossy(g.sequence(id))
        ));
    }
    for e in g.edge_ids() {
        let edge = g.edge(e);
        out.push(format!(
            "edge\t{label}\t{stage}\t{}\t{}\t{}\t{}",
            index[edge.source], index[edge.target], edge.multiplicity, edge.is_ref
        ));
    }
    let show = |v: Option<usize>| v.map_or("-".to_string(), |v| index[v].to_string());
    out.push(format!(
        "summary\t{label}\t{stage}\tvertices={}\tedges={}\trefsource={}\trefsink={}",
        g.vertex_count(),
        g.edge_count(),
        show(g.reference_source_vertex()),
        show(g.reference_sink_vertex())
    ));
}

fn error(out: &mut Vec<String>, label: &str, stage: &str, e: &SeqGraphError) {
    out.push(format!(
        "error\t{label}\t{stage}\tIllegalStateException: {e:?}"
    ));
}

fn cleanup(out: &mut Vec<String>, label: &str, seq: &mut SeqGraph) {
    seq.zip_linear_chains();
    print(out, label, "zipped", seq);
    seq.remove_singleton_orphan_vertices();
    seq.remove_vertices_not_connected_to_ref_regardless_of_edge_direction();
    print(out, label, "pruned", seq);
    if let Err(e) = seq.simplify_graph() {
        return error(out, label, "cleanup", &e);
    }
    print(out, label, "merged", seq);
    if seq.reference_source_vertex().is_none() || seq.reference_sink_vertex().is_none() {
        out.push(format!("status\t{label}\tJUST_ASSEMBLED_REFERENCE"));
        return;
    }
    if let Err(e) = seq.remove_paths_not_connected_to_ref() {
        return error(out, label, "cleanup", &e);
    }
    if let Err(e) = seq.simplify_graph() {
        return error(out, label, "cleanup", &e);
    }
    if seq.vertex_count() == 1 {
        let complete = seq.vertex_ids().next().unwrap();
        let dummy = seq.add_vertex(b"");
        seq.add_edge(complete, dummy, true, 0);
    }
    print(out, label, "final", seq);
    out.push(format!("status\t{label}\tASSEMBLED_SOME_VARIATION"));
}

fn from_reads(out: &mut Vec<String>, label: &str, k: usize, reads: &[(String, usize)]) {
    out.push(format!("case\t{label}"));
    let mut rtg = ReadThreadingGraph::new(k, 10, 1);
    rtg.add_sequence(REF.as_bytes(), true).unwrap();
    for (bases, copies) in reads {
        for _ in 0..*copies {
            rtg.add_read("s1", bases.as_bytes(), &vec![30; bases.len()])
                .unwrap();
        }
    }
    rtg.build_graph_if_necessary().unwrap();
    prune_low_weight_chains(&mut rtg, 2);
    rtg.remove_paths_not_connected_to_ref().unwrap();
    let mut seq = SeqGraph::from_read_threading_graph(&rtg);
    print(out, label, "seq", &seq);
    seq.clean_non_ref_paths();
    print(out, label, "cleaned", &seq);
    cleanup(out, label, &mut seq);
}

fn simplify(out: &mut Vec<String>, label: &str, g: &mut SeqGraph) {
    out.push(format!("case\t{label}"));
    print(out, label, "built", g);
    let zipped = g.zip_linear_chains();
    out.push(format!("result\t{label}\tzip\t{zipped}"));
    print(out, label, "zipped", g);
    match g.simplify_graph() {
        Ok(()) => print(out, label, "simplified", g),
        Err(e) => error(out, label, "simplify", &e),
    }
}

fn by_hand(out: &mut Vec<String>) {
    let tail = "GGGGCCCCAAAATTTT";
    let mut g = SeqGraph::new();
    let top = g.add_vertex(b"ACGTACGTAC");
    let t1 = g.add_vertex(format!("A{tail}").as_bytes());
    let t2 = g.add_vertex(format!("CC{tail}").as_bytes());
    let t3 = g.add_vertex(format!("TTT{tail}").as_bytes());
    g.add_edge(top, t1, true, 5);
    g.add_edge(top, t2, false, 2);
    g.add_edge(top, t3, false, 1);
    simplify(out, "tails", &mut g);

    let mut g = SeqGraph::new();
    let top = g.add_vertex(b"ACGTACGTAC");
    let t1 = g.add_vertex(b"AGG");
    let t2 = g.add_vertex(b"CGG");
    g.add_edge(top, t1, true, 5);
    g.add_edge(top, t2, false, 2);
    simplify(out, "short-tails", &mut g);

    let mut g = SeqGraph::new();
    let s1 = g.add_vertex(b"ACGTTT");
    let s2 = g.add_vertex(b"CGATTT");
    let m1 = g.add_vertex(b"AAACCC");
    let m2 = g.add_vertex(b"GGACCC");
    let bottom = g.add_vertex(b"TTGCA");
    g.add_edge(s1, m1, true, 4);
    g.add_edge(s2, m2, false, 3);
    g.add_edge(m1, bottom, true, 4);
    g.add_edge(m2, bottom, false, 3);
    simplify(out, "common-suffix", &mut g);

    let mut g = SeqGraph::new();
    let a = g.add_vertex(b"ACGT");
    let b = g.add_vertex(b"TTGG");
    let p1 = g.add_vertex(b"CCA");
    let p2 = g.add_vertex(b"CCA");
    let bottom = g.add_vertex(b"GAGA");
    g.add_edge(a, p1, true, 2);
    g.add_edge(b, p2, false, 2);
    g.add_edge(p1, bottom, true, 2);
    g.add_edge(p2, bottom, false, 2);
    simplify(out, "shared-sequence", &mut g);

    let mut g = SeqGraph::new();
    let a = g.add_vertex(b"ACGT");
    let p1 = g.add_vertex(b"CCA");
    let p2 = g.add_vertex(b"CCA");
    let bottom = g.add_vertex(b"GAGA");
    g.add_edge(a, p1, true, 2);
    g.add_edge(a, p2, false, 3);
    g.add_edge(p1, bottom, true, 2);
    g.add_edge(p2, bottom, false, 3);
    simplify(out, "shared-sequence-one-source", &mut g);

    let mut g = SeqGraph::new();
    let top = g.add_vertex(b"ACGTA");
    let x1 = g.add_vertex(b"CCGG");
    let x2 = g.add_vertex(b"CCAGG");
    let x3 = g.add_vertex(b"CCTTGG");
    let bottom = g.add_vertex(b"TACGT");
    g.add_edge(top, x1, true, 3);
    g.add_edge(top, x2, false, 2);
    g.add_edge(top, x3, false, 1);
    g.add_edge(x1, bottom, true, 3);
    g.add_edge(x2, bottom, false, 2);
    g.add_edge(x3, bottom, false, 1);
    simplify(out, "diamond-prefix-suffix", &mut g);

    let mut g = SeqGraph::new();
    let top = g.add_vertex(b"ACGTA");
    let x1 = g.add_vertex(b"C");
    let x2 = g.add_vertex(b"G");
    let bottom = g.add_vertex(b"TACGT");
    g.add_edge(top, x1, true, 3);
    g.add_edge(top, x2, false, 2);
    g.add_edge(x1, bottom, true, 3);
    g.add_edge(x2, bottom, false, 2);
    simplify(out, "diamond-nothing-shared", &mut g);

    let mut g = SeqGraph::new();
    let a = g.add_vertex(b"AAA");
    let b = g.add_vertex(b"CCC");
    let c = g.add_vertex(b"GGG");
    let d = g.add_vertex(b"TTT");
    g.add_edge(a, b, true, 1);
    g.add_edge(b, c, true, 1);
    g.add_edge(c, d, false, 1);
    simplify(out, "zip-ref-boundary", &mut g);

    let mut g = SeqGraph::new();
    let a = g.add_vertex(b"AAA");
    let b = g.add_vertex(b"CCC");
    let c = g.add_vertex(b"GGG");
    g.add_edge(a, b, true, 1);
    g.add_edge(b, b, false, 1);
    g.add_edge(b, c, true, 1);
    simplify(out, "self-loop", &mut g);
}

/// A case: its label, the kmer size, and its reads with their counts.
type Case<'a> = (&'a str, usize, Vec<(String, usize)>);

#[test]
fn every_stage_of_every_case_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/seq_graph.txt.gz"),
    );
    let expected: Vec<&str> = golden.lines().filter(|l| !l.starts_with('#')).collect();

    let snp = format!("{}T{}", &REF[..20], &REF[21..]);
    let snp2 = format!("{}A{}", &REF[..12], &REF[13..]);
    let snp_b = format!("{}C{}", &REF[..20], &REF[21..]);
    let ins = format!("{}GG{}", &REF[..20], &REF[20..]);
    let long_ins = format!("{}TTACCGGT{}", &REF[..20], &REF[20..]);
    let del = format!("{}{}", &REF[..18], &REF[21..]);
    let mnp = format!("{}TT{}", &REF[..19], &REF[21..]);
    let before = format!("TTCCAAGGTT{}", &REF[..20]);
    let r = REF.to_string();
    let mut actual = Vec::new();
    let cases: Vec<Case> = vec![
        ("ref-only", 10, vec![]),
        ("snp", 10, vec![(r.clone(), 3), (snp.clone(), 3)]),
        ("two-snps-apart", 10, vec![(snp.clone(), 3), (snp2, 3)]),
        ("triallelic-snp", 10, vec![(snp.clone(), 3), (snp_b, 3)]),
        ("insertion", 10, vec![(ins, 3)]),
        ("long-insertion", 10, vec![(long_ins, 3)]),
        ("deletion", 10, vec![(del.clone(), 3)]),
        ("snp-and-deletion", 10, vec![(snp.clone(), 3), (del, 3)]),
        ("mnp", 10, vec![(mnp, 3)]),
        ("prefix-before-reference", 10, vec![(before, 3)]),
        ("matching-reads-only", 10, vec![(r, 5)]),
        ("k5-snp", 5, vec![(snp.clone(), 3)]),
        ("k25-snp", 25, vec![(snp, 3)]),
    ];
    for (label, k, reads) in &cases {
        from_reads(&mut actual, label, *k, reads);
    }
    by_hand(&mut actual);

    for (i, (want, got)) in expected.iter().zip(&actual).enumerate() {
        assert_eq!(got, want, "line {i}");
    }
    assert_eq!(actual.len(), expected.len(), "line count");
}
