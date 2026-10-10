//! `GraphBasedKBestHaplotypeFinder`, against the golden
//! `tools/readfilter-conformance/KBestHaplotypeDump.java` wrote in the pinned container.
//!
//! The graphs are the dump's: threaded from reads and taken through `cleanupSeqGraph`, or built by
//! hand. Each is searched from its reference source to its reference sink and from all of its
//! sources to all of its sinks, for several k; scores are compared by their bits.

use gatk_corpus as corpus;
use gatk_engine::kbest_haplotype::find_best_haplotypes;
use gatk_engine::read_threading_graph::{prune_low_weight_chains, ReadThreadingGraph};
use gatk_engine::seq_graph::SeqGraph;
use gatk_engine::tsv_table::java_double_to_string;

const REF: &str = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";
const KS: [usize; 4] = [1, 2, 3, 128];

fn from_reads(k: usize, reads: &[(String, usize)]) -> SeqGraph {
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
    seq.clean_non_ref_paths();
    seq.zip_linear_chains();
    seq.remove_singleton_orphan_vertices();
    seq.remove_vertices_not_connected_to_ref_regardless_of_edge_direction();
    seq.simplify_graph().unwrap();
    seq.remove_paths_not_connected_to_ref().unwrap();
    seq.simplify_graph().unwrap();
    if seq.vertex_count() == 1 {
        let complete = seq.vertex_ids().next().unwrap();
        let dummy = seq.add_vertex(b"");
        seq.add_edge(complete, dummy, true, 0);
    }
    seq
}

fn search(out: &mut Vec<String>, label: &str, g: &SeqGraph) {
    out.push(format!("case\t{label}"));
    let mut index = vec![usize::MAX; g.vertex_ids().max().map_or(0, |m| m + 1)];
    for (position, id) in g.vertex_ids().enumerate() {
        index[id] = position;
        out.push(format!(
            "vertex\t{label}\t{position}\t{}",
            String::from_utf8_lossy(g.sequence(id))
        ));
    }
    for e in g.edge_ids() {
        let edge = g.edge(e);
        out.push(format!(
            "edge\t{label}\t{}\t{}\t{}\t{}",
            index[edge.source], index[edge.target], edge.multiplicity, edge.is_ref
        ));
    }
    for mode in ["refpath", "all"] {
        for k in KS {
            let (sources, sinks) = if mode == "refpath" {
                (
                    g.reference_source_vertex().into_iter().collect::<Vec<_>>(),
                    g.reference_sink_vertex().into_iter().collect::<Vec<_>>(),
                )
            } else {
                (g.sources(), g.sinks())
            };
            match find_best_haplotypes(g, &sources, &sinks, k) {
                Ok(found) => {
                    for (rank, h) in found.iter().enumerate() {
                        let path: Vec<String> =
                            h.vertices.iter().map(|&v| index[v].to_string()).collect();
                        out.push(format!(
                            "haplotype\t{label}\t{mode}\t{k}\t{rank}\t{}\t{:016x},{}\t{}\t{}",
                            String::from_utf8_lossy(h.bases()),
                            h.score.to_bits(),
                            java_double_to_string(h.score),
                            h.is_reference,
                            path.join(",")
                        ));
                    }
                }
                Err(e) => out.push(format!(
                    "error\t{label}\t{mode}\t{k}\tIllegalStateException: {e:?}"
                )),
            }
        }
    }
}

fn by_hand(out: &mut Vec<String>) {
    let diamond = |m1: i32, m2: i32, ref_mult: i32| {
        let mut g = SeqGraph::new();
        let top = g.add_vertex(b"ACGT");
        let x = g.add_vertex(b"A");
        let y = g.add_vertex(b"C");
        let bottom = g.add_vertex(b"TTTT");
        g.add_edge(top, x, true, ref_mult);
        g.add_edge(top, y, false, m1);
        g.add_edge(x, bottom, true, ref_mult);
        g.add_edge(y, bottom, false, m2);
        g
    };
    search(out, "tie-on-score", &diamond(5, 5, 5));
    search(out, "zero-multiplicity", &diamond(0, 0, 5));
    search(out, "all-zero", &diamond(0, 0, 0));

    let mut g = SeqGraph::new();
    let a = g.add_vertex(b"AAAA");
    let b1 = g.add_vertex(b"C");
    let b2 = g.add_vertex(b"G");
    let m = g.add_vertex(b"TTTT");
    let c1 = g.add_vertex(b"A");
    let c2 = g.add_vertex(b"C");
    let c3 = g.add_vertex(b"G");
    let z = g.add_vertex(b"GGGG");
    g.add_edge(a, b1, true, 6);
    g.add_edge(a, b2, false, 3);
    g.add_edge(b1, m, true, 6);
    g.add_edge(b2, m, false, 3);
    g.add_edge(m, c1, true, 4);
    g.add_edge(m, c2, false, 3);
    g.add_edge(m, c3, false, 2);
    g.add_edge(c1, z, true, 4);
    g.add_edge(c2, z, false, 3);
    g.add_edge(c3, z, false, 2);
    search(out, "six-haplotypes", &g);

    let mut g = SeqGraph::new();
    let s = g.add_vertex(b"AAAA");
    let p = g.add_vertex(b"CC");
    let q = g.add_vertex(b"GG");
    let t = g.add_vertex(b"TTTT");
    g.add_edge(s, p, true, 4);
    g.add_edge(p, q, true, 4);
    g.add_edge(q, p, false, 1);
    g.add_edge(s, q, false, 2);
    g.add_edge(q, t, true, 4);
    search(out, "cycle", &g);

    let mut g = SeqGraph::new();
    let s1 = g.add_vertex(b"AAAA");
    let s2 = g.add_vertex(b"CCCC");
    let mid = g.add_vertex(b"GG");
    let t1 = g.add_vertex(b"TTTT");
    let t2 = g.add_vertex(b"TATA");
    g.add_edge(s1, mid, true, 4);
    g.add_edge(s2, mid, false, 2);
    g.add_edge(mid, t1, true, 4);
    g.add_edge(mid, t2, false, 2);
    search(out, "two-sources-two-sinks", &g);
}

/// The sign bit of a NaN is the FPU's choice, not the algorithm's: `-Infinity - -Infinity` is
/// `0xfff8...` on x86-64, where the golden was made, and `0x7ff8...` on Apple silicon (htsjdk-rs
/// decision 0032, "the sign of a NaN is not a measurement"). Both are compared as NaN.
fn nan_sign_masked(line: &str) -> String {
    line.replace("fff8000000000000,NaN", "NaN,NaN")
        .replace("7ff8000000000000,NaN", "NaN,NaN")
}

/// A case: its label, the kmer size, and its reads with their counts.
type Case<'a> = (&'a str, usize, Vec<(String, usize)>);

#[test]
fn the_k_best_haplotypes_match_the_reference_bit_for_bit() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/kbest_haplotype.txt.gz"),
    );
    let expected: Vec<&str> = golden.lines().filter(|l| !l.starts_with('#')).collect();

    let snp = format!("{}T{}", &REF[..20], &REF[21..]);
    let snp2 = format!("{}A{}", &REF[..12], &REF[13..]);
    let snp_b = format!("{}C{}", &REF[..20], &REF[21..]);
    let ins = format!("{}GG{}", &REF[..20], &REF[20..]);
    let del = format!("{}{}", &REF[..18], &REF[21..]);
    let r = REF.to_string();
    let cases: Vec<Case> = vec![
        ("ref-only", 10, vec![]),
        ("snp", 10, vec![(r.clone(), 3), (snp.clone(), 3)]),
        ("snp-minor", 10, vec![(r, 9), (snp.clone(), 2)]),
        ("two-snps-apart", 10, vec![(snp.clone(), 3), (snp2, 4)]),
        ("triallelic-snp", 10, vec![(snp.clone(), 3), (snp_b, 3)]),
        ("insertion", 10, vec![(ins, 3)]),
        ("snp-and-deletion", 10, vec![(snp.clone(), 3), (del, 5)]),
        ("k5-snp", 5, vec![(snp, 3)]),
    ];
    let mut actual = Vec::new();
    for (label, k, reads) in &cases {
        search(&mut actual, label, &from_reads(*k, reads));
    }
    by_hand(&mut actual);
    for (i, (want, got)) in expected.iter().zip(&actual).enumerate() {
        assert_eq!(nan_sign_masked(got), nan_sign_masked(want), "line {i}");
    }
    assert_eq!(actual.len(), expected.len(), "line count");
}
