//! Dangling-end recovery on the read threading graph, against the golden
//! `tools/readfilter-conformance/DanglingEndDump.java` wrote in the pinned container.
//!
//! Each case is built and chain-pruned as the dump does, then printed after the tails are
//! recovered, after the heads are, and after the removal of paths not connected to the reference.
//! A vertex's index is its position in the current vertex set, as in the dump.

use gatk_corpus as corpus;
use gatk_engine::read_threading_graph::{prune_low_weight_chains, GraphError, ReadThreadingGraph};
use gatk_engine::smith_waterman::{SmithWatermanJavaAligner, STANDARD_NGS};

const REF: &str = "ATGCGTACCTGAGTCAAGCTTGGACTCTAGAGCAATCGGTACCATTGACGTAGCTAGGCA";

struct Case {
    label: &'static str,
    k: usize,
    prune_factor: i32,
    min_branch: i32,
    recover_all: bool,
    min_matching: i32,
    reads: Vec<(&'static str, String, usize)>,
}

/// `DanglingEndDump.alt`: a base that is not the reference's.
fn alt(c: u8) -> char {
    match c {
        b'A' => 'C',
        b'C' => 'G',
        b'G' => 'T',
        _ => 'A',
    }
}

fn mutate(s: &str, at: usize) -> String {
    format!("{}{}{}", &s[..at], alt(s.as_bytes()[at]), &s[at + 1..])
}

fn case(
    label: &'static str,
    k: usize,
    prune_factor: i32,
    min_branch: i32,
    recover_all: bool,
    min_matching: i32,
    reads: Vec<(&'static str, String, usize)>,
) -> Case {
    Case {
        label,
        k,
        prune_factor,
        min_branch,
        recover_all,
        min_matching,
        reads,
    }
}

/// `DanglingEndDump.tails`: the defaults, k of 10.
fn tails(label: &'static str, reads: Vec<(&'static str, String, usize)>) -> Case {
    case(label, 10, 2, 4, false, -1, reads)
}

fn one(bases: &str) -> Vec<(&'static str, String, usize)> {
    vec![("s1", bases.to_string(), 3)]
}

fn cases() -> Vec<Case> {
    let r = REF;
    let b = r.as_bytes();
    let tail_snp4 = format!("{}{}{}", &r[..40], alt(b[40]), &r[41..45]);
    let tail_snp1 = format!("{}{}{}", &r[..40], alt(b[40]), &r[41..42]);
    let tail_snp8 = format!("{}{}{}", &r[..40], alt(b[40]), &r[41..49]);
    let tail_ins = format!("{}GA{}", &r[..40], &r[40..46]);
    let tail_del = format!("{}{}", &r[..40], &r[43..49]);
    let tail_mnp = format!("{}{}{}{}", &r[..40], alt(b[40]), alt(b[41]), &r[42..47]);
    let tail_junk = format!("{}TTTTTTTT", &r[..40]);
    let tail_other = format!("{}{}GGGG", &r[..40], alt(b[40]));
    let tail_near_start = format!("{}{}{}", &r[..12], alt(b[12]), &r[13..17]);
    let tail_at_source = format!("{}{}{}", &r[..10], alt(b[10]), &r[11..16]);
    let head2 = mutate(&r[20..], 2);
    let head4 = mutate(&r[20..], 4);
    let head6 = mutate(&r[20..], 6);
    let head9 = mutate(&r[20..], 9);
    let head_ins = format!("{}GA{}", &r[20..26], &r[26..]);
    let head_del = format!("{}{}", &r[20..26], &r[29..]);
    let head_junk = format!("TTTTTTTT{}", &r[28..]);
    let head_long = mutate(&r[5..], 3);
    let pair = |a: &str, an: usize, b: &str, bn: usize| {
        vec![("s1", a.to_string(), an), ("s1", b.to_string(), bn)]
    };

    vec![
        case("ref-only", 10, 2, 4, false, -1, vec![]),
        tails("tail-snp-4", one(&tail_snp4)),
        tails("tail-snp-1", one(&tail_snp1)),
        case(
            "tail-snp-1-minbranch-0",
            10,
            2,
            0,
            false,
            -1,
            one(&tail_snp1),
        ),
        tails("tail-snp-8", one(&tail_snp8)),
        tails("tail-insertion", one(&tail_ins)),
        tails("tail-deletion", one(&tail_del)),
        tails("tail-mnp", one(&tail_mnp)),
        tails("tail-junk", one(&tail_junk)),
        case("tail-junk-recover-all", 10, 2, 4, true, -1, one(&tail_junk)),
        tails("tail-near-start", one(&tail_near_start)),
        tails("tail-lca-is-ref-source", one(&tail_at_source)),
        case(
            "tail-snp-4-minmatching-3",
            10,
            2,
            4,
            false,
            3,
            one(&tail_snp4),
        ),
        case(
            "tail-snp-4-minmatching-5",
            10,
            2,
            4,
            false,
            5,
            one(&tail_snp4),
        ),
        case(
            "tail-deletion-minmatching-2",
            10,
            2,
            4,
            false,
            2,
            one(&tail_del),
        ),
        tails("tail-fork", pair(&tail_snp4, 3, &tail_other, 3)),
        case(
            "tail-fork-recover-all",
            10,
            2,
            4,
            true,
            -1,
            pair(&tail_snp4, 3, &tail_other, 3),
        ),
        tails("tail-light-end", pair(&tail_snp4, 3, &tail_snp8, 1)),
        case(
            "tail-light-end-factor-0",
            10,
            0,
            4,
            false,
            -1,
            pair(&tail_snp4, 3, &tail_snp8, 1),
        ),
        tails("head-2", one(&head2)),
        tails("head-4", one(&head4)),
        tails("head-6", one(&head6)),
        tails("head-9", one(&head9)),
        tails("head-insertion", one(&head_ins)),
        tails("head-deletion", one(&head_del)),
        tails("head-junk", one(&head_junk)),
        tails("head-long", one(&head_long)),
        case("head-2-minbranch-0", 10, 2, 0, false, -1, one(&head2)),
        case("head-4-minmatching-3", 10, 2, 4, false, 3, one(&head4)),
        case("head-6-minmatching-3", 10, 2, 4, false, 3, one(&head6)),
        case("head-9-minmatching-3", 10, 2, 4, false, 3, one(&head9)),
        case("head-9-minmatching-0", 10, 2, 4, false, 0, one(&head9)),
        case(
            "head-insertion-minmatching-3",
            10,
            2,
            4,
            false,
            3,
            one(&head_ins),
        ),
        case(
            "head-deletion-minmatching-3",
            10,
            2,
            4,
            false,
            3,
            one(&head_del),
        ),
        case("head-junk-recover-all", 10, 2, 4, true, -1, one(&head_junk)),
        case("head-6-recover-all", 10, 2, 4, true, -1, one(&head6)),
        tails("head-and-tail", pair(&head6, 3, &tail_snp4, 3)),
        case(
            "head-and-tail-minmatching-3",
            10,
            2,
            4,
            false,
            3,
            pair(&head6, 3, &tail_snp4, 3),
        ),
        tails(
            "two-samples",
            vec![("s1", head9.clone(), 2), ("s2", tail_ins.clone(), 2)],
        ),
        case("k5-tail-snp-4", 5, 2, 4, false, -1, one(&tail_snp4)),
        case("k5-head-6", 5, 2, 4, false, -1, one(&head6)),
        case(
            "negative-prune-factor",
            10,
            -1,
            4,
            false,
            -1,
            one(&tail_snp4),
        ),
        case("negative-min-branch", 10, 2, -1, false, -1, one(&tail_snp4)),
    ]
}

fn print(out: &mut Vec<String>, label: &str, stage: &str, graph: &ReadThreadingGraph) {
    let mut index = vec![None; graph.vertex_ids().max().map_or(0, |m| m + 1)];
    for (position, id) in graph.vertex_ids().enumerate() {
        index[id] = Some(position);
    }
    let at = |id: usize| index.get(id).copied().flatten();
    for id in graph.vertex_ids() {
        out.push(format!(
            "vertex\t{label}\t{stage}\t{}\t{}",
            at(id).unwrap(),
            String::from_utf8_lossy(&graph.vertex(id).sequence)
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

/// The simple name of the exception the reference throws.
fn describe(error: &GraphError) -> &'static str {
    match error {
        GraphError::IllegalArgument(_) => "IllegalArgumentException",
        GraphError::IndexOutOfBounds => "IndexOutOfBoundsException",
        GraphError::NoSuchElement => "NoSuchElementException",
        _ => "IllegalStateException",
    }
}

fn render(c: &Case) -> Vec<String> {
    let mut out = vec![format!(
        "case\t{}\tk={}\tprunefactor={}\tminbranch={}\trecoverall={}\tminmatching={}",
        c.label, c.k, c.prune_factor, c.min_branch, c.recover_all, c.min_matching
    )];
    let mut graph = ReadThreadingGraph::new(c.k, 10, 1);
    graph.set_min_matching_bases_to_dangling_end_recovery(c.min_matching);
    graph.add_sequence(REF.as_bytes(), true).unwrap();
    for (sample, bases, copies) in &c.reads {
        for _ in 0..*copies {
            graph
                .add_read(sample, bases.as_bytes(), &vec![30; bases.len()])
                .unwrap();
        }
    }
    graph.build_graph_if_necessary().unwrap();
    prune_low_weight_chains(&mut graph, c.prune_factor.max(0));
    print(&mut out, c.label, "pruned", &graph);

    let aligner = SmithWatermanJavaAligner;
    let stages = (|| {
        graph.recover_dangling_tails(
            c.prune_factor,
            c.min_branch,
            c.recover_all,
            &aligner,
            &STANDARD_NGS,
        )?;
        print(&mut out, c.label, "tails", &graph);
        graph.recover_dangling_heads(
            c.prune_factor,
            c.min_branch,
            c.recover_all,
            &aligner,
            &STANDARD_NGS,
        )?;
        print(&mut out, c.label, "heads", &graph);
        graph.remove_paths_not_connected_to_ref()?;
        print(&mut out, c.label, "connected", &graph);
        Ok::<(), GraphError>(())
    })();
    if let Err(error) = stages {
        out.push(format!("error\t{}\t{}", c.label, describe(&error)));
    }
    out
}

#[test]
fn every_case_recovers_dangling_ends_as_the_reference_does() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dangling_ends.txt.gz"),
    );
    let expected: Vec<&str> = golden.lines().filter(|l| !l.starts_with('#')).collect();
    let actual: Vec<String> = cases().iter().flat_map(render).collect();
    for (i, (want, got)) in expected.iter().zip(&actual).enumerate() {
        assert_eq!(got, want, "line {i}");
    }
    assert_eq!(actual.len(), expected.len(), "line count");
}
