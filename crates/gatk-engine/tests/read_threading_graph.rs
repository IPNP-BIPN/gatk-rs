//! `ReadThreadingGraph` after `buildGraphIfNecessary`, against the golden
//! `tools/readfilter-conformance/ReadThreadingGraphDump.java` wrote in the pinned container.
//!
//! The cases are the dump's, built the same way: the reference through `addSequence`, then each
//! read through `addRead` under its sample, then the build. The comparison is line for line, so
//! the order of vertices, edges and k-mer map entries is compared as well as their contents: that
//! order is what the threading loop and every later stage of the assembler read.

use gatk_corpus as corpus;
use gatk_engine::read_threading_graph::{GraphError, ReadThreadingGraph};

const REF: &str = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";

struct Read {
    sample: &'static str,
    bases: String,
    quals: Option<Vec<u8>>,
}

struct Case {
    label: &'static str,
    k: usize,
    min_base_quality: u8,
    pruning_samples: usize,
    only_at_existing: bool,
    through_branches: bool,
    reference: String,
    reads: Vec<Read>,
}

fn read(sample: &'static str, bases: &str) -> Read {
    Read {
        sample,
        bases: bases.to_string(),
        quals: None,
    }
}

fn quals(length: usize, value: u8, low: &[usize]) -> Vec<u8> {
    let mut q = vec![value; length];
    for &i in low {
        q[i] = 2;
    }
    q
}

#[allow(clippy::too_many_arguments)]
fn case(
    label: &'static str,
    k: usize,
    min_base_quality: u8,
    pruning_samples: usize,
    only_at_existing: bool,
    through_branches: bool,
    reference: &str,
    reads: Vec<Read>,
) -> Case {
    Case {
        label,
        k,
        min_base_quality,
        pruning_samples,
        only_at_existing,
        through_branches,
        reference: reference.to_string(),
        reads,
    }
}

fn cases() -> Vec<Case> {
    let snp = format!("{}T{}", &REF[..20], &REF[21..]);
    let ins = format!("{}GG{}", &REF[..20], &REF[20..]);
    let del = format!("{}{}", &REF[..18], &REF[21..]);
    let repeat = "ACGTACGTACGTACGTTGCATGTCGCATG";
    let low = Read {
        sample: "s1",
        bases: REF.to_string(),
        quals: Some(quals(REF.len(), 30, &[3, 20, 21])),
    };
    let low_again = Read {
        sample: "s1",
        bases: REF.to_string(),
        quals: Some(quals(REF.len(), 30, &[3, 20, 21])),
    };
    vec![
        case("ref-only", 5, 10, 1, false, false, REF, vec![]),
        case(
            "ref-and-matching-reads",
            5,
            10,
            1,
            false,
            false,
            REF,
            vec![
                read("s1", REF),
                read("s1", &REF[5..35]),
                read("s1", &REF[10..]),
            ],
        ),
        case(
            "snp-bubble",
            5,
            10,
            1,
            false,
            false,
            REF,
            vec![read("s1", &snp), read("s1", &snp)],
        ),
        case(
            "insertion-bubble",
            5,
            10,
            1,
            false,
            false,
            REF,
            vec![read("s1", &ins)],
        ),
        case(
            "deletion-bubble",
            5,
            10,
            1,
            false,
            false,
            REF,
            vec![read("s1", &del)],
        ),
        case(
            "read-split-at-n",
            5,
            10,
            1,
            false,
            false,
            REF,
            vec![read("s1", &format!("{}N{}", &REF[..15], &REF[16..]))],
        ),
        case(
            "read-split-at-low-quality",
            5,
            10,
            1,
            false,
            false,
            REF,
            vec![low],
        ),
        case(
            "min-quality-zero",
            5,
            0,
            1,
            false,
            false,
            REF,
            vec![low_again],
        ),
        case(
            "repeat-in-reference",
            4,
            10,
            1,
            false,
            false,
            repeat,
            vec![read("s1", &repeat[4..])],
        ),
        case(
            "low-complexity",
            3,
            10,
            1,
            false,
            false,
            "AAAAAAAAAACAAAAAAAAAA",
            vec![read("s1", "AAAAAAAACAAAAAAA")],
        ),
        case(
            "read-without-a-start",
            3,
            10,
            1,
            false,
            false,
            REF,
            vec![read("s1", "ATATATATATAT")],
        ),
        case(
            "read-before-reference",
            5,
            10,
            1,
            false,
            false,
            &REF[10..],
            vec![read("s1", &REF[..25])],
        ),
        case(
            "rearranged-read",
            5,
            10,
            1,
            false,
            false,
            REF,
            vec![read("s1", &format!("{}{}", &REF[25..], &REF[..15]))],
        ),
        case(
            "start-only-at-existing",
            5,
            10,
            1,
            true,
            false,
            &REF[10..],
            vec![read("s1", &REF[..25])],
        ),
        case(
            "two-samples",
            5,
            10,
            2,
            false,
            false,
            REF,
            vec![
                read("s1", &snp),
                read("s1", &snp),
                read("s2", &snp),
                read("s2", REF),
            ],
        ),
        case(
            "two-samples-one-pruning-sample",
            5,
            10,
            1,
            false,
            false,
            REF,
            vec![
                read("s1", &snp),
                read("s1", &snp),
                read("s2", &snp),
                read("s2", REF),
            ],
        ),
        case(
            "counts-through-branches",
            5,
            10,
            1,
            false,
            true,
            REF,
            vec![read("s1", &snp), read("s1", &snp[12..])],
        ),
        case(
            "counts-stop-at-branches",
            5,
            10,
            1,
            false,
            false,
            REF,
            vec![read("s1", &snp), read("s1", &snp[12..])],
        ),
        case(
            "k10-snp",
            10,
            10,
            1,
            false,
            false,
            REF,
            vec![read("s1", &snp), read("s1", REF)],
        ),
        case(
            "k25-snp",
            25,
            10,
            1,
            false,
            false,
            REF,
            vec![read("s1", &snp)],
        ),
    ]
}

fn error_line(label: &str, error: &GraphError) -> String {
    format!("error\t{label}\t{error:?}")
}

fn render(c: &Case) -> Vec<String> {
    let mut out = vec![format!(
        "case\t{}\tk={}\tminbq={}\tpruning={}\texisting={}\tbranches={}",
        c.label, c.k, c.min_base_quality, c.pruning_samples, c.only_at_existing, c.through_branches
    )];
    let mut graph = ReadThreadingGraph::new(c.k, c.min_base_quality, c.pruning_samples);
    graph.set_threading_start_only_at_existing_vertex(c.only_at_existing);
    graph.set_increase_counts_through_branches(c.through_branches);
    let built = graph
        .add_sequence(c.reference.as_bytes(), true)
        .and_then(|_| {
            for r in &c.reads {
                let q = r.quals.clone().unwrap_or_else(|| vec![30; r.bases.len()]);
                graph.add_read(r.sample, r.bases.as_bytes(), &q)?;
            }
            Ok(())
        })
        .and_then(|_| graph.build_graph_if_necessary());
    if let Err(error) = built {
        out.push(error_line(c.label, &error));
        return out;
    }

    let label = c.label;
    for (i, v) in graph.vertices().iter().enumerate() {
        out.push(format!(
            "vertex\t{label}\t{i}\t{}\t{}",
            String::from_utf8_lossy(&v.sequence),
            v.additional_info
        ));
    }
    for e in graph.edges() {
        out.push(format!(
            "edge\t{label}\t{}\t{}\t{}\t{}\t{}",
            e.source,
            e.target,
            e.multiplicity,
            e.pruning_multiplicity(),
            e.is_ref
        ));
    }
    let path: Vec<String> = graph
        .reference_path()
        .unwrap_or(&[])
        .iter()
        .map(|v| v.to_string())
        .collect();
    out.push(format!("refpath\t{label}\t{}", path.join(",")));
    for (kmer, vertex) in graph.kmer_to_vertex() {
        out.push(format!(
            "kmer\t{label}\t{}\t{vertex}",
            String::from_utf8_lossy(kmer)
        ));
    }
    let mut non_unique: Vec<String> = graph
        .non_unique_kmers()
        .iter()
        .map(|k| String::from_utf8_lossy(k).into_owned())
        .collect();
    non_unique.sort();
    out.push(format!("nonunique\t{label}\t{}", non_unique.join(",")));
    let show = |v: Option<usize>| v.map_or("-".to_string(), |v| v.to_string());
    out.push(format!(
        "summary\t{label}\tvertices={}\tedges={}\tcycles={}\tlowquality={}\trefsource={}\trefsink={}",
        graph.vertices().len(),
        graph.edges().len(),
        graph.has_cycles(),
        graph.is_low_quality_graph(),
        show(graph.reference_source_vertex()),
        show(graph.reference_sink_vertex())
    ));
    out
}

#[test]
fn every_case_builds_the_reference_graph_line_for_line() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/read_threading_graph.txt.gz"),
    );
    let expected: Vec<&str> = golden.lines().filter(|l| !l.starts_with('#')).collect();
    let actual: Vec<String> = cases().iter().flat_map(render).collect();
    for (i, (want, got)) in expected.iter().zip(&actual).enumerate() {
        assert_eq!(got, want, "line {i}");
    }
    assert_eq!(actual.len(), expected.len(), "line count");
}
