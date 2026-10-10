//! `Mutect2Engine.logLikelihoodRatio` and `AdaptiveChainPruner` on the read threading graph,
//! against the golden `tools/readfilter-conformance/AdaptiveChainPrunerDump.java` wrote in the
//! pinned container.
//!
//! The log likelihood ratios and every chain's log odds are compared by their bits, because the
//! pruner compares them against thresholds and an ulp can move a chain across one.

use gatk_corpus as corpus;
use gatk_engine::mutect_engine::log_likelihood_ratio_from_error;
use gatk_engine::read_threading_graph::{
    chain_log_odds, find_all_chains, prune_adaptive, AdaptivePruning, ReadThreadingGraph,
};
use gatk_engine::tsv_table::java_double_to_string;

const REF: &str = "ACGTTGCATGTCGCATGATGCATGAGAGCTCAGTCTAGGC";

fn render(value: f64) -> String {
    format!("{:016x},{}", value.to_bits(), java_double_to_string(value))
}

/// `ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD` and the seeding
/// one: `MathUtils.log10ToLog(1.0)` and `log10ToLog(4.0)`, a factor of the natural logarithm of ten.
fn lod() -> f64 {
    jmath::math::log(10.0)
}

fn seed() -> f64 {
    4.0 * jmath::math::log(10.0)
}

struct Case {
    label: &'static str,
    k: usize,
    params: AdaptivePruning,
    reads: Vec<(String, usize)>,
}

fn params(initial: f64, threshold: f64, seeding: f64, max: i32) -> AdaptivePruning {
    AdaptivePruning {
        initial_error_probability: initial,
        log_odds_threshold: threshold,
        seeding_log_odds_threshold: seeding,
        max_unpruned_variants: max,
    }
}

fn cases() -> Vec<Case> {
    let snp = format!("{}T{}", &REF[..20], &REF[21..]);
    let snp2 = format!("{}A{}", &REF[..12], &REF[13..]);
    let tail = format!("{}TTTTTTTT", &REF[..32]);
    let r = REF.to_string();
    let d = params(0.001, lod(), seed(), 100);
    let c = |label, k, params, reads| Case {
        label,
        k,
        params,
        reads,
    };
    vec![
        c("real-bubble", 10, d, vec![(r.clone(), 8), (snp.clone(), 2)]),
        c(
            "error-bubble",
            10,
            d,
            vec![(r.clone(), 20), (snp.clone(), 1)],
        ),
        c(
            "two-bubbles",
            10,
            d,
            vec![(r.clone(), 10), (snp.clone(), 5), (snp2.clone(), 1)],
        ),
        c(
            "dangling-tail",
            10,
            d,
            vec![(r.clone(), 10), (tail.clone(), 1)],
        ),
        c(
            "heavy-dangling-tail",
            10,
            d,
            vec![(r.clone(), 4), (tail, 4)],
        ),
        c(
            "max-variants-zero",
            10,
            params(0.001, lod(), seed(), 0),
            vec![(r.clone(), 8), (snp.clone(), 4), (snp2.clone(), 4)],
        ),
        c(
            "max-variants-one",
            10,
            params(0.001, lod(), seed(), 1),
            vec![(r.clone(), 8), (snp.clone(), 4), (snp2, 4)],
        ),
        c(
            "high-threshold",
            10,
            params(0.001, 50.0, 60.0, 100),
            vec![(r.clone(), 8), (snp.clone(), 3)],
        ),
        c(
            "zero-threshold",
            10,
            params(0.001, 0.0, 0.0, 100),
            vec![(r.clone(), 8), (snp.clone(), 1)],
        ),
        c(
            "initial-error-rate-high",
            10,
            params(0.2, lod(), seed(), 100),
            vec![(r.clone(), 8), (snp.clone(), 2)],
        ),
        c("ref-only", 10, d, vec![]),
        c("k5-real-bubble", 5, d, vec![(r, 6), (snp, 3)]),
    ]
}

fn llr_lines() -> Vec<String> {
    let mut out = Vec::new();
    for reference in [0, 1, 2, 5, 10, 30, 100] {
        for alt in [1, 2, 3, 10, 40] {
            for (p, shown) in [
                (0.001, "0.001"),
                (0.01, "0.01"),
                (0.05, "0.05"),
                (0.2, "0.2"),
                (0.5, "0.5"),
                (0.0, "0.0"),
                (1.0, "1.0"),
                (f64::NAN, "NaN"),
                (1.5, "1.5"),
            ] {
                let value = log_likelihood_ratio_from_error(reference, alt, p)
                    .map_or("exception:IllegalArgumentException".to_string(), render);
                out.push(format!("llr\t{reference}\t{alt}\t{shown}\t{value}"));
            }
        }
    }
    out
}

fn positions(graph: &ReadThreadingGraph) -> Vec<Option<usize>> {
    let mut index = vec![None; graph.vertex_ids().max().map_or(0, |m| m + 1)];
    for (position, id) in graph.vertex_ids().enumerate() {
        index[id] = Some(position);
    }
    index
}

fn render_case(c: &Case) -> Vec<String> {
    let p = &c.params;
    let mut out = vec![format!(
        "case\t{}\tk={}\tinitial={}\tthreshold={}\tseeding={}\tmaxvariants={}",
        c.label,
        c.k,
        java_double_to_string(p.initial_error_probability),
        java_double_to_string(p.log_odds_threshold),
        java_double_to_string(p.seeding_log_odds_threshold),
        p.max_unpruned_variants
    )];
    let mut graph = ReadThreadingGraph::new(c.k, 10, 1);
    graph.add_sequence(REF.as_bytes(), true).unwrap();
    for (bases, copies) in &c.reads {
        for _ in 0..*copies {
            graph
                .add_read("s1", bases.as_bytes(), &vec![30; bases.len()])
                .unwrap();
        }
    }
    graph.build_graph_if_necessary().unwrap();

    let chains = find_all_chains(&graph);
    if !chains.is_empty() {
        let index = positions(&graph);
        let first = gatk_engine::read_threading_graph::likely_error_chains(
            &graph,
            &chains,
            p.initial_error_probability,
            p,
        )
        .unwrap();
        let pass = |out: &mut Vec<String>, name: &str, rate: f64, errors: &[usize]| {
            for (i, chain) in chains.iter().enumerate() {
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
                let (left, right) = chain_log_odds(&graph, chain, rate).unwrap();
                out.push(format!(
                    "chain\t{}\t{name}\t{}\t{}\t{}\t{}\t{}",
                    c.label,
                    vertices.join(","),
                    multiplicities.join(","),
                    render(left),
                    render(right),
                    errors.contains(&i)
                ));
            }
        };
        pass(&mut out, "initial", p.initial_error_probability, &first);
        let error_count: i32 = first
            .iter()
            .map(|&i| graph.edge(*chains[i].edges.last().unwrap()).multiplicity)
            .sum();
        let total: i32 = chains
            .iter()
            .map(|chain| {
                chain
                    .edges
                    .iter()
                    .map(|&e| graph.edge(e).multiplicity)
                    .sum::<i32>()
            })
            .sum();
        let rate = f64::from(error_count) / f64::from(total);
        out.push(format!("errorrate\t{}\t{}", c.label, render(rate)));
        let second =
            gatk_engine::read_threading_graph::likely_error_chains(&graph, &chains, rate, p)
                .unwrap();
        pass(&mut out, "estimated", rate, &second);
    }

    prune_adaptive(&mut graph, p).unwrap();
    let index = positions(&graph);
    for id in graph.vertex_ids() {
        out.push(format!(
            "vertex\t{}\tpruned\t{}\t{}",
            c.label,
            index[id].unwrap(),
            String::from_utf8_lossy(&graph.vertex(id).sequence)
        ));
    }
    for id in graph.edge_ids() {
        let e = graph.edge(id);
        out.push(format!(
            "edge\t{}\tpruned\t{}\t{}\t{}\t{}",
            c.label,
            index[e.source].unwrap(),
            index[e.target].unwrap(),
            e.multiplicity,
            e.is_ref
        ));
    }
    out.push(format!(
        "summary\t{}\tpruned\tvertices={}\tedges={}",
        c.label,
        graph.vertex_count(),
        graph.edge_count()
    ));
    out
}

#[test]
fn the_ratio_and_every_chain_decision_match_the_reference_bit_for_bit() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/adaptive_chain_pruner.txt.gz"),
    );
    let expected: Vec<&str> = golden.lines().filter(|l| !l.starts_with('#')).collect();
    let mut actual = llr_lines();
    for c in cases() {
        actual.extend(render_case(&c));
    }
    for (i, (want, got)) in expected.iter().zip(&actual).enumerate() {
        assert_eq!(got, want, "line {i}");
    }
    assert_eq!(actual.len(), expected.len(), "line count");
}
