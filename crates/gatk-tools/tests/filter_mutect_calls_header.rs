//! The header half of the `filter-mutect-calls` conformance, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/FilterMutectCallsDump.java`.
//!
//! The tool rewrites the header rather than appending to it: Mutect2's `filtering_status` line is
//! dropped and replaced under the same key, every `##FILTER` in `MUTECT_FILTER_NAMES` is added
//! whether or not its filter runs, `AS_FilterStatus` and `STRQ` arrive as `##INFO`, and the input's
//! own `##FILTER` lines survive. `VCFHeader` writes its metadata sorted, so the output's order is the
//! strings' order and not the order they were added in.

use gatk_corpus as corpus;
use gatk_tools::filter_mutect_calls::output_header_lines;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/filter_mutect_calls.txt.gz"),
    )
}

/// The one `##FILTER` line the input header carries, which the output keeps.
const INPUT_FILTER_LINES: [&str; 1] = ["##FILTER=<ID=LowQual,Description=\"Low quality\">"];

#[test]
fn every_header_row_matches_the_golden() {
    let text = golden();
    let mut by_run: std::collections::BTreeMap<String, Vec<String>> =
        std::collections::BTreeMap::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.splitn(3, '\t');
        if fields.next() != Some("header") {
            continue;
        }
        let run = fields.next().expect("a run").to_string();
        // The dump escapes the line; the header lines carry no escape-worthy character but the
        // quotes, so the payload is the line as written.
        by_run
            .entry(run)
            .or_default()
            .push(fields.next().expect("a line").to_string());
    }
    assert_eq!(by_run.len(), 5, "the runs that produced a header");

    let inputs: Vec<String> = INPUT_FILTER_LINES.iter().map(|l| l.to_string()).collect();
    let ours = output_header_lines(&inputs);
    for (run, expected) in &by_run {
        assert_eq!(&ours, expected, "the header of the {run} run");
    }
}
