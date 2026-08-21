//! Conformance for `ReferenceBlockConcordance` against GATK 4.6.2.0, compared as the three
//! metrics files of every run.
//!
//! Golden from `tools/readfilter-conformance/ReferenceBlockConcordanceDump.java`, which carries
//! each run's two GVCFs and its three outputs.
//!
//! # What this suite is for
//!
//!  * **the three histograms**, bins, order and values;
//!  * **the string sort**, so `1,80` precedes `100,20` and `50,40` comes last;
//!  * **the concordance histogram counting bases**, 50, 20 and 1 for the three overlapping pairs;
//!  * **a filtered block being walked** and a variant site not;
//!  * **a block on one side only** still reaching that side's histogram;
//!  * **and the multi-sample refusal**, which the hom-ref filter lets through and the length
//!    extraction raises.
//!
//! The walk itself is `gatk_engine::concordance_walker`, already oracle-backed under
//! `concordance-walker`; this suite drives it and compares what the tool builds on top.

use gatk_corpus as corpus;
use gatk_engine::concordance_walker::{concordance, ConcordanceRecord};
use gatk_tools::reference_block_concordance::{accumulate, write_histogram, Block, MultiSample};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/reference_block_concordance.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

/// The `histogram\t<label>\t<which>=` rows.
fn histogram(text: &str, label: &str, which: &str) -> String {
    let prefix = format!("histogram\t{label}\t{which}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries histogram/{label}/{which}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

/// The blocks of one GVCF, filtered the way the tool filters: hom-ref on genotype 0.
fn blocks(gvcf: &str) -> Vec<Block> {
    let mut found = Vec::new();
    for line in gvcf.lines() {
        if line.starts_with('#') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        let start: i32 = columns[1].parse().expect("a position");
        let end = columns[7]
            .split(';')
            .find_map(|field| field.strip_prefix("END="))
            .and_then(|value| value.parse().ok())
            .unwrap_or(start);
        let format: Vec<&str> = columns[8].split(':').collect();
        let genotypes = columns.len() - 9;
        let first: Vec<&str> = columns[9].split(':').collect();
        let call = first[format.iter().position(|f| *f == "GT").expect("a GT")];
        let gq: i32 = first[format.iter().position(|f| *f == "GQ").expect("a GQ")]
            .parse()
            .expect("a GQ");
        let is_hom_ref = call.split(['/', '|']).all(|allele| allele == "0");
        if !is_hom_ref {
            continue;
        }
        found.push(Block {
            contig: columns[0].to_string(),
            start,
            end,
            gq,
            is_hom_ref,
            genotypes,
            rendered: String::new(),
        });
    }
    found
}

/// The blocks as the walker reads them, which needs only three fields.
struct Step(Block);

impl ConcordanceRecord for Step {
    fn contig(&self) -> &str {
        &self.0.contig
    }
    fn start(&self) -> i32 {
        self.0.start
    }
    fn is_filtered(&self) -> bool {
        // The tool's filters read the genotype alone, so nothing here is ever "filtered" as far
        // as the walker is concerned.
        false
    }
}

fn run(text: &str, label: &str) -> (String, String, String) {
    let truth = blocks(&value(text, "truth", label));
    let eval = blocks(&value(text, "eval", label));
    let truth_steps: Vec<Step> = truth.iter().cloned().map(Step).collect();
    let eval_steps: Vec<Step> = eval.iter().cloned().map(Step).collect();
    // `areVariantsAtSameLocusConcordant` is `true` for this tool.
    let steps = concordance(&truth_steps, &eval_steps, &["chr1".to_string()], |_, _| {
        true
    });
    let pairs: Vec<(Option<usize>, Option<usize>)> =
        steps.iter().map(|step| (step.truth, step.eval)).collect();
    let histograms = accumulate(&truth, &eval, &pairs).expect("a run the tool allows");
    // The two headers the engine writes, which the golden masks.
    let headers = vec!["MASKED".to_string(), "MASKED".to_string()];
    (
        write_histogram(&histograms.truth_blocks, &headers),
        write_histogram(&histograms.eval_blocks, &headers),
        write_histogram(&histograms.confidence_concordance, &headers),
    )
}

/// The golden's files carry the masked header lines; this rebuilds the same masking so the two
/// can be compared as whole files.
fn masked(file: &str) -> String {
    file.lines()
        .map(|line| {
            if line.starts_with("# ") {
                "# MASKED"
            } else {
                line
            }
        })
        .collect::<Vec<&str>>()
        .join("\n")
        + "\n"
}

#[test]
fn every_histogram_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in ["blocks", "identical", "empty-truth"] {
        let (truth, eval, concordance) = run(&text, label);
        assert_eq!(
            masked(&truth),
            histogram(&text, label, "truth-blocks"),
            "{label}: truth blocks"
        );
        assert_eq!(
            masked(&eval),
            histogram(&text, label, "eval-blocks"),
            "{label}: eval blocks"
        );
        assert_eq!(
            masked(&concordance),
            histogram(&text, label, "concordance"),
            "{label}: concordance"
        );
        compared += 1;
    }
    assert_eq!(compared, 3, "the golden's runs");
}

/// The bins are sorted as strings, which is neither the lengths' order nor the file's.
#[test]
fn the_bins_are_sorted_as_strings() {
    let text = golden();
    let file = histogram(&text, "blocks", "truth-blocks");
    let bins: Vec<&str> = file
        .lines()
        .skip_while(|line| !line.starts_with("BIN"))
        .skip(1)
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').next().expect("a bin"))
        .collect();
    assert_eq!(bins, vec!["1,80", "100,20", "50,40", "50,60"]);
}

/// The concordance histogram counts bases: fifty for the first overlap, twenty for the second and
/// one for the pair that shares a single base.
#[test]
fn the_concordance_histogram_counts_bases() {
    let text = golden();
    let rows: Vec<(String, String)> = histogram(&text, "blocks", "concordance")
        .lines()
        .skip_while(|line| !line.starts_with("BIN"))
        .skip(1)
        .filter(|line| !line.is_empty())
        .map(|line| {
            let mut fields = line.split('\t');
            (
                fields.next().expect("a bin").to_string(),
                fields.next().expect("a value").to_string(),
            )
        })
        .collect();
    assert_eq!(
        rows,
        vec![
            ("20,30".to_string(), "50".to_string()),
            ("40,40".to_string(), "20".to_string()),
            ("60,70".to_string(), "1".to_string()),
        ]
    );
}

#[test]
fn a_multi_sample_record_is_refused_after_the_filter() {
    let text = golden();
    let rendered = refusal(&text, "two-samples")
        .split_once('"')
        .expect("a quoted record")
        .1
        .trim_end_matches("\".")
        .to_string();
    let error = MultiSample { rendered };
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "two-samples")
    );

    // And the block itself passes the hom-ref filter first, which is why the refusal happens at
    // all rather than the record being dropped.
    let parsed = blocks(&value(&text, "truth", "two-samples"));
    assert_eq!(parsed.len(), 1);
    assert!(parsed[0].is_hom_ref);
    assert_eq!(parsed[0].genotypes, 2);
    let steps = vec![(Some(0usize), Some(0usize))];
    assert!(accumulate(&parsed, &parsed, &steps).is_err());
}
