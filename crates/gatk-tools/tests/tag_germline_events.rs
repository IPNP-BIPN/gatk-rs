//! Conformance for `TagGermlineEvents` against GATK 4.6.2.0, compared as the whole output file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/TagGermlineEventsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the tagging is an OR over two different lists**: `shifted` tags through the reciprocal
//!    overlap of the merged tumour run while its breakpoints are out of reach of the padding, and
//!    `reciprocal-default` is the same pair failing at the default threshold;
//!  * **only non-neutral merged normal regions are considered**, which `all-neutral` shows;
//!  * **the per-segment filter is a third test**, whose intersection comparison is strict;
//!  * **the default tag is `0`**, so an untagged segment carries a value;
//!  * **the tag column lands where its name sorts**, which `other-call-column` shows by putting
//!    `POSSIBLE_GERMLINE` before `call_state`;
//!  * **and an empty call on either side is a refusal**, the normal being checked first.

use gatk_corpus as corpus;
use gatk_tools::annotated_interval::read;
use gatk_tools::tag_germline_events::{
    tag_tumour_segments, TagError, DEFAULT_PADDING_IN_BP, DEFAULT_RECIPROCAL_THRESHOLD,
    GERMLINE_TAG_HEADER,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/tag_germline_events.txt.gz"),
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

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

fn dictionary() -> Vec<String> {
    vec!["chr1".to_string(), "chr2".to_string()]
}

const TUMOUR: &str = "CONTIG\tSTART\tEND\tCALL\tMEAN_LOG2_COPY_RATIO\tNUM_POINTS\n\
                      chr1\t1\t100\t+\t0.7\t10\n\
                      chr1\t101\t200\t+\t0.7\t10\n\
                      chr1\t201\t300\t0\t0.0\t10\n\
                      chr2\t1\t100\t-\t-0.8\t10\n";

/// The two files, the call column, the padding and the threshold of each labelled run.
fn run(label: &str) -> (&'static str, &'static str, &'static str, i32, f64) {
    let padding = DEFAULT_PADDING_IN_BP;
    let threshold = DEFAULT_RECIPROCAL_THRESHOLD;
    match label {
        "exact-match" => (
            TUMOUR,
            "CONTIG\tSTART\tEND\tCALL\n\
             chr1\t1\t100\t+\n\
             chr1\t101\t200\t+\n\
             chr1\t201\t300\t0\n\
             chr2\t1\t100\t0\n",
            "CALL",
            padding,
            threshold,
        ),
        "shifted" | "shifted-padded" => (
            TUMOUR,
            "CONTIG\tSTART\tEND\tCALL\n\
             chr1\t20\t180\t+\n\
             chr1\t201\t300\t0\n\
             chr2\t1\t100\t0\n",
            "CALL",
            if label == "shifted-padded" {
                25
            } else {
                padding
            },
            threshold,
        ),
        "all-neutral" => (
            TUMOUR,
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t300\t0\nchr2\t1\t100\t0\n",
            "CALL",
            padding,
            threshold,
        ),
        "second-contig" => (
            TUMOUR,
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t300\t0\nchr2\t1\t100\t-\n",
            "CALL",
            padding,
            threshold,
        ),
        "normal-contains" => (
            TUMOUR,
            "CONTIG\tSTART\tEND\tCALL\n\
             chr1\t1\t200\t+\n\
             chr1\t201\t300\t0\n\
             chr2\t1\t100\t0\n",
            "CALL",
            padding,
            threshold,
        ),
        "reciprocal" | "reciprocal-default" | "zero-threshold" => (
            TUMOUR,
            "CONTIG\tSTART\tEND\tCALL\n\
             chr1\t40\t240\t+\n\
             chr1\t241\t300\t0\n\
             chr2\t1\t100\t0\n",
            "CALL",
            padding,
            match label {
                "reciprocal" => 0.5,
                "zero-threshold" => 0.0,
                _ => threshold,
            },
        ),
        "other-call-column" => (
            "CONTIG\tSTART\tEND\tcall_state\nchr1\t1\t100\t+\nchr1\t101\t200\t0\n",
            "CONTIG\tSTART\tEND\tcall_state\nchr1\t1\t100\t+\nchr1\t101\t200\t0\n",
            "call_state",
            padding,
            threshold,
        ),
        "empty-tumour-call" => (
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t\n",
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n",
            "CALL",
            padding,
            threshold,
        ),
        "empty-normal-call" => (
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n",
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t\n",
            "CALL",
            padding,
            threshold,
        ),
        "negative-padding" => (
            TUMOUR,
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n",
            "CALL",
            -1,
            threshold,
        ),
        "threshold-above-one" => (
            TUMOUR,
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n",
            "CALL",
            padding,
            1.5,
        ),
        "overlapping-tumour" => (
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\nchr1\t50\t150\t+\n",
            "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n",
            "CALL",
            padding,
            threshold,
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// One run, tagged and written as the tool writes it.
fn tagged(label: &str) -> Result<String, TagError> {
    let (tumour, normal, call, padding, threshold) = run(label);
    let mut collection = read(tumour).expect("a tumour file the codec accepts");
    let normal = read(normal).expect("a normal file the codec accepts");
    let records = tag_tumour_segments(
        &collection.records,
        &normal.records,
        call,
        &dictionary(),
        GERMLINE_TAG_HEADER,
        padding,
        threshold,
    )?;
    collection.records = records;
    if !collection
        .annotations
        .iter()
        .any(|name| name == GERMLINE_TAG_HEADER)
    {
        collection.annotations.push(GERMLINE_TAG_HEADER.to_string());
        collection.annotations.sort();
    }
    Ok(collection.write())
}

#[test]
fn every_run_writes_what_the_reference_writes() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "exact-match",
        "shifted",
        "shifted-padded",
        "all-neutral",
        "second-contig",
        "normal-contains",
        "reciprocal",
        "reciprocal-default",
        "zero-threshold",
        "other-call-column",
    ] {
        let ours = tagged(label).expect("a run that finishes");
        assert_eq!(ours, value(&text, "tagged", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 10, "the golden's outputs");
}

#[test]
fn the_five_refusals_carry_the_references_messages() {
    let text = golden();
    let mut refused = 0;
    for label in [
        "empty-tumour-call",
        "empty-normal-call",
        "negative-padding",
        "threshold-above-one",
        "overlapping-tumour",
    ] {
        let error = tagged(label).expect_err("a run the tagger refuses");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
        refused += 1;
    }
    assert_eq!(refused, 5, "the golden's refusals");
}
