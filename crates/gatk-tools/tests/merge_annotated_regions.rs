//! Conformance for `MergeAnnotatedRegions` against GATK 4.6.2.0, compared as the whole output file
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/MergeAnnotatedRegionsDump.java`.
//!
//! # What this suite is for
//!
//!  * **abutting regions are not merged**, `IntervalUtils.overlaps` being a real overlap;
//!  * **a conflicting annotation is split on the separator, deduplicated, sorted and rejoined**,
//!    so `b`, `a`, `b` gives `a__b` and a value that already carries the separator is split first;
//!  * **an annotation missing from one side passes through unchanged**;
//!  * **the annotation columns come out alphabetical** and the locatable ones renamed;
//!  * **`Position` names both the start and the end**, so a file with one coordinate column parses
//!    as one-base regions;
//!  * **the rows are sorted by the dictionary first**, which is what lets a chain merge;
//!  * **an unknown contig passes through** rather than being refused;
//!  * **and a file of no rows throws**, the annotations being read off the first record.

use gatk_corpus as corpus;
use gatk_tools::annotated_interval::{merge_regions, read, CollectionError, DEFAULT_SEPARATOR};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/merge_annotated_regions.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> Option<String> {
    let prefix = format!("{kind}\t{label}=");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .map(unescape)
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .map(unescape)
        .unwrap_or_else(|| panic!("the golden carries error/{label}"))
}

fn dictionary() -> Vec<String> {
    vec!["chr1".to_string(), "chr2".to_string()]
}

/// The input file of each labelled run, as the dump wrote it.
fn input(label: &str) -> &'static str {
    match label {
        "plain" => {
            "CONTIG\tSTART\tEND\tname\tvalue\tcall\n\
             chr1\t1\t100\tone\t0.5\t+\n\
             chr1\t50\t150\ttwo\t0.5\t+\n\
             chr1\t300\t400\tthree\t1.5\t-\n"
        }
        "abutting" => "CONTIG\tSTART\tEND\tname\nchr1\t1\t100\ta\nchr1\t101\t200\tb\n",
        "chain" => {
            "CONTIG\tSTART\tEND\tname\n\
             chr1\t1\t100\tb\n\
             chr1\t50\t150\ta\n\
             chr1\t120\t200\tb\n"
        }
        "skipped-overlap" => {
            "CONTIG\tSTART\tEND\tname\n\
             chr1\t1\t100\ta\n\
             chr1\t120\t130\tb\n\
             chr1\t90\t200\tc\n"
        }
        "missing-annotation" => {
            "CONTIG\tSTART\tEND\tname\tvalue\nchr1\t1\t100\ta\t\nchr1\t50\t150\tb\t7\n"
        }
        "separator-in-value" => {
            "CONTIG\tSTART\tEND\tname\nchr1\t1\t100\tb__c\nchr1\t50\t150\ta__b\n"
        }
        "unsorted" => {
            "CONTIG\tSTART\tEND\tname\n\
             chr2\t1\t100\tz\n\
             chr1\t200\t240\ty\n\
             chr1\t1\t100\tx\n"
        }
        "other-column-names" => {
            "Chromosome\tStart_Position\tEnd_Position\tname\nchr1\t1\t100\ta\nchr1\t50\t150\tb\n"
        }
        "position-column" => "chrom\tPosition\tname\nchr1\t10\ta\nchr1\t10\tb\nchr1\t20\tc\n",
        "comments" => "#a note\n#another\nCONTIG\tSTART\tEND\tname\nchr1\t1\t100\ta\n",
        "sam-header" => {
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:240\nCONTIG\tSTART\tEND\tname\nchr1\t1\t100\ta\n"
        }
        "no-rows" => "CONTIG\tSTART\tEND\tname\n",
        "unknown-contig" => "CONTIG\tSTART\tEND\tname\nchrX\t1\t100\ta\n",
        "no-locatable-columns" => "name\tvalue\na\t1\n",
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_run_writes_what_the_reference_writes() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "plain",
        "abutting",
        "chain",
        "skipped-overlap",
        "missing-annotation",
        "separator-in-value",
        "unsorted",
        "other-column-names",
        "position-column",
        "comments",
        "sam-header",
        "unknown-contig",
    ] {
        let mut collection = read(input(label)).expect("a file the codec accepts");
        collection.records = merge_regions(&collection.records, &dictionary(), DEFAULT_SEPARATOR);
        assert_eq!(
            collection.write(),
            value(&text, "merged", label).expect("the golden carries the output"),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 12, "the golden's outputs");
}

#[test]
fn the_two_refusals_carry_the_references_messages() {
    let text = golden();
    for (label, expected) in [
        ("no-rows", CollectionError::NoRecords),
        ("no-locatable-columns", CollectionError::NoLocatableColumns),
    ] {
        let error = read(input(label)).expect_err("a file the collection refuses");
        assert_eq!(error, expected);
        // The tribble message names the whole contig list and then the input's URI, which is the
        // container's working directory: the reader here is given text rather than a path, so the
        // source is supplied to the message.
        let source = "file:///work/merge-annotated-regions-dump/no-locatable-columns.seg";
        let line = format!(
            "{}:{}",
            error.java_class(),
            error.message_with_source(if label == "no-locatable-columns" {
                source
            } else {
                ""
            })
        );
        assert_eq!(line, refusal(&text, label), "{label}");
    }
}
