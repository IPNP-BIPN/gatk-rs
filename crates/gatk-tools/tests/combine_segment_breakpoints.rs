//! Conformance for `CombineSegmentBreakpoints` against GATK 4.6.2.0, compared as the whole output
//! file of every run.
//!
//! Golden from `tools/readfilter-conformance/CombineSegmentBreakpointsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the breakpoints are pooled and sorted with starts before ends**, which `single-base` is
//!    there for;
//!  * **a piece between two segments is shrunk at both ends and dropped when they were adjacent**,
//!    which `disjoint` and `adjacent` separate;
//!  * **a stretch only one file covers is still a row**, with empty strings from the other side;
//!  * **a column both files carry is suffixed with that file's label** and one only a single file
//!    carries keeps its name;
//!  * **and the output columns are sorted alphabetically**, so `labels` puts `CALL_normal` before
//!    `CALL_tumour` whatever order the files were given in.
//!
//! # What is compared, and what is supplied
//!
//! Every byte of every output. The SAM header lines are supplied to the writer rather than derived
//! here: `@HD VN:1.6 GO:none SO:coordinate` and the merge of the two inputs' dictionaries with the
//! reference's are htsjdk's `SAMFileHeader`, not this tool's, and the golden's `M5` and `UR` are
//! the harness's mask.

use gatk_corpus as corpus;
use gatk_tools::annotated_interval::read;
use gatk_tools::combine_segment_breakpoints::{combine, DEFAULT_LABELS};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/combine_segment_breakpoints.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, label: &str) -> String {
    let prefix = format!("combined\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {label}")),
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

/// The header lines the merged `SAMFileHeader` writes, which are htsjdk's rather than this tool's.
fn header_lines(only_first_contig: bool) -> Vec<String> {
    let mut lines = vec!["@HD\tVN:1.6\tGO:none\tSO:coordinate".to_string()];
    if only_first_contig {
        lines.push("@SQ\tSN:chr1\tLN:240".to_string());
    } else {
        for contig in ["chr1", "chr2"] {
            lines.push(format!(
                "@SQ\tSN:{contig}\tLN:240\tM5:<masked>\tUR:<masked>"
            ));
        }
    }
    lines
}

const FIRST: &str = "CONTIG\tSTART\tEND\tCALL\tMEAN\n\
                     chr1\t1\t100\t+\t0.5\n\
                     chr1\t101\t200\t0\t0.0\n\
                     chr2\t1\t100\t-\t-0.5\n";

const SECOND: &str = "CONTIG\tSTART\tEND\tCALL\tNAME\n\
                      chr1\t50\t150\t+\tsecond-one\n\
                      chr1\t180\t260\t-\tsecond-two\n\
                      chr2\t200\t240\t+\tsecond-three\n";

/// The two files, the labels and the columns of interest of each run.
fn run(label: &str) -> (&'static str, &'static str, [&'static str; 2], Vec<String>) {
    let columns = |names: &[&str]| names.iter().map(|n| (*n).to_string()).collect();
    match label {
        "default" => (
            FIRST,
            SECOND,
            DEFAULT_LABELS,
            columns(&["CALL", "MEAN", "NAME"]),
        ),
        "labels" => (
            FIRST,
            SECOND,
            ["tumour", "normal"],
            columns(&["CALL", "MEAN", "NAME"]),
        ),
        "column-of-interest" => (FIRST, SECOND, DEFAULT_LABELS, columns(&["CALL"])),
        "columns-from-both" => (FIRST, SECOND, DEFAULT_LABELS, columns(&["MEAN", "NAME"])),
        "identical" => (FIRST, FIRST, DEFAULT_LABELS, columns(&["CALL", "MEAN"])),
        "adjacent" => (
            "CONTIG\tSTART\tEND\tA\nchr1\t1\t100\ta1\nchr1\t101\t200\ta2\n",
            "CONTIG\tSTART\tEND\tB\nchr1\t1\t200\tb1\n",
            DEFAULT_LABELS,
            columns(&["A", "B"]),
        ),
        "disjoint" => (
            "CONTIG\tSTART\tEND\tA\nchr1\t1\t50\ta1\n",
            "CONTIG\tSTART\tEND\tB\nchr1\t150\t200\tb1\n",
            DEFAULT_LABELS,
            columns(&["A", "B"]),
        ),
        "single-base" => (
            "CONTIG\tSTART\tEND\tA\nchr1\t100\t100\ta1\n",
            "CONTIG\tSTART\tEND\tB\nchr1\t1\t200\tb1\n",
            DEFAULT_LABELS,
            columns(&["A", "B"]),
        ),
        "one-sam-header" => (
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:240\nCONTIG\tSTART\tEND\tA\nchr1\t1\t100\ta1\n",
            "CONTIG\tSTART\tEND\tB\nchr1\t1\t100\tb1\n",
            DEFAULT_LABELS,
            columns(&["A", "B"]),
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_run_writes_what_the_reference_writes() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "default",
        "labels",
        "column-of-interest",
        "columns-from-both",
        "identical",
        "adjacent",
        "disjoint",
        "single-base",
        "one-sam-header",
    ] {
        let (first, second, labels, columns) = run(label);
        let first = read(first).expect("a file the codec accepts");
        let second = read(second).expect("a file the codec accepts");
        let records = combine(
            &first.records,
            &second.records,
            &dictionary(),
            labels,
            &columns,
        );
        let annotations: Vec<String> = records
            .first()
            .map(|record| record.annotations.keys().cloned().collect())
            .unwrap_or_default();
        let mut collection = first;
        collection.header_lines = header_lines(label == "one-sam-header");
        collection.comments.clear();
        collection.annotations = annotations;
        collection.records = records;
        assert_eq!(collection.write(), value(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 9, "the golden's outputs");
}

#[test]
fn the_three_refusals_are_the_tools_and_the_collections() {
    let text = golden();
    // The tool's own check on the columns of interest, and the collection's two, which the
    // `merge-annotated-regions` suite already pins. All three are asserted as the golden's text.
    assert_eq!(
        refusal(&text, "missing-column"),
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput:Bad input: Some columns \
         of interest specified by the user were not seen in any input files: ABSENT"
    );
    assert_eq!(
        refusal(&text, "empty-file"),
        "java.lang.IndexOutOfBoundsException:Index 0 out of bounds for length 0"
    );
    assert_eq!(
        refusal(&text, "overlapping-input"),
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput:Bad input: Overlap \
         detected in input:  chr1:1-100 overlapped chr1:1-100, chr1:50-150"
    );
}
