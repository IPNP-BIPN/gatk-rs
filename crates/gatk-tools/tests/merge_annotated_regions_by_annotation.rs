//! Conformance for `MergeAnnotatedRegionsByAnnotation` against GATK 4.6.2.0, compared as the whole
//! output file of every run.
//!
//! Golden from `tools/readfilter-conformance/MergeAnnotatedRegionsByAnnotationDump.java`.
//!
//! # What this suite is for
//!
//!  * **abutting regions are one apart, not zero**, so `zero-distance` merges nothing;
//!  * **the comparison is against the region built so far**, so `default-distance` walks a chain
//!    of hops and merges three regions spanning twelve hundred bases;
//!  * **the annotations not named are still reconciled** with the separator, which is what turns
//!    `a`, `b`, `c` into `a__b__c`;
//!  * **the three output column names are arguments**;
//!  * **and this tool writes through the writer**, so no SAM header reaches the output, not even
//!    the `@CO` lines its sibling writes.

use gatk_corpus as corpus;
use gatk_tools::annotated_interval::{
    merge_regions_by_annotation, read, write_without_header, CollectionError, DEFAULT_SEPARATOR,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/merge_annotated_regions_by_annotation.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, label: &str) -> String {
    let prefix = format!("merged\t{label}=");
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

const SEGMENTS: &str = "CONTIG\tSTART\tEND\tCALL\tNAME\n\
                        chr1\t1\t100\t+\ta\n\
                        chr1\t101\t200\t+\tb\n\
                        chr1\t1201\t1300\t+\tc\n\
                        chr1\t1301\t1400\t-\td\n\
                        chr2\t1\t100\t+\te\n\
                        chr2\t50000\t50100\t+\tf\n";

/// The input, the annotations matched on, the distance and the three column names of each run.
fn run(label: &str) -> (&'static str, Vec<String>, i64, [&'static str; 3]) {
    let default_columns = ["CONTIG", "START", "END"];
    let call = vec!["CALL".to_string()];
    match label {
        "default-distance" => (SEGMENTS, call, 1_000_000, default_columns),
        "zero-distance" => (SEGMENTS, call, 0, default_columns),
        "middle-distance" => (SEGMENTS, call, 1001, default_columns),
        "match-on-name" => (SEGMENTS, vec!["NAME".to_string()], 1_000_000, default_columns),
        "match-on-both" => (
            SEGMENTS,
            vec!["CALL".to_string(), "NAME".to_string()],
            1_000_000,
            default_columns,
        ),
        "renamed-columns" => (
            SEGMENTS,
            call,
            1_000_000,
            ["chrom", "chromStart", "chromEnd"],
        ),
        "sam-header" => (
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:240\nCONTIG\tSTART\tEND\tCALL\tNAME\nchr1\t1\t100\t+\ta\n",
            call,
            1_000_000,
            default_columns,
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_run_writes_what_the_reference_writes() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "default-distance",
        "zero-distance",
        "middle-distance",
        "match-on-name",
        "match-on-both",
        "renamed-columns",
        "sam-header",
    ] {
        let (input, annotations, distance, columns) = run(label);
        let collection = read(input).expect("a file the codec accepts");
        let merged = merge_regions_by_annotation(
            &collection.records,
            &dictionary(),
            &annotations,
            DEFAULT_SEPARATOR,
            distance,
        );
        let ours = write_without_header(&merged, columns[0], columns[1], columns[2])
            .expect("a run with regions");
        assert_eq!(ours, value(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 7, "the golden's outputs");
}

#[test]
fn the_two_refusals_are_the_tools_own_argument_checks() {
    let text = golden();
    // Both are the tool's argument validation rather than anything the merge does, so they are
    // asserted as the constants the golden carries.
    assert_eq!(
        refusal(&text, "missing-annotation"),
        "java.lang.IllegalArgumentException:Input file did not have all of the specified \
         annotations.  Missing annotations were: ABSENT"
    );
    assert_eq!(
        refusal(&text, "negative-distance"),
        "java.lang.IllegalArgumentException:Cannot have a negative value for distance."
    );
    // The writer's own refusal, which no run of the golden reaches: a merge of nothing has no
    // first region to take the columns from.
    assert_eq!(
        write_without_header(&[], "CONTIG", "START", "END"),
        Err(CollectionError::NoRecords)
    );
}
