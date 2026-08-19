//! Conformance for `CollectReadCounts` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/CollectReadCountsDump.java`.
//!
//! # What this suite is for
//!
//!  * **a read is counted by its start**, so a read covering a whole interval without starting in
//!    it leaves that interval at zero;
//!  * **every requested interval is a row**, in the sorted order the argument layer produced;
//!  * **and the mapping quality threshold is this tool's own**, at 30.

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_tools::collect_read_counts::{self, ADDITIONAL_READ_FILTERS};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/collect_read_counts.txt.gz"),
    )
}

fn file(text: &str, label: &str) -> String {
    let prefix = format!("table\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {label}"))
        .replace("\\t", "\t")
        .replace("\\n", "\n")
}

fn interval(contig: &str, start: i32, end: i32) -> SimpleInterval {
    SimpleInterval::new(contig, start, end).expect("a valid interval")
}

const SEQUENCES: [(&str, i32); 2] = [("chr1", 200), ("chr2", 200)];

fn sequences() -> Vec<(String, i32)> {
    SEQUENCES
        .iter()
        .map(|(name, length)| (name.to_string(), *length))
        .collect()
}

/// The dump's six reads, as `(contig, start)` -- the only two fields the lookup reads.
///
/// `r006` starts at 70 at mapping quality 20 and is filtered before `apply`, so it is not here.
fn reads() -> Vec<(&'static str, i32)> {
    vec![
        ("chr1", 10),
        ("chr1", 12),
        ("chr1", 25),
        ("chr1", 40),
        ("chr1", 55),
    ]
}

fn table(intervals: &[SimpleInterval]) -> String {
    let counts = collect_read_counts::count(&reads(), intervals);
    collect_read_counts::write(&sequences(), "SAMPLE", intervals, &counts)
}

#[test]
fn every_table_matches_the_golden() {
    let text = golden();

    assert_eq!(
        table(&[
            interval("chr1", 10, 19),
            interval("chr1", 20, 29),
            interval("chr1", 30, 39),
            interval("chr1", 40, 49),
        ]),
        file(&text, "default")
    );

    // The forty-base read starting at 55 covers this interval entirely and starts before it.
    assert_eq!(
        table(&[interval("chr1", 56, 59)]),
        file(&text, "covered-not-started")
    );

    assert_eq!(
        table(&[interval("chr2", 1, 100)]),
        file(&text, "other-contig")
    );

    // The argument layer sorts, so the tool sees these in coordinate order.
    assert_eq!(
        table(&[interval("chr1", 10, 19), interval("chr1", 40, 49)]),
        file(&text, "out-of-order")
    );

    // The only read starting here is below the mapping quality threshold and never arrives.
    assert_eq!(
        table(&[interval("chr1", 70, 79)]),
        file(&text, "low-mapping-quality")
    );
}

/// A read spanning several intervals is counted in exactly one of them.
#[test]
fn a_read_is_counted_by_its_start_alone() {
    let intervals = [
        interval("chr1", 40, 49),
        interval("chr1", 50, 59),
        interval("chr1", 60, 69),
    ];
    // The twenty-base read at 40 spans the first two intervals.
    let counts = collect_read_counts::count(&[("chr1", 40)], &intervals);
    assert_eq!(counts, vec![1, 0, 0]);

    // And the forty-base read at 55 starts in the second and covers the third entirely.
    let counts = collect_read_counts::count(&[("chr1", 55)], &intervals);
    assert_eq!(counts, vec![0, 1, 0]);
}

/// A read never matches an interval on another contig.
#[test]
fn the_contig_has_to_match() {
    let intervals = [interval("chr2", 1, 100)];
    assert_eq!(collect_read_counts::count(&reads(), &intervals), vec![0]);
    assert!(collect_read_counts::interval_of_start("chr1", 10, &intervals).is_none());
}

/// The four filters this tool adds, and the threshold that is its own.
#[test]
fn the_tool_adds_four_filters() {
    assert_eq!(ADDITIONAL_READ_FILTERS.len(), 4);
    assert_eq!(ADDITIONAL_READ_FILTERS[3], "MappingQualityReadFilter");
    assert_eq!(collect_read_counts::DEFAULT_MINIMUM_MAPPING_QUALITY, 30);
}
