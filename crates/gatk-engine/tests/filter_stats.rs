//! Conformance for the Mutect filtering stats file against GATK 4.6.2.0, compared as the whole file
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/FilterStatsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the two roundings disagree**, two decimals in the columns and three in the metadata;
//!  * **both turn NaN into zero**, so a run with no passing call writes `0.0` rather than `NaN`;
//!  * **infinity saturates** to `9.223372036854776E15`;
//!  * **a negative rounds the other way**, the rounding being `floor(x + 0.5)`;
//!  * **and the rows are the caller's**, unfiltered and unsorted.

use gatk_corpus as corpus;
use gatk_engine::filtering_stats::{column, metadata, write_summary, FilterStats};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/filter_stats.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

/// The dump's own escaping, undone in one pass.
///
/// `ReferenceQueryDump.escape` escapes the backslash **first** and then the tab and the newline, so
/// undoing it with two independent replacements turns a filter name holding `\"` into one holding
/// `\\"`. This golden is the first one to carry a backslash at all, which is why it is the first
/// one that needs a real unescaper.
fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// The whole file the reference wrote for one run.
fn expected(text: &str, label: &str) -> String {
    unescape(
        rows(text, "table")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no table {label}"))[1],
    )
}

fn stats(name: &str, values: [f64; 4]) -> FilterStats {
    FilterStats {
        filter_name: name.to_string(),
        false_positive_count: values[0],
        false_discovery_rate: values[1],
        false_negative_count: values[2],
        false_negative_rate: values[3],
    }
}

fn pairs(list: &[(&str, &str)]) -> Vec<(String, String)> {
    list.iter()
        .map(|(key, value)| (key.to_string(), value.to_string()))
        .collect()
}

/// The seven runs of the dump, rebuilt.
fn ours(label: &str) -> String {
    let clustering = pairs(&[("clustering", "1"), ("other", "two")]);
    match label {
        "baseline" => write_summary(
            &[
                stats("weak_evidence", [3.0, 0.15, 2.0, 0.08]),
                stats("slippage", [1.0, 0.05, 4.0, 0.16]),
            ],
            &clustering,
            0.234,
            20.0,
            25.0,
            3.0,
            5.0,
        ),
        "rounding" => write_summary(
            &[stats("halves", [1.005, 2.675, 0.12345, 0.005])],
            &[],
            0.123456,
            3.0,
            2.0,
            1.0,
            1.0,
        ),
        "no-calls" => write_summary(
            &[stats("weak_evidence", [0.0, 0.0, 0.0, 0.0])],
            &[],
            0.5,
            0.0,
            0.0,
            0.0,
            0.0,
        ),
        "infinite-fdr" => write_summary(
            &[stats("weak_evidence", [1.0, 1.0, 0.0, 0.0])],
            &[],
            0.5,
            0.0,
            1.0,
            1.0,
            0.0,
        ),
        "negatives" => write_summary(
            &[stats("odd", [-1.005, -0.125, -2.5, -0.0])],
            &[],
            -0.0005,
            4.0,
            2.0,
            1.0,
            1.0,
        ),
        "no-rows" => write_summary(&[], &clustering, 0.1, 10.0, 9.0, 1.0, 2.0),
        "awkward-names" => write_summary(
            &[
                stats("has\ttab", [1.0, 0.1, 1.0, 0.1]),
                stats("has\"quote", [1.0, 0.1, 1.0, 0.1]),
                stats("has,comma", [1.0, 0.1, 1.0, 0.1]),
            ],
            &[],
            0.2,
            10.0,
            9.0,
            1.0,
            2.0,
        ),
        other => panic!("no run {other}"),
    }
}

#[test]
fn every_file_matches_the_golden_byte_for_byte() {
    let text = golden();
    let tables = rows(&text, "table");
    assert_eq!(tables.len(), 7, "one file per run");
    for row in tables {
        let label = row[0];
        assert_eq!(ours(label), expected(&text, label), "{label}");
    }
}

#[test]
fn the_two_roundings_disagree() {
    let text = golden();
    // The same value written twice: 0.123456 as a threshold and 0.12345 as a column.
    let rounding = expected(&text, "rounding");
    assert!(
        rounding.contains("#<METADATA>threshold=0.123\n"),
        "{rounding}"
    );
    assert!(rounding.contains("\t0.12\t"), "{rounding}");
    assert_eq!(column(0.123456), "0.12");
    assert_eq!(metadata(0.123456), "0.123");
    // And a half goes up in the columns.
    assert!(
        rounding.contains("halves\t1.01\t2.68\t0.12\t0.01\n"),
        "{rounding}"
    );
}

#[test]
fn a_division_by_zero_is_written_as_zero_or_as_a_very_large_number() {
    let text = golden();
    // 0/0 three times.
    let none = expected(&text, "no-calls");
    assert!(none.contains("#<METADATA>fdr=0.0\n"), "{none}");
    assert!(none.contains("#<METADATA>sensitivity=0.0\n"), "{none}");
    // And a division by zero that is not 0/0.
    let infinite = expected(&text, "infinite-fdr");
    assert!(
        infinite.contains("#<METADATA>fdr=9.223372036854776E15\n"),
        "{infinite}"
    );
}

#[test]
fn a_negative_rounds_the_other_way() {
    let text = golden();
    let negatives = expected(&text, "negatives");
    // -1.005 comes out -1.0 where 1.005 comes out 1.01, and -0.125 comes out -0.12.
    assert!(
        negatives.contains("odd\t-1.0\t-0.12\t-2.5\t0.0\n"),
        "{negatives}"
    );
    // And a small negative threshold loses its sign entirely.
    assert!(
        negatives.contains("#<METADATA>threshold=0.0\n"),
        "{negatives}"
    );
    assert_eq!(column(-1.005), "-1.0");
    assert_eq!(column(1.005), "1.01");
}

#[test]
fn the_rows_are_the_callers_and_a_file_may_have_none() {
    let text = golden();
    let none = expected(&text, "no-rows");
    assert!(none.ends_with("filter\tFP\tFDR\tFN\tFNR\n"), "{none}");
    // The clustering pairs still come first, and the three of the writer's own after them.
    assert!(
        none.starts_with("#<METADATA>clustering=1\n#<METADATA>other=two\n#<METADATA>threshold=")
    );
    // Two rows in the caller's order, neither of them dropped.
    let baseline = expected(&text, "baseline");
    let names: Vec<&str> = baseline
        .lines()
        .skip_while(|line| line.starts_with('#'))
        .skip(1)
        .filter_map(|line| line.split('\t').next())
        .collect();
    assert_eq!(names, vec!["weak_evidence", "slippage"]);
}

#[test]
fn a_name_that_needs_quoting_gets_it() {
    let text = golden();
    let awkward = expected(&text, "awkward-names");
    assert!(awkward.contains("\"has\ttab\"\t"), "{awkward}");
    assert!(awkward.contains("has,comma\t"), "{awkward}");
    // Whatever the writer does with a quote, the port does the same: the whole file was compared.
    assert_eq!(ours("awkward-names"), awkward);
}
