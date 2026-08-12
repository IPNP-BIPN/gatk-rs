//! Conformance for `GATKReport` against GATK 4.6.2.0, compared as **text, space by space**.
//!
//! Golden from `tools/readfilter-conformance/GATKReportDump.java`. This format is decided by
//! spaces, and a diff cannot show them, so the golden carries every line twice: once as written and
//! once with its spaces replaced by underscores. The port is compared against the second.
//!
//! # What this suite is for
//!
//! The file format both BQSR tools are written in, settled before either tool that uses it:
//!
//!  * **the column width is the widest formatted value**, name included, and a `%.4f` column is as
//!    wide as its rendering rather than as its `toString`;
//!  * **columns are separated by exactly two spaces**, with the padding inside each column, so a
//!    left-aligned last column carries trailing spaces;
//!  * **alignment is right until a value asks for left**, and five renderings are exempt;
//!  * **two escapes in `writeRow`**: an untyped column renders a double `%.8f`, and a non-finite
//!    double leaves its own format entirely;
//!  * **the header carries the formats and not the widths**, which is why a parse and a second
//!    writing reproduce the first.

use gatk_corpus as corpus;
use gatk_engine::gatk_report::{Alignment, Report, Sorting, Table, Value};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gatk_report.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// The lines of one labelled report, with the spaces still shown as underscores.
fn lines_of(text: &str, label: &str) -> Vec<String> {
    let mut lines: Vec<(usize, String)> = rows(text, "line")
        .into_iter()
        .filter(|row| row[0] == label)
        .map(|row| (row[1].parse::<usize>().unwrap(), row[2].to_string()))
        .collect();
    lines.sort_by_key(|(n, _)| *n);
    lines.into_iter().map(|(_, line)| line).collect()
}

/// The same rendering the golden uses, so the two can be compared at all.
fn underscored(text: &str) -> Vec<String> {
    text.split('\n')
        .map(|line| line.replace(' ', "_"))
        .collect()
}

/// One column per data type, with the two `writeRow` escapes in it.
fn types() -> Report {
    let mut table = Table::new("Types", "one column per data type", Sorting::DoNotSort);
    for (name, format) in [
        ("Name", "%s"),
        ("Count", "%d"),
        ("Rate", "%.4f"),
        ("Flag", "%b"),
        ("Letter", "%c"),
        ("Untyped", ""),
    ] {
        table.add_column(name, format);
    }
    let rows = [
        (
            "0",
            Value::Str("short".into()),
            1i64,
            Value::Double(0.5),
            true,
            'A',
            Value::Double(0.5),
        ),
        (
            "1",
            Value::Str("a considerably longer value".into()),
            1234567,
            Value::Double(0.123456789),
            false,
            'z',
            Value::Double(1.0 / 3.0),
        ),
        (
            "2",
            Value::Null,
            0,
            Value::Double(f64::NAN),
            true,
            'x',
            Value::Double(f64::INFINITY),
        ),
        (
            "3",
            Value::Str("neg".into()),
            -42,
            Value::Double(f64::NEG_INFINITY),
            false,
            'y',
            Value::Double(f64::NAN),
        ),
    ];
    for (key, name, count, rate, flag, letter, untyped) in rows {
        table.set(key, "Name", name);
        table.set(key, "Count", Value::Int(count));
        table.set(key, "Rate", rate);
        table.set(key, "Flag", Value::Bool(flag));
        table.set(key, "Letter", Value::Char(letter));
        table.set(key, "Untyped", untyped);
    }
    let mut report = Report::new();
    report.add_table(table);
    report
}

/// Numeric throughout, so nothing ever asks for a left alignment.
fn numeric() -> Report {
    let mut table = Table::new("Numeric", "right aligned throughout", Sorting::DoNotSort);
    table.add_column("Quality", "%d");
    table.add_column("EmpiricalQuality", "%.4f");
    for i in 0..3i64 {
        let key = i.to_string();
        table.set(&key, "Quality", Value::Int(10 + i * 10));
        table.set(
            &key,
            "EmpiricalQuality",
            Value::Double(10.0 + i as f64 * 10.0 + 0.5),
        );
    }
    let mut report = Report::new();
    report.add_table(table);
    report
}

/// Two tables in one report, which is the shape of a recalibration report.
fn two() -> Report {
    let mut first = Table::new("First", "the first table", Sorting::DoNotSort);
    first.add_column("Argument", "%s");
    first.set("0", "Argument", Value::Str("value".into()));
    let mut second = Table::new("Second", "the second table", Sorting::DoNotSort);
    second.add_column("Key", "%s");
    second.add_column("Value", "%d");
    second.set("0", "Key", Value::Str("k".into()));
    second.set("0", "Value", Value::Int(7));
    let mut report = Report::new();
    report.add_table(first);
    report.add_table(second);
    report
}

/// The same rows added out of order, under one sorting.
fn sorted(sorting: Sorting) -> Report {
    let mut table = Table::new("Sorted", "rows added out of order", sorting);
    table.add_column("RowKey", "%s");
    table.add_column("Value", "%d");
    for (key, value) in [("bbb", 2i64), ("aaa", 1), ("ccc", 3)] {
        table.set(key, "RowKey", Value::Str(key.into()));
        table.set(key, "Value", Value::Int(value));
    }
    let mut report = Report::new();
    report.add_table(table);
    report
}

fn report_for(label: &str) -> Report {
    match label {
        "types" => types(),
        "numeric" => numeric(),
        "two" => two(),
        "sort-SORT_BY_COLUMN" => sorted(Sorting::SortByColumn),
        "sort-SORT_BY_ROW" => sorted(Sorting::SortByRow),
        "sort-DO_NOT_SORT" => sorted(Sorting::DoNotSort),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_report_is_written_space_for_space() {
    let text = golden();
    let labels: Vec<String> = rows(&text, "roundtrip")
        .into_iter()
        .map(|row| row[0].to_string())
        .collect();
    assert_eq!(labels.len(), 6, "six reports");

    for label in &labels {
        let expected = lines_of(&text, label);
        let ours = underscored(&report_for(label).write());
        assert_eq!(ours.len(), expected.len(), "{label}: line count");
        for (n, (ours, theirs)) in ours.iter().zip(&expected).enumerate() {
            assert_eq!(ours, theirs, "{label}: line {n}");
        }
    }
    println!(
        "gatk-report: {} reports written space for space",
        labels.len()
    );
}

/// The widths and alignments the format computed, asserted directly rather than through the bytes.
#[test]
fn the_widths_and_alignments_are_the_references() {
    let text = golden();
    for row in rows(&text, "width") {
        let (label, name, width, alignment) = (row[0], row[1], row[2], row[3]);
        let report = report_for(label);
        let column = report
            .tables
            .iter()
            .flat_map(|table| &table.columns)
            .find(|column| column.name == name)
            .unwrap_or_else(|| panic!("{label}: no column {name}"));
        assert_eq!(column.width().to_string(), width, "{label}/{name}: width");
        let ours = match column.alignment() {
            Alignment::Left => "LEFT",
            Alignment::Right => "RIGHT",
        };
        assert_eq!(ours, alignment, "{label}/{name}: alignment");
    }
}

/// The header carries the formats and not the widths, which is what makes a re-write reproducible.
#[test]
fn the_reference_round_tripped_every_report() {
    let text = golden();
    let roundtrips = rows(&text, "roundtrip");
    assert_eq!(roundtrips.len(), 6);
    for row in roundtrips {
        assert_eq!(row[1], "true", "{}", row[0]);
    }
    // And the version line the reference wrote.
    let version = text
        .lines()
        .find_map(|line| line.strip_prefix("version\t"))
        .expect("the golden lost its version row");
    assert_eq!(version, "#:GATKReport.v1.1:1");
}
