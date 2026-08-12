//! Conformance for GATK's tsv table format against 4.6.2.0, compared as written text and as the
//! values read back.
//!
//! Golden from `tools/readfilter-conformance/TsvTableDump.java`.
//!
//! # What this suite is for
//!
//!  * **metadata is a tagged comment**, and an untagged one that looks like metadata is not;
//!  * **the map is filled by the comment hook**, so a subclass that overrides it loses every pair,
//!    which the golden holds as two readings of the same file;
//!  * **a value is quoted only when it has to be**, and a comma is not one of the reasons;
//!  * **a double is `String.valueOf`'s spelling**;
//!  * **and the three refusals are three different things**, of which one is not even a
//!    `UserException`.

use gatk_corpus as corpus;
use gatk_engine::tsv_table::{quote_if_needed, write_table, Table, TableError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/tsv_table.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// The reverse of the dump's `escape`, which turned `\` into `\\` before it touched tabs and
/// newlines. Scanning once is the only correct reverse: replacing `\\t` first would turn a real
/// backslash followed by a `t` into a tab.
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

/// The file each label was written as, from the golden itself.
fn written(text: &str, label: &str) -> String {
    let row = rows(text, "written")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no file {label}"));
    unescape(row[1])
}

/// The rows each written label was built from, mirroring the dump's tables.
fn source_rows(label: &str) -> Vec<Vec<String>> {
    let row = |contig: &str, position: &str, value: &str, note: &str| {
        vec![
            contig.to_string(),
            position.to_string(),
            value.to_string(),
            note.to_string(),
        ]
    };
    match label {
        "plain" => vec![
            row("chr1", "100", "1.0", "ordinary"),
            row("chr2", "200", "0.5", "another"),
        ],
        "quoting" => vec![
            row("chr1", "1", "1.0", "has a space"),
            row("chr1", "2", "1.0", "has\ta tab"),
            row("chr1", "3", "1.0", "has a \"quote\""),
            row("chr1", "4", "1.0", "has a \\backslash"),
            row("chr1", "5", "1.0", "has a , comma"),
            row("chr1", "6", "1.0", ""),
            row("chr1", "7", "1.0", "#starts with a comment prefix"),
        ],
        // `String.valueOf` for each of the doubles the dump wrote.
        "doubles" => vec![
            row("chr1", "1", "1.0", "integral"),
            row("chr1", "2", "0.1", "tenth"),
            row("chr1", "3", "0.3333333333333333", "third"),
            row("chr1", "4", "1.0E-7", "small"),
            row("chr1", "5", "1.0E21", "large"),
            row("chr1", "6", "NaN", "nan"),
            row("chr1", "7", "Infinity", "infinity"),
        ],
        "empty" | "no-header" => vec![],
        other => panic!("no written case {other}"),
    }
}

#[test]
fn every_written_table_is_the_reference() {
    let text = golden();
    for label in ["plain", "quoting", "doubles", "empty", "no-header"] {
        let ours = write_table(
            &["contig", "position", "value", "note"],
            &source_rows(label),
            &[("sample", "s1")],
        );
        assert_eq!(ours, written(&text, label), "written/{label}");
    }
}

#[test]
fn every_read_row_is_the_reference() {
    let text = golden();
    for label in ["plain", "quoting", "doubles"] {
        let table = Table::parse(&written(&text, label), "x").expect("the table parses");
        let expected: Vec<String> = rows(&text, "read")
            .into_iter()
            .filter(|row| row[0] == label)
            .map(|row| unescape(row[2]))
            .collect();
        let ours: Vec<String> = table.rows.iter().map(|row| row.join(",")).collect();
        assert_eq!(ours, expected, "read/{label}");
    }
}

/// The tagged comment fills the map; the untagged one does not.
#[test]
fn the_metadata_tag_is_what_makes_metadata() {
    let text = golden();

    // The reader that overrides the hook saw no metadata at all.
    for row in rows(&text, "metadata") {
        assert_eq!(row.get(1).copied().unwrap_or(""), "", "metadata/{}", row[0]);
    }
    // The one that did not saw the pair.
    for row in rows(&text, "metadatadefault") {
        assert_eq!(row[1], "sample=s1", "metadatadefault/{}", row[0]);
    }

    // And the port collects it whichever way it is read, from the tagged line only.
    let table = Table::parse(&written(&text, "plain"), "x").expect("parses");
    assert_eq!(table.metadata.get("sample").map(String::as_str), Some("s1"));

    let hand_written = Table::parse(&written(&text, "short-row"), "x");
    // The hand-written file's `#sample=s1` is a comment, not metadata, and the row is still short.
    assert!(hand_written.is_err());
}

#[test]
fn every_refusal_is_the_reference() {
    let text = golden();
    let expected = |label: &str| -> String {
        rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no refusal {label}"))[1]
            .to_string()
    };

    // A row of the wrong width, in both directions, is the same message.
    for label in ["short-row", "long-row"] {
        let file = written(&text, label);
        let error = Table::parse(&file, &format!("tsvtable-dump/{label}.table"))
            .expect_err("this table is refused");
        assert_eq!(
            format!("{}:Bad input: {}", error.java_class(), error.message()),
            expected(label),
            "error/{label}"
        );
    }

    // A missing column is a plain IllegalArgumentException, raised when it is asked for.
    let table = Table::parse(&written(&text, "missing-column"), "x").expect("the table parses");
    let error = table.get(&table.rows[0], "note").unwrap_err();
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        expected("missing-column")
    );

    // And an empty numeric field only fails when it is asked for as a number.
    let table = Table::parse(&written(&text, "empty-number"), "x").expect("the table parses");
    assert_eq!(table.get(&table.rows[0], "position").expect("a value"), "");
    let error = table
        .get_int(
            &table.rows[0],
            "position",
            "tsvtable-dump/empty-number.table",
            3,
        )
        .unwrap_err();
    assert!(matches!(error, TableError::NotAnInteger { .. }));
    assert_eq!(
        format!("{}:Bad input: {}", error.java_class(), error.message()),
        expected("empty-number")
    );
}

/// A comma forces nothing, which is the difference between this and a CSV writer.
#[test]
fn a_comma_is_not_a_reason_to_quote() {
    assert_eq!(quote_if_needed("has a , comma"), "has a , comma");
    assert_eq!(quote_if_needed("has\ta tab"), "\"has\ta tab\"");
}
