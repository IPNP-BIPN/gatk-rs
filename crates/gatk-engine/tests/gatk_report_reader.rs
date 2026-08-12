//! Conformance for **reading** a `GATKReport` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/GATKReportReaderDump.java`. The writing side has its
//! own suite; this is the direction `ApplyBQSR` uses, and what it settles is which Java type each
//! cell comes back as, because that decides which branch `RecalibrationReport`'s `asLong`,
//! `asDouble` and `decodeByte` take.
//!
//! # What this suite is for
//!
//!  * **the columns are cut at the positions of the header line's words**, and no width is declared
//!    anywhere in the file;
//!  * **a hand-edited file is cut in the wrong places**, which is a `NumberFormatException` on a
//!    merged value and a `StringIndexOutOfBoundsException` on a short line;
//!  * **a `%d` column parses to `Long` and not `Integer`**;
//!  * **an empty format is written `%s` and read back as a string**, so an untyped double comes back
//!    as the characters it was rendered with;
//!  * **a parsed table is `DO_NOT_SORT`** whatever it was written with, which is what makes a parse
//!    and a second writing reproduce the first;
//!  * **`getReadGroups` sorts and deduplicates**.

use gatk_corpus as corpus;
use gatk_engine::gatk_report::{
    split_fixed_width, word_starts, DataType, Report, ReportReadError, Value,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/gatk_report_reader.txt.gz"),
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

/// The Java class the reference's parser produced for a cell, from the value this port produced.
fn java_class(value: &Value) -> &'static str {
    match value {
        Value::Null => "null",
        Value::Int(_) => "Long",
        Value::Double(_) => "Double",
        Value::Bool(_) => "Boolean",
        Value::Char(_) => "Character",
        Value::Str(_) => "String",
    }
}

/// `String.valueOf(value)` on a parsed cell, which is how the golden writes it.
fn value_text(value: &Value) -> String {
    match value {
        Value::Null => "null".to_string(),
        Value::Int(number) => number.to_string(),
        Value::Double(number) => {
            if number.is_nan() {
                "NaN".to_string()
            } else if *number == number.trunc() && number.abs() < 1e7 {
                format!("{number:.1}")
            } else {
                format!("{number}")
            }
        }
        Value::Bool(flag) => flag.to_string(),
        Value::Char(character) => character.to_string(),
        Value::Str(text) => text.clone(),
    }
}

/// The reference's `GATKReportDataType.toString()`, which is its format-matching pattern and not
/// its name.
fn data_type_pattern(data_type: DataType) -> &'static str {
    match data_type {
        DataType::Boolean => "%[Bb]",
        DataType::Character => "%[Cc]",
        DataType::Decimal => "%.*[EeFf]",
        DataType::Integer => "%[Dd]",
        // An empty format becomes `%s`, and its column reports the string pattern.
        DataType::String | DataType::Unknown => "%[Ss]",
    }
}

/// The reports the dump parsed, keyed by label. The four hand-written ones are literal; the three
/// the reference wrote are recovered from the golden's own round-trip claim by rebuilding them here.
fn hand_written(label: &str) -> Option<String> {
    Some(match label {
        "no-description" => "#:GATKReport.v1.1:1\n\
             #:GATKTable:2:1:%s:%d:;\n\
             #:GATKTable:Bare\n\
             Key  Value\n\
             k        7\n\
             \n"
        .to_string(),
        "no-rows" => "#:GATKReport.v1.1:1\n\
             #:GATKTable:2:0:%s:%d:;\n\
             #:GATKTable:Empty:nothing in it\n\
             Key  Value\n\
             \n"
        .to_string(),
        "ragged" => "#:GATKReport.v1.1:1\n\
             #:GATKTable:2:2:%s:%d:;\n\
             #:GATKTable:Ragged:values wider than the header\n\
             Key  Value\n\
             kkkkkk  7\n\
             k  8\n\
             \n"
        .to_string(),
        _ => return None,
    })
}

/// The three reports the reference wrote, rebuilt with this port's writer.
///
/// The writer is already conformant (the gatk-report suite), so writing them here and parsing the
/// result is the same text the reference parsed. The golden's `roundtrip` rows are what says so.
fn written(label: &str) -> Option<String> {
    use gatk_engine::gatk_report::{Sorting, Table};

    let mut report = Report::new();
    match label {
        "types" => {
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
            table.set("0", "Name", Value::Str("short".into()));
            table.set("0", "Count", Value::Int(1));
            table.set("0", "Rate", Value::Double(0.5));
            table.set("0", "Flag", Value::Bool(true));
            table.set("0", "Letter", Value::Char('A'));
            table.set("0", "Untyped", Value::Double(0.5));
            table.set("1", "Name", Value::Str("null".into()));
            table.set("1", "Count", Value::Int(1234567));
            table.set("1", "Rate", Value::Double(f64::NAN));
            table.set("1", "Flag", Value::Bool(false));
            table.set("1", "Letter", Value::Char('z'));
            table.set("1", "Untyped", Value::Str("text".into()));
            report.add_table(table);
        }
        "two" => {
            let mut first = Table::new("First", "the first table", Sorting::DoNotSort);
            first.add_column("Argument", "%s");
            first.set("0", "Argument", Value::Str("value".into()));
            report.add_table(first);
            let mut second = Table::new("Second", "the second table", Sorting::DoNotSort);
            second.add_column("Key", "%s");
            second.add_column("Value", "%d");
            second.set("0", "Key", Value::Str("k".into()));
            second.set("0", "Value", Value::Int(7));
            report.add_table(second);
        }
        "sorted" => {
            let mut table = Table::new("Sorted", "rows out of order", Sorting::SortByColumn);
            table.add_column("RowKey", "%s");
            table.add_column("Value", "%d");
            for key in ["bbb", "aaa", "ccc"] {
                table.set(key, "RowKey", Value::Str(key.into()));
                table.set(key, "Value", Value::Int(key.len() as i64));
            }
            report.add_table(table);
        }
        _ => return None,
    }
    Some(report.write())
}

fn source(label: &str) -> String {
    hand_written(label)
        .or_else(|| written(label))
        .unwrap_or_else(|| panic!("{label} is in the golden but not configured here"))
}

/// The fixed-width split, which is where every parsed column boundary comes from.
#[test]
fn the_fixed_width_split_is_the_reference() {
    let text = golden();
    for row in rows(&text, "starts") {
        let line = row[0].replace('_', " ");
        let ours: Vec<String> = word_starts(&line)
            .into_iter()
            .map(|start| start.to_string())
            .collect();
        let expected = row.get(1).copied().unwrap_or("");
        assert_eq!(ours.join(","), expected, "getWordStarts({:?})", row[0]);
    }

    let header = "Alpha  Beta  Gamma";
    let header_starts = word_starts(header);
    for row in rows(&text, "split") {
        let (line, starts) = match row[0].strip_prefix("header:") {
            Some(data) => (data.replace('_', " "), header_starts.clone()),
            None => {
                let line = row[0].replace('_', " ");
                let starts = word_starts(&line);
                (line, starts)
            }
        };
        let expected = row.get(1).copied().unwrap_or("");
        match split_fixed_width(&line, &starts) {
            Ok(fields) => assert_eq!(fields.join("|"), expected, "split({:?})", row[0]),
            Err(error) => assert_eq!(
                format!("E:StringIndexOutOfBoundsException:{}", error.message()),
                expected,
                "split({:?})",
                row[0]
            ),
        }
    }
}

/// Every report the dump parsed: its tables, its columns' recovered types, and every cell's class.
#[test]
fn every_parsed_report_is_the_reference() {
    let text = golden();
    let labels: Vec<String> = {
        let mut seen: Vec<String> = Vec::new();
        for row in rows(&text, "report") {
            if !seen.iter().any(|label| label == row[0]) {
                seen.push(row[0].to_string());
            }
        }
        seen
    };
    assert_eq!(labels.len(), 5, "five reports parsed successfully");

    let mut cells = 0;
    for label in &labels {
        let report = Report::parse(&source(label))
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));

        let expected_report = rows(&text, "report")
            .into_iter()
            .find(|row| row[0] == *label)
            .unwrap();
        assert_eq!(expected_report[1], "v1.1", "{label}: version");
        assert_eq!(
            report.tables.len().to_string(),
            expected_report[2],
            "{label}: table count"
        );

        for row in rows(&text, "table")
            .into_iter()
            .filter(|row| row[0] == *label)
        {
            let index: usize = row[1].parse().unwrap();
            let table = &report.tables[index];
            assert_eq!(table.name, row[2], "{label}: table {index} name");
            assert_eq!(
                table.description, row[3],
                "{label}: table {index} description"
            );
            assert_eq!(table.rows.len().to_string(), row[4], "{label}: rows");
            assert_eq!(table.columns.len().to_string(), row[5], "{label}: columns");
        }

        for row in rows(&text, "column")
            .into_iter()
            .filter(|row| row[0] == *label)
        {
            let table = report.table_named(row[1]).unwrap();
            let index: usize = row[2].parse().unwrap();
            let column = &table.columns[index];
            assert_eq!(column.name, row[3], "{label}/{}: column name", row[1]);
            assert_eq!(column.format, row[4], "{label}/{}: format", row[1]);
            assert_eq!(
                data_type_pattern(column.data_type),
                row[5],
                "{label}/{}: data type",
                row[1]
            );
        }

        for row in rows(&text, "cell")
            .into_iter()
            .filter(|row| row[0] == *label)
        {
            let table = report.table_named(row[1]).unwrap();
            let index: usize = row[2].parse().unwrap();
            let column = table
                .columns
                .iter()
                .position(|column| column.name == row[3])
                .unwrap();
            let value = &table.rows[index][column];
            assert_eq!(
                java_class(value),
                row[4],
                "{label}/{}/{}: class",
                row[1],
                row[3]
            );
            assert_eq!(
                value_text(value),
                row[5],
                "{label}/{}/{}: value",
                row[1],
                row[3]
            );
            cells += 1;
        }
    }
    println!("gatk-report-reader: {cells} cells compared");
}

/// A parse and a second writing, which is what makes a gathered report reproducible.
#[test]
fn the_round_trips_are_the_references() {
    let text = golden();
    for row in rows(&text, "roundtrip") {
        let label = row[0];
        let original = source(label);
        let report = Report::parse(&original).unwrap();
        assert_eq!(
            (report.write() == original).to_string(),
            row[1],
            "{label}: round trip"
        );
    }
}

/// `getReadGroups`, which reads one named table and sorts what it finds.
#[test]
fn the_read_groups_are_sorted_and_deduplicated() {
    let text = golden();
    let recal = "#:GATKReport.v1.1:1\n\
         #:GATKTable:3:4:%s:%s:%d:;\n\
         #:GATKTable:RecalTable0:\n\
         ReadGroup  EventType  EmpiricalQuality\n\
         zebra      M                        30\n\
         alpha      M                        29\n\
         alpha      I                        45\n\
         middle     D                        45\n\
         \n";
    let report = Report::parse(recal).unwrap();
    let expected = rows(&text, "readgroups")
        .into_iter()
        .find(|row| row[0] == "recal")
        .unwrap()[1]
        .to_string();
    assert_eq!(report.read_groups().unwrap().join(","), expected);

    // Four rows and three read groups: the set is what removes the repeat.
    let table = report.table_named("RecalTable0").unwrap();
    let expected_table = rows(&text, "table")
        .into_iter()
        .find(|row| row[0] == "recal")
        .unwrap();
    assert_eq!(table.rows.len().to_string(), expected_table[4]);
    assert_eq!(table.columns.len().to_string(), expected_table[5]);
    assert_eq!(table.description, expected_table[3]);
}

/// Every input the reader refuses, worded as the reference words it.
#[test]
fn the_refusals_are_worded_like_the_reference() {
    let text = golden();
    let message = |what: &str| -> String {
        rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == what)
            .unwrap_or_else(|| panic!("no error row {what}"))[2]
            .to_string()
    };

    assert_eq!(
        Report::parse("").unwrap_err().message(),
        message("empty-stream")
    );
    assert_eq!(
        Report::parse("#:GATKReport.v0.1:1\n")
            .unwrap_err()
            .message(),
        message("legacy-version")
    );
    assert_eq!(
        Report::parse("#:GATKReport.v9.9:1\n")
            .unwrap_err()
            .message(),
        message("no-such-version")
    );
    assert_eq!(
        Report::parse("hello\n").unwrap_err().message(),
        message("not-a-report")
    );
    assert_eq!(
        Report::parse(
            "#:GATKReport.v1.1:2\n\
             #:GATKTable:1:1:%s:;\n\
             #:GATKTable:One:\n\
             Key\n\
             k\n\
             \n"
        )
        .unwrap_err()
        .message(),
        message("too-few-tables")
    );
    assert_eq!(
        Report::parse(
            "#:GATKReport.v1.1:1\n\
             #:GATKTable:1:3:%s:;\n\
             #:GATKTable:One:\n\
             Key\n\
             k\n\
             \n"
        )
        .unwrap_err()
        .message(),
        message("too-few-rows")
    );
    assert_eq!(
        Report::parse(
            "#:GATKReport.v1.1:1\n\
             #:GATKTable:1:1:%d:;\n\
             #:GATKTable:One:\n\
             Value\n\
             abc\n\
             \n"
        )
        .unwrap_err()
        .message(),
        message("unparseable-integer")
    );
    // The ragged file, whose merged value a `%d` column refuses.
    assert_eq!(
        Report::parse(&source("ragged")).unwrap_err().message(),
        rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == "read@ragged")
            .unwrap()[2]
    );

    let report = Report::parse(
        "#:GATKReport.v1.1:1\n\
         #:GATKTable:1:1:%s:;\n\
         #:GATKTable:One:\n\
         Key\n\
         k\n\
         \n",
    )
    .unwrap();
    assert_eq!(
        report.table_named("Nonesuch").unwrap_err().message(),
        message("unknown-table")
    );
    assert_eq!(
        report.table_named("Nonesuch").unwrap_err(),
        ReportReadError::NoSuchTable("Nonesuch".to_string())
    );
}
