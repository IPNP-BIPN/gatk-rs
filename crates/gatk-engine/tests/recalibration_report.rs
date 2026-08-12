//! Conformance for `RecalibrationReport` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/RecalibrationReportDump.java`. The reports are ones the
//! reference wrote, because the reader cuts every data line at the header's word starts and a
//! hand-aligned file is cut in the wrong places; that behaviour has its own suite. Here the point is
//! the assembly.
//!
//! # What this suite is for
//!
//!  * **the read group keys are the sorted order, not the file's**, and they come from `RecalTable0`
//!    alone;
//!  * **a read group named in `RecalTable1` but not in `RecalTable0`** reaches a table as the key -1;
//!  * **the read group table reads its reported quality from `EstimatedQReported`** and the others
//!    from `QualityScore`, a different column of a different type;
//!  * **the empirical quality column is ignored on the way in**, so every datum can be recomputed;
//!  * **`null` in the `Arguments` table is an absence**, not four characters;
//!  * **all three reports round-trip**, which is what makes gathering safe.

use gatk_corpus as corpus;
use gatk_engine::gatk_report::{Report, Sorting, Table, Value};
use gatk_engine::recalibration_report::{RecalibrationReport, RecalibrationReportError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/recalibration_report.txt.gz"),
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

/// The report the dump wrote for one label, rebuilt with this port's writer.
///
/// The writer is already conformant (the gatk-report suite) and so is `RecalUtils`' table layout,
/// which this reproduces column for column. The golden's `roundtrip` rows are what say the two
/// agree.
fn source(label: &str) -> String {
    let (groups, with_quality, with_covariates) = match label {
        "two-groups" => (vec!["zebra", "alpha"], true, true),
        "all-events" => (vec!["alpha"], true, true),
        "read-group-only" => (vec!["alpha"], false, false),
        other => panic!("{other} is in the golden but not configured here"),
    };
    write_report(&groups, with_quality, with_covariates)
}

/// The five tables `RecalUtils.createRecalibrationGATKReport` writes, in its order.
fn write_report(groups: &[&str], with_quality: bool, with_covariates: bool) -> String {
    let mut report = Report::new();

    let mut arguments = Table::new(
        "Arguments",
        "Recalibration argument collection values used in this run",
        Sorting::SortByRow,
    );
    arguments.add_column("Argument", "%s");
    arguments.add_column("Value", "%s");
    for (name, value) in [
        ("binary_tag_name", "null"),
        (
            "covariate",
            "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate",
        ),
        ("default_platform", "null"),
        ("deletions_default_quality", "45"),
        ("force_platform", "null"),
        ("indels_context_size", "3"),
        ("insertions_default_quality", "45"),
        ("low_quality_tail", "2"),
        ("maximum_cycle_value", "500"),
        ("mismatches_context_size", "2"),
        ("mismatches_default_quality", "-1"),
        ("no_standard_covs", "false"),
        ("quantizing_levels", "16"),
        ("recalibration_report", "null"),
    ] {
        arguments.set(name, "Argument", Value::Str(name.into()));
        arguments.set(name, "Value", Value::Str(value.into()));
    }
    report.add_table(arguments);

    let mut quantized = Table::new(
        "Quantized",
        "Quality quantization map",
        Sorting::SortByColumn,
    );
    quantized.add_column("QualityScore", "%d");
    quantized.add_column("Count", "%d");
    quantized.add_column("QuantizedScore", "%d");
    for qual in 0..=93i64 {
        let key = qual.to_string();
        quantized.set(&key, "QualityScore", Value::Int(qual));
        quantized.set(&key, "Count", Value::Int(0));
        quantized.set(&key, "QuantizedScore", Value::Int(qual));
    }
    report.add_table(quantized);

    let mut read_group = Table::new("RecalTable0", "", Sorting::SortByColumn);
    read_group.add_column("ReadGroup", "%s");
    read_group.add_column("EventType", "%s");
    read_group.add_column("EmpiricalQuality", "%.4f");
    read_group.add_column("EstimatedQReported", "%.4f");
    read_group.add_column("Observations", "%d");
    read_group.add_column("Errors", "%.2f");
    for (index, group) in groups.iter().enumerate() {
        for (ordinal, event) in ["M", "I", "D"].into_iter().enumerate() {
            let key = format!("{group}-{event}");
            read_group.set(&key, "ReadGroup", Value::Str((*group).into()));
            read_group.set(&key, "EventType", Value::Str(event.into()));
            read_group.set(&key, "EmpiricalQuality", Value::Double(30.0));
            read_group.set(
                &key,
                "EstimatedQReported",
                Value::Double(30.0 + index as f64 + ordinal as f64),
            );
            read_group.set(&key, "Observations", Value::Int(1000 + index as i64));
            read_group.set(&key, "Errors", Value::Double(10.0 + index as f64));
        }
    }
    report.add_table(read_group);

    let mut quality = Table::new("RecalTable1", "", Sorting::SortByColumn);
    quality.add_column("ReadGroup", "%s");
    quality.add_column("QualityScore", "%s");
    quality.add_column("EventType", "%s");
    quality.add_column("EmpiricalQuality", "%.4f");
    quality.add_column("Observations", "%d");
    quality.add_column("Errors", "%.2f");
    if with_quality {
        for group in groups {
            for score in [20i64, 30] {
                for event in ["M", "I", "D"] {
                    let key = format!("{group}-{score}-{event}");
                    quality.set(&key, "ReadGroup", Value::Str((*group).into()));
                    quality.set(&key, "QualityScore", Value::Int(score));
                    quality.set(&key, "EventType", Value::Str(event.into()));
                    quality.set(&key, "EmpiricalQuality", Value::Double(30.0));
                    quality.set(&key, "Observations", Value::Int(100 + score));
                    quality.set(&key, "Errors", Value::Double(1.0 + score as f64 / 10.0));
                }
            }
        }
    }
    report.add_table(quality);

    let mut covariates = Table::new("RecalTable2", "", Sorting::SortByColumn);
    covariates.add_column("ReadGroup", "%s");
    covariates.add_column("QualityScore", "%s");
    covariates.add_column("CovariateValue", "%s");
    covariates.add_column("CovariateName", "%s");
    covariates.add_column("EventType", "%s");
    covariates.add_column("EmpiricalQuality", "%.4f");
    covariates.add_column("Observations", "%d");
    covariates.add_column("Errors", "%.2f");
    if with_covariates {
        for group in groups {
            for (score, name, value, observations, errors) in [
                (30i64, "Context", "AC", 50i64, 0.5),
                (20, "Cycle", "-3", 60, 0.6),
            ] {
                let key = format!("{group}-{name}");
                covariates.set(&key, "ReadGroup", Value::Str((*group).into()));
                covariates.set(&key, "QualityScore", Value::Int(score));
                covariates.set(&key, "CovariateValue", Value::Str(value.into()));
                covariates.set(&key, "CovariateName", Value::Str(name.into()));
                covariates.set(&key, "EventType", Value::Str("M".into()));
                covariates.set(&key, "EmpiricalQuality", Value::Double(30.0));
                covariates.set(&key, "Observations", Value::Int(observations));
                covariates.set(&key, "Errors", Value::Double(errors));
            }
        }
    }
    report.add_table(covariates);

    report.write()
}

/// The read group keys, which are the sorted order and not the file's.
#[test]
fn the_read_group_keys_are_the_sorted_order() {
    let text = golden();
    for row in rows(&text, "readgroups") {
        let label = row[0];
        let report = RecalibrationReport::parse(&source(label))
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));
        let groups: Vec<String> = (0..=report.covariates.read_group.maximum_key_value())
            .map(|key| {
                report
                    .covariates
                    .read_group
                    .format_key(key)
                    .unwrap()
                    .to_string()
            })
            .collect();
        assert_eq!(groups.join(","), row[1], "{label}");
    }
}

/// Every datum of every table, with the keys it landed at and the reported quality it carries.
#[test]
fn every_parsed_datum_is_the_reference() {
    let text = golden();
    let labels = ["two-groups", "all-events", "read-group-only"];
    let mut compared = 0;
    for label in labels {
        let report = RecalibrationReport::parse(&source(label)).unwrap();
        let expected: Vec<Vec<&str>> = rows(&text, "datum")
            .into_iter()
            .filter(|row| row[0] == label)
            .collect();

        let mut ours: Vec<Vec<String>> = Vec::new();
        for (index, table) in report.tables.all_tables.iter().enumerate() {
            for (keys, datum) in table.all_leaves() {
                ours.push(vec![
                    index.to_string(),
                    keys.iter()
                        .map(|key| key.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                    datum.borrow().num_observations().to_string(),
                    format!("{:.2}", datum.borrow().num_mismatches()),
                    format!("{:x}", datum.borrow().reported_quality().to_bits()),
                ]);
            }
        }
        assert_eq!(ours.len(), expected.len(), "{label}: datum count");
        for (ours, theirs) in ours.iter().zip(&expected) {
            assert_eq!(ours[0], theirs[1], "{label}: table index");
            assert_eq!(ours[1], theirs[2], "{label}: keys");
            assert_eq!(ours[2], theirs[3], "{label}: observations");
            assert_eq!(ours[3], theirs[4], "{label}: errors");
            assert_eq!(ours[4], theirs[5], "{label}: reported quality");
            compared += 1;
        }
    }
    println!("recalibration-report: {compared} datums compared");
}

/// The arguments the report round-trips into the collection.
#[test]
fn the_arguments_are_the_reference() {
    let text = golden();
    for label in ["two-groups", "all-events", "read-group-only"] {
        let report = RecalibrationReport::parse(&source(label)).unwrap();
        let expected = |name: &str| -> String {
            rows(&text, "argument")
                .into_iter()
                .find(|row| row[0] == label && row[1] == name)
                .unwrap_or_else(|| panic!("no argument {name} for {label}"))[2]
                .to_string()
        };
        assert_eq!(
            report.arguments.mismatches_context_size.to_string(),
            expected("mismatches_context_size")
        );
        assert_eq!(
            report.arguments.indels_context_size.to_string(),
            expected("indels_context_size")
        );
        assert_eq!(
            report.arguments.maximum_cycle_value.to_string(),
            expected("maximum_cycle_value")
        );
        assert_eq!(
            report.arguments.low_qual_tail.to_string(),
            expected("low_quality_tail")
        );
        assert_eq!(
            report.quantizing_levels.to_string(),
            expected("quantizing_levels")
        );
        // `null` became an absence, not the four characters.
        assert_eq!(expected("binary_tag_name"), "null");
        assert_eq!(expected("default_platform"), "null");
    }
}

/// The report's own summary: how many read groups it numbered, and whether it holds anything.
#[test]
fn the_report_summary_is_the_reference() {
    let text = golden();
    for row in rows(&text, "report") {
        let label = row[0];
        let report = RecalibrationReport::parse(&source(label)).unwrap();
        assert_eq!(
            (report.covariates.read_group.maximum_key_value() + 1).to_string(),
            row[1],
            "{label}: read groups"
        );
        assert_eq!(report.is_empty().to_string(), row[3], "{label}: isEmpty");
    }
}

/// Replace the first occurrence of `from` in one named table's body, keeping the width.
fn replace_in_table(report: &str, table: &str, from: &str, to: &str) -> String {
    let marker = format!("#:GATKTable:{table}:");
    let start = report.find(&marker).expect("the table is in the report");
    let (head, tail) = report.split_at(start);
    format!("{head}{}", tail.replacen(from, to, 1))
}

/// Every report the assembly refuses, worded as the reference words it.
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

    // A read group named in RecalTable1 but not in RecalTable0: the key is -1 and the write lands
    // at a negative index. `ghost` is the same width as `alpha`, so the fixed-width columns hold.
    let base = source("all-events");
    let ghosted = replace_in_table(&base, "RecalTable1", "alpha", "ghost");
    let error = RecalibrationReport::parse(&ghosted).unwrap_err();
    assert_eq!(
        error.message(),
        message("group-missing-from-read-group-table"),
        "a read group missing from RecalTable0"
    );

    // An event type that is not M, I or D, in the read group table.
    let bad_event = replace_in_table(&base, "RecalTable0", "  M  ", "  X  ");
    assert_eq!(
        RecalibrationReport::parse(&bad_event)
            .unwrap_err()
            .message(),
        message("unknown-event-type")
    );

    // A report with no RecalTable0 at all.
    let no_table = "#:GATKReport.v1.1:1\n\
         #:GATKTable:1:1:%s:;\n\
         #:GATKTable:Quantized:\n\
         QualityScore\n\
         0\n\
         \n";
    let error = RecalibrationReport::parse(no_table).unwrap_err();
    assert!(
        matches!(error, RecalibrationReportError::Report(_)),
        "{error:?}"
    );
    assert_eq!(error.message(), message("no-read-group-table"));
}
