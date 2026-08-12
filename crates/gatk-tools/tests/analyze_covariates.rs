//! Conformance for `AnalyzeCovariates` against GATK 4.6.2.0, compared as the whole intermediate
//! csv and as the refusals.
//!
//! Golden from `tools/readfilter-conformance/AnalyzeCovariatesDump.java`.
//!
//! # What this suite is for
//!
//!  * **the quality score table is filed one past the last covariate index**, and that is what
//!    names the row `QualityScore`;
//!  * **the optional covariates drop the reported quality from their key**, so a context row is
//!    summed over it;
//!  * **the read group table never reaches the csv**;
//!  * **the modes come out in the map's order**, not the command line's;
//!  * **and the consistency check fires on a mismatch context and not on an indel context**.

use gatk_corpus as corpus;
use gatk_engine::recalibration_report::RecalibrationReport;
use gatk_tools::analyze_covariates::{analyze_covariates, AnalyzeCovariatesError, RoleReport};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/analyze_covariates.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// The reverse of the dump's `escape`, scanning once so a real backslash is never read as a tab.
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

fn labelled(text: &str, kind: &str, label: &str) -> String {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no {kind} row for {label}"))
        .get(1)
        .copied()
        .unwrap_or("")
        .to_string()
}

/// One recalibration report, parsed from the text the reference wrote.
fn report(text: &str, label: &str) -> RecalibrationReport {
    RecalibrationReport::parse(&unescape(&labelled(text, "report", label)))
        .unwrap_or_else(|error| panic!("report {label}: {}", error.message()))
}

/// Every csv the golden holds, each from the reports its label was built from.
#[test]
fn every_csv_is_the_reference() {
    let text = golden();
    let plain = report(&text, "plain");
    let groups = report(&text, "two-groups");
    let after = report(&text, "after");
    let unchecked = report(&text, "unchecked-argument");

    let cases: Vec<(&str, Vec<RoleReport<'_>>)> = vec![
        (
            "one-report",
            vec![RoleReport {
                role: "BQSR",
                report: &plain,
            }],
        ),
        (
            "two-groups",
            vec![RoleReport {
                role: "BQSR",
                report: &groups,
            }],
        ),
        (
            "before-after",
            vec![
                RoleReport {
                    role: "Before",
                    report: &plain,
                },
                RoleReport {
                    role: "After",
                    report: &after,
                },
            ],
        ),
        // The same two, given the other way round: the csv must not change.
        (
            "after-before",
            vec![
                RoleReport {
                    role: "After",
                    report: &after,
                },
                RoleReport {
                    role: "Before",
                    report: &plain,
                },
            ],
        ),
        (
            "all-three",
            vec![
                RoleReport {
                    role: "BQSR",
                    report: &plain,
                },
                RoleReport {
                    role: "Before",
                    report: &plain,
                },
                RoleReport {
                    role: "After",
                    report: &after,
                },
            ],
        ),
        // An argument the check never looks at, so the two are combined.
        (
            "unchecked-argument",
            vec![
                RoleReport {
                    role: "Before",
                    report: &plain,
                },
                RoleReport {
                    role: "After",
                    report: &unchecked,
                },
            ],
        ),
    ];

    for (label, reports) in &cases {
        let ours = analyze_covariates(reports, true)
            .unwrap_or_else(|error| panic!("csv/{label}: {}", error.message()));
        assert_eq!(
            ours,
            unescape(&labelled(&text, "csv", label)),
            "csv/{label}"
        );
    }
}

/// The read group table has data in every report and no row of its own in any csv.
#[test]
fn no_row_of_the_csv_comes_from_the_read_group_table() {
    let text = golden();
    let plain = report(&text, "plain");
    assert!(
        !plain.tables.read_group_table().all_leaves().is_empty(),
        "the read group table is not empty to begin with"
    );
    let csv = analyze_covariates(
        &[RoleReport {
            role: "BQSR",
            report: &plain,
        }],
        true,
    )
    .expect("a csv");
    for line in csv.lines().skip(1) {
        let name = line.split(',').nth(2).expect("a covariate name");
        assert_ne!(name, "ReadGroup", "{line}");
    }
}

/// The context rows are summed over the reported quality, which is the dropped key.
#[test]
fn a_context_row_is_the_sum_over_the_reported_qualities() {
    let text = golden();
    let csv = unescape(&labelled(&text, "csv", "one-report"));
    let context = csv
        .lines()
        .find(|line| line.contains(",Context,"))
        .expect("a context row");
    // 50 + 20 and 50 + 30, one row.
    assert_eq!(
        context,
        "alpha,AC,Context,Base Substitution,150,1.00,22.00,22.84,-0.84,BQSR"
    );
}

#[test]
fn every_refusal_is_the_reference() {
    let text = golden();
    let expected = |label: &str| labelled(&text, "error", label);
    let plain = report(&text, "plain");
    let different = report(&text, "different-arguments");

    let no_report = analyze_covariates(&[], true).expect_err("no report at all");
    assert_eq!(
        format!("{}:{}", no_report.java_class(), no_report.message()),
        expected("no-report")
    );

    let no_output = analyze_covariates(
        &[RoleReport {
            role: "BQSR",
            report: &plain,
        }],
        false,
    )
    .expect_err("no output requested");
    assert_eq!(
        format!("{}:{}", no_output.java_class(), no_output.message()),
        expected("no-output")
    );

    let inconsistent = analyze_covariates(
        &[
            RoleReport {
                role: "Before",
                report: &plain,
            },
            RoleReport {
                role: "After",
                report: &different,
            },
        ],
        true,
    )
    .expect_err("the mismatch context differs");
    assert!(matches!(
        inconsistent,
        AnalyzeCovariatesError::IncompatibleParameters(_)
    ));
    assert_eq!(
        format!("{}:{}", inconsistent.java_class(), inconsistent.message()),
        expected("inconsistent")
    );
}
