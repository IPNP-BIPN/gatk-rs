//! Conformance for `GroundTruthScorer` against GATK 4.6.2.0, compared as the CSV and the report of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/GroundTruthScorerDump.java`.
//!
//! The flow-based scoring is not measured or ported. What is compared is the report's shape, the
//! percentile columns, the NaN an emptied flow produces, and the refusals.
//!
//! # What this suite is for
//!
//!  * **the report being four tables in a fixed order**;
//!  * **`--quality-percentiles` naming columns of the first**;
//!  * **`--exclude-zero-flows` emptying a row rather than dropping it**;
//!  * **`--add-mean-call` adding two columns**;
//!  * **the score threshold being compared against a negative default**;
//!  * **an empty CSV still carrying its header**;
//!  * **`--gt-no-output` and the threshold being told apart by the report**;
//!  * **the four-level table's row count**;
//!  * **and a read that is not flow-based being refused by name and position.**

use gatk_corpus as corpus;
use gatk_tools::ground_truth_scorer::{
    bin_to_base, bin_to_deviation, deviation_to_bin, four_level_row_count, four_level_table_name,
    keeps_read, keeps_row, one_level_table_name, percentile_columns, phred, two_level_table_name,
    Accumulator, Arguments, PercentileReport, BASE_VALUE_MAX, DEFAULT_FLOW_ORDER,
    DEFAULT_QUALITY_PERCENTILES, HMER_VALUE_MAX, MEAN_CALL_COLUMNS,
    NORMALIZED_SCORE_THRESHOLD_DEFAULT, PERCENTILE_FIXED_COLUMNS, PERCENTILE_TABLE_NAME,
    QUAL_VALUE_MAX,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/ground_truth_scorer.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{name}")),
    )
}

/// The table names one report declares, in order.
fn table_names(text: &str, label: &str) -> Vec<String> {
    section(text, "report", label)
        .lines()
        .filter_map(|line| line.strip_prefix("#:GATKTable:"))
        .filter(|rest| !rest.starts_with(|c: char| c.is_ascii_digit()))
        .map(|rest| rest.split(':').next().expect("a name").to_string())
        .collect()
}

/// The declared width and height of each table, from its `#:GATKTable:<cols>:<rows>:` line.
fn table_shapes(text: &str, label: &str) -> Vec<(usize, usize)> {
    section(text, "report", label)
        .lines()
        .filter_map(|line| line.strip_prefix("#:GATKTable:"))
        .filter(|rest| rest.starts_with(|c: char| c.is_ascii_digit()))
        .map(|rest| {
            let mut parts = rest.split(':');
            (
                parts.next().expect("a width").parse().expect("a number"),
                parts.next().expect("a height").parse().expect("a number"),
            )
        })
        .collect()
}

/// One run's CSV, as its header and its rows.
fn csv(text: &str, label: &str) -> (Vec<String>, Vec<String>) {
    let file = section(text, "out", label);
    let mut lines = file.lines().filter(|line| !line.is_empty());
    let header = lines
        .next()
        .expect("a header")
        .split(',')
        .map(str::to_string)
        .collect();
    (header, lines.map(str::to_string).collect())
}

/// The first table's data rows, split on whitespace.
fn percentile_rows(text: &str, label: &str) -> Vec<Vec<String>> {
    section(text, "report", label)
        .lines()
        .skip_while(|line| !line.trim_start().starts_with("flow"))
        .skip(1)
        .take_while(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| line.split_whitespace().map(str::to_string).collect())
        .collect()
}

/// Four tables, in the order the tool builds them.
#[test]
fn the_report_is_four_tables() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "default",
        "percentiles",
        "exclude-zero-flows",
        "mean-call",
        "use-softclipped",
        "no-output",
        "threshold-mild",
        "threshold-above-everything",
    ] {
        assert_eq!(
            table_names(&text, label),
            vec![
                PERCENTILE_TABLE_NAME.to_string(),
                one_level_table_name("qual"),
                two_level_table_name("qual", "hmer"),
                four_level_table_name("qual", "hmer", "deviation", "base"),
            ],
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 8, "the runs the port reproduces");
    // The names are built from the column names and not written out.
    assert_eq!(one_level_table_name("qual"), "qualReport");
    assert_eq!(two_level_table_name("qual", "hmer"), "qual_hmerReport");
    assert_eq!(
        four_level_table_name("qual", "hmer", "deviation", "base"),
        "qual_hmer_deviation_base_Report"
    );
}

/// Seven fixed columns plus one per percentile.
#[test]
fn the_percentile_columns_are_named_by_the_argument() {
    let text = golden();
    assert_eq!(PERCENTILE_FIXED_COLUMNS.len(), 7);
    // The default five make twelve columns.
    assert_eq!(
        percentile_columns(DEFAULT_QUALITY_PERCENTILES),
        vec![
            "flow", "count", "min", "max", "mean", "median", "std", "p10", "p25", "p50", "p75",
            "p90"
        ]
    );
    assert_eq!(table_shapes(&text, "default")[0].0, 12);
    // Three make ten.
    assert_eq!(percentile_columns("5,50,95").len(), 10);
    assert_eq!(table_shapes(&text, "percentiles")[0].0, 10);
    // And the header line names them.
    let header: Vec<String> = section(&text, "report", "percentiles")
        .lines()
        .find(|line| line.trim_start().starts_with("flow"))
        .expect("a header")
        .split_whitespace()
        .map(str::to_string)
        .collect();
    assert_eq!(header, percentile_columns("5,50,95"));
}

/// The row stays with a count of zero and every statistic comes out NaN.
#[test]
fn excluding_the_zero_flows_empties_a_row_rather_than_dropping_it() {
    let text = golden();
    // Both runs have the same number of rows.
    assert_eq!(
        table_shapes(&text, "default")[0].1,
        table_shapes(&text, "exclude-zero-flows")[0].1
    );
    let plain = percentile_rows(&text, "default");
    let excluded = percentile_rows(&text, "exclude-zero-flows");
    assert_eq!(plain.len(), excluded.len());
    // The first flow is emptied: a count of zero and NaN for everything after it.
    assert_ne!(plain[0][1], "0");
    assert_eq!(excluded[0][1], "0");
    for cell in &excluded[0][2..] {
        assert_eq!(cell, "NaN");
    }
    // Which is what an empty series gives.
    let empty = PercentileReport::default();
    let row = empty.row(0, DEFAULT_QUALITY_PERCENTILES);
    assert_eq!(row[1], 0.0);
    assert!(row[2..].iter().all(|value| value.is_nan()));
    // A series with one probability in it holds its phred, not the probability.
    let mut one = PercentileReport::default();
    one.add_probability(0.012);
    let row = one.row(0, "50");
    assert_eq!(row[1], 1.0);
    assert!((row[2] - (-10.0 * 0.012f64.log10())).abs() < 1e-12);
}

/// Two more columns, and the ones already there unchanged.
#[test]
fn the_mean_call_adds_two_columns() {
    let text = golden();
    let (plain, plain_rows) = csv(&text, "default");
    let (extended, extended_rows) = csv(&text, "mean-call");
    assert_eq!(plain.len(), 15);
    assert_eq!(extended.len(), 17);
    // The first fifteen are the same, in the same order.
    assert_eq!(extended[..15], plain[..]);
    assert_eq!(extended[15..], MEAN_CALL_COLUMNS);
    // And the same reads are written either way.
    assert_eq!(plain_rows.len(), extended_rows.len());
    let name = |row: &str| row.split(',').next().expect("a name").to_string();
    assert_eq!(
        plain_rows.iter().map(|r| name(r)).collect::<Vec<_>>(),
        extended_rows.iter().map(|r| name(r)).collect::<Vec<_>>()
    );
}

/// The default is negative, so a positive threshold is already above every score.
#[test]
fn the_score_threshold_is_compared_against_a_negative_default() {
    let text = golden();
    assert_eq!(NORMALIZED_SCORE_THRESHOLD_DEFAULT, -0.1);
    let arguments = Arguments::default();
    assert!(keeps_read(0.0, &arguments));
    assert!(keeps_read(-0.1, &arguments), "the boundary is kept");
    assert!(!keeps_read(-0.2, &arguments));
    // A threshold of 0.1 is above every score this fixture produces, so both threshold runs are
    // empty rather than one being mild.
    let (_, mild) = csv(&text, "threshold-mild");
    let (_, above) = csv(&text, "threshold-above-everything");
    assert!(mild.is_empty());
    assert!(above.is_empty());
    let raised = Arguments {
        normalized_score_threshold: 0.1,
        ..Arguments::default()
    };
    assert!(!keeps_read(0.0, &raised));
}

/// The two arguments empty the CSV alike and differ in the report.
#[test]
fn an_empty_csv_still_carries_its_header() {
    let text = golden();
    for label in ["threshold-mild", "threshold-above-everything", "no-output"] {
        let (header, rows) = csv(&text, label);
        assert_eq!(header.len(), 15, "{label}");
        assert!(rows.is_empty(), "{label}");
    }
    // The threshold empties the report's first table as well; --gt-no-output does not.
    assert_eq!(table_shapes(&text, "threshold-mild")[0].1, 0);
    assert_eq!(
        table_shapes(&text, "no-output")[0].1,
        table_shapes(&text, "default")[0].1
    );
    // So the two are told apart by the report and not by the CSV.
    assert_eq!(
        section(&text, "out", "no-output"),
        section(&text, "out", "threshold-mild")
    );
    assert_ne!(
        section(&text, "report", "no-output"),
        section(&text, "report", "threshold-mild")
    );
}

/// Its size is the product of the allocated dimensions, not of what was observed.
#[test]
fn the_four_level_table_is_allocated_rather_than_grown() {
    assert_eq!(QUAL_VALUE_MAX, 60);
    assert_eq!(HMER_VALUE_MAX, 100);
    assert_eq!(BASE_VALUE_MAX, 3);
    assert_eq!(deviation_to_bin(101), 202);
    assert_eq!(four_level_row_count(), 61 * 101 * 202 * 4);
    assert!(four_level_row_count() > 4_900_000);
    // The golden's own runs all omit the zeros, so their tables are a handful of rows.
    let text = golden();
    let shapes = table_shapes(&text, "default");
    assert!(shapes[3].1 < 1000, "{:?}", shapes[3]);
}

/// The origin row is always written; the zeros elsewhere are what the option drops.
#[test]
fn omitting_the_zeros_keeps_the_origin() {
    assert!(keeps_row(true, &[0], 0));
    assert!(keeps_row(true, &[0, 0], 0));
    assert!(keeps_row(true, &[0, 0, 0, 0], 0));
    assert!(!keeps_row(true, &[1], 0));
    assert!(!keeps_row(true, &[0, 1], 0));
    assert!(keeps_row(true, &[1], 5));
    // With the option off, every row is written whatever it holds.
    assert!(keeps_row(false, &[7], 0));
}

/// A rate of zero over a non-empty cell reports the threshold's phred, not zero.
#[test]
fn a_cell_with_no_errors_reports_the_thresholds_phred() {
    let mut cell = Accumulator::default();
    for _ in 0..10 {
        cell.add(true);
    }
    assert_eq!(cell.count(), 10);
    assert_eq!(cell.false_rate(), 0.0);
    // With a threshold, the phred is the threshold's.
    let expected = (-10.0 * 0.003f64.log10()).ceil() as i64;
    assert_eq!(phred(0.0, 10, 0.003), expected);
    assert_eq!(expected, 26);
    // With no threshold, or over an empty cell, it is zero.
    assert_eq!(phred(0.0, 10, 0.0), 0);
    assert_eq!(phred(0.0, 0, 0.003), 0);
    // And an ordinary rate is its own phred, rounded up.
    assert_eq!(phred(0.1, 10, 0.003), 10);
    assert_eq!(phred(0.5, 10, 0.003), 4);
    // An empty cell's rate is zero rather than a NaN.
    assert_eq!(Accumulator::default().false_rate(), 0.0);
}

/// `0,-1,1,-2,2` become `0,1,2,3,4`, and the string form is not the number's.
#[test]
fn the_deviation_folds_into_a_bin() {
    assert_eq!(deviation_to_bin(0), 0);
    assert_eq!(deviation_to_bin(-1), 1);
    assert_eq!(deviation_to_bin(1), 2);
    assert_eq!(deviation_to_bin(-2), 3);
    assert_eq!(deviation_to_bin(2), 4);
    // The string form signs the positives and not the negatives or zero.
    assert_eq!(bin_to_deviation(0), "0");
    assert_eq!(bin_to_deviation(1), "-1");
    assert_eq!(bin_to_deviation(2), "+1");
    assert_eq!(bin_to_deviation(3), "-2");
    assert_eq!(bin_to_deviation(4), "+2");
    // The base bin indexes the flow order.
    assert_eq!(DEFAULT_FLOW_ORDER, "TGCA");
    assert_eq!(bin_to_base(0), 'T');
    assert_eq!(bin_to_base(3), 'A');
}

/// Both files move when the soft-clipped bases are scored.
#[test]
fn scoring_the_soft_clipped_bases_changes_both_files() {
    let text = golden();
    assert_ne!(
        section(&text, "out", "default"),
        section(&text, "out", "use-softclipped")
    );
    assert_ne!(
        section(&text, "report", "default"),
        section(&text, "report", "use-softclipped")
    );
    // It also changes WHICH reads are written, and it LOSES one: scoring the clipped bases takes
    // the soft-clipped read out of the file rather than adding anything to it.
    let (_, plain) = csv(&text, "default");
    let (_, clipped) = csv(&text, "use-softclipped");
    assert_eq!(plain.len(), 5);
    assert_eq!(clipped.len(), 4);
    let name = |row: &str| row.split(',').next().expect("a name").to_string();
    let plain_names: Vec<String> = plain.iter().map(|r| name(r)).collect();
    let clipped_names: Vec<String> = clipped.iter().map(|r| name(r)).collect();
    // Every read the clipped run wrote is in the default too.
    for name in &clipped_names {
        assert!(plain_names.contains(name), "{name}");
    }
    // And the one it lost is the soft-clipped read itself.
    let lost: Vec<&String> = plain_names
        .iter()
        .filter(|name| !clipped_names.contains(name))
        .collect();
    assert_eq!(lost, vec!["r-clipped"]);
    assert!(!Arguments::default().use_softclipped_bases);
}

/// Refused by name and position, the tool having nothing to score it with.
#[test]
fn a_read_that_is_not_flow_based_is_refused() {
    let text = golden();
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix("error\tnot-flow-based\t"))
        .expect("its refusal");
    let (class, message) = row.split_once(':').expect("a class and a message");
    assert_eq!(class, "java.lang.IllegalArgumentException");
    assert!(
        message.starts_with("read must be flow based: "),
        "{message}"
    );
    // It names the read and where it aligned.
    assert!(message.contains("r-plain"), "{message}");
    assert!(message.contains("chr1:1000-1099"), "{message}");
    // The plain BAM's own read group carries no flow order, which is what makes it plain.
    let sam = section(&text, "sam", "plain");
    assert!(!sam.is_empty());
}
