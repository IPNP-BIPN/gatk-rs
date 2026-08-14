//! Conformance for `Concordance`'s filter-analysis table against GATK 4.6.2.0, compared as the
//! whole table, the summary beside it, and both crashes.
//!
//! Golden from `tools/readfilter-conformance/ConcordanceFilterAnalysisDump.java`.
//!
//! # What this suite is for
//!
//!  * **the counting happens without the flag that asks for it**, the guard being
//!    `(flag && FTN) || FFN`, which is only visible as a crash: the same undeclared filter is fatal
//!    at a truth locus with no `--filter-analysis` and harmless on a record standing alone;
//!  * **a record with two filters increments neither unique column**;
//!  * **a declared filter nothing carries still gets a row**;
//!  * **and the row order is a `HashMap`'s**, not the header's.
//!
//! The traversal itself is [`gatk_engine::concordance_walker`], measured with its own suite.

use gatk_corpus as corpus;
use gatk_engine::concordance_walker::{concordance, ConcordanceRecord};
use gatk_tools::concordance::{
    truth_variant_filter, variants_at_same_locus_are_concordant, FilterAnalysis,
    FilterAnalysisError, Summary,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/concordance_filter_analysis.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// One record, decoded as far as the walker and the two tables look at it.
struct Record {
    contig: String,
    start: i32,
    reference: String,
    alternates: Vec<String>,
    /// Every FILTER the record carries, which is empty for `PASS` and for `.`.
    filters: Vec<String>,
}

impl ConcordanceRecord for Record {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn is_filtered(&self) -> bool {
        !self.filters.is_empty()
    }
}

impl Record {
    /// `isSNP()`, which is all this fixture carries.
    fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternates.iter().all(|alt| alt.len() == 1)
    }
}

/// The whole text of one input.
fn input(text: &str, label: &str) -> String {
    unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    )
}

fn records(whole: &str) -> Vec<Record> {
    whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            Record {
                contig: field[0].to_string(),
                start: field[1].parse().expect("a position"),
                reference: field[3].to_string(),
                alternates: field[4].split(',').map(|alt| alt.to_string()).collect(),
                filters: match field[6] {
                    "." | "PASS" => Vec::new(),
                    list => list.split(';').map(|name| name.to_string()).collect(),
                },
            }
        })
        .collect()
}

/// `evalHeader.getFilterLines()`, in the order the header declares them.
fn declared_filters(whole: &str) -> Vec<String> {
    whole
        .lines()
        .filter_map(|line| line.strip_prefix("##FILTER=<ID="))
        .filter_map(|rest| rest.split(',').next())
        .map(|id| id.to_string())
        .collect()
}

/// Every run: its truth file, its eval file, and whether `--filter-analysis` was given.
fn setup(run: &str) -> (&'static str, &'static str, bool) {
    match run {
        "baseline" => ("truth", "eval", true),
        "ftn-undeclared-no-flag" => ("ghost-truth", "ghost-eval", false),
        "ftn-undeclared-with-flag" => ("ghost-truth", "ghost-eval", true),
        "ffn-undeclared-no-flag" => ("ghost-locus-truth", "ghost-locus-eval", false),
        other => panic!("no setup for {other}"),
    }
}

const DICTIONARY: [&str; 1] = ["chr1"];

/// The whole traversal, and both outputs of it.
fn run(text: &str, run: &str) -> Result<(String, Option<String>), FilterAnalysisError> {
    let (truth_file, eval_file, requested) = setup(run);
    let dictionary: Vec<String> = DICTIONARY.iter().map(|name| name.to_string()).collect();
    let eval_text = input(text, eval_file);
    let truth: Vec<Record> = records(&input(text, truth_file))
        .into_iter()
        // No symbolic record in this fixture, so the truth filter is the filtered test alone.
        .filter(|record| truth_variant_filter(record.is_filtered(), false))
        .collect();
    let eval = records(&eval_text);

    let mut summary = Summary::default();
    let mut analysis = FilterAnalysis::new(&declared_filters(&eval_text));
    for step in concordance(&truth, &eval, &dictionary, |truth, eval| {
        variants_at_same_locus_are_concordant(
            &truth.reference,
            &truth.alternates,
            &eval.reference,
            &eval.alternates,
        )
    }) {
        let record = match step.truth {
            Some(index) => &truth[index],
            None => &eval[step.eval.expect("a step has one side or the other")],
        };
        summary.add(step.state, record.is_snp());
        let filters = match step.eval {
            Some(index) => eval[index].filters.clone(),
            None => Vec::new(),
        };
        analysis.apply(step.state, &filters, requested)?;
    }
    let table = if requested {
        Some(analysis.table()?)
    } else {
        None
    };
    Ok((summary.table(), table))
}

fn expected(text: &str, kind: &str, label: &str) -> Option<String> {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == label)
        .map(|row| unescape(row[1]))
}

#[test]
fn every_output_matches_the_golden_byte_for_byte() {
    let text = golden();
    for label in ["baseline", "ftn-undeclared-no-flag"] {
        let (summary, table) = run(&text, label).expect("this run reaches the end");
        assert_eq!(
            Some(summary),
            expected(&text, "summary", label),
            "the summary of {label}"
        );
        assert_eq!(
            table,
            expected(&text, "table", label),
            "the table of {label}"
        );
    }
}

#[test]
fn the_counting_happens_without_the_flag_that_asks_for_it() {
    let text = golden();
    // A filtered false negative with no flag: the reference crashes, and it wrote no table at all.
    let error = run(&text, "ffn-undeclared-no-flag").expect_err("the undeclared filter is a null");
    assert!(expected(&text, "table", "ffn-undeclared-no-flag").is_none());
    // The same undeclared filter on a record standing alone, with no flag, reaches the end.
    assert!(run(&text, "ftn-undeclared-no-flag").is_ok());
    assert!(expected(&text, "summary", "ftn-undeclared-no-flag").is_some());
    // And with the flag, the very same file crashes.
    let with_flag = run(&text, "ftn-undeclared-with-flag").expect_err("now it is reached");

    for (label, error) in [
        ("ffn-undeclared-no-flag", error),
        ("ftn-undeclared-with-flag", with_flag),
    ] {
        let row = rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no refusal {label}"));
        let (class, message) = row[1].split_once(':').expect("class and message");
        assert_eq!(error.class(), class);
        assert_eq!(error.message(), message);
    }
}

#[test]
fn a_record_with_two_filters_increments_neither_unique_column() {
    let text = golden();
    let table = expected(&text, "table", "baseline").expect("the baseline table");
    // `weak` is alone on one false negative and one true negative and shares the other two.
    assert!(
        table.lines().any(|line| line == "weak\t1\t2\t1\t1"),
        "{table}"
    );
    // `shallow` is never alone.
    assert!(
        table.lines().any(|line| line == "shallow\t1\t1\t0\t0"),
        "{table}"
    );
}

#[test]
fn every_declared_filter_has_a_row_in_the_hash_maps_order() {
    let text = golden();
    let eval = input(&text, "eval");
    assert_eq!(
        declared_filters(&eval),
        vec!["weak", "shallow", "noisy", "unused"]
    );
    let table = expected(&text, "table", "baseline").expect("the baseline table");
    let order: Vec<&str> = table
        .lines()
        .skip(1)
        .filter_map(|line| line.split('\t').next())
        .collect();
    assert_eq!(order, vec!["shallow", "unused", "noisy", "weak"]);
    // The one nothing carries is a row of zeroes rather than an absence.
    assert!(table.lines().any(|line| line == "unused\t0\t0\t0\t0"));
}

#[test]
fn a_filtered_state_is_still_counted_by_the_summary() {
    let text = golden();
    let summary = expected(&text, "summary", "baseline").expect("the baseline summary");
    // Two filtered false negatives in the FN column, and the two filtered true negatives nowhere.
    assert_eq!(
        summary.lines().nth(1).expect("the SNP row"),
        "SNP\t1\t0\t2\t0.333\t1.0"
    );
}
