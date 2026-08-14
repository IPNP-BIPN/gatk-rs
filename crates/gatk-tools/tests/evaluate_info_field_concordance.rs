//! Conformance for `EvaluateInfoFieldConcordance` against GATK 4.6.2.0, compared as the whole
//! summary table of every run and as both refusals.
//!
//! Golden from `tools/readfilter-conformance/EvaluateInfoFieldConcordanceDump.java`.
//!
//! # What this suite is for
//!
//!  * **a true positive whose record lacks the key is counted but contributes no delta**;
//!  * **the variance is the cancelling form**, and reports `0.39999999999999997` where the algebra
//!    says 0.4;
//!  * **an empty bucket is two NaN columns**;
//!  * **and a missing key in either header is refused before any record**.
//!
//! The traversal itself is [`gatk_engine::concordance_walker`], measured with its own suite: this
//! one replays it to decide which records reach the arithmetic.

use gatk_corpus as corpus;
use gatk_engine::concordance_walker::{concordance, ConcordanceRecord, ConcordanceState};
use gatk_tools::evaluate_info_field_concordance::{check_keys, Concordance, EvalType};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/evaluate_info_field_concordance.txt.gz"),
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

/// One record, decoded as far as the walker and the arithmetic look at it.
struct Record {
    contig: String,
    start: i32,
    reference: String,
    alternates: Vec<String>,
    filtered: bool,
    info: Vec<(String, String)>,
}

impl ConcordanceRecord for Record {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn is_filtered(&self) -> bool {
        self.filtered
    }
}

impl Record {
    fn attribute(&self, key: &str) -> Option<f64> {
        self.info
            .iter()
            .find(|(name, _)| name == key)
            .and_then(|(_, value)| value.parse().ok())
    }

    /// `isSNP()` and `isIndel()` for a biallelic record, which is all this fixture carries.
    fn eval_type(&self) -> EvalType {
        let alternate = &self.alternates[0];
        if alternate.len() == self.reference.len() {
            if self.reference.len() == 1 {
                EvalType::Snp
            } else {
                EvalType::Other
            }
        } else {
            EvalType::Indel
        }
    }
}

/// The records of one input, minus the ones this walker's filters drop on both sides.
fn records(text: &str, label: &str) -> Vec<Record> {
    let whole = unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    );
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
                filtered: !matches!(field[6], "." | "PASS"),
                info: field[7]
                    .split(';')
                    .filter_map(|entry| entry.split_once('='))
                    .map(|(key, value)| (key.to_string(), value.to_string()))
                    .collect(),
            }
        })
        // Both of this walker's filters: filtered records dropped on either side.
        .filter(|record| !record.is_filtered())
        .collect()
}

/// The rule this tool declares: the same reference allele, and truth's first alternate among eval's.
fn concordant(truth: &Record, eval: &Record) -> bool {
    truth.reference == eval.reference && eval.alternates.contains(&truth.alternates[0])
}

/// Every run that produced a table: its two files and the two keys.
fn setup(run: &str) -> (&'static str, &'static str, &'static str, &'static str) {
    match run {
        "baseline" => ("truth", "eval", "SCORE", "SCORE"),
        "spread" => ("spread", "spread-eval", "SCORE", "SCORE"),
        "nothing-agrees" => ("nothing-agrees", "nothing-agrees-eval", "SCORE", "SCORE"),
        "different-keys" => ("truth", "eval", "SCORE", "OTHER"),
        other => panic!("no setup for {other}"),
    }
}

const DICTIONARY: [&str; 1] = ["chr1"];

/// The whole traversal and the arithmetic on top of it.
fn table(text: &str, run: &str) -> String {
    let (truth_file, eval_file, eval_key, truth_key) = setup(run);
    let dictionary: Vec<String> = DICTIONARY.iter().map(|name| name.to_string()).collect();
    let truth = records(text, truth_file);
    let eval = records(text, eval_file);

    let mut concordance_summary = Concordance::default();
    for step in concordance(&truth, &eval, &dictionary, concordant) {
        if step.state != ConcordanceState::TruePositive {
            continue;
        }
        let eval_record = &eval[step.eval.expect("a true positive has both")];
        let truth_record = &truth[step.truth.expect("a true positive has both")];
        concordance_summary.add(
            eval_record.eval_type(),
            eval_record.attribute(eval_key),
            truth_record.attribute(truth_key),
        );
    }
    concordance_summary.table(eval_key, truth_key)
}

fn expected(text: &str, run: &str) -> String {
    unescape(
        rows(text, "table")
            .into_iter()
            .find(|row| row[0] == run)
            .unwrap_or_else(|| panic!("no table {run}"))[1],
    )
}

#[test]
fn every_table_matches_the_golden_byte_for_byte() {
    let text = golden();
    let tables = rows(&text, "table");
    assert_eq!(tables.len(), 4, "four of the six runs write a table");
    for row in tables {
        let run = row[0];
        assert_eq!(
            table(&text, run),
            expected(&text, run),
            "the table of {run}"
        );
    }
}

#[test]
fn the_variance_does_not_quite_cancel() {
    let text = golden();
    let snp = expected(&text, "baseline")
        .lines()
        .nth(1)
        .expect("the SNP row")
        .to_string();
    assert!(
        snp.ends_with("\t0.8\t0.39999999999999997"),
        "the algebra says 0.4: {snp}"
    );
}

#[test]
fn an_empty_bucket_is_two_nan_columns() {
    let text = golden();
    // The spread run has no indel among its true positives.
    let indel = expected(&text, "spread")
        .lines()
        .nth(2)
        .expect("the INDEL row")
        .to_string();
    assert!(indel.ends_with("\tNaN\tNaN"), "{indel}");
    // And a run where nothing agrees has two of them.
    let nothing = expected(&text, "nothing-agrees");
    assert_eq!(
        nothing
            .lines()
            .filter(|line| line.ends_with("\tNaN\tNaN"))
            .count(),
        2
    );
}

#[test]
fn both_refusals_carry_the_references_class_and_words() {
    let text = golden();
    let refusal = |label: &str| -> (String, String) {
        let row = rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no refusal {label}"));
        let (class, message) = row[1].split_once(':').expect("class and message");
        (class.to_string(), message.to_string())
    };

    let (class, message) = refusal("missing-eval-key");
    let missing_eval = check_keys(
        false,
        "SCORE",
        "evaluateinfofieldconcordance-dump/no-key.vcf",
        true,
        "SCORE",
        "evaluateinfofieldconcordance-dump/truth.vcf",
    )
    .expect_err("the eval header lacks it");
    assert_eq!(missing_eval.class(), class);
    assert_eq!(missing_eval.message(), message);

    let (class, message) = refusal("missing-truth-key");
    let missing_truth = check_keys(
        true,
        "SCORE",
        "evaluateinfofieldconcordance-dump/eval.vcf",
        false,
        "SCORE",
        "evaluateinfofieldconcordance-dump/no-key.vcf",
    )
    .expect_err("the truth header lacks it");
    assert_eq!(missing_truth.class(), class);
    assert_eq!(missing_truth.message(), message);
}
