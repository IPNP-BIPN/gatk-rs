//! Conformance for `AbstractConcordanceWalker.ConcordanceIterator` against GATK 4.6.2.0, compared
//! as the whole sequence of steps of every run.
//!
//! Golden from `tools/readfilter-conformance/ConcordanceWalkerDump.java`, two probe walkers run
//! through the real command line.
//!
//! # What this suite is for
//!
//!  * **a same-locus disagreement advances truth alone**, so one disagreement produces three steps;
//!  * **a filtered eval record is labelled by what truth has**;
//!  * **two of the five states are unreachable** for a walker dropping filtered records on both
//!    sides;
//!  * **and the order is the dictionary's**.

use gatk_engine::concordance_walker::{concordance, ConcordanceRecord, ConcordanceState};

fn golden() -> String {
    gatk_corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/concordance_walker.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// One record, as the iterator and the dump's own rendering see it.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Record {
    contig: String,
    start: i32,
    reference: String,
    alternates: Vec<String>,
    filters: String,
}

impl ConcordanceRecord for Record {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn is_filtered(&self) -> bool {
        self.filters != "PASS"
    }
}

impl Record {
    /// The dump's `render`: `<contig>:<start> <ref>/<alts> <filters>`.
    fn render(&self) -> String {
        format!(
            "{}:{} {}/{} {}",
            self.contig,
            self.start,
            self.reference,
            self.alternates.join(","),
            self.filters
        )
    }
}

/// The records of one input file, in file order.
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
                filters: match field[6] {
                    "." => "PASS".to_string(),
                    other => other.to_string(),
                },
            }
        })
        .collect()
}

/// The rule both probe walkers use: the same reference allele, and truth's first alternate among
/// eval's.
fn concordant(truth: &Record, eval: &Record) -> bool {
    truth.reference == eval.reference && eval.alternates.contains(&truth.alternates[0])
}

/// Every run: its two files and the filters the probe applied to each side.
fn setup(run: &str) -> (&'static str, &'static str, bool) {
    match run {
        // The base class's own: filtered truth dropped, every eval record kept.
        "default-filters" => ("truth", "eval", false),
        // EvaluateInfoFieldConcordance's: filtered records dropped on both sides.
        "filtered-dropped-both-sides" => ("truth", "eval", true),
        "across-contigs" => ("across-contigs", "second-contig-only", false),
        other => panic!("no setup for {other}"),
    }
}

const DICTIONARY: [&str; 2] = ["chr1", "chr2"];

/// The steps the golden recorded for one run, as `(state, truth, eval)` renderings.
fn expected(text: &str, run: &str) -> Vec<(String, String, String)> {
    rows(text, "state")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| (row[2].to_string(), row[3].to_string(), row[4].to_string()))
        .collect()
}

/// The same, replayed through the port.
fn mine(text: &str, run: &str) -> Vec<(String, String, String)> {
    let (truth_file, eval_file, drop_filtered_eval) = setup(run);
    let dictionary: Vec<String> = DICTIONARY.iter().map(|name| name.to_string()).collect();

    let truth: Vec<Record> = records(text, truth_file)
        .into_iter()
        .filter(|record| !record.is_filtered())
        .collect();
    let eval: Vec<Record> = records(text, eval_file)
        .into_iter()
        .filter(|record| !drop_filtered_eval || !record.is_filtered())
        .collect();

    concordance(&truth, &eval, &dictionary, concordant)
        .into_iter()
        .map(|step| {
            (
                step.state.name().to_string(),
                step.truth.map_or("-".to_string(), |at| truth[at].render()),
                step.eval.map_or("-".to_string(), |at| eval[at].render()),
            )
        })
        .collect()
}

#[test]
fn every_run_matches_the_golden_step_for_step() {
    let text = golden();
    for run in [
        "default-filters",
        "filtered-dropped-both-sides",
        "across-contigs",
    ] {
        assert_eq!(mine(&text, run), expected(&text, run), "the steps of {run}");
    }
}

#[test]
fn one_disagreement_produces_three_steps() {
    let text = golden();
    let steps = expected(&text, "default-filters");
    // The disagreement at 200 and what follows it.
    assert_eq!(
        &steps[5..8],
        &[
            (
                "FALSE_NEGATIVE".to_string(),
                "chr1:200 AT/A PASS".to_string(),
                "-".to_string()
            ),
            (
                "FALSE_POSITIVE".to_string(),
                "-".to_string(),
                "chr1:200 A/C PASS".to_string()
            ),
            (
                "FALSE_NEGATIVE".to_string(),
                "chr1:210 A/G PASS".to_string(),
                "-".to_string()
            ),
        ]
    );
}

#[test]
fn a_filtered_eval_record_is_labelled_by_what_truth_has() {
    let text = golden();
    let steps = expected(&text, "default-filters");
    let alone = steps
        .iter()
        .find(|(_, _, eval)| eval.starts_with("chr1:130"))
        .expect("the eval-only filtered record");
    assert_eq!(alone.0, ConcordanceState::FilteredTrueNegative.name());

    let paired = steps
        .iter()
        .find(|(_, truth, _)| truth.starts_with("chr1:140"))
        .expect("the filtered record at a truth locus");
    assert_eq!(paired.0, ConcordanceState::FilteredFalseNegative.name());
}

#[test]
fn dropping_filtered_records_on_both_sides_removes_two_states() {
    let text = golden();
    let states: Vec<String> = expected(&text, "filtered-dropped-both-sides")
        .into_iter()
        .map(|(state, _, _)| state)
        .collect();
    assert!(!states.iter().any(|state| state.starts_with("FILTERED_")));
    // And the truth record whose eval was filtered becomes a plain false negative.
    assert!(expected(&text, "filtered-dropped-both-sides").iter().any(
        |(state, truth, eval)| state == "FALSE_NEGATIVE"
            && truth.starts_with("chr1:140")
            && eval == "-"
    ));
}

#[test]
fn a_filtered_truth_record_never_reaches_the_iterator() {
    let text = golden();
    let at_300 = expected(&text, "default-filters")
        .into_iter()
        .find(|(_, _, eval)| eval.starts_with("chr1:300"))
        .expect("the eval record at 300");
    assert_eq!(at_300.0, "FALSE_POSITIVE");
    assert_eq!(at_300.1, "-", "its truth record was filtered away");
}
