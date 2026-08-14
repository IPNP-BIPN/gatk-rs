//! Conformance for `Concordance`'s summary table against GATK 4.6.2.0, compared as the whole table
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/ConcordanceDump.java`.
//!
//! # What this suite is for
//!
//!  * **an empty callset reports a recall of `0.0`**, the rounding turning the division's NaN into
//!    zero, where [`gatk_tools::evaluate_info_field_concordance`] writes `NaN` for the same `0/0`;
//!  * **a filtered eval record alone leaves no trace**, while a filtered one at a truth locus
//!    reaches the FN column;
//!  * **an MNP and a symbolic record share the INDEL row**, and the symbolic record survives on the
//!    eval side only;
//!  * **agreement needs the same number of alternates but only truth's first**.
//!
//! The traversal itself is [`gatk_engine::concordance_walker`], measured with its own suite: this
//! one replays it to decide which state each step is in.

use gatk_corpus as corpus;
use gatk_engine::concordance_walker::{concordance, ConcordanceRecord};
use gatk_tools::concordance::{
    truth_variant_filter, variants_at_same_locus_are_concordant, Summary,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/concordance.txt.gz"),
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

/// One record, decoded as far as the walker and the summary look at it.
struct Record {
    contig: String,
    start: i32,
    reference: String,
    alternates: Vec<String>,
    filtered: bool,
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
    /// `hasSymbolicAlleles() || isSV()`, which for this fixture is the angle brackets.
    fn is_symbolic_or_sv(&self) -> bool {
        self.alternates.iter().any(|alt| alt.starts_with('<'))
    }

    /// `isSNP()`: one base against one base, for every alternate.
    fn is_snp(&self) -> bool {
        self.reference.len() == 1
            && !self.alternates.is_empty()
            && self
                .alternates
                .iter()
                .all(|alt| alt.len() == 1 && !alt.starts_with('<'))
    }
}

/// The records of one input, in file order and undropped.
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
            }
        })
        .collect()
}

/// Every run: its truth file and its eval file.
fn setup(run: &str) -> (&'static str, &'static str) {
    match run {
        "baseline" => ("truth", "eval"),
        "filtered" => ("filtered-truth", "filtered-eval"),
        "empty" => ("empty-truth", "empty-eval"),
        "types" => ("types-truth", "types-eval"),
        "alleles" => ("alleles-truth", "alleles-eval"),
        other => panic!("no setup for {other}"),
    }
}

const DICTIONARY: [&str; 1] = ["chr1"];

/// The whole traversal and the table on top of it.
fn table(text: &str, run: &str) -> String {
    let (truth_file, eval_file) = setup(run);
    let dictionary: Vec<String> = DICTIONARY.iter().map(|name| name.to_string()).collect();
    // The truth filter drops filtered and symbolic records; the eval side has no filter at all.
    let truth: Vec<Record> = records(text, truth_file)
        .into_iter()
        .filter(|record| truth_variant_filter(record.is_filtered(), record.is_symbolic_or_sv()))
        .collect();
    let eval = records(text, eval_file);

    let mut summary = Summary::default();
    for step in concordance(&truth, &eval, &dictionary, |truth, eval| {
        variants_at_same_locus_are_concordant(
            &truth.reference,
            &truth.alternates,
            &eval.reference,
            &eval.alternates,
        )
    }) {
        // `getTruthIfPresentElseEval()`.
        let record = match step.truth {
            Some(index) => &truth[index],
            None => &eval[step.eval.expect("a step has one side or the other")],
        };
        summary.add(step.state, record.is_snp());
    }
    summary.table()
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
    assert_eq!(tables.len(), 5, "one table per run");
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
fn an_empty_callset_reports_zero_rather_than_nan() {
    let text = golden();
    let empty = expected(&text, "empty");
    for row in empty.lines().skip(1) {
        assert!(row.ends_with("\t0\t0\t0\t0.0\t0.0"), "{row}");
    }
    // The division itself is NaN: only the rounding makes it zero.
    assert!(Summary::default().snp.sensitivity().is_nan());
}

#[test]
fn a_filtered_eval_record_alone_leaves_no_trace() {
    let text = golden();
    // Three eval records, one of them unmatched and filtered, and the precision is still one.
    assert_eq!(records(&text, "filtered-eval").len(), 3);
    let snp = expected(&text, "filtered")
        .lines()
        .nth(1)
        .expect("the SNP row")
        .to_string();
    assert_eq!(snp, "SNP\t1\t0\t1\t0.5\t1.0", "{snp}");
}

#[test]
fn an_mnp_and_a_symbolic_record_share_the_indel_row() {
    let text = golden();
    let types = expected(&text, "types");
    // The multi-allelic SNP is the whole SNP row; the MNP and the symbolic record are the other.
    assert_eq!(
        types.lines().nth(1).expect("the SNP row"),
        "SNP\t1\t0\t0\t1.0\t1.0"
    );
    // The symbolic record is in both files and comes out a false positive, the truth side having
    // dropped it.
    assert_eq!(
        types.lines().nth(2).expect("the INDEL row"),
        "INDEL\t1\t1\t0\t1.0\t0.5"
    );
}

#[test]
fn the_rounding_is_half_up_after_an_ulp_is_added() {
    let text = golden();
    let snp = expected(&text, "baseline")
        .lines()
        .nth(1)
        .expect("the SNP row")
        .to_string();
    assert!(snp.ends_with("\t0.667\t0.667"), "two thirds, twice: {snp}");
}
