//! Conformance for `Concordance`'s three annotated VCFs against GATK 4.6.2.0, compared as which
//! record reached which file with which STATUS, and as the headers the three are written against.
//!
//! Golden from `tools/readfilter-conformance/ConcordanceAnnotatedVcfsDump.java`.
//!
//! # What this suite is for
//!
//!  * **one step writes two different records**, a true positive being truth's record in `-tpfn`
//!    and eval's in `-tpfp`;
//!  * **`-tpfn` is written against the truth header**, so its sample column is the truth file's;
//!  * **a filtered false negative is `FFN` in both files it reaches**, and keeps its FILTER;
//!  * **the two header-building orders reach nothing**, the two runs asking for one output each
//!    writing the same header lines.
//!
//! The traversal itself is [`gatk_engine::concordance_walker`], measured with its own suite: this
//! one replays it to decide which step is in which state.

use gatk_corpus as corpus;
use gatk_engine::concordance_walker::{concordance, ConcordanceRecord};
use gatk_tools::concordance::{
    truth_variant_filter, variants_at_same_locus_are_concordant, writes, AnnotatedVcf, Side,
    TRUTH_STATUS_HEADER_LINE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/concordance_annotated_vcfs.txt.gz"),
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

/// One record, decoded as far as the walker and the routing look at it.
struct Record {
    contig: String,
    start: i32,
    identifier: String,
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
                identifier: field[2].to_string(),
                reference: field[3].to_string(),
                alternates: field[4].split(',').map(|alt| alt.to_string()).collect(),
                filtered: !matches!(field[6], "." | "PASS"),
            }
        })
        .collect()
}

/// Every line of one written VCF, in file order.
fn written(text: &str, label: &str) -> Vec<String> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let rest = line
                .strip_prefix("vcfline\t")
                .or_else(|| line.strip_prefix("commandline\t"))?;
            let (row, body) = rest.split_once('\t')?;
            (row == label).then(|| unescape(body))
        })
        .collect()
}

/// The `ID` and the `STATUS` of every record line of one written VCF.
fn identified(text: &str, label: &str) -> Vec<(String, String)> {
    written(text, label)
        .into_iter()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let status = field[7]
                .split(';')
                .find_map(|entry| entry.strip_prefix("STATUS="))
                .expect("every written record is annotated")
                .to_string();
            (field[2].to_string(), status)
        })
        .collect()
}

const DICTIONARY: [&str; 1] = ["chr1"];

/// What the traversal says each file should hold: the source record's ID and the STATUS beside it.
fn routed(text: &str, file: AnnotatedVcf) -> Vec<(String, String)> {
    let dictionary: Vec<String> = DICTIONARY.iter().map(|name| name.to_string()).collect();
    let truth: Vec<Record> = records(&input(text, "truth"))
        .into_iter()
        .filter(|record| truth_variant_filter(record.is_filtered(), false))
        .collect();
    let eval = records(&input(text, "eval"));

    let mut lines = Vec::new();
    for step in concordance(&truth, &eval, &dictionary, |truth, eval| {
        variants_at_same_locus_are_concordant(
            &truth.reference,
            &truth.alternates,
            &eval.reference,
            &eval.alternates,
        )
    }) {
        for (target, side) in writes(step.state) {
            if target != file {
                continue;
            }
            let record = match side {
                Side::Truth => &truth[step.truth.expect("a truth-side write has one")],
                Side::Eval => &eval[step.eval.expect("an eval-side write has one")],
            };
            lines.push((
                record.identifier.clone(),
                step.state.abbreviation().to_string(),
            ));
        }
    }
    lines
}

#[test]
fn one_step_writes_two_different_records() {
    let text = golden();
    // chr1:100 is a true positive, and the two files carry the two sides of it under one status.
    assert!(identified(&text, "all-three-tpfn").contains(&("truth_tp".into(), "TP".into())));
    assert!(identified(&text, "all-three-tpfp").contains(&("eval_tp".into(), "TP".into())));
    // And the two lines agree on nothing else.
    let truth_line = written(&text, "all-three-tpfn")
        .into_iter()
        .find(|line| line.contains("truth_tp"))
        .expect("the truth side");
    let eval_line = written(&text, "all-three-tpfp")
        .into_iter()
        .find(|line| line.contains("eval_tp"))
        .expect("the eval side");
    assert!(truth_line.contains("TRUTHONLY=1") && truth_line.contains("\t11\t"));
    assert!(eval_line.contains("EVALONLY=11") && eval_line.contains("\t44\t"));
}

#[test]
fn every_file_holds_what_the_routing_says_it_does() {
    let text = golden();
    for (file, label) in [
        (
            AnnotatedVcf::TruePositivesAndFalseNegatives,
            "all-three-tpfn",
        ),
        (
            AnnotatedVcf::TruePositivesAndFalsePositives,
            "all-three-tpfp",
        ),
        (
            AnnotatedVcf::FilteredTrueNegativesAndFalseNegatives,
            "all-three-ftnfn",
        ),
    ] {
        assert_eq!(routed(&text, file), identified(&text, label), "{label}");
    }
}

#[test]
fn a_filtered_false_negative_is_ffn_in_both_files_and_keeps_its_filter() {
    let text = golden();
    assert!(identified(&text, "all-three-tpfn").contains(&("truth_ffn".into(), "FFN".into())));
    assert!(identified(&text, "all-three-ftnfn").contains(&("eval_ffn".into(), "FFN".into())));
    // The FILTER column of the eval copy is the one that made it a filtered state.
    let filtered = written(&text, "all-three-ftnfn")
        .into_iter()
        .find(|line| line.contains("eval_ffn"))
        .expect("the eval copy");
    assert_eq!(filtered.split('\t').nth(6), Some("weak"));
}

#[test]
fn only_the_first_file_is_written_against_the_truth_header() {
    let text = golden();
    let samples = |label: &str| -> String {
        written(&text, label)
            .into_iter()
            .find(|line| line.starts_with("#CHROM"))
            .expect("a column header")
            .split('\t')
            .nth(9)
            .expect("one sample")
            .to_string()
    };
    assert_eq!(samples("all-three-tpfn"), "truthsample");
    assert_eq!(samples("all-three-tpfp"), "evalsample");
    assert_eq!(samples("all-three-ftnfn"), "evalsample");
    assert_eq!(
        AnnotatedVcf::TruePositivesAndFalseNegatives.header(),
        Side::Truth
    );
    // The truth file declares the truth file's INFO key and not the eval file's.
    let tpfn = written(&text, "all-three-tpfn");
    assert!(tpfn.iter().any(|line| line.contains("ID=TRUTHONLY")));
    assert!(!tpfn.iter().any(|line| line.contains("ID=EVALONLY")));
}

#[test]
fn the_two_header_building_orders_reach_nothing() {
    let text = golden();
    let header = |label: &str| -> Vec<String> {
        written(&text, label)
            .into_iter()
            .filter(|line| line.starts_with("##") && !line.starts_with("##GATKCommandLine"))
            .collect()
    };
    // -tpfp adds the default lines and then STATUS; -ftnfn adds STATUS and then the default lines.
    assert_eq!(header("tpfp-alone-tpfp"), header("ftnfn-alone-ftnfn"));
    assert!(header("tpfp-alone-tpfp").contains(&TRUTH_STATUS_HEADER_LINE.to_string()));
    // The default tool lines do reach the file, unlike the sibling tool's.
    assert!(header("tpfp-alone-tpfp").contains(&"##source=Concordance".to_string()));
}
