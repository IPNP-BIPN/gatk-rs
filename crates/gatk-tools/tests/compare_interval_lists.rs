//! Conformance for `CompareIntervalLists` against GATK 4.6.2.0, compared as the line every pair of
//! files leaves behind.
//!
//! Golden from `tools/readfilter-conformance/CompareIntervalListsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the walk is not symmetric**: a master that outlasts the test dies with a
//!    `NoSuchElementException` while a test that outlasts the master is reported;
//!  * **a test interval wider than the master's is equal**, and the same pair reversed throws;
//!  * **each file is ALL-merged first**, so abutting or overlapping intervals in one file and the
//!    single interval they add up to in the other are equal;
//!  * **and equality is printed while inequality is thrown**.
//!
//! # What is compared, and what is not
//!
//! The three rows whose answer is a parse refusal are not replayed through this port: they are the
//! interval argument collection's, already pinned by its own suite, and one of them names a file by
//! URI. They are asserted as the constants the golden carries, so that a change in either message
//! is still a failure here.

use gatk_corpus as corpus;
use gatk_engine::interval::{MergingRule, SimpleInterval};
use gatk_engine::interval_args::{load_intervals, SetRule};
use gatk_tools::compare_interval_lists::{equate_intervals, Comparison};
use htsjdk_bam::header::{SamHeader, SequenceRecord};

const CONTIG_LENGTH: i32 = 240;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/compare_interval_lists.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, label: &str) -> String {
    let prefix = format!("compare\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {label}")),
    )
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    for contig in ["chr1", "chr2"] {
        header
            .sequences
            .push(SequenceRecord::new(contig, CONTIG_LENGTH));
    }
    header
}

/// One file, sorted and ALL-merged the way `getGenomeLocs` does it.
fn file(queries: &[&str]) -> Vec<SimpleInterval> {
    let queries: Vec<String> = queries.iter().map(|q| (*q).to_string()).collect();
    let (loaded, _) = load_intervals(&queries, &header(), SetRule::Union, MergingRule::All, 0)
        .expect("intervals the dictionary allows");
    loaded
}

/// The two files of each labelled run, as the dump wrote them.
fn run(label: &str) -> (Vec<&'static str>, Vec<&'static str>) {
    match label {
        "equal" => (vec!["chr1:1-100"], vec!["chr1:1-100"]),
        "equal-two-contigs" => (
            vec!["chr1:1-100", "chr2:1-50"],
            vec!["chr1:1-100", "chr2:1-50"],
        ),
        "abutting-merged" => (vec!["chr1:1-50", "chr1:51-100"], vec!["chr1:1-100"]),
        "overlapping-merged" => (vec!["chr1:1-60", "chr1:40-100"], vec!["chr1:1-100"]),
        "test-wider" => (vec!["chr1:10-20"], vec!["chr1:1-100"]),
        "master-wider" => (vec!["chr1:1-100"], vec!["chr1:10-20"]),
        "test-longer" => (vec!["chr1:1-100"], vec!["chr1:1-100", "chr2:1-50"]),
        "master-longer" => (vec!["chr1:1-100", "chr2:1-50"], vec!["chr1:1-100"]),
        "disjoint" => (vec!["chr1:1-100"], vec!["chr1:201-240"]),
        "different-contigs" => (vec!["chr1:1-100"], vec!["chr2:1-100"]),
        "test-inside-master" => (vec!["chr1:1-100"], vec!["chr1:40-60"]),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_pair_answers_what_the_reference_answers() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "equal",
        "equal-two-contigs",
        "abutting-merged",
        "overlapping-merged",
        "test-wider",
        "master-wider",
        "test-longer",
        "master-longer",
        "disjoint",
        "different-contigs",
        "test-inside-master",
    ] {
        let (first, second) = run(label);
        let ours = equate_intervals(&file(&first), &file(&second));
        assert_eq!(ours.line(), value(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 11, "the golden's comparable pairs");

    // The asymmetry, stated as itself rather than inferred from two rows.
    assert_eq!(
        equate_intervals(&file(&["chr1:10-20"]), &file(&["chr1:1-100"])),
        Comparison::Equal
    );
    assert_eq!(
        equate_intervals(&file(&["chr1:1-100"]), &file(&["chr1:10-20"])),
        Comparison::TestExhausted
    );
}

#[test]
fn the_three_parse_refusals_are_the_argument_collections() {
    let text = golden();
    for (label, message) in [
        (
            "empty-master",
            "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile:File \
             file:<masked>/empty-master-1.list is malformed: It contains no intervals.",
        ),
        (
            "empty-test",
            "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile:File \
             file:<masked>/empty-test-2.list is malformed: It contains no intervals.",
        ),
        (
            "off-the-end",
            "org.broadinstitute.hellbender.exceptions.UserException$MalformedGenomeLoc:Badly \
             formed genome unclippedLoc: Parameters to GenomeLocParser are incorrect:The genome \
             loc coordinates 1-1000 exceed the contig size (240)",
        ),
    ] {
        assert_eq!(value(&text, label), message, "{label}");
    }
}
