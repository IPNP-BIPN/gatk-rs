//! Conformance for `VcfToIntervalList` against Picard 3.4.0, compared as the whole output file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/VcfToIntervalListDump.java`, which carries each run's
//! input as well as its output.
//!
//! # What this suite is for
//!
//!  * **the merging is a stream and never sorts**, so an unsorted file comes out unsorted and two
//!    records that would have merged do not when a third stands between them;
//!  * **abutting intervals merge and intervals one base apart do not**;
//!  * **an unnamed record is `interval-<n>` counted only over unnamed records**, and counted after
//!    the filtering, so `INCLUDE_FILTERED` renumbers everything after a filtered one;
//!  * **`INCLUDE_FILTERED` changes which intervals merge**, a filtered record between two others
//!    being the bridge that joins them;
//!  * **a `PASS` filter is not a filter**;
//!  * **and a header with no contig lines is a null dictionary**, which the codec walks into.
//!
//! # What the golden does not pin down
//!
//! `CONCAT_ALL` rebuilds a merged interval from the group's minimum start, while `USE_FIRST` keeps
//! the start of the group's first member. The two differ only when a record overlaps one that
//! started later, which no run here produces. The port follows the Java on both branches.

use gatk_corpus as corpus;
use gatk_tools::vcf_to_interval_list::{convert, IdMethod};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/vcf_to_interval_list.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

/// The arguments of each run.
fn arguments(label: &str) -> (bool, IdMethod) {
    match label {
        "defaults" | "unsorted" | "no-records" | "no-contigs" => (false, IdMethod::ConcatAll),
        "include-filtered" => (true, IdMethod::ConcatAll),
        "use-first" => (false, IdMethod::UseFirst),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_interval_list_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "defaults",
        "include-filtered",
        "unsorted",
        "no-records",
        "use-first",
    ] {
        let (include_filtered, method) = arguments(label);
        let ours = convert(&value(&text, "input", label), include_filtered, method)
            .expect("a run the tool allows");
        assert_eq!(ours, value(&text, "list", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 5, "the golden's outputs");
}

#[test]
fn a_header_with_no_contigs_walks_into_a_null_dictionary() {
    let text = golden();
    let (include_filtered, method) = arguments("no-contigs");
    let error = convert(
        &value(&text, "input", "no-contigs"),
        include_filtered,
        method,
    )
    .expect_err("the null dictionary");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-contigs")
    );
}

/// The two behaviours the whole-file comparison would let a wrong port hide: an interval named
/// after a record that is not the one before it, and a merge that a filtered record made.
#[test]
fn the_numbering_and_the_bridge_move_with_include_filtered() {
    let text = golden();
    let names = |label: &str| -> Vec<String> {
        value(&text, "list", label)
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(|line| line.rsplit('\t').next().expect("a name column").to_string())
            .collect()
    };
    let default_names = names("defaults");
    let filtered_names = names("include-filtered");
    // The unnamed record at 400 is the fourth unnamed one by default and the fifth when the
    // filtered record before it is kept.
    assert!(default_names.contains(&"interval-4|rs9;rs10".to_string()));
    assert!(filtered_names.contains(&"interval-5|rs9;rs10".to_string()));
    // And the filtered record at 81 joins the two around it, which are separate by default.
    assert!(default_names.contains(&"rs7".to_string()));
    assert!(filtered_names.contains(&"rs7|bridge|rs8".to_string()));
}
