//! Conformance for `MultiFeatureWalker` against GATK 4.6.2.0, through `ExampleMultiFeatureWalker`,
//! which prints every feature it is handed so the merge order is the output.
//!
//! Golden from `tools/readfilter-conformance/MultiFeatureWalkerDump.java`.
//!
//! # What this suite is for
//!
//!  * **the heap not being stable**, so three files each holding the same interval come out first,
//!    third, second;
//!  * **the comparison being contig index, then start, then end**, the index coming from the
//!    dictionary;
//!  * **a contig the dictionary does not name being misdiagnosed** as an unsorted input;
//!  * **the sort check firing late and against the same input**;
//!  * **the larger dictionary winning**, and the two refusals when the smaller is not a subset of
//!    it in the same order;
//!  * **and no dictionary at all being refused.**
//!
//! Two runs of the golden are deliberately not compared here, because neither is this class's
//! doing: `unknown-contig` is refused by the engine's own dictionary comparison before the walk
//! starts, and `unsorted-index` is refused by `IndexFeatureFile`. They are measured so that the
//! reachability of the walker's own checks is on the record.

use gatk_corpus as corpus;
use gatk_engine::multi_feature_walker::{
    choose_dictionary, merge, DictSource, Located, WalkerError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/multi_feature_walker.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// `features\t<label>=...`, the whole of what the tool printed.
fn features(text: &str, label: &str) -> String {
    let prefix = format!("features\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries features/{label}")),
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

/// One `.rd.txt` row as the walker holds it: the file is zero-based, the walker one-based, and
/// `DepthEvidence.toString()` prints the one-based start.
fn bin(contig: &str, start: i32, end: i32, count: i32) -> Located {
    Located {
        contig: contig.to_string(),
        start: start + 1,
        end,
        text: format!("{contig}\t{}\t{end}\t{count}", start + 1),
    }
}

fn dictionary(source: &str, contigs: &[&str]) -> DictSource {
    DictSource {
        contigs: contigs.iter().map(|name| name.to_string()).collect(),
        source: source.to_string(),
    }
}

fn printed(features: &[Located]) -> String {
    features
        .iter()
        .map(|feature| format!("{}\n", feature.text))
        .collect()
}

fn full() -> DictSource {
    dictionary("sequence-dictionary", &["chr1", "chr2", "chr3"])
}

#[test]
fn every_merged_stream_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, inputs) in [
        (
            "interleaved",
            vec![
                vec![bin("chr1", 100, 200, 1), bin("chr1", 300, 400, 3)],
                vec![bin("chr1", 200, 300, 2), bin("chr1", 400, 500, 4)],
            ],
        ),
        (
            "tie-two",
            vec![
                vec![bin("chr1", 100, 200, 1)],
                vec![bin("chr1", 100, 200, 2)],
            ],
        ),
        (
            "tie-three",
            vec![
                vec![bin("chr1", 100, 200, 1)],
                vec![bin("chr1", 100, 200, 2)],
                vec![bin("chr1", 100, 200, 3)],
            ],
        ),
        (
            "same-start",
            vec![
                vec![bin("chr1", 100, 500, 1)],
                vec![bin("chr1", 100, 200, 2)],
            ],
        ),
        (
            "two-contigs",
            vec![
                vec![bin("chr1", 100, 200, 1), bin("chr3", 100, 200, 3)],
                vec![bin("chr2", 100, 200, 2)],
            ],
        ),
    ] {
        let merged = merge(&inputs, &full()).expect("a walk that is not refused");
        assert_eq!(printed(&merged), features(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 5, "the golden's merged streams");
}

/// Three files each holding the same interval come out first, third, second. Not the order they
/// were named in, and not an order anything in the comparator decides.
#[test]
fn the_heap_is_not_stable() {
    let text = golden();
    let merged = merge(
        &[
            vec![bin("chr1", 100, 200, 1)],
            vec![bin("chr1", 100, 200, 2)],
            vec![bin("chr1", 100, 200, 3)],
        ],
        &full(),
    )
    .expect("a walk that is not refused");
    let counts: Vec<&str> = merged
        .iter()
        .map(|feature| feature.text.rsplit('\t').next().expect("a count"))
        .collect();
    assert_eq!(counts, vec!["1", "3", "2"]);
    assert_eq!(printed(&merged), features(&text, "tie-three"));
}

/// A contig the dictionary does not name has index -1 and sorts before every named contig, so a
/// sorted file is reported as unsorted and the message names neither the dictionary nor the
/// contig that is missing from it.
#[test]
fn an_unnamed_contig_is_reported_as_an_unsorted_input() {
    let text = golden();
    let error = merge(
        &[vec![bin("chr1", 100, 200, 1), bin("chr2", 100, 200, 2)]],
        &dictionary("sequence-dictionary", &["chr1", "chr3"]),
    )
    .expect_err("a refused walk");
    assert_eq!(
        error,
        WalkerError::NotSorted {
            contig: "chr2".to_string(),
            start: 101
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "partly-unknown-contig")
    );
}

/// The check compares the replacement against the entry just returned, from the same input, so it
/// fires when the next record is drawn and names the new feature's locus.
#[test]
fn the_sort_check_names_the_new_feature() {
    let text = golden();
    let error = merge(
        &[
            vec![bin("chr1", 300, 400, 3), bin("chr1", 100, 200, 1)],
            vec![bin("chr1", 200, 300, 2)],
        ],
        &full(),
    )
    .expect_err("a refused walk");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "unsorted")
    );
}

#[test]
fn the_dictionary_refusals_match_the_golden() {
    let text = golden();

    let error = choose_dictionary(None, None).expect_err("no dictionary");
    assert_eq!(error, WalkerError::NoDictionary);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-dictionary")
    );

    // A master dictionary holding a contig the reference does not.
    let error = choose_dictionary(
        Some(dictionary("sequence-dictionary", &["chr1", "chrX"])),
        Some(dictionary("reference", &["chr1", "chr2", "chr3"])),
    )
    .expect_err("an absent contig");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "dict-contig-absent")
    );

    // The same contigs in another order. Both are three long, and on equal sizes the reference is
    // offered as the new dictionary and is therefore the one called small.
    let error = choose_dictionary(
        Some(dictionary("sequence-dictionary", &["chr2", "chr1", "chr3"])),
        Some(dictionary("reference", &["chr1", "chr2", "chr3"])),
    )
    .expect_err("contigs out of order");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "dict-out-of-order")
    );
}

/// The larger dictionary wins when the smaller is a subset of it in the same order.
#[test]
fn the_more_comprehensive_dictionary_wins() {
    let chosen = choose_dictionary(
        Some(dictionary("sequence-dictionary", &["chr1", "chr3"])),
        Some(dictionary("reference", &["chr1", "chr2", "chr3"])),
    )
    .expect("a compatible pair");
    assert_eq!(chosen.source, "reference");
    assert_eq!(chosen.contigs, vec!["chr1", "chr2", "chr3"]);
}
