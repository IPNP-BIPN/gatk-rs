//! Conformance for `AddCommentsToBam` against Picard 3.4.0, compared as the whole output file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/AddCommentsToBamDump.java`, whose fixtures and outputs
//! travel as base64.
//!
//! # What this suite is for
//!
//!  * **the rewritten bytes**, which are the reference's exactly, for every accepted run;
//!  * **the comments landing after the file's own**, in the order given;
//!  * **the `@CO` prefix being added once**, not twice, for a comment that already carries it;
//!  * **a tab surviving inside a comment**, which lets one forge extra header fields;
//!  * **the sam refusal reading the name**, so a BAM named `.sam` is refused and a sam named
//!    `.bam` gets past it and fails in the copy;
//!  * **and no comment at all still being a rewrite**.
//!
//! The `md5` and `indexed` runs write files beside the output, which the reference records and the
//! port does not produce; their outputs are compared like the rest, which is the point: neither
//! changes a byte of the BAM.

use gatk_corpus as corpus;
use gatk_tools::add_comments_to_bam::{add_comments, is_refused_by_name, CommentError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/add_comments_to_bam.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .to_string()
}

fn bytes(text: &str, kind: &str, label: &str) -> Vec<u8> {
    corpus::decode_base64(&row(text, kind, label))
}

fn sam(text: &str, label: &str) -> String {
    let prefix = format!("sam\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries sam/{label}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    unescape(&row(text, "error", label))
}

const PATH: &str = "<dir>/reads.bam";

#[test]
fn every_rewritten_file_matches_the_golden() {
    let text = golden();
    let input = bytes(&text, "fixture", "reads");
    let mut compared = 0;
    for (label, comments) in [
        ("one-comment", vec!["a comment"]),
        ("two-comments", vec!["first", "second"]),
        ("no-comment", vec![]),
        ("already-prefixed", vec!["@CO\talready prefixed"]),
        ("with-tab", vec!["key\tvalue"]),
        // Neither of these changes the output's own bytes.
        ("md5", vec!["a comment"]),
        ("indexed", vec!["a comment"]),
    ] {
        let comments: Vec<String> = comments.iter().map(|c| c.to_string()).collect();
        let ours = add_comments(&input, PATH, &comments).expect("a run the tool allows");
        assert_eq!(ours, bytes(&text, "output", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 7, "the golden's outputs");
}

/// The prefix is added once, which is what the measurement corrected and what the upstream fix
/// made true of the header type this port uses.
#[test]
fn the_comment_prefix_is_added_once() {
    let text = golden();
    let comments: Vec<String> = vec!["@CO\talready prefixed".to_string()];
    let ours = add_comments(&bytes(&text, "fixture", "reads"), PATH, &comments)
        .expect("a run the tool allows");
    assert_eq!(ours, bytes(&text, "output", "already-prefixed"));
    let written = sam(&text, "already-prefixed");
    let comment_lines: Vec<&str> = written
        .lines()
        .filter(|line| line.starts_with("@CO"))
        .collect();
    assert_eq!(
        comment_lines,
        vec!["@CO\tan existing comment", "@CO\talready prefixed"]
    );
}

/// A tab inside a comment reaches the file whole, so the comment shows up as two fields.
#[test]
fn a_tab_survives_inside_a_comment() {
    let text = golden();
    let written = sam(&text, "with-tab");
    assert!(written.lines().any(|line| line == "@CO\tkey\tvalue"));
}

#[test]
fn the_sam_refusal_reads_the_name() {
    let text = golden();
    assert!(is_refused_by_name("<dir>/misnamed.sam"));
    assert!(!is_refused_by_name("<dir>/really-sam.bam"));

    // A BAM named `.sam`: refused by the tool.
    let error = add_comments(
        &bytes(&text, "fixture", "reads"),
        "<dir>/misnamed.sam",
        &["a comment".to_string()],
    )
    .expect_err("the suffix refusal");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "named-sam")
    );

    // A sam named `.bam`: past that check, and refused by the copy instead.
    let error = add_comments(
        &bytes(&text, "fixture", "really-sam"),
        "<dir>/really-sam.bam",
        &["a comment".to_string()],
    )
    .expect_err("the copy's refusal");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "sam-named-bam")
    );
}

/// The tool's own newline check, which the command line can never reach: Picard's parser refuses
/// the argument first, and the golden records that instead.
#[test]
fn a_newline_is_refused_by_the_tool_and_by_the_parser_before_it() {
    let text = golden();
    let error = add_comments(
        &bytes(&text, "fixture", "reads"),
        PATH,
        &["first line\nsecond line".to_string()],
    )
    .expect_err("the tool's own check");
    assert_eq!(error, CommentError::ContainsNewline);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        "picard.PicardException:Comments can not contain a new line"
    );
    // What the reference actually answered from the command line. The character in the message is
    // a REAL newline, not the two characters of an escape, which is why this literal carries one.
    assert_eq!(
        refusal(&text, "with-newline"),
        "java.lang.IllegalArgumentException:Supplied String contains illegal character '\n'."
    );
}
