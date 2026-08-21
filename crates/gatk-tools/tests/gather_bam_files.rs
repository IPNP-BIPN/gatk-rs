//! Conformance for `GatherBamFiles` against Picard 3.4.0, compared as the whole output file of
//! every run that takes the block copying path.
//!
//! Golden from `tools/readfilter-conformance/GatherBamFilesDump.java`, whose fixtures and outputs
//! all travel as base64.
//!
//! # What this suite is for
//!
//!  * **the gathered bytes**, which are the reference's exactly, terminator included;
//!  * **the header being the first file's**, so a record referencing a read group the header never
//!    declares comes through untouched;
//!  * **an empty shard and a list file changing nothing**, both producing the same bytes as
//!    naming the files directly;
//!  * **the order never being checked**, so shards gathered backwards give a file whose header
//!    lies about its own sort order;
//!  * **the path choice**, which reads the files rather than their names;
//!  * **and the `.md5`**, which is the digest of those same bytes as hex.
//!
//! # What this port does not do
//!
//! The record-by-record gather, which runs when any input is not a BAM. Its bytes are in the
//! golden under `with-sam`; reproducing them needs a BAM writer driven from a sam reader, which is
//! a brick of its own. The suite asserts that the choice is made correctly and that the two paths
//! disagree, not the second path's bytes.

use gatk_corpus as corpus;
use gatk_tools::gather_bam_files::{gather, is_bam_file, md5_file, unroll, use_block_copying};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/gather_bam_files.txt.gz"),
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

/// The `sam\t<label>=` rows, which are the output read back as text.
fn sam(text: &str, label: &str) -> String {
    let prefix = format!("sam\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries sam/{label}")),
    )
}

fn shards(text: &str, labels: &[&str]) -> Vec<Vec<u8>> {
    labels
        .iter()
        .map(|label| bytes(text, "fixture", label))
        .collect()
}

fn gathered(text: &str, labels: &[&str]) -> Vec<u8> {
    let held = shards(text, labels);
    let inputs: Vec<&[u8]> = held.iter().map(|shard| shard.as_slice()).collect();
    gather(&inputs).expect("a run the tool allows")
}

#[test]
fn every_block_copied_output_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, inputs) in [
        ("two-shards", vec!["first", "second"]),
        ("single", vec!["first"]),
        ("with-empty", vec!["first", "empty", "second"]),
        ("other-read-group", vec!["first", "other-rg"]),
        ("out-of-order", vec!["second", "earlier"]),
        // The index and the digest are written beside the output and change none of its bytes.
        ("indexed", vec!["first", "second"]),
        ("md5", vec!["first", "second"]),
        // A list file names the same two shards.
        ("unrolled", vec!["first", "second"]),
    ] {
        assert_eq!(
            gathered(&text, &inputs),
            bytes(&text, "output", label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 8, "the golden's block copied outputs");
}

/// The header is the first file's, so the second shard's read group is gone from the header while
/// its records still name it.
#[test]
fn a_record_can_reference_a_read_group_the_header_never_declares() {
    let text = golden();
    let output = sam(&text, "other-read-group");
    let groups: Vec<&str> = output
        .lines()
        .filter(|line| line.starts_with("@RG"))
        .collect();
    assert_eq!(groups.len(), 1);
    assert!(groups[0].contains("ID:rg1"));
    assert!(output.contains("RG:Z:rg2"));
}

/// Nothing checks the order, so the header keeps saying `coordinate` over records that are not.
#[test]
fn the_order_is_never_checked() {
    let text = golden();
    let output = sam(&text, "out-of-order");
    assert!(output
        .lines()
        .next()
        .expect("a header")
        .contains("SO:coordinate"));
    let starts: Vec<i32> = output
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| {
            line.split('\t')
                .nth(3)
                .expect("a position")
                .parse()
                .expect("a number")
        })
        .collect();
    assert_eq!(starts, vec![300, 400, 10, 20]);
}

/// The choice reads the files: the sam fixture is not a BAM however it is named, and the run that
/// includes it produces different bytes.
#[test]
fn one_sam_sends_the_whole_run_down_the_other_path() {
    let text = golden();
    assert!(is_bam_file(&bytes(&text, "fixture", "first")));
    assert!(!is_bam_file(&bytes(&text, "fixture", "second-sam")));
    let held = shards(&text, &["first", "second-sam"]);
    let inputs: Vec<&[u8]> = held.iter().map(|shard| shard.as_slice()).collect();
    assert!(!use_block_copying(&inputs));
    // And the reference's own output for that run is not the block copy's.
    assert_ne!(
        bytes(&text, "output", "with-sam"),
        bytes(&text, "output", "two-shards")
    );
}

#[test]
fn the_md5_is_the_digest_of_those_same_bytes() {
    let text = golden();
    assert_eq!(
        md5_file(&gathered(&text, &["first", "second"])),
        row(&text, "md5", "md5")
    );
}

/// `IOUtil.unrollFiles` replaces a list by the paths it names and leaves a BAM alone.
#[test]
fn a_list_of_paths_is_unrolled() {
    let entries = vec![
        "/tmp/first.bam".to_string(),
        "/tmp/a.bam\n/tmp/b.bam\n".to_string(),
    ];
    assert_eq!(
        unroll(&entries),
        vec!["/tmp/first.bam", "/tmp/a.bam", "/tmp/b.bam"]
    );
}
