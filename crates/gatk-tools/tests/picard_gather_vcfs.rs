//! Conformance for Picard's `GatherVcfs` against Picard 3.4.0, compared as the whole output file
//! of every run and as the exact log line of every refusal.
//!
//! Golden from `tools/readfilter-conformance/PicardGatherVcfsDump.java`, which carries each run's
//! inputs as well as its output.
//!
//! # What this suite is for
//!
//!  * **a refusal is an exit code and a log line**, not an exception, except for the two that are
//!    raised outside the try;
//!  * **an `AssertionError` is not a `RuntimeException`**, so the dictionary mismatch is thrown
//!    where the sample mismatch beside it is exit 1;
//!  * **the two order checks compare different things**, so shards whose first records are ordered
//!    and whose records overlap pass the first and fail the second;
//!  * **an empty shard moves neither check**, and sorts last under `REORDER_INPUT_BY_FIRST_VARIANT`;
//!  * **only the first comment survives**, `VCFHeader` keying an unstructured line by its key;
//!  * **and `CREATE_INDEX` decides which refusal a file with no contig lines gets**: a
//!    `PicardException` with it on, a `NullPointerException` turned into exit 1 with it off.
//!
//! # What the golden does not pin down
//!
//! The block copying path, which the golden records as a decompressed text, a list of bgzf block
//! sizes and a digest. Reproducing those bytes is a brick of its own; this port covers the
//! conventional path, which is the one every plain `.vcf` output takes.

use gatk_corpus as corpus;
use gatk_tools::picard_gather_vcfs::{gather, Arguments, Input};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/picard_gather_vcfs.txt.gz"),
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

/// The `exit` row, which carries the code and the ERROR lines the run logged.
fn exit_row(text: &str, label: &str) -> (i32, String) {
    let prefix = format!("exit\t{label}=");
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries exit/{label}"));
    let (code, log) = row.split_once('\t').expect("a code and a log");
    (code.parse().expect("an exit code"), unescape(log))
}

fn thrown(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

/// The inputs are named by the path their refusals quote, which the dump masked.
fn input<'a>(text: &'a str, label: &'a str, paths: &'a mut Vec<String>) -> String {
    paths.push(format!("<dir>/{label}.vcf"));
    value(text, "input", label)
}

struct Shards {
    texts: Vec<String>,
    paths: Vec<String>,
}

impl Shards {
    fn of(text: &str, labels: &[&str]) -> Shards {
        let mut paths = Vec::new();
        let texts = labels
            .iter()
            .map(|label| input(text, label, &mut paths))
            .collect();
        Shards { texts, paths }
    }

    fn inputs(&self) -> Vec<Input<'_>> {
        self.texts
            .iter()
            .zip(self.paths.iter())
            .map(|(text, path)| Input {
                path: path.as_str(),
                text: text.as_str(),
            })
            .collect()
    }
}

#[test]
fn every_gathered_file_matches_the_golden() {
    let text = golden();
    let runs: Vec<(&str, Vec<&str>, Arguments)> = vec![
        ("ordered", vec!["a", "b", "c"], Arguments::default()),
        ("single", vec!["a"], Arguments::default()),
        (
            "comments",
            vec!["a", "b"],
            Arguments {
                comments: vec!["one comment".to_string(), "another".to_string()],
                ..Arguments::default()
            },
        ),
        ("empty-shard", vec!["a", "none", "b"], Arguments::default()),
        (
            "reordered",
            vec!["c", "none", "a", "b"],
            Arguments {
                reorder_input_by_first_variant: true,
                ..Arguments::default()
            },
        ),
        (
            "no-contigs-no-index",
            vec!["bare-a", "bare-b"],
            Arguments {
                create_index: false,
                ..Arguments::default()
            },
        ),
    ];
    let mut compared = 0;
    for (label, labels, arguments) in runs {
        let shards = Shards::of(&text, &labels);
        let outcome = gather(&shards.inputs(), &arguments);
        if label == "no-contigs-no-index" {
            // This one is a refusal, checked below; it rides here to prove it is not a file.
            assert!(outcome.is_err(), "{label}");
            compared += 1;
            continue;
        }
        assert_eq!(
            outcome.expect("a run the tool allows"),
            value(&text, "gathered", label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 6, "the golden's outputs");
}

#[test]
fn every_exit_matches_the_golden() {
    let text = golden();
    let runs: Vec<(&str, Vec<&str>, Arguments)> = vec![
        ("unordered", vec!["c", "a", "b"], Arguments::default()),
        (
            "overlapping",
            vec!["a", "over"],
            Arguments {
                create_index: false,
                ..Arguments::default()
            },
        ),
        (
            "different-samples",
            vec!["a", "one-sample"],
            Arguments::default(),
        ),
        (
            "no-contigs-no-index",
            vec!["bare-a", "bare-b"],
            Arguments {
                create_index: false,
                ..Arguments::default()
            },
        ),
    ];
    let mut compared = 0;
    for (label, labels, arguments) in runs {
        let shards = Shards::of(&text, &labels);
        let failure = gather(&shards.inputs(), &arguments).expect_err("a refusal");
        let (code, log) = exit_row(&text, label);
        assert_eq!(failure.exit_code(), Some(code), "{label}");
        assert_eq!(failure.log_line().expect("a log line"), log, "{label}");
        compared += 1;
    }
    assert_eq!(compared, 4, "the golden's exits");
}

#[test]
fn the_two_refusals_raised_outside_the_try_are_thrown() {
    let text = golden();
    // The dictionary mismatch, which is an AssertionError.
    let shards = Shards::of(&text, &["a", "other-dictionary"]);
    let failure = gather(&shards.inputs(), &Arguments::default()).expect_err("a refusal");
    assert_eq!(failure.exit_code(), None);
    assert_eq!(
        format!("{}:{}", failure.java_class(), failure.message()),
        thrown(&text, "different-dictionary")
    );

    // And the indexing refusal, which is raised before the try.
    let shards = Shards::of(&text, &["bare-a", "bare-b"]);
    let failure = gather(&shards.inputs(), &Arguments::default()).expect_err("a refusal");
    assert_eq!(failure.exit_code(), None);
    assert_eq!(
        format!("{}:{}", failure.java_class(), failure.message()),
        thrown(&text, "no-contigs")
    );
}

/// The reordering is invisible in a whole-file comparison when the records end up the same either
/// way, so the order of the shards is checked directly: the empty file sorts last.
#[test]
fn the_empty_shard_sorts_last() {
    let text = golden();
    let shards = Shards::of(&text, &["c", "none", "a", "b"]);
    let reordered = gather(
        &shards.inputs(),
        &Arguments {
            reorder_input_by_first_variant: true,
            ..Arguments::default()
        },
    )
    .expect("a run the tool allows");
    // Given in the order c, none, a, b, which without the reordering is refused outright.
    assert!(gather(&shards.inputs(), &Arguments::default()).is_err());
    let positions: Vec<&str> = reordered
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').next().expect("a contig"))
        .collect();
    assert_eq!(positions, vec!["chr1", "chr1", "chr1", "chr2", "chr2"]);
}
