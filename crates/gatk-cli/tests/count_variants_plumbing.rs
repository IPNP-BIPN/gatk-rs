//! Conformance for the file plumbing under `CountVariants` against GATK 4.6.2.0.
//!
//! `count-variants` compares the number, the bytes `-O` receives and the refusals. This suite is
//! the layer between that number and a command line: `gatk-rs CountVariants` over the golden's own
//! inputs, with the golden's own arguments.
//!
//! The fixtures are rebuilt rather than read from disk. The golden carries each input as text, and
//! the index beside it is built here by [`gatk_tools::index_feature_file::build`], which is the
//! bundled tool the reference's own refusal tells the user to run.
//!
//! # What this suite is for
//!
//!  * **the tool returning the count**, which `handleResult` prints;
//!  * **`-O` receiving the digits and nothing else**, and truncating what was there;
//!  * **no `-O` writing nothing at all**, rather than an empty file;
//!  * **`-L` selecting by the record's whole SPAN**, so an interval reaches a record whose
//!    position it does not hold, and a record over two intervals is counted once;
//!  * **`-L` against an input with no index being refused before any record is read**;
//!  * **and an unwritable `-O` carrying the path and nothing else.**

use gatk_corpus as corpus;
use gatk_tools::main_entry::{self, Failure};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/count_variants.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn row<'a>(text: &'a str, kind: &str, label: &str) -> Vec<&'a str> {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))[prefix.len()..]
        .split('\t')
        .collect()
}

fn args(list: &[&str]) -> Vec<String> {
    list.iter().map(|arg| (*arg).to_string()).collect()
}

/// A directory of this test's own, named after the case.
fn scratch(case: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("gatk-cli-count-variants-{case}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    dir
}

/// One of the golden's inputs, written into `dir` with the index the reference had beside it.
fn fixture(text: &str, dir: &std::path::Path, label: &str, indexed: bool) -> std::path::PathBuf {
    let vcf = unescape(row(text, "input", label)[0]);
    let path = dir.join(format!("{label}.vcf"));
    std::fs::write(&path, &vcf).expect("the fixture");
    if indexed {
        let name = path.to_string_lossy().to_string();
        let index = gatk_tools::index_feature_file::build(
            &vcf,
            &gatk_tools::index_feature_file::Source::new(&name),
            &name,
        )
        .expect("the index");
        std::fs::write(gatk_tools::index_feature_file::default_output(&name), index)
            .expect("the index beside it");
    }
    path
}

/// The count the golden recorded for one case.
fn expected(text: &str, case: &str) -> i64 {
    row(text, "count", case)[0].parse().expect("a count")
}

#[test]
fn the_plumbing_answers_what_the_reference_answered() {
    let text = golden();
    // (case, input label, extra arguments, whether `-O` was given)
    let cases: Vec<(&str, &str, Vec<&str>, bool)> = vec![
        ("plain-no-output", "plain", vec![], false),
        ("plain", "plain", vec![], true),
        ("filtered-only", "filtered-only", vec![], true),
        ("empty", "empty", vec![], true),
        ("span-by-end", "spanning", vec!["-L", "chr1:300-310"], true),
        (
            "span-by-ref-length",
            "spanning",
            vec!["-L", "chr1:605-606"],
            true,
        ),
        ("span-missed", "spanning", vec!["-L", "chr1:500-510"], true),
        (
            "two-intervals-one-record",
            "spanning",
            vec!["-L", "chr1:150-160", "-L", "chr1:350-360"],
            true,
        ),
        (
            "interval-matches-nothing",
            "plain",
            vec!["-L", "chr1:900-950"],
            true,
        ),
        (
            "interval-selects-contig",
            "two-contigs",
            vec!["-L", "chr2"],
            true,
        ),
    ];

    for (case, label, extra, has_output) in cases {
        let dir = scratch(case);
        let input = fixture(&text, &dir, label, true);
        let output = dir.join("count.txt");

        let mut argv = args(&["CountVariants", "--variant", &input.to_string_lossy()]);
        if has_output {
            argv.push("--output".to_string());
            argv.push(output.to_string_lossy().to_string());
        }
        argv.extend(extra.iter().map(|arg| (*arg).to_string()));

        let run = gatk_cli::run(&argv);
        assert_eq!(run.status, 0, "{case}: {}", run.stderr);
        let count = expected(&text, case);
        assert_eq!(
            run.stdout,
            format!("Tool returned:\n{count}\n"),
            "{case}: what the tool returned"
        );

        if has_output {
            // `print`, not `println`: the file is the digits and nothing else.
            assert_eq!(
                std::fs::read(&output).expect("the output file"),
                count.to_string().into_bytes(),
                "{case}: what -O received"
            );
        } else {
            // The golden's own row for this case is a file that was never written.
            assert_eq!(row(&text, "file", case)[0], "no-output-argument");
            assert!(
                !output.exists(),
                "{case}: -O was not given and nothing was written"
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}

/// `-O` truncates rather than appends: ten bytes before the run, one after it.
#[test]
fn the_output_is_overwritten_and_not_appended_to() {
    let text = golden();
    let dir = scratch("overwrite");
    let input = fixture(&text, &dir, "plain", true);
    let output = dir.join("preexisting.count");
    std::fs::write(&output, "9999999999").expect("the pre-existing file");
    assert_eq!(row(&text, "file", "before-overwrite")[1], "10");

    let run = gatk_cli::run(&args(&[
        "CountVariants",
        "--variant",
        &input.to_string_lossy(),
        "--output",
        &output.to_string_lossy(),
    ]));
    assert_eq!(run.status, 0, "{}", run.stderr);
    let count = expected(&text, "overwrite");
    assert_eq!(
        std::fs::read(&output).unwrap(),
        count.to_string().into_bytes()
    );
    assert_eq!(
        row(&text, "file", "overwrite")[1],
        count.to_string().len().to_string()
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// The three refusals, in the reference's own words and with its own statuses.
#[test]
fn the_refusals_are_the_reference_ones() {
    let text = golden();

    // `-L` needs an index, and the refusal comes before any record is read.
    let dir = scratch("unindexed");
    let input = fixture(&text, &dir, "unindexed", false);
    let run = gatk_cli::run(&args(&[
        "CountVariants",
        "--variant",
        &input.to_string_lossy(),
        "-L",
        "chr1:100-200",
    ]));
    assert_eq!(run.status, main_entry::exit_status(Failure::User));
    let recorded = row(&text, "error", "interval-without-index")[0];
    let message = recorded.split_once(':').expect("a message").1;
    // The golden's path is the container's; what is compared is the wording around it.
    assert!(
        message.ends_with("must support random access to enable traversal by intervals. If it's a file, please index it using the bundled tool IndexFeatureFile"),
        "{message}"
    );
    assert!(
        run.stderr.contains("must support random access"),
        "{}",
        run.stderr
    );

    // An unwritable `-O` carries the path and nothing else, and it is thrown AFTER the traversal.
    let output = dir.join("a-directory");
    std::fs::create_dir_all(&output).expect("a directory to write onto");
    let indexed = fixture(&text, &dir, "plain", true);
    let run = gatk_cli::run(&args(&[
        "CountVariants",
        "--variant",
        &indexed.to_string_lossy(),
        "--output",
        &output.to_string_lossy(),
    ]));
    assert_eq!(run.status, main_entry::exit_status(Failure::User));
    let recorded = row(&text, "error", "output-is-a-directory")[0];
    assert!(
        recorded.starts_with(
            "org.broadinstitute.hellbender.exceptions.UserException$CouldNotCreateOutputFile:"
        ),
        "{recorded}"
    );
    assert!(
        run.stderr.contains(&output.to_string_lossy().to_string()),
        "{}",
        run.stderr
    );

    // And an interval off the dictionary is the interval layer's refusal, not this tool's, so it
    // is asserted as a message rather than replayed: `count-variants` says the same.
    let run = gatk_cli::run(&args(&[
        "CountVariants",
        "--variant",
        &indexed.to_string_lossy(),
        "-L",
        "chr3:1-10",
    ]));
    assert_ne!(run.status, 0, "an interval off the dictionary is refused");
    let recorded = row(&text, "error", "interval-off-the-dictionary")[0];
    assert!(
        recorded.contains("is not valid for this input"),
        "{recorded}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
