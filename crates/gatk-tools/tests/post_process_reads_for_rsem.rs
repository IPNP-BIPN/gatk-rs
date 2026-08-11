//! Conformance for `PostProcessReadsForRSEM` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/PostProcessReadsForRSEMDump.java`. One run finishes and
//! four crash or refuse, so most of this suite is messages rather than bytes: the tool's two null
//! dereferences are as much a part of what it does as the file it writes.
//!
//! # What this suite is for
//!
//!  * **a fourth `getDefaultReadFilters` pattern**: the whole list replaced by a single filter that
//!    is not `Wellformed`, so a supplementary alignment never reaches the tool;
//!  * **two of its own null guards dereference null**, and the two shapes of the second one print
//!    different JVM messages because the comparison that follows evaluates its left operand first;
//!  * **the output order is primary pair, then each secondary pair**, first-of-pair before
//!    second-of-pair, which is not the order the query-name group arrived in;
//!  * **three reasons drop both reads**, and a failing primary takes its secondary alignments with
//!    it while a failing secondary drops only itself.

use gatk_corpus as corpus;
use gatk_engine::reads::ReadsDataSource;
use gatk_tools::post_process_reads_for_rsem::{
    self as tool, RsemError, READ1_IS_NULL, READ1_READS_IS_NULL, READ2_READS_IS_NULL,
};
use gatk_tools::sam_output::Options;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/post_process_reads_for_rsem.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn of_run<'a>(text: &'a str, kind: &str, label: &str) -> Vec<Vec<&'a str>> {
    rows(text, kind)
        .into_iter()
        .filter(|row| row[0] == label)
        .collect()
}

/// Which fixture each labelled run was given.
fn fixture_for(label: &str) -> &'static str {
    match label {
        "plain" => "plain",
        "nofirst" => "no_first",
        "onesided" => "one_sided",
        "othersided" => "other_sided",
        "coordinate" => "coordinate",
        other => panic!("{other} is in the golden but not configured here"),
    }
}

fn install(text: &str, dir: &std::path::Path) {
    std::fs::create_dir_all(dir).expect("a scratch directory");
    for row in rows(text, "fixture") {
        std::fs::write(
            dir.join(format!("{}.bam", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture bam");
    }
}

fn run(
    text: &str,
    dir: &std::path::Path,
    label: &str,
) -> Result<(Vec<u8>, Option<Vec<u8>>), RsemError> {
    let source = ReadsDataSource::open_unindexed(&dir.join(format!("{}.bam", fixture_for(label))))
        .expect("the fixture opens with no index");
    let command_line = of_run(text, "commandline", label)
        .first()
        .map(|row| row.get(1).copied().unwrap_or(""))
        .unwrap_or("");
    let options = Options {
        command_line,
        ..Options::default()
    };
    tool::post_process_reads_for_rsem(&source, &options).expect("the source reads")
}

#[test]
fn the_output_file_is_byte_identical() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-rsem-{}", std::process::id()));
    install(&text, &dir);

    let outputs = rows(&text, "output");
    assert_eq!(outputs.len(), 1, "one run finishes and four do not");

    let (label, expected_base64) = (outputs[0][0], outputs[0][1]);
    let (ours, our_index) = run(&text, &dir, label).expect("this label does not refuse");

    let expected = corpus::decode_base64(expected_base64);
    assert_eq!(ours.len(), expected.len(), "{label}: output length differs");
    if ours != expected {
        let at = ours
            .iter()
            .zip(&expected)
            .position(|(a, b)| a != b)
            .unwrap_or(0);
        panic!("{label}: first byte difference at offset {at}");
    }
    assert!(our_index.is_none(), "a queryname output has no index");
    assert_eq!(
        of_run(&text, "index", label).first().map(|row| row[1]),
        Some("absent")
    );

    println!("post-process-reads-for-rsem: the output is byte-identical");
}

/// The finding: two guards written to handle null dereference it.
#[test]
fn the_null_guards_dereference_null_and_name_two_different_lists() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-rsem-npe-{}", std::process::id()));
    install(&text, &dir);

    let message = |label: &str| -> String {
        of_run(&text, "error", label)
            .first()
            .map(|row| row[1].to_string())
            .unwrap_or_else(|| panic!("the golden lost the {label} crash"))
    };

    for (label, expected) in [
        ("nofirst", READ1_IS_NULL),
        ("onesided", READ2_READS_IS_NULL),
        ("othersided", READ1_READS_IS_NULL),
    ] {
        assert_eq!(
            message(label),
            format!("java.lang.NullPointerException:{expected}"),
            "{label}"
        );
        assert_eq!(
            run(&text, &dir, label).expect_err("the reference crashes here"),
            RsemError::NullDereference(expected),
            "{label}"
        );
    }

    // The two secondary shapes are mirror images and print different names, because
    // `read1Reads.size() != read2Reads.size()` evaluates its left operand first.
    assert_ne!(message("onesided"), message("othersided"));
}

/// The sort order is checked on the header, and the check is the tool's own.
#[test]
fn a_coordinate_sorted_input_is_refused() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-rsem-so-{}", std::process::id()));
    install(&text, &dir);

    let message = of_run(&text, "error", "coordinate")
        .first()
        .map(|row| row[1].to_string())
        .expect("the golden lost the coordinate refusal");
    assert_eq!(
        message,
        format!(
            "org.broadinstitute.hellbender.exceptions.UserException:{}",
            RsemError::NotQueryNameSorted.message()
        )
    );
    assert_eq!(
        run(&text, &dir, "coordinate").expect_err("refused"),
        RsemError::NotQueryNameSorted
    );
}

/// What survived, in what order, which is the whole point of the tool.
#[test]
fn the_output_is_reordered_and_most_of_the_input_is_gone() {
    let text = golden();
    let inputs: Vec<String> = rows(&text, "reads")
        .iter()
        .filter(|row| row[0] == "in:plain")
        .map(|row| format!("{}:{}", row[1], row[4]))
        .collect();
    let outputs: Vec<String> = of_run(&text, "reads", "plain")
        .iter()
        .map(|row| format!("{}:{}", row[1], row[4]))
        .collect();

    assert_eq!(inputs.len(), 20);
    assert_eq!(outputs.len(), 8, "twelve reads did not survive");

    // p1's primary pair, then its secondary pair, first-of-pair before second-of-pair each time.
    assert_eq!(
        outputs,
        vec![
            "p1:100", "p1:300", "p1:500", "p1:700", // primary then secondary
            "p6:2100", "p6:2300", // the supplementary alignment was filtered out
            "p8:3100", "p8:3300", // the chimeric secondary pair dropped, the primary survived
        ]
    );

    // p2 unmapped mate, p3 chimeric, p4 two cigar elements, p5 a single element that is not M,
    // p7 no second-of-pair: none of them reaches the output.
    for gone in ["p2", "p3", "p4", "p5", "p7"] {
        assert!(
            inputs.iter().any(|r| r.starts_with(gone)),
            "{gone} was in the input"
        );
        assert!(
            !outputs.iter().any(|r| r.starts_with(gone)),
            "{gone} should be gone"
        );
    }
}

/// One filter, and it is not `Wellformed`.
#[test]
fn the_filter_chain_is_a_single_filter() {
    let text = golden();
    let filters: Vec<&str> = rows(&text, "filters").iter().map(|row| row[0]).collect();
    assert_eq!(filters, vec!["NotSupplementaryAlignmentReadFilter"]);
    assert_eq!(
        tool::DEFAULT_READ_FILTERS.to_vec(),
        vec!["NotSupplementaryAlignmentReadFilter"]
    );
}
