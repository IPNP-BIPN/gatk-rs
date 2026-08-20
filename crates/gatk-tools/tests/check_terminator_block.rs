//! Conformance for `CheckTerminatorBlock` against Picard 3.4.0, compared as the termination answer
//! and the exit code of every file.
//!
//! Golden from `tools/readfilter-conformance/CheckTerminatorBlockDump.java`, whose nine files it
//! carries as base64, so this replays the reference's own bytes rather than bytes a Rust writer
//! produced.
//!
//! # What this suite is for
//!
//!  * **a file with no terminator passes**, its last block being healthy;
//!  * **a file one byte shorter than that is defective**, because the last block is then short;
//!  * **a file shorter than the terminator is defective** without being searched;
//!  * **a file that is not gzip at all is defective** rather than an exception;
//!  * **a corrupt payload is not noticed**, the check never decompressing anything;
//!  * **and a terminator with one byte wrong falls through to the backwards search**, which finds
//!    the block healthy.

use gatk_corpus as corpus;
use gatk_tools::check_terminator_block::{check, termination_name};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/check_terminator_block.txt.gz"),
    )
}

fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .map(|rest| rest.trim_start_matches(['=', '\t']).to_string())
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
}

#[test]
fn every_file_answers_what_the_reference_answers() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "complete",
        "no-terminator",
        "truncated-block",
        "terminator-only",
        "too-short",
        "empty",
        "not-gzip",
        "corrupt-payload",
        "wrong-terminator",
    ] {
        let bytes = corpus::decode_base64(&row(&text, "fixture", label));
        let (termination, code) = check(&bytes);
        assert_eq!(
            termination_name(termination),
            row(&text, "termination", label),
            "{label}: the termination"
        );
        assert_eq!(
            code.to_string(),
            row(&text, "exit", label),
            "{label}: the exit code"
        );
        compared += 1;
    }
    assert_eq!(compared, 9, "the golden's files");
}
