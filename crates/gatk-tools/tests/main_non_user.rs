//! Conformance for `handleNonUserException` against GATK 4.6.2.0, compared as the line it writes.
//!
//! Golden from `tools/readfilter-conformance/MainNonUserDump.java`. `main-entry` measured the
//! OTHER handler and left this one alone, so the port had a single banner and printed it for both:
//! a failure the reference reports as `java.lang.IllegalArgumentException: Dictionary cannot have
//! size zero` came back as `A USER ERROR has occurred: Dictionary cannot have size zero`, with the
//! right message, the right status and the wrong wrapper (#1020).
//!
//! # What this suite is for
//!
//!  * **the handler being `printStackTrace` and nothing else**: no banner, no prefix, no notice
//!    about a system property;
//!  * **the first line being `Throwable.toString()`**, which is the class's BINARY name, so a
//!    nested class carries a `$` where its source carries a dot;
//!  * **a null message being the class alone**, with no trailing colon, where an empty message is
//!    the colon and nothing after it;
//!  * **and the `Error` overload being the same handler**, so an `OutOfMemoryError` differs from
//!    an exception only in the status `mainEntry` exits with.
//!
//! The frames the reference prints under that line are its own stack, and this port has none to
//! print. That is a boundary rather than an omission, and it is stated rather than hidden: the
//! golden's `shape` row says everything after the first line was a frame, and this test asserts
//! that the port's report is the message and stops.
//!
//! While the suite is `golden-pending` the dump is named by `MAIN_NON_USER_DUMP`.

use gatk_tools::main_entry::non_user_exception_report;

/// The throwables the harness constructs, which are this test's inputs where the golden is its
/// answer. The rendering is what is compared; the class and the message are what produced it.
const CASES: &[(&str, &str, Option<&str>)] = &[
    (
        "illegal-argument",
        "java.lang.IllegalArgumentException",
        Some("Dictionary cannot have size zero"),
    ),
    (
        "sam-format",
        "htsjdk.samtools.SAMFormatException",
        Some("Error parsing text SAM file. Not enough fields; Line 1"),
    ),
    ("no-message", "java.lang.IllegalStateException", None),
    ("empty-message", "java.lang.IllegalStateException", Some("")),
    (
        "nested-class",
        "MainNonUserDump$Nested",
        Some("thrown from a class inside another"),
    ),
    (
        "gatk-nested",
        "org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException",
        Some("a branch that was supposed to be unreachable"),
    ),
    (
        "multi-line-message",
        "java.lang.IllegalArgumentException",
        Some("first line\nsecond line"),
    ),
    (
        "out-of-memory",
        "java.lang.OutOfMemoryError",
        Some("Java heap space"),
    ),
];

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(dump: &str, kind: &str, name: &str) -> String {
    let prefix = format!("{kind}\t{name}\t");
    dump.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
        .unwrap_or_else(|| panic!("{kind}/{name} is not in the dump"))
}

#[test]
fn the_non_user_handler_prints_the_class_and_the_message() {
    let dump = match std::env::var("MAIN_NON_USER_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by MAIN_NON_USER_DUMP"),
        Err(_) => {
            println!(
                "skipped: the main-non-user golden is still pending. Run the suite and point \
                 MAIN_NON_USER_DUMP at tools/conformance/pending/main-non-user.MainNonUserDump.txt"
            );
            return;
        }
    };

    for (case, exception, message) in CASES {
        let report = non_user_exception_report(exception, *message);
        let first = report.lines().next().unwrap_or_default();
        assert_eq!(first, field(&dump, "non-user", case), "{case}");

        // The message is the whole of what follows, and the port adds nothing to it: no banner
        // above, no property notice below. The reference's frames are the difference, and the
        // golden's own row is what says they were frames and nothing else.
        let expected_lines = message.map(|text| text.lines().count().max(1)).unwrap_or(1);
        assert_eq!(
            report.lines().count(),
            expected_lines,
            "{case}: extra lines"
        );
        assert!(report.ends_with('\n'), "{case}: the handler ends its line");
        assert!(
            field(&dump, "shape", case).contains("stdout=0"),
            "{case}: the handler writes to stderr alone"
        );
    }
}

/// The two handlers write different things, which is the whole of the finding behind #1020.
#[test]
fn the_two_handlers_do_not_agree() {
    let user = gatk_tools::main_entry::user_exception_report("Dictionary cannot have size zero");
    let non_user = non_user_exception_report(
        "java.lang.IllegalArgumentException",
        Some("Dictionary cannot have size zero"),
    );
    assert!(user.contains("A USER ERROR has occurred: "));
    assert!(!non_user.contains("A USER ERROR has occurred: "));
    assert!(!non_user.contains(gatk_tools::main_entry::BANNER));
    assert!(!non_user.contains("system property"));
    assert_eq!(
        non_user,
        "java.lang.IllegalArgumentException: Dictionary cannot have size zero\n"
    );
}
