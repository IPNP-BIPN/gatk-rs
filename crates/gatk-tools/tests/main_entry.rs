//! Conformance for `Main`'s entry against GATK 4.6.2.0, compared as the stream each path writes
//! to, the exit status each exception carries, and what the two handlers print.
//!
//! Golden from `tools/readfilter-conformance/MainEntryDump.java`. `mainEntry` ends in
//! `System.exit`, so the dump measures the pair rather than watching a process die: the five
//! constants read off the class, and the exception each path throws.
//!
//! # What this suite is for
//!
//!  * **no arguments printing the usage to stdout and returning nothing**;
//!  * **`-h` and `--help` being that path, and only as the first argument**;
//!  * **an unknown name printing the same usage to stderr and then refusing**;
//!  * **`--version` being scanned over every argument**;
//!  * **the version's first line going to stdout whatever stream it is handed**;
//!  * **the five statuses, and a name that does not resolve exiting differently from one that
//!    does**;
//!  * **a Picard tool's non-zero return arriving wrapped**;
//!  * **`handleResult` being silent for null**;
//!  * **and the decorated user error naming the property that prints a stack trace.**

use gatk_corpus as corpus;
use gatk_tools::main_entry::{
    decorated_exception_message, exit_status, route, tool_arguments, tool_returned,
    user_exception_report, version_lines, Failure, Route, Stream, ANY_OTHER_EXCEPTION_EXIT_VALUE,
    COMMANDLINE_EXCEPTION_EXIT_VALUE, OUT_OF_MEMORY_EXIT_VALUE, PICARD_TOOL_EXCEPTION,
    USER_EXCEPTION_EXIT_VALUE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/main_entry.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

fn words(text: &str) -> Vec<String> {
    text.split(' ').map(str::to_string).collect()
}

/// One case's `shape` row, as the two `<stream>=<lines>,<first line>` halves it holds.
fn shape(text: &str, case: &str) -> (String, String) {
    let row = field(text, "shape", case).unwrap_or_else(|| panic!("shape/{case}"));
    let parts = words(&row);
    (parts[0].clone(), parts[1..].join(" "))
}

fn line_count(half: &str) -> usize {
    let rest = half.split_once('=').expect("a count").1;
    rest.split_once(',')
        .expect("a first line")
        .0
        .parse()
        .expect("a number")
}

/// The five statuses the golden read off the class.
#[test]
fn the_five_statuses_are_the_ported_ones() {
    let text = golden();
    for (name, ours) in [
        (
            "COMMANDLINE_EXCEPTION_EXIT_VALUE",
            COMMANDLINE_EXCEPTION_EXIT_VALUE,
        ),
        ("USER_EXCEPTION_EXIT_VALUE", USER_EXCEPTION_EXIT_VALUE),
        ("PICARD_TOOL_EXCEPTION", PICARD_TOOL_EXCEPTION),
        (
            "ANY_OTHER_EXCEPTION_EXIT_VALUE",
            ANY_OTHER_EXCEPTION_EXIT_VALUE,
        ),
        ("OUT_OF_MEMORY_EXIT_VALUE", OUT_OF_MEMORY_EXIT_VALUE),
    ] {
        let theirs: i32 = field(&text, "status", name)
            .unwrap_or_else(|| panic!("status/{name}"))
            .parse()
            .expect("a number");
        assert_eq!(ours, theirs, "{name}");
    }
    // Which the catch blocks map to in the reference's own order.
    assert_eq!(exit_status(Failure::CommandLine), 1);
    assert_eq!(exit_status(Failure::User), 2);
    assert_eq!(exit_status(Failure::PicardNonZero), 4);
    assert_eq!(exit_status(Failure::Other), 3);
    assert_eq!(exit_status(Failure::OutOfMemory), 137);
}

/// No arguments is not an error: the usage goes to stdout and the run returns nothing.
#[test]
fn no_arguments_prints_the_usage_to_stdout() {
    let text = golden();
    for case in ["no-arguments", "dash-h", "long-help"] {
        let (out, err) = shape(&text, case);
        assert!(line_count(&out) > 300, "{case} printed the usage");
        assert_eq!(line_count(&err), 0, "{case} said nothing on stderr");
        assert_eq!(
            field(&text, "result", case).as_deref(),
            Some("null"),
            "{case}"
        );
    }
    assert_eq!(route(&[]), Route::Usage(Stream::Stdout));
    assert_eq!(route(&["-h".to_string()]), Route::Usage(Stream::Stdout));
    assert_eq!(route(&["--help".to_string()]), Route::Usage(Stream::Stdout));
}

/// The help test looks at the first argument alone, so a tool followed by `-h` resolves the tool.
#[test]
fn help_after_a_tool_resolves_the_tool() {
    let text = golden();
    let (out, err) = shape(&text, "help-after-a-tool");
    assert_eq!(line_count(&out), 0);
    // The usage is the tool's own, on stderr, and the run returns a code rather than null.
    assert!(err.contains("USAGE: CountReads [arguments]"), "{err}");
    assert_eq!(
        field(&text, "result", "help-after-a-tool").as_deref(),
        Some("0")
    );
    let args = vec!["CountReads".to_string(), "-h".to_string()];
    assert_eq!(
        route(&args),
        Route::Tool {
            name: "CountReads".to_string()
        }
    );
    // And the tool never sees its own name.
    assert_eq!(tool_arguments(&args), ["-h".to_string()]);
}

/// `--version` is scanned over every argument, not only the first.
#[test]
fn the_version_is_scanned_over_every_argument() {
    let text = golden();
    for case in ["version", "version-short", "version-after-a-tool"] {
        let (out, err) = shape(&text, case);
        assert_eq!(line_count(&out), 4, "{case}");
        assert_eq!(line_count(&err), 0, "{case}");
        assert_eq!(
            field(&text, "result", case).as_deref(),
            Some("null"),
            "{case}"
        );
    }
    assert_eq!(route(&["--version".to_string()]), Route::Version);
    assert_eq!(route(&["-version".to_string()]), Route::Version);
    assert_eq!(
        route(&["CountReads".to_string(), "--version".to_string()]),
        Route::Version
    );
    // A tool that resolves is still not run, which is what makes this a scan and not a fallback.
    assert_ne!(
        route(&["CountReads".to_string(), "--version".to_string()]),
        Route::Tool {
            name: "CountReads".to_string()
        }
    );
}

/// The version's first line goes to stdout whatever stream the printer is handed.
#[test]
fn the_version_is_torn_in_half_by_a_stream_that_is_not_stdout() {
    let text = golden();
    let first = field(&text, "out", "version-to-stdout").expect("the first line");
    let rest = field(&text, "out", "version-to-elsewhere").expect("the rest");
    assert_eq!(first, "The Genome Analysis Toolkit (GATK) v4.6.2.0\n");
    assert_eq!(rest, "HTSJDK Version: 4.2.0\nPicard Version: 3.4.0\n");
    let (ours_first, ours_rest) = version_lines(
        "The Genome Analysis Toolkit (GATK)",
        "4.6.2.0",
        "4.2.0",
        "3.4.0",
    );
    assert_eq!(ours_first, first);
    assert_eq!(ours_rest, rest);
    // The three versions are the reference's own pins, which is what the golden holds them for.
    assert!(rest.contains("4.2.0") && rest.contains("3.4.0"));
}

/// An unknown name prints the same usage to stderr and then refuses.
#[test]
fn an_unknown_name_prints_the_usage_to_stderr_and_refuses() {
    let text = golden();
    let (out, err) = shape(&text, "unknown-name");
    let (help_out, _) = shape(&text, "no-arguments");
    assert_eq!(line_count(&out), 0);
    // The very same usage, down to its line count, on the other stream.
    assert_eq!(line_count(&err), line_count(&help_out));
    let message = field(&text, "error", "unknown-name").expect("the refusal");
    assert!(message.starts_with("org.broadinstitute.hellbender.exceptions.UserException:"));
    assert!(message.contains("'NoSuchToolAtAll' is not a valid command."));
    match route(&["NoSuchToolAtAll".to_string()]) {
        Route::Unknown { message: ours } => {
            assert!(
                ours.starts_with("'NoSuchToolAtAll' is not a valid command.\n"),
                "{ours}"
            );
        }
        other => panic!("{other:?}"),
    }
}

/// A deprecated name is a user exception like any other, carrying the notice.
#[test]
fn a_deprecated_name_carries_its_notice() {
    let text = golden();
    let message = field(&text, "error", "deprecated-name").expect("the notice");
    assert!(message.starts_with("org.broadinstitute.hellbender.exceptions.UserException:"));
    assert!(message.contains(
        "IndelRealigner is no longer included in GATK as of version 4.0.0.0. \
         Please use GATK3 to run this tool"
    ));
    match route(&["IndelRealigner".to_string()]) {
        Route::Unknown { message: ours } => assert_eq!(
            ours,
            "IndelRealigner is no longer included in GATK as of version 4.0.0.0. \
             Please use GATK3 to run this tool"
        ),
        other => panic!("{other:?}"),
    }
}

/// A name that does not resolve exits 2; a name that does but is given nothing exits 1.
#[test]
fn the_two_refusals_exit_differently() {
    let text = golden();
    let unknown = field(&text, "error", "unknown-name").expect("the refusal");
    let missing = field(&text, "error", "gatk-tool-no-arguments").expect("the parse failure");
    assert!(unknown.contains("UserException"));
    assert!(missing.starts_with("org.broadinstitute.barclay.argparser.CommandLineException"));
    assert!(missing.contains("Argument 'input' is required"));
    assert_eq!(exit_status(Failure::User), 2);
    assert_eq!(exit_status(Failure::CommandLine), 1);
    assert_ne!(
        exit_status(Failure::User),
        exit_status(Failure::CommandLine)
    );
    // And the parse failure said nothing on either stream: the tool's usage is Barclay's to print.
    let (out, err) = shape(&text, "gatk-tool-no-arguments");
    assert_eq!((line_count(&out), line_count(&err)), (0, 0));
}

/// A Picard tool is wrapped, and its non-zero return arrives as an exception carrying the code.
#[test]
fn a_picard_tools_code_arrives_wrapped() {
    let text = golden();
    let message = field(&text, "error", "picard-tool-no-arguments").expect("the wrapper");
    assert!(
        message.starts_with("org.broadinstitute.hellbender.exceptions.PicardNonZeroExitException:")
    );
    // The code is the exception's own field and not its message, which is why the message is null.
    assert!(message.contains(":null code=1"), "{message}");
    assert_eq!(exit_status(Failure::PicardNonZero), 4);
    // Picard printed its own usage on the way out, so the wrapper is not silent.
    let (_, err) = shape(&text, "picard-tool-no-arguments");
    assert!(err.contains("USAGE: MarkDuplicates [arguments]"), "{err}");
}

/// `handleResult` is silent for null and prints the value under a line of its own otherwise.
#[test]
fn a_tool_that_returns_nothing_is_silent() {
    let text = golden();
    let printed = field(&text, "out", "handlers").expect("the handler output");
    // Three calls, one of them null, and only two lines' worth of output.
    assert_eq!(printed, "Tool returned:\n0\nTool returned:\na string\n");
    assert_eq!(tool_returned(None), None);
    assert_eq!(
        tool_returned(Some("0")),
        Some("Tool returned:\n0".to_string())
    );
    assert_eq!(
        tool_returned(Some("a string")),
        Some("Tool returned:\na string".to_string())
    );
}

/// The decorated user error, banner and property line both.
#[test]
fn the_user_error_names_the_property_that_prints_a_stack_trace() {
    let text = golden();
    let printed = field(&text, "err", "handlers").expect("the handler output");
    assert_eq!(
        printed,
        user_exception_report("the message a refusal carries")
    );
    // The banner is printed twice, with a blank line on each side of the message.
    let decorated = decorated_exception_message("A USER ERROR has occurred: ", "boom");
    assert_eq!(decorated.lines().count(), 5);
    assert_eq!(decorated.lines().next(), decorated.lines().nth(4));
    assert_eq!(decorated.lines().next().expect("a banner").len(), 71);
    assert!(printed.contains("GATK_STACKTRACE_ON_USER_EXCEPTION"));
}
