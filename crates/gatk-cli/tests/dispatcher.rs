//! Conformance for the dispatcher against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/MainEntryDump.java`, which ran `Main.instanceMain`
//! over eleven command lines and recorded, for each, which stream it wrote to, how many lines it
//! wrote, what the first of them was, what the run returned and what exception it ended in.
//!
//! # What this suite is for
//!
//!  * **the paths that write nothing being silent, and the ones that refuse writing the usage**;
//!  * **the version being four lines whatever else is on the command line**;
//!  * **a name that does not resolve costing status two, and one that went out costing the same**;
//!  * **and the boundary: a tool this port cannot run yet refuses in its OWN words, which is not
//!    a claim about the reference and is stated here so that it cannot be mistaken for one.**

use gatk_corpus as corpus;
use gatk_tools::main_entry::{self, Failure};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/main_entry.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, name: &str) -> String {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
        .unwrap_or_else(|| panic!("{kind}/{name}"))
}

/// The golden writes a shape as `out=<lines>,<first line> err=<lines>,<first line>`.
fn shape(text: &str, case: &str) -> ((usize, String), (usize, String)) {
    let row = field(text, "shape", case);
    let (out, err) = row.split_once(" err=").expect("the two streams");
    let read = |part: &str| {
        let (count, first) = part.split_once(',').expect("a count and a line");
        (count.parse::<usize>().expect("a count"), first.to_string())
    };
    (
        read(out.strip_prefix("out=").expect("the stdout shape")),
        read(err),
    )
}

/// A stream as the golden counts it: `split("\n", -1)`, so a trailing newline is a last line.
fn counted(text: &str) -> (usize, String) {
    if text.is_empty() {
        return (0, String::new());
    }
    let lines: Vec<&str> = text.split('\n').collect();
    (lines.len(), lines[0].to_string())
}

fn args(list: &[&str]) -> Vec<String> {
    list.iter().map(|arg| (*arg).to_string()).collect()
}

/// The five statuses are the golden's, and they are read off the constants and not off a process.
#[test]
fn the_statuses_are_the_reference_ones() {
    let text = golden();
    let status = |name: &str| {
        field(&text, "status", name)
            .parse::<i32>()
            .expect("a number")
    };
    assert_eq!(
        main_entry::exit_status(Failure::CommandLine),
        status("COMMANDLINE_EXCEPTION_EXIT_VALUE")
    );
    assert_eq!(
        main_entry::exit_status(Failure::User),
        status("USER_EXCEPTION_EXIT_VALUE")
    );
    assert_eq!(
        main_entry::exit_status(Failure::PicardNonZero),
        status("PICARD_TOOL_EXCEPTION")
    );
    assert_eq!(
        main_entry::exit_status(Failure::Other),
        status("ANY_OTHER_EXCEPTION_EXIT_VALUE")
    );
    assert_eq!(
        main_entry::exit_status(Failure::OutOfMemory),
        status("OUT_OF_MEMORY_EXIT_VALUE")
    );
}

/// Asking for the usage writes it to stdout, returns nothing and is not an error.
#[test]
fn asking_for_the_usage_is_not_an_error() {
    let text = golden();
    for case in ["no-arguments", "dash-h", "long-help"] {
        let (out, err) = shape(&text, case);
        let written = gatk_cli::run(&args(match case {
            "no-arguments" => &[],
            "dash-h" => &["-h"],
            _ => &["--help"],
        }));
        assert_eq!(counted(&written.stdout).1, out.1, "{case}");
        assert_eq!(counted(&written.stderr), (err.0, err.1), "{case}");
        assert_eq!(written.status, 0, "{case}");
        // `handleResult(null)` prints nothing, which is why the run that reaches it with no
        // program at all is silent past the usage.
        assert_eq!(field(&text, "result", case), "null");
        // The listing under that first line is the port's own: the names resolve, the one-line
        // summaries are not measured, and the golden is not asked to agree with the count.
        assert!(counted(&written.stdout).0 > 1, "{case}");
        assert_ne!(counted(&written.stdout).0, out.0, "{case}");
    }
}

/// The version is four lines, wherever on the command line it was asked for.
#[test]
fn the_version_is_four_lines_wherever_it_is_asked_for() {
    let text = golden();
    for (case, argv) in [
        ("version", vec!["--version"]),
        ("version-short", vec!["-version"]),
        ("version-after-a-tool", vec!["CountReads", "--version"]),
    ] {
        let (out, err) = shape(&text, case);
        let written = gatk_cli::run(&args(&argv));
        assert_eq!(counted(&written.stdout), (out.0, out.1), "{case}");
        assert_eq!(counted(&written.stderr), (err.0, err.1), "{case}");
        assert_eq!(written.status, 0, "{case}");
    }
    // And the two halves are the two the reference splits it into, which is what makes printing
    // the version to any other stream tear it in half.
    let (first, rest) = main_entry::version_lines(
        gatk_cli::TOOLKIT_NAME,
        gatk_cli::TOOLKIT_VERSION,
        gatk_cli::HTSJDK_VERSION,
        gatk_cli::PICARD_VERSION,
    );
    assert_eq!(first, field(&text, "out", "version-to-stdout"));
    assert_eq!(rest, field(&text, "out", "version-to-elsewhere"));
}

/// A name that resolves to nothing is refused with the usage and the message both.
#[test]
fn a_name_that_does_not_resolve_is_refused_with_the_usage() {
    let text = golden();
    for (case, name) in [
        ("unknown-name", "NoSuchToolAtAll"),
        ("deprecated-name", "IndelRealigner"),
    ] {
        let (out, err) = shape(&text, case);
        let written = gatk_cli::run(&args(&[name]));
        assert_eq!(counted(&written.stdout), (out.0, out.1), "{case}");
        // The refusal goes to stderr, starting with the usage's own first line.
        assert_eq!(counted(&written.stderr).1, err.1, "{case}");
        assert_eq!(
            written.status,
            main_entry::exit_status(Failure::User),
            "{case}"
        );
        // And the message under it is the reference's, word for word.
        let message = field(&text, "error", case)
            .split_once("UserException:")
            .expect("the exception")
            .1
            .to_string();
        assert!(written.stderr.contains(&message), "{case}");
    }
}

/// The two handlers write what the golden says they write.
#[test]
fn the_handlers_are_the_reference_ones() {
    let text = golden();
    let mut printed = String::new();
    for result in [None, Some("0"), Some("a string")] {
        if let Some(line) = main_entry::tool_returned(result) {
            printed.push_str(&line);
            printed.push('\n');
        }
    }
    assert_eq!(printed, field(&text, "out", "handlers"));
    assert_eq!(
        main_entry::user_exception_report("the message a refusal carries"),
        field(&text, "err", "handlers")
    );
}

/// The first tool a command line reaches: the parse is the reference's, decision for decision.
#[test]
fn a_tool_that_is_no_walker_parses() {
    let declarations = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/tool_declarations.txt.gz"),
    );
    let outcome = |case: &str| {
        let prefix = format!("parse\tIndexFeatureFile\t{case}\t");
        declarations
            .lines()
            .find(|line| line.starts_with(&prefix))
            .map(|line| unescape(&line[prefix.len()..]))
            .unwrap_or_else(|| panic!("parse/{case}"))
    };
    // The golden's own three command lines, through the binary rather than through the parser.
    assert!(outcome("no-arguments").contains("Argument input was missing"));
    let refused = gatk_cli::run(&args(&["IndexFeatureFile"]));
    assert!(
        refused
            .stderr
            .contains("Argument input was missing: Argument 'input' is required"),
        "{}",
        refused.stderr
    );
    // A parse that fails is a CommandLineException, which is status one and not two.
    assert_eq!(
        refused.status,
        main_entry::exit_status(Failure::CommandLine)
    );
    // `mainEntry` prints the tool's own usage above that message, and so does the port: the usage
    // is composed from the same declarations the parser is built from, and it is compared against
    // the usage golden below.
    assert!(refused
        .stderr
        .starts_with("USAGE: IndexFeatureFile [arguments]"));
    // A command line the reference accepts is accepted, and the port then refuses it for its own
    // reason, which is that it cannot RUN the tool.
    assert!(outcome("input-only").ends_with("ok"));
    let accepted = gatk_cli::run(&args(&["IndexFeatureFile", "-I", "/dev/null"]));
    assert!(!accepted.stderr.contains("Argument input was missing"));
    assert!(accepted.stderr.contains("this port does not carry yet"));
    assert_eq!(accepted.status, main_entry::exit_status(Failure::User));
    // And an argument the tool does not declare is refused by the parser, as the golden says.
    assert!(outcome("an-interval").contains("not a recognized option"));
    let unknown = gatk_cli::run(&args(&[
        "IndexFeatureFile",
        "-I",
        "/dev/null",
        "-L",
        "chr1",
    ]));
    assert!(
        unknown.stderr.contains("is not a recognized option"),
        "{}",
        unknown.stderr
    );
    assert_eq!(
        unknown.status,
        main_entry::exit_status(Failure::CommandLine)
    );
}

/// A tool asked for its own help answers with its usage, byte for byte.
#[test]
fn a_tool_asked_for_help_answers_with_its_usage() {
    let text = golden();
    let usage_golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/usage_text.txt.gz"),
    );
    let prefix = "usage\tIndexFeatureFile\t";
    let expected = usage_golden
        .lines()
        .find(|line| line.starts_with(prefix))
        .map(|line| unescape(&line[prefix.len()..]))
        .expect("the tool's usage");
    for flag in ["-h", "--help"] {
        let written = gatk_cli::run(&args(&["IndexFeatureFile", flag]));
        // The tool's usage goes to STDERR and the run is a success, which is what the golden's
        // `help-after-a-tool` shape says of the tool it measured.
        assert_eq!(written.stderr, expected, "{flag}");
        assert!(written.stdout.is_empty(), "{flag}");
        assert_eq!(written.status, 0, "{flag}");
    }
    let (out, err) = shape(&text, "help-after-a-tool");
    assert_eq!(out.0, 0);
    assert!(err.1.starts_with("USAGE: CountReads [arguments]"));
    // The tool the golden measured is a walker, whose usage carries the conditional blocks one
    // per read filter the descriptor discovered. The port has no descriptor, so it lays out no
    // usage for it and the dispatcher falls through to its own refusal.
    assert_eq!(gatk_cli::tool_usage("CountReads"), None);
    let walker = gatk_cli::run(&args(&["CountReads", "-h"]));
    assert!(walker.stderr.contains("this port does not carry yet"));
}

/// A tool this port cannot run yet refuses in its own words, and says so.
#[test]
fn a_tool_this_port_cannot_run_refuses_in_its_own_words() {
    let text = golden();
    let written = gatk_cli::run(&args(&["CountReads"]));
    // The reference gets past the dispatch here and is refused by the PARSER, which is status one
    // and a CommandLineException. This port has no declarations to parse against yet, so it stops
    // one step earlier, and its refusal carries neither that status nor that message.
    let reference = field(&text, "error", "gatk-tool-no-arguments");
    assert!(reference.contains("Argument input was missing"));
    assert!(!written.stderr.contains("Argument input was missing"));
    assert_ne!(
        written.status,
        main_entry::exit_status(Failure::CommandLine)
    );
    assert!(written.stderr.contains("this port does not carry yet"));
    // What IS the reference's is the step before: the name resolved, so the dispatcher did not
    // refuse it the way it refuses a name that does not.
    assert!(!written.stderr.contains("is not a valid command"));
    assert!(gatk_cli::runner("CountReads").is_none());
}
