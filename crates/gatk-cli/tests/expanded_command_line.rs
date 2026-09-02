//! Conformance for `getCommandLine()` against GATK 4.6.2.0, over command lines the dump parsed.
//!
//! Golden from `tools/argument-conformance/ExpandedCommandLineDump.java`. A BAM's `@PG` carries
//! this string in `CL` and a VCF's `##GATKCommandLine` carries it in `CommandLine`, so it is a
//! byte of the output for every tool that writes either.
//!
//! # What this suite is for
//!
//!  * **the two groups**: the arguments the user set, in the parser's declaration order, then
//!    those that were not set but have a non-null default, in that same order;
//!  * **the long form always**, so a short alias is expanded;
//!  * **`--name=value` not being a form the parser accepts at all**;
//!  * **a collection being one pair per element**;
//!  * **an argument set to its own default moving to the FIRST group**, which changes where it
//!    appears;
//!  * **a flag given with no value printing one**;
//!  * **and a default of `null` being omitted where an EMPTY default is printed with a trailing
//!    space.**
//!
//! While the suite is `golden-pending` the dump is named by `EXPANDED_COMMAND_LINE_DUMP`.

use gatk_cli::command_line::expanded;

/// The command lines the harness ran, in its own order.
const CASES: &[(&str, &str, &[&str])] = &[
    (
        "PrintReads",
        "input-and-output",
        &["-I", "/dev/null", "-O", "/dev/null"],
    ),
    (
        "PrintReads",
        "long-names",
        &["--input", "/dev/null", "--output", "/dev/null"],
    ),
    (
        "PrintReads",
        "equals-form",
        &["--input=/dev/null", "--output=/dev/null"],
    ),
    (
        "PrintReads",
        "two-inputs",
        &["-I", "/dev/null", "-I", "/dev/null", "-O", "/dev/null"],
    ),
    (
        "PrintReads",
        "default-set-explicitly",
        &[
            "-I",
            "/dev/null",
            "-O",
            "/dev/null",
            "--create-output-bam-index",
            "true",
        ],
    ),
    (
        "PrintReads",
        "flag-without-a-value",
        &[
            "-I",
            "/dev/null",
            "-O",
            "/dev/null",
            "--add-output-sam-program-record",
        ],
    ),
    (
        "PrintReads",
        "with-an-interval",
        &["-I", "/dev/null", "-O", "/dev/null", "-L", "chr1:1-100"],
    ),
    ("CountReads", "one-input", &["-I", "/dev/null"]),
    ("CountVariants", "one-variant", &["-V", "/dev/null"]),
    ("IndexFeatureFile", "one-input", &["-I", "/dev/null"]),
    (
        "CreateHadoopBamSplittingIndex",
        "one-input",
        &["-I", "/dev/null"],
    ),
];

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t")
        .replace("\\n", "\n")
        .replace("\\\\", "\\")
}

fn field(dump: &str, kind: &str, tool: &str, label: &str) -> Option<String> {
    let prefix = format!("{kind}\t{tool}\t{label}\t");
    dump.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

#[test]
fn every_command_line_expands_as_the_reference_expands_it() {
    let dump = match std::env::var("EXPANDED_COMMAND_LINE_DUMP") {
        Ok(path) => {
            std::fs::read_to_string(path).expect("the dump named by EXPANDED_COMMAND_LINE_DUMP")
        }
        Err(_) => {
            println!(
                "skipped: the expanded-command-line golden is still pending. Run the suite and \
                 point EXPANDED_COMMAND_LINE_DUMP at \
                 tools/conformance/pending/expanded-command-line.ExpandedCommandLineDump.txt"
            );
            return;
        }
    };

    for (tool, label, argv) in CASES {
        let args: Vec<String> = argv.iter().map(|arg| (*arg).to_string()).collect();
        let parsed = gatk_cli::parse_for(tool, &args);
        match parsed {
            Ok(parser) => {
                let ours = expanded(tool, &parser);
                assert_eq!(
                    ours,
                    field(&dump, "line", tool, label)
                        .unwrap_or_else(|| panic!("{tool}/{label}: the golden refused")),
                    "{tool}/{label}"
                );
            }
            Err(message) => {
                let theirs = field(&dump, "error", tool, label).unwrap_or_else(|| {
                    panic!("{tool}/{label}: the port refused, the golden did not")
                });
                // The golden carries the exception class and this carries the message, so what is
                // compared is that the reference refused with THIS wording.
                assert!(
                    theirs.contains(&message),
                    "{tool}/{label}: {theirs} vs {message}"
                );
            }
        }
    }
}
