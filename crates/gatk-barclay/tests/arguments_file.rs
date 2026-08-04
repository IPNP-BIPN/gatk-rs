//! Conformance for Barclay's `--arguments_file`, against the oracle.
//!
//! Golden from `tools/argument-conformance/BarclayArgumentsFileDump.java`.
//!
//! # What the golden settles
//!
//! ```text
//! field  file-then-command-line  plain-collection  [a, b, cli]
//! field  command-line-then-file  plain-collection  [a, b, cli]
//! result scalar-in-both          ... Argument 'plain-scalar' cannot be specified more than once.
//! field  self-referencing-file      plain-collection  [self]
//! field  mutually-referencing-files plain-collection  [pong, ping]
//! ```
//!
//! The first two rows are the finding: the file's arguments come first **wherever
//! `--arguments_file` sat on the command line**, because the expansion is what the original argv
//! is appended to. So a scalar given in both places is a duplicate rather than an override, and a
//! collection is ordered by where the value came from rather than by where it was written.
//!
//! The last two are the recursion's bound. It is a set of file names, not a depth limit, and every
//! file *named* goes into it — including one skipped for already being there.

use std::collections::HashMap;

use gatk_barclay::{
    Annotation, Definition, Error, FileSource, IoError, Parser, Value, ValueClass,
    ARGUMENTS_FILE_FULLNAME,
};
use gatk_corpus as corpus;

/// The nine files the dump wrote, by the path it wrote them to.
struct Fixtures(HashMap<&'static str, &'static str>);

impl Fixtures {
    fn new() -> Self {
        let mut files = HashMap::new();
        files.insert("fixtures/scalar.txt", "--plain-scalar fromfile\n");
        files.insert(
            "fixtures/collection.txt",
            "--plain-collection a --plain-collection b\n",
        );
        files.insert(
            "fixtures/messy.txt",
            "# a comment\n\n   --plain-collection    one\t\t--plain-collection two   \n\n",
        );
        files.insert("fixtures/flag.txt", "--flag\n");
        files.insert(
            "fixtures/outer.txt",
            "--arguments_file fixtures/inner.txt\n--plain-collection outer\n",
        );
        files.insert("fixtures/inner.txt", "--plain-collection inner\n");
        files.insert(
            "fixtures/self.txt",
            "--arguments_file fixtures/self.txt\n--plain-collection self\n",
        );
        files.insert(
            "fixtures/ping.txt",
            "--arguments_file fixtures/pong.txt\n--plain-collection ping\n",
        );
        files.insert(
            "fixtures/pong.txt",
            "--arguments_file fixtures/ping.txt\n--plain-collection pong\n",
        );
        Fixtures(files)
    }
}

impl FileSource for Fixtures {
    fn read(&self, path: &str) -> Result<String, IoError> {
        self.0.get(path).map(|text| text.to_string()).ok_or(IoError)
    }
}

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/arguments_file.txt.gz"),
    )
}

/// `BarclayArgumentsFileDump.Args`, in declaration order: the special collection first, since it
/// is the first field, and its own `arguments_file` is therefore definition zero.
fn definitions() -> Vec<Definition> {
    vec![
        Definition::new(
            Annotation {
                full_name: ARGUMENTS_FILE_FULLNAME,
                doc: "read one or more arguments files and add them to the command line",
                optional: true,
                ..Annotation::default()
            },
            "ARGUMENTS_FILE",
            ValueClass::Text,
            true,
            false,
            Value::List(Vec::new()),
        ),
        Definition::new(
            Annotation {
                full_name: "plain-scalar",
                doc: "a scalar",
                optional: true,
                ..Annotation::default()
            },
            "plainScalar",
            ValueClass::Text,
            false,
            false,
            Value::Null,
        ),
        Definition::new(
            Annotation {
                full_name: "plain-collection",
                doc: "a collection",
                optional: true,
                ..Annotation::default()
            },
            "plainCollection",
            ValueClass::Text,
            true,
            false,
            Value::List(Vec::new()),
        ),
        Definition::new(
            Annotation {
                full_name: "flag",
                doc: "a flag",
                optional: true,
                ..Annotation::default()
            },
            "flag",
            ValueClass::Boolean,
            false,
            true,
            Value::Bool(false),
        ),
    ]
}

fn render_error(error: &Error) -> String {
    format!("E:{}:{}", error.class, error.message.replace('\n', "\\n"))
}

#[test]
fn every_expansion_is_the_one_the_reference_performs() {
    let text = golden();
    let files = Fixtures::new();

    let mut produced: Vec<String> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        expected.push(line);

        let Some(rest) = line.strip_prefix("case\t") else {
            continue;
        };
        let mut parts = rest.splitn(2, '\t');
        let label = parts.next().expect("a label");
        let argv_text = parts.next().unwrap_or("");
        let argv: Vec<&str> = if argv_text.is_empty() {
            Vec::new()
        } else {
            argv_text.split(' ').collect()
        };

        produced.push(format!("case\t{label}\t{argv_text}"));

        let mut parser = Parser::new(definitions());
        match parser.parse_arguments_with(&argv, &files) {
            Ok(()) => {
                produced.push(format!("result\t{label}\tok"));
                for name in ["plain-scalar", "plain-collection", "flag"] {
                    let value = parser.value_of(name).expect("a declared field");
                    produced.push(format!(
                        "field\t{label}\t{name}\t{}",
                        value.to_java_string()
                    ));
                }
                // The dump reports `ARGUMENTS_FILE` last, and it holds whatever the *final* pass
                // parsed rather than the union of every pass.
                let value = parser
                    .value_of(ARGUMENTS_FILE_FULLNAME)
                    .expect("a declared field");
                produced.push(format!(
                    "field\t{label}\targuments_file\t{}",
                    value.to_java_string()
                ));
            }
            Err(error) => produced.push(format!("result\t{label}\t{}", render_error(&error))),
        }
    }

    assert_eq!(produced.len(), expected.len(), "row count");
    for (index, (produced, oracle)) in produced.iter().zip(expected.iter()).enumerate() {
        assert_eq!(produced, oracle, "row {index}");
    }
}

/// The file's arguments are prepended, wherever the argument that named the file sat.
#[test]
fn the_files_arguments_come_first_either_way() {
    let files = Fixtures::new();
    for argv in [
        vec![
            "--arguments_file",
            "fixtures/collection.txt",
            "--plain-collection",
            "cli",
        ],
        vec![
            "--plain-collection",
            "cli",
            "--arguments_file",
            "fixtures/collection.txt",
        ],
    ] {
        let mut parser = Parser::new(definitions());
        parser
            .parse_arguments_with(&argv, &files)
            .expect("the expansion parses");
        assert_eq!(
            parser
                .value_of("plain-collection")
                .unwrap()
                .to_java_string(),
            "[a, b, cli]"
        );
    }
}

/// The recursion is bounded by a set of names, so a cycle of any length is read once round.
#[test]
fn a_cycle_is_read_once() {
    let files = Fixtures::new();

    let mut parser = Parser::new(definitions());
    parser
        .parse_arguments_with(&["--arguments_file", "fixtures/self.txt"], &files)
        .expect("a self-referencing file terminates");
    assert_eq!(
        parser
            .value_of("plain-collection")
            .unwrap()
            .to_java_string(),
        "[self]"
    );

    let mut parser = Parser::new(definitions());
    parser
        .parse_arguments_with(&["--arguments_file", "fixtures/ping.txt"], &files)
        .expect("a two-file cycle terminates");
    // `pong` before `ping`: the second pass prepends pong.txt's expansion to a command line that
    // already began with ping.txt's.
    assert_eq!(
        parser
            .value_of("plain-collection")
            .unwrap()
            .to_java_string(),
        "[pong, ping]"
    );
}
