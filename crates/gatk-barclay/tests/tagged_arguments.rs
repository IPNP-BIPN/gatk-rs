//! Conformance for Barclay's tagged arguments and its collection-file expansion, against the
//! oracle.
//!
//! Golden from `tools/argument-conformance/BarclayTaggedArgumentDump.java`.
//!
//! # What the golden settles
//!
//! ```text
//! result  same-option-tag-and-value-twice  ... "tagged-collection:one a.bam" was duplicated ...
//! result  same-tag-different-values        ok
//! result  tag-on-plain-scalar              The argument: "/plain-scalar" does not accept tags ...
//! field   expand-list       plain-collection  [first, second, third]
//! field   no-expand-scalar  plain-scalar      fixtures/values.list
//! tag     expand-tagged     tagged-collection 2  one
//! ```
//!
//! The first pair is the surrogate key doing double duty as a uniqueness test: it is built from the
//! option string **and** the value, so the same tag twice with the same value collides and the same
//! tag twice with different values does not. The third row is an argument whose short name is
//! empty reporting itself as `/plain-scalar`, because the message is built from
//! `getShortName() + "/" + getFullName()` with no guard. The last three are expansion: a
//! collection value ending in `.list` becomes the file's lines, the identical value on a scalar
//! stays a path, and a tag is written onto **every** value the file produced.

use std::collections::HashMap;

use gatk_barclay::{Annotation, Definition, Error, FileSource, IoError, Parser, Value, ValueClass};
use gatk_corpus as corpus;

/// The expansion files the dump wrote, by the path it wrote them to.
///
/// Given to the parser through [`FileSource`] rather than written to disk: these are the bytes the
/// reference read, and a file rebuilt here to match a description of them would be a second
/// fixture.
struct Fixtures(HashMap<&'static str, &'static str>);

impl Fixtures {
    fn new() -> Self {
        let body = "first\n\n# a comment\n  second  \nthird\n";
        let mut files = HashMap::new();
        files.insert("fixtures/values.list", body);
        files.insert("fixtures/values.args", body);
        files.insert("fixtures/values.txt", body);
        files.insert("fixtures/looks-like-a-dict.list", "@HD\tVN:1.6\nchr1\n");
        files.insert("fixtures/empty.list", "\n# only a comment\n");
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
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/tagged_arguments.txt.gz"),
    )
}

/// The five fields of `BarclayTaggedArgumentDump.Args`, in declaration order.
fn definitions() -> Vec<Definition> {
    vec![
        Definition::new(
            Annotation {
                full_name: "tagged-collection",
                doc: "a collection that accepts tags",
                optional: true,
                ..Annotation::default()
            },
            "taggedCollection",
            ValueClass::Tagged,
            true,
            false,
            Value::List(Vec::new()),
        ),
        Definition::new(
            Annotation {
                full_name: "tagged-scalar",
                short_name: "T",
                doc: "a scalar that accepts tags",
                optional: true,
                ..Annotation::default()
            },
            "taggedScalar",
            ValueClass::Tagged,
            false,
            false,
            Value::Null,
        ),
        Definition::new(
            Annotation {
                full_name: "plain-collection",
                doc: "a collection of plain strings",
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
                full_name: "plain-scalar",
                doc: "a scalar that does not accept tags",
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
                full_name: "no-expansion",
                doc: "a collection that refuses expansion",
                optional: true,
                suppress_file_expansion: true,
                ..Annotation::default()
            },
            "noExpansion",
            ValueClass::Text,
            true,
            false,
            Value::List(Vec::new()),
        ),
    ]
}

const REPORTED: [&str; 5] = [
    "tagged-collection",
    "tagged-scalar",
    "plain-collection",
    "plain-scalar",
    "no-expansion",
];

/// `tag\t<label>\t<field>\t<index>\t<tag>\t<attributes>`, the shape the dump reports a populated
/// tag in. `String.valueOf(null)` is the four characters, and the attributes are joined with `,`
/// in insertion order.
fn tag_row(label: &str, field: &str, index: usize, value: &Value) -> Option<String> {
    let Value::Tagged {
        tag, attributes, ..
    } = value
    else {
        return None;
    };
    let rendered: Vec<String> = attributes
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect();
    Some(format!(
        "tag\t{label}\t{field}\t{index}\t{}\t{}",
        tag.clone().unwrap_or_else(|| "null".to_string()),
        rendered.join(",")
    ))
}

fn render_error(error: &Error) -> String {
    // The dump escapes the one message that embeds a newline, so the port does too.
    format!("E:{}:{}", error.class, error.message.replace('\n', "\\n"))
}

#[test]
fn every_tag_and_every_expansion_is_read_the_way_the_reference_reads_it() {
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
                for name in REPORTED {
                    let value = parser.value_of(name).expect("a declared field");
                    produced.push(format!(
                        "field\t{label}\t{name}\t{}",
                        value.to_java_string()
                    ));
                }
                // The tag rows come after every field row, collection first, then the scalar.
                if let Some(Value::List(values)) = parser.value_of("tagged-collection") {
                    for (index, value) in values.iter().enumerate() {
                        if let Some(row) = tag_row(label, "tagged-collection", index, value) {
                            produced.push(row);
                        }
                    }
                }
                if let Some(value @ Value::Tagged { .. }) = parser.value_of("tagged-scalar") {
                    if let Some(row) = tag_row(label, "tagged-scalar", 0, value) {
                        produced.push(row);
                    }
                }
            }
            Err(error) => produced.push(format!("result\t{label}\t{}", render_error(&error))),
        }
    }

    assert_eq!(produced.len(), expected.len(), "row count");
    for (index, (produced, oracle)) in produced.iter().zip(expected.iter()).enumerate() {
        assert_eq!(produced, oracle, "row {index}");
    }
}

/// The surrogate key is the option string **and** the value, so it collides only on both.
#[test]
fn the_surrogate_key_is_what_makes_a_duplicate() {
    let files = Fixtures::new();

    let mut same = Parser::new(definitions());
    let error = same
        .parse_arguments_with(
            &[
                "--tagged-collection:one",
                "a.bam",
                "--tagged-collection:one",
                "a.bam",
            ],
            &files,
        )
        .expect_err("the same option, tag and value twice is a duplicate");
    assert!(error.message.contains("was duplicated on the command line"));

    let mut different = Parser::new(definitions());
    different
        .parse_arguments_with(
            &[
                "--tagged-collection:one",
                "a.bam",
                "--tagged-collection:one",
                "b.bam",
            ],
            &files,
        )
        .expect("the same tag with two values is two values");
}

/// Expansion is a collection-only mechanism, and `suppressFileExpansion` turns it off.
#[test]
fn only_a_collection_expands() {
    let files = Fixtures::new();

    let mut collection = Parser::new(definitions());
    collection
        .parse_arguments_with(&["--plain-collection", "fixtures/values.list"], &files)
        .expect("a collection expands");
    assert_eq!(
        collection
            .value_of("plain-collection")
            .unwrap()
            .to_java_string(),
        "[first, second, third]"
    );

    let mut scalar = Parser::new(definitions());
    scalar
        .parse_arguments_with(&["--plain-scalar", "fixtures/values.list"], &files)
        .expect("a scalar takes the path");
    assert_eq!(
        scalar.value_of("plain-scalar").unwrap().to_java_string(),
        "fixtures/values.list"
    );

    let mut suppressed = Parser::new(definitions());
    suppressed
        .parse_arguments_with(&["--no-expansion", "fixtures/values.list"], &files)
        .expect("expansion is suppressed");
    assert_eq!(
        suppressed
            .value_of("no-expansion")
            .unwrap()
            .to_java_string(),
        "[fixtures/values.list]"
    );
}
