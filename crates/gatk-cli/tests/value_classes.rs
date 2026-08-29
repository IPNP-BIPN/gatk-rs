//! Conformance for the four value classes against GATK 4.6.2.0 and Barclay 5.0.0.
//!
//! Golden from `tools/argument-conformance/ToolArgumentValueClassDump.java`, which parsed
//! thirty-four command lines against one target carrying a field of each class.
//!
//! # What this suite is for
//!
//!  * **a class built from a string accepting every string, so a bad path is not a bad value**;
//!  * **a tag being a tag on one class and an error on another**;
//!  * **`Float.valueOf`'s grammar, which is not `str::parse`'s in either direction**;
//!  * **and the default rendering, which is what decides optionality.**

use gatk_barclay::{java_float, java_long, Annotation, Definition, Parser, Value, ValueClass};
use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/tool_argument_value_classes.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn row(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

fn outcome(text: &str, case: &str) -> String {
    row(text, "result", case).unwrap_or_else(|| panic!("result/{case}"))
}

fn field(text: &str, case: &str, name: &str) -> String {
    let prefix = format!("field\t{case}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
        .unwrap_or_else(|| panic!("field/{case}/{name}"))
}

/// One definition of the given class, named the way the dump's target names it.
fn definition(long_name: &'static str, short_name: &'static str, class: ValueClass) -> Definition {
    Definition::new(
        Annotation {
            full_name: long_name,
            short_name,
            doc: "a value",
            optional: true,
            ..Annotation::default()
        },
        long_name,
        class,
        false,
        false,
        Value::Null,
    )
}

fn path() -> Definition {
    definition(
        "path",
        "P",
        ValueClass::Constructed {
            simple_name: "GATKPath",
            taggable: true,
        },
    )
}

fn file() -> Definition {
    definition(
        "file",
        "F",
        ValueClass::Constructed {
            simple_name: "File",
            taggable: false,
        },
    )
}

/// A class built from a string accepts every string, so no bad value exists for it.
#[test]
fn a_constructed_class_accepts_every_string() {
    let text = golden();
    for case in [
        "a-path",
        "a-uri",
        "a-relative-path",
        "a-path-that-is-not-there",
        "an-empty-path",
        "a-file",
        "a-file-that-is-not-there",
        "an-empty-file",
    ] {
        assert_eq!(outcome(&text, case), "ok", "{case}");
    }
    // And the value is the string it was given rather than a resolved path.
    assert_eq!(field(&text, "a-uri", "path"), "gs://bucket/reads.bam");
    assert_eq!(field(&text, "a-relative-path", "path"), "reads.bam");
    let mut parser = Parser::new(vec![path(), file()]);
    assert!(parser
        .parse_arguments(&["--path", "gs://bucket/reads.bam", "--file", ""])
        .is_ok());
    assert_eq!(
        parser.definitions()[0].value,
        Value::Tagged {
            value: "gs://bucket/reads.bam".to_string(),
            tag: None,
            attributes: Vec::new(),
        }
    );
    assert_eq!(parser.definitions()[1].value, Value::Str(String::new()));
}

/// A tag is a tag on the taggable class and an error on the one that is not.
#[test]
fn a_tag_is_a_tag_on_one_class_and_an_error_on_the_other() {
    let text = golden();
    assert_eq!(outcome(&text, "a-tagged-path"), "ok");
    // The tag is part of the field's own rendering, not only of its accessors.
    assert_eq!(
        field(&text, "a-tagged-path", "path"),
        "/tmp/reads.bam:name,key=value"
    );
    let mut parser = Parser::new(vec![path()]);
    assert!(parser
        .parse_arguments(&["--path:name,key=value", "/tmp/reads.bam"])
        .is_ok());
    assert_eq!(
        parser.definitions()[0].value,
        Value::Tagged {
            value: "/tmp/reads.bam".to_string(),
            tag: Some("name".to_string()),
            attributes: vec![("key".to_string(), "value".to_string())],
        }
    );
    // The same command line against the untaggable class is refused, and the message names the
    // argument as `shortName/fullName`.
    let refusal = outcome(&text, "a-tagged-file");
    assert!(
        refusal.contains("The argument: \"F/file\" does not accept tags: \"name\""),
        "{refusal}"
    );
    let mut parser = Parser::new(vec![file()]);
    let error = parser
        .parse_arguments(&["--file:name", "/tmp/reads.bam"])
        .expect_err("the refusal");
    assert!(refusal.ends_with(&error.message), "{}", error.message);
}

/// `Float.valueOf` is not `str::parse`, and the golden holds both sides of the difference.
#[test]
fn the_float_grammar_is_javas() {
    let text = golden();
    let accepted = [
        ("a-fraction", "0.25", 0.25f32),
        ("a-fraction-in-exponent-form", "1e3", 1000.0),
        ("a-fraction-with-a-float-suffix", "1.5f", 1.5),
        ("a-fraction-with-a-double-suffix", "1.5d", 1.5),
        ("a-fraction-in-hexadecimal", "0x1p3", 8.0),
        ("a-fraction-with-a-leading-space", " 1.5", 1.5),
        ("a-fraction-with-a-trailing-space", "1.5 ", 1.5),
        ("a-fraction-with-a-leading-plus", "+1.5", 1.5),
    ];
    for (case, spelling, value) in accepted {
        assert_eq!(outcome(&text, case), "ok", "{case}");
        assert_eq!(java_float(spelling), Some(value), "{spelling}");
        // And the reference's own rendering of the field agrees with the port's value.
        let written: f32 = field(&text, case, "fraction").parse().expect("a number");
        assert_eq!(written, value, "{case}");
    }
    // The two capitalised words, and the value that overflows into one of them.
    assert_eq!(java_float("Infinity"), Some(f32::INFINITY));
    assert_eq!(java_float("-Infinity"), Some(f32::NEG_INFINITY));
    assert!(java_float("NaN").expect("a value").is_nan());
    assert_eq!(
        field(&text, "a-fraction-out-of-a-floats-range", "fraction"),
        "Infinity"
    );
    assert_eq!(java_float("1e40"), Some(f32::INFINITY));
    // The three refusals, two of which Rust's own parse accepts.
    for (case, spelling) in [
        ("a-fraction-spelled-inf", "inf"),
        ("a-fraction-spelled-nan-in-lower-case", "nan"),
        ("a-fraction-with-an-underscore", "1_000"),
        ("a-fraction-that-is-not-a-number", "abc"),
    ] {
        let refusal = outcome(&text, case);
        assert!(
            refusal.contains(&format!(
                "Failure constructing 'Float' from the string '{spelling}'."
            )),
            "{refusal}"
        );
        assert_eq!(java_float(spelling), None, "{spelling}");
    }
    assert!("inf".parse::<f32>().is_ok());
    assert!("nan".parse::<f32>().is_ok());
    // A double's message names Double where this one names Float, which is why they are two
    // classes and not one.
    let mut parser = Parser::new(vec![definition("fraction", "FR", ValueClass::Float)]);
    let error = parser
        .parse_arguments(&["--fraction", "abc"])
        .expect_err("the refusal");
    assert!(error.message.contains("'Float'"), "{}", error.message);
}

/// A `Long` is neither an `Integer` with a wider range nor a `Float` with a narrower grammar.
#[test]
fn the_long_grammar_is_its_own() {
    let text = golden();
    // The range: an int cannot hold this and a long can, and one past a long is refused.
    assert_eq!(outcome(&text, "a-count-past-an-ints-limit"), "ok");
    assert_eq!(java_long("2147483648"), Some(2147483648));
    assert_eq!(
        field(&text, "a-count-at-a-longs-limit", "count"),
        i64::MAX.to_string()
    );
    assert_eq!(java_long("9223372036854775807"), Some(i64::MAX));
    let past = outcome(&text, "a-count-past-a-longs-limit");
    assert!(
        past.contains("Failure constructing 'Long' from the string '9223372036854775808'."),
        "{past}"
    );
    assert_eq!(java_long("9223372036854775808"), None);
    // The grammar: a leading plus is taken, and a leading space and a hexadecimal literal are not,
    // which is where it parts company with the float.
    assert_eq!(outcome(&text, "a-count-with-a-plus"), "ok");
    assert_eq!(java_long("+42"), Some(42));
    for (case, spelling) in [
        ("a-count-with-a-space", " 42"),
        ("a-count-in-hexadecimal", "0x2a"),
        ("a-count-that-is-not-a-number", "abc"),
    ] {
        assert!(
            outcome(&text, case).contains("Failure constructing 'Long'"),
            "{case}"
        );
        assert_eq!(java_long(spelling), None, "{spelling}");
    }
    // The float takes both of the first two, which is what makes them two grammars.
    assert_eq!(java_float(" 1.5"), Some(1.5));
    assert_eq!(java_float("0x1p3"), Some(8.0));
    // And a refusal names the class it failed to build, which is how the two are told apart.
    let mut parser = Parser::new(vec![definition("count", "C", ValueClass::Long)]);
    let error = parser
        .parse_arguments(&["--count", "abc"])
        .expect_err("the refusal");
    assert!(error.message.contains("'Long'"), "{}", error.message);
}

/// The default rendering is `String.valueOf(field)`, which is what decides optionality.
#[test]
fn the_default_rendering_is_the_fields_own() {
    let text = golden();
    let written = |name: &str| {
        let prefix = format!("default\t{name}\t");
        text.lines()
            .find(|line| line.starts_with(&prefix))
            .map(|line| unescape(&line[prefix.len()..]))
            .unwrap_or_else(|| panic!("default/{name}"))
    };
    assert_eq!(written("path"), "null");
    assert_eq!(written("file"), "null");
    assert_eq!(written("initialised-file"), "already/here");
    assert_eq!(written("primitive-fraction"), "0.5");
    // An empty collection renders as null, which is the reference's own collapse of the two.
    assert_eq!(written("paths"), "null");
    let uninitialised = path();
    assert_eq!(uninitialised.default_value_as_string(), "null");
    let initialised = Definition::new(
        Annotation {
            full_name: "initialised-file",
            doc: "a value",
            optional: true,
            ..Annotation::default()
        },
        "initialised-file",
        ValueClass::Constructed {
            simple_name: "File",
            taggable: false,
        },
        false,
        false,
        Value::Str("already/here".to_string()),
    );
    assert_eq!(initialised.default_value_as_string(), "already/here");
    assert!(initialised.is_optional());
}

/// A feature input's tag is its name, which is what a walker looks the feature up by.
#[test]
fn a_feature_inputs_tag_is_its_name() {
    let text = golden();
    assert_eq!(field(&text, "a-feature", "feature-name"), "/tmp/sites.vcf");
    assert_eq!(field(&text, "a-feature", "feature-path"), "/tmp/sites.vcf");
    assert_eq!(field(&text, "a-tagged-feature", "feature-name"), "known");
    assert_eq!(
        field(&text, "a-tagged-feature", "feature-path"),
        "/tmp/sites.vcf"
    );
    // The port carries the tag and the value apart, which is what the two answers are.
    let mut parser = Parser::new(vec![definition(
        "feature",
        "",
        ValueClass::Constructed {
            simple_name: "FeatureInput",
            taggable: true,
        },
    )]);
    assert!(parser
        .parse_arguments(&["--feature:known", "/tmp/sites.vcf"])
        .is_ok());
    assert_eq!(
        parser.definitions()[0].value,
        Value::Tagged {
            value: "/tmp/sites.vcf".to_string(),
            tag: Some("known".to_string()),
            attributes: Vec::new(),
        }
    );
}
