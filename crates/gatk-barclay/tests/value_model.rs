//! Conformance for Barclay's argument value model, against the oracle.
//!
//! Golden from `tools/argument-conformance/BarclayValueModelDump.java`, the real
//! `CommandLineArgumentParser` over the same thirteen `@Argument` fields declared below.
//!
//! # What the golden settles
//!
//! ```text
//! result  scalar-equals-syntax   E:CommandLineException:Can't parse option name containing an embedded '=' ...
//! result  bounded-below          ... has a bad value: 0. allowed range [1, 10].
//! result  bounded-null           ... has a bad value: null. allowed range [1.0, 10.0].
//! result  recommended-far-out    ok
//! result  flag-bad-value         ... Positional arguments were provided ',maybe}' ...
//! field   collection-one         collection  [a]        the declared default is gone
//! ```
//!
//! The two `bounded-` rows are the same argument and the same bounds, formatted two ways, because
//! the formatter asks whether the *value* is an `Integer` and a null is not. `recommended-far-out`
//! is a value far outside a recommended range that produces no error and no warning, which is what
//! makes that branch nearly unreachable.

use gatk_barclay::{Annotation, Definition, Error, Parser, Value, ValueClass};
use gatk_corpus as corpus;

/// The enum `BarclayValueModelDump.Mode`, whose constants appear in the message a bad value gets.
const MODE: ValueClass = ValueClass::Enum {
    simple_name: "Mode",
    constants: &["FAST", "SLOW"],
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/value_model.txt.gz"),
    )
}

/// The thirteen fields of `BarclayValueModelDump.Args`, in declaration order, which is the order
/// values are propagated in and the order `validateArgumentValues` reports a missing required
/// argument in.
fn definitions() -> Vec<Definition> {
    vec![
        // Uninitialised and not marked optional: required.
        Definition::new(
            Annotation {
                full_name: "required-string",
                short_name: "R",
                doc: "required, uninitialised",
                ..Annotation::default()
            },
            "requiredString",
            ValueClass::Text,
            false,
            false,
            Value::Null,
        ),
        // Initialised to an **empty** collection, which `convertDefaultValueToString` maps back to
        // "null": required, despite being initialised.
        Definition::new(
            Annotation {
                full_name: "required-collection",
                doc: "required: an empty collection reads as null",
                ..Annotation::default()
            },
            "requiredCollection",
            ValueClass::Text,
            true,
            false,
            Value::List(Vec::new()),
        ),
        Definition::new(
            Annotation {
                full_name: "optional-string",
                short_name: "S",
                doc: "optional scalar",
                optional: true,
                ..Annotation::default()
            },
            "optionalString",
            ValueClass::Text,
            false,
            false,
            Value::Null,
        ),
        // Declared without `optional = true`, initialised, therefore optional.
        Definition::new(
            Annotation {
                full_name: "defaulted-required",
                doc: "declared required, initialised, therefore optional",
                ..Annotation::default()
            },
            "defaultedRequired",
            ValueClass::Text,
            false,
            false,
            Value::Str("preset".to_string()),
        ),
        Definition::new(
            Annotation {
                full_name: "flag",
                doc: "a boolean, which may appear with no value",
                optional: true,
                ..Annotation::default()
            },
            "flag",
            ValueClass::Boolean,
            false,
            true,
            Value::Bool(false),
        ),
        Definition::new(
            Annotation {
                full_name: "bounded-int",
                doc: "a hard range",
                optional: true,
                min_value: 1.0,
                max_value: 10.0,
                ..Annotation::default()
            },
            "boundedInt",
            ValueClass::Integer,
            false,
            false,
            Value::Null,
        ),
        // A primitive field: boxed to `Integer` for the range check, and primitive for the null
        // check, which is the only place the difference shows.
        Definition::new(
            Annotation {
                full_name: "primitive-int",
                doc: "a primitive field, so not nullable",
                optional: true,
                ..Annotation::default()
            },
            "primitiveInt",
            ValueClass::Integer,
            false,
            true,
            Value::Int(7),
        ),
        Definition::new(
            Annotation {
                full_name: "recommended-int",
                doc: "a recommended range and no hard range",
                optional: true,
                min_recommended_value: 5.0,
                max_recommended_value: 8.0,
                ..Annotation::default()
            },
            "recommendedInt",
            ValueClass::Integer,
            false,
            false,
            Value::Null,
        ),
        Definition::new(
            Annotation {
                full_name: "bounded-double",
                doc: "a minimum only, so the message names one bound",
                optional: true,
                min_value: 0.0,
                ..Annotation::default()
            },
            "boundedDouble",
            ValueClass::Double,
            false,
            false,
            Value::Null,
        ),
        // Initialised to a **non-empty** collection: optional, and the contents are discarded as
        // soon as the argument is named.
        Definition::new(
            Annotation {
                full_name: "collection",
                doc: "an optional collection",
                optional: true,
                ..Annotation::default()
            },
            "collection",
            ValueClass::Text,
            true,
            false,
            Value::List(vec![Value::Str("declared".to_string())]),
        ),
        Definition::new(
            Annotation {
                full_name: "mutex-a",
                doc: "mutex with b",
                optional: true,
                mutex: &["mutex-b"],
                ..Annotation::default()
            },
            "mutexA",
            ValueClass::Text,
            false,
            false,
            Value::Null,
        ),
        Definition::new(
            Annotation {
                full_name: "mutex-b",
                doc: "mutex with a",
                optional: true,
                mutex: &["mutex-a"],
                ..Annotation::default()
            },
            "mutexB",
            ValueClass::Text,
            false,
            false,
            Value::Null,
        ),
        Definition::new(
            Annotation {
                full_name: "enum-arg",
                doc: "an enum, built by Enum.valueOf",
                optional: true,
                ..Annotation::default()
            },
            "enumArg",
            MODE,
            false,
            false,
            Value::Null,
        ),
    ]
}

/// The long names the dump prints a `field` row for, in its order.
const REPORTED: [&str; 13] = [
    "required-string",
    "required-collection",
    "optional-string",
    "defaulted-required",
    "flag",
    "bounded-int",
    "primitive-int",
    "recommended-int",
    "bounded-double",
    "collection",
    "mutex-a",
    "mutex-b",
    "enum-arg",
];

/// `E:<class>:<message>`, the shape the dump reports a rejection in.
fn render_error(error: &Error) -> String {
    format!("E:{}:{}", error.class, error.message)
}

#[test]
fn every_vector_is_read_the_way_the_reference_reads_it() {
    let text = golden();

    // The cases come out of the golden's own `case` rows, so the port is run over exactly the
    // command lines the reference was run over rather than a transcription of them. The expected
    // *answers* are the `result` and `field` rows, which is what is being compared.
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
        // The dump joins argv with a single space, and no fixture value here contains one. The
        // `scalar-empty-value` case ends in an empty argument, which `split` keeps only because the
        // trailing tab is preserved by the row above.
        let argv: Vec<&str> = if argv_text.is_empty() {
            Vec::new()
        } else {
            argv_text.split(' ').collect()
        };

        produced.push(format!("case\t{label}\t{argv_text}"));

        let mut parser = Parser::new(definitions());
        match parser.parse_arguments(&argv) {
            Ok(()) => {
                produced.push(format!("result\t{label}\tok"));
                for name in REPORTED {
                    let value = parser.value_of(name).expect("a declared field");
                    produced.push(format!(
                        "field\t{label}\t{name}\t{}",
                        value.to_java_string()
                    ));
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

/// `isOptional` is the annotation **or** an initialised default, and an empty collection is not an
/// initialised default.
#[test]
fn an_initialised_field_is_optional_whatever_the_annotation_says() {
    let definitions = definitions();
    let by_name = |name: &str| {
        definitions
            .iter()
            .find(|definition| definition.long_name() == name)
            .expect("a declared field")
    };

    assert!(!by_name("required-string").is_optional());
    // Declared exactly like `required-string` bar the initialiser.
    assert!(by_name("defaulted-required").is_optional());
    // Initialised, but to an empty collection, which renders as "null".
    assert!(!by_name("required-collection").is_optional());
    assert!(by_name("collection").is_optional());
    assert_eq!(
        by_name("required-collection").default_value_as_string(),
        "null"
    );
    assert_eq!(
        by_name("collection").default_value_as_string(),
        "[declared]"
    );
}

/// The recommended range is checked against the hard bounds, so a value far outside it is not out
/// of range at all.
#[test]
fn the_recommended_range_is_checked_against_the_hard_one() {
    let definitions = definitions();
    let recommended = definitions
        .iter()
        .find(|definition| definition.long_name() == "recommended-int")
        .expect("a declared field");

    assert!(recommended.has_recommended_range());
    assert!(!recommended.has_bounded_range());
    // 100 is far outside [5, 8] and inside the hard range, which is unbounded.
    assert_eq!(
        recommended.check_argument_range(&Value::Int(100)),
        Ok(false)
    );
    // The only value the warning can fire for.
    assert_eq!(recommended.check_argument_range(&Value::Null), Ok(true));
}
