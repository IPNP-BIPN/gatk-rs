//! Conformance for Barclay's `@ArgumentCollection` flattening, against the oracle.
//!
//! Golden from `tools/argument-conformance/BarclayArgumentCollectionDump.java`.
//!
//! # What the golden settles
//!
//! ```text
//! defs  derived  0  derived-required     the subclass's own field, first
//! defs  derived  1  middle-before        the nested collection, at the point the field appears
//! defs  derived  2  inner-one            and its own nesting, depth-first
//! defs  derived  5  derived-last
//! defs  derived  6  base-required        the superclass, last
//! result  nothing-given  ... Argument 'derived-required' is required
//! result  clashing-aliases  ... clash has already been used.
//! ```
//!
//! The order is the whole finding. `getAllFields` adds a class's own declared fields and *then*
//! climbs to its superclass, so a base class's required argument is checked after the subclass's,
//! and a user who omits both is told about the subclass's one. A nested collection is spliced in
//! where its field sits, not appended.

use gatk_barclay::{
    create_argument_definitions, Annotation, ClassDecl, Definition, Error, FieldDecl, Parser,
    Value, ValueClass,
};
use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/argument_collections.txt.gz"),
    )
}

/// One `@Argument String` field, which is every field this dump declares.
fn text(full_name: &'static str, field_name: &'static str, optional: bool) -> FieldDecl {
    FieldDecl::Argument(Box::new(Definition::new(
        Annotation {
            full_name,
            optional,
            doc: "",
            ..Annotation::default()
        },
        field_name,
        ValueClass::Text,
        false,
        false,
        Value::Null,
    )))
}

/// `Inner`, the innermost collection.
fn inner() -> ClassDecl {
    ClassDecl::new(
        "BarclayArgumentCollectionDump$Inner",
        vec![
            text("inner-one", "innerOne", true),
            text("inner-two", "innerTwo", true),
        ],
    )
}

/// `Middle`, which holds `Inner` **between** two of its own arguments.
fn middle() -> ClassDecl {
    ClassDecl::new(
        "BarclayArgumentCollectionDump$Middle",
        vec![
            text("middle-before", "middleBefore", true),
            FieldDecl::Collection(inner()),
            text("middle-after", "middleAfter", true),
        ],
    )
}

/// `Derived extends Base`.
fn derived() -> ClassDecl {
    ClassDecl::new(
        "BarclayArgumentCollectionDump$Derived",
        vec![
            text("derived-required", "derivedRequired", false),
            FieldDecl::Collection(middle()),
            text("derived-last", "derivedLast", true),
        ],
    )
    .extending(ClassDecl::new(
        "BarclayArgumentCollectionDump$Base",
        vec![
            text("base-required", "baseRequired", false),
            text("base-optional", "baseOptional", true),
        ],
    ))
}

fn clashing() -> ClassDecl {
    ClassDecl::new(
        "BarclayArgumentCollectionDump$Clashing",
        vec![
            FieldDecl::Collection(ClassDecl::new(
                "BarclayArgumentCollectionDump$ClashingA",
                vec![text("clash", "clash", true)],
            )),
            FieldDecl::Collection(ClassDecl::new(
                "BarclayArgumentCollectionDump$ClashingB",
                vec![text("clash", "clash", true)],
            )),
        ],
    )
}

fn uninitialised() -> ClassDecl {
    ClassDecl::new(
        "BarclayArgumentCollectionDump$Uninitialised",
        vec![FieldDecl::UninitialisedCollection {
            field_name: "inner",
        }],
    )
}

fn both_annotations() -> ClassDecl {
    ClassDecl::new(
        "BarclayArgumentCollectionDump$BothAnnotations",
        vec![FieldDecl::BothAnnotations { field_name: "both" }],
    )
}

/// An outer object and a nested one declaring the same name: a construction failure, not shadowing.
fn shadowing() -> ClassDecl {
    ClassDecl::new(
        "BarclayArgumentCollectionDump$Shadowing",
        vec![
            text("derived-last", "derivedLast", true),
            FieldDecl::Collection(ClassDecl::new(
                "BarclayArgumentCollectionDump$ShadowingInner",
                vec![text("derived-last", "derivedLast", true)],
            )),
        ],
    )
}

fn render_error(error: &Error) -> String {
    format!("E:{}:{}", error.class, error.message.replace('\n', "\\n"))
}

/// The declaration each label was run over.
fn class_for(label: &str) -> ClassDecl {
    match label {
        "clashing-aliases" => clashing(),
        "uninitialised-collection" => uninitialised(),
        "both-annotations" => both_annotations(),
        "shadowing" => shadowing(),
        _ => derived(),
    }
}

#[test]
fn the_flattening_is_the_one_the_reference_builds() {
    let text = golden();

    let mut produced: Vec<String> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        expected.push(line);

        if let Some(rest) = line.strip_prefix("defs\t") {
            // One `defs` block per label, emitted once when its first row is seen.
            let mut parts = rest.splitn(3, '\t');
            let label = parts.next().expect("a label");
            let index: usize = parts
                .next()
                .expect("an index")
                .parse()
                .expect("a numeric index");
            if index == 0 {
                let definitions =
                    create_argument_definitions(&class_for(label)).expect("the class builds");
                for (position, definition) in definitions.iter().enumerate() {
                    produced.push(format!(
                        "defs\t{label}\t{position}\t{}",
                        definition.long_name()
                    ));
                }
            }
            continue;
        }

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

        // The construction failures happen here, before any command line is looked at, which is
        // why `Parser::from_class` returns a `Result` at all.
        match Parser::from_class(&class_for(label)) {
            Err(error) => produced.push(format!("result\t{label}\t{}", render_error(&error))),
            Ok(mut parser) => match parser.parse_arguments(&argv) {
                Ok(()) => {
                    produced.push(format!("result\t{label}\tok"));
                    // The dump reports every definition, in the flattened order.
                    let rows: Vec<String> = parser
                        .definitions()
                        .iter()
                        .map(|definition| {
                            format!(
                                "field\t{label}\t{}\t{}",
                                definition.long_name(),
                                definition.value.to_java_string()
                            )
                        })
                        .collect();
                    produced.extend(rows);
                }
                Err(error) => produced.push(format!("result\t{label}\t{}", render_error(&error))),
            },
        }
    }

    assert_eq!(produced.len(), expected.len(), "row count");
    for (index, (produced, oracle)) in produced.iter().zip(expected.iter()).enumerate() {
        assert_eq!(produced, oracle, "row {index}");
    }
}

/// The order, stated as the property rather than as a list of rows.
#[test]
fn a_subclass_is_registered_before_the_class_it_extends() {
    let definitions = create_argument_definitions(&derived()).expect("the class builds");
    let names: Vec<&str> = definitions
        .iter()
        .map(|definition| definition.long_name())
        .collect();

    // The subclass's own field first, the superclass's last.
    assert_eq!(names.first(), Some(&"derived-required"));
    assert_eq!(names.last(), Some(&"base-optional"));

    // The nested collection is spliced in where its field sits, not appended: `middle-before`
    // through `middle-after` sit between the two fields the subclass declares around it.
    assert_eq!(
        names,
        vec![
            "derived-required",
            "middle-before",
            "inner-one",
            "inner-two",
            "middle-after",
            "derived-last",
            "base-required",
            "base-optional",
        ]
    );
}
