//! Conformance for the enum-valued arguments against GATK 4.6.2.0.
//!
//! Golden from `tools/argument-conformance/ToolArgumentEnumDump.java`, which asked each ported tool
//! for its parser, took every argument whose underlying field class is an enum, and printed the
//! type's constants in declaration order.
//!
//! # What this suite is for
//!
//!  * **the constants being in declaration order and not a sorted one**;
//!  * **the conversion being `Enum.valueOf` and therefore case sensitive**;
//!  * **the refusal listing every constant, which is why the list and not the count is carried**;
//!  * **and a `ClpEnum`'s per-constant documentation belonging to the usage and not the refusal.**

use gatk_corpus as corpus;
use gatk_tools::tool_declarations::{
    declarations, enum_type, ENUM_TYPES, {self as decl},
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/tool_argument_enums.txt.gz"),
    )
}

fn rows(text: &str, kind: &str) -> Vec<String> {
    let prefix = format!("{kind}\t");
    text.lines()
        .filter(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
        .collect()
}

/// Every type the golden holds is in the table, with its constants in the golden's order.
#[test]
fn the_constants_are_in_declaration_order() {
    let text = golden();
    let listed = rows(&text, "enum");
    assert_eq!(listed.len(), ENUM_TYPES.len());
    for row in listed {
        let (name, constants) = row.split_once('\t').expect("a name and its constants");
        let type_ = enum_type(name).unwrap_or_else(|| panic!("{name}"));
        let written: Vec<&str> = constants.split(',').collect();
        assert_eq!(type_.constants, written.as_slice(), "{name}");
        // Declaration order, which is not the sorted one for most of them.
        let mut sorted = written.clone();
        sorted.sort_unstable();
        if name == "LogLevel" || name == "Mode" || name == "ValidationStringency" {
            assert_ne!(type_.constants, sorted.as_slice(), "{name}");
        }
    }
}

/// Every enum-valued argument points at a type the table carries, and its default is a constant.
#[test]
fn an_arguments_default_is_one_of_the_constants() {
    let text = golden();
    for row in rows(&text, "arg") {
        let mut parts = row.split('\t');
        let tool = parts.next().expect("a tool");
        let name = parts.next().expect("a long name");
        let body = parts.next().expect("a type and a default");
        let (type_name, default) = body.split_once('|').expect("a type and a default");
        let list = declarations(tool).unwrap_or_else(|| panic!("{tool}"));
        let declaration = list
            .iter()
            .find(|declaration| declaration.long_name == name)
            .unwrap_or_else(|| panic!("{tool}/{name}"));
        assert_eq!(declaration.type_name, type_name, "{tool}/{name}");
        let type_ = enum_type(type_name).unwrap_or_else(|| panic!("{type_name}"));
        // An unset enum argument has no default at all; a set one holds a constant.
        if default == "null" {
            assert_eq!(declaration.default, None, "{tool}/{name}");
        } else {
            assert!(type_.constants.contains(&default), "{tool}/{name}");
            assert_eq!(declaration.default, Some(default), "{tool}/{name}");
        }
    }
}

/// A `ClpEnum` documents its constants, and that documentation is the usage text's.
#[test]
fn a_clp_enum_documents_its_constants() {
    let text = golden();
    let documented: Vec<String> = rows(&text, "clp");
    // TWO types in this corpus implement it now: `Mode` was alone until `SplitIntervals` brought
    // `IntervalListScatterMode`, whose five constants each carry a sentence of their own. Every
    // constant of a documented type is documented, and the golden's rows are exactly those.
    let implementing: Vec<&str> = ENUM_TYPES
        .iter()
        .filter(|type_| !type_.docs.is_empty())
        .map(|type_| type_.name)
        .collect();
    assert_eq!(implementing, vec!["IntervalListScatterMode", "Mode"]);
    let expected: usize = ENUM_TYPES
        .iter()
        .filter(|type_| !type_.docs.is_empty())
        .map(|type_| {
            assert_eq!(type_.docs.len(), type_.constants.len(), "{}", type_.name);
            type_.docs.len()
        })
        .sum();
    assert_eq!(documented.len(), expected);
    for row in documented {
        let (name, body) = row.split_once('\t').expect("a type and a body");
        let (constant, doc) = body.split_once('=').expect("a constant and its doc");
        let type_ = enum_type(name).unwrap_or_else(|| panic!("{name}"));
        let written = type_
            .docs
            .iter()
            .find(|(written, _)| *written == constant)
            .unwrap_or_else(|| panic!("{name}/{constant}"));
        assert_eq!(written.1, doc, "{name}/{constant}");
    }
    // And every type that does NOT implement it carries no documentation at all, which is what
    // keeps the two kinds apart.
    for type_ in ENUM_TYPES {
        if !implementing.contains(&type_.name) {
            assert!(type_.docs.is_empty(), "{}", type_.name);
        }
    }
}

/// The refusal a value outside the type produces lists every constant, in that order.
#[test]
fn the_refusal_lists_every_constant() {
    let text = golden();
    let outcome = |case: &str| {
        rows(&text, "parse")
            .into_iter()
            .find(|row| row.starts_with(&format!("CountReads\t{case}\t")))
            .map(|row| row.rsplit_once('\t').expect("an outcome").1.to_string())
            .unwrap_or_else(|| panic!("{case}"))
    };
    assert_eq!(outcome("a-constant"), "ok");
    assert_eq!(outcome("the-other-constant"), "ok");
    // The conversion is `Enum.valueOf`, so the lower-case spelling of a constant is not one.
    let lower = outcome("lower-case");
    assert!(lower.contains("'union' is not a valid value"), "{lower}");
    // And the message lists the type's constants, in the table's own order.
    let rule = enum_type("IntervalSetRule").expect("IntervalSetRule");
    for constant in rule.constants {
        assert!(lower.contains(constant), "{lower}");
    }
    let position = |needle: &str| lower.find(needle).expect("a constant");
    assert!(position(rule.constants[0]) < position(rule.constants[1]));
    // An empty value is refused the same way rather than falling back to the default.
    assert!(outcome("an-empty-value").contains("is not a valid value"));
    // A second type, so the shape is not one type's own.
    assert_eq!(outcome("a-stringency"), "ok");
    let stringency = outcome("not-a-stringency");
    for constant in enum_type("ValidationStringency").expect("it").constants {
        assert!(stringency.contains(constant), "{stringency}");
    }
    // The table is reachable from a declaration, which is what a parser needs.
    let declaration = decl::COUNTREADS
        .iter()
        .find(|declaration| declaration.long_name == "interval-set-rule")
        .expect("the set rule");
    assert!(enum_type(declaration.type_name).is_some());
}
