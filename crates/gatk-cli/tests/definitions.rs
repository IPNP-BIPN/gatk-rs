//! Conformance for the declaration-to-definition map against GATK 4.6.2.0.
//!
//! The declarations are the reference's own parser reporting itself, and the definitions are what
//! the ported parser consumes. This suite checks the map between them, and it checks the SIZE of
//! the gap: four classes the reference names have no measured conversion, and every argument of
//! those classes is a definition this port declines to build.
//!
//! # What this suite is for
//!
//!  * **a definition carrying what the declaration carries**;
//!  * **the default rendering surviving the round trip**;
//!  * **an enum reaching its constants**;
//!  * **and the gap being named rather than papered over: no tool parses yet, and the reason is
//!    one class each of them declares.**

use gatk_barclay::{Parser, ValueClass};
use gatk_cli::definitions::{
    definition, definitions, missing, unconvertible, UNCONVERTIBLE_CLASSES,
};
use gatk_tools::tool_declarations::{declarations, Declaration, COUNTREADS, INDEXFEATUREFILE};

fn find<'a>(list: &'a [Declaration], name: &str) -> &'a Declaration {
    list.iter()
        .find(|declaration| declaration.long_name == name)
        .unwrap_or_else(|| panic!("{name}"))
}

/// A definition carries the declaration's own answers, not a second reading of them.
#[test]
fn a_definition_carries_the_declaration() {
    let input = find(COUNTREADS, "interval-set-rule");
    let built = definition(input).expect("a definition");
    assert_eq!(built.long_name(), "interval-set-rule");
    // `getArgumentAliases()` is the short name first and then the long one, which is what an
    // error message names the argument by.
    assert_eq!(built.alias_display_string(), "isr/interval-set-rule");
    // The default survives the round trip through the value, which is what decides optionality.
    assert_eq!(built.default_value_as_string(), "UNION");
    assert!(built.is_optional());
    // An enum reaches its constants, which is what the conversion and the refusal both need.
    assert!(matches!(
        built.class,
        ValueClass::Enum {
            simple_name: "IntervalSetRule",
            ..
        }
    ));
    // A required argument stays required, and a collection stays a collection.
    let filters = find(COUNTREADS, "read-filter");
    let built = definition(filters).expect("a definition");
    assert!(built.is_collection);
    // An initialised-but-empty collection renders as null, which is the reference's own collapse.
    assert_eq!(built.default_value_as_string(), "null");
    // A flag is a primitive boolean whose default is false.
    let flag = definition(find(COUNTREADS, "disable-tool-default-read-filters")).expect("a flag");
    assert!(flag.is_flag());
    assert!(flag.field_is_primitive);
    assert_eq!(flag.default_value_as_string(), "false");
}

/// The gap is four classes, and it is counted against the declarations rather than written down.
#[test]
fn the_gap_is_named_and_counted() {
    for tool in [
        "CountReads",
        "CountVariants",
        "PrintReads",
        "ApplyBQSR",
        "SelectVariants",
        "IndexFeatureFile",
        "GatherVcfsCloud",
    ] {
        let list = declarations(tool).unwrap_or_else(|| panic!("{tool}"));
        let built = definitions(list);
        let skipped = missing(list);
        assert_eq!(built.len() + skipped.len(), list.len(), "{tool}");
        // Every skipped argument is skipped for its class and for no other reason.
        for name in &skipped {
            let declaration = find(list, name);
            assert!(unconvertible(declaration), "{tool}/{name}");
        }
        // And every one of them is one of the four classes, which is what makes the gap a list of
        // measurements to take rather than an open question.
        for declaration in list.iter().filter(|d| unconvertible(d)) {
            assert!(
                UNCONVERTIBLE_CLASSES.contains(&declaration.type_name),
                "{tool}"
            );
            assert!(definition(declaration).is_none(), "{tool}");
        }
        // No tool is parseable yet, and the reason is the same for all seven: each declares a
        // path, and the conversion a path goes through is not measured.
        assert!(!gatk_cli::parseable(tool), "{tool}");
        assert!(skipped.iter().any(|name| {
            let declaration = find(list, name);
            declaration.type_name == "GATKPath"
        }));
    }
    // The smallest tool is the clearest: fourteen arguments, of which three are paths and one is
    // a `File`, which is a different class with a different message and not the same gap twice.
    assert_eq!(INDEXFEATUREFILE.len(), 14);
    assert_eq!(
        missing(INDEXFEATUREFILE),
        ["input", "output", "tmp-dir", "arguments_file"]
    );
}

/// The definitions that ARE built parse a command line the way the reference's parser does.
#[test]
fn the_definitions_that_are_built_parse() {
    // A parser over the convertible half of CountReads answers about those arguments. It does not
    // accept a command line, and the reason is worth stating: the plugin-controlled arguments the
    // read-filter descriptor owns read as REQUIRED in the declarations, and the reference trims
    // them before the required check runs because no filter selected them. That trim is the
    // descriptor's, and it is not wired here, so the parser asks for an argument the reference
    // never asks for.
    let mut parser = Parser::new(definitions(COUNTREADS));
    let untrimmed = parser
        .parse_arguments(&["--interval-set-rule", "INTERSECTION"])
        .expect_err("the plugin arguments are still required");
    let named = COUNTREADS
        .iter()
        .find(|declaration| untrimmed.message.contains(declaration.long_name))
        .expect("the refusal names an argument");
    assert!(named.controlled_by.is_some(), "{}", untrimmed.message);
    assert!(named.required, "{}", untrimmed.message);
    let mut parser = Parser::new(definitions(COUNTREADS));
    let refusal = parser
        .parse_arguments(&["--interval-set-rule", "union"])
        .expect_err("the refusal");
    assert!(
        refusal.message.contains("'union' is not a valid value"),
        "{}",
        refusal.message
    );
    // And the refusal lists the constants, which came from the second golden.
    assert!(refusal.message.contains("UNION"), "{}", refusal.message);
    assert!(
        refusal.message.contains("INTERSECTION"),
        "{}",
        refusal.message
    );
    // The port refuses to hand a whole command line to that parser, because a third of the tool's
    // arguments are missing from it: `parse_failure` says nothing rather than something wrong.
    assert_eq!(
        gatk_cli::parse_failure(
            "CountReads",
            &["--input".to_string(), "/dev/null".to_string()]
        ),
        None
    );
}
