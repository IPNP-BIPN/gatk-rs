//! Conformance for the declaration-to-definition map against GATK 4.6.2.0.
//!
//! The declarations are the reference's own parser reporting itself, and the definitions are what
//! the ported parser consumes. This suite checks the map between them, and it checks what is left
//! of the gap: every class converts now, and what still stops a walker is the plugin trim.
//!
//! # What this suite is for
//!
//!  * **a definition carrying what the declaration carries**;
//!  * **the default rendering surviving the round trip**;
//!  * **an enum reaching its constants**;
//!  * **every declared class converting, so the gap is counted at zero rather than described**;
//!  * **and the one gap that is left being named: a walker's plugin-controlled arguments read as
//!    required, and the trim that removes them is the descriptor's.**

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

/// Every class converts, and what is left is the plugin trim.
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
        // Nothing is skipped any more, and it is counted against the declarations rather than
        // written down: every class the seven tools name has a measured conversion.
        assert!(skipped.is_empty(), "{tool}: {skipped:?}");
        assert!(!list.iter().any(unconvertible), "{tool}");
        assert!(UNCONVERTIBLE_CLASSES.is_empty());
        // What decides parseability now is the plugin trim. A walker's descriptor owns arguments
        // that read as required, and the reference removes them before the required check because
        // no filter selected them; the port has no descriptor, so it declines to parse rather
        // than asking for an argument the reference never asks for.
        let controlled = list
            .iter()
            .any(|declaration| declaration.controlled_by.is_some() && declaration.required);
        assert_eq!(gatk_cli::parseable(tool), !controlled, "{tool}");
    }
    // The smallest tool is the clearest: fourteen arguments, every one of them convertible, and no
    // plugin descriptor at all, so it is the first tool this port can hand a command line to.
    assert_eq!(INDEXFEATUREFILE.len(), 14);
    assert!(missing(INDEXFEATUREFILE).is_empty());
    assert!(gatk_cli::parseable("IndexFeatureFile"));
    assert!(!gatk_cli::parseable("CountReads"));
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
        .parse_arguments(&[
            "--input",
            "/dev/null",
            "--interval-set-rule",
            "INTERSECTION",
        ])
        .expect_err("the plugin arguments are still required");
    let named = COUNTREADS
        .iter()
        .find(|declaration| untrimmed.message.contains(declaration.long_name))
        .expect("the refusal names an argument");
    assert!(named.controlled_by.is_some(), "{}", untrimmed.message);
    assert!(named.required, "{}", untrimmed.message);
    let mut parser = Parser::new(definitions(COUNTREADS));
    let refusal = parser
        .parse_arguments(&["--input", "/dev/null", "--interval-set-rule", "union"])
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
