//! Conformance for the declaration-to-definition map against GATK 4.6.2.0.
//!
//! The declarations are the reference's own parser reporting itself, and the definitions are what
//! the ported parser consumes. This suite checks the map between them, and it checks that the map
//! is now total: every class converts, and every controlled argument carries the filter that
//! declared it, so a walker's command line reaches the parser.
//!
//! # What this suite is for
//!
//!  * **a definition carrying what the declaration carries**;
//!  * **the default rendering surviving the round trip**;
//!  * **an enum reaching its constants**;
//!  * **every declared class converting, so the gap is counted at zero rather than described**;
//!  * **and a controlled argument carrying its owner, which is what the trim runs over.** The trim
//!    itself is `plugin_ownership.rs`'s suite.

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

/// Every class converts, and every controlled argument reaches the trim with its owner.
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
        // Nothing is skipped for the seven the golden started with: every class they name has a
        // measured conversion. The two tools added later are checked below, one of them naming a
        // class that is not measured yet.
        assert!(skipped.is_empty(), "{tool}: {skipped:?}");
        assert!(!list.iter().any(unconvertible), "{tool}");
        assert!(UNCONVERTIBLE_CLASSES.is_empty());
        // Every tool parses now. A walker's descriptor owns arguments that read as required, and
        // the trim removes them before the required check because no filter selected them; the
        // port runs that trim over the measured ownership table, so the twelve required
        // controlled arguments no longer stop a command line.
        assert!(gatk_cli::parseable(tool), "{tool}");
        for declaration in list.iter().filter(|d| d.controlled_by.is_some()) {
            assert!(
                definition(declaration)
                    .expect("a definition")
                    .controlled_by
                    .is_some(),
                "{tool}: {}",
                declaration.long_name
            );
        }
    }
    // A class the seven did not name: `CreateHadoopBamSplittingIndex` declares a `Long`, which
    // was the last class in the nine tools' declarations without a measured conversion. It has
    // one now, so nothing is skipped and the tool parses.
    let spark = declarations("CreateHadoopBamSplittingIndex").expect("its declarations");
    assert_eq!(find(spark, "splitting-index-granularity").type_name, "Long");
    assert!(missing(spark).is_empty());
    assert!(gatk_cli::parseable("CreateHadoopBamSplittingIndex"));
    assert!(!spark
        .iter()
        .any(|declaration| declaration.controlled_by.is_some()));
    // Where the other new tool declares nothing exotic and does parse.
    assert!(gatk_cli::parseable("PrintBGZFBlockInformation"));
    assert!(
        missing(declarations("PrintBGZFBlockInformation").expect("its declarations")).is_empty()
    );

    // The smallest tool is the clearest: fourteen arguments, every one of them convertible, and no
    // plugin descriptor at all, so it is the first tool this port can hand a command line to.
    assert_eq!(INDEXFEATUREFILE.len(), 14);
    assert!(missing(INDEXFEATUREFILE).is_empty());
    assert!(gatk_cli::parseable("IndexFeatureFile"));
    assert!(!INDEXFEATUREFILE
        .iter()
        .any(|declaration| declaration.controlled_by.is_some()));
}

/// The definitions that ARE built parse a command line the way the reference's parser does.
#[test]
fn the_definitions_that_are_built_parse() {
    // A parser over CountReads' definitions accepts the command line the reference accepts. The
    // twelve plugin-controlled arguments that read as REQUIRED in the declarations are removed by
    // the trim before the required check, because no filter selected them, and the trim now knows
    // which filter declared each of them.
    let mut parser = Parser::new(definitions(COUNTREADS));
    parser
        .parse_arguments(&[
            "--input",
            "/dev/null",
            "--interval-set-rule",
            "INTERSECTION",
        ])
        .expect("the plugin arguments are trimmed");
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
    // And a whole command line through the dispatcher's own parser is accepted: `parse_failure`
    // says nothing because there is nothing to refuse.
    assert_eq!(
        gatk_cli::parse_failure(
            "CountReads",
            &["--input".to_string(), "/dev/null".to_string()]
        ),
        None
    );
}
