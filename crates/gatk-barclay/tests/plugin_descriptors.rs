//! Conformance for Barclay's plugin descriptors, against the oracle.
//!
//! Golden from `tools/argument-conformance/BarclayPluginDescriptorDump.java`, which runs GATK's
//! own `GATKReadFilterPluginDescriptor` over GATK's own read filters. `--read-filter` **is** this
//! mechanism, so a stand-in descriptor would have measured the stand-in.
//!
//! # What the golden settles
//!
//! ```text
//! defs  discovered  1  ambig-filter-frac        AmbiguousBaseReadFilter
//! defs  discovered  2  minimum-mapping-quality  MappingQualityReadFilter
//! result  no-filter                 ok
//! result  dependent-without-filter  ... Argument "/ambig-filter-frac" is only valid when the
//!                                       argument "AmbiguousBaseReadFilter" is specified
//! result  unknown-filter-name       ... Unrecognized read filter name: NoSuchReadFilter
//! ```
//!
//! Two arguments of two read filters are in the parser on **every** run, whether or not anybody
//! asked for those filters. `validatePluginArgumentValues` runs before the required check and
//! removes each one that nobody set and whose filter nobody named — so they are not merely
//! optional, they are absent, and a required argument of an unselected filter would not fire at
//! all.
//!
//! Given without its filter, the same argument is an error naming the filter **class**, not
//! `--read-filter`, and building the argument's own name from `getShortName() + "/" +
//! getLongName()` with no guard, which is where the leading slash comes from.
//!
//! The last row is not Barclay's: `validateAndResolvePlugins` is the descriptor's own hook, and
//! GATK uses it to refuse a name that matches no filter.

use gatk_barclay::{
    Annotation, Definition, Error, Parser, PluginControl, PluginResolution, Value, ValueClass,
};
use gatk_corpus as corpus;

/// `ReadFilterArgumentDefinitions.READ_FILTER_LONG_NAME`.
const SELECTOR: &str = "read-filter";

/// The read filters this suite names. Not the whole library: the descriptor registers every
/// filter's arguments at once, and that count is a property of the library rather than of this
/// mechanism.
const KNOWN: [&str; 2] = ["AmbiguousBaseReadFilter", "MappingQualityReadFilter"];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/plugin_descriptors.txt.gz"),
    )
}

/// The four arguments the dump reports on, in its order.
fn definitions() -> Vec<Definition> {
    let mut selector = Definition::new(
        Annotation {
            full_name: SELECTOR,
            doc: "the read filters to apply",
            optional: true,
            ..Annotation::default()
        },
        "userReadFilterNames",
        ValueClass::Text,
        true,
        false,
        Value::List(Vec::new()),
    );
    selector.controlled_by = None;

    // `AmbiguousBaseReadFilter.maxAmbiguousBaseFraction`, whose default is 0.05.
    let mut ambiguous = Definition::new(
        Annotation {
            full_name: "ambig-filter-frac",
            doc: "threshold fraction of ambiguous bases",
            optional: true,
            ..Annotation::default()
        },
        "maxAmbiguousBaseFraction",
        ValueClass::Double,
        false,
        false,
        Value::Double(0.05),
    );
    ambiguous.controlled_by = Some(PluginControl {
        predecessor: "AmbiguousBaseReadFilter",
        selector: SELECTOR,
    });

    // `MappingQualityReadFilter.minMappingQualityScore`, whose default is 10.
    let mut mapping_quality = Definition::new(
        Annotation {
            full_name: "minimum-mapping-quality",
            doc: "minimum mapping quality to keep",
            optional: true,
            ..Annotation::default()
        },
        "minMappingQualityScore",
        ValueClass::Integer,
        false,
        true,
        Value::Int(10),
    );
    mapping_quality.controlled_by = Some(PluginControl {
        predecessor: "MappingQualityReadFilter",
        selector: SELECTOR,
    });

    let tool = Definition::new(
        Annotation {
            full_name: "tool-arg",
            doc: "the tool's own",
            optional: true,
            ..Annotation::default()
        },
        "toolArg",
        ValueClass::Text,
        false,
        false,
        Value::Null,
    );

    vec![selector, ambiguous, mapping_quality, tool]
}

fn parser() -> Parser {
    Parser::new(definitions()).with_plugin_resolution(PluginResolution {
        selector: SELECTOR,
        known: KNOWN.iter().map(|name| name.to_string()).collect(),
        // `GATKReadFilterPluginDescriptor.validateAndResolvePlugins`.
        unrecognized_prefix: "Unrecognized read filter name: ",
    })
}

/// The dump's `defs` rows carry which object declared each argument.
fn owner(definition: &Definition) -> String {
    match &definition.controlled_by {
        Some(control) => control.predecessor.to_string(),
        None => "tool".to_string(),
    }
}

fn render_error(error: &Error) -> String {
    format!("E:{}:{}", error.class, error.message.replace('\n', "\\n"))
}

#[test]
fn every_plugin_decision_is_the_one_the_reference_makes() {
    let text = golden();

    let mut produced: Vec<String> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        expected.push(line);

        if let Some(rest) = line.strip_prefix("defs\t") {
            let mut parts = rest.splitn(2, '\t');
            let label = parts.next().expect("a label");
            let index: usize = parts
                .next()
                .and_then(|rest| rest.split('\t').next().map(str::to_string))
                .expect("an index")
                .parse()
                .expect("a numeric index");
            if index == 0 {
                for (position, definition) in definitions().iter().enumerate() {
                    produced.push(format!(
                        "defs\t{label}\t{position}\t{}\t{}",
                        definition.long_name(),
                        owner(definition)
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

        let mut parser = parser();
        match parser.parse_arguments(&argv) {
            Ok(()) => {
                produced.push(format!("result\t{label}\tok"));
                // Every reported argument, including the controlled ones the trim dropped: the
                // dump reads the definition's value, which the trim does not touch.
                for definition in parser.definitions() {
                    produced.push(format!(
                        "field\t{label}\t{}\t{}",
                        definition.long_name(),
                        definition.value.to_java_string()
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

/// A controlled argument is *removed* rather than made optional, which is only visible from the
/// fact that setting it without its filter is an error while leaving it alone is not.
#[test]
fn a_dangling_dependent_argument_names_the_filter_class() {
    let mut untouched = parser();
    untouched
        .parse_arguments(&[])
        .expect("an unselected filter's argument is absent, not merely unset");

    let mut dangling = parser();
    let error = dangling
        .parse_arguments(&["--ambig-filter-frac", "0.1"])
        .expect_err("setting it without its filter is refused");
    assert_eq!(
        error.message,
        "Argument \"/ambig-filter-frac\" is only valid when the argument \
         \"AmbiguousBaseReadFilter\" is specified"
    );

    let mut selected = parser();
    selected
        .parse_arguments(&[
            "--read-filter",
            "AmbiguousBaseReadFilter",
            "--ambig-filter-frac",
            "0.1",
        ])
        .expect("naming the filter makes its argument usable");
}

/// The unknown-name refusal is the descriptor's, not the argument layer's.
#[test]
fn an_unknown_filter_name_is_refused_by_the_descriptor() {
    let mut parser = parser();
    let error = parser
        .parse_arguments(&["--read-filter", "NoSuchReadFilter"])
        .expect_err("a name matching no filter is refused");
    assert_eq!(
        error.message,
        "Unrecognized read filter name: NoSuchReadFilter"
    );
}
