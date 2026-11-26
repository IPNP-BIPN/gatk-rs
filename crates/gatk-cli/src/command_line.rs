//! `CommandLineArgumentParser.getCommandLine()`: the string a tool records in every file it writes.
//!
//! A BAM's `@PG` carries it in `CL` and a VCF's `##GATKCommandLine` carries it in `CommandLine`,
//! so it is a byte of the output for every record-transform and variant-transform tool there is.
//! `sam_output::Options` takes it as an INPUT because no port could invent it, and that is what
//! kept `PrintReads` from having a runner.
//!
//! # Two groups, and neither is the order the user typed
//!
//! The arguments the user SET come first, in the parser's own declaration order, and then the ones
//! that were not set but have a non-null default, in that same order. So `-I` then `-O` on the
//! command line comes back as `--output` then `--input` when the tool declares the output first.
//!
//! # Everything is the long form
//!
//! A short alias is expanded and a value is separated from its name by a space. (`--name=value` is
//! not a form the parser accepts at all: it refuses with `Can't parse option name containing an
//! embedded '='`.)
//!
//! # A collection is one pair per element
//!
//! Two `-I` become two `--input` pairs rather than one with a list after it.
//!
//! # A plugin's own argument is not a default
//!
//! The read-filter descriptor registers every filter it discovers, so an unselected one's
//! arguments sit in the parser with their defaults. None of them reaches the line: it stops at
//! `--disable-tool-default-read-filters`, which is the descriptor's own argument rather than a
//! plugin's. One that was SET is printed, because setting it is what selects the plugin.
//!
//! # A default of null is omitted, and an EMPTY default is not
//!
//! The filter is on the string `null`, so a collection whose default renders that way disappears
//! while `--gcs-project-for-requester-pays` is printed with nothing after it.
//!
//! Ported from `org.broadinstitute.barclay.argparser.CommandLineArgumentParser.getCommandLine` and
//! `NamedArgumentDefinition.getCommandLineDisplayString`.

use gatk_barclay::{Parser, Value};

/// `NULL_ARGUMENT_STRING`, which is what a default has to differ from to be printed.
pub const NULL_ARGUMENT_STRING: &str = "null";

/// The values one argument contributes, each of which becomes its own `--name value` pair.
fn display_values(value: &Value) -> Vec<String> {
    match value {
        Value::Null => Vec::new(),
        Value::List(values) => values.iter().map(display_one).collect(),
        other => vec![display_one(other)],
    }
}

fn display_one(value: &Value) -> String {
    match value {
        // A tagged value prints its own string; the tag itself belongs to the NAME, and none of
        // the tools that reach this carry one on a set argument.
        Value::Tagged { value, .. } => value.clone(),
        other => other.to_java_string(),
    }
}

/// What one definition's field prints: a `float` is held widened and prints as the `float` it is,
/// so `--min-category-reads-percentage` says 0.05 where the widened double says 0.05000000074505806.
fn definition_values(definition: &gatk_barclay::Definition) -> Vec<String> {
    match (&definition.class, &definition.value) {
        (gatk_barclay::ValueClass::Float, Value::Double(value)) => {
            vec![jmath::to_string::float_to_string(*value as f32)]
        }
        // `Set<MetricAccumulationLevel>` starts as an `EnumSet` and stays one: the levels print in
        // declaration order, once each, whatever order they were asked in.
        (gatk_barclay::ValueClass::Enum { constants, .. }, Value::List(_))
            if definition.long_name() == "metric-accumulation-level" =>
        {
            let asked = display_values(&definition.value);
            constants
                .iter()
                .filter(|constant| asked.iter().any(|level| level == *constant))
                .map(|constant| constant.to_string())
                .collect()
        }
        (_, value) => display_values(value),
    }
}

/// Whether any value of an argument was given a tag on the command line (`--mask:BED file`).
fn has_tag(value: &Value) -> bool {
    match value {
        Value::Tagged { tag, .. } => tag.is_some(),
        Value::List(values) => values.iter().any(has_tag),
        _ => false,
    }
}

/// Each value of a tagged argument as Barclay prints it: the tag joins the NAME, and a
/// `FeatureInput` the user named prints `tag:absolutePath` (`FeatureInput.toString`).
///
/// shortcut: a tag's `key=value` attributes are not printed, as no measured line carries one;
/// upgrade when one does.
fn tagged_pairs(definition: &gatk_barclay::Definition) -> Vec<(String, String)> {
    let values = match &definition.value {
        Value::List(values) => values.clone(),
        other => vec![other.clone()],
    };
    let feature_input = matches!(
        definition.class,
        gatk_barclay::ValueClass::Constructed {
            simple_name: "FeatureInput",
            ..
        }
    );
    values
        .iter()
        .map(|value| match value {
            Value::Tagged {
                value,
                tag: Some(tag),
                ..
            } => {
                let shown = if feature_input {
                    let absolute = if value.contains("://") {
                        value.clone()
                    } else {
                        std::path::absolute(value)
                            .map(|path| path.display().to_string())
                            .unwrap_or_else(|_| value.clone())
                    };
                    format!("{tag}:{absolute}")
                } else {
                    value.clone()
                };
                (format!("{}:{tag}", definition.long_name()), shown)
            }
            other => (definition.long_name().to_string(), display_one(other)),
        })
        .collect()
}

/// `getCommandLine()`: the class's simple name, then the set arguments, then the defaulted ones.
///
/// `class_name` is `callerArguments.getClass().getSimpleName()`, which for every tool here is the
/// tool's own name.
pub fn expanded(class_name: &str, parser: &Parser) -> String {
    let mut line = String::from(class_name);

    let push = |line: &mut String, name: &str, values: Vec<String>| {
        for value in values {
            line.push_str(" --");
            line.push_str(name);
            line.push(' ');
            line.push_str(&value);
        }
    };

    for definition in parser.definitions() {
        if definition.has_been_set() && has_tag(&definition.value) {
            for (name, value) in tagged_pairs(definition) {
                push(&mut line, &name, vec![value]);
            }
            continue;
        }
        if definition.has_been_set() {
            let values = definition_values(definition);
            if values.is_empty() {
                // A collection SET to `null` is emptied and still counts as set: its display string
                // is empty and the join still puts a space for it, so the line carries two spaces
                // where the argument would be. Measured on `ComposeSTRTableFile --intervals null`.
                line.push(' ');
            } else {
                push(&mut line, definition.long_name(), values);
            }
        }
    }
    for definition in parser.definitions() {
        if definition.has_been_set() {
            continue;
        }
        if definition.default_value_as_string() == NULL_ARGUMENT_STRING {
            continue;
        }
        // A PLUGIN's own argument reaches the line when its plugin is SELECTED, and a tool's
        // DEFAULT filter counts as selected. The descriptor registers every read filter it
        // discovers, so an unselected one's arguments sit in the parser with their defaults and
        // reach no command line the reference records; `PrintDistantMates` is the first tool here
        // whose default filter has an argument of its own, and the reference prints
        // `--mate-too-distant-length 1000` for it.
        if let Some(control) = &definition.controlled_by {
            if !parser.plugin_is_selected(control) {
                continue;
            }
        }
        // The DEFAULT is what is printed, and the parser left it in the value, so the two agree.
        // A default that renders as the empty string still prints its name and a trailing space.
        let values = definition_values(definition);
        if values.is_empty() {
            line.push_str(" --");
            line.push_str(definition.long_name());
            line.push(' ');
        } else {
            push(&mut line, definition.long_name(), values);
        }
    }
    line
}
