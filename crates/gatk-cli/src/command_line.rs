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
        if definition.has_been_set() {
            push(
                &mut line,
                definition.long_name(),
                display_values(&definition.value),
            );
        }
    }
    for definition in parser.definitions() {
        if definition.has_been_set() {
            continue;
        }
        if definition.default_value_as_string() == NULL_ARGUMENT_STRING {
            continue;
        }
        // A PLUGIN's own argument does not reach the line unless somebody set it. The descriptor
        // registers every read filter it discovers, so an unselected one's arguments are in the
        // parser with their defaults and in no command line the reference records: the line stops
        // at `--disable-tool-default-read-filters`, which is the descriptor's own.
        if definition.controlled_by.is_some() {
            continue;
        }
        // The DEFAULT is what is printed, and the parser left it in the value, so the two agree.
        // A default that renders as the empty string still prints its name and a trailing space.
        let values = display_values(&definition.value);
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
