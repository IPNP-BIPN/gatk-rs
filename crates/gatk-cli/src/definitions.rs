//! A tool's declarations, as Barclay definitions the parser can be built from.
//!
//! [`gatk_tools::tool_declarations`] is what the reference's own parser reports: a long name, its
//! aliases, whether it is required, whether it is a collection, its default as the parser renders
//! it, its type, its bounds and its documentation. [`gatk_barclay::Definition`] is what the ported
//! parser consumes. This module is the map between them, and it is deliberately partial.
//!
//! # Which types convert
//!
//! Every class the seven tools declare converts, the four that used to be missing having been
//! measured in the `tool-argument-value-classes` golden: `Float` through Java's own grammar rather
//! than `str::parse`, and `File`, `GATKPath` and `FeatureInput` through a constructor that accepts
//! every string, which is why a bad path is not a bad value. What separates the last three is the
//! name a refusal would carry and whether a tag on the argument is a tag or an error.
//!
//! [`UNCONVERTIBLE_CLASSES`] is therefore empty. It stays as a named place for the next class the
//! declarations grow, and the test beside this file asserts the emptiness against the declarations
//! rather than against a number written down here.

use gatk_barclay::{Annotation, Definition, PluginControl, Value, ValueClass};
use gatk_tools::plugin_ownership;
use gatk_tools::tool_declarations::{enum_type, Declaration};

/// The classes the declarations name and this module does not convert, which is none of them.
pub const UNCONVERTIBLE_CLASSES: [&str; 0] = [];

/// Whether an argument's type is one of those.
pub fn unconvertible(declaration: &Declaration) -> bool {
    UNCONVERTIBLE_CLASSES.contains(&declaration.type_name)
}

/// The `ValueClass` a declared type converts through, if it is one this port has measured.
pub fn value_class(type_name: &str) -> Option<ValueClass> {
    match type_name {
        "Integer" => Some(ValueClass::Integer),
        "Double" => Some(ValueClass::Double),
        "String" => Some(ValueClass::Text),
        "Boolean" => Some(ValueClass::Boolean),
        // The three classes built from a string. `GATKPath` and `FeatureInput` implement
        // `TaggedArgument` and `File` does not, which is the difference the same command line
        // turns into a tag on one and a refusal on the other.
        "GATKPath" => Some(ValueClass::Constructed {
            simple_name: "GATKPath",
            taggable: true,
        }),
        "FeatureInput" => Some(ValueClass::Constructed {
            simple_name: "FeatureInput",
            taggable: true,
        }),
        "File" => Some(ValueClass::Constructed {
            simple_name: "File",
            taggable: false,
        }),
        "Float" => Some(ValueClass::Float),
        "Long" => Some(ValueClass::Long),
        name => enum_type(name).map(|type_| ValueClass::Enum {
            simple_name: type_.name,
            constants: type_.constants,
        }),
    }
}

/// The value a constructed instance holds, rebuilt from the rendering the golden carries.
///
/// The rendering is `String.valueOf` on the field, so it is lossy in exactly one direction: it
/// cannot say whether a collection was empty or null, both of which render as `null`. The
/// reference's own `convertDefaultValueToString` collapses them the same way, so a rebuilt empty
/// list and a rebuilt null render identically and the definition is not changed by the choice.
pub fn initial_value(declaration: &Declaration, class: &ValueClass) -> Value {
    if declaration.collection {
        return Value::List(Vec::new());
    }
    let Some(default) = declaration.default else {
        return Value::Null;
    };
    match class {
        ValueClass::Integer => default
            .parse::<i32>()
            .map(Value::Int)
            .unwrap_or(Value::Null),
        ValueClass::Double => default
            .parse::<f64>()
            .map(Value::Double)
            .unwrap_or(Value::Null),
        ValueClass::Boolean => match default {
            "true" => Value::Bool(true),
            "false" => Value::Bool(false),
            _ => Value::Null,
        },
        ValueClass::Enum { .. } => Value::Enum(default.to_string()),
        ValueClass::Float => gatk_barclay::java_float(default)
            .map(|value| Value::Double(f64::from(value)))
            .unwrap_or(Value::Null),
        ValueClass::Long => gatk_barclay::java_long(default)
            .map(Value::Int64)
            .unwrap_or(Value::Null),
        ValueClass::Text | ValueClass::Tagged | ValueClass::Constructed { .. } => {
            Value::Str(default.to_string())
        }
    }
}

/// One declaration as a definition, or nothing when its class is not one this port converts.
pub fn definition(declaration: &Declaration) -> Option<Definition> {
    let class = value_class(declaration.type_name)?;
    // `getArgumentAliases()` is the short name first and then the long one, so an argument with
    // two aliases has a short name and one with a single alias does not.
    let short_name = match declaration.aliases {
        [short, long] if *long == declaration.long_name => *short,
        _ => "",
    };
    let annotation = Annotation {
        full_name: declaration.long_name,
        short_name,
        doc: declaration.doc,
        // `isOptional()` is the annotation OR an initialised field, and what the golden carries is
        // the answer to that question rather than the annotation's own flag. Handing the answer
        // back as the flag reproduces the answer, which is what the parser reads.
        optional: !declaration.required,
        mutex: declaration.mutex,
        min_value: declaration.min_value,
        max_value: declaration.max_value,
        min_recommended_value: declaration.min_recommended_value,
        max_recommended_value: declaration.max_recommended_value,
        suppress_file_expansion: false,
    };
    let initial = initial_value(declaration, &class);
    let mut built = Definition::new(
        annotation,
        declaration.long_name,
        class,
        declaration.collection,
        declaration.primitive,
        initial,
    );
    // An argument a descriptor controls carries the filter that declared it, which the trim needs
    // and the declarations do not have: they record the descriptor. The table comes from the
    // `plugin-argument-ownership` golden, and an argument the golden does not name is left as the
    // tool's own rather than guessed at.
    if declaration.controlled_by.is_some() {
        if let Some(owner) = plugin_ownership::owner(declaration.long_name) {
            built.controlled_by = Some(PluginControl {
                predecessor: owner,
                selector: plugin_ownership::SELECTOR,
            });
        }
    }
    Some(built)
}

/// Every definition a tool's declarations produce, in the parser's own order.
///
/// The list is shorter than the declarations by the arguments whose class is not converted, and
/// [`missing`] says which those are for a given tool.
pub fn definitions(declarations: &[Declaration]) -> Vec<Definition> {
    declarations.iter().filter_map(definition).collect()
}

/// The long names a tool declares that this module cannot build a definition for.
pub fn missing(declarations: &[Declaration]) -> Vec<&'static str> {
    declarations
        .iter()
        .filter(|declaration| definition(declaration).is_none())
        .map(|declaration| declaration.long_name)
        .collect()
}
