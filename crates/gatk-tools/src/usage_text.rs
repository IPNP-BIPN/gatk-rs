//! The usage text Barclay renders for a tool.
//!
//! `gatk <Tool> -h` prints a rendering of the same annotations the argument definitions come from.
//! It is not prose: it has a column layout, a wrapping rule, an ordering and a set of trailer
//! sentences, and all of them are byte-comparable.
//!
//! The documentation strings themselves are data from the annotations and are not held here: this
//! module takes them and lays them out.
//!
//! Ported from `org.broadinstitute.barclay.argparser.CommandLineArgumentParser` (Barclay 5.0.0).

/// `CommandLineArgumentParser.ARGUMENT_COLUMN_WIDTH`: where an argument's description starts.
use crate::tool_declarations::Declaration;

pub const ARGUMENT_COLUMN_WIDTH: usize = 30;
/// `CommandLineArgumentParser.DESCRIPTION_COLUMN_WIDTH`: how wide that description may be.
pub const DESCRIPTION_COLUMN_WIDTH: usize = 90;

/// The heading of the block a plugin descriptor's conditional arguments sit under.
pub const CONDITIONAL_TITLE_PREFIX: &str = "Conditional Arguments for ";

/// The three section titles, in the order they are printed.
pub const REQUIRED_TITLE: &str = "Required Arguments:";
pub const OPTIONAL_TITLE: &str = "Optional Arguments:";
pub const ADVANCED_TITLE: &str = "Advanced Arguments:";

/// One argument, as the usage needs it: the declaration's names and type, and the text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entry {
    /// The long name, without its dashes.
    pub long_name: String,
    /// The other aliases, without their dash, in the order the parser reports them.
    pub short_names: Vec<String>,
    /// The type as the usage writes it, without its angle brackets.
    pub type_name: String,
    /// The documentation and the trailer sentences, already joined: this module wraps it and does
    /// not compose it.
    pub description: String,
}

/// The description column of one entry, composed the way `getArgumentDescription` composes it.
///
/// The order is fixed and every piece carries its own trailing space, which is why the rendered
/// text has two spaces after the documentation and one after everything else:
///
///  1. the documentation, followed by TWO spaces, and only if there is any;
///  2. for a collection, how many times it may be given, which is a different sentence for a
///     required one than for an optional one;
///  3. either `Default value: <rendering>. ` or `Required. `, never both;
///  4. the possible values, which exist for a boolean and for an enum and for nothing else;
///  5. and the mutually exclusive arguments, introduced by a sentence with no leading space of its
///     own, so the piece before it supplies one.
pub fn description(
    doc: &str,
    collection: bool,
    optional: bool,
    default: &str,
    possible_values: Option<&str>,
    mutex: &[(&str, &str)],
) -> String {
    let mut text = String::new();
    if !doc.is_empty() {
        text.push_str(doc);
        text.push_str("  ");
    }
    if collection {
        text.push_str(if optional {
            "This argument may be specified 0 or more times. "
        } else {
            "This argument must be specified at least once. "
        });
    }
    if optional {
        text.push_str(&format!("Default value: {default}. "));
    } else {
        text.push_str("Required. ");
    }
    if let Some(values) = possible_values {
        text.push_str(values);
    }
    if !mutex.is_empty() {
        text.push_str(" Cannot be used in conjunction with argument(s)");
        for (field, short) in mutex {
            text.push(' ');
            text.push_str(field);
            if !short.is_empty() {
                text.push_str(&format!(" ({short})"));
            }
        }
    }
    text
}

/// `getOptionsAsDisplayString`, which answers for a boolean and for an enum and for nothing else.
///
/// The suffix carries its own trailing space, so a description that ends in the possible values
/// ends in a space like every other.
pub fn possible_values(type_name: &str, constants: Option<&[&str]>) -> Option<String> {
    if type_name == "Boolean" {
        return Some("Possible values: {true, false} ".to_string());
    }
    constants.map(|constants| format!("Possible values: {{{}}} ", constants.join(", ")))
}

/// The `--long,-short <Type>` column of one entry.
pub fn name_column(entry: &Entry) -> String {
    let mut names = format!("--{}", entry.long_name);
    for short in &entry.short_names {
        names.push_str(&format!(",-{short}"));
    }
    format!("{names} <{}>", entry.type_name)
}

/// Barclay's wrapping: greedy over spaces, and exactly ONE space is consumed at a break.
///
/// A description that ends in a space therefore keeps it, and a double space (which is what a
/// trailer sentence sits behind) leaves one of the two at the end of the line it broke on. That is
/// why a required argument's first line ends `format ` and its second reads `Required. `.
pub fn wrap(text: &str, width: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut line = String::new();
    for token in text.split(' ') {
        let candidate = if line.is_empty() {
            token.to_string()
        } else {
            format!("{line} {token}")
        };
        if candidate.len() <= width || line.is_empty() {
            line = candidate;
        } else {
            lines.push(line);
            line = token.to_string();
        }
    }
    lines.push(line);
    lines
}

/// One entry, laid out.
///
/// The name column is padded to [`ARGUMENT_COLUMN_WIDTH`] when it fits; when it does not, it takes
/// a line of its own and the description starts on the next one, indented the same.
pub fn entry_lines(entry: &Entry) -> Vec<String> {
    let names = name_column(entry);
    let wrapped = wrap(&entry.description, DESCRIPTION_COLUMN_WIDTH);
    let indent = " ".repeat(ARGUMENT_COLUMN_WIDTH);
    let mut lines = Vec::new();
    // A name column exactly as wide as the column takes no padding and keeps the description on
    // its own line: `--ambig-filter-bases <Integer>Threshold number of ambiguous bases` is thirty
    // characters and then the text, with nothing between them.
    if names.len() <= ARGUMENT_COLUMN_WIDTH {
        let mut first = names.clone();
        first.push_str(&" ".repeat(ARGUMENT_COLUMN_WIDTH - names.len()));
        first.push_str(wrapped.first().map(String::as_str).unwrap_or(""));
        lines.push(first);
        for continuation in wrapped.iter().skip(1) {
            lines.push(format!("{indent}{continuation}"));
        }
    } else {
        lines.push(names);
        for continuation in &wrapped {
            lines.push(format!("{indent}{continuation}"));
        }
    }
    lines
}

/// The banner an experimental or a beta tool prints ABOVE its usage line.
///
/// Two blank lines, the banner, and a blank line, and only then the `USAGE:` line, so the header
/// of such a tool does not start where a reader would expect. The two strings are the reference's
/// own.
pub const EXPERIMENTAL_BANNER: &str = "**EXPERIMENTAL FEATURE - USE AT YOUR OWN RISK**";
pub const BETA_BANNER: &str = "**BETA FEATURE - WORK IN PROGRESS**";

/// Which banner a tool prints, if any.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Maturity {
    Released,
    Beta,
    Experimental,
}

impl Maturity {
    pub fn banner(self) -> Option<&'static str> {
        match self {
            Maturity::Released => None,
            Maturity::Beta => Some(BETA_BANNER),
            Maturity::Experimental => Some(EXPERIMENTAL_BANNER),
        }
    }
}

/// The header: the banner if there is one, the usage line, the summary, and the version.
///
/// The summary is the tool's own `summary` annotation and it is NOT wrapped: a long one is
/// printed as the annotation wrote it, over as many lines as it already had.
pub fn header(tool: &str, summary: &str, version: &str, maturity: Maturity) -> Vec<String> {
    let mut lines = Vec::new();
    if let Some(banner) = maturity.banner() {
        lines.push(String::new());
        lines.push(String::new());
        lines.push(banner.to_string());
        lines.push(String::new());
    }
    lines.push(format!("USAGE: {tool} [arguments]"));
    lines.push(String::new());
    lines.extend(summary.split('\n').map(str::to_string));
    lines.push(format!("Version:{version}"));
    lines.push(String::new());
    lines.push(String::new());
    lines
}

/// One declaration as the usage prints it.
///
/// The type name is the underlying field's, which is what the name column carries; the possible
/// values exist for a boolean and for an enum; and the mutex sentence names the other argument by
/// its FIELD name rather than by its long name, which is why the whole list is passed in.
pub fn entry_for(declaration: &Declaration, all: &[Declaration]) -> Entry {
    let constants =
        crate::tool_declarations::enum_type(declaration.type_name).map(|type_| type_.constants);
    let mutex: Vec<(&str, &str)> = declaration
        .mutex
        .iter()
        .map(|name| {
            let other = all.iter().find(|other| other.long_name == *name);
            let short = other.and_then(short_name).unwrap_or("");
            (*name, short)
        })
        .collect();
    Entry {
        long_name: declaration.long_name.to_string(),
        short_names: short_name(declaration)
            .map(str::to_string)
            .into_iter()
            .collect(),
        type_name: declaration.type_name.to_string(),
        description: description(
            declaration.doc,
            declaration.collection,
            !declaration.required,
            declaration.default.unwrap_or("null"),
            possible_values(declaration.type_name, constants).as_deref(),
            &mutex,
        ),
    }
}

/// The short name the name column prints, which is the alias that is not the long name.
///
/// `getArgumentUsage` prints it only when it differs from the long name, so an argument whose two
/// aliases are the same word prints one.
pub fn short_name(declaration: &Declaration) -> Option<&'static str> {
    match declaration.aliases {
        [short, long] if *long == declaration.long_name && *short != declaration.long_name => {
            Some(short)
        }
        _ => None,
    }
}

/// The three sections a tool's own arguments fall into.
///
/// A hidden argument is printed in none of them, an advanced one in its own, and everything else
/// in required or optional by its own declaration. An argument a plugin descriptor controls is in
/// none of the three either: it belongs to a conditional block.
pub fn sections(list: &[Declaration]) -> (Vec<Entry>, Vec<Entry>, Vec<Entry>) {
    let mut required = Vec::new();
    let mut optional = Vec::new();
    let mut advanced = Vec::new();
    for declaration in list {
        if declaration.hidden || declaration.controlled_by.is_some() {
            continue;
        }
        let entry = entry_for(declaration, list);
        if declaration.advanced {
            advanced.push(entry);
        } else if declaration.required {
            required.push(entry);
        } else {
            optional.push(entry);
        }
    }
    (required, optional, advanced)
}

/// One plugin descriptor's conditional arguments: the arguments a tool only accepts once a plugin
/// has been named.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Conditional {
    /// The descriptor's own name, as the heading writes it.
    pub descriptor: String,
    /// One group per plugin, in the order the reference prints them, each with its arguments.
    pub groups: Vec<(String, Vec<Entry>)>,
}

/// `Valid only if "<plugin>" is specified:`, the line each group opens on.
pub fn conditional_predicate(plugin: &str) -> String {
    format!("Valid only if \"{plugin}\" is specified:")
}

/// The whole text: the header, each section that has an entry in it, and the conditional blocks.
///
/// The entries of a section are ordered by long name, case ignored, which is not the order the
/// parser reports its definitions in. A blank line follows the title and each entry, and a second
/// one closes the section. A conditional group is tighter: its predicate line is followed at once
/// by its first entry.
#[allow(clippy::too_many_arguments)]
pub fn render(
    tool: &str,
    summary: &str,
    version: &str,
    maturity: Maturity,
    required: &[Entry],
    optional: &[Entry],
    advanced: &[Entry],
    conditional: &[Conditional],
) -> String {
    let mut lines = header(tool, summary, version, maturity);
    for (title, entries) in [
        (REQUIRED_TITLE, required),
        (OPTIONAL_TITLE, optional),
        (ADVANCED_TITLE, advanced),
    ] {
        if entries.is_empty() {
            continue;
        }
        lines.push(title.to_string());
        lines.push(String::new());
        let mut sorted: Vec<&Entry> = entries.iter().collect();
        sorted.sort_by(|a, b| {
            a.long_name
                .to_lowercase()
                .cmp(&b.long_name.to_lowercase())
                .then_with(|| a.long_name.cmp(&b.long_name))
        });
        for entry in sorted {
            lines.extend(entry_lines(entry));
            lines.push(String::new());
        }
        lines.push(String::new());
    }
    // The last section's closing blank line is the one the first conditional heading follows, so
    // the heading replaces it rather than adding to it.
    if !conditional.is_empty() {
        lines.pop();
    }
    for block in conditional {
        lines.push(format!("{CONDITIONAL_TITLE_PREFIX}{}:", block.descriptor));
        lines.push(String::new());
        for (plugin, entries) in &block.groups {
            lines.push(conditional_predicate(plugin));
            for entry in entries {
                lines.extend(entry_lines(entry));
                lines.push(String::new());
            }
        }
    }
    // The text ends on a blank line of its own, after the last entry's.
    if !conditional.is_empty() {
        lines.push(String::new());
    }
    lines.join("\n")
}
