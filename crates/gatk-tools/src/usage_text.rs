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

/// The section titles, in the order they are printed.
pub const POSITIONAL_TITLE: &str = "Positional Arguments:";
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
    // A line break in the description starts a paragraph of its own, which is how a `ClpEnum`'s
    // constants come out one per line.
    text.split('\n')
        .flat_map(|paragraph| wrap_paragraph(paragraph, width))
        .collect()
}

fn wrap_paragraph(text: &str, width: usize) -> Vec<String> {
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
    // A description that ends in a space whose break fell at the very end leaves nothing to put on
    // the last line, and the reference emits no line rather than an indent with nothing after it.
    // That is why `--platform-filter-name`'s entry ends `Required.` with no trailing space in the
    // usage golden while `--library`'s, which did not wrap there, ends `Required. ` with one.
    if !line.is_empty() || lines.is_empty() {
        lines.push(line);
    }
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
    entry_for_in(declaration, all, None)
}

/// The same, for a named tool, whose plugin descriptors answer for some of its arguments and
/// which may deprecate some.
///
/// `--read-filter` and `--disable-read-filter` are `String`s with no enum behind them, and their
/// possible values come from the descriptor: the whole catalogue for the first and the TOOL'S OWN
/// DEFAULTS for the second; the annotation descriptor's three the same way. All of them, and the
/// deprecations, are measured per tool in the `usage-text` golden.
pub fn entry_for_in(declaration: &Declaration, all: &[Declaration], tool: Option<&str>) -> Entry {
    // `--read-filter`, `--disable-read-filter`, `--annotation` and their siblings are `String`s
    // with no enum behind them: their possible values are the descriptor's, per tool.
    let descriptor_values = tool
        .and_then(|tool| crate::tool_declarations::descriptor_values(tool, declaration.long_name));
    let enum_type = crate::tool_declarations::enum_type(declaration.type_class);
    let values = match (descriptor_values, enum_type) {
        // A descriptor that lists nothing, `--annotations-to-exclude` with no default annotation,
        // allows anything, and says so with no trailing space.
        (Some([]), _) => Some("Any value allowed".to_string()),
        (Some(values), _) => possible_values(declaration.type_name, Some(values)),
        // A `ClpEnum` documents each constant, and prints one per line with no `Possible values`
        // and no trailing space, so a mutex sentence follows the last constant on its line.
        (None, Some(type_)) if !type_.docs.is_empty() => Some(
            type_
                .docs
                .iter()
                .map(|(constant, doc)| format!("{constant} ({doc})"))
                .collect::<Vec<_>>()
                .join("\n"),
        ),
        (None, type_) => possible_values(declaration.type_name, type_.map(|type_| type_.constants)),
    };
    let mutex: Vec<(&str, &str)> = declaration
        .mutex
        .iter()
        .map(|name| {
            let other = all.iter().find(|other| other.long_name == *name);
            // Unlike the name column, the sentence prints the short name even when it repeats
            // the long one: `DO_NOT_VALIDATE_FILTERED (do-not-validate-filtered-records)`.
            let short = match other.map(|other| (other.aliases, other.long_name)) {
                Some(([short, long], name)) if *long == name => *short,
                _ => "",
            };
            // The sentence names the target's FIELD, which is the long name only where the field
            // is named after the argument.
            let printed = crate::plugin_ownership::mutex_field_name(name).unwrap_or(name);
            (printed, short)
        })
        .collect();
    // A deprecated argument says so in front of its documentation, with the reason; which ones
    // are is the tool's and is measured in the usage golden.
    let doc = match tool
        .and_then(|tool| crate::tool_declarations::deprecation(tool, declaration.long_name))
    {
        Some(detail) => format!(
            "This argument is DEPRECATED ({detail}). {}",
            declaration.doc
        ),
        None => declaration.doc.to_string(),
    };
    // `COMPRESSION_LEVEL`'s default is htsjdk's `Defaults.COMPRESSION_LEVEL`, read from
    // `samjdk.compression_level`, which `Main` has set from `GATKConfig` before the tool exists.
    let default = match declaration.long_name {
        "COMPRESSION_LEVEL" => crate::gatk_config::default(crate::gatk_config::COMPRESSION_LEVEL),
        _ => declaration.default,
    }
    .unwrap_or("null");
    Entry {
        long_name: declaration.long_name.to_string(),
        short_names: short_name(declaration)
            .map(str::to_string)
            .into_iter()
            .collect(),
        type_name: declaration.type_name.to_string(),
        description: description(
            &doc,
            declaration.collection,
            !declaration.required,
            default,
            values.as_deref(),
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
    sections_for(None, list)
}

/// The same, for a named tool, whose descriptor answers for two of its arguments.
pub fn sections_for(
    tool: Option<&str>,
    list: &[Declaration],
) -> (Vec<Entry>, Vec<Entry>, Vec<Entry>) {
    let mut required = Vec::new();
    let mut optional = Vec::new();
    let mut advanced = Vec::new();
    for declaration in list {
        if declaration.hidden || declaration.controlled_by.is_some() {
            continue;
        }
        let entry = entry_for_in(declaration, list, tool);
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

/// The conditional blocks a tool's controlled arguments fall into, one group per plugin.
///
/// One block per descriptor, `annotation` before `readFilter`. The read filters' groups are in the
/// ownership table's order, which is the golden's: by the owner's class name and then by the
/// argument's long name. A plugin that declares no argument has no group, which is why fifty-six
/// filters make twenty-two blocks.
pub fn conditionals(list: &[Declaration]) -> Vec<Conditional> {
    let mut blocks = Vec::new();
    // The descriptors are in the order of their names, which puts `annotation` before
    // `readFilter`. An annotation's arguments are grouped by the annotation that declares them,
    // the annotations and the arguments under each in the order of their names.
    let mut annotations: Vec<(&str, &Declaration)> = list
        .iter()
        .filter(|declaration| !declaration.hidden)
        .filter_map(|declaration| {
            crate::plugin_ownership::annotation_owner(declaration.long_name)
                .filter(|_| declaration.controlled_by == Some(ANNOTATION_DESCRIPTOR_CLASS))
                .map(|owner| (owner, declaration))
        })
        .collect();
    annotations.sort_by(|a, b| (a.0, a.1.long_name).cmp(&(b.0, b.1.long_name)));
    let mut groups: Vec<(String, Vec<Entry>)> = Vec::new();
    for (owner, declaration) in annotations {
        let rendered = entry_for(declaration, list);
        match groups.last_mut() {
            Some((last, entries)) if last == owner => entries.push(rendered),
            _ => groups.push((owner.to_string(), vec![rendered])),
        }
    }
    if !groups.is_empty() {
        blocks.push(Conditional {
            descriptor: ANNOTATION_DESCRIPTOR.to_string(),
            groups,
        });
    }
    let mut groups: Vec<(String, Vec<Entry>)> = Vec::new();
    for entry in crate::plugin_ownership::OWNERSHIP.iter() {
        let Some(declaration) = list
            .iter()
            .find(|declaration| declaration.long_name == entry.long_name)
        else {
            continue;
        };
        // `VariantEval` declares a `--sample` of its own, which is not `SampleReadFilter`'s.
        if declaration.hidden || declaration.controlled_by.is_none() {
            continue;
        }
        let rendered = entry_for(declaration, list);
        match groups.iter_mut().find(|(owner, _)| owner == entry.owner) {
            Some((_, entries)) => entries.push(rendered),
            None => groups.push((entry.owner.to_string(), vec![rendered])),
        }
    }
    if !groups.is_empty() {
        blocks.push(Conditional {
            descriptor: crate::plugin_ownership::DESCRIPTOR.to_string(),
            groups,
        });
    }
    blocks
}

/// The annotation descriptor, by the class a declaration names and the name a heading prints.
pub const ANNOTATION_DESCRIPTOR_CLASS: &str = "GATKAnnotationPluginDescriptor";
pub const ANNOTATION_DESCRIPTOR: &str = "annotation";

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
    render_with_positional(
        tool,
        summary,
        version,
        maturity,
        &[],
        required,
        optional,
        advanced,
        conditional,
    )
}

/// The same, with the positional arguments a tool declares in a section of their own, first.
///
/// The header's second blank line opens the first of the positional and the required sections;
/// a tool with neither goes straight from one blank line to its optional arguments.
#[allow(clippy::too_many_arguments)]
pub fn render_with_positional(
    tool: &str,
    summary: &str,
    version: &str,
    maturity: Maturity,
    positional: &[Entry],
    required: &[Entry],
    optional: &[Entry],
    advanced: &[Entry],
    conditional: &[Conditional],
) -> String {
    let mut lines = header(tool, summary, version, maturity);
    if positional.is_empty() && required.is_empty() {
        lines.pop();
    }
    for (title, entries) in [
        (POSITIONAL_TITLE, positional),
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
