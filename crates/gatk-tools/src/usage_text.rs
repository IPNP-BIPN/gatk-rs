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

/// The header: the usage line, the summary, and the version the jar was built from.
///
/// The summary is the tool's own `oneLineSummary` and it is NOT wrapped: a long one is printed as
/// the annotation wrote it, over as many lines as it already had.
pub fn header(tool: &str, summary: &str, version: &str) -> Vec<String> {
    let mut lines = vec![format!("USAGE: {tool} [arguments]"), String::new()];
    lines.extend(summary.split('\n').map(str::to_string));
    lines.push(format!("Version:{version}"));
    lines.push(String::new());
    lines.push(String::new());
    lines
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
pub fn render(
    tool: &str,
    summary: &str,
    version: &str,
    required: &[Entry],
    optional: &[Entry],
    advanced: &[Entry],
    conditional: &[Conditional],
) -> String {
    let mut lines = header(tool, summary, version);
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
