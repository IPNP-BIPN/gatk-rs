//! Conformance for the usage text against GATK 4.6.2.0 and Barclay 5.0.0.
//!
//! Golden from `tools/argument-conformance/UsageTextDump.java`, which holds three tools' whole
//! usage. The test takes each entry apart at the column boundary, hands the pieces back to the
//! port, and compares what it lays out against the bytes the reference wrote.
//!
//! The documentation strings are the annotations' and are not ported; the layout is.
//!
//! # What this suite is for
//!
//!  * **the two column widths**;
//!  * **the wrapping consuming exactly one space at a break**;
//!  * **a name column that does not fit taking a line of its own**;
//!  * **the sections, their titles and their blank lines**;
//!  * **the entries being ordered by long name**;
//!  * **and the header carrying the version.**

use gatk_corpus as corpus;
use gatk_tools::usage_text::{
    conditional_predicate, entry_lines, header, name_column, render, wrap, Conditional, Entry,
    ADVANCED_TITLE, ARGUMENT_COLUMN_WIDTH, CONDITIONAL_TITLE_PREFIX, DESCRIPTION_COLUMN_WIDTH,
    OPTIONAL_TITLE, REQUIRED_TITLE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/usage_text.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn usage(text: &str, tool: &str) -> Vec<String> {
    let prefix = format!("usage\t{tool}\t");
    let body = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
        .unwrap_or_else(|| panic!("usage/{tool}"));
    body.split('\n').map(str::to_string).collect()
}

/// One entry as the golden wrote it: the lines it occupies, and the pieces the port needs.
struct Written {
    lines: Vec<String>,
    entry: Entry,
}

/// One section: its title and the entries under it.
type Section = (String, Vec<Written>);
/// One conditional block: the descriptor's name, and one group per plugin.
type Block = (String, Vec<Section>);
/// A parsed usage: its header lines, its sections, and its conditional blocks.
type Parsed = (Vec<String>, Vec<Section>, Vec<Block>);

fn parse(lines: &[String]) -> Parsed {
    let indent = " ".repeat(ARGUMENT_COLUMN_WIDTH);
    let mut sections: Vec<Section> = Vec::new();
    let mut conditional: Vec<Block> = Vec::new();
    let mut header = Vec::new();
    let mut index = 0;
    let entry_at = |index: &mut usize| -> Written {
        let mut owned = vec![lines[*index].clone()];
        *index += 1;
        while *index < lines.len() && lines[*index].starts_with(&indent) {
            owned.push(lines[*index].clone());
            *index += 1;
        }
        written(&owned)
    };
    while index < lines.len() && !lines[index].ends_with("Arguments:") {
        header.push(lines[index].clone());
        index += 1;
    }
    while index < lines.len() {
        if lines[index].is_empty() {
            index += 1;
            continue;
        }
        if lines[index].starts_with(CONDITIONAL_TITLE_PREFIX) {
            let descriptor = lines[index]
                .trim_start_matches(CONDITIONAL_TITLE_PREFIX)
                .trim_end_matches(':')
                .to_string();
            index += 2;
            let mut groups: Vec<Section> = Vec::new();
            while index < lines.len() && lines[index].starts_with("Valid only if ") {
                let plugin = lines[index]
                    .trim_start_matches("Valid only if \"")
                    .trim_end_matches("\" is specified:")
                    .to_string();
                index += 1;
                let mut entries = Vec::new();
                while index < lines.len() && lines[index].starts_with("--") {
                    entries.push(entry_at(&mut index));
                    index += 1; // the blank line after the entry
                }
                groups.push((plugin, entries));
            }
            conditional.push((descriptor, groups));
            continue;
        }
        let title = lines[index].clone();
        index += 2; // the title and the blank line under it
        let mut entries = Vec::new();
        while index < lines.len() && lines[index].starts_with("--") {
            entries.push(entry_at(&mut index));
            index += 1; // the blank line after the entry
        }
        sections.push((title, entries));
        // A section is closed by a second blank line, unless what follows is a conditional block,
        // which the last entry's own blank line is enough to separate.
        if index < lines.len() && lines[index].is_empty() {
            index += 1;
        }
    }
    (header, sections, conditional)
}

/// The pieces of one written entry: its names, its type and its description, put back together.
///
/// The description is rebuilt by joining the wrapped segments with ONE space, which is exactly
/// what the wrapping consumed, so the join is lossless.
fn written(lines: &[String]) -> Written {
    let names_line = &lines[0];
    // The names end at the type's closing bracket. When that lands before the column boundary the
    // description follows on the same line, and `--ambig-filter-bases <Integer>` ends exactly ON
    // it, with the text immediately after and no space between.
    let close = names_line.find('>').expect("a type");
    let (names, first) = if close < ARGUMENT_COLUMN_WIDTH {
        (
            &names_line[..=close],
            names_line[ARGUMENT_COLUMN_WIDTH..].to_string(),
        )
    } else {
        (names_line.as_str(), String::new())
    };
    let mut segments = Vec::new();
    if !first.is_empty() {
        segments.push(first);
    }
    for line in &lines[1..] {
        segments.push(line[ARGUMENT_COLUMN_WIDTH..].to_string());
    }
    let (declaration, type_name) = names.split_once(" <").expect("a type");
    let mut aliases = declaration.split(',');
    let long_name = aliases
        .next()
        .expect("a long name")
        .trim_start_matches("--")
        .to_string();
    Written {
        lines: lines.to_vec(),
        entry: Entry {
            long_name,
            short_names: aliases
                .map(|a| a.trim_start_matches('-').to_string())
                .collect(),
            type_name: type_name.trim_end_matches('>').to_string(),
            description: segments.join(" "),
        },
    }
}

/// Every entry of every tool, laid out again, comes back byte for byte.
#[test]
fn the_port_lays_out_every_entry_as_written() {
    let text = golden();
    let mut checked = 0;
    for tool in ["CountReads", "IndexFeatureFile", "GatherVcfsCloud"] {
        let (_, sections, conditional) = parse(&usage(&text, tool));
        assert!(!sections.is_empty(), "{tool}");
        for (_, entries) in sections {
            for written in entries {
                assert_eq!(entry_lines(&written.entry), written.lines, "{tool}");
                checked += 1;
            }
        }
        for (_, groups) in conditional {
            for (_, entries) in groups {
                for written in entries {
                    assert_eq!(entry_lines(&written.entry), written.lines, "{tool}");
                    checked += 1;
                }
            }
        }
    }
    // Over a hundred entries between the three tools, the conditional ones included.
    assert!(checked > 100, "{checked} entries laid out");
}

/// The whole text comes back too, header and sections and blank lines.
#[test]
fn the_port_renders_the_whole_text() {
    let text = golden();
    for tool in ["CountReads", "IndexFeatureFile", "GatherVcfsCloud"] {
        let lines = usage(&text, tool);
        let (head, sections, conditional) = parse(&lines);
        // The summary is everything between the blank line and the version, which is more than
        // one line for a tool whose one-line summary is not one line.
        let version_at = head
            .iter()
            .position(|l| l.starts_with("Version:"))
            .expect("a version");
        let summary = head[2..version_at].join("\n");
        let version = head[version_at].trim_start_matches("Version:").to_string();
        let section = |title: &str| -> Vec<Entry> {
            sections
                .iter()
                .find(|(name, _)| name == title)
                .map(|(_, entries)| entries.iter().map(|w| w.entry.clone()).collect())
                .unwrap_or_default()
        };
        let blocks: Vec<Conditional> = conditional
            .iter()
            .map(|(descriptor, groups)| Conditional {
                descriptor: descriptor.clone(),
                groups: groups
                    .iter()
                    .map(|(plugin, entries)| {
                        (
                            plugin.clone(),
                            entries.iter().map(|w| w.entry.clone()).collect(),
                        )
                    })
                    .collect(),
            })
            .collect();
        let ours = render(
            tool,
            &summary,
            &version,
            &section(REQUIRED_TITLE),
            &section(OPTIONAL_TITLE),
            &section(ADVANCED_TITLE),
            &blocks,
        );
        // Compare line by line: a three-hundred-line panic message says nothing, and the index
        // of the first difference says everything.
        let theirs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let mine: Vec<&str> = ours.split('\n').collect();
        for (index, (a, b)) in mine.iter().zip(&theirs).enumerate() {
            assert_eq!(a, b, "{tool} line {index}");
        }
        assert_eq!(mine.len(), theirs.len(), "{tool} line count");
    }
}

/// The wrapping consumes exactly one space at a break, which is why a trailer keeps one.
#[test]
fn the_wrapping_consumes_one_space() {
    assert_eq!(wrap("one two three", 7), vec!["one two", "three"]);
    // A double space leaves one behind on the line it broke on.
    assert_eq!(wrap("aaa  bbb", 5), vec!["aaa ", "bbb"]);
    // A word longer than the width is not broken.
    assert_eq!(wrap("aaaaaaaa bb", 4), vec!["aaaaaaaa", "bb"]);
    // And that is the rule the golden's required argument shows.
    let text = golden();
    let (_, sections, _) = parse(&usage(&text, "IndexFeatureFile"));
    let (title, required) = &sections[0];
    assert_eq!(title, REQUIRED_TITLE);
    assert_eq!(required.len(), 1);
    assert!(required[0].lines[0].ends_with("format "));
    assert_eq!(required[0].lines[1].trim(), "Required.");
    assert!(required[0].lines[1].ends_with("Required. "));
}

/// A name column that does not fit takes a line of its own.
#[test]
fn a_long_name_column_takes_its_own_line() {
    let short = Entry {
        long_name: "input".to_string(),
        short_names: vec!["I".to_string()],
        type_name: "GATKPath".to_string(),
        description: "a doc  Required. ".to_string(),
    };
    assert_eq!(name_column(&short), "--input,-I <GATKPath>");
    assert!(entry_lines(&short)[0].starts_with("--input,-I <GATKPath>         a doc"));
    let long = Entry {
        long_name: "use-jdk-deflater".to_string(),
        short_names: vec!["jdk-deflater".to_string()],
        type_name: "Boolean".to_string(),
        description: "a doc ".to_string(),
    };
    let lines = entry_lines(&long);
    assert_eq!(lines[0], "--use-jdk-deflater,-jdk-deflater <Boolean>");
    assert_eq!(
        lines[1],
        format!("{}a doc ", " ".repeat(ARGUMENT_COLUMN_WIDTH))
    );
    // Which is what the golden does with that very argument.
    let text = golden();
    let (_, sections, _) = parse(&usage(&text, "IndexFeatureFile"));
    let written = sections
        .iter()
        .flat_map(|(_, entries)| entries)
        .find(|w| w.entry.long_name == "use-jdk-deflater")
        .expect("the argument");
    assert_eq!(
        written.lines[0],
        "--use-jdk-deflater,-jdk-deflater <Boolean>"
    );
}

/// The sections, their order, and the header's version.
#[test]
fn the_sections_and_the_header() {
    let text = golden();
    let lines = usage(&text, "CountReads");
    let (head, sections, conditional) = parse(&lines);
    assert_eq!(head[0], "USAGE: CountReads [arguments]");
    assert_eq!(head[3], "Version:4.6.2.0");
    assert_eq!(header("CountReads", &head[2], "4.6.2.0"), head);
    // A tool whose one-line summary is not one line keeps its own line breaks.
    let (other, _, _) = parse(&usage(&text, "GatherVcfsCloud"));
    let version_at = other
        .iter()
        .position(|l| l.starts_with("Version:"))
        .expect("a version");
    assert!(version_at > 3, "the summary runs over more than one line");
    assert_eq!(
        header(
            "GatherVcfsCloud",
            &other[2..version_at].join("\n"),
            "4.6.2.0"
        ),
        other
    );
    let titles: Vec<&str> = sections.iter().map(|(title, _)| title.as_str()).collect();
    assert_eq!(titles, vec![REQUIRED_TITLE, OPTIONAL_TITLE, ADVANCED_TITLE]);
    // Each section's entries are ordered by long name, which the parser's own order is not.
    for (title, entries) in &sections {
        // The order ignores case, which is why `--output` comes before `--QUIET`.
        let names: Vec<&str> = entries.iter().map(|w| w.entry.long_name.as_str()).collect();
        let mut sorted = names.clone();
        sorted.sort_by_key(|name| name.to_lowercase());
        assert_eq!(names, sorted, "{title}");
    }
    // A tool that is no walker has no read-filter arguments to print, and no conditional block
    // either: the block exists because a plugin descriptor does.
    let (_, plain, plain_conditional) = parse(&usage(&text, "IndexFeatureFile"));
    assert!(!plain
        .iter()
        .flat_map(|(_, entries)| entries)
        .any(|w| w.entry.long_name == "read-filter"));
    assert!(plain_conditional.is_empty());
    assert!(sections
        .iter()
        .flat_map(|(_, entries)| entries)
        .any(|w| w.entry.long_name == "read-filter"));
    // The walker's block is the read-filter descriptor's, and each group opens on its predicate.
    assert_eq!(conditional.len(), 1);
    assert_eq!(conditional[0].0, "readFilter");
    let (plugin, entries) = &conditional[0].1[0];
    assert_eq!(plugin, "AmbiguousBaseReadFilter");
    assert_eq!(entries[0].entry.long_name, "ambig-filter-bases");
    assert_eq!(
        conditional_predicate(plugin),
        "Valid only if \"AmbiguousBaseReadFilter\" is specified:"
    );
    assert_eq!(DESCRIPTION_COLUMN_WIDTH, 90);
}

/// The whole usage of a tool that is no walker, composed from its declarations.
///
/// Every other test here takes the golden apart and hands the pieces back. This one starts from
/// the declarations alone: the documentation, the type, the default, the flags and the enum
/// constants, all of them measured, and asks the port to write the file the reference wrote.
#[test]
fn a_tools_usage_is_composed_from_its_declarations() {
    use gatk_tools::tool_declarations::INDEXFEATUREFILE;
    use gatk_tools::usage_text::sections;

    let text = golden();
    let expected = usage(&text, "IndexFeatureFile").join("\n");
    let (required, optional, advanced) = sections(INDEXFEATUREFILE);
    // Every one of the fourteen declarations is printed, none being hidden or plugin-controlled.
    assert_eq!(
        required.len() + optional.len() + advanced.len(),
        INDEXFEATUREFILE.len()
    );
    let written = render(
        "IndexFeatureFile",
        "Creates an index for a feature file, e.g. VCF or BED file.",
        "4.6.2.0",
        &required,
        &optional,
        &advanced,
        &[],
    );
    assert_eq!(written, expected);
}
