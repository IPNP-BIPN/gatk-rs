//! Conformance for the plugin trim against GATK 4.6.2.0.
//!
//! The rule was measured in `barclay-plugin-descriptors` and ported into
//! [`gatk_barclay::Parser`]: a controlled argument whose plugin nobody selected and nobody set
//! leaves the definition list before the required check. This suite is about the TABLE that rule
//! runs over, measured in `plugin-argument-ownership`, and about what the two together do to a
//! walker's command line.
//!
//! # What this suite is for
//!
//!  * **the ownership table matching the golden row for row**, owner, short name and required;
//!  * **the twelve required controlled arguments, which are what makes the trim load-bearing**;
//!  * **the command lines: a required argument of nobody dropped when absent, refused with the
//!    CLASS named when given, and still required once its filter is named**;
//!  * **a tool default counting as selected**;
//!  * **and a walker being parseable now, while its usage still is not.**

use gatk_barclay::Parser;
use gatk_corpus as corpus;
use gatk_tools::plugin_ownership::{self, OWNERSHIP};
use gatk_tools::tool_declarations::{Declaration, COUNTREADS};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/plugin_argument_ownership.txt.gz"),
    )
}

fn rows(text: &str, kind: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|line| line.starts_with(&format!("{kind}\t")))
        .map(|line| line.split('\t').skip(1).map(str::to_string).collect())
        .collect()
}

/// The table is the golden's, row for row and in its order.
#[test]
fn the_ownership_table_is_the_goldens() {
    let text = golden();
    let recorded = rows(&text, "owner");
    assert_eq!(recorded.len(), OWNERSHIP.len());
    for (row, entry) in recorded.iter().zip(OWNERSHIP.iter()) {
        assert_eq!(row[0], entry.owner);
        assert_eq!(row[1], entry.long_name);
        // The golden writes a bare dash where there is no short name.
        assert_eq!(
            row[2],
            if entry.short_name.is_empty() {
                "-"
            } else {
                entry.short_name
            }
        );
        assert_eq!(row[3] == "required", entry.required, "{}", entry.long_name);
    }
    // Twelve of the twenty-eight are required, which is the count the trim has to remove from a
    // plain command line before the required check sees it.
    assert_eq!(OWNERSHIP.iter().filter(|entry| entry.required).count(), 12);
    // The names do not name their filters, which is why the table had to be measured.
    assert_eq!(
        plugin_ownership::owner("library"),
        Some("LibraryReadFilter")
    );
    assert_eq!(
        plugin_ownership::owner("black-listed-lanes"),
        Some("PlatformUnitReadFilter")
    );
    assert_eq!(plugin_ownership::owner("tool-arg"), None);
}

/// Every argument CountReads declares as controlled is one the table names.
#[test]
fn the_table_covers_what_the_declarations_control() {
    let controlled: Vec<&Declaration> = COUNTREADS
        .iter()
        .filter(|declaration| declaration.controlled_by.is_some())
        .collect();
    assert_eq!(controlled.len(), OWNERSHIP.len());
    for declaration in &controlled {
        let entry = plugin_ownership::ownership(declaration.long_name)
            .unwrap_or_else(|| panic!("{}", declaration.long_name));
        // The declarations and the ownership dump agree about which of them are required, which
        // they were measured by two different programs to answer.
        assert_eq!(entry.required, declaration.required, "{}", entry.long_name);
    }
    // And the argument that selects them is the tool's own rather than a plugin's.
    let selector = COUNTREADS
        .iter()
        .find(|declaration| declaration.long_name == plugin_ownership::SELECTOR)
        .expect("the selector");
    assert!(selector.controlled_by.is_none());
    assert!(selector.collection);
}

/// A parser over the shape the dump measured: the controlled arguments and the ones that select
/// them, without the rest of a tool's surface.
fn measured_parser() -> Parser {
    let list: Vec<Declaration> = COUNTREADS
        .iter()
        .filter(|declaration| {
            declaration.controlled_by.is_some()
                || declaration.long_name == plugin_ownership::SELECTOR
                || declaration.long_name.starts_with("disable-")
        })
        .cloned()
        .collect();
    Parser::new(gatk_cli::definitions::definitions(&list))
}

fn result(text: &str, label: &str) -> String {
    text.lines()
        .find(|line| line.starts_with(&format!("result\t{label}\t")))
        .map(|line| line.rsplit('\t').next().expect("a result").to_string())
        .unwrap_or_else(|| panic!("{label}"))
}

/// The command lines the golden recorded, replayed through the ported parser.
#[test]
fn the_trim_answers_the_goldens_command_lines() {
    let text = golden();

    // Nothing named: twelve required arguments in the list and none of them fires.
    assert_eq!(result(&text, "nothing-named"), "ok");
    assert!(measured_parser().parse_arguments(&[]).is_ok());

    // A required argument of a filter nobody named, given. The message names the CLASS, and the
    // argument's own name is the short name, a slash and the long name with no guard.
    let recorded = result(&text, "a-required-argument-of-nobody");
    let refusal = measured_parser()
        .parse_arguments(&["--library", "lib1"])
        .expect_err("the refusal");
    assert!(
        recorded.ends_with(&refusal.message),
        "{recorded} / {}",
        refusal.message
    );
    assert!(refusal.message.contains("\"library/library\""));
    assert!(refusal.message.contains("\"LibraryReadFilter\""));
    // The same for one with no short name, which reports a leading slash.
    let recorded = result(&text, "a-default-filter-argument");
    let refusal = measured_parser()
        .parse_arguments(&["--ambig-filter-frac", "0.1"])
        .expect_err("the refusal");
    assert!(
        recorded.ends_with(&refusal.message),
        "{recorded} / {}",
        refusal.message
    );
    assert!(refusal.message.contains("\"/ambig-filter-frac\""));

    // Its filter named, so the argument is allowed.
    assert_eq!(result(&text, "its-filter-named"), "ok");
    assert!(measured_parser()
        .parse_arguments(&["--read-filter", "LibraryReadFilter", "--library", "lib1"])
        .is_ok());

    // Named without it, and the required check does fire: the trim keeps a selected plugin's
    // required argument in the list.
    let recorded = result(&text, "its-filter-named-without-it");
    let refusal = measured_parser()
        .parse_arguments(&["--read-filter", "LibraryReadFilter"])
        .expect_err("the missing argument");
    assert!(
        recorded.ends_with(&refusal.message),
        "{recorded} / {}",
        refusal.message
    );
    assert_eq!(
        refusal.message,
        "Argument library was missing: Argument 'library' is required"
    );

    // A second one, so the first is not a special case.
    assert_eq!(result(&text, "a-second-required-argument"), "ok");
    assert!(measured_parser()
        .parse_arguments(&[
            "--read-filter",
            "PlatformReadFilter",
            "--platform-filter-name",
            "ILLUMINA"
        ])
        .is_ok());
    let recorded = result(&text, "a-second-required-argument-missing");
    let refusal = measured_parser()
        .parse_arguments(&["--read-filter", "PlatformReadFilter"])
        .expect_err("the missing argument");
    assert!(
        recorded.ends_with(&refusal.message),
        "{recorded} / {}",
        refusal.message
    );
}

/// A tool default counts as selected, which is what the `allowed` rows say.
#[test]
fn a_tool_default_counts_as_selected() {
    let text = golden();
    // The rows the port can answer are the ones whose filter owns an argument: the predicate is
    // read off a command line that sets one. `WellformedReadFilter` and `MappedReadFilter` declare
    // none, which is why the dump asked the descriptor directly and this asks the parser.
    for row in rows(&text, "allowed") {
        let Some(argument) = OWNERSHIP.iter().find(|entry| entry.owner == row[1]) else {
            continue;
        };
        let defaults: Vec<String> = if row[0] == "no-defaults" {
            Vec::new()
        } else {
            vec![
                "WellformedReadFilter".to_string(),
                "MappedReadFilter".to_string(),
            ]
        };
        let allowed = measured_parser()
            .with_default_plugins(defaults)
            .parse_arguments(&[&format!("--{}", argument.long_name), "lib1"])
            .is_ok();
        assert_eq!(allowed, row[2] == "true", "{} {}", row[0], row[1]);
    }
    // The same filter as a default is allowed, with no `--read-filter` on the command line: the
    // descriptor does not distinguish a default from a name the command line gave it.
    assert!(measured_parser()
        .with_default_plugins(vec!["AmbiguousBaseReadFilter".to_string()])
        .parse_arguments(&["--ambig-filter-frac", "0.1"])
        .is_ok());
    assert_eq!(result(&text, "a-default-filters-own-argument"), "ok");
    assert!(measured_parser()
        .parse_arguments(&[
            "--read-filter",
            "AmbiguousBaseReadFilter",
            "--ambig-filter-frac",
            "0.1"
        ])
        .is_ok());
}

/// A walker parses now, and its usage is composed too.
#[test]
fn a_walker_parses_and_so_does_its_usage() {
    assert!(gatk_cli::parseable("CountReads"));
    assert!(gatk_cli::parseable("PrintReads"));
    assert!(gatk_cli::parseable("IndexFeatureFile"));
    // A whole command line, which before the table was refused for twelve arguments the reference
    // never asks for.
    assert_eq!(
        gatk_cli::parse_failure(
            "CountReads",
            &["--input".to_string(), "/dev/null".to_string()]
        ),
        None
    );
    // And the usage, whose conditional blocks the same table orders. `walker_usage.rs` compares it
    // against the golden; what is asserted here is that the two questions have one answer now.
    assert!(gatk_cli::usage_composable("CountReads"));
    assert!(gatk_cli::usage_composable("IndexFeatureFile"));
    let usage = gatk_cli::tool_usage("CountReads").expect("a walker's usage");
    assert!(usage.contains("Conditional Arguments for readFilter:"));
    assert!(gatk_cli::tool_usage("IndexFeatureFile").is_some());
}
