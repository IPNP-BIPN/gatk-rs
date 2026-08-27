//! Conformance for `Main`'s dispatch against GATK 4.6.2.0, compared as the message a tool name
//! that resolves to nothing produces.
//!
//! Golden from `tools/readfilter-conformance/MainDispatchDump.java`, which asks the reference's
//! own `getUnknownCommandMessage` over a fixed catalogue of five classes rather than its whole
//! class path: what is pinned is the search and not the catalogue.
//!
//! # What this suite is for
//!
//!  * **a deprecated tool short-circuiting the search**;
//!  * **a prefix scoring zero whatever its length**;
//!  * **a substring scoring zero only from five characters**;
//!  * **the distance's weights, a deletion costing four against an insertion's one**;
//!  * **the floor being seven**;
//!  * **`this` becoming `one of these` at two matches**;
//!  * **the suggestions running together with no separator**;
//!  * **every tool scoring zero suppressing the suggestion**;
//!  * **and a name that resolves being a refusal rather than an answer.**

use gatk_corpus as corpus;
use gatk_tools::main_dispatch::{
    distance, levenshtein_distance, suggested_alternate_command, tool_deprecation_info,
    unknown_command_message, DELETION_COST, DEPRECATED_TOOLS, HELP_SIMILARITY_FLOOR,
    INSERTION_COST, MINIMUM_SUBSTRING_LENGTH, SUBSTITUTION_COST, SWAP_COST,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/main_dispatch.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// The catalogue the dump scored against, which it prints itself.
fn catalogue(text: &str) -> Vec<String> {
    field(text, "tools", "catalogue")
        .expect("the catalogue is in the golden")
        .split(',')
        .map(str::to_string)
        .collect()
}

/// The command each case asked about, which the message quotes back.
const COMMANDS: &[(&str, &str)] = &[
    ("deprecated-short-circuits", "IndelRealigner"),
    ("prefix-one-character", "P"),
    ("prefix-four-characters", "Prin"),
    ("substring-four-characters", "ount"),
    ("substring-five-characters", "ountR"),
    ("one-insertion", "PrintReadss"),
    ("one-deletion", "PrntReads"),
    ("one-substitution", "PrintReadz"),
    ("two-too-many", "PrintReadsxy"),
    ("two-too-few", "PrtReads"),
    ("nothing-close", "Zzzzzzzzzzzzzzzzzzzz"),
    ("every-tool-scores-zero", ""),
    ("prefix-beats-a-neighbour", "CountRead"),
    ("single-tool-zero", "Print"),
];

/// Every case's message is what the port reaches, over the catalogue the golden names.
#[test]
fn every_case_reaches_the_same_message() {
    let text = golden();
    let catalogue = catalogue(&text);
    for (case, command) in COMMANDS {
        let classes = if *case == "single-tool-zero" {
            vec!["PrintReads".to_string()]
        } else {
            catalogue.clone()
        };
        let expected = field(&text, "message", case).unwrap_or_else(|| panic!("{case}"));
        assert_eq!(
            unknown_command_message(&classes, command).expect("no refusal"),
            expected,
            "{case}"
        );
    }
    // And the two-tool catalogue, which the dump names separately.
    let two = vec!["CountReads".to_string(), "CountBases".to_string()];
    assert_eq!(
        unknown_command_message(&two, "ount").expect("no refusal"),
        field(&text, "message", "substring-of-two").expect("substring-of-two")
    );
}

/// Every deprecated tool's notice is the registry's own sentence, and a live tool has none.
#[test]
fn the_deprecation_notices_are_the_registrys_own() {
    let text = golden();
    for (tool, _, _) in DEPRECATED_TOOLS {
        assert_eq!(
            tool_deprecation_info(tool),
            field(&text, "deprecated", tool),
            "{tool}"
        );
    }
    assert_eq!(tool_deprecation_info("PrintReads"), None);
    assert_eq!(
        field(&text, "deprecated", "PrintReads").as_deref(),
        Some("null")
    );
    // The notice short-circuits the search, so no suggestion follows it.
    let message = field(&text, "message", "deprecated-short-circuits").expect("the message");
    assert_eq!(
        message,
        tool_deprecation_info("IndelRealigner").expect("a notice")
    );
    assert!(!message.contains("Did you mean"));
}

/// A prefix scores zero whatever its length, while a substring needs five characters.
#[test]
fn a_substring_needs_five_characters() {
    assert_eq!(distance("P", "PrintReads"), Ok(0));
    assert_eq!(distance("Prin", "PrintReads"), Ok(0));
    assert_eq!(distance("ountR", "CountReads"), Ok(0));
    assert_eq!(MINIMUM_SUBSTRING_LENGTH, 5);
    // Four characters fall back on the distance, which still finds it: six insertions.
    assert_eq!(distance("ount", "CountReads"), Ok(6));
    assert_eq!(distance("ount", "CountBases"), Ok(6));
}

/// Dropping a character from the command costs four and adding one costs one, which is what puts
/// `PrintReadsxy` over the floor and `PrtReads` well under it.
#[test]
fn a_deletion_costs_four_and_an_insertion_one() {
    assert_eq!(SWAP_COST, 0);
    assert_eq!(SUBSTITUTION_COST, 2);
    assert_eq!(INSERTION_COST, 1);
    assert_eq!(DELETION_COST, 4);
    assert_eq!(HELP_SIMILARITY_FLOOR, 7);
    assert_eq!(distance("PrintReadsxy", "PrintReads"), Ok(8));
    assert_eq!(distance("PrtReads", "PrintReads"), Ok(2));
    let over = distance("PrintReadsxy", "PrintReads").expect("a distance");
    let under = distance("PrtReads", "PrintReads").expect("a distance");
    assert!(over >= HELP_SIMILARITY_FLOOR);
    assert!(under < HELP_SIMILARITY_FLOOR);
    // Which is exactly what the two messages show.
    let text = golden();
    let catalogue = catalogue(&text);
    assert!(!suggested_alternate_command(&catalogue, "PrintReadsxy")
        .expect("no refusal")
        .contains("Did you mean"));
    assert!(suggested_alternate_command(&catalogue, "PrtReads")
        .expect("no refusal")
        .contains("        PrintReads"));
}

/// Two matches at the best distance change the question and print both names with no separator.
#[test]
fn two_matches_run_together_on_one_line() {
    let two = vec!["CountReads".to_string(), "CountBases".to_string()];
    let message = suggested_alternate_command(&two, "ount").expect("no refusal");
    assert!(message.contains("Did you mean one of these?"), "{message}");
    assert!(
        message.contains("        CountReads        CountBases"),
        "{message}"
    );
    assert!(message.ends_with("CountBases"), "{message}");
    // One match asks the other question.
    let one = vec!["CountReads".to_string(), "VariantFiltration".to_string()];
    let message = suggested_alternate_command(&one, "ountR").expect("no refusal");
    assert!(message.contains("Did you mean this?"), "{message}");
}

/// When every tool scores zero the suggestion is suppressed rather than every tool being listed,
/// which a one-tool catalogue the command prefixes reaches as surely as the empty command does.
#[test]
fn every_tool_scoring_zero_suppresses_the_suggestion() {
    let text = golden();
    let catalogue = catalogue(&text);
    let empty = suggested_alternate_command(&catalogue, "").expect("no refusal");
    assert_eq!(empty, "'' is not a valid command.\n");
    assert!(catalogue.iter().all(|name| distance("", name) == Ok(0)));
    let one = vec!["PrintReads".to_string()];
    assert_eq!(
        suggested_alternate_command(&one, "Print").expect("no refusal"),
        "'Print' is not a valid command.\n"
    );
    // Two tools of which only one scores zero is not that case.
    let two = vec!["PrintReads".to_string(), "VariantFiltration".to_string()];
    assert!(suggested_alternate_command(&two, "Print")
        .expect("no refusal")
        .contains("        PrintReads"));
}

/// A name that resolves is a refusal, the search being reached only once resolution has failed.
#[test]
fn a_name_that_resolves_is_refused() {
    let text = golden();
    let catalogue = catalogue(&text);
    assert_eq!(
        suggested_alternate_command(&catalogue, "PrintReads"),
        Err("Command matches: PrintReads".to_string())
    );
    let error = field(&text, "error", "name-matches-a-tool").expect("the refusal");
    assert_eq!(
        error,
        "java.lang.RuntimeException:Command matches: PrintReads"
    );
}

/// The message always opens on the same line and always ends it, suggestion or not.
#[test]
fn the_message_always_opens_on_the_same_line() {
    let text = golden();
    for (case, command) in COMMANDS {
        if *case == "deprecated-short-circuits" {
            continue;
        }
        let message = field(&text, "message", case).unwrap_or_else(|| panic!("{case}"));
        assert!(
            message.starts_with(&format!("'{command}' is not a valid command.\n")),
            "{case}: {message}"
        );
    }
}

/// The distance itself, whose swap is a transposition of two adjacent characters.
#[test]
fn the_distance_prices_a_transposition_apart() {
    // With the dispatcher's weights a transposition is free and a substitution is not.
    assert_eq!(levenshtein_distance("ab", "ba", 0, 2, 1, 4), 0);
    assert_eq!(levenshtein_distance("ab", "cb", 0, 2, 1, 4), 2);
    // The plain distance counts each edit as one.
    assert_eq!(levenshtein_distance("kitten", "sitting", 1, 1, 1, 1), 3);
    assert_eq!(levenshtein_distance("same", "same", 0, 2, 1, 4), 0);
    // The first row's LAST cell is left at zero rather than seeded, so an empty first string
    // answers zero whatever the second is. An empty second one costs a deletion per character.
    assert_eq!(levenshtein_distance("", "abc", 0, 2, 1, 4), 0);
    assert_eq!(levenshtein_distance("abc", "", 0, 2, 1, 4), 12);
}
