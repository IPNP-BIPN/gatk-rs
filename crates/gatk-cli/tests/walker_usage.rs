//! Conformance for a walker's composed usage against GATK 4.6.2.0.
//!
//! `usage-text` measures the whole text `CountReads -h` prints, conditional blocks and all. The
//! port composes it from the declarations, the ownership table (`plugin-argument-ownership`), the
//! catalogue and defaults (`read-filter-catalogue`) and the mutex field names
//! (`mutex-target-names`), and it is the reference's text byte for byte.
//!
//! # What this suite is for
//!
//!  * **the conditional blocks: one per read filter that declares an argument, in the ownership
//!    table's order, with the arguments under their owner**;
//!  * **the two arguments the descriptor answers for: `--read-filter` prints the whole catalogue
//!    and `--disable-read-filter` prints the TOOL'S OWN defaults**;
//!  * **the wrapping, which drops a line that would hold nothing but the indent**;
//!  * **and the mutex sentence, which names the target definition's FIELD rather than the long
//!    name the declarations carry.**

use gatk_corpus as corpus;

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(character) = chars.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn golden(tool: &str) -> String {
    let text = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/usage_text.txt.gz"),
    );
    let prefix = format!("usage\t{tool}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
        .unwrap_or_else(|| panic!("{tool}"))
}

/// A walker's usage, byte for byte.
///
/// It was four lines short until the mutex field names were measured: the sentence names the
/// target definition's FIELD, and both the declarations golden and the annotation hold the long
/// name it resolves to. With `mutex-target-names` frozen the composition is the reference's.
#[test]
fn a_walkers_usage_is_the_goldens() {
    for tool in [
        "CountReads",
        "PrintBGZFBlockInformation",
        "IndexFeatureFile",
        "GatherVcfsCloud",
    ] {
        let expected = golden(tool);
        let produced = gatk_cli::composed_usage(tool).unwrap_or_else(|| panic!("{tool}"));
        assert_eq!(produced, expected, "{tool}");
        // Which is what the dispatcher answers `-h` with now, for a walker as for the rest.
        assert_eq!(
            gatk_cli::tool_usage(tool).as_deref(),
            Some(expected.as_str()),
            "{tool}"
        );
    }
    // The one that carries the sentence, so the field name is asserted and not merely reached.
    let usage = gatk_cli::composed_usage("CountReads").expect("the composition");
    assert!(usage.contains("argument(s) maxAmbiguousBaseFraction"));
    assert!(!usage.contains("argument(s) ambig-filter-frac"));
    assert_eq!(
        gatk_tools::plugin_ownership::mutex_field_name("ambig-filter-frac"),
        Some("maxAmbiguousBaseFraction")
    );
    // An argument whose field is named after it keeps its long name, which is why the sentence
    // looked right everywhere it had been compared before.
    assert_eq!(
        gatk_tools::plugin_ownership::mutex_field_name("intervals"),
        None
    );
}

/// The blocks themselves: their heading, their order, and what is under each.
#[test]
fn the_conditional_blocks_are_one_per_filter_that_declares_an_argument() {
    let produced = gatk_cli::composed_usage("CountReads").expect("the composition");
    assert!(produced.contains("Conditional Arguments for readFilter:"));
    let headings: Vec<&str> = produced
        .lines()
        .filter(|line| line.starts_with("Valid only if "))
        .collect();
    // Twenty-eight arguments over twenty filters, which is what the ownership table holds.
    let owners: Vec<&str> = {
        let mut seen: Vec<&str> = Vec::new();
        for entry in gatk_tools::plugin_ownership::OWNERSHIP.iter() {
            if !seen.contains(&entry.owner) {
                seen.push(entry.owner);
            }
        }
        seen
    };
    assert_eq!(headings.len(), owners.len());
    for (heading, owner) in headings.iter().zip(&owners) {
        assert_eq!(*heading, format!("Valid only if \"{owner}\" is specified:"));
    }
    // A filter that declares nothing has no block, so the catalogue is longer than the blocks.
    assert!(gatk_tools::plugin_ownership::CATALOGUE.len() > owners.len());
}

/// The two arguments the descriptor answers for, whose possible values are not an enum's.
#[test]
fn the_descriptor_answers_for_two_arguments() {
    let produced = gatk_cli::composed_usage("CountReads").expect("the composition");
    // `--disable-read-filter` lists the TOOL's defaults, which is one filter for a read walker.
    assert!(produced.contains("Possible values: {WellformedReadFilter} "));
    // `--read-filter` lists the whole catalogue, whose first and last names bracket it.
    assert!(produced.contains("Possible values: {AlignmentAgreesWithHeaderReadFilter,"));
    assert!(produced.contains("WellformedReadFilter} "));
    assert_eq!(
        gatk_tools::plugin_ownership::default_filters("CountReads"),
        Some(&["WellformedReadFilter"][..])
    );
    // A tool that is no walker builds no descriptor, so it has neither argument to answer for.
    assert_eq!(
        gatk_tools::plugin_ownership::default_filters("IndexFeatureFile"),
        None
    );
    let plain = gatk_cli::composed_usage("IndexFeatureFile").expect("the composition");
    assert!(!plain.contains("Conditional Arguments"));
    assert_eq!(gatk_cli::tool_usage("IndexFeatureFile"), Some(plain));
}
