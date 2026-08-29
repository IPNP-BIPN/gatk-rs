//! Conformance for a walker's composed usage against GATK 4.6.2.0.
//!
//! `usage-text` measures the whole text `CountReads -h` prints, conditional blocks and all. The
//! port composes it from the declarations, the ownership table (`plugin-argument-ownership`) and
//! the catalogue and defaults (`read-filter-catalogue`), and this suite says how close that is:
//! two hundred and ninety-three of the two hundred and ninety-seven lines, exactly.
//!
//! # What this suite is for
//!
//!  * **the conditional blocks: one per read filter that declares an argument, in the ownership
//!    table's order, with the arguments under their owner**;
//!  * **the two arguments the descriptor answers for: `--read-filter` prints the whole catalogue
//!    and `--disable-read-filter` prints the TOOL'S OWN defaults**;
//!  * **the wrapping, which drops a line that would hold nothing but the indent**;
//!  * **and the one gap, counted rather than described: four lines whose mutex sentence names the
//!    other argument by its FIELD name, which no golden carries yet.**

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

/// The four lines the port cannot write yet, and why.
///
/// `Cannot be used in conjunction with argument(s) maxAmbiguousBaseFraction` names the FIELD the
/// annotation named, and the declarations golden measured `getMutexTargetList()`, which is the
/// long name that field resolves to. Four of the twenty-eight controlled arguments carry the
/// sentence, and each contributes one line.
const UNMEASURED_MUTEX_FIELDS: [&str; 4] = [
    "maxAmbiguousBaseFraction",
    "maxAmbiguousBases",
    "maximumSoftClippedRatio",
    "maximumLeadingTrailingSoftClippedRatio",
];

/// A walker's usage, line for line, with the gap counted.
#[test]
fn a_walkers_usage_is_the_goldens_but_for_the_mutex_fields() {
    let expected = golden("CountReads");
    let produced = gatk_cli::composed_usage("CountReads").expect("the composition");
    let expected_lines: Vec<&str> = expected.lines().collect();
    let produced_lines: Vec<&str> = produced.lines().collect();
    assert_eq!(expected_lines.len(), 297);
    assert_eq!(produced_lines.len(), expected_lines.len());

    let mut differing = 0;
    for (index, (left, right)) in expected_lines.iter().zip(&produced_lines).enumerate() {
        if left == right {
            continue;
        }
        differing += 1;
        // Every difference is the same difference: the reference names a field, the port names the
        // long name that field resolves to.
        let field = UNMEASURED_MUTEX_FIELDS
            .iter()
            .find(|field| left.contains(**field))
            .unwrap_or_else(|| panic!("line {index}: {left:?} vs {right:?}"));
        assert!(left.contains(*field), "line {index}");
        assert!(!right.contains(*field), "line {index}");
        assert!(right.contains("argument(s) "), "line {index}: {right:?}");
    }
    assert_eq!(differing, UNMEASURED_MUTEX_FIELDS.len());
    // Which is why the dispatcher does not answer `-h` with it.
    assert!(gatk_cli::tool_usage("CountReads").is_none());
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
