//! Conformance for the two names of a mutex target against GATK 4.6.2.0.
//!
//! Golden from `tools/argument-conformance/MutexTargetsDump.java`. A mutex target is declared by
//! one name and printed by another: `getMutexTargetList()` and the annotation's own `mutex()` hold
//! the LONG name, and the usage prints the target definition's FIELD name.
//!
//! # What this suite is for
//!
//!  * **the map between the two names, for every argument that has a mutex**;
//!  * **the sentence the usage builds from the second**;
//!  * **and the declarations agreeing about the first, which is what made the difference
//!    invisible until a walker's usage was composed.**

use gatk_corpus as corpus;
use gatk_tools::plugin_ownership::mutex_field_name;
use gatk_tools::tool_declarations::{declarations, COUNTREADS};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/mutex_target_names.txt.gz"),
    )
}

fn rows(text: &str, kind: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|line| line.starts_with(&format!("{kind}\t")))
        .map(|line| line.split('\t').skip(1).map(str::to_string).collect())
        .collect()
}

/// Every row's field name is the port's, and every row's long name is the declarations'.
#[test]
fn the_two_names_are_the_goldens() {
    let text = golden();
    let recorded = rows(&text, "mutex");
    assert!(!recorded.is_empty());
    for row in &recorded {
        let (tool, argument, target, field, source) = (&row[0], &row[1], &row[2], &row[3], &row[4]);
        // The annotation named the target by its long name, which is what the declarations carry.
        assert_eq!(source, "annotated", "{tool}/{argument}");
        let list = declarations(tool).unwrap_or_else(|| panic!("{tool}"));
        let declared = list
            .iter()
            .find(|declaration| declaration.long_name == *argument)
            .unwrap_or_else(|| panic!("{tool}/{argument}"));
        assert!(
            declared.mutex.contains(&target.as_str()),
            "{tool}/{argument}"
        );
        // And the name the usage prints is the target's field, which is what the port holds.
        assert_eq!(mutex_field_name(target), Some(field.as_str()), "{target}");
        assert_ne!(target, field, "{target}");
    }
}

/// The sentence, which is the line a usage carries.
#[test]
fn the_sentence_is_built_from_the_field_name() {
    let text = golden();
    let sentences = rows(&text, "sentence");
    assert!(!sentences.is_empty());
    let composed = gatk_cli::composed_usage("CountReads").expect("the composition");
    let flattened = composed.split_whitespace().collect::<Vec<&str>>().join(" ");
    for row in &sentences {
        if row[0] != "CountReads" {
            continue;
        }
        assert!(flattened.contains(&row[1]), "{}", row[1]);
    }
    // Four mutex arguments, all of them a read filter's, and none of the tool's own.
    let mutexed: Vec<&str> = COUNTREADS
        .iter()
        .filter(|declaration| !declaration.mutex.is_empty())
        .map(|declaration| declaration.long_name)
        .collect();
    assert_eq!(mutexed.len(), 4);
    for name in mutexed {
        assert!(mutex_field_name(name).is_some(), "{name}");
    }
}
