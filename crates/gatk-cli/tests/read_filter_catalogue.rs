//! Conformance for the read filter catalogue and each tool's defaults against GATK 4.6.2.0.
//!
//! Golden from `tools/argument-conformance/ReadFilterCatalogueDump.java`: the filter library once,
//! then the nine declared tools with the descriptor they built, their defaults, and the possible
//! values of the two arguments that print them.
//!
//! # What this suite is for
//!
//!  * **the catalogue, which is the library and not the ownership table**;
//!  * **the defaults, which are per tool and count as SELECTED**;
//!  * **a tool that is no walker having no descriptor at all, which is not an empty list**;
//!  * **and the two arguments the usage reads those two sets through.**

use gatk_corpus as corpus;
use gatk_tools::plugin_ownership::{self, CATALOGUE, DESCRIPTOR, OWNERSHIP};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/read_filter_catalogue.txt.gz"),
    )
}

fn rows(text: &str, kind: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|line| line.starts_with(&format!("{kind}\t")))
        .map(|line| line.split('\t').skip(1).map(str::to_string).collect())
        .collect()
}

/// The catalogue is the golden's, name for name and in its order.
#[test]
fn the_catalogue_is_the_goldens() {
    let text = golden();
    let row = &rows(&text, "catalogue")[0];
    assert_eq!(row[0].parse::<usize>().expect("a count"), CATALOGUE.len());
    let names: Vec<&str> = row[1].split(' ').collect();
    assert_eq!(names, CATALOGUE.to_vec());
    // It is the library rather than the ownership table: every owner is in it, and most of it
    // declares no argument at all.
    for entry in OWNERSHIP.iter() {
        assert!(CATALOGUE.contains(&entry.owner), "{}", entry.owner);
    }
    let owners = OWNERSHIP
        .iter()
        .map(|entry| entry.owner)
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(owners.len(), 20);
    assert!(CATALOGUE.len() > owners.len());
}

/// The defaults are the golden's, per tool, and a tool with no descriptor has none.
#[test]
fn the_defaults_are_the_goldens() {
    let text = golden();
    for row in rows(&text, "descriptor") {
        let tool = &row[0];
        match row[1].as_str() {
            "none" => assert_eq!(plugin_ownership::default_filters(tool), None, "{tool}"),
            name => {
                assert_eq!(name, DESCRIPTOR, "{tool}");
                assert!(plugin_ownership::default_filters(tool).is_some(), "{tool}");
            }
        }
    }
    for row in rows(&text, "defaults") {
        let expected: Vec<&str> = row[1].split(' ').collect();
        let ported =
            plugin_ownership::default_filters(&row[0]).unwrap_or_else(|| panic!("{}", row[0]));
        assert_eq!(ported.to_vec(), expected, "{}", row[0]);
    }
    // Five walkers and four tools that are not, which is the whole of the declared set.
    assert_eq!(rows(&text, "defaults").len(), 5);
    assert_eq!(rows(&text, "descriptor").len(), 9);
}

/// The two arguments the usage reads the sets through, whose values are the descriptor's.
#[test]
fn the_two_arguments_carry_the_two_sets() {
    let text = golden();
    for row in rows(&text, "allowed") {
        let (tool, argument, values) = (&row[0], &row[1], &row[2]);
        let expected: Vec<&str> = values.split(' ').collect();
        let ported: Vec<&str> = match argument.as_str() {
            "read-filter" => CATALOGUE.to_vec(),
            "disable-read-filter" => plugin_ownership::default_filters(tool)
                .expect("a walker's defaults")
                .to_vec(),
            other => panic!("{other}"),
        };
        assert_eq!(ported, expected, "{tool}/{argument}");
    }
}
