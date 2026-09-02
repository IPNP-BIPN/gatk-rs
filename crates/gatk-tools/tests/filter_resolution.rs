//! Conformance for `GATKReadFilterPluginDescriptor`'s resolution against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/FilterResolutionDump.java`. The runner read
//! `--read-filter` and `--disable-tool-default-read-filters` and ignored the other two arguments
//! the descriptor owns, so `--disable-read-filter` and `--inverted-read-filter` changed nothing in
//! the port and change the filters, and therefore the count, in the reference.
//!
//! # What this suite is for
//!
//!  * **the order being defaults, then enabled, then inverted**, with disabling applied to the
//!    defaults only and before anything is added;
//!  * **an enabled filter already among the defaults not being added twice**;
//!  * **an inverted filter being a different filter**, appended rather than replacing anything;
//!  * **inverting a tool default being a refusal of its own** unless the defaults were disabled;
//!  * **enabling or disabling the same filter twice being two different refusals**;
//!  * **an unknown name being one refusal on `--disable-read-filter` and a different one, with a
//!    different exception class, on the other two**;
//!  * **and the "enabled and inverted" refusal listing the empty set**, because the reference
//!    formats the wrong variable. That is its behaviour, and the port reproduces it.
//!
//! While the suite is `golden-pending` the dump is named by `FILTER_RESOLUTION_DUMP`.

use gatk_tools::filter_resolution::resolve;
use gatk_tools::plugin_ownership::CATALOGUE;

/// The defaults the harness gave its imaginary tool, in the order it declares them.
const DEFAULTS: &[&str] = &[
    "MappedReadFilter",
    "NotDuplicateReadFilter",
    "PrimaryLineReadFilter",
];

/// One command line: enabled, disabled, inverted, and whether the defaults were dropped.
type Case = (
    &'static str,
    &'static [&'static str],
    &'static [&'static str],
    &'static [&'static str],
    bool,
);

const CASES: &[Case] = &[
    ("defaults-only", &[], &[], &[], false),
    ("one-enabled", &["GoodCigarReadFilter"], &[], &[], false),
    (
        "two-enabled",
        &["GoodCigarReadFilter", "FirstOfPairReadFilter"],
        &[],
        &[],
        false,
    ),
    (
        "enabled-is-a-default",
        &["MappedReadFilter"],
        &[],
        &[],
        false,
    ),
    ("disable-a-default", &[], &["MappedReadFilter"], &[], false),
    (
        "disable-two-defaults",
        &[],
        &["MappedReadFilter", "PrimaryLineReadFilter"],
        &[],
        false,
    ),
    (
        "disable-a-non-default",
        &[],
        &["GoodCigarReadFilter"],
        &[],
        false,
    ),
    ("disable-unknown", &[], &["NoSuchReadFilter"], &[], false),
    ("enabled-unknown", &["NoSuchReadFilter"], &[], &[], false),
    ("inverted-unknown", &[], &[], &["NoSuchReadFilter"], false),
    (
        "invert-a-non-default",
        &[],
        &[],
        &["GoodCigarReadFilter"],
        false,
    ),
    ("invert-a-default", &[], &[], &["MappedReadFilter"], false),
    (
        "enabled-and-inverted",
        &["GoodCigarReadFilter"],
        &[],
        &["GoodCigarReadFilter"],
        false,
    ),
    ("no-defaults", &[], &[], &[], true),
    (
        "no-defaults-one-enabled",
        &["GoodCigarReadFilter"],
        &[],
        &[],
        true,
    ),
    (
        "no-defaults-and-disable",
        &[],
        &["MappedReadFilter"],
        &[],
        true,
    ),
    (
        "enabled-twice",
        &["GoodCigarReadFilter", "GoodCigarReadFilter"],
        &[],
        &[],
        false,
    ),
    (
        "disabled-twice",
        &[],
        &["MappedReadFilter", "MappedReadFilter"],
        &[],
        false,
    ),
    (
        "enabled-and-disabled",
        &["GoodCigarReadFilter"],
        &["GoodCigarReadFilter"],
        &[],
        false,
    ),
];

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_string()).collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t")
        .replace("\\n", "\n")
        .replace("\\\\", "\\")
}

fn field(dump: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    dump.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

#[test]
fn every_command_line_resolves_to_the_reference_filters() {
    let dump = match std::env::var("FILTER_RESOLUTION_DUMP") {
        Ok(path) => {
            std::fs::read_to_string(path).expect("the dump named by FILTER_RESOLUTION_DUMP")
        }
        Err(_) => {
            println!(
                "skipped: the filter-resolution golden is still pending. Run the suite and point \
                 FILTER_RESOLUTION_DUMP at \
                 tools/conformance/pending/filter-resolution.FilterResolutionDump.txt"
            );
            return;
        }
    };

    for (case, enabled, disabled, inverted, disable_defaults) in CASES {
        let answer = resolve(
            DEFAULTS,
            &CATALOGUE,
            &strings(enabled),
            &strings(disabled),
            &strings(inverted),
            *disable_defaults,
        );
        match answer {
            Ok(filters) => {
                // A negated filter is a `ReadFilterNegate` and prints as one: the list carries the
                // wrapper's class name, not the wrapped filter's.
                let rendered = if filters.is_empty() {
                    "(empty)".to_string()
                } else {
                    filters
                        .iter()
                        .map(|filter| {
                            if filter.negated {
                                "ReadFilterNegate".to_string()
                            } else {
                                filter.name.clone()
                            }
                        })
                        .collect::<Vec<_>>()
                        .join(" ")
                };
                assert_eq!(
                    rendered,
                    field(&dump, "filters", case)
                        .unwrap_or_else(|| panic!("{case}: the golden refused, the port did not")),
                    "{case}"
                );
            }
            Err(error) => {
                assert_eq!(
                    format!("{}: {}", error.java_class(), error.message()),
                    field(&dump, "error", case)
                        .unwrap_or_else(|| panic!("{case}: the port refused, the golden did not")),
                    "{case}"
                );
            }
        }
    }

    let cases: std::collections::BTreeSet<&str> = dump
        .lines()
        .filter_map(|line| line.split('\t').nth(1))
        .collect();
    assert_eq!(
        cases.len(),
        CASES.len(),
        "the dump carries a case this test does not"
    );
}
