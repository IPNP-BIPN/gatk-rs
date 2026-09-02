//! Conformance for the main usage listing against GATK 4.6.2.0, line for line.
//!
//! Golden from `tools/readfilter-conformance/MainUsageDump.java`. `main-entry` measured the
//! listing's SHAPE, because the other three hundred and seventy-two lines are tool names and
//! summaries no golden carried; this carries them, which is what C.1 was waiting on (#818).
//!
//! The dump holds the DECLARATIONS as well as the rendering, and that is what makes this a test of
//! the layout rather than a copy of it: the port renders from the annotations and the comparison is
//! against what the reference printed.
//!
//! # What this suite is for
//!
//!  * **the padding being decided by the SIMPLE name and printed with the DISPLAY name**, so a
//!    Picard tool is tested at 45 characters and printed at 54;
//!  * **the long branch writing four spaces where the short one pads to 45**;
//!  * **a group description never being truncated**, `%-45s` being a minimum and not a width;
//!  * **the two orderings being different comparators**, groups by heading and tools by display
//!    name, so the Picard suffix takes part in the sort;
//!  * **every escape being a byte of the output**;
//!  * **and the same text reaching stdout for no arguments and stderr for a name that does not
//!    resolve.**
//!
//! While the suite is `golden-pending` the dump is named by `MAIN_USAGE_DUMP`.

use gatk_tools::main_usage::usage;

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t")
        .replace("\\n", "\n")
        .replace("\\\\", "\\")
}

/// The lines one stream carried, in index order.
fn lines(dump: &str, stream: &str) -> Vec<String> {
    let prefix = format!("line\t{stream}\t");
    let mut rows: Vec<(usize, String)> = dump
        .lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .map(|rest| {
            let (index, text) = rest.split_once('\t').expect("an index and a line");
            (index.parse().expect("an index"), unescape(text))
        })
        .collect();
    rows.sort_by_key(|(index, _)| *index);
    rows.into_iter().map(|(_, text)| text).collect()
}

fn count(dump: &str, stream: &str) -> usize {
    let prefix = format!("count\t{stream}\tlines\t");
    dump.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .expect("a count")
        .parse()
        .expect("a number")
}

#[test]
fn the_listing_is_the_reference_one_line_for_line() {
    let dump = match std::env::var("MAIN_USAGE_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by MAIN_USAGE_DUMP"),
        Err(_) => {
            println!(
                "skipped: the main-usage golden is still pending. Run the suite and point \
                 MAIN_USAGE_DUMP at tools/conformance/pending/main-usage.MainUsageDump.txt"
            );
            return;
        }
    };

    // `getCommandLineName()` is empty for GATK, which is why the first line has two spaces.
    let ours: Vec<String> = usage("").split('\n').map(str::to_string).collect();
    let theirs = lines(&dump, "no-arguments-out");
    assert_eq!(theirs.len(), count(&dump, "no-arguments-out"));

    for (index, (ours, theirs)) in ours.iter().zip(theirs.iter()).enumerate() {
        assert_eq!(ours, theirs, "line {index}");
    }
    assert_eq!(
        ours.len(),
        theirs.len(),
        "the listing is a different length"
    );

    // The other stream carries nothing on that path, and the same text on the other one: a port
    // that built the listing for one of them only would pass a suite that looked at one.
    assert_eq!(count(&dump, "no-arguments-err"), 0);
    assert_eq!(count(&dump, "unknown-name-out"), 0);
    assert_eq!(lines(&dump, "unknown-name-err"), theirs);
}

/// The declarations the port carries are the ones the dump measured, and nothing has been dropped.
#[test]
fn every_tool_and_group_the_reference_declares_is_here() {
    let dump = match std::env::var("MAIN_USAGE_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by MAIN_USAGE_DUMP"),
        Err(_) => return,
    };
    use gatk_tools::main_usage_catalogue::{GROUPS, TOOLS};

    let tools: Vec<&str> = dump
        .lines()
        .filter_map(|line| line.strip_prefix("tool\t"))
        .map(|rest| rest.split('\t').next().expect("a display name"))
        .collect();
    assert_eq!(tools.len(), TOOLS.len());
    for (row, tool) in tools.iter().zip(TOOLS) {
        assert_eq!(*row, tool.display_name);
    }

    let groups: Vec<&str> = dump
        .lines()
        .filter_map(|line| line.strip_prefix("group\t"))
        .map(|rest| rest.split('\t').next().expect("a class"))
        .collect();
    assert_eq!(groups.len(), GROUPS.len());
    for (row, group) in groups.iter().zip(GROUPS) {
        assert_eq!(*row, group.class);
    }
}
