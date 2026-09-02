//! `Main.printUsage`: the listing `gatk` answers with when it is given nothing.
//!
//! Three hundred and seventy-three lines, and none of them is plain text: every one carries the
//! ANSI escapes Barclay writes, and those are bytes of the output like any other. The data is
//! [`crate::main_usage_catalogue`], generated from the reference's own annotations; what is here
//! is the layout, and the layout is where a port goes wrong.
//!
//! # The padding is decided by a name the line does not show
//!
//! `getDisplaySummaryForTool` tests `toolClass.getSimpleName().length() >= 45`, and then prints
//! `toolDisplayName`, which is the simple name PLUS `" (Picard)"` for a Picard tool. So a Picard
//! tool whose simple name is 44 characters is padded to 45 while its printed name is 53, and a
//! GATK tool of the same printed length is not. The test and the value are two different strings.
//!
//! And the two branches are not the same layout with a different width: the long one writes the
//! name and then FOUR spaces, where the short one pads the name to 45 and writes none.
//!
//! # The group heading is two padded fields, and the second one is never truncated
//!
//! `String.format("%s%-48s %-45s%s\n", ...)` pads and does not cut: `%-45s` on a description
//! longer than 45 characters leaves it whole. A port that formatted with a width instead of a
//! minimum would truncate every long group description.
//!
//! # The order is two different comparators
//!
//! Groups sort by `getName()`, which is the heading rather than the class; tools sort by
//! `toolDisplayName`, which is the name WITH the Picard suffix, so `SortSam (Picard)` and
//! `SortSamSpark` interleave by that string and not by their class names.
//!
//! # The trailing lines
//!
//! The builder ends with a rule, and `println` then adds the newline that makes the last line
//! blank. The listing therefore ends with a rule, a `KNRM` on its own line and an empty line, and
//! a port that stopped at the rule is two lines short.
//!
//! Ported from `org.broadinstitute.hellbender.Main.printUsage` and
//! `Main.getDisplaySummaryForTool`.

use crate::main_usage_catalogue::{Group, Maturity, Tool, GROUPS, TOOLS};

/// The escapes `Main` writes, under the names it gives them.
pub const KNRM: &str = "\u{1b}[0m";
pub const RED: &str = "\u{1b}[31m";
pub const GREEN: &str = "\u{1b}[32m";
pub const CYAN: &str = "\u{1b}[36m";
pub const WHITE: &str = "\u{1b}[37m";
pub const BOLDRED: &str = "\u{1b}[1m\u{1b}[31m";

/// The rule between groups: eighty-six dashes.
pub const RULE: &str =
    "--------------------------------------------------------------------------------------";

/// `String.format("%-<width>s", value)`: pad on the right, never cut.
fn pad(value: &str, width: usize) -> String {
    if value.len() >= width {
        value.to_string()
    } else {
        format!("{value}{}", " ".repeat(width - value.len()))
    }
}

/// `getDisplaySummaryForTool`, which is one line and two branches.
pub fn tool_line(tool: &Tool) -> String {
    let summary = match tool.maturity {
        Maturity::Experimental => {
            format!("{RED}(EXPERIMENTAL Tool) {CYAN}{}", tool.one_line_summary)
        }
        Maturity::Beta => format!("{RED}(BETA Tool) {CYAN}{}", tool.one_line_summary),
        Maturity::Released => format!("{CYAN}{}", tool.one_line_summary),
    };
    // The test is on the SIMPLE name and the value printed is the DISPLAY name, which for a Picard
    // tool is nine characters longer.
    if tool.simple_name.len() >= 45 {
        format!("{GREEN}    {}    {summary}{KNRM}\n", tool.display_name)
    } else {
        format!("{GREEN}    {}{summary}{KNRM}\n", pad(tool.display_name, 45))
    }
}

/// `printUsage`'s group heading: the name with a colon, then the description.
pub fn group_heading(group: &Group) -> String {
    format!(
        "{WHITE}{RULE}\n{KNRM}{RED}{} {}{KNRM}\n",
        pad(&format!("{}:", group.name), 48),
        pad(group.description, 45)
    )
}

/// The whole listing, as `printUsage` builds it and `println` closes it.
///
/// `command_line_name` is `getCommandLineName()`, which GATK leaves EMPTY, so the first line reads
/// `USAGE:  <program name>` with two spaces where a name would have been.
pub fn usage(command_line_name: &str) -> String {
    let mut out = format!(
        "{BOLDRED}USAGE: {command_line_name} {GREEN}<program name>{BOLDRED} [-h]\n\n{KNRM}\
         {BOLDRED}Available Programs:\n{KNRM}"
    );

    // Groups by their NAME, which is the heading and not the class.
    let mut groups: Vec<&Group> = GROUPS.iter().collect();
    groups.sort_by(|a, b| a.name.cmp(b.name));

    for group in groups {
        out.push_str(&group_heading(group));
        // Tools by their DISPLAY name, so the Picard suffix takes part in the ordering.
        let mut tools: Vec<&Tool> = TOOLS.iter().filter(|t| t.group == group.class).collect();
        tools.sort_by(|a, b| a.display_name.cmp(b.display_name));
        for tool in tools {
            out.push_str(&tool_line(tool));
        }
        out.push('\n');
    }
    out.push_str(&format!("{WHITE}{RULE}\n{KNRM}"));
    // `println(builder.toString())`, whose newline is what makes the final line blank.
    out.push('\n');
    out
}
