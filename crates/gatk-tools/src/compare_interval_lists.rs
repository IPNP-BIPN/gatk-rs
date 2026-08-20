//! `CompareIntervalLists`, ported from
//! `org.broadinstitute.hellbender.tools.CompareIntervalLists` and `IntervalUtils.equateIntervals`
//! (GATK 4.6.2.0).
//!
//! Two interval files, each sorted and `ALL`-merged on its own, walked against one another until
//! one of them runs out.
//!
//! # The walk is not symmetric
//!
//! ```java
//! while ( ! master.isEmpty() ) {
//!     final GenomeLoc masterHead = master.pop();
//!     final GenomeLoc testHead = test.pop();
//! ```
//!
//! The loop condition reads the master and the body pops both. A master that outlasts the test
//! therefore pops an empty list and dies with a `NoSuchElementException`, while a test that
//! outlasts the master falls out of the loop and is reported as a message. Swapping the two files
//! is not a no-op: the same pair can be equal one way and a crash the other.
//!
//! # A wider test interval is equal
//!
//! Only the master's remainder is pushed back onto the master:
//!
//! ```java
//! reverse(masterHead.subtract(testHead)).forEach(master::push);
//! ```
//!
//! `subtract` of a wider test leaves nothing, so the test's overhang is never examined and never
//! compared with anything. `chr1:10-20` as master against `chr1:1-100` as test is reported equal.
//!
//! # The remainder is pushed in reverse
//!
//! `subtract` returns the AFTER piece before the BEFORE piece, and `reverse` sorts them descending
//! so that pushing leaves the earlier one on top. A test interval strictly inside the master
//! therefore leaves the master's two remaining pieces in coordinate order, which is what makes the
//! next comparison meaningful.

use gatk_engine::interval::SimpleInterval;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK CompareIntervalLists";

/// What the comparison answers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Comparison {
    /// `null` from `equateIntervals`, which the tool prints as `Intervals are equal`.
    Equal,
    /// A difference, which the tool wraps in a `UserException`.
    Different(String),
    /// `test.pop()` on an empty list, which is a `NoSuchElementException` with no message.
    TestExhausted,
}

impl Comparison {
    /// The line the tool leaves behind: its own output, or the exception and its message.
    pub fn line(&self) -> String {
        match self {
            Comparison::Equal => "Intervals are equal".to_string(),
            Comparison::Different(difference) => format!(
                "org.broadinstitute.hellbender.exceptions.UserException:Intervals are not equal: \n{difference}"
            ),
            Comparison::TestExhausted => "java.util.NoSuchElementException:null".to_string(),
        }
    }
}

/// `GenomeLoc.toString()` for a mapped interval, which is what the difference messages carry.
fn display(interval: &SimpleInterval) -> String {
    format!("{}:{}-{}", interval.contig, interval.start, interval.end)
}

fn overlaps(left: &SimpleInterval, right: &SimpleInterval) -> bool {
    left.contig == right.contig && left.start <= right.end && right.start <= left.end
}

fn contains(outer: &SimpleInterval, inner: &SimpleInterval) -> bool {
    outer.contig == inner.contig && outer.start <= inner.start && outer.end >= inner.end
}

/// `GenomeLoc.subtract`, which is only ever called on overlapping intervals here.
///
/// The order of the two pieces is the reference's: the AFTER piece first, then the BEFORE piece.
fn subtract(this: &SimpleInterval, that: &SimpleInterval) -> Vec<SimpleInterval> {
    if this == that {
        return Vec::new();
    }
    if contains(this, that) {
        let mut pieces = Vec::new();
        if this.end - (that.end + 1) >= 0 {
            pieces.push(SimpleInterval {
                contig: this.contig.clone(),
                start: that.end + 1,
                end: this.end,
            });
        }
        if (that.start - 1) - this.start >= 0 {
            pieces.push(SimpleInterval {
                contig: this.contig.clone(),
                start: this.start,
                end: that.start - 1,
            });
        }
        return pieces;
    }
    if contains(that, this) {
        return Vec::new();
    }
    let piece = if that.start < this.start {
        SimpleInterval {
            contig: this.contig.clone(),
            start: that.end + 1,
            end: this.end,
        }
    } else {
        SimpleInterval {
            contig: this.contig.clone(),
            start: this.start,
            end: that.start - 1,
        }
    };
    vec![piece]
}

/// `IntervalUtils.equateIntervals(master, test)`.
///
/// Both lists are expected already sorted and `ALL`-merged, which is what the tool's own
/// `getGenomeLocs` does to each file before it calls this.
pub fn equate_intervals(master: &[SimpleInterval], test: &[SimpleInterval]) -> Comparison {
    let mut master: std::collections::VecDeque<SimpleInterval> = master.iter().cloned().collect();
    let mut test: std::collections::VecDeque<SimpleInterval> = test.iter().cloned().collect();

    while let Some(master_head) = master.pop_front() {
        let Some(test_head) = test.pop_front() else {
            return Comparison::TestExhausted;
        };
        if !overlaps(&test_head, &master_head) {
            return Comparison::Different(format!(
                "Incompatible locs detected masterHead={}, testHead={}",
                display(&master_head),
                display(&test_head)
            ));
        }
        // `reverse` sorts the remainder descending, and each `push` puts one on the front, so the
        // earliest piece ends up on top.
        let mut remainder = subtract(&master_head, &test_head);
        remainder.sort_by(|left, right| {
            (&right.contig, right.start, right.end).cmp(&(&left.contig, left.start, left.end))
        });
        for piece in remainder {
            master.push_front(piece);
        }
    }

    match test.front() {
        None => Comparison::Equal,
        Some(first) => Comparison::Different(format!(
            "Remaining elements found in test: first={}",
            display(first)
        )),
    }
}
