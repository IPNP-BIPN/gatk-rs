//! `CountVariants`, ported from `org.broadinstitute.hellbender.tools.walkers.CountVariants`
//! (GATK 4.6.2.0).
//!
//! The tool's own body is three lines: a counter, an increment in `apply`, and `out.print(count)`
//! in `onTraversalSuccess`. Everything there is to be identical about therefore sits on either side
//! of it: the traversal that decides what `apply` is called on, and the output collection that
//! decides what reaches a file.
//!
//! # The count reaches no stream without `-O`
//!
//! The class documentation says "The tool prints the count to standard output (and can optionally
//! write it to a file)". It does not. `OptionalTextOutputArgumentCollection.print` is
//!
//! ```java
//! public void print(Object value) {
//!     if (output != null) { Files.write(output.toPath(), value.toString().getBytes()); }
//! }
//! ```
//!
//! so with no `-O` the count is written nowhere at all. It survives as the traversal's return
//! value, a `java.lang.Long`, and in the log line the engine prints on the way out. The golden
//! records the return value for exactly that reason: it is the only place a run without `-O` puts
//! the number.
//!
//! # The file has no trailing newline, and it is truncated rather than appended
//!
//! `onTraversalSuccess` calls `print` and not `println`, so a count of 5 is one byte. `Files.write`
//! with no options is `CREATE`, `TRUNCATE_EXISTING`, `WRITE`, so the same one byte replaces a
//! ten-byte file rather than following it.
//!
//! # Every row counts
//!
//! There is no variant filter on this walker, so a `FILTER` column, a symbolic allele, a
//! multi-allelic site and a second record at a position already seen each count once. The count is
//! rows, not variant alleles, which is what the tool's own summary says.
//!
//! # A record is selected by its whole span, not by its position
//!
//! Traversal by `-L` is a Tribble query, and the span a query matches is the decoded record's:
//! `END` when the record carries one, otherwise `start + len(REF) - 1`. A record at `chr1:100` with
//! `END=400` is counted by `-L chr1:300-310`, which its `POS` never reaches, and a deletion whose
//! `REF` is ten bases long is counted by an interval over its sixth base.
//!
//! A record spanning two intervals is still counted once: [`gatk_engine::variant_source`] carries
//! the one-interval memory of `FeatureIntervalIterator.featureIsNovel`, which drops a feature that
//! overlaps the previous interval.
//!
//! # The two refusals come from opposite ends of the run
//!
//! `-L` against a file with no index is refused by `setIntervalsForTraversal` before a record is
//! read. `-O` onto a path that cannot be written is refused after the whole traversal has already
//! run, and its message is **the path alone**: `new UserException.CouldNotCreateOutputFile(output
//! .toString(), e)` selects the `(String message, Exception e)` overload, so the filename is the
//! whole message and the `IOException` behind it is dropped from the text. None of the "Couldn't
//! write file %s because %s" overloads is the one this call site reaches.

use std::path::Path;

use gatk_engine::interval::SimpleInterval;
use gatk_engine::variant_source::{traverse, Located};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK CountVariants";

/// What this tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CountVariantsError {
    /// `-L` against an input with no index, thrown before any record is read.
    IntervalsWithoutRandomAccess { path: String },
    /// `-O` onto a path that cannot be written, thrown after the traversal.
    CouldNotCreateOutputFile { path: String },
}

impl CountVariantsError {
    /// The Java class the reference throws.
    pub fn class(&self) -> &'static str {
        match self {
            CountVariantsError::IntervalsWithoutRandomAccess { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            CountVariantsError::CouldNotCreateOutputFile { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$CouldNotCreateOutputFile"
            }
        }
    }

    /// The message it carries, which for the output file is the path and nothing else.
    pub fn message(&self) -> String {
        match self {
            CountVariantsError::IntervalsWithoutRandomAccess { path } => format!(
                "Input {path} must support random access to enable traversal by intervals. \
                 If it's a file, please index it using the bundled tool IndexFeatureFile"
            ),
            CountVariantsError::CouldNotCreateOutputFile { path } => path.clone(),
        }
    }
}

/// `CountVariants.apply`, counted over the traversal: one per record handed to it.
///
/// `intervals` is the merged list the argument layer produces, and `indexed` is whether the input
/// supports random access, which only matters when there are intervals to traverse by.
pub fn count<T: Located>(
    features: &[T],
    intervals: Option<&[SimpleInterval]>,
    indexed: bool,
    path: &str,
) -> Result<i64, CountVariantsError> {
    let restricted = gatk_engine::variant_source::intervals_for_traversal(intervals);
    if restricted.is_some() && !indexed {
        return Err(CountVariantsError::IntervalsWithoutRandomAccess {
            path: path.to_string(),
        });
    }
    Ok(traverse(features, intervals).len() as i64)
}

/// `OptionalTextOutputArgumentCollection.print(count)`: the bytes the file holds afterwards.
///
/// `Long.toString`, with no trailing newline. A run with no `-O` writes nothing, which is what the
/// `None` return stands for rather than an empty file.
pub fn output_bytes(count: i64, has_output: bool) -> Option<Vec<u8>> {
    if has_output {
        Some(count.to_string().into_bytes())
    } else {
        None
    }
}

/// The same, written where `-O` points, truncating whatever was there.
///
/// The error carries the path as its whole message, as the reference's overload does.
pub fn write_output(output: Option<&Path>, count: i64) -> Result<(), CountVariantsError> {
    let Some(output) = output else {
        return Ok(());
    };
    std::fs::write(output, count.to_string().as_bytes()).map_err(|_| {
        CountVariantsError::CouldNotCreateOutputFile {
            path: output.display().to_string(),
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Variant {
        contig: &'static str,
        start: i32,
        stop: i32,
    }

    impl Located for Variant {
        fn contig(&self) -> &str {
            self.contig
        }
        fn start(&self) -> i32 {
            self.start
        }
        fn stop(&self) -> i32 {
            self.stop
        }
    }

    /// The `spanning` fixture: an `END` block from 100 to 400, and a ten-base deletion at 600.
    const SPANNING: &[Variant] = &[
        Variant {
            contig: "chr1",
            start: 100,
            stop: 400,
        },
        Variant {
            contig: "chr1",
            start: 600,
            stop: 609,
        },
    ];

    fn interval(contig: &str, start: i32, end: i32) -> SimpleInterval {
        SimpleInterval::new(contig, start, end).expect("a valid interval")
    }

    #[test]
    fn the_count_is_rows_and_the_file_has_no_trailing_newline() {
        assert_eq!(output_bytes(5, true), Some(b"5".to_vec()));
        assert_eq!(output_bytes(5, true).unwrap().len(), 1);
        assert_eq!(output_bytes(0, true), Some(b"0".to_vec()));
    }

    #[test]
    fn without_an_output_argument_the_count_is_written_nowhere() {
        assert_eq!(output_bytes(5, false), None);
    }

    #[test]
    fn a_record_is_matched_by_its_span_and_not_by_its_position() {
        let by_end = [interval("chr1", 300, 310)];
        assert_eq!(
            count(SPANNING, Some(&by_end), true, "spanning.vcf").unwrap(),
            1
        );

        let by_ref_length = [interval("chr1", 605, 606)];
        assert_eq!(
            count(SPANNING, Some(&by_ref_length), true, "spanning.vcf").unwrap(),
            1
        );

        let between_them = [interval("chr1", 500, 510)];
        assert_eq!(
            count(SPANNING, Some(&between_them), true, "spanning.vcf").unwrap(),
            0
        );
    }

    #[test]
    fn one_record_over_two_intervals_is_counted_once() {
        let two = [interval("chr1", 150, 160), interval("chr1", 350, 360)];
        assert_eq!(
            count(SPANNING, Some(&two), true, "spanning.vcf").unwrap(),
            1
        );
    }

    #[test]
    fn intervals_without_an_index_are_refused_before_any_record() {
        let one = [interval("chr1", 100, 200)];
        let refused = count(
            SPANNING,
            Some(&one),
            false,
            "countvariants-dump/unindexed.vcf",
        )
        .expect_err("no index");
        assert_eq!(
            refused.message(),
            "Input countvariants-dump/unindexed.vcf must support random access to enable \
             traversal by intervals. If it's a file, please index it using the bundled tool \
             IndexFeatureFile"
        );
    }

    #[test]
    fn no_intervals_needs_no_index() {
        assert_eq!(count(SPANNING, None, false, "spanning.vcf").unwrap(), 2);
        assert_eq!(
            count(SPANNING, Some(&[]), false, "spanning.vcf").unwrap(),
            2
        );
    }

    #[test]
    fn the_output_file_refusal_is_the_path_alone() {
        let refused = CountVariantsError::CouldNotCreateOutputFile {
            path: "countvariants-dump/.".to_string(),
        };
        assert_eq!(refused.message(), "countvariants-dump/.");
        assert_eq!(
            refused.class(),
            "org.broadinstitute.hellbender.exceptions.UserException$CouldNotCreateOutputFile"
        );
    }
}
