//! `RenameSampleInVcf.doWork`, over [`gatk_tools::rename_sample_in_vcf`].
//!
//! The library does the renaming and its two `IllegalArgumentException`s. What is here is the
//! file handling around it and one exception the library cannot tell apart: `OLD_SAMPLE_NAME`
//! given for a file with no sample at all, where the reference never reaches its own message
//! because building it calls `get(0)` on an empty list:
//!
//! ```java
//! if (OLD_SAMPLE_NAME != null && !OLD_SAMPLE_NAME.equals(header.getGenotypeSamples().get(0))) {
//! ```
//!
//! # Reading this file without reading Rust
//!
//! - `?` returns the error to the caller at once, Java's `throws`;
//! - `old.as_deref()` lends the `Option<String>` as an `Option<&str>`, a value or nothing,
//!   without copying the string;
//! - `match error { Pattern if condition => ..., other => ... }` picks the first arm whose
//!   shape and guard fit; `other` catches everything else, and the compiler checks no case is
//!   left out.

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::rename_sample_in_vcf::{rename, RenameError};

use super::{assert_readable, assert_writable, thrown, vcf_text, write_vcf};
use crate::runners::{argument, flag};

/// `RenameSampleInVcf.doWork`: INPUT readable, OUTPUT writable, then the rename and the write.
pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let text = vcf_text(&input)?;
    let old = argument(parser, "OLD_SAMPLE_NAME");
    let written = rename(
        &text,
        &argument(parser, "NEW_SAMPLE_NAME").unwrap_or_default(),
        old.as_deref(),
    )
    .map_err(|error| match error {
        // `getGenotypeSamples().get(0)` of a file with no sample at all. The library reports an
        // empty `contained` both for this and for a sample whose name is empty, so the header
        // line is looked at again to tell the two apart.
        RenameError::UnexpectedSample { contained }
            if contained.is_empty() && no_samples(&text) =>
        {
            thrown(
                "java.lang.IndexOutOfBoundsException",
                "Index 0 out of bounds for length 0".to_string(),
            )
        }
        other => thrown(other.java_class(), other.message()),
    })?;
    write_vcf(parser, &output, &written, flag(parser, "CREATE_INDEX"))?;
    Ok(0)
}

/// Whether the `#CHROM` line names no sample: nine columns or fewer (`#CHROM` through `FORMAT`),
/// or no such line at all.
///
/// `is_none_or` is true for a missing line and otherwise asks the closure about the line.
fn no_samples(text: &str) -> bool {
    text.lines()
        .find(|line| line.starts_with("#CHROM"))
        .is_none_or(|line| line.split('\t').count() <= 9)
}
