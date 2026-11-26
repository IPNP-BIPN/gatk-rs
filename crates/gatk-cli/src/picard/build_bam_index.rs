//! `BuildBamIndex.doWork`, over [`gatk_tools::build_bam_index`].
//!
//! With no `OUTPUT` the index is named after the input's FILE NAME and resolved against the
//! working directory, not beside the input:
//!
//! ```java
//! final String baseFileName = inputPath.getFileName().toString();
//! ```
//!
//! ```java
//! OUTPUT = new File(baseFileName.substring(0, index) + FileExtensions.BAI_INDEX);
//! ```
//!
//! # Reading this file without reading Rust
//!
//! - `Result<i32, Thrown>` is "the exit code, or the exception the reference would throw";
//!   `?` after a call hands that exception straight back to the caller, Java's `throws`;
//! - `.map_err(|error| ...)` rewrites the error a call failed with and leaves a success alone;
//!   `|error| ...` is a closure, Java's lambda `error -> ...`;
//! - `Option<String>` is a value or nothing; `.unwrap_or_default()` turns nothing into `""` and
//!   `.unwrap_or_else(|| ...)` computes the fallback only when it is needed.

use gatk_barclay::Parser;
use gatk_tools::build_bam_index::{build, default_output};
use gatk_tools::main_entry::{Thrown, PORT_FAILURE};

use super::{assert_writable, thrown};
use crate::runners::argument;

/// `BuildBamIndex.doWork`: name the output, check it is writable, then build and write the BAI.
///
/// INPUT is not asserted readable first; the reference only asserts the OUTPUT and lets
/// `SamReaderFactory.open` refuse a missing input. `build` answers the two `SAMException`s the
/// reference throws after opening (not a BAM, not coordinate sorted) under their own class and
/// message.
pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    // No OUTPUT: the input's last path component, never its directory, so the index lands in
    // the working directory. `default_output` swaps a `.bam` suffix for `.bai` or appends `.bai`.
    let output = argument(parser, "OUTPUT").unwrap_or_else(|| {
        let name = std::path::Path::new(&input)
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default();
        default_output(&name)
    });
    assert_writable(&output)?;
    // `SamReaderFactory.open`, which is where a missing input is refused, as the stream's own.
    let file = std::fs::read(&input).map_err(|_| {
        thrown(
            "htsjdk.samtools.util.RuntimeIOException",
            format!(
                "java.io.FileNotFoundException: {} (No such file or directory)",
                super::absolute(&input)
            ),
        )
    })?;
    let bai = build(&file).map_err(|error| thrown(error.java_class(), error.message()))?;
    std::fs::write(&output, bai)
        .map_err(|error| Thrown::non_user(PORT_FAILURE, format!("{output}: {error}")))?;
    Ok(0)
}
