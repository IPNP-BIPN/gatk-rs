//! `BuildBamIndex.doWork`, over [`gatk_tools::build_bam_index`].
//!
//! With no `OUTPUT` the index is named after the input's FILE NAME and resolved against the
//! working directory, not beside the input.

use gatk_barclay::Parser;
use gatk_tools::build_bam_index::{build, default_output};
use gatk_tools::main_entry::{Thrown, PORT_FAILURE};

use super::{assert_writable, thrown};
use crate::runners::argument;

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
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
