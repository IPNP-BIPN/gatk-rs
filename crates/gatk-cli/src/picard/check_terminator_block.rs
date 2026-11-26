//! `CheckTerminatorBlock.doWork`, over [`gatk_tools::check_terminator_block`].
//!
//! The termination's name goes to stderr and the exit code is 100 for `DEFECTIVE`, 0 otherwise:
//!
//! ```java
//! final FileTermination term = BlockCompressedInputStream.checkTermination(INPUT);
//! System.err.println(term.name());
//! ```
//!
//! # Reading this file without reading Rust
//!
//! `?` after a call returns its error to the caller at once, Java's `throws`; `.map_err(|_| ...)`
//! replaces whatever error the read failed with by the one the reference throws (`|_| ...` is a
//! closure that ignores its argument). `let (a, b) = f()` unpacks the pair `check` returns.

use gatk_barclay::Parser;
use gatk_tools::check_terminator_block::{check, termination_name};
use gatk_tools::main_entry::Thrown;

use super::{absolute, assert_readable, thrown};
use crate::runners::argument;

/// `CheckTerminatorBlock.doWork`: the readable check, the termination, its name, its exit code.
pub(super) fn run(parser: &Parser, stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    assert_readable(&input)?;
    // `checkTermination` reads the last block's worth of the file and nothing more. The port reads
    // the whole file and `check` looks at its tail; an `IOException` there becomes the reference's
    // own `PicardException`, quoted from its `catch`:
    // `throw new PicardException("Exception reading terminator block of file: "`
    // `+ INPUT.getAbsolutePath());`
    let data = std::fs::read(&input).map_err(|_| {
        thrown(
            "picard.PicardException",
            format!(
                "Exception reading terminator block of file: {}",
                absolute(&input)
            ),
        )
    })?;
    let (termination, code) = check(&data);
    // `System.err.println(term.name())`: the enum constant's name, then a newline.
    stderr.push_str(termination_name(termination));
    stderr.push('\n');
    Ok(code)
}
