//! `CheckTerminatorBlock.doWork`, over [`gatk_tools::check_terminator_block`].

use gatk_barclay::Parser;
use gatk_tools::check_terminator_block::{check, termination_name};
use gatk_tools::main_entry::Thrown;

use super::{absolute, assert_readable, thrown};
use crate::runners::argument;

pub(super) fn run(parser: &Parser, stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    assert_readable(&input)?;
    // `checkTermination` reads the last block's worth of the file and nothing more.
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
    stderr.push_str(termination_name(termination));
    stderr.push('\n');
    Ok(code)
}
