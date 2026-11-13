//! `VcfToIntervalList.doWork`, over [`gatk_tools::vcf_to_interval_list`].

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::vcf_to_interval_list::{convert, IdMethod};

use super::{assert_readable, assert_writable, thrown, vcf_text};
use crate::runners::{argument, flag, scalar, write_file};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let method = match scalar(parser, "VARIANT_ID_METHOD").as_deref() {
        Some("USE_FIRST") => IdMethod::UseFirst,
        _ => IdMethod::ConcatAll,
    };
    let list = convert(&vcf_text(&input)?, flag(parser, "INCLUDE_FILTERED"), method)
        .map_err(|error| thrown(error.java_class(), error.message()))?;
    write_file(&output, list.as_bytes())?;
    Ok(0)
}
