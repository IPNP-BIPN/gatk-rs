//! `MakeSitesOnlyVcf.doWork`, over [`gatk_tools::make_sites_only_vcf`].

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;

use super::{assert_readable, assert_writable, thrown, vcf_text, write_vcf};
use crate::runners::{argument, arguments, flag};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let text = vcf_text(&input)?;
    let written = gatk_tools::make_sites_only_vcf::make_sites_only(
        &text,
        &arguments(parser, "SAMPLE"),
        flag(parser, "CREATE_INDEX"),
    )
    .map_err(|error| thrown(error.java_class(), error.message()))?;
    write_vcf(parser, &output, &written, flag(parser, "CREATE_INDEX"))?;
    Ok(0)
}
