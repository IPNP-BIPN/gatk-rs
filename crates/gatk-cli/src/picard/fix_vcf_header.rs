//! `FixVcfHeader.doWork`, after `customCommandLineValidation`.

use gatk_barclay::Parser;
use gatk_tools::fix_vcf_header::{fix, fix_with_header};
use gatk_tools::main_entry::Thrown;

use super::vcf_format_converter::output_type;
use super::{assert_readable, thrown, vcf_text, write_vcf};
use crate::runners::{argument, flag, scalar};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    let header = argument(parser, "HEADER").filter(|header| header != "null");
    let check_first_n_records = scalar(parser, "CHECK_FIRST_N_RECORDS")
        .and_then(|value| value.parse::<i32>().ok())
        .unwrap_or(-1);
    assert_readable(&input)?;
    if let Some(header) = &header {
        assert_readable(header)?;
    }
    let text = vcf_text(&input)?;
    let fixed = match &header {
        Some(header) => fix_with_header(
            &text,
            &vcf_text(header)?,
            flag(parser, "ENFORCE_SAME_SAMPLES"),
        ),
        None => fix(&text, check_first_n_records),
    }
    .map_err(|error| thrown(error.java_class(), error.message()))?;

    output_type(&output)?;
    write_indexed_vcf(parser, &output, &fixed)?;
    Ok(0)
}

/// `write_vcf` for a writer built with htsjdk's default options, which hold `INDEX_ON_THE_FLY`:
/// Picard's `CREATE_INDEX` reaches only `SAMFileWriterFactory`, so the VCF is indexed whatever it
/// says.
pub(super) fn write_indexed_vcf(parser: &Parser, output: &str, text: &str) -> Result<(), Thrown> {
    write_vcf(parser, output, text, true)
}
