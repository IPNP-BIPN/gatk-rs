//! `VcfFormatConverter.doWork`, over [`gatk_tools::vcf_format_converter`].
//!
//! The output's format is the writer builder's guess from its name: `.vcf` is text, `.vcf.gz` and
//! `.vcf.bgz` are BGZF, and anything else but `.bcf` is the builder's refusal. `REQUIRE_INDEX`
//! looks for the index the reader would open: a `.tbi` beside a block-compressed name, else `.idx`.

use gatk_barclay::Parser;
use gatk_tools::main_entry::{Thrown, PORT_LIMITATION};
use gatk_tools::vcf_format_converter::{convert, Arguments, Input};

use super::{
    absolute, assert_readable, assert_writable, thrown, vcf_lines, vcf_text, write_vcf, VcfStream,
};
use crate::runners::{argument, flag};
use crate::vcf_stream::{Body, Records};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    if is_bcf(&input) || is_bcf(&output) {
        return Err(bcf_limitation());
    }
    let tabix = [".gz", ".gzip", ".bgz", ".bgzf"]
        .iter()
        .any(|extension| input.to_lowercase().ends_with(extension))
        && std::path::Path::new(&format!("{input}.tbi")).exists();
    let indexed = tabix || std::path::Path::new(&format!("{input}.idx")).exists();
    let (require_index, create_index) =
        (flag(parser, "REQUIRE_INDEX"), flag(parser, "CREATE_INDEX"));
    if (indexed || !require_index) && streamed(parser, &input, &output, create_index).is_some() {
        return Ok(0);
    }
    let text = vcf_text(&input)?;
    let path = absolute(&input);
    let converted = convert(
        &Input {
            path: &path,
            text: &text,
            indexed,
        },
        &Arguments {
            require_index,
            create_index,
        },
    )
    .map_err(|error| thrown(error.java_class(), error.message()))?;
    output_type(&output)?;
    write_vcf(parser, &output, &converted, create_index)?;
    Ok(0)
}

/// `convert` and [`write_vcf`] a record at a time. `None` on any refusal, having written nothing,
/// so the run above refuses in its own words and order.
fn streamed(parser: &Parser, input: &str, output: &str, create_index: bool) -> Option<()> {
    let mut lines = vcf_lines(input)?;
    let (text, header_lines) = lines.read_header().ok()?;
    let file = htsjdk_vcf::reader::read_vcf(&text).ok()?;
    let has_dictionary = file
        .header
        .lines
        .iter()
        .any(|line| matches!(line, htsjdk_vcf::header::HeaderLine::Contig { .. }));
    if create_index && !has_dictionary {
        return None;
    }
    output_type(output).ok()?;
    let mut out = VcfStream::create(parser, output, &file.header, create_index, None)?;
    for record in Records::new(Body::Stream(lines), &file, header_lines) {
        out.write(&record.ok()?)?;
    }
    out.finish(parser)
}

/// `VariantContextWriterBuilder.build()`'s refusal of a name it cannot infer a format from.
pub(super) fn output_type(output: &str) -> Result<(), Thrown> {
    if is_bcf(output) {
        return Err(bcf_limitation());
    }
    if output.ends_with(".vcf") || output.ends_with(".vcf.gz") || output.ends_with(".vcf.bgz") {
        return Ok(());
    }
    Err(thrown(
        "java.lang.IllegalArgumentException",
        "Output format type is not set, or could not be inferred from the output path. If a \
         path was used, does it have a valid VCF extension (.vcf, .vcf.gz, .vcf.bgz, .bcf)?"
            .to_string(),
    ))
}

pub(super) fn bcf_limitation() -> Thrown {
    Thrown::non_user(
        PORT_LIMITATION,
        "BCF is a format this port does not read or write yet. This message is the port's own and not GATK's.",
    )
}

pub(super) fn is_bcf(path: &str) -> bool {
    path.to_lowercase().ends_with(".bcf")
}
