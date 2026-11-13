//! `VcfFormatConverter.doWork`, over [`gatk_tools::vcf_format_converter`].
//!
//! The output's format is the writer builder's guess from its name: `.vcf` is text, `.vcf.gz` and
//! `.vcf.bgz` are BGZF, and anything else but `.bcf` is the builder's refusal. `REQUIRE_INDEX`
//! looks for the index the reader would open: a `.tbi` beside a block-compressed name, else `.idx`.

use gatk_barclay::Parser;
use gatk_tools::main_entry::{Thrown, PORT_LIMITATION};
use gatk_tools::vcf_format_converter::{convert, Arguments, Input};

use super::{absolute, assert_readable, assert_writable, thrown, vcf_text, write_vcf};
use crate::runners::{argument, flag};

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
    let text = vcf_text(&input)?;
    let path = absolute(&input);
    let converted = convert(
        &Input {
            path: &path,
            text: &text,
            indexed,
        },
        &Arguments {
            require_index: flag(parser, "REQUIRE_INDEX"),
            create_index: flag(parser, "CREATE_INDEX"),
        },
    )
    .map_err(|error| thrown(error.java_class(), error.message()))?;
    output_type(&output)?;
    write_vcf(parser, &output, &converted, flag(parser, "CREATE_INDEX"))?;
    Ok(0)
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
