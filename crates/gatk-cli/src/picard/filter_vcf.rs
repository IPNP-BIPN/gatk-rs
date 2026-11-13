//! `FilterVcf.doWork`.

use gatk_barclay::Parser;
use gatk_tools::filter_vcf::{filter_to, Thresholds};
use gatk_tools::main_entry::{Thrown, PORT_LIMITATION};
use htsjdk_vcf::reader::{read_vcf, ReadError, ReadFailure};

use super::fix_vcf_header::write_indexed_vcf;
use super::vcf_format_converter::{is_bcf, output_type};
use super::{assert_readable, assert_writable, thrown, vcf_text};
use crate::runners::{argument, scalar};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let defaults = Thresholds::default();
    let number = |name: &str, default: f64| {
        scalar(parser, name)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(default)
    };
    let integer = |name: &str, default: i32| {
        scalar(parser, name)
            .and_then(|value| value.parse::<i32>().ok())
            .unwrap_or(default)
    };
    let thresholds = Thresholds {
        min_ab: number("MIN_AB", defaults.min_ab),
        max_fs: number("MAX_FS", defaults.max_fs),
        min_qd: number("MIN_QD", defaults.min_qd),
        min_gq: integer("MIN_GQ", defaults.min_gq),
        min_dp: integer("MIN_DP", defaults.min_dp),
    };
    let file_name = output.rsplit('/').next().unwrap_or_default();
    let vcf_or_bcf = file_name.ends_with(".vcf") || file_name.ends_with(".bcf");
    if is_bcf(&output) {
        return Err(super::vcf_format_converter::bcf_limitation());
    }
    let text = vcf_text(&input)?;
    // The script filter is built once `new VCFFileReader` has read the header, before the first
    // record is read and before the dictionary check: only a header failure comes first.
    if argument(parser, "JAVASCRIPT_FILE").is_some_and(|file| file != "null")
        && !matches!(
            read_vcf(&text),
            Err(ReadFailure {
                error: ReadError::Header(_) | ReadError::HeaderLine(_),
                ..
            })
        )
    {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "FilterVcf's JAVASCRIPT_FILE runs a Nashorn script, which this port does not. This \
             message is the port's own and not GATK's.",
        ));
    }
    let filtered = filter_to(&text, &thresholds, vcf_or_bcf)
        .map_err(|error| thrown(error.java_class(), error.message()))?;
    output_type(&output)?;
    // `new VariantContextWriterBuilder()` keeps the default `INDEX_ON_THE_FLY`.
    write_indexed_vcf(parser, &output, &filtered)?;
    Ok(0)
}
