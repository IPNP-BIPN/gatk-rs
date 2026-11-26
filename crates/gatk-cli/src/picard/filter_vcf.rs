//! `FilterVcf.doWork`.

use gatk_barclay::Parser;
use gatk_tools::filter_vcf::{filter_to, ScriptFilter, Thresholds};
use gatk_tools::javascript_variant_filter::{JavascriptVariantFilter, SCRIPT_ERROR};
use gatk_tools::main_entry::Thrown;
use htsjdk_vcf::reader::{read_vcf, ReadError, ReadFailure};

use super::vcf_format_converter::{is_bcf, output_type};
use super::{assert_readable, assert_writable, thrown, vcf_text, write_vcf};
use crate::runners::{argument, flag, scalar};

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
    let header_fails = || {
        matches!(
            read_vcf(&text),
            Err(ReadFailure {
                error: ReadError::Header(_) | ReadError::HeaderLine(_),
                ..
            })
        )
    };
    let mut script = match argument(parser, "JAVASCRIPT_FILE").filter(|file| file != "null") {
        Some(file) if !header_fails() => Some(script_filter(&file)?),
        _ => None,
    };
    let filtered = filter_to(&text, &thresholds, vcf_or_bcf, script.as_mut())
        .map_err(|error| thrown(error.java_class(), error.message()))?;
    output_type(&output)?;
    // Picard 3.5.0 sets or unsets `INDEX_ON_THE_FLY` from `CREATE_INDEX` (b976447cf); 3.4.0 kept
    // the builder's default and indexed whatever it said.
    write_vcf(parser, &output, &filtered, flag(parser, "CREATE_INDEX"))?;
    Ok(0)
}

/// `new VariantContextJavascriptFilter(JAVASCRIPT_FILE, header)`: the file read, the script
/// compiled, and the filter named after `IOUtil.basename`, or `JSFILTER` when that is empty.
fn script_filter(file: &str) -> Result<ScriptFilter, Thrown> {
    // `new FileReader(scriptFile)`, whose IOException Picard wraps.
    let bytes = std::fs::read(file).map_err(|_| {
        thrown(
            "picard.PicardException",
            "javascript-related error".to_string(),
        )
    })?;
    let script = JavascriptVariantFilter::new(&String::from_utf8_lossy(&bytes)).map_err(|_| {
        thrown(
            "htsjdk.samtools.util.RuntimeScriptException",
            SCRIPT_ERROR.to_string(),
        )
    })?;
    let file_name = file.rsplit('/').next().unwrap_or_default();
    let name = match file_name.rfind('.') {
        Some(dot) if dot > 0 => &file_name[..dot],
        _ => file_name,
    };
    Ok(ScriptFilter {
        name: if name.is_empty() { "JSFILTER" } else { name }.to_string(),
        // shortcut: `File.toString()` also collapses `//` and drops a trailing `/`; a path given
        // that way is written as given.
        path: file.to_string(),
        script,
    })
}
