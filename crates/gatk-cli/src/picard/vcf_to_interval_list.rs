//! `VcfToIntervalList.doWork`, over [`gatk_tools::vcf_to_interval_list`].

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::vcf_to_interval_list::{
    convert, dictionary, list_header, to_interval, write_interval, IdMethod, Merger,
};
use std::io::{BufRead, Write};

use super::{assert_readable, assert_writable, thrown, vcf_text};
use crate::runners::{argument, flag, scalar, write_file};
use crate::vcf_stream::{Body, LineReader, Records};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let method = match scalar(parser, "VARIANT_ID_METHOD").as_deref() {
        Some("USE_FIRST") => IdMethod::UseFirst,
        _ => IdMethod::ConcatAll,
    };
    let include_filtered = flag(parser, "INCLUDE_FILTERED");
    if streamed(&input, &output, include_filtered, method).is_some() {
        return Ok(0);
    }
    let list = convert(&vcf_text(&input)?, include_filtered, method)
        .map_err(|error| thrown(error.java_class(), error.message()))?;
    write_file(&output, list.as_bytes())?;
    Ok(0)
}

/// [`convert`] a record at a time, written beside `output` and renamed onto it at the end.
///
/// `None` when anything is refused or fails: the scratch file is removed and nothing is written,
/// so the whole-text path refuses in its own words. The text is the one [`vcf_text`] reads: a
/// file starting with the gzip magic is inflated by the reader `decompress_all` runs serially,
/// which refuses the same files, and both are decoded lossily, line by line here.
fn streamed(input: &str, output: &str, include_filtered: bool, method: IdMethod) -> Option<()> {
    let mut file = std::io::BufReader::new(std::fs::File::open(input).ok()?);
    let reader: Box<dyn BufRead> = if file.fill_buf().ok()?.starts_with(&[0x1f, 0x8b]) {
        Box::new(std::io::BufReader::new(htsjdk_bgzf::read::BgzfReader::new(
            file,
        )))
    } else {
        Box::new(file)
    };
    let mut lines = LineReader::new(reader, false);
    let (header, header_lines) = lines.read_header().ok()?;
    let vcf = htsjdk_vcf::reader::read_vcf(&header).ok()?;
    let sequences = dictionary(&vcf.header)?;

    let partial = format!("{output}.partial");
    let written = (|| {
        let mut out = std::io::BufWriter::new(std::fs::File::create(&partial).ok()?);
        out.write_all(list_header(&sequences).as_bytes()).ok()?;
        let (mut count, mut merger) = (0, Merger::new(method.concatenates()));
        for record in Records::new(Body::Stream(lines), &vcf, header_lines) {
            let record = record.ok()?;
            if let Some(interval) = to_interval(&record, include_filtered, &mut count) {
                if let Some(merged) = merger.push(interval) {
                    write_interval(&mut out, &merged).ok()?;
                }
            }
        }
        if let Some(merged) = merger.finish() {
            write_interval(&mut out, &merged).ok()?;
        }
        out.into_inner().ok()?;
        std::fs::rename(&partial, output).ok()
    })();
    if written.is_none() {
        let _ = std::fs::remove_file(&partial);
    }
    written
}
