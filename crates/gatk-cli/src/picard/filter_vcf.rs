//! `FilterVcf.doWork`.

use gatk_barclay::Parser;
use gatk_tools::filter_vcf::{filter_stream, filter_to, FilterError, Thresholds};
use gatk_tools::main_entry::{Thrown, PORT_FAILURE, PORT_LIMITATION};
use htsjdk_vcf::reader::{read_vcf, ReadError, ReadFailure};

use super::vcf_format_converter::{is_bcf, output_type};
use super::{assert_readable, assert_writable, thrown, vcf_text, write_vcf};
use crate::runners::{argument, flag, scalar, StreamedVcf};

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
    let javascript = argument(parser, "JAVASCRIPT_FILE").is_some_and(|file| file != "null");
    if !javascript && output.ends_with(".vcf") {
        if let Some(lines) = PlainLines::open(&input)? {
            return streamed(
                lines,
                &input,
                &output,
                &thresholds,
                flag(parser, "CREATE_INDEX"),
            );
        }
    }
    let text = vcf_text(&input)?;
    // The script filter is built once `new VCFFileReader` has read the header, before the first
    // record is read and before the dictionary check: only a header failure comes first.
    if javascript
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
    // Picard 3.5.0 sets or unsets `INDEX_ON_THE_FLY` from `CREATE_INDEX` (b976447cf); 3.4.0 kept
    // the builder's default and indexed whatever it said.
    write_vcf(parser, &output, &filtered, flag(parser, "CREATE_INDEX"))?;
    Ok(0)
}

/// A plain-text `.vcf` input to a `.vcf` output, a record at a time: what the whole-file path
/// writes, without the file, its records or the output ever held whole.
fn streamed(
    mut lines: PlainLines,
    input: &str,
    output: &str,
    thresholds: &Thresholds,
    create_index: bool,
) -> Result<i32, Thrown> {
    let mut sink = StreamedVcf::create(output, create_index)?;
    let mut write_error = None;
    let result = filter_stream(&mut lines, thresholds, true, |text| {
        sink.write(text).map_err(|error| {
            write_error = Some(error);
            FilterError::Vcf(String::new(), String::new())
        })
    });
    let failure = if let Some(error) = lines.error.take() {
        Thrown::non_user(PORT_FAILURE, format!("{input}: {error}"))
    } else if let Some(error) = write_error {
        Thrown::non_user(PORT_FAILURE, format!("could not write {output}: {error}"))
    } else if let Err(error) = result {
        thrown(error.java_class(), error.message())
    } else {
        sink.finish()?;
        return Ok(0);
    };
    sink.abandon();
    Err(failure)
}

/// The lines of a file that is not gzipped, as `vcf_text(path).lines()` gives them: decoded
/// lossily, without their `\n` or `\r\n`. A read failure ends the lines and is kept in `error`.
struct PlainLines {
    reader: std::io::BufReader<std::fs::File>,
    line: Vec<u8>,
    error: Option<std::io::Error>,
}

impl PlainLines {
    /// `None` for a gzipped file, which [`vcf_text`] inflates whole.
    fn open(path: &str) -> Result<Option<Self>, Thrown> {
        use std::io::BufRead;
        let failed =
            |error: std::io::Error| Thrown::non_user(PORT_FAILURE, format!("{path}: {error}"));
        let file = std::fs::File::open(path).map_err(failed)?;
        let mut reader = std::io::BufReader::with_capacity(1 << 20, file);
        if reader
            .fill_buf()
            .map_err(failed)?
            .starts_with(&[0x1f, 0x8b])
        {
            return Ok(None);
        }
        Ok(Some(PlainLines {
            reader,
            line: Vec::new(),
            error: None,
        }))
    }
}

impl Iterator for PlainLines {
    type Item = String;

    fn next(&mut self) -> Option<String> {
        use std::io::BufRead;
        self.line.clear();
        match self.reader.read_until(b'\n', &mut self.line) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(error) => {
                self.error = Some(error);
                return None;
            }
        }
        // `str::lines`: a `\r` goes only with the `\n` after it.
        if self.line.ends_with(b"\n") {
            self.line.pop();
            if self.line.ends_with(b"\r") {
                self.line.pop();
            }
        }
        Some(String::from_utf8_lossy(&self.line).into_owned())
    }
}
