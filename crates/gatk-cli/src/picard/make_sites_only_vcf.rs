//! `MakeSitesOnlyVcf.doWork`, over [`gatk_tools::make_sites_only_vcf`].

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;

use super::{assert_readable, assert_writable, thrown, vcf_lines, vcf_text, write_vcf, VcfStream};
use crate::runners::{argument, arguments, flag};
use crate::vcf_stream::{Body, Records};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    assert_readable(&input)?;
    assert_writable(&output)?;
    let samples = arguments(parser, "SAMPLE");
    let create_index = flag(parser, "CREATE_INDEX");
    if streamed(parser, &input, &output, &samples, create_index).is_some() {
        return Ok(0);
    }
    let text = vcf_text(&input)?;
    let written = gatk_tools::make_sites_only_vcf::make_sites_only(&text, &samples, create_index)
        .map_err(|error| thrown(error.java_class(), error.message()))?;
    write_vcf(parser, &output, &written, create_index)?;
    Ok(0)
}

/// `make_sites_only` and [`write_vcf`] a record at a time. `None` on any refusal, having written
/// nothing, so the whole-text path refuses in its own words.
fn streamed(
    parser: &Parser,
    input: &str,
    output: &str,
    samples: &[String],
    create_index: bool,
) -> Option<()> {
    let mut lines = vcf_lines(input)?;
    let (header_text, header_lines) = lines.read_header().ok()?;
    let file = htsjdk_vcf::reader::read_vcf(&header_text).ok()?;
    let has_dictionary = file
        .header
        .lines
        .iter()
        .any(|line| matches!(line, htsjdk_vcf::header::HeaderLine::Contig { .. }));
    if create_index && !has_dictionary {
        return None;
    }
    let wanted = gatk_tools::make_sites_only_vcf::wanted_samples(samples);
    let mut header = file.header.clone();
    header.samples = wanted.clone();
    let mut out = VcfStream::create(parser, output, &header, create_index, None)?;
    for record in Records::new(Body::Stream(lines), &file, header_lines) {
        let mut record = record.ok()?;
        record
            .genotypes
            .retain(|genotype| wanted.contains(&genotype.sample_name));
        out.write(&record)?;
    }
    out.finish(parser)
}
