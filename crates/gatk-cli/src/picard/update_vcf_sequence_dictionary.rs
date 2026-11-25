//! `UpdateVcfSequenceDictionary.doWork`, over [`gatk_tools::picard_update_vcf_dictionary`].
//!
//! `SAMSequenceDictionaryExtractor` reads the dictionary by the file's kind: a FASTA's is the
//! `.dict` beside it, a VCF's is its contig lines (none at all is `null`), and anything else is a
//! SAM header's `@SQ` lines. A `null` dictionary is refused by the builder when it is to write a
//! Tribble index, and by `setSequenceDictionary` otherwise.

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::picard_update_vcf_dictionary::{read_dictionary, update, updated_header, Sequence};

use super::{assert_readable, assert_writable, thrown, vcf_lines, vcf_text, write_vcf, VcfStream};
use crate::runners::{argument, flag};
use crate::vcf_stream::{Body, Records};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let input = argument(parser, "INPUT").unwrap_or_default();
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    let source = argument(parser, "SEQUENCE_DICTIONARY").unwrap_or_default();
    assert_readable(&input)?;
    assert_readable(&source)?;
    assert_writable(&output)?;

    let dictionary = extract_dictionary(&source)?;
    if let Some(sequences) = &dictionary {
        if streamed(parser, &input, &output, sequences).is_some() {
            return Ok(0);
        }
    }
    let text = vcf_text(&input)?;
    let compressed = output.ends_with(".gz") || output.ends_with(".bgz");
    if flag(parser, "CREATE_INDEX") && dictionary.is_none() && !compressed {
        return Err(thrown(
            "java.lang.IllegalArgumentException",
            "A reference dictionary is required for creating Tribble indices on the fly"
                .to_string(),
        ));
    }
    let Some(dictionary) = dictionary else {
        return Err(thrown(
            "java.lang.NullPointerException",
            "Cannot invoke \"htsjdk.samtools.SAMSequenceDictionary.getSequences()\" because \
             \"dictionary\" is null"
                .to_string(),
        ));
    };
    let updated =
        update(&text, &dictionary).map_err(|error| thrown(error.java_class(), error.message()))?;
    write_vcf(parser, &output, &updated, flag(parser, "CREATE_INDEX"))?;
    Ok(0)
}

/// `update` and [`write_vcf`] a record at a time. `None` on any refusal, having written nothing,
/// so the run above refuses in its own words and order.
fn streamed(parser: &Parser, input: &str, output: &str, dictionary: &[Sequence]) -> Option<()> {
    let mut lines = vcf_lines(input)?;
    let (text, header_lines) = lines.read_header().ok()?;
    let file = htsjdk_vcf::reader::read_vcf(&text).ok()?;
    let header = updated_header(&file.header, dictionary);
    let mut out = VcfStream::create(parser, output, &header, flag(parser, "CREATE_INDEX"), None)?;
    for record in Records::new(Body::Stream(lines), &file, header_lines) {
        out.write(&record.ok()?)?;
    }
    out.finish(parser)
}

/// `SAMSequenceDictionaryExtractor.extractDictionary`, `None` for a VCF with no contig lines.
fn extract_dictionary(path: &str) -> Result<Option<Vec<Sequence>>, Thrown> {
    let read = |path: &str| {
        std::fs::read_to_string(path)
            .map_err(|error| thrown("java.io.UncheckedIOException", format!("{path}: {error}")))
    };
    let lower = path.to_lowercase();
    for extension in [".fasta", ".fa", ".fna"] {
        if let Some(stem) = lower.strip_suffix(extension) {
            let beside = format!("{}.dict", &path[..stem.len()]);
            return Ok(Some(read_dictionary(&read(&beside)?)));
        }
    }
    if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
        let file = htsjdk_vcf::reader::read_vcf(&vcf_text(path)?)
            .map_err(|failure| thrown(failure.error.class(), failure.error.message()))?;
        return Ok(
            gatk_tools::vcf_to_interval_list::dictionary(&file.header).map(|sequences| {
                sequences
                    .into_iter()
                    .map(|sequence| Sequence {
                        name: sequence.name,
                        length: sequence.length,
                        assembly: sequence.assembly,
                    })
                    .collect()
            }),
        );
    }
    Ok(Some(read_dictionary(&read(path)?)))
}
