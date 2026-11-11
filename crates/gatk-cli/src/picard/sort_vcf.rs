//! `SortVcf.doWork`, over [`gatk_tools::sort_vcf`].
//!
//! The library sorts texts; what is here is what the tool does with paths first: every input and
//! the dictionary asserted readable, the dictionary read through `SamReaderFactory`, and handed to
//! each input that has none of its own. An input that declares one keeps it, and the given one is
//! only checked against it.

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::picard_gather_vcfs::{assert_same_dictionary, dictionary, Sequence};
use gatk_tools::sort_vcf::{sort, SortError};

use super::{absolute, assert_readable, thrown, vcf_text, with_contigs, write_vcf};
use crate::runners::{argument, arguments, dictionary_from_any, flag};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let inputs = arguments(parser, "INPUT");
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    for input in &inputs {
        assert_readable(input)?;
    }
    let dictionary = match argument(parser, "SEQUENCE_DICTIONARY") {
        None => None,
        Some(path) => {
            assert_readable(&path)?;
            Some(dictionary_from_any(&path)?)
        }
    };

    let mut texts = Vec::new();
    for input in &inputs {
        let text = vcf_text(input)?;
        let has_contigs = text
            .lines()
            .take_while(|line| line.starts_with("##"))
            .any(|line| line.starts_with("##contig="));
        texts.push(match (&dictionary, has_contigs) {
            // `header.setSequenceDictionary(samSequenceDictionary)`.
            (Some(sequences), false) => with_contigs(&text, sequences),
            // `samSequenceDictionary.assertSameDictionary(dict)`, rethrown as the cause of an
            // `IllegalArgumentException`.
            (Some(sequences), true) => {
                if let Err(message) = assert_same_dictionary(&given(sequences), &own(&text)) {
                    return Err(thrown(
                        "java.lang.IllegalArgumentException",
                        format!("java.lang.AssertionError: {message}"),
                    ));
                }
                text
            }
            _ => text,
        });
    }

    let sorted = sort(&texts).map_err(|error| {
        // The library masks the path its message names; the tool prints the input's absolute one.
        let culprit = match error {
            SortError::MissingDictionary => {
                texts.iter().position(|text| !text.contains("\n##contig="))
            }
            SortError::DifferentSamples => {
                let first = samples(&texts[0]);
                texts.iter().position(|text| samples(text) != first)
            }
            _ => None,
        };
        let message = match culprit {
            Some(index) => error
                .message()
                .replace("<masked>", &absolute(&inputs[index])),
            None => error.message(),
        };
        thrown(error.java_class(), message)
    })?;
    write_vcf(parser, &output, &sorted, flag(parser, "CREATE_INDEX"))?;
    Ok(0)
}

/// `getSampleNamesInOrder()`, which is sorted.
fn samples(text: &str) -> Vec<String> {
    let mut names: Vec<String> = text
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .map(|line| line.split('\t').skip(9).map(str::to_string).collect())
        .unwrap_or_default();
    names.sort();
    names
}

/// The input's own dictionary, from its contig lines.
fn own(text: &str) -> Vec<Sequence> {
    htsjdk_vcf::reader::read_vcf(text)
        .ok()
        .and_then(|file| dictionary(&file.header))
        .unwrap_or_default()
}

/// The given dictionary, as `SAMSequenceRecord`s print.
fn given(sequences: &[htsjdk_bam::header::SequenceRecord]) -> Vec<Sequence> {
    sequences
        .iter()
        .enumerate()
        .map(|(index, record)| Sequence {
            name: record.name.clone(),
            length: i64::from(record.length),
            index,
            assembly: record.attributes.get("AS").map(str::to_string),
        })
        .collect()
}
