//! `MergeVcfs.doWork`, over [`gatk_tools::merge_vcfs`].
//!
//! The library merges texts and checks their headers in input order; what is here is the order
//! the tool meets its files in. The dictionary is opened first, the inputs are unrolled (a name
//! that is not a VCF is a list of them), and each input is asserted readable only when its turn
//! comes, so a header refusal of an earlier input wins over a missing later one. The library masks
//! the path a header refusal names; the culprit is the shortest prefix of the inputs that the
//! library refuses.

use gatk_barclay::Parser;
use gatk_tools::main_entry::Thrown;
use gatk_tools::merge_vcfs::{merge, MergeError};

use super::{
    absolute, assert_readable, thrown, unroll_paths, vcf_text, with_contigs,
    write_vcf_with_dictionary,
};
use crate::runners::{argument, arguments, dictionary_from_any, flag};

pub(super) fn run(parser: &Parser, _stderr: &mut String) -> Result<i32, Thrown> {
    let output = argument(parser, "OUTPUT").unwrap_or_default();
    // `SamReaderFactory.open(SEQUENCE_DICTIONARY)`, with no readability check in front of it.
    let mut dictionary = match argument(parser, "SEQUENCE_DICTIONARY") {
        None => None,
        Some(path) if !std::path::Path::new(&path).exists() => {
            return Err(thrown(
                "htsjdk.samtools.util.RuntimeIOException",
                format!("java.io.FileNotFoundException: {path} (No such file or directory)"),
            ))
        }
        Some(path) => Some(dictionary_from_any(&path)?),
    };
    let paths = unroll_paths(&arguments(parser, "INPUT"))?;

    let mut texts = Vec::new();
    for path in &paths {
        if let Err(unreadable) = assert_readable(path) {
            // The inputs before this one were opened, and checked, first.
            return Err(header_refusal(&texts, &paths).unwrap_or(unreadable));
        }
        let text = vcf_text(path)?;
        let has_contigs = text
            .lines()
            .take_while(|line| line.starts_with("##"))
            .any(|line| line.starts_with("##contig="));
        // An input with no contig lines takes the dictionary given, or else the first input's:
        // `if (sequenceDictionary == null) sequenceDictionary = fileHeader.getSequenceDictionary()`.
        let text = match &dictionary {
            Some(sequences) if !has_contigs => with_contigs(&text, sequences),
            _ => text,
        };
        if dictionary.is_none() && has_contigs {
            dictionary = Some(sequences_of(&text));
        }
        texts.push(text);
    }

    let merged = match merge(&texts, &arguments(parser, "COMMENT")) {
        Ok(merged) => merged,
        Err(error) => {
            return Err(header_refusal(&texts, &paths)
                .unwrap_or_else(|| thrown(error.java_class(), error.message())))
        }
    };
    // `setReferenceDictionary(sequenceDictionary)`: SEQUENCE_DICTIONARY when given, which is what
    // the index's DICT properties list even where the merged header's contigs differ.
    write_vcf_with_dictionary(
        parser,
        &output,
        &merged,
        flag(parser, "CREATE_INDEX"),
        dictionary.as_deref(),
    )?;
    Ok(0)
}

/// `VCFHeader.getSequenceDictionary()`: each contig line's ID, length and assembly.
fn sequences_of(text: &str) -> Vec<htsjdk_bam::header::SequenceRecord> {
    text.lines()
        .take_while(|line| line.starts_with("##"))
        .filter_map(|line| line.strip_prefix("##contig=<"))
        .map(|body| {
            let field = |key: &str| {
                body.trim_end_matches('>')
                    .split(',')
                    .find_map(|pair| pair.strip_prefix(key))
                    .map(str::to_string)
            };
            let mut record = htsjdk_bam::header::SequenceRecord::new(
                &field("ID=").unwrap_or_default(),
                field("length=")
                    .and_then(|value| value.parse().ok())
                    .unwrap_or(0),
            );
            if let Some(assembly) = field("assembly=") {
                record.attributes.set("AS", &assembly);
            }
            record
        })
        .collect()
}

/// The first header refusal over the inputs read so far, naming the input that caused it.
fn header_refusal(texts: &[String], paths: &[String]) -> Option<Thrown> {
    (1..=texts.len()).find_map(|end| match merge(&texts[..end], &[]) {
        Err(
            error @ (MergeError::MissingDictionary
            | MergeError::IncompatibleContigs
            | MergeError::DifferentSamples),
        ) => Some(thrown(
            error.java_class(),
            error
                .message()
                .replace("<masked>", &absolute(&paths[end - 1])),
        )),
        _ => None,
    })
}
