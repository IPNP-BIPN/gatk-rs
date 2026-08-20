//! `UpdateVcfSequenceDictionary`, ported from `picard.vcf.UpdateVcfSequenceDictionary`
//! (Picard 3.4.0).
//!
//! A VCF whose contig lines are replaced by a dictionary read from elsewhere.
//!
//! # Not the GATK tool of almost the same name
//!
//! [`crate::update_vcf_sequence_dictionary`] is GATK's, and the two agree on very little. GATK's
//! refuses an input that already has a dictionary unless `--replace` is given, refuses a record on
//! a contig the dictionary lacks, and refuses a record that runs past a sequence's end. Picard's
//! does none of that: it replaces the contig lines and writes every record through.
//!
//! # So the output can declare fewer contigs than its records use
//!
//! A contig the input had and the dictionary lacks is gone from the header, and its records are
//! still written. A dictionary with no sequences at all is accepted, leaving a file with no contig
//! lines and records on two contigs.
//!
//! # Only part of a `.dict` survives into a contig line
//!
//! `setSequenceDictionary` builds one `##contig` line per sequence, and the line carries the ID,
//! the length and `assembly` when the sequence had `AS`. A sequence's `M5` and `UR` do not reach
//! the VCF, so a round trip through this tool loses them whichever side they started on.

use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::vcf_file::write_vcf;

/// One sequence of the dictionary, as far as a contig line reads it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub length: i64,
    /// `AS`, which becomes `assembly=`. `M5` and `UR` have nowhere to go and are not modelled.
    pub assembly: Option<String>,
}

/// What the tool refuses, which is only what the reader and the writer refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UpdateError {
    Vcf(String, String),
}

impl UpdateError {
    pub fn java_class(&self) -> &str {
        match self {
            UpdateError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            UpdateError::Vcf(_, message) => message.clone(),
        }
    }
}

/// `SAMSequenceDictionaryExtractor.extractDictionary` over a `.dict`'s `@SQ` lines.
pub fn read_dictionary(text: &str) -> Vec<Sequence> {
    text.lines()
        .filter(|line| line.starts_with("@SQ\t"))
        .map(|line| {
            let field = |tag: &str| {
                line.split('\t')
                    .find_map(|part| part.strip_prefix(&format!("{tag}:")))
                    .map(str::to_string)
            };
            Sequence {
                name: field("SN").unwrap_or_default(),
                length: field("LN")
                    .and_then(|value| value.parse().ok())
                    .unwrap_or(0),
                assembly: field("AS"),
            }
        })
        .collect()
}

/// The `##contig` line one sequence writes.
pub fn contig_line(sequence: &Sequence, index: i32) -> HeaderLine {
    let mut fields = vec![
        ("ID".to_string(), sequence.name.clone()),
        ("length".to_string(), sequence.length.to_string()),
    ];
    if let Some(assembly) = &sequence.assembly {
        fields.push(("assembly".to_string(), assembly.clone()));
    }
    HeaderLine::Contig { index, fields }
}

/// `doWork()`: the whole run, text in and text out.
pub fn update(input: &str, dictionary: &[Sequence]) -> Result<String, UpdateError> {
    let file = read_vcf(input).map_err(|failure| {
        UpdateError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    let mut header = file.header.clone();
    // `setSequenceDictionary` replaces the contig lines and touches nothing else.
    header
        .lines
        .retain(|line| !matches!(line, HeaderLine::Contig { .. }));
    for (index, sequence) in dictionary.iter().enumerate() {
        header.lines.push(contig_line(sequence, index as i32));
    }
    write_vcf(&header, &file.records).map_err(|error| {
        UpdateError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!("{error:?}"),
        )
    })
}
