//! `FixVcfHeader`, ported from `picard.vcf.FixVcfHeader` (Picard 3.4.0).
//!
//! A VCF whose header does not declare everything its records use, written back out with the
//! missing declarations invented, or with a header taken wholesale from another file.
//!
//! # An invented line says as little as it can
//!
//! ```java
//! new VCFInfoHeaderLine(id, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,
//!     "Missing description: this INFO line was added by Picard's FixVCFHeader")
//! ```
//!
//! Always `Number=.` and always `Type=String`, whatever the value in the record looked like. A
//! fixed header therefore describes its keys LESS precisely than the data does, and a downstream
//! reader that trusts the header reads an integer as a string.
//!
//! # The standard FORMAT lines arrive whether or not the file uses them
//!
//! `addStandardFormatLines(headerLines, false, Genotype.PRIMARY_KEYS)` puts GT, AD, DP, GQ, PL and
//! FT into the rebuilt header. A file carrying nothing but `GT` comes out declaring six FORMAT
//! lines, five of which nothing in it uses.
//!
//! # A record limit can leave the header wrong
//!
//! `CHECK_FIRST_N_RECORDS` stops the search, and the writer is built with
//! `ALLOW_MISSING_FIELDS_IN_HEADER` unset, so a key that first appears in a later record is not
//! declared and the write then refuses AT that record. The fixing tool fails on the file it was
//! asked to fix.
//!
//! # A replacement header replaces everything but the samples
//!
//! Nothing of the input's own header survives. `ENFORCE_SAME_SAMPLES` is on by default and compares
//! the two sample lists pairwise, naming the index of the first that differs; with it off the
//! input's samples are kept and the header file's are discarded, so a sites-only header still
//! writes the input's columns.

use htsjdk_vcf::encoder::EncodeError;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::standard_header_lines::{standard_format_line, PRIMARY_KEYS};
use htsjdk_vcf::vcf_file::write_vcf;

/// The descriptions the tool invents, one per kind of line.
pub const FILTER_DESCRIPTION: &str =
    "Missing description: this FILTER line was added by Picard's FixVCFHeader";
pub const INFO_DESCRIPTION: &str =
    "Missing description: this INFO line was added by Picard's FixVCFHeader";
pub const FORMAT_DESCRIPTION: &str =
    "Missing description: this FORMAT line was added by Picard's FixVCFHeader";

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FixError {
    /// `enforceSameSamples`, when the two lists are different lengths.
    DifferentSampleCount,
    /// `enforceSameSamples`, naming the first index that differs.
    SampleMismatch {
        index: usize,
        reader: String,
        input: String,
    },
    /// A key the fixing did not declare, met by the writer at the record that uses it.
    MissingFromHeader {
        key: String,
        field: String,
        contig: String,
        start: i64,
    },
    /// The reader refused.
    Vcf(String, String),
}

impl FixError {
    pub fn java_class(&self) -> &str {
        match self {
            FixError::DifferentSampleCount | FixError::SampleMismatch { .. } => {
                "picard.PicardException"
            }
            FixError::MissingFromHeader { .. } => "java.lang.IllegalStateException",
            FixError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            FixError::DifferentSampleCount => {
                "The input VCF had a different # of samples than the input VCF header.".to_string()
            }
            FixError::SampleMismatch {
                index,
                reader,
                input,
            } => format!("Mismatch in the {index}th sample: '{reader}' != '{input}'"),
            FixError::MissingFromHeader {
                key,
                field,
                contig,
                start,
            } => format!(
                "Key {key} found in VariantContext field {field} at {contig}:{start} but this key \
                 isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers \
                 by default."
            ),
            FixError::Vcf(_, message) => message.clone(),
        }
    }
}

/// The `##FILTER` line the tool invents, which is a filter line and not a compound one: a FILTER
/// carries an ID and a description and nothing else.
fn invented_filter(id: &str) -> HeaderLine {
    HeaderLine::Filter {
        id: id.to_string(),
        description: FILTER_DESCRIPTION.to_string(),
    }
}

fn invented(key: &str, id: &str, description: &str) -> HeaderLine {
    HeaderLine::Compound {
        key: key.to_string(),
        id: id.to_string(),
        number: Cardinality::Unbounded,
        line_type: LineType::String,
        description: description.to_string(),
        extra: Vec::new(),
    }
}

fn has_line(header: &VcfHeader, key: &str, id: &str) -> bool {
    header.lines.iter().any(|line| match line {
        HeaderLine::Compound {
            key: line_key,
            id: line_id,
            ..
        } => line_key == key && line_id == id,
        HeaderLine::Filter { id: line_id, .. } => key == "FILTER" && line_id == id,
        HeaderLine::Structured {
            key: line_key,
            fields,
        } => {
            line_key == key
                && fields
                    .iter()
                    .any(|(name, value)| name == "ID" && value == id)
        }
        _ => false,
    })
}

/// `doWork()` with no `HEADER`: the header rebuilt from the records under it.
///
/// `check_first_n_records` is the tool's `-1` for "every record".
pub fn fix(input: &str, check_first_n_records: i32) -> Result<String, FixError> {
    let file = read_vcf(input).map_err(|failure| {
        FixError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    let existing = file.header.clone();
    let mut added: Vec<HeaderLine> = Vec::new();
    let mut seen: Vec<(String, String)> = Vec::new();

    for (index, record) in file.records.iter().enumerate() {
        if check_first_n_records > 0 && index as i32 >= check_first_n_records {
            break;
        }
        // `getFilters()` is empty for a passing record, so PASS invents nothing.
        for filter in record.filters.iter().flatten() {
            let key = ("FILTER".to_string(), filter.clone());
            if !has_line(&existing, "FILTER", filter) && !seen.contains(&key) {
                added.push(invented_filter(filter));
                seen.push(key);
            }
        }
        for (id, _) in &record.attributes {
            let key = ("INFO".to_string(), id.clone());
            if !has_line(&existing, "INFO", id) && !seen.contains(&key) {
                added.push(invented("INFO", id, INFO_DESCRIPTION));
                seen.push(key);
            }
        }
        for genotype in &record.genotypes {
            // `getExtendedAttributes()`, which is everything but the typed fields.
            for (id, _) in &genotype.extended {
                let key = ("FORMAT".to_string(), id.clone());
                if !has_line(&existing, "FORMAT", id) && !seen.contains(&key) {
                    added.push(invented("FORMAT", id, FORMAT_DESCRIPTION));
                    seen.push(key);
                }
            }
        }
    }

    let mut header = existing.clone();
    // `addStandardFormatLines(headerLines, false, PRIMARY_KEYS)`, which the file need not use.
    for id in PRIMARY_KEYS {
        if !has_line(&header, "FORMAT", id) {
            if let Some(line) = standard_format_line(id) {
                header.lines.push(line);
            }
        }
    }
    header.lines.extend(added);
    write(&header, &file, existing.samples.clone())
}

/// `doWork()` with a `HEADER`: the replacement header, and the samples one side or the other's.
pub fn fix_with_header(
    input: &str,
    replacement: &str,
    enforce_same_samples: bool,
) -> Result<String, FixError> {
    let file = read_vcf(input).map_err(|failure| {
        FixError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    let header_file = read_vcf(replacement).map_err(|failure| {
        FixError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;

    let mut reader_samples = file.header.samples.clone();
    reader_samples.sort();
    let mut input_samples = header_file.header.samples.clone();
    input_samples.sort();

    let samples = if enforce_same_samples {
        if reader_samples.len() != input_samples.len() {
            return Err(FixError::DifferentSampleCount);
        }
        for (index, (reader, given)) in reader_samples.iter().zip(input_samples.iter()).enumerate()
        {
            if reader != given {
                return Err(FixError::SampleMismatch {
                    index,
                    reader: reader.clone(),
                    input: given.clone(),
                });
            }
        }
        input_samples
    } else {
        // `new VCFHeader(inputHeader.getMetaDataInInputOrder(), existingHeader.getSampleNamesInOrder())`.
        reader_samples
    };

    let header = VcfHeader {
        lines: header_file.header.lines.clone(),
        samples: samples.clone(),
    };
    write(&header, &file, samples)
}

fn write(
    header: &VcfHeader,
    file: &htsjdk_vcf::reader::VcfFile,
    samples: Vec<String>,
) -> Result<String, FixError> {
    let mut header = header.clone();
    header.samples = samples;
    write_vcf(&header, &file.records).map_err(|error| match error {
        EncodeError::MissingFromHeader {
            key,
            field,
            contig,
            start,
        } => FixError::MissingFromHeader {
            key,
            field: field.to_string(),
            contig,
            start,
        },
        other => FixError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!("{other:?}"),
        ),
    })
}
