//! `RenameSampleInVcf`, ported from `picard.vcf.RenameSampleInVcf` (Picard 3.4.0).
//!
//! A single-sample VCF written back out with its one sample column renamed.
//!
//! # The tool is two checks and a list of one name
//!
//! ```java
//! if (header.getGenotypeSamples().size() > 1) { throw ... }
//! if (OLD_SAMPLE_NAME != null && !OLD_SAMPLE_NAME.equals(header.getGenotypeSamples().get(0))) { throw ... }
//! final VCFHeader outHeader = new VCFHeader(header.getMetaDataInInputOrder(), makeList(NEW_SAMPLE_NAME));
//! ```
//!
//! Everything else in the output is the writer's. The records are not touched at all: they are read
//! and handed straight back, so what changes in them is whatever the encoder spells differently
//! from the input, which is mostly the QUAL column.
//!
//! # The genotypes are carried through untouched, and this port has to re-key them
//!
//! htsjdk leaves a record's genotypes LAZY after reading: the encoder finds unparsed text and
//! appends it verbatim rather than looking each sample up by name.
//!
//! ```java
//! if (gc.isLazyWithData() && ((LazyGenotypesContext) gc).getUnparsedGenotypeData() instanceof String) {
//!     vcfOutput.append(((LazyGenotypesContext) gc).getUnparsedGenotypeData().toString());
//! ```
//!
//! That is the only reason a rename keeps its data: the writer would otherwise ask the record for
//! a sample called by the NEW name and find nothing, writing `./.` on every row. This port decodes
//! its genotypes eagerly, so it renames the genotype's own sample as well as the header's, which
//! is the same bytes out for the same bytes in.
//!
//! # A sites-only VCF is renamed rather than refused
//!
//! The check is `size() > 1`, and a file with NO samples passes it. The output then declares one
//! sample, and every record, having no genotype for it, is written with `GT` and `./.`. That is a
//! different file in a way the tool never mentions.
//!
//! # The new name is not validated
//!
//! A name with a space, a name that is a number, and the name the file already had are all
//! accepted, and go into the column header as they are.

use htsjdk_vcf::encoder::EncodeError;
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::vcf_file::write_vcf;

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RenameError {
    /// More than one sample column.
    NotSingleSample,
    /// `OLD_SAMPLE_NAME` given and not the file's, which names the sample that was there.
    UnexpectedSample { contained: String },
    /// The reader or the writer refused, which includes a record using an undeclared INFO key.
    Vcf(String, String),
}

impl RenameError {
    pub fn java_class(&self) -> &str {
        match self {
            RenameError::NotSingleSample | RenameError::UnexpectedSample { .. } => {
                "java.lang.IllegalArgumentException"
            }
            RenameError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            RenameError::NotSingleSample => "Input VCF must be single-sample.".to_string(),
            RenameError::UnexpectedSample { contained } => {
                format!("Input VCF did not contain expected sample. Contained: {contained}")
            }
            RenameError::Vcf(_, message) => message.clone(),
        }
    }
}

/// `doWork()`: the whole run, text in and text out.
///
/// `old_name` is the optional `OLD_SAMPLE_NAME`, checked against the first sample only.
pub fn rename(input: &str, new_name: &str, old_name: Option<&str>) -> Result<String, RenameError> {
    let file = read_vcf(input).map_err(|failure| {
        RenameError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    if file.header.samples.len() > 1 {
        return Err(RenameError::NotSingleSample);
    }
    if let Some(old) = old_name {
        // `getGenotypeSamples().get(0)` on a sites-only file would throw; no measured run reaches
        // that, and the check is written as the reference's comparison against the first sample.
        let first = file.header.samples.first().cloned().unwrap_or_default();
        if old != first {
            return Err(RenameError::UnexpectedSample { contained: first });
        }
    }
    let mut header = file.header.clone();
    header.samples = vec![new_name.to_string()];
    // The reference does this by not decoding at all; see the note above.
    let mut records = file.records.clone();
    for record in &mut records {
        for genotype in &mut record.genotypes {
            genotype.sample_name = new_name.to_string();
        }
    }
    write_vcf(&header, &records).map_err(|error| match error {
        // The writer's own refusal, which a record using an undeclared INFO key reaches. The
        // reference's message names the key, the field, and the record it was on.
        EncodeError::MissingFromHeader {
            key,
            field,
            contig,
            start,
        } => RenameError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!(
                "Key {key} found in VariantContext field {field} at {contig}:{start} but this key \
                 isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers \
                 by default."
            ),
        ),
        other => RenameError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!("{other:?}"),
        ),
    })
}
