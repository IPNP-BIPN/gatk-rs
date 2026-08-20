//! `MakeSitesOnlyVcf`, ported from `picard.vcf.MakeSitesOnlyVcf` (Picard 3.4.0).
//!
//! A VCF with its genotype columns dropped, or kept for a named subset.
//!
//! # The header is built from the names asked for, not the names the file has
//!
//! ```java
//! final VCFHeader header = new VCFHeader(inputVcfHeader.getMetaDataInInputOrder(), SAMPLE);
//! ```
//!
//! `SAMPLE` is a `TreeSet`, so the output's columns are ALPHABETICAL however the user typed them
//! and whatever order the input had. And a name the input never carried is not dropped: it becomes
//! a column of its own, with `./.` on every record, because the record's subset has no genotype for
//! it and the writer fills the gap.
//!
//! # The annotations are not recomputed
//!
//! `subsetToSamplesWithOriginalAnnotations` keeps the record's INFO fields as they were and resets
//! the alleles from the original record. So a one-sample output still carries the whole file's `AC`
//! and `AN`, and an ALT no remaining genotype calls still appears in the column. The tool's name
//! says sites-only; what it does not say is that the numbers describe a file that is no longer
//! there.
//!
//! # The default is no samples at all
//!
//! Which leaves eight columns: no FORMAT column and no sample columns, rather than an empty FORMAT.
//! `CREATE_INDEX` is on by default for this tool, so an input whose header declares no contigs is a
//! refusal rather than an unindexed output.

use htsjdk_vcf::encoder::EncodeError;
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::vcf_file::write_vcf;

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SitesOnlyError {
    /// `CREATE_INDEX` with no dictionary to index against.
    NoSequenceDictionary,
    /// The reader or the writer refused.
    Vcf(String, String),
}

impl SitesOnlyError {
    pub fn java_class(&self) -> &str {
        match self {
            SitesOnlyError::NoSequenceDictionary => "picard.PicardException",
            SitesOnlyError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            SitesOnlyError::NoSequenceDictionary => "A sequence dictionary must be available \
                 (either through the input file or by setting it explicitly) when creating \
                 indexed output."
                .to_string(),
            SitesOnlyError::Vcf(_, message) => message.clone(),
        }
    }
}

/// `doWork()`: the whole run, text in and text out.
///
/// `samples` is the `SAMPLE` argument; empty means sites-only. `create_index` is on by default in
/// this tool's constructor, which is why it is a parameter here rather than an assumption.
pub fn make_sites_only(
    input: &str,
    samples: &[String],
    create_index: bool,
) -> Result<String, SitesOnlyError> {
    let file = read_vcf(input).map_err(|failure| {
        SitesOnlyError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    // `inputVcfHeader.getSequenceDictionary()`, which is the contig lines and nothing else.
    let has_dictionary = file
        .header
        .lines
        .iter()
        .any(|line| matches!(line, htsjdk_vcf::header::HeaderLine::Contig { .. }));
    if create_index && !has_dictionary {
        return Err(SitesOnlyError::NoSequenceDictionary);
    }

    // The `TreeSet`: sorted, and de-duplicated.
    let mut wanted: Vec<String> = samples.to_vec();
    wanted.sort();
    wanted.dedup();

    let mut header = file.header.clone();
    header.samples = wanted.clone();
    let mut records = file.records.clone();
    for record in &mut records {
        // `subsetToSamples`, which keeps the genotypes whose sample was asked for and adds nothing
        // for a sample the record does not have.
        record
            .genotypes
            .retain(|genotype| wanted.contains(&genotype.sample_name));
    }
    write_vcf(&header, &records).map_err(|error| match error {
        EncodeError::MissingFromHeader {
            key,
            field,
            contig,
            start,
        } => SitesOnlyError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!(
                "Key {key} found in VariantContext field {field} at {contig}:{start} but this key \
                 isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers \
                 by default."
            ),
        ),
        other => SitesOnlyError::Vcf(
            "java.lang.IllegalStateException".to_string(),
            format!("{other:?}"),
        ),
    })
}
