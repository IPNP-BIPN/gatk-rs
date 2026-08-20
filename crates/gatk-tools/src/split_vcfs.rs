//! `SplitVcfs`, ported from `picard.vcf.SplitVcfs` (Picard 3.4.0), with the record typing of
//! `htsjdk.variant.variantcontext.VariantContext.determineType` under it.
//!
//! One VCF in, two out: the indels in one file and the SNPs in the other.
//!
//! # The type of a record is a pairwise comparison that collapses to MIXED
//!
//! ```java
//! Type biallelicType = typeOfBiallelicVariant(REF, allele);
//! if (type == null) { type = biallelicType; }
//! else if (biallelicType != type) { return Type.MIXED; }
//! ```
//!
//! Each alternate is typed against the reference on its own, and the record's type is that common
//! type or `MIXED`. A record with no alternate at all is `NO_VARIATION`.
//!
//! A SPANNING DELETION IS ONE BASE LONG, so `A -> C,*` types both alternates SNP and the record
//! goes to the SNP file. That is not an obvious reading of a star, and the golden is what says it.
//!
//! # Three of the six types go nowhere
//!
//! The tool writes indels to one file and SNPs to the other; a MIXED, MNP, SYMBOLIC or NO_VARIATION
//! record is counted and dropped. `STRICT` is ON BY DEFAULT, so that count is normally a refusal
//! instead, and the refusal comes after the earlier records have already been written.
//!
//! # Both files carry the whole input header
//!
//! Samples and all, however few records they end up holding. `CREATE_INDEX` is on by default, so an
//! input whose header declares no contigs is a refusal before either file is opened.

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::VariantContext;
use htsjdk_vcf::vcf_file::write_vcf;

/// `VariantContext.Type`, as `determineType` produces it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

impl RecordType {
    /// The enum constant's own name, which the STRICT refusal interpolates.
    pub fn name(&self) -> &'static str {
        match self {
            RecordType::NoVariation => "NO_VARIATION",
            RecordType::Snp => "SNP",
            RecordType::Mnp => "MNP",
            RecordType::Indel => "INDEL",
            RecordType::Symbolic => "SYMBOLIC",
            RecordType::Mixed => "MIXED",
        }
    }
}

/// `typeOfBiallelicVariant`, one alternate against the reference.
fn biallelic_type(reference: &Allele, alternate: &Allele) -> RecordType {
    if alternate.is_symbolic() {
        return RecordType::Symbolic;
    }
    if reference.len() == alternate.len() {
        if alternate.len() == 1 {
            return RecordType::Snp;
        }
        return RecordType::Mnp;
    }
    RecordType::Indel
}

/// `determineType`: the common type of every alternate, or MIXED.
pub fn record_type(record: &VariantContext) -> RecordType {
    if record.alleles.len() <= 1 {
        return RecordType::NoVariation;
    }
    let reference = &record.alleles[0];
    let mut common: Option<RecordType> = None;
    for allele in &record.alleles[1..] {
        let this = biallelic_type(reference, allele);
        match common {
            None => common = Some(this),
            Some(seen) if seen != this => return RecordType::Mixed,
            Some(_) => {}
        }
    }
    common.unwrap_or(RecordType::NoVariation)
}

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplitError {
    /// `CREATE_INDEX` with no dictionary to index against.
    NoSequenceDictionary,
    /// `STRICT` and a record that is neither a SNP nor an indel.
    UnexpectedType(RecordType),
    /// The reader or the writer refused.
    Vcf(String, String),
}

impl SplitError {
    pub fn java_class(&self) -> &str {
        match self {
            SplitError::NoSequenceDictionary => "picard.PicardException",
            SplitError::UnexpectedType(_) => "java.lang.IllegalStateException",
            SplitError::Vcf(class, _) => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            SplitError::NoSequenceDictionary => "A sequence dictionary must be available \
                 (either through the input file or by setting it explicitly) when creating \
                 indexed output."
                .to_string(),
            SplitError::UnexpectedType(kind) => {
                format!("Found a record with type {}", kind.name())
            }
            SplitError::Vcf(_, message) => message.clone(),
        }
    }
}

/// The two files a run writes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Split {
    pub snps: String,
    pub indels: String,
}

/// `doWork()`: the whole run.
///
/// `strict` is the tool's `STRICT`, which defaults to true; `create_index` is `CREATE_INDEX`, which
/// this tool also leaves on.
pub fn split(input: &str, strict: bool, create_index: bool) -> Result<Split, SplitError> {
    let file = read_vcf(input).map_err(|failure| {
        SplitError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    let has_dictionary = file
        .header
        .lines
        .iter()
        .any(|line| matches!(line, HeaderLine::Contig { .. }));
    if create_index && !has_dictionary {
        return Err(SplitError::NoSequenceDictionary);
    }

    let mut snps = Vec::new();
    let mut indels = Vec::new();
    for record in &file.records {
        // The indel test runs first, though no record can answer both.
        match record_type(record) {
            RecordType::Indel => indels.push(record.clone()),
            RecordType::Snp => snps.push(record.clone()),
            other => {
                if strict {
                    return Err(SplitError::UnexpectedType(other));
                }
            }
        }
    }

    let encode = |records: &[VariantContext]| {
        write_vcf(&file.header, records).map_err(|error| {
            SplitError::Vcf(
                "java.lang.IllegalStateException".to_string(),
                format!("{error:?}"),
            )
        })
    };
    Ok(Split {
        snps: encode(&snps)?,
        indels: encode(&indels)?,
    })
}
