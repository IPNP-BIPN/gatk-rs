//! `SplitVcfs`, ported from `picard.vcf.SplitVcfs` (Picard 3.5.0; `SplitVcfs.java` is unchanged
//! since 3.4.0), with the record typing of
//! `htsjdk.variant.variantcontext.VariantContext.determineType` (htsjdk 5.0.0) under it.
//!
//! One VCF in, two out: the indels in one file and the SNPs in the other.
//!
//! # The type of a record is a pairwise comparison that collapses to MIXED
//!
//! ```java
//! Type biallelicType = typeOfBiallelicVariant(REF, allele);
//! // [...]
//! if (type == null) {
//!     type = biallelicType;
//! }
//! // [...]
//! else if (biallelicType != type) {
//!     return Type.MIXED;
//! }
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
//! instead, and the refusal comes after the earlier records have already been written: an index
//! refusing one of those comes first ([`split_stream`]).
//!
//! # Both files carry the whole input header
//!
//! Samples and all, however few records they end up holding. `CREATE_INDEX` is on by default, so an
//! input whose header declares no contigs is a refusal before either file is opened.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a value or nothing" (`Some(x)` or `None`), Java's nullable reference with the
//!   null check enforced by the compiler; `Result<T, E>` is "a value or an error" (`Ok` or `Err`),
//!   how this codebase models a Java exception;
//! - `?` after a call returns its error to the caller at once, like an exception propagating;
//! - `match` must cover every case of an enum, so adding a type would not compile until each
//!   `match` decides what to do with it;
//! - `&[T]` is a borrowed, read-only view of a list; `&x` lends `x` without copying it, and
//!   `.clone()` is the explicit copy;
//! - a function's last expression is its return value; `return` is only written for early exits.

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::VariantContext;
use htsjdk_vcf::vcf_file::write_vcf;

/// `VariantContext.Type`, as `determineType` produces it. Each Rust variant is the Java constant
/// of the same name ([`RecordType::name`] spells it).
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
///
/// ```java
/// if (allele.isSymbolic()) return Type.SYMBOLIC;
///
/// if (ref.length() == allele.length()) {
///     if (allele.length() == 1) return Type.SNP;
///     else return Type.MNP;
/// }
/// ```
///
/// Anything else is `INDEL`: the reference compares lengths only, with no prefix test, so a
/// same-length pair is a SNP or an MNP however many of its bases differ. The Java's first check,
/// a symbolic reference allele refused with `IllegalStateException`, is not carried here.
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
///
/// One allele is `NO_VARIATION`, as in the Java's `case 1`. The Java refuses zero alleles with an
/// `IllegalStateException`; `<= 1` folds that case into `NO_VARIATION` instead. The first
/// alternate that disagrees with the type seen so far returns `MIXED` at once, so later alternates
/// are never typed.
pub fn record_type(record: &VariantContext) -> RecordType {
    if record.alleles.len() <= 1 {
        return RecordType::NoVariation;
    }
    let reference = &record.alleles[0];
    // `None` until the first alternate is typed: the Java's `type == null`.
    let mut common: Option<RecordType> = None;
    // `[1..]` is every allele after the reference, the Java's `if (allele == REF) continue;`.
    for allele in &record.alleles[1..] {
        let this = biallelic_type(reference, allele);
        // `Some(seen) if seen != this` is a match arm with a guard: it is taken only when the
        // condition after `if` holds, so an equal type falls through to `Some(_)` and is kept.
        match common {
            None => common = Some(this),
            Some(seen) if seen != this => return RecordType::Mixed,
            Some(_) => {}
        }
    }
    // Unreachable fallback: two or more alleles always type at least one alternate.
    common.unwrap_or(RecordType::NoVariation)
}

/// What the tool refuses: each variant is one Java exception, named by [`SplitError::java_class`]
/// and worded by [`SplitError::message`].
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
    /// The fully qualified class of the exception the Java throws for this refusal. The STRICT
    /// refusal is a bare `IllegalStateException` (`SplitVcfs.java` line 114), not a
    /// `PicardException`.
    pub fn java_class(&self) -> &str {
        match self {
            SplitError::NoSequenceDictionary => "picard.PicardException",
            SplitError::UnexpectedType(_) => "java.lang.IllegalStateException",
            SplitError::Vcf(class, _) => class,
        }
    }

    /// The exception's message, copied from `doWork()`; the trailing `\` of each string line joins
    /// it to the next with the leading spaces dropped, so the text is the Java's single literal.
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

/// The two files a run writes, each as the whole VCF text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Split {
    pub snps: String,
    pub indels: String,
    /// The type of each record written, SNP or INDEL, in input order: the order the two writers
    /// were handed them, which is what an index they share sees.
    pub routed: Vec<RecordType>,
    /// The record `STRICT` stopped at, after every record in `routed` was written.
    pub unexpected: Option<RecordType>,
}

/// `doWork()`: the whole run.
///
/// `strict` is the tool's `STRICT`, which defaults to true; `create_index` is `CREATE_INDEX`, which
/// this tool also leaves on:
///
/// ```java
/// public Boolean STRICT = true;
/// ```
///
/// ```java
/// public SplitVcfs() {
///     this.CREATE_INDEX = true;
/// }
/// ```
///
/// A run [`split_stream`] stopped at a `STRICT` record becomes the error here, and the records
/// written before it are dropped with the [`Split`].
pub fn split(input: &str, strict: bool, create_index: bool) -> Result<Split, SplitError> {
    let split = split_stream(input, strict, create_index)?;
    match split.unexpected {
        Some(kind) => Err(SplitError::UnexpectedType(kind)),
        None => Ok(split),
    }
}

/// [`split`] as Picard streams it: a `STRICT` refusal stops the run at its record, and the records
/// before it are written all the same, so what their indexes refuse comes first.
pub fn split_stream(input: &str, strict: bool, create_index: bool) -> Result<Split, SplitError> {
    // `map_err` turns the reader's error into this tool's, and `?` returns it to the caller.
    let file = read_vcf(input).map_err(|failure| {
        SplitError::Vcf(failure.error.class().to_string(), failure.error.message())
    })?;
    let has_dictionary = file
        .header
        .lines
        .iter()
        .any(|line| matches!(line, HeaderLine::Contig { .. }));
    // `fileHeader.getSequenceDictionary()` is null exactly when the header has no contig line;
    // `SEQUENCE_DICTIONARY`, which would supply one instead, is not an input of this function.
    if create_index && !has_dictionary {
        return Err(SplitError::NoSequenceDictionary);
    }

    let mut snps = Vec::new();
    let mut indels = Vec::new();
    let mut routed = Vec::new();
    let mut unexpected = None;
    for record in &file.records {
        // The indel test runs first, though no record can answer both:
        // `if (context.isIndel()) indelWriter.add(context);`
        let kind = record_type(record);
        match kind {
            RecordType::Indel => indels.push(record.clone()),
            RecordType::Snp => snps.push(record.clone()),
            other => {
                // `STRICT` throws at this record; without it the record is only counted (for a
                // debug log line) and skipped, which `continue` does.
                if strict {
                    unexpected = Some(other);
                    break;
                }
                continue;
            }
        }
        routed.push(kind);
    }

    // A closure (Java's lambda) that writes the whole input header and the given records.
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
        routed,
        unexpected,
    })
}
