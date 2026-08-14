//! `ValidateVariants`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants` (GATK 4.6.2.0).
//!
//! The tool writes nothing. Its whole output is whether it threw and what it said, so the message
//! is the only thing there is to be identical about.
//!
//! # Excluding one type puts REF back in the set
//!
//! `ALL` is the default and the one value that does not go through the concrete set. Excluding
//! anything at all builds `CONCRETE_TYPES` (`REF`, `IDS`, `ALLELES`, `CHR_COUNTS`) and removes the
//! exclusions from it, so `--validation-type-to-exclude ALLELES` on a run with no reference is a
//! `MissingReference` refusal: the argument that turns one check off turns another on.
//!
//! # What `ALL` actually checks
//!
//! Two things without a reference, four with one and a dbSNP file. The reference-base check needs a
//! reference and the ID check needs the dbSNP IDs, so a plain run tests less than the name says.
//!
//! # `--validate-GVCF` is three checks, not two
//!
//! The `<NON_REF>` allele in every record, the records being ordered, and the file covering the
//! whole reference. The last counts every uncovered locus and names the first gap, so a two-record
//! GVCF over a 1900-base reference is refused for the 1898 loci it does not describe. It runs at
//! the end of the traversal, so the per-record checks fire first.

/// The checks, as `ValidationType` names them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationType {
    Ref,
    Ids,
    Alleles,
    ChrCounts,
}

/// `CONCRETE_TYPES`, in the enum's own order, which is what an exclusion is subtracted from.
pub const CONCRETE_TYPES: [ValidationType; 4] = [
    ValidationType::Ref,
    ValidationType::Ids,
    ValidationType::Alleles,
    ValidationType::ChrCounts,
];

/// The arguments that decide which checks run.
#[derive(Debug, Clone, Default)]
pub struct Arguments {
    pub types_to_exclude: Vec<ValidationType>,
    pub do_not_validate_filtered_records: bool,
    pub warn_on_errors: bool,
    pub validate_gvcf: bool,
    pub has_reference: bool,
    pub has_dbsnp: bool,
}

/// As much of a record as the validator reads.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub contig: String,
    pub start: i32,
    pub reference: String,
    pub alternates: Vec<String>,
    /// Empty for `.` or `PASS`.
    pub filters: Vec<String>,
    /// The `AC` entries as the record declares them, one per alternate.
    pub allele_counts: Vec<i32>,
    /// The `AN` entry, where the record has one.
    pub allele_number: Option<i32>,
    /// One call per sample, as allele indices.
    pub genotypes: Vec<Vec<Option<usize>>>,
    /// The QUAL column, which the GVCF messages print as `%.2f`.
    pub qual: Option<f64>,
    /// The type name the messages carry, which `determineType` decides.
    pub variant_type: String,
    /// Every INFO field, in the SORTED order the message prints them in.
    pub attributes: Vec<(String, String)>,
}

/// What the tool refuses with, each carrying the reference's own words.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValidationError {
    /// `REF` asked for with no reference, before any record is read.
    MissingReference,
    /// One of the strict checks, which names the input file and the type that was asked for.
    Strict { input: String, detail: String },
    /// A GVCF record without `<NON_REF>`.
    MissingNonRef { record: String },
    /// A GVCF record behind one already traversed.
    OutOfOrder { record: String },
    /// The gaps a GVCF leaves, counted at the end of the traversal.
    NotCovering { loci: i64, first_gap: String },
}

impl ValidationError {
    pub fn message(&self) -> String {
        match self {
            ValidationError::MissingReference => {
                "Validation type REF was selected but no reference was provided.  A reference is \
                 specified with the -R command line argument."
                    .to_string()
            }
            ValidationError::Strict { input, detail } => {
                format!("Input {input} fails strict validation of type ALL: {detail}")
            }
            ValidationError::MissingNonRef { record } => format!(
                "In a GVCF all records must contain a <NON_REF> allele. Offending record: {record}"
            ),
            // The reference's own grammar, "must ordered", is kept.
            ValidationError::OutOfOrder { record } => format!(
                "In a GVCF all records must ordered. Record: {record} covers a position previously \
                 traversed."
            ),
            ValidationError::NotCovering { loci, first_gap } => format!(
                "A GVCF must cover the entire region. Found {loci} loci with no VariantContext \
                 covering it. The first uncovered segment is:{first_gap}"
            ),
        }
    }

    pub fn java_class(&self) -> &'static str {
        match self {
            ValidationError::MissingReference => {
                "org.broadinstitute.hellbender.exceptions.UserException$MissingReference"
            }
            ValidationError::Strict { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$FailsStrictValidation"
            }
            ValidationError::NotCovering { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$ValidationFailure"
            }
            _ => "org.broadinstitute.hellbender.exceptions.UserException",
        }
    }
}

/// The checks a run makes, which is `ALL` unless something was excluded.
///
/// Returns the refusal `ALL` cannot make: asking for `REF` with no reference is decided here,
/// before any record is read, and it is reached by excluding anything at all.
pub fn types_to_apply(arguments: &Arguments) -> Result<Vec<ValidationType>, ValidationError> {
    if arguments.types_to_exclude.is_empty() {
        // `ALL` on its own: the reference-base check runs only where there is a reference, and the
        // ID check only where there are IDs, so the set is decided by what the run has rather than
        // by what it asked for.
        let mut types = vec![ValidationType::Alleles, ValidationType::ChrCounts];
        if arguments.has_reference {
            types.insert(0, ValidationType::Ref);
        }
        if arguments.has_dbsnp {
            types.push(ValidationType::Ids);
        }
        return Ok(types);
    }
    let types: Vec<ValidationType> = CONCRETE_TYPES
        .iter()
        .filter(|kind| !arguments.types_to_exclude.contains(kind))
        .copied()
        .collect();
    if types.contains(&ValidationType::Ref) && !arguments.has_reference {
        return Err(ValidationError::MissingReference);
    }
    Ok(types)
}

/// One record against the checks the run makes.
///
/// `reference_base` is what the reference holds at the record's start, where the run has one.
pub fn validate_record(
    record: &Record,
    input: &str,
    types: &[ValidationType],
    reference_base: Option<&str>,
    arguments: &Arguments,
) -> Result<(), ValidationError> {
    if arguments.do_not_validate_filtered_records && !record.filters.is_empty() {
        return Ok(());
    }

    if arguments.validate_gvcf && !record.alternates.iter().any(|alt| alt == "<NON_REF>") {
        return Err(ValidationError::MissingNonRef {
            record: rendered(record),
        });
    }

    let strict = |detail: String| ValidationError::Strict {
        input: input.to_string(),
        detail,
    };

    for kind in types {
        match kind {
            ValidationType::Ref => {
                if let Some(observed) = reference_base {
                    if !observed.eq_ignore_ascii_case(&record.reference) {
                        return Err(strict(format!(
                            "the REF allele is incorrect for the record at position {}:{}, fasta \
                             says {observed} vs. VCF says {}",
                            record.contig, record.start, record.reference
                        )));
                    }
                }
            }
            ValidationType::Alleles => {
                // The check is about the GENOTYPES rather than the ALT column: an alternate no
                // sample calls is the failure. SYMBOLIC alleles are excluded from both sides,
                // because a GVCF's `<NON_REF>` is expected never to be called, so the check is
                // about the plain alternates alone.
                let symbolic =
                    |index: usize| index > 0 && record.alternates[index - 1].starts_with('<');
                let observed: Vec<usize> = record
                    .genotypes
                    .iter()
                    .flatten()
                    .flatten()
                    .copied()
                    .filter(|index| !symbolic(*index))
                    .collect();
                let unused = (1..record.alternates.len() + 1)
                    .filter(|index| !symbolic(*index))
                    .any(|index| !observed.contains(&index));
                if !record.genotypes.is_empty() && unused {
                    return Err(strict(format!(
                        "one or more of the ALT allele(s) for the record at position {}:{} are not \
                         observed at all in the sample genotypes",
                        record.contig, record.start
                    )));
                }
            }
            ValidationType::ChrCounts => {
                if record.genotypes.is_empty() {
                    continue;
                }
                let called: Vec<usize> = record
                    .genotypes
                    .iter()
                    .flatten()
                    .flatten()
                    .copied()
                    .collect();
                for (index, declared) in record.allele_counts.iter().enumerate() {
                    let counted =
                        called.iter().filter(|allele| **allele == index + 1).count() as i32;
                    if *declared != counted {
                        return Err(strict(format!(
                            "the Allele Count (AC) tag is incorrect for the record at position \
                             {}:{}, {declared} vs. {counted}",
                            record.contig, record.start
                        )));
                    }
                }
                if let Some(declared) = record.allele_number {
                    let counted = called.len() as i32;
                    if declared != counted {
                        return Err(strict(format!(
                            "the Allele Number (AN) tag is incorrect for the record at position \
                             {}:{}, {declared} vs. {counted}",
                            record.contig, record.start
                        )));
                    }
                }
            }
            ValidationType::Ids => {}
        }
    }
    Ok(())
}

/// `validateVariantsOrder`: the previous start, reset at a contig change.
#[derive(Debug, Default)]
pub struct OrderCheck {
    previous_contig: Option<String>,
    previous_start: i32,
}

impl OrderCheck {
    pub fn new() -> OrderCheck {
        OrderCheck {
            previous_contig: None,
            previous_start: -1,
        }
    }

    pub fn check(&mut self, record: &Record) -> Result<(), ValidationError> {
        if self.previous_contig.as_deref() != Some(record.contig.as_str()) {
            self.previous_contig = Some(record.contig.clone());
            self.previous_start = -1;
        }
        if self.previous_start > -1 && record.start < self.previous_start {
            return Err(ValidationError::OutOfOrder {
                record: rendered(record),
            });
        }
        self.previous_start = record.start;
        Ok(())
    }
}

/// `toStringWithoutGenotypes`, which the two GVCF messages quote whole.
///
/// `String.format("[VC %s @ %s Q%s of type=%s alleles=%s attr=%s filters=%s", ...)`, and the
/// format string has no closing bracket: the message ends after the filters, mid-structure. The
/// source is `Unknown` for a file the walker opened, the position collapses to one number when the
/// record spans one base, the quality is `%.2f`, the alleles are SORTED with the reference first
/// and the attributes are a sorted `key=value` map.
fn rendered(record: &Record) -> String {
    let mut alleles: Vec<String> = record.alternates.clone();
    alleles.sort();
    let alleles = std::iter::once(format!("{}*", record.reference))
        .chain(alleles)
        .collect::<Vec<_>>()
        .join(", ");
    let attributes = record
        .attributes
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect::<Vec<_>>()
        .join(", ");
    let quality = match record.qual {
        Some(qual) => format!("{qual:.2}"),
        None => ".".to_string(),
    };
    format!(
        "[VC Unknown @ {}:{} Q{quality} of type={} alleles=[{alleles}] attr={{{attributes}}} filters={}",
        record.contig,
        record.start,
        record.variant_type,
        record.filters.join(",")
    )
}
