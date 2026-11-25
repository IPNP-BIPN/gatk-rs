//! Ported from `org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup`
//! (GATK 4.7.0.0; both Java files are unchanged since 4.6.2.0):
//! `BasicSomaticShortMutationValidator` and `BasicValidationResult`.
//!
//! Whether a discovery call is confirmed by a validation pileup, and the table that answer is
//! written as.
//!
//! # Only the first alternate, and only a diploid genotype
//!
//! The validator reads `getAllele(0)` and `getAllele(1)` and nothing else, so a multiallelic call
//! is validated on its first alternate alone. A genotype of any other ploidy is not refused: the
//! type of the variant is computed from the second allele BEFORE the ploidy that was just computed
//! is consulted, so a haploid genotype throws out of the list access. The golden carries that
//! throw, and this port carries it as a refusal rather than pretending it is a `false`.
//!
//! # The four ways the answer is nothing
//!
//! A genotype that cannot be validated, a null validation pileup, and then two more after the work
//! has started: a NaN alternate ratio, and a discovery pileup of no reads at all. All four answer
//! `None`, and the tool writes no row for them.
//!
//! ```java
//! if (!isAbleToValidateGenotype(genotype, referenceAllele) || validationNormalPileup == null){
//!     return null;
//! }
//! ...
//! if (Double.isNaN(maxAltRatioSeenInNormalValidation) || discoveryTumorTotalCount == 0) {
//!     return null;
//! }
//! ```
//!
//! # Reading this file without reading Rust
//!
//! - `Result<Option<T>, ValidationError>` has three outcomes: `Ok(Some(row))` is a result,
//!   `Ok(None)` is the reference's `return null`, and `Err(..)` is a throw. Keeping `null` apart
//!   from a throw is the point of the nesting.
//! - `?` after a `Result` returns its error to the caller at once; `.map_err(..)` first wraps a
//!   lower layer's error in this file's type; `.ok_or(e)` turns a missing value into `Err(e)`.
//! - `&[T]` and `&T` are borrowed, read-only views; the caller keeps ownership. `.clone()` is an
//!   explicit copy, which Java never needs because it shares references.

use crate::read_pileup::ReadPileup;
use crate::somatic_validation_power::{
    calculate_max_alt_ratio, calculate_min_count_for_signal, calculate_num_reads_supporting_allele,
    calculate_power, PowerError,
};
use crate::tsv_table::java_double_to_string;
use crate::variant_context_utils::{is_complex_indel, type_of_variant, Allele, VariantType};

/// `VALIDATABLE_TYPES`: `Arrays.asList(VariantContext.Type.SNP, VariantContext.Type.MNP,
/// VariantContext.Type.INDEL)`. Only `contains` is asked of it, so the order does not matter.
pub const VALIDATABLE_TYPES: [VariantType; 3] =
    [VariantType::Snp, VariantType::Mnp, VariantType::Indel];

/// As much of a `Genotype` as the validator reads.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidationGenotype {
    /// The called alleles, as alleles rather than indices, because that is what the validator
    /// compares against the record's reference.
    pub alleles: Vec<Allele>,
    /// `AD`, whose absence and whose length are both tested.
    pub ad: Option<Vec<i32>>,
    /// `getFilters`, which is null for an unfiltered genotype and is concatenated onto the
    /// record's filters without a separator.
    pub filters: Option<String>,
}

/// What the validator refuses, one variant per throw site, kept as a value so the caller can
/// print the class and message the reference would have.
#[derive(Debug, Clone, PartialEq)]
pub enum ValidationError {
    /// The list access `getAllele(1)` on a genotype of ploidy one.
    AlleleIndexOutOfBounds { index: usize, size: usize },
    /// A count check, or the power calculation underneath.
    Power(PowerError),
    /// `ParamUtils.isPositiveOrZero` on one of the three counts.
    NegativeCount(&'static str),
}

impl ValidationError {
    /// The exception class the reference throws for this refusal.
    pub fn java_class(&self) -> &'static str {
        match self {
            ValidationError::AlleleIndexOutOfBounds { .. } => "java.lang.IndexOutOfBoundsException",
            ValidationError::Power(error) => error.java_class(),
            ValidationError::NegativeCount(_) => "java.lang.IllegalArgumentException",
        }
    }

    /// The exception's message.
    pub fn message(&self) -> String {
        match self {
            ValidationError::AlleleIndexOutOfBounds { index, size } => {
                format!("Index: {index} Size: {size}")
            }
            ValidationError::Power(error) => error.message(),
            ValidationError::NegativeCount(message) => (*message).to_string(),
        }
    }
}

/// `isAbleToValidateGenotype(genotype, referenceAllele)`.
///
/// Five conditions, and the reference computes all five before testing any of them, which is why
/// the type computation reaches the second allele of a genotype that has only one:
///
/// ```java
/// final boolean isDiploid = genotype.getAlleles().size() == 2;
/// final boolean doesGenotypeHaveReference = genotype.getAllele(0).equals(referenceAllele);
/// final boolean isReferenceNotSymbolic = !referenceAllele.isSymbolic();
/// final VariantContext.Type variantType =
///     GATKVariantContextUtils.typeOfVariant(genotype.getAllele(0), genotype.getAllele(1));
/// ```
///
/// (`BasicSomaticShortMutationValidator.java` lines 38 to 41, the last wrapped at a space.) The
/// port keeps that order, so the refusals come out in the reference's order too: an empty
/// genotype fails at index 0 and a haploid one at index 1. The `logger.info` line the reference
/// writes for an unvalidatable genotype is not reproduced.
pub fn is_able_to_validate_genotype(
    genotype: &ValidationGenotype,
    reference: &Allele,
) -> Result<bool, ValidationError> {
    let is_diploid = genotype.alleles.len() == 2;
    // `.first()` is an `Option`: `None` for an empty list. The `unwrap_or(false)` is never what
    // decides the answer, because an empty list is refused just below, and the reference's own
    // `getAllele(0)` throws at this line instead. Both report index 0.
    let does_genotype_have_reference = genotype
        .alleles
        .first()
        .map(|allele| allele == reference)
        .unwrap_or(false);
    let is_reference_not_symbolic = !reference.is_symbolic();
    let first = genotype
        .alleles
        .first()
        .ok_or(ValidationError::AlleleIndexOutOfBounds {
            index: 0,
            size: genotype.alleles.len(),
        })?;
    let second = genotype
        .alleles
        .get(1)
        .ok_or(ValidationError::AlleleIndexOutOfBounds {
            index: 1,
            size: genotype.alleles.len(),
        })?;
    // `typeOfVariant` can throw for a symbolic or spanning-deletion first allele, and that throw
    // escapes `isAbleToValidateGenotype` in the reference as it escapes here.
    let variant_type = type_of_variant(first, second)
        .map_err(|error| ValidationError::Power(PowerError::Allele(error)))?;
    let is_validatable_variant_type =
        VALIDATABLE_TYPES.contains(&variant_type) && !is_complex_indel(first, second);
    // `genotype.hasAD() && (genotype.getAD().length == 2)`.
    let has_known_coverage = genotype
        .ad
        .as_ref()
        .map(|depths| depths.len() == 2)
        .unwrap_or(false);
    Ok(is_diploid
        && does_genotype_have_reference
        && is_validatable_variant_type
        && has_known_coverage
        && is_reference_not_symbolic)
}

/// `BasicValidationResult`: one validated call, and the row it is written as. The fields are the
/// reference's getters, one to one.
#[derive(Debug, Clone, PartialEq)]
pub struct BasicValidationResult {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub reference: Allele,
    pub alternate: Allele,
    pub minimum_validation_read_count: i32,
    pub is_enough_validation_reads: bool,
    pub is_out_of_noise_floor: bool,
    pub power: f64,
    pub validation_alt_count: i32,
    pub validation_ref_count: i32,
    pub discovery_alt_count: i32,
    pub discovery_ref_count: i32,
    pub filters: String,
    pub num_alt_supporting_reads_in_normal: i64,
}

/// `BasicValidationResultTableColumn.COLUMNS`, in the order the writer composes them: each
/// enum constant's `toString`, which is the column name it was built with.
pub const COLUMNS: [&str; 15] = [
    "CONTIG",
    "START",
    "END",
    "ref_allele",
    "alt_allele",
    "t_alt_count",
    "t_ref_count",
    "tv_alt_count",
    "tv_ref_count",
    "min_val_count",
    "power",
    "validated",
    "sufficient_tv_alt_coverage",
    "discovery_vcf_filter",
    "num_alt_reads_in_validation_normal",
];

impl BasicValidationResult {
    /// The record's own line: the discovery counts first, then the validation ones.
    ///
    /// `validated` is `isOutOfNoiseFloor` and `sufficient_tv_alt_coverage` is
    /// `isEnoughValidationReads`, which is the opposite pairing to the one the column names
    /// suggest reading left to right.
    ///
    /// The order is `composeLine`'s, which is also `COLUMNS`'; it is the reverse of the
    /// constructor, which takes `isEnoughValidationReads` before `isOutOfNoiseFloor`.
    ///
    /// Each value is the text `DataLine.set` makes of it: `Integer.toString`, `Long.toString`,
    /// `Boolean.toString`, and for the power `Double.toString`, because `set(int, double)` ends in
    /// `return set(index, Double.toString(value));` whichever branch it took before that
    /// (`DataLine.java` line 287), so a whole-number power is still written `1.0`. The alleles are
    /// `getBaseString()`, their bases as text.
    pub fn line(&self) -> String {
        [
            self.contig.clone(),
            self.start.to_string(),
            self.end.to_string(),
            String::from_utf8_lossy(&self.reference.bases).into_owned(),
            String::from_utf8_lossy(&self.alternate.bases).into_owned(),
            self.discovery_alt_count.to_string(),
            self.discovery_ref_count.to_string(),
            self.validation_alt_count.to_string(),
            self.validation_ref_count.to_string(),
            self.minimum_validation_read_count.to_string(),
            java_double_to_string(self.power),
            self.is_out_of_noise_floor.to_string(),
            self.is_enough_validation_reads.to_string(),
            self.filters.clone(),
            self.num_alt_supporting_reads_in_normal.to_string(),
        ]
        .join("\t")
    }
}

/// `BasicValidationResult.write`: the header, then one line per record, each ended by a newline
/// and separated by tabs.
pub fn write_table(records: &[BasicValidationResult]) -> String {
    let mut text = COLUMNS.join("\t");
    text.push('\n');
    for record in records {
        text.push_str(&record.line());
        text.push('\n');
    }
    text
}

/// `calculateBasicValidationResult`.
///
/// `validation_normal_pileup` is `None` for the null the reference tests for, which is a `None`
/// answer and not a refusal.
///
/// The reference's interval becomes `contig`, `start` and `end`, which are copied into the result
/// unchanged.
#[allow(clippy::too_many_arguments, reason = "the reference's own signature")]
pub fn calculate_basic_validation_result(
    genotype: &ValidationGenotype,
    reference: &Allele,
    validation_normal_pileup: Option<&ReadPileup<'_>>,
    validation_tumor_alt_count: i32,
    validation_tumor_total_count: i32,
    minimum_base_quality: i32,
    contig: &str,
    start: i32,
    end: i32,
    filters: &str,
) -> Result<Option<BasicValidationResult>, ValidationError> {
    if !is_able_to_validate_genotype(genotype, reference)? {
        return Ok(None);
    }
    // The second half of the reference's first `if`: a null pileup is `return null`, tested
    // after the genotype (so a genotype refusal is reported even when the pileup is also null).
    let Some(pileup) = validation_normal_pileup else {
        return Ok(None);
    };
    // The three `ParamUtils.isPositiveOrZero` calls, in the reference's order and with its text.
    if validation_tumor_alt_count < 0 {
        return Err(ValidationError::NegativeCount(
            "Validation alt count must be >= 0",
        ));
    }
    if validation_tumor_total_count < 0 {
        return Err(ValidationError::NegativeCount(
            "Validation total count must be >= 0",
        ));
    }
    if minimum_base_quality < 0 {
        return Err(ValidationError::NegativeCount(
            "Minimum base quality cutoff must be >= 0",
        ));
    }
    let max_alt_ratio = calculate_max_alt_ratio(pileup, reference, minimum_base_quality)
        .map_err(ValidationError::Power)?;
    // `genotype.getAD()`. `expect` would panic on `None`, which cannot happen: a validatable
    // genotype has an AD of exactly two entries.
    let depths = genotype.ad.as_ref().expect("the genotype was validatable");
    let discovery_alt_count = depths[1];
    let discovery_total_count = depths[0] + discovery_alt_count;
    // The second pair of `null` answers, tested after the ratio and the discovery counts have
    // been computed, as in the reference.
    if max_alt_ratio.is_nan() || discovery_total_count == 0 {
        return Ok(None);
    }
    let minimum_count = calculate_min_count_for_signal(validation_tumor_total_count, max_alt_ratio)
        .map_err(ValidationError::Power)?;
    let power = calculate_power(
        validation_tumor_total_count,
        discovery_alt_count,
        discovery_total_count,
        minimum_count,
    )
    .map_err(ValidationError::Power)?;
    // `genotype.getFilters() == null ? "" : genotype.getFilters()`: `unwrap_or_default` gives
    // the empty string for `None`.
    let genotype_filters = genotype.filters.clone().unwrap_or_default();
    let supporting = calculate_num_reads_supporting_allele(
        pileup,
        &genotype.alleles[0],
        &genotype.alleles[1],
        minimum_base_quality,
    )
    .map_err(ValidationError::Power)?;
    Ok(Some(BasicValidationResult {
        contig: contig.to_string(),
        start,
        end,
        reference: reference.clone(),
        alternate: genotype.alleles[1].clone(),
        minimum_validation_read_count: minimum_count,
        // `validationTumorAltCount >= 2`, a literal two in the reference and not
        // `MINIMUM_NUM_READS_FOR_SIGNAL_COUNT`, though the two are equal.
        is_enough_validation_reads: validation_tumor_alt_count >= 2,
        is_out_of_noise_floor: validation_tumor_alt_count >= minimum_count,
        power,
        validation_alt_count: validation_tumor_alt_count,
        validation_ref_count: validation_tumor_total_count - validation_tumor_alt_count,
        discovery_alt_count,
        discovery_ref_count: discovery_total_count - discovery_alt_count,
        // `filters + genotypeFilters`: concatenated with no separator.
        filters: format!("{filters}{genotype_filters}"),
        num_alt_supporting_reads_in_normal: supporting,
    }))
}
