//! Ported from `org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup`
//! (GATK 4.6.2.0): `BasicSomaticShortMutationValidator` and `BasicValidationResult`.
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
//! # The three ways the answer is nothing
//!
//! A genotype that cannot be validated, a null validation pileup, and then two more after the work
//! has started: a NaN alternate ratio, and a discovery pileup of no reads at all. All four answer
//! `None`, and the tool writes no row for them.

use crate::read_pileup::ReadPileup;
use crate::somatic_validation_power::{
    calculate_max_alt_ratio, calculate_min_count_for_signal, calculate_num_reads_supporting_allele,
    calculate_power, PowerError,
};
use crate::tsv_table::java_double_to_string;
use crate::variant_context_utils::{is_complex_indel, type_of_variant, Allele, VariantType};

/// `VALIDATABLE_TYPES`.
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

/// What the validator refuses.
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
    pub fn java_class(&self) -> &'static str {
        match self {
            ValidationError::AlleleIndexOutOfBounds { .. } => "java.lang.IndexOutOfBoundsException",
            ValidationError::Power(error) => error.java_class(),
            ValidationError::NegativeCount(_) => "java.lang.IllegalArgumentException",
        }
    }

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
/// the type computation reaches the second allele of a genotype that has only one.
pub fn is_able_to_validate_genotype(
    genotype: &ValidationGenotype,
    reference: &Allele,
) -> Result<bool, ValidationError> {
    let is_diploid = genotype.alleles.len() == 2;
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
    let variant_type = type_of_variant(first, second)
        .map_err(|error| ValidationError::Power(PowerError::Allele(error)))?;
    let is_validatable_variant_type =
        VALIDATABLE_TYPES.contains(&variant_type) && !is_complex_indel(first, second);
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

/// One validated call, and the row it is written as.
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

/// `BasicValidationResultTableColumn.COLUMNS`, in the order the writer composes them.
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

/// `BasicValidationResult.write`: the header, then one line per record.
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
    let Some(pileup) = validation_normal_pileup else {
        return Ok(None);
    };
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
    let depths = genotype.ad.as_ref().expect("the genotype was validatable");
    let discovery_alt_count = depths[1];
    let discovery_total_count = depths[0] + discovery_alt_count;
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
        is_enough_validation_reads: validation_tumor_alt_count >= 2,
        is_out_of_noise_floor: validation_tumor_alt_count >= minimum_count,
        power,
        validation_alt_count: validation_tumor_alt_count,
        validation_ref_count: validation_tumor_total_count - validation_tumor_alt_count,
        discovery_alt_count,
        discovery_ref_count: discovery_total_count - discovery_alt_count,
        filters: format!("{filters}{genotype_filters}"),
        num_alt_supporting_reads_in_normal: supporting,
    }))
}
