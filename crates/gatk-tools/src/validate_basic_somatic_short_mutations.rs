//! `ValidateBasicSomaticShortMutations`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.ValidateBasicSomaticShortMutations`
//! (GATK 4.6.2.0).
//!
//! Each discovery call is asked whether a separate tumour-normal pair of bams confirms it. The
//! counting is [`gatk_engine::allele_pileup_counter`] and the arithmetic is
//! [`gatk_engine::somatic_validation_power`]; what is here is the walker's own decisions.
//!
//! # The validation depth is the counter's, not the pileup's
//!
//! ```java
//! final int validationTumorAltCount = validationTumorAllelicCounts.getOrDefault(altAllele, new MutableInt(0)).intValue();
//! final int validationTumorRefCount = validationTumorAllelicCounts.get(referenceAllele).intValue();
//! ```
//!
//! Both counts come out of the map, so the total the power is computed at is the reference plus the
//! first alternate and nothing else: a read carrying a second alternate, or none, lowers the depth
//! rather than raising it. The alternate count has a default and the reference count does not,
//! which is safe only because the counter always holds the reference key.
//!
//! # A null result is skipped and then dereferenced
//!
//! ```java
//! if (basicValidationResult != null) {
//!     results.add(basicValidationResult);
//! }
//! final boolean normalArtifact = basicValidationResult.getNumAltSupportingReadsInNormal() > maxValidationNormalCount;
//! ```
//!
//! The append is guarded and the next line is not. A genotype that is validatable but whose result
//! is null, which is a discovery pileup of no reads at all or a NaN alternate ratio, therefore ends
//! the run with a null pointer. This port carries that as a refusal rather than repairing it: the
//! golden's `zero-ad` row is the tool dying, and a port that answered anything else would disagree
//! with the reference on a real input.
//!
//! # An artifact is powered whatever its power is
//!
//! `normalArtifact || power > minPower`, so a record with alternate reads in the validation normal
//! counts as a false positive even when there was no power to validate it at all. The judgment and
//! the table disagree on purpose for such a record: `validated` in the table is
//! `isOutOfNoiseFloor`, which does not know about the normal, and `JUDGMENT` in the VCF does.

use gatk_engine::allele_pileup_counter::AllelePileupCounter;
use gatk_engine::basic_somatic_short_mutation_validator::{
    calculate_basic_validation_result, is_able_to_validate_genotype, BasicValidationResult,
    ValidationError, ValidationGenotype,
};
use gatk_engine::read_pileup::ReadPileup;
use gatk_engine::variant_context_utils::Allele;
use htsjdk_vcf::variant::format_vcf_double;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK ValidateBasicSomaticShortMutations";

/// `DEFAULT_MIN_BQ_CUTOFF`.
pub const DEFAULT_MIN_BQ_CUTOFF: i32 = 20;
/// `DEFAULT_MIN_POWER`.
pub const DEFAULT_MIN_POWER: f64 = 0.9;
/// `DEFAULT_MAX_VALIDATION_NORMAL_COUNT`.
pub const DEFAULT_MAX_VALIDATION_NORMAL_COUNT: i32 = 1;

/// `POWER_INFO_FIELD_KEY`, `VALIDATION_AD_INFO_FIELD_KEY` and `JUDGMENT_INFO_FIELD_KEY`.
pub const POWER_KEY: &str = "POWER";
pub const VALIDATION_AD_KEY: &str = "VAL_AD";
pub const JUDGMENT_KEY: &str = "JUDGMENT";

/// The three INFO lines the tool adds to the annotated VCF's header.
pub const HEADER_LINES: [&str; 3] = [
    "##INFO=<ID=JUDGMENT,Number=1,Type=String,Description=\"Validation judgment: validated, unvalidated, or skipped.\">",
    "##INFO=<ID=POWER,Number=1,Type=Float,Description=\"Power to validate variant in validation bam.\">",
    "##INFO=<ID=VAL_AD,Number=A,Type=Integer,Description=\"Ref and alt allele count in validation bam.\">",
];

/// `Judgment`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Judgment {
    Validated,
    Unvalidated,
    Skipped,
}

impl Judgment {
    /// The enum constant's name, which is what the attribute is written as.
    pub fn name(&self) -> &'static str {
        match self {
            Judgment::Validated => "VALIDATED",
            Judgment::Unvalidated => "UNVALIDATED",
            Judgment::Skipped => "SKIPPED",
        }
    }
}

/// The arguments the walker reads, with the reference's defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub discovery_sample: String,
    pub validation_case_name: String,
    pub validation_control_name: String,
    pub min_power: f64,
    pub max_validation_normal_count: i32,
    pub min_bq_cutoff: i32,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            discovery_sample: String::new(),
            validation_case_name: String::new(),
            validation_control_name: String::new(),
            min_power: DEFAULT_MIN_POWER,
            max_validation_normal_count: DEFAULT_MAX_VALIDATION_NORMAL_COUNT,
            min_bq_cutoff: DEFAULT_MIN_BQ_CUTOFF,
        }
    }
}

/// What one `apply` produced.
#[derive(Debug, Clone, PartialEq)]
pub struct Applied {
    /// The judgment written to the annotated VCF.
    pub judgment: Judgment,
    /// `POWER`, absent for a skipped record.
    pub power: Option<f64>,
    /// `VAL_AD`, which is the reference count then the alternate count.
    pub validation_ad: Option<(i32, i32)>,
    /// The row the validation table gets, absent for a skipped record.
    pub result: Option<BasicValidationResult>,
    /// Whether the record counted as validated, which is not the table's `validated` column.
    pub validated: bool,
    /// Whether the record counted towards the false positives when it was not validated.
    pub powered: bool,
}

impl Applied {
    /// The INFO field the annotated VCF carries, in the writer's sorted order.
    pub fn info(&self) -> String {
        match (self.power, self.validation_ad) {
            (Some(power), Some((reference, alternate))) => format!(
                "{JUDGMENT_KEY}={};{POWER_KEY}={};{VALIDATION_AD_KEY}={reference},{alternate}",
                self.judgment.name(),
                format_vcf_double(power)
            ),
            _ => format!("{JUDGMENT_KEY}={}", self.judgment.name()),
        }
    }
}

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum ToolError {
    /// The null pointer of a validatable genotype whose result is null.
    NullResult,
    /// The layers underneath refused.
    Validation(ValidationError),
    /// The counter refused, which the constructor's own checks reach.
    Counter(gatk_engine::allele_pileup_counter::CounterError),
}

impl ToolError {
    pub fn java_class(&self) -> &'static str {
        match self {
            ToolError::NullResult => "java.lang.NullPointerException",
            ToolError::Validation(error) => error.java_class(),
            ToolError::Counter(error) => error.java_class(),
        }
    }

    /// The message, the null pointer's helpful text included, which names the expression and the
    /// variable rather than the line.
    pub fn message(&self) -> String {
        match self {
            ToolError::NullResult => "Cannot invoke \"org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BasicValidationResult.getNumAltSupportingReadsInNormal()\" because \"basicValidationResult\" is null".to_string(),
            ToolError::Validation(error) => error.message(),
            ToolError::Counter(error) => error.message(),
        }
    }
}

/// `apply`, for one record.
///
/// `case_pileup` and `control_pileup` are the two samples' pileups at the record's start. A control
/// of `None` is the record being skipped in silence: nothing is written, not even a judgment. A
/// case of `None` is a null pileup handed to the counter, which counts nothing and leaves the
/// validation depth at zero.
#[allow(clippy::too_many_arguments, reason = "the walker's own inputs")]
pub fn apply(
    contig: &str,
    start: i32,
    end: i32,
    reference: &Allele,
    alternates: &[Allele],
    filters: &[String],
    genotype: &ValidationGenotype,
    case_pileup: Option<&ReadPileup<'_>>,
    control_pileup: Option<&ReadPileup<'_>>,
    arguments: &Arguments,
) -> Result<Option<Applied>, ToolError> {
    if reference.is_symbolic() {
        // Unreachable from a VCF, whose reference allele is bases by construction.
        return Ok(None);
    }
    if !is_able_to_validate_genotype(genotype, reference).map_err(ToolError::Validation)? {
        return Ok(Some(Applied {
            judgment: Judgment::Skipped,
            power: None,
            validation_ad: None,
            result: None,
            validated: false,
            powered: false,
        }));
    }
    let Some(control_pileup) = control_pileup else {
        return Ok(None);
    };
    let alternate = &genotype.alleles[1];

    let mut counter = AllelePileupCounter::new(reference, alternates, arguments.min_bq_cutoff)
        .map_err(ToolError::Counter)?;
    if let Some(pileup) = case_pileup {
        counter.add_pileup(pileup).map_err(|error| {
            ToolError::Validation(ValidationError::Power(
                gatk_engine::somatic_validation_power::PowerError::Allele(error),
            ))
        })?;
    }
    let validation_alt_count = counter.count(alternate).unwrap_or(0);
    let validation_ref_count = counter
        .count(reference)
        .expect("the counter always holds the reference key");

    // `getFilters().stream().sorted().collect(joining(";"))`, which is empty for a PASS record.
    let mut sorted = filters.to_vec();
    sorted.sort();
    let filter_string = sorted.join(";");

    let result = calculate_basic_validation_result(
        genotype,
        reference,
        Some(control_pileup),
        validation_alt_count,
        validation_ref_count + validation_alt_count,
        arguments.min_bq_cutoff,
        contig,
        start,
        end,
        &filter_string,
    )
    .map_err(ToolError::Validation)?;

    // The reference reads the result here without a null check, one line after the check it made.
    let Some(result) = result else {
        return Err(ToolError::NullResult);
    };
    let normal_artifact = result.num_alt_supporting_reads_in_normal
        > i64::from(arguments.max_validation_normal_count);
    let validated = !normal_artifact && result.is_out_of_noise_floor;
    let powered = normal_artifact || result.power > arguments.min_power;
    Ok(Some(Applied {
        judgment: if validated {
            Judgment::Validated
        } else {
            Judgment::Unvalidated
        },
        power: Some(result.power),
        validation_ad: Some((validation_ref_count, validation_alt_count)),
        result: Some(result),
        validated,
        powered,
    }))
}

/// `onTraversalSuccess`'s summary, which counts a record by its own type rather than its genotype's.
///
/// A record is a true positive when it validated and a false positive when it did not but was
/// powered; a record that is neither is counted nowhere, and the false negatives are the literal
/// zero the tool passes in.
pub fn count_towards_summary(
    summary: &mut crate::concordance::Summary,
    applied: &Applied,
    is_snp: bool,
) {
    if applied.validated {
        summary.add(
            gatk_engine::concordance_walker::ConcordanceState::TruePositive,
            is_snp,
        );
    } else if applied.powered {
        summary.add(
            gatk_engine::concordance_walker::ConcordanceState::FalsePositive,
            is_snp,
        );
    }
}

/// The annotated VCF's metadata lines: the input's, the tool's three, and `##source`, sorted.
///
/// `##fileformat` stays first, as the writer puts it, and everything else is compared as strings.
/// The lines this is given are the ones htsjdk holds after reading the input, which is not always
/// what the input file said: htsjdk replaces a `FT` declaration with its own reserved one.
pub fn header_lines(input_lines: &[String]) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    for line in input_lines
        .iter()
        .filter(|line| !line.starts_with("##fileformat"))
    {
        if !lines.contains(line) {
            lines.push(line.clone());
        }
    }
    for line in HEADER_LINES {
        let line = line.to_string();
        if !lines.contains(&line) {
            lines.push(line);
        }
    }
    let source = format!("##source={}", TOOL_NAME.trim_start_matches("GATK "));
    if !lines.contains(&source) {
        lines.push(source);
    }
    lines.sort();
    let mut all = vec!["##fileformat=VCFv4.2".to_string()];
    all.extend(lines);
    all
}
