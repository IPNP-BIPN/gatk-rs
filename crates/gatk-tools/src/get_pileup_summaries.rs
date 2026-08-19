//! Ported from `org.broadinstitute.hellbender.tools.walkers.contamination.GetPileupSummaries`
//! (GATK 4.6.2.0).
//!
//! The tool that writes what [`crate::calculate_contamination`] reads. Both ends of that chain are
//! now pinned by goldens, which is worth more than either alone: a table this tool writes is a
//! table that tool can read, and both agree with the reference.
//!
//! # Eleven filters, and one of them is far stricter than anything else here
//!
//! `getDefaultReadFilters` returns eleven, replacing the walker's two. Ten are library filters; the
//! eleventh is `MappingQualityReadFilter` at **50**, which is well above the 10 or 20 the rest of
//! the port meets. A read at 30 is wellformed, mapped, primary and useless to this tool.
//!
//! # Only the first variant at a locus, and only if it is a biallelic SNP in range
//!
//! `featureContext.getValues(variants).get(0)`: the first record, whatever else overlaps. It has to
//! be biallelic AND a SNP AND have an allele frequency strictly inside the bounds, which are
//! exclusive at both ends -- a site at exactly the default 0.01 is not summarised.
//!
//! # The two refusals happen at opposite ends of the run
//!
//! A population VCF whose HEADER has no `AF` is refused in `onTraversalStart`, before a read is
//! touched. A VCF whose header declares `AF` but whose records never carry one is refused in
//! `onTraversalSuccess`, after the whole traversal has run and produced an empty table. Same
//! exception class, and a port that checked both at the start would refuse a run the reference
//! completes.

use gatk_engine::pileup_summary::{self, PileupSummary};
use gatk_readfilter::{self as filters, with_header, Parameterized};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::{Value, VariantContext};

/// `DEFAULT_MIN_POPULATION_AF`.
pub const DEFAULT_MIN_POPULATION_AF: f64 = 0.01;
/// `DEFAULT_MAX_POPULATION_AF`.
pub const DEFAULT_MAX_POPULATION_AF: f64 = 0.2;
/// `DEFAULT_MINIMUM_MAPPING_QUALITY`, which is this tool's own and not the library's.
pub const DEFAULT_MINIMUM_MAPPING_QUALITY: i32 = 50;

/// The eleven filters `getDefaultReadFilters` returns, in its order.
pub const DEFAULT_READ_FILTERS: [&str; 11] = [
    "MappingQualityReadFilter",
    "MappingQualityAvailableReadFilter",
    "MappingQualityNotZeroReadFilter",
    "MappedReadFilter",
    "PrimaryLineReadFilter",
    "NotDuplicateReadFilter",
    "PassesVendorQualityCheckReadFilter",
    "NonZeroReferenceLengthAlignmentReadFilter",
    "MateOnSameContigOrNoMappedMateReadFilter",
    "GoodCigarReadFilter",
    "WellformedReadFilter",
];

/// The conjunction those eleven names make.
pub fn default_read_filter(read: &BamRecord, header: &SamHeader) -> bool {
    Parameterized::MappingQuality {
        min: DEFAULT_MINIMUM_MAPPING_QUALITY,
        max: None,
    }
    .test(read)
        && filters::mapping_quality_available(read)
        && filters::mapping_quality_not_zero(read)
        && filters::mapped(read)
        && filters::primary_line(read)
        && filters::not_duplicate(read)
        && filters::passes_vendor_quality_check(read)
        && filters::non_zero_reference_length_alignment(read)
        && filters::mate_on_same_contig_or_no_mapped_mate(read)
        && filters::good_cigar(read)
        && with_header::wellformed(read, header)
}

/// What the tool refuses, both of them `UserException.BadInput`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SummariesError {
    /// `onTraversalStart`: the header declares no `AF`.
    HeaderWithoutAlleleFrequency,
    /// `onTraversalSuccess`: no record carried one.
    NoRecordWithAlleleFrequency,
}

impl SummariesError {
    /// The exception class, which is the same for both.
    pub fn java_class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    }

    /// The message, with the `Bad input: ` prefix the exception adds.
    pub fn message(&self) -> String {
        match self {
            SummariesError::HeaderWithoutAlleleFrequency => {
                "Bad input: Population vcf does not have an allele frequency (AF) info field in \
                 its header."
                    .to_string()
            }
            SummariesError::NoRecordWithAlleleFrequency => {
                "Bad input: No variants in population vcf had an allele frequency (AF) field."
                    .to_string()
            }
        }
    }
}

/// `getAttributeAsDouble(AF, ...)` for this tool: the first value, or nothing.
fn allele_frequency(variant: &VariantContext) -> Option<f64> {
    variant
        .attributes
        .iter()
        .find(|(key, _)| key == "AF")
        .and_then(|(_, value)| first_double(value))
}

/// The first number of an INFO value, however it arrived.
///
/// `getAttributeAsDouble` reads a list by taking its first element and a string by parsing it, so a
/// record whose `AF` is `0.1,0.2` answers 0.1 rather than refusing.
fn first_double(value: &Value) -> Option<f64> {
    match value {
        Value::Double(number) => Some(*number),
        Value::Int(number) => Some(*number as f64),
        Value::Str(text) => text.split(',').next().and_then(|first| first.parse().ok()),
        Value::List(values) => values.first().and_then(first_double),
        Value::Missing | Value::Bool(_) => None,
    }
}

/// `isBiallelic() && isSNP()`.
fn is_biallelic_snp(variant: &VariantContext) -> bool {
    variant.alleles.len() == 2
        && crate::remove_nearby_indels::variant_type(variant)
            == crate::remove_nearby_indels::VariantType::Snp
}

/// One locus's summary, or nothing when the site is not summarised.
///
/// `base_counts` is `ReadPileup.getBaseCounts()` over the reads that survived the filters.
pub fn summarise(
    variant: &VariantContext,
    base_counts: [i32; 4],
    minimum: f64,
    maximum: f64,
) -> Option<PileupSummary> {
    if !is_biallelic_snp(variant) {
        return None;
    }
    // The bounds are strict at both ends, so a site at exactly either one is excluded.
    let frequency = allele_frequency(variant)?;
    if !(minimum < frequency && frequency < maximum) {
        return None;
    }
    PileupSummary::from_base_counts(
        &variant.contig,
        variant.start as i32,
        frequency,
        variant.alleles[0].base_string().as_bytes()[0],
        variant.alleles[1].base_string().as_bytes()[0],
        base_counts,
    )
}

/// `doWork`: the table, or the refusal.
///
/// `sites` is what the traversal produced: the variants overlapping each locus, in order, with the
/// base counts of the filtered pileup at it. `header_declares_allele_frequency` is the check
/// `onTraversalStart` makes against the VCF header.
pub fn run(
    sites: &[(Vec<VariantContext>, [i32; 4])],
    header_declares_allele_frequency: bool,
    sample: &str,
    minimum: f64,
    maximum: f64,
) -> Result<String, SummariesError> {
    if !header_declares_allele_frequency {
        return Err(SummariesError::HeaderWithoutAlleleFrequency);
    }

    let mut summaries = Vec::new();
    let mut saw_with = false;
    let mut saw_without = false;

    for (variants, base_counts) in sites {
        // `getValues(variants)` empty is a locus with no record, which returns before anything
        // else is looked at.
        let Some(variant) = variants.first() else {
            continue;
        };
        // The two flags are set by `alleleFrequencyInRange`, which is only reached for a biallelic
        // SNP: a triallelic record with no AF sets neither.
        if is_biallelic_snp(variant) {
            match allele_frequency(variant) {
                Some(_) => saw_with = true,
                None => saw_without = true,
            }
        }
        if let Some(summary) = summarise(variant, *base_counts, minimum, maximum) {
            summaries.push(summary);
        }
    }

    // The second check, after the traversal rather than before it.
    if saw_without && !saw_with {
        return Err(SummariesError::NoRecordWithAlleleFrequency);
    }

    Ok(pileup_summary::write_to_file(sample, &summaries))
}
