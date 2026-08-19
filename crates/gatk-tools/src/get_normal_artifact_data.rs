//! Ported from `org.broadinstitute.hellbender.tools.walkers.mutect.GetNormalArtifactData`
//! (GATK 4.6.2.0).
//!
//! The training data for the normal-artifact filter this port already has
//! ([`gatk_engine::normal_artifact_filter`]). A locus walker over two samples, two rejection rules,
//! a binomial p-value and a random downsample.
//!
//! # The downsample is random, seeded, and shared
//!
//! `Utils.getRandomGenerator()` is one `java.util.Random(47382911)` for the whole process, and this
//! tool draws from it once per candidate locus. So which loci survive depends on HOW MANY DRAWS
//! CAME BEFORE, and the draw happens BEFORE the last rejection rule -- a locus rejected for having
//! too much tumour support has already consumed a number. A port that drew later, or from its own
//! generator, keeps a different set of loci from the same data.
//!
//! # The allele is chosen in the normal and counted in the tumour
//!
//! `maxElementIndex` over the normal's six counts picks the allele; the tumour's count of THAT
//! allele is what goes in the table, whatever the tumour's own most common base is.
//!
//! # Six slots, not four
//!
//! Indices 0 to 3 are `ACGT`, index 4 is "this element is before an insertion" and index 5 is
//! "before a deletion start". An index at or above 4 makes the record an `INDEL`, which is how an
//! indel becomes an allele without ever being a base.

use gatk_engine::read_pileup::ReadPileup;
use gatk_readfilter::{self as filters, with_header, Parameterized};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// `Mutect2Engine.READ_QUALITY_FILTER_THRESHOLD`.
pub const READ_QUALITY_FILTER_THRESHOLD: i32 = 20;
/// `Mutect2Engine.MIN_READ_LENGTH`, which is why a short read never reaches `apply`.
pub const MIN_READ_LENGTH: i32 = 30;
/// `--error-prob`'s default.
pub const DEFAULT_ERROR_PROBABILITY: f64 = 0.001;
/// The floor `Math.max(1 - tumorPValue, 0.05)` puts under the keep probability.
pub const MINIMUM_DOWNSAMPLE_PROBABILITY: f64 = 0.05;

/// `Mutect2Engine.makeStandardMutect2ReadFilters`, in its order.
pub const STANDARD_MUTECT2_READ_FILTERS: [&str; 12] = [
    "MappingQualityReadFilter",
    "MappingQualityAvailableReadFilter",
    "MappingQualityNotZeroReadFilter",
    "MappedReadFilter",
    "NotSecondaryAlignmentReadFilter",
    "NotDuplicateReadFilter",
    "PassesVendorQualityCheckReadFilter",
    "NonChimericOriginalAlignmentReadFilter",
    "NonZeroReferenceLengthAlignmentReadFilter",
    "ReadLengthReadFilter",
    "GoodCigarReadFilter",
    "WellformedReadFilter",
];

/// The conjunction those twelve names make.
pub fn standard_mutect2_read_filter(read: &BamRecord, header: &SamHeader) -> bool {
    Parameterized::MappingQuality {
        min: READ_QUALITY_FILTER_THRESHOLD,
        max: None,
    }
    .test(read)
        && filters::mapping_quality_available(read)
        && filters::mapping_quality_not_zero(read)
        && filters::mapped(read)
        && filters::not_secondary_alignment(read)
        && filters::not_duplicate(read)
        && filters::passes_vendor_quality_check(read)
        && filters::non_chimeric_original_alignment(read)
        && filters::non_zero_reference_length_alignment(read)
        && Parameterized::ReadLength {
            min: MIN_READ_LENGTH,
            max: i32::MAX,
        }
        .test(read)
        && filters::good_cigar(read)
        && with_header::wellformed(read, header)
}

/// One row of the table.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalArtifactRecord {
    /// The normal's count of the chosen allele.
    pub normal_alt_count: i32,
    /// The normal pileup's size.
    pub normal_depth: i32,
    /// The tumour's count of the SAME allele.
    pub tumor_alt_count: i32,
    /// The tumour pileup's size.
    pub tumor_depth: i32,
    /// The probability this locus was kept with.
    pub downsampling: f64,
    /// `SNV` or `INDEL`, decided by which of the six slots won.
    pub kind: &'static str,
}

/// `getBaseCounts`: four bases and two indel slots.
///
/// A deletion is skipped; an element before an insertion or before a deletion start is counted in
/// slot 4 or 5 and NOT as its base, so an indel-flanking base never contributes to `ACGT`. Only a
/// base that differs from the reference is counted at all, which is what makes these counts
/// alternate counts rather than depths.
pub fn base_counts(pileup: &ReadPileup, reference_base: u8) -> [f64; 6] {
    let mut counts = [0.0; 6];
    for element in &pileup.elements {
        if element.is_deletion() {
            continue;
        } else if element.is_before_insertion() {
            counts[4] += 1.0;
        } else if element.is_before_deletion_start() {
            counts[5] += 1.0;
        } else if element.base() != reference_base {
            let index = gatk_engine::base_utils::simple_base_to_base_index(element.base());
            if index != -1 {
                counts[index as usize] += 1.0;
            }
        }
    }
    counts
}

/// What one locus produced, and what it consumed.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// No record, and no draw from the generator.
    SkippedBeforeDraw,
    /// No record, but a number was drawn.
    SkippedAfterDraw,
    /// A record, and a number was drawn.
    Kept(Box<NormalArtifactRecord>),
}

/// `apply` for one locus.
///
/// `draw` is the next `nextDouble()` from the shared generator, and the caller decides when to
/// take it -- which is the point: this function says exactly when the reference would have.
pub fn apply(
    normal: &ReadPileup,
    tumor: &ReadPileup,
    reference_base: u8,
    error_probability: f64,
    draw: impl FnOnce() -> f64,
) -> Outcome {
    let normal_counts = base_counts(normal, reference_base);
    let best_allele = gatk_engine::math_utils::max_element_index(&normal_counts, 0, 6);
    let normal_alt_count = normal_counts[best_allele] as i32;
    // No alternate at all, or so much of one that this is a variant rather than an artefact.
    if normal_alt_count == 0 || f64::from(normal_alt_count) > 0.2 * normal.size() as f64 {
        return Outcome::SkippedBeforeDraw;
    }

    let tumor_alt_count = base_counts(tumor, reference_base)[best_allele] as i32;
    // `1 - BinomialDistribution(tumorDepth, errorProb).cumulativeProbability(tumorAlt - 1)`, which
    // for a count of zero is `1 - 0` because the distribution answers zero below its support.
    let cumulative = jmath::binomial::cumulative_probability(
        tumor.size() as i32,
        error_probability,
        tumor_alt_count - 1,
    )
    .expect("the binomial parameters are valid");
    let tumor_p_value = 1.0 - cumulative;
    let downsampling = (1.0 - tumor_p_value).max(MINIMUM_DOWNSAMPLE_PROBABILITY);

    // The draw happens here, BEFORE the last rejection rule, so a locus rejected below has still
    // consumed a number.
    if draw() > downsampling {
        return Outcome::SkippedAfterDraw;
    }
    if f64::from(tumor_alt_count) > 0.5 * tumor.size() as f64 {
        return Outcome::SkippedAfterDraw;
    }

    Outcome::Kept(Box::new(NormalArtifactRecord {
        normal_alt_count,
        normal_depth: normal.size() as i32,
        tumor_alt_count,
        tumor_depth: tumor.size() as i32,
        downsampling,
        kind: if best_allele < 4 { "SNV" } else { "INDEL" },
    }))
}

/// `NormalArtifactRecord.writeToFile`'s column names.
pub const COLUMNS: [&str; 6] = [
    "normal_alt",
    "normal_dp",
    "tumor_alt",
    "tumor_dp",
    "downsampling",
    "type",
];

/// The whole table: the header, then one line per record.
pub fn write(records: &[NormalArtifactRecord]) -> String {
    let mut text = COLUMNS.join("\t");
    text.push('\n');
    for record in records {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            record.normal_alt_count,
            record.normal_depth,
            record.tumor_alt_count,
            record.tumor_depth,
            gatk_engine::tsv_table::java_double_to_string(record.downsampling),
            record.kind
        ));
    }
    text
}
