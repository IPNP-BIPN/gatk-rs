//! `MethylationTypeCaller`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.MethylationTypeCaller` (GATK 4.6.2.0).
//!
//! A `LocusWalker` over bisulfite-sequenced reads, and the first tool here whose output is a **VCF**
//! rather than a BAM. Bisulfite treatment turns an unmethylated C into a T, so at every reference C
//! the tool counts the reads that stayed C against those that became T; on the other strand the same
//! event reads as a G that became an A.
//!
//! # The strand is chosen by the reference base, not by the reads
//!
//! ```java
//! if (referenceBase == (byte)'C') {
//!     final ReadPileup forwardBasePileup = alignmentContext.stratify(FORWARD).getBasePileup();
//! } else if (referenceBase == (byte)'G') {
//!     final ReadPileup reverseBasePileup = alignmentContext.stratify(REVERSE).getBasePileup();
//! } else { return; }
//! ```
//!
//! A C is counted from forward reads alone and a G from reverse reads alone. A C covered only by
//! reverse reads therefore produces **no record**, while those same reads still count towards `DP`
//! at a neighbouring site. `DP` is the whole pileup, both strands and every base, where the two
//! counts are one strand and two bases: the three numbers of a record are not meant to add up.
//!
//! # The two contexts are read on their own strands
//!
//! The forward context is `getBases(0, 2)`, three bases starting at the site. The reverse context is
//! `getBases(2, 0)` **reverse complemented**, three bases ending at the site read backwards, so both
//! begin with the C of the site as that strand sees it. Near a contig edge the window is trimmed
//! rather than padded, so a site two bases from the end gets a two-base context and the last base of
//! the contig gets one.

use gatk_engine::base_utils;
use gatk_engine::context::ReferenceContext;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::locus_iterator::AlignmentContext;
use gatk_engine::read;
use gatk_engine::reads::ReadsDataSource;
use gatk_engine::reference::{ReferenceError, ReferenceFileSource};
use htsjdk_bam::header::SamHeader;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Value, VariantContext};

/// `GATKVCFConstants.UNCONVERTED_BASE_COVERAGE_KEY`.
pub const UNCONVERTED_BASE_COVERAGE_KEY: &str = "UNCONVERTED_BASE_COV";
/// `GATKVCFConstants.CONVERTED_BASE_COVERAGE_KEY`.
pub const CONVERTED_BASE_COVERAGE_KEY: &str = "CONVERTED_BASE_COV";
/// `GATKVCFConstants.METHYLATION_REFERENCE_CONTEXT_KEY`.
pub const METHYLATION_REFERENCE_CONTEXT_KEY: &str = "REFERENCE_CONTEXT";
/// `VCFConstants.DEPTH_KEY`.
pub const DEPTH_KEY: &str = "DP";
/// `VCFConstants.GENOTYPE_KEY`.
pub const GENOTYPE_KEY: &str = "GT";

/// `REFERENCE_CONTEXT_LENGTH`.
pub const REFERENCE_CONTEXT_LENGTH: i32 = 2;

/// What this tool can refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MethylationError {
    /// `createMethylationHeader`: the reads header was null.
    NoReadsHeader,
    /// The reference could not be read.
    Reference(String),
}

impl MethylationError {
    /// The message the reference carries.
    pub fn message(&self) -> String {
        match self {
            MethylationError::NoReadsHeader => {
                "Error writing header, getHeaderForReads() returns null".to_string()
            }
            MethylationError::Reference(text) => text.clone(),
        }
    }
}

/// `createMethylationHeader(header, headerLines)`.
///
/// The samples are the read groups' `SM` values **sorted then deduplicated**, in that order: sorting
/// first is what makes the deduplication remove neighbours rather than every later repeat.
///
/// `default_lines` is what `getDefaultToolVCFHeaderLines` returned. It is a parameter because it
/// carries the run's own date, which no golden can hold: a caller that wants a comparable file
/// passes nothing, which is what `--add-output-vcf-command-line false` does.
pub fn create_methylation_header(header: &SamHeader, default_lines: Vec<HeaderLine>) -> VcfHeader {
    let mut lines = default_lines;
    lines.push(HeaderLine::info(
        UNCONVERTED_BASE_COVERAGE_KEY,
        Cardinality::Fixed(1),
        LineType::Integer,
        // The trailing space is the reference's own.
        "Count of reads supporting methylation that are unconverted ",
    ));
    lines.push(HeaderLine::info(
        CONVERTED_BASE_COVERAGE_KEY,
        Cardinality::Fixed(1),
        LineType::Integer,
        "Count of reads supporting methylation that are converted ",
    ));
    lines.push(HeaderLine::info(
        METHYLATION_REFERENCE_CONTEXT_KEY,
        Cardinality::Fixed(1),
        LineType::String,
        "Forward Strand Reference context",
    ));
    // `VCFStandardHeaderLines.getInfoLine(DEPTH_KEY)` and `getFormatLine(GENOTYPE_KEY)`.
    lines.push(HeaderLine::info(
        DEPTH_KEY,
        Cardinality::Fixed(1),
        LineType::Integer,
        "Approximate read depth; some reads may have been filtered",
    ));
    lines.push(HeaderLine::format(
        GENOTYPE_KEY,
        Cardinality::Fixed(1),
        LineType::String,
        "Genotype",
    ));

    let mut samples: Vec<String> = header
        .read_groups
        .iter()
        .filter_map(|group| group.attributes.get("SM").map(|s| s.to_string()))
        .collect();
    samples.sort();
    samples.dedup();

    VcfHeader { lines, samples }
}

/// `MethylationTypeCaller`: the whole tool, from the reads to the VCF text.
///
/// `default_lines` is what `getDefaultToolVCFHeaderLines` would have returned. It stays a parameter
/// because the default set carries the run's own date, which is why the comparable runs of the
/// suite pass `--add-output-vcf-command-line false` and hand an empty list here.
pub fn methylation_type_caller(
    source: &ReadsDataSource,
    reference: &mut ReferenceFileSource,
    intervals: Option<&[SimpleInterval]>,
    default_lines: Vec<HeaderLine>,
) -> Result<String, MethylationError> {
    let header = source.header().clone();
    let reads = crate::read_walker::traverse(source, intervals.unwrap_or(&[]), &|_| true)
        .map_err(|error| MethylationError::Reference(format!("{error:?}")))?;

    let filter = crate::locus_walker::default_filter(&header);
    let applied = crate::locus_walker::traverse(
        &reads,
        &header,
        Some(reference),
        intervals,
        crate::locus_walker::Options::default(),
        &filter,
    )
    .map_err(|error| MethylationError::Reference(format!("{error:?}")))?;

    let mut records = Vec::new();
    for one in &applied {
        if let Some(record) = apply(&one.context, reference, &one.reference)? {
            records.push(record);
        }
    }

    let vcf_header = create_methylation_header(&header, default_lines);
    htsjdk_vcf::vcf_file::write_vcf(&vcf_header, &records)
        .map_err(|error| MethylationError::Reference(format!("{error:?}")))
}

/// `apply`: the record one locus produces, or nothing.
///
/// Returns `None` at a reference base that is neither C nor G, and at a site whose methylated
/// coverage is zero. Both are the reference's own early exits, and the second is why a site covered
/// only by reads showing some third base writes nothing at all.
pub fn apply(
    context: &AlignmentContext<'_>,
    reference: &mut ReferenceFileSource,
    reference_context: &ReferenceContext,
) -> Result<Option<VariantContext>, MethylationError> {
    let bases = reference_context
        .bases_of(reference, &site(context))
        .map_err(|error| MethylationError::Reference(format!("{error:?}")))?;
    let Some(&reference_base) = bases.first() else {
        return Ok(None);
    };

    let (alt, unconverted, converted, context_bases) = match reference_base {
        b'C' => {
            let forward = context
                .pileup
                .filtered(|element| !read::is_reverse_strand(element.read));
            let counts = forward.base_counts();
            let unconverted = counts[base_index(b'C')];
            let converted = counts[base_index(b'T')];
            let bases = if unconverted + converted > 0 {
                // `getBases(0, REFERENCE_CONTEXT_LENGTH)`: the site and the two bases after it.
                Some(window(
                    reference,
                    &site(context),
                    0,
                    REFERENCE_CONTEXT_LENGTH,
                )?)
            } else {
                None
            };
            (b'T', unconverted, converted, bases)
        }
        b'G' => {
            let reverse = context
                .pileup
                .filtered(|element| read::is_reverse_strand(element.read));
            let counts = reverse.base_counts();
            let unconverted = counts[base_index(b'G')];
            let converted = counts[base_index(b'A')];
            let bases = if unconverted + converted > 0 {
                // `getBases(REFERENCE_CONTEXT_LENGTH, 0)` reverse complemented: the two bases
                // before the site and the site, read on the other strand.
                let forward = window(reference, &site(context), REFERENCE_CONTEXT_LENGTH, 0)?;
                Some(base_utils::simple_reverse_complement(&forward))
            } else {
                None
            };
            (b'A', unconverted, converted, bases)
        }
        // Neither C nor G: the reference strand does not support methylation.
        _ => return Ok(None),
    };

    let Some(context_bases) = context_bases else {
        return Ok(None);
    };

    let alleles = vec![
        Allele::create(&[reference_base], true).expect("a single base is an allele"),
        Allele::create(&[alt], false).expect("a single base is an allele"),
    ];
    let mut variant = VariantContext::new(&context.contig, context.position as i64, alleles);
    variant.id = ".".to_string();
    // `vcb.unfiltered()`: filters were never applied, which prints `.` rather than `PASS`.
    variant.filters = None;
    variant.attributes = vec![
        (
            UNCONVERTED_BASE_COVERAGE_KEY.to_string(),
            Value::Int(unconverted as i64),
        ),
        (
            CONVERTED_BASE_COVERAGE_KEY.to_string(),
            Value::Int(converted as i64),
        ),
        (
            METHYLATION_REFERENCE_CONTEXT_KEY.to_string(),
            Value::Str(String::from_utf8_lossy(&context_bases).into_owned()),
        ),
        // `alignmentContext.size()`: the whole pileup, not the strand that was counted.
        (
            DEPTH_KEY.to_string(),
            Value::Int(context.pileup.size() as i64),
        ),
    ];
    // `vcb.noGenotypes()`. The encoder still writes a `./.` column for every sample of the header.
    variant.genotypes = Vec::new();
    Ok(Some(variant))
}

/// The one-base interval the locus is.
fn site(context: &AlignmentContext<'_>) -> SimpleInterval {
    SimpleInterval {
        contig: context.contig.clone(),
        start: context.position,
        end: context.position,
    }
}

/// `ReferenceContext.getBases(leading, trailing)` for a context whose window is the site itself.
///
/// Trimmed to the contig rather than padded, which is what gives a site near the end a shorter
/// context than one in the middle.
fn window(
    reference: &mut ReferenceFileSource,
    site: &SimpleInterval,
    leading: i32,
    trailing: i32,
) -> Result<Vec<u8>, MethylationError> {
    let length = reference
        .sequence_length(&site.contig)
        .ok_or_else(|| MethylationError::Reference(format!("unknown contig {}", site.contig)))?
        as i32;
    let start = (site.start - leading).max(1);
    let end = (site.end + trailing).min(length);
    reference
        .query(&site.contig, start, end)
        .map_err(|error: ReferenceError| MethylationError::Reference(format!("{error:?}")))
}

/// `BaseUtils.simpleBaseToBaseIndex`, as an index into the four counts.
fn base_index(base: u8) -> usize {
    // The four bases this tool counts are all in the table, so the -1 an unknown base returns is
    // unreachable here; it would index out of the counts array rather than be silently ignored.
    let index = base_utils::simple_base_to_base_index(base);
    assert!(index >= 0, "A, C, G and T are the four counted bases");
    index as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::ReadGroup;

    fn header_with(samples: &[&str]) -> SamHeader {
        let mut header = SamHeader::default();
        for (index, sample) in samples.iter().enumerate() {
            let mut group = ReadGroup::new(&format!("rg{index}"));
            group.attributes.set("SM", sample);
            header.read_groups.push(group);
        }
        header
    }

    #[test]
    fn the_samples_are_sorted_and_then_deduplicated() {
        let header = header_with(&["s2", "s1", "s2"]);
        let vcf = create_methylation_header(&header, Vec::new());
        assert_eq!(vcf.samples, vec!["s1".to_string(), "s2".to_string()]);
    }

    #[test]
    fn the_header_carries_five_lines_of_its_own() {
        let header = header_with(&["s1"]);
        let vcf = create_methylation_header(&header, Vec::new());
        assert_eq!(vcf.lines.len(), 5);
        let rendered = vcf.write();
        assert!(rendered.contains("##INFO=<ID=UNCONVERTED_BASE_COV,"));
        assert!(rendered.contains("##FORMAT=<ID=GT,"));
        // The description's trailing space is the reference's own and survives.
        assert!(rendered.contains("that are unconverted \">"));
    }
}
